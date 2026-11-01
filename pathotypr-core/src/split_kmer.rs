//! Implements a dynamic variant detection engine for FASTQ-based genotyping.
//!
//! This module creates diagnostic k-mers for SNPs and MNVs based on a
//! user-provided marker file, then scans FASTQ reads to count evidence
//! for reference vs. alternate alleles.
//!
//! **Note:** Indels are intentionally skipped in the FASTQ workflow because
//! short reads across repetitive regions (e.g., PE/PPE in M. tuberculosis)
//! produce unreliable k-mer matches and high false positive rates. Indels
//! are only reliably detected via the FASTA assembly classify workflow,
//! which uses full-length alignment.

use crate::errors::CancellationToken;
use anyhow::Result;
use fxhash::FxHashMap;
use log::{debug, info, warn};
use needletail::Sequence;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

/// Default length of the diagnostic k-mers used for identification.
pub const DEFAULT_MARKER_KMER_LEN: usize = 31;

// --- Bloom Filter for fast k-mer rejection ---
// 128KB bloom filter (1M bits) — fits in L2 cache (512KB/core on modern CPUs).
// With ~32k entries and 2 hash functions, false positive rate is ~0.004%.
const BLOOM_BITS: usize = 1_048_576;
const BLOOM_WORDS: usize = BLOOM_BITS / 64;

pub type BloomFilter = Box<[u64; BLOOM_WORDS]>;

fn bloom_new() -> BloomFilter {
    // Box to avoid 32KB on the stack — length is guaranteed by BLOOM_WORDS
    vec![0u64; BLOOM_WORDS]
        .into_boxed_slice()
        .try_into()
        .expect("bloom filter size mismatch (this is a compile-time invariant)")
}

fn bloom_insert(filter: &mut [u64; BLOOM_WORDS], key: u64) {
    let h1 = key.wrapping_mul(0x9E3779B97F4A7C15) as usize % BLOOM_BITS;
    let h2 = key.wrapping_mul(0x517CC1B727220A95) as usize % BLOOM_BITS;
    filter[h1 / 64] |= 1u64 << (h1 % 64);
    filter[h2 / 64] |= 1u64 << (h2 % 64);
}

#[inline(always)]
fn bloom_may_contain(filter: &[u64; BLOOM_WORDS], key: u64) -> bool {
    let h1 = key.wrapping_mul(0x9E3779B97F4A7C15) as usize % BLOOM_BITS;
    let h2 = key.wrapping_mul(0x517CC1B727220A95) as usize % BLOOM_BITS;
    (filter[h1 / 64] & (1u64 << (h1 % 64)) != 0)
        && (filter[h2 / 64] & (1u64 << (h2 % 64)) != 0)
}

#[derive(Debug, Clone)]
pub enum MarkerHitList {
    Single((usize, bool)),
    Multi(Vec<(usize, bool)>),
}

pub type MarkerIndex = FxHashMap<u64, MarkerHitList>;

/// Represents a single marker, capable of describing all variant types.
/// MODIFIED: The `lineages` field is now a vector of strings to hold the hierarchy.
#[derive(Debug, Clone)]
pub struct Marker {
    pub pos: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub lineages: Vec<String>, // e.g., ["L2", "L2.2", "L2.2.1"]
    pub alt_kmer: Vec<u8>,
    pub ref_kmer: Vec<u8>,
    /// Whether the reference genome already carries this marker's ALT k-mer
    /// somewhere else. Reads from that copy count as the alternate allele, so
    /// the counts cannot settle the marker: only an anchor can.
    pub alt_in_reference: bool,
    /// Whether the reference carries this marker's window more than once.
    /// Reads from the copy count as the reference allele and drown the real
    /// locus, so again only an anchor can settle it.
    pub ref_repeated: bool,
    /// 0-based start of the reference window both k-mers were cut from. A read
    /// that matches either k-mer matches it exactly, so it reads every base of
    /// `kmer_start .. kmer_start + k`.
    pub kmer_start: usize,
    pub annotations: Vec<String>,
}

// --- Private Helper Functions ---

fn hash_kmer(kmer: &[u8]) -> Option<u64> {
    let mut val: u64 = 0;
    for &base in kmer {
        val = (val << 2)
            | match base {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => return None,
            };
    }
    Some(val)
}

fn hash_reverse_complement(kmer: &[u8]) -> Option<u64> {
    let mut val: u64 = 0;
    for &base in kmer.iter().rev() {
        val = (val << 2)
            | match base {
                b'A' | b'a' => 3,
                b'C' | b'c' => 2,
                b'G' | b'g' => 1,
                b'T' | b't' => 0,
                _ => return None,
            };
    }
    Some(val)
}

fn read_ref_sequence(path: &str) -> Result<Vec<u8>> {
    let mut reader = needletail::parse_fastx_file(path)
        .map_err(|e| anyhow::anyhow!("Cannot open reference FASTA '{}': {}", path, e))?;

    let record = reader
        .next()
        .ok_or_else(|| anyhow::anyhow!("No records found in reference FASTA '{}'.", path))?
        .map_err(|e| anyhow::anyhow!("Error reading reference FASTA '{}': {}", path, e))?;

    let sequence = record.normalize(true).into_owned();

    if reader.next().is_some() {
        anyhow::bail!(
            "Reference FASTA '{}' contains multiple records; provide a single-record FASTA.",
            path
        );
    }

    Ok(sequence)
}

/// Builds a diagnostic k-mer for a small variant (SNP, MNV, small indel).
fn build_small_variant_kmer(
    pos0: usize,
    ref_allele: &[u8],
    alt_allele: &[u8],
    ref_seq: &[u8],
    kmer_len: usize,
) -> Option<(Vec<u8>, Vec<u8>, usize)> {
    let ref_len_allele = ref_allele.len();
    let alt_len_allele = alt_allele.len();
    let max_allele_len = ref_len_allele.max(alt_len_allele);
    let flank_len = (kmer_len - max_allele_len) / 2;
    let kmer_start_pos = pos0.saturating_sub(flank_len);
    let left_flank_len = pos0 - kmer_start_pos;
    let right_flank_start_ref = pos0 + ref_len_allele;

    // Compute separate right flank lengths so each k-mer is exactly kmer_len.
    // For indels, the shorter allele needs more right flank to compensate.
    let ref_right_needed = kmer_len.saturating_sub(left_flank_len + ref_len_allele);
    let alt_right_needed = kmer_len.saturating_sub(left_flank_len + alt_len_allele);
    let max_right_needed = ref_right_needed.max(alt_right_needed);

    if right_flank_start_ref + max_right_needed > ref_seq.len() {
        return None;
    }

    let mut ref_kmer = Vec::with_capacity(kmer_len);
    ref_kmer.extend_from_slice(&ref_seq[kmer_start_pos..pos0]);
    ref_kmer.extend_from_slice(ref_allele);
    ref_kmer.extend_from_slice(
        &ref_seq[right_flank_start_ref..right_flank_start_ref + ref_right_needed],
    );

    let mut alt_kmer = Vec::with_capacity(kmer_len);
    alt_kmer.extend_from_slice(&ref_seq[kmer_start_pos..pos0]);
    alt_kmer.extend_from_slice(alt_allele);
    alt_kmer.extend_from_slice(
        &ref_seq[right_flank_start_ref..right_flank_start_ref + alt_right_needed],
    );

    Some((ref_kmer, alt_kmer, kmer_start_pos))
}

// --- Public API ---

pub fn build_markers(
    ref_fasta: &str,
    tsv_markers: &str,
    cancel_token: &Option<CancellationToken>,
    kmer_len: usize,
) -> Result<Vec<Marker>> {
    if let Some(token) = cancel_token {
        if token.is_cancelled() {
            return Err(anyhow::anyhow!("Task cancelled by user"));
        }
    }

    info!("  Reading reference genome...");
    let ref_seq = read_ref_sequence(ref_fasta)?;

    info!("  Parsing marker file and generating diagnostic k-mers...");
    let mut markers = Vec::new();
    let mut indels_skipped = 0usize;
    // Markers the reference cannot give a probe for: the window would run past
    // the end of the sequence, or the allele is longer than k. They are
    // reported rather than left to vanish between the panel and the results.
    let mut without_kmer = 0usize;
    let file = File::open(tsv_markers)
        .map_err(|e| anyhow::anyhow!("Failed to open marker file '{}': {}", tsv_markers, e))?;
    let mut line_count = 0usize;
    // Marker rows, which is not the file's line count: the header, blank lines
    // and comments are not markers, and the number is reported to the user
    // beside how many of them got a k-mer.
    let mut marker_lines = 0usize;
    let mut header_validated = false;
    for line in BufReader::new(file).lines() {
        line_count += 1;
        if line_count % 4096 == 0 {
            if let Some(token) = cancel_token {
                if token.is_cancelled() {
                    return Err(anyhow::anyhow!("Task cancelled by user"));
                }
            }
        }

        let line_str = line?;
        if line_str.trim().is_empty() || line_str.starts_with('#') {
            continue;
        }
        // Validate header row: expect at least position, ref, alt, lineage columns
        if !header_validated {
            let lower = line_str.to_lowercase();
            if lower.contains("pos") || lower.contains("lineage") || lower.contains("ref") {
                info!("  Detected header row, skipping: {}", line_str.chars().take(80).collect::<String>());
                header_validated = true;
                continue;
            }
            header_validated = true; // no header, proceed with data
        }
        marker_lines += 1;
        let fields: Vec<&str> = line_str.split('\t').collect();
        if fields.len() < 4 {
            warn!(
                "Skipping marker line with fewer than 4 columns: {}",
                line_str
            );
            continue;
        }

        let pos0: usize = match fields[0].parse::<usize>() {
            Ok(p) => p.saturating_sub(1),
            Err(_) => {
                warn!("Skipping marker with invalid position: {}", fields[0]);
                continue;
            }
        };
        let ref_allele_str = fields[1].to_string();
        let alt_allele_str = fields[2].to_string();

        // MODIFICATION: Read multiple lineage columns
        let mut lineage_cols: Vec<String> = Vec::new();
        let mut annotation_cols: Vec<String> = Vec::new();
        let mut reading_lineages = true;

        for field in fields[3..].iter() {
            if field.trim().is_empty() {
                reading_lineages = false; // Stop reading lineages after the first empty cell
                continue;
            }
            if reading_lineages {
                lineage_cols.push(field.to_string());
            } else {
                annotation_cols.push(field.to_string());
            }
        }

        if lineage_cols.is_empty() {
            warn!(
                "Skipping marker at pos {} due to no lineage information.",
                pos0 + 1
            );
            continue;
        }

        let ref_allele = ref_allele_str.as_bytes();
        let alt_allele = alt_allele_str.as_bytes();

        let max_allele_len = ref_allele.len().max(alt_allele.len());

        // FASTQ mode: only SNPs and MNVs (ref_len == alt_len).
        // Indels produce unreliable k-mer matches in short reads, especially
        // across repetitive regions. Use the FASTA classify workflow for indels.
        if ref_allele.len() != alt_allele.len() {
            indels_skipped += 1;
            debug!(
                "Skipping indel at pos {} (ref={}bp, alt={}bp) — not supported in FASTQ mode",
                pos0 + 1,
                ref_allele.len(),
                alt_allele.len()
            );
            continue;
        }

        if max_allele_len < kmer_len {
            if let Some((ref_kmer, alt_kmer, kmer_start)) =
                build_small_variant_kmer(pos0, ref_allele, alt_allele, &ref_seq, kmer_len)
            {
                markers.push(Marker {
                    pos: pos0,
                    lineages: lineage_cols,
                    ref_allele: ref_allele_str,
                    alt_allele: alt_allele_str,
                    ref_kmer,
                    alt_kmer,
                    // Filled by one pass over the reference once the panel is
                    // read: see mark_reference_copies.
                    alt_in_reference: false,
                    ref_repeated: false,
                    kmer_start,
                    annotations: annotation_cols,
                });
            } else {
                without_kmer += 1;
                debug!(
                    "No k-mer for the variant at pos {} (the window runs past the end of the reference)",
                    pos0 + 1
                );
            }
        } else {
            without_kmer += 1;
            debug!(
                "Skipping variant at pos {} (allele len {} >= k-mer len {})",
                pos0 + 1,
                max_allele_len,
                kmer_len
            );
        }
    }

    if let Some(token) = cancel_token {
        if token.is_cancelled() {
            return Err(anyhow::anyhow!("Task cancelled by user"));
        }
    }

    mark_reference_copies(&mut markers, &ref_seq, kmer_len);

    info!(
        "  Built {} diagnostic k-mers from {} marker lines{}.",
        markers.len(),
        marker_lines,
        {
            let mut notes = Vec::new();
            if indels_skipped > 0 {
                notes.push(format!(
                    "{} indels skipped — use FASTA classify for indel detection",
                    indels_skipped
                ));
            }
            if without_kmer > 0 {
                notes.push(format!(
                    "{} markers the reference gives no k-mer for (run with -v to list them)",
                    without_kmer
                ));
            }
            if notes.is_empty() {
                String::new()
            } else {
                format!(" ({})", notes.join("; "))
            }
        }
    );

    Ok(markers)
}

pub fn scan_fastq(
    fastqs: &[String],
    markers: &[Marker],
    cancel_token: &Option<CancellationToken>,
    kmer_len: usize,
) -> Result<Vec<[u32; 2]>> {
    if let Some(token) = cancel_token {
        if token.is_cancelled() {
            return Err(anyhow::anyhow!("Task cancelled by user"));
        }
    }

    info!(
        "  Indexing {} markers (including reverse complements)...",
        markers.len()
    );
    let (index, bloom) = build_marker_index(markers);
    scan_fastq_with_index(fastqs, &index, &bloom, markers.len(), cancel_token, kmer_len)
}

fn insert_marker_hit(index: &mut MarkerIndex, hash: u64, hit: (usize, bool)) {
    match index.entry(hash) {
        std::collections::hash_map::Entry::Vacant(entry) => {
            entry.insert(MarkerHitList::Single(hit));
        }
        std::collections::hash_map::Entry::Occupied(mut entry) => {
            let slot = entry.get_mut();
            match slot {
                MarkerHitList::Single(first) => {
                    let existing = *first;
                    *slot = MarkerHitList::Multi(vec![existing, hit]);
                }
                MarkerHitList::Multi(hits) => {
                    hits.push(hit);
                }
            }
        }
    }
}


/// Longest allele an anchor reads. The marker panels reach eleven bases at
/// rpoB 450; anything longer keeps the k-mer answer.
const MAX_ANCHOR_ALLELE: usize = 32;

// --- Anchors: reading the allele from one side ---

/// Which side of the variant an anchor sits on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AnchorSide {
    Left,
    Right,
}

/// One marker reachable from an anchor, and which way round a read carries it.
#[derive(Clone, Copy, Debug)]
struct AnchorRef {
    marker: usize,
    side: AnchorSide,
    /// The read holds the anchor as the reference spells it, rather than its
    /// reverse complement. Which side of the match the allele sits on, and
    /// whether it has to be complemented, both follow from this.
    forward: bool,
}

/// Reads carrying each allele of a marker, as read from one of its anchors.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AnchorTally {
    pub reference: u32,
    pub alternate: u32,
    /// Reads whose bases at the locus are neither of the marker's alleles.
    pub other: u32,
}

impl AnchorTally {
    /// Reads that read the locus at all, whatever they found there.
    pub fn depth(&self) -> u32 {
        self.reference
            .saturating_add(self.alternate)
            .saturating_add(self.other)
    }
}

/// The `k-1` bases on each side of a variant, indexed on their own.
///
/// A k-mer spans the variant and both its flanks, so one other change anywhere
/// in it stops the match and the marker can be read neither way. Each flank on
/// its own survives a change in the other one, and the allele is then read by
/// counting bases from the anchor instead of being recognised whole. That reads
/// alleles nobody catalogued, which an index of expected k-mers cannot do.
///
/// An anchor is one base shorter than a k-mer but says nothing about where it
/// matched, so it is kept only when it occurs **exactly once** in the
/// reference, counting both strands: a repeated anchor cannot say which copy a
/// read came from. An anchor that is its own reverse complement is dropped too,
/// because a read matching it could be either way round and the allele would be
/// on either side.
///
/// Only substitutions get anchors. With an indel the allele's length is what is
/// in question, so bases counted from one side cannot settle which allele is
/// there; split-fastq skips indels anyway.
pub struct AnchorIndex {
    len: usize,
    lookup: FxHashMap<u64, Vec<AnchorRef>>,
    bloom: BloomFilter,
    /// Markers with at least one anchor that is unique in the reference.
    readable: Vec<bool>,
}

impl AnchorIndex {
    pub fn is_empty(&self) -> bool {
        self.lookup.is_empty()
    }

    /// Bases per anchor.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether any anchor of this marker can say what the locus holds. False
    /// means every anchor it has is repeated in the reference, so a read
    /// matching one would not say which copy it came from.
    pub fn can_read(&self, marker_id: usize) -> bool {
        self.readable.get(marker_id).copied().unwrap_or(false)
    }
}

/// Records one anchor of one marker, pooling markers that share the anchor.
fn push_anchor(
    bases: &[u8],
    anchor_ref: AnchorRef,
    slot_of: &mut FxHashMap<u64, usize>,
    slot_refs: &mut Vec<Vec<AnchorRef>>,
    slot_hashes: &mut Vec<(u64, u64)>,
) {
    let (Some(forward), Some(reverse)) = (hash_kmer(bases), hash_reverse_complement(bases)) else {
        return;
    };
    let slot = *slot_of.entry(forward).or_insert_with(|| {
        slot_refs.push(Vec::new());
        slot_hashes.push((forward, reverse));
        slot_refs.len() - 1
    });
    slot_refs[slot].push(anchor_ref);
}

/// Builds the anchors of the given markers from the reference.
pub fn build_anchor_index(
    reference_path: &str,
    markers: &[Marker],
    wanted: &[usize],
    kmer_len: usize,
) -> Result<AnchorIndex> {
    let ref_seq = read_ref_sequence(reference_path)?;
    Ok(build_anchor_index_from_sequence(&ref_seq, markers, wanted, kmer_len))
}

/// The anchors of the given markers, read off a reference already in memory.
pub fn build_anchor_index_from_sequence(
    ref_seq: &[u8],
    markers: &[Marker],
    wanted: &[usize],
    kmer_len: usize,
) -> AnchorIndex {
    let anchor_len = kmer_len.saturating_sub(1);
    let mut index = AnchorIndex {
        len: anchor_len,
        lookup: FxHashMap::default(),
        bloom: bloom_new(),
        readable: vec![false; markers.len()],
    };
    if anchor_len == 0 || wanted.is_empty() {
        return index;
    }

    let seq_len = ref_seq.len();
    let mut slot_of: FxHashMap<u64, usize> = FxHashMap::default();
    let mut slot_refs: Vec<Vec<AnchorRef>> = Vec::new();
    let mut slot_hashes: Vec<(u64, u64)> = Vec::new();

    for &marker_id in wanted {
        let Some(marker) = markers.get(marker_id) else {
            continue;
        };
        let allele_len = marker.ref_allele.len();
        if allele_len != marker.alt_allele.len() {
            continue;
        }
        if marker.pos >= anchor_len {
            push_anchor(
                &ref_seq[marker.pos - anchor_len..marker.pos],
                AnchorRef { marker: marker_id, side: AnchorSide::Left, forward: true },
                &mut slot_of,
                &mut slot_refs,
                &mut slot_hashes,
            );
        }
        let right_start = marker.pos + allele_len;
        if right_start + anchor_len <= seq_len {
            push_anchor(
                &ref_seq[right_start..right_start + anchor_len],
                AnchorRef { marker: marker_id, side: AnchorSide::Right, forward: true },
                &mut slot_of,
                &mut slot_refs,
                &mut slot_hashes,
            );
        }
    }

    // How often the reference carries each anchor, on either strand. One pass
    // over the genome answers it for all of them at once.
    let mut slot_by_hash: FxHashMap<u64, usize> = FxHashMap::default();
    for (slot, &(forward, reverse)) in slot_hashes.iter().enumerate() {
        slot_by_hash.insert(forward, slot);
        slot_by_hash.insert(reverse, slot);
    }
    let mut occurrences = vec![0u32; slot_refs.len()];
    for (_, bitkmer, _) in ref_seq.bit_kmers(anchor_len as u8, false) {
        if let Some(&slot) = slot_by_hash.get(&bitkmer.0) {
            occurrences[slot] = occurrences[slot].saturating_add(1);
        }
    }

    let mut repeated = 0usize;
    for (slot, refs) in slot_refs.into_iter().enumerate() {
        let (forward, reverse) = slot_hashes[slot];
        if occurrences[slot] != 1 || forward == reverse {
            repeated += 1;
            continue;
        }
        for anchor_ref in &refs {
            index.readable[anchor_ref.marker] = true;
        }
        bloom_insert(&mut index.bloom, forward);
        bloom_insert(&mut index.bloom, reverse);
        index
            .lookup
            .entry(reverse)
            .or_default()
            .extend(refs.iter().map(|r| AnchorRef { forward: false, ..*r }));
        index.lookup.entry(forward).or_default().extend(refs);
    }
    if repeated > 0 {
        debug!(
            "{} of the anchors are not unique in the reference and cannot read their locus",
            repeated
        );
    }
    index
}

/// Copies into `out` the bases a read holds where a marker's allele should be,
/// found by counting from an anchor match at `at`.
///
/// From a left anchor the allele begins where the anchor ends; from a right
/// anchor it ends where the anchor begins. A read carrying the reverse
/// complement of the anchor holds the allele on the other side of the match,
/// and reverse complemented.
fn allele_from_anchor(
    read: &[u8],
    at: usize,
    anchor_len: usize,
    side: AnchorSide,
    forward: bool,
    out: &mut [u8],
) -> bool {
    let allele_len = out.len();
    let after = matches!((side, forward), (AnchorSide::Left, true) | (AnchorSide::Right, false));
    let bases = if after {
        let start = at + anchor_len;
        let end = start + allele_len;
        if end > read.len() {
            return false;
        }
        &read[start..end]
    } else {
        if at < allele_len {
            return false;
        }
        &read[at - allele_len..at]
    };
    if forward {
        out.copy_from_slice(bases);
    } else {
        for (slot, &base) in out.iter_mut().zip(bases.iter().rev()) {
            *slot = match base {
                b'A' | b'a' => b'T',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                b'T' | b't' => b'A',
                other => other,
            };
        }
    }
    true
}

/// Adds what one read says at every locus it anchors, into `tallies`.
///
/// `allele` is a scratch buffer, passed in so the hot loop allocates nothing.
fn tally_anchors_in_read(
    read: &[u8],
    markers: &[Marker],
    index: &AnchorIndex,
    allele: &mut [u8; MAX_ANCHOR_ALLELE],
    tallies: &mut [[AnchorTally; 2]],
) {
    for (at, bitkmer, _) in read.bit_kmers(index.len as u8, false) {
        if !bloom_may_contain(&index.bloom, bitkmer.0) {
            continue;
        }
        let Some(refs) = index.lookup.get(&bitkmer.0) else {
            continue;
        };
        for anchor_ref in refs {
            let marker = &markers[anchor_ref.marker];
            let allele_len = marker.ref_allele.len();
            if allele_len == 0 || allele_len > MAX_ANCHOR_ALLELE {
                continue;
            }
            let out = &mut allele[..allele_len];
            if !allele_from_anchor(read, at, index.len, anchor_ref.side, anchor_ref.forward, out) {
                continue;
            }
            let slot = match anchor_ref.side {
                AnchorSide::Left => 0,
                AnchorSide::Right => 1,
            };
            let tally = &mut tallies[anchor_ref.marker][slot];
            if out.eq_ignore_ascii_case(marker.alt_allele.as_bytes()) {
                tally.alternate = tally.alternate.saturating_add(1);
            } else if out.eq_ignore_ascii_case(marker.ref_allele.as_bytes()) {
                tally.reference = tally.reference.saturating_add(1);
            } else {
                tally.other = tally.other.saturating_add(1);
            }
        }
    }
}

/// Counts, per marker and per anchor, which allele the reads hold at its locus.
pub fn scan_fastq_anchors(
    fastqs: &[String],
    markers: &[Marker],
    index: &AnchorIndex,
    cancel_token: &Option<CancellationToken>,
) -> Result<Vec<[AnchorTally; 2]>> {
    let num_markers = markers.len();
    let mut tallies = vec![[AnchorTally::default(); 2]; num_markers];
    if index.is_empty() {
        return Ok(tallies);
    }
    if let Some(token) = cancel_token {
        if token.is_cancelled() {
            return Err(anyhow::anyhow!("Task cancelled by user"));
        }
    }

    info!("  Reading the unread loci from their anchors...");
    const BATCH_SIZE: usize = 8192;
    let cancelled = AtomicBool::new(false);
    let parse_error: Arc<Mutex<Option<String>>> = Arc::new(Mutex::new(None));

    for path in fastqs {
        if cancelled.load(Ordering::Relaxed) {
            break;
        }
        if let Some(token) = cancel_token {
            if token.is_cancelled() {
                cancelled.store(true, Ordering::Relaxed);
                break;
            }
        }
        let mut rdr = match needletail::parse_fastx_file(path) {
            Ok(reader) => reader,
            Err(e) => {
                let mut first_error = parse_error.lock().unwrap_or_else(|e| e.into_inner());
                if first_error.is_none() {
                    *first_error = Some(format!("Could not open or parse FASTQ '{}': {}", path, e));
                }
                cancelled.store(true, Ordering::Relaxed);
                break;
            }
        };
        let mut batch = ReadBatch::with_capacity(BATCH_SIZE, 150);
        let mut done = false;
        loop {
            batch.clear();
            while batch.len() < BATCH_SIZE {
                if let Some(record) = rdr.next() {
                    match record {
                        // seq(), not sequence(): the newlines of a wrapped FASTA
                        // would cut every k-mer that spans one.
                        Ok(record) => batch.push(&record.seq()),
                        Err(e) => {
                            let mut first_error =
                                parse_error.lock().unwrap_or_else(|e| e.into_inner());
                            if first_error.is_none() {
                                *first_error =
                                    Some(format!("Malformed record in '{}': {}", path, e));
                            }
                            cancelled.store(true, Ordering::Relaxed);
                            done = true;
                            break;
                        }
                    }
                } else {
                    done = true;
                    break;
                }
            }
            if batch.is_empty() {
                break;
            }
            if let Some(token) = cancel_token {
                if token.is_cancelled() {
                    cancelled.store(true, Ordering::Relaxed);
                    break;
                }
            }

            let chunk_size = (batch.len() / (rayon::current_num_threads() * 2)).max(64);
            let chunk_ranges = batch.chunk_ranges(chunk_size);
            let batch_tallies = chunk_ranges
                .par_iter()
                .fold(
                    || vec![[AnchorTally::default(); 2]; num_markers],
                    |mut local, &(start, end)| {
                        let mut allele = [0u8; MAX_ANCHOR_ALLELE];
                        for i in start..end {
                            let seq = batch.get(i);
                            tally_anchors_in_read(seq, markers, index, &mut allele, &mut local);
                        }
                        local
                    },
                )
                .reduce(
                    || vec![[AnchorTally::default(); 2]; num_markers],
                    |mut acc, local| {
                        for (slot, value) in acc.iter_mut().zip(local.iter()) {
                            for side in 0..2 {
                                slot[side].reference = slot[side]
                                    .reference
                                    .saturating_add(value[side].reference);
                                slot[side].alternate = slot[side]
                                    .alternate
                                    .saturating_add(value[side].alternate);
                                slot[side].other =
                                    slot[side].other.saturating_add(value[side].other);
                            }
                        }
                        acc
                    },
                );
            for (slot, value) in tallies.iter_mut().zip(batch_tallies.iter()) {
                for side in 0..2 {
                    slot[side].reference =
                        slot[side].reference.saturating_add(value[side].reference);
                    slot[side].alternate =
                        slot[side].alternate.saturating_add(value[side].alternate);
                    slot[side].other = slot[side].other.saturating_add(value[side].other);
                }
            }
            if done {
                break;
            }
        }
    }

    if let Some(message) = parse_error.lock().unwrap_or_else(|e| e.into_inner()).take() {
        return Err(anyhow::anyhow!(message));
    }
    if cancelled.load(Ordering::Relaxed) {
        return Err(anyhow::anyhow!("Task cancelled by user"));
    }
    Ok(tallies)
}

/// Records which markers have a k-mer the reference genome itself carries.
///
/// A k-mer says what it matched and never where, so reads from a paralogous
/// copy count towards the wrong marker. Of the WHO catalogue markers this
/// module indexes, five have an ALT k-mer that already sits elsewhere in H37Rv
/// (three of them in PPE35, one in ndh, one upstream of dnaA), and sixty-seven
/// have a second copy of their reference window. Marking them here
/// is what lets the caller refuse to settle them on read counts alone and ask
/// the anchors, which do say where they matched.
fn mark_reference_copies(markers: &mut [Marker], ref_seq: &[u8], kmer_len: usize) {
    if kmer_len == 0 || ref_seq.len() < kmer_len || markers.is_empty() {
        return;
    }
    let mut alt_of: FxHashMap<u64, Vec<usize>> = FxHashMap::default();
    let mut ref_of: FxHashMap<u64, Vec<usize>> = FxHashMap::default();
    for (id, marker) in markers.iter().enumerate() {
        for (word, map) in [(&marker.alt_kmer, &mut alt_of), (&marker.ref_kmer, &mut ref_of)] {
            if let Some(hash) = hash_kmer(word) {
                map.entry(hash).or_default().push(id);
            }
            if let Some(hash) = hash_reverse_complement(word) {
                // A palindrome hashes the same both ways; recording it once
                // keeps one occurrence from counting as two.
                let entry = map.entry(hash).or_default();
                if entry.last() != Some(&id) {
                    entry.push(id);
                }
            }
        }
    }

    let mut alt_hits = vec![0u32; markers.len()];
    let mut ref_hits = vec![0u32; markers.len()];
    for (_, bitkmer, _) in ref_seq.bit_kmers(kmer_len as u8, false) {
        if let Some(ids) = alt_of.get(&bitkmer.0) {
            for &id in ids {
                alt_hits[id] = alt_hits[id].saturating_add(1);
            }
        }
        if let Some(ids) = ref_of.get(&bitkmer.0) {
            for &id in ids {
                ref_hits[id] = ref_hits[id].saturating_add(1);
            }
        }
    }

    let (mut with_alt_copy, mut with_ref_copy) = (0usize, 0usize);
    for (id, marker) in markers.iter_mut().enumerate() {
        marker.alt_in_reference = alt_hits[id] > 0;
        marker.ref_repeated = ref_hits[id] > 1;
        if marker.alt_in_reference {
            with_alt_copy += 1;
            debug!(
                "The reference already carries the ALT k-mer of the marker at pos {} ({}>{}), so only its anchors can settle it",
                marker.pos + 1, marker.ref_allele, marker.alt_allele
            );
        }
        if marker.ref_repeated {
            with_ref_copy += 1;
        }
    }
    if with_alt_copy > 0 || with_ref_copy > 0 {
        info!(
            "  {} markers have an ALT k-mer the reference already carries and {} a repeated reference window; read counts cannot settle those, so their anchors do.",
            with_alt_copy, with_ref_copy
        );
    }
}

pub fn build_marker_index(markers: &[Marker]) -> (MarkerIndex, BloomFilter) {
    let mut index: MarkerIndex = FxHashMap::default();
    index.reserve(markers.len().saturating_mul(4));
    let mut bloom = bloom_new();
    for (id, marker) in markers.iter().enumerate() {
        if let Some(ref_hash) = hash_kmer(&marker.ref_kmer) {
            insert_marker_hit(&mut index, ref_hash, (id, false));
            bloom_insert(&mut bloom, ref_hash);
        }
        if let Some(alt_hash) = hash_kmer(&marker.alt_kmer) {
            insert_marker_hit(&mut index, alt_hash, (id, true));
            bloom_insert(&mut bloom, alt_hash);
        }
        if let Some(rc_ref_hash) = hash_reverse_complement(&marker.ref_kmer) {
            insert_marker_hit(&mut index, rc_ref_hash, (id, false));
            bloom_insert(&mut bloom, rc_ref_hash);
        }
        if let Some(rc_alt_hash) = hash_reverse_complement(&marker.alt_kmer) {
            insert_marker_hit(&mut index, rc_alt_hash, (id, true));
            bloom_insert(&mut bloom, rc_alt_hash);
        }
    }
    (index, bloom)
}

/// Contiguous buffer for read sequences — avoids per-read allocations.
/// Sequences are stored back-to-back in `data`, with `offsets` marking each start.
struct ReadBatch {
    data: Vec<u8>,
    offsets: Vec<usize>,
}

impl ReadBatch {
    fn with_capacity(num_reads: usize, avg_read_len: usize) -> Self {
        Self {
            data: Vec::with_capacity(num_reads * avg_read_len),
            offsets: Vec::with_capacity(num_reads),
        }
    }

    fn push(&mut self, seq: &[u8]) {
        self.offsets.push(self.data.len());
        self.data.extend_from_slice(seq);
    }

    fn clear(&mut self) {
        self.data.clear();
        self.offsets.clear();
    }

    fn len(&self) -> usize {
        self.offsets.len()
    }

    fn get(&self, i: usize) -> &[u8] {
        let start = self.offsets[i];
        let end = self.offsets.get(i + 1).copied().unwrap_or(self.data.len());
        &self.data[start..end]
    }

    fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    /// Returns ranges for chunked parallel processing.
    fn chunk_ranges(&self, chunk_size: usize) -> Vec<(usize, usize)> {
        (0..self.len())
            .step_by(chunk_size)
            .map(|start| (start, (start + chunk_size).min(self.len())))
            .collect()
    }
}

pub fn scan_fastq_with_index(
    fastqs: &[String],
    index: &MarkerIndex,
    bloom: &[u64; BLOOM_WORDS],
    num_markers: usize,
    cancel_token: &Option<CancellationToken>,
    kmer_len: usize,
) -> Result<Vec<[u32; 2]>> {
    if let Some(token) = cancel_token {
        if token.is_cancelled() {
            return Err(anyhow::anyhow!("Task cancelled by user"));
        }
    }

    // Fast stop flag shared across threads
    let cancelled = AtomicBool::new(false);
    let parse_error: Arc<Mutex<Option<String>>> = Arc::new(Mutex::new(None));

    info!("  Scanning FASTQ reads in parallel...");
    const BATCH_SIZE: usize = 8192;

    let mut counts = vec![[0u32; 2]; num_markers];

    for path in fastqs {
        // Check cancellation at file level
        if cancelled.load(Ordering::Relaxed) {
            break;
        }
        if let Some(token) = cancel_token {
            if token.is_cancelled() {
                cancelled.store(true, Ordering::Relaxed);
                break;
            }
        }

        let mut rdr = match needletail::parse_fastx_file(path) {
            Ok(reader) => reader,
            Err(e) => {
                let mut first_error = parse_error.lock().unwrap_or_else(|e| e.into_inner());
                if first_error.is_none() {
                    *first_error =
                        Some(format!("Could not open or parse FASTQ '{}': {}", path, e));
                }
                cancelled.store(true, Ordering::Relaxed);
                break;
            }
        };

        // Contiguous buffer: avoids per-read Vec<u8> allocations.
        // Buffer memory is reused across batches (clear doesn't free).
        let mut batch = ReadBatch::with_capacity(BATCH_SIZE, 150);
        let mut record_count = 0u64;
        let mut done = false;

        loop {
            // Fill batch from sequential reader
            batch.clear();
            while batch.len() < BATCH_SIZE {
                if let Some(record) = rdr.next() {
                    let record = match record {
                        Ok(rec) => rec,
                        Err(e) => {
                            let mut first_error = parse_error.lock().unwrap_or_else(|e| e.into_inner());
                            if first_error.is_none() {
                                *first_error =
                                    Some(format!("Invalid FASTQ record in '{}': {}", path, e));
                            }
                            cancelled.store(true, Ordering::Relaxed);
                            done = true;
                            break;
                        }
                    };
                    // seq() drops the line breaks a wrapped FASTA carries inside
                    // its sequence; sequence() keeps them, and every k-mer that
                    // spans one was skipped as holding an invalid base. Reads in
                    // FASTQ have none, so for them nothing changes.
                    batch.push(&record.seq());

                    record_count += 1;
                    if record_count & 8191 == 0 {
                        if cancelled.load(Ordering::Relaxed) {
                            done = true;
                            break;
                        }
                        if let Some(token) = cancel_token {
                            if token.is_cancelled() {
                                cancelled.store(true, Ordering::Relaxed);
                                done = true;
                                break;
                            }
                        }
                    }
                } else {
                    done = true;
                    break;
                }
            }

            if !batch.is_empty() {
                // Process batch in parallel: each rayon thread gets thread-local counts.
                let chunk_size = (batch.len() / (rayon::current_num_threads() * 2)).max(64);
                let chunk_ranges = batch.chunk_ranges(chunk_size);
                let batch_counts = chunk_ranges
                    .par_iter()
                    .fold(
                        || vec![[0u32; 2]; num_markers],
                        |mut local, &(start, end)| {
                            for i in start..end {
                                let seq = batch.get(i);
                                // canonical=false: index already contains both forward and RC k-mers
                                for (_, bitkmer, _) in
                                    seq.bit_kmers(kmer_len as u8, false)
                                {
                                    // Bloom filter pre-check: reject ~99% of misses
                                    // without touching the hash table.
                                    if !bloom_may_contain(bloom, bitkmer.0) {
                                        continue;
                                    }
                                    if let Some(hits) = index.get(&bitkmer.0) {
                                        match hits {
                                            MarkerHitList::Single((marker_id, is_alt)) => {
                                                local[*marker_id][*is_alt as usize] += 1;
                                            }
                                            MarkerHitList::Multi(hit_list) => {
                                                for &(marker_id, is_alt) in hit_list {
                                                    local[marker_id][is_alt as usize] += 1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            local
                        },
                    )
                    .reduce(
                        || vec![[0u32; 2]; num_markers],
                        |mut acc, local| {
                            for i in 0..num_markers {
                                acc[i][0] += local[i][0];
                                acc[i][1] += local[i][1];
                            }
                            acc
                        },
                    );

                // Merge batch counts into global counts
                for i in 0..num_markers {
                    counts[i][0] += batch_counts[i][0];
                    counts[i][1] += batch_counts[i][1];
                }
            }

            if done {
                break;
            }
        }
    }

    if let Some(message) = parse_error.lock().unwrap_or_else(|e| e.into_inner()).clone() {
        return Err(anyhow::anyhow!(message));
    }
    if cancelled.load(Ordering::Relaxed) {
        return Err(anyhow::anyhow!("Task cancelled by user"));
    }
    if let Some(token) = cancel_token {
        if token.is_cancelled() {
            return Err(anyhow::anyhow!("Task cancelled by user"));
        }
    }

    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::errors::CancellationToken;

    #[test]
    fn scan_fastq_respects_cancellation_before_work() {
        let token = CancellationToken::new();
        token.cancel();
        let result = scan_fastq(&Vec::new(), &Vec::new(), &Some(token), DEFAULT_MARKER_KMER_LEN);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .to_lowercase()
            .contains("cancel"));
    }

    /// The four ways a read can carry an anchor. A read is either orientation,
    /// and from each side the allele sits on a different side of the match, so
    /// getting any of these backwards would read the wrong bases and still
    /// return an answer.
    #[test]
    fn an_allele_is_read_from_either_anchor_in_either_orientation() {
        const A: usize = 10;
        let left = b"ACGTACGTAC";
        let right = b"GGTTGGTTGG";
        let mut forward = Vec::new();
        forward.extend_from_slice(left);
        forward.extend_from_slice(b"TG");
        forward.extend_from_slice(right);

        let mut allele = [0u8; 2];
        assert!(allele_from_anchor(&forward, 0, A, AnchorSide::Left, true, &mut allele));
        assert_eq!(&allele, b"TG", "left anchor, read as the reference spells it");
        assert!(allele_from_anchor(&forward, 12, A, AnchorSide::Right, true, &mut allele));
        assert_eq!(&allele, b"TG", "right anchor, same read");

        // The same stretch sequenced the other way round: the read holds the
        // reverse complement, so the allele is on the other side of the match
        // and has to be complemented back.
        let reverse: Vec<u8> = forward
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => other,
            })
            .collect();
        assert!(allele_from_anchor(&reverse, 12, A, AnchorSide::Left, false, &mut allele));
        assert_eq!(&allele, b"TG", "left anchor, read the other way round");
        assert!(allele_from_anchor(&reverse, 0, A, AnchorSide::Right, false, &mut allele));
        assert_eq!(&allele, b"TG", "right anchor, read the other way round");
    }

    /// A read that stops before the allele has not read it, and must not be
    /// counted as having read the reference.
    #[test]
    fn a_read_that_stops_short_reads_no_allele() {
        const A: usize = 10;
        let mut allele = [0u8; 2];
        // The anchor is the whole read: nothing follows it.
        assert!(!allele_from_anchor(b"ACGTACGTAC", 0, A, AnchorSide::Left, true, &mut allele));
        // One base short of the allele.
        assert!(!allele_from_anchor(b"ACGTACGTACT", 0, A, AnchorSide::Left, true, &mut allele));
        assert!(allele_from_anchor(b"ACGTACGTACTG", 0, A, AnchorSide::Left, true, &mut allele));
        // The anchor starts the read, so nothing precedes it.
        assert!(!allele_from_anchor(b"ACGTACGTAC", 0, A, AnchorSide::Right, true, &mut allele));
    }

    /// A reference with no repeat in it, long enough for a 31-base window plus
    /// a 30-base anchor on each side.
    fn test_reference() -> Vec<u8> {
        concat!(
            "ACGTTGCAAG", "GCTTAACCGG", "ATCGATTCAG", "CTAGCCATGG", "TACGTTAACG",
            "GCATTGCAGT", "TTGACCGATC", "AGGCTTATCC", "GATTCCAGGT", "CCTAAGGTTC",
            "AGTCCATGAC", "GGTTACCTGA", "TCAGGCATTC", "CCAAGTTGCA", "TGCCAATGGT",
        )
        .as_bytes()
        .to_vec()
    }

    fn one_marker(reference: &[u8], pos0: usize, alt: &str, kmer_len: usize) -> Marker {
        let ref_allele = String::from_utf8(reference[pos0..pos0 + alt.len()].to_vec()).unwrap();
        let (ref_kmer, alt_kmer, kmer_start) = build_small_variant_kmer(
            pos0,
            ref_allele.as_bytes(),
            alt.as_bytes(),
            reference,
            kmer_len,
        )
        .expect("the window fits");
        Marker {
            pos: pos0,
            ref_allele,
            alt_allele: alt.to_string(),
            lineages: vec!["RIF".to_string(), "rpoB_X1Y".to_string()],
            alt_kmer,
            ref_kmer,
            alt_in_reference: false,
            ref_repeated: false,
            kmer_start,
            annotations: Vec::new(),
        }
    }

    /// The case this exists for: the read carries the variant *and* another
    /// change seven bases along, inside the same window. Neither k-mer matches,
    /// so the k-mers alone would report a marker nobody read. The anchor on the
    /// far side of the neighbour still matches, and reads the variant.
    #[test]
    fn an_anchor_reads_a_variant_a_neighbouring_change_hides() {
        let reference = test_reference();
        let k = 31;
        let pos0 = 75;
        let marker = one_marker(&reference, pos0, "A", k);
        let markers = vec![marker];
        let index = build_anchor_index_from_sequence(&reference, &markers, &[0], k);
        assert!(index.can_read(0), "both anchors are unique in this reference");

        // The sample: the variant, plus a neighbour inside the same window.
        let mut sample = reference.clone();
        sample[pos0] = b'A';
        sample[pos0 + 7] = if reference[pos0 + 7] == b'G' { b'C' } else { b'G' };
        // Neither k-mer of the marker is in this read.
        let read = &sample[pos0 - 40..pos0 + 40];
        assert!(
            !read
                .windows(k)
                .any(|w| w == markers[0].alt_kmer.as_slice() || w == markers[0].ref_kmer.as_slice()),
            "the neighbour must break both k-mers, or the test proves nothing"
        );

        let mut allele = [0u8; MAX_ANCHOR_ALLELE];
        let mut tallies = vec![[AnchorTally::default(); 2]; 1];
        tally_anchors_in_read(read, &markers, &index, &mut allele, &mut tallies);
        let [left, right] = tallies[0];
        assert_eq!(left.alternate, 1, "the left anchor is clear of the neighbour");
        assert_eq!(left.reference, 0);
        assert_eq!(right.depth(), 0, "the neighbour sits under the right anchor");
    }

    /// The same shape without the variant: the anchor has to read the reference
    /// and leave the marker absent. A one-base slip in the allele offset would
    /// read the neighbour's base and call something that is not there.
    #[test]
    fn an_anchor_reads_the_reference_when_only_the_neighbour_changed() {
        let reference = test_reference();
        let k = 31;
        let pos0 = 75;
        let markers = vec![one_marker(&reference, pos0, "A", k)];
        let index = build_anchor_index_from_sequence(&reference, &markers, &[0], k);

        let mut sample = reference.clone();
        sample[pos0 + 7] = if reference[pos0 + 7] == b'G' { b'C' } else { b'G' };
        let read = &sample[pos0 - 40..pos0 + 40];

        let mut allele = [0u8; MAX_ANCHOR_ALLELE];
        let mut tallies = vec![[AnchorTally::default(); 2]; 1];
        tally_anchors_in_read(read, &markers, &index, &mut allele, &mut tallies);
        assert_eq!(tallies[0][0].reference, 1, "the left anchor reads the reference base");
        assert_eq!(tallies[0][0].alternate, 0, "and must not read the neighbour as the allele");
    }

    /// A read sequenced the other way round reads the same allele.
    #[test]
    fn an_anchor_reads_the_allele_from_a_reverse_complemented_read() {
        let reference = test_reference();
        let k = 31;
        let pos0 = 75;
        let markers = vec![one_marker(&reference, pos0, "A", k)];
        let index = build_anchor_index_from_sequence(&reference, &markers, &[0], k);

        let mut sample = reference.clone();
        sample[pos0] = b'A';
        let forward = &sample[pos0 - 40..pos0 + 40];
        let reverse: Vec<u8> = forward
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => other,
            })
            .collect();

        let mut allele = [0u8; MAX_ANCHOR_ALLELE];
        let mut tallies = vec![[AnchorTally::default(); 2]; 1];
        tally_anchors_in_read(&reverse, &markers, &index, &mut allele, &mut tallies);
        assert_eq!(
            tallies[0][0].alternate + tallies[0][1].alternate,
            2,
            "both anchors are in this read, each reading the variant once"
        );
        assert_eq!(tallies[0][0].reference + tallies[0][1].reference, 0);
    }

    /// An anchor that the reference carries twice says nothing about which copy
    /// a read came from, so the marker keeps the k-mer answer.
    #[test]
    fn a_repeated_anchor_is_not_used() {
        let k = 11;
        let anchor_len = k - 1;
        // The 10 bases before the variant occur again further along.
        let left = "ACGTTGCAAG";
        let reference = format!("{left}G{left}TTTTTTTTTTTT{left}CCCCCCCCCCCC").into_bytes();
        let pos0 = left.len();
        assert_eq!(reference[pos0], b'G');
        let markers = vec![one_marker(&reference, pos0, "A", k)];
        let index = build_anchor_index_from_sequence(&reference, &markers, &[0], k);
        assert_eq!(index.len(), anchor_len);
        assert!(
            !index.can_read(0),
            "the left anchor is repeated and the right one runs off the end"
        );
    }

    /// A wrapped FASTA carries newlines inside its sequence. Read raw, every
    /// k-mer spanning one held an invalid base and was skipped, so a variant
    /// whose k-mer crossed a line break was never seen.
    #[test]
    fn a_kmer_across_a_line_break_of_a_wrapped_fasta_is_counted() {
        let reference = test_reference();
        let k = 31;
        let pos0 = 75;
        let markers = vec![one_marker(&reference, pos0, "A", k)];
        let (index, bloom) = build_marker_index(&markers);

        let mut sample = reference.clone();
        sample[pos0] = b'A';
        // Lines of 70 bases put the variant at column 76 of the sequence, so
        // its k-mer starts on the first line and ends on the second.
        let mut wrapped = b">sample\n".to_vec();
        for line in sample.chunks(70) {
            wrapped.extend_from_slice(line);
            wrapped.push(b'\n');
        }
        let dir = std::env::temp_dir().join(format!("pathotypr-wrap-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("wrapped.fasta");
        std::fs::write(&path, &wrapped).unwrap();

        let counts = scan_fastq_with_index(
            &[path.to_str().unwrap().to_string()],
            &index,
            &bloom,
            markers.len(),
            &None,
            k,
        )
        .unwrap();
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(counts[0][1], 1, "the alternate k-mer spans the line break and must be counted");
    }

    /// The counts are per allele, so a read carrying neither is neither.
    #[test]
    fn anchor_depth_adds_up_the_three_outcomes() {
        let tally = AnchorTally { reference: 3, alternate: 7, other: 1 };
        assert_eq!(tally.depth(), 11);
        assert_eq!(AnchorTally::default().depth(), 0);
    }
}
