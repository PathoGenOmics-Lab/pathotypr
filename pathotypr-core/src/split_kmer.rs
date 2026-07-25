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
) -> Option<(Vec<u8>, Vec<u8>)> {
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

    Some((ref_kmer, alt_kmer))
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
    let file = File::open(tsv_markers)
        .map_err(|e| anyhow::anyhow!("Failed to open marker file '{}': {}", tsv_markers, e))?;
    let mut line_count = 0usize;
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
            if let Some((ref_kmer, alt_kmer)) =
                build_small_variant_kmer(pos0, ref_allele, alt_allele, &ref_seq, kmer_len)
            {
                markers.push(Marker {
                    pos: pos0,
                    lineages: lineage_cols,
                    ref_allele: ref_allele_str,
                    alt_allele: alt_allele_str,
                    ref_kmer,
                    alt_kmer,
                    annotations: annotation_cols,
                });
            }
        } else {
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

    info!(
        "  Built {} diagnostic k-mers from {} marker lines{}.",
        markers.len(),
        line_count,
        if indels_skipped > 0 {
            format!(" ({} indels skipped — use FASTA classify for indel detection)", indels_skipped)
        } else {
            String::new()
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
                    batch.push(record.sequence());

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
}
