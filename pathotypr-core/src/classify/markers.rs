//! Marker variant parsing, k-mer generation, and genome scanning.
//!
//! Reads marker definitions from TSV, generates diagnostic k-mers against
//! a reference, and scans query genomes for marker matches.

use csv::ReaderBuilder;
use log::{debug, info, trace, warn};
use needletail::Sequence;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::errors::AppResult;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A single marker variant parsed from the markers TSV.
#[derive(Clone, Debug)]
pub struct MarkerVariant {
    pub pos: usize,          // 1-based genomic position
    pub ref_allele: String,
    pub alt_allele: String,
    pub lineage: String,
    pub gene: Option<String>,
    pub mutation: Option<String>,
}

const MAX_ENCODED_MARKER_K: usize = 32;

#[derive(Clone, Debug)]
pub struct MarkerKmerEntry {
    pub kmer: String,
    /// The same window with the reference allele where the alternative one
    /// sits, or `None` when it cannot be built or looked up. Without it the
    /// marker can be seen present but never confirmed absent.
    pub ref_kmer: Option<String>,
    /// The flank immediately before the variant, and the one immediately after
    /// it, both `anchor_len` bases. Only built for substitutions, where the two
    /// alleles have the same length and the allele can therefore be read off a
    /// genome by counting bases from one anchor. See [`AnchorIndex`].
    pub left_anchor: Option<String>,
    pub right_anchor: Option<String>,
    pub ref_position: usize,
    pub lineage: String,
    pub ref_allele_len: usize,
    pub alt_allele_len: usize,
    pub left_flank_len: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub gene: Option<String>,
    pub mutation: Option<String>,
}

/// A k-mer lookup, 2-bit encoded when every k-mer of the panel fits in a u64.
#[derive(Clone, Debug)]
enum KmerLookup {
    Encoded(FxHashMap<u64, usize>),
    Text(FxHashMap<String, usize>),
}

/// Every marker k-mer of a panel, indexed both ways.
///
/// Each marker has an alternative k-mer, which says the variant is there, and
/// where it can be built a reference k-mer, which says it is not. A genome
/// that shows neither has some other change inside the window, and the marker
/// is then undetermined rather than absent.
///
/// Markers whose reference k-mers are identical share a *site*: same position,
/// same reference allele, same window. Reading any one allele of a site reads
/// the position for all of them, so a marker is only undetermined when nothing
/// at its site was found.
#[derive(Clone, Debug)]
pub struct MarkerIndex {
    /// Sorted by reference position, then by alternative k-mer, so every walk
    /// over the panel comes out in the same order.
    markers: Vec<MarkerKmerEntry>,
    /// Alternative k-mer to marker id.
    alt: KmerLookup,
    /// Reference k-mer to site id.
    reference: KmerLookup,
    /// The site each marker belongs to; `None` when it has no reference k-mer.
    marker_site: Vec<Option<usize>>,
    site_count: usize,
    /// The half-length anchors, used only to resolve markers the k-mers left
    /// undetermined.
    anchors: AnchorIndex,
}

// ---------------------------------------------------------------------------
// Anchors
// ---------------------------------------------------------------------------

/// The two flanks of every substitution marker, indexed on their own.
///
/// A k-mer spans the variant and both its flanks, so a single other change
/// anywhere in it stops the match. Each flank on its own survives a change in
/// the other one, and the allele can then be read off the genome by counting
/// bases from the anchor instead of being recognised as a whole. That reads
/// alleles nobody catalogued, which an index of expected k-mers cannot do.
///
/// Only substitutions get anchors. With an indel the allele's length is what
/// is in question, so bases counted from one anchor cannot say which allele is
/// there; those markers keep the k-mer answer.
///
/// An anchor is half a k-mer and therefore much less specific, so an allele is
/// read only from an anchor that occurs **exactly once** in the genome. A
/// repeated anchor says nothing about which copy was read.
#[derive(Clone, Debug, Default)]
pub struct AnchorIndex {
    /// Bases per anchor; 0 when there are none.
    len: usize,
    left: Option<KmerLookup>,
    right: Option<KmerLookup>,
    /// Markers reachable from each slot, and the longest allele to read there.
    left_markers: Vec<Vec<usize>>,
    right_markers: Vec<Vec<usize>>,
    left_span: Vec<usize>,
    right_span: Vec<usize>,
    /// Slot of each marker, when it has anchors.
    marker_left: Vec<Option<usize>>,
    marker_right: Vec<Option<usize>>,
}

/// Where an allele read from the genome came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AnchorEvidence {
    Left,
    Right,
    Both,
}

impl AnchorEvidence {
    pub fn as_str(self) -> &'static str {
        match self {
            AnchorEvidence::Left => "left_anchor",
            AnchorEvidence::Right => "right_anchor",
            AnchorEvidence::Both => "both_anchors",
        }
    }
}

/// What the anchors say about one marker.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AnchorRead {
    /// No anchor of this marker occurs exactly once, so nothing was read.
    Unusable,
    /// The two anchors disagree, so neither is trusted.
    Conflict,
    /// The bases found where the allele should be.
    Allele(String, AnchorEvidence),
}

/// Anchor occurrences accumulated over every contig and strand of one genome.
#[derive(Clone, Debug)]
pub struct AnchorScan {
    left_count: Vec<u32>,
    right_count: Vec<u32>,
    left_bases: Vec<Option<Vec<u8>>>,
    right_bases: Vec<Option<Vec<u8>>>,
}

impl AnchorScan {
    pub fn new(index: &AnchorIndex) -> Self {
        Self {
            left_count: vec![0; index.left_markers.len()],
            right_count: vec![0; index.right_markers.len()],
            left_bases: vec![None; index.left_markers.len()],
            right_bases: vec![None; index.right_markers.len()],
        }
    }

    /// The allele of a marker, read from whichever anchors are usable.
    pub fn allele_of(&self, index: &AnchorIndex, marker_id: usize, allele_len: usize) -> AnchorRead {
        let from_left = index
            .marker_left
            .get(marker_id)
            .copied()
            .flatten()
            .filter(|&slot| self.left_count[slot] == 1)
            .and_then(|slot| self.left_bases[slot].as_ref())
            .filter(|bases| bases.len() >= allele_len)
            // the allele starts where the left anchor ends
            .map(|bases| bases[..allele_len].to_vec());
        let from_right = index
            .marker_right
            .get(marker_id)
            .copied()
            .flatten()
            .filter(|&slot| self.right_count[slot] == 1)
            .and_then(|slot| self.right_bases[slot].as_ref())
            .filter(|bases| bases.len() >= allele_len)
            // and ends where the right anchor begins
            .map(|bases| bases[bases.len() - allele_len..].to_vec());
        let (bases, evidence) = match (from_left, from_right) {
            (None, None) => return AnchorRead::Unusable,
            (Some(l), Some(r)) if l != r => return AnchorRead::Conflict,
            (Some(l), Some(_)) => (l, AnchorEvidence::Both),
            (Some(l), None) => (l, AnchorEvidence::Left),
            (None, Some(r)) => (r, AnchorEvidence::Right),
        };
        AnchorRead::Allele(String::from_utf8_lossy(&bases).into_owned(), evidence)
    }
}

impl AnchorIndex {
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn len(&self) -> usize {
        self.len
    }
}

/// Records every anchor occurrence of one sequence into `scan`.
///
/// Called once per contig and strand, so the counts add up over the whole
/// genome and a locus present twice is seen as repeated.
pub fn scan_anchors(genome_sequence: &str, index: &AnchorIndex, scan: &mut AnchorScan) {
    let (left, right) = match (index.left.as_ref(), index.right.as_ref()) {
        (Some(l), Some(r)) => (l, r),
        _ => return,
    };
    let a = index.len;
    let bytes = genome_sequence.as_bytes();
    let seq_len = bytes.len();
    if a == 0 || seq_len < a {
        return;
    }

    let mut record = |slot_left: Option<usize>, slot_right: Option<usize>, start: usize| {
        if let Some(slot) = slot_left {
            scan.left_count[slot] = scan.left_count[slot].saturating_add(1);
            if scan.left_bases[slot].is_none() {
                let from = start + a;
                let to = (from + index.left_span[slot]).min(seq_len);
                if from <= to {
                    scan.left_bases[slot] = Some(bytes[from..to].to_vec());
                }
            }
        }
        if let Some(slot) = slot_right {
            scan.right_count[slot] = scan.right_count[slot].saturating_add(1);
            if scan.right_bases[slot].is_none() {
                let span = index.right_span[slot];
                let from = start.saturating_sub(span);
                if start >= span {
                    scan.right_bases[slot] = Some(bytes[from..start].to_vec());
                }
            }
        }
    };

    match (left, right) {
        (KmerLookup::Encoded(left), KmerLookup::Encoded(right)) => {
            for (start, bitkmer, _) in bytes.bit_kmers(a as u8, false) {
                let l = left.get(&bitkmer.0).copied();
                let r = right.get(&bitkmer.0).copied();
                if l.is_some() || r.is_some() {
                    record(l, r, start);
                }
            }
        }
        (KmerLookup::Text(left), KmerLookup::Text(right)) => {
            for start in 0..=(seq_len - a) {
                let window = &genome_sequence[start..start + a];
                let l = left.get(window).copied();
                let r = right.get(window).copied();
                if l.is_some() || r.is_some() {
                    record(l, r, start);
                }
            }
        }
        _ => unreachable!("anchor lookups built with different storage"),
    }
}

impl MarkerIndex {
    pub fn len(&self) -> usize {
        self.markers.len()
    }

    pub fn is_empty(&self) -> bool {
        self.markers.is_empty()
    }

    pub fn storage_name(&self) -> &'static str {
        match self.alt {
            KmerLookup::Encoded(_) => "u64-encoded",
            KmerLookup::Text(_) => "string",
        }
    }

    /// Markers in index order; a marker id is its position in this slice.
    pub fn markers(&self) -> &[MarkerKmerEntry] {
        &self.markers
    }

    /// The site of a marker, or `None` when it has no reference k-mer.
    pub fn site_of(&self, marker_id: usize) -> Option<usize> {
        self.marker_site.get(marker_id).copied().flatten()
    }

    pub fn site_count(&self) -> usize {
        self.site_count
    }

    pub fn anchors(&self) -> &AnchorIndex {
        &self.anchors
    }
}

/// Result of a marker k-mer match against a genome sequence.
#[derive(Clone, Debug)]
pub struct MarkerMatch {
    /// Position of the marker in [`MarkerIndex::markers`].
    pub marker_id: usize,
    pub kmer: String,
    pub genome_position: usize,
    pub ref_position: usize,
    pub lineage: String,
    pub ref_allele_len: usize,
    pub alt_allele_len: usize,
    pub left_flank_len: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub gene: Option<String>,
    pub mutation: Option<String>,
}

// ---------------------------------------------------------------------------
// K-mer encoding
// ---------------------------------------------------------------------------

fn encode_base_2bit(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn encode_kmer_2bit(kmer: &[u8]) -> Option<u64> {
    let mut encoded = 0u64;
    for &base in kmer {
        encoded = (encoded << 2) | encode_base_2bit(base)?;
    }
    Some(encoded)
}

// ---------------------------------------------------------------------------
// Marker TSV parsing
// ---------------------------------------------------------------------------

/// Reads marker definitions from a TSV file.
pub fn get_positions(tsv_file: &str) -> AppResult<Vec<MarkerVariant>> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .from_path(tsv_file)?;
    let mut markers = Vec::new();
    for result in rdr.records() {
        let record = result?;
        if record.len() < 4 {
            continue;
        }

        if record.get(0).map(|s| s.trim().starts_with('#')).unwrap_or(false) {
            continue;
        }

        let pos: usize = match record[0].parse() {
            Ok(p) => p,
            Err(_) => {
                debug!("Skipping non-numeric marker position row: {:?}", record);
                continue;
            }
        };
        let ref_allele = record[1].trim();
        let alt_allele = record[2].trim();
        if alt_allele.is_empty() {
            warn!("Skipping marker at position {} because ALT allele is empty.", pos);
            continue;
        }
        if ref_allele.is_empty() {
            warn!("Skipping marker at position {} because REF allele is empty.", pos);
            continue;
        }

        let cols_from_3: Vec<&str> = record.iter().skip(3).collect();
        let lineage_path = cols_from_3
            .iter()
            .take_while(|&&s| !s.trim().is_empty())
            .copied()
            .collect::<Vec<&str>>()
            .join(";");
        if lineage_path.is_empty() {
            warn!("Skipping marker at position {} because lineage path is empty.", pos);
            continue;
        }

        let empty_idx = cols_from_3.iter().position(|s| s.trim().is_empty());
        let (gene, mutation) = if let Some(idx) = empty_idx {
            let gene = cols_from_3.get(idx + 1)
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string());
            let mutation = cols_from_3.get(idx + 2)
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string());
            (gene, mutation)
        } else {
            (None, None)
        };

        markers.push(MarkerVariant {
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            lineage: lineage_path,
            gene,
            mutation,
        });
    }
    Ok(markers)
}

// ---------------------------------------------------------------------------
// K-mer generation
// ---------------------------------------------------------------------------

/// Generate marker k-mers from variants against a reference sequence.
pub fn generate_markerkmer(
    markers: &[MarkerVariant],
    ref_seq: &str,
    k: usize,
    min_flank_bases: usize,
) -> MarkerIndex {
    let mut text_index: FxHashMap<String, MarkerKmerEntry> = FxHashMap::default();
    let mut encoded_index: FxHashMap<u64, MarkerKmerEntry> = FxHashMap::default();
    let mut use_encoded_index = k > 0 && k <= MAX_ENCODED_MARKER_K;
    let mut fallback_reason: Option<&'static str> = if use_encoded_index {
        None
    } else {
        Some("k-mer size is out of range for u64 encoding")
    };

    let effective_max_allele = k.saturating_sub(2 * min_flank_bases);
    // One base short of a k-mer, entirely on one side of the variant. Half that
    // length tolerates a neighbour slightly further away, but 4^15 is only a few
    // hundred times the size of a bacterial genome, so some anchors then occur
    // somewhere else by chance; when the real one is broken and the chance one
    // is not, the allele gets read from the wrong locus. At k-1 that stops
    // happening, and a single neighbour on either side still leaves one anchor
    // standing, which is the point of having two.
    let anchor_len = k - 1;
    let mut skipped_flank: usize = 0;

    let seq_bytes = ref_seq.as_bytes();
    let seq_len = seq_bytes.len();

    for marker in markers {
        let pos = marker.pos;
        let ref_allele = marker.ref_allele.as_bytes();
        let alt_allele = marker.alt_allele.as_bytes();
        let ref_len = ref_allele.len();
        let alt_len = alt_allele.len();
        let max_allele_len = ref_len.max(alt_len);

        if min_flank_bases > 0 && max_allele_len > effective_max_allele {
            skipped_flank += 1;
            debug!(
                "Skipping marker at pos {} (allele len {} > max {} for min_flank {})",
                pos, max_allele_len, effective_max_allele, min_flank_bases
            );
            continue;
        }

        if max_allele_len >= k {
            warn!(
                "Skipping marker at position {} because allele length ({}) >= k-mer size ({}).",
                pos, max_allele_len, k
            );
            continue;
        }

        let pos0 = pos.saturating_sub(1);
        let flank_len = (k - max_allele_len) / 2;
        let kmer_start = pos0.saturating_sub(flank_len);
        let left_flank_len = pos0 - kmer_start;
        let right_flank_start = pos0 + ref_len;
        let alt_right_needed = k.saturating_sub(left_flank_len + alt_len);

        if right_flank_start + alt_right_needed > seq_len {
            warn!("Skipping marker at position {} because k-mer window exceeds reference length.", pos);
            continue;
        }
        if kmer_start + k > seq_len {
            continue;
        }

        let mut kmer_bytes = Vec::with_capacity(k);
        kmer_bytes.extend_from_slice(&seq_bytes[kmer_start..pos0]);
        kmer_bytes.extend_from_slice(alt_allele);
        kmer_bytes.extend_from_slice(&seq_bytes[right_flank_start..right_flank_start + alt_right_needed]);

        debug_assert_eq!(kmer_bytes.len(), k,
            "K-mer length mismatch at position {}: expected {}, got {} (ref={}, alt={})",
            pos, k, kmer_bytes.len(), marker.ref_allele, marker.alt_allele
        );
        if kmer_bytes.len() != k {
            warn!("Skipping marker at position {} due to k-mer length mismatch: expected {}, got {}.", pos, k, kmer_bytes.len());
            continue;
        }

        let new_kmer = String::from_utf8_lossy(&kmer_bytes).into_owned();
        trace!("Generated marker k-mer '{}' for position {} (ref={}, alt={})", new_kmer, pos, marker.ref_allele, marker.alt_allele);

        // The reference k-mer uses the same left flank and the reference
        // allele from the marker file, the way split-fastq builds its pair. An
        // insertion needs more right flank for the shorter reference allele,
        // which can run off the end of the reference.
        let ref_right_needed = k.saturating_sub(left_flank_len + ref_len);
        let ref_kmer = (right_flank_start + ref_right_needed <= seq_len).then(|| {
            let mut ref_bytes = Vec::with_capacity(k);
            ref_bytes.extend_from_slice(&seq_bytes[kmer_start..pos0]);
            ref_bytes.extend_from_slice(ref_allele);
            ref_bytes.extend_from_slice(&seq_bytes[right_flank_start..right_flank_start + ref_right_needed]);
            String::from_utf8_lossy(&ref_bytes).into_owned()
        }).filter(|r| r.len() == k);

        // Anchors are half a k-mer on each side of the variant, the same length
        // for every marker so one pass over the genome reads them all. Only
        // substitutions get them: with an indel the allele's own length is what
        // is in question, so counting bases from an anchor cannot resolve it.
        let (left_anchor, right_anchor) = if ref_len == alt_len
            && anchor_len > 0
            && pos0 >= anchor_len
            && right_flank_start + anchor_len <= seq_len
        {
            (
                Some(String::from_utf8_lossy(&seq_bytes[pos0 - anchor_len..pos0]).into_owned()),
                Some(String::from_utf8_lossy(&seq_bytes[right_flank_start..right_flank_start + anchor_len]).into_owned()),
            )
        } else {
            (None, None)
        };

        let entry = MarkerKmerEntry {
            kmer: new_kmer.clone(),
            ref_kmer,
            left_anchor,
            right_anchor,
            ref_position: pos,
            lineage: marker.lineage.clone(),
            ref_allele_len: ref_len,
            alt_allele_len: alt_len,
            left_flank_len,
            ref_allele: marker.ref_allele.clone(),
            alt_allele: marker.alt_allele.clone(),
            gene: marker.gene.clone(),
            mutation: marker.mutation.clone(),
        };

        if use_encoded_index {
            if let Some(encoded_kmer) = encode_kmer_2bit(new_kmer.as_bytes()) {
                if let Some(prev) = encoded_index.insert(encoded_kmer, entry) {
                    warn!("Duplicate marker k-mer at position {}: overwrites previous marker at position {} (lineage '{}' replaces '{}')",
                        pos, prev.ref_position, marker.lineage, prev.lineage);
                }
            } else {
                use_encoded_index = false;
                fallback_reason = Some("marker contains non-ACGT nucleotide(s) not encodable as u64");
                for (_, e) in encoded_index.drain() {
                    text_index.insert(e.kmer.clone(), e);
                }
                text_index.insert(new_kmer, entry);
            }
        } else if let Some(prev) = text_index.insert(new_kmer, entry) {
            warn!("Duplicate marker k-mer at position {}: overwrites previous marker at position {} (lineage '{}' replaces '{}')",
                pos, prev.ref_position, marker.lineage, prev.lineage);
        }
    }

    if skipped_flank > 0 {
        info!("Skipped {} markers with allele length > {} (min_flank_bases={}, k={}).",
            skipped_flank, effective_max_allele, min_flank_bases, k);
    }

    let entries: Vec<MarkerKmerEntry> = if use_encoded_index {
        encoded_index.into_values().collect()
    } else {
        if let Some(reason) = fallback_reason {
            info!("Using string marker index fallback: {}", reason);
        }
        text_index.into_values().collect()
    };
    build_index(entries, use_encoded_index, anchor_len)
}

/// Turns the final set of markers into the two-way index.
///
/// Built from the markers that survived duplicate alternative k-mers, so the
/// overwrite rule above decides which marker owns a k-mer, as it always has.
fn build_index(mut markers: Vec<MarkerKmerEntry>, encoded: bool, anchor_len: usize) -> MarkerIndex {
    markers.sort_by(|a, b| a.ref_position.cmp(&b.ref_position).then_with(|| a.kmer.cmp(&b.kmer)));

    let mut alt_encoded: FxHashMap<u64, usize> = FxHashMap::default();
    let mut alt_text: FxHashMap<String, usize> = FxHashMap::default();
    let mut ref_encoded: FxHashMap<u64, usize> = FxHashMap::default();
    let mut ref_text: FxHashMap<String, usize> = FxHashMap::default();
    let mut marker_site = Vec::with_capacity(markers.len());
    let mut site_count = 0usize;
    let mut without_reference = 0usize;
    let mut indistinguishable = 0usize;

    for (id, marker) in markers.iter_mut().enumerate() {
        if encoded {
            if let Some(key) = encode_kmer_2bit(marker.kmer.as_bytes()) {
                alt_encoded.insert(key, id);
            }
            // The flanks were encodable for the alternative k-mer, so only a
            // reference allele outside A, C, G and T can fail here.
            if marker.ref_kmer.as_deref().is_some_and(|r| encode_kmer_2bit(r.as_bytes()).is_none()) {
                marker.ref_kmer = None;
            }
        } else {
            alt_text.insert(marker.kmer.clone(), id);
        }

        let site = match marker.ref_kmer.as_deref() {
            None => {
                without_reference += 1;
                None
            }
            Some(ref_kmer) => {
                if ref_kmer == marker.kmer {
                    indistinguishable += 1;
                }
                let next = site_count;
                let site = if encoded {
                    let key = encode_kmer_2bit(ref_kmer.as_bytes())
                        .expect("non-encodable reference k-mers were dropped above");
                    *ref_encoded.entry(key).or_insert(next)
                } else {
                    *ref_text.entry(ref_kmer.to_string()).or_insert(next)
                };
                if site == next {
                    site_count += 1;
                }
                Some(site)
            }
        };
        marker_site.push(site);
    }

    if without_reference > 0 {
        warn!(
            "{} markers have no reference k-mer (the window runs past the end of the reference, or the \
             reference allele holds a base other than A, C, G or T). They can be found present but never \
             confirmed absent, so a genome without them lists them as undetermined.",
            without_reference
        );
    }
    if indistinguishable > 0 {
        warn!(
            "{} markers have the same k-mer for both alleles, so the two cannot be told apart. A genome \
             carrying either is reported with both alleles found.",
            indistinguishable
        );
    }
    debug!("Indexed {} markers at {} reference sites.", markers.len(), site_count);

    let (alt, reference) = if encoded {
        (KmerLookup::Encoded(alt_encoded), KmerLookup::Encoded(ref_encoded))
    } else {
        (KmerLookup::Text(alt_text), KmerLookup::Text(ref_text))
    };
    let anchors = build_anchor_index(&markers, encoded, anchor_len);
    MarkerIndex { markers, alt, reference, marker_site, site_count, anchors }
}

/// Collects the anchors of every substitution marker into their own lookups.
fn build_anchor_index(markers: &[MarkerKmerEntry], encoded: bool, anchor_len: usize) -> AnchorIndex {
    let mut index = AnchorIndex {
        len: anchor_len,
        marker_left: vec![None; markers.len()],
        marker_right: vec![None; markers.len()],
        ..Default::default()
    };
    if anchor_len == 0 {
        return index;
    }
    let mut left_encoded: FxHashMap<u64, usize> = FxHashMap::default();
    let mut left_text: FxHashMap<String, usize> = FxHashMap::default();
    let mut right_encoded: FxHashMap<u64, usize> = FxHashMap::default();
    let mut right_text: FxHashMap<String, usize> = FxHashMap::default();

    for (id, marker) in markers.iter().enumerate() {
        let allele_len = marker.alt_allele_len;
        for side in [true, false] {
            let anchor = if side { &marker.left_anchor } else { &marker.right_anchor };
            let Some(anchor) = anchor.as_deref() else { continue };
            let (slots, spans) = if side {
                (&mut index.left_markers, &mut index.left_span)
            } else {
                (&mut index.right_markers, &mut index.right_span)
            };
            let next = slots.len();
            let slot = if encoded {
                let Some(key) = encode_kmer_2bit(anchor.as_bytes()) else { continue };
                let map = if side { &mut left_encoded } else { &mut right_encoded };
                *map.entry(key).or_insert(next)
            } else {
                let map = if side { &mut left_text } else { &mut right_text };
                *map.entry(anchor.to_string()).or_insert(next)
            };
            if slot == next {
                slots.push(Vec::new());
                spans.push(0);
            }
            slots[slot].push(id);
            spans[slot] = spans[slot].max(allele_len);
            if side {
                index.marker_left[id] = Some(slot);
            } else {
                index.marker_right[id] = Some(slot);
            }
        }
    }

    if index.left_markers.is_empty() && index.right_markers.is_empty() {
        index.len = 0;
        return index;
    }
    let (left, right) = if encoded {
        (KmerLookup::Encoded(left_encoded), KmerLookup::Encoded(right_encoded))
    } else {
        (KmerLookup::Text(left_text), KmerLookup::Text(right_text))
    };
    index.left = Some(left);
    index.right = Some(right);
    debug!(
        "Indexed {}-base anchors for {} markers at {} left and {} right slots.",
        anchor_len,
        index.marker_left.iter().filter(|s| s.is_some()).count(),
        index.left_markers.len(),
        index.right_markers.len()
    );
    index
}

// ---------------------------------------------------------------------------
// Genome scanning
// ---------------------------------------------------------------------------

/// What one pass over one sequence found.
#[derive(Clone, Debug, Default)]
pub struct MarkerScan {
    /// Markers whose alternative k-mer occurs, once each, in sequence order.
    pub matches: Vec<MarkerMatch>,
    /// Sites whose reference k-mer occurs, once each.
    pub reference_sites: Vec<usize>,
}

/// Records what one k-mer of the genome turned out to be.
fn record_hit(
    marker_index: &MarkerIndex,
    marker: Option<usize>,
    site: Option<usize>,
    genome_position: usize,
    seen_markers: &mut FxHashSet<usize>,
    seen_sites: &mut FxHashSet<usize>,
    scan: &mut MarkerScan,
) {
    if let Some(id) = marker {
        if seen_markers.insert(id) {
            let entry = &marker_index.markers[id];
            trace!("Found marker {} at genome position {}", entry.kmer, genome_position);
            scan.matches.push(MarkerMatch {
                marker_id: id,
                kmer: entry.kmer.clone(),
                genome_position,
                ref_position: entry.ref_position,
                lineage: entry.lineage.clone(),
                ref_allele_len: entry.ref_allele_len,
                alt_allele_len: entry.alt_allele_len,
                left_flank_len: entry.left_flank_len,
                ref_allele: entry.ref_allele.clone(),
                alt_allele: entry.alt_allele.clone(),
                gene: entry.gene.clone(),
                mutation: entry.mutation.clone(),
            });
        }
    }
    if let Some(site) = site {
        if seen_sites.insert(site) {
            scan.reference_sites.push(site);
        }
    }
}

/// Scan a genome sequence for both k-mers of every marker in a single pass.
pub fn scan_markers(
    genome_sequence: &str,
    marker_index: &MarkerIndex,
    k: usize,
) -> MarkerScan {
    debug!("Scanning genome for marker k-mer matches.");
    let mut scan = MarkerScan::default();

    let seq_len = genome_sequence.len();
    if k == 0 || seq_len < k {
        return scan;
    }

    let mut seen_markers: FxHashSet<usize> = FxHashSet::default();
    let mut seen_sites: FxHashSet<usize> = FxHashSet::default();
    match (&marker_index.alt, &marker_index.reference) {
        (KmerLookup::Encoded(alt), KmerLookup::Encoded(reference)) => {
            for (genome_position, bitkmer, _) in genome_sequence.as_bytes().bit_kmers(k as u8, false) {
                let marker = alt.get(&bitkmer.0).copied();
                let site = reference.get(&bitkmer.0).copied();
                if marker.is_some() || site.is_some() {
                    record_hit(marker_index, marker, site, genome_position, &mut seen_markers, &mut seen_sites, &mut scan);
                }
            }
        }
        (KmerLookup::Text(alt), KmerLookup::Text(reference)) => {
            for genome_position in 0..=(seq_len - k) {
                let kmer = &genome_sequence[genome_position..genome_position + k];
                let marker = alt.get(kmer).copied();
                let site = reference.get(kmer).copied();
                if marker.is_some() || site.is_some() {
                    record_hit(marker_index, marker, site, genome_position, &mut seen_markers, &mut seen_sites, &mut scan);
                }
            }
        }
        // `build_index` always builds both lookups with the same storage.
        _ => unreachable!("marker lookups built with different storage"),
    }
    debug!(
        "Found {} alternative and {} reference k-mers in genome.",
        scan.matches.len(),
        scan.reference_sites.len()
    );
    scan
}

/// Scan a genome sequence for marker k-mer matches (alternative alleles only).
pub fn find_markers(
    genome_sequence: &str,
    marker_index: &MarkerIndex,
    k: usize,
) -> Vec<MarkerMatch> {
    scan_markers(genome_sequence, marker_index, k).matches
}
