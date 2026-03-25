//! Marker variant parsing, k-mer generation, and genome scanning.
//!
//! Reads marker definitions from TSV, generates diagnostic k-mers against
//! a reference, and scans query genomes for marker matches.

use csv::ReaderBuilder;
use log::{debug, info, trace, warn};
use needletail::Sequence;
use rustc_hash::FxHashMap;
use std::collections::HashMap;

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

#[derive(Clone, Debug)]
pub enum MarkerIndex {
    Encoded(FxHashMap<u64, MarkerKmerEntry>),
    Text(FxHashMap<String, MarkerKmerEntry>),
}

impl MarkerIndex {
    pub fn len(&self) -> usize {
        match self {
            MarkerIndex::Encoded(map) => map.len(),
            MarkerIndex::Text(map) => map.len(),
        }
    }

    pub fn storage_name(&self) -> &'static str {
        match self {
            MarkerIndex::Encoded(_) => "u64-encoded",
            MarkerIndex::Text(_) => "string",
        }
    }
}

/// Result of a marker k-mer match against a genome sequence.
#[derive(Clone, Debug)]
pub struct MarkerMatch {
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

        let entry = MarkerKmerEntry {
            kmer: new_kmer.clone(),
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

    if use_encoded_index {
        MarkerIndex::Encoded(encoded_index)
    } else {
        if let Some(reason) = fallback_reason {
            info!("Using string marker index fallback: {}", reason);
        }
        MarkerIndex::Text(text_index)
    }
}

// ---------------------------------------------------------------------------
// Genome scanning
// ---------------------------------------------------------------------------

/// Scan a genome sequence for marker k-mer matches.
pub fn find_markers(
    genome_sequence: &str,
    marker_index: &MarkerIndex,
    k: usize,
) -> Vec<MarkerMatch> {
    debug!("Scanning genome for marker k-mer matches.");
    let mut matched_markers: HashMap<String, MarkerMatch> = HashMap::new();

    let seq_len = genome_sequence.len();
    if k == 0 || seq_len < k {
        return Vec::new();
    }

    match marker_index {
        MarkerIndex::Encoded(encoded_index) => {
            for (genome_position, bitkmer, _) in genome_sequence.as_bytes().bit_kmers(k as u8, false) {
                if let Some(entry) = encoded_index.get(&bitkmer.0) {
                    trace!("Found marker {} at genome position {}", entry.kmer, genome_position);
                    matched_markers.entry(entry.kmer.clone()).or_insert(MarkerMatch {
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
        }
        MarkerIndex::Text(text_index) => {
            for genome_position in 0..=(seq_len - k) {
                let kmer = &genome_sequence[genome_position..genome_position + k];
                if let Some(entry) = text_index.get(kmer) {
                    trace!("Found marker {} at genome position {}", entry.kmer, genome_position);
                    matched_markers.entry(entry.kmer.clone()).or_insert(MarkerMatch {
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
        }
    }
    debug!("Found {} matched markers in genome.", matched_markers.len());
    matched_markers.into_values().collect()
}
