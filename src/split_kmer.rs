//! Implements a dynamic variant detection engine.
//!
//! This module dynamically creates diagnostic k-mers for SNPs, MNVs, and Indels
//! (both small and large) based on a user-provided marker file. It then scans
//! FASTQ reads to count evidence for reference vs. alternate alleles.

use anyhow::Result;
use bio::io::fasta::{self, FastaRead};
use fxhash::FxHashMap;
use log::{info, warn};
use needletail::kmer;
use needletail::Sequence;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// The total length of the diagnostic k-mers used for identification.
const MARKER_KMER_LEN: usize = 31;

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
        val = (val << 2) | match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
    }
    Some(val)
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut rc_seq: Vec<u8> = Vec::with_capacity(seq.len());
    for &base in seq.iter().rev() {
        rc_seq.push(match base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        });
    }
    rc_seq
}

fn read_ref_sequence(path: &str) -> Result<Vec<u8>> {
    let mut reader = fasta::Reader::from_file(path)?;
    let mut record = fasta::Record::new();
    reader.read(&mut record)?;
    Ok(record.seq().to_ascii_uppercase())
}

/// Builds a diagnostic k-mer for a small variant (SNP, MNV, small indel).
fn build_small_variant_kmer(
    pos0: usize,
    ref_allele: &[u8],
    alt_allele: &[u8],
    ref_seq: &[u8],
) -> Option<(Vec<u8>, Vec<u8>)> {
    let ref_len_allele = ref_allele.len();
    let max_allele_len = ref_len_allele.max(alt_allele.len());
    let flank_len = (MARKER_KMER_LEN - max_allele_len) / 2;
    let kmer_start_pos = pos0.saturating_sub(flank_len);
    let right_flank_start_ref = pos0 + ref_len_allele;
    let needed_right_flank = MARKER_KMER_LEN.saturating_sub(pos0.saturating_sub(kmer_start_pos) + max_allele_len);

    if kmer_start_pos + flank_len + max_allele_len + needed_right_flank > ref_seq.len() {
        return None;
    }

    let mut ref_kmer = Vec::with_capacity(MARKER_KMER_LEN);
    ref_kmer.extend_from_slice(&ref_seq[kmer_start_pos..pos0]);
    ref_kmer.extend_from_slice(ref_allele);
    ref_kmer.extend_from_slice(&ref_seq[right_flank_start_ref..right_flank_start_ref + needed_right_flank]);
    ref_kmer.resize(MARKER_KMER_LEN, b'N');

    let mut alt_kmer = Vec::with_capacity(MARKER_KMER_LEN);
    alt_kmer.extend_from_slice(&ref_seq[kmer_start_pos..pos0]);
    alt_kmer.extend_from_slice(alt_allele);
    alt_kmer.extend_from_slice(&ref_seq[right_flank_start_ref..right_flank_start_ref + needed_right_flank]);
    alt_kmer.resize(MARKER_KMER_LEN, b'N');

    Some((ref_kmer, alt_kmer))
}

/// Builds diagnostic k-mer(s) for a large structural variant (insertion or deletion).
fn build_large_variant_kmers(
    pos0: usize,
    ref_allele: &[u8],
    alt_allele: &[u8],
    ref_seq: &[u8],
    base_marker: &Marker,
) -> Vec<Marker> {
    warn!(
        "Handling large variant at pos {} as a structural variant.",
        pos0 + 1
    );
    let k_half = MARKER_KMER_LEN / 2;
    let mut sv_markers = Vec::new();

    if ref_allele.len() > alt_allele.len() {
        // Large Deletion
        if pos0 >= k_half {
            let right_flank_start = pos0 + ref_allele.len();
            if right_flank_start + k_half <= ref_seq.len() {
                let mut junction_kmer = Vec::with_capacity(MARKER_KMER_LEN);
                junction_kmer.extend_from_slice(&ref_seq[pos0 - k_half..pos0]);
                junction_kmer.extend_from_slice(alt_allele);
                let needed_right_len = MARKER_KMER_LEN.saturating_sub(junction_kmer.len());
                if right_flank_start + needed_right_len <= ref_seq.len() {
                    junction_kmer.extend_from_slice(
                        &ref_seq[right_flank_start..right_flank_start + needed_right_len],
                    );
                    junction_kmer.resize(MARKER_KMER_LEN, b'N');
                    sv_markers.push(Marker {
                        alt_kmer: junction_kmer,
                        ref_kmer: vec![b'N'; MARKER_KMER_LEN], // Ref kmer is meaningless for a junction
                        ..base_marker.clone()
                    });
                }
            }
        }
    } else {
        // Large Insertion
        // Left Junction
        if pos0 >= k_half {
            let mut left_junction_kmer = Vec::with_capacity(MARKER_KMER_LEN);
            left_junction_kmer.extend_from_slice(&ref_seq[pos0 - k_half..pos0]);
            left_junction_kmer.extend_from_slice(&alt_allele[0..MARKER_KMER_LEN - k_half]);
            sv_markers.push(Marker {
                alt_kmer: left_junction_kmer,
                ref_kmer: vec![b'N'; MARKER_KMER_LEN],
                ..base_marker.clone()
            });
        }
        // Right Junction
        let right_flank_start = pos0 + ref_allele.len();
        if right_flank_start + k_half <= ref_seq.len() {
            let mut right_junction_kmer = Vec::with_capacity(MARKER_KMER_LEN);
            right_junction_kmer
                .extend_from_slice(&alt_allele[alt_allele.len() - (MARKER_KMER_LEN - k_half)..]);
            right_junction_kmer
                .extend_from_slice(&ref_seq[right_flank_start..right_flank_start + k_half]);
            sv_markers.push(Marker {
                alt_kmer: right_junction_kmer,
                ref_kmer: vec![b'N'; MARKER_KMER_LEN],
                ..base_marker.clone()
            });
        }
    }
    sv_markers
}

// --- Public API ---

pub fn build_markers(ref_fasta: &str, tsv_markers: &str) -> Result<Vec<Marker>> {
    info!("  Reading reference genome...");
    let ref_seq = read_ref_sequence(ref_fasta)?;

    info!("  Parsing marker file and generating diagnostic k-mers...");
    let mut markers = Vec::new();
    let file = File::open(tsv_markers)?;
    for line in BufReader::new(file).lines() {
        let line_str = line?;
        if line_str.trim().is_empty() || line_str.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line_str.split('\t').collect();
        if fields.len() < 4 {
            warn!("Skipping marker line with fewer than 4 columns: {}", line_str);
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
            warn!("Skipping marker at pos {} due to no lineage information.", pos0 + 1);
            continue;
        }

        let ref_allele = ref_allele_str.as_bytes();
        let alt_allele = alt_allele_str.as_bytes();

        let max_allele_len = ref_allele.len().max(alt_allele.len());

        if max_allele_len < MARKER_KMER_LEN {
            // --- Logic for Small Variants (SNPs, Indels, MNVs) ---
            if let Some((ref_kmer, alt_kmer)) =
                build_small_variant_kmer(pos0, ref_allele, alt_allele, &ref_seq)
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
            // --- Logic for Large Variants (SVs) ---
            let base_marker = Marker {
                pos: pos0,
                lineages: lineage_cols,
                ref_allele: ref_allele_str.clone(),
                alt_allele: alt_allele_str.clone(),
                ref_kmer: Vec::new(), // Placeholder
                alt_kmer: Vec::new(), // Placeholder
                annotations: annotation_cols,
            };
            let sv_markers =
                build_large_variant_kmers(pos0, ref_allele, alt_allele, &ref_seq, &base_marker);
            markers.extend(sv_markers);
        }
    }
    Ok(markers)
}

pub fn scan_fastq(fastqs: &[String], markers: &[Marker]) -> Result<Vec<[u32; 2]>> {
    info!(
        "  Indexing {} markers (including reverse complements)...",
        markers.len()
    );
    let mut index: FxHashMap<u64, Vec<(usize, bool)>> = FxHashMap::default();
    for (id, marker) in markers.iter().enumerate() {
        if let Some(ref_hash) = hash_kmer(&marker.ref_kmer) {
            index.entry(ref_hash).or_default().push((id, false));
        }
        if let Some(alt_hash) = hash_kmer(&marker.alt_kmer) {
            index.entry(alt_hash).or_default().push((id, true));
        }
        let rc_ref = reverse_complement(&marker.ref_kmer);
        if let Some(rc_ref_hash) = hash_kmer(&rc_ref) {
            index.entry(rc_ref_hash).or_default().push((id, false));
        }
        let rc_alt = reverse_complement(&marker.alt_kmer);
        if let Some(rc_alt_hash) = hash_kmer(&rc_alt) {
            index.entry(rc_alt_hash).or_default().push((id, true));
        }
    }
    let counts = vec![[0u32; 2]; markers.len()];
    let counts_mutex = std::sync::Mutex::new(counts);
    info!("  Scanning FASTQ reads in parallel...");
    fastqs.par_iter().for_each(|path| {
        if let Ok(mut rdr) = needletail::parse_fastx_file(path) {
            while let Some(Ok(record)) = rdr.next() {
                for kmer_slice in kmer::Kmers::new(record.sequence(), MARKER_KMER_LEN as u8) {
                    if let Some(kmer_hash) = hash_kmer(kmer_slice) {
                        if let Some(hits) = index.get(&kmer_hash) {
                            let mut guard = counts_mutex.lock().unwrap();
                            for &(marker_id, is_alt) in hits {
                                if is_alt {
                                    guard[marker_id][1] += 1;
                                } else {
                                    guard[marker_id][0] += 1;
                                }
                            }
                        }
                    }
                }
            }
        } else {
            warn!("Could not open or parse FASTQ file: {}. Skipping.", path);
        }
    });
    Ok(counts_mutex.into_inner().unwrap())
}
