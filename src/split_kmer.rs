//! Implements a split k-mer analysis engine.
//!
//! This technique is used to identify single nucleotide polymorphisms (SNPs)
//! directly from raw sequencing reads without alignment. A marker is defined
//! by a 'left' k-mer, a 'right' k-mer, and the central base that differs.
//! This module provides the logic to build these markers and efficiently
//! scan FASTQ files to count them.

use anyhow::{anyhow, Context}; // Import anyhow features
use fxhash::FxHashMap;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

/// The length of the left and right arms of the split k-mer.
/// A total k-mer of size `2 * K + 1` is scanned.
pub const K: usize = 15;

/// Bitmask to ensure encoded k-mers fit within `2 * K` bits.
const MASK: u32 = (1u32 << (2 * K)) - 1;

/// Encodes a k-mer sequence into a compact `u32` representation (2 bits per base).
///
/// Returns `None` if the k-mer contains invalid characters (e.g., 'N').
fn encode(kmer: &[u8]) -> Option<u32> {
    if kmer.len() != K {
        return None;
    }
    let mut val: u32 = 0;
    for &base in kmer {
        val = (val << 2)
            | match base {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => return None, // Contains 'N' or other non-ACGT base
            };
    }
    Some(val & MASK)
}

/// Represents a single split k-mer marker.
#[derive(Clone, Debug)]
pub struct Marker {
    /// The encoded left part of the k-mer.
    pub left: u32,
    /// The encoded right part of the k-mer.
    pub right: u32,
    /// The base at this position in the reference genome.
    pub ref_base: u8,
    /// The alternate base that defines the marker.
    pub alt_base: u8,
    /// The 0-based position of the SNP in the reference genome.
    pub pos: usize,
    /// The lineage associated with this marker.
    pub lineage: String,
}

/// Constructs a list of `Marker` objects from a reference FASTA and a TSV file.
///
/// # Arguments
/// * `ref_fasta` - Path to the reference genome FASTA file.
/// * `tsv_markers` - Path to the TSV file defining markers (pos, alt, lineage).
pub fn build_markers(ref_fasta: &str, tsv_markers: &str) -> anyhow::Result<Vec<Marker>> {
    // Read the reference sequence.
    let mut rdr = parse_fastx_file(ref_fasta)
        .with_context(|| format!("Failed to open reference FASTA: {}", ref_fasta))?;
    let rec = rdr
        .next()
        .ok_or_else(|| anyhow!("Reference FASTA file is empty"))??;
    let ref_seq = rec.seq().to_vec();
    let ref_len = ref_seq.len();

    let mut markers = Vec::new();
    let file = File::open(tsv_markers)
        .with_context(|| format!("Failed to open markers TSV: {}", tsv_markers))?;
    for line in BufReader::new(file).lines() {
        let line_str = line?;
        if line_str.trim().is_empty() || line_str.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line_str.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        // Parse 1-based position from the file.
        let pos1: usize = match fields[0].parse() {
            Ok(v) => v,
            Err(_) => {
                log::warn!("Skipping non-numeric line in markers TSV: {}", fields[0]);
                continue;
            }
        };

        let pos0 = pos1 - 1; // Convert to 0-based index.
        if pos0 < K || pos0 + K >= ref_len {
            log::warn!(
                "Marker position {} is too close to the sequence ends. Skipping.",
                pos1
            );
            continue;
        }

        // Extract and encode the left and right k-mer arms, returning a proper error on failure.
        let left_kmer = encode(&ref_seq[pos0 - K..pos0])
            .ok_or_else(|| anyhow!("Invalid base found in left k-mer flank for marker at position {}", pos1))?;
        let right_kmer = encode(&ref_seq[pos0 + 1..pos0 + 1 + K])
            .ok_or_else(|| anyhow!("Invalid base found in right k-mer flank for marker at position {}", pos1))?;

        markers.push(Marker {
            left: left_kmer,
            right: right_kmer,
            ref_base: ref_seq[pos0],
            alt_base: fields[1].as_bytes()[0],
            pos: pos0,
            lineage: fields[2].to_string(),
        });
    }
    Ok(markers)
}

/// Scans FASTQ files to count occurrences of split k-mer markers.
///
/// This is the performance-critical part of the analysis. It uses `rayon` for
/// parallelism and `FxHashMap` for fast lookups.
pub fn scan_fastq(fastqs: &[String], markers: &[Marker]) -> anyhow::Result<Vec<[u32; 4]>> {
    // Create an index mapping k-mer hashes to the markers they belong to.
    // The boolean indicates if it's a left arm (true) or right arm (false).
    let mut index: FxHashMap<u32, Vec<(usize, bool)>> = FxHashMap::default();
    for (id, marker) in markers.iter().enumerate() {
        index.entry(marker.left).or_default().push((id, true));
        index.entry(marker.right).or_default().push((id, false));
    }

    // A thread-safe structure to hold the base counts (A,C,G,T) for each marker.
    let counts = vec![[0u32; 4]; markers.len()];
    let counts_mutex = std::sync::Mutex::new(counts);

    // Process each FASTQ file in parallel.
    fastqs.par_iter().try_for_each(|path| -> anyhow::Result<()> {
        let mut reader = parse_fastx_file(path)
            .with_context(|| format!("Failed to open FASTQ file: {}", path))?;
        while let Some(record) = reader.next() {
            let rec = record.with_context(|| format!("Failed to parse record in {}", path))?;
            let seq = rec.seq();
            let seq_len = seq.len();
            if seq_len < 2 * K + 1 {
                continue;
            }

            // Pre-calculate hashes for all K-sized windows in the read.
            let mut hashes = Vec::with_capacity(seq_len - K + 1);
            for window in seq.windows(K) {
                hashes.push(encode(window).unwrap_or(u32::MAX));
            }

            // Slide a (2*K+1)-sized window across the read to find marker pairs.
            for i in 0..=seq_len - (2 * K + 1) {
                let left_hash = hashes[i];
                let right_hash = hashes[i + K + 1];
                if left_hash == u32::MAX || right_hash == u32::MAX {
                    continue;
                }

                // Check if the left and right hashes exist in our marker index.
                if let (Some(left_hits), Some(right_hits)) =
                    (index.get(&left_hash), index.get(&right_hash))
                {
                    // Iterate through potential matches.
                    for &(id1, is_left1) in left_hits {
                        if !is_left1 {
                            continue;
                        } // Must be a left arm.
                        for &(id2, is_left2) in right_hits {
                            // If they belong to the same marker and form a valid pair...
                            if id1 == id2 && !is_left2 {
                                let central_base = seq[i + K];
                                let mut guard = counts_mutex.lock().unwrap();
                                // ...increment the count for the observed central base.
                                match central_base {
                                    b'A' | b'a' => guard[id1][0] += 1,
                                    b'C' | b'c' => guard[id1][1] += 1,
                                    b'G' | b'g' => guard[id1][2] += 1,
                                    b'T' | b't' => guard[id1][3] += 1,
                                    _ => {}
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    })?;

    // Return the final counts.
    Ok(counts_mutex.into_inner().unwrap())
}
