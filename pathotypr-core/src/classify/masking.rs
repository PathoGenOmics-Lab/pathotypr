//! FASTA sequence masking at marker positions.
//!
//! Replaces marker-affected positions with 'N' for phylogenetic analysis
//! that needs to exclude DR-related positions.

use log::info;
use needletail::Sequence;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};

use super::markers::MarkerVariant;
use crate::errors::{AppError, AppResult};

/// Collect all unique genomic positions (0-based ranges) that should be masked.
pub fn collect_mask_positions(markers: &[MarkerVariant]) -> Vec<(usize, usize)> {
    let mut positions: HashSet<(usize, usize)> = HashSet::new();
    for m in markers {
        let start = m.pos.saturating_sub(1);
        let end = start + m.ref_allele.len();
        positions.insert((start, end));
    }
    let mut sorted: Vec<(usize, usize)> = positions.into_iter().collect();
    sorted.sort_unstable();
    // Merge overlapping ranges
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (s, e) in sorted {
        if let Some(last) = merged.last_mut() {
            if s <= last.1 {
                last.1 = last.1.max(e);
                continue;
            }
        }
        merged.push((s, e));
    }
    merged
}

/// Apply mask positions to a sequence, replacing masked bases with 'N'.
pub fn mask_sequence(seq: &[u8], mask_ranges: &[(usize, usize)]) -> Vec<u8> {
    let mut masked = seq.to_vec();
    for &(start, end) in mask_ranges {
        let end = end.min(masked.len());
        if start < masked.len() {
            for b in &mut masked[start..end] {
                *b = b'N';
            }
        }
    }
    masked
}

/// Write a masked multi-FASTA file from input genomes.
pub fn write_masked_fasta(
    fasta_path: &str,
    output_path: &str,
    mask_ranges: &[(usize, usize)],
) -> AppResult<()> {
    let mut reader = needletail::parse_fastx_file(fasta_path).map_err(|e| {
        AppError::Generic(format!("Failed to open FASTA file {}: {}", fasta_path, e))
    })?;
    let mut out = BufWriter::new(File::create(output_path)?);
    while let Some(record) = reader.next() {
        let record = record.map_err(|e| {
            AppError::Generic(format!("Failed to read FASTA record: {}", e))
        })?;
        let id = String::from_utf8_lossy(record.id());
        let seq = record.normalize(true);
        let masked = mask_sequence(&seq, mask_ranges);
        writeln!(out, ">{}", id)?;
        for chunk in masked.chunks(80) {
            out.write_all(chunk)?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    info!("Masked FASTA written to: {}", output_path);
    Ok(())
}
