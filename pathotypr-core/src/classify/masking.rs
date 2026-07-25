//! FASTA sequence masking at marker positions.
//!
//! Replaces marker-affected positions with 'N' for phylogenetic analysis
//! that needs to exclude DR-related positions.

use log::{info, warn};
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
///
/// The mask ranges are reference coordinates, so they may only be applied to
/// sequences that live in the reference coordinate system — a reference-based
/// consensus genome, or an alignment whose records all have the reference
/// length. The contigs of a de-novo assembly bear no relation to those
/// coordinates, and masking one would silently replace arbitrary positions with
/// `N`, quietly corrupting a file that is typically fed straight into a
/// phylogeny.
///
/// Returns `Ok(false)` (writing nothing) when the input is not in reference
/// coordinates, so the caller can report it rather than emit a wrong file.
pub fn write_masked_fasta(
    fasta_path: &str,
    output_path: &str,
    mask_ranges: &[(usize, usize)],
    ref_len: usize,
) -> AppResult<bool> {
    // First pass: confirm every record is in the reference coordinate system.
    {
        let mut probe = needletail::parse_fastx_file(fasta_path).map_err(|e| {
            AppError::Generic(format!("Failed to open FASTA file {}: {}", fasta_path, e))
        })?;
        let mut n_records = 0usize;
        while let Some(record) = probe.next() {
            let record = record.map_err(|e| {
                AppError::Generic(format!("Failed to read FASTA record: {}", e))
            })?;
            n_records += 1;
            let len = record.seq().len();
            if len != ref_len {
                warn!(
                    "Not writing a masked FASTA for '{}': record '{}' is {} bp but the reference is \
                     {} bp, so the marker positions do not apply to it. Masking expects a \
                     reference-coordinate sequence (e.g. a consensus genome), not a de-novo assembly.",
                    fasta_path,
                    String::from_utf8_lossy(record.id()),
                    len,
                    ref_len
                );
                return Ok(false);
            }
        }
        if n_records == 0 {
            warn!("Not writing a masked FASTA for '{}': it contains no records.", fasta_path);
            return Ok(false);
        }
    }

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
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn write_fasta(dir: &std::path::Path, name: &str, records: &[(&str, &str)]) -> String {
        let path = dir.join(name);
        let mut f = File::create(&path).unwrap();
        for (id, seq) in records {
            writeln!(f, ">{}", id).unwrap();
            writeln!(f, "{}", seq).unwrap();
        }
        path.to_string_lossy().to_string()
    }

    #[test]
    fn masks_a_reference_length_consensus() {
        let dir = tempfile::tempdir().unwrap();
        let reference = "ACGTACGTACGTACGTACGT"; // 20 bp
        let input = write_fasta(dir.path(), "consensus.fasta", &[("s1", reference)]);
        let out = dir.path().join("out.fasta").to_string_lossy().to_string();

        // Mask the 5th base (0-based range 4..5).
        let written = write_masked_fasta(&input, &out, &[(4, 5)], reference.len()).unwrap();
        assert!(written, "a reference-length record must be masked");

        let produced = std::fs::read_to_string(&out).unwrap();
        let seq: String = produced.lines().filter(|l| !l.starts_with('>')).collect();
        assert_eq!(seq.len(), reference.len());
        assert_eq!(&seq[4..5], "N", "the marker position must be masked");
        assert_eq!(&seq[..4], &reference[..4], "other bases must be untouched");
    }

    #[test]
    fn refuses_to_mask_a_de_novo_assembly() {
        let dir = tempfile::tempdir().unwrap();
        let reference = "ACGTACGTACGTACGTACGT"; // 20 bp
        // Contigs of a draft assembly bear no relation to reference coordinates.
        let input = write_fasta(
            dir.path(),
            "assembly.fasta",
            &[("contig_1", "ACGTACGTAC"), ("contig_2", "GGGGTTTTCCCCAAAA")],
        );
        let out = dir.path().join("out.fasta").to_string_lossy().to_string();

        let written = write_masked_fasta(&input, &out, &[(4, 5)], reference.len()).unwrap();
        assert!(!written, "a de-novo assembly must not be masked");
        assert!(
            !std::path::Path::new(&out).exists(),
            "no misleading masked file may be produced"
        );
    }
}
