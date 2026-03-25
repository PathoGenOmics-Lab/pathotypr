//! FASTA file reading utilities.
//!
//! Uses needletail for fast FASTA parsing (~3× faster than line-by-line
//! BufReader for large files). Progress is reported via indicatif when
//! the file size is available.

use indicatif::{ProgressBar, ProgressStyle};
use log::{info, warn};
use needletail::parse_fastx_file;
use std::fs;

use crate::errors::{AppError, AppResult};

// ---------------------------------------------------------------------------
// Header parsing
// ---------------------------------------------------------------------------

/// Extracts the lineage label from a FASTA header.
///
/// Takes the first whitespace-delimited token as the lineage.
pub fn parse_lineage_from_header(header: &str) -> String {
    header
        .split_whitespace()
        .next()
        .unwrap_or(header)
        .to_string()
}

/// Creates a progress bar for file reading, or a hidden spinner if size unknown.
fn make_read_progress(path: &str) -> ProgressBar {
    let file_size = fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    if file_size > 0 {
        let pb = ProgressBar::new(file_size);
        let style = ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")
            .unwrap_or_else(|_| ProgressStyle::default_bar())
            .progress_chars("#>-");
        pb.set_style(style);
        pb
    } else {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .unwrap_or_else(|_| ProgressStyle::default_spinner()),
        );
        pb
    }
}

// ---------------------------------------------------------------------------
// FASTA reader (labeled, for training)
// ---------------------------------------------------------------------------

/// Reads sequences, lineages, and full headers from a labeled FASTA file.
///
/// Returns `(sequences, lineages, headers)`. Entries with empty sequences
/// are skipped with a warning.
pub fn read_fasta(path: &str) -> AppResult<(Vec<String>, Vec<String>, Vec<String>)> {
    let pb = make_read_progress(path);
    let mut reader = parse_fastx_file(path)
        .map_err(|e| AppError::Parsing(format!("Failed to open FASTA {}: {}", path, e)))?;

    let mut sequences = Vec::new();
    let mut lineages = Vec::new();
    let mut headers = Vec::new();
    let mut bytes_read: u64 = 0;

    while let Some(record) = reader.next() {
        let record = record
            .map_err(|e| AppError::Parsing(format!("Error reading FASTA record: {}", e)))?;

        let header = std::str::from_utf8(record.id())
            .map_err(|e| AppError::Parsing(format!("Non-UTF8 header: {}", e)))?
            .to_string();

        let norm_seq = record.seq();
        let seq = String::from_utf8(norm_seq.to_vec())
            .map_err(|e| AppError::Parsing(format!("Non-UTF8 sequence in {}: {}", header, e)))?;

        // Approximate bytes: header line + sequence + newlines
        bytes_read += header.len() as u64 + norm_seq.len() as u64 + 2;
        pb.set_position(bytes_read);

        if seq.is_empty() {
            warn!("Skipping entry with empty sequence: {}", header);
            continue;
        }

        lineages.push(parse_lineage_from_header(&header));
        headers.push(header);
        sequences.push(seq);
    }

    pb.finish_and_clear();
    info!("  Read {} sequences from {}", sequences.len(), path);
    Ok((sequences, lineages, headers))
}

// ---------------------------------------------------------------------------
// General-purpose reader (no lineage parsing)
// ---------------------------------------------------------------------------

/// Reads sequences and their full headers from a FASTA file.
///
/// Returns `Vec<(header, sequence)>`. Entries with empty sequences are
/// skipped with a warning. Used by `predict` and any workflow that doesn't
/// need lineage labels parsed from headers.
pub fn read_fasta_records(path: &str) -> AppResult<Vec<(String, String)>> {
    let pb = make_read_progress(path);
    let mut reader = parse_fastx_file(path)
        .map_err(|e| AppError::Parsing(format!("Failed to open FASTA {}: {}", path, e)))?;

    let mut records = Vec::new();
    let mut bytes_read: u64 = 0;

    while let Some(record) = reader.next() {
        let record = record
            .map_err(|e| AppError::Parsing(format!("Error reading FASTA record: {}", e)))?;

        let header = std::str::from_utf8(record.id())
            .map_err(|e| AppError::Parsing(format!("Non-UTF8 header: {}", e)))?
            .to_string();

        let norm_seq = record.seq();
        let seq = String::from_utf8(norm_seq.to_vec())
            .map_err(|e| AppError::Parsing(format!("Non-UTF8 sequence in {}: {}", header, e)))?;

        bytes_read += header.len() as u64 + norm_seq.len() as u64 + 2;
        pb.set_position(bytes_read);

        if seq.is_empty() {
            warn!("Skipping entry with empty sequence: {}", header);
            continue;
        }

        records.push((header, seq));
    }

    pb.finish_and_clear();
    info!("  Read {} records from {}", records.len(), path);
    Ok(records)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_lineage_preserves_underscores() {
        assert_eq!(parse_lineage_from_header("L4_1_2"), "L4_1_2");
        assert_eq!(parse_lineage_from_header("L4_1_2 extra info"), "L4_1_2");
    }
}
