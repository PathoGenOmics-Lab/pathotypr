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
#[derive(Debug, Clone)]
pub struct Marker {
    pub pos: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub lineage: String,
    pub alt_kmer: Vec<u8>,
    pub ref_kmer: Vec<u8>,
    pub annotations: Vec<String>,
}

fn hash_kmer(kmer: &[u8]) -> Option<u64> {
    let mut val: u64 = 0;
    for &base in kmer { val = (val << 2) | match base { b'A' | b'a' => 0, b'C' | b'c' => 1, b'G' | b'g' => 2, b'T' | b't' => 3, _ => return None }; }
    Some(val)
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut rc_seq: Vec<u8> = Vec::with_capacity(seq.len());
    for &base in seq.iter().rev() { rc_seq.push(match base { b'A' | b'a' => b'T', b'C' | b'c' => b'G', b'G' | b'g' => b'C', b'T' | b't' => b'A', _ => b'N' }); }
    rc_seq
}

fn read_ref_sequence(path: &str) -> Result<Vec<u8>> {
    let mut reader = fasta::Reader::from_file(path)?;
    let mut record = fasta::Record::new();
    reader.read(&mut record)?;
    Ok(record.seq().to_ascii_uppercase())
}

pub fn build_markers(ref_fasta: &str, tsv_markers: &str) -> Result<Vec<Marker>> {
    info!("  Reading reference genome...");
    let ref_seq = read_ref_sequence(ref_fasta)?;
    let ref_len = ref_seq.len();

    info!("  Parsing marker file and generating diagnostic k-mers...");
    let mut markers = Vec::new();
    let file = File::open(tsv_markers)?;
    for line in BufReader::new(file).lines() {
        let line_str = line?;
        if line_str.trim().is_empty() || line_str.starts_with('#') { continue; }
        let fields: Vec<&str> = line_str.split('\t').collect();
        if fields.len() < 4 { continue; }

        let pos0: usize = match fields[0].parse::<usize>() { Ok(p) => p.saturating_sub(1), Err(_) => { warn!("Skipping marker with invalid position: {}", fields[0]); continue; } };
        let ref_allele_str = fields[1].to_string();
        let alt_allele_str = fields[2].to_string();
        let lineage = fields[3].to_string();
        let annotations = if fields.len() > 4 { fields[4..].iter().map(|s| s.to_string()).collect() } else { Vec::new() };
        
        let ref_allele = ref_allele_str.as_bytes();
        let alt_allele = alt_allele_str.as_bytes();
        let ref_len_allele = ref_allele.len();
        let alt_len_allele = alt_allele.len();

        let max_allele_len = ref_len_allele.max(alt_len_allele);

        if max_allele_len < MARKER_KMER_LEN {
            // --- Logic for Small Variants (SNPs, Indels, MNVs) ---
            let flank_len = (MARKER_KMER_LEN - max_allele_len) / 2;
            let kmer_start_pos = pos0.saturating_sub(flank_len);
            let right_flank_start_ref = pos0 + ref_len_allele;
            let needed_right_flank = MARKER_KMER_LEN.saturating_sub(kmer_start_pos.saturating_sub(pos0) + max_allele_len);
            if kmer_start_pos + flank_len + max_allele_len + needed_right_flank > ref_len { continue; }

            let mut ref_kmer = Vec::with_capacity(MARKER_KMER_LEN);
            ref_kmer.extend_from_slice(&ref_seq[kmer_start_pos..pos0]);
            ref_kmer.extend_from_slice(ref_allele);
            ref_kmer.extend_from_slice(&ref_seq[right_flank_start_ref .. right_flank_start_ref + needed_right_flank]);
            ref_kmer.resize(MARKER_KMER_LEN, b'N');

            let mut alt_kmer = Vec::with_capacity(MARKER_KMER_LEN);
            alt_kmer.extend_from_slice(&ref_seq[kmer_start_pos..pos0]);
            alt_kmer.extend_from_slice(alt_allele);
            alt_kmer.extend_from_slice(&ref_seq[right_flank_start_ref .. right_flank_start_ref + needed_right_flank]);
            alt_kmer.resize(MARKER_KMER_LEN, b'N');

            markers.push(Marker { pos: pos0, lineage, ref_allele: ref_allele_str, alt_allele: alt_allele_str, ref_kmer, alt_kmer, annotations });
        } else {
            // --- NEW LOGIC for Large Variants (SVs) ---
            warn!("Handling large variant at pos {} as a structural variant.", pos0 + 1);
            let k_half = MARKER_KMER_LEN / 2;
            let dummy_ref_kmer = vec![b'N'; MARKER_KMER_LEN];

            if ref_len_allele > alt_len_allele { // Large Deletion or complex delins
                let mut junction_kmer = Vec::with_capacity(MARKER_KMER_LEN);
                if pos0 < k_half { continue; }
                let right_flank_start = pos0 + ref_len_allele;
                if right_flank_start + k_half > ref_len { continue; }
                
                junction_kmer.extend_from_slice(&ref_seq[pos0 - k_half .. pos0]);
                junction_kmer.extend_from_slice(alt_allele);
                let needed_right_len = MARKER_KMER_LEN.saturating_sub(junction_kmer.len());
                if right_flank_start + needed_right_len > ref_len { continue; }
                junction_kmer.extend_from_slice(&ref_seq[right_flank_start .. right_flank_start + needed_right_len]);
                junction_kmer.resize(MARKER_KMER_LEN, b'N');

                markers.push(Marker { pos: pos0, lineage, ref_allele: ref_allele_str, alt_allele: alt_allele_str, ref_kmer: dummy_ref_kmer, alt_kmer: junction_kmer, annotations });
            } else { // Large Insertion
                // Left Junction
                if pos0 >= k_half {
                    let mut left_junction_kmer = Vec::with_capacity(MARKER_KMER_LEN);
                    left_junction_kmer.extend_from_slice(&ref_seq[pos0 - k_half .. pos0]);
                    left_junction_kmer.extend_from_slice(&alt_allele[0..MARKER_KMER_LEN - k_half]);
                    markers.push(Marker { pos: pos0, lineage: lineage.clone(), ref_allele: ref_allele_str.clone(), alt_allele: alt_allele_str.clone(), ref_kmer: dummy_ref_kmer.clone(), alt_kmer: left_junction_kmer, annotations: annotations.clone() });
                }
                // Right Junction
                let right_flank_start = pos0 + ref_len_allele;
                if right_flank_start + k_half <= ref_len {
                    let mut right_junction_kmer = Vec::with_capacity(MARKER_KMER_LEN);
                    right_junction_kmer.extend_from_slice(&alt_allele[alt_len_allele - (MARKER_KMER_LEN - k_half) ..]);
                    right_junction_kmer.extend_from_slice(&ref_seq[right_flank_start .. right_flank_start + k_half]);
                    markers.push(Marker { pos: pos0, lineage, ref_allele: ref_allele_str, alt_allele: alt_allele_str, ref_kmer: dummy_ref_kmer, alt_kmer: right_junction_kmer, annotations });
                }
            }
        }
    }
    Ok(markers)
}

pub fn scan_fastq(fastqs: &[String], markers: &[Marker]) -> Result<Vec<[u32; 2]>> {
    info!("  Indexing {} markers (including reverse complements)...", markers.len());
    let mut index: FxHashMap<u64, Vec<(usize, bool)>> = FxHashMap::default();
    for (id, marker) in markers.iter().enumerate() {
        if let Some(ref_hash) = hash_kmer(&marker.ref_kmer) { index.entry(ref_hash).or_default().push((id, false)); }
        if let Some(alt_hash) = hash_kmer(&marker.alt_kmer) { index.entry(alt_hash).or_default().push((id, true)); }
        let rc_ref = reverse_complement(&marker.ref_kmer);
        if let Some(rc_ref_hash) = hash_kmer(&rc_ref) { index.entry(rc_ref_hash).or_default().push((id, false)); }
        let rc_alt = reverse_complement(&marker.alt_kmer);
        if let Some(rc_alt_hash) = hash_kmer(&rc_alt) { index.entry(rc_alt_hash).or_default().push((id, true)); }
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
                                if is_alt { guard[marker_id][1] += 1; } 
                                else { guard[marker_id][0] += 1; }
                            }
                        }
                    }
                }
            }
        } else { warn!("Could not open or parse FASTQ file: {}. Skipping.", path); }
    });
    Ok(counts_mutex.into_inner().unwrap())
}
