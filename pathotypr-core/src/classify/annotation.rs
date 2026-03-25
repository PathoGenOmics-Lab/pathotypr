//! GFF annotation parsing and variant-to-amino-acid translation.
//!
//! Handles CDS feature extraction from GFF3 files and translates
//! genomic variants (SNPs, MNVs) to amino acid changes.

use log::{debug, warn};
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::OnceLock;

use crate::errors::AppResult;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Gene {
    pub id: String,
    pub strand: Strand,
}

// ---------------------------------------------------------------------------
// Genetic code
// ---------------------------------------------------------------------------

fn genetic_code_3_letter() -> &'static HashMap<&'static str, &'static str> {
    static CODE: OnceLock<HashMap<&'static str, &'static str>> = OnceLock::new();
    CODE.get_or_init(|| {
        [
            ("GCA", "Ala"), ("GCC", "Ala"), ("GCG", "Ala"), ("GCT", "Ala"),
            ("CGA", "Arg"), ("CGC", "Arg"), ("CGG", "Arg"), ("CGT", "Arg"),
            ("AGA", "Arg"), ("AGG", "Arg"),
            ("AAC", "Asn"), ("AAT", "Asn"),
            ("GAC", "Asp"), ("GAT", "Asp"),
            ("TGC", "Cys"), ("TGT", "Cys"),
            ("GAA", "Glu"), ("GAG", "Glu"),
            ("CAA", "Gln"), ("CAG", "Gln"),
            ("GGA", "Gly"), ("GGC", "Gly"), ("GGG", "Gly"), ("GGT", "Gly"),
            ("CAC", "His"), ("CAT", "His"),
            ("ATA", "Ile"), ("ATC", "Ile"), ("ATT", "Ile"),
            ("CTA", "Leu"), ("CTC", "Leu"), ("CTG", "Leu"), ("CTT", "Leu"),
            ("TTA", "Leu"), ("TTG", "Leu"),
            ("AAA", "Lys"), ("AAG", "Lys"),
            ("ATG", "Met"),
            ("TTC", "Phe"), ("TTT", "Phe"),
            ("CCA", "Pro"), ("CCC", "Pro"), ("CCG", "Pro"), ("CCT", "Pro"),
            ("TCA", "Ser"), ("TCC", "Ser"), ("TCG", "Ser"), ("TCT", "Ser"),
            ("AGC", "Ser"), ("AGT", "Ser"),
            ("ACA", "Thr"), ("ACC", "Thr"), ("ACG", "Thr"), ("ACT", "Thr"),
            ("TGG", "Trp"),
            ("TAC", "Tyr"), ("TAT", "Tyr"),
            ("GTA", "Val"), ("GTC", "Val"), ("GTG", "Val"), ("GTT", "Val"),
            ("TAA", "Stp"), ("TAG", "Stp"), ("TGA", "Stp"),
        ]
        .iter()
        .cloned()
        .collect()
    })
}

// ---------------------------------------------------------------------------
// GFF parsing
// ---------------------------------------------------------------------------

/// Parse a GFF3 file and build an interval tree of CDS features.
pub fn parse_gff_and_build_tree(gff_file: &str) -> AppResult<Lapper<usize, Gene>> {
    debug!("Parsing GFF file: {}", gff_file);
    let file = File::open(gff_file)?;
    let reader = BufReader::new(file);
    let mut intervals = Vec::new();

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                warn!("Skipping unreadable GFF line: {}", e);
                continue;
            }
        };
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut fields = line.splitn(9, '\t');
        let _seqid = match fields.next() { Some(f) => f, None => continue };
        let _source = match fields.next() { Some(f) => f, None => continue };
        let feature_type = match fields.next() { Some(f) => f, None => continue };

        if feature_type != "CDS" {
            continue;
        }

        let start_str = match fields.next() { Some(f) => f, None => continue };
        let end_str = match fields.next() { Some(f) => f, None => continue };
        let _score = match fields.next() { Some(f) => f, None => continue };
        let strand_str = match fields.next() { Some(f) => f, None => continue };
        let _phase = match fields.next() { Some(f) => f, None => continue };
        let attributes_str = match fields.next() { Some(f) => f, None => continue };

        let start: usize = match start_str.parse() {
            Ok(v) => v,
            Err(_) => { warn!("Skipping CDS with invalid start: {}", start_str); continue; }
        };
        let end: usize = match end_str.parse() {
            Ok(v) => v,
            Err(_) => { warn!("Skipping CDS with invalid end: {}", end_str); continue; }
        };

        let strand = match strand_str {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::Unknown,
        };

        let mut gene_name: Option<&str> = None;
        for attr in attributes_str.split(';') {
            let attr = attr.trim();
            if let Some((key, val)) = attr.split_once('=') {
                match key {
                    "gene" => { gene_name = Some(val); break; }
                    "locus_tag" if gene_name.is_none() => { gene_name = Some(val); }
                    "Name" if gene_name.is_none() => { gene_name = Some(val); }
                    "ID" if gene_name.is_none() => { gene_name = Some(val); }
                    _ => {}
                }
            }
        }

        let gene = Gene {
            id: gene_name.unwrap_or("Unknown").to_string(),
            strand,
        };
        intervals.push(Interval {
            start: start - 1,
            stop: end,
            val: gene,
        });
    }

    let num_intervals = intervals.len();
    let lapper = Lapper::new(intervals);
    debug!("Successfully parsed {} CDS features from GFF.", num_intervals);
    Ok(lapper)
}

// ---------------------------------------------------------------------------
// Complement
// ---------------------------------------------------------------------------

pub(crate) fn complement_byte(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => b'N',
    }
}

// ---------------------------------------------------------------------------
// Variant translation
// ---------------------------------------------------------------------------

/// Translate a variant (SNP or MNV) to amino acid change(s).
/// Handles multi-base substitutions that may span codon boundaries.
pub(crate) fn translate_variant_info(
    gene_interval: &Interval<usize, Gene>,
    ref_pos_0based: usize,
    alt_bases: &str,
    ref_allele_len: usize,
    ref_seq: &str,
    ref_seq_rc: &str,
) -> (Option<String>, Option<String>, Option<String>) {
    let gene = &gene_interval.val;
    let code = genetic_code_3_letter();

    // For indel-like MNVs (ref_len != alt_len), report as frameshift
    if ref_allele_len != alt_bases.len() {
        let first_aa_pos = match gene.strand {
            Strand::Forward => {
                if ref_pos_0based < gene_interval.start { return (None, None, None); }
                (ref_pos_0based - gene_interval.start) / 3 + 1
            }
            Strand::Reverse => {
                let last_ref_pos = ref_pos_0based + ref_allele_len - 1;
                let pos_rc = ref_seq.len() - last_ref_pos - 1;
                let start_rc = ref_seq.len() - gene_interval.stop;
                if pos_rc < start_rc { return (None, None, None); }
                (pos_rc - start_rc) / 3 + 1
            }
            _ => return (None, None, None),
        };
        return (
            Some(gene.id.clone()),
            Some(first_aa_pos.to_string()),
            Some(format!("frameshift_pos{}", first_aa_pos)),
        );
    }

    let variant_len = ref_allele_len;

    match gene.strand {
        Strand::Forward => {
            if ref_pos_0based < gene_interval.start {
                return (None, None, None);
            }
            let offset_in_gene = ref_pos_0based - gene_interval.start;
            let first_codon_idx = offset_in_gene / 3;
            let last_affected_pos = ref_pos_0based + variant_len - 1;
            let clamped_last = last_affected_pos.min(gene_interval.stop.saturating_sub(1));
            let last_offset = clamped_last - gene_interval.start;
            let last_codon_idx = last_offset / 3;

            let ref_bytes = ref_seq.as_bytes();
            let alt_bytes = alt_bases.as_bytes();
            let mut changes = Vec::new();

            for codon_idx in first_codon_idx..=last_codon_idx {
                let codon_start = gene_interval.start + codon_idx * 3;
                if codon_start + 3 > ref_bytes.len() { break; }

                let ref_codon: [u8; 3] = [ref_bytes[codon_start], ref_bytes[codon_start + 1], ref_bytes[codon_start + 2]];
                let mut alt_codon = ref_codon;

                for i in 0..3 {
                    let genome_pos = codon_start + i;
                    if genome_pos >= ref_pos_0based && genome_pos < ref_pos_0based + variant_len {
                        let alt_idx = genome_pos - ref_pos_0based;
                        alt_codon[i] = alt_bytes[alt_idx];
                    }
                }

                let ref_codon_str = String::from_utf8_lossy(&ref_codon);
                let alt_codon_str = String::from_utf8_lossy(&alt_codon);
                let ref_aa = code.get(ref_codon_str.as_ref()).unwrap_or(&"???");
                let alt_aa = code.get(alt_codon_str.as_ref()).unwrap_or(&"???");
                let aa_pos = codon_idx + 1;
                let syn = if ref_aa == alt_aa { " syn" } else { "" };
                changes.push((aa_pos, format!("{}{}{}({}>{}){}",
                    ref_aa, aa_pos, alt_aa, ref_codon_str, alt_codon_str, syn)));
            }

            if changes.is_empty() {
                return (None, None, None);
            }
            let first_aa = changes[0].0;
            let change_str = changes.into_iter().map(|(_, s)| s).collect::<Vec<_>>().join(";");
            (Some(gene.id.clone()), Some(first_aa.to_string()), Some(change_str))
        }
        Strand::Reverse => {
            let last_ref_pos = ref_pos_0based + variant_len - 1;
            let start_rc = ref_seq.len() - gene_interval.stop;
            let first_pos_rc = ref_seq.len() - last_ref_pos - 1;

            if first_pos_rc < start_rc {
                return (None, None, None);
            }

            let rc_bytes = ref_seq_rc.as_bytes();
            let offset_in_gene_rc = first_pos_rc - start_rc;
            let first_codon_idx = offset_in_gene_rc / 3;
            let last_pos_rc = ref_seq.len() - ref_pos_0based - 1;
            let end_rc = ref_seq.len() - gene_interval.start - 1;
            let clamped_last_rc = last_pos_rc.min(end_rc);
            let last_offset_rc = clamped_last_rc - start_rc;
            let last_codon_idx = last_offset_rc / 3;

            let alt_bytes_orig = alt_bases.as_bytes();
            let alt_rc: Vec<u8> = alt_bytes_orig.iter().rev().map(|&b| complement_byte(b)).collect();

            let mut changes = Vec::new();

            for codon_idx in first_codon_idx..=last_codon_idx {
                let codon_start = start_rc + codon_idx * 3;
                if codon_start + 3 > rc_bytes.len() { break; }

                let ref_codon: [u8; 3] = [rc_bytes[codon_start], rc_bytes[codon_start + 1], rc_bytes[codon_start + 2]];
                let mut alt_codon = ref_codon;

                for i in 0..3 {
                    let rc_pos = codon_start + i;
                    if rc_pos >= first_pos_rc && rc_pos <= last_pos_rc {
                        let alt_idx = rc_pos - first_pos_rc;
                        if alt_idx < alt_rc.len() {
                            alt_codon[i] = alt_rc[alt_idx];
                        }
                    }
                }

                let ref_codon_str = String::from_utf8_lossy(&ref_codon);
                let alt_codon_str = String::from_utf8_lossy(&alt_codon);
                let ref_aa = code.get(ref_codon_str.as_ref()).unwrap_or(&"???");
                let alt_aa = code.get(alt_codon_str.as_ref()).unwrap_or(&"???");
                let aa_pos = codon_idx + 1;
                let syn = if ref_aa == alt_aa { " syn" } else { "" };
                changes.push((aa_pos, format!("{}{}{}({}>{}){}",
                    ref_aa, aa_pos, alt_aa, ref_codon_str, alt_codon_str, syn)));
            }

            if changes.is_empty() {
                return (None, None, None);
            }
            let first_aa = changes[0].0;
            let change_str = changes.into_iter().map(|(_, s)| s).collect::<Vec<_>>().join(";");
            (Some(gene.id.clone()), Some(first_aa.to_string()), Some(change_str))
        }
        _ => (None, None, None),
    }
}
