//! Feature hashing vectorizer for k-mer-based ML.
//!
//! Maps k-mers to fixed-size buckets via bitmask (hashing trick).
//! Stateless — the same hash function is used at train and predict time.

use log::info;
use needletail::Sequence;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Feature hasher
// ---------------------------------------------------------------------------

/// Stateless feature hasher that maps k-mers to fixed-size buckets via bitmask.
#[derive(Serialize, Deserialize, Debug)]
pub struct FeatureHasher {
    /// Number of hash buckets (must be a power of 2).
    pub num_buckets: usize,
}

impl FeatureHasher {
    /// Creates a new feature hasher with the given number of buckets.
    ///
    /// # Panics
    /// Panics if `num_buckets` is not a power of 2.
    pub fn new(num_buckets: usize) -> Self {
        assert!(
            num_buckets > 0 && num_buckets.is_power_of_two(),
            "num_buckets must be a power of 2, got {}",
            num_buckets
        );
        Self { num_buckets }
    }

    /// Returns the number of features (= number of buckets).
    pub fn num_features(&self) -> usize {
        self.num_buckets
    }

    /// Transforms sequences into sparse feature vectors using the hashing trick.
    ///
    /// Returns rows sorted by feature index for efficient binary search during
    /// tree prediction — callers no longer need to sort manually.
    pub fn transform_sparse<S: AsRef<[u8]> + Sync>(
        &self,
        sequences: &[S],
        k: usize,
    ) -> Vec<Vec<(usize, f32)>> {
        let mask = self.num_buckets - 1;
        sequences
            .par_iter()
            .map(|seq| {
                let seq_bytes = seq.as_ref();
                let mut counts: FxHashMap<usize, f32> = FxHashMap::default();
                counts.reserve(8192);
                for (_, bitkmer, _) in seq_bytes.bit_kmers(k as u8, true) {
                    let idx = (bitkmer.0 as usize) & mask;
                    *counts.entry(idx).or_insert(0.0) += 1.0;
                }
                let mut sparse: Vec<(usize, f32)> = counts.into_iter().collect();
                sparse.sort_unstable_by_key(|&(f, _)| f);
                sparse
            })
            .collect()
    }

    /// Reverse-map bucket indices to the k-mer sequences that hash into them.
    pub fn reverse_map_buckets<S: AsRef<[u8]> + Sync>(
        &self,
        sequences: &[S],
        k: usize,
        target_buckets: &FxHashSet<usize>,
    ) -> FxHashMap<usize, Vec<String>> {
        info!(
            "  Reverse-mapping {} important buckets to k-mer sequences...",
            target_buckets.len()
        );
        let mask = self.num_buckets - 1;

        let per_thread: Vec<FxHashMap<usize, FxHashSet<u64>>> = sequences
            .par_iter()
            .map(|seq| {
                let seq_bytes = seq.as_ref();
                let mut local: FxHashMap<usize, FxHashSet<u64>> = FxHashMap::default();
                for (_, bitkmer, _) in seq_bytes.bit_kmers(k as u8, true) {
                    let idx = (bitkmer.0 as usize) & mask;
                    if target_buckets.contains(&idx) {
                        local.entry(idx).or_default().insert(bitkmer.0);
                    }
                }
                local
            })
            .collect();

        let mut merged: FxHashMap<usize, FxHashSet<u64>> = FxHashMap::default();
        for local in per_thread {
            for (bucket, kmers) in local {
                merged.entry(bucket).or_default().extend(kmers);
            }
        }

        merged
            .into_iter()
            .map(|(bucket, kmer_set)| {
                let strings: Vec<String> = kmer_set
                    .into_iter()
                    .map(|bits| bitkmer_to_string(bits, k))
                    .collect();
                (bucket, strings)
            })
            .collect()
    }

    /// Reverse-map target buckets WITH genomic coordinates.
    ///
    /// For each k-mer that hashes into a target bucket, records every occurrence
    /// (sequence index + 1-based position) across all training sequences.
    pub fn reverse_map_with_coords<S: AsRef<[u8]> + Sync>(
        &self,
        sequences: &[S],
        k: usize,
        target_buckets: &FxHashSet<usize>,
    ) -> Vec<KmerCoord> {
        info!(
            "  Mapping {} important buckets to genomic coordinates...",
            target_buckets.len()
        );
        let mask = self.num_buckets - 1;

        let per_seq: Vec<Vec<KmerCoord>> = sequences
            .par_iter()
            .enumerate()
            .map(|(seq_idx, seq)| {
                let seq_bytes = seq.as_ref();
                let mut hits = Vec::new();
                for (pos, bitkmer, _) in seq_bytes.bit_kmers(k as u8, true) {
                    let idx = (bitkmer.0 as usize) & mask;
                    if target_buckets.contains(&idx) {
                        hits.push(KmerCoord {
                            bucket: idx,
                            kmer: bitkmer_to_string(bitkmer.0, k),
                            seq_index: seq_idx,
                            position: pos + 1, // 1-based
                        });
                    }
                }
                hits
            })
            .collect();

        per_seq.into_iter().flatten().collect()
    }
}

// ---------------------------------------------------------------------------
// Supporting types
// ---------------------------------------------------------------------------

/// A genomic coordinate hit for a discriminant k-mer.
#[derive(Debug)]
pub struct KmerCoord {
    pub bucket: usize,
    pub kmer: String,
    pub seq_index: usize,
    /// 1-based position in the sequence.
    pub position: usize,
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Convert a 2-bit encoded k-mer (from needletail) back to a DNA string.
///
/// Needletail encoding: A=0, C=1, G=2, T=3.
fn bitkmer_to_string(bits: u64, k: usize) -> String {
    let mut s = vec![b'N'; k];
    for i in 0..k {
        let base_bits = (bits >> (2 * (k - 1 - i))) & 0b11;
        s[i] = match base_bits {
            0b00 => b'A',
            0b01 => b'C',
            0b10 => b'G',
            0b11 => b'T',
            _ => b'N',
        };
    }
    String::from_utf8(s).unwrap_or_else(|_| "N".repeat(k))
}
