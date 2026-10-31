use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::hint::black_box;
use std::path::Path;

use pathotypr_core::classify;
use pathotypr_core::{FeatureHasher, DEFAULT_HASH_BUCKETS};

// Paths to test data. Set PATHOTYPR_TEST_DIR to override.
// Defaults to TEST_PATHOTYPR at the workspace root, which is the directory
// benches/run_benchmarks.sh looks in for the real data set.
// For synthetic data: bash pathotypr-core/benches/run_benchmarks.sh --synthetic
const DEFAULT_TEST_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../TEST_PATHOTYPR");

/// Marker k-mer length, matching the split-fastq default.
const MARKER_K: usize = pathotypr_core::split_kmer::DEFAULT_MARKER_KMER_LEN;

/// Flanking bases kept on each side of an allele inside a marker k-mer. There
/// is no library constant for it: this mirrors the `--min-flank-bases` default
/// of `pathotypr classify`, so the timings describe the configuration people
/// actually run.
const MIN_FLANK_BASES: usize = 10;

fn get_test_dir() -> String {
    std::env::var("PATHOTYPR_TEST_DIR").unwrap_or_else(|_| DEFAULT_TEST_DIR.to_string())
}

fn try_test_path(name: &str) -> Option<String> {
    let base = get_test_dir();
    let p = format!("{}/{}", base, name);
    if Path::new(&p).exists() {
        Some(p)
    } else {
        None
    }
}

/// Try multiple possible filenames, return the first that exists.
fn try_test_paths(names: &[&str]) -> Option<String> {
    for name in names {
        if let Some(p) = try_test_path(name) {
            return Some(p);
        }
    }
    None
}

// ===========================================================================
// Group 1: I/O, FASTA and GFF parsing
// ===========================================================================
fn bench_io(c: &mut Criterion) {
    let mut group = c.benchmark_group("io");
    group.sample_size(20);

    let ref_path = try_test_paths(&["MTB_ancestor_reference.fasta", "synthetic_reference.fasta"])
        .expect("No reference FASTA found");
    group.bench_function("get_ref (reference)", |b| {
        b.iter(|| classify::get_ref(black_box(&ref_path)).unwrap())
    });

    let sample_path = try_test_paths(&["sequence_L4.fasta", "synthetic_sample_L4.fasta"])
        .expect("No sample FASTA found");
    group.bench_function("get_ref (sample)", |b| {
        b.iter(|| classify::get_ref(black_box(&sample_path)).unwrap())
    });

    group.bench_function("get_genomes_from_fasta", |b| {
        b.iter(|| classify::get_genomes_from_fasta(black_box(&sample_path)).unwrap())
    });

    let markers_path = try_test_paths(&["c_set.txt", "synthetic_markers.tsv"])
        .expect("No markers file found");
    group.bench_function("get_positions (markers)", |b| {
        b.iter(|| classify::get_positions(black_box(&markers_path)).unwrap())
    });

    if let Some(gff_path) = try_test_paths(&["sequence_L4.gff3", "synthetic_annotation.gff3"]) {
        group.bench_function("parse_gff_and_build_tree", |b| {
            b.iter(|| classify::parse_gff_and_build_tree(black_box(&gff_path)).unwrap())
        });
    }

    group.finish();
}

// ===========================================================================
// Group 2: Marker generation and k-mer indexing
// ===========================================================================
fn bench_marker_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("marker_generation");
    group.sample_size(20);

    let ref_path = try_test_paths(&["MTB_ancestor_reference.fasta", "synthetic_reference.fasta"])
        .expect("No reference FASTA found");
    let markers_path = try_test_paths(&["c_set.txt", "synthetic_markers.tsv"])
        .expect("No markers file found");

    let ref_seq = classify::get_ref(&ref_path).unwrap();
    let markers = classify::get_positions(&markers_path).unwrap();

    for &k in &[21, MARKER_K] {
        group.bench_with_input(
            BenchmarkId::new("generate_markerkmer", format!("k={}", k)),
            &k,
            |b, &k| {
                b.iter(|| {
                    classify::generate_markerkmer(
                        black_box(&markers),
                        black_box(&ref_seq),
                        k,
                        MIN_FLANK_BASES,
                    )
                })
            },
        );
    }

    group.finish();
}

// ===========================================================================
// Group 3: K-mer matching (find_markers)
// ===========================================================================
fn bench_find_markers(c: &mut Criterion) {
    let mut group = c.benchmark_group("find_markers");
    group.sample_size(10);

    let ref_path = try_test_paths(&["MTB_ancestor_reference.fasta", "synthetic_reference.fasta"])
        .expect("No reference FASTA found");
    let markers_path = try_test_paths(&["c_set.txt", "synthetic_markers.tsv"])
        .expect("No markers file found");

    let ref_seq = classify::get_ref(&ref_path).unwrap();
    let markers = classify::get_positions(&markers_path).unwrap();
    let k = MARKER_K;
    let marker_index = classify::generate_markerkmer(&markers, &ref_seq, k, MIN_FLANK_BASES);

    group.bench_function(format!("find_markers (reference, k={})", k), |b| {
        b.iter(|| classify::find_markers(black_box(&ref_seq), black_box(&marker_index), k))
    });

    let sample_path = try_test_paths(&["sequence_L4.fasta", "synthetic_sample_L4.fasta"])
        .expect("No sample FASTA found");
    let sample_seq = classify::get_ref(&sample_path).unwrap();
    group.bench_function(format!("find_markers (sample, k={})", k), |b| {
        b.iter(|| classify::find_markers(black_box(&sample_seq), black_box(&marker_index), k))
    });

    group.finish();
}

// ===========================================================================
// Group 4: Reverse complement
// ===========================================================================
fn bench_reverse_complement(c: &mut Criterion) {
    let mut group = c.benchmark_group("reverse_complement");

    let ref_path = try_test_paths(&["MTB_ancestor_reference.fasta", "synthetic_reference.fasta"])
        .expect("No reference FASTA found");
    let ref_seq = classify::get_ref(&ref_path).unwrap();

    // The input is whichever reference the test directory holds, so the name
    // carries no size: 4.3 MB for the MTB ancestor, far less for the synthetic
    // one, and a figure in the name would be wrong for one of them.
    group.bench_function("reverse_complement_sequence", |b| {
        b.iter(|| classify::reverse_complement_sequence(black_box(&ref_seq)))
    });

    group.finish();
}

// ===========================================================================
// Group 5: Full genome analysis (analyze_genome_seq)
// ===========================================================================
fn bench_analyze_genome(c: &mut Criterion) {
    let mut group = c.benchmark_group("analyze_genome");
    group.sample_size(10);

    let ref_path = try_test_paths(&["MTB_ancestor_reference.fasta", "synthetic_reference.fasta"])
        .expect("No reference FASTA found");
    let markers_path = try_test_paths(&["c_set.txt", "synthetic_markers.tsv"])
        .expect("No markers file found");

    let ref_seq = classify::get_ref(&ref_path).unwrap();
    let ref_seq_rc = classify::reverse_complement_sequence(&ref_seq);
    let markers = classify::get_positions(&markers_path).unwrap();
    let k = MARKER_K;
    let marker_index = classify::generate_markerkmer(&markers, &ref_seq, k, MIN_FLANK_BASES);

    let gff_path = try_test_paths(&["sequence_L4.gff3", "synthetic_annotation.gff3"]);
    let annotations = gff_path
        .as_ref()
        .map(|p| classify::parse_gff_and_build_tree(p).unwrap());

    let sample_path = try_test_paths(&["sequence_L4.fasta", "synthetic_sample_L4.fasta"])
        .expect("No sample FASTA found");
    let sample_seq = classify::get_ref(&sample_path).unwrap();

    group.bench_function("analyze_genome_seq (no GFF)", |b| {
        b.iter(|| {
            classify::analyze_genome_seq(
                black_box("test_sample"),
                black_box(&sample_seq),
                black_box(&marker_index),
                black_box(&None),
                k,
                black_box(&ref_seq),
                black_box(&ref_seq_rc),
            )
        })
    });

    group.bench_function("analyze_genome_seq (with GFF)", |b| {
        b.iter(|| {
            classify::analyze_genome_seq(
                black_box("test_sample"),
                black_box(&sample_seq),
                black_box(&marker_index),
                black_box(&annotations),
                k,
                black_box(&ref_seq),
                black_box(&ref_seq_rc),
            )
        })
    });

    group.finish();
}

// ===========================================================================
// Group 6: split_kmer, marker building and FASTQ scanning
// ===========================================================================
fn bench_split_kmer(c: &mut Criterion) {
    let mut group = c.benchmark_group("split_kmer");
    group.sample_size(10);

    let ref_path = try_test_paths(&["MTB_ancestor_reference.fasta", "synthetic_reference.fasta"])
        .expect("No reference FASTA found");
    let markers_path = try_test_paths(&["c_set.txt", "synthetic_markers.tsv"])
        .expect("No markers file found");

    group.bench_function("build_markers", |b| {
        b.iter(|| {
            pathotypr_core::split_kmer::build_markers(
                black_box(&ref_path),
                black_box(&markers_path),
                &None,
                MARKER_K,
            )
            .unwrap()
        })
    });

    let markers =
        pathotypr_core::split_kmer::build_markers(&ref_path, &markers_path, &None, MARKER_K)
            .unwrap();
    group.bench_function("build_marker_index", |b| {
        b.iter(|| pathotypr_core::split_kmer::build_marker_index(black_box(&markers)))
    });

    // Try real FASTQ first, then synthetic
    if let Some(fastq_path) =
        try_test_paths(&["SRR36978519.fastq.gz", "synthetic_reads_R1.fastq.gz"])
    {
        let (index, bloom) = pathotypr_core::split_kmer::build_marker_index(&markers);
        group.sample_size(10);
        group.measurement_time(std::time::Duration::from_secs(60));
        group.bench_function("scan_fastq_with_index", |b| {
            b.iter(|| {
                pathotypr_core::split_kmer::scan_fastq_with_index(
                    black_box(&[fastq_path.clone()]),
                    black_box(&index),
                    black_box(&bloom),
                    markers.len(),
                    &None,
                    MARKER_K,
                )
                .unwrap()
            })
        });
    }

    group.finish();
}

// ===========================================================================
// Group 7: Feature hashing (vectorizer.rs)
// ===========================================================================
fn bench_vectorizer(c: &mut Criterion) {
    let mut group = c.benchmark_group("vectorizer");
    group.sample_size(10);

    let sample_path = try_test_paths(&["sequence_L4.fasta", "synthetic_sample_L4.fasta"])
        .expect("No sample FASTA found");
    let sample_seq = classify::get_ref(&sample_path).unwrap();

    // This group used to time a CountVectorizer built from a 10,000 k-mer
    // vocabulary. Feature hashing replaced it and needs no vocabulary pass, so
    // the setup is one constant and the measurement is the transform alone,
    // which is where train and predict spend their vectorizing time.
    let hasher = FeatureHasher::new(DEFAULT_HASH_BUCKETS);
    let sequences = vec![sample_seq.as_str()];
    let k = MARKER_K;

    group.bench_function(
        format!("transform_sparse (1 seq, {} buckets)", DEFAULT_HASH_BUCKETS),
        |b| b.iter(|| hasher.transform_sparse(black_box(&sequences), k)),
    );

    group.finish();
}

// ===========================================================================
// Group 8: train, FASTA reading and data preparation
// ===========================================================================
fn bench_train(c: &mut Criterion) {
    let mut group = c.benchmark_group("train");
    group.sample_size(10);

    // Try training FASTA (multi-seq) or single sample
    let training_path = try_test_paths(&["synthetic_training.fasta", "sequence_L4.fasta"])
        .expect("No training FASTA found");

    group.bench_function("read_fasta", |b| {
        b.iter(|| pathotypr_core::read_fasta(black_box(&training_path)).unwrap())
    });

    // Collect several sequences for the vectorizing step
    let mut sequences = Vec::new();
    for name in &[
        "sequence_L4.fasta",
        "synthetic_sample_L4.fasta",
        "MTB_ancestor_reference.fasta",
        "synthetic_reference.fasta",
        "L7_DE0090.fasta",
        "E1ASM0009.circlator.fasta",
        "E1ASM0057.circlator.fasta",
    ] {
        if let Some(p) = try_test_path(name) {
            sequences.push(classify::get_ref(&p).unwrap());
        }
        if sequences.len() >= 5 {
            break;
        }
    }

    // prepare_data needs sequences plus labels with at least two classes
    let labels: Vec<String> = (0..sequences.len())
        .map(|i| {
            if i % 2 == 0 {
                "L4".to_string()
            } else {
                "L7".to_string()
            }
        })
        .collect();

    let n_seqs = sequences.len();
    group.bench_function(format!("prepare_data ({} seqs, k=21)", n_seqs), |b| {
        b.iter(|| {
            pathotypr_core::train::prepare_data(
                black_box(&sequences),
                black_box(&labels),
                21,
                &None,
            )
            .unwrap()
        })
    });

    group.finish();
}

// ===========================================================================
// Group 9: predict, model loading, FASTA reading and prediction
// ===========================================================================
fn bench_predict(c: &mut Criterion) {
    let mut group = c.benchmark_group("predict");
    group.sample_size(10);

    let sample_path = try_test_paths(&["sequence_L4.fasta", "synthetic_sample_L4.fasta"])
        .expect("No sample FASTA found");
    group.bench_function("read_fasta_for_prediction", |b| {
        b.iter(|| {
            pathotypr_core::predict::read_fasta_for_prediction(black_box(&sample_path)).unwrap()
        })
    });

    // A trained model is neither part of the real test directory nor something
    // the synthetic generator produces, so the three model benchmarks run only
    // when one is handed to them: PATHOTYPR_BENCH_MODEL, or a model under one
    // of the names below inside the test directory. Build one with
    // `pathotypr train`.
    let model_path = std::env::var("PATHOTYPR_BENCH_MODEL")
        .ok()
        .filter(|p| Path::new(p).exists())
        .or_else(|| try_test_paths(&["test_model.pathotypr", "test_model.pathotypr.zst"]));

    if let Some(model_path) = model_path {
        if let Ok(bundle) = pathotypr_core::predict::load_model_bundle(&model_path) {
            group.bench_function("load_model_bundle", |b| {
                b.iter(|| {
                    pathotypr_core::predict::load_model_bundle(black_box(&model_path)).unwrap()
                })
            });

            let records = pathotypr_core::predict::read_fasta_for_prediction(&sample_path).unwrap();
            let sequences: Vec<&str> = records.iter().map(|(_, seq)| seq.as_str()).collect();
            let k = bundle.config.kmer_size;

            group.bench_function("transform_sparse (predict)", |b| {
                b.iter(|| bundle.vectorizer.transform_sparse(black_box(&sequences), k))
            });

            let x_sparse = bundle.vectorizer.transform_sparse(&sequences, k);
            let n_classes = bundle.label_encoder.int_to_label.len();

            // The trees vote straight off the sparse row, so this is the loop
            // predict::run takes per record. The dense matrix this benchmark
            // used to build belonged to the smartcore trees that
            // SparseDecisionTree replaced.
            group.bench_function(
                format!("forest_predict ({} trees)", bundle.trees.len()),
                |b| {
                    b.iter(|| {
                        for sparse_row in &x_sparse {
                            let mut votes = vec![0u16; n_classes];
                            for tree in &bundle.trees {
                                let pred = tree.predict_one(sparse_row);
                                if pred < n_classes {
                                    votes[pred] += 1;
                                }
                            }
                            black_box(&votes);
                        }
                    })
                },
            );
        }
    }

    group.finish();
}

// ===========================================================================
// Group 10: match, reference loading
// ===========================================================================
fn bench_match(c: &mut Criterion) {
    let mut group = c.benchmark_group("match");
    group.sample_size(10);

    let sample_path = try_test_paths(&["sequence_L4.fasta", "synthetic_sample_L4.fasta"])
        .expect("No sample FASTA found");
    group.bench_function("read_references_from_multifasta", |b| {
        b.iter(|| {
            pathotypr_core::r#match::read_references_from_multifasta(black_box(&sample_path))
                .unwrap()
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_io,
    bench_marker_generation,
    bench_find_markers,
    bench_reverse_complement,
    bench_analyze_genome,
    bench_split_kmer,
    bench_vectorizer,
    bench_train,
    bench_predict,
    bench_match,
);
criterion_main!(benches);
