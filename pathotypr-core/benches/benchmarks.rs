use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use std::path::Path;

use pathotypr_core::classify;

// Paths to test data — set PATHOTYPR_TEST_DIR env var to override.
// Defaults to TEST_PATHOTYPR alongside the workspace root.
// For synthetic data: PATHOTYPR_TEST_DIR=.../TEST_SYNTHETIC cargo bench
const DEFAULT_TEST_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../TEST_PATHOTYPR"
);

fn get_test_dir() -> String {
    std::env::var("PATHOTYPR_TEST_DIR").unwrap_or_else(|_| DEFAULT_TEST_DIR.to_string())
}

fn test_path(name: &str) -> String {
    let base = get_test_dir();
    let p = format!("{}/{}", base, name);
    assert!(Path::new(&p).exists(), "Test file not found: {}", p);
    p
}

fn try_test_path(name: &str) -> Option<String> {
    let base = get_test_dir();
    let p = format!("{}/{}", base, name);
    if Path::new(&p).exists() { Some(p) } else { None }
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
// Group 1: I/O — FASTA / GFF parsing
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
// Group 2: Marker generation + k-mer indexing
// ===========================================================================
fn bench_marker_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("marker_generation");
    group.sample_size(20);

    let ref_path = try_test_paths(&["MTB_ancestor_reference.fasta", "synthetic_reference.fasta"])
        .expect("No reference FASTA found");
    let markers_path = try_test_paths(&["c_set.txt", "synthetic_markers.tsv"])
        .expect("No markers file found");

    let ref_seq = classify::get_ref(&ref_path).unwrap();
    let (reference_positions, markers_lineage) = classify::get_positions(&markers_path).unwrap();

    for &k in &[21, 31] {
        group.bench_with_input(
            BenchmarkId::new("generate_markerkmer", format!("k={}", k)),
            &k,
            |b, &k| {
                b.iter(|| {
                    classify::generate_markerkmer(
                        black_box(&reference_positions),
                        black_box(&ref_seq),
                        black_box(&markers_lineage),
                        k,
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
    let (reference_positions, markers_lineage) = classify::get_positions(&markers_path).unwrap();
    let k = 31;
    let marker_index =
        classify::generate_markerkmer(&reference_positions, &ref_seq, &markers_lineage, k);

    group.bench_function("find_markers (reference, k=31)", |b| {
        b.iter(|| classify::find_markers(black_box(&ref_seq), black_box(&marker_index), k))
    });

    let sample_path = try_test_paths(&["sequence_L4.fasta", "synthetic_sample_L4.fasta"])
        .expect("No sample FASTA found");
    let sample_seq = classify::get_ref(&sample_path).unwrap();
    group.bench_function("find_markers (sample, k=31)", |b| {
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

    group.bench_function("reverse_complement_sequence (4.3 MB)", |b| {
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
    let (reference_positions, markers_lineage) = classify::get_positions(&markers_path).unwrap();
    let k = 31;
    let marker_index =
        classify::generate_markerkmer(&reference_positions, &ref_seq, &markers_lineage, k);

    let gff_path = try_test_paths(&["sequence_L4.gff3", "synthetic_annotation.gff3"]);
    let annotations = gff_path.as_ref().map(|p| classify::parse_gff_and_build_tree(p).unwrap());

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
// Group 6: split_kmer — marker building + FASTQ scanning
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
            )
            .unwrap()
        })
    });

    let markers =
        pathotypr_core::split_kmer::build_markers(&ref_path, &markers_path, &None).unwrap();
    group.bench_function("build_marker_index", |b| {
        b.iter(|| pathotypr_core::split_kmer::build_marker_index(black_box(&markers)))
    });

    // Try real FASTQ first, then synthetic
    if let Some(fastq_path) = try_test_paths(&[
        "SRR36978519.fastq.gz",
        "synthetic_reads_R1.fastq.gz",
    ]) {
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
                )
                .unwrap()
            })
        });
    }

    group.finish();
}

// ===========================================================================
// Group 7: CountVectorizer (common.rs)
// ===========================================================================
fn bench_vectorizer(c: &mut Criterion) {
    let mut group = c.benchmark_group("vectorizer");
    group.sample_size(10);

    let sample_path = try_test_paths(&["sequence_L4.fasta", "synthetic_sample_L4.fasta"])
        .expect("No sample FASTA found");
    let sample_seq = classify::get_ref(&sample_path).unwrap();

    use needletail::Sequence;
    use rustc_hash::FxHashMap;
    let k = 31u8;
    let mut vocab = FxHashMap::default();
    let mut idx = 0usize;
    for (_, bitkmer, _) in sample_seq.as_bytes().bit_kmers(k, true) {
        vocab.entry(bitkmer.0).or_insert_with(|| {
            let i = idx;
            idx += 1;
            i
        });
        if idx >= 10000 {
            break;
        }
    }
    let vectorizer = pathotypr_core::CountVectorizer {
        vocabulary: vocab,
        num_features: idx,
    };

    let sequences: Vec<&[u8]> = vec![sample_seq.as_bytes()];
    group.bench_function("transform_sparse (1 seq, 10k features)", |b| {
        b.iter(|| vectorizer.transform_sparse(black_box(&sequences), k as usize))
    });

    group.finish();
}

// ===========================================================================
// Group 8: train — FASTA reading, k-mer counting, data preparation
// ===========================================================================
fn bench_train(c: &mut Criterion) {
    let mut group = c.benchmark_group("train");
    group.sample_size(10);

    // Try training FASTA (multi-seq) or single sample
    let training_path = try_test_paths(&["synthetic_training.fasta", "sequence_L4.fasta"])
        .expect("No training FASTA found");

    group.bench_function("read_fasta", |b| {
        b.iter(|| pathotypr_core::train::read_fasta(black_box(&training_path)).unwrap())
    });

    // Collect multiple sequences for parallel k-mer counting
    let mut sequences = Vec::new();
    for name in &[
        "sequence_L4.fasta", "synthetic_sample_L4.fasta",
        "MTB_ancestor_reference.fasta", "synthetic_reference.fasta",
        "L7_DE0090.fasta", "E1ASM0009.circlator.fasta", "E1ASM0057.circlator.fasta",
    ] {
        if let Some(p) = try_test_path(name) {
            sequences.push(classify::get_ref(&p).unwrap());
        }
        if sequences.len() >= 5 { break; }
    }

    let n_seqs = sequences.len();
    group.bench_function(
        &format!("parallel_kmer_discovery ({} seqs, k=21)", n_seqs),
        |b| {
            b.iter(|| {
                pathotypr_core::train::parallel_kmer_discovery(
                    black_box(&sequences),
                    21,
                    &None,
                )
                .unwrap()
            })
        },
    );

    // prepare_data — needs sequences + labels with >=2 classes
    let labels: Vec<String> = (0..sequences.len())
        .map(|i| if i % 2 == 0 { "L4".to_string() } else { "L7".to_string() })
        .collect();

    group.bench_function(
        &format!("prepare_data ({} seqs, k=21)", n_seqs),
        |b| {
            b.iter(|| {
                pathotypr_core::train::prepare_data(
                    black_box(&sequences),
                    black_box(&labels),
                    21,
                    &None,
                )
                .unwrap()
            })
        },
    );

    group.finish();
}

// ===========================================================================
// Group 9: predict — model loading, FASTA reading, prediction
// ===========================================================================
fn bench_predict(c: &mut Criterion) {
    let mut group = c.benchmark_group("predict");
    group.sample_size(10);

    let sample_path = try_test_paths(&["sequence_L4.fasta", "synthetic_sample_L4.fasta"])
        .expect("No sample FASTA found");
    group.bench_function("read_fasta_for_prediction", |b| {
        b.iter(|| {
            pathotypr_core::predict::read_fasta_for_prediction(black_box(&sample_path))
                .unwrap()
        })
    });

    // load_model_bundle — only if compatible model exists
    if let Some(model_path) = try_test_path("test_bundle.bin.gz") {
        if let Ok(bundle) = pathotypr_core::predict::load_model_bundle(&model_path) {
            group.bench_function("load_model_bundle", |b| {
                b.iter(|| {
                    pathotypr_core::predict::load_model_bundle(black_box(&model_path)).unwrap()
                })
            });

            let records =
                pathotypr_core::predict::read_fasta_for_prediction(&sample_path).unwrap();
            let sequences: Vec<&str> = records.iter().map(|(_, seq)| seq.as_str()).collect();

            group.bench_function("transform_sparse (predict)", |b| {
                b.iter(|| {
                    bundle
                        .vectorizer
                        .transform_sparse(black_box(&sequences), bundle.config.kmer_size)
                })
            });

            let x_sparse =
                bundle
                    .vectorizer
                    .transform_sparse(&sequences, bundle.config.kmer_size);
            let n_features = bundle.vectorizer.num_features;
            let n_classes = bundle.label_encoder.int_to_label.len();

            group.bench_function("forest_predict (100 trees)", |b| {
                b.iter(|| {
                    use smartcore::linalg::basic::matrix::DenseMatrix;
                    for sparse_row in &x_sparse {
                        let mut dense = vec![0.0_f32; n_features];
                        for &(idx, val) in sparse_row {
                            if idx < n_features {
                                dense[idx] = val;
                            }
                        }
                        let x = DenseMatrix::new(1, n_features, dense, false).unwrap();
                        let mut votes = vec![0u16; n_classes];
                        for tree in &bundle.trees {
                            let pred = tree.predict(&x).unwrap();
                            if pred[0] < n_classes {
                                votes[pred[0]] += 1;
                            }
                        }
                        black_box(&votes);
                    }
                })
            });
        }
    }

    group.finish();
}

// ===========================================================================
// Group 10: match — reference loading
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
