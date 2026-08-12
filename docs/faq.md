# FAQ & troubleshooting

Answers to the questions that come up most often when running `pathotypr`, grouped by workflow. If your issue isn't here, check the per-command guides ([`train`](train.md), [`predict`](predict.md), [`classify`](classify.md), [`split-fastq`](split-fastq.md), [`match`](match.md)) or open an issue on [GitHub](https://github.com/PathoGenOmics-Lab/pathotypr).

!!! tip "First thing to try: turn up the logs"
    Almost every "why did it do that?" question is answered faster by re-running with `-v` (debug) or `-vv` (trace). See [Seeing more detail](#seeing-more-detail-v-vv) at the bottom of this page.

## Installation & first launch

!!! question "The desktop app won't open — macOS says it's from an *unidentified developer*, or Windows SmartScreen blocks it"
    The GUI is **not** signed with a paid developer certificate, so on the **first launch** the OS shows a warning. This is expected and only happens once.

    === "macOS"

        Right-click (or Control-click) the app in **Applications** → **Open** → click **Open** in the dialog. After the first time it launches normally.

    === "Windows"

        If **SmartScreen** appears, click **More info** → **Run anyway**.

    Full per-platform details are on the [Installation](installation.md#desktop-gui-pre-built) page.

## K-mer size

!!! question "My run stopped with `k-mer size must be between 1 and 31`"
    Every command that accepts a k-mer size validates it to the range **1–31** (the k-mer must fit in a 2-bit-packed `u64`). A value of `0` or anything above `31` aborts the run with:

    ```text
    k-mer size must be between 1 and 31, got 41
    ```

    Pick a value inside 1–31 and re-run.

Defaults and flags differ per command:

| Command | Flag | Default | Notes |
|---|---|---|---|
| [`train`](train.md) | `-k, --kmer-size` | `21` | Baked into the model bundle |
| [`predict`](predict.md) | — | — | Inherited from the model; not a flag |
| [`classify`](classify.md) | `--kmer-size` | `31` | Long form only |
| [`split-fastq`](split-fastq.md) | `-k, --kmer-size` | `31` | See odd/11–31 note below |
| [`match`](match.md) | `-k, --kmer-size` | `31` | Longer k = more specific |

!!! tip "split-fastq: prefer odd k between 11 and 31"
    For `split-fastq`, the `-k` help specifies that diagnostic k-mers should be **odd and between 11 and 31**, so the variant sits centred with symmetric flanking bases on each side. Keep the default `31`, or choose an odd value in that range (e.g. `21`, `27`). Values outside 1–31 still error as above.

## Input data

!!! question "Does `classify` work on multi-contig draft assemblies?"
    **Yes — in every input mode.** Draft assemblies from bacterial WGS routinely contain many contigs, and `classify` reads *every* record in a query FASTA and scans all of them for markers. The single-record restriction applies **only** to the reference genome.

!!! warning "The *reference* (`-r`) must be a single record"
    Both `classify` and `split-fastq` build markers against a **single-record** reference FASTA. A multi-record reference stops the run:

    ```text
    Reference FASTA 'ref.fasta' contains multiple records; provide a single-record FASTA.
    ```

    Your **query** genomes have no such limit — only the reference does. See [Input formats → Reference FASTA](input-formats.md#reference-fasta).

!!! tip "How multi-contig results are grouped"
    - **`-l, --input-list`** — all contigs in each listed FASTA are aggregated under that one sample name. This is the natural mode for *one draft assembly = one sample*.
    - **`-i, --input`** and **`--input-files`** — each FASTA record is reported as its own genome (identified by its sequence ID; `--input-files` also prefixes the file name).

    If you want a multi-contig assembly scored as a single sample, list it with `--input-list`.

## Paired-end reads (`split-fastq`)

!!! question "My paired-end FASTQs aren't being paired, or the wrong files got paired together"
    By default `split-fastq` **auto-detects** read pairs from file names. It recognises these suffix conventions:

    | Pattern | Example |
    |---|---|
    | `_R1_001` / `_R2_001` | `sample_R1_001.fastq.gz`, `sample_R2_001.fastq.gz` |
    | `_R1` / `_R2` | `sample_R1.fq.gz`, `sample_R2.fq.gz` |
    | `_1` / `_2` | `sample_1.fastq.gz`, `sample_2.fastq.gz` |
    | `.1` / `.2` | `sample.1.fastq.gz`, `sample.2.fastq.gz` |

    If your files don't follow one of these conventions, take control of pairing explicitly:

    === "Let it auto-detect (default)"

        No flags needed — just pass the files:

        ```bash
        pathotypr split-fastq \
          -m markers.tsv -r reference.fasta \
          -i sample_R1.fastq.gz -i sample_R2.fastq.gz \
          -o run
        ```

    === "Force manual pairing (--paired)"

        Files are grouped **sequentially in pairs**, so list R1 then its R2, and so on. An odd number of files is an error (`--paired requires an even number of input files.`):

        ```bash
        pathotypr split-fastq \
          -m markers.tsv -r reference.fasta \
          -i a_fwd.fq.gz -i a_rev.fq.gz \
          -i b_fwd.fq.gz -i b_rev.fq.gz \
          --paired -o run
        ```

    === "Treat every file as single-end (--no-auto-paired)"

        Disables auto-detection; each file becomes its own single-end sample:

        ```bash
        pathotypr split-fastq \
          -m markers.tsv -r reference.fasta \
          -i lane1.fq.gz -i lane2.fq.gz \
          --no-auto-paired -o run
        ```

    === "Name the samples yourself (--input-list)"

        The most robust option for irregular names — you assign the sample name and its file(s) in a TSV (one FASTQ = single-end, two = paired):

        ```bash
        pathotypr split-fastq \
          -m markers.tsv -r reference.fasta \
          -l samples.tsv -o run
        ```

    !!! note "Precedence"
        `--input-list` wins if provided; otherwise `--paired` forces sequential pairing; otherwise auto-detection runs unless `--no-auto-paired` turns it off. See the [Sample list format](input-formats.md#sample-list-tsv-split-fastq-and-match-input-list).

## Performance & memory

!!! question "A run used a lot of RAM — which command is the heaviest, and how do I cap it?"
    `match` is typically the most memory-demanding workflow: it holds the entire query k-mer table in memory while scoring references. Crucially, it scores references in **streaming batches**, so peak memory stays **bounded and does *not* grow with the total number of references** (1 or 500+, memory stays flat).

    What the batch memory *does* scale with is the **thread count**: each batch holds `num_threads` references (clamped to 4–64). To lower peak RAM, reduce `--threads`:

    ```bash
    pathotypr match \
      -l samples.tsv -r references.fasta \
      -o match.tsv --threads 4
    ```

    Typical usage lands around a few hundred MB. Noise filtering (`--min-kmer-count`, default `2`) also trims the query k-mer table once unique k-mers exceed 100K. See [Reference matching → Memory model](algorithms/reference-matching.md#memory-model).

## Output

!!! question "How do I get Excel (`.xlsx`) output instead of just TSV?"
    Add the `--excel` flag. The command still writes its normal TSV files and additionally produces `.xlsx` versions alongside them:

    ```bash
    pathotypr classify \
      -m markers.tsv -r reference.fasta \
      -i sample.fasta -o run --excel
    ```

    `--excel` is available on [`classify`](classify.md), [`split-fastq`](split-fastq.md), [`predict`](predict.md), and [`match`](match.md).

    !!! note
        [`train`](train.md) does not have `--excel` — its importance reports are written as TSV only.

## Reference data

!!! question "Where do the MTBC markers and the pre-trained model live?"
    Ready-to-use marker panels and a pre-trained Random Forest model for the *Mycobacterium tuberculosis* complex (MTBC) are published on Zenodo under the concept DOI [10.5281/zenodo.19210043](https://doi.org/10.5281/zenodo.19210043), which always resolves to the newest version:

    | File | Description |
    |---|---|
    | `pathotypr_lineage_markers_*.tsv` | 3,707 lineage SNPs (L1–L10, A1–A4) |
    | `pathotypr_dr_markers_ancestor_*.tsv` | DR mutations from the WHO catalogue (2nd edition, 2023), ancestor coordinates |
    | `pathotypr_dr_markers_H37Rv_*.tsv` | The same catalogue in H37Rv coordinates |
    | `pathotypr_rf_model_*.pathotypr` | Pre-trained RF model (k=31, 100 trees) |

    Filenames carry the catalogue version, so they change between releases; the desktop app resolves the newest deposit on its own. Download links are on the [Installation](installation.md#mtbc-marker-files-pre-trained-model) page. Pass the marker TSV to `-m` (for `classify`/`split-fastq`) and the model to `-m/--model` (for `predict`).

!!! question "Can I use pathotypr for an organism other than *M. tuberculosis*?"
    Technically yes — nothing is hard-coded to one species. The marker panel and
    reference you supply define the organism, and `train` builds a model from
    whatever labelled genomes you give it.

    But pathotypr has **only been developed and validated on the MTBC**. The
    published panels, the pre-trained model and every benchmark on this site are
    MTBC, and no other species has been tested. If you point it at another
    organism, treat the output as exploratory: validate the calls against a
    truth set you trust before using them for anything that matters. See the
    [marker format](marker_format.md) reference for how to write your own panel.

## Seeing more detail (`-v` / `-vv`)

!!! question "How do I make a run print more about what it's doing?"
    Every subcommand accepts the global, repeatable verbosity flag. It works before or after the subcommand:

    | Flag | Log level |
    |---|---|
    | *(none)* | info |
    | `-v` | debug |
    | `-vv` (or more) | trace |

    ```bash
    # both forms are equivalent
    pathotypr -vv split-fastq -m markers.tsv -r reference.fasta -i reads.fq.gz -o run
    pathotypr split-fastq -m markers.tsv -r reference.fasta -i reads.fq.gz -o run -vv
    ```

    Debug and trace surface per-sample, per-batch, and (for `train`) per-fold detail that is otherwise hidden. `-h/--help` is available on every subcommand; `-V/--version` works only on the top-level command (`pathotypr --version`).

## See also

- [Installation](installation.md) — downloads, Bioconda, building from source, and the first-launch note.
- [Input formats](input-formats.md) — marker TSV, sample lists, reference and GFF3 requirements.
- [Desktop GUI](gui.md) — the same workflows with drag-and-drop and progress bars.
