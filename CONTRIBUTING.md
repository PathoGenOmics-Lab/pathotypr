# Contributing to pathotypr

Thank you for your interest in contributing to pathotypr. 🧬

## Where to go

| I want to | Go here |
|---|---|
| Ask how to run something | [Q&A discussions](https://github.com/PathoGenOmics-Lab/pathotypr/discussions/new?category=q-a) |
| Report a crash or an error | [Bug report](https://github.com/PathoGenOmics-Lab/pathotypr/issues/new?template=bug_report.yml) |
| Report a call I do not believe | [Unexpected results](https://github.com/PathoGenOmics-Lab/pathotypr/issues/new?template=unexpected_results.yml) |
| Suggest a capability | [Feature request](https://github.com/PathoGenOmics-Lab/pathotypr/issues/new?template=feature_request.yml) |
| Fix or question the docs | [Documentation problem](https://github.com/PathoGenOmics-Lab/pathotypr/issues/new?template=documentation.yml) |
| Propose a marker panel | [Marker proposals](https://github.com/PathoGenOmics-Lab/pathotypr/discussions/new?category=marker-proposals) |
| Send code | [Open a pull request](https://github.com/PathoGenOmics-Lab/pathotypr/pulls) |

Everything about using the tool lives on the
[documentation site](https://pathogenomics-lab.github.io/pathotypr/).

## Proposing new marker sets

pathotypr is organism-agnostic: the panel you supply defines what is typed.
Community marker proposals are currently accepted for the **Mycobacterium
tuberculosis complex (MTBC)** only, including *M. tuberculosis*, *M. bovis*,
*M. caprae*, *M. africanum*, *M. microti*, *M. pinnipedii*, *M. orygis*, and
related species.

### How to propose

1. **Prepare your markers** as a TSV in the format described in the
   [marker format reference](https://pathogenomics-lab.github.io/pathotypr/marker_format/).
2. **Test them locally** with `pathotypr classify` on a validation set, using the
   same reference genome the positions were called against.
3. **Open a discussion** in the
   [Marker proposals](https://github.com/PathoGenOmics-Lab/pathotypr/discussions/new?category=marker-proposals)
   category. A form will guide you through what is needed.
4. **Wait for automated validation.** A bot parses the TSV preview you paste and
   comments with what it found or what is malformed. Editing the discussion
   re-runs it.
5. **A maintainer reviews the full set.**

### Marker quality guidelines

| What we look for | Why |
|---|---|
| At least 20 markers per lineage | Robust classification with noise tolerance |
| Tested on at least 50 genomes per lineage | Reduces false positives |
| At least 99% concordance with known lineages | Ensures accuracy |
| Documented reference genome | Reproducibility |
| Published evidence (paper, preprint) | Scientific rigour |

### Marker file checklist

The full rules, including how the lineage hierarchy is split from the annotation
columns, are in the
[marker format reference](https://pathogenomics-lab.github.io/pathotypr/marker_format/).
The short version:

- [ ] Tab-separated, UTF-8, Unix line endings.
- [ ] Column order is `position`, `REF`, `ALT`, one or more lineage levels, then
      optionally an **empty separator column** followed by `gene` and `mutation`.
- [ ] One lineage level per column, root to leaf. Do not pack a hierarchy into a
      single semicolon-separated cell: `classify` tolerates it but `split-fastq`
      reads it as one literal level, and the two would then disagree.
- [ ] Positions are 1-based on the exact reference you pass with `-r`, which must
      be a single-record FASTA.
- [ ] Alleles are plain `A`, `C`, `G`, `T`. SNPs and equal-length MNVs work
      everywhere; indels are accepted by `classify` and skipped by `split-fastq`.
- [ ] Allele lengths fit the k-mer window: under the defaults that means 11 bp or
      shorter for `classify`.
- [ ] Each lineage has at least 5 markers, ideally 20 or more.
- [ ] A dry run with `-vv` reports no unexpected skips.

## Code contributions

1. Fork the repository.
2. Create a branch (`git checkout -b feature/my-feature`).
3. Make your changes and add tests.
4. Run `cargo test` and `cargo check`.
5. Open a pull request.

### Code style

- Rust: follow `rustfmt` defaults.
- Tests: unit tests for new functions, integration tests for new CLI commands.
- Docs: update the relevant pages under `docs/`. Every page on the site has an
  edit button that opens a pull request for you.

Every pull request is built and tested by CI. All checks must pass before it can
be merged.

## License

By contributing, you agree that your contributions will be licensed under the
AGPL-3.0 license.
