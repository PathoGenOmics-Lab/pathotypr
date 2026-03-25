# Contributing to pathotypr

Thank you for your interest in contributing to pathotypr! 🧬

## Proposing New Marker Sets

pathotypr is organism-agnostic — it works with any set of lineage-defining markers. We currently accept community marker proposals for the **Mycobacterium tuberculosis complex (MTBC)** only, including *M. tuberculosis*, *M. bovis*, *M. caprae*, *M. africanum*, *M. microti*, *M. pinnipedii*, *M. orygis*, and related species.

### How to Propose

1. **Prepare your markers** in TSV format ([see format docs](docs/marker_format.md))
2. **Test locally** with `pathotypr classify` on a validation set
3. **Open a Discussion** in the [Marker Proposals](https://github.com/mycolega/pathotypr/discussions/categories/ideas) category
4. **Fill in the template** with organism, reference, evidence, and a TSV preview
5. **Wait for automated validation** — a bot will check your TSV format
6. **A maintainer will review** the full marker set

### Marker Quality Guidelines

| What we look for | Why |
|---|---|
| ≥ 20 markers per lineage | Robust classification with noise tolerance |
| Tested on ≥ 50 genomes per lineage | Reduces false positives |
| ≥ 99% concordance with known lineages | Ensures accuracy |
| Documented reference genome | Reproducibility |
| Published evidence (paper, preprint) | Scientific rigor |

### Marker File Checklist

- [ ] Tab-separated, UTF-8 encoded
- [ ] Header: `position	ref	alt	lineage`
- [ ] 1-based positions relative to the reference
- [ ] Single nucleotide variants only (no indels)
- [ ] Each lineage has ≥ 5 markers (ideally ≥ 20)
- [ ] Tested on independent validation genomes

## Reporting Bugs

Use [Issues](https://github.com/mycolega/pathotypr/issues) with:
- pathotypr version (`pathotypr --version`)
- OS and architecture
- Input data description
- Error message or unexpected output
- Steps to reproduce

## Feature Requests

Open a [Discussion](https://github.com/mycolega/pathotypr/discussions/categories/ideas) describing:
- What you'd like to do
- Why it would be useful
- Any related tools or papers

## Code Contributions

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Make changes and add tests
4. Run `cargo test` and `cargo check`
5. Open a Pull Request

### Code Style

- Rust: follow `rustfmt` defaults
- Tests: add unit tests for new functions, integration tests for new CLI commands
- Docs: update relevant markdown files in `docs/`

## License

By contributing, you agree that your contributions will be licensed under the AGPL-3.0 license.
