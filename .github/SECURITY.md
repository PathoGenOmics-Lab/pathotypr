# Security policy

## Reporting a vulnerability

Please report privately through
[GitHub security advisories](https://github.com/PathoGenOmics-Lab/pathotypr/security/advisories/new)
rather than opening a public issue. That keeps the details between us until
there is a fix to point people at.

If you would rather not use GitHub, write to paula.ruiz.rodriguez@csic.es.

Please include what you found, how to reproduce it, and which version and
platform you saw it on. If you have a proof of concept, a minimal one is more
useful than a complete one.

You can expect an acknowledgement within a week. If the report is valid you will
be credited in the advisory unless you ask otherwise.

## What is in scope

pathotypr is a command line tool and a desktop application that reads files you
give it and, in the desktop app, downloads marker panels and models over HTTPS.
The interesting surfaces are therefore:

- **Parsing.** A crafted FASTA, FASTQ, TSV, GFF or model bundle that causes a
  crash, an out-of-bounds read, unbounded memory growth, or code execution.
- **The model bundle.** `predict` deserializes a file produced by `train`. A
  bundle that does something other than load a model when opened is in scope.
- **Downloads in the desktop app.** Anything that lets a download escape the URL
  validation or the SSRF guard, write outside the chosen directory, or reach a
  local or reserved address.
- **Paths.** Output filenames derived from input data that escape the output
  directory.

## What is not in scope

- **A wrong scientific call.** A marker panel that misclassifies, or a model that
  predicts badly, is a correctness problem rather than a security one. Please
  use the [Unexpected results](https://github.com/PathoGenOmics-Lab/pathotypr/issues/new?template=unexpected_results.yml)
  form, which asks the questions needed to chase it down.
- **Resource use on enormous but well-formed input.** Genomics files are large by
  nature; that a 200 GB FASTQ needs time and memory is expected behaviour.
- **Advisories against a dependency that this project cannot reach.** These are
  still worth reporting, and are still fixed, but through a normal issue or pull
  request rather than a private advisory.

## Supported versions

Fixes land on the latest release. There are no long-term support branches, so
the answer to "which version is patched" is always the newest one.

Dependency updates arrive monthly and grouped, and GitHub security advisories
against a dependency open their own pull request the day they are published,
independently of that schedule.
