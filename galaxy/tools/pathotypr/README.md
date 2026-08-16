# Galaxy wrappers for pathotypr

Five tools, one per subcommand, sharing `macros.xml`.

| Tool | Wraps | Outputs |
|---|---|---|
| `pathotypr_train` | `train` | model, and optionally the two feature-importance reports |
| `pathotypr_predict` | `predict` | predictions TSV |
| `pathotypr_classify` | `classify` | per-marker hits, summary, optionally masked FASTA |
| `pathotypr_split_fastq` | `split-fastq` | per-marker calls, summary |
| `pathotypr_match` | `match` | best-match report |

## Running the tests

```bash
python3 -m venv .venv && .venv/bin/pip install planemo
.venv/bin/planemo lint tools/pathotypr/*.xml
.venv/bin/planemo test tools/pathotypr/
```

`planemo test` resolves `pathotypr` from Bioconda. On a machine where that
package has no build for the local platform, put a `pathotypr` binary on `PATH`
and add `--no_dependency_resolution`; that still exercises the command lines,
the outputs and the assertions, and leaves only the dependency resolution to CI.

## The test data

Everything under `test-data/` is synthetic and generated deterministically:
a 2 kb reference, eight markers, two samples, paired reads at about 18x, and a
twelve-genome training set in two classes. 64 KB in total.

It is not filler. The generator asserts that **every marker k-mer occurs exactly
once in the reference**, in both orientations, because a marker that matched in
two places would make the tests pass for the wrong reason. The expected calls
follow from how the samples were built:

- `sample1` carries the L2 to L2.1 to L2.1.1 path plus the RIF marker, so
  `classify` calls L2 and reports `rpoB` / `S450L`.
- `sample2` carries L3 to L3.2 plus the INH marker.
- `split-fastq` on `sample1`'s reads produces the **same summary** as `classify`
  on its assembly, which is the property the marker format promises.
- `predict` assigns the two held-out queries to L2 and L3.
- `match` picks `ref_A`, which is `sample1`'s own genome.

## Three things the wrappers have to work around

1. **Sample names come from filenames.** Galaxy datasets arrive as
   `dataset_NNN.dat`, so `classify` and `split-fastq` would name their outputs
   after a Galaxy id. Each input is symlinked to a stable name first.
2. **`match` echoes input paths into its report.** Same fix, otherwise the first
   column is full of absolute paths that differ between instances.
3. **`--min-alt-percent` is an integer.** A Galaxy `float` parameter renders
   `95.0`, which the tool rejects outright. It is declared as an integer, so
   fractional thresholds such as 99.5% are not available.

## Submitting to the IUC

The directory layout matches `galaxyproject/tools-iuc`, so submission is a copy:

```bash
cp -r galaxy/tools/pathotypr <tools-iuc>/tools/pathotypr
```

then a pull request against that repository. Their CI runs the same planemo
lint and tests on Linux, where the Bioconda package resolves.

## Not done yet

Marker panels and models come from the history. A **data manager** plus a
`.loc` table would let an administrator install the published MTBC panels once
so users pick them from a dropdown, which is how Galaxy normally handles
reference data. That is a separate piece of work and was left out deliberately
rather than half-built: shipping a data table with nothing to populate it moves
the burden to administrators without helping anyone.
