# Building the drug-resistance marker catalogue

These scripts regenerate `pathotypr_dr_markers_*.tsv` from the WHO catalogue of mutations
in *Mycobacterium tuberculosis* complex and their association with drug resistance
(second edition, WHO/UCN/TB/2023.7) and from the mutations curated from the literature by
TBProfiler. The released files are published on Zenodo; this directory documents how they
are produced so the catalogue can be audited and rebuilt when the WHO catalogue is updated.

The lineage marker panel and the pre-trained Random Forest model are produced by a
different workflow and are not covered here.

## Inputs

Place these in this directory before running:

| File | Source |
|---|---|
| `WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt` | WHO catalogue, genomic coordinates of each variant |
| `WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt` | WHO catalogue master file, one row per variant and drug |
| `tbdb_mutations.csv` | TBProfiler `mutations.csv` |
| `tbdb_genome.gff` | TBProfiler genome annotation |
| `sequence_L4.fasta` | H37Rv reference (NC_000962.3 = AL123456.3) |
| `MTB_ancestor_reference.fasta` | Inferred MTBC ancestor reference |

## Order

```bash
python3 convert_who_markers.py        # WHO catalogue  -> both coordinate frames
python3 resolve_tbdb_mutations.py     # appends the TBProfiler-only mutations (H37Rv)
python3 tbdb_markers_to_ancestor.py N # ports those N rows to ancestor coordinates
```

`resolve_tbdb_mutations.py` prints how many markers it appended; pass that number to
`tbdb_markers_to_ancestor.py`. It writes only to the H37Rv file, so without the third step
the ancestor catalogue would be missing the TBProfiler entries.

Outputs are `pathotypr_dr_markers_H37Rv.tsv` and `pathotypr_dr_markers_ancestor.tsv`.

## How a variant gets its drug

The WHO catalogue grades every variant–drug pair separately, so the same variant can be an
established marker for one drug and of uncertain significance for another. `rrs` n.1401A>G,
for instance, is grade 1 for amikacin, kanamycin and capreomycin but grade 3 for
streptomycin. Drug labels and grades are therefore read per variant from the master file,
never inferred from the gene.

A variant graded 1 or 2 for more than one drug carries a composite label listing them, as
the catalogue has always done for `BDQ_CFZ` and `FQ`: `AMI_KAN_CAP`, `INH_ETH`, `AMI_KAN`.
The grade of a composite label is the strongest grade among its drugs. One row per drug
would be equivalent for `split-fastq`, whose index holds several markers per k-mer, but
`classify` keeps a single marker per diagnostic k-mer, so the extra rows would overwrite
each other and only the last drug would ever be reported.

A variant that reaches neither grade 1 nor 2 for any drug is written once under the label
of its gene, taken from `GENE_DRUG` in `convert_who_markers.py`. Candidate genes such as
`glpK` were screened by the catalogue against most drugs and sit at grade 3–5 for all of
them; labelling those per drug would repeat the same uninformative variant once per
gene–drug pair. `GENE_DRUG` also decides which genes are covered at all.
