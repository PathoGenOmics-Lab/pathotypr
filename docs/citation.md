# Citation

**How to cite pathotypr, the DOIs for the software and marker panels, and the licence.**

If pathotypr contributed to your results, please cite the preprint. Citing the
software DOI as well lets readers retrieve the exact marker panels and model you
ran against.

## Cite the paper

> Ruiz-Rodriguez P, Coscollá M. **Pathotypr: harmonised MTBC lineage assignment
> and resistance-associated variant detection for genomic surveillance.**
> *bioRxiv* (2026). doi: [10.64898/2026.03.24.714002](https://doi.org/10.64898/2026.03.24.714002)

=== "BibTeX"

    ```bibtex
    @article{ruizrodriguez2026pathotypr,
      title   = {Pathotypr: harmonised MTBC lineage assignment and
                 resistance-associated variant detection for genomic surveillance},
      author  = {Ruiz-Rodriguez, Paula and Coscoll{\'a}, Mireia},
      journal = {bioRxiv},
      year    = {2026},
      doi     = {10.64898/2026.03.24.714002}
    }
    ```

=== "RIS"

    ```text
    TY  - JOUR
    TI  - Pathotypr: harmonised MTBC lineage assignment and resistance-associated
          variant detection for genomic surveillance
    AU  - Ruiz-Rodriguez, Paula
    AU  - Coscollá, Mireia
    JO  - bioRxiv
    PY  - 2026
    DO  - 10.64898/2026.03.24.714002
    ER  -
    ```

=== "APA"

    ```text
    Ruiz-Rodriguez, P., & Coscollá, M. (2026). Pathotypr: harmonised MTBC lineage
    assignment and resistance-associated variant detection for genomic
    surveillance. bioRxiv. https://doi.org/10.64898/2026.03.24.714002
    ```

## Cite the software and markers

The archived release, including the marker panels and the trained model, has its
own DOI:

[10.5281/zenodo.19210043](https://doi.org/10.5281/zenodo.19210043)

!!! tip "Record the version you ran"
    Marker panels change between releases, so a lineage or resistance call is
    only reproducible alongside the version that produced it. `pathotypr
    --version` prints it, and every output file records it in its header.

The repository also ships a [`CITATION.cff`](https://github.com/PathoGenOmics-Lab/pathotypr/blob/main/CITATION.cff)
file, so GitHub's **Cite this repository** button always returns an up-to-date
reference.

## Authors

| | |
|---|---|
| **Paula Ruiz-Rodriguez** | [ORCID 0000-0003-0727-5974](https://orcid.org/0000-0003-0727-5974) |
| **Mireia Coscollá** | [ORCID 0000-0003-2078-1032](https://orcid.org/0000-0003-2078-1032) |

Institute for Integrative Systems Biology (I²SysBio), CSIC and University of
Valencia.

## Licence

pathotypr is released under the
[GNU Affero General Public License v3.0](https://github.com/PathoGenOmics-Lab/pathotypr/blob/main/LICENSE).
You may use, modify and redistribute it, including in a service offered over a
network, provided derivative works stay under the same licence and their source
remains available.

## See also

- [Benchmarks](benchmarks.md) for the performance figures reported in the paper.
- [Marker format](marker_format.md) for how to describe a panel you curated yourself.
