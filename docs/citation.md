# Citation

**How to cite pathotypr, the DOIs for the software and marker panels, and the licence.**

If pathotypr contributed to your results, please cite the paper. Citing the
software DOI as well lets readers retrieve the exact marker panels and model you
ran against.

## Cite the paper

> Ruiz-Rodriguez P, Coscollá M. **Alignment-Free Tool for Assemblies and Reads in
> Genomic Surveillance of Tuberculosis.** *Emerg Infect Dis*. 2026;32(11):1913-1918.
> doi: [10.3201/eid3211.261074](https://doi.org/10.3201/eid3211.261074)

=== "BibTeX"

    ```bibtex
    @article{ruizrodriguez2026pathotypr,
      title   = {Alignment-Free Tool for Assemblies and Reads in Genomic
                 Surveillance of Tuberculosis},
      author  = {Ruiz-Rodriguez, Paula and Coscoll{\'a}, Mireia},
      journal = {Emerging Infectious Diseases},
      year    = {2026},
      volume  = {32},
      number  = {11},
      pages   = {1913--1918},
      doi     = {10.3201/eid3211.261074}
    }
    ```

=== "RIS"

    ```text
    TY  - JOUR
    TI  - Alignment-Free Tool for Assemblies and Reads in Genomic Surveillance of Tuberculosis
    AU  - Ruiz-Rodriguez, Paula
    AU  - Coscollá, Mireia
    JO  - Emerging Infectious Diseases
    PY  - 2026
    VL  - 32
    IS  - 11
    SP  - 1913
    EP  - 1918
    DO  - 10.3201/eid3211.261074
    ER  -
    ```

=== "APA"

    ```text
    Ruiz-Rodriguez, P., & Coscollá, M. (2026). Alignment-free tool for assemblies
    and reads in genomic surveillance of tuberculosis. Emerging Infectious
    Diseases, 32(11), 1913-1918. https://doi.org/10.3201/eid3211.261074
    ```

!!! note "Already cited the preprint?"
    The same work was first posted on bioRxiv as *Pathotypr: harmonised MTBC
    lineage assignment and resistance-associated variant detection for genomic
    surveillance*, doi:
    [10.64898/2026.03.24.714002](https://doi.org/10.64898/2026.03.24.714002).
    Please cite the published article above in new work.

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
| **Mireia Coscollá** | [ORCID 0000-0003-0752-0538](https://orcid.org/0000-0003-0752-0538) |

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
