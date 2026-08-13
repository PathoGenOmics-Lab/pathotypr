#!/usr/bin/env python3
"""
Convert WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt to Pathotypr marker TSV format.

Input:  WHO catalogue genomic coordinates (variant, chromosome, position, ref, alt)
Output: Pathotypr-compatible marker TSV for DR classification

Reference: sequence_L4.fasta (H37Rv AL123456.3 = NC_000962.3, same coordinates)

Format:
  pos  ref  alt  DRUG  DRUG-R  gene_mutation  grade  gene  mutation

The catalogue grades every variant x drug pair separately, so a variant can be graded 1
for one drug and 3 for another. Drug labels and grades are therefore taken per variant
from the master file, and a variant graded 1-2 for several drugs gets a composite label
such as AMI_KAN_CAP, in the same way the marker files already use BDQ_CFZ and FQ.
GENE_DRUG below decides which genes are covered and labels the variants that carry no
grade 1-2 evidence.

Deriving the label from the gene instead misreports variants whose gene serves more than
one drug: rrs n.1401A>G is grade 1 for amikacin, kanamycin and capreomycin but only
grade 3 for streptomycin, and the inhA promoter variants are graded for both isoniazid
and ethionamide.
"""

import csv
import re
import sys
from collections import defaultdict

WHO_COORDS_FILE = "WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"
WHO_MASTER_FILE = "WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt"
OUTPUT_L4 = "pathotypr_dr_markers_H37Rv.tsv"
OUTPUT_ANC = "pathotypr_dr_markers_ancestor.tsv"
REF_L4 = "sequence_L4.fasta"
REF_ANC = "MTB_ancestor_reference.fasta"

# Gene → Drug fallback, used only for variants that are absent from the master file
# (their drug comes from the catalogue itself; see load_variant_drug_grades)
GENE_DRUG = {
    # First-line
    "rpoB": "RIF",
    "katG": "INH", "katG_LoF": "INH",
    "inhA": "INH", "fabG1": "INH", "fabG1-inhA": "INH",
    "embB": "EMB", "embA": "EMB", "embC": "EMB",
    "embA-embB": "EMB", "embC-embA": "EMB",
    "pncA": "PZA", "pncA_LoF": "PZA",
    "panD": "PZA", "panD_LoF": "PZA",
    "rpsL": "STR", "rrs": "STR",
    "gid": "STR", "gid_LoF": "STR",

    # Fluoroquinolones
    "gyrA": "FQ", "gyrB": "FQ",

    # Second-line injectables
    "eis": "KAN",  # eis promoter
    "tlyA": "CAP", "tlyA_LoF": "CAP",

    # Ethionamide
    "ethA": "ETH", "ethA_LoF": "ETH",
    "ethR": "ETH",

    # New/repurposed drugs
    "rplC": "LZD",
    "rrl": "LZD",  # rrl 23S rRNA - also BDQ/CFZ cross?
    "Rv0678": "BDQ_CFZ", "Rv0678_LoF": "BDQ_CFZ",
    "mmpR5": "BDQ_CFZ", "mmpR5_LoF": "BDQ_CFZ",
    "mmpS5": "BDQ_CFZ", "mmpS5_LoF": "BDQ_CFZ",
    "mmpL5": "BDQ_CFZ", "mmpL5_LoF": "BDQ_CFZ",
    "atpE": "BDQ",
    "pepQ": "BDQ_CFZ", "pepQ_LoF": "BDQ_CFZ",
    "Rv1979c": "BDQ_CFZ", "Rv1979c_LoF": "BDQ_CFZ",
    "ddn": "DLM", "ddn_LoF": "DLM",
    "fbiA": "DLM", "fbiA_LoF": "DLM",
    "fbiB": "DLM", "fbiB_LoF": "DLM",
    "fbiC": "DLM", "fbiC_LoF": "DLM",
    "fgd1": "DLM", "fgd1_LoF": "DLM",

    # Additional genes in WHO catalogue
    "rpoA": "RIF",  # compensatory
    "rpoC": "RIF",  # compensatory
    "ndh": "INH",
    "ahpC": "INH",
    "kasA": "INH",
    "Rv2752c": "ETH",  # hypothetical

    # Previously excluded genes — now included
    "dnaA": "OTHER",        # replication-associated
    "PPE35": "OTHER",       # PE/PPE family
    "aftB": "EMB",          # cell wall (arabinosyltransferase, EMB target pathway)
    "mtrB": "OTHER",        # two-component system
    "mshA": "ETH",          # mycothiol biosynthesis (ETH activation)
    "glpK": "OTHER",        # glycerol kinase
    "clpC1": "OTHER",       # protease (pyrazinamide-related, experimental)
    "Rv2477c": "OTHER",     # hypothetical
    "whiB6": "OTHER",       # transcription factor
    "ubiA": "EMB",          # decaprenylphosphoryl-5-phosphoribose oxidase (EMB pathway)
}

# WHO drug name → Pathotypr label. Levofloxacin and moxifloxacin stay collapsed as FQ,
# and bedaquiline+clofazimine as BDQ_CFZ when a variant is graded for both, so that the
# label vocabulary of previous marker files is preserved.
WHO_DRUG_ABBR = {
    "Rifampicin": "RIF", "Isoniazid": "INH", "Ethambutol": "EMB",
    "Pyrazinamide": "PZA", "Streptomycin": "STR", "Ethionamide": "ETH",
    "Amikacin": "AMI", "Kanamycin": "KAN", "Capreomycin": "CAP",
    "Linezolid": "LZD", "Delamanid": "DLM",
    "Levofloxacin": "LEV", "Moxifloxacin": "MXF",
    "Bedaquiline": "BDQ", "Clofazimine": "CFZ",
}

GRADE_RANK = {"1)": 0, "2)": 1, "3)": 2, "4)": 3, "5)": 4}


def strongest(grades):
    """Return the lowest-numbered (strongest) grade of an iterable."""
    return sorted(grades, key=lambda g: GRADE_RANK.get(g[:2], 9))[0]


LABEL_ORDER = ["RIF", "INH", "EMB", "PZA", "STR", "FQ", "ETH", "AMI", "KAN", "CAP",
               "LZD", "BDQ", "CFZ", "DLM"]


def collapse_labels(drug_grades):
    """{WHO drug: grade} → (label, grade) for the drugs graded 1-2, or None.

    A variant associated with several drugs gets one composite label rather than one row
    per drug, following the convention the marker files already use for BDQ_CFZ and FQ.
    One row per drug would be equivalent for split-fastq, whose index holds several
    markers per k-mer, but classify keeps a single marker per k-mer, so the extra rows
    would overwrite each other and only the last drug would ever be reported.

    The composite carries the strongest grade of its member drugs, again as BDQ_CFZ and
    FQ already do.
    """
    abbr = {WHO_DRUG_ABBR[d]: g for d, g in drug_grades.items()
            if d in WHO_DRUG_ABBR and GRADE_RANK.get(g[:2], 9) <= 1}
    if not abbr:
        return None

    fq = [abbr.pop(d) for d in ("LEV", "MXF") if d in abbr]
    if fq:
        abbr["FQ"] = strongest(fq)

    parts = sorted(abbr, key=lambda l: LABEL_ORDER.index(l) if l in LABEL_ORDER else 99)
    return "_".join(parts), strongest(abbr.values())


def load_reference(fasta_path):
    """Load a FASTA reference sequence, return as uppercase string."""
    seq_lines = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq_lines.append(line.strip().upper())
    return "".join(seq_lines)


def get_ref_bases(ref_seq, pos_1based, length):
    """Extract bases from reference at 1-based position."""
    start = pos_1based - 1  # Convert to 0-based
    return ref_seq[start:start + length]


def parse_variant_name(variant):
    """Extract gene and mutation from WHO variant name.

    Examples:
        rpoB_p.Ser450Leu -> (rpoB, S450L)
        katG_p.Ser315Thr -> (katG, S315T)
        inhA_c.-15C>T    -> (fabG1-inhA, c-15t)
        rpoB_c.1349C>T   -> (rpoB, c.1349C>T)
        pncA_LoF_p.Trp68* -> (pncA_LoF, W68*)
    """
    # Split gene from mutation
    parts = variant.split("_", 1)
    if len(parts) < 2:
        return variant, variant

    gene_raw = parts[0]
    mut_part = parts[1]

    # Handle LoF genes
    if mut_part.startswith("LoF_"):
        gene_raw = gene_raw + "_LoF"
        mut_part = mut_part[4:]  # Remove "LoF_"

    # Three-letter to one-letter amino acid mapping
    aa3to1 = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    }

    if mut_part.startswith("p."):
        # Protein change: p.Ser450Leu
        prot = mut_part[2:]
        # Match: AminoAcid Position AminoAcid_or_*
        m = re.match(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|del|ins.*)", prot)
        if m:
            ref_aa = aa3to1.get(m.group(1), m.group(1))
            pos = m.group(2)
            alt_raw = m.group(3)
            alt_aa = aa3to1.get(alt_raw, alt_raw)
            mutation = f"{ref_aa}{pos}{alt_aa}"
        else:
            mutation = prot
    elif mut_part.startswith("c.") or mut_part.startswith("n."):
        # Coding/non-coding DNA change: c.-15C>T or n.-15C>T
        mutation = mut_part
    elif mut_part.startswith("r."):
        # RNA level
        mutation = mut_part
    else:
        mutation = mut_part

    return gene_raw, mutation


def get_display_gene(gene_raw):
    """Clean gene name for display (remove _LoF suffix)."""
    return gene_raw.replace("_LoF", "")


def load_variant_drug_grades(master_file):
    """Load FINAL CONFIDENCE GRADING per variant AND drug from the WHO master file.

    The master file has one row per variant+drug, so a variant graded for several drugs
    carries a separate grade for each. Keeping only the first grade found would attach
    the grade of one drug to the label of another.
    """
    grades = defaultdict(dict)
    with open(master_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            variant = row.get("variant", "").strip()
            drug = row.get("drug", "").strip()
            grade = row.get("FINAL CONFIDENCE GRADING", "").strip()
            if variant and drug and grade:
                grades[variant][drug] = grade
    return grades


def main():
    # Load reference sequences
    print("Loading references...")
    ref_l4 = load_reference(REF_L4)
    ref_anc = load_reference(REF_ANC)
    print(f"  L4 (H37Rv):  {len(ref_l4):,} bp")
    print(f"  Ancestor:    {len(ref_anc):,} bp")

    # Count differences between references
    diffs = sum(1 for a, b in zip(ref_l4, ref_anc) if a != b)
    print(f"  Differences: {diffs:,} positions")
    print()

    # Load grades from master file
    print("Loading confidence grades from master file...")
    variant_grades = load_variant_drug_grades(WHO_MASTER_FILE)
    n_multi = sum(1 for dg in variant_grades.values()
                  if (r := collapse_labels(dg)) is not None and "_" in r[0])
    print(f"  Grades loaded for {len(variant_grades):,} variants")
    print(f"  Of these, {n_multi:,} are graded 1-2 for more than one drug")

    grade_counts = defaultdict(int)
    for dg in variant_grades.values():
        grade_counts[strongest(dg.values())] += 1
    for g, n in sorted(grade_counts.items()):
        print(f"    {g}: {n:,}")
    print()

    # Parse WHO genomic coordinates file (ALL variants, no grade filter)
    variants = []
    skipped_no_drug = defaultdict(int)

    with open(WHO_COORDS_FILE) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            variant_name = row["variant"]
            pos = int(row["position"])
            ref = row["reference_nucleotide"]
            alt = row["alternative_nucleotide"]

            gene_raw, mutation = parse_variant_name(variant_name)
            display_gene = get_display_gene(gene_raw)

            # GENE_DRUG remains the gene whitelist: it decides which genes are covered.
            # What it no longer decides is the drug of each variant, which now comes from
            # the catalogue. Letting every catalogued gene in would widen the panel and
            # pack more markers into the same regions, which costs FASTA markers to k-mer
            # collisions in classify.
            gene_drug = GENE_DRUG.get(gene_raw) or GENE_DRUG.get(display_gene)
            if gene_drug is None:
                skipped_no_drug[gene_raw] += 1
                continue

            # Drug label and grade come from the catalogue, per variant
            resolved = collapse_labels(variant_grades.get(variant_name, {}))

            if resolved is not None:
                drug_label, grade = resolved
            else:
                # No drug reaches grade 1-2: keep the gene-level label, with the strongest
                # grade the variant carries. Candidate genes such as glpK were screened
                # against most drugs and sit at grade 3-5 for all of them, so labelling
                # them by drug would say nothing the gene does not already say.
                all_grades = variant_grades.get(variant_name, {}).values()
                drug_label = gene_drug
                grade = strongest(all_grades) if all_grades else "Unknown"

            variants.append({
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "drug": drug_label,
                "resistance": f"{drug_label}-R",
                "marker_name": f"{display_gene}_{mutation}",
                "gene": display_gene,
                "mutation": mutation,
                "variant_name": variant_name,
                "grade": grade,
            })

    # Deduplicate: same pos + ref + alt → keep the strongest grade
    seen = {}
    for v in variants:
        key = (v["pos"], v["ref"], v["alt"])
        kept = seen.get(key)
        if kept is None or GRADE_RANK.get(v["grade"][:2], 9) < GRADE_RANK.get(kept["grade"][:2], 9):
            seen[key] = v
    unique_variants = list(seen.values())

    # Sort by drug then position
    drug_order = ["RIF", "INH", "EMB", "PZA", "STR", "FQ", "ETH", "AMI",
                  "KAN", "CAP", "LZD", "BDQ", "CFZ", "BDQ_CFZ", "DLM", "OTHER"]

    def sort_key(v):
        d = v["drug"]
        try:
            return (drug_order.index(d), v["pos"])
        except ValueError:
            return (999, v["pos"])

    unique_variants.sort(key=sort_key)

    # Write both output files
    def write_markers(output_path, variants, ref_seq, ref_name):
        """Write marker TSV, replacing REF column with bases from the given reference."""
        n_ref_mismatch = 0
        n_ref_equals_alt = 0
        written = 0

        with open(output_path, "w") as f:
            f.write("#pos\tref\talt\tdrug\tresistance\tmarker_name\tgrade\tgene\tmutation\n")
            for v in variants:
                who_ref = v["ref"]
                ref_from_seq = get_ref_bases(ref_seq, v["pos"], len(who_ref))

                # Use the reference base from the actual FASTA
                actual_ref = ref_from_seq if ref_from_seq else who_ref

                # Skip if ref == alt in this reference (variant is the reference state)
                if actual_ref == v["alt"]:
                    n_ref_equals_alt += 1
                    continue

                if actual_ref != who_ref:
                    n_ref_mismatch += 1

                f.write(f"{v['pos']}\t{actual_ref}\t{v['alt']}\t"
                        f"{v['drug']}\t{v['resistance']}\t{v['marker_name']}\t"
                        f"{v['grade']}\t{v['gene']}\t{v['mutation']}\n")
                written += 1

        print(f"  {ref_name}:")
        print(f"    Written:          {written:,}")
        print(f"    REF≠WHO (changed):{n_ref_mismatch:>7,}")
        print(f"    REF=ALT (skipped):{n_ref_equals_alt:>7,}")
        return written

    print("Writing marker files...")
    n_l4 = write_markers(OUTPUT_L4, unique_variants, ref_l4, "L4 (H37Rv)")
    n_anc = write_markers(OUTPUT_ANC, unique_variants, ref_anc, "Ancestor")
    print()

    # Summary
    print(f"Total WHO variants parsed:     {len(variants)}")
    print(f"After deduplication:           {len(unique_variants)}")
    print()

    # By drug
    drug_counts = defaultdict(int)
    for v in unique_variants:
        drug_counts[v["drug"]] += 1

    print("Markers per drug:")
    for d in drug_order:
        if d in drug_counts:
            print(f"  {d:12s} {drug_counts[d]:>6,}")
    print()

    # By variant type
    snp = sum(1 for v in unique_variants if len(v["ref"]) == 1 and len(v["alt"]) == 1)
    mnv = sum(1 for v in unique_variants if len(v["ref"]) > 1 or len(v["alt"]) > 1)
    print(f"SNPs:  {snp:>6,}")
    print(f"MNVs:  {mnv:>6,}")
    print()

    # By grade
    out_grade_counts = defaultdict(int)
    for v in unique_variants:
        out_grade_counts[v["grade"]] += 1
    print("Markers per confidence grade:")
    for g, n in sorted(out_grade_counts.items()):
        print(f"  {g:40s} {n:>6,}")
    print()

    # Skipped
    if skipped_no_drug:
        print(f"\nSkipped (no drug mapping): {sum(skipped_no_drug.values()):,}")
        for g, n in sorted(skipped_no_drug.items(), key=lambda x: -x[1])[:10]:
            print(f"  {g}: {n:,}")

    print(f"\nOutput L4:       {OUTPUT_L4}")
    print(f"Output Ancestor: {OUTPUT_ANC}")


if __name__ == "__main__":
    main()
