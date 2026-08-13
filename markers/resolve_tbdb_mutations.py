#!/usr/bin/env python3
"""
Resolve TBProfiler-exclusive mutations (source=tbdb) to genomic coordinates
using the GFF annotation and H37Rv reference, then append to pathotypr_dr_markers.tsv.

Mutation types handled:
  - p.Ser450Leu  → protein: codon_pos = (AA_pos-1)*3+1 in CDS → genomic
  - c.1210C>T    → coding DNA position in CDS → genomic
  - c.-15C>T     → promoter: upstream of CDS start → genomic
  - n.1401A>G    → non-coding (rRNA): position within gene → genomic
  - c.*NNN       → downstream of stop codon
  - del/ins      → deletions, insertions (complex)
"""

import csv
import re
import sys
from collections import defaultdict

# --- File paths (run from markers_dr/) ---
TBDB_MUTATIONS = "tbdb_mutations.csv"
GFF_FILE = "tbdb_genome.gff"
REF_FASTA = "sequence_L4.fasta"
OUTPUT_MARKERS = "pathotypr_dr_markers_H37Rv.tsv"

# Three-letter to one-letter amino acid
AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

# Codon table (standard genetic code)
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C",
              "a": "t", "t": "a", "c": "g", "g": "c", "N": "N"}


def rev_comp(seq):
    return "".join(COMPLEMENT.get(b, b) for b in reversed(seq))


def load_reference(fasta_path):
    seq_lines = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq_lines.append(line.strip().upper())
    return "".join(seq_lines)


def load_gff_genes(gff_path):
    """Load gene coordinates from GFF. Returns dict: gene_name -> {start, end, strand, rv_id}"""
    genes = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type != "gene":
                continue

            start = int(parts[3])  # 1-based
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]

            name = None
            rv_id = None
            for attr in attrs.split(";"):
                if attr.startswith("Name="):
                    name = attr.split("=")[1]
                elif attr.startswith("gene_id="):
                    rv_id = attr.split("=")[1]

            if name:
                genes[name] = {"start": start, "end": end, "strand": strand, "rv_id": rv_id}
            if rv_id and rv_id not in genes:
                genes[rv_id] = {"start": start, "end": end, "strand": strand, "rv_id": rv_id}

    return genes


def get_ref_base(ref_seq, pos_1based, length=1):
    return ref_seq[pos_1based - 1: pos_1based - 1 + length]


def cds_to_genomic(gene_info, cds_pos):
    """Convert CDS position (1-based) to genomic position (1-based)."""
    if gene_info["strand"] == "+":
        return gene_info["start"] + cds_pos - 1
    else:
        return gene_info["end"] - cds_pos + 1


def promoter_to_genomic(gene_info, offset):
    """Convert promoter position (negative, e.g. -15) to genomic position."""
    if gene_info["strand"] == "+":
        return gene_info["start"] + offset  # offset is negative
    else:
        return gene_info["end"] - offset  # offset is negative, so this adds


def noncoding_to_genomic(gene_info, nc_pos):
    """Convert non-coding position to genomic (same as CDS for rRNA genes)."""
    if gene_info["strand"] == "+":
        return gene_info["start"] + nc_pos - 1
    else:
        return gene_info["end"] - nc_pos + 1


def resolve_protein_mutation(gene_info, ref_seq, mutation):
    """Resolve p.Ser450Leu to genomic pos, ref, alt.

    Returns list of (pos, ref, alt) tuples (one per codon base that changes).
    """
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|del|ins.*)", mutation)
    if not m:
        return None

    ref_aa_3 = m.group(1)
    aa_pos = int(m.group(2))
    alt_aa_3 = m.group(3)

    ref_aa = AA3TO1.get(ref_aa_3, ref_aa_3)
    alt_aa = AA3TO1.get(alt_aa_3, alt_aa_3)

    # Stop codon or deletion — skip (can't resolve to simple SNP)
    if alt_aa in ("*", "del") or alt_aa_3.startswith("ins"):
        # For stop codons, try to find the codon change
        if alt_aa == "*":
            pass  # continue below
        else:
            return None

    # CDS position of first base in codon
    codon_cds_start = (aa_pos - 1) * 3 + 1

    # Get the 3 genomic positions for this codon
    results = []
    if gene_info["strand"] == "+":
        gpos1 = gene_info["start"] + codon_cds_start - 1
        ref_codon = ref_seq[gpos1 - 1: gpos1 + 2]
    else:
        gpos3 = gene_info["end"] - codon_cds_start + 1  # 3rd base of codon in fwd
        gpos1 = gpos3 - 2
        ref_codon = rev_comp(ref_seq[gpos1 - 1: gpos1 + 2])

    # Verify ref amino acid
    ref_aa_check = CODON_TABLE.get(ref_codon, "?")
    if ref_aa_check != ref_aa:
        # Mismatch — might be annotation issue, skip
        return None

    # Find all possible alt codons that encode alt_aa
    if alt_aa == "*":
        target_codons = [c for c, aa in CODON_TABLE.items() if aa == "*"]
    else:
        target_codons = [c for c, aa in CODON_TABLE.items() if aa == alt_aa]

    # Find the alt codon with minimum changes from ref_codon
    best_codon = None
    best_changes = 4
    for tc in target_codons:
        changes = sum(1 for a, b in zip(ref_codon, tc) if a != b)
        if changes < best_changes:
            best_changes = changes
            best_codon = tc

    if best_codon is None:
        return None

    # Generate individual SNP changes
    for i in range(3):
        if ref_codon[i] != best_codon[i]:
            if gene_info["strand"] == "+":
                gpos = gene_info["start"] + codon_cds_start - 1 + i
                ref_base = ref_codon[i]
                alt_base = best_codon[i]
            else:
                gpos = gene_info["end"] - codon_cds_start + 1 - i
                ref_base = COMPLEMENT[ref_codon[i]]
                alt_base = COMPLEMENT[best_codon[i]]

            results.append((gpos, ref_base, alt_base))

    return results if results else None


def resolve_coding_snp(gene_info, ref_seq, mutation):
    """Resolve c.1210C>T or c.-15C>T to genomic pos, ref, alt."""
    # Promoter: c.-15C>T
    m = re.match(r"c\.(-?\d+)([ACGT])>([ACGT])", mutation)
    if not m:
        return None

    pos = int(m.group(1))
    ref_base = m.group(2)
    alt_base = m.group(3)

    if pos < 0:
        # Promoter position
        gpos = promoter_to_genomic(gene_info, pos)
    else:
        # CDS position
        gpos = cds_to_genomic(gene_info, pos)

    # For reverse strand, complement the bases
    if gene_info["strand"] == "-":
        ref_base = COMPLEMENT.get(ref_base, ref_base)
        alt_base = COMPLEMENT.get(alt_base, alt_base)

    # Verify against reference
    actual_ref = get_ref_base(ref_seq, gpos)
    if actual_ref != ref_base:
        # Try without complementing (some entries are already in genomic coords)
        if actual_ref == m.group(2):
            ref_base = m.group(2)
            alt_base = m.group(3)
        else:
            return None

    return [(gpos, ref_base, alt_base)]


def resolve_noncoding_snp(gene_info, ref_seq, mutation):
    """Resolve n.1401A>G to genomic pos, ref, alt."""
    m = re.match(r"n\.(\d+)([ACGT])>([ACGT])", mutation)
    if not m:
        return None

    pos = int(m.group(1))
    ref_base = m.group(2)
    alt_base = m.group(3)

    gpos = noncoding_to_genomic(gene_info, pos)

    if gene_info["strand"] == "-":
        ref_base = COMPLEMENT.get(ref_base, ref_base)
        alt_base = COMPLEMENT.get(alt_base, alt_base)

    actual_ref = get_ref_base(ref_seq, gpos)
    if actual_ref != ref_base:
        if actual_ref == m.group(2):
            ref_base = m.group(2)
            alt_base = m.group(3)
        else:
            return None

    return [(gpos, ref_base, alt_base)]


# Drug name normalization
DRUG_MAP = {
    "rifampicin": "RIF",
    "isoniazid": "INH",
    "ethambutol": "EMB",
    "pyrazinamide": "PZA",
    "streptomycin": "STR",
    "levofloxacin": "FQ",
    "moxifloxacin": "FQ",
    "fluoroquinolones": "FQ",
    "ethionamide": "ETH",
    "kanamycin": "KAN",
    "amikacin": "AMI",
    "capreomycin": "CAP",
    "linezolid": "LZD",
    "bedaquiline": "BDQ",
    "clofazimine": "CFZ",
    "delamanid": "DLM",
    "para-aminosalicylic_acid": "PAS",
    "cycloserine": "DCS",
}


def main():
    print("Loading reference...")
    ref_seq = load_reference(REF_FASTA)
    print(f"  H37Rv: {len(ref_seq):,} bp")

    print("Loading GFF annotations...")
    genes = load_gff_genes(GFF_FILE)
    print(f"  Genes loaded: {len(genes)}")

    # Load existing markers to avoid duplicates
    existing = set()
    with open(OUTPUT_MARKERS) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                existing.add((parts[0], parts[1], parts[2]))  # pos, ref, alt

    print(f"  Existing markers: {len(existing):,}")

    # Parse tbdb-exclusive mutations
    tbdb_mutations = []
    with open(TBDB_MUTATIONS) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["source"] != "tbdb":
                continue
            tbdb_mutations.append(row)

    print(f"  TBProfiler-exclusive mutations: {len(tbdb_mutations)}")
    print()

    # Resolve each mutation
    resolved = 0
    skipped_no_gene = 0
    skipped_no_resolve = 0
    skipped_duplicate = 0
    new_markers = []

    for mut in tbdb_mutations:
        gene_name = mut["Gene"]
        mutation = mut["Mutation"]
        drug_raw = mut["drug"]
        confidence = mut.get("confidence", "")

        # Map drug name
        drug = DRUG_MAP.get(drug_raw, "OTHER")

        # Find gene in GFF
        gene_info = genes.get(gene_name)
        if gene_info is None:
            skipped_no_gene += 1
            continue

        # Resolve mutation to genomic coordinates
        positions = None
        if mutation.startswith("p."):
            positions = resolve_protein_mutation(gene_info, ref_seq, mutation)
        elif mutation.startswith("c."):
            positions = resolve_coding_snp(gene_info, ref_seq, mutation)
        elif mutation.startswith("n."):
            positions = resolve_noncoding_snp(gene_info, ref_seq, mutation)

        if positions is None:
            skipped_no_resolve += 1
            continue

        for gpos, ref_base, alt_base in positions:
            key = (str(gpos), ref_base, alt_base)
            if key in existing:
                skipped_duplicate += 1
                continue

            existing.add(key)

            # Map confidence to WHO-like grade
            grade = "TBProfiler_literature"
            if confidence:
                if "Assoc w R" in confidence and "Not" not in confidence:
                    if "Interim" in confidence:
                        grade = "2) Assoc w R - Interim"
                    else:
                        grade = "1) Assoc w R"
                elif "Uncertain" in confidence:
                    grade = "3) Uncertain significance"
                elif "Not assoc" in confidence:
                    if "Interim" in confidence:
                        grade = "4) Not assoc w R - Interim"
                    else:
                        grade = "5) Not assoc w R"

            # Short mutation label
            aa_m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)", mutation)
            if aa_m:
                r1 = AA3TO1.get(aa_m.group(1), aa_m.group(1))
                a1 = AA3TO1.get(aa_m.group(3), aa_m.group(3))
                short_mut = f"{r1}{aa_m.group(2)}{a1}"
            else:
                short_mut = mutation

            marker_name = f"{gene_name}_{short_mut}"
            resistance = f"{drug}-R"

            new_markers.append({
                "pos": gpos,
                "ref": ref_base,
                "alt": alt_base,
                "drug": drug,
                "resistance": resistance,
                "marker_name": marker_name,
                "grade": grade,
                "gene": gene_name,
                "mutation": short_mut,
            })
            resolved += 1

    print(f"Results:")
    print(f"  Resolved to genomic coords: {resolved}")
    print(f"  Skipped (gene not in GFF):  {skipped_no_gene}")
    print(f"  Skipped (can't resolve):    {skipped_no_resolve}")
    print(f"  Skipped (duplicate):        {skipped_duplicate}")
    print()

    # Append to existing markers file
    if new_markers:
        with open(OUTPUT_MARKERS, "a") as f:
            for m in new_markers:
                f.write(f"{m['pos']}\t{m['ref']}\t{m['alt']}\t"
                        f"{m['drug']}\t{m['resistance']}\t{m['marker_name']}\t"
                        f"{m['grade']}\t{m['gene']}\t{m['mutation']}\n")
        print(f"Appended {len(new_markers)} new markers to {OUTPUT_MARKERS}")

        # Summary by drug
        drug_counts = defaultdict(int)
        for m in new_markers:
            drug_counts[m["drug"]] += 1
        print("\nNew markers per drug:")
        for d, n in sorted(drug_counts.items(), key=lambda x: -x[1]):
            print(f"  {d:12s} {n:>4}")
    else:
        print("No new markers to add.")


if __name__ == "__main__":
    main()
