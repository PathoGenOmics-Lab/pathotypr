#!/usr/bin/env python3
"""
Port the TBProfiler-literature markers appended by resolve_tbdb_mutations.py from H37Rv
to inferred-ancestor coordinates.

convert_who_markers.py writes both coordinate frames, but resolve_tbdb_mutations.py only
appends to the H37Rv file, so without this step the ancestor marker file loses the 402
TBProfiler-exclusive entries. Positions are shared between the two references; only the
reference base differs, and a marker whose ancestral base already equals the alternate
allele is dropped, exactly as convert_who_markers.py does.

Run after resolve_tbdb_mutations.py:
    python3 tbdb_markers_to_ancestor.py <n_appended>
"""
import sys
from pathlib import Path

L4_MARKERS = "pathotypr_dr_markers_H37Rv.tsv"
ANC_MARKERS = "pathotypr_dr_markers_ancestor.tsv"
REF_ANC = "MTB_ancestor_reference.fasta"


def load_reference(path):
    with open(path) as f:
        return "".join(l.strip().upper() for l in f if not l.startswith(">"))


def main(n_appended):
    ref = load_reference(REF_ANC)
    print(f"Ancestor reference: {len(ref):,} bp")

    rows = Path(L4_MARKERS).read_text().rstrip("\n").split("\n")
    tail = rows[-n_appended:]
    print(f"Porting the last {len(tail)} rows of {L4_MARKERS}")

    already = set()
    for line in Path(ANC_MARKERS).read_text().rstrip("\n").split("\n"):
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) >= 4:
            already.add((p[0], p[1], p[2], p[3]))

    written = skipped_same = skipped_dup = changed = 0
    with open(ANC_MARKERS, "a") as out:
        for line in tail:
            p = line.split("\t")
            if len(p) < 9:
                continue
            pos, l4_ref, alt = p[0], p[1], p[2]
            anc_ref = ref[int(pos) - 1: int(pos) - 1 + len(l4_ref)] or l4_ref
            if anc_ref == alt:
                skipped_same += 1
                continue
            if (pos, anc_ref, alt, p[3]) in already:
                skipped_dup += 1
                continue
            if anc_ref != l4_ref:
                changed += 1
            out.write("\t".join([pos, anc_ref, alt] + p[3:]) + "\n")
            written += 1

    print(f"  written:                {written}")
    print(f"  ref base differed:      {changed}")
    print(f"  skipped (ref == alt):   {skipped_same}")
    print(f"  skipped (duplicate):    {skipped_dup}")


if __name__ == "__main__":
    main(int(sys.argv[1]))
