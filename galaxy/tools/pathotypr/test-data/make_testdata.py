#!/usr/bin/env python3
"""Generate the smallest test data that still exercises every pathotypr module.

Everything is deterministic: a fixed seed, and every derived file is built from
the same reference so the expected calls can be reasoned about rather than
copied out of a run.

Size matters here. tools-iuc carries this data forever, so the reference is a
couple of kb rather than a genome.
"""
import gzip
import os
import random
import sys

random.seed(20260816)

OUT = sys.argv[1] if len(sys.argv) > 1 else "test-data"
os.makedirs(OUT, exist_ok=True)

REF_LEN = 2000
K = 31            # default --kmer-size
FLANK = 10        # default --min-flank-bases
READ_LEN = 100
DEPTH = 20        # comfortably above the default --min-depth of 10


def w(name, text):
    with open(os.path.join(OUT, name), "w") as fh:
        fh.write(text)
    return os.path.join(OUT, name)


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# ---------------------------------------------------------------- reference
# Random sequence, but every k-mer window around a marker has to be unique or
# the diagnostic k-mers match in more than one place. Random 2 kb over a 31-mer
# window is unique with overwhelming probability; asserted below rather than
# assumed.
ref = "".join(random.choice("ACGT") for _ in range(REF_LEN))

# Marker positions are 1-based. Keep them clear of the ends so a centred 31-mer
# with 10 flanking bases always fits.
MARKERS = [
    # pos, alt, lineage levels, gene, mutation
    (200,  "L1",),
    (400,  "L2",),
    (600,  "L2", "L2.1"),
    (800,  "L2", "L2.1", "L2.1.1"),
    (1000, "L3",),
    (1200, "L3", "L3.2"),
]
ANNOTATED = [
    (1400, "RIF", "rpoB", "S450L"),
    (1600, "INH", "katG", "Ser315Thr"),
]


def alt_of(base):
    return {"A": "G", "G": "A", "C": "T", "T": "C"}[base]


rows = ["#position\tref\talt\tlevel1\tlevel2\tlevel3"]
marker_alts = {}
for m in MARKERS:
    pos, levels = m[0], list(m[1:])
    r = ref[pos - 1]
    a = alt_of(r)
    marker_alts[pos] = (r, a)
    rows.append("\t".join([str(pos), r, a] + levels))
for pos, lineage, gene, mut in ANNOTATED:
    r = ref[pos - 1]
    a = alt_of(r)
    marker_alts[pos] = (r, a)
    # The empty cell after the lineage is what makes gene and mutation be read
    # as annotations rather than as two more lineage levels.
    rows.append("\t".join([str(pos), r, a, lineage, "", gene, mut]))

# Every marker k-mer must be unique in the reference, otherwise a hit is
# ambiguous and the test data would be testing the wrong thing.
half = K // 2
for pos in marker_alts:
    kmer = ref[pos - 1 - half: pos + half]
    assert len(kmer) == K, f"marker {pos} too close to the end"
    assert ref.count(kmer) == 1 and ref.count(revcomp(kmer)) == 0, (
        f"marker k-mer at {pos} is not unique")
    assert pos - 1 - half >= FLANK and pos + half <= REF_LEN - FLANK

w("reference.fasta", ">MTB_test_ref synthetic 2 kb reference\n" +
  "\n".join(ref[i:i + 60] for i in range(0, len(ref), 60)) + "\n")
w("markers.tsv", "\n".join(rows) + "\n")


def genome_with(alt_positions, name, desc):
    """The reference with the ALT allele substituted at the given positions."""
    g = list(ref)
    for pos in alt_positions:
        g[pos - 1] = marker_alts[pos][1]
    seq = "".join(g)
    w(name, f">{desc}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    return seq


# sample1 carries the full L2 -> L2.1 -> L2.1.1 path plus the RIF marker.
s1 = genome_with([400, 600, 800, 1400], "sample1.fasta", "sample1")
# sample2 carries L3 -> L3.2 and the INH marker.
s2 = genome_with([1000, 1200, 1600], "sample2.fasta", "sample2")

# ------------------------------------------------------------------- FASTQ
def reads(seq, path, depth=DEPTH, read_len=READ_LEN, paired=True):
    """Tile the genome so every position, markers included, is covered `depth`
    times. Tiling rather than random sampling keeps the file small and the
    coverage guaranteed."""
    step = max(1, (read_len * 2) // depth)
    r1, r2 = [], []
    n = 0
    for start in range(0, len(seq) - read_len * 2, step):
        f = seq[start:start + read_len]
        rc = revcomp(seq[start + read_len: start + read_len * 2])
        if len(f) < read_len or len(rc) < read_len:
            continue
        n += 1
        q = "I" * read_len
        r1.append(f"@read{n}/1\n{f}\n+\n{q}\n")
        r2.append(f"@read{n}/2\n{rc}\n+\n{q}\n")
    if paired:
        with gzip.open(path + "_R1.fastq.gz", "wt") as fh:
            fh.write("".join(r1))
        with gzip.open(path + "_R2.fastq.gz", "wt") as fh:
            fh.write("".join(r2))
    return n


n1 = reads(s1, os.path.join(OUT, "sample1"))

# ---------------------------------------------------------------- training
# train takes the label from the first whitespace-separated token of the
# header. Several genomes per class, each a mutated copy, so the forest has
# something to separate.
def mutate(seq, n_changes, rng):
    g = list(seq)
    for _ in range(n_changes):
        i = rng.randrange(len(g))
        g[i] = alt_of(g[i])
    return "".join(g)


rng = random.Random(7)
train_records = []
for i in range(6):
    train_records.append((f"L2 train_L2_{i}", mutate(s1, 20, rng)))
for i in range(6):
    train_records.append((f"L3 train_L3_{i}", mutate(s2, 20, rng)))

w("training.fasta", "".join(
    f">{h}\n" + "\n".join(s[j:j + 60] for j in range(0, len(s), 60)) + "\n"
    for h, s in train_records))

# query set for predict: two unseen genomes, one per class
w("query.fasta", "".join(
    f">{h}\n" + "\n".join(s[j:j + 60] for j in range(0, len(s), 60)) + "\n"
    for h, s in [("query_a", mutate(s1, 20, rng)), ("query_b", mutate(s2, 20, rng))]))

# ------------------------------------------------------------- match refs
w("references.fasta",
  ">ref_A\n" + "\n".join(s1[i:i + 60] for i in range(0, len(s1), 60)) + "\n" +
  ">ref_B\n" + "\n".join(s2[i:i + 60] for i in range(0, len(s2), 60)) + "\n")

print(f"reference   {REF_LEN} bp")
print(f"markers     {len(marker_alts)} ({len(MARKERS)} lineage, {len(ANNOTATED)} annotated)")
print(f"reads       {n1} pairs, ~{n1 * READ_LEN * 2 / REF_LEN:.0f}x")
print(f"training    {len(train_records)} records, 2 classes")
print()
for f in sorted(os.listdir(OUT)):
    print(f"  {os.path.getsize(os.path.join(OUT, f)):>8} B  {f}")
