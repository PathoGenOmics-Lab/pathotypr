"use strict";

const test = require("node:test");
const assert = require("node:assert");
const { validateMarkerTsv, renderComment } = require("./validate-marker-tsv.js");

// Every case below is taken from docs/marker_format.md, so the validator and the
// documentation cannot drift apart without a test failing.
const fence = (tsv) => "some prose\n\n```text\n" + tsv + "\n```\n\nmore prose";

test("the documented lineage example parses cleanly", () => {
  const r = validateMarkerTsv(
    fence(
      [
        "#position\tref\talt\tlevel1\tlevel2\tlevel3",
        "615938\tA\tG\tL1",
        "1799921\tC\tG\tL2",
        "801959\tA\tG\tL2\tL2.2",
        "2831482\tC\tT\tL2\tL2.2\tL2.2.1",
      ].join("\n")
    )
  );
  assert.deepStrictEqual(r.errors, []);
  assert.strictEqual(r.stats.validRows, 4);
  assert.strictEqual(r.stats.snps, 4);
  assert.strictEqual(r.stats.maxDepth, 3);
  assert.deepStrictEqual([...r.stats.lineages].sort(), ["L1", "L2", "L2.2", "L2.2.1"]);
});

test("the documented annotation example is recognised as annotated", () => {
  const r = validateMarkerTsv(
    fence(
      [
        "#position\tref\talt\tlineage\t\tgene\tmutation",
        "761155\tC\tT\tRIF\t\trpoB\tS450L",
        "2155168\tG\tA\tINH\t\tkatG\tSer315Thr",
      ].join("\n")
    )
  );
  assert.deepStrictEqual(r.errors, []);
  assert.strictEqual(r.stats.annotated, 2);
  assert.strictEqual(r.stats.maxDepth, 1);
  assert.deepStrictEqual([...r.stats.lineages].sort(), ["INH", "RIF"]);
});

test("an equal-length MNV is an MNV and not a SNP", () => {
  const r = validateMarkerTsv(
    fence(["#position\tref\talt\tlineage", "3820\tAC\tGT\tL4"].join("\n"))
  );
  assert.deepStrictEqual(r.errors, []);
  assert.strictEqual(r.stats.mnvs, 1);
  assert.strictEqual(r.stats.snps, 0);
  assert.strictEqual(r.stats.indels, 0);
});

test("a multi-base allele is not rejected as an invalid base", () => {
  // The previous validator tested `['A','C','G','T'].includes(ref)`, which fails
  // for every MNV and indel the format explicitly supports.
  const r = validateMarkerTsv(
    fence(["#position\tref\talt\tlineage", "3820\tACGT\tGTCA\tL4"].join("\n"))
  );
  assert.deepStrictEqual(r.errors, []);
});

test("an indel parses but warns that split-fastq skips it", () => {
  const r = validateMarkerTsv(
    fence(["#position\tref\talt\tlineage", "4247431\tCC\tC\tL1"].join("\n"))
  );
  assert.deepStrictEqual(r.errors, []);
  assert.strictEqual(r.stats.indels, 1);
  assert.match(r.warnings.join(" "), /split-fastq skips them/);
});

test("each documented invalid row is reported", () => {
  const cases = [
    ["615938\tA\tG", /at least 4 are required/],
    ["foo\tA\tG\tL1", /not an integer/],
    ["615938\t\tG\tL1", /empty REF/],
    ["615938\tA\t\tL1", /empty ALT/],
    ["615938\tA\tG\t", /no non-empty lineage level/],
  ];
  for (const [row, pattern] of cases) {
    const r = validateMarkerTsv(fence("#position\tref\talt\tlineage\n" + row));
    assert.match(r.errors.join(" "), pattern, `row: ${JSON.stringify(row)}`);
  }
});

test("a non-ACGT allele is reported", () => {
  const r = validateMarkerTsv(
    fence(["#position\tref\talt\tlineage", "615938\tN\tG\tL1"].join("\n"))
  );
  assert.match(r.errors.join(" "), /other than A, C, G and T/);
});

test("annotations without the empty separator column are flagged", () => {
  // docs/marker_format.md calls this the mistake most people make: the trailing
  // cells are swallowed into the hierarchy and gene/mutation are never set.
  const r = validateMarkerTsv(
    fence(
      [
        "#position\tref\talt\tlineage\tgene\tmutation",
        "761155\tC\tT\tRIF\trpoB\tS450L",
      ].join("\n")
    )
  );
  assert.deepStrictEqual(r.errors, []);
  assert.match(r.warnings.join(" "), /empty separator column/);
  assert.strictEqual(r.stats.maxDepth, 3);
});

test("a space-separated preview is caught", () => {
  const r = validateMarkerTsv(
    fence(["position ref alt lineage", "615938 A G L1"].join("\n"))
  );
  assert.match(r.errors.join(" "), /no tab characters/);
});

test("a semicolon-packed hierarchy warns about the two engines disagreeing", () => {
  const r = validateMarkerTsv(
    fence(["#position\tref\talt\tlineage", "615938\tA\tG\tL2;L2.2;L2.2.1"].join("\n"))
  );
  assert.deepStrictEqual(r.errors, []);
  assert.match(r.warnings.join(" "), /split-fastq reads it as a single literal level/);
});

test("an allele longer than the k-mer window is an error", () => {
  const long = "A".repeat(35);
  const r = validateMarkerTsv(
    fence(["#position\tref\talt\tlineage", `615938\t${long}\t${long}\tL1`].join("\n"))
  );
  assert.match(r.errors.join(" "), /does not fit inside the default k-mer window/);
});

test("an allele over the classify flank limit warns but is not an error", () => {
  const mid = "A".repeat(15);
  const r = validateMarkerTsv(
    fence(["#position\tref\talt\tlineage", `615938\t${mid}\t${mid}\tL1`].join("\n"))
  );
  assert.deepStrictEqual(r.errors, []);
  assert.match(r.warnings.join(" "), /--min-flank-bases/);
});

test("a header without a leading hash is skipped, not reported as a bad row", () => {
  const r = validateMarkerTsv(
    fence(["position\tref\talt\tlineage", "615938\tA\tG\tL1"].join("\n"))
  );
  assert.deepStrictEqual(r.errors, []);
  assert.strictEqual(r.stats.headers, 1);
  assert.strictEqual(r.stats.validRows, 1);
});

test("duplicate rows are counted once and reported", () => {
  const r = validateMarkerTsv(
    fence(
      [
        "#position\tref\talt\tlineage",
        "615938\tA\tG\tL1",
        "615938\tA\tG\tL1",
      ].join("\n")
    )
  );
  assert.strictEqual(r.stats.duplicates, 1);
  assert.match(r.warnings.join(" "), /repeat a position/);
});

test("a body with no code fence is reported rather than crashing", () => {
  const r = validateMarkerTsv("I forgot to paste anything");
  assert.match(r.errors.join(" "), /No TSV preview found/);
  assert.strictEqual(r.stats, null);
});

test("carriage returns from a Windows paste do not break parsing", () => {
  const r = validateMarkerTsv(
    "```text\r\n#position\tref\talt\tlineage\r\n615938\tA\tG\tL1\r\n```"
  );
  assert.deepStrictEqual(r.errors, []);
  assert.strictEqual(r.stats.validRows, 1);
});

test("the rendered comment reflects success and failure", () => {
  const ok = renderComment(
    validateMarkerTsv(fence("#position\tref\talt\tlineage\n615938\tA\tG\tL1"))
  );
  assert.match(ok, /parses cleanly/);
  assert.doesNotMatch(ok, /\u2014/, "no em-dashes in user-facing output");

  const bad = renderComment(validateMarkerTsv("nothing here"));
  assert.match(bad, /has problems/);
  assert.doesNotMatch(bad, /\u2014/, "no em-dashes in user-facing output");
});
