"use strict";

/**
 * Validates a marker TSV preview pasted into a discussion body.
 *
 * The rules mirror docs/marker_format.md exactly. They are deliberately kept in
 * one pure function with no GitHub API access so they can be tested directly:
 * see validate-marker-tsv.test.js.
 */

// Default --kmer-size. An allele must fit strictly inside the window.
const DEFAULT_K = 31;
// classify additionally requires max(len(REF), len(ALT)) <= k - 2 * min_flank_bases,
// which is 11 at the default --min-flank-bases of 10.
const CLASSIFY_MAX_ALLELE = 11;

// How many individual row problems to list before collapsing the rest.
const MAX_LISTED = 8;

const DOCS = "https://pathogenomics-lab.github.io/pathotypr/marker_format/";

// Hidden marker so an edited discussion updates the existing comment instead of
// posting another one underneath it.
const COMMENT_MARKER = "<!-- pathotypr-marker-validation -->";

function extractPreview(body) {
  const fenced = body.match(/```(?:text|tsv|tab|plain)?[ \t]*\r?\n([\s\S]*?)```/);
  return fenced ? fenced[1].replace(/\r\n/g, "\n") : null;
}

function looksLikeHeader(line) {
  const l = line.toLowerCase();
  return l.includes("pos") || l.includes("lineage") || l.includes("ref");
}

function validateMarkerTsv(body) {
  const errors = [];
  const warnings = [];
  const notes = [];

  const preview = extractPreview(body || "");
  if (preview === null) {
    errors.push(
      "No TSV preview found. Paste your first rows inside a fenced code block so they can be checked."
    );
    return { errors, warnings, notes, stats: null };
  }

  const lines = preview.split("\n");
  const stats = {
    linesSeen: 0,
    comments: 0,
    headers: 0,
    dataRows: 0,
    validRows: 0,
    snps: 0,
    mnvs: 0,
    indels: 0,
    annotated: 0,
    maxDepth: 0,
    lineages: new Set(),
    duplicates: 0,
  };

  const rowErrors = [];
  const rowWarnings = [];
  const seen = new Set();
  let sawTab = false;
  let sawAnyContent = false;
  let headerConsumed = false;
  let headerCols = null;

  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    const lineNo = i + 1;

    if (line.trim() === "") continue;
    stats.linesSeen++;
    sawAnyContent = true;

    if (line.startsWith("#")) {
      stats.comments++;
      headerConsumed = true;
      if (line.includes("\t")) sawTab = true;
      if (headerCols === null) headerCols = line.replace(/^#/, "").split("\t");
      continue;
    }

    if (line.includes("\t")) sawTab = true;

    const cols = line.split("\t");
    const posRaw = cols[0].trim();

    // A first field that is not an integer is skipped by both parsers. On the
    // first such line that reads like column names, that is just a header.
    if (!/^\d+$/.test(posRaw)) {
      if (!headerConsumed && looksLikeHeader(line)) {
        stats.headers++;
        headerConsumed = true;
        if (headerCols === null) headerCols = cols;
        continue;
      }
      rowErrors.push(
        `line ${lineNo}: position "${posRaw}" is not an integer, so this row would be skipped`
      );
      continue;
    }
    headerConsumed = true;

    const pos = parseInt(posRaw, 10);
    if (pos < 1) {
      rowErrors.push(`line ${lineNo}: position ${pos} is not 1-based`);
      continue;
    }

    stats.dataRows++;

    if (cols.length < 4) {
      rowErrors.push(
        `line ${lineNo}: ${cols.length} columns, at least 4 are required (position, REF, ALT and one lineage level)`
      );
      continue;
    }

    const ref = cols[1].trim();
    const alt = cols[2].trim();

    if (ref === "") {
      rowErrors.push(`line ${lineNo}: empty REF allele`);
      continue;
    }
    if (alt === "") {
      rowErrors.push(`line ${lineNo}: empty ALT allele`);
      continue;
    }
    if (!/^[ACGT]+$/i.test(ref)) {
      rowErrors.push(
        `line ${lineNo}: REF "${ref}" contains characters other than A, C, G and T`
      );
      continue;
    }
    if (!/^[ACGT]+$/i.test(alt)) {
      rowErrors.push(
        `line ${lineNo}: ALT "${alt}" contains characters other than A, C, G and T`
      );
      continue;
    }

    // Lineage levels run from column 4 until the first empty cell. Whatever
    // follows that empty cell is annotation: gene, then mutation.
    const levels = [];
    let cut = cols.length;
    for (let c = 3; c < cols.length; c++) {
      if (cols[c].trim() === "") {
        cut = c;
        break;
      }
      levels.push(cols[c].trim());
    }

    if (levels.length === 0) {
      rowErrors.push(`line ${lineNo}: no non-empty lineage level`);
      continue;
    }

    const annotations = cols.slice(cut + 1).filter((c) => c.trim() !== "");
    if (annotations.length > 0) stats.annotated++;

    for (const level of levels) {
      if (level.includes(";")) {
        rowWarnings.push(
          `line ${lineNo}: "${level}" packs a hierarchy into one cell; classify tolerates this but split-fastq reads it as a single literal level, so the two would disagree`
        );
      }
      stats.lineages.add(level);
    }
    stats.maxDepth = Math.max(stats.maxDepth, levels.length);

    const longest = Math.max(ref.length, alt.length);
    if (longest >= DEFAULT_K) {
      rowErrors.push(
        `line ${lineNo}: longest allele is ${longest} bp, which does not fit inside the default k-mer window of ${DEFAULT_K}`
      );
      continue;
    }
    if (longest > CLASSIFY_MAX_ALLELE) {
      rowWarnings.push(
        `line ${lineNo}: longest allele is ${longest} bp, above the ${CLASSIFY_MAX_ALLELE} bp that classify allows at default --min-flank-bases`
      );
    }

    if (ref.length !== alt.length) {
      stats.indels++;
    } else if (ref.length > 1) {
      stats.mnvs++;
    } else {
      stats.snps++;
    }

    const key = `${pos}|${ref.toUpperCase()}|${alt.toUpperCase()}`;
    if (seen.has(key)) {
      stats.duplicates++;
    } else {
      seen.add(key);
    }

    stats.validRows++;
  }

  if (!sawAnyContent) {
    errors.push("The TSV preview is empty.");
    return { errors, warnings, notes, stats: null };
  }

  if (!sawTab) {
    errors.push(
      "The preview contains no tab characters. The marker format is tab-separated, and a space-separated file will not parse. Copy the rows straight out of your file rather than retyping them."
    );
  }

  // The mistake docs/marker_format.md calls out as the most common one: gene and
  // mutation columns are only recognised when an empty cell separates them from
  // the last lineage level. Without it they are silently swallowed into the
  // hierarchy, and nothing in the run reports it.
  if (headerCols && stats.validRows > 0 && stats.annotated === 0) {
    const named = headerCols.map((c) => c.trim().toLowerCase());
    if (named.includes("gene") || named.includes("mutation")) {
      warnings.push(
        "The header names a gene or mutation column, but no row has the empty separator column that marks where the lineage hierarchy ends. Without it those cells are read as extra lineage levels and the annotations are never set."
      );
    }
  }

  // Aggregate rather than repeating the same sentence per row.
  if (stats.indels > 0) {
    warnings.push(
      `${stats.indels} row(s) have REF and ALT of different lengths. classify accepts these; split-fastq skips them.`
    );
  }
  if (stats.duplicates > 0) {
    warnings.push(
      `${stats.duplicates} row(s) repeat a position with the same REF and ALT.`
    );
  }

  const collapse = (list, label) => {
    if (list.length <= MAX_LISTED) return list;
    const head = list.slice(0, MAX_LISTED);
    head.push(`and ${list.length - MAX_LISTED} more ${label}`);
    return head;
  };

  errors.push(...collapse(rowErrors, "errors"));
  warnings.push(...collapse(rowWarnings, "warnings"));

  if (stats.validRows === 0 && rowErrors.length === 0) {
    errors.push("No marker rows were found in the preview.");
  }

  return { errors, warnings, notes, stats };
}

function renderComment(result, scopeNote) {
  const { errors, warnings, stats } = result;
  const out = [COMMENT_MARKER];

  if (errors.length === 0) {
    out.push("### ✅ The TSV preview parses cleanly");
  } else {
    out.push("### ❌ The TSV preview has problems");
  }
  out.push("");

  if (errors.length > 0) {
    out.push("**Errors**");
    out.push("");
    for (const e of errors) out.push(`- ${e}`);
    out.push("");
  }

  if (warnings.length > 0) {
    out.push("**Worth checking**");
    out.push("");
    for (const w of warnings) out.push(`- ${w}`);
    out.push("");
  }

  if (stats) {
    const lineages = [...stats.lineages].sort();
    const shown =
      lineages.length > 12
        ? `${lineages.slice(0, 12).join(", ")} and ${lineages.length - 12} more`
        : lineages.join(", ");
    out.push("| Metric | Value |");
    out.push("|---|---|");
    out.push(`| Rows that parsed | ${stats.validRows} of ${stats.dataRows} |`);
    out.push(`| SNPs | ${stats.snps} |`);
    out.push(`| MNVs | ${stats.mnvs} |`);
    out.push(`| Indels | ${stats.indels} |`);
    out.push(`| Rows with gene and mutation | ${stats.annotated} |`);
    out.push(`| Deepest lineage hierarchy | ${stats.maxDepth} level(s) |`);
    out.push(`| Distinct lineage levels | ${stats.lineages.size}${shown ? ` (${shown})` : ""} |`);
    out.push("");
  }

  if (scopeNote) {
    out.push(`> ${scopeNote}`);
    out.push("");
  }

  if (errors.length === 0) {
    out.push(
      "This checks the preview only. A maintainer will review the full marker set, so please make sure the complete file is attached or linked above."
    );
  } else {
    out.push(
      "Edit the discussion to run these checks again."
    );
  }

  out.push("");
  out.push("---");
  out.push(`*Automated check against the [marker format reference](${DOCS}).*`);

  return out.join("\n");
}

module.exports = {
  validateMarkerTsv,
  renderComment,
  extractPreview,
  DOCS,
  COMMENT_MARKER,
};
