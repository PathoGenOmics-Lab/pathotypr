// ============================================================================
// Drug-resistance profile, lineage levels and allele-fraction insights
// ============================================================================
//
// The genotyping output is the same shape whether the marker panel describes
// lineages or drug resistance: every marker carries a `lineage_path` built from
// the marker TSV columns after REF/ALT, joined with ';'.
//
// For the WHO drug-resistance panel those columns are
//   drug ; resistance ; marker_name ; WHO grade ; gene ; mutation
// e.g. "RIF;RIF-R;rpoB_S450L;1) Assoc w R;rpoB;S450L"
//
// so the path is not a lineage at all. This module recognises that shape and
// turns it into a clinical resistance profile, keeps the genuine lineage paths
// as a proper hierarchy, and summarises allele fractions.

import { escapeHtml } from './utils.js';

// ---------------------------------------------------------------------------
// WHO confidence grading
// ---------------------------------------------------------------------------

/** Verdict a single marker contributes, derived from its WHO grade. */
export const GRADE = {
  RESISTANT: 'resistant',   // WHO 1) and 2) — associated with resistance
  UNCERTAIN: 'uncertain',   // WHO 3) — uncertain significance
  NOT_ASSOC: 'not_assoc',   // WHO 4) and 5) — not associated with resistance
  UNGRADED: 'ungraded'      // Unknown / literature-derived
};

const GRADE_RANK = {
  [GRADE.RESISTANT]: 3,
  [GRADE.UNCERTAIN]: 2,
  [GRADE.UNGRADED]: 1,
  [GRADE.NOT_ASSOC]: 0
};

const GRADE_LABEL = {
  [GRADE.RESISTANT]: 'Resistant',
  [GRADE.UNCERTAIN]: 'Uncertain',
  [GRADE.NOT_ASSOC]: 'Not associated',
  [GRADE.UNGRADED]: 'Ungraded marker'
};

/**
 * Classify a WHO grade string into a verdict.
 * Accepts the catalogue's "N) text" form as well as bare text.
 */
export function classifyGrade(rawGrade) {
  const grade = String(rawGrade || '').trim();
  if (!grade) return GRADE.UNGRADED;
  const leading = grade.match(/^\s*([1-5])\s*\)/);
  if (leading) {
    const n = Number(leading[1]);
    if (n === 1 || n === 2) return GRADE.RESISTANT;
    if (n === 3) return GRADE.UNCERTAIN;
    return GRADE.NOT_ASSOC;
  }
  const lower = grade.toLowerCase();
  if (lower.includes('not assoc')) return GRADE.NOT_ASSOC;
  if (lower.includes('assoc')) return GRADE.RESISTANT;
  if (lower.includes('uncertain')) return GRADE.UNCERTAIN;
  return GRADE.UNGRADED;
}

// ---------------------------------------------------------------------------
// Drugs: labels and clinical ordering
// ---------------------------------------------------------------------------

const DRUG_INFO = {
  RIF:     { label: 'Rifampicin',      group: 'First line' },
  INH:     { label: 'Isoniazid',       group: 'First line' },
  PZA:     { label: 'Pyrazinamide',    group: 'First line' },
  EMB:     { label: 'Ethambutol',      group: 'First line' },
  FQ:      { label: 'Fluoroquinolones', group: 'Group A' },
  BDQ:     { label: 'Bedaquiline',     group: 'Group A' },
  LZD:     { label: 'Linezolid',       group: 'Group A' },
  BDQ_CFZ: { label: 'Bedaquiline / Clofazimine', group: 'Group A/B' },
  CFZ:     { label: 'Clofazimine',     group: 'Group B' },
  DLM:     { label: 'Delamanid',       group: 'Group C' },
  ETH:     { label: 'Ethionamide',     group: 'Group C' },
  AMK:     { label: 'Amikacin',        group: 'Injectables' },
  KAN:     { label: 'Kanamycin',       group: 'Injectables' },
  CAP:     { label: 'Capreomycin',     group: 'Injectables' },
  STR:     { label: 'Streptomycin',    group: 'Injectables' },
  PAS:     { label: 'Para-aminosalicylic acid', group: 'Group C' },
  OTHER:   { label: 'Other',           group: 'Other' }
};

// Clinical reading order: first line, then Group A, companions, injectables.
const DRUG_ORDER = [
  'RIF', 'INH', 'PZA', 'EMB',
  'FQ', 'BDQ', 'LZD', 'BDQ_CFZ', 'CFZ',
  'DLM', 'ETH', 'PAS',
  'AMK', 'KAN', 'CAP', 'STR',
  'OTHER'
];

export function drugLabel(code) {
  return DRUG_INFO[code]?.label || code;
}

export function drugGroup(code) {
  return DRUG_INFO[code]?.group || 'Other';
}

function drugRank(code) {
  const idx = DRUG_ORDER.indexOf(code);
  return idx === -1 ? DRUG_ORDER.length : idx;
}

// ---------------------------------------------------------------------------
// Marker path parsing
// ---------------------------------------------------------------------------

/** A grade cell looks like "1) Assoc w R", "Unknown" or "TBProfiler_literature". */
function looksLikeGrade(value) {
  const v = String(value || '').trim();
  if (!v) return false;
  return /^[1-5]\s*\)/.test(v) || /^unknown$/i.test(v) || /literature/i.test(v);
}

/**
 * Parse one marker's `lineage_path` into its resistance components.
 *
 * Returns `null` when the path does not have the drug-resistance shape, which
 * is how a genuine lineage path (L2;L2.2;L2.2.1) is told apart.
 */
export function parseDrPath(lineagePath) {
  const parts = String(lineagePath || '').split(';').map(p => p.trim());
  if (parts.length < 4) return null;

  const [drug, resistance, markerName, grade, gene, mutation] = parts;
  if (!drug) return null;
  // The grade column is the reliable signal: no lineage level ever looks like
  // "1) Assoc w R".
  if (!looksLikeGrade(grade)) return null;

  return {
    drug: drug.toUpperCase(),
    resistance: resistance || '',
    markerName: markerName || '',
    grade: grade || '',
    verdict: classifyGrade(grade),
    gene: gene || '',
    mutation: mutation || markerName || ''
  };
}

/**
 * Decide whether a set of marker paths comes from a drug-resistance panel.
 * A single stray match is not enough — require a clear majority so a lineage
 * panel with an odd annotation is never mistaken for a DR panel.
 */
export function isDrPanel(lineagePaths) {
  const paths = (lineagePaths || []).filter(Boolean);
  if (paths.length === 0) return false;
  let drLike = 0;
  for (const p of paths) if (parseDrPath(p)) drLike += 1;
  return drLike >= Math.max(1, Math.ceil(paths.length * 0.5));
}

// ---------------------------------------------------------------------------
// Resistance profile
// ---------------------------------------------------------------------------

/**
 * Build the per-drug resistance profile.
 *
 * `records` are the called variants: `{ lineagePath, altFraction, coverage }`.
 * Every row in the genotyping output is already a *called* variant (the core
 * only writes markers passing --min-depth and --min-alt-percent), so presence
 * in the input means the mutation was detected.
 *
 * A drug's verdict is the strongest evidence found for it: one WHO grade 1/2
 * mutation makes it resistant regardless of how many uncertain ones there are.
 */
export function buildDrProfile(records) {
  const byDrug = new Map();

  for (const rec of records || []) {
    const parsed = parseDrPath(rec?.lineagePath);
    if (!parsed) continue;

    if (!byDrug.has(parsed.drug)) {
      byDrug.set(parsed.drug, { drug: parsed.drug, mutations: [] });
    }
    byDrug.get(parsed.drug).mutations.push({
      ...parsed,
      altFraction: Number.isFinite(rec?.altFraction) ? rec.altFraction : null,
      coverage: Number.isFinite(rec?.coverage) ? rec.coverage : null
    });
  }

  const profile = [];
  for (const entry of byDrug.values()) {
    // Strongest first, then by depth so the best-supported call leads.
    entry.mutations.sort((a, b) => {
      const rank = GRADE_RANK[b.verdict] - GRADE_RANK[a.verdict];
      if (rank !== 0) return rank;
      return (b.coverage || 0) - (a.coverage || 0);
    });
    const verdict = entry.mutations[0]?.verdict || GRADE.NOT_ASSOC;
    profile.push({
      drug: entry.drug,
      label: drugLabel(entry.drug),
      group: drugGroup(entry.drug),
      verdict,
      verdictLabel: GRADE_LABEL[verdict],
      mutations: entry.mutations,
      // Only grade 1/2 mutations actually drive a resistance call.
      resistantCount: entry.mutations.filter(m => m.verdict === GRADE.RESISTANT).length
    });
  }

  profile.sort((a, b) => drugRank(a.drug) - drugRank(b.drug) || a.drug.localeCompare(b.drug));
  return profile;
}

/**
 * One-line verdict for a sample, for the KPI header and the single-sample call.
 * Without this the header would present a raw marker path
 * ("FQ;FQ-R;gyrA_D94G;1) Assoc w R;gyrA;D94G") as if it were a lineage call.
 */
export function summariseDrProfile(profile) {
  const list = profile || [];
  const resistant = list.filter(p => p.verdict === GRADE.RESISTANT).map(p => p.drug);
  if (resistant.length > 0) return `Resistant: ${resistant.join(', ')}`;
  const uncertain = list.filter(p => p.verdict === GRADE.UNCERTAIN).map(p => p.drug);
  if (uncertain.length > 0) return `No resistance markers · uncertain: ${uncertain.join(', ')}`;
  return 'No resistance markers';
}

/** Label used when tallying a resistance run per drug. */
export function drCallLabel(entry) {
  return `${entry.drug} · ${entry.verdictLabel}`;
}

// ---------------------------------------------------------------------------
// Lineage levels
// ---------------------------------------------------------------------------

/**
 * Turn `lineage:count` entries into the nested hierarchy they describe.
 *
 * The backend reports every prefix of a path (L2, L2;L2.2, L2;L2.2;L2.2.1), so
 * a flat list hides the structure. This rebuilds the tree and keeps each
 * level's own support count.
 */
export function buildLineageLevels(entries) {
  const root = { name: '', label: '', count: 0, depth: 0, children: new Map() };

  for (const entry of entries || []) {
    const path = String(entry?.lineage || '').trim();
    if (!path) continue;
    const count = Number.isFinite(entry?.count) ? entry.count : 0;
    const parts = path.split(';').map(p => p.trim()).filter(Boolean);

    let node = root;
    parts.forEach((part, depth) => {
      if (!node.children.has(part)) {
        node.children.set(part, {
          name: part,
          label: part,
          count: 0,
          depth: depth + 1,
          children: new Map()
        });
      }
      node = node.children.get(part);
    });
    // The count belongs to the deepest level of this entry's own path.
    node.count = Math.max(node.count, count);
  }

  const toArray = (node) => {
    const children = [...node.children.values()]
      .map(toArray)
      .sort((a, b) => b.count - a.count || a.name.localeCompare(b.name));
    return { ...node, children };
  };

  return toArray(root).children;
}

// ---------------------------------------------------------------------------
// Mixed infections: fixed markers on incompatible lineage branches
// ---------------------------------------------------------------------------

/**
 * Two lineage paths are compatible when one is an ancestor of the other (or
 * they are equal): `L2` and `L2;L2.2` describe one strain seen at different
 * depths. Paths that diverge — `L2` vs `L4`, or `L2;L2.2` vs `L2;L2.1` — cannot
 * both come from a single clone.
 */
export function isCompatiblePath(a, b) {
  const pa = String(a || '').split(';').map(s => s.trim()).filter(Boolean);
  const pb = String(b || '').split(';').map(s => s.trim()).filter(Boolean);
  const n = Math.min(pa.length, pb.length);
  if (n === 0) return false;
  for (let i = 0; i < n; i += 1) {
    if (pa[i] !== pb[i]) return false;
  }
  return true;
}

/**
 * Resolve the observed lineage paths into the deepest branches they support,
 * and measure how much evidence is unique to each.
 *
 * Because the core only writes variants at or above `--min-alt-percent`, every
 * marker in the file is effectively *fixed*. A single clone can only fix
 * markers along one root-to-tip path, so fixed markers on two divergent
 * branches are the signature of a mixed infection — and unlike intermediate
 * allele fractions, this signal survives the default 95% threshold.
 *
 * `total` counts every marker compatible with the branch (including the
 * ancestral ones it shares with its siblings); `specific` counts only markers
 * that no other branch can explain, which is the evidence that actually
 * establishes a second strain.
 */
export function buildLineageBranches(entries) {
  const counts = new Map();
  for (const entry of entries || []) {
    const path = String(entry?.lineage || '').trim();
    if (!path) continue;
    const count = Number.isFinite(entry?.count) ? entry.count : 0;
    counts.set(path, (counts.get(path) || 0) + count);
  }
  if (counts.size === 0) return [];

  const all = [...counts.keys()];
  // A branch tip is a path no other observed path extends.
  const tips = all.filter(p => !all.some(q => q !== p && q.startsWith(`${p};`)));

  const branches = tips.map(branch => {
    let total = 0;
    let specific = 0;
    for (const [path, count] of counts) {
      if (!isCompatiblePath(path, branch)) continue;
      total += count;
      const explainedElsewhere = tips.some(other => other !== branch && isCompatiblePath(path, other));
      if (!explainedElsewhere) specific += count;
    }
    return { branch, total, specific };
  });

  branches.sort((a, b) => b.specific - a.specific || b.total - a.total || a.branch.localeCompare(b.branch));
  return branches;
}

/**
 * Decide whether the observed branches look like more than one strain.
 *
 * Requires at least two branches carrying evidence of their own: a branch with
 * no specific markers is just an ancestor of another and proves nothing.
 * `minSpecific` guards against a single stray marker calling a mixture.
 */
export function detectMixedLineages(branches, minSpecific = 2) {
  const supported = (branches || []).filter(b => b.specific >= minSpecific);
  const totalSpecific = supported.reduce((sum, b) => sum + b.specific, 0);
  return {
    mixed: supported.length >= 2,
    branches: supported,
    totalSpecific,
    shares: supported.map(b => ({
      ...b,
      share: totalSpecific > 0 ? (b.specific / totalSpecific) * 100 : 0
    }))
  };
}

// ---------------------------------------------------------------------------
// Allele fractions
// ---------------------------------------------------------------------------

/**
 * Bin allele fractions (percentages) for the mixed-infection histogram.
 *
 * A clonal sample is strongly bimodal near 0% and 100%; markers sitting at
 * intermediate fractions are the signature of more than one strain.
 */
export function buildAlleleHistogram(fractions, binCount = 10) {
  const bins = Math.max(2, Math.floor(binCount));
  const width = 100 / bins;
  const out = Array.from({ length: bins }, (_, i) => ({
    from: i * width,
    to: (i + 1) * width,
    count: 0
  }));

  let total = 0;
  let intermediate = 0;
  for (const raw of fractions || []) {
    const value = Number(raw);
    if (!Number.isFinite(value)) continue;
    const clamped = Math.min(100, Math.max(0, value));
    const idx = Math.min(bins - 1, Math.floor(clamped / width));
    out[idx].count += 1;
    total += 1;
    // 20–80% is the band where a clonal sample should have very little weight.
    if (clamped >= 20 && clamped <= 80) intermediate += 1;
  }

  return {
    bins: out,
    total,
    intermediate,
    intermediateRatio: total > 0 ? intermediate / total : 0,
    maxCount: out.reduce((m, b) => Math.max(m, b.count), 0)
  };
}

// ---------------------------------------------------------------------------
// Rendering
// ---------------------------------------------------------------------------

function formatPct(value) {
  if (!Number.isFinite(value)) return '—';
  return `${value.toFixed(value >= 10 ? 0 : 1)}%`;
}

function formatDepth(value) {
  if (!Number.isFinite(value) || value <= 0) return '—';
  return `${Math.round(value)}×`;
}

/** Render the per-drug resistance profile. */
export function renderDrProfileHtml(profile) {
  if (!Array.isArray(profile) || profile.length === 0) return '';

  const cards = profile.map(entry => {
    const top = entry.mutations[0];
    const detail = top
      ? `${escapeHtml(top.gene || top.markerName || '—')} ${escapeHtml(top.mutation || '')}`.trim()
      : '';
    const extra = entry.mutations.length > 1
      ? `<span class="dr-more">+${entry.mutations.length - 1} more</span>`
      : '';

    const rows = entry.mutations.map(m => `
      <tr class="dr-mut dr-mut--${m.verdict}">
        <td class="dr-mut-gene">${escapeHtml(m.gene || '—')}</td>
        <td class="dr-mut-change">${escapeHtml(m.mutation || m.markerName || '—')}</td>
        <td class="dr-mut-grade">${escapeHtml(m.grade || '—')}</td>
        <td class="dr-mut-depth">${formatDepth(m.coverage)}</td>
        <td class="dr-mut-frac">${formatPct(m.altFraction)}</td>
      </tr>
    `).join('');

    return `
      <details class="dr-card dr-card--${entry.verdict}">
        <summary class="dr-card-head">
          <span class="dr-drug" title="${escapeHtml(entry.label)}">${escapeHtml(entry.drug)}</span>
          <span class="dr-verdict">${escapeHtml(entry.verdictLabel)}</span>
          <span class="dr-top">${detail}</span>
          ${extra}
        </summary>
        <div class="dr-card-body">
          <table class="dr-mut-table">
            <thead>
              <tr><th>Gene</th><th>Change</th><th>WHO grade</th><th>Depth</th><th>Allele fraction</th></tr>
            </thead>
            <tbody>${rows}</tbody>
          </table>
        </div>
      </details>
    `;
  }).join('');

  const resistant = profile.filter(p => p.verdict === GRADE.RESISTANT).map(p => p.drug);
  const headline = resistant.length > 0
    ? `Resistance markers found for ${resistant.length} drug${resistant.length > 1 ? 's' : ''}: ${escapeHtml(resistant.join(', '))}`
    : 'No WHO grade 1/2 resistance markers detected';

  return `
    <section class="viz-dr-profile">
      <header class="viz-dr-header">
        <h4>Drug resistance profile</h4>
        <p class="viz-dr-headline">${headline}</p>
      </header>
      <div class="dr-cards">${cards}</div>
      <p class="viz-dr-note">
        Verdicts follow the WHO catalogue grading: grades 1 and 2 are associated with
        resistance, grade 3 is of uncertain significance, grades 4 and 5 are not
        associated. Only drugs with at least one detected marker are listed.
      </p>
    </section>
  `;
}

/** Render the lineage hierarchy as nested levels. */
export function renderLineageLevelsHtml(levels) {
  if (!Array.isArray(levels) || levels.length === 0) return '';

  const renderNodes = (nodes) => nodes.map(node => {
    const children = node.children?.length ? renderNodes(node.children) : '';
    const count = node.count > 0 ? `<span class="lvl-count">${node.count}</span>` : '';
    return `
      <li class="lvl-node lvl-depth-${node.depth}">
        <span class="lvl-name">${escapeHtml(node.label)}</span>${count}
        ${children ? `<ul class="lvl-children">${children}</ul>` : ''}
      </li>
    `;
  }).join('');

  return `
    <section class="viz-lineage-levels">
      <header><h4>Lineage levels</h4></header>
      <ul class="lvl-tree">${renderNodes(levels)}</ul>
    </section>
  `;
}

/**
 * Render the lineage-branch support, which is the primary mixed-infection
 * signal: fixed markers that belong to divergent branches.
 */
export function renderMixedLineagesHtml(detection, options = {}) {
  const all = detection?.shares || [];
  if (all.length === 0) return '';

  // A genuine mixture is two or three strains. A long tail of branches usually
  // means marker noise rather than that many co-infecting strains, so show the
  // best-supported ones and account for the rest instead of listing them all.
  const maxBranches = Math.max(2, options.maxBranches || 6);
  const shares = all.slice(0, maxBranches);
  const hidden = all.length - shares.length;

  const rows = shares.map(entry => `
    <li class="mx-branch${detection.mixed ? '' : ' mx-branch--single'}">
      <div class="mx-branch-head">
        <span class="mx-branch-name">${escapeHtml(entry.branch)}</span>
        <span class="mx-branch-share">${entry.share.toFixed(0)}%</span>
      </div>
      <div class="mx-bar"><div class="mx-bar-fill" style="width:${entry.share.toFixed(1)}%"></div></div>
      <div class="mx-branch-meta">
        <strong>${entry.specific}</strong> marker${entry.specific === 1 ? '' : 's'} unique to this branch
        <span class="mx-branch-total">· ${entry.total} compatible in total</span>
      </div>
    </li>
  `).join('');

  const verdict = detection.mixed
    ? `<p class="mx-verdict mx-verdict--mixed">
         Mixed infection likely — fixed markers support ${all.length} incompatible lineages
       </p>`
    : `<p class="mx-verdict mx-verdict--clonal">
         Consistent with a single strain — all fixed markers lie on one lineage path
       </p>`;

  const overflow = hidden > 0
    ? `<li class="mx-more">+${hidden} further branch${hidden === 1 ? '' : 'es'} with at least
         two unique markers, not shown.</li>`
    : '';

  const manyBranches = all.length > 3
    ? ` Note that ${all.length} branches is more than a co-infection usually produces: a genuine
        mixture is normally two or three strains, so check the marker panel and depth before
        reading this as that many strains.`
    : '';

  const note = detection.mixed
    ? `A single clone can only fix markers along one root-to-tip path, so markers fixed on
       divergent branches point to more than one strain in the sample. Only markers that no
       other branch can explain are counted, so shared ancestral markers never inflate the
       call.${manyBranches}`
    : `Every observed marker lies on one root-to-tip path, so no second strain is implied.`;

  return `
    <section class="viz-mixed-lineages">
      <header><h4>Lineage support</h4></header>
      ${verdict}
      <ul class="mx-branches">${rows}${overflow}</ul>
      <p class="mx-note">${note}</p>
    </section>
  `;
}

/**
 * Render the allele-fraction histogram.
 *
 * `minAltPercent` is the threshold the run used. The core only writes variants
 * at or above it, so when it is high the intermediate fractions that reveal a
 * mixed infection were filtered out before reaching this file — say so rather
 * than let an empty middle read as "clonal".
 */
export function renderAlleleHistogramHtml(hist, options = {}) {
  if (!hist || hist.total === 0) return '';
  const { minAltPercent = null } = options;

  const bars = hist.bins.map(bin => {
    const pct = hist.maxCount > 0 ? (bin.count / hist.maxCount) * 100 : 0;
    const mid = (bin.from + bin.to) / 2;
    const isIntermediate = mid >= 20 && mid <= 80;
    return `
      <div class="af-bin${isIntermediate ? ' af-bin--intermediate' : ''}"
           title="${bin.from.toFixed(0)}–${bin.to.toFixed(0)}%: ${bin.count} marker(s)">
        ${bin.count > 0 ? `<div class="af-bar" style="height:${pct.toFixed(1)}%"></div>` : ''}
        <span class="af-tick">${bin.from.toFixed(0)}</span>
      </div>
    `;
  }).join('');

  const ratio = (hist.intermediateRatio * 100).toFixed(0);
  let note;
  if (minAltPercent !== null && Number(minAltPercent) >= 20 && hist.intermediate === 0) {
    note = `No variant below ${escapeHtml(String(Math.round(Number(minAltPercent))))}% alternate allele
            fraction appears in this file, so it was filtered by <code>--min-alt-percent</code> before
            being written. Intermediate fractions are the signature of a mixed infection, so this view
            cannot rule one out — re-run with a lower <code>--min-alt-percent</code> (for example 10)
            to look for mixtures.`;
  } else if (hist.intermediate > 0) {
    note = `${hist.intermediate} of ${hist.total} markers (${ratio}%) fall between 20% and 80%.
            A clonal sample is strongly bimodal near 0% and 100%, so weight in the middle
            suggests more than one strain.`;
  } else {
    note = `All ${hist.total} markers sit outside the 20–80% band, which is what a single
            clonal strain looks like.`;
  }

  return `
    <section class="viz-allele-fraction">
      <header><h4>Allele fraction distribution</h4></header>
      <div class="af-chart">${bars}</div>
      <div class="af-axis"><span>0%</span><span>alternate allele fraction</span><span>100%</span></div>
      <p class="af-note">${note}</p>
    </section>
  `;
}
