// ============================================================================
// Batch charts: resistance matrix and lineage composition
// ============================================================================
//
// Two charts for a run with more than one sample, where a per-sample card list
// stops being readable:
//
//   * Resistance matrix — samples x drugs. Each cell is a *state*, so it uses
//     the fixed status palette, and because status colour must never carry
//     meaning alone every cell also has a glyph, a legend and a table view.
//   * Lineage composition — the share of branch-specific markers per sample, as
//     a stacked bar. That is identity, so it uses the categorical slots.
//
// Palettes are the validated defaults, checked with the data-viz validator
// against this app's own surfaces (#ffffff light, #192734 dark):
//   categorical 1-3  PASS both modes (light aqua is sub-3:1 -> labels + table)
//   status           used as documented; warning/serious are sub-3:1 on light
//                    by design, mitigated by the glyph + label pairing.

import { escapeHtml } from './utils.js';
import { parseDrPath, GRADE, drugLabel, buildDrProfile } from './dr-insights.js';
import { buildLineageBranches } from './dr-insights.js';

// Drug columns in clinical reading order; only those present are shown.
const DRUG_ORDER = [
  'RIF', 'INH', 'PZA', 'EMB',
  'FQ', 'BDQ', 'LZD', 'BDQ_CFZ', 'CFZ',
  'DLM', 'ETH', 'PAS',
  'AMK', 'KAN', 'CAP', 'STR',
  'OTHER'
];

// Status roles. The glyph is the secondary encoding that lets the matrix be
// read without colour at all.
const STATE = {
  [GRADE.RESISTANT]: { key: 'resistant', label: 'Resistant', glyph: '●', rank: 3 },
  [GRADE.UNCERTAIN]: { key: 'uncertain', label: 'Uncertain', glyph: '◐', rank: 2 },
  [GRADE.UNGRADED]:  { key: 'ungraded',  label: 'Ungraded',  glyph: '◍', rank: 1 },
  [GRADE.NOT_ASSOC]: { key: 'notassoc',  label: 'Not associated', glyph: '○', rank: 0 }
};
const EMPTY_STATE = { key: 'none', label: 'No marker detected', glyph: '·', rank: -1 };

/**
 * Build the samples x drugs matrix from the loaded rows.
 * Returns null unless this is a resistance panel with more than one sample —
 * for a single sample the per-drug cards are the better form than a one-row grid.
 */
export function buildResistanceMatrix(records) {
  const bySample = new Map();
  for (const rec of records || []) {
    if (!parseDrPath(rec?.lineagePath)) continue;
    const sample = rec.sample || 'Sample';
    if (!bySample.has(sample)) bySample.set(sample, []);
    bySample.get(sample).push(rec);
  }
  if (bySample.size < 2) return null;

  const samples = [...bySample.keys()].sort();
  const present = new Set();
  const cells = new Map();

  for (const sample of samples) {
    const profile = buildDrProfile(bySample.get(sample));
    for (const entry of profile) {
      present.add(entry.drug);
      cells.set(`${sample}|${entry.drug}`, {
        state: STATE[entry.verdict] || EMPTY_STATE,
        top: entry.mutations[0],
        count: entry.mutations.length
      });
    }
  }

  const drugs = DRUG_ORDER.filter(d => present.has(d))
    .concat([...present].filter(d => !DRUG_ORDER.includes(d)).sort());

  return { samples, drugs, cells };
}

/** Escape a value for use inside a double-quoted attribute. */
function attr(value) {
  return escapeHtml(String(value ?? ''));
}

/**
 * Render the resistance matrix.
 * Ships a legend and a table view, so nothing depends on colour alone.
 */
export function renderResistanceMatrixHtml(matrix) {
  if (!matrix || matrix.samples.length === 0 || matrix.drugs.length === 0) return '';
  const { samples, drugs, cells } = matrix;

  const head = drugs.map(d =>
    `<th scope="col" title="${attr(drugLabel(d))}"><span>${escapeHtml(d)}</span></th>`
  ).join('');

  const body = samples.map(sample => {
    const tds = drugs.map(drug => {
      const cell = cells.get(`${sample}|${drug}`);
      const state = cell?.state || EMPTY_STATE;
      const detail = cell?.top
        ? `${cell.top.gene || ''} ${cell.top.mutation || ''}`.trim()
        : 'no marker detected';
      const extra = cell && cell.count > 1 ? ` (+${cell.count - 1} more)` : '';
      const tip = `${sample} · ${drugLabel(drug)}: ${state.label}${cell ? ` — ${detail}${extra}` : ''}`;
      return `<td class="rm-cell rm-cell--${state.key}" tabindex="0"
                  aria-label="${attr(tip)}" title="${attr(tip)}">
                <span class="rm-glyph" aria-hidden="true">${state.glyph}</span>
              </td>`;
    }).join('');
    return `<tr><th scope="row" class="rm-sample" title="${attr(sample)}">${escapeHtml(sample)}</th>${tds}</tr>`;
  }).join('');

  const legend = [STATE[GRADE.RESISTANT], STATE[GRADE.UNCERTAIN], STATE[GRADE.NOT_ASSOC], EMPTY_STATE]
    .map(s => `<span class="rm-key rm-key--${s.key}">
                 <span class="rm-glyph" aria-hidden="true">${s.glyph}</span>${escapeHtml(s.label)}
               </span>`).join('');

  // Table view: the WCAG-clean twin, values reachable without colour or hover.
  const tableRows = samples.map(sample => {
    const resistant = drugs.filter(d => cells.get(`${sample}|${d}`)?.state.key === 'resistant');
    const uncertain = drugs.filter(d => cells.get(`${sample}|${d}`)?.state.key === 'uncertain');
    return `<tr>
      <td>${escapeHtml(sample)}</td>
      <td>${resistant.length ? escapeHtml(resistant.join(', ')) : '—'}</td>
      <td>${uncertain.length ? escapeHtml(uncertain.join(', ')) : '—'}</td>
    </tr>`;
  }).join('');

  return `
    <section class="viz-chart viz-resistance-matrix">
      <header class="viz-chart-head">
        <h4>Resistance matrix</h4>
        <p class="viz-chart-sub">${samples.length} samples × ${drugs.length} drugs · WHO catalogue grading</p>
      </header>

      <div class="rm-scroll">
        <table class="rm-grid">
          <thead><tr><th scope="col" class="rm-corner"></th>${head}</tr></thead>
          <tbody>${body}</tbody>
        </table>
      </div>

      <div class="viz-legend">${legend}</div>

      <details class="viz-table-view">
        <summary>Table view</summary>
        <table class="viz-data-table">
          <thead><tr><th>Sample</th><th>Resistant</th><th>Uncertain</th></tr></thead>
          <tbody>${tableRows}</tbody>
        </table>
      </details>
    </section>
  `;
}

/**
 * Render lineage composition as one stacked bar per sample.
 *
 * Identity, so the categorical slots apply. Capped at three branches plus
 * "Other": past that the hues stop being safely distinguishable, and a genuine
 * mixture is two or three strains anyway.
 */
export function renderLineageCompositionHtml(perSample) {
  const entries = (perSample || []).filter(s => s.branches.length > 0);
  if (entries.length === 0) return '';

  const MAX_SLOTS = 3;
  // Colour follows the branch, not its rank within a sample, so a lineage keeps
  // its hue across every bar. Slots go to the branches that actually lead a
  // sample — ranking by global marker count would leave a sample whose lineage
  // is rare in the batch painted entirely as "Other".
  const leads = new Map();
  const totals = new Map();
  entries.forEach(s => {
    const lead = s.branches[0];
    if (lead) leads.set(lead.branch, (leads.get(lead.branch) || 0) + 1);
    s.branches.forEach(b => totals.set(b.branch, (totals.get(b.branch) || 0) + b.specific));
  });
  const ranked = new Map();
  [...leads.entries()]
    .sort((a, b) => b[1] - a[1] || (totals.get(b[0]) || 0) - (totals.get(a[0]) || 0) || a[0].localeCompare(b[0]))
    .forEach(([branch], i) => { if (i < MAX_SLOTS) ranked.set(branch, i + 1); });

  const bars = entries.map(sample => {
    const total = sample.branches.reduce((sum, b) => sum + b.specific, 0) || 1;
    const named = sample.branches.filter(b => ranked.has(b.branch));
    // Fold the whole tail into a single "Other" segment: one gray block per
    // sample, not a shredded row of them.
    const tail = sample.branches
      .filter(b => !ranked.has(b.branch))
      .reduce((sum, b) => sum + b.specific, 0);
    const segments = named.map(b => ({
      slot: ranked.get(b.branch),
      pct: (b.specific / total) * 100,
      label: b.branch,
      specific: b.specific
    })).sort((a, b) => a.slot - b.slot);
    if (tail > 0) {
      segments.push({ slot: 0, pct: (tail / total) * 100, label: 'Other lineages', specific: tail });
    }

    const fills = segments.map(seg => `
      <div class="lc-seg lc-seg--${seg.slot}" style="width:${seg.pct.toFixed(2)}%"
           tabindex="0" title="${attr(`${sample.sample} · ${seg.label}: ${seg.specific} unique markers (${seg.pct.toFixed(0)}%)`)}">
        ${seg.pct >= 18 ? `<span class="lc-seg-label">${escapeHtml(seg.label)}</span>` : ''}
      </div>
    `).join('');

    return `
      <div class="lc-row">
        <span class="lc-name" title="${attr(sample.sample)}">${escapeHtml(sample.sample)}</span>
        <div class="lc-bar">${fills}</div>
      </div>
    `;
  }).join('');

  const legend = [...ranked.entries()]
    .sort((a, b) => a[1] - b[1])
    .map(([branch, slot]) =>
      `<span class="lc-key"><i class="lc-swatch lc-seg--${slot}"></i>${escapeHtml(branch)}</span>`)
    .concat(totals.size > ranked.size
      ? ['<span class="lc-key"><i class="lc-swatch lc-seg--0"></i>Other lineages</span>'] : [])
    .join('');

  const tableRows = entries.map(s => `
    <tr>
      <td>${escapeHtml(s.sample)}</td>
      <td>${escapeHtml(s.branches.map(b => `${b.branch} (${b.specific})`).join(', '))}</td>
    </tr>`).join('');

  return `
    <section class="viz-chart viz-lineage-composition">
      <header class="viz-chart-head">
        <h4>Lineage composition</h4>
        <p class="viz-chart-sub">Share of branch-specific markers per sample</p>
      </header>
      <div class="lc-rows">${bars}</div>
      <div class="viz-legend">${legend}</div>
      <details class="viz-table-view">
        <summary>Table view</summary>
        <table class="viz-data-table">
          <thead><tr><th>Sample</th><th>Branches (unique markers)</th></tr></thead>
          <tbody>${tableRows}</tbody>
        </table>
      </details>
    </section>
  `;
}

/** Per-sample lineage branches, ready for the composition chart. */
export function buildLineageComposition(recordsBySample) {
  const out = [];
  for (const [sample, entries] of recordsBySample) {
    const branches = buildLineageBranches(entries).filter(b => b.specific > 0);
    if (branches.length > 0) out.push({ sample, branches });
  }
  return out.sort((a, b) => a.sample.localeCompare(b.sample));
}
