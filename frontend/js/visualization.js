// ============================================================================
// Visualization Module - Charts and Graphs
// ============================================================================

import { lineageColors, fallbackColors, TIMING } from './config.js';
import { getPanelResultsData, initToolCharts, setToolChart, destroyToolCharts, getToolCharts } from './state.js';
import { logMessage } from './console.js';
import { escapeHtml } from './utils.js';
import { readFastaRange, readTextFile } from './tauri.js';
import {
  isDrPanel,
  buildDrProfile,
  buildLineageLevels,
  buildLineageBranches,
  detectMixedLineages,
  buildAlleleHistogram,
  renderDrProfileHtml,
  renderLineageLevelsHtml,
  renderMixedLineagesHtml,
  renderAlleleHistogramHtml,
  summariseDrProfile,
  drCallLabel
} from './dr-insights.js';
import {
  buildResistanceMatrix,
  renderResistanceMatrixHtml,
  buildLineageComposition,
  renderLineageCompositionHtml
} from './genotype-charts.js';

// Store which tools have active visualizations and their column names
const activeVisualizations = {};
const visualizationColumns = {
  classify: 'major_lineage',
  predict: 'prediction',
  splitfq: 'major_lineage',
  match: 'best_match'
};
const VIZ_PREFERENCES_KEY = 'pathotypr-viz-preferences';
const DEFAULT_VIZ_PREFS = {
  topN: '10',
  metric: 'count', // count | percent
  groupOther: true
};
const TOOL_QUALITY_METRICS = {
  predict: [
    { tokens: ['confidence', 'probability', 'prob', 'score'], label: 'Prediction Confidence', warn: 0.7, danger: 0.5 }
  ],
  classify: [
    { tokens: ['confidence', 'support', 'probability', 'score'], label: 'Classification Confidence', warn: 0.7, danger: 0.5 },
    { tokens: ['depth', 'coverage', 'dp'], label: 'Read Depth', warn: 50, danger: 10 },
    { tokens: ['alt_percent', 'alt%', 'vaf', 'allele_fraction'], label: 'Alt Allele (%)', warn: 80, danger: 50 }
  ],
  splitfq: [
    { tokens: ['alt_fraction', 'alt_percent', 'alt%', 'vaf', 'allele_fraction'], label: 'Alt Fraction (%)', warn: 80, danger: 50 },
    { tokens: ['depth', 'coverage', 'dp'], label: 'Read Depth', warn: 50, danger: 10 },
    { tokens: ['confidence', 'support', 'score'], label: 'Genotype Confidence', warn: 0.7, danger: 0.5 }
  ],
  match: [
    { tokens: ['shared_kmer_fraction', 'identity', 'ani', 'similarity', 'score', 'fraction'], label: 'Match Score', warn: 0.5, danger: 0.2 },
    { tokens: ['coverage', 'cov'], label: 'Coverage', warn: 0.5, danger: 0.2 }
  ]
};
const GENOMIC_TRACK_TOOLS = new Set(['classify']);
const TRACK_STATE = {};
const TRACK_PARSED_DATA_CACHE = new WeakMap();
const TRACK_RENDER_CACHE = {};
const TRACK_FASTA_WINDOW_CACHE = new Map();
const TRACK_MAX_POINTS = 350;
const TRACK_SEQUENCE_FETCH_MAX_SPAN = 1400;
const TRACK_SEQUENCE_SHOW_MAX_SPAN = 900;
const TRACK_SEQUENCE_MAX_RENDER_BASES = 1400;
const TRACK_SEQUENCE_LETTER_MAX_SPAN = 900;
const TRACK_SEQUENCE_FETCH_PADDING_FRACTION = 0.25;
const TRACK_SEQUENCE_FETCH_PADDING_MIN = 40;
const TRACK_SEQUENCE_FETCH_PADDING_MAX = 260;
const TRACK_SEQUENCE_FETCH_DEBOUNCE_MS = 140;
const TRACK_FETCH_TIMERS = {};
const TRACK_GFF_CACHE = new Map();
const TRACK_FASTA_WINDOW_CACHE_MAX_ENTRIES = 160;
const TRACK_GFF_TYPE_PRIORITY = Object.freeze({
  cds: 0,
  gene: 1,
  mrna: 2,
  transcript: 2,
  orf: 3
});
const TRACK_EXPANDED_ANCHORS = {};
const SPLITFQ_TRACK_STATE = {};
const SPLITFQ_TRACK_MAX_POINTS = 900;
const SPLITFQ_TRACK_MIN_SPAN_FACTOR = 0.02;
const SPLITFQ_TRACK_RENDER_THROTTLES = {};
const MATCH_TRACK_STATE = {};
const MATCH_TRACK_MAX_POINTS = 900;
const MATCH_TRACK_MIN_SPAN_FACTOR = 0.08;
const MATCH_TRACK_RENDER_THROTTLES = {};
const TRACK_RECORD_SELECTION_EVENT = 'pathotypr:track-record-selected';
let trackExpandEscapeBound = false;
let visualizationPreferences = loadVisualizationPreferences();

function normalizeTrackSourceRowIndex(rawRowIndex) {
  const value = Number(rawRowIndex);
  return Number.isFinite(value) ? Math.trunc(value) : Number.NaN;
}

function getTrackRenderCache(toolId) {
  if (!TRACK_RENDER_CACHE[toolId]) {
    TRACK_RENDER_CACHE[toolId] = {
      sampleViewCache: new Map(),
      controlsKey: '',
      svgDomKey: '',
      svgInfoCache: new Map()
    };
  }
  return TRACK_RENDER_CACHE[toolId];
}

function clearTrackRenderCache(toolId) {
  delete TRACK_RENDER_CACHE[toolId];
}

function buildTrackFastaWindowCacheKey(fastaPath, start, end, recordName = null) {
  const safePath = String(fastaPath || '').trim();
  const safeRecord = String(recordName || '').trim();
  const from = Math.max(1, Math.floor(Number(start) || 1));
  const to = Math.max(from, Math.ceil(Number(end) || from));
  return `${safePath}::${safeRecord}::${from}-${to}`;
}

function setTrackFastaWindowCacheValue(key, value) {
  TRACK_FASTA_WINDOW_CACHE.set(key, value);
  if (TRACK_FASTA_WINDOW_CACHE.size > TRACK_FASTA_WINDOW_CACHE_MAX_ENTRIES) {
    const oldestKey = TRACK_FASTA_WINDOW_CACHE.keys().next().value;
    TRACK_FASTA_WINDOW_CACHE.delete(oldestKey);
  }
}

async function readFastaRangeCached(fastaPath, start, end, recordName = null) {
  const key = buildTrackFastaWindowCacheKey(fastaPath, start, end, recordName);
  const cached = TRACK_FASTA_WINDOW_CACHE.get(key);
  if (cached?.status === 'ready') {
    // Refresh key order for LRU.
    TRACK_FASTA_WINDOW_CACHE.delete(key);
    setTrackFastaWindowCacheValue(key, cached);
    return cached.data;
  }
  if (cached?.status === 'loading' && cached.promise) {
    return cached.promise;
  }

  const promise = readFastaRange(fastaPath, start, end, recordName);
  setTrackFastaWindowCacheValue(key, { status: 'loading', promise });
  try {
    const response = await promise;
    setTrackFastaWindowCacheValue(key, { status: 'ready', data: response });
    return response;
  } catch (err) {
    TRACK_FASTA_WINDOW_CACHE.delete(key);
    throw err;
  }
}

function emitTrackRecordSelectionEvent(toolId, record) {
  if (!toolId || !record) return;
  const sourceRowIndex = normalizeTrackSourceRowIndex(record.sourceRowIndex);
  if (!Number.isFinite(sourceRowIndex)) return;

  const trackKey = String(record.trackKey || '').trim();
  try {
    document.dispatchEvent(new CustomEvent(TRACK_RECORD_SELECTION_EVENT, {
      detail: {
        toolId,
        sourceRowIndex,
        trackKey
      }
    }));
  } catch {
    // Non-browser environments are intentionally ignored.
  }
}

function loadVisualizationPreferences() {
  try {
    const raw = localStorage.getItem(VIZ_PREFERENCES_KEY);
    if (!raw) return {};
    const parsed = JSON.parse(raw);
    return parsed && typeof parsed === 'object' ? parsed : {};
  } catch {
    return {};
  }
}

function saveVisualizationPreferences() {
  localStorage.setItem(VIZ_PREFERENCES_KEY, JSON.stringify(visualizationPreferences));
}

function getVisualizationPreferences(toolId) {
  return {
    ...DEFAULT_VIZ_PREFS,
    ...(visualizationPreferences[toolId] || {})
  };
}

function setVisualizationPreferences(toolId, patch) {
  visualizationPreferences[toolId] = {
    ...getVisualizationPreferences(toolId),
    ...patch
  };
  saveVisualizationPreferences();
}

function hasChartLibrary() {
  return typeof window.Chart === 'function';
}

function formatPlotDecimal(value, digits = 1) {
  const num = Number(value);
  if (!Number.isFinite(num)) return (0).toFixed(digits);
  return num.toFixed(digits).replace(',', '.');
}

function normalizeHeader(value) {
  return String(value || '')
    .trim()
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '_');
}

function normalizeValue(value) {
  return String(value || '').trim();
}

function isUnknownOutcome(value) {
  const normalized = normalizeValue(value).toLowerCase();
  return (
    normalized === '' ||
    normalized === 'n/a' ||
    normalized === 'na' ||
    normalized === 'unknown' ||
    normalized === 'unresolved' ||
    normalized === 'unassigned' ||
    normalized === 'none' ||
    normalized === '-' ||
    normalized === 'null'
  );
}

function parseNumericValue(raw) {
  const value = String(raw ?? '').trim();
  if (!value) return null;
  const normalized = value.replace(',', '.');
  const parsed = Number.parseFloat(normalized);
  return Number.isFinite(parsed) ? parsed : null;
}

function quantile(values, q) {
  if (!Array.isArray(values) || values.length === 0) return null;
  const sorted = [...values].sort((a, b) => a - b);
  const pos = (sorted.length - 1) * q;
  const base = Math.floor(pos);
  const rest = pos - base;
  if (sorted[base + 1] !== undefined) {
    return sorted[base] + rest * (sorted[base + 1] - sorted[base]);
  }
  return sorted[base];
}

function findHeaderByTokens(headers, tokens) {
  const normalizedTokens = tokens.map(token => normalizeHeader(token));
  return headers.findIndex(header => {
    const normalized = normalizeHeader(header);
    return normalizedTokens.some(token => normalized.includes(token));
  });
}

function findSampleColumnIndex(headers, labelColIdx = -1) {
  const candidateTokens = ['sample', 'genome', 'file', 'isolate', 'strain', 'name'];
  const idx = headers.findIndex((header, colIdx) => {
    if (colIdx === labelColIdx) return false;
    const normalized = normalizeHeader(header);
    return candidateTokens.some(token => normalized.includes(token));
  });
  return idx;
}

function estimateSampleCount(data, labelColIdx = -1) {
  const headers = Array.isArray(data?.headers) ? data.headers : [];
  const rows = Array.isArray(data?.rows) ? data.rows : [];
  if (rows.length === 0) return 0;

  const sampleColIdx = findSampleColumnIndex(headers, labelColIdx);
  if (sampleColIdx === -1) return rows.length;

  const unique = new Set();
  rows.forEach(row => {
    const value = normalizeValue(row[sampleColIdx]);
    if (value) unique.add(value);
  });
  return unique.size || rows.length;
}

function findQualityMetricSummary(toolId, data) {
  const headers = Array.isArray(data?.headers) ? data.headers : [];
  const rows = Array.isArray(data?.rows) ? data.rows : [];
  if (headers.length === 0 || rows.length === 0) return null;

  const candidates = TOOL_QUALITY_METRICS[toolId] || [];
  for (const candidate of candidates) {
    const columnIdx = findHeaderByTokens(headers, candidate.tokens);
    if (columnIdx === -1) continue;

    const values = rows
      .map(row => parseNumericValue(row[columnIdx]))
      .filter(v => v !== null);

    if (values.length < Math.max(3, Math.ceil(rows.length * 0.2))) continue;

    const median = quantile(values, 0.5);
    const p90 = quantile(values, 0.9);
    const min = Math.min(...values);
    const max = Math.max(...values);
    return {
      label: candidate.label,
      median,
      p90,
      min,
      max,
      count: values.length,
      warn: candidate.warn,
      danger: candidate.danger
    };
  }

  return null;
}

/**
 * Setup visualization for all tools
 */
export function setupAllVisualizations() {
  if (!hasChartLibrary()) {
    logMessage('Chart.js not loaded. Visualization is disabled.', 'warning');
  }

  setupToolVisualization('classify', 'Lineage Distribution', 'major_lineage');
  setupToolVisualization('predict', 'Prediction Distribution', 'prediction');
  setupToolVisualization('splitfq', 'Genotype Distribution', 'major_lineage');
  setupToolVisualization('match', 'Reference Distribution', 'best_match');
}

function ensureVisualizationControls(vizPanel, toolId) {
  if (!vizPanel || vizPanel.querySelector('.visualization-controls')) return;
  const prefs = getVisualizationPreferences(toolId);
  const controls = document.createElement('div');
  controls.className = 'visualization-controls';
  controls.innerHTML = `
    <label>
      Metric
      <select class="viz-metric-select" data-tool="${toolId}">
        <option value="count" ${prefs.metric === 'count' ? 'selected' : ''}>Count</option>
        <option value="percent" ${prefs.metric === 'percent' ? 'selected' : ''}>Percent</option>
      </select>
    </label>
    <label>
      Top
      <select class="viz-topn-select" data-tool="${toolId}">
        <option value="5" ${prefs.topN === '5' ? 'selected' : ''}>Top 5</option>
        <option value="10" ${prefs.topN === '10' ? 'selected' : ''}>Top 10</option>
        <option value="20" ${prefs.topN === '20' ? 'selected' : ''}>Top 20</option>
        <option value="all" ${prefs.topN === 'all' ? 'selected' : ''}>All</option>
      </select>
    </label>
    <label class="viz-group-other">
      <input type="checkbox" class="viz-other-checkbox" data-tool="${toolId}" ${prefs.groupOther ? 'checked' : ''}>
      Group remainder as "Other"
    </label>
  `;
  const content = vizPanel.querySelector('.visualization-content');
  if (content) content.prepend(controls);
}

function ensureVisualizationInsights(vizPanel, toolId) {
  if (!vizPanel || vizPanel.querySelector('.viz-insights')) return;

  const insights = document.createElement('section');
  insights.className = 'viz-insights';
  insights.innerHTML = `
    <div class="viz-insights-grid" id="${toolId}-viz-insights-grid"></div>
    <div class="viz-outcome-breakdown" id="${toolId}-viz-outcome-breakdown"></div>
    <section class="viz-genomic-track hidden" id="${toolId}-viz-genomic-track">
      <div class="viz-track-header">
        <h4>Genomic Track</h4>
        <div class="viz-track-controls" id="${toolId}-viz-track-controls"></div>
      </div>
      <div class="viz-track-body" id="${toolId}-viz-track-body"></div>
    </section>
  `;

  const content = vizPanel.querySelector('.visualization-content');
  if (content) content.appendChild(insights);
}

function ensureSingleSampleSummary(vizPanel, toolId) {
  if (!vizPanel) return;
  const content = vizPanel.querySelector('.visualization-content');
  const chartContainer = content?.querySelector('.chart-container');
  if (!content || !chartContainer || vizPanel.querySelector('.viz-single-sample-summary')) return;

  const summary = document.createElement('section');
  summary.className = 'viz-single-sample-summary hidden';
  summary.id = `${toolId}-viz-single-sample-summary`;
  content.insertBefore(summary, chartContainer);
}

function renderSingleSampleSummary(toolId, summary) {
  const vizPanel = document.getElementById(`${toolId}-visualization`);
  const summaryEl = document.getElementById(`${toolId}-viz-single-sample-summary`);
  const chartContainer = vizPanel?.querySelector('.chart-container');
  const legend = document.getElementById(`${toolId}-chart-legend`);
  if (!vizPanel || !summaryEl || !chartContainer || !legend) return;

  if (!summary) {
    vizPanel.classList.remove('single-sample-mode');
    summaryEl.classList.add('hidden');
    summaryEl.replaceChildren();
    chartContainer.classList.remove('hidden');
    legend.classList.remove('hidden');
    return;
  }

  vizPanel.classList.add('single-sample-mode');
  summaryEl.classList.remove('hidden');
  chartContainer.classList.add('hidden');
  legend.classList.add('hidden');

  summaryEl.innerHTML = `
    <div class="viz-single-main">
      <span class="viz-single-label">Single Sample Mode</span>
      <strong class="viz-single-sample-name">${escapeHtml(summary.sampleName)}</strong>
      <p class="viz-single-call">
        Primary call: <span>${escapeHtml(summary.outcomeLabel)}</span>
      </p>
    </div>
    <div class="viz-single-metrics">
      <article class="viz-single-metric">
        <span>Assigned</span>
        <strong>${escapeHtml(String(summary.assignedCount))}</strong>
      </article>
      <article class="viz-single-metric">
        <span>Rows</span>
        <strong>${escapeHtml(String(summary.rowCount))}</strong>
      </article>
      <article class="viz-single-metric">
        <span>${escapeHtml(summary.qualityLabel)}</span>
        <strong>${escapeHtml(summary.qualityValue)}</strong>
      </article>
      <article class="viz-single-metric">
        <span>Unresolved</span>
        <strong>${escapeHtml(String(summary.unresolvedCount))}</strong>
      </article>
    </div>
  `;
}

function clearChartEmptyState(toolId) {
  ['donut', 'bar'].forEach(chartType => {
    const canvas = document.getElementById(`${toolId}-${chartType}-chart`);
    const wrapper = canvas?.closest('.chart-wrapper');
    if (!wrapper) return;
    wrapper.classList.remove('chart-empty');
    wrapper.querySelector('.chart-empty-message')?.remove();
  });
}

function setChartEmptyState(toolId, message) {
  ['donut', 'bar'].forEach(chartType => {
    const canvas = document.getElementById(`${toolId}-${chartType}-chart`);
    const wrapper = canvas?.closest('.chart-wrapper');
    if (!wrapper) return;
    wrapper.classList.add('chart-empty');
    let empty = wrapper.querySelector('.chart-empty-message');
    if (!empty) {
      empty = document.createElement('div');
      empty.className = 'chart-empty-message';
      wrapper.appendChild(empty);
    }
    empty.textContent = message;
  });
}

function formatIntegerForViz(value) {
  const num = Number(value);
  if (!Number.isFinite(num)) return '0';
  return Math.round(num).toLocaleString('en-US');
}

function resetChartCanvases(toolId) {
  ['donut', 'bar'].forEach(chartType => {
    const existingCanvas = document.getElementById(`${toolId}-${chartType}-chart`);
    if (!existingCanvas || !existingCanvas.parentElement) return;

    const replacement = document.createElement('canvas');
    replacement.id = existingCanvas.id;
    const ariaLabel = existingCanvas.getAttribute('aria-label');
    if (ariaLabel) {
      replacement.setAttribute('aria-label', ariaLabel);
    }
    existingCanvas.parentElement.replaceChild(replacement, existingCanvas);
  });
}

function parseSplitFastqLineageCounts(rawValue) {
  const text = String(rawValue || '').trim();
  if (!text) return [];

  const counts = new Map();
  // Entries are separated by whitespace (both backends join with " ").
  // ';' must NOT be treated as a separator: it is the delimiter *inside* a
  // nested lineage path (e.g. "L2;L2.2:5"), so splitting on it shredded the
  // path and mislabelled the lineage.
  text.split(/[\s,]+/).forEach(entry => {
    const part = String(entry || '').trim();
    if (!part) return;

    const sep = part.lastIndexOf(':');
    if (sep <= 0 || sep >= part.length - 1) return;
    const lineage = part.slice(0, sep).trim();
    const count = parseNumericValue(part.slice(sep + 1));
    if (!lineage || count === null) return;

    counts.set(lineage, (counts.get(lineage) || 0) + count);
  });

  return [...counts.entries()]
    .map(([lineage, count]) => ({ lineage, count }))
    .sort((a, b) => b.count - a.count || a.lineage.localeCompare(b.lineage, undefined, { sensitivity: 'base' }));
}

function parseSplitFastqPosition(rawValue, fallbackValue) {
  const text = normalizeValue(rawValue);
  if (!text) return { value: fallbackValue, inferred: true };

  const direct = parseNumericValue(text);
  if (direct !== null) return { value: direct, inferred: false };

  const matched = text.match(/-?\d+(?:\.\d+)?/);
  if (matched) {
    const parsed = Number.parseFloat(matched[0]);
    if (Number.isFinite(parsed)) return { value: parsed, inferred: false };
  }

  return { value: fallbackValue, inferred: true };
}

function getSplitFastqSupportColor(supportPct) {
  if (!Number.isFinite(supportPct)) return '#9ca3af';
  if (supportPct >= 90) return '#16a34a';
  if (supportPct >= 70) return '#f59e0b';
  return '#dc2626';
}

function splitfqTrackFormatAxisValue(value, useBpAxis) {
  return useBpAxis ? formatTrackBp(value) : formatIntegerForViz(value);
}

function splitfqTrackBuildModelKey(model) {
  const firstId = model.trackItems?.[0]?.id || '';
  const lastId = model.trackItems?.[model.trackItems.length - 1]?.id || '';
  return [
    model.mode || 'summary',
    model.trackItems?.length || 0,
    Math.round(model.trackDomainStart || 0),
    Math.round(model.trackDomainEnd || 0),
    firstId,
    lastId
  ].join('|');
}

function splitfqTrackGetDomain(items) {
  if (!Array.isArray(items) || items.length === 0) {
    return { start: 0, end: 1 };
  }
  let min = Number.POSITIVE_INFINITY;
  let max = Number.NEGATIVE_INFINITY;
  items.forEach(item => {
    const x = Number(item?.x);
    if (!Number.isFinite(x)) return;
    if (x < min) min = x;
    if (x > max) max = x;
  });
  if (!Number.isFinite(min) || !Number.isFinite(max)) {
    return { start: 0, end: 1 };
  }
  if (max <= min) max = min + 1;
  return { start: min, end: max };
}

function splitfqTrackDefaultSpan(mode, domainStart, domainEnd, itemCount) {
  const fullSpan = Math.max(1, domainEnd - domainStart);
  if (!Number.isFinite(fullSpan) || fullSpan <= 0) return 1;
  if (mode === 'detailed') {
    return Math.min(fullSpan, Math.max(60, fullSpan * 0.3));
  }
  const samplesWindow = Math.max(8, Math.min(itemCount, 40));
  return Math.min(fullSpan, Math.max(4, samplesWindow));
}

function splitfqTrackSetWindow(state, start, end, domainStart, domainEnd) {
  const minBound = Number.isFinite(domainStart) ? domainStart : 0;
  let maxBound = Number.isFinite(domainEnd) ? domainEnd : minBound + 1;
  if (maxBound <= minBound) maxBound = minBound + 1;

  const fullSpan = Math.max(1, maxBound - minBound);
  const minSpan = Math.min(fullSpan, Math.max(1, fullSpan * SPLITFQ_TRACK_MIN_SPAN_FACTOR));

  let windowStart = Number.isFinite(start) ? start : minBound;
  let windowEnd = Number.isFinite(end) ? end : maxBound;
  if (windowEnd < windowStart) {
    const tmp = windowStart;
    windowStart = windowEnd;
    windowEnd = tmp;
  }

  if ((windowEnd - windowStart) < minSpan) {
    const center = (windowStart + windowEnd) / 2;
    windowStart = center - (minSpan / 2);
    windowEnd = center + (minSpan / 2);
  }

  if (windowStart < minBound) {
    windowEnd += (minBound - windowStart);
    windowStart = minBound;
  }
  if (windowEnd > maxBound) {
    windowStart -= (windowEnd - maxBound);
    windowEnd = maxBound;
  }

  if (windowStart < minBound) windowStart = minBound;
  if ((windowEnd - windowStart) < minSpan) {
    windowEnd = Math.min(maxBound, windowStart + minSpan);
    windowStart = Math.max(minBound, windowEnd - minSpan);
  }

  state.windowStart = windowStart;
  state.windowEnd = windowEnd;
}

function splitfqTrackZoom(state, factor, center, domainStart, domainEnd) {
  const currentSpan = Math.max(1, state.windowEnd - state.windowStart);
  const targetSpan = currentSpan * factor;
  const xCenter = Number.isFinite(center)
    ? center
    : (state.windowStart + state.windowEnd) / 2;
  splitfqTrackSetWindow(
    state,
    xCenter - targetSpan / 2,
    xCenter + targetSpan / 2,
    domainStart,
    domainEnd
  );
}

function splitfqTrackPan(state, fraction, domainStart, domainEnd) {
  const span = Math.max(1, state.windowEnd - state.windowStart);
  const shift = span * fraction;
  splitfqTrackSetWindow(
    state,
    state.windowStart + shift,
    state.windowEnd + shift,
    domainStart,
    domainEnd
  );
}

function splitfqTrackDownsampleItems(items, maxPoints, selectedItemId = '') {
  if (!Array.isArray(items) || items.length <= maxPoints) return items || [];
  const stride = Math.ceil(items.length / maxPoints);
  const sampled = items.filter((_, index) => index % stride === 0);
  if (selectedItemId && !sampled.some(item => item.id === selectedItemId)) {
    const selected = items.find(item => item.id === selectedItemId);
    if (selected) sampled.push(selected);
  }
  return sampled.sort((a, b) => a.x - b.x || a.id.localeCompare(b.id));
}

function splitfqTrackItemSearchText(item) {
  return [
    item.sample,
    item.label,
    item.lineage,
    item.mode === 'detailed' ? `pos ${formatTrackBp(item.x)}` : `sample ${formatIntegerForViz(item.x)}`,
    Number.isFinite(item.supportPct) ? `${formatPlotDecimal(item.supportPct, 1)}%` : ''
  ].filter(Boolean).join(' ').toLowerCase();
}

function splitfqTrackFindMatches(items, query) {
  const normalized = String(query || '').trim().toLowerCase();
  if (!normalized) return [];
  return items.filter(item => splitfqTrackItemSearchText(item).includes(normalized));
}

function splitfqTrackBuildSearchOptions(items, limit = 180) {
  if (!Array.isArray(items) || items.length === 0) return '';
  const unique = new Set();
  const options = [];
  for (const item of items) {
    const label = item.mode === 'detailed'
      ? `${item.sample} · ${item.label}`
      : `${item.sample} · ${item.lineage || 'Unresolved'}`;
    const value = String(label || '').trim();
    if (!value || unique.has(value)) continue;
    unique.add(value);
    options.push(value);
    if (options.length >= limit) break;
  }
  return options.map(value => `<option value="${escapeHtml(value)}"></option>`).join('');
}

function splitfqTrackSearchStatus(state) {
  const query = String(state.searchQuery || '').trim();
  const matches = Array.isArray(state.searchMatchIds) ? state.searchMatchIds.length : 0;
  if (!query) return 'Type to search by sample, marker, lineage, or position.';
  if (matches === 0) return 'No matching points.';
  const index = Math.max(0, Number(state.searchMatchIndex || 0)) + 1;
  return `${index}/${matches} match${matches === 1 ? '' : 'es'}`;
}

function splitfqTrackFindItemById(items, itemId) {
  if (!itemId) return null;
  return items.find(item => item.id === itemId) || null;
}

function splitfqTrackFocusItem(state, item, domainStart, domainEnd, zoom = true) {
  if (!item) return;
  state.selectedItemId = item.id;
  if (!zoom) {
    if (item.x < state.windowStart || item.x > state.windowEnd) {
      const span = Math.max(1, state.windowEnd - state.windowStart);
      splitfqTrackSetWindow(
        state,
        item.x - span / 2,
        item.x + span / 2,
        domainStart,
        domainEnd
      );
    }
    return;
  }
  const span = Math.max(1, state.windowEnd - state.windowStart);
  const targetSpan = span * 0.45;
  splitfqTrackSetWindow(
    state,
    item.x - targetSpan / 2,
    item.x + targetSpan / 2,
    domainStart,
    domainEnd
  );
}

function splitfqTrackDomainValueFromEvent(event, svgEl, state, useOverview = false) {
  if (!svgEl) return (state.windowStart + state.windowEnd) / 2;

  const viewWidth = Number(svgEl.dataset.viewWidth || 980);
  const left = useOverview
    ? Number(svgEl.dataset.overviewLeft || 0)
    : Number(svgEl.dataset.plotLeft || 0);
  const width = useOverview
    ? Number(svgEl.dataset.overviewWidth || 1)
    : Number(svgEl.dataset.plotWidth || 1);
  const domainStart = useOverview ? state.domainStart : state.windowStart;
  const domainEnd = useOverview ? state.domainEnd : state.windowEnd;
  const domainSpan = Math.max(1, domainEnd - domainStart);

  const xView = getSvgViewXFromClient(svgEl, event.clientX, viewWidth / 2);
  const relPlot = Math.max(0, Math.min(1, (xView - left) / Math.max(1, width)));
  return domainStart + relPlot * domainSpan;
}

function getSvgViewXFromClient(svgEl, clientX, fallback = 0) {
  if (!svgEl || !Number.isFinite(clientX)) return fallback;

  try {
    if (typeof svgEl.createSVGPoint === 'function' && typeof svgEl.getScreenCTM === 'function') {
      const point = svgEl.createSVGPoint();
      point.x = clientX;
      point.y = 0;
      const ctm = svgEl.getScreenCTM();
      if (ctm && typeof ctm.inverse === 'function') {
        const transformed = point.matrixTransform(ctm.inverse());
        if (Number.isFinite(transformed.x)) return transformed.x;
      }
    }
  } catch {
    // Keep fallback conversion below.
  }

  const rect = svgEl.getBoundingClientRect();
  const viewWidth = Number(svgEl.dataset.viewWidth || 980);
  if (rect.width <= 0) return fallback;
  const rel = Math.max(0, Math.min(1, (clientX - rect.left) / rect.width));
  return rel * viewWidth;
}

function buildSplitfqTrackSvg(model, state, items) {
  const readLaneEnabled = model.mode === 'detailed';
  const width = 980;
  const height = readLaneEnabled ? 372 : 316;
  const marginLeft = 58;
  const marginRight = 18;
  const marginTop = 20;
  const marginBottom = readLaneEnabled ? 136 : 86;
  const plotWidth = width - marginLeft - marginRight;
  const plotHeight = height - marginTop - marginBottom;
  const plotBottom = marginTop + plotHeight;
  const readLaneTop = plotBottom + 12;
  const readLaneHeight = readLaneEnabled ? 40 : 0;
  const readLaneBottom = readLaneTop + readLaneHeight;
  const overviewY = readLaneEnabled ? (readLaneBottom + 24) : (plotBottom + 26);
  const overviewHeight = 8;

  const windowSpan = Math.max(1, state.windowEnd - state.windowStart);
  const fullSpan = Math.max(1, state.domainEnd - state.domainStart);
  const xForWindow = x => marginLeft + ((x - state.windowStart) / windowSpan) * plotWidth;
  const yForSupport = support => marginTop + (1 - (support / 100)) * plotHeight;
  const xForOverview = x => marginLeft + ((x - state.domainStart) / fullSpan) * plotWidth;

  const visibleItems = items
    .filter(item => item.x >= state.windowStart && item.x <= state.windowEnd)
    .sort((a, b) => a.x - b.x || a.id.localeCompare(b.id));
  const renderItems = splitfqTrackDownsampleItems(visibleItems, SPLITFQ_TRACK_MAX_POINTS, state.selectedItemId);
  const supports = visibleItems
    .map(item => Number(item.supportPct))
    .filter(value => Number.isFinite(value))
    .map(value => Math.max(0, Math.min(100, value)));
  const medianSupport = supports.length > 0 ? quantile(supports, 0.5) : null;
  const selectedItem = splitfqTrackFindItemById(items, state.selectedItemId);
  const matchIds = new Set(Array.isArray(state.searchMatchIds) ? state.searchMatchIds : []);

  const yTicks = [0, 25, 50, 75, 100];
  const yTickSvg = yTicks.map(value => {
    const y = yForSupport(value);
    return `
      <line x1="${marginLeft}" y1="${y}" x2="${marginLeft + plotWidth}" y2="${y}" class="viz-splitfq-track-grid"/>
      <text x="${marginLeft - 8}" y="${y + 4}" text-anchor="end" class="viz-splitfq-track-axis-label">${value}%</text>
    `;
  }).join('');

  const xTickCount = 5;
  const xTickLabelY = readLaneEnabled ? (readLaneBottom + 16) : (plotBottom + 19);
  const xTickSvg = Array.from({ length: xTickCount + 1 }, (_, index) => {
    const value = state.windowStart + (windowSpan * index / xTickCount);
    const x = xForWindow(value);
    const label = splitfqTrackFormatAxisValue(Math.round(value), model.trackUseBpAxis);
    return `
      <line x1="${x}" y1="${plotBottom}" x2="${x}" y2="${plotBottom + 6}" class="viz-splitfq-track-axis-tick"/>
      <text x="${x}" y="${xTickLabelY}" text-anchor="middle" class="viz-splitfq-track-axis-label">${escapeSvg(label)}</text>
    `;
  }).join('');

  const pointSvg = renderItems.map(item => {
    const support = Math.max(0, Math.min(100, Number(item.supportPct) || 0));
    const x = xForWindow(item.x);
    const y = yForSupport(support);
    const color = getSplitFastqSupportColor(support);
    const selected = item.id === state.selectedItemId;
    const matched = matchIds.has(item.id);
    const title = model.mode === 'detailed'
      ? `${item.sample} · ${item.label}\nALT ${formatPlotDecimal(support, 1)}% · ${formatIntegerForViz(item.altCount || 0)} ALT / ${formatIntegerForViz(item.coverage || 0)} total`
      : `${item.sample} · ${item.lineage || 'Unresolved'}\nSupport ${formatPlotDecimal(support, 1)}%`;
    return `
      <line x1="${x}" y1="${plotBottom}" x2="${x}" y2="${y}" class="viz-splitfq-point-stem${selected ? ' is-selected' : ''}" stroke="${color}"/>
      <circle cx="${x}" cy="${y}" r="${selected ? 4.5 : 3.2}" class="viz-splitfq-point${selected ? ' is-selected' : ''}${matched ? ' is-match' : ''}" fill="${color}"/>
      <circle
        cx="${x}"
        cy="${y}"
        r="8"
        class="viz-splitfq-point-hit"
        data-item-id="${escapeHtml(item.id)}"
      >
        <title>${escapeSvg(title)}</title>
      </circle>
    `;
  }).join('');

  const medianSvg = Number.isFinite(medianSupport)
    ? `
      <line
        x1="${marginLeft}"
        y1="${yForSupport(medianSupport)}"
        x2="${marginLeft + plotWidth}"
        y2="${yForSupport(medianSupport)}"
        class="viz-splitfq-track-median"
      />
      <text x="${marginLeft + plotWidth - 4}" y="${yForSupport(medianSupport) - 6}" text-anchor="end" class="viz-splitfq-track-median-label">
        median ${formatPlotDecimal(medianSupport, 1)}%
      </text>
    `
    : '';

  const readSupportLaneSvg = readLaneEnabled
    ? (() => {
      const bars = splitfqTrackDownsampleItems(renderItems, 560, state.selectedItemId);
      if (bars.length === 0) return '';
      const maxCoverage = Math.max(1, ...bars.map(item => Math.max(0, Number(item.coverage) || 0)));
      const widthPerPoint = plotWidth / Math.max(10, bars.length || 1);
      const barWidth = Math.max(1.4, Math.min(4.8, widthPerPoint * 0.55));
      const midCoverage = maxCoverage / 2;

      const panelGuides = `
        <rect x="${marginLeft}" y="${readLaneTop}" width="${plotWidth}" height="${readLaneHeight}" class="viz-splitfq-track-read-panel-bg"/>
        <line x1="${marginLeft}" y1="${readLaneBottom}" x2="${marginLeft + plotWidth}" y2="${readLaneBottom}" class="viz-splitfq-track-read-axis"/>
        <line x1="${marginLeft}" y1="${readLaneTop}" x2="${marginLeft + plotWidth}" y2="${readLaneTop}" class="viz-splitfq-track-read-grid"/>
        <line x1="${marginLeft}" y1="${readLaneTop + (readLaneHeight / 2)}" x2="${marginLeft + plotWidth}" y2="${readLaneTop + (readLaneHeight / 2)}" class="viz-splitfq-track-read-grid"/>
        <text x="${marginLeft - 8}" y="${readLaneTop + 4}" text-anchor="end" class="viz-splitfq-track-axis-label">${formatIntegerForViz(maxCoverage)}</text>
        <text x="${marginLeft - 8}" y="${readLaneTop + (readLaneHeight / 2) + 4}" text-anchor="end" class="viz-splitfq-track-axis-label">${formatIntegerForViz(midCoverage)}</text>
        <text x="${marginLeft - 8}" y="${readLaneBottom + 4}" text-anchor="end" class="viz-splitfq-track-axis-label">0</text>
        <text x="${marginLeft}" y="${readLaneTop - 4}" class="viz-splitfq-track-axis-label">SNP-support reads (stacked ALT/REF)</text>
      `;

      const laneBars = bars.map(item => {
        const coverage = Math.max(0, Number(item.coverage) || 0);
        if (coverage <= 0) return '';
        const altCountRaw = Math.max(0, Number(item.altCount) || 0);
        const altCount = Math.min(coverage, altCountRaw);
        const x = xForWindow(item.x) - (barWidth / 2);
        const scaledHeight = Math.max(1, (coverage / maxCoverage) * readLaneHeight);
        const y = readLaneBottom - scaledHeight;
        const altHeight = scaledHeight * (altCount / coverage);
        const refHeight = Math.max(0, scaledHeight - altHeight);
        const selectedClass = item.id === state.selectedItemId ? ' is-selected' : '';
        const hitWidth = Math.max(6, barWidth + 2);
        const hitX = x + (barWidth / 2) - (hitWidth / 2);
        const title = `${item.sample} · ${item.label}\nALT ${formatIntegerForViz(altCount)} / REF ${formatIntegerForViz(coverage - altCount)} · total ${formatIntegerForViz(coverage)}`;

        return `
          <rect x="${x}" y="${y}" width="${barWidth}" height="${scaledHeight}" rx="0.7" class="viz-splitfq-track-read-bg${selectedClass}"/>
          <rect x="${x}" y="${y}" width="${barWidth}" height="${refHeight}" rx="0.7" class="viz-splitfq-track-read-ref${selectedClass}"/>
          <rect x="${x}" y="${y + refHeight}" width="${barWidth}" height="${Math.max(0.9, altHeight)}" rx="0.7" class="viz-splitfq-track-read-alt${selectedClass}"/>
          <rect
            x="${hitX}"
            y="${readLaneTop - 3}"
            width="${hitWidth}"
            height="${readLaneHeight + 6}"
            class="viz-splitfq-point-hit viz-splitfq-read-hit"
            data-item-id="${escapeHtml(item.id)}"
          >
            <title>${escapeSvg(title)}</title>
          </rect>
        `;
      }).join('');

      return `
        ${panelGuides}
        ${laneBars}
      `;
    })()
    : '';

  const overviewPoints = splitfqTrackDownsampleItems(items, 650, state.selectedItemId)
    .map(item => {
      const x = xForOverview(item.x);
      const support = Math.max(0, Math.min(100, Number(item.supportPct) || 0));
      const color = getSplitFastqSupportColor(support);
      return `<rect x="${x}" y="${overviewY}" width="1.5" height="${overviewHeight}" fill="${color}" opacity="0.65"/>`;
    })
    .join('');

  const overviewWindowX = xForOverview(state.windowStart);
  const overviewWindowW = Math.max(2, (windowSpan / fullSpan) * plotWidth);
  const overviewSvg = `
    <rect x="${marginLeft}" y="${overviewY}" width="${plotWidth}" height="${overviewHeight}" class="viz-splitfq-overview-bg"/>
    ${overviewPoints}
    <rect
      x="${overviewWindowX}"
      y="${overviewY - 1}"
      width="${overviewWindowW}"
      height="${overviewHeight + 2}"
      class="viz-splitfq-overview-window"
    />
    <rect
      x="${marginLeft}"
      y="${overviewY - 5}"
      width="${plotWidth}"
      height="${overviewHeight + 10}"
      class="viz-splitfq-overview-hit"
      data-splitfq-role="overview-hit"
    />
    <text x="${marginLeft}" y="${overviewY + 20}" class="viz-splitfq-track-axis-label">${escapeSvg(splitfqTrackFormatAxisValue(state.domainStart, model.trackUseBpAxis))}</text>
    <text x="${marginLeft + plotWidth}" y="${overviewY + 20}" text-anchor="end" class="viz-splitfq-track-axis-label">${escapeSvg(splitfqTrackFormatAxisValue(state.domainEnd, model.trackUseBpAxis))}</text>
  `;

  const svg = `
    <svg
      class="viz-splitfq-track-svg"
      viewBox="0 0 ${width} ${height}"
      preserveAspectRatio="xMidYMid meet"
      style="aspect-ratio: ${width} / ${height};"
      role="img"
      aria-label="Interactive split fastq evidence track"
      data-view-width="${width}"
      data-view-height="${height}"
      data-plot-left="${marginLeft}"
      data-plot-width="${plotWidth}"
      data-overview-left="${marginLeft}"
      data-overview-width="${plotWidth}"
    >
      <rect x="${marginLeft}" y="${marginTop}" width="${plotWidth}" height="${plotHeight}" class="viz-splitfq-track-plot-bg"/>
      ${yTickSvg}
      ${medianSvg}
      <line x1="${marginLeft}" y1="${plotBottom}" x2="${marginLeft + plotWidth}" y2="${plotBottom}" class="viz-splitfq-track-axis"/>
      ${readSupportLaneSvg}
      ${xTickSvg}
      ${pointSvg}
      <text x="${marginLeft + plotWidth / 2}" y="${height - 6}" text-anchor="middle" class="viz-splitfq-track-axis-title">
        ${escapeSvg(model.trackAxisLabel)}
      </text>
      <text
        x="${14}"
        y="${marginTop + plotHeight / 2}"
        text-anchor="middle"
        class="viz-splitfq-track-axis-title"
        transform="rotate(-90 14 ${marginTop + plotHeight / 2})"
      >
        ALT support (%)
      </text>
      ${overviewSvg}
    </svg>
  `;

  return {
    svg,
    visibleCount: visibleItems.length,
    renderedCount: renderItems.length,
    selectedItem,
    readLaneEnabled
  };
}

function buildSplitfqTrackDetailHtml(model, selectedItem, visibleCount) {
  if (!selectedItem) {
    return `
      <div class="viz-splitfq-track-detail-empty">
        Click a point to focus it. ${visibleCount > 0 ? 'Use drag, wheel zoom, and search to navigate.' : 'No points in current window.'}
      </div>
    `;
  }

  if (model.mode === 'detailed') {
    const coverage = Math.max(0, Number(selectedItem.coverage) || 0);
    const altCount = Math.min(coverage, Math.max(0, Number(selectedItem.altCount) || 0));
    const refCount = Math.max(0, Number(selectedItem.refCount) || Math.max(0, coverage - altCount));
    const altPct = coverage > 0 ? (altCount / coverage) * 100 : 0;
    const refPct = coverage > 0 ? Math.max(0, 100 - altPct) : 0;
    return `
      <div class="viz-splitfq-track-detail-row">
        <span class="viz-splitfq-track-detail-label">Sample</span>
        <strong>${escapeHtml(selectedItem.sample || 'n/a')}</strong>
      </div>
      <div class="viz-splitfq-track-detail-row">
        <span class="viz-splitfq-track-detail-label">Marker</span>
        <strong>${escapeHtml(selectedItem.label || 'n/a')}</strong>
      </div>
      <div class="viz-splitfq-track-detail-row">
        <span class="viz-splitfq-track-detail-label">ALT Support</span>
        <strong>${formatPlotDecimal(selectedItem.supportPct || 0, 1)}%</strong>
      </div>
      <div class="viz-splitfq-track-detail-row">
        <span class="viz-splitfq-track-detail-label">Reads</span>
        <strong>${formatIntegerForViz(altCount)} ALT / ${formatIntegerForViz(coverage)} total</strong>
      </div>
      <div class="viz-splitfq-track-readbar-wrap">
        <span class="viz-splitfq-track-detail-label">SNP read support</span>
        <div class="viz-splitfq-track-readbar" role="img" aria-label="Read support bar for selected SNP">
          <span class="viz-splitfq-track-readbar-ref" style="width:${formatPlotDecimal(refPct, 2)}%"></span>
          <span class="viz-splitfq-track-readbar-alt" style="width:${formatPlotDecimal(altPct, 2)}%"></span>
        </div>
        <span class="viz-splitfq-track-readbar-values">
          REF ${formatIntegerForViz(refCount)} (${formatPlotDecimal(refPct, 1)}%) · ALT ${formatIntegerForViz(altCount)} (${formatPlotDecimal(altPct, 1)}%)
        </span>
      </div>
    `;
  }

  return `
    <div class="viz-splitfq-track-detail-row">
      <span class="viz-splitfq-track-detail-label">Sample</span>
      <strong>${escapeHtml(selectedItem.sample || 'n/a')}</strong>
    </div>
    <div class="viz-splitfq-track-detail-row">
      <span class="viz-splitfq-track-detail-label">Major lineage</span>
      <strong>${escapeHtml(selectedItem.lineage || 'Unresolved')}</strong>
    </div>
    <div class="viz-splitfq-track-detail-row">
      <span class="viz-splitfq-track-detail-label">Major support</span>
      <strong>${formatPlotDecimal(selectedItem.supportPct || 0, 1)}%</strong>
    </div>
    <div class="viz-splitfq-track-detail-row">
      <span class="viz-splitfq-track-detail-label">Evidence</span>
      <strong>${formatIntegerForViz(selectedItem.majorCount || 0)} / ${formatIntegerForViz(selectedItem.totalEvidence || 0)} markers</strong>
    </div>
  `;
}

function buildSplitFastqDetailedModel(headers, rows) {
  const sampleIdx = findSampleColumnIndex(headers, -1);
  const refCountIdx = findHeaderIndexByTokens(headers, ['ref_count', 'reference_count', 'ref_reads']);
  const altCountIdx = findHeaderIndexByTokens(headers, ['alt_count', 'alternate_count', 'alt_reads', 'mut_count']);
  const altFractionIdx = findHeaderIndexByTokens(headers, ['alt_fraction', 'alt_percent', 'alt_pct', 'vaf', 'allele_fraction']);
  const markerIdx = findHeaderIndexByTokens(headers, ['pos', 'position', 'marker', 'marker_id', 'locus', 'site']);
  const lineageIdx = findHeaderIndexByTokens(headers, ['lineage_path', 'lineage', 'major_lineage']);

  const hasReadEvidence = refCountIdx !== -1 || altCountIdx !== -1 || altFractionIdx !== -1;
  if (!hasReadEvidence) return null;

  const sampleStats = new Map();
  const trackItems = [];
  const altFractions = [];
  const markerRows = [];
  let lociWithReadSupport = 0;
  let totalCoverage = 0;
  let totalAltReads = 0;
  let numericPositionCount = 0;

  rows.forEach((row, rowIdx) => {
    const sampleName = sampleIdx !== -1
      ? (normalizeValue(row[sampleIdx]) || `Sample ${rowIdx + 1}`)
      : 'All samples';
    const refCount = refCountIdx !== -1 ? parseNumericValue(row[refCountIdx]) : null;
    const altCount = altCountIdx !== -1 ? parseNumericValue(row[altCountIdx]) : null;
    const coverage = Math.max(0, (refCount || 0) + (altCount || 0));

    let altFraction = altFractionIdx !== -1 ? parseNumericValue(row[altFractionIdx]) : null;
    if (altFraction === null && coverage > 0 && altCount !== null) {
      altFraction = (altCount / coverage) * 100;
    }

    if (coverage <= 0 && altFraction === null) return;

    const markerValue = markerIdx !== -1 ? normalizeValue(row[markerIdx]) : '';
    const markerLabel = markerValue || `Locus ${rowIdx + 1}`;
    const positionInfo = parseSplitFastqPosition(markerValue, rowIdx + 1);
    if (!positionInfo.inferred) numericPositionCount += 1;

    lociWithReadSupport += 1;
    totalCoverage += coverage;
    if (altCount !== null) totalAltReads += altCount;

    if (altFraction !== null) {
      altFractions.push(altFraction);
      markerRows.push({
        marker: markerLabel,
        sample: sampleName,
        altFraction,
        coverage
      });
    }

    if (!sampleStats.has(sampleName)) {
      sampleStats.set(sampleName, {
        sample: sampleName,
        loci: 0,
        coverage: 0,
        altReads: 0,
        altFractions: []
      });
    }

    const stats = sampleStats.get(sampleName);
    stats.loci += 1;
    stats.coverage += coverage;
    if (altCount !== null) stats.altReads += altCount;
    if (altFraction !== null) stats.altFractions.push(altFraction);

    const supportPct = Number.isFinite(altFraction)
      ? altFraction
      : (coverage > 0 && altCount !== null ? (altCount / coverage) * 100 : null);
    if (!Number.isFinite(supportPct)) return;

    trackItems.push({
      id: `splitfq-d-${rowIdx}`,
      mode: 'detailed',
      x: positionInfo.value,
      sample: sampleName,
      label: markerLabel,
      lineage: lineageIdx !== -1 ? normalizeValue(row[lineageIdx]) : '',
      supportPct,
      coverage,
      refCount: refCount || 0,
      altCount: altCount || 0,
      totalEvidence: coverage,
      majorCount: altCount || 0
    });
  });

  if (lociWithReadSupport === 0 || trackItems.length === 0) return null;

  const weightedAlt = totalCoverage > 0 ? (totalAltReads / totalCoverage) * 100 : null;
  const medianAlt = quantile(altFractions, 0.5);
  const highAltLoci = altFractions.filter(value => value >= 95).length;
  const coverages = trackItems.map(item => item.coverage).filter(c => c > 0);
  const medianDepth = coverages.length > 0 ? quantile(coverages, 0.5) : null;
  const lowDepthCount = trackItems.filter(item => item.coverage > 0 && item.coverage < 10).length;

  const sampleRows = [...sampleStats.values()]
    .map(stats => {
      const supportPct = stats.coverage > 0
        ? (stats.altReads / stats.coverage) * 100
        : quantile(stats.altFractions, 0.5);
      return {
        sample: stats.sample,
        loci: stats.loci,
        coverage: stats.coverage,
        supportPct
      };
    })
    .sort((a, b) => {
      const aSupport = Number.isFinite(a.supportPct) ? a.supportPct : -1;
      const bSupport = Number.isFinite(b.supportPct) ? b.supportPct : -1;
      return bSupport - aSupport || b.coverage - a.coverage || a.sample.localeCompare(b.sample, undefined, { sensitivity: 'base' });
    });

  const topMarkers = markerRows
    .sort((a, b) => b.altFraction - a.altFraction || b.coverage - a.coverage)
    .slice(0, 6);
  const sampleOptions = [...new Set(trackItems.map(item => item.sample))]
    .sort((a, b) => a.localeCompare(b, undefined, { sensitivity: 'base' }));
  const domain = splitfqTrackGetDomain(trackItems);
  const useBpAxis = numericPositionCount >= Math.max(2, Math.round(trackItems.length * 0.6));

  trackItems.sort((a, b) => a.x - b.x || a.id.localeCompare(b.id));

  return {
    mode: 'detailed',
    heading: 'Read Support Track',
    caption: 'Interactive marker-level support from ref/alt read counts.',
    kpis: [
      {
        label: 'Loci with Reads',
        value: formatIntegerForViz(lociWithReadSupport),
        sub: `${formatIntegerForViz(rows.length)} rows parsed`
      },
      {
        label: 'Weighted ALT',
        value: weightedAlt === null ? 'n/a' : `${formatPlotDecimal(weightedAlt, 1)}%`,
        sub: `${formatIntegerForViz(totalCoverage)} total reads`
      },
      {
        label: 'Median ALT',
        value: medianAlt === null ? 'n/a' : `${formatPlotDecimal(medianAlt, 1)}%`,
        sub: `${formatIntegerForViz(altFractions.length)} loci with ALT%`
      },
      {
        label: 'High ALT Loci',
        value: formatIntegerForViz(highAltLoci),
        sub: 'ALT >= 95%'
      },
      {
        label: 'Median Depth',
        value: medianDepth === null ? 'n/a' : formatIntegerForViz(Math.round(medianDepth)),
        sub: `${formatIntegerForViz(coverages.length)} loci with reads`
      },
      {
        label: 'Low-Depth Markers',
        value: formatIntegerForViz(lowDepthCount),
        sub: 'coverage < 10×'
      }
    ],
    sampleRows: sampleRows.slice(0, 8),
    markerRows: topMarkers,
    topLineages: [],
    totalEvidence: totalCoverage,
    sampleOptions,
    trackItems,
    trackDomainStart: domain.start,
    trackDomainEnd: domain.end,
    trackAxisLabel: useBpAxis ? 'Marker position (bp)' : 'Marker index (result order)',
    trackUseBpAxis: useBpAxis
  };
}

function buildSplitFastqSummaryModel(headers, rows, primaryColumn) {
  const lineageCountIdx = findHeaderIndexByTokens(headers, [
    'lineage:count',
    'lineage_count',
    'lineage_counts',
    'lineage_call_counts'
  ]);
  let majorIdx = findHeaderIndexByTokens(headers, ['major_lineage', primaryColumn || 'major_lineage'], { exactOnly: true });
  if (majorIdx === -1 && primaryColumn) {
    majorIdx = findHeaderIndexByTokens(headers, [primaryColumn]);
  }
  const sampleIdx = findSampleColumnIndex(headers, majorIdx);

  if (lineageCountIdx === -1 && majorIdx === -1) return null;

  const lineageTotals = new Map();
  const sampleRows = [];
  const supportValues = [];
  const trackItems = [];
  let lowSupportCount = 0;
  let multihitCount = 0;

  rows.forEach((row, rowIdx) => {
    const sampleName = sampleIdx !== -1
      ? (normalizeValue(row[sampleIdx]) || `Sample ${rowIdx + 1}`)
      : `Sample ${rowIdx + 1}`;
    const lineageEntries = lineageCountIdx !== -1 ? parseSplitFastqLineageCounts(row[lineageCountIdx]) : [];
    lineageEntries.forEach(entry => {
      lineageTotals.set(entry.lineage, (lineageTotals.get(entry.lineage) || 0) + entry.count);
    });

    let majorLineage = majorIdx !== -1 ? normalizeValue(row[majorIdx]) : '';
    if (!majorLineage && lineageEntries.length > 0) {
      majorLineage = lineageEntries[0].lineage;
    }
    if (!majorLineage) majorLineage = 'Unresolved';

    if (lineageEntries.length === 0 && !isUnknownOutcome(majorLineage)) {
      lineageTotals.set(majorLineage, (lineageTotals.get(majorLineage) || 0) + 1);
    }

    const totalEvidence = lineageEntries.reduce((sum, entry) => sum + entry.count, 0);
    const matched = lineageEntries.find(entry => normalizeHeader(entry.lineage) === normalizeHeader(majorLineage));
    const majorCount = matched
      ? matched.count
      : (lineageEntries.length > 0 ? lineageEntries[0].count : null);
    const supportPct = totalEvidence > 0 && majorCount !== null
      ? (majorCount / totalEvidence) * 100
      : null;

    if (lineageEntries.length > 1) multihitCount += 1;
    if (Number.isFinite(supportPct)) {
      supportValues.push(supportPct);
      if (supportPct < 60) lowSupportCount += 1;
    }

    sampleRows.push({
      sample: sampleName,
      majorLineage,
      majorCount,
      totalEvidence,
      supportPct
    });

    trackItems.push({
      id: `splitfq-s-${rowIdx}`,
      mode: 'summary',
      x: rowIdx + 1,
      sample: sampleName,
      label: sampleName,
      lineage: majorLineage,
      supportPct: Number.isFinite(supportPct) ? supportPct : 0,
      coverage: totalEvidence,
      refCount: 0,
      altCount: majorCount || 0,
      totalEvidence,
      majorCount: majorCount || 0
    });
  });

  if (sampleRows.length === 0 || trackItems.length === 0) return null;

  const medianSupport = quantile(supportValues, 0.5);
  const topLineages = [...lineageTotals.entries()]
    .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0], undefined, { sensitivity: 'base' }))
    .slice(0, 5);
  const rankedRows = [...sampleRows]
    .sort((a, b) => {
      const aSupport = Number.isFinite(a.supportPct) ? a.supportPct : 101;
      const bSupport = Number.isFinite(b.supportPct) ? b.supportPct : 101;
      return aSupport - bSupport || b.totalEvidence - a.totalEvidence || a.sample.localeCompare(b.sample, undefined, { sensitivity: 'base' });
    })
    .slice(0, 8);
  const totalEvidence = [...lineageTotals.values()].reduce((acc, value) => acc + value, 0);

  return {
    mode: 'summary',
    heading: 'Lineage Evidence Track',
    caption: 'Interactive confidence track derived from lineage:count evidence.',
    kpis: [
      {
        label: 'Samples',
        value: formatIntegerForViz(sampleRows.length),
        sub: `${formatIntegerForViz(rows.length)} rows`
      },
      {
        label: 'Median Major Support',
        value: medianSupport === null ? 'n/a' : `${formatPlotDecimal(medianSupport, 1)}%`,
        sub: `${formatIntegerForViz(supportValues.length)} samples with marker counts`
      },
      {
        label: 'Low Support',
        value: formatIntegerForViz(lowSupportCount),
        sub: 'major support < 60%'
      },
      {
        label: 'Multi-lineage',
        value: formatIntegerForViz(multihitCount),
        sub: 'samples with >1 candidate'
      }
    ],
    sampleRows: rankedRows,
    markerRows: [],
    topLineages,
    totalEvidence,
    sampleOptions: [...new Set(trackItems.map(item => item.sample))]
      .sort((a, b) => a.localeCompare(b, undefined, { sensitivity: 'base' })),
    trackItems,
    trackDomainStart: 1,
    trackDomainEnd: Math.max(2, trackItems.length),
    trackAxisLabel: 'Input row # (1-based)',
    trackUseBpAxis: false
  };
}

function buildSplitFastqInsightModel(data, primaryColumn) {
  const headers = Array.isArray(data?.headers) ? data.headers : [];
  const rows = Array.isArray(data?.rows) ? data.rows : [];
  if (headers.length === 0 || rows.length === 0) return null;

  const detailed = buildSplitFastqDetailedModel(headers, rows);
  if (detailed) return detailed;
  return buildSplitFastqSummaryModel(headers, rows, primaryColumn);
}

function ensureSplitfqTrackState(toolId, model) {
  const modelKey = splitfqTrackBuildModelKey(model);
  const optionSet = new Set(['__all__', ...(Array.isArray(model.sampleOptions) ? model.sampleOptions : [])]);
  const existing = SPLITFQ_TRACK_STATE[toolId];

  if (!existing || existing.modelKey !== modelKey) {
    const sampleFilter = existing && optionSet.has(existing.sampleFilter)
      ? existing.sampleFilter
      : '__all__';
    SPLITFQ_TRACK_STATE[toolId] = {
      modelKey,
      sampleFilter,
      activeFilter: '',
      domainStart: model.trackDomainStart,
      domainEnd: model.trackDomainEnd,
      windowStart: model.trackDomainStart,
      windowEnd: model.trackDomainEnd,
      selectedItemId: '',
      searchQuery: existing?.searchQuery || '',
      searchMatchIds: [],
      searchMatchIndex: -1
    };
  } else if (!optionSet.has(existing.sampleFilter)) {
    existing.sampleFilter = '__all__';
  }

  return SPLITFQ_TRACK_STATE[toolId];
}

function scheduleSplitfqTrackRender(toolId, data, primaryColumn) {
  if (SPLITFQ_TRACK_RENDER_THROTTLES[toolId]) return;
  const schedule = (typeof window !== 'undefined' && typeof window.requestAnimationFrame === 'function')
    ? window.requestAnimationFrame.bind(window)
    : (cb => setTimeout(cb, 16));
  SPLITFQ_TRACK_RENDER_THROTTLES[toolId] = schedule(() => {
    delete SPLITFQ_TRACK_RENDER_THROTTLES[toolId];
    const breakdownEl = document.getElementById(`${toolId}-viz-outcome-breakdown`);
    if (!breakdownEl) return;
    renderSplitFastqLightweightInsights(toolId, data, primaryColumn, breakdownEl);
    renderClassifyGenotypingInsights(toolId, data, breakdownEl);
  });
}

function attachSplitfqTrackInteractions(toolId, data, primaryColumn, model, state, items, domain, section) {
  if (!section) return;

  const sampleSelect = section.querySelector('.viz-splitfq-track-sample-select');
  const searchInput = section.querySelector('.viz-splitfq-track-search-input');
  const startInput = section.querySelector('.viz-splitfq-track-start-input');
  const endInput = section.querySelector('.viz-splitfq-track-end-input');
  const svgEl = section.querySelector('.viz-splitfq-track-svg');
  const findBtn = section.querySelector('[data-splitfq-action="find"]');
  const applyRangeBtn = section.querySelector('[data-splitfq-action="apply-range"]');

  sampleSelect?.addEventListener('change', () => {
    state.sampleFilter = sampleSelect.value || '__all__';
    state.activeFilter = '';
    state.selectedItemId = '';
    state.searchMatchIndex = -1;
    scheduleSplitfqTrackRender(toolId, data, primaryColumn);
  });

  searchInput?.addEventListener('input', () => {
    state.searchQuery = searchInput.value || '';
    state.searchMatchIndex = 0;
    scheduleSplitfqTrackRender(toolId, data, primaryColumn);
  });
  searchInput?.addEventListener('keydown', event => {
    if (event.key !== 'Enter') return;
    event.preventDefault();
    findBtn?.click();
  });
  [startInput, endInput].forEach(input => {
    input?.addEventListener('keydown', event => {
      if (event.key !== 'Enter') return;
      event.preventDefault();
      applyRangeBtn?.click();
    });
  });

  const moveToMatch = (direction = 1, fromCurrent = true) => {
    if (!Array.isArray(state.searchMatchIds) || state.searchMatchIds.length === 0) return;
    let index = Math.max(0, Number(state.searchMatchIndex || 0));
    if (fromCurrent) {
      index = (index + direction + state.searchMatchIds.length) % state.searchMatchIds.length;
    }
    state.searchMatchIndex = index;
    const matchId = state.searchMatchIds[index];
    const matchItem = splitfqTrackFindItemById(items, matchId);
    if (matchItem) splitfqTrackFocusItem(state, matchItem, domain.start, domain.end, true);
  };

  section.querySelectorAll('[data-splitfq-action]').forEach(button => {
    button.addEventListener('click', () => {
      const action = button.dataset.splitfqAction;
      if (!action) return;

      switch (action) {
        case 'zoom-in':
          splitfqTrackZoom(state, 0.7, null, domain.start, domain.end);
          break;
        case 'zoom-out':
          splitfqTrackZoom(state, 1.35, null, domain.start, domain.end);
          break;
        case 'pan-left':
          splitfqTrackPan(state, -0.2, domain.start, domain.end);
          break;
        case 'pan-right':
          splitfqTrackPan(state, 0.2, domain.start, domain.end);
          break;
        case 'reset': {
          const span = splitfqTrackDefaultSpan(model.mode, domain.start, domain.end, items.length);
          splitfqTrackSetWindow(state, domain.start, domain.start + span, domain.start, domain.end);
          state.selectedItemId = '';
          break;
        }
        case 'apply-range': {
          const start = Number.parseFloat(startInput?.value || '');
          const end = Number.parseFloat(endInput?.value || '');
          if (Number.isFinite(start) && Number.isFinite(end)) {
            splitfqTrackSetWindow(state, start, end, domain.start, domain.end);
          }
          break;
        }
        case 'find': {
          if (state.searchMatchIds.length === 0 && String(state.searchQuery || '').trim()) {
            state.searchMatchIds = splitfqTrackFindMatches(items, state.searchQuery).map(item => item.id);
            state.searchMatchIndex = 0;
          }
          moveToMatch(0, false);
          break;
        }
        case 'search-next':
          moveToMatch(1, true);
          break;
        case 'search-prev':
          moveToMatch(-1, true);
          break;
        default:
          return;
      }

      scheduleSplitfqTrackRender(toolId, data, primaryColumn);
    });
  });

  if (!svgEl) return;

  svgEl.addEventListener('click', event => {
    const pointHit = event.target.closest('.viz-splitfq-point-hit');
    if (pointHit) {
      const itemId = pointHit.getAttribute('data-item-id') || '';
      const item = splitfqTrackFindItemById(items, itemId);
      if (!item) return;
      splitfqTrackFocusItem(state, item, domain.start, domain.end, true);
      scheduleSplitfqTrackRender(toolId, data, primaryColumn);
      return;
    }

    const role = event.target.getAttribute('data-splitfq-role');
    if (role === 'overview-hit') {
      const center = splitfqTrackDomainValueFromEvent(event, svgEl, state, true);
      const span = Math.max(1, state.windowEnd - state.windowStart);
      splitfqTrackSetWindow(
        state,
        center - span / 2,
        center + span / 2,
        domain.start,
        domain.end
      );
      scheduleSplitfqTrackRender(toolId, data, primaryColumn);
    }
  });

  svgEl.addEventListener('wheel', event => {
    event.preventDefault();
    const center = splitfqTrackDomainValueFromEvent(event, svgEl, state, false);
    splitfqTrackZoom(
      state,
      event.deltaY < 0 ? 0.78 : 1.25,
      center,
      domain.start,
      domain.end
    );
    scheduleSplitfqTrackRender(toolId, data, primaryColumn);
  }, { passive: false });

  svgEl.addEventListener('mousedown', event => {
    if (event.button !== 0) return;
    if (event.target.closest('.viz-splitfq-point-hit')) return;
    if (event.target.getAttribute('data-splitfq-role') === 'overview-hit') return;

    const viewWidth = Number(svgEl.dataset.viewWidth || 980);
    const plotWidth = Math.max(1, Number(svgEl.dataset.plotWidth || 1));
    const startViewX = getSvgViewXFromClient(svgEl, event.clientX, viewWidth / 2);
    const spanAtStart = Math.max(1, state.windowEnd - state.windowStart);
    const drag = {
      startViewX,
      startWindowStart: state.windowStart,
      startWindowEnd: state.windowEnd,
      plotWidth,
      spanAtStart
    };

    document.body.classList.add('track-dragging');
    event.preventDefault();

    const onMove = moveEvent => {
      const currentViewX = getSvgViewXFromClient(svgEl, moveEvent.clientX, drag.startViewX);
      const deltaView = currentViewX - drag.startViewX;
      const deltaDomain = (deltaView / drag.plotWidth) * drag.spanAtStart;
      splitfqTrackSetWindow(
        state,
        drag.startWindowStart - deltaDomain,
        drag.startWindowEnd - deltaDomain,
        domain.start,
        domain.end
      );
      scheduleSplitfqTrackRender(toolId, data, primaryColumn);
    };

    const onUp = () => {
      document.body.classList.remove('track-dragging');
      window.removeEventListener('mousemove', onMove);
      window.removeEventListener('mouseup', onUp);
    };

    window.addEventListener('mousemove', onMove);
    window.addEventListener('mouseup', onUp);
  });
}

function renderPredictConfidenceInsights(toolId, data, breakdownEl) {
  if (toolId !== 'predict' || !breakdownEl) return;
  const headers = Array.isArray(data?.headers) ? data.headers : [];
  const rows = Array.isArray(data?.rows) ? data.rows : [];
  if (rows.length === 0) return;

  const confIdx = findHeaderByTokens(headers, ['confidence', 'probability', 'score']);
  if (confIdx === -1) return;

  const values = rows.map(row => parseNumericValue(row[confIdx])).filter(v => v !== null);
  if (values.length === 0) return;

  // Confidence histogram bins
  const bins = [
    { label: 'Very Low (<0.5)', min: -Infinity, max: 0.5, color: '#ef4444', count: 0 },
    { label: 'Low (0.5–0.7)', min: 0.5, max: 0.7, color: '#f59e0b', count: 0 },
    { label: 'Good (0.7–0.9)', min: 0.7, max: 0.9, color: '#eab308', count: 0 },
    { label: 'High (≥0.9)', min: 0.9, max: Infinity, color: '#22c55e', count: 0 }
  ];

  for (const v of values) {
    for (const bin of bins) {
      if (v >= bin.min && v < bin.max) { bin.count++; break; }
    }
  }
  // Handle edge case: value === Infinity or exactly max
  if (bins[3].count === 0) {
    const exactMax = values.filter(v => v >= 0.9).length;
    if (exactMax > 0) bins[3].count = exactMax;
  }

  const maxCount = Math.max(1, ...bins.map(b => b.count));
  const lowCount = bins[0].count + bins[1].count;
  const meanConf = values.reduce((a, b) => a + b, 0) / values.length;

  const existingSection = breakdownEl.querySelector('.viz-predict-confidence');
  const html = `
    <div class="viz-predict-confidence" style="margin-top:12px;">
      <div class="viz-outcome-header">
        <h4>Confidence Distribution</h4>
        <span>Mean: ${formatPlotDecimal(meanConf, 3)} · ${lowCount > 0 ? `${lowCount} low-confidence` : 'All confident'}</span>
      </div>
      <div class="viz-outcome-list">
        ${bins.map(bin => {
          const widthPct = (bin.count / maxCount) * 100;
          const sharePct = values.length > 0 ? (bin.count / values.length) * 100 : 0;
          return `
            <div class="viz-outcome-item">
              <div class="viz-outcome-meta">
                <span class="viz-outcome-label">${escapeHtml(bin.label)}</span>
                <span class="viz-outcome-value">${bin.count} (${formatPlotDecimal(sharePct, 1)}%)</span>
              </div>
              <div class="viz-outcome-track">
                <span class="viz-outcome-fill" style="width:${formatPlotDecimal(widthPct, 2)}%;background:${bin.color}"></span>
              </div>
            </div>
          `;
        }).join('')}
      </div>
    </div>
  `;

  if (existingSection) {
    existingSection.outerHTML = html;
  } else {
    breakdownEl.insertAdjacentHTML('beforeend', html);
  }
}

/**
 * Build the resistance / lineage-level / allele-fraction panels straight from
 * the loaded table, so they do not depend on the insight model's internal shape.
 *
 * The same output format carries both lineage and drug-resistance panels: the
 * marker's trailing columns are joined into `lineage_path`. For the WHO
 * resistance catalogue that path is drug;resistance;marker;grade;gene;mutation,
 * which is detected and rendered as a clinical profile instead of a lineage.
 */
function buildGenotypingInsightSections(data) {
  const headers = Array.isArray(data?.headers) ? data.headers : [];
  const rows = Array.isArray(data?.rows) ? data.rows : [];
  if (headers.length === 0 || rows.length === 0) return '';

  const lineageIdx = findHeaderIndexByTokens(headers, ['lineage_path', 'lineage', 'major_lineage']);
  const lineageCountIdx = findHeaderIndexByTokens(headers, ['lineage:count', 'lineage_count']);
  const fracIdx = findHeaderIndexByTokens(headers, ['alt_fraction', 'alt_percent', 'alt_pct', 'vaf', 'allele_fraction']);
  const refIdx = findHeaderIndexByTokens(headers, ['ref_count', 'reference_count', 'ref_reads']);
  const altIdx = findHeaderIndexByTokens(headers, ['alt_count', 'alternate_count', 'alt_reads', 'mut_count']);

  // Summary tables carry no per-marker evidence, only the lineage:count cell.
  if (lineageIdx === -1 && lineageCountIdx === -1) return '';

  const sampleIdx = findSampleColumnIndex(headers, lineageIdx);
  const records = [];
  const fractions = [];
  if (lineageIdx !== -1) {
    rows.forEach(row => {
      const lineagePath = normalizeValue(row[lineageIdx]);
      if (!lineagePath) return;
      const sample = sampleIdx !== -1 ? (normalizeValue(row[sampleIdx]) || 'Sample') : 'Sample';
      const refCount = refIdx !== -1 ? parseNumericValue(row[refIdx]) : null;
      const altCount = altIdx !== -1 ? parseNumericValue(row[altIdx]) : null;
      const coverage = (refCount || 0) + (altCount || 0);
      let altFraction = fracIdx !== -1 ? parseNumericValue(row[fracIdx]) : null;
      if (altFraction === null && coverage > 0 && altCount !== null) {
        altFraction = (altCount / coverage) * 100;
      }
      records.push({ sample, lineagePath, altFraction, coverage: coverage > 0 ? coverage : null });
      if (Number.isFinite(altFraction)) fractions.push(altFraction);
    });
  }

  let html = '';

  if (records.length > 0 && isDrPanel(records.map(r => r.lineagePath))) {
    // With several samples a matrix reads far better than a list of cards, and
    // the cards would silently merge every sample's mutations into one profile.
    // buildResistanceMatrix returns null for a single sample, where the cards
    // are the right form and a one-row grid would not be.
    const matrix = buildResistanceMatrix(records);
    if (matrix) {
      html += renderResistanceMatrixHtml(matrix);
    } else {
      html += renderDrProfileHtml(buildDrProfile(records));
    }
  } else {
    // Genuine lineage panel: rebuild the hierarchy the flat paths describe.
    const counts = new Map();
    if (records.length > 0) {
      records.forEach(r => counts.set(r.lineagePath, (counts.get(r.lineagePath) || 0) + 1));
    } else if (lineageCountIdx !== -1) {
      rows.forEach(row => {
        parseSplitFastqLineageCounts(row[lineageCountIdx]).forEach(entry => {
          counts.set(entry.lineage, (counts.get(entry.lineage) || 0) + entry.count);
        });
      });
    }
    if (records.length > 0) {
      const bySample = new Map();
      records.forEach(r => {
        if (!bySample.has(r.sample)) bySample.set(r.sample, new Map());
        const m = bySample.get(r.sample);
        m.set(r.lineagePath, (m.get(r.lineagePath) || 0) + 1);
      });
      const multiSample = bySample.size > 1;
      if (multiSample) {
        const perSample = buildLineageComposition(
          new Map([...bySample].map(([sample, m]) =>
            [sample, [...m].map(([lineage, count]) => ({ lineage, count }))]))
        );
        html += renderLineageCompositionHtml(perSample);
      }
    }
    if (counts.size > 0) {
      const entries = [...counts].map(([lineage, count]) => ({ lineage, count }));
      // The mixed-infection verdict describes a single sample; pooling a batch
      // would count every sample's lineage as a co-infecting strain. For a
      // batch the composition chart above already carries that story.
      if (!(records.length > 0 && new Set(records.map(r => r.sample)).size > 1)) {
        html += renderMixedLineagesHtml(detectMixedLineages(buildLineageBranches(entries)));
      }
      html += renderLineageLevelsHtml(buildLineageLevels(entries));
    }
  }

  if (fractions.length > 0) {
    // The core only writes variants at or above --min-alt-percent, so the
    // lowest fraction present tells us where the data was truncated.
    const observedMin = fractions.reduce((m, v) => Math.min(m, v), Infinity);
    html += renderAlleleHistogramHtml(buildAlleleHistogram(fractions, 10), {
      minAltPercent: Number.isFinite(observedMin) ? observedMin : null
    });
  }

  return html;
}

/**
 * Classify runs the very same marker panels as split-fastq — including the WHO
 * resistance catalogue — so it gets the same resistance and lineage panels.
 * Assembly output carries no read depth, so the allele-fraction panel simply
 * does not appear for it.
 */
function renderClassifyGenotypingInsights(toolId, data, breakdownEl) {
  if (toolId !== 'classify' || !breakdownEl) return;
  const existing = breakdownEl.querySelector('.viz-classify-genotyping');
  const html = buildGenotypingInsightSections(data);
  if (!html) {
    existing?.remove();
    return;
  }
  const section = existing || document.createElement('section');
  section.className = 'viz-classify-genotyping';
  section.innerHTML = html;
  if (!existing) breakdownEl.appendChild(section);
}

function renderSplitFastqLightweightInsights(toolId, data, primaryColumn, breakdownEl) {
  if (toolId !== 'splitfq' || !breakdownEl) return;
  const model = buildSplitFastqInsightModel(data, primaryColumn);
  const existingSection = breakdownEl.querySelector('.viz-splitfq-lightweight');
  if (!model) {
    existingSection?.remove();
    return;
  }

  const state = ensureSplitfqTrackState(toolId, model);
  let items = Array.isArray(model.trackItems) ? [...model.trackItems] : [];
  if (state.sampleFilter && state.sampleFilter !== '__all__') {
    items = items.filter(item => item.sample === state.sampleFilter);
  }
  if (items.length === 0) {
    state.sampleFilter = '__all__';
    items = Array.isArray(model.trackItems) ? [...model.trackItems] : [];
  }
  items.sort((a, b) => a.x - b.x || a.id.localeCompare(b.id));

  const domain = splitfqTrackGetDomain(items);
  const domainChanged = state.domainStart !== domain.start || state.domainEnd !== domain.end;
  const filterChanged = state.activeFilter !== state.sampleFilter;
  if (domainChanged || filterChanged) {
    const span = splitfqTrackDefaultSpan(model.mode, domain.start, domain.end, items.length);
    splitfqTrackSetWindow(state, domain.start, domain.start + span, domain.start, domain.end);
    if (filterChanged) state.selectedItemId = '';
  } else {
    splitfqTrackSetWindow(state, state.windowStart, state.windowEnd, domain.start, domain.end);
  }
  state.domainStart = domain.start;
  state.domainEnd = domain.end;
  state.activeFilter = state.sampleFilter;

  const query = String(state.searchQuery || '').trim();
  if (!query) {
    state.searchMatchIds = [];
    state.searchMatchIndex = -1;
  } else {
    state.searchMatchIds = splitfqTrackFindMatches(items, query).map(item => item.id);
    if (state.searchMatchIds.length === 0) {
      state.searchMatchIndex = -1;
    } else if (state.searchMatchIndex < 0 || state.searchMatchIndex >= state.searchMatchIds.length) {
      state.searchMatchIndex = 0;
    }
    const selectedIndex = state.searchMatchIds.indexOf(state.selectedItemId);
    if (selectedIndex !== -1) state.searchMatchIndex = selectedIndex;
  }

  if (state.selectedItemId && !splitfqTrackFindItemById(items, state.selectedItemId)) {
    state.selectedItemId = '';
  }

  const trackInfo = buildSplitfqTrackSvg(model, state, items);
  const topLineages = Array.isArray(model.topLineages) ? model.topLineages : [];
  const sampleOptions = Array.isArray(model.sampleOptions) ? model.sampleOptions : [];
  const searchStatus = splitfqTrackSearchStatus(state);
  const searchOptions = splitfqTrackBuildSearchOptions(items);
  const windowStartValue = Math.round(state.windowStart);
  const windowEndValue = Math.round(state.windowEnd);
  const selectedDetail = buildSplitfqTrackDetailHtml(model, trackInfo.selectedItem, trackInfo.visibleCount);
  const splitfqWindowLabel = model.trackUseBpAxis
    ? 'Window start/end (bp)'
    : (model.mode === 'summary' ? 'Window start/end (input row #)' : 'Window start/end (marker index)');
  const splitfqXAxisHint = model.trackUseBpAxis
    ? 'X-axis: marker position in reference genome (bp).'
    : (model.mode === 'summary'
      ? 'X-axis: input row number from the loaded results table (1-based, stable after sorting/filtering).'
      : 'X-axis: marker index (1-based, order in result rows).');

  const genotypingSections = buildGenotypingInsightSections(data);

  const section = existingSection || document.createElement('section');
  section.className = 'viz-splitfq-lightweight';
  section.innerHTML = `
    <div class="viz-splitfq-head">
      <h5>${escapeHtml(model.heading || 'Split FASTQ Evidence')}</h5>
      <span>${escapeHtml(model.caption || '')}</span>
    </div>
    ${genotypingSections}
    <div class="viz-splitfq-kpis">
      ${(model.kpis || []).map(kpi => `
        <article class="viz-splitfq-kpi">
          <span class="viz-splitfq-kpi-label">${escapeHtml(kpi.label || '')}</span>
          <strong class="viz-splitfq-kpi-value">${escapeHtml(String(kpi.value || '0'))}</strong>
          <span class="viz-splitfq-kpi-sub">${escapeHtml(kpi.sub || '')}</span>
        </article>
      `).join('')}
    </div>
    ${topLineages.length > 0 ? `
      <div class="viz-splitfq-lineages">
        ${topLineages.map(([lineage, value]) => {
          const share = model.totalEvidence > 0 ? (value / model.totalEvidence) * 100 : 0;
          return `
            <span class="viz-splitfq-lineage-chip" title="${escapeHtml(lineage)}">
              <span>${escapeHtml(lineage)}</span>
              <strong>${formatIntegerForViz(value)}</strong>
              <em>${formatPlotDecimal(share, 1)}%</em>
            </span>
          `;
        }).join('')}
      </div>
    ` : ''}
    <div class="viz-splitfq-track-shell">
      <div class="viz-splitfq-track-controls">
        <div class="viz-splitfq-track-toolbar-top">
          <label class="viz-splitfq-track-field">
            <span class="viz-splitfq-track-field-label">Sample</span>
            <select class="viz-splitfq-track-sample-select">
              <option value="__all__" ${state.sampleFilter === '__all__' ? 'selected' : ''}>All samples</option>
              ${sampleOptions.map(sample => `
                <option value="${escapeHtml(sample)}" ${sample === state.sampleFilter ? 'selected' : ''}>${escapeHtml(sample)}</option>
              `).join('')}
            </select>
          </label>
          <div class="viz-splitfq-track-field viz-splitfq-track-search">
            <span class="viz-splitfq-track-field-label">Find SNP/read support</span>
            <div class="viz-splitfq-track-search-wrap">
              <input
                type="text"
                class="viz-splitfq-track-search-input"
                placeholder="Sample, marker, lineage, or position"
                value="${escapeHtml(state.searchQuery || '')}"
                list="${toolId}-viz-splitfq-search-options"
                autocomplete="off"
              />
              <datalist id="${toolId}-viz-splitfq-search-options">
                ${searchOptions}
              </datalist>
              <div class="viz-splitfq-track-search-actions">
                <button type="button" class="viz-track-btn viz-track-btn-search" data-splitfq-action="find">
                  <span class="viz-track-btn-icon" aria-hidden="true">⌕</span>
                  <span class="viz-track-btn-label">Find</span>
                </button>
                <button type="button" class="viz-track-btn viz-track-btn-search-nav" data-splitfq-action="search-prev" ${state.searchMatchIds.length > 0 ? '' : 'disabled'}>
                  <span class="viz-track-btn-icon" aria-hidden="true">‹</span>
                  <span class="viz-track-btn-label">Prev</span>
                </button>
                <button type="button" class="viz-track-btn viz-track-btn-search-nav" data-splitfq-action="search-next" ${state.searchMatchIds.length > 0 ? '' : 'disabled'}>
                  <span class="viz-track-btn-icon" aria-hidden="true">›</span>
                  <span class="viz-track-btn-label">Next</span>
                </button>
              </div>
            </div>
            <span class="viz-splitfq-track-search-status" aria-live="polite">${escapeHtml(searchStatus)}</span>
          </div>
        </div>
        <div class="viz-splitfq-track-toolbar-bottom">
          <label class="viz-splitfq-track-field viz-splitfq-track-window">
            <span class="viz-splitfq-track-field-label">${escapeHtml(splitfqWindowLabel)}</span>
            <div class="viz-splitfq-track-window-wrap">
              <input
                type="number"
                class="viz-splitfq-track-start-input"
                value="${windowStartValue}"
                step="1"
                min="${Math.floor(domain.start)}"
                max="${Math.ceil(domain.end)}"
              />
              <span class="viz-splitfq-track-window-sep">to</span>
              <input
                type="number"
                class="viz-splitfq-track-end-input"
                value="${windowEndValue}"
                step="1"
                min="${Math.floor(domain.start)}"
                max="${Math.ceil(domain.end)}"
              />
              <button type="button" class="viz-track-btn viz-track-btn-primary" data-splitfq-action="apply-range">
                <span class="viz-track-btn-icon" aria-hidden="true">✓</span>
                <span class="viz-track-btn-label">Apply</span>
              </button>
            </div>
          </label>
          <div class="viz-splitfq-track-nav" role="group" aria-label="Split fastq track navigation">
            <div class="viz-track-btn-group">
              <button type="button" class="viz-track-btn" data-splitfq-action="zoom-in">
                <span class="viz-track-btn-icon" aria-hidden="true">＋</span>
                <span class="viz-track-btn-label">Zoom In</span>
              </button>
              <button type="button" class="viz-track-btn" data-splitfq-action="zoom-out">
                <span class="viz-track-btn-icon" aria-hidden="true">－</span>
                <span class="viz-track-btn-label">Zoom Out</span>
              </button>
              <button type="button" class="viz-track-btn" data-splitfq-action="pan-left">
                <span class="viz-track-btn-icon" aria-hidden="true">←</span>
                <span class="viz-track-btn-label">Pan Left</span>
              </button>
              <button type="button" class="viz-track-btn" data-splitfq-action="pan-right">
                <span class="viz-track-btn-icon" aria-hidden="true">→</span>
                <span class="viz-track-btn-label">Pan Right</span>
              </button>
            </div>
            <button type="button" class="viz-track-btn" data-splitfq-action="reset">
              <span class="viz-track-btn-icon" aria-hidden="true">↺</span>
              <span class="viz-track-btn-label">Reset</span>
            </button>
          </div>
        </div>
      </div>
      <div class="viz-splitfq-track-meta">
        <span>${trackInfo.visibleCount} point(s) in window · ${items.length} in current filter</span>
        <span>${trackInfo.renderedCount < trackInfo.visibleCount ? `rendered ${trackInfo.renderedCount} for performance` : 'full resolution render'}</span>
        <span>${escapeHtml(splitfqXAxisHint)}</span>
        ${trackInfo.readLaneEnabled ? '<span>Read panel: stacked bars ALT/REF, height = total supporting reads</span>' : ''}
        <span>Controls: drag to pan, wheel to zoom, click point to focus</span>
      </div>
      <div class="viz-splitfq-track-canvas">
        ${trackInfo.svg}
      </div>
      <div class="viz-splitfq-track-detail">
        ${selectedDetail}
      </div>
    </div>
  `;

  if (!existingSection) {
    breakdownEl.appendChild(section);
  }

  attachSplitfqTrackInteractions(toolId, data, primaryColumn, model, state, items, domain, section);
}

function ensureMatchTrackState(toolId, model) {
  const modelKey = matchTrackBuildModelKey(model);
  const optionSet = new Set(['__all__', ...(Array.isArray(model.referenceOptions) ? model.referenceOptions : [])]);
  const existing = MATCH_TRACK_STATE[toolId];

  if (!existing || existing.modelKey !== modelKey) {
    const referenceFilter = existing && optionSet.has(existing.referenceFilter)
      ? existing.referenceFilter
      : '__all__';
    MATCH_TRACK_STATE[toolId] = {
      modelKey,
      referenceFilter,
      activeReferenceFilter: '',
      domainStart: model.trackDomainStart,
      domainEnd: model.trackDomainEnd,
      windowStart: model.trackDomainStart,
      windowEnd: model.trackDomainEnd,
      selectedPointId: '',
      searchQuery: existing?.searchQuery || '',
      searchMatchIds: [],
      searchMatchIndex: -1
    };
  } else if (!optionSet.has(existing.referenceFilter)) {
    existing.referenceFilter = '__all__';
  }

  return MATCH_TRACK_STATE[toolId];
}

function scheduleMatchTrackRender(toolId, data, primaryColumn) {
  if (MATCH_TRACK_RENDER_THROTTLES[toolId]) return;
  const schedule = (typeof window !== 'undefined' && typeof window.requestAnimationFrame === 'function')
    ? window.requestAnimationFrame.bind(window)
    : (cb => setTimeout(cb, 16));
  MATCH_TRACK_RENDER_THROTTLES[toolId] = schedule(() => {
    delete MATCH_TRACK_RENDER_THROTTLES[toolId];
    const breakdownEl = document.getElementById(`${toolId}-viz-outcome-breakdown`);
    if (!breakdownEl) return;
    renderMatchRefTrackInsights(toolId, data, primaryColumn, breakdownEl);
  });
}

function matchTrackGetDomain(points) {
  if (!Array.isArray(points) || points.length === 0) return { start: 0, end: 1 };
  let min = Number.POSITIVE_INFINITY;
  let max = Number.NEGATIVE_INFINITY;
  points.forEach(point => {
    const value = Number(point?.x);
    if (!Number.isFinite(value)) return;
    if (value < min) min = value;
    if (value > max) max = value;
  });
  if (!Number.isFinite(min) || !Number.isFinite(max)) {
    return { start: 0, end: 1 };
  }
  if (max <= min) max = min + 1;
  return { start: min, end: max };
}

function matchTrackSetWindow(state, start, end, domainStart, domainEnd) {
  const minBound = Number.isFinite(domainStart) ? domainStart : 0;
  let maxBound = Number.isFinite(domainEnd) ? domainEnd : minBound + 1;
  if (maxBound <= minBound) maxBound = minBound + 1;

  const fullSpan = Math.max(1, maxBound - minBound);
  const minSpan = Math.min(fullSpan, Math.max(1, fullSpan * MATCH_TRACK_MIN_SPAN_FACTOR));

  let windowStart = Number.isFinite(start) ? start : minBound;
  let windowEnd = Number.isFinite(end) ? end : maxBound;
  if (windowEnd < windowStart) {
    const tmp = windowStart;
    windowStart = windowEnd;
    windowEnd = tmp;
  }

  if ((windowEnd - windowStart) < minSpan) {
    const center = (windowStart + windowEnd) / 2;
    windowStart = center - (minSpan / 2);
    windowEnd = center + (minSpan / 2);
  }

  if (windowStart < minBound) {
    windowEnd += (minBound - windowStart);
    windowStart = minBound;
  }
  if (windowEnd > maxBound) {
    windowStart -= (windowEnd - maxBound);
    windowEnd = maxBound;
  }
  if (windowStart < minBound) windowStart = minBound;

  if ((windowEnd - windowStart) < minSpan) {
    windowEnd = Math.min(maxBound, windowStart + minSpan);
    windowStart = Math.max(minBound, windowEnd - minSpan);
  }

  state.windowStart = windowStart;
  state.windowEnd = windowEnd;
}

function matchTrackZoom(state, factor, center, domainStart, domainEnd) {
  const currentSpan = Math.max(1, state.windowEnd - state.windowStart);
  const targetSpan = currentSpan * factor;
  const anchor = Number.isFinite(center)
    ? center
    : (state.windowStart + state.windowEnd) / 2;
  matchTrackSetWindow(
    state,
    anchor - targetSpan / 2,
    anchor + targetSpan / 2,
    domainStart,
    domainEnd
  );
}

function matchTrackPan(state, fraction, domainStart, domainEnd) {
  const span = Math.max(1, state.windowEnd - state.windowStart);
  const shift = span * fraction;
  matchTrackSetWindow(
    state,
    state.windowStart + shift,
    state.windowEnd + shift,
    domainStart,
    domainEnd
  );
}

function matchTrackDownsamplePoints(points, maxPoints, selectedPointId = '') {
  if (!Array.isArray(points) || points.length <= maxPoints) return points || [];
  const stride = Math.ceil(points.length / maxPoints);
  const sampled = points.filter((_, index) => index % stride === 0);
  if (selectedPointId && !sampled.some(point => point.id === selectedPointId)) {
    const selected = points.find(point => point.id === selectedPointId);
    if (selected) sampled.push(selected);
  }
  return sampled.sort((a, b) => a.x - b.x || a.id.localeCompare(b.id));
}

function matchTrackFindItemById(points, pointId) {
  if (!pointId) return null;
  return points.find(point => point.id === pointId) || null;
}

function matchTrackItemSearchText(item) {
  return [
    item.sampleLabel,
    item.reference,
    item.query,
    Number.isFinite(item.x) ? `${Math.round(item.x)}` : '',
    Number.isFinite(item.score) ? `${formatPlotDecimal(item.score * 100, 2)}%` : '',
    `input row ${Number.isFinite(item.x) ? Math.round(item.x) : ''}`
  ].filter(Boolean).join(' ').toLowerCase();
}

function matchTrackFindMatches(points, query) {
  const normalized = String(query || '').trim().toLowerCase();
  if (!normalized) return [];
  return points.filter(point => matchTrackItemSearchText(point).includes(normalized));
}

function matchTrackBuildSearchOptions(points, limit = 180) {
  if (!Array.isArray(points) || points.length === 0) return '';
  const unique = new Set();
  const options = [];
  for (const point of points) {
    const value = `sample ${point.sampleLabel} · ${point.reference || 'Unassigned'} · ${formatPlotDecimal(point.score * 100, 2)}% · input row ${Math.round(point.x)}`;
    if (unique.has(value)) continue;
    unique.add(value);
    options.push(`<option value="${escapeHtml(value)}"></option>`);
    if (options.length >= limit) break;
  }
  return options.join('');
}

function matchTrackSearchStatus(state) {
  const query = String(state.searchQuery || '').trim();
  const matches = Array.isArray(state.searchMatchIds) ? state.searchMatchIds.length : 0;
  if (!query) return 'Type sample, reference, input row, or score.';
  if (matches === 0) return 'No matching points.';
  const index = Math.max(0, Number(state.searchMatchIndex || 0)) + 1;
  return `${index}/${matches} match${matches === 1 ? '' : 'es'}`;
}

function matchTrackDomainValueFromEvent(event, svgEl, state, useOverview = false) {
  if (!svgEl) return (state.windowStart + state.windowEnd) / 2;
  const viewWidth = Number(svgEl.dataset.viewWidth || 980);
  const plotLeft = Number(svgEl.dataset.plotLeft || 0);
  const plotWidth = Number(svgEl.dataset.plotWidth || 1);
  const x = getSvgViewXFromClient(svgEl, event.clientX, viewWidth / 2);
  const rel = Math.max(0, Math.min(1, (x - plotLeft) / Math.max(1, plotWidth)));
  // The overview strip spans the full domain (same plot geometry), so a click
  // on it must map through the full domain, not the current zoom window.
  const domainStart = useOverview ? state.domainStart : state.windowStart;
  const domainEnd = useOverview ? state.domainEnd : state.windowEnd;
  return domainStart + rel * (domainEnd - domainStart);
}

function matchTrackFocusItem(state, point, domainStart, domainEnd, zoom = true) {
  if (!point) return;
  state.selectedPointId = point.id;
  if (!zoom) {
    if (point.x < state.windowStart || point.x > state.windowEnd) {
      const span = Math.max(1, state.windowEnd - state.windowStart);
      matchTrackSetWindow(state, point.x - span / 2, point.x + span / 2, domainStart, domainEnd);
    }
    return;
  }
  const span = Math.max(1, state.windowEnd - state.windowStart);
  const targetSpan = span * 0.45;
  matchTrackSetWindow(state, point.x - targetSpan / 2, point.x + targetSpan / 2, domainStart, domainEnd);
}

function formatReferenceLabel(label) {
  const trimmed = normalizeValue(label);
  return trimmed || 'Unassigned';
}

function resolveMatchTrackReferenceColor(point, colorMap) {
  if (!colorMap) return '#64748b';
  return colorMap[point.reference] || '#64748b';
}

function buildMatchTrackModel(data, primaryColumn) {
  const headers = Array.isArray(data?.headers) ? data.headers : [];
  const rows = Array.isArray(data?.rows) ? data.rows : [];
  if (headers.length === 0 || rows.length === 0) return null;

  const queryIdx = findHeaderIndexByTokens(headers, [
    'query_files',
    'query_file',
    'query',
    'sample',
    'file',
    'genome',
    'input'
  ]);
  const referenceIdx = findHeaderIndexByTokens(headers, [
    'best_match_reference',
    'best_reference',
    'best_match',
    'reference',
    'match_reference'
  ]);
  const scoreIdx = findHeaderIndexByTokens(headers, [
    'shared_kmer_fraction',
    'match_score',
    'identity',
    'similarity',
    'ani',
    'coverage'
  ]);
  const fallbackReferenceIdx = findHeaderIndexByTokens(headers, [primaryColumn]);
  const resolvedReferenceIdx = referenceIdx !== -1 ? referenceIdx : fallbackReferenceIdx;
  if (queryIdx === -1 && resolvedReferenceIdx === -1 && scoreIdx === -1) return null;

  const referenceCounts = new Map();
  const points = [];
  let minScore = Number.POSITIVE_INFINITY;
  let maxScore = Number.NEGATIVE_INFINITY;
  let missingScoreRows = 0;

  rows.forEach((row, idx) => {
    const rawQuery = queryIdx !== -1 ? normalizeValue(row[queryIdx]) : '';
    const firstSample = rawQuery
      .split(',')
      .map(value => normalizeValue(value))
      .find(value => value);
    const sampleLabel = firstSample || `Sample ${idx + 1}`;
    const reference = formatReferenceLabel(resolvedReferenceIdx !== -1 ? row[resolvedReferenceIdx] : '');
    let score = parseNumericValue(scoreIdx !== -1 ? row[scoreIdx] : null);
    if (Number.isFinite(score)) {
      if (score > 1.2 && score <= 100) score = score / 100;
      score = Math.max(0, Math.min(1, score));
    } else {
      score = 0;
      missingScoreRows += 1;
    }

    minScore = Math.min(minScore, score);
    maxScore = Math.max(maxScore, score);
    referenceCounts.set(reference, (referenceCounts.get(reference) || 0) + 1);
    points.push({
      id: `match-${idx}`,
      x: idx + 1,
      query: rawQuery || sampleLabel,
      sampleLabel,
      reference,
      score
    });
  });

  if (points.length === 0) return null;
  const references = [...referenceCounts.entries()]
    .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0], undefined, { sensitivity: 'base' }));
  const bestPoint = points.reduce((best, point) => (point.score > best.score ? point : best), points[0]);
  const meanScore = points.reduce((acc, point) => acc + point.score, 0) / points.length;
  const referenceNames = [...referenceCounts.keys()];
  const referenceColors = {};
  generateChartColors(referenceNames).forEach((color, index) => {
    referenceColors[referenceNames[index]] = color;
  });

  return {
    heading: 'Reference Match Timeline',
    caption: 'Interactive best-reference match score by input row.',
    kpis: [
      {
        label: 'Input Rows',
        value: formatIntegerForViz(points.length),
        sub: `${missingScoreRows > 0 ? `${missingScoreRows} rows with fallback score` : 'all rows with score'}`
      },
      {
        label: 'Best Reference',
        value: bestPoint.reference,
        sub: `${formatPlotDecimal(bestPoint.score * 100, 2)}% best match`
      },
      {
        label: 'Mean Match Score',
        value: `${formatPlotDecimal(meanScore * 100, 2)}%`,
        sub: `${formatPlotDecimal((minScore || 0) * 100, 2)}-${formatPlotDecimal((maxScore || 0) * 100, 2)}%`
      },
      {
        label: 'Distinct References',
        value: formatIntegerForViz(references.length),
        sub: `top: ${references[0]?.[0] || 'n/a'}`
      }
    ],
    topReferences: references.slice(0, 7),
    trackItems: points,
    trackDomainStart: 1,
    trackDomainEnd: Math.max(2, points.length),
    trackReferenceColors: referenceColors,
    referenceOptions: referenceNames,
    referenceCounts,
    scoreMin: Number.isFinite(minScore) ? minScore : 0,
    scoreMax: Number.isFinite(maxScore) ? maxScore : 1
  };
}

function buildMatchTrackDetailHtml(model, selectedPoint) {
  if (!selectedPoint) {
    return `
      <div class="viz-match-track-detail-empty">
        Click a point to inspect it. Drag, wheel and zoom to navigate.
      </div>
    `;
  }

  return `
    <div class="viz-match-track-detail-row">
      <span class="viz-match-track-detail-label">Sample / Query</span>
      <strong>${escapeHtml(selectedPoint.query || selectedPoint.sampleLabel || 'Unknown')}</strong>
    </div>
    <div class="viz-match-track-detail-row">
      <span class="viz-match-track-detail-label">Reference</span>
      <strong>${escapeHtml(selectedPoint.reference || 'Unassigned')}</strong>
    </div>
    <div class="viz-match-track-detail-row">
      <span class="viz-match-track-detail-label">Match Score</span>
      <strong>${formatPlotDecimal(selectedPoint.score * 100, 2)}%</strong>
    </div>
    <div class="viz-match-track-detail-row">
      <span class="viz-match-track-detail-label">Input row # (1-based)</span>
      <strong>${formatIntegerForViz(selectedPoint.x)}</strong>
    </div>
  `;
}

function buildMatchTrackSvg(model, state, points, domain) {
  const width = 980;
  const height = 340;
  const marginLeft = 68;
  const marginRight = 18;
  const marginTop = 20;
  const marginBottom = 68;
  const plotWidth = width - marginLeft - marginRight;
  const plotHeight = height - marginTop - marginBottom;
  const plotBottom = marginTop + plotHeight;
  const overviewY = plotBottom + 14;
  const overviewHeight = 10;
  const overviewWindowHeight = 18;

  const visibleItems = points
    .filter(point => point.x >= state.windowStart && point.x <= state.windowEnd)
    .sort((a, b) => a.x - b.x || a.id.localeCompare(b.id));
  const renderedPoints = matchTrackDownsamplePoints(visibleItems, MATCH_TRACK_MAX_POINTS, state.selectedPointId);
  const selectedPoint = matchTrackFindItemById(points, state.selectedPointId);
  const matches = new Set(Array.isArray(state.searchMatchIds) ? state.searchMatchIds : []);
  const windowSpan = Math.max(1, state.windowEnd - state.windowStart);
  const xForWindow = value => marginLeft + ((value - state.windowStart) / windowSpan) * plotWidth;
  const xForDomain = value => marginLeft + ((value - domain.start) / Math.max(1, domain.end - domain.start)) * plotWidth;
  const yForScore = score => marginTop + (1 - score) * plotHeight;
  const yTicks = [0, 0.2, 0.4, 0.6, 0.8, 1];
  const yTickSvg = yTicks.map(value => {
    const y = yForScore(value);
    return `
      <line x1="${marginLeft}" y1="${y}" x2="${marginLeft + plotWidth}" y2="${y}" class="viz-match-track-grid"/>
      <text x="${marginLeft - 8}" y="${y + 4}" text-anchor="end" class="viz-match-track-axis-label">${formatPlotDecimal(value * 100, 0)}%</text>
    `;
  }).join('');

  const xTickCount = 6;
  const xTickSvg = Array.from({ length: xTickCount + 1 }, (_, index) => {
    const value = state.windowStart + (windowSpan * index / xTickCount);
    const x = xForWindow(value);
    return `
      <line x1="${x}" y1="${plotBottom}" x2="${x}" y2="${plotBottom + 6}" class="viz-match-track-axis-tick"/>
      <text x="${x}" y="${plotBottom + 20}" text-anchor="middle" class="viz-match-track-axis-label">${formatIntegerForViz(Math.round(value))}</text>
    `;
  }).join('');

  const pointSvg = renderedPoints.map(point => {
    const x = xForWindow(point.x);
    const y = yForScore(point.score);
    const isSelected = point.id === state.selectedPointId;
    const isMatched = matches.has(point.id);
    const selectedClass = isSelected ? ' is-selected' : '';
    const matchedClass = isMatched ? ' is-match' : '';
    const color = resolveMatchTrackReferenceColor(point, model.trackReferenceColors);
    const hitRadius = isSelected ? 8 : 7;
    return `
      <line x1="${x}" y1="${plotBottom}" x2="${x}" y2="${y}" class="viz-match-track-stem${selectedClass}" stroke="${color}"/>
      <circle cx="${x}" cy="${y}" r="${isSelected ? 4.8 : 3.6}" class="viz-match-point${selectedClass}${matchedClass}" fill="${color}"/>
      <circle
        cx="${x}"
        cy="${y}"
        r="${hitRadius}"
        class="viz-match-point-hit"
        data-match-id="${escapeHtml(point.id)}"
      >
        <title>${escapeSvg(`${point.sampleLabel}: ${point.reference} · ${formatPlotDecimal(point.score * 100, 2)}% · row ${Math.round(point.x)}`)}</title>
      </circle>
    `;
  }).join('');

  const overviewLine = points
    .map(point => {
      const x = xForDomain(point.x);
      const y = overviewY + (overviewHeight / 2) - (point.score * (overviewHeight / 2));
      const color = resolveMatchTrackReferenceColor(point, model.trackReferenceColors);
      return `<line x1="${x}" x2="${x}" y1="${overviewY}" y2="${overviewY + overviewHeight}" class="viz-match-overview-line" stroke="${color}" stroke-opacity="0.5"/>`;
    }).join('');

  const overviewWindowWidth = Math.max(2, (windowSpan / Math.max(1, domain.end - domain.start)) * plotWidth);
  const overviewWindowX = xForDomain(state.windowStart);
  const windowSummary = `Input row window ${formatIntegerForViz(state.windowStart)}-${formatIntegerForViz(state.windowEnd)} / ${formatIntegerForViz(Math.round(domain.start))}-${formatIntegerForViz(Math.round(domain.end))} (1-based)`;

  const svg = `
    <svg
      class="viz-match-track-svg"
      viewBox="0 0 ${width} ${height}"
      preserveAspectRatio="xMidYMid meet"
      style="aspect-ratio: ${width} / ${height};"
      role="img"
      aria-label="Interactive match reference timeline"
      data-view-width="${width}"
      data-view-height="${height}"
      data-plot-left="${marginLeft}"
      data-plot-width="${plotWidth}"
    >
      <rect x="${marginLeft}" y="${marginTop}" width="${plotWidth}" height="${plotHeight}" class="viz-match-track-plot-bg"/>
      ${yTickSvg}
      <line x1="${marginLeft}" y1="${plotBottom}" x2="${marginLeft + plotWidth}" y2="${plotBottom}" class="viz-match-track-axis"/>
      ${xTickSvg}
      <rect
        x="${overviewWindowX}"
        y="${overviewY - 1}"
        width="${overviewWindowWidth}"
        height="${overviewHeight + 2}"
        class="viz-match-overview-window"
      />
      <rect
        x="${marginLeft}"
        y="${overviewY - 4}"
        width="${plotWidth}"
        height="${overviewHeight + 8}"
        class="viz-match-overview-bg"
      />
      ${overviewLine}
      <rect
        x="${marginLeft}"
        y="${overviewY - 5}"
        width="${plotWidth}"
        height="${overviewHeight + 10}"
        class="viz-match-overview-hit"
        data-match-action="overview-hit"
      />
      ${pointSvg}
      <text x="${marginLeft + plotWidth / 2}" y="${height - 7}" text-anchor="middle" class="viz-match-track-axis-title">Input row # (1-based)</text>
      <text
        x="15"
        y="${marginTop + plotHeight / 2}"
        text-anchor="middle"
        class="viz-match-track-axis-title"
        transform="rotate(-90 15 ${marginTop + plotHeight / 2})"
      >
        Match score
      </text>
    </svg>
  `;

  return {
    svg,
    visibleCount: visibleItems.length,
    renderedCount: renderedPoints.length,
    selectedPoint,
    overviewWindowLabel: windowSummary
  };
}

function matchTrackBuildModelKey(model) {
  const firstId = model?.trackItems?.[0]?.id || '';
  const lastId = model?.trackItems?.[model.trackItems.length - 1]?.id || '';
  return [
    'match-v1',
    model?.trackItems?.length || 0,
    Math.round(model?.trackDomainStart || 0),
    Math.round(model?.trackDomainEnd || 0),
    firstId,
    lastId
  ].join('|');
}

function attachMatchTrackInteractions(toolId, data, primaryColumn, model, state, points, domain, section) {
  if (!section) return;

  const referenceSelect = section.querySelector('.viz-match-track-reference-select');
  const searchInput = section.querySelector('.viz-match-track-search-input');
  const startInput = section.querySelector('.viz-match-track-start-input');
  const endInput = section.querySelector('.viz-match-track-end-input');
  const svgEl = section.querySelector('.viz-match-track-svg');
  const applyRangeBtn = section.querySelector('[data-match-action="apply-range"]');

  referenceSelect?.addEventListener('change', () => {
    state.referenceFilter = referenceSelect.value || '__all__';
    state.activeReferenceFilter = '';
    state.selectedPointId = '';
    state.searchMatchIndex = -1;
    scheduleMatchTrackRender(toolId, data, primaryColumn);
  });

  searchInput?.addEventListener('input', () => {
    state.searchQuery = searchInput.value || '';
    state.searchMatchIndex = 0;
    scheduleMatchTrackRender(toolId, data, primaryColumn);
  });
  searchInput?.addEventListener('keydown', event => {
    if (event.key !== 'Enter') return;
    event.preventDefault();
    section.querySelector('[data-match-action="find"]')?.click();
  });
  [startInput, endInput].forEach(input => {
    input?.addEventListener('keydown', event => {
      if (event.key !== 'Enter') return;
      event.preventDefault();
      applyRangeBtn?.click();
    });
  });

  const moveToMatch = (direction = 1, fromCurrent = true) => {
    if (!Array.isArray(state.searchMatchIds) || state.searchMatchIds.length === 0) return;
    let index = Math.max(0, Number(state.searchMatchIndex || 0));
    if (fromCurrent) {
      index = (index + direction + state.searchMatchIds.length) % state.searchMatchIds.length;
    }
    state.searchMatchIndex = index;
    const matchId = state.searchMatchIds[index];
    const matchPoint = matchTrackFindItemById(points, matchId);
    if (matchPoint) matchTrackFocusItem(state, matchPoint, domain.start, domain.end, true);
  };

  section.querySelectorAll('[data-match-action]').forEach(button => {
    button.addEventListener('click', () => {
      const action = button.dataset.matchAction;
      if (!action) return;

      switch (action) {
        case 'zoom-in':
          matchTrackZoom(state, 0.7, null, domain.start, domain.end);
          break;
        case 'zoom-out':
          matchTrackZoom(state, 1.35, null, domain.start, domain.end);
          break;
        case 'pan-left':
          matchTrackPan(state, -0.2, domain.start, domain.end);
          break;
        case 'pan-right':
          matchTrackPan(state, 0.2, domain.start, domain.end);
          break;
        case 'reset': {
          const span = Math.max(2, domain.end - domain.start);
          matchTrackSetWindow(state, domain.start, domain.start + span, domain.start, domain.end);
          state.selectedPointId = '';
          break;
        }
        case 'apply-range': {
          const start = Number.parseFloat(startInput?.value || '');
          const end = Number.parseFloat(endInput?.value || '');
          if (Number.isFinite(start) && Number.isFinite(end)) {
            matchTrackSetWindow(state, start, end, domain.start, domain.end);
          }
          break;
        }
        case 'find':
          if (state.searchMatchIds.length === 0 && String(state.searchQuery || '').trim()) {
            state.searchMatchIds = matchTrackFindMatches(points, state.searchQuery).map(point => point.id);
            state.searchMatchIndex = 0;
          }
          moveToMatch(0, false);
          break;
        case 'search-next':
          moveToMatch(1, true);
          break;
        case 'search-prev':
          moveToMatch(-1, true);
          break;
        default:
          return;
      }

      scheduleMatchTrackRender(toolId, data, primaryColumn);
    });
  });

  if (!svgEl) return;

  svgEl.addEventListener('click', event => {
    const pointHit = event.target.closest('.viz-match-point-hit');
    if (pointHit) {
      const pointId = pointHit.getAttribute('data-match-id') || '';
      const matchPoint = matchTrackFindItemById(points, pointId);
      if (!matchPoint) return;
      matchTrackFocusItem(state, matchPoint, domain.start, domain.end, true);
      scheduleMatchTrackRender(toolId, data, primaryColumn);
      return;
    }
    if (event.target.closest('.viz-match-overview-hit')) {
      const center = matchTrackDomainValueFromEvent(event, svgEl, state, true);
      const span = Math.max(1, state.windowEnd - state.windowStart);
      matchTrackSetWindow(state, center - span / 2, center + span / 2, domain.start, domain.end);
      scheduleMatchTrackRender(toolId, data, primaryColumn);
    }
  });

  svgEl.addEventListener('wheel', event => {
    event.preventDefault();
    const center = matchTrackDomainValueFromEvent(event, svgEl, state);
    matchTrackZoom(state, event.deltaY < 0 ? 0.78 : 1.25, center, domain.start, domain.end);
    scheduleMatchTrackRender(toolId, data, primaryColumn);
  }, { passive: false });

  svgEl.addEventListener('mousedown', event => {
    if (event.button !== 0) return;
    if (event.target.closest('.viz-match-point-hit')) return;
    if (event.target.closest('.viz-match-overview-hit')) return;

    const viewWidth = Number(svgEl.dataset.viewWidth || 980);
    const plotWidth = Math.max(1, Number(svgEl.dataset.plotWidth || 1));
    const startViewX = getSvgViewXFromClient(svgEl, event.clientX, viewWidth / 2);
    const drag = {
      startViewX,
      startWindowStart: state.windowStart,
      startWindowEnd: state.windowEnd,
      plotWidth,
      spanAtStart: Math.max(1, state.windowEnd - state.windowStart)
    };

    document.body.classList.add('track-dragging');
    event.preventDefault();

    const onMove = moveEvent => {
      const currentViewX = getSvgViewXFromClient(svgEl, moveEvent.clientX, drag.startViewX);
      const deltaView = currentViewX - drag.startViewX;
      const deltaDomain = (deltaView / drag.plotWidth) * drag.spanAtStart;
      matchTrackSetWindow(
        state,
        drag.startWindowStart - deltaDomain,
        drag.startWindowEnd - deltaDomain,
        domain.start,
        domain.end
      );
      scheduleMatchTrackRender(toolId, data, primaryColumn);
    };

    const onUp = () => {
      document.body.classList.remove('track-dragging');
      window.removeEventListener('mousemove', onMove);
      window.removeEventListener('mouseup', onUp);
    };

    window.addEventListener('mousemove', onMove);
    window.addEventListener('mouseup', onUp);
  });
}

function renderMatchRefTrackInsights(toolId, data, primaryColumn, breakdownEl) {
  if (toolId !== 'match' || !breakdownEl) return;
  const model = buildMatchTrackModel(data, primaryColumn);
  const existingSection = breakdownEl.querySelector('.viz-match-lightweight');
  if (!model) {
    existingSection?.remove();
    return;
  }

  const state = ensureMatchTrackState(toolId, model);
  let points = Array.isArray(model.trackItems) ? [...model.trackItems] : [];
  if (state.referenceFilter && state.referenceFilter !== '__all__') {
    points = points.filter(point => point.reference === state.referenceFilter);
  }
  if (points.length === 0) {
    state.referenceFilter = '__all__';
    points = Array.isArray(model.trackItems) ? [...model.trackItems] : [];
  }
  points.sort((a, b) => a.x - b.x || a.id.localeCompare(b.id));

  const domain = matchTrackGetDomain(points);
  const domainChanged = state.domainStart !== domain.start || state.domainEnd !== domain.end;
  const filterChanged = state.activeReferenceFilter !== state.referenceFilter;
  if (domainChanged || filterChanged) {
    matchTrackSetWindow(state, domain.start, domain.end, domain.start, domain.end);
    if (filterChanged) state.selectedPointId = '';
  } else {
    matchTrackSetWindow(state, state.windowStart, state.windowEnd, domain.start, domain.end);
  }
  state.domainStart = domain.start;
  state.domainEnd = domain.end;
  state.activeReferenceFilter = state.referenceFilter;

  const query = String(state.searchQuery || '').trim();
  if (!query) {
    state.searchMatchIds = [];
    state.searchMatchIndex = -1;
  } else {
    state.searchMatchIds = matchTrackFindMatches(points, query).map(point => point.id);
    if (state.searchMatchIds.length === 0) {
      state.searchMatchIndex = -1;
    } else if (state.searchMatchIndex < 0 || state.searchMatchIndex >= state.searchMatchIds.length) {
      state.searchMatchIndex = 0;
    }
    const selectedIndex = state.searchMatchIds.indexOf(state.selectedPointId);
    if (selectedIndex !== -1) state.searchMatchIndex = selectedIndex;
  }

  if (state.selectedPointId && !matchTrackFindItemById(points, state.selectedPointId)) {
    state.selectedPointId = '';
  }

  const trackInfo = buildMatchTrackSvg(model, state, points, domain);
  const referenceOptions = Array.isArray(model.referenceOptions) ? model.referenceOptions : [];
  const searchStatus = matchTrackSearchStatus(state);
  const searchOptions = matchTrackBuildSearchOptions(points);
  const windowStartValue = Math.round(state.windowStart);
  const windowEndValue = Math.round(state.windowEnd);
  const selectedDetail = buildMatchTrackDetailHtml(model, trackInfo.selectedPoint);
  const matchWindowLabel = 'Window start/end (input row #)';

  const section = existingSection || document.createElement('section');
  section.className = 'viz-match-lightweight';
  section.innerHTML = `
    <div class="viz-match-head">
      <h5>${escapeHtml(model.heading || 'Match Timeline')}</h5>
      <span>${escapeHtml(model.caption || '')}</span>
    </div>
    <div class="viz-match-kpis">
      ${(model.kpis || []).map(kpi => `
        <article class="viz-match-kpi">
          <span class="viz-match-kpi-label">${escapeHtml(kpi.label || '')}</span>
          <strong class="viz-match-kpi-value">${escapeHtml(String(kpi.value || '0'))}</strong>
          <span class="viz-match-kpi-sub">${escapeHtml(kpi.sub || '')}</span>
        </article>
      `).join('')}
    </div>
    ${model.topReferences.length > 0 ? `
      <div class="viz-match-reference-chips">
        ${model.topReferences.map(([reference, value]) => {
          const share = points.length > 0 ? (value / model.trackItems.length) * 100 : 0;
          return `
            <span class="viz-match-reference-chip" title="${escapeHtml(reference)}">
              <span>${escapeHtml(reference)}</span>
              <strong>${formatIntegerForViz(value)}</strong>
              <em>${formatPlotDecimal(share, 1)}%</em>
            </span>
          `;
        }).join('')}
      </div>
    ` : ''}
    <div class="viz-match-track-shell">
      <div class="viz-match-track-controls">
        <div class="viz-match-track-toolbar-top">
          <label class="viz-match-track-field">
            <span class="viz-match-track-field-label">Reference</span>
            <select class="viz-match-track-reference-select">
              <option value="__all__" ${state.referenceFilter === '__all__' ? 'selected' : ''}>All references</option>
              ${referenceOptions.map(reference => `
                <option value="${escapeHtml(reference)}" ${reference === state.referenceFilter ? 'selected' : ''}>
                  ${escapeHtml(reference)}
                </option>
              `).join('')}
            </select>
          </label>
          <div class="viz-match-track-field viz-match-track-search">
            <span class="viz-match-track-field-label">Find match</span>
            <div class="viz-match-track-search-wrap">
              <input
                type="text"
                class="viz-match-track-search-input"
                placeholder="Sample, reference, input row, score..."
                value="${escapeHtml(state.searchQuery || '')}"
                list="${toolId}-viz-match-search-options"
                autocomplete="off"
              />
              <datalist id="${toolId}-viz-match-search-options">
                ${searchOptions}
              </datalist>
              <div class="viz-match-track-search-actions">
                <button type="button" class="viz-track-btn viz-track-btn-search" data-match-action="find">
                  <span class="viz-track-btn-icon" aria-hidden="true">⌕</span>
                  <span class="viz-track-btn-label">Find</span>
                </button>
                <button type="button" class="viz-track-btn viz-track-btn-search-nav" data-match-action="search-prev" ${state.searchMatchIds.length > 0 ? '' : 'disabled'}>
                  <span class="viz-track-btn-icon" aria-hidden="true">‹</span>
                  <span class="viz-track-btn-label">Prev</span>
                </button>
                <button type="button" class="viz-track-btn viz-track-btn-search-nav" data-match-action="search-next" ${state.searchMatchIds.length > 0 ? '' : 'disabled'}>
                  <span class="viz-track-btn-icon" aria-hidden="true">›</span>
                  <span class="viz-track-btn-label">Next</span>
                </button>
              </div>
            </div>
            <span class="viz-match-track-search-status" aria-live="polite">${escapeHtml(searchStatus)}</span>
          </div>
        </div>
        <div class="viz-match-track-toolbar-bottom">
          <label class="viz-match-track-field viz-match-track-window">
            <span class="viz-match-track-field-label">${matchWindowLabel}</span>
            <div class="viz-match-track-window-wrap">
              <input
                type="number"
                class="viz-match-track-start-input"
                value="${windowStartValue}"
                step="1"
                min="${Math.floor(domain.start)}"
                max="${Math.ceil(domain.end)}"
              />
              <span class="viz-match-track-window-sep">to</span>
              <input
                type="number"
                class="viz-match-track-end-input"
                value="${windowEndValue}"
                step="1"
                min="${Math.floor(domain.start)}"
                max="${Math.ceil(domain.end)}"
              />
              <button type="button" class="viz-track-btn viz-track-btn-primary" data-match-action="apply-range">
                <span class="viz-track-btn-icon" aria-hidden="true">✓</span>
                <span class="viz-track-btn-label">Apply</span>
              </button>
            </div>
          </label>
          <div class="viz-match-track-nav" role="group" aria-label="Match timeline navigation">
            <div class="viz-track-btn-group">
              <button type="button" class="viz-track-btn" data-match-action="zoom-in">
                <span class="viz-track-btn-icon" aria-hidden="true">＋</span>
                <span class="viz-track-btn-label">Zoom In</span>
              </button>
              <button type="button" class="viz-track-btn" data-match-action="zoom-out">
                <span class="viz-track-btn-icon" aria-hidden="true">－</span>
                <span class="viz-track-btn-label">Zoom Out</span>
              </button>
              <button type="button" class="viz-track-btn" data-match-action="pan-left">
                <span class="viz-track-btn-icon" aria-hidden="true">←</span>
                <span class="viz-track-btn-label">Pan Left</span>
              </button>
              <button type="button" class="viz-track-btn" data-match-action="pan-right">
                <span class="viz-track-btn-icon" aria-hidden="true">→</span>
                <span class="viz-track-btn-label">Pan Right</span>
              </button>
            </div>
            <button type="button" class="viz-track-btn" data-match-action="reset">
              <span class="viz-track-btn-icon" aria-hidden="true">↺</span>
              <span class="viz-track-btn-label">Reset</span>
            </button>
          </div>
        </div>
      </div>
      <div class="viz-match-track-meta">
        <span>${trackInfo.visibleCount} point(s) in window · ${points.length} filtered</span>
        <span>${trackInfo.renderedCount < trackInfo.visibleCount ? `rendered ${trackInfo.renderedCount} for performance` : 'full resolution render'}</span>
        <span>X-axis: input row number from loaded results (1-based).</span>
        <span>Row #1 = first loaded row; row numbers stay fixed even if table is sorted or filtered.</span>
        <span>${trackInfo.overviewWindowLabel || ''}</span>
        <span>Controls: drag to pan, wheel to zoom, click point to focus</span>
      </div>
      <div class="viz-match-track-canvas">
        ${trackInfo.svg}
      </div>
      <div class="viz-match-track-detail">
        ${selectedDetail}
      </div>
    </div>
  `;

  if (!existingSection) {
    breakdownEl.appendChild(section);
  }

  attachMatchTrackInteractions(toolId, data, primaryColumn, model, state, points, domain, section);
}

function renderToolInsights(toolId, data, counts, primaryColumn) {
  const gridEl = document.getElementById(`${toolId}-viz-insights-grid`);
  const breakdownEl = document.getElementById(`${toolId}-viz-outcome-breakdown`);
  if (!gridEl || !breakdownEl) return;

  const headers = Array.isArray(data?.headers) ? data.headers : [];
  const rows = Array.isArray(data?.rows) ? data.rows : [];
  if (rows.length === 0) {
    gridEl.innerHTML = '';
    breakdownEl.innerHTML = '<p class="viz-outcome-empty">No rows available for insight analysis.</p>';
    return;
  }

  const labelColIdx = resolvePrimaryColumnIndex(headers, primaryColumn);
  const sampleCount = estimateSampleCount(data, labelColIdx);
  // On a resistance panel one sample contributes a count per drug, so summing
  // the tally would exceed the number of samples ("600% assigned"). What is
  // assigned there is the sample itself: it carries at least one marker.
  const drCalls = buildDrCallSummary(data, primaryColumn);
  const assignedCount = drCalls
    ? drCalls.sampleSummaries.size
    : Object.values(counts).reduce((acc, value) => acc + value, 0);
  const unknownCount = Math.max(sampleCount - assignedCount, 0);
  const distinctCount = Object.keys(counts).length;
  const topEntry = Object.entries(counts)
    .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0], undefined, { sensitivity: 'base' }))[0];
  const topLabel = topEntry?.[0] || 'N/A';
  const topCount = topEntry?.[1] || 0;
  const topShare = sampleCount > 0 ? (topCount / sampleCount) * 100 : 0;

  const quality = findQualityMetricSummary(toolId, data);
  const cards = [
    {
      label: 'Samples',
      value: String(sampleCount),
      sub: `${rows.length} result rows`
    },
    {
      label: 'Assigned',
      value: `${assignedCount}`,
      sub: `${formatPlotDecimal(sampleCount > 0 ? (assignedCount / sampleCount) * 100 : 0, 1)}% assigned`
    },
    {
      label: 'Distinct Calls',
      value: String(distinctCount),
      sub: unknownCount > 0 ? `${unknownCount} unresolved` : 'No unresolved samples'
    },
    {
      label: 'Top Call',
      value: topLabel,
      sub: `${topCount} samples (${formatPlotDecimal(topShare, 1)}%)`
    }
  ];

  if (quality) {
    let qualityClass = 'quality-ok';
    if (quality.danger != null && quality.median < quality.danger) {
      qualityClass = 'quality-danger';
    } else if (quality.warn != null && quality.median < quality.warn) {
      qualityClass = 'quality-warn';
    }
    cards.push({
      label: quality.label,
      value: formatPlotDecimal(quality.median, 2),
      sub: `P90 ${formatPlotDecimal(quality.p90, 2)} · n=${quality.count}`,
      extraClass: qualityClass
    });
  }

  gridEl.innerHTML = cards.map(card => `
    <article class="viz-kpi-card${card.extraClass ? ' ' + card.extraClass : ''}">
      <span class="viz-kpi-label">${escapeHtml(card.label)}</span>
      <strong class="viz-kpi-value">${escapeHtml(String(card.value))}</strong>
      <span class="viz-kpi-sub">${escapeHtml(card.sub)}</span>
    </article>
  `).join('');

  if (sampleCount === 1) {
    breakdownEl.innerHTML = `
      <div class="viz-outcome-header">
        <h4>Single-sample outcome</h4>
        <span>Distribution charts hidden for clarity</span>
      </div>
      <div class="viz-outcome-list">
        <div class="viz-outcome-item">
          <div class="viz-outcome-meta">
            <span class="viz-outcome-label">${escapeHtml(topLabel)}</span>
            <span class="viz-outcome-value">${topCount} sample</span>
          </div>
          <div class="viz-outcome-track">
            <span class="viz-outcome-fill" style="width:100%"></span>
          </div>
        </div>
      </div>
    `;
    renderMatchRefTrackInsights(toolId, data, primaryColumn, breakdownEl);
    renderSplitFastqLightweightInsights(toolId, data, primaryColumn, breakdownEl);
    renderClassifyGenotypingInsights(toolId, data, breakdownEl);
    return;
  }

  const sortedEntries = Object.entries(counts)
    .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0], undefined, { sensitivity: 'base' }));
  const maxValue = Math.max(1, ...sortedEntries.map(([, value]) => value), unknownCount);
  const topEntries = sortedEntries.slice(0, 6);
  if (unknownCount > 0) topEntries.push(['Unresolved', unknownCount]);

  if (topEntries.length === 0) {
    breakdownEl.innerHTML = '<p class="viz-outcome-empty">No categorical outcomes to display.</p>';
  } else {
    breakdownEl.innerHTML = `
      <div class="viz-outcome-header">
        <h4>Outcome Breakdown</h4>
        <span>Top ${Math.min(topEntries.length, 6)} + unresolved</span>
      </div>
      <div class="viz-outcome-list">
        ${topEntries.map(([label, value]) => {
          const widthPct = (value / maxValue) * 100;
          const sharePct = sampleCount > 0 ? (value / sampleCount) * 100 : 0;
          return `
            <div class="viz-outcome-item">
              <div class="viz-outcome-meta">
                <span class="viz-outcome-label">${escapeHtml(label)}</span>
                <span class="viz-outcome-value">${value} (${formatPlotDecimal(sharePct, 1)}%)</span>
              </div>
              <div class="viz-outcome-track">
                <span class="viz-outcome-fill" style="width:${formatPlotDecimal(widthPct, 2)}%"></span>
              </div>
            </div>
          `;
        }).join('')}
      </div>
    `;
  }

  renderMatchRefTrackInsights(toolId, data, primaryColumn, breakdownEl);
  renderSplitFastqLightweightInsights(toolId, data, primaryColumn, breakdownEl);
  renderClassifyGenotypingInsights(toolId, data, breakdownEl);
  renderPredictConfidenceInsights(toolId, data, breakdownEl);
}

function clearToolInsights(toolId) {
  const gridEl = document.getElementById(`${toolId}-viz-insights-grid`);
  const breakdownEl = document.getElementById(`${toolId}-viz-outcome-breakdown`);
  if (gridEl) gridEl.replaceChildren();
  if (breakdownEl) breakdownEl.replaceChildren();
  if (SPLITFQ_TRACK_RENDER_THROTTLES[toolId]) {
    if (typeof window !== 'undefined' && typeof window.cancelAnimationFrame === 'function') {
      window.cancelAnimationFrame(SPLITFQ_TRACK_RENDER_THROTTLES[toolId]);
    } else {
      clearTimeout(SPLITFQ_TRACK_RENDER_THROTTLES[toolId]);
    }
    delete SPLITFQ_TRACK_RENDER_THROTTLES[toolId];
  }
  delete SPLITFQ_TRACK_STATE[toolId];
  if (MATCH_TRACK_RENDER_THROTTLES[toolId]) {
    if (typeof window !== 'undefined' && typeof window.cancelAnimationFrame === 'function') {
      window.cancelAnimationFrame(MATCH_TRACK_RENDER_THROTTLES[toolId]);
    } else {
      clearTimeout(MATCH_TRACK_RENDER_THROTTLES[toolId]);
    }
    delete MATCH_TRACK_RENDER_THROTTLES[toolId];
  }
  delete MATCH_TRACK_STATE[toolId];
}

function clearGenomicTrack(toolId) {
  setTrackExpanded(toolId, false);
  const section = document.getElementById(`${toolId}-viz-genomic-track`);
  const bodyEl = document.getElementById(`${toolId}-viz-track-body`);
  const controlsEl = document.getElementById(`${toolId}-viz-track-controls`);
  if (bodyEl) bodyEl.replaceChildren();
  if (controlsEl) controlsEl.replaceChildren();
  if (section) section.classList.add('hidden');
  Object.keys(TRACK_FETCH_TIMERS).forEach(key => {
    if (key === toolId || key.startsWith(`${toolId}:`)) {
      clearTimeout(TRACK_FETCH_TIMERS[key]);
      delete TRACK_FETCH_TIMERS[key];
    }
  });
  delete TRACK_STATE[toolId];
  clearTrackRenderCache(toolId);
}

function getTrackToolIdFromSection(section) {
  const sectionId = String(section?.id || '');
  if (!sectionId.endsWith('-viz-genomic-track')) return '';
  return sectionId.slice(0, -'-viz-genomic-track'.length);
}

function ensureTrackExpandAnchor(toolId, section) {
  const existingAnchor = TRACK_EXPANDED_ANCHORS[toolId];
  if (existingAnchor?.isConnected) return existingAnchor;
  if (!section?.parentNode) return null;

  const anchor = document.createElement('span');
  anchor.className = 'viz-track-expand-anchor';
  anchor.setAttribute('aria-hidden', 'true');
  section.parentNode.insertBefore(anchor, section);
  TRACK_EXPANDED_ANCHORS[toolId] = anchor;
  return anchor;
}

function moveTrackSectionToBody(toolId, section) {
  if (!section?.isConnected || section.parentElement === document.body) return;
  const anchor = ensureTrackExpandAnchor(toolId, section);
  if (!anchor?.isConnected) return;
  document.body.appendChild(section);
}

function restoreTrackSectionFromBody(toolId, section) {
  const anchor = TRACK_EXPANDED_ANCHORS[toolId];
  if (anchor?.parentNode) {
    anchor.parentNode.insertBefore(section, anchor);
    anchor.remove();
  }
  delete TRACK_EXPANDED_ANCHORS[toolId];
}

function collapseExpandedTrackSection(section) {
  if (!section) return;
  const otherToolId = getTrackToolIdFromSection(section);
  section.classList.remove('expanded');
  if (otherToolId) {
    restoreTrackSectionFromBody(otherToolId, section);
  }
}

function ensureTrackExpandEscapeHandler() {
  if (trackExpandEscapeBound) return;
  document.addEventListener('keydown', event => {
    if (event.key !== 'Escape') return;
    const expanded = document.querySelectorAll('.viz-genomic-track.expanded');
    if (expanded.length === 0) return;
    expanded.forEach(section => collapseExpandedTrackSection(section));
    syncTrackExpandedBodyState();
  });
  trackExpandEscapeBound = true;
}

function isTrackExpanded(toolId) {
  const section = document.getElementById(`${toolId}-viz-genomic-track`);
  return Boolean(section?.classList.contains('expanded'));
}

function syncTrackExpandedBodyState() {
  const hasExpandedTrack = document.querySelector('.viz-genomic-track.expanded');
  document.body.classList.toggle('track-expanded-open', Boolean(hasExpandedTrack));
}

function setTrackExpanded(toolId, expanded) {
  ensureTrackExpandEscapeHandler();
  const section = document.getElementById(`${toolId}-viz-genomic-track`);
  if (!section) {
    syncTrackExpandedBodyState();
    return;
  }
  const shouldExpand = Boolean(expanded);
  if (shouldExpand) {
    document.querySelectorAll('.viz-genomic-track.expanded').forEach(node => {
      if (node !== section) collapseExpandedTrackSection(node);
    });
    moveTrackSectionToBody(toolId, section);
  } else {
    restoreTrackSectionFromBody(toolId, section);
  }
  section.classList.toggle('expanded', shouldExpand);
  syncTrackExpandedBodyState();
}

function escapeSvg(value) {
  return String(value || '')
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
}

function normalizeNtBase(base) {
  const upper = String(base || '').trim().toUpperCase();
  return /^[ACGTN]$/.test(upper) ? upper : 'N';
}

function getNtBaseColor(base) {
  return NT_BASE_COLORS[normalizeNtBase(base)] || NT_BASE_COLORS.N;
}

function findHeaderIndexByTokens(headers, tokens, options = {}) {
  const { exactOnly = false } = options;
  const normalizedHeaders = headers.map(normalizeHeader);
  for (const token of tokens) {
    const normalizedToken = normalizeHeader(token);
    let index = normalizedHeaders.findIndex(header => header === normalizedToken);
    if (index !== -1) return index;
    if (exactOnly) continue;
    if (normalizedToken.length <= 3) continue;
    index = normalizedHeaders.findIndex(header => header.includes(normalizedToken));
    if (index !== -1) return index;
  }
  return -1;
}

function parseMutationTrackData(data) {
  const headers = Array.isArray(data?.headers) ? data.headers : [];
  const rows = Array.isArray(data?.rows) ? data.rows : [];
  if (headers.length === 0 || rows.length === 0) {
    return { records: [], samples: [], reason: 'No tabular results available.' };
  }
  const normalizedHeaders = headers.map(header => normalizeHeader(header));

  const refPosIdx = (() => {
    const prioritized = findHeaderIndexByTokens(headers, [
      'SNPreference',
      'reference_pos',
      'reference_position',
      'ref_pos',
      'refposition'
    ]);
    if (prioritized !== -1) return prioritized;
    return findHeaderIndexByTokens(headers, ['pos', 'position'], { exactOnly: true });
  })();

  const genomePosIdx = (() => {
    const prioritized = findHeaderIndexByTokens(headers, [
      'SNPgenome',
      'genome_pos',
      'sample_pos',
      'query_pos',
      'k-merpos',
      'k_merpos',
      'kmer_pos',
      'kmerpos'
    ]);
    if (prioritized !== -1) return prioritized;
    return -1;
  })();

  if (refPosIdx === -1 && genomePosIdx === -1) {
    return {
      records: [],
      samples: [],
      reason: 'No mutation position columns were found in this result file.'
    };
  }

  const genomePosHeader = genomePosIdx >= 0 ? normalizedHeaders[genomePosIdx] : '';
  const genomePosUsesKmerPos = (
    genomePosHeader === 'k_merpos' ||
    genomePosHeader === 'kmer_pos' ||
    genomePosHeader === 'kmerpos'
  );

  const sampleIdx = findHeaderIndexByTokens(headers, [
    'genome',
    'sample',
    'sample_name',
    'isolate',
    'file'
  ]);
  const geneIdx = findHeaderIndexByTokens(headers, ['gene', 'gene_id', 'locus_tag']);
  const geneStartIdx = findHeaderIndexByTokens(headers, [
    'gene_start',
    'gff_start',
    'cds_start'
  ]);
  const geneEndIdx = findHeaderIndexByTokens(headers, [
    'gene_end',
    'gff_end',
    'cds_end'
  ]);
  const lineageIdx = findHeaderIndexByTokens(headers, [
    'lineage_path',
    'lineage',
    'major_lineage',
    'prediction',
    'best_match'
  ]);
  const refAlleleIdx = findHeaderIndexByTokens(headers, ['ref_allele', 'ref', 'reference_allele']);
  const altAlleleIdx = findHeaderIndexByTokens(headers, ['alt_allele', 'alt', 'alternate_allele']);
  const kmerIdx = findHeaderIndexByTokens(headers, ['k-mer', 'kmer', 'marker_kmer']);
  const aaPosIdx = findHeaderIndexByTokens(headers, ['aa_pos', 'amino_acid_pos']);
  const aaChangeIdx = findHeaderIndexByTokens(headers, ['aa_change', 'amino_acid_change']);

  const records = [];
  for (let i = 0; i < rows.length; i++) {
    const row = rows[i];
    const sourceRowIndex = normalizeTrackSourceRowIndex(row?.__sourceRowIndex ?? i);
    const refPosRaw = refPosIdx >= 0 ? parseNumericValue(row[refPosIdx]) : null;
    const genomePosRaw = genomePosIdx >= 0 ? parseNumericValue(row[genomePosIdx]) : null;
    if (refPosRaw === null && genomePosRaw === null) continue;
    const resolvedRefPosRaw = refPosRaw !== null ? refPosRaw : genomePosRaw;
    const resolvedGenomePosRaw = genomePosRaw !== null ? genomePosRaw : refPosRaw;
    if (!Number.isFinite(resolvedRefPosRaw) || !Number.isFinite(resolvedGenomePosRaw)) continue;
    const refPos = Math.max(1, Math.round(resolvedRefPosRaw));

    const sampleRaw = sampleIdx >= 0 ? normalizeValue(row[sampleIdx]) : '';
    const sample = sampleRaw || 'All samples';
    const gene = geneIdx >= 0 ? normalizeValue(row[geneIdx]) : '';
    const geneStartRaw = geneStartIdx >= 0 ? parseNumericValue(row[geneStartIdx]) : null;
    const geneEndRaw = geneEndIdx >= 0 ? parseNumericValue(row[geneEndIdx]) : null;
    const geneStart = geneStartRaw === null ? null : Math.max(1, Math.round(geneStartRaw));
    const geneEnd = geneEndRaw === null ? null : Math.max(1, Math.round(geneEndRaw));
    const lineage = lineageIdx >= 0 ? normalizeValue(row[lineageIdx]) : '';
    let refAllele = refAlleleIdx >= 0 ? normalizeValue(row[refAlleleIdx]) : '';
    let altAllele = altAlleleIdx >= 0 ? normalizeValue(row[altAlleleIdx]) : '';
    let altAlleleSource = altAllele ? 'column' : 'unknown';
    const kmer = kmerIdx >= 0 ? normalizeValue(row[kmerIdx]) : '';

    // SNPgenome and k-merPOS already arrive 1-based from classify
    // (format_marker_match emits genome_position + 1 / variant_start + 1), so
    // no extra offset is applied here.
    let normalizedGenomePosRaw = resolvedGenomePosRaw;
    if (genomePosRaw !== null && genomePosUsesKmerPos && kmer) {
      const center = Math.floor(kmer.length / 2);
      if (Number.isFinite(center) && center > 0) {
        normalizedGenomePosRaw += center;
      }
    }

    if (!altAllele && kmer) {
      const center = Math.floor(kmer.length / 2);
      const base = (kmer[center] || '').toUpperCase();
      if (/^[ACGTN]$/.test(base)) {
        altAllele = base;
        altAlleleSource = 'kmer';
      }
    }
    if (!refAllele && /^[ACGTN]$/.test(altAllele)) refAllele = '?';
    const aaPos = aaPosIdx >= 0 ? normalizeValue(row[aaPosIdx]) : '';
    const aaChange = aaChangeIdx >= 0 ? normalizeValue(row[aaChangeIdx]) : '';

    records.push({
      trackKey: `${sample}::${sourceRowIndex}`,
      sourceRowIndex,
      sample,
      refPos,
      genomePos: Math.max(1, Math.round(normalizedGenomePosRaw)),
      gene,
      geneStart,
      geneEnd,
      lineage,
      refAllele,
      altAllele,
      altAlleleSource,
      kmer,
      aaPos,
      aaChange
    });
  }

  if (records.length === 0) {
    return {
      records: [],
      samples: [],
      reason: 'No mutation rows with valid positions were found.'
    };
  }

  records.sort((a, b) => a.refPos - b.refPos || a.genomePos - b.genomePos);
  const samples = Array.from(new Set(records.map(record => record.sample))).sort((a, b) =>
    a.localeCompare(b, undefined, { sensitivity: 'base' })
  );
  return { records, samples, reason: '' };
}

function getCachedMutationTrackData(data) {
  if (!data || typeof data !== 'object') {
    return {
      ...parseMutationTrackData(data),
      recordsBySample: new Map()
    };
  }

  const headersLength = Array.isArray(data.headers) ? data.headers.length : 0;
  const rows = Array.isArray(data.rows) ? data.rows : [];
  const rowsLength = rows.length;
  const firstRowWidth = rowsLength > 0 && Array.isArray(rows[0]) ? rows[0].length : 0;
  const lastRowWidth = rowsLength > 0 && Array.isArray(rows[rowsLength - 1]) ? rows[rowsLength - 1].length : 0;
  const signature = `${headersLength}:${rowsLength}:${firstRowWidth}:${lastRowWidth}`;

  const cached = TRACK_PARSED_DATA_CACHE.get(data);
  if (cached?.signature === signature && cached.parsed) return cached.parsed;

  const parsed = parseMutationTrackData(data);
  const recordsBySample = new Map();
  if (Array.isArray(parsed.records)) {
    parsed.records.forEach(record => {
      const sample = String(record?.sample || '').trim() || 'All samples';
      if (!recordsBySample.has(sample)) {
        recordsBySample.set(sample, []);
      }
      recordsBySample.get(sample).push(record);
    });
  }

  const parsedPayload = {
    ...parsed,
    recordsBySample
  };
  TRACK_PARSED_DATA_CACHE.set(data, {
    signature,
    parsed: parsedPayload
  });
  return parsedPayload;
}

function summarizeTrackFeatureCollection(features) {
  if (!Array.isArray(features) || features.length === 0) return '0';
  const first = features[0] || {};
  const last = features[features.length - 1] || {};
  return [
    features.length,
    Number(first.start) || 0,
    Number(first.end) || 0,
    Number(last.start) || 0,
    Number(last.end) || 0
  ].join(':');
}

function buildTrackSampleRecordSignature(sampleRecords) {
  if (!Array.isArray(sampleRecords) || sampleRecords.length === 0) return '0';
  const first = sampleRecords[0] || {};
  const last = sampleRecords[sampleRecords.length - 1] || {};
  return [
    sampleRecords.length,
    String(first.trackKey || ''),
    String(last.trackKey || '')
  ].join(':');
}

function buildTrackSampleViewCacheKey(sample, sampleRecords, trackState) {
  const safeSample = String(sample || '').trim() || 'All samples';
  const gffSourceKey = String(trackState?.gffSourceKey || '');
  const selectedField = getResolvedTrackGffLabelField(trackState);
  const gffSignature = summarizeTrackFeatureCollection(
    Array.isArray(trackState?.gffFeatures) ? trackState.gffFeatures : []
  );
  const sampleSignature = buildTrackSampleRecordSignature(sampleRecords);
  return `${safeSample}::${sampleSignature}::${gffSourceKey}::${selectedField}::${gffSignature}`;
}

function getCachedTrackSampleView(toolId, sample, sampleRecords, trackState) {
  const cache = getTrackRenderCache(toolId);
  const cacheKey = buildTrackSampleViewCacheKey(sample, sampleRecords, trackState);
  const cachedView = cache.sampleViewCache.get(cacheKey);
  if (cachedView && cachedView.sampleRecords === sampleRecords) {
    return cachedView;
  }

  const displayRecords = buildTrackDisplayRecords(sampleRecords, trackState);
  const recordByKey = new Map(displayRecords.map(record => [record.trackKey, record]));
  const geneSegmentData = buildTrackGeneSegments(displayRecords);
  const view = {
    cacheKey,
    sampleRecords,
    displayRecords,
    recordByKey,
    geneSegments: geneSegmentData.segments,
    geneMode: geneSegmentData.mode,
    overviewRecords: downsampleTrackPoints(displayRecords),
    mutationSearchOptions: buildTrackMutationSearchOptions(displayRecords),
    mutationSelectOptions: buildTrackMutationSelectOptions(displayRecords)
  };

  cache.sampleViewCache.set(cacheKey, view);
  if (cache.sampleViewCache.size > 24) {
    const oldest = cache.sampleViewCache.keys().next().value;
    cache.sampleViewCache.delete(oldest);
  }
  return view;
}

function patchTrackControls(controlsEl, {
  windowStartBp,
  windowEndBp,
  searchStatusText,
  searchHasMatches,
  searchQuery,
  selectedKey
}) {
  const startInput = controlsEl?.querySelector('.viz-track-window-start-input');
  if (startInput) startInput.value = String(windowStartBp);
  const endInput = controlsEl?.querySelector('.viz-track-window-end-input');
  if (endInput) endInput.value = String(windowEndBp);

  const searchStatusEl = controlsEl?.querySelector('.viz-track-search-status');
  if (searchStatusEl) searchStatusEl.textContent = String(searchStatusText || 'No search');
  const prevBtn = controlsEl?.querySelector('[data-track-action="search-prev"]');
  if (prevBtn) prevBtn.disabled = !searchHasMatches;
  const nextBtn = controlsEl?.querySelector('[data-track-action="search-next"]');
  if (nextBtn) nextBtn.disabled = !searchHasMatches;

  const searchInput = controlsEl?.querySelector('.viz-track-search-input');
  if (searchInput && searchInput.value !== String(searchQuery || '')) {
    searchInput.value = String(searchQuery || '');
  }

  const mutationSelectEl = controlsEl?.querySelector('.viz-track-mutation-select');
  if (mutationSelectEl) {
    mutationSelectEl.value = String(selectedKey || '');
  }
}

function buildTrackSvgDomKey(sampleView, trackState, fullSampleSequence, fullReferenceSequence) {
  const compact = value => {
    const num = Number(value);
    return Number.isFinite(num) ? Number(num.toFixed(3)) : 0;
  };
  const sampleSeqKey = fullSampleSequence
    ? `${fullSampleSequence.start}:${fullSampleSequence.end}:${String(fullSampleSequence.sequence || '').length}`
    : 'none';
  const referenceSeqKey = fullReferenceSequence
    ? `${fullReferenceSequence.start}:${fullReferenceSequence.end}:${String(fullReferenceSequence.sequence || '').length}`
    : 'none';
  return [
    sampleView?.cacheKey || '',
    String(trackState?.sequenceSourceKey || ''),
    String(trackState?.referenceSequenceSourceKey || ''),
    compact(trackState?.domainStart),
    compact(trackState?.domainEnd),
    compact(trackState?.windowStart),
    compact(trackState?.windowEnd),
    compact(trackState?.referenceSequenceTotalLength),
    String(trackState?.selectedKey || ''),
    String(trackState?.aaFocusGene || ''),
    sampleSeqKey,
    referenceSeqKey
  ].join('::');
}

function ensureTrackBodyShell(bodyEl) {
  if (!bodyEl) return { metaEl: null, svgHostEl: null, detailHostEl: null };
  let metaEl = bodyEl.querySelector('.viz-track-meta');
  let svgHostEl = bodyEl.querySelector('.viz-track-svg-host');
  let detailHostEl = bodyEl.querySelector('.viz-track-detail-host');
  if (metaEl && svgHostEl && detailHostEl) {
    return { metaEl, svgHostEl, detailHostEl };
  }
  bodyEl.innerHTML = `
    <div class="viz-track-meta"></div>
    <div class="viz-track-svg-host"></div>
    <div class="viz-track-detail-host"></div>
  `;
  metaEl = bodyEl.querySelector('.viz-track-meta');
  svgHostEl = bodyEl.querySelector('.viz-track-svg-host');
  detailHostEl = bodyEl.querySelector('.viz-track-detail-host');
  return { metaEl, svgHostEl, detailHostEl };
}

function getCachedTrackSvgInfo(renderCache, svgKey, buildFn) {
  const cache = renderCache?.svgInfoCache;
  if (!(cache instanceof Map)) {
    return buildFn();
  }
  const cached = cache.get(svgKey);
  if (cached) {
    // Refresh LRU order.
    cache.delete(svgKey);
    cache.set(svgKey, cached);
    return cached;
  }
  const built = buildFn();
  cache.set(svgKey, built);
  if (cache.size > 18) {
    const oldest = cache.keys().next().value;
    cache.delete(oldest);
  }
  return built;
}

function downsampleTrackPoints(records) {
  if (records.length <= TRACK_MAX_POINTS) return records;
  const sampled = [];
  const step = records.length / TRACK_MAX_POINTS;
  for (let i = 0; i < TRACK_MAX_POINTS; i++) {
    sampled.push(records[Math.floor(i * step)]);
  }
  return sampled;
}

const TRACK_LAYOUT = Object.freeze({
  width: 1120,
  height: 380,
  xStart: 150,
  xEnd: 1060,
  yRef: 66,
  yGenome: 140,
  yGene: 214,
  yAa: 252,
  yOverview: 332
});
const TRACK_PAN_STEP = 0.28;
const TRACK_ZOOM_IN_FACTOR = 0.5;
const TRACK_ZOOM_OUT_FACTOR = 1.75;
const TRACK_CLICK_FOCUS_ZOOM_FACTOR = 0.34;
const TRACK_SEARCH_FOCUS_ZOOM_FACTOR = 0.28;
const TRACK_MIN_WINDOW_BP = 6;
const TRACK_MIN_WINDOW_FRACTION = 0;
const TRACK_MIN_WINDOW_MAX_BP = 12;
const TRACK_SEQUENCE_LETTER_MIN_PX = 6.2;
const TRACK_REFERENCE_POINT_LABEL_MIN_GAP_PX = 34;

const NT_BASE_COLORS = Object.freeze({
  A: '#22c55e',
  C: '#3b82f6',
  G: '#f59e0b',
  T: '#ef4444',
  N: '#64748b'
});
const AA_THREE_TO_ONE = Object.freeze({
  ALA: 'A',
  ARG: 'R',
  ASN: 'N',
  ASP: 'D',
  CYS: 'C',
  GLN: 'Q',
  GLU: 'E',
  GLY: 'G',
  HIS: 'H',
  ILE: 'I',
  LEU: 'L',
  LYS: 'K',
  MET: 'M',
  PHE: 'F',
  PRO: 'P',
  SER: 'S',
  THR: 'T',
  TRP: 'W',
  TYR: 'Y',
  VAL: 'V',
  TER: '*',
  STOP: '*'
});
const AA_COLORS = Object.freeze({
  A: '#22c55e',
  R: '#2563eb',
  N: '#60a5fa',
  D: '#ef4444',
  C: '#14b8a6',
  Q: '#6366f1',
  E: '#f97316',
  G: '#f59e0b',
  H: '#a855f7',
  I: '#10b981',
  L: '#059669',
  K: '#3b82f6',
  M: '#84cc16',
  F: '#e11d48',
  P: '#8b5cf6',
  S: '#0ea5e9',
  T: '#06b6d4',
  W: '#7c3aed',
  Y: '#d946ef',
  V: '#16a34a',
  '*': '#64748b',
  X: '#94a3b8'
});
const CODON_TO_AA = Object.freeze({
  TTT: 'F', TTC: 'F', TTA: 'L', TTG: 'L',
  TCT: 'S', TCC: 'S', TCA: 'S', TCG: 'S',
  TAT: 'Y', TAC: 'Y', TAA: '*', TAG: '*',
  TGT: 'C', TGC: 'C', TGA: '*', TGG: 'W',
  CTT: 'L', CTC: 'L', CTA: 'L', CTG: 'L',
  CCT: 'P', CCC: 'P', CCA: 'P', CCG: 'P',
  CAT: 'H', CAC: 'H', CAA: 'Q', CAG: 'Q',
  CGT: 'R', CGC: 'R', CGA: 'R', CGG: 'R',
  ATT: 'I', ATC: 'I', ATA: 'I', ATG: 'M',
  ACT: 'T', ACC: 'T', ACA: 'T', ACG: 'T',
  AAT: 'N', AAC: 'N', AAA: 'K', AAG: 'K',
  AGT: 'S', AGC: 'S', AGA: 'R', AGG: 'R',
  GTT: 'V', GTC: 'V', GTA: 'V', GTG: 'V',
  GCT: 'A', GCC: 'A', GCA: 'A', GCG: 'A',
  GAT: 'D', GAC: 'D', GAA: 'E', GAG: 'E',
  GGT: 'G', GGC: 'G', GGA: 'G', GGG: 'G'
});
const NT_COMPLEMENT = Object.freeze({
  A: 'T',
  C: 'G',
  G: 'C',
  T: 'A',
  N: 'N'
});

function formatTrackBp(value) {
  const rounded = Math.round(Number(value) || 0);
  return rounded.toLocaleString('en-US');
}

function getTrackEffectiveSequenceOrientation(trackState) {
  return 'forward';
}

function getTrackSequenceSourceStateKey(source, orientation) {
  return `${source?.sourceKey || ''}::forward`;
}

function isTrackInformativeNt(value) {
  return /^[ACGTN]$/.test(String(value || '').trim().toUpperCase());
}

function getTrackExpectedSampleAllele(record) {
  const altAllele = String(record?.altAllele || '').trim().toUpperCase();
  if (isTrackInformativeNt(altAllele)) return altAllele;
  const markerKmer = String(record?.kmer || '').trim().toUpperCase();
  if (markerKmer) {
    const center = Math.floor(markerKmer.length / 2);
    const centerBase = markerKmer[center] || '';
    if (isTrackInformativeNt(centerBase)) return centerBase;
  }
  return '';
}

function getTrackAaSymbol(aaChange) {
  const text = normalizeValue(aaChange).toUpperCase();
  if (!text) return 'X';
  const chunks = text.match(/[A-Z\*]{1,5}/g);
  if (!chunks || chunks.length === 0) return 'X';
  const rawToken = chunks[chunks.length - 1];
  if (rawToken === '*' || rawToken === 'X') return rawToken;
  if (rawToken.length === 1) return rawToken;
  const mapped = AA_THREE_TO_ONE[rawToken.slice(0, 3)];
  if (mapped) return mapped;
  if (rawToken === 'STOP' || rawToken === 'TER') return '*';
  return rawToken[0] || 'X';
}

function getTrackAaColor(aaChange) {
  const symbol = getTrackAaSymbol(aaChange);
  return AA_COLORS[symbol] || AA_COLORS.X;
}

function getTrackComplementBase(base) {
  const normalized = normalizeNtBase(base);
  return NT_COMPLEMENT[normalized] || 'N';
}

function getTrackTranslatedAaFromCodon(codon) {
  const normalized = String(codon || '').trim().toUpperCase();
  if (normalized.length !== 3) return 'X';
  if (!/^[ACGTN]{3}$/.test(normalized)) return 'X';
  if (normalized.includes('N')) return 'X';
  return CODON_TO_AA[normalized] || 'X';
}

function getTrackAaFromSequenceCodon(sequenceData, codonStart, codonEnd, strand = '') {
  if (!sequenceData || !Number.isFinite(codonStart) || !Number.isFinite(codonEnd)) return '';
  const min = Math.min(codonStart, codonEnd);
  const max = Math.max(codonStart, codonEnd);
  const bases = [];
  for (let pos = min; pos <= max; pos++) {
    const base = getTrackSequenceBaseAtPosition(sequenceData, pos);
    if (!base) return '';
    bases.push(base);
  }
  if (bases.length !== 3) return '';
  let codon = bases.join('');
  if (strand === '-') {
    codon = bases.reverse().map(getTrackComplementBase).join('');
  }
  return getTrackTranslatedAaFromCodon(codon);
}

function getTrackSequenceBaseAtPosition(sequenceData, pos) {
  if (!sequenceData || !Number.isFinite(pos)) return '';
  const start = Number(sequenceData.start);
  const end = Number(sequenceData.end);
  const seq = String(sequenceData.sequence || '').toUpperCase();
  if (!Number.isFinite(start) || !Number.isFinite(end) || end < start || !seq) return '';
  const roundedPos = Math.round(pos);
  if (roundedPos < start || roundedPos > end) return '';
  const idx = roundedPos - start;
  if (idx < 0 || idx >= seq.length) return '';
  const nt = seq[idx];
  return isTrackInformativeNt(nt) ? nt : '';
}

function getTrackPathBasename(path) {
  const value = String(path || '').trim();
  if (!value) return '';
  const slash = Math.max(value.lastIndexOf('/'), value.lastIndexOf('\\'));
  return slash >= 0 ? value.slice(slash + 1) : value;
}

function getTrackPathStem(path) {
  const base = getTrackPathBasename(path).toLowerCase();
  if (!base) return '';
  return base
    .replace(/\.(fasta|fa|fna|fas)\.gz$/i, '')
    .replace(/\.(fasta|fa|fna|fas)$/i, '')
    .replace(/\.[^.]+$/, '');
}

function getTrackGffPathStem(path) {
  const base = getTrackPathBasename(path).toLowerCase();
  if (!base) return '';
  return base
    .replace(/\.(gff3?|gtf)\.gz$/i, '')
    .replace(/\.(gff3?|gtf)$/i, '')
    .replace(/\.[^.]+$/, '');
}

function parseTrackBracketSample(sample) {
  const match = String(sample || '').match(/^\[(.+?)\]\s+(.+)$/);
  if (!match) return null;
  return {
    fileName: match[1].trim(),
    recordName: match[2].trim()
  };
}

function resolveSampleMappedPath(sampleToPath, sample, bracketSample = null) {
  if (!sampleToPath || typeof sampleToPath !== 'object') return '';
  const sampleName = String(sample || '').trim();
  if (!sampleName) return '';

  const direct = String(sampleToPath[sampleName] || '').trim();
  if (direct) return direct;

  const lowerSample = sampleName.toLowerCase();
  const lowerRecord = String(bracketSample?.recordName || '').trim().toLowerCase();
  const lowerFile = String(bracketSample?.fileName || '').trim().toLowerCase();
  for (const [key, value] of Object.entries(sampleToPath)) {
    const mappedPath = String(value || '').trim();
    if (!mappedPath) continue;
    const normalizedKey = String(key || '').trim().toLowerCase();
    if (!normalizedKey) continue;
    if (normalizedKey === lowerSample) return mappedPath;
    if (lowerRecord && normalizedKey === lowerRecord) return mappedPath;
    if (lowerFile) {
      if (normalizedKey === lowerFile) return mappedPath;
      if (getTrackPathBasename(normalizedKey) === lowerFile) return mappedPath;
    }
  }

  return '';
}

function resolveSampleMappedFastaPath(sampleToFasta, sample, bracketSample = null) {
  return resolveSampleMappedPath(sampleToFasta, sample, bracketSample);
}

function resolveTrackSequenceSource(data, sample) {
  const ctx = data?.context?.trackContext;
  if (!ctx || !sample) return null;
  const fastaFiles = Array.isArray(ctx.fastaFiles) ? ctx.fastaFiles : [];
  const sampleToFasta = ctx.sampleToFasta && typeof ctx.sampleToFasta === 'object'
    ? ctx.sampleToFasta
    : {};

  const bracketSample = parseTrackBracketSample(sample);
  let fastaPath = '';
  let recordName = null;

  fastaPath = resolveSampleMappedFastaPath(sampleToFasta, sample, bracketSample);
  if (bracketSample?.recordName) {
    recordName = bracketSample.recordName;
  }

  if (!fastaPath && bracketSample && fastaFiles.length > 0) {
    const wanted = bracketSample.fileName.toLowerCase();
    const matched = fastaFiles.find(path => getTrackPathBasename(path).toLowerCase() === wanted);
    if (matched) {
      fastaPath = matched;
      recordName = bracketSample.recordName || null;
    }
  }

  if (!fastaPath && fastaFiles.length === 1) {
    fastaPath = fastaFiles[0];
    recordName = bracketSample?.recordName || null;
  }

  if (!fastaPath && fastaFiles.length > 1) {
    const normalizedSample = String(sample || '').trim().toLowerCase();
    const normalizedRecord = String(bracketSample?.recordName || '').trim().toLowerCase();
    const wantedNames = [normalizedSample, normalizedRecord].filter(Boolean);
    const matched = fastaFiles.find(path => {
      const base = getTrackPathBasename(path).toLowerCase();
      const stem = getTrackPathStem(path);
      return wantedNames.some(name => name === base || name === stem);
    });
    if (matched) {
      fastaPath = matched;
    }
  }

  if (!fastaPath) {
    const uniqueMappedPaths = Array.from(
      new Set(
        Object.values(sampleToFasta)
          .map(value => String(value || '').trim())
          .filter(Boolean)
      )
    );
    if (uniqueMappedPaths.length === 1) {
      fastaPath = uniqueMappedPaths[0];
    }
  }

  if (!fastaPath) return null;
  return {
    fastaPath,
    recordName: recordName || null,
    sourceKey: `${fastaPath}::${recordName || ''}`
  };
}

function resolveTrackGffSource(data, sample) {
  const ctx = data?.context?.trackContext;
  if (!ctx || !sample) return null;
  const sampleToGff = ctx.sampleToGff && typeof ctx.sampleToGff === 'object'
    ? ctx.sampleToGff
    : {};
  const gffFiles = Array.isArray(ctx.gffFiles) ? ctx.gffFiles : [];
  const singleGffPath = String(ctx.gffPath || '').trim();
  const bracketSample = parseTrackBracketSample(sample);
  let gffPath = resolveSampleMappedPath(sampleToGff, sample, bracketSample);

  if (!gffPath && bracketSample && gffFiles.length > 0) {
    const wanted = bracketSample.fileName.toLowerCase();
    const matched = gffFiles.find(path => getTrackPathBasename(path).toLowerCase() === wanted);
    if (matched) {
      gffPath = matched;
    }
  }

  if (!gffPath && gffFiles.length === 1) {
    gffPath = gffFiles[0];
  }

  if (!gffPath && gffFiles.length > 1) {
    const normalizedSample = String(sample || '').trim().toLowerCase();
    const normalizedRecord = String(bracketSample?.recordName || '').trim().toLowerCase();
    const wantedNames = [normalizedSample, normalizedRecord].filter(Boolean);
    const matched = gffFiles.find(path => {
      const base = getTrackPathBasename(path).toLowerCase();
      const stem = getTrackGffPathStem(path);
      return wantedNames.some(name => name === base || name === stem);
    });
    if (matched) {
      gffPath = matched;
    }
  }

  if (!gffPath) {
    const uniqueMappedPaths = Array.from(
      new Set(
        Object.values(sampleToGff)
          .map(value => String(value || '').trim())
          .filter(Boolean)
      )
    );
    if (uniqueMappedPaths.length === 1) {
      gffPath = uniqueMappedPaths[0];
    }
  }

  if (!gffPath && singleGffPath) {
    gffPath = singleGffPath;
  }

  if (!gffPath) return null;
  return {
    gffPath,
    sourceKey: gffPath
  };
}

function resolveReferenceTrackSequenceSource(data) {
  const referencePath = String(data?.context?.trackContext?.referencePath || '').trim();
  if (!referencePath) return null;
  return {
    fastaPath: referencePath,
    recordName: null,
    sourceKey: `${referencePath}::`
  };
}

function getTrackGffTypePriority(type) {
  const normalized = String(type || '').trim().toLowerCase();
  return Number.isFinite(TRACK_GFF_TYPE_PRIORITY[normalized]) ? TRACK_GFF_TYPE_PRIORITY[normalized] : 9;
}

function parseTrackGffAttributes(rawAttributes) {
  const attributes = {};
  const text = String(rawAttributes || '').trim();
  if (!text || text === '.') return attributes;

  text.split(';').forEach(part => {
    const token = String(part || '').trim();
    if (!token) return;

    let key = '';
    let value = '';
    const equalIndex = token.indexOf('=');
    if (equalIndex >= 0) {
      key = token.slice(0, equalIndex).trim();
      value = token.slice(equalIndex + 1).trim();
    } else {
      const spaceIndex = token.indexOf(' ');
      if (spaceIndex >= 0) {
        key = token.slice(0, spaceIndex).trim();
        value = token.slice(spaceIndex + 1).trim();
      } else {
        key = token.trim();
      }
    }
    if (!key) return;
    let clean = value.replace(/^"+|"+$/g, '').trim();
    if (clean.includes(',')) clean = clean.split(',')[0].trim();
    try {
      clean = decodeURIComponent(clean);
    } catch {
      // Leave as-is when decoding fails.
    }
    attributes[key] = clean;
  });

  return attributes;
}

function sortTrackGffAttributeKeys(keys) {
  const preferredOrder = ['locus_tag', 'gene', 'name', 'id', 'gene_id', 'product', 'parent'];
  return [...keys].sort((a, b) => {
    const aNorm = normalizeHeader(a);
    const bNorm = normalizeHeader(b);
    const aPref = preferredOrder.indexOf(aNorm);
    const bPref = preferredOrder.indexOf(bNorm);
    if (aPref !== -1 || bPref !== -1) {
      if (aPref === -1) return 1;
      if (bPref === -1) return -1;
      if (aPref !== bPref) return aPref - bPref;
    }
    return String(a).localeCompare(String(b), undefined, { sensitivity: 'base' });
  });
}

function chooseDefaultTrackGffLabelField(attributeKeys) {
  if (!Array.isArray(attributeKeys) || attributeKeys.length === 0) return '';
  const preferredOrder = ['locus_tag', 'gene', 'name', 'id', 'gene_id', 'product'];
  for (const preferred of preferredOrder) {
    const match = attributeKeys.find(key => normalizeHeader(key) === preferred);
    if (match) return match;
  }
  return attributeKeys[0];
}

function getTrackGffAttributeValue(feature, fieldName) {
  if (!feature || !fieldName) return '';
  const attributes = feature.attributes && typeof feature.attributes === 'object'
    ? feature.attributes
    : {};
  const direct = normalizeValue(attributes[fieldName]);
  if (direct) return direct;
  const wanted = normalizeHeader(fieldName);
  for (const [key, value] of Object.entries(attributes)) {
    if (normalizeHeader(key) !== wanted) continue;
    const normalized = normalizeValue(value);
    if (normalized) return normalized;
  }
  return '';
}

function getTrackGffFeatureLabel(feature, selectedField = '') {
  if (!feature) return '';
  const picked = selectedField ? getTrackGffAttributeValue(feature, selectedField) : '';
  if (picked) return picked;
  const fallbackKeys = ['locus_tag', 'gene', 'Name', 'ID', 'gene_id', 'product'];
  for (const key of fallbackKeys) {
    const value = getTrackGffAttributeValue(feature, key);
    if (value) return value;
  }
  const attributes = feature.attributes && typeof feature.attributes === 'object'
    ? feature.attributes
    : {};
  for (const value of Object.values(attributes)) {
    const normalized = normalizeValue(value);
    if (normalized) return normalized;
  }
  return '';
}

function parseTrackGffContent(content) {
  if (content === null || content === undefined) {
    throw new Error('GFF file could not be read.');
  }
  const lines = String(content).split(/\r?\n/);
  const allFeatures = [];
  for (const line of lines) {
    if (!line || line.startsWith('#')) continue;
    const cols = line.split('\t');
    if (cols.length < 9) continue;
    const type = String(cols[2] || '').trim();
    const start = Number.parseInt(cols[3], 10);
    const end = Number.parseInt(cols[4], 10);
    if (!Number.isFinite(start) || !Number.isFinite(end)) continue;
    const min = Math.max(1, Math.min(start, end));
    const max = Math.max(1, Math.max(start, end));
    const strandRaw = String(cols[6] || '').trim();
    const strand = strandRaw === '+' || strandRaw === '-' ? strandRaw : '';
    const phaseParsed = Number.parseInt(String(cols[7] || '').trim(), 10);
    const phase = (phaseParsed === 0 || phaseParsed === 1 || phaseParsed === 2) ? phaseParsed : null;
    const attributes = parseTrackGffAttributes(cols[8]);
    allFeatures.push({
      type,
      start: min,
      end: max,
      strand,
      phase,
      attributes
    });
  }

  if (allFeatures.length === 0) {
    return { features: [], attributeKeys: [] };
  }

  const filtered = allFeatures.filter(feature => getTrackGffTypePriority(feature.type) < 9);
  const features = (filtered.length > 0 ? filtered : allFeatures)
    .sort((a, b) => a.start - b.start || a.end - b.end);

  const keySet = new Set();
  features.forEach(feature => {
    Object.keys(feature.attributes || {}).forEach(key => {
      const normalized = normalizeValue(key);
      if (normalized) keySet.add(normalized);
    });
  });

  return {
    features,
    attributeKeys: sortTrackGffAttributeKeys(Array.from(keySet))
  };
}

async function loadTrackGffSource(source) {
  if (!source?.sourceKey || !source?.gffPath) {
    return { features: [], attributeKeys: [] };
  }

  const cached = TRACK_GFF_CACHE.get(source.sourceKey);
  if (cached?.status === 'ready') return cached.data;
  if (cached?.status === 'loading' && cached.promise) return cached.promise;
  if (cached?.status === 'error') throw new Error(cached.error || 'GFF cache error');

  const promise = (async () => {
    const content = await readTextFile(source.gffPath);
    const parsed = parseTrackGffContent(content);
    return parsed;
  })();

  TRACK_GFF_CACHE.set(source.sourceKey, { status: 'loading', promise });
  try {
    const parsed = await promise;
    TRACK_GFF_CACHE.set(source.sourceKey, { status: 'ready', data: parsed });
    return parsed;
  } catch (err) {
    TRACK_GFF_CACHE.set(source.sourceKey, { status: 'error', error: String(err?.message || err || '') });
    throw err;
  }
}

function hasTrackGffForSource(trackState, source) {
  if (!trackState || !source) return false;
  return (
    trackState.gffSourceKey === source.sourceKey &&
    Array.isArray(trackState.gffFeatures)
  );
}

function getResolvedTrackGffLabelField(trackState) {
  if (!trackState) return '';
  const keys = Array.isArray(trackState.gffAttributeKeys) ? trackState.gffAttributeKeys : [];
  if (keys.length === 0) return '';
  const selected = String(trackState.gffLabelField || '').trim();
  if (selected && keys.includes(selected)) return selected;
  return chooseDefaultTrackGffLabelField(keys);
}

function mapTrackRecordsToGffFeatures(sampleRecords, gffFeatures) {
  const mappings = new Map();
  if (!Array.isArray(sampleRecords) || sampleRecords.length === 0) return mappings;
  if (!Array.isArray(gffFeatures) || gffFeatures.length === 0) return mappings;

  const sortedRecords = [...sampleRecords].sort((a, b) => a.genomePos - b.genomePos);
  let featureIndex = 0;
  let active = [];
  for (const record of sortedRecords) {
    const pos = Math.round(Number(record.genomePos));
    if (!Number.isFinite(pos)) continue;
    while (featureIndex < gffFeatures.length && gffFeatures[featureIndex].start <= pos) {
      active.push(gffFeatures[featureIndex]);
      featureIndex += 1;
    }
    if (active.length > 0) {
      active = active.filter(feature => feature.end >= pos);
    }
    let best = null;
    let bestPriority = Number.POSITIVE_INFINITY;
    let bestSpan = Number.POSITIVE_INFINITY;
    active.forEach(feature => {
      if (feature.start > pos || feature.end < pos) return;
      const priority = getTrackGffTypePriority(feature.type);
      const span = Math.max(1, feature.end - feature.start);
      if (priority < bestPriority || (priority === bestPriority && span < bestSpan)) {
        best = feature;
        bestPriority = priority;
        bestSpan = span;
      }
    });
    if (best && record.trackKey) {
      mappings.set(record.trackKey, best);
    }
  }

  return mappings;
}

function getTrackRecordGeneLabel(record) {
  return normalizeValue(record?.displayGene || record?.gene);
}

function getTrackRecordGeneStart(record) {
  // Test finiteness on the raw value: Number(null) and Number('') are 0, which
  // would defeat the null guard and read a missing gene start as position 0.
  const raw = record?.displayGeneStart;
  if (Number.isFinite(raw)) return raw;
  return Number.isFinite(record?.geneStart) ? record.geneStart : null;
}

function getTrackRecordGeneEnd(record) {
  const raw = record?.displayGeneEnd;
  if (Number.isFinite(raw)) return raw;
  return Number.isFinite(record?.geneEnd) ? record.geneEnd : null;
}

function getTrackRecordGeneStrand(record) {
  const strand = normalizeValue(record?.displayGeneStrand);
  if (strand === '+' || strand === '-') return strand;
  return '';
}

function parseTrackAaPosition(record) {
  const text = normalizeValue(record?.aaPos);
  if (!text) return null;
  const match = text.match(/\d+/);
  if (!match) return null;
  const parsed = Number.parseInt(match[0], 10);
  return Number.isFinite(parsed) && parsed > 0 ? parsed : null;
}

function resolveTrackAaGeneContext(sampleRecords, visibleRaw, trackState) {
  if (!Array.isArray(sampleRecords) || sampleRecords.length === 0) return null;
  const selectedRecord = sampleRecords.find(record => record.trackKey === trackState?.selectedKey) || null;
  const selectedGeneLabel = getTrackRecordGeneLabel(selectedRecord);
  const focusedGeneLabel = normalizeValue(trackState?.aaFocusGene);
  let geneLabel = selectedGeneLabel || focusedGeneLabel;

  if (!geneLabel) {
    const visibleGeneLabels = Array.from(new Set(
      (Array.isArray(visibleRaw) ? visibleRaw : [])
        .map(getTrackRecordGeneLabel)
        .filter(label => label && !isUnknownOutcome(label))
    ));
    if (visibleGeneLabels.length >= 1) {
      [geneLabel] = visibleGeneLabels;
    }
  }
  if (!geneLabel) return null;

  const geneRecords = sampleRecords.filter(record => getTrackRecordGeneLabel(record) === geneLabel);
  if (geneRecords.length === 0) return null;

  const starts = geneRecords.map(getTrackRecordGeneStart).filter(Number.isFinite);
  const ends = geneRecords.map(getTrackRecordGeneEnd).filter(Number.isFinite);
  if (starts.length === 0 || ends.length === 0) return null;
  const geneMin = Math.min(...starts, ...ends);
  const geneMax = Math.max(...starts, ...ends);

  const strandCounts = new Map();
  geneRecords.forEach(record => {
    const strand = normalizeValue(record.displayGeneStrand);
    if (!strand) return;
    strandCounts.set(strand, (strandCounts.get(strand) || 0) + 1);
  });
  const strand = Array.from(strandCounts.entries())
    .sort((a, b) => b[1] - a[1])[0]?.[0] || '';

  const phaseCounts = new Map();
  geneRecords.forEach(record => {
    const phaseRaw = Number(record.displayGenePhase);
    if (phaseRaw !== 0 && phaseRaw !== 1 && phaseRaw !== 2) return;
    phaseCounts.set(phaseRaw, (phaseCounts.get(phaseRaw) || 0) + 1);
  });
  const phase = Array.from(phaseCounts.entries())
    .sort((a, b) => b[1] - a[1])[0]?.[0] ?? 0;

  return {
    label: geneLabel,
    start: geneMin,
    end: geneMax,
    strand: strand === '-' ? '-' : '+',
    phase
  };
}

function buildTrackGeneCodonSeries(geneContext) {
  if (!geneContext) return [];
  const start = Math.min(geneContext.start, geneContext.end);
  const end = Math.max(geneContext.start, geneContext.end);
  const phase = (geneContext.phase === 0 || geneContext.phase === 1 || geneContext.phase === 2)
    ? geneContext.phase
    : 0;

  const codons = [];
  if (geneContext.strand === '-') {
    let codonEnd = end - phase;
    let aaIndex = 1;
    while ((codonEnd - 2) >= start) {
      codons.push({ start: codonEnd - 2, end: codonEnd, aaIndex });
      codonEnd -= 3;
      aaIndex += 1;
    }
  } else {
    let codonStart = start + phase;
    let aaIndex = 1;
    while ((codonStart + 2) <= end) {
      codons.push({ start: codonStart, end: codonStart + 2, aaIndex });
      codonStart += 3;
      aaIndex += 1;
    }
  }
  return codons;
}

function buildTrackDisplayRecords(sampleRecords, trackState) {
  if (!Array.isArray(sampleRecords) || sampleRecords.length === 0) return [];
  const gffFeatures = Array.isArray(trackState?.gffFeatures) ? trackState.gffFeatures : [];
  const selectedField = getResolvedTrackGffLabelField(trackState);
  const mappedFeatures = mapTrackRecordsToGffFeatures(sampleRecords, gffFeatures);

  return sampleRecords.map(record => {
    const matchedFeature = mappedFeatures.get(record.trackKey) || null;
    const gffLabel = getTrackGffFeatureLabel(matchedFeature, selectedField);
    return {
      ...record,
      displayGene: gffLabel || record.gene || '',
      displayGeneStart: matchedFeature ? matchedFeature.start : record.geneStart,
      displayGeneEnd: matchedFeature ? matchedFeature.end : record.geneEnd,
      displayGeneType: matchedFeature?.type || '',
      displayGeneStrand: matchedFeature?.strand || '',
      displayGenePhase: Number.isFinite(matchedFeature?.phase) ? matchedFeature.phase : null
    };
  });
}

function getTrackDomain(records, totalLength = null) {
  const values = [1];
  records.forEach(record => {
    values.push(record.genomePos);
    if (record.geneStart !== null) values.push(record.geneStart);
    if (record.geneEnd !== null) values.push(record.geneEnd);
  });
  const knownTotal = Number.isFinite(totalLength) ? Number(totalLength) : null;
  if (knownTotal !== null && knownTotal > 0) values.push(knownTotal);
  const domainStart = 1;
  const domainEndRaw = Math.max(...values);
  const domainEnd = domainEndRaw <= domainStart ? domainStart + 1 : domainEndRaw;
  return { domainStart, domainEnd };
}

function getTrackMinWindow(domainStart, domainEnd) {
  const span = Math.max(1, domainEnd - domainStart);
  const scaled = Math.max(TRACK_MIN_WINDOW_BP, span * TRACK_MIN_WINDOW_FRACTION);
  return Math.min(TRACK_MIN_WINDOW_MAX_BP, scaled);
}

function clampTrackWindow(windowStart, windowEnd, domainStart, domainEnd) {
  const minWindow = getTrackMinWindow(domainStart, domainEnd);
  const domainSpan = Math.max(minWindow, domainEnd - domainStart);
  let start = Number.isFinite(windowStart) ? windowStart : domainStart;
  let end = Number.isFinite(windowEnd) ? windowEnd : domainEnd;
  if (end < start) [start, end] = [end, start];
  if (end - start < minWindow) {
    const center = (start + end) / 2;
    start = center - minWindow / 2;
    end = center + minWindow / 2;
  }
  if (start < domainStart) {
    const shift = domainStart - start;
    start += shift;
    end += shift;
  }
  if (end > domainEnd) {
    const shift = end - domainEnd;
    start -= shift;
    end -= shift;
  }
  if (start < domainStart) start = domainStart;
  if (end > domainEnd) end = domainEnd;
  if (end - start < minWindow) {
    start = domainStart;
    end = Math.min(domainEnd, domainStart + Math.min(domainSpan, minWindow));
  }
  return { windowStart: start, windowEnd: end };
}

function zoomTrackWindow(windowStart, windowEnd, factor, domainStart, domainEnd, anchor = null) {
  const currentSpan = Math.max(getTrackMinWindow(domainStart, domainEnd), windowEnd - windowStart);
  const nextSpan = Math.max(getTrackMinWindow(domainStart, domainEnd), Math.min(domainEnd - domainStart, currentSpan * factor));
  const center = Number.isFinite(anchor) ? anchor : (windowStart + windowEnd) / 2;
  return clampTrackWindow(center - (nextSpan / 2), center + (nextSpan / 2), domainStart, domainEnd);
}

function panTrackWindow(windowStart, windowEnd, direction, domainStart, domainEnd) {
  const span = windowEnd - windowStart;
  const offset = span * TRACK_PAN_STEP * direction;
  return clampTrackWindow(windowStart + offset, windowEnd + offset, domainStart, domainEnd);
}

function getNextWindowFromAction(action, current, controlsEl) {
  if (!current) return null;
  if (action === 'zoom-in') {
    return zoomTrackWindow(current.windowStart, current.windowEnd, TRACK_ZOOM_IN_FACTOR, current.domainStart, current.domainEnd);
  }
  if (action === 'zoom-out') {
    return zoomTrackWindow(current.windowStart, current.windowEnd, TRACK_ZOOM_OUT_FACTOR, current.domainStart, current.domainEnd);
  }
  if (action === 'pan-left') {
    return panTrackWindow(current.windowStart, current.windowEnd, -1, current.domainStart, current.domainEnd);
  }
  if (action === 'pan-right') {
    return panTrackWindow(current.windowStart, current.windowEnd, 1, current.domainStart, current.domainEnd);
  }
  if (action === 'reset') {
    return { windowStart: current.domainStart, windowEnd: current.domainEnd };
  }
  if (action === 'apply-range') {
    const startInput = controlsEl?.querySelector('.viz-track-window-start-input');
    const endInput = controlsEl?.querySelector('.viz-track-window-end-input');
    const startValue = parseNumericValue(startInput?.value || '');
    const endValue = parseNumericValue(endInput?.value || '');
    if (startValue === null && endValue === null) return null;

    let nextStart = Number.isFinite(startValue) ? Math.floor(startValue) : Math.floor(current.windowStart);
    let nextEnd = Number.isFinite(endValue) ? Math.ceil(endValue) : Math.ceil(current.windowEnd);

    if (nextStart > nextEnd) {
      const temp = nextStart;
      nextStart = nextEnd;
      nextEnd = temp;
    }

    if (nextEnd === nextStart) {
      const minSpan = getTrackMinWindow(current.domainStart, current.domainEnd);
      nextEnd = nextStart + minSpan;
    }

    return clampTrackWindow(nextStart, nextEnd, current.domainStart, current.domainEnd);
  }
  return null;
}

function hasTrackSequenceForWindow(trackState, source, windowStart, windowEnd) {
  if (!trackState || !source || !trackState.sequence) return false;
  const effectiveOrientation = getTrackEffectiveSequenceOrientation(trackState);
  const expectedSourceKey = getTrackSequenceSourceStateKey(source, effectiveOrientation);
  if (trackState.sequenceSourceKey !== expectedSourceKey) return false;
  if (!Number.isFinite(trackState.sequenceStart) || !Number.isFinite(trackState.sequenceEnd)) return false;
  return (
    trackState.sequenceStart <= Math.floor(windowStart) &&
    trackState.sequenceEnd >= Math.ceil(windowEnd)
  );
}

function hasReferenceTrackSequenceForWindow(trackState, source, windowStart, windowEnd) {
  if (!trackState || !source || !trackState.referenceSequence) return false;
  if (trackState.referenceSequenceSourceKey !== source.sourceKey) return false;
  if (!Number.isFinite(trackState.referenceSequenceStart) || !Number.isFinite(trackState.referenceSequenceEnd)) return false;
  return (
    trackState.referenceSequenceStart <= Math.floor(windowStart) &&
    trackState.referenceSequenceEnd >= Math.ceil(windowEnd)
  );
}

function buildTrackSequenceRequestForRange(
  windowStart,
  windowEnd,
  source
) {
  if (!source) return null;
  const span = Math.max(1, windowEnd - windowStart);
  if (span > TRACK_SEQUENCE_FETCH_MAX_SPAN) return null;

  const padding = Math.min(
    TRACK_SEQUENCE_FETCH_PADDING_MAX,
    Math.max(TRACK_SEQUENCE_FETCH_PADDING_MIN, Math.round(span * TRACK_SEQUENCE_FETCH_PADDING_FRACTION))
  );
  const targetStart = Math.max(1, Math.floor(windowStart) - padding);
  const targetEnd = Math.max(targetStart, Math.ceil(windowEnd) + padding);
  const requestStart = targetStart;
  const requestEnd = targetEnd;
  return {
    requestStart,
    requestEnd,
    targetStart,
    targetEnd,
    orientation: 'forward',
    mappedReverse: false,
    requestKey: `${source.sourceKey}:forward:${targetStart}:${targetEnd}:${requestStart}:${requestEnd}`
  };
}

function buildTrackSequenceRequest(trackState, source) {
  if (!trackState) return null;
  return buildTrackSequenceRequestForRange(
    trackState.windowStart,
    trackState.windowEnd,
    source
  );
}

async function fetchTrackSequence(toolId, data, source) {
  const current = TRACK_STATE[toolId];
  if (!current || current.sample !== source.sample) return;

  const request = buildTrackSequenceRequest(current, source);
  if (!request) return;
  const { requestStart, requestEnd, requestKey } = request;

  TRACK_STATE[toolId] = {
    ...current,
    sequenceLoading: true,
    sequenceRequestKey: requestKey,
    sequenceError: ''
  };

  const response = await readFastaRangeCached(
    source.fastaPath,
    requestStart,
    requestEnd,
    source.recordName
  );
  const latest = TRACK_STATE[toolId];
  if (!latest || latest.sample !== source.sample || latest.sequenceRequestKey !== requestKey) return;

  if (!response || !response.sequence) {
    TRACK_STATE[toolId] = {
      ...latest,
      sequenceLoading: false,
      sequenceError: 'Could not load sample FASTA sequence for this window.'
    };
    renderGenomicTrack(toolId, data);
    return;
  }

  const responseTotalLength = Number.isFinite(response.total_length)
    ? Number(response.total_length)
    : (Number.isFinite(current.sequenceTotalLength) ? Number(current.sequenceTotalLength) : null);
  const sequenceStart = Number.isFinite(response.start) ? response.start : request.requestStart;
  const sequenceEnd = Number.isFinite(response.end) ? response.end : request.requestEnd;
  const sequence = String(response.sequence || '').toUpperCase();

  TRACK_STATE[toolId] = {
    ...latest,
    sequenceOrientation: 'forward',
    sequence,
    sequenceStart,
    sequenceEnd,
    sequenceTotalLength: Number.isFinite(responseTotalLength) ? responseTotalLength : (latest.sequenceTotalLength || null),
    sequenceSourceKey: getTrackSequenceSourceStateKey(source, 'forward'),
    sequenceRecordName: response.record_name || source.recordName || '',
    sequenceLoading: false,
    sequenceRequestKey: '',
    sequenceError: ''
  };
  renderGenomicTrack(toolId, data);
}

async function fetchTrackSequenceLengthProbe(toolId, data, source) {
  const current = TRACK_STATE[toolId];
  if (!current || current.sample !== source.sample) return;
  const requestKey = `${source.sourceKey}:probe-length`;
  TRACK_STATE[toolId] = {
    ...current,
    sequenceLoading: true,
    sequenceRequestKey: requestKey,
    sequenceError: ''
  };

  const response = await readFastaRangeCached(source.fastaPath, 1, 1, source.recordName);
  const latest = TRACK_STATE[toolId];
  if (!latest || latest.sample !== source.sample || latest.sequenceRequestKey !== requestKey) return;

  if (!response || !Number.isFinite(response.total_length)) {
    TRACK_STATE[toolId] = {
      ...latest,
      sequenceLoading: false,
      sequenceRequestKey: '',
      sequenceError: 'Could not determine sample genome length from FASTA.'
    };
    renderGenomicTrack(toolId, data);
    return;
  }

  TRACK_STATE[toolId] = {
    ...latest,
    sequenceTotalLength: response.total_length,
    sequenceLoading: false,
    sequenceRequestKey: '',
    sequenceError: ''
  };
  renderGenomicTrack(toolId, data);
}

async function fetchTrackReferenceSequence(toolId, data, source, sample, windowStart, windowEnd) {
  const current = TRACK_STATE[toolId];
  if (!current || current.sample !== sample) return;

  const request = buildTrackSequenceRequestForRange(windowStart, windowEnd, source);
  if (!request) return;
  const { requestStart, requestEnd, requestKey } = request;

  TRACK_STATE[toolId] = {
    ...current,
    referenceSequenceLoading: true,
    referenceSequenceRequestKey: requestKey,
    referenceSequenceError: ''
  };

  const response = await readFastaRangeCached(
    source.fastaPath,
    requestStart,
    requestEnd,
    source.recordName
  );
  const latest = TRACK_STATE[toolId];
  if (!latest || latest.sample !== sample || latest.referenceSequenceRequestKey !== requestKey) return;

  if (!response || !response.sequence) {
    TRACK_STATE[toolId] = {
      ...latest,
      referenceSequenceLoading: false,
      referenceSequenceError: 'Could not load reference FASTA sequence for this window.'
    };
    renderGenomicTrack(toolId, data);
    return;
  }

  TRACK_STATE[toolId] = {
    ...latest,
    referenceSequence: String(response.sequence || '').toUpperCase(),
    referenceSequenceStart: Number.isFinite(response.start) ? response.start : requestStart,
    referenceSequenceEnd: Number.isFinite(response.end) ? response.end : requestEnd,
    referenceSequenceTotalLength: Number.isFinite(response.total_length)
      ? response.total_length
      : (latest.referenceSequenceTotalLength || null),
    referenceSequenceSourceKey: source.sourceKey,
    referenceSequenceRecordName: response.record_name || source.recordName || '',
    referenceSequenceLoading: false,
    referenceSequenceRequestKey: '',
    referenceSequenceError: ''
  };
  renderGenomicTrack(toolId, data);
}

async function fetchTrackReferenceLengthProbe(toolId, data, source, sample) {
  const current = TRACK_STATE[toolId];
  if (!current || current.sample !== sample) return;
  const requestKey = `${source.sourceKey}:probe-length`;
  TRACK_STATE[toolId] = {
    ...current,
    referenceSequenceLoading: true,
    referenceSequenceRequestKey: requestKey,
    referenceSequenceError: ''
  };

  const response = await readFastaRangeCached(source.fastaPath, 1, 1, source.recordName);
  const latest = TRACK_STATE[toolId];
  if (!latest || latest.sample !== sample || latest.referenceSequenceRequestKey !== requestKey) return;

  if (!response || !Number.isFinite(response.total_length)) {
    TRACK_STATE[toolId] = {
      ...latest,
      referenceSequenceLoading: false,
      referenceSequenceRequestKey: '',
      referenceSequenceError: 'Could not determine reference genome length from FASTA.'
    };
    renderGenomicTrack(toolId, data);
    return;
  }

  TRACK_STATE[toolId] = {
    ...latest,
    referenceSequenceTotalLength: response.total_length,
    referenceSequenceLoading: false,
    referenceSequenceRequestKey: '',
    referenceSequenceError: ''
  };
  renderGenomicTrack(toolId, data);
}

async function fetchTrackGffAnnotations(toolId, data, source, sample) {
  const current = TRACK_STATE[toolId];
  if (!current || current.sample !== sample) return;
  const requestKey = `${source.sourceKey}:gff`;
  TRACK_STATE[toolId] = {
    ...current,
    gffSourceKey: source.sourceKey,
    gffLoading: true,
    gffRequestKey: requestKey,
    gffError: ''
  };

  try {
    const parsed = await loadTrackGffSource(source);
    const latest = TRACK_STATE[toolId];
    if (!latest || latest.sample !== sample || latest.gffRequestKey !== requestKey) return;
    const attributeKeys = Array.isArray(parsed.attributeKeys) ? parsed.attributeKeys : [];
    const selectedField = String(latest.gffLabelField || '').trim();
    TRACK_STATE[toolId] = {
      ...latest,
      gffSourceKey: source.sourceKey,
      gffFeatures: Array.isArray(parsed.features) ? parsed.features : [],
      gffAttributeKeys: attributeKeys,
      gffLabelField: selectedField && attributeKeys.includes(selectedField) ? selectedField : '',
      gffLoading: false,
      gffRequestKey: '',
      gffError: ''
    };
  } catch (err) {
    const latest = TRACK_STATE[toolId];
    if (!latest || latest.sample !== sample || latest.gffRequestKey !== requestKey) return;
    TRACK_STATE[toolId] = {
      ...latest,
      gffFeatures: [],
      gffAttributeKeys: [],
      gffLoading: false,
      gffRequestKey: '',
      gffError: 'Could not load sample GFF annotations for this sample.'
    };
  }
  renderGenomicTrack(toolId, data);
}

function scheduleTrackSequenceFetch(toolId, data, source) {
  const timerKey = `${toolId}:sample`;
  if (TRACK_FETCH_TIMERS[timerKey]) {
    clearTimeout(TRACK_FETCH_TIMERS[timerKey]);
  }
  TRACK_FETCH_TIMERS[timerKey] = setTimeout(() => {
    delete TRACK_FETCH_TIMERS[timerKey];
    fetchTrackSequence(toolId, data, source);
  }, TRACK_SEQUENCE_FETCH_DEBOUNCE_MS);
}

function scheduleTrackSequenceLengthProbe(toolId, data, source) {
  const timerKey = `${toolId}:sample-probe`;
  if (TRACK_FETCH_TIMERS[timerKey]) {
    clearTimeout(TRACK_FETCH_TIMERS[timerKey]);
  }
  TRACK_FETCH_TIMERS[timerKey] = setTimeout(() => {
    delete TRACK_FETCH_TIMERS[timerKey];
    fetchTrackSequenceLengthProbe(toolId, data, source);
  }, TRACK_SEQUENCE_FETCH_DEBOUNCE_MS);
}

function scheduleTrackReferenceSequenceFetch(toolId, data, source, sample, windowStart, windowEnd) {
  const timerKey = `${toolId}:reference`;
  if (TRACK_FETCH_TIMERS[timerKey]) {
    clearTimeout(TRACK_FETCH_TIMERS[timerKey]);
  }
  TRACK_FETCH_TIMERS[timerKey] = setTimeout(() => {
    delete TRACK_FETCH_TIMERS[timerKey];
    fetchTrackReferenceSequence(toolId, data, source, sample, windowStart, windowEnd);
  }, TRACK_SEQUENCE_FETCH_DEBOUNCE_MS);
}

function scheduleTrackReferenceLengthProbe(toolId, data, source, sample) {
  const timerKey = `${toolId}:reference-probe`;
  if (TRACK_FETCH_TIMERS[timerKey]) {
    clearTimeout(TRACK_FETCH_TIMERS[timerKey]);
  }
  TRACK_FETCH_TIMERS[timerKey] = setTimeout(() => {
    delete TRACK_FETCH_TIMERS[timerKey];
    fetchTrackReferenceLengthProbe(toolId, data, source, sample);
  }, TRACK_SEQUENCE_FETCH_DEBOUNCE_MS);
}

function scheduleTrackGffFetch(toolId, data, source, sample) {
  const timerKey = `${toolId}:gff`;
  if (TRACK_FETCH_TIMERS[timerKey]) {
    clearTimeout(TRACK_FETCH_TIMERS[timerKey]);
  }
  TRACK_FETCH_TIMERS[timerKey] = setTimeout(() => {
    delete TRACK_FETCH_TIMERS[timerKey];
    fetchTrackGffAnnotations(toolId, data, source, sample);
  }, TRACK_SEQUENCE_FETCH_DEBOUNCE_MS);
}

function buildTrackGeneSegments(records) {
  const map = new Map();
  let mode = 'inferred';
  records.forEach(record => {
    const geneLabel = getTrackRecordGeneLabel(record);
    if (!geneLabel || isUnknownOutcome(geneLabel)) return;
    const geneStart = getTrackRecordGeneStart(record);
    const geneEnd = getTrackRecordGeneEnd(record);
    if (geneStart === null || geneEnd === null) return;
    const min = Math.min(geneStart, geneEnd);
    const max = Math.max(geneStart, geneEnd);
    const strand = getTrackRecordGeneStrand(record);
    const prev = map.get(geneLabel);
    if (!prev) {
      map.set(geneLabel, {
        gene: geneLabel,
        min,
        max,
        count: 1,
        strandPlus: strand === '+' ? 1 : 0,
        strandMinus: strand === '-' ? 1 : 0
      });
    } else {
      prev.min = Math.min(prev.min, min);
      prev.max = Math.max(prev.max, max);
      prev.count += 1;
      if (strand === '+') prev.strandPlus += 1;
      if (strand === '-') prev.strandMinus += 1;
    }
  });
  if (map.size === 0) {
    records.forEach(record => {
      const geneLabel = getTrackRecordGeneLabel(record);
      if (!geneLabel || isUnknownOutcome(geneLabel)) return;
      const strand = getTrackRecordGeneStrand(record);
      const prev = map.get(geneLabel);
      if (!prev) {
        map.set(geneLabel, {
          gene: geneLabel,
          min: record.genomePos,
          max: record.genomePos,
          count: 1,
          strandPlus: strand === '+' ? 1 : 0,
          strandMinus: strand === '-' ? 1 : 0
        });
      } else {
        prev.min = Math.min(prev.min, record.genomePos);
        prev.max = Math.max(prev.max, record.genomePos);
        prev.count += 1;
        if (strand === '+') prev.strandPlus += 1;
        if (strand === '-') prev.strandMinus += 1;
      }
    });
  } else {
    mode = 'gff';
  }
  const segments = Array.from(map.values())
    .map(segment => ({
      ...segment,
      strand: (segment.strandPlus > 0 || segment.strandMinus > 0)
        ? (segment.strandPlus >= segment.strandMinus ? '+' : '-')
        : ''
    }))
    .sort((a, b) => b.count - a.count || a.gene.localeCompare(b.gene, undefined, { sensitivity: 'base' }));
  return { segments, mode };
}

function ensureTrackState(toolId, sample, sampleRecords) {
  const current = TRACK_STATE[toolId] || {};
  const sampleChanged = current.sample !== sample;
  const autoFollowDomain = sampleChanged
    ? true
    : (typeof current.autoFollowDomain === 'boolean' ? current.autoFollowDomain : true);
  const sequenceTotalLength = sampleChanged ? null : current.sequenceTotalLength;
  const { domainStart, domainEnd } = getTrackDomain(sampleRecords, sequenceTotalLength);
  const previousDomainEnd = Number.isFinite(current.domainEnd) ? current.domainEnd : null;
  let windowStart = current.windowStart;
  let windowEnd = current.windowEnd;
  if (sampleChanged || !Number.isFinite(windowStart) || !Number.isFinite(windowEnd)) {
    windowStart = domainStart;
    windowEnd = domainEnd;
  } else if (
    autoFollowDomain &&
    previousDomainEnd !== null &&
    domainEnd > previousDomainEnd &&
    Math.abs(windowStart - domainStart) < 2 &&
    Math.abs(windowEnd - previousDomainEnd) < 3
  ) {
    // Auto-extend to full domain only when user is still in "whole-genome view".
    windowStart = domainStart;
    windowEnd = domainEnd;
  }
  ({ windowStart, windowEnd } = clampTrackWindow(windowStart, windowEnd, domainStart, domainEnd));
  const selectedKey = sampleRecords.some(record => record.trackKey === current.selectedKey) ? current.selectedKey : '';
  const searchQuery = sampleChanged ? '' : String(current.searchQuery || '');
  const validRecordKeys = new Set(sampleRecords.map(record => record.trackKey));
  const searchMatchKeys = sampleChanged
    ? []
    : (Array.isArray(current.searchMatchKeys)
      ? current.searchMatchKeys.filter(key => validRecordKeys.has(key))
      : []);
  const rawSearchMatchIndex = Number(current.searchMatchIndex);
  const searchMatchIndex = searchMatchKeys.length === 0
    ? -1
    : Math.min(
        Math.max(Number.isFinite(rawSearchMatchIndex) ? Math.trunc(rawSearchMatchIndex) : 0, 0),
        searchMatchKeys.length - 1
      );
  const next = {
    sample,
    domainStart,
    domainEnd,
    windowStart,
    windowEnd,
    autoFollowDomain,
    selectedKey,
    searchQuery,
    searchMatchKeys,
    searchMatchIndex,
    sequenceOrientation: 'forward',
    sequence: sampleChanged ? '' : (current.sequence || ''),
    sequenceStart: sampleChanged ? null : current.sequenceStart,
    sequenceEnd: sampleChanged ? null : current.sequenceEnd,
    sequenceTotalLength: sampleChanged ? null : current.sequenceTotalLength,
    sequenceSourceKey: sampleChanged ? '' : (current.sequenceSourceKey || ''),
    sequenceRecordName: sampleChanged ? '' : (current.sequenceRecordName || ''),
    sequenceLoading: sampleChanged ? false : Boolean(current.sequenceLoading),
    sequenceRequestKey: sampleChanged ? '' : (current.sequenceRequestKey || ''),
    sequenceError: sampleChanged ? '' : (current.sequenceError || ''),
    referenceSequence: sampleChanged ? '' : (current.referenceSequence || ''),
    referenceSequenceStart: sampleChanged ? null : current.referenceSequenceStart,
    referenceSequenceEnd: sampleChanged ? null : current.referenceSequenceEnd,
    referenceSequenceTotalLength: sampleChanged ? null : current.referenceSequenceTotalLength,
    referenceSequenceSourceKey: sampleChanged ? '' : (current.referenceSequenceSourceKey || ''),
    referenceSequenceRecordName: sampleChanged ? '' : (current.referenceSequenceRecordName || ''),
    referenceSequenceLoading: sampleChanged ? false : Boolean(current.referenceSequenceLoading),
    referenceSequenceRequestKey: sampleChanged ? '' : (current.referenceSequenceRequestKey || ''),
    referenceSequenceError: sampleChanged ? '' : (current.referenceSequenceError || ''),
    gffFeatures: sampleChanged
      ? []
      : (Array.isArray(current.gffFeatures) ? current.gffFeatures : []),
    gffAttributeKeys: sampleChanged
      ? []
      : (Array.isArray(current.gffAttributeKeys) ? current.gffAttributeKeys : []),
    gffSourceKey: sampleChanged ? '' : (current.gffSourceKey || ''),
    gffLoading: sampleChanged ? false : Boolean(current.gffLoading),
    gffRequestKey: sampleChanged ? '' : (current.gffRequestKey || ''),
    gffError: sampleChanged ? '' : (current.gffError || ''),
    gffLabelField: sampleChanged ? '' : (current.gffLabelField || ''),
    aaFocusGene: sampleChanged ? '' : (current.aaFocusGene || '')
  };
  TRACK_STATE[toolId] = next;
  return next;
}

function buildNtBoxesFromSequence(sequenceData, toX, rangeStart, rangeEnd, y, span, rowClass, forcedPositions = []) {
  if (!sequenceData?.sequence) {
    return { svg: '', showNt: false, source: '', lettersVisible: false };
  }
  const seqStart = Number(sequenceData.start);
  const seqEnd = Number(sequenceData.end);
  if (!Number.isFinite(seqStart) || !Number.isFinite(seqEnd) || seqEnd < seqStart) {
    return { svg: '', showNt: false, source: '', lettersVisible: false };
  }
  const seq = String(sequenceData.sequence || '');
  if (!seq) return { svg: '', showNt: false, source: '', lettersVisible: false };

  const minPos = Math.max(Math.floor(rangeStart), seqStart);
  const maxPos = Math.min(Math.ceil(rangeEnd), seqEnd);
  if (maxPos < minPos) return { svg: '', showNt: false, source: '', lettersVisible: false };

  const baseStep = span <= 120 ? 1 : span <= 240 ? 2 : span <= 420 ? 3 : span <= 700 ? 4 : 6;
  const densityStep = Math.max(1, Math.ceil((maxPos - minPos + 1) / TRACK_SEQUENCE_MAX_RENDER_BASES));
  const step = Math.max(baseStep, densityStep);
  const pixelPerBp = (TRACK_LAYOUT.xEnd - TRACK_LAYOUT.xStart) / Math.max(1, span);

  const renderPositions = new Set();
  for (let pos = minPos; pos <= maxPos; pos += step) {
    renderPositions.add(Math.round(pos));
  }
  forcedPositions.forEach(pos => {
    const roundedPos = Math.round(Number(pos));
    if (!Number.isFinite(roundedPos)) return;
    if (roundedPos < minPos || roundedPos > maxPos) return;
    renderPositions.add(roundedPos);
  });

  const sortedPositions = Array.from(renderPositions).sort((a, b) => a - b);
  let minBpDistance = Number.POSITIVE_INFINITY;
  for (let i = 1; i < sortedPositions.length; i++) {
    const diff = sortedPositions[i] - sortedPositions[i - 1];
    if (diff > 0 && diff < minBpDistance) {
      minBpDistance = diff;
    }
  }
  const effectiveBpStep = Number.isFinite(minBpDistance) ? Math.max(1, minBpDistance) : step;
  const cellPx = pixelPerBp * effectiveBpStep;
  const showLetters = span <= TRACK_SEQUENCE_LETTER_MAX_SPAN && cellPx >= TRACK_SEQUENCE_LETTER_MIN_PX;
  const boxWidth = Math.max(2.2, Math.min(16, cellPx - (showLetters ? 1.2 : 0.8)));
  const boxHeight = showLetters ? 12 : 9;
  const yTop = y + (showLetters ? 8 : 7);

  const boxes = [];
  sortedPositions.forEach(pos => {
    const idx = pos - seqStart;
    if (idx < 0 || idx >= seq.length) return;
    const base = normalizeNtBase(seq[idx]);
    const xCenter = toX(pos);
    const x = xCenter - (boxWidth / 2);
    boxes.push(`
      <g class="track-nt-cell ${rowClass}">
        <rect x="${formatPlotDecimal(x, 2)}" y="${formatPlotDecimal(yTop, 2)}" width="${formatPlotDecimal(boxWidth, 2)}" height="${formatPlotDecimal(boxHeight, 2)}" rx="1.8" class="track-nt-box track-nt-box-${base}" style="fill:${getNtBaseColor(base)}"/>
        ${showLetters ? `<text x="${formatPlotDecimal(xCenter, 2)}" y="${formatPlotDecimal(yTop + boxHeight - 2.2, 2)}" text-anchor="middle" class="track-nt-letter">${base}</text>` : ''}
      </g>
    `);
  });

  if (boxes.length === 0) return { svg: '', showNt: false, source: '', lettersVisible: false };
  return {
    svg: boxes.join(''),
    showNt: true,
    source: 'full',
    lettersVisible: showLetters
  };
}

function buildNtBoxesFromMutations(records, posAccessor, baseAccessor, toX, y, span, rowClass) {
  if (!Array.isArray(records) || records.length === 0 || span > 420) {
    return { svg: '', showNt: false, source: '', lettersVisible: false };
  }

  const unique = new Map();
  records.forEach(record => {
    const pos = Number(posAccessor(record));
    if (!Number.isFinite(pos)) return;
    const base = normalizeNtBase(baseAccessor(record));
    const roundedPos = Math.round(pos);
    if (!unique.has(roundedPos)) {
      unique.set(roundedPos, { pos, base });
    }
  });
  const points = Array.from(unique.values())
    .sort((a, b) => a.pos - b.pos)
    .slice(0, TRACK_SEQUENCE_MAX_RENDER_BASES);
  if (points.length === 0) return { svg: '', showNt: false, source: '', lettersVisible: false };

  const pixelPerBp = (TRACK_LAYOUT.xEnd - TRACK_LAYOUT.xStart) / Math.max(1, span);
  let minBpDistance = Number.POSITIVE_INFINITY;
  for (let i = 1; i < points.length; i++) {
    const diff = points[i].pos - points[i - 1].pos;
    if (diff > 0 && diff < minBpDistance) {
      minBpDistance = diff;
    }
  }
  const effectiveBpStep = Number.isFinite(minBpDistance) ? Math.max(1, minBpDistance) : 1;
  const cellPx = pixelPerBp * effectiveBpStep;
  const showLetters = span <= TRACK_SEQUENCE_LETTER_MAX_SPAN && cellPx >= TRACK_SEQUENCE_LETTER_MIN_PX;
  const boxWidth = Math.max(2.1, Math.min(10, cellPx - (showLetters ? 1.2 : 0.8)));
  const boxHeight = showLetters ? 11 : 8.5;
  const yTop = y + (showLetters ? 8 : 7);

  const svg = points.map(point => {
    const xCenter = toX(point.pos);
    const x = xCenter - (boxWidth / 2);
    return `
      <g class="track-nt-cell ${rowClass}">
        <rect x="${formatPlotDecimal(x, 2)}" y="${formatPlotDecimal(yTop, 2)}" width="${formatPlotDecimal(boxWidth, 2)}" height="${formatPlotDecimal(boxHeight, 2)}" rx="1.8" class="track-nt-box track-nt-box-${point.base}" style="fill:${getNtBaseColor(point.base)}"/>
        ${showLetters ? `<text x="${formatPlotDecimal(xCenter, 2)}" y="${formatPlotDecimal(yTop + boxHeight - 2, 2)}" text-anchor="middle" class="track-nt-letter">${point.base}</text>` : ''}
      </g>
    `;
  }).join('');

  return {
    svg,
    showNt: true,
    source: 'mutations',
    lettersVisible: showLetters
  };
}

function buildNtOverlay({
  visibleRecords,
  toWindowX,
  toRefX,
  windowStart,
  windowEnd,
  refStart,
  refEnd,
  yGenome,
  yRef,
  fullSampleSequence = null,
  fullReferenceSequence = null
}) {
  const span = Math.max(1, windowEnd - windowStart);
  if (span > 2500) {
    return {
      grid: '',
      sampleSvg: '',
      referenceSvg: '',
      sampleShowNt: false,
      sampleSource: '',
      sampleLettersVisible: false,
      referenceShowNt: false,
      referenceSource: '',
      referenceLettersVisible: false
    };
  }

  const step = span <= 120 ? 1 : span <= 300 ? 5 : span <= 900 ? 10 : 25;
  let grid = '';
  const maxTicks = 160;
  let ticks = 0;
  for (
    let pos = Math.ceil(windowStart / step) * step;
    pos <= windowEnd && ticks < maxTicks;
    pos += step
  ) {
    const x = toWindowX(pos);
    grid += `
      <line x1="${formatPlotDecimal(x, 2)}" y1="${yGenome - 18}" x2="${formatPlotDecimal(x, 2)}" y2="${yGenome + 18}" class="track-grid-line"/>
    `;
    ticks += 1;
  }

  let sampleLayer = { svg: '', showNt: false, source: '', lettersVisible: false };
  if (fullSampleSequence && span <= TRACK_SEQUENCE_SHOW_MAX_SPAN) {
    sampleLayer = buildNtBoxesFromSequence(
      fullSampleSequence,
      toWindowX,
      windowStart,
      windowEnd,
      yGenome,
      span,
      'track-nt-row-sample',
      visibleRecords.map(record => record.genomePos)
    );
  }
  if (!sampleLayer.showNt) {
    sampleLayer = buildNtBoxesFromMutations(
      visibleRecords,
      record => record.genomePos,
      record => record.altAllele || record.refAllele || '',
      toWindowX,
      yGenome,
      span,
      'track-nt-row-sample'
    );
  }

  const refSpan = Math.max(1, refEnd - refStart);
  let referenceLayer = { svg: '', showNt: false, source: '', lettersVisible: false };
  if (fullReferenceSequence && refSpan <= TRACK_SEQUENCE_SHOW_MAX_SPAN) {
    referenceLayer = buildNtBoxesFromSequence(
      fullReferenceSequence,
      toRefX,
      refStart,
      refEnd,
      yRef,
      refSpan,
      'track-nt-row-reference',
      visibleRecords.map(record => record.refPos)
    );
  }
  if (!referenceLayer.showNt) {
    referenceLayer = buildNtBoxesFromMutations(
      visibleRecords,
      record => record.refPos,
      record => record.refAllele || record.altAllele || '',
      toRefX,
      yRef,
      refSpan,
      'track-nt-row-reference'
    );
  }

  return {
    grid,
    sampleSvg: sampleLayer.svg,
    referenceSvg: referenceLayer.svg,
    sampleShowNt: sampleLayer.showNt,
    sampleSource: sampleLayer.source,
    sampleLettersVisible: sampleLayer.lettersVisible,
    referenceShowNt: referenceLayer.showNt,
    referenceSource: referenceLayer.source,
    referenceLettersVisible: referenceLayer.lettersVisible
  };
}

function buildTrackAminoAcidAnnotations(
  sampleRecords,
  visibleRaw,
  trackState,
  toWindowX,
  yAa,
  fullSampleSequence = null
) {
  const geneContext = resolveTrackAaGeneContext(sampleRecords, visibleRaw, trackState);
  if (!geneContext) {
    return { svg: '', visibleAaCount: 0, totalAaCount: 0, geneLabel: '' };
  }

  const codonSeries = buildTrackGeneCodonSeries(geneContext);
  if (codonSeries.length === 0) {
    return { svg: '', visibleAaCount: 0, totalAaCount: 0, geneLabel: geneContext.label };
  }

  const visibleCodons = codonSeries.filter(codon =>
    codon.end >= trackState.windowStart && codon.start <= trackState.windowEnd
  );
  if (visibleCodons.length === 0) {
    return { svg: '', visibleAaCount: 0, totalAaCount: codonSeries.length, geneLabel: geneContext.label };
  }

  const mutationByAaIndex = new Map();
  sampleRecords.forEach(record => {
    if (getTrackRecordGeneLabel(record) !== geneContext.label) return;
    const aaIndex = parseTrackAaPosition(record);
    if (!aaIndex) return;
    if (mutationByAaIndex.has(aaIndex)) return;
    mutationByAaIndex.set(aaIndex, record);
  });

  const boxHeight = 7.6;
  const blockY = yAa + 4;

  const svg = visibleCodons.map(codon => {
    const codonStart = Math.max(trackState.windowStart, codon.start);
    const codonEnd = Math.min(trackState.windowEnd, codon.end);
    const x1 = toWindowX(codonStart - 0.5);
    const x2 = toWindowX(codonEnd + 0.5);
    const blockWidth = Math.max(4.2, x2 - x1);
    const centerX = x1 + (blockWidth / 2);
    const codonCellWidth = blockWidth / 3;
    const showLetter = codonCellWidth >= TRACK_SEQUENCE_LETTER_MIN_PX;

    const mutationRecord = mutationByAaIndex.get(codon.aaIndex) || null;
    const aaFromSeq = getTrackAaFromSequenceCodon(
      fullSampleSequence,
      codon.start,
      codon.end,
      geneContext.strand
    );
    const aaSymbol = aaFromSeq || getTrackAaSymbol(mutationRecord?.aaChange);
    const color = aaFromSeq ? (AA_COLORS[aaFromSeq] || AA_COLORS.X) : getTrackAaColor(mutationRecord?.aaChange);

    const tooltip = [
      `Gene ${geneContext.label}`,
      `AA ${codon.aaIndex}${aaSymbol ? ` (${aaSymbol})` : ''}`,
      `Codon ${codon.start}-${codon.end}`,
      mutationRecord?.aaChange ? `Mutation ${mutationRecord.aaChange}` : null
    ].filter(Boolean).join(' | ');

    const codonBoxes = [0, 1, 2].map(offset => {
      const cellX = x1 + (offset * codonCellWidth) + 0.16;
      const cellWidth = Math.max(1.2, codonCellWidth - 0.32);
      return `<rect x="${formatPlotDecimal(cellX, 2)}" y="${formatPlotDecimal(blockY, 2)}" width="${formatPlotDecimal(cellWidth, 2)}" height="${formatPlotDecimal(boxHeight, 2)}" rx="1.6" class="track-aa-codon-nt" style="fill:${color}"/>`;
    }).join('');

    const aaIndexText = String(codon.aaIndex);
    const minAaIndexLabelWidth = Math.max(10, aaIndexText.length * 4 + 4);
    const showAaIndexLabel = Boolean(mutationRecord) && blockWidth >= minAaIndexLabelWidth;

    return `
      <g class="track-aa-point">
        <title>${escapeSvg(tooltip)}</title>
        ${codonBoxes}
        ${showLetter && aaSymbol ? `<text x="${formatPlotDecimal(centerX, 2)}" y="${formatPlotDecimal(blockY + boxHeight - 1.7, 2)}" text-anchor="middle" class="track-aa-label">${escapeSvg(aaSymbol)}</text>` : ''}
        ${showAaIndexLabel ? `<text x="${formatPlotDecimal(centerX, 2)}" y="${formatPlotDecimal(blockY + boxHeight + 8.3, 2)}" text-anchor="middle" class="track-aa-pos-label">${escapeSvg(aaIndexText)}</text>` : ''}
      </g>
    `;
  }).join('');

  return {
    svg,
    visibleAaCount: visibleCodons.length,
    totalAaCount: codonSeries.length,
    geneLabel: geneContext.label
  };
}

function buildGenomicTrackSvg(
  sampleRecords,
  trackState,
  fullSampleSequence = null,
  fullReferenceSequence = null,
  precomputed = null
) {
  const { width, height, xStart, xEnd, yRef, yGenome, yGene, yAa, yOverview } = TRACK_LAYOUT;
  const domainSpan = Math.max(1, trackState.domainEnd - trackState.domainStart);
  const windowSpan = Math.max(1, trackState.windowEnd - trackState.windowStart);
  const visibleRaw = sampleRecords.filter(record =>
    record.genomePos >= trackState.windowStart && record.genomePos <= trackState.windowEnd
  );
  const visibleRecords = downsampleTrackPoints(visibleRaw);
  const refSource = sampleRecords;
  const refValues = refSource.map(record => record.refPos);
  const refAxisStart = 1;
  const refAxisMaxFromData = Math.max(1, ...refValues);
  const refAxisEndRaw = Number.isFinite(trackState.referenceSequenceTotalLength)
    ? Math.max(refAxisMaxFromData, Number(trackState.referenceSequenceTotalLength))
    : refAxisMaxFromData;
  const refAxisEnd = refAxisEndRaw <= refAxisStart ? refAxisStart + 1 : refAxisEndRaw;
  const refSpan = Math.max(1, refAxisEnd - refAxisStart);
  const refOverlaySource = visibleRaw.length > 0 ? visibleRaw : sampleRecords;
  const refOverlayValues = refOverlaySource.map(record => record.refPos);
  const refOverlayStart = Math.max(1, Math.floor(Math.min(...refOverlayValues)));
  const refOverlayEnd = Math.max(refOverlayStart, Math.ceil(Math.max(...refOverlayValues)));

  const toWindowX = pos =>
    xStart + ((pos - trackState.windowStart) / windowSpan) * (xEnd - xStart);
  const toDomainX = pos => xStart + ((pos - trackState.domainStart) / domainSpan) * (xEnd - xStart);
  const toRefX = pos => xStart + ((pos - refAxisStart) / refSpan) * (xEnd - xStart);
  const ntOverlay = buildNtOverlay({
    visibleRecords: visibleRaw,
    toWindowX,
    toRefX,
    windowStart: trackState.windowStart,
    windowEnd: trackState.windowEnd,
    refStart: refOverlayStart,
    refEnd: refOverlayEnd,
    yGenome,
    yRef,
    fullSampleSequence,
    fullReferenceSequence
  });
  const showReferenceGenomePosLabels = Boolean(ntOverlay.referenceShowNt);
  const referenceLabelKeys = new Set();
  if (showReferenceGenomePosLabels) {
    const selectedKey = String(trackState.selectedKey || '');
    let lastLabelX = Number.NEGATIVE_INFINITY;
    visibleRecords.forEach(record => {
      const key = String(record.trackKey || '');
      if (!key) return;
      const refX = toRefX(record.refPos);
      const isSelected = selectedKey && key === selectedKey;
      if (isSelected || (refX - lastLabelX) >= TRACK_REFERENCE_POINT_LABEL_MIN_GAP_PX) {
        referenceLabelKeys.add(key);
        lastLabelX = refX;
      }
    });
    if (selectedKey) referenceLabelKeys.add(selectedKey);
  }

  const selectedRecord = precomputed?.recordByKey instanceof Map
    ? (precomputed.recordByKey.get(trackState.selectedKey) || null)
    : (sampleRecords.find(record => record.trackKey === trackState.selectedKey) || null);
  const selectedGene = getTrackRecordGeneLabel(selectedRecord);
  const connectors = visibleRecords.map(record => {
    const refX = toRefX(record.refPos);
    const genomeX = toWindowX(record.genomePos);
    const mutationLabel = `${record.refAllele || '?'}>${record.altAllele || '?'}`;
    const isSelected = record.trackKey === trackState.selectedKey;
    const meta = [
      `Ref ${Math.round(record.refPos)}`,
      `Sample ${Math.round(record.genomePos)}`,
      mutationLabel !== '?>?' ? mutationLabel : null,
      getTrackRecordGeneLabel(record) ? `Gene ${getTrackRecordGeneLabel(record)}` : null,
      record.aaPos ? `AA pos ${record.aaPos}` : null,
      record.aaChange ? `AA ${record.aaChange}` : null,
      record.lineage ? record.lineage : null
    ].filter(Boolean).join(' | ');
    const selectedClass = isSelected ? ' is-selected' : '';
    const showPointPositionLabels = referenceLabelKeys.has(String(record.trackKey || ''));
    const referencePosLabelSvg = showPointPositionLabels
      ? `<text x="${formatPlotDecimal(refX, 2)}" y="${yRef - 24}" text-anchor="middle" class="track-point-genome-label${selectedClass}">${escapeSvg(formatTrackBp(record.refPos))}</text>`
      : '';
    const samplePosLabelY = ntOverlay.sampleShowNt ? (yGenome - 12) : (yGenome - 10);
    const samplePosLabelSvg = showPointPositionLabels
      ? `<text x="${formatPlotDecimal(genomeX, 2)}" y="${formatPlotDecimal(samplePosLabelY, 2)}" text-anchor="middle" class="track-point-sample-label${selectedClass}">${escapeSvg(formatTrackBp(record.genomePos))}</text>`
      : '';
    return `
      <g class="track-point${selectedClass}" data-track-key="${escapeSvg(record.trackKey)}">
        <title>${escapeSvg(meta)}</title>
        ${referencePosLabelSvg}
        ${samplePosLabelSvg}
        <line x1="${formatPlotDecimal(refX, 2)}" y1="${yRef}" x2="${formatPlotDecimal(genomeX, 2)}" y2="${yGenome}" class="track-link${selectedClass}"/>
        <circle cx="${formatPlotDecimal(refX, 2)}" cy="${yRef}" r="4.2" class="track-dot track-dot-ref track-mutation-hit${selectedClass}" data-track-key="${escapeSvg(record.trackKey)}" data-genome-pos="${formatPlotDecimal(record.genomePos, 3)}"/>
        <circle cx="${formatPlotDecimal(genomeX, 2)}" cy="${yGenome}" r="4.2" class="track-dot track-dot-genome track-mutation-hit${selectedClass}" data-track-key="${escapeSvg(record.trackKey)}" data-genome-pos="${formatPlotDecimal(record.genomePos, 3)}"/>
      </g>
    `;
  }).join('');

  const segmentData = (
    precomputed &&
    Array.isArray(precomputed.geneSegments) &&
    typeof precomputed.geneMode === 'string'
  )
    ? { segments: precomputed.geneSegments, mode: precomputed.geneMode }
    : buildTrackGeneSegments(sampleRecords);
  const { segments, mode } = segmentData;
  const visibleGenes = segments
    .filter(segment => segment.max >= trackState.windowStart && segment.min <= trackState.windowEnd)
    .slice(0, 30);
  const geneSegments = visibleGenes.map((segment, index) => {
    const clippedStart = Math.max(trackState.windowStart, segment.min);
    const clippedEnd = Math.min(trackState.windowEnd, segment.max);
    const x1 = toWindowX(clippedStart);
    const x2 = toWindowX(clippedEnd);
    const barWidth = Math.max(5, x2 - x1);
    const row = index % 3;
    const y = yGene - 18 + (row * 12);
    const activeClass = selectedGene && selectedGene === segment.gene ? ' is-active' : '';
    const strand = segment.strand === '+' || segment.strand === '-' ? segment.strand : '';
    const strandArrow = strand === '+' ? '→' : strand === '-' ? '←' : '';
    const strandX = strand === '+'
      ? x1 + Math.max(6, barWidth - 5)
      : x1 + 5;
    const strandAnchor = strand === '+' ? 'end' : 'start';
    const geneLabelText = strand ? `${segment.gene} (${strand})` : segment.gene;
    const titleText = strand
      ? `${segment.gene}: ${segment.count} mutation(s) · strand ${strand}`
      : `${segment.gene}: ${segment.count} mutation(s)`;
    return `
      <g class="track-gene">
        <title>${escapeSvg(titleText)}<\/title>
        <rect x="${formatPlotDecimal(x1, 2)}" y="${formatPlotDecimal(y, 2)}" width="${formatPlotDecimal(barWidth, 2)}" height="9.5" rx="4" class="track-gene-box track-gene-hit${activeClass}" data-gene-start="${formatPlotDecimal(segment.min, 0)}" data-gene-end="${formatPlotDecimal(segment.max, 0)}" data-gene-label="${escapeSvg(segment.gene)}" data-gene-strand="${escapeSvg(strand)}"/>
        ${strandArrow && barWidth > 12 ? `<text x="${formatPlotDecimal(strandX, 2)}" y="${formatPlotDecimal(y + 7.2, 2)}" text-anchor="${strandAnchor}" class="track-gene-strand">${strandArrow}</text>` : ''}
        ${barWidth > 26 ? `<text x="${formatPlotDecimal(x1 + (barWidth / 2), 2)}" y="${formatPlotDecimal(y - 2, 2)}" text-anchor="middle" class="track-gene-label">${escapeSvg(geneLabelText)}</text>` : ''}
      </g>
    `;
  }).join('');
  const aaAnnotations = buildTrackAminoAcidAnnotations(
    sampleRecords,
    visibleRaw,
    trackState,
    toWindowX,
    yAa,
    fullSampleSequence
  );

  const overviewSource = Array.isArray(precomputed?.overviewRecords)
    ? precomputed.overviewRecords
    : downsampleTrackPoints(sampleRecords);
  const overviewTicks = overviewSource.map(record => {
    const x = toDomainX(record.genomePos);
    return `<line x1="${formatPlotDecimal(x, 2)}" y1="${yOverview - 6}" x2="${formatPlotDecimal(x, 2)}" y2="${yOverview + 6}" class="track-overview-tick"/>`;
  }).join('');
  const overviewStart = toDomainX(trackState.windowStart);
  const overviewEnd = toDomainX(trackState.windowEnd);
  const overviewWidth = Math.max(6, overviewEnd - overviewStart);
  const sampleWindowTickY = ntOverlay.sampleShowNt ? (yGenome - 30) : (yGenome - 22);

  const svg = `
    <svg class="viz-track-svg" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" role="img" aria-label="Interactive genomic browser track">
      <line x1="${xStart}" y1="${yRef}" x2="${xEnd}" y2="${yRef}" class="track-axis"/>
      <line x1="${xStart}" y1="${yGenome}" x2="${xEnd}" y2="${yGenome}" class="track-axis"/>
      <line x1="${xStart}" y1="${yGene}" x2="${xEnd}" y2="${yGene}" class="track-axis track-axis-gene"/>
      <line x1="${xStart}" y1="${yAa}" x2="${xEnd}" y2="${yAa}" class="track-axis track-axis-aa"/>
      <text x="20" y="${yRef + 4}" class="track-label">Reference (mapped)</text>
      <text x="20" y="${yGenome + 4}" class="track-label">Sample Genome</text>
      <text x="20" y="${yGene + 4}" class="track-label">Genes (Sample GFF)</text>
      <text x="20" y="${yAa + 4}" class="track-label track-label-aa">Amino acids (gene)</text>
      <text x="${xStart}" y="${yRef - 12}" class="track-tick">${formatTrackBp(refAxisStart)}</text>
      <text x="${xEnd}" y="${yRef - 12}" text-anchor="end" class="track-tick">${formatTrackBp(refAxisEnd)}</text>
      <text x="${xStart}" y="${sampleWindowTickY}" class="track-tick">${formatTrackBp(trackState.windowStart)}</text>
      <text x="${xEnd}" y="${sampleWindowTickY}" text-anchor="end" class="track-tick">${formatTrackBp(trackState.windowEnd)}</text>
      ${ntOverlay.grid}
      ${connectors}
      ${ntOverlay.referenceSvg}
      ${ntOverlay.sampleSvg}
      ${geneSegments}
      ${aaAnnotations.svg}
      <line x1="${xStart}" y1="${yOverview}" x2="${xEnd}" y2="${yOverview}" class="track-overview-axis"/>
      ${overviewTicks}
      <rect x="${xStart}" y="${yOverview - 11}" width="${xEnd - xStart}" height="22" class="track-overview-hit"/>
      <rect x="${formatPlotDecimal(overviewStart, 2)}" y="${yOverview - 8}" width="${formatPlotDecimal(overviewWidth, 2)}" height="16" rx="4" class="track-overview-window"/>
      <text x="20" y="${yOverview + 3}" class="track-label track-overview-label">Overview</text>
      <text x="${xStart}" y="${yOverview + 22}" class="track-tick">${formatTrackBp(trackState.domainStart)}</text>
      <text x="${xEnd}" y="${yOverview + 22}" text-anchor="end" class="track-tick">${formatTrackBp(trackState.domainEnd)}</text>
    </svg>
  `;

  return {
    svg,
    totalMutations: sampleRecords.length,
    inWindowMutations: visibleRaw.length,
    renderedMutations: visibleRecords.length,
    visibleGeneCount: visibleGenes.length,
    visibleAaCount: aaAnnotations.visibleAaCount,
    totalAaCount: aaAnnotations.totalAaCount,
    aaGeneLabel: aaAnnotations.geneLabel,
    geneMode: mode,
    sampleShowNt: ntOverlay.sampleShowNt,
    sampleNtSource: ntOverlay.sampleSource,
    sampleNtLettersVisible: ntOverlay.sampleLettersVisible,
    referenceShowNt: ntOverlay.referenceShowNt,
    referenceNtSource: ntOverlay.referenceSource,
    referenceNtLettersVisible: ntOverlay.referenceLettersVisible,
    selectedRecord
  };
}

function buildTrackDetailHtml(trackInfo, trackState) {
  if (!trackInfo.selectedRecord) {
    return `
      <div class="viz-track-detail">
        <span class="viz-track-detail-title">Selection</span>
        <span>Click a mutation point to inspect details and auto-zoom. Click a gene block to jump to that region.</span>
      </div>
    `;
  }
  const record = trackInfo.selectedRecord;
  const mutation = `${record.refAllele || '?'}>${record.altAllele || '?'}`;
  const sampleSequenceData = (
    trackState?.sequence &&
    Number.isFinite(trackState?.sequenceStart) &&
    Number.isFinite(trackState?.sequenceEnd)
  ) ? {
      start: trackState.sequenceStart,
      end: trackState.sequenceEnd,
      sequence: trackState.sequence
    } : null;
  const referenceSequenceData = (
    trackState?.referenceSequence &&
    Number.isFinite(trackState?.referenceSequenceStart) &&
    Number.isFinite(trackState?.referenceSequenceEnd)
  ) ? {
      start: trackState.referenceSequenceStart,
      end: trackState.referenceSequenceEnd,
      sequence: trackState.referenceSequence
    } : null;
  const observedSampleNt = getTrackSequenceBaseAtPosition(sampleSequenceData, record.genomePos);
  const observedReferenceNt = getTrackSequenceBaseAtPosition(referenceSequenceData, record.refPos);
  const expectedSampleNt = getTrackExpectedSampleAllele(record);
  const expectedReferenceNt = String(record.refAllele || '').trim().toUpperCase();
  const diagnostics = [];
  if (observedSampleNt) {
    diagnostics.push(
      expectedSampleNt
        ? `Sample FASTA base ${observedSampleNt} vs expected ${expectedSampleNt}${observedSampleNt === expectedSampleNt ? ' (match)' : ' (mismatch)'}`
        : `Sample FASTA base ${observedSampleNt}`
    );
  }
  if (observedReferenceNt) {
    diagnostics.push(
      isTrackInformativeNt(expectedReferenceNt)
        ? `Reference FASTA base ${observedReferenceNt} vs expected ${expectedReferenceNt}${observedReferenceNt === expectedReferenceNt ? ' (match)' : ' (mismatch)'}`
        : `Reference FASTA base ${observedReferenceNt}`
    );
  }
  if (record.altAlleleSource === 'kmer') {
    diagnostics.push('Expected ALT was inferred from marker k-mer center.');
  }
  const details = [
    `Ref ${formatTrackBp(record.refPos)}`,
    `Sample ${formatTrackBp(record.genomePos)}`,
    mutation !== '?>?' ? mutation : null,
    getTrackRecordGeneLabel(record) ? `Gene ${getTrackRecordGeneLabel(record)}` : null,
    getTrackRecordGeneStrand(record) ? `Strand ${getTrackRecordGeneStrand(record)}` : null,
    record.aaPos ? `AA pos ${record.aaPos}` : null,
    record.aaChange ? `AA ${record.aaChange}` : null,
    record.lineage ? `Lineage ${record.lineage}` : null
  ].filter(Boolean).join(' · ');
  return `
    <div class="viz-track-detail">
      <span class="viz-track-detail-title">Selected Mutation</span>
      <span>${escapeHtml(details)}</span>
      ${diagnostics.length > 0 ? `<span class="viz-track-detail-window">${escapeHtml(diagnostics.join(' · '))}</span>` : ''}
      <span class="viz-track-detail-window">Sample window ${formatTrackBp(trackState.windowStart)} - ${formatTrackBp(trackState.windowEnd)}</span>
    </div>
  `;
}

function buildTrackMutationSearchOptions(sampleRecords, maxOptions = 1200) {
  if (!Array.isArray(sampleRecords) || sampleRecords.length === 0) return '';
  const byMutationKey = new Map();
  for (const record of sampleRecords) {
    const genomePos = Math.round(record.genomePos);
    if (!Number.isFinite(genomePos)) continue;
    const mutation = `${record.refAllele || '?'}>${record.altAllele || '?'}`;
    const geneLabel = getTrackRecordGeneLabel(record) || 'intergenic';
    const optionKey = `${genomePos}|${geneLabel}|${mutation}`;
    if (byMutationKey.has(optionKey)) continue;
    byMutationKey.set(optionKey, record);
    if (byMutationKey.size >= maxOptions) break;
  }
  return Array.from(byMutationKey.values())
    .sort((a, b) => a.genomePos - b.genomePos)
    .map(record => {
      const genomePos = Math.round(record.genomePos);
      const refPos = Math.round(record.refPos);
      const mutation = `${record.refAllele || '?'}>${record.altAllele || '?'}`;
      const geneLabel = getTrackRecordGeneLabel(record) || 'intergenic';
      const value = `${genomePos} ${geneLabel} ${mutation}`;
      const desc = `sample ${genomePos} · ref ${refPos} · ${geneLabel} · ${mutation}`;
      return `<option value="${escapeHtml(value)}" label="${escapeHtml(desc)}">${escapeHtml(desc)}</option>`;
    })
    .join('');
}

function buildTrackMutationSelectOptions(sampleRecords, maxOptions = 2500) {
  if (!Array.isArray(sampleRecords) || sampleRecords.length === 0) return '';
  const unique = new Map();
  for (const record of sampleRecords) {
    if (!record?.trackKey) continue;
    if (unique.has(record.trackKey)) continue;
    unique.set(record.trackKey, record);
    if (unique.size >= maxOptions) break;
  }
  return Array.from(unique.values())
    .sort((a, b) => a.genomePos - b.genomePos || a.refPos - b.refPos)
    .map(record => {
      const genomePos = Math.round(record.genomePos);
      const refPos = Math.round(record.refPos);
      const mutation = `${record.refAllele || '?'}>${record.altAllele || '?'}`;
      const geneLabel = getTrackRecordGeneLabel(record) || 'intergenic';
      const text = `sample ${genomePos} · ref ${refPos} · ${geneLabel} · ${mutation}`;
      return `<option value="${escapeHtml(record.trackKey)}">${escapeHtml(text)}</option>`;
    })
    .join('');
}

function findTrackRecordsByQuery(sampleRecords, query, maxMatches = 2500) {
  const text = String(query || '').trim();
  if (!text || !Array.isArray(sampleRecords) || sampleRecords.length === 0) return [];
  const numeric = parseNumericValue(text);
  if (numeric !== null) {
    const rounded = Math.round(numeric);
    return sampleRecords.filter(record =>
      Math.round(record.genomePos) === rounded || Math.round(record.refPos) === rounded
    )
      .sort((a, b) => a.genomePos - b.genomePos || a.refPos - b.refPos)
      .slice(0, maxMatches);
  }

  const tokens = text.toLowerCase().split(/\s+/).filter(Boolean);
  if (tokens.length === 0) return [];

  const scored = [];
  sampleRecords.forEach(record => {
    const mutation = `${record.refAllele || '?'}>${record.altAllele || '?'}`.toLowerCase();
    const searchable = [
      String(Math.round(record.genomePos)),
      String(Math.round(record.refPos)),
      normalizeValue(getTrackRecordGeneLabel(record)).toLowerCase(),
      normalizeValue(record.aaPos).toLowerCase(),
      normalizeValue(record.aaChange).toLowerCase(),
      normalizeValue(record.lineage).toLowerCase(),
      mutation
    ];
    let score = 0;
    tokens.forEach(token => {
      if (!token) return;
      if (searchable.some(field => field === token)) score += 14;
      else if (searchable.some(field => field.startsWith(token))) score += 8;
      else if (searchable.some(field => field.includes(token))) score += 4;
    });
    if (score > 0) {
      scored.push({ record, score });
    }
  });

  scored.sort((a, b) =>
    b.score - a.score ||
    a.record.genomePos - b.record.genomePos ||
    a.record.refPos - b.record.refPos
  );
  return scored.slice(0, maxMatches).map(entry => entry.record);
}

function focusTrackRecord(
  toolId,
  data,
  record,
  zoomFactor = TRACK_CLICK_FOCUS_ZOOM_FACTOR,
  searchQuery = '',
  searchMeta = null
) {
  const current = TRACK_STATE[toolId];
  if (!current || !record) return false;
  const focusPos = Number(record.genomePos);
  const nextWindow = Number.isFinite(focusPos)
    ? zoomTrackWindow(
        current.windowStart,
        current.windowEnd,
        zoomFactor,
        current.domainStart,
        current.domainEnd,
        focusPos
      )
    : { windowStart: current.windowStart, windowEnd: current.windowEnd };
  const nextState = {
    ...current,
    selectedKey: record.trackKey,
    aaFocusGene: getTrackRecordGeneLabel(record) || current.aaFocusGene || '',
    searchQuery: String(searchQuery || current.searchQuery || ''),
    autoFollowDomain: false,
    ...nextWindow
  };
  if (searchMeta && Array.isArray(searchMeta.searchMatchKeys)) {
    const matchKeys = searchMeta.searchMatchKeys.filter(key => typeof key === 'string' && key.length > 0);
    const rawIndex = Number(searchMeta.searchMatchIndex);
    const clampedIndex = matchKeys.length === 0
      ? -1
      : Math.min(
          Math.max(Number.isFinite(rawIndex) ? Math.trunc(rawIndex) : 0, 0),
          matchKeys.length - 1
        );
    nextState.searchMatchKeys = matchKeys;
    nextState.searchMatchIndex = clampedIndex;
  }
  TRACK_STATE[toolId] = nextState;
  renderGenomicTrack(toolId, data);
  emitTrackRecordSelectionEvent(toolId, record);
  return true;
}

export function focusGenomicTrackRecordBySourceRow(toolId, data, sourceRowIndex) {
  const targetSourceRow = normalizeTrackSourceRowIndex(sourceRowIndex);
  if (!Number.isFinite(targetSourceRow)) return false;

  const { records, reason } = getCachedMutationTrackData(data);
  if (records.length === 0) {
    if (reason) {
      logMessage(reason, 'warning');
    }
    return false;
  }

  const record = records.find(entry => normalizeTrackSourceRowIndex(entry?.sourceRowIndex) === targetSourceRow);
  if (!record) return false;

  const sample = String(record.sample || '').trim() || 'All samples';
  const sampleRecords = records.filter(entry => String(entry?.sample || '').trim() === sample);
  if (sampleRecords.length === 0) return false;

  const nextState = ensureTrackState(toolId, sample, sampleRecords);
  TRACK_STATE[toolId] = {
    ...nextState,
    searchMatchKeys: [],
    searchMatchIndex: -1,
    autoFollowDomain: false
  };

  return focusTrackRecord(toolId, data, record, TRACK_CLICK_FOCUS_ZOOM_FACTOR);
}

function hasSameTrackSearchMatches(trackState, nextKeys) {
  if (!trackState || !Array.isArray(trackState.searchMatchKeys)) return false;
  if (trackState.searchMatchKeys.length !== nextKeys.length) return false;
  for (let i = 0; i < nextKeys.length; i++) {
    if (trackState.searchMatchKeys[i] !== nextKeys[i]) return false;
  }
  return true;
}

function applyTrackMutationSearch(toolId, data, sampleRecords, query, direction = 'first') {
  const text = String(query || '').trim();
  const current = TRACK_STATE[toolId];
  if (!current || !text) return false;

  const matches = findTrackRecordsByQuery(sampleRecords, text);
  if (matches.length === 0) {
    TRACK_STATE[toolId] = {
      ...current,
      searchQuery: text,
      searchMatchKeys: [],
      searchMatchIndex: -1
    };
    renderGenomicTrack(toolId, data);
    logMessage(`No mutation match found for "${text}" in this sample.`, 'warning');
    return false;
  }

  const matchKeys = matches.map(record => record.trackKey);
  const sameQuery = String(current.searchQuery || '').trim().toLowerCase() === text.toLowerCase();
  const sameMatches = hasSameTrackSearchMatches(current, matchKeys) && sameQuery;
  const currentIndex = Number(current.searchMatchIndex);
  let nextIndex = 0;
  if (direction === 'next' && sameMatches && Number.isFinite(currentIndex) && currentIndex >= 0) {
    nextIndex = (Math.trunc(currentIndex) + 1) % matchKeys.length;
  } else if (direction === 'prev' && sameMatches && Number.isFinite(currentIndex) && currentIndex >= 0) {
    nextIndex = (Math.trunc(currentIndex) - 1 + matchKeys.length) % matchKeys.length;
  }

  return focusTrackRecord(
    toolId,
    data,
    matches[nextIndex],
    TRACK_SEARCH_FOCUS_ZOOM_FACTOR,
    text,
    {
      searchMatchKeys: matchKeys,
      searchMatchIndex: nextIndex
    }
  );
}

function getTrackSearchStatusText(trackState) {
  const query = String(trackState?.searchQuery || '').trim();
  const matchCount = Array.isArray(trackState?.searchMatchKeys) ? trackState.searchMatchKeys.length : 0;
  if (!query) return 'No search';
  if (matchCount <= 0) return 'No matches';
  const rawIndex = Number(trackState?.searchMatchIndex);
  const currentIndex = Number.isFinite(rawIndex) ? Math.trunc(rawIndex) : 0;
  const clampedIndex = Math.min(Math.max(currentIndex, 0), matchCount - 1);
  return `Match ${clampedIndex + 1}/${matchCount}`;
}

function attachTrackInteractions(toolId, data, sampleRecords = []) {
  const controlsEl = document.getElementById(`${toolId}-viz-track-controls`);
  const bodyEl = document.getElementById(`${toolId}-viz-track-body`);
  if (!controlsEl || !bodyEl) return;
  let ignoreClickUntil = 0;
  let dragState = null;
  const recordByKey = new Map(sampleRecords.map(record => [record.trackKey, record]));
  if (controlsEl.dataset.trackControlsBound !== '1') {
    controlsEl.dataset.trackControlsBound = '1';

    const sampleSelect = controlsEl.querySelector('.viz-track-sample-select');
    sampleSelect?.addEventListener('change', () => {
      const nextSample = sampleSelect.value;
      TRACK_STATE[toolId] = { sample: nextSample };
      renderGenomicTrack(toolId, data);
    });

    const mutationSelect = controlsEl.querySelector('.viz-track-mutation-select');
    mutationSelect?.addEventListener('change', () => {
      const key = String(mutationSelect.value || '').trim();
      if (!key) return;
      const match = sampleRecords.find(record => record.trackKey === key);
      if (!match) return;
      focusTrackRecord(toolId, data, match, TRACK_SEARCH_FOCUS_ZOOM_FACTOR);
    });

    const gffFieldSelect = controlsEl.querySelector('.viz-track-gff-field-select');
    gffFieldSelect?.addEventListener('change', () => {
      const current = TRACK_STATE[toolId];
      if (!current) return;
      const selectedField = String(gffFieldSelect.value || '').trim();
      TRACK_STATE[toolId] = {
        ...current,
        gffLabelField: selectedField === '__auto__' ? '' : selectedField,
        searchMatchKeys: [],
        searchMatchIndex: -1
      };
      renderGenomicTrack(toolId, data);
    });

    const searchInput = controlsEl.querySelector('.viz-track-search-input');
    const runMutationSearch = (direction = 'first') => {
      const fallbackQuery = TRACK_STATE[toolId]?.searchQuery || '';
      const query = String(searchInput?.value || fallbackQuery).trim();
      if (!query) return;
      applyTrackMutationSearch(toolId, data, sampleRecords, query, direction);
    };
    searchInput?.addEventListener('input', () => {
      const current = TRACK_STATE[toolId];
      if (!current) return;
      TRACK_STATE[toolId] = {
        ...current,
        searchQuery: String(searchInput.value || '').trim(),
        searchMatchKeys: [],
        searchMatchIndex: -1
      };
    });
    searchInput?.addEventListener('keydown', event => {
      if (event.key !== 'Enter') return;
      event.preventDefault();
      const current = TRACK_STATE[toolId];
      const query = String(searchInput?.value || '').trim();
      const canCycle = Boolean(
        current &&
        query &&
        query.toLowerCase() === String(current.searchQuery || '').trim().toLowerCase() &&
        Array.isArray(current.searchMatchKeys) &&
        current.searchMatchKeys.length > 0
      );
      runMutationSearch(canCycle ? 'next' : 'first');
    });

    controlsEl.querySelectorAll('.viz-track-btn[data-track-action]').forEach(btn => {
      btn.addEventListener('click', () => {
        const action = btn.dataset.trackAction;
        if (action === 'toggle-expand') {
          setTrackExpanded(toolId, !isTrackExpanded(toolId));
          renderGenomicTrack(toolId, data);
          return;
        }
        if (action === 'find-mutation') {
          runMutationSearch('first');
          return;
        }
        if (action === 'search-prev') {
          runMutationSearch('prev');
          return;
        }
        if (action === 'search-next') {
          runMutationSearch('next');
          return;
        }
        const current = TRACK_STATE[toolId];
        if (!current) return;
        const nextWindow = getNextWindowFromAction(action, current, controlsEl);
        if (!nextWindow) return;
        TRACK_STATE[toolId] = {
          ...current,
          ...nextWindow,
          autoFollowDomain: action === 'reset'
        };
        renderGenomicTrack(toolId, data);
      });
    });

    controlsEl.querySelectorAll('.viz-track-window-start-input, .viz-track-window-end-input').forEach(input => {
      input.addEventListener('keydown', (event) => {
        if (event.key !== 'Enter') return;
        event.preventDefault();
        const current = TRACK_STATE[toolId];
        if (!current) return;
        const nextWindow = getNextWindowFromAction('apply-range', current, controlsEl);
        if (!nextWindow) return;
        TRACK_STATE[toolId] = { ...current, ...nextWindow, autoFollowDomain: false };
        renderGenomicTrack(toolId, data);
      });
    });
  }

  const svg = bodyEl.querySelector('.viz-track-svg');
  if (!svg) return;
  if (svg.dataset.trackInteractionsBound === '1') return;
  svg.dataset.trackInteractionsBound = '1';
  const wheelState = {
    raf: 0,
    factor: 1,
    anchor: null
  };

  const stopDrag = () => {
    window.removeEventListener('mousemove', onDragMove);
    document.body.classList.remove('track-dragging');
    if (dragState?.raf) cancelAnimationFrame(dragState.raf);
    const moved = Boolean(dragState?.moved);
    dragState = null;
    if (moved) {
      ignoreClickUntil = Date.now() + 180;
      renderGenomicTrack(toolId, data);
    }
  };

  const flushWheelZoom = () => {
    wheelState.raf = 0;
    const current = TRACK_STATE[toolId];
    if (!current) {
      wheelState.factor = 1;
      wheelState.anchor = null;
      return;
    }
    const factor = Math.min(6, Math.max(0.12, wheelState.factor));
    if (!Number.isFinite(factor) || factor === 1) {
      wheelState.factor = 1;
      wheelState.anchor = null;
      return;
    }
    const anchor = Number.isFinite(wheelState.anchor)
      ? wheelState.anchor
      : (current.windowStart + current.windowEnd) / 2;
    const nextWindow = zoomTrackWindow(
      current.windowStart,
      current.windowEnd,
      factor,
      current.domainStart,
      current.domainEnd,
      anchor
    );
    TRACK_STATE[toolId] = { ...current, ...nextWindow, autoFollowDomain: false };
    wheelState.factor = 1;
    wheelState.anchor = null;
    renderGenomicTrack(toolId, data);
  };

  const onDragMove = (event) => {
    if (!dragState) return;
    const current = TRACK_STATE[toolId];
    if (!current) return;
    const deltaPx = event.clientX - dragState.startClientX;
    if (Math.abs(deltaPx) > 2) dragState.moved = true;
    let nextWindow;
    if (dragState.mode === 'overview') {
      const domainSpan = Math.max(1, dragState.initialDomainEnd - dragState.initialDomainStart);
      const deltaBp = (deltaPx / dragState.plotWidthPx) * domainSpan;
      nextWindow = clampTrackWindow(
        dragState.initialWindowStart + deltaBp,
        dragState.initialWindowEnd + deltaBp,
        dragState.initialDomainStart,
        dragState.initialDomainEnd
      );
    } else {
      const windowSpan = Math.max(1, dragState.initialWindowEnd - dragState.initialWindowStart);
      const deltaBp = -(deltaPx / dragState.plotWidthPx) * windowSpan;
      nextWindow = clampTrackWindow(
        dragState.initialWindowStart + deltaBp,
        dragState.initialWindowEnd + deltaBp,
        dragState.initialDomainStart,
        dragState.initialDomainEnd
      );
    }
    TRACK_STATE[toolId] = { ...current, ...nextWindow, autoFollowDomain: false };
    if (!dragState.raf) {
      dragState.raf = requestAnimationFrame(() => {
        if (!dragState) return;
        dragState.raf = 0;
        renderGenomicTrack(toolId, data);
      });
    }
  };

  svg.addEventListener('mousedown', (event) => {
    if (event.button !== 0) return;
    const current = TRACK_STATE[toolId];
    if (!current) return;
    const target = event.target;
    if (target.closest('.track-mutation-hit') || target.closest('.track-gene-hit') || target.closest('.track-overview-hit')) {
      return;
    }
    const rect = svg.getBoundingClientRect();
    if (!rect.width || !rect.height) return;
    const svgX = ((event.clientX - rect.left) / rect.width) * TRACK_LAYOUT.width;
    const svgY = ((event.clientY - rect.top) / rect.height) * TRACK_LAYOUT.height;
    const plotWidthPx = rect.width * ((TRACK_LAYOUT.xEnd - TRACK_LAYOUT.xStart) / TRACK_LAYOUT.width);
    if (!plotWidthPx) return;

    const onOverviewWindow = Boolean(target.closest('.track-overview-window'));
    const onMainTrack = (
      svgX >= TRACK_LAYOUT.xStart &&
      svgX <= TRACK_LAYOUT.xEnd &&
      svgY >= TRACK_LAYOUT.yRef - 24 &&
      svgY <= TRACK_LAYOUT.yAa + 44
    );
    if (!onOverviewWindow && !onMainTrack) return;

    dragState = {
      mode: onOverviewWindow ? 'overview' : 'main',
      startClientX: event.clientX,
      plotWidthPx,
      initialWindowStart: current.windowStart,
      initialWindowEnd: current.windowEnd,
      initialDomainStart: current.domainStart,
      initialDomainEnd: current.domainEnd,
      moved: false,
      raf: 0
    };
    document.body.classList.add('track-dragging');
    window.addEventListener('mousemove', onDragMove);
    window.addEventListener('mouseup', stopDrag, { once: true });
    event.preventDefault();
  });

  svg.addEventListener('wheel', (event) => {
    event.preventDefault();
    const current = TRACK_STATE[toolId];
    if (!current) return;
    const rect = svg.getBoundingClientRect();
    if (!rect.width) return;
    const svgX = ((event.clientX - rect.left) / rect.width) * TRACK_LAYOUT.width;
    if (svgX < TRACK_LAYOUT.xStart || svgX > TRACK_LAYOUT.xEnd) return;
    const ratio = (svgX - TRACK_LAYOUT.xStart) / (TRACK_LAYOUT.xEnd - TRACK_LAYOUT.xStart);
    const anchor = current.windowStart + (ratio * (current.windowEnd - current.windowStart));
    const eventFactor = event.deltaY < 0 ? 0.78 : 1.32;
    wheelState.factor *= eventFactor;
    wheelState.anchor = anchor;
    if (!wheelState.raf) {
      wheelState.raf = requestAnimationFrame(flushWheelZoom);
    }
  }, { passive: false });

  svg.addEventListener('click', (event) => {
    if (Date.now() < ignoreClickUntil) return;
    const mutationHit = event.target.closest('.track-mutation-hit');
    if (mutationHit?.dataset.trackKey) {
      const key = mutationHit.dataset.trackKey;
      const match = recordByKey.get(key);
      if (match) {
        focusTrackRecord(toolId, data, match, TRACK_CLICK_FOCUS_ZOOM_FACTOR);
      }
      return;
    }

    const geneHit = event.target.closest('.track-gene-hit');
    if (geneHit?.dataset.geneStart && geneHit.dataset.geneEnd) {
      const current = TRACK_STATE[toolId];
      if (!current) return;
      const start = Number.parseFloat(geneHit.dataset.geneStart);
      const end = Number.parseFloat(geneHit.dataset.geneEnd);
      const geneLabel = normalizeValue(geneHit.dataset.geneLabel);
      if (!Number.isFinite(start) || !Number.isFinite(end)) return;
      const margin = Math.max(40, Math.abs(end - start) * 0.25);
      const nextWindow = clampTrackWindow(start - margin, end + margin, current.domainStart, current.domainEnd);
      TRACK_STATE[toolId] = {
        ...current,
        ...nextWindow,
        aaFocusGene: geneLabel || current.aaFocusGene || '',
        autoFollowDomain: false
      };
      renderGenomicTrack(toolId, data);
      return;
    }

    const overviewHit = event.target.closest('.track-overview-hit');
    if (overviewHit) {
      const current = TRACK_STATE[toolId];
      if (!current) return;
      const rect = svg.getBoundingClientRect();
      if (!rect.width) return;
      const svgX = ((event.clientX - rect.left) / rect.width) * TRACK_LAYOUT.width;
      const ratio = (svgX - TRACK_LAYOUT.xStart) / (TRACK_LAYOUT.xEnd - TRACK_LAYOUT.xStart);
      const clampedRatio = Math.min(Math.max(ratio, 0), 1);
      const center = current.domainStart + (clampedRatio * (current.domainEnd - current.domainStart));
      const span = current.windowEnd - current.windowStart;
      const nextWindow = clampTrackWindow(center - (span / 2), center + (span / 2), current.domainStart, current.domainEnd);
      TRACK_STATE[toolId] = { ...current, ...nextWindow, autoFollowDomain: false };
      renderGenomicTrack(toolId, data);
    }
  });
}

function renderGenomicTrack(toolId, data) {
  if (!GENOMIC_TRACK_TOOLS.has(toolId)) {
    clearGenomicTrack(toolId);
    return;
  }

  const section = document.getElementById(`${toolId}-viz-genomic-track`);
  const bodyEl = document.getElementById(`${toolId}-viz-track-body`);
  const controlsEl = document.getElementById(`${toolId}-viz-track-controls`);
  if (!section || !bodyEl || !controlsEl) return;

  const renderCache = getTrackRenderCache(toolId);
  const { records, samples, recordsBySample, reason } = getCachedMutationTrackData(data);
  if (records.length === 0) {
    section.classList.remove('hidden');
    controlsEl.replaceChildren();
    controlsEl.dataset.trackControlsBound = '';
    renderCache.controlsKey = '';
    renderCache.svgDomKey = '';
    bodyEl.innerHTML = `<p class="viz-track-empty">${escapeHtml(reason || 'No mutation track data available.')}</p>`;
    return;
  }

  const selectedSample = TRACK_STATE[toolId]?.sample && samples.includes(TRACK_STATE[toolId].sample)
    ? TRACK_STATE[toolId].sample
    : samples[0];
  const sampleRecords = recordsBySample.get(selectedSample) || [];
  if (sampleRecords.length === 0) {
    section.classList.remove('hidden');
    controlsEl.replaceChildren();
    controlsEl.dataset.trackControlsBound = '';
    renderCache.controlsKey = '';
    renderCache.svgDomKey = '';
    bodyEl.innerHTML = '<p class="viz-track-empty">No mutation points found for the selected sample.</p>';
    return;
  }

  let trackState = ensureTrackState(toolId, selectedSample, sampleRecords);
  const sequenceSource = resolveTrackSequenceSource(data, selectedSample);
  const referenceSequenceSource = resolveReferenceTrackSequenceSource(data);
  const gffSource = resolveTrackGffSource(data, selectedSample);
  if (!gffSource && (
    trackState.gffSourceKey ||
    trackState.gffLoading ||
    trackState.gffError ||
    (Array.isArray(trackState.gffFeatures) && trackState.gffFeatures.length > 0)
  )) {
    TRACK_STATE[toolId] = {
      ...trackState,
      gffFeatures: [],
      gffAttributeKeys: [],
      gffSourceKey: '',
      gffLoading: false,
      gffRequestKey: '',
      gffError: '',
      gffLabelField: ''
    };
    trackState = TRACK_STATE[toolId];
  }
  if (gffSource && trackState.gffSourceKey && trackState.gffSourceKey !== gffSource.sourceKey) {
    TRACK_STATE[toolId] = {
      ...trackState,
      gffFeatures: [],
      gffAttributeKeys: [],
      gffSourceKey: '',
      gffLoading: false,
      gffRequestKey: '',
      gffError: ''
    };
    trackState = TRACK_STATE[toolId];
  }

  const sampleView = getCachedTrackSampleView(toolId, selectedSample, sampleRecords, trackState);
  const sampleRecordsForTrack = sampleView.displayRecords;
  const windowSpan = Math.max(1, trackState.windowEnd - trackState.windowStart);
  const refWindowRecords = sampleRecordsForTrack.filter(record =>
    record.genomePos >= trackState.windowStart && record.genomePos <= trackState.windowEnd
  );
  const refRangeSource = refWindowRecords.length > 0 ? refWindowRecords : sampleRecordsForTrack;
  const refValues = refRangeSource.map(record => record.refPos);
  const refWindowStart = Math.max(1, Math.floor(Math.min(...refValues)));
  const refWindowEnd = Math.ceil(Math.max(...refValues));
  const refWindowSpan = Math.max(1, refWindowEnd - refWindowStart);

  const hasFullSequence = hasTrackSequenceForWindow(
    trackState,
    sequenceSource,
    trackState.windowStart,
    trackState.windowEnd
  );
  const hasFullReferenceSequence = hasReferenceTrackSequenceForWindow(
    trackState,
    referenceSequenceSource,
    refWindowStart,
    refWindowEnd
  );
  const sequenceRequest = buildTrackSequenceRequest(trackState, sequenceSource);
  const referenceSequenceRequest = buildTrackSequenceRequestForRange(
    refWindowStart,
    refWindowEnd,
    referenceSequenceSource
  );
  const shouldRetryFailedRequest = Boolean(
    sequenceRequest &&
    trackState.sequenceError &&
    trackState.sequenceRequestKey !== sequenceRequest.requestKey
  );
  const shouldRetryFailedReferenceRequest = Boolean(
    referenceSequenceRequest &&
    trackState.referenceSequenceError &&
    trackState.referenceSequenceRequestKey !== referenceSequenceRequest.requestKey
  );
  const shouldScheduleSequenceFetch = Boolean(
    sequenceSource &&
    sequenceRequest &&
    !hasFullSequence &&
    !trackState.sequenceLoading &&
    (!trackState.sequenceError || shouldRetryFailedRequest)
  );
  const shouldProbeSequenceLength = Boolean(
    sequenceSource &&
    !Number.isFinite(trackState.sequenceTotalLength) &&
    !trackState.sequenceLoading &&
    !trackState.sequenceError &&
    !sequenceRequest
  );
  const shouldScheduleReferenceSequenceFetch = Boolean(
    referenceSequenceSource &&
    referenceSequenceRequest &&
    !hasFullReferenceSequence &&
    !trackState.referenceSequenceLoading &&
    (!trackState.referenceSequenceError || shouldRetryFailedReferenceRequest)
  );
  const shouldProbeReferenceLength = Boolean(
    referenceSequenceSource &&
    !Number.isFinite(trackState.referenceSequenceTotalLength) &&
    !trackState.referenceSequenceLoading &&
    !trackState.referenceSequenceError &&
    !referenceSequenceRequest
  );
  const hasResolvedGff = hasTrackGffForSource(trackState, gffSource);
  const shouldScheduleGffFetch = Boolean(
    gffSource &&
    !hasResolvedGff &&
    !trackState.gffLoading &&
    !trackState.gffError
  );
  if (
    (!sequenceSource || windowSpan > TRACK_SEQUENCE_FETCH_MAX_SPAN || hasFullSequence) &&
    TRACK_FETCH_TIMERS[`${toolId}:sample`]
  ) {
    clearTimeout(TRACK_FETCH_TIMERS[`${toolId}:sample`]);
    delete TRACK_FETCH_TIMERS[`${toolId}:sample`];
  }
  if (
    (!sequenceSource || Number.isFinite(trackState.sequenceTotalLength)) &&
    TRACK_FETCH_TIMERS[`${toolId}:sample-probe`]
  ) {
    clearTimeout(TRACK_FETCH_TIMERS[`${toolId}:sample-probe`]);
    delete TRACK_FETCH_TIMERS[`${toolId}:sample-probe`];
  }
  if (
    (!referenceSequenceSource || refWindowSpan > TRACK_SEQUENCE_FETCH_MAX_SPAN || hasFullReferenceSequence) &&
    TRACK_FETCH_TIMERS[`${toolId}:reference`]
  ) {
    clearTimeout(TRACK_FETCH_TIMERS[`${toolId}:reference`]);
    delete TRACK_FETCH_TIMERS[`${toolId}:reference`];
  }
  if (
    (!referenceSequenceSource || Number.isFinite(trackState.referenceSequenceTotalLength)) &&
    TRACK_FETCH_TIMERS[`${toolId}:reference-probe`]
  ) {
    clearTimeout(TRACK_FETCH_TIMERS[`${toolId}:reference-probe`]);
    delete TRACK_FETCH_TIMERS[`${toolId}:reference-probe`];
  }
  if (
    (!gffSource || hasResolvedGff || trackState.gffError) &&
    TRACK_FETCH_TIMERS[`${toolId}:gff`]
  ) {
    clearTimeout(TRACK_FETCH_TIMERS[`${toolId}:gff`]);
    delete TRACK_FETCH_TIMERS[`${toolId}:gff`];
  }
  if (
    shouldScheduleSequenceFetch &&
    !document.body.classList.contains('track-dragging')
  ) {
    scheduleTrackSequenceFetch(toolId, data, {
      ...sequenceSource,
      sample: selectedSample
    });
  }
  if (
    shouldProbeSequenceLength &&
    !document.body.classList.contains('track-dragging')
  ) {
    scheduleTrackSequenceLengthProbe(toolId, data, { ...sequenceSource, sample: selectedSample });
  }
  if (
    shouldScheduleReferenceSequenceFetch &&
    !document.body.classList.contains('track-dragging')
  ) {
    scheduleTrackReferenceSequenceFetch(
      toolId,
      data,
      referenceSequenceSource,
      selectedSample,
      refWindowStart,
      refWindowEnd
    );
  }
  if (
    shouldProbeReferenceLength &&
    !document.body.classList.contains('track-dragging')
  ) {
    scheduleTrackReferenceLengthProbe(
      toolId,
      data,
      referenceSequenceSource,
      selectedSample
    );
  }
  if (
    shouldScheduleGffFetch &&
    !document.body.classList.contains('track-dragging')
  ) {
    scheduleTrackGffFetch(toolId, data, gffSource, selectedSample);
  }
  const fullSampleSequence = hasFullSequence
    ? {
        start: trackState.sequenceStart,
        end: trackState.sequenceEnd,
        sequence: trackState.sequence
      }
    : null;
  const fullReferenceSequence = hasFullReferenceSequence
    ? {
        start: trackState.referenceSequenceStart,
        end: trackState.referenceSequenceEnd,
        sequence: trackState.referenceSequence
      }
    : null;
  const svgDomKey = buildTrackSvgDomKey(
    sampleView,
    trackState,
    fullSampleSequence,
    fullReferenceSequence
  );
  const trackInfo = getCachedTrackSvgInfo(
    renderCache,
    svgDomKey,
    () => buildGenomicTrackSvg(
      sampleRecordsForTrack,
      trackState,
      fullSampleSequence,
      fullReferenceSequence,
      sampleView
    )
  );
  const windowStartBp = Math.round(trackState.windowStart);
  const windowEndBp = Math.round(trackState.windowEnd);
  const windowSpanBp = Math.max(1, Math.round(trackState.windowEnd - trackState.windowStart));
  const trackSearchOptions = sampleView.mutationSearchOptions;
  const trackMutationSelectOptions = sampleView.mutationSelectOptions;
  const trackSearchQuery = String(trackState.searchQuery || '');
  const trackSelectedKey = sampleRecordsForTrack.some(record => record.trackKey === trackState.selectedKey)
    ? trackState.selectedKey
    : '';
  const trackSearchHasMatches = Array.isArray(trackState.searchMatchKeys) && trackState.searchMatchKeys.length > 0;
  const trackSearchStatusText = getTrackSearchStatusText(trackState);
  const trackExpanded = isTrackExpanded(toolId);
  const gffAttributeKeys = Array.isArray(trackState.gffAttributeKeys) ? trackState.gffAttributeKeys : [];
  const gffResolvedField = getResolvedTrackGffLabelField(trackState);
  const gffSelectedField = String(trackState.gffLabelField || '').trim();
  const gffSelectedValue = gffSelectedField && gffAttributeKeys.includes(gffSelectedField)
    ? gffSelectedField
    : '__auto__';
  const gffFieldOptions = gffAttributeKeys.map(field =>
    `<option value="${escapeHtml(field)}" ${gffSelectedValue === field ? 'selected' : ''}>${escapeHtml(field)}</option>`
  ).join('');
  const sampleNtStatus = (() => {
    if (trackInfo.sampleNtSource === 'full') {
      return trackInfo.sampleNtLettersVisible
        ? 'Sample NT boxes with letters visible'
        : 'Sample NT boxes visible (zoom in to show letters)';
    }
    if (trackState.sequenceLoading && windowSpan <= TRACK_SEQUENCE_FETCH_MAX_SPAN) return 'Loading sample FASTA bases...';
    if (trackState.sequenceError && windowSpan <= TRACK_SEQUENCE_FETCH_MAX_SPAN) return 'Sample FASTA bases unavailable for this sample';
    if (!sequenceSource && windowSpan <= TRACK_SEQUENCE_FETCH_MAX_SPAN) return 'Sample FASTA source not mapped';
    if (trackInfo.sampleNtSource === 'mutations') {
      return trackInfo.sampleNtLettersVisible
        ? 'Sample mutation NT boxes visible'
        : 'Sample mutation NT boxes visible (zoom in to show letters)';
    }
    return 'Zoom in to display sample NT boxes';
  })();
  const referenceNtStatus = (() => {
    if (trackInfo.referenceNtSource === 'full') {
      return trackInfo.referenceNtLettersVisible
        ? 'Reference NT boxes with letters visible'
        : 'Reference NT boxes visible (zoom in to show letters)';
    }
    if (trackState.referenceSequenceLoading && refWindowSpan <= TRACK_SEQUENCE_FETCH_MAX_SPAN) return 'Loading reference FASTA bases...';
    if (trackState.referenceSequenceError && refWindowSpan <= TRACK_SEQUENCE_FETCH_MAX_SPAN) return 'Reference FASTA bases unavailable for this window';
    if (!referenceSequenceSource && refWindowSpan <= TRACK_SEQUENCE_FETCH_MAX_SPAN) return 'Reference FASTA source not mapped';
    if (trackInfo.referenceNtSource === 'mutations') {
      return trackInfo.referenceNtLettersVisible
        ? 'Reference mutation NT boxes visible'
        : 'Reference mutation NT boxes visible (zoom in to show letters)';
    }
    return 'Zoom in to display reference NT boxes';
  })();

  const controlsKey = [
    selectedSample,
    samples.join('|'),
    sampleView.cacheKey,
    gffSelectedValue,
    gffAttributeKeys.join('|'),
    trackExpanded ? '1' : '0'
  ].join('::');

  let didRebuildControls = false;
  if (renderCache.controlsKey !== controlsKey) {
    controlsEl.innerHTML = `
      <label>
        Genome sample
        <select class="viz-track-sample-select" ${samples.length <= 1 ? 'disabled' : ''}>
          ${samples.map(sample => `
            <option value="${escapeHtml(sample)}" ${sample === selectedSample ? 'selected' : ''}>
              ${escapeHtml(sample)}
            </option>
          `).join('')}
        </select>
      </label>
      <label>
        GFF label field
        <select class="viz-track-gff-field-select" ${gffAttributeKeys.length > 0 ? '' : 'disabled'}>
          <option value="__auto__" ${gffSelectedValue === '__auto__' ? 'selected' : ''}>
            Auto (${escapeHtml(gffResolvedField || 'best available')})
          </option>
          ${gffFieldOptions}
        </select>
      </label>
      <label>
        Preselect mutation
        <select class="viz-track-mutation-select">
          <option value="">Choose a mutation...</option>
          ${trackMutationSelectOptions}
        </select>
      </label>
      <label class="viz-track-search">
        Find mutation
        <div class="viz-track-search-wrap">
          <input
            type="text"
            class="viz-track-search-input"
            placeholder="Position, gene/label, AA change, ref>alt"
            value="${escapeHtml(trackSearchQuery)}"
            list="${toolId}-viz-track-search-options"
            autocomplete="off"
            aria-label="Search mutation and zoom to it"
          />
          <datalist id="${toolId}-viz-track-search-options">
            ${trackSearchOptions}
          </datalist>
          <button type="button" class="viz-track-btn viz-track-btn-search" data-track-action="find-mutation" title="Find and zoom to mutation">
            <span class="viz-track-btn-icon" aria-hidden="true">⌕</span>
            <span class="viz-track-btn-label">Find</span>
          </button>
          <button
            type="button"
            class="viz-track-btn viz-track-btn-search-nav"
            data-track-action="search-prev"
            title="Previous match"
            ${trackSearchHasMatches ? '' : 'disabled'}
          >
            <span class="viz-track-btn-icon" aria-hidden="true">‹</span>
            <span class="viz-track-btn-label">Prev</span>
          </button>
          <button
            type="button"
            class="viz-track-btn viz-track-btn-search-nav"
            data-track-action="search-next"
            title="Next match"
            ${trackSearchHasMatches ? '' : 'disabled'}
          >
            <span class="viz-track-btn-icon" aria-hidden="true">›</span>
            <span class="viz-track-btn-label">Next</span>
          </button>
          <span class="viz-track-search-status" aria-live="polite">${escapeHtml(trackSearchStatusText)}</span>
        </div>
      </label>
      <div class="viz-track-bottom-row">
        <label class="viz-track-window-range">
          Window start/end (bp)
          <div class="viz-track-window-range-wrap">
            <input
              type="number"
              class="viz-track-window-start-input"
              step="1"
              inputmode="numeric"
              value="${windowStartBp}"
              min="${Math.floor(trackState.domainStart)}"
              max="${Math.ceil(trackState.domainEnd)}"
              title="Window start in base pairs"
            />
            <span class="viz-track-window-sep">to</span>
            <input
              type="number"
              class="viz-track-window-end-input"
              step="1"
              inputmode="numeric"
              value="${windowEndBp}"
              min="${Math.floor(trackState.domainStart)}"
              max="${Math.ceil(trackState.domainEnd)}"
              title="Window end in base pairs"
            />
            <button type="button" class="viz-track-btn viz-track-btn-primary" data-track-action="apply-range" title="Apply start/end range">
              <span class="viz-track-btn-icon" aria-hidden="true">✓</span>
              <span class="viz-track-btn-label">Apply</span>
            </button>
          </div>
        </label>
        <div class="viz-track-nav" role="group" aria-label="Track navigation commands">
          <div class="viz-track-btn-group">
            <button type="button" class="viz-track-btn" data-track-action="zoom-in" title="Zoom in">
              <span class="viz-track-btn-icon" aria-hidden="true">＋</span>
              <span class="viz-track-btn-label">Zoom In</span>
            </button>
            <button type="button" class="viz-track-btn" data-track-action="zoom-out" title="Zoom out">
              <span class="viz-track-btn-icon" aria-hidden="true">－</span>
              <span class="viz-track-btn-label">Zoom Out</span>
            </button>
            <button type="button" class="viz-track-btn" data-track-action="pan-left" title="Pan left">
              <span class="viz-track-btn-icon" aria-hidden="true">←</span>
              <span class="viz-track-btn-label">Pan Left</span>
            </button>
            <button type="button" class="viz-track-btn" data-track-action="pan-right" title="Pan right">
              <span class="viz-track-btn-icon" aria-hidden="true">→</span>
              <span class="viz-track-btn-label">Pan Right</span>
            </button>
          </div>
          <button type="button" class="viz-track-btn" data-track-action="reset" title="Reset view">
            <span class="viz-track-btn-icon" aria-hidden="true">↺</span>
            <span class="viz-track-btn-label">Reset</span>
          </button>
          <button type="button" class="viz-track-btn viz-track-btn-emphasis" data-track-action="toggle-expand" title="${trackExpanded ? 'Reduce genomic track' : 'Expand genomic track'}">
            <span class="viz-track-btn-icon" aria-hidden="true">${trackExpanded ? '⤡' : '⤢'}</span>
            <span class="viz-track-btn-label">${trackExpanded ? 'Reduce' : 'Expand'}</span>
          </button>
        </div>
      </div>
      <div class="visually-hidden" aria-live="polite">
        Current window ${formatTrackBp(windowStartBp)} - ${formatTrackBp(windowEndBp)} bp (${formatTrackBp(windowSpanBp)} bp span).
      </div>
    `;
    controlsEl.dataset.trackControlsBound = '';
    renderCache.controlsKey = controlsKey;
    didRebuildControls = true;
  }

  patchTrackControls(controlsEl, {
    windowStartBp,
    windowEndBp,
    searchStatusText: trackSearchStatusText,
    searchHasMatches: trackSearchHasMatches,
    searchQuery: trackSearchQuery,
    selectedKey: trackSelectedKey
  });

  const { metaEl, svgHostEl, detailHostEl } = ensureTrackBodyShell(bodyEl);
  const metaHtml = `
    <span><strong>${trackInfo.inWindowMutations}</strong> mutation(s) in view · ${trackInfo.totalMutations} total</span>
    <span>Genome range 1 - ${formatTrackBp(trackState.domainEnd)}</span>
    <span>Sample window ${formatTrackBp(trackState.windowStart)} - ${formatTrackBp(trackState.windowEnd)}</span>
    <span>${trackInfo.geneMode === 'gff'
      ? `${trackInfo.visibleGeneCount} GFF gene interval(s) in view${gffResolvedField ? ` · field ${gffResolvedField}` : ''}`
      : `${trackInfo.visibleGeneCount} inferred gene interval(s) in view`}</span>
    <span>${trackInfo.visibleAaCount || 0}/${trackInfo.totalAaCount || 0} amino-acid codons in view${trackInfo.aaGeneLabel ? ` · gene ${trackInfo.aaGeneLabel}` : ''}</span>
    <span>NT: ${sampleNtStatus} · ${referenceNtStatus}</span>
  `;
  if (metaEl) metaEl.innerHTML = metaHtml;

  const detailHtml = buildTrackDetailHtml(trackInfo, trackState);
  if (detailHostEl) detailHostEl.innerHTML = detailHtml;

  let didReplaceSvg = false;
  if (svgHostEl && renderCache.svgDomKey !== svgDomKey) {
    svgHostEl.innerHTML = trackInfo.svg;
    renderCache.svgDomKey = svgDomKey;
    didReplaceSvg = true;
  }
  section.classList.remove('hidden');
  if (didReplaceSvg || didRebuildControls) {
    attachTrackInteractions(toolId, data, sampleRecordsForTrack);
  }
}

/**
 * Setup visualization for a specific tool
 */
function setupToolVisualization(toolId, title, columnName) {
  const vizBtn = document.getElementById(`${toolId}-visualize-btn`);
  const vizPanel = document.getElementById(`${toolId}-visualization`);
  const closeBtn = vizPanel?.querySelector('.btn-close-viz');
  const fullscreenBtn = vizPanel?.querySelector('.btn-fullscreen-viz');
  ensureVisualizationControls(vizPanel, toolId);
  ensureVisualizationInsights(vizPanel, toolId);
  ensureSingleSampleSummary(vizPanel, toolId);

  // Initialize charts storage
  initToolCharts(toolId);
  activeVisualizations[toolId] = { columnName, visible: false, defaultTitle: title };

  if (vizBtn) {
    vizBtn.addEventListener('click', () => {
      if (vizBtn.dataset.ready !== 'true') return;
      showToolVisualization(toolId, columnName);
    });
  }

  if (closeBtn) {
    closeBtn.addEventListener('click', () => {
      setTrackExpanded(toolId, false);
      vizPanel.classList.add('hidden');
      vizPanel.classList.remove('fullscreen');
      vizBtn?.classList.remove('hidden');
      const expandIcon = fullscreenBtn?.querySelector('.icon-expand');
      const collapseIcon = fullscreenBtn?.querySelector('.icon-collapse');
      expandIcon?.classList.remove('hidden');
      collapseIcon?.classList.add('hidden');
      if (activeVisualizations[toolId]) {
        activeVisualizations[toolId].visible = false;
      }
    });
  }

  if (fullscreenBtn) {
    fullscreenBtn.addEventListener('click', () => {
      handleFullscreenToggle(vizPanel, fullscreenBtn, toolId, columnName);
    });
  }

  const metricSelect = vizPanel?.querySelector('.viz-metric-select');
  if (metricSelect) {
    metricSelect.addEventListener('change', () => {
      setVisualizationPreferences(toolId, { metric: metricSelect.value });
      refreshToolVisualization(toolId, columnName);
    });
  }

  const topNSelect = vizPanel?.querySelector('.viz-topn-select');
  if (topNSelect) {
    topNSelect.addEventListener('change', () => {
      setVisualizationPreferences(toolId, { topN: topNSelect.value });
      refreshToolVisualization(toolId, columnName);
    });
  }

  const groupOtherCheckbox = vizPanel?.querySelector('.viz-other-checkbox');
  if (groupOtherCheckbox) {
    groupOtherCheckbox.addEventListener('change', () => {
      setVisualizationPreferences(toolId, { groupOther: groupOtherCheckbox.checked });
      refreshToolVisualization(toolId, columnName);
    });
  }
}

/**
 * Handle fullscreen toggle for visualization
 */
function handleFullscreenToggle(vizPanel, fullscreenBtn, toolId, columnName) {
  const expandIcon = fullscreenBtn.querySelector('.icon-expand');
  const collapseIcon = fullscreenBtn.querySelector('.icon-collapse');

  if (vizPanel.classList.contains('fullscreen')) {
    vizPanel.classList.remove('fullscreen');
    expandIcon?.classList.remove('hidden');
    collapseIcon?.classList.add('hidden');

    // Recreate charts at original size
    const data = getPanelResultsData(`${toolId}-results`);
    if (data) {
      setTimeout(() => {
        renderVisualization(toolId, columnName);
      }, TIMING.CHART_RESIZE);
    }
  } else {
    vizPanel.classList.add('fullscreen');
    expandIcon?.classList.add('hidden');
    collapseIcon?.classList.remove('hidden');

    // Recreate charts at fullscreen size
    const data = getPanelResultsData(`${toolId}-results`);
    if (data) {
      setTimeout(() => {
        renderVisualization(toolId, columnName);
      }, TIMING.CHART_RESIZE);
    }
  }
}

/**
 * Show visualization for a tool
 */
function showToolVisualization(toolId, columnName) {
  const vizPanel = document.getElementById(`${toolId}-visualization`);
  const vizBtn = document.getElementById(`${toolId}-visualize-btn`);
  const data = getPanelResultsData(`${toolId}-results`);

  if (!data || !vizPanel) return;

  if (!hasChartLibrary()) {
    logMessage('Chart.js is unavailable. Cannot render visualization.', 'error');
    return;
  }

  if (!Array.isArray(data.rows) || data.rows.length === 0) {
    logMessage('No data found for visualization', 'warning');
    return;
  }

  // Store active visualization info for theme refresh
  activeVisualizations[toolId] = {
    ...activeVisualizations[toolId],
    columnName,
    visible: true
  };

  vizBtn?.classList.add('hidden');
  vizPanel.classList.remove('hidden');

  renderVisualization(toolId, columnName);

  setTimeout(() => {
    vizPanel.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  }, TIMING.SCROLL_SETTLE);
}

/**
 * Parse tool data for visualization
 * For tools with per-marker detail rows (classify, splitfq), groups by genome
 * and counts unique genomes per lineage instead of individual marker rows.
 */
function resolvePrimaryColumnIndex(headers, primaryColumn) {
  // Find the primary column
  let colIdx = headers.findIndex(h =>
    h.toLowerCase() === primaryColumn.toLowerCase()
  );

  // Try other common variations
  if (colIdx === -1) {
    const variations = [
      'classification', 'lineage', 'prediction', 'best_match',
      'match', 'result', 'genotype', 'type', 'category'
    ];
    colIdx = headers.findIndex(h =>
      variations.some(v => h.toLowerCase().includes(v))
    );
  }

  // Use last column if not found
  if (colIdx === -1) colIdx = headers.length - 1;
  return colIdx;
}

function parseToolData(data, primaryColumn) {
  const counts = {};
  const headers = data.headers;
  const rows = data.rows;

  const colIdx = resolvePrimaryColumnIndex(headers, primaryColumn);

  // Check if this data has a genome column (detailed results with multiple rows per genome)
  const genomeColIdx = headers.findIndex(h =>
    h.toLowerCase() === 'genome' || h.toLowerCase() === 'file'
  );

  if (genomeColIdx !== -1 && genomeColIdx !== colIdx) {
    // Group rows by genome and determine majority lineage per genome
    const genomeLineages = {};
    rows.forEach(row => {
      const genome = row[genomeColIdx]?.trim();
      const value = row[colIdx]?.trim();
      if (!genome || isUnknownOutcome(value)) return;
      if (!genomeLineages[genome]) genomeLineages[genome] = {};
      genomeLineages[genome][value] = (genomeLineages[genome][value] || 0) + 1;
    });

    // For each genome, pick the majority lineage, then count genomes per lineage
    for (const genome of Object.keys(genomeLineages)) {
      const lineageCounts = genomeLineages[genome];
      const majorLineage = Object.entries(lineageCounts)
        .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0], undefined, { sensitivity: 'base' }))[0]?.[0];
      if (!isUnknownOutcome(majorLineage)) {
        counts[majorLineage] = (counts[majorLineage] || 0) + 1;
      }
    }
  } else {
    // Simple mode: one row per sample (summary files, predict, match)
    rows.forEach(row => {
      const value = row[colIdx]?.trim();
      if (!isUnknownOutcome(value)) {
        counts[value] = (counts[value] || 0) + 1;
      }
    });
  }

  return counts;
}

/**
 * Re-tally a WHO resistance run by drug.
 *
 * For a resistance panel the label column holds the whole marker path
 * (drug;resistance;marker;grade;gene;mutation), so the generic tally produces a
 * distribution over raw paths and the KPI header ends up presenting one of them
 * as the sample's "call". Group by sample instead and describe each sample by
 * the drugs it carries markers for, which is what the cards, the chart and the
 * single-sample outcome should all be about.
 *
 * Returns `null` for a genuine lineage panel, leaving the normal tally alone.
 */
function buildDrCallSummary(data, primaryColumn) {
  const headers = Array.isArray(data?.headers) ? data.headers : [];
  const rows = Array.isArray(data?.rows) ? data.rows : [];
  if (headers.length === 0 || rows.length === 0) return null;

  const labelColIdx = resolvePrimaryColumnIndex(headers, primaryColumn);
  if (labelColIdx === -1) return null;

  const paths = rows.map(row => normalizeValue(row[labelColIdx])).filter(Boolean);
  if (!isDrPanel(paths)) return null;

  const sampleColIdx = findSampleColumnIndex(headers, labelColIdx);
  const bySample = new Map();
  rows.forEach(row => {
    const lineagePath = normalizeValue(row[labelColIdx]);
    if (!lineagePath) return;
    const sample = sampleColIdx !== -1
      ? (normalizeValue(row[sampleColIdx]) || 'Sample')
      : 'Sample';
    if (!bySample.has(sample)) bySample.set(sample, []);
    bySample.get(sample).push({ lineagePath });
  });
  if (bySample.size === 0) return null;

  const counts = {};
  const sampleSummaries = new Map();
  for (const [sample, records] of bySample) {
    const profile = buildDrProfile(records);
    sampleSummaries.set(sample, summariseDrProfile(profile));
    // One count per sample per drug, so a multi-sample run answers
    // "how many samples carry markers for this drug".
    profile.forEach(entry => {
      const label = drCallLabel(entry);
      counts[label] = (counts[label] || 0) + 1;
    });
  }

  return { counts, sampleSummaries };
}

/** Retitle the panel: "Lineage Distribution" is wrong for a resistance run. */
function setVisualizationTitle(toolId, title) {
  const heading = document.querySelector(`#${toolId}-visualization .visualization-title h3`);
  if (!heading) return;
  const next = title || activeVisualizations[toolId]?.defaultTitle;
  if (next && heading.textContent !== next) heading.textContent = next;
}

function renderVisualization(toolId, columnName) {
  const data = getPanelResultsData(`${toolId}-results`);
  if (!data) return;

  const drCalls = buildDrCallSummary(data, columnName);
  const counts = drCalls ? drCalls.counts : parseToolData(data, columnName);
  setVisualizationTitle(toolId, drCalls ? 'Resistance profile' : null);

  const singleSampleSummary = buildSingleSampleVisualizationSummary(toolId, data, counts, columnName);
  if (singleSampleSummary && drCalls) {
    // Describe the sample by its resistance verdict rather than by whichever
    // marker path happened to sort first.
    singleSampleSummary.outcomeLabel =
      drCalls.sampleSummaries.get(singleSampleSummary.sampleName)
      || summariseDrProfile([]);
  }
  renderSingleSampleSummary(toolId, singleSampleSummary);

  if (singleSampleSummary) {
    destroyToolCharts(toolId);
    resetChartCanvases(toolId);
    clearChartEmptyState(toolId);
    createToolLegend(toolId, [], [], [], 0);
  } else {
    createToolCharts(toolId, counts, getVisualizationPreferences(toolId));
  }

  renderToolInsights(toolId, data, counts, columnName);
  renderGenomicTrack(toolId, data);
}

function buildSingleSampleVisualizationSummary(toolId, data, counts, primaryColumn) {
  const headers = Array.isArray(data?.headers) ? data.headers : [];
  const rows = Array.isArray(data?.rows) ? data.rows : [];
  if (rows.length === 0 || headers.length === 0) return null;

  const labelColIdx = resolvePrimaryColumnIndex(headers, primaryColumn);
  const sampleCount = estimateSampleCount(data, labelColIdx);
  if (sampleCount !== 1) return null;

  const sampleColIdx = findSampleColumnIndex(headers, labelColIdx);
  const sampleName = sampleColIdx === -1
    ? 'Single sample'
    : (normalizeValue(rows.find(row => normalizeValue(row?.[sampleColIdx]))?.[sampleColIdx]) || 'Single sample');

  const sortedEntries = Object.entries(counts)
    .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0], undefined, { sensitivity: 'base' }));
  const topEntry = sortedEntries[0] || null;
  const assignedCount = topEntry?.[1] || 0;
  const unresolvedCount = Math.max(sampleCount - assignedCount, 0);
  const quality = findQualityMetricSummary(toolId, data);

  return {
    sampleName,
    outcomeLabel: topEntry?.[0] || 'Unresolved',
    assignedCount,
    unresolvedCount,
    rowCount: rows.length,
    qualityLabel: quality?.label || 'Quality',
    qualityValue: quality ? formatPlotDecimal(quality.median, 2) : 'N/A'
  };
}

/**
 * Create charts for a tool
 */
function createToolCharts(toolId, dataCounts, preferences = DEFAULT_VIZ_PREFS) {
  // Destroy existing charts
  destroyToolCharts(toolId);
  resetChartCanvases(toolId);

  if (!hasChartLibrary()) return;
  clearChartEmptyState(toolId);

  const prefs = {
    ...DEFAULT_VIZ_PREFS,
    ...preferences
  };

  // Keep chart and legend ordering deterministic and aligned.
  const entries = Object.entries(dataCounts)
    .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0], undefined, { sensitivity: 'base' }));
  if (entries.length === 0) {
    setChartEmptyState(toolId, 'No assigned categories to chart');
    createToolLegend(toolId, [], [], [], 0);
    return;
  }
  const overallTotal = entries.reduce((acc, [, value]) => acc + value, 0);

  const requestedTopN = prefs.topN === 'all'
    ? Number.POSITIVE_INFINITY
    : Number.parseInt(prefs.topN, 10);
  const topN = Number.isFinite(requestedTopN) ? Math.max(1, requestedTopN) : Number.POSITIVE_INFINITY;

  let chartEntries = entries;
  if (Number.isFinite(topN) && entries.length > topN) {
    const topEntries = entries.slice(0, topN);
    if (prefs.groupOther) {
      const otherCount = entries.slice(topN).reduce((acc, [, value]) => acc + value, 0);
      topEntries.push(['Other', otherCount]);
    }
    chartEntries = topEntries;
  }

  const labels = chartEntries.map(([label]) => label);
  const countValues = chartEntries.map(([, value]) => value);
  const total = countValues.reduce((a, b) => a + b, 0);
  const metricValues = prefs.metric === 'percent'
    ? countValues.map(value => (overallTotal > 0 ? (value / overallTotal) * 100 : 0))
    : countValues;

  const colors = generateChartColors(labels);

  // Get theme-aware colors
  const isDark = document.documentElement.getAttribute('data-theme') === 'dark';
  const textColor = isDark ? '#e7e9ea' : '#18181b';
  const gridColor = isDark ? 'rgba(255, 255, 255, 0.15)' : 'rgba(0, 0, 0, 0.1)';
  const axisColor = isDark ? '#8b98a5' : '#71717a';

  // Create donut chart
  const donutCtx = document.getElementById(`${toolId}-donut-chart`)?.getContext('2d');
  if (donutCtx) {
    const donutChart = new Chart(donutCtx, {
      type: 'doughnut',
      data: {
        labels: labels,
        datasets: [{
          data: countValues,
          backgroundColor: colors,
          borderColor: 'transparent',
          borderWidth: 0,
          hoverOffset: 8
        }]
      },
      options: {
        locale: 'en-US',
        responsive: true,
        maintainAspectRatio: true,
        cutout: '60%',
        plugins: {
          legend: { display: false },
          tooltip: {
            backgroundColor: 'rgba(0,0,0,0.85)',
            titleColor: '#fff',
            bodyColor: '#fff',
            titleFont: { size: 13, weight: '600' },
            bodyFont: { size: 12 },
            padding: 12,
            cornerRadius: 8,
            callbacks: {
              label: (ctx) => {
                const count = countValues[ctx.dataIndex] || 0;
                const pct = formatPlotDecimal(overallTotal > 0 ? (count / overallTotal) * 100 : 0, 1);
                return ` ${ctx.label}: ${count} (${pct}%)`;
              }
            }
          }
        }
      }
    });
    setToolChart(toolId, 'donut', donutChart);
  }

  // Create horizontal bar chart
  const barCtx = document.getElementById(`${toolId}-bar-chart`)?.getContext('2d');
  if (barCtx) {
    const barChart = new Chart(barCtx, {
      type: 'bar',
      data: {
        labels: labels,
        datasets: [{
          data: metricValues,
          backgroundColor: colors,
          borderRadius: 6,
          borderSkipped: false
        }]
      },
      options: {
        locale: 'en-US',
        indexAxis: 'y',
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          legend: { display: false },
          tooltip: {
            backgroundColor: 'rgba(0,0,0,0.85)',
            titleColor: '#fff',
            bodyColor: '#fff',
            titleFont: { size: 13, weight: '600' },
            bodyFont: { size: 12 },
            padding: 12,
            cornerRadius: 8,
            callbacks: {
              label: (ctx) => {
                const count = countValues[ctx.dataIndex] || 0;
                const pct = formatPlotDecimal(overallTotal > 0 ? (count / overallTotal) * 100 : 0, 1);
                const valueLabel = prefs.metric === 'percent'
                  ? `Value: ${formatPlotDecimal(ctx.raw, 1)}%`
                  : `Count: ${count}`;
                return ` ${valueLabel} (${pct}% of total)`;
              }
            }
          }
        },
        scales: {
          x: {
            beginAtZero: true,
            border: { color: axisColor },
            grid: { color: gridColor },
            ticks: {
              color: textColor,
              font: { size: 11 },
              callback: (value) => {
                if (prefs.metric === 'percent') return `${formatPlotDecimal(value, 1)}%`;
                return String(value).replace(',', '.');
              }
            },
            title: {
              display: true,
              text: prefs.metric === 'percent' ? 'Percent (%)' : 'Genomes',
              color: axisColor,
              font: { size: 12 }
            }
          },
          y: {
            border: { color: axisColor },
            grid: { display: false },
            ticks: { color: textColor, font: { size: 12, weight: '500' } }
          }
        }
      }
    });
    setToolChart(toolId, 'bar', barChart);
  }

  // Create legend
  createToolLegend(toolId, labels, countValues, colors, overallTotal);
}

/**
 * Generate colors for chart labels
 */
function generateChartColors(labels) {
  let fallbackIdx = 0;

  return labels.map(label => {
    const upperLabel = label.toUpperCase();

    // Check for lineage color matches (longest first to handle L10 before L1)
    for (const [lineage, color] of Object.entries(lineageColors).sort((a, b) => b[0].length - a[0].length)) {
      if (upperLabel.includes(lineage.toUpperCase()) ||
          upperLabel.startsWith(lineage.toUpperCase()) ||
          upperLabel.split(';').some(part => part.trim().toUpperCase().startsWith(lineage.toUpperCase()))) {
        return color;
      }
    }

    // Fallback color
    const color = fallbackColors[fallbackIdx % fallbackColors.length];
    fallbackIdx++;
    return color;
  });
}

/**
 * Create chart legend
 */
function createToolLegend(toolId, labels, values, colors, total) {
  const legendEl = document.getElementById(`${toolId}-chart-legend`);
  if (!legendEl) return;

  legendEl.replaceChildren();
  labels.forEach((label, i) => {
    const pct = formatPlotDecimal(total > 0 ? (values[i] / total) * 100 : 0, 1);
    const item = document.createElement('div');
    item.className = 'legend-item';
    item.dataset.index = String(i);
    item.dataset.tool = toolId;

    const colorEl = document.createElement('span');
    colorEl.className = 'legend-color';
    colorEl.style.backgroundColor = colors[i];
    item.appendChild(colorEl);

    const labelEl = document.createElement('span');
    labelEl.className = 'legend-label';
    labelEl.textContent = label;
    item.appendChild(labelEl);

    const countEl = document.createElement('span');
    countEl.className = 'legend-count';
    countEl.textContent = `(${values[i]} - ${pct}%)`;
    item.appendChild(countEl);

    legendEl.appendChild(item);
  });

  // Add hover interactivity
  legendEl.querySelectorAll('.legend-item').forEach(item => {
    item.addEventListener('mouseenter', () => {
      const idx = parseInt(item.dataset.index, 10);
      const tool = item.dataset.tool;
      highlightChartSegment(tool, idx);
    });
    item.addEventListener('mouseleave', () => {
      const tool = item.dataset.tool;
      resetChartHighlight(tool);
    });
  });
}

/**
 * Highlight chart segment on hover
 */
function highlightChartSegment(toolId, index) {
  const charts = getToolCharts(toolId);
  if (charts?.donut) {
    charts.donut.setActiveElements([{ datasetIndex: 0, index }]);
    charts.donut.update();
  }
  if (charts?.bar) {
    charts.bar.setActiveElements([{ datasetIndex: 0, index }]);
    charts.bar.update();
  }
}

/**
 * Reset chart highlight
 */
function resetChartHighlight(toolId) {
  const charts = getToolCharts(toolId);
  if (charts?.donut) {
    charts.donut.setActiveElements([]);
    charts.donut.update();
  }
  if (charts?.bar) {
    charts.bar.setActiveElements([]);
    charts.bar.update();
  }
}

/**
 * Refresh all active charts for theme change
 * Called when the theme is toggled to update chart colors
 */
export function refreshChartsForTheme() {
  for (const [toolId, vizInfo] of Object.entries(activeVisualizations)) {
    if (!vizInfo.visible) continue;
    refreshToolVisualization(toolId, vizInfo.columnName);
  }
}

export function refreshToolChartsLayout(toolId, options = {}) {
  const { forceRecreate = false } = options;
  const vizPanel = document.getElementById(`${toolId}-visualization`);
  if (!vizPanel || vizPanel.classList.contains('hidden')) return false;

  const charts = getToolCharts(toolId);
  const chartTypes = [];
  if (charts?.donut) chartTypes.push(charts.donut);
  if (charts?.bar) chartTypes.push(charts.bar);

  if (chartTypes.length === 0) {
    if (vizPanel.classList.contains('single-sample-mode')) {
      return true;
    }
    const vizInfo = activeVisualizations[toolId];
    const columnName = vizInfo?.columnName || visualizationColumns[toolId];
    if (columnName) {
      refreshToolVisualization(toolId, columnName);
      return false;
    }
    return false;
  }

  let shouldRecreate = forceRecreate;
  chartTypes.forEach(chart => {
    if (!chart || !chart.canvas || !chart.canvas.isConnected) {
      shouldRecreate = true;
      return;
    }
    const parent = chart.canvas.parentElement;
    const canvas = chart.canvas;
    const parentRect = parent?.getBoundingClientRect();
    const canvasRect = canvas?.getBoundingClientRect();
    if (
      !parent ||
      !canvasRect ||
      parentRect.width <= 0 ||
      parentRect.height <= 0 ||
      canvasRect.width <= 0 ||
      canvasRect.height <= 0
    ) {
      shouldRecreate = true;
      return;
    }
    if (!Number.isFinite(chart.width) || chart.width <= 0 || !Number.isFinite(chart.height) || chart.height <= 0) {
      shouldRecreate = true;
      return;
    }
    try {
      if (shouldRecreate) return;
      const width = Math.max(1, Math.floor(canvasRect.width));
      const height = Math.max(1, Math.floor(canvasRect.height));
      if (chart.canvas) {
        chart.resize(width, height);
      } else {
        chart.resize();
      }
      chart.update('none');
    } catch {
      // Ignore resize/update errors for already-disposed charts.
      shouldRecreate = true;
    }
  });

  if (shouldRecreate) {
    const vizInfo = activeVisualizations[toolId];
    const columnName = vizInfo?.columnName || visualizationColumns[toolId];
    if (columnName) {
      refreshToolVisualization(toolId, columnName);
      return false;
    }
    return false;
  }

  return true;
}

/**
 * Refresh a specific tool visualization after new results are loaded.
 * If the panel is currently open, charts are rebuilt from the latest table data.
 */
export function refreshVisualizationForTool(toolId) {
  const vizInfo = activeVisualizations[toolId];
  const columnName = vizInfo?.columnName || visualizationColumns[toolId];
  if (!columnName || !vizInfo?.visible) return;
  refreshToolVisualization(toolId, columnName);
}

/**
 * Reset visualization state for a tool before running a new analysis.
 */
export function resetVisualizationForTool(toolId) {
  const vizPanel = document.getElementById(`${toolId}-visualization`);
  const vizBtn = document.getElementById(`${toolId}-visualize-btn`);
  const legendEl = document.getElementById(`${toolId}-chart-legend`);
  const singleSampleSummaryEl = document.getElementById(`${toolId}-viz-single-sample-summary`);

  setTrackExpanded(toolId, false);
  destroyToolCharts(toolId);

  if (activeVisualizations[toolId]) {
    activeVisualizations[toolId].visible = false;
  }

  if (vizPanel) {
    vizPanel.classList.add('hidden');
    vizPanel.classList.remove('fullscreen');
    vizPanel.classList.remove('single-sample-mode');
    // See resetInlineResults: the tab bar is only rebuilt on multi-marker runs,
    // so it must be torn down here or a later run shows stale tabs.
    vizPanel.querySelector('.result-tabs')?.remove();
    const chartContainer = vizPanel.querySelector('.chart-container');
    chartContainer?.classList.remove('hidden');
  }

  const fullscreenBtn = vizPanel?.querySelector('.btn-fullscreen-viz');
  const expandIcon = fullscreenBtn?.querySelector('.icon-expand');
  const collapseIcon = fullscreenBtn?.querySelector('.icon-collapse');
  expandIcon?.classList.remove('hidden');
  collapseIcon?.classList.add('hidden');

  if (legendEl) legendEl.replaceChildren();
  legendEl?.classList.remove('hidden');
  if (singleSampleSummaryEl) {
    singleSampleSummaryEl.classList.add('hidden');
    singleSampleSummaryEl.replaceChildren();
  }
  clearToolInsights(toolId);
  clearGenomicTrack(toolId);
  clearChartEmptyState(toolId);

  if (vizBtn) {
    vizBtn.classList.add('hidden');
    vizBtn.dataset.ready = 'false';
  }
}

function refreshToolVisualization(toolId, columnName) {
  const vizPanel = document.getElementById(`${toolId}-visualization`);
  if (!vizPanel || vizPanel.classList.contains('hidden')) {
    if (activeVisualizations[toolId]) {
      activeVisualizations[toolId].visible = false;
    }
    return;
  }

  const data = getPanelResultsData(`${toolId}-results`);
  if (!data) return;

  renderVisualization(toolId, columnName);
}

/**
 * Render a summary card for training results with KPI tiles.
 * Called from forms.js when train completes with report data.
 */
export function renderTrainSummaryCard(toolId, report) {
  const container = document.getElementById('train-summary-card');
  if (!container || !report) return;

  const isCv = report.cv_mean_accuracy_pct != null;
  const accuracyValue = isCv
    ? `${Number(report.cv_mean_accuracy_pct).toFixed(2)}% ± ${Number(report.cv_std_accuracy_pct || 0).toFixed(2)}%`
    : `${Number(report.accuracy_pct || 0).toFixed(2)}%`;
  const accuracySub = isCv
    ? `${(report.cv_fold_accuracies || []).length}-fold stratified CV`
    : `${report.n_test || 0} test samples`;

  const oobCard = report.oob_accuracy_pct != null
    ? [{ label: 'OOB Accuracy', value: `${Number(report.oob_accuracy_pct).toFixed(2)}%`, sub: 'out-of-bag estimate' }]
    : [];

  const cards = [
    { label: 'Accuracy', value: accuracyValue, sub: accuracySub },
    ...oobCard,
    { label: 'Training Samples', value: String(report.n_training || 0), sub: `${report.n_classes || 0} classes` },
    { label: 'Features', value: Number(report.n_features || 0).toLocaleString(), sub: `k=${report.kmer_size || '?'}` },
    { label: 'Trees', value: String(report.n_trees || 0), sub: 'Random Forest ensemble' },
    { label: 'Model Size', value: `${Number(report.model_size_mb || 0).toFixed(2)} MB`, sub: 'zstd compressed' },
    { label: 'Time', value: `${Number(report.training_time_secs || 0).toFixed(1)}s`, sub: 'total training time' }
  ];

  const classListHtml = Array.isArray(report.class_names) && report.class_names.length > 0
    ? `<div class="train-summary-classes">
         <span class="train-summary-classes-label">Classes detected:</span>
         <span class="train-summary-classes-list">${report.class_names.map(c => escapeHtml(c)).join(', ')}</span>
       </div>`
    : '';

  container.innerHTML = `
    <h3 class="train-summary-title">Training Report</h3>
    <div class="viz-kpi-grid">
      ${cards.map(card => `
        <article class="viz-kpi-card">
          <span class="viz-kpi-label">${escapeHtml(card.label)}</span>
          <strong class="viz-kpi-value">${escapeHtml(String(card.value))}</strong>
          <span class="viz-kpi-sub">${escapeHtml(card.sub)}</span>
        </article>
      `).join('')}
    </div>
    ${classListHtml}
  `;
  container.classList.remove('hidden');
}
