// ============================================================================
// Results Module - Inline Results Viewer
// ============================================================================

import { TIMING } from './config.js';
import {
  setLastOutputPath,
  setPanelResultsData,
  clearPanelResultsData,
  getPanelResultsData
} from './state.js';
import { escapeHtml, getFileName } from './utils.js';
import { openFileLocation } from './tauri.js';
import { logMessage } from './console.js';
import {
  focusGenomicTrackRecordBySourceRow,
  refreshToolChartsLayout,
  refreshVisualizationForTool
} from './visualization.js';

// Debounce timers per viewer for filter input
const filterTimers = {};
const chartRefreshTimers = new Map();
const chartRefreshGeneration = new Map();
const DEFAULT_PAGE_SIZE = 100;
const PAGE_SIZE_OPTIONS = [50, 100, 250, 500];
const TRACK_RECORD_SELECTION_EVENT = 'pathotypr:track-record-selected';
const CHART_READY_RETRIES = 6;
const QUICK_FILTER_MAX_OPTIONS = 5;
const QUICK_FILTER_UNKNOWN_VALUES = new Set([
  '',
  'n/a',
  'na',
  'unknown',
  'unresolved',
  'unassigned',
  'none',
  '-',
  'null',
  '.'
]);

function normalizeSourceRowIndex(rawRowIndex) {
  const value = Number(rawRowIndex);
  return Number.isFinite(value) ? Math.trunc(value) : Number.NaN;
}

function getTrackViewerIdFromViewer(viewerId) {
  return String(viewerId || '').replace(/-results$/, '');
}

function applySelectedRowPage(viewer, data) {
  if (!data || !Array.isArray(data.filteredRows)) return;
  const sourceRowIndex = normalizeSourceRowIndex(data.selectedSourceRowIndex);
  if (!Number.isFinite(sourceRowIndex)) return;

  const pageSize = Math.max(1, Number.parseInt(data.pageSize || DEFAULT_PAGE_SIZE, 10));
  const filteredIndex = data.filteredRows.findIndex(
    row => normalizeSourceRowIndex(row?.__sourceRowIndex) === sourceRowIndex
  );
  if (filteredIndex < 0) return;

  const totalPages = Math.max(1, Math.ceil(data.filteredRows.length / pageSize));
  data.page = Math.min(Math.max(1, Math.floor(filteredIndex / pageSize) + 1), totalPages);
}

function ensureRowsHaveSourceIndex(rows) {
  if (!Array.isArray(rows)) return [];
  return rows.map((row, rowIndex) => {
    if (!Array.isArray(row)) return row;
    if (!Number.isFinite(normalizeSourceRowIndex(row.__sourceRowIndex))) {
      row.__sourceRowIndex = rowIndex;
    }
    return row;
  });
}

function normalizeHeader(value) {
  return String(value || '')
    .trim()
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '_');
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

function normalizeNtValue(value) {
  const text = String(value || '').trim().toUpperCase();
  if (!text || text === '.') return '';
  if (/^[ACGTN]$/.test(text)) return text;
  return text;
}

function inferAltNtFromKmer(value) {
  const text = String(value || '').trim().toUpperCase();
  if (!text) return '';
  const center = Math.floor(text.length / 2);
  const base = text[center] || '';
  return /^[ACGTN]$/.test(base) ? base : '';
}

function augmentNucleotideColumns(headers, rows) {
  if (!Array.isArray(headers) || !Array.isArray(rows) || headers.length === 0) {
    return { headers, rows };
  }

  const refIdx = findHeaderIndexByTokens(headers, [
    'ref_nt',
    'reference_nt',
    'ref_allele',
    'reference_allele'
  ]);
  const altIdx = findHeaderIndexByTokens(headers, [
    'mut_nt',
    'mutation_nt',
    'alt_nt',
    'alt_allele',
    'alternate_allele'
  ]);
  const kmerIdx = findHeaderIndexByTokens(headers, ['k-mer', 'kmer', 'marker_kmer']);

  if (refIdx === -1 && altIdx === -1 && kmerIdx === -1) {
    return { headers, rows };
  }

  const nextHeaders = [...headers];
  let refOutIdx = findHeaderIndexByTokens(nextHeaders, ['ref_nt', 'reference_nt'], { exactOnly: true });
  if (refOutIdx === -1) {
    refOutIdx = nextHeaders.length;
    nextHeaders.push('Ref NT');
  }
  let altOutIdx = findHeaderIndexByTokens(nextHeaders, ['mut_nt', 'mutation_nt', 'alt_nt'], { exactOnly: true });
  if (altOutIdx === -1) {
    altOutIdx = nextHeaders.length;
    nextHeaders.push('Mut NT');
  }

  const nextRows = rows.map(row => {
    const next = [...row];
    if (Object.prototype.hasOwnProperty.call(row, '__sourceRowIndex')) {
      next.__sourceRowIndex = row.__sourceRowIndex;
    }
    const refNt = refIdx >= 0 ? normalizeNtValue(next[refIdx]) : '';
    let altNt = altIdx >= 0 ? normalizeNtValue(next[altIdx]) : '';
    if (!altNt && kmerIdx >= 0) altNt = inferAltNtFromKmer(next[kmerIdx]);

    next[refOutIdx] = refNt;
    next[altOutIdx] = altNt;
    return next;
  });

  return { headers: nextHeaders, rows: nextRows };
}

function parseNumericCell(value) {
  const normalized = String(value ?? '').trim().replace(',', '.');
  if (!normalized) return Number.NaN;
  const parsed = Number.parseFloat(normalized);
  return Number.isFinite(parsed) ? parsed : Number.NaN;
}

function isUnknownOutcomeValue(value) {
  return QUICK_FILTER_UNKNOWN_VALUES.has(String(value ?? '').trim().toLowerCase());
}

function buildQuickFilterOptions(data) {
  const headers = Array.isArray(data?.headers) ? data.headers : [];
  const rows = Array.isArray(data?.rows) ? data.rows : [];
  const options = [{
    id: 'all',
    label: 'All',
    title: 'Show all rows',
    count: rows.length,
    filter: null
  }];

  if (headers.length === 0 || rows.length === 0) return options;

  const categoryCol = findHeaderIndexByTokens(headers, [
    'major_lineage',
    'lineage',
    'prediction',
    'best_match',
    'classification',
    'genotype',
    'result',
    'call'
  ]);

  if (categoryCol !== -1) {
    const counts = new Map();
    rows.forEach(row => {
      const value = String(row?.[categoryCol] ?? '').trim();
      if (!value) return;
      counts.set(value, (counts.get(value) || 0) + 1);
    });

    Array.from(counts.entries())
      .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0], undefined, { sensitivity: 'base' }))
      .slice(0, QUICK_FILTER_MAX_OPTIONS)
      .forEach(([value, count], idx) => {
        options.push({
          id: `category-${idx}`,
          label: value.length > 24 ? `${value.slice(0, 23)}…` : value,
          title: `${value} (${count} rows)`,
          count,
          filter: {
            type: 'category',
            col: categoryCol,
            value
          }
        });
      });

    const unresolvedCount = rows.reduce((acc, row) => (
      isUnknownOutcomeValue(row?.[categoryCol]) ? acc + 1 : acc
    ), 0);
    if (unresolvedCount > 0) {
      options.push({
        id: 'category-unresolved',
        label: 'Unresolved',
        title: `${unresolvedCount} unresolved rows`,
        count: unresolvedCount,
        filter: {
          type: 'unknown',
          col: categoryCol
        }
      });
    }
  }

  const confidenceCol = findHeaderIndexByTokens(headers, [
    'confidence',
    'probability',
    'posterior',
    'score',
    'support',
    'identity',
    'ani'
  ]);

  if (confidenceCol !== -1) {
    const numericValues = rows
      .map(row => parseNumericCell(row?.[confidenceCol]))
      .filter(value => Number.isFinite(value));

    if (numericValues.length >= 3) {
      const max = Math.max(...numericValues);
      const min = Math.min(...numericValues);
      const normalizedConfidence = min >= 0 && max <= 1.1;
      const highThreshold = normalizedConfidence ? 0.9 : 90;
      const lowThreshold = normalizedConfidence ? 0.5 : 50;
      const formatter = normalizedConfidence
        ? (value) => value.toFixed(2)
        : (value) => String(Math.round(value));

      options.push({
        id: 'confidence-high',
        label: `Conf ≥ ${formatter(highThreshold)}`,
        title: `Confidence >= ${formatter(highThreshold)}`,
        filter: {
          type: 'numeric',
          col: confidenceCol,
          op: '>=',
          value: highThreshold
        }
      });

      options.push({
        id: 'confidence-low',
        label: `Conf < ${formatter(lowThreshold)}`,
        title: `Confidence < ${formatter(lowThreshold)}`,
        filter: {
          type: 'numeric',
          col: confidenceCol,
          op: '<',
          value: lowThreshold
        }
      });
    }
  }

  return options;
}

function matchesQuickFilter(row, quickFilter) {
  if (!quickFilter || typeof quickFilter !== 'object') return true;
  const colIndex = Number.parseInt(quickFilter.col, 10);
  if (!Number.isInteger(colIndex) || colIndex < 0 || colIndex >= row.length) return true;
  const value = row[colIndex];

  if (quickFilter.type === 'category') {
    const target = String(quickFilter.value ?? '').trim().toLowerCase();
    return String(value ?? '').trim().toLowerCase() === target;
  }

  if (quickFilter.type === 'unknown') {
    return isUnknownOutcomeValue(value);
  }

  if (quickFilter.type === 'numeric') {
    const left = parseNumericCell(value);
    const right = Number(quickFilter.value);
    if (!Number.isFinite(left) || !Number.isFinite(right)) return false;
    switch (quickFilter.op) {
      case '>': return left > right;
      case '<': return left < right;
      case '>=': return left >= right;
      case '<=': return left <= right;
      case '=': return left === right;
      default: return true;
    }
  }

  return true;
}

function ensureResultsEnhancements(viewer) {
  const filterRow = viewer.querySelector('.results-viewer-filter');
  if (filterRow && !viewer.querySelector('.results-quick-filters')) {
    const quickFilters = document.createElement('div');
    quickFilters.className = 'results-quick-filters hidden';
    quickFilters.innerHTML = `
      <span class="results-quick-filters-label">Quick filters</span>
      <div class="results-filter-chips" role="group" aria-label="Quick filters"></div>
    `;
    filterRow.insertAdjacentElement('afterend', quickFilters);
  }

  const tableContainer = viewer.querySelector('.results-viewer-table');
  if (tableContainer && !viewer.querySelector('.results-viewer-layout')) {
    const layout = document.createElement('div');
    layout.className = 'results-viewer-layout';
    tableContainer.insertAdjacentElement('beforebegin', layout);
    layout.appendChild(tableContainer);

    const detailPanel = document.createElement('aside');
    detailPanel.className = 'results-row-detail hidden';
    detailPanel.innerHTML = `
      <div class="results-row-detail-header">
        <h4>Row details</h4>
        <div class="results-row-detail-actions">
          <button type="button" class="btn-results-action btn-results-action-mini" data-action="copy-selected-row">
            Copy row
          </button>
          <button type="button" class="results-row-detail-close" aria-label="Close details">&times;</button>
        </div>
      </div>
      <div class="results-row-detail-body"></div>
    `;
    layout.appendChild(detailPanel);
  }
}

function renderQuickFilterChips(viewer, data) {
  const container = viewer.querySelector('.results-quick-filters');
  const chipsEl = container?.querySelector('.results-filter-chips');
  if (!container || !chipsEl) return;

  const options = Array.isArray(data?.quickFilterOptions) ? data.quickFilterOptions : [];
  if (options.length <= 1) {
    chipsEl.replaceChildren();
    container.classList.add('hidden');
    return;
  }

  const activeId = String(data?.quickFilterId || 'all');
  chipsEl.innerHTML = options.map(option => {
    const isActive = option.id === activeId;
    const count = Number.isFinite(option.count) ? `<span class="results-filter-chip-count">${option.count}</span>` : '';
    const title = option.title ? ` title="${escapeHtml(option.title)}"` : '';
    return `
      <button
        type="button"
        class="results-filter-chip${isActive ? ' is-active' : ''}"
        data-filter-id="${escapeHtml(option.id)}"${title}
      >
        <span>${escapeHtml(option.label)}</span>
        ${count}
      </button>
    `;
  }).join('');
  container.classList.remove('hidden');
}

function renderSelectedRowDetail(viewer, data) {
  const detailPanel = viewer.querySelector('.results-row-detail');
  const detailBody = detailPanel?.querySelector('.results-row-detail-body');
  const layout = viewer.querySelector('.results-viewer-layout');
  if (!detailPanel || !detailBody || !data) return;

  const selectedSourceRowIndex = normalizeSourceRowIndex(data.selectedSourceRowIndex);
  if (!Number.isFinite(selectedSourceRowIndex)) {
    detailPanel.classList.add('hidden');
    detailBody.replaceChildren();
    layout?.classList.remove('has-detail');
    return;
  }

  const selectedRow = Array.isArray(data.rows)
    ? data.rows.find(row => normalizeSourceRowIndex(row?.__sourceRowIndex) === selectedSourceRowIndex)
    : null;

  if (!selectedRow || !Array.isArray(data.headers)) {
    detailPanel.classList.add('hidden');
    detailBody.replaceChildren();
    layout?.classList.remove('has-detail');
    return;
  }

  const visibleInCurrentFilters = Array.isArray(data.filteredRows)
    ? data.filteredRows.some(row => normalizeSourceRowIndex(row?.__sourceRowIndex) === selectedSourceRowIndex)
    : false;

  const populatedFields = data.headers
    .map((header, idx) => ({ header, value: selectedRow[idx] ?? '' }))
    .filter(item => String(item.value ?? '').trim() !== '');
  const shownFields = populatedFields.slice(0, 16);

  detailBody.innerHTML = `
    <p class="results-row-detail-meta">
      Source row #${selectedSourceRowIndex + 1}
      <span>${visibleInCurrentFilters ? 'Visible in current filters' : 'Filtered out in current view'}</span>
    </p>
    <dl class="results-row-detail-grid">
      ${shownFields.map(item => `
        <div class="results-row-detail-item">
          <dt>${escapeHtml(item.header)}</dt>
          <dd>${escapeHtml(String(item.value))}</dd>
        </div>
      `).join('')}
    </dl>
    ${populatedFields.length > shownFields.length
      ? `<p class="results-row-detail-note">Showing ${shownFields.length} of ${populatedFields.length} non-empty fields.</p>`
      : ''}
  `;
  detailPanel.classList.remove('hidden');
  layout?.classList.add('has-detail');
}

function clearAllFilters(viewer, data) {
  const filterInput = viewer.querySelector('.results-filter-input');
  if (filterInput) filterInput.value = '';
  const thead = viewer.querySelector('thead');
  thead?.querySelectorAll('.col-filter').forEach(input => {
    input.value = '';
  });
  thead?.querySelectorAll('.col-filter-op').forEach(select => {
    select.value = '=';
  });
  data.colFilters = {};
  data.quickFilterId = 'all';
  data.quickFilter = null;
  data.page = 1;
}

function debouncedApplyFilters(viewer, delay = TIMING.FILTER_DEBOUNCE) {
  const id = viewer.id;
  if (filterTimers[id]) clearTimeout(filterTimers[id]);
  filterTimers[id] = setTimeout(() => {
    applyFilters(viewer, { forceRecreate: false });
    filterTimers[id] = null;
  }, delay);
}

function requestChartsLayoutRefresh(viewer, options = {}) {
  const panelId = getTrackViewerIdFromViewer(viewer?.id);
  if (!panelId) return;

  const forceRecreate = options?.forceRecreate === true;
  const generation = (chartRefreshGeneration.get(panelId) || 0) + 1;
  chartRefreshGeneration.set(panelId, generation);

  if (chartRefreshTimers.has(panelId)) {
    clearTimeout(chartRefreshTimers.get(panelId));
    chartRefreshTimers.delete(panelId);
  }

  const scheduleRefresh = (delayMs, callback) => {
    chartRefreshTimers.set(
      panelId,
      setTimeout(() => {
        if (chartRefreshGeneration.get(panelId) !== generation) return;
        callback();
      }, delayMs)
    );
  };

  const isReady = () => {
    const donutCanvas = document.getElementById(`${panelId}-donut-chart`);
    const barCanvas = document.getElementById(`${panelId}-bar-chart`);
    const hasUsableCanvas = (canvas) => {
      if (!canvas || !canvas.isConnected) return false;
      const canvasRect = canvas.getBoundingClientRect();
      if (canvasRect.width <= 0 || canvasRect.height <= 0) return false;
      const wrapper = canvas.closest('.chart-wrapper');
      if (!wrapper) return false;
      const wrapperRect = wrapper.getBoundingClientRect();
      return wrapperRect.width > 0 && wrapperRect.height > 0;
    };
    return hasUsableCanvas(donutCanvas) && hasUsableCanvas(barCanvas);
  };

  const runRefresh = (attempt = 0) => {
    if (chartRefreshGeneration.get(panelId) !== generation) return;

    const vizPanel = document.getElementById(`${panelId}-visualization`);
    if (!vizPanel || vizPanel.classList.contains('hidden')) return;

    if (!forceRecreate && attempt === 0) {
      const success = refreshToolChartsLayout(panelId, options);
      if (success) {
        chartRefreshTimers.delete(panelId);
        return;
      }
      scheduleRefresh(TIMING.CHART_RESIZE, () => runRefresh(1));
      return;
    }

    if (attempt === 0) {
      refreshVisualizationForTool(panelId);
    }

    requestAnimationFrame(() => {
      if (chartRefreshGeneration.get(panelId) !== generation) return;

      if (isReady()) {
        const success = refreshToolChartsLayout(panelId, { forceRecreate: false });
        if (!success && attempt < CHART_READY_RETRIES) {
          scheduleRefresh(TIMING.CHART_RESIZE, () => runRefresh(attempt + 1));
          return;
        }

        chartRefreshTimers.delete(panelId);
        requestAnimationFrame(() => {
          window.dispatchEvent(new Event('resize'));
        });
        return;
      }

      if (attempt < CHART_READY_RETRIES) {
        const delay = TIMING.CHART_RESIZE * Math.max(1, attempt + 1);
        scheduleRefresh(delay, () => runRefresh(attempt + 1));
        return;
      }

      refreshToolChartsLayout(panelId, { forceRecreate: true });
      chartRefreshTimers.delete(panelId);
      requestAnimationFrame(() => {
        window.dispatchEvent(new Event('resize'));
      });
    });
  };

  if (!forceRecreate) {
    scheduleRefresh(TIMING.CHART_RESIZE, () => runRefresh(0));
    return;
  }

  runRefresh(0);
}

function ensureResultsControls(viewer) {
  const filterRow = viewer.querySelector('.results-viewer-filter');
  if (!filterRow || filterRow.querySelector('.results-filter-actions')) return;

  const actions = document.createElement('div');
  actions.className = 'results-filter-actions';
  actions.innerHTML = `
    <button type="button" class="btn-results-action" data-action="copy-filtered" title="Copy filtered rows as TSV">
      Copy Filtered TSV
    </button>
    <button type="button" class="btn-results-action" data-action="reset-filters" title="Reset search and all filters">
      Reset Filters
    </button>
    <div class="results-pagination">
      <label class="results-page-size-label">
        Rows:
        <select class="results-page-size">
          ${PAGE_SIZE_OPTIONS.map(size => `<option value="${size}">${size}</option>`).join('')}
        </select>
      </label>
      <button type="button" class="btn-results-page" data-direction="-1" title="Previous page">&lsaquo;</button>
      <span class="results-page-info">Page 1 / 1</span>
      <button type="button" class="btn-results-page" data-direction="1" title="Next page">&rsaquo;</button>
    </div>
  `;
  filterRow.appendChild(actions);
}

/**
 * Setup results modal and inline viewers
 */
export function setupResultsModal() {
  // Setup inline results viewers for each panel
  setupInlineResultsViewer('predict-results');
  setupInlineResultsViewer('classify-results');
  setupInlineResultsViewer('splitfq-results');
  setupInlineResultsViewer('match-results');
}

/**
 * Setup individual inline results viewer
 */
function setupInlineResultsViewer(viewerId) {
  const viewer = document.getElementById(viewerId);
  if (!viewer) return;
  ensureResultsControls(viewer);
  ensureResultsEnhancements(viewer);

  const panelId = viewerId.replace('-results', '');
  const trackToolId = getTrackViewerIdFromViewer(viewerId);
  const closeBtn = viewer.querySelector('.btn-close-results');
  const openFolderBtn = viewer.querySelector('.btn-open-folder');
  const fullscreenBtn = viewer.querySelector('.btn-fullscreen-results');
  const filterInput = viewer.querySelector('.results-filter-input');
  const showResultsBtn = document.getElementById(`${panelId}-show-results`);

  if (closeBtn) {
    closeBtn.addEventListener('click', () => {
      viewer.classList.add('hidden');
      viewer.classList.remove('fullscreen');
      if (showResultsBtn && viewer.dataset.outputPath) {
        showResultsBtn.classList.remove('hidden');
      }
    });
  }

  if (openFolderBtn) {
    openFolderBtn.addEventListener('click', () => {
      const path = viewer.dataset.outputPath;
      if (path) openOutputFolder(path);
    });
  }

  if (fullscreenBtn) {
    fullscreenBtn.addEventListener('click', () => {
      const expandIcon = fullscreenBtn.querySelector('.icon-expand');
      const collapseIcon = fullscreenBtn.querySelector('.icon-collapse');

      if (viewer.classList.contains('fullscreen')) {
        viewer.classList.add('exiting-fullscreen');
        expandIcon?.classList.remove('hidden');
        collapseIcon?.classList.add('hidden');

        setTimeout(() => {
          viewer.classList.remove('fullscreen', 'exiting-fullscreen');
        }, TIMING.FULLSCREEN_EXIT);
      } else {
        viewer.classList.add('fullscreen');
        expandIcon?.classList.add('hidden');
        collapseIcon?.classList.remove('hidden');
      }
    });
  }

  if (filterInput) {
    filterInput.addEventListener('input', () => {
      debouncedApplyFilters(viewer);
    });
  }

  const pageSizeSelect = viewer.querySelector('.results-page-size');
  if (pageSizeSelect) {
    pageSizeSelect.value = String(DEFAULT_PAGE_SIZE);
    pageSizeSelect.addEventListener('change', () => {
      const data = getPanelResultsData(viewerId);
      if (!data) return;
      const parsed = Number.parseInt(pageSizeSelect.value, 10);
      data.pageSize = Number.isFinite(parsed) ? parsed : DEFAULT_PAGE_SIZE;
      data.page = 1;
      applyFilters(viewer, { forceRecreate: false });
    });
  }

  viewer.querySelectorAll('.btn-results-page').forEach(btn => {
    btn.addEventListener('click', () => {
      const direction = Number.parseInt(btn.dataset.direction || '0', 10);
      const data = getPanelResultsData(viewerId);
      if (!data || !direction) return;
      data.page = Math.max(1, (data.page || 1) + direction);
      applyFilters(viewer, { forceRecreate: false });
    });
  });

  const copyFilteredBtn = viewer.querySelector('.btn-results-action[data-action="copy-filtered"]');
  if (copyFilteredBtn) {
    copyFilteredBtn.addEventListener('click', async () => {
      const data = getPanelResultsData(viewerId);
      if (!data) return;
      await copyFilteredRowsAsTsv(data);
    });
  }

  const resetFiltersBtn = viewer.querySelector('.btn-results-action[data-action="reset-filters"]');
  if (resetFiltersBtn) {
    resetFiltersBtn.addEventListener('click', () => {
      const data = getPanelResultsData(viewerId);
      if (!data) return;
      clearAllFilters(viewer, data);
      applyFilters(viewer, { forceRecreate: false });
    });
  }

  const quickFilters = viewer.querySelector('.results-quick-filters');
  if (quickFilters && !quickFilters.dataset.bound) {
    quickFilters.dataset.bound = 'true';
    quickFilters.addEventListener('click', (event) => {
      const chip = event.target.closest('.results-filter-chip[data-filter-id]');
      if (!chip) return;
      const data = getPanelResultsData(viewerId);
      if (!data) return;

      const nextId = String(chip.dataset.filterId || 'all');
      const option = (Array.isArray(data.quickFilterOptions) ? data.quickFilterOptions : [])
        .find(entry => entry.id === nextId);
      data.quickFilterId = nextId;
      data.quickFilter = option?.filter || null;
      data.page = 1;
      applyFilters(viewer, { forceRecreate: false });
    });
  }

  const detailPanel = viewer.querySelector('.results-row-detail');
  if (detailPanel && !detailPanel.dataset.bound) {
    detailPanel.dataset.bound = 'true';
    detailPanel.addEventListener('click', async (event) => {
      const closeBtn = event.target.closest('.results-row-detail-close');
      if (closeBtn) {
        const data = getPanelResultsData(viewerId);
        if (!data) return;
        data.selectedSourceRowIndex = -1;
        applyFilters(viewer, { forceRecreate: false });
        return;
      }

      const copySelectedBtn = event.target.closest('[data-action="copy-selected-row"]');
      if (copySelectedBtn) {
        const data = getPanelResultsData(viewerId);
        if (!data) return;
        const selectedSourceRowIndex = normalizeSourceRowIndex(data.selectedSourceRowIndex);
        if (!Number.isFinite(selectedSourceRowIndex)) return;
        const selectedRow = Array.isArray(data.rows)
          ? data.rows.find(row => normalizeSourceRowIndex(row?.__sourceRowIndex) === selectedSourceRowIndex)
          : null;
        if (!selectedRow) return;
        await copyRowAsTsv(selectedRow);
      }
    });
  }

  const tbody = viewer.querySelector('tbody');
  if (tbody) {
    tbody.addEventListener('click', async (event) => {
      const row = event.target.closest('tr[data-row-index]');
      if (!row) return;
      const index = Number.parseInt(row.dataset.rowIndex || '-1', 10);
      const data = getPanelResultsData(viewerId);
      if (!data || index < 0 || !Array.isArray(data.filteredRows) || !data.filteredRows[index]) return;
      const sourceRowIndex = normalizeSourceRowIndex(row.dataset.sourceRowIndex);
      if (Number.isFinite(sourceRowIndex)) {
        data.selectedSourceRowIndex = sourceRowIndex;
        applySelectedRowPage(viewer, data);
        if (trackToolId === 'classify') {
          focusGenomicTrackRecordBySourceRow(trackToolId, data, sourceRowIndex);
        }
        applyFilters(viewer, { forceRecreate: false });
      }
      await copyRowAsTsv(data.filteredRows[index]);
    });
  }

  const onTrackRecordSelection = (event) => {
    const detail = event?.detail || {};
    if (detail.toolId !== trackToolId) return;
    const sourceRowIndex = normalizeSourceRowIndex(detail.sourceRowIndex);
    if (!Number.isFinite(sourceRowIndex)) return;
    const data = getPanelResultsData(viewerId);
    if (!data) return;
    data.selectedSourceRowIndex = sourceRowIndex;
    applySelectedRowPage(viewer, data);
    applyFilters(viewer, { forceRecreate: false });
  };
  if (!viewer.dataset.trackSelectionBound) {
    viewer.dataset.trackSelectionBound = 'true';
    document.addEventListener(TRACK_RECORD_SELECTION_EVENT, onTrackRecordSelection);
  }

  if (showResultsBtn) {
    showResultsBtn.addEventListener('click', () => {
      viewer.classList.remove('hidden');
      showResultsBtn.classList.add('hidden');
      setTimeout(() => {
        viewer.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
      }, TIMING.SCROLL_SETTLE);
    });
  }
}

/**
 * Show inline results for a panel
 * @param {string} panelId - The panel ID
 * @param {string} outputPath - Path to the output file
 * @param {string} tsvContent - TSV content to display
 * @param {string|string[]} [inputFiles] - Optional input file(s) to show in File column
 */
export function showInlineResults(panelId, outputPath, tsvContent, inputFiles = null, context = null) {
  const viewerId = `${panelId}-results`;
  const viewer = document.getElementById(viewerId);
  const showResultsBtn = document.getElementById(`${panelId}-show-results`);

  if (!viewer) {
    console.warn('Results viewer not found:', viewerId);
    logMessage('Results viewer not found for panel: ' + panelId, 'warning');
    return;
  }

  // Store output path
  viewer.dataset.outputPath = outputPath;
  setLastOutputPath(outputPath);

  // Update path display
  const pathEl = viewer.querySelector('.results-viewer-path');
  if (pathEl) {
    const fileName = getFileName(outputPath);
    pathEl.textContent = fileName;
    pathEl.title = outputPath;
  }

  // Reset filter and fullscreen
  const filterInput = viewer.querySelector('.results-filter-input');
  if (filterInput) filterInput.value = '';
  viewer.classList.remove('fullscreen');

  const trimmedContent = String(tsvContent || '').trim();
  if (!trimmedContent) {
    console.warn('[showInlineResults] Empty TSV content');
    logMessage('Results file is empty.', 'warning');
    return;
  }

  // Parse and display TSV
  const lines = trimmedContent
    .split(/\r?\n/)
    .map(line => line.replace(/\r$/, ''))
    .filter(line => line.trim() !== '');
  if (lines.length === 0) {
    console.warn('[showInlineResults] No lines in TSV content');
    return;
  }

  // Parse headers and rows from TSV
  let headers = lines[0].split('\t');
  let rows = lines.slice(1).map(line => line.split('\t'));

  // Add "File" column with input filename if provided
  if (inputFiles) {
    // Normalize to array and extract just filenames
    const files = Array.isArray(inputFiles) ? inputFiles : [inputFiles];
    const fileNames = files.map(f => getFileName(f));

    // Add File column at the beginning
    headers = ['File', ...headers];
    rows = rows.map((row, idx) => {
      // Use corresponding filename or first one if single file for all rows
      const fileName = fileNames.length === 1 ? fileNames[0] : (fileNames[idx] || fileNames[0]);
      return [fileName, ...row];
    });
  }

  rows = ensureRowsHaveSourceIndex(rows);

  ({ headers, rows } = augmentNucleotideColumns(headers, rows));

  // Detect numeric columns
  const numericCols = headers.map((_, colIdx) => {
    if (rows.length === 0) return false;
    let numericCount = 0;
    rows.forEach(row => {
      const val = row[colIdx];
      if (val !== undefined && val !== '' && !isNaN(parseFloat(val))) {
        numericCount++;
      }
    });
    return numericCount / rows.length > 0.5;
  });

  // Store data for filtering and sorting
  const panelData = {
    headers,
    rows,
    filteredRows: rows,
    sortCol: -1,
    sortDir: 'asc',
    colFilters: {},
    numericCols,
    page: 1,
    pageSize: DEFAULT_PAGE_SIZE,
    selectedSourceRowIndex: -1,
    quickFilterOptions: [],
    quickFilterId: 'all',
    quickFilter: null,
    context: context && typeof context === 'object' ? context : null
  };
  panelData.quickFilterOptions = buildQuickFilterOptions(panelData);
  setPanelResultsData(viewerId, panelData);
  renderQuickFilterChips(viewer, panelData);
  renderSelectedRowDetail(viewer, panelData);

  const pageSizeSelect = viewer.querySelector('.results-page-size');
  if (pageSizeSelect) pageSizeSelect.value = String(DEFAULT_PAGE_SIZE);

  // Populate table
  const thead = viewer.querySelector('thead');
  const tbody = viewer.querySelector('tbody');

  if (thead) {
    thead.innerHTML = `
      <tr>
        ${headers.map((h, i) => `
          <th scope="col" data-col="${i}">
            <div class="th-content">
              <span>${escapeHtml(h)}</span>
              <svg class="sort-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" aria-hidden="true">
                <path d="M12 5v14M5 12l7-7 7 7"/>
              </svg>
            </div>
            ${numericCols[i] ? `
              <div class="numeric-filter" data-col="${i}">
                <select class="col-filter-op" data-col="${i}">
                  <option value="=">=</option>
                  <option value=">">&gt;</option>
                  <option value="<">&lt;</option>
                  <option value=">=">&ge;</option>
                  <option value="<=">&le;</option>
                </select>
                <input type="number" class="col-filter col-filter-num" placeholder="0" data-col="${i}" step="any">
              </div>
            ` : `
              <input type="text" class="col-filter" placeholder="Filter..." data-col="${i}">
            `}
          </th>
        `).join('')}
      </tr>
    `;

    // Setup sorting
    thead.querySelectorAll('th').forEach(th => {
      th.addEventListener('click', (e) => {
        if (e.target.classList.contains('col-filter')) return;
        const col = parseInt(th.dataset.col);
        sortTable(viewer, col);
      });
    });

    // Setup column filters
    thead.querySelectorAll('.col-filter').forEach(input => {
      input.addEventListener('input', (e) => {
        e.stopPropagation();
        const col = parseInt(input.dataset.col);
        const data = getPanelResultsData(viewerId);

        if (input.classList.contains('col-filter-num')) {
          const opSelect = thead.querySelector(`.col-filter-op[data-col="${col}"]`);
          const op = opSelect?.value || '=';
          const numVal = input.value !== '' ? parseFloat(input.value) : null;
          data.colFilters[col] = { type: 'numeric', op, value: numVal };
        } else {
          data.colFilters[col] = { type: 'text', value: input.value.toLowerCase() };
        }
        debouncedApplyFilters(viewer);
      });
      input.addEventListener('click', (e) => e.stopPropagation());
    });

    // Setup numeric filter operator changes
    thead.querySelectorAll('.col-filter-op').forEach(select => {
      select.addEventListener('change', (e) => {
        e.stopPropagation();
        const col = parseInt(select.dataset.col);
        const data = getPanelResultsData(viewerId);
        const numInput = thead.querySelector(`.col-filter-num[data-col="${col}"]`);
        const numVal = numInput && numInput.value !== '' ? parseFloat(numInput.value) : null;
        data.colFilters[col] = { type: 'numeric', op: select.value, value: numVal };
        debouncedApplyFilters(viewer);
      });
      select.addEventListener('click', (e) => e.stopPropagation());
    });
  }

  applyFilters(viewer, { forceRecreate: false });

  // Hide show results button, show viewer
  if (showResultsBtn) showResultsBtn.classList.add('hidden');
  viewer.classList.remove('hidden');

  // Show visualization button
  const vizBtn = document.getElementById(`${panelId}-visualize-btn`);
  if (vizBtn) {
    vizBtn.classList.remove('hidden');
    vizBtn.dataset.ready = 'true';
  }

  // Scroll to results
  setTimeout(() => {
    viewer.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  }, TIMING.SCROLL_SETTLE);
}

/**
 * Render table body with rows
 */
function renderTableBody(viewer, rows, globalFilter = '', rowOffset = 0, emptyMessage = 'No matching results') {
  const tbody = viewer.querySelector('tbody');
  if (!tbody) return;
  if (rows.length === 0) {
    const colCount = viewer.querySelector('thead tr')?.children.length || 1;
    tbody.innerHTML = `<tr class="empty-state-row"><td colspan="${colCount}">${emptyMessage}</td></tr>`;
    return;
  }
  const lowerFilter = globalFilter.toLowerCase();
  const data = getPanelResultsData(viewer.id);
  const selectedSourceRowIndex = normalizeSourceRowIndex(data?.selectedSourceRowIndex);

  // Detect confidence column for row-level quality highlighting (predict panel).
  const confColIdx = detectConfidenceColumn(viewer, data);

  tbody.innerHTML = rows.map((row, idx) => {
    const classes = [];
    if (selectedSourceRowIndex === normalizeSourceRowIndex(row?.__sourceRowIndex)) {
      classes.push('row-selected');
    }
    if (confColIdx >= 0) {
      const confVal = Number.parseFloat(row[confColIdx]);
      if (Number.isFinite(confVal)) {
        if (confVal < 0.5) classes.push('row-very-low-confidence');
        else if (confVal < 0.7) classes.push('row-low-confidence');
      }
    }
    return `<tr
      data-row-index="${rowOffset + idx}"
      data-source-row-index="${normalizeSourceRowIndex(row?.__sourceRowIndex)}"
      class="${classes.join(' ')}"
    >${row.map(cell => {
      const escaped = escapeHtml(cell ?? '');
      const matchClass = lowerFilter && String(cell || '').toLowerCase().includes(lowerFilter)
        ? ' class="filter-match"'
        : '';
      return `<td${matchClass}>${escaped}</td>`;
    }).join('')}</tr>`;
  }).join('');
}

/** Detect the Confidence column index for predict-panel result tables. Returns -1 if not applicable. */
function detectConfidenceColumn(viewer, data) {
  const viewerId = String(viewer?.id || '');
  if (!viewerId.startsWith('predict')) return -1;
  const headers = Array.isArray(data?.headers) ? data.headers : [];
  const tokens = ['confidence', 'probability', 'score'];
  return headers.findIndex(h => {
    const norm = String(h || '').toLowerCase().replace(/[_\- ]/g, '');
    return tokens.some(t => norm.includes(t));
  });
}

function sortRows(rows, data) {
  if (data.sortCol < 0) return [...rows];
  const colIndex = data.sortCol;
  return [...rows].sort((a, b) => {
    const valA = a[colIndex] || '';
    const valB = b[colIndex] || '';

    const numA = Number.parseFloat(valA);
    const numB = Number.parseFloat(valB);
    if (!Number.isNaN(numA) && !Number.isNaN(numB)) {
      return data.sortDir === 'asc' ? numA - numB : numB - numA;
    }

    const cmp = String(valA).localeCompare(String(valB), undefined, { numeric: true, sensitivity: 'base' });
    return data.sortDir === 'asc' ? cmp : -cmp;
  });
}

function rowMatchesFilters(row, globalFilter, colFilters, quickFilter) {
  const lowerGlobal = globalFilter.toLowerCase();
  const matchesGlobal = !lowerGlobal || row.some(cell =>
    String(cell || '').toLowerCase().includes(lowerGlobal)
  );
  if (!matchesGlobal) return false;
  if (!matchesQuickFilter(row, quickFilter)) return false;

  for (let colIdx = 0; colIdx < row.length; colIdx++) {
    const colFilter = colFilters[colIdx];
    if (!colFilter) continue;

    const text = String(row[colIdx] || '').toLowerCase();
    if (colFilter.type === 'text' && colFilter.value) {
      if (!text.includes(colFilter.value)) return false;
      continue;
    }

    if (colFilter.type === 'numeric' && colFilter.value !== null) {
      const cellNum = Number.parseFloat(row[colIdx]);
      if (Number.isNaN(cellNum)) return false;
      switch (colFilter.op) {
        case '=':
          if (cellNum !== colFilter.value) return false;
          break;
        case '>':
          if (!(cellNum > colFilter.value)) return false;
          break;
        case '<':
          if (!(cellNum < colFilter.value)) return false;
          break;
        case '>=':
          if (!(cellNum >= colFilter.value)) return false;
          break;
        case '<=':
          if (!(cellNum <= colFilter.value)) return false;
          break;
        default:
          break;
      }
    }
  }

  return true;
}

/**
 * Sort table by column
 */
function sortTable(viewer, colIndex) {
  const viewerId = viewer.id;
  const data = getPanelResultsData(viewerId);
  if (!data) return;

  // Toggle sort direction
  if (data.sortCol === colIndex) {
    data.sortDir = data.sortDir === 'asc' ? 'desc' : 'asc';
  } else {
    data.sortCol = colIndex;
    data.sortDir = 'asc';
  }
  data.page = 1;

  // Update header styles
  viewer.querySelectorAll('th').forEach(th => {
    th.classList.remove('sort-asc', 'sort-desc');
    if (parseInt(th.dataset.col) === colIndex) {
      th.classList.add(data.sortDir === 'asc' ? 'sort-asc' : 'sort-desc');
    }
  });
  applyFilters(viewer, { forceRecreate: false });
}

/**
 * Apply all filters to table
 */
function applyFilters(viewer, options = {}) {
  const viewerId = viewer.id;
  const data = getPanelResultsData(viewerId);
  if (!data) return;

  const globalFilter = viewer.querySelector('.results-filter-input')?.value.toLowerCase() || '';
  const sortedRows = sortRows(data.rows, data);
  const filteredRows = sortedRows.filter(row =>
    rowMatchesFilters(row, globalFilter, data.colFilters, data.quickFilter)
  );
  data.filteredRows = filteredRows;

  const pageSize = Math.max(1, Number.parseInt(data.pageSize || DEFAULT_PAGE_SIZE, 10));
  const totalPages = Math.max(1, Math.ceil(filteredRows.length / pageSize));
  data.page = Math.min(Math.max(1, data.page || 1), totalPages);

  const startIndex = (data.page - 1) * pageSize;
  const endIndex = startIndex + pageSize;
  const pageRows = filteredRows.slice(startIndex, endIndex);

  const emptyMessage = data.rows.length > 0
    ? 'No rows match the current filters'
    : 'No matching results';
  renderTableBody(viewer, pageRows, globalFilter, startIndex, emptyMessage);
  updateInlineRowCount(viewer, pageRows.length, filteredRows.length, data.rows.length, data.page, pageSize);
  updatePaginationControls(viewer, data.page, totalPages, filteredRows.length);
  renderQuickFilterChips(viewer, data);
  renderSelectedRowDetail(viewer, data);
  requestChartsLayoutRefresh(viewer, options);
}

/**
 * Update row count display
 */
function updateInlineRowCount(viewer, visiblePageRows, filteredTotal, totalRows, page, pageSize) {
  const countEl = viewer.querySelector('.results-row-count');
  if (!countEl) return;

  if (filteredTotal === 0) {
    countEl.textContent = totalRows === 0 ? '0 rows' : `0 / ${totalRows} rows`;
    return;
  }

  const start = (page - 1) * pageSize + 1;
  const end = start + visiblePageRows - 1;
  if (filteredTotal === totalRows) {
    countEl.textContent = `${start}-${end} of ${totalRows}`;
  } else {
    countEl.textContent = `${start}-${end} of ${filteredTotal} (filtered from ${totalRows})`;
  }
}

function updatePaginationControls(viewer, page, totalPages, filteredCount) {
  const pageInfo = viewer.querySelector('.results-page-info');
  const prevBtn = viewer.querySelector('.btn-results-page[data-direction="-1"]');
  const nextBtn = viewer.querySelector('.btn-results-page[data-direction="1"]');
  const pagination = viewer.querySelector('.results-pagination');
  if (!pageInfo || !prevBtn || !nextBtn || !pagination) return;

  pageInfo.textContent = `Page ${page} / ${totalPages}`;
  prevBtn.disabled = page <= 1;
  nextBtn.disabled = page >= totalPages;
  pagination.classList.toggle('hidden', filteredCount <= 0);
}

async function copyRowAsTsv(row) {
  try {
    await navigator.clipboard.writeText(row.join('\t'));
    logMessage('Row copied to clipboard (TSV).', 'info');
  } catch {
    logMessage('Could not copy row to clipboard.', 'warning');
  }
}

async function copyFilteredRowsAsTsv(data) {
  const rows = Array.isArray(data.filteredRows) ? data.filteredRows : [];
  if (rows.length === 0) {
    logMessage('No filtered rows to copy.', 'warning');
    return;
  }
  const lines = [data.headers.join('\t'), ...rows.map(row => row.join('\t'))];
  try {
    await navigator.clipboard.writeText(lines.join('\n'));
    logMessage(`Copied ${rows.length} filtered rows as TSV.`, 'success');
  } catch {
    logMessage('Could not copy filtered rows to clipboard.', 'warning');
  }
}

/**
 * Open folder containing output file
 */
export async function openOutputFolder(path) {
  if (!path) return;
  const success = await openFileLocation(path);
  if (!success) {
    try {
      await navigator.clipboard.writeText(path);
      logMessage('Path copied to clipboard', 'info');
    } catch {
      logMessage('Could not open folder or copy path to clipboard.', 'warning');
    }
  }
}

/**
 * Reset/hide inline results for a panel before starting a new run.
 */
export function resetInlineResults(panelId) {
  const viewerId = `${panelId}-results`;
  const viewer = document.getElementById(viewerId);
  const showResultsBtn = document.getElementById(`${panelId}-show-results`);
  const vizBtn = document.getElementById(`${panelId}-visualize-btn`);

  clearPanelResultsData(viewerId);

  if (viewer) {
    viewer.classList.add('hidden');
    viewer.classList.remove('fullscreen');
    delete viewer.dataset.outputPath;

    // Drop the multi-marker-set tab bar. It is only rebuilt on multi-marker
    // runs, so leaving it in place made a later single-marker run show stale
    // tabs pointing at the previous run's output sets.
    viewer.querySelector('.result-tabs')?.remove();

    const pathEl = viewer.querySelector('.results-viewer-path');
    if (pathEl) {
      pathEl.textContent = '';
      pathEl.removeAttribute('title');
    }

    const filterInput = viewer.querySelector('.results-filter-input');
    if (filterInput) filterInput.value = '';

    const countEl = viewer.querySelector('.results-row-count');
    if (countEl) countEl.textContent = '';

    const pageInfo = viewer.querySelector('.results-page-info');
    if (pageInfo) pageInfo.textContent = 'Page 1 / 1';

    const pageSize = viewer.querySelector('.results-page-size');
    if (pageSize) pageSize.value = String(DEFAULT_PAGE_SIZE);

    viewer.querySelectorAll('.btn-results-page').forEach(btn => {
      btn.disabled = true;
    });

    const thead = viewer.querySelector('thead');
    const tbody = viewer.querySelector('tbody');
    if (thead) thead.replaceChildren();
    if (tbody) tbody.replaceChildren();

    const quickFilters = viewer.querySelector('.results-quick-filters');
    const quickFilterChips = quickFilters?.querySelector('.results-filter-chips');
    if (quickFilterChips) quickFilterChips.replaceChildren();
    quickFilters?.classList.add('hidden');

    const detailPanel = viewer.querySelector('.results-row-detail');
    const detailBody = detailPanel?.querySelector('.results-row-detail-body');
    if (detailBody) detailBody.replaceChildren();
    detailPanel?.classList.add('hidden');
    viewer.querySelector('.results-viewer-layout')?.classList.remove('has-detail');
  }

  if (showResultsBtn) showResultsBtn.classList.add('hidden');
  if (vizBtn) {
    vizBtn.classList.add('hidden');
    vizBtn.dataset.ready = 'false';
  }
}

/**
 * Close all results viewers (legacy support)
 */
export function closeResultsModal() {
  document.querySelectorAll('.results-viewer').forEach(v => v.classList.add('hidden'));
}
