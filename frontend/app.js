// ============================================================================
// Pathotypr GUI - Enhanced Frontend Application (Tauri v2)
// Results Modal, File Chips, Stepped Progress, Smart Validation
// ============================================================================
// NOTE: This file is legacy and kept for reference only.
// Active frontend modules are under frontend/js/ and loaded via js/main.js in index.html.

// Global state
let currentDropTarget = null;
let isProcessing = false;
let consoleMessageCount = 0;
let resultsChart = null;
let lastOutputPath = null;
let progressInterval = null; // Store reference to progress interval for cleanup

// Step definitions for each command
const processSteps = {
  train: ['Loading data', 'Extracting k-mers', 'Training model', 'Saving'],
  predict: ['Loading model', 'Processing sequences', 'Classifying'],
  classify: ['Loading markers', 'Building k-mers', 'Classifying genomes'],
  'split-fastq': ['Indexing reference', 'Processing reads', 'Genotyping'],
  match: ['Indexing references', 'Processing reads', 'Matching']
};

// Wait for Tauri to be ready
document.addEventListener('DOMContentLoaded', async () => {
  await new Promise(resolve => setTimeout(resolve, 100));
  initApp();
});

function initApp() {
  console.log('Initializing Pathotypr GUI...');
  console.log('Tauri available:', typeof window.__TAURI__ !== 'undefined');

  if (typeof window.__TAURI__ === 'undefined') {
    console.error('Tauri API not available!');
    logMessage('Error: Tauri API not loaded. Please restart the application.', 'error');
    return;
  }

  setupNavigation();
  setupQuickActions();
  setupDropzones();
  setupTauriDragDrop();
  setupSaveButtons();
  setupRadioHandlers();
  setupForms();
  setupThemeToggle();
  setupRunButton();
  setupConsoleDrawer();
  setupResultsModal();
  setupDemoButton();
  setupDefaultButtons();
  document.getElementById('btn-clear-output').addEventListener('click', clearOutput);
  getVersion();
  loadSavedTheme();
  loadSavedFormState();
  logMessage('Pathotypr GUI ready', 'success');
}

// ============================================================================
// Tauri Native Drag & Drop
// ============================================================================

async function setupTauriDragDrop() {
  try {
    const { getCurrentWindow } = window.__TAURI__.window;
    const appWindow = getCurrentWindow();

    await appWindow.onDragDropEvent((event) => {
      if (event.payload.type === 'over') {
        document.querySelectorAll('.dropzone').forEach(dz => {
          dz.classList.add('drag-hover-global');
        });

        const pos = event.payload.position;
        if (pos) {
          const elementUnderCursor = document.elementFromPoint(pos.x, pos.y);
          const dropzone = elementUnderCursor?.closest('.dropzone');

          document.querySelectorAll('.dropzone').forEach(dz => {
            dz.classList.remove('drag-over');
          });

          if (dropzone) {
            currentDropTarget = dropzone;
            dropzone.classList.add('drag-over');
          }
        }
      } else if (event.payload.type === 'drop') {
        document.querySelectorAll('.dropzone').forEach(dz => {
          dz.classList.remove('drag-hover-global', 'drag-over');
        });

        const paths = event.payload.paths;
        if (paths && paths.length > 0) {
          const pos = event.payload.position;
          let targetDropzone = currentDropTarget;

          if (pos) {
            const elementUnderCursor = document.elementFromPoint(pos.x, pos.y);
            targetDropzone = elementUnderCursor?.closest('.dropzone') || currentDropTarget;
          }

          if (targetDropzone) {
            const targetId = targetDropzone.dataset.target;
            const isMultiple = targetDropzone.dataset.multiple === 'true';
            const validExtensions = targetDropzone.dataset.extensions?.split(',') || [];

            // Validate extensions
            const invalidFiles = paths.filter(p => !validateFileExtension(p, validExtensions));
            if (invalidFiles.length > 0) {
              showExtensionError(targetDropzone, validExtensions);
              return;
            }

            if (isMultiple) {
              setDropzoneFiles(targetDropzone, targetId, paths);
            } else {
              setDropzoneFile(targetDropzone, targetId, paths[0]);
            }

            targetDropzone.classList.add('drop-received');
            setTimeout(() => targetDropzone.classList.remove('drop-received'), 500);
            logMessage(`File loaded: ${paths[0].split(/[/\\]/).pop()}`, 'success');
          } else {
            logMessage(`Dropped ${paths.length} file(s). Please drop directly on a file input field.`, 'warning');
          }
        }
      } else if (event.payload.type === 'leave' || event.payload.type === 'cancel') {
        document.querySelectorAll('.dropzone').forEach(dz => {
          dz.classList.remove('drag-hover-global', 'drag-over');
        });
        currentDropTarget = null;
      }
    });

    logMessage('Drag & drop enabled', 'info');
  } catch (err) {
    console.warn('Tauri drag-drop setup failed:', err);
  }
}

// ============================================================================
// Smart Validation
// ============================================================================

function validateFileExtension(filePath, validExtensions) {
  if (!validExtensions || validExtensions.length === 0) return true;
  const fileName = filePath.toLowerCase();
  return validExtensions.some(ext => fileName.endsWith(ext.toLowerCase().trim()));
}

function showExtensionError(dropzone, validExtensions) {
  dropzone.classList.add('extension-error');
  logMessage(`Invalid file type. Expected: ${validExtensions.join(', ')}`, 'error');

  setTimeout(() => {
    dropzone.classList.remove('extension-error');
  }, 2000);
}

function formatFileSize(bytes) {
  if (bytes === 0) return '0 B';
  const k = 1024;
  const sizes = ['B', 'KB', 'MB', 'GB'];
  const i = Math.floor(Math.log(bytes) / Math.log(k));
  return parseFloat((bytes / Math.pow(k, i)).toFixed(1)) + ' ' + sizes[i];
}

// ============================================================================
// Console Drawer
// ============================================================================

function setupConsoleDrawer() {
  const drawer = document.getElementById('console-drawer');
  const backdrop = document.getElementById('drawer-backdrop');
  const toggleBtn = document.getElementById('btn-console-toggle');
  const closeBtn = document.getElementById('btn-close-drawer');
  const handle = document.getElementById('drawer-handle');

  if (toggleBtn) toggleBtn.addEventListener('click', toggleConsoleDrawer);
  if (closeBtn) closeBtn.addEventListener('click', closeConsoleDrawer);
  if (backdrop) backdrop.addEventListener('click', closeConsoleDrawer);

  // Drawer resize
  if (handle) {
    let isResizing = false;
    let startY, startHeight;

    handle.addEventListener('mousedown', (e) => {
      isResizing = true;
      startY = e.clientY;
      startHeight = drawer.offsetHeight;
      document.body.style.cursor = 'row-resize';
      e.preventDefault();
    });

    document.addEventListener('mousemove', (e) => {
      if (!isResizing) return;
      const diff = startY - e.clientY;
      const newHeight = Math.min(Math.max(startHeight + diff, 150), window.innerHeight * 0.7);
      drawer.style.height = `${newHeight}px`;
    });

    document.addEventListener('mouseup', () => {
      if (isResizing) {
        isResizing = false;
        document.body.style.cursor = '';
        localStorage.setItem('pathotypr-drawer-height', drawer.style.height);
      }
    });

    const savedHeight = localStorage.getItem('pathotypr-drawer-height');
    if (savedHeight) drawer.style.height = savedHeight;
  }

  document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape' && drawer.classList.contains('open')) {
      closeConsoleDrawer();
    }
  });
}

function toggleConsoleDrawer() {
  const drawer = document.getElementById('console-drawer');
  const backdrop = document.getElementById('drawer-backdrop');
  if (drawer.classList.contains('open')) {
    closeConsoleDrawer();
  } else {
    drawer.classList.add('open');
    backdrop.classList.add('visible');
    resetConsoleBadge();
  }
}

function closeConsoleDrawer() {
  document.getElementById('console-drawer').classList.remove('open');
  document.getElementById('drawer-backdrop').classList.remove('visible');
}

function openConsoleDrawer() {
  document.getElementById('console-drawer').classList.add('open');
  document.getElementById('drawer-backdrop').classList.add('visible');
  resetConsoleBadge();
}

function updateConsoleBadge() {
  const badge = document.getElementById('console-badge');
  const drawer = document.getElementById('console-drawer');
  if (!drawer.classList.contains('open')) {
    consoleMessageCount++;
    badge.textContent = consoleMessageCount > 99 ? '99+' : consoleMessageCount;
    badge.classList.add('has-messages');
  }
}

function resetConsoleBadge() {
  const badge = document.getElementById('console-badge');
  consoleMessageCount = 0;
  badge.classList.remove('has-messages');
  badge.textContent = '';
}

// ============================================================================
// Inline Results Viewer
// ============================================================================

// Store current results data for each panel
const panelResultsData = {};

function setupResultsModal() {
  // Setup inline results viewers for each panel
  setupInlineResultsViewer('predict-results');
  setupInlineResultsViewer('classify-results');
  setupInlineResultsViewer('splitfq-results');
  setupInlineResultsViewer('match-results');

  // Setup visualization for all tools
  setupToolVisualization('classify', 'Lineage Distribution', 'major_lineage');
  setupToolVisualization('predict', 'Prediction Distribution', 'prediction');
  setupToolVisualization('splitfq', 'Genotype Distribution', 'major_lineage');
  setupToolVisualization('match', 'Reference Distribution', 'best_match');
}

function setupInlineResultsViewer(viewerId) {
  const viewer = document.getElementById(viewerId);
  if (!viewer) return;

  const panelId = viewerId.replace('-results', '');
  const closeBtn = viewer.querySelector('.btn-close-results');
  const openFolderBtn = viewer.querySelector('.btn-open-folder');
  const fullscreenBtn = viewer.querySelector('.btn-fullscreen-results');
  const filterInput = viewer.querySelector('.results-filter-input');
  const showResultsBtn = document.getElementById(`${panelId}-show-results`);

  if (closeBtn) {
    closeBtn.addEventListener('click', () => {
      viewer.classList.add('hidden');
      viewer.classList.remove('fullscreen');
      // Show the "View Results" button
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
        // Exit fullscreen with animation
        viewer.classList.add('exiting-fullscreen');
        expandIcon?.classList.remove('hidden');
        collapseIcon?.classList.add('hidden');

        // Wait for animation to complete
        setTimeout(() => {
          viewer.classList.remove('fullscreen', 'exiting-fullscreen');
        }, 400);
      } else {
        // Enter fullscreen
        viewer.classList.add('fullscreen');
        expandIcon?.classList.add('hidden');
        collapseIcon?.classList.remove('hidden');
      }
    });
  }

  if (filterInput) {
    filterInput.addEventListener('input', (e) => {
      filterInlineResults(viewer, e.target.value);
    });
  }

  // Setup show results button
  if (showResultsBtn) {
    showResultsBtn.addEventListener('click', () => {
      viewer.classList.remove('hidden');
      showResultsBtn.classList.add('hidden');
      setTimeout(() => {
        viewer.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
      }, 100);
    });
  }
}

function showInlineResults(panelId, outputPath, tsvContent) {
  const viewerId = `${panelId}-results`;
  const viewer = document.getElementById(viewerId);
  const showResultsBtn = document.getElementById(`${panelId}-show-results`);
  console.log('[showInlineResults] Panel:', panelId, 'Viewer ID:', viewerId, 'Found:', !!viewer);

  if (!viewer) {
    console.warn('Results viewer not found:', viewerId);
    logMessage('Results viewer not found for panel: ' + panelId, 'warning');
    return;
  }

  // Store output path
  viewer.dataset.outputPath = outputPath;
  lastOutputPath = outputPath;

  // Update path display
  const pathEl = viewer.querySelector('.results-viewer-path');
  if (pathEl) {
    const fileName = outputPath.split(/[/\\]/).pop();
    pathEl.textContent = fileName;
    pathEl.title = outputPath;
  }

  // Reset filter and fullscreen
  const filterInput = viewer.querySelector('.results-filter-input');
  if (filterInput) filterInput.value = '';
  viewer.classList.remove('fullscreen');

  // Parse and display TSV
  const lines = tsvContent.trim().split('\n');
  console.log('[showInlineResults] TSV lines:', lines.length);

  if (lines.length === 0) {
    console.warn('[showInlineResults] No lines in TSV content');
    return;
  }

  const headers = lines[0].split('\t');
  const rows = lines.slice(1).map(line => line.split('\t'));
  console.log('[showInlineResults] Headers:', headers.length, 'Rows:', rows.length);

  // Detect which columns are numeric (if >50% of values are numbers)
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
  console.log('[showInlineResults] Numeric columns:', numericCols);

  // Store for filtering and sorting
  panelResultsData[viewerId] = { headers, rows, sortCol: -1, sortDir: 'asc', colFilters: {}, numericCols };

  // Populate table with sortable headers and column filters
  const thead = viewer.querySelector('thead');
  const tbody = viewer.querySelector('tbody');

  thead.innerHTML = `
    <tr>
      ${headers.map((h, i) => `
        <th data-col="${i}">
          <div class="th-content">
            <span>${escapeHtml(h)}</span>
            <svg class="sort-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
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

  renderTableBody(viewer, rows);

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
      const data = panelResultsData[viewerId];

      if (input.classList.contains('col-filter-num')) {
        // Numeric filter - store value and operator
        const opSelect = thead.querySelector(`.col-filter-op[data-col="${col}"]`);
        const op = opSelect?.value || '=';
        const numVal = input.value !== '' ? parseFloat(input.value) : null;
        data.colFilters[col] = { type: 'numeric', op, value: numVal };
      } else {
        // Text filter
        data.colFilters[col] = { type: 'text', value: input.value.toLowerCase() };
      }
      applyFilters(viewer);
    });
    input.addEventListener('click', (e) => e.stopPropagation());
  });

  // Setup numeric filter operator changes
  thead.querySelectorAll('.col-filter-op').forEach(select => {
    select.addEventListener('change', (e) => {
      e.stopPropagation();
      const col = parseInt(select.dataset.col);
      const data = panelResultsData[viewerId];
      const numInput = thead.querySelector(`.col-filter-num[data-col="${col}"]`);
      const numVal = numInput && numInput.value !== '' ? parseFloat(numInput.value) : null;
      data.colFilters[col] = { type: 'numeric', op: select.value, value: numVal };
      applyFilters(viewer);
    });
    select.addEventListener('click', (e) => e.stopPropagation());
  });

  // Update row count
  updateInlineRowCount(viewer, rows.length, rows.length);

  // Hide show results button, show viewer
  if (showResultsBtn) showResultsBtn.classList.add('hidden');
  viewer.classList.remove('hidden');

  // Show visualization button for all panels with results
  const vizBtn = document.getElementById(`${panelId}-visualize-btn`);
  if (vizBtn) {
    vizBtn.classList.remove('hidden');
    // Store data for visualization
    vizBtn.dataset.ready = 'true';
  }

  // Scroll to results
  setTimeout(() => {
    viewer.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  }, 100);
}

function renderTableBody(viewer, rows) {
  const tbody = viewer.querySelector('tbody');
  tbody.innerHTML = rows.map((row, idx) =>
    `<tr data-row-index="${idx}">${row.map(cell => `<td>${escapeHtml(cell)}</td>`).join('')}</tr>`
  ).join('');
}

function sortTable(viewer, colIndex) {
  const viewerId = viewer.id;
  const data = panelResultsData[viewerId];
  if (!data) return;

  // Toggle sort direction
  if (data.sortCol === colIndex) {
    data.sortDir = data.sortDir === 'asc' ? 'desc' : 'asc';
  } else {
    data.sortCol = colIndex;
    data.sortDir = 'asc';
  }

  // Sort rows
  const sortedRows = [...data.rows].sort((a, b) => {
    const valA = a[colIndex] || '';
    const valB = b[colIndex] || '';

    // Try numeric sort first
    const numA = parseFloat(valA);
    const numB = parseFloat(valB);
    if (!isNaN(numA) && !isNaN(numB)) {
      return data.sortDir === 'asc' ? numA - numB : numB - numA;
    }

    // Fall back to string sort
    const cmp = valA.localeCompare(valB, undefined, { numeric: true, sensitivity: 'base' });
    return data.sortDir === 'asc' ? cmp : -cmp;
  });

  // Update header styles
  viewer.querySelectorAll('th').forEach(th => {
    th.classList.remove('sort-asc', 'sort-desc');
    if (parseInt(th.dataset.col) === colIndex) {
      th.classList.add(data.sortDir === 'asc' ? 'sort-asc' : 'sort-desc');
    }
  });

  // Re-render and apply filters
  renderTableBody(viewer, sortedRows);
  applyFilters(viewer);
}

function applyFilters(viewer) {
  const viewerId = viewer.id;
  const data = panelResultsData[viewerId];
  if (!data) return;

  const globalFilter = viewer.querySelector('.results-filter-input')?.value.toLowerCase() || '';
  const tbody = viewer.querySelector('tbody');
  const rows = tbody.querySelectorAll('tr');
  let visibleCount = 0;

  rows.forEach((row) => {
    const cells = row.querySelectorAll('td');
    let matchesGlobal = !globalFilter;
    let matchesColFilters = true;

    cells.forEach((cell, colIdx) => {
      const text = cell.textContent.toLowerCase();
      cell.classList.remove('filter-match');

      // Check global filter
      if (globalFilter && text.includes(globalFilter)) {
        matchesGlobal = true;
        cell.classList.add('filter-match');
      }

      // Check column filter
      const colFilter = data.colFilters[colIdx];
      if (colFilter) {
        if (colFilter.type === 'numeric') {
          // Numeric filter with comparison operator
          if (colFilter.value !== null) {
            const cellNum = parseFloat(cell.textContent);
            if (!isNaN(cellNum)) {
              let matches = false;
              switch (colFilter.op) {
                case '=': matches = cellNum === colFilter.value; break;
                case '>': matches = cellNum > colFilter.value; break;
                case '<': matches = cellNum < colFilter.value; break;
                case '>=': matches = cellNum >= colFilter.value; break;
                case '<=': matches = cellNum <= colFilter.value; break;
              }
              if (!matches) matchesColFilters = false;
            } else {
              // Non-numeric cell doesn't match numeric filter
              matchesColFilters = false;
            }
          }
        } else if (colFilter.type === 'text') {
          // Text filter
          if (colFilter.value && !text.includes(colFilter.value)) {
            matchesColFilters = false;
          }
        }
      }
    });

    if (matchesGlobal && matchesColFilters) {
      row.classList.remove('filtered-out');
      visibleCount++;
    } else {
      row.classList.add('filtered-out');
    }
  });

  updateInlineRowCount(viewer, visibleCount, data.rows.length);
}

function filterInlineResults(viewer, searchTerm) {
  // Use the unified filter function
  applyFilters(viewer);
}

function updateInlineRowCount(viewer, visible, total) {
  const countEl = viewer.querySelector('.results-row-count');
  if (countEl) {
    if (visible === total) {
      countEl.textContent = `${total} rows`;
    } else {
      countEl.textContent = `${visible} / ${total}`;
    }
  }
}

// ============================================================================
// Tool Visualization (Generic)
// ============================================================================

// Store charts for each tool
const toolCharts = {};

function setupToolVisualization(toolId, title, columnName) {
  const vizBtn = document.getElementById(`${toolId}-visualize-btn`);
  const vizPanel = document.getElementById(`${toolId}-visualization`);
  const closeBtn = vizPanel?.querySelector('.btn-close-viz');
  const fullscreenBtn = vizPanel?.querySelector('.btn-fullscreen-viz');

  // Initialize charts storage for this tool
  toolCharts[toolId] = { donut: null, bar: null };

  if (vizBtn) {
    vizBtn.addEventListener('click', () => {
      if (vizBtn.dataset.ready !== 'true') return;
      showToolVisualization(toolId, columnName);
    });
  }

  if (closeBtn) {
    closeBtn.addEventListener('click', () => {
      vizPanel.classList.add('hidden');
      vizPanel.classList.remove('fullscreen');
      vizBtn?.classList.remove('hidden');
    });
  }

  if (fullscreenBtn) {
    fullscreenBtn.addEventListener('click', () => {
      const expandIcon = fullscreenBtn.querySelector('.icon-expand');
      const collapseIcon = fullscreenBtn.querySelector('.icon-collapse');

      if (vizPanel.classList.contains('fullscreen')) {
        // Exit fullscreen - recreate charts at original size
        vizPanel.classList.remove('fullscreen');
        expandIcon?.classList.remove('hidden');
        collapseIcon?.classList.add('hidden');

        // Recreate charts after DOM settles
        const data = panelResultsData[`${toolId}-results`];
        if (data) {
          setTimeout(() => {
            const counts = parseToolData(data, columnName);
            createToolCharts(toolId, counts);
          }, 100);
        }
      } else {
        // Enter fullscreen
        vizPanel.classList.add('fullscreen');
        expandIcon?.classList.add('hidden');
        collapseIcon?.classList.remove('hidden');

        // Recreate charts at fullscreen size
        const data = panelResultsData[`${toolId}-results`];
        if (data) {
          setTimeout(() => {
            const counts = parseToolData(data, columnName);
            createToolCharts(toolId, counts);
          }, 100);
        }
      }
    });
  }
}

function showToolVisualization(toolId, columnName) {
  const vizPanel = document.getElementById(`${toolId}-visualization`);
  const vizBtn = document.getElementById(`${toolId}-visualize-btn`);
  const data = panelResultsData[`${toolId}-results`];

  if (!data || !vizPanel) return;

  // Parse the data to extract distribution
  const counts = parseToolData(data, columnName);

  if (Object.keys(counts).length === 0) {
    logMessage('No data found for visualization', 'warning');
    return;
  }

  // Hide button, show panel
  vizBtn?.classList.add('hidden');
  vizPanel.classList.remove('hidden');

  // Create charts
  createToolCharts(toolId, counts);

  // Scroll to visualization
  setTimeout(() => {
    vizPanel.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  }, 100);
}

function parseToolData(data, primaryColumn) {
  const counts = {};
  const headers = data.headers;
  const rows = data.rows;

  // Find the primary column specifically
  let colIdx = headers.findIndex(h =>
    h.toLowerCase() === primaryColumn.toLowerCase()
  );

  // If not found, try other common variations based on context
  if (colIdx === -1) {
    const variations = [
      'classification', 'lineage', 'prediction', 'best_match',
      'match', 'result', 'genotype', 'type', 'category'
    ];
    colIdx = headers.findIndex(h =>
      variations.some(v => h.toLowerCase().includes(v))
    );
  }

  // If still not found, use the last column
  if (colIdx === -1) colIdx = headers.length - 1;

  console.log('[parseToolData] Using column:', headers[colIdx], 'at index:', colIdx);

  rows.forEach(row => {
    const value = row[colIdx]?.trim();
    if (value && value !== '' && value !== 'N/A' && value !== 'Unknown') {
      counts[value] = (counts[value] || 0) + 1;
    }
  });

  return counts;
}

function createToolCharts(toolId, dataCounts) {
  // Destroy existing charts for this tool
  if (toolCharts[toolId]?.donut) toolCharts[toolId].donut.destroy();
  if (toolCharts[toolId]?.bar) toolCharts[toolId].bar.destroy();

  const labels = Object.keys(dataCounts);
  const values = Object.values(dataCounts);
  const total = values.reduce((a, b) => a + b, 0);

  // Generate colors based on labels
  const colors = generateChartColors(labels);

  // Get theme-aware text color
  const textColor = getComputedStyle(document.documentElement).getPropertyValue('--text-primary').trim() || '#18181b';
  const gridColor = getComputedStyle(document.documentElement).getPropertyValue('--border').trim() || '#e4e4e7';

  // Create donut chart
  const donutCtx = document.getElementById(`${toolId}-donut-chart`)?.getContext('2d');
  if (donutCtx) {
    toolCharts[toolId].donut = new Chart(donutCtx, {
      type: 'doughnut',
      data: {
        labels: labels,
        datasets: [{
          data: values,
          backgroundColor: colors,
          borderColor: 'transparent',
          borderWidth: 0,
          hoverOffset: 8
        }]
      },
      options: {
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
                const pct = ((ctx.raw / total) * 100).toFixed(1);
                return ` ${ctx.label}: ${ctx.raw} (${pct}%)`;
              }
            }
          }
        }
      }
    });
  }

  // Create horizontal bar chart
  const barCtx = document.getElementById(`${toolId}-bar-chart`)?.getContext('2d');
  if (barCtx) {
    // Sort by value descending
    const sorted = labels.map((l, i) => ({ label: l, value: values[i], color: colors[i] }))
      .sort((a, b) => b.value - a.value);

    toolCharts[toolId].bar = new Chart(barCtx, {
      type: 'bar',
      data: {
        labels: sorted.map(d => d.label),
        datasets: [{
          data: sorted.map(d => d.value),
          backgroundColor: sorted.map(d => d.color),
          borderRadius: 6,
          borderSkipped: false
        }]
      },
      options: {
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
                const pct = ((ctx.raw / total) * 100).toFixed(1);
                return ` Count: ${ctx.raw} (${pct}%)`;
              }
            }
          }
        },
        scales: {
          x: {
            beginAtZero: true,
            grid: {
              color: gridColor + '40'
            },
            ticks: {
              color: textColor,
              font: { size: 11 }
            }
          },
          y: {
            grid: { display: false },
            ticks: {
              color: textColor,
              font: { size: 12, weight: '500' }
            }
          }
        }
      }
    });
  }

  // Create legend
  createToolLegend(toolId, labels, values, colors, total);
}

// Lineage-specific colors
const lineageColors = {
  'A1': '#d1ae00',
  'A2': '#8ef5c8',
  'A3': '#73c2ff',
  'A4': '#ff9cdb',
  'L1': '#ff3091',
  'L2': '#001aff',
  'L3': '#8a0bd2',
  'L4': '#ff0000',
  'L5': '#995200',
  'L6': '#1eb040',
  'L7': '#fbff00',
  'L8': '#ff9d00',
  'L9': '#37ff30',
  'L10': '#8fbda1'
};

const fallbackColors = [
  '#14b8a6', '#8b5cf6', '#f59e0b', '#ef4444', '#3b82f6',
  '#ec4899', '#10b981', '#6366f1', '#f97316', '#06b6d4',
  '#84cc16', '#a855f7', '#22c55e', '#e11d48', '#0ea5e9'
];

function generateChartColors(labels) {
  let fallbackIdx = 0;

  return labels.map(label => {
    // Check if label contains any of the lineage names
    const upperLabel = label.toUpperCase();

    // Check for exact matches first (L10 before L1)
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

function createToolLegend(toolId, labels, values, colors, total) {
  const legendEl = document.getElementById(`${toolId}-chart-legend`);
  if (!legendEl) return;

  legendEl.innerHTML = labels.map((label, i) => {
    const pct = ((values[i] / total) * 100).toFixed(1);
    return `
      <div class="legend-item" data-index="${i}" data-tool="${toolId}">
        <span class="legend-color" style="background: ${colors[i]}"></span>
        <span class="legend-label">${label}</span>
        <span class="legend-count">(${values[i]} - ${pct}%)</span>
      </div>
    `;
  }).join('');

  // Add hover interactivity
  legendEl.querySelectorAll('.legend-item').forEach(item => {
    item.addEventListener('mouseenter', () => {
      const idx = parseInt(item.dataset.index);
      const tool = item.dataset.tool;
      highlightChartSegment(tool, idx);
    });
    item.addEventListener('mouseleave', () => {
      const tool = item.dataset.tool;
      resetChartHighlight(tool);
    });
  });
}

function highlightChartSegment(toolId, index) {
  const charts = toolCharts[toolId];
  if (charts?.donut) {
    charts.donut.setActiveElements([{ datasetIndex: 0, index }]);
    charts.donut.update();
  }
  if (charts?.bar) {
    charts.bar.setActiveElements([{ datasetIndex: 0, index }]);
    charts.bar.update();
  }
}

function resetChartHighlight(toolId) {
  const charts = toolCharts[toolId];
  if (charts?.donut) {
    charts.donut.setActiveElements([]);
    charts.donut.update();
  }
  if (charts?.bar) {
    charts.bar.setActiveElements([]);
    charts.bar.update();
  }
}

function resizeCharts(toolId) {
  const charts = toolCharts[toolId];
  if (charts?.donut) {
    charts.donut.resize();
    charts.donut.update('none');
  }
  if (charts?.bar) {
    charts.bar.resize();
    charts.bar.update('none');
  }
}

function escapeHtml(text) {
  const div = document.createElement('div');
  div.textContent = text;
  return div.innerHTML;
}

async function openOutputFolder(path) {
  const targetPath = path || lastOutputPath;
  if (!targetPath) return;
  try {
    const { invoke } = window.__TAURI__.core;
    await invoke('open_file_location', { path: targetPath });
  } catch (err) {
    // Fallback: copy path to clipboard
    navigator.clipboard.writeText(targetPath);
    logMessage('Path copied to clipboard', 'info');
  }
}

// Legacy modal functions (kept for compatibility)
function openResultsModal(outputPath, tsvContent) {
  // Determine which panel we're in based on active panel
  const activePanel = document.querySelector('.panel.active');
  if (!activePanel) return;

  const panelId = activePanel.id.replace('panel-', '');
  showInlineResults(panelId, outputPath, tsvContent);
}

function closeResultsModal() {
  // Close all inline results viewers
  document.querySelectorAll('.results-viewer').forEach(v => v.classList.add('hidden'));
}

// ============================================================================
// Theme Toggle
// ============================================================================

function setupThemeToggle() {
  const themeToggle = document.getElementById('theme-toggle');
  if (themeToggle) themeToggle.addEventListener('click', toggleTheme);
}

function toggleTheme() {
  const html = document.documentElement;
  const currentTheme = html.getAttribute('data-theme');
  const newTheme = currentTheme === 'dark' ? 'light' : 'dark';
  html.setAttribute('data-theme', newTheme);
  updateThemeUI(newTheme);
  localStorage.setItem('pathotypr-theme', newTheme);
}

function updateThemeUI(theme) {
  const iconEl = document.getElementById('theme-icon');
  const textEl = document.getElementById('theme-text');

  if (theme === 'dark') {
    iconEl.innerHTML = `<path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/>`;
    textEl.textContent = 'Dark Mode';
  } else {
    iconEl.innerHTML = `<circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/><line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/><line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/>`;
    textEl.textContent = 'Light Mode';
  }
}

function loadSavedTheme() {
  const savedTheme = localStorage.getItem('pathotypr-theme') || 'light';
  document.documentElement.setAttribute('data-theme', savedTheme);
  updateThemeUI(savedTheme);
}

// ============================================================================
// Default Parameters
// ============================================================================

function setupDefaultButtons() {
  document.querySelectorAll('.btn-defaults').forEach(btn => {
    btn.addEventListener('click', () => {
      const tool = btn.dataset.tool;
      setDefaultParameters(tool);
    });
  });
}

function setDefaultParameters(tool) {
  const defaults = {
    train: {
      'train-kmer': 21,
      'train-split': 0.2,
      'train-threads': ''  // Auto
    },
    predict: {
      'predict-threads': ''  // Auto
    },
    classify: {
      'classify-kmer': 21,
      'classify-threads': '',  // Auto
      'classify-nested': false
    },
    splitfq: {
      'splitfq-min-depth': 10,
      'splitfq-min-alt': 95,
      'splitfq-threads': '',  // Auto
      'splitfq-paired': false,
      'splitfq-nested': false
    },
    match: {
      'match-kmer': 21,
      'match-threads': ''  // Auto
    }
  };

  const toolDefaults = defaults[tool];
  if (!toolDefaults) return;

  Object.entries(toolDefaults).forEach(([id, value]) => {
    const element = document.getElementById(id);
    if (!element) return;

    if (element.type === 'checkbox') {
      element.checked = value;
    } else {
      element.value = value;
    }
  });

  logMessage(`Default parameters set for ${tool}`, 'info');
}

// ============================================================================
// Demo Data
// ============================================================================

function setupDemoButton() {
  const btn = document.getElementById('btn-load-demo');
  if (btn) btn.addEventListener('click', loadDemoData);
}

function loadDemoData() {
  // Fill Train form with demo values
  document.getElementById('train-kmer').value = 21;
  document.getElementById('train-split').value = 0.2;
  document.getElementById('train-threads').value = 4;

  // Fill Classify form
  document.getElementById('classify-kmer').value = 21;
  document.getElementById('classify-threads').value = 4;
  document.getElementById('classify-nested').checked = true;

  // Fill Split-FASTQ form
  document.getElementById('splitfq-min-depth').value = 10;
  document.getElementById('splitfq-min-alt').value = 95;
  document.getElementById('splitfq-threads').value = 4;

  // Fill Match form
  document.getElementById('match-kmer').value = 31;
  document.getElementById('match-threads').value = 4;

  logMessage('Demo parameters loaded. Select your input files to run.', 'success');
  saveFormState();
}

// ============================================================================
// Form State Persistence
// ============================================================================

function saveFormState() {
  const state = {
    'train-kmer': document.getElementById('train-kmer')?.value,
    'train-split': document.getElementById('train-split')?.value,
    'train-threads': document.getElementById('train-threads')?.value,
    'classify-kmer': document.getElementById('classify-kmer')?.value,
    'classify-threads': document.getElementById('classify-threads')?.value,
    'classify-nested': document.getElementById('classify-nested')?.checked,
    'splitfq-min-depth': document.getElementById('splitfq-min-depth')?.value,
    'splitfq-min-alt': document.getElementById('splitfq-min-alt')?.value,
    'splitfq-threads': document.getElementById('splitfq-threads')?.value,
    'splitfq-paired': document.getElementById('splitfq-paired')?.checked,
    'splitfq-nested': document.getElementById('splitfq-nested')?.checked,
    'match-kmer': document.getElementById('match-kmer')?.value,
    'match-threads': document.getElementById('match-threads')?.value
  };
  localStorage.setItem('pathotypr-form-state', JSON.stringify(state));
}

function loadSavedFormState() {
  try {
    const state = JSON.parse(localStorage.getItem('pathotypr-form-state'));
    if (!state) return;

    Object.entries(state).forEach(([id, value]) => {
      const el = document.getElementById(id);
      if (el) {
        if (el.type === 'checkbox') {
          el.checked = value;
        } else {
          el.value = value || el.defaultValue;
        }
      }
    });
  } catch (e) {
    console.warn('Failed to load form state:', e);
  }
}

// Save form state on input change
document.addEventListener('change', (e) => {
  if (e.target.closest('form')) {
    saveFormState();
  }
});

// ============================================================================
// Run Button Handler
// ============================================================================

function setupRunButton() {
  const runBtn = document.getElementById('btn-run');
  if (runBtn) runBtn.addEventListener('click', handleRunClick);
}

async function handleRunClick() {
  console.log('[Run] Run button clicked - isProcessing:', isProcessing);

  if (isProcessing) {
    // Cancel the running task
    console.log('[Run] Attempting to cancel running task');
    logMessage('Cancelling task...', 'warning');
    try {
      const cancelled = await tauriInvoke('cancel_task');
      if (cancelled) {
        logMessage('Cancel signal sent. Waiting for task to stop...', 'info');
      }
    } catch (err) {
      console.error('[Run] Cancel error:', err);
    }
    return;
  }

  const activePanel = document.querySelector('.panel.active');
  if (!activePanel) return;

  if (activePanel.id === 'panel-home') {
    logMessage('Select a tool from the sidebar or quick actions to run', 'warning');
    return;
  }

  const form = activePanel.querySelector('form');
  if (form) {
    console.log('[Run] Dispatching submit event to form');
    form.dispatchEvent(new Event('submit', { cancelable: true, bubbles: true }));
  }
}

// ============================================================================
// Tauri API Helpers
// ============================================================================

async function tauriInvoke(cmd, args = {}) {
  try {
    return await window.__TAURI__.core.invoke(cmd, args);
  } catch (err) {
    console.error(`Tauri invoke error (${cmd}):`, err);
    throw err;
  }
}

async function openFileDialog(options = {}) {
  try {
    return await window.__TAURI__.dialog.open(options);
  } catch (err) {
    console.error('Open dialog error:', err);
    throw err;
  }
}

async function saveFileDialog(options = {}) {
  try {
    return await window.__TAURI__.dialog.save(options);
  } catch (err) {
    console.error('Save dialog error:', err);
    throw err;
  }
}

async function readTextFile(path) {
  try {
    // Use our custom Tauri command to read files (more reliable than fs plugin)
    const { invoke } = window.__TAURI__.core;
    return await invoke('read_text_file', { path });
  } catch (err) {
    console.error('[readTextFile] Error reading file:', path, err);
    return null;
  }
}

async function getFileMetadata(path) {
  try {
    const { stat } = window.__TAURI__.fs;
    return await stat(path);
  } catch (err) {
    return null;
  }
}

// ============================================================================
// Navigation
// ============================================================================

const panelInfo = {
  'home': {
    iconPath: '<path d="M3 9l9-7 9 7v11a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2z"/><polyline points="9 22 9 12 15 12 15 22"/>',
    title: 'Home',
    description: 'Welcome to Pathotypr - Select a tool to get started'
  },
  'train': {
    iconPath: '<path d="M12 20V10"/><path d="M18 20V4"/><path d="M6 20v-4"/>',
    title: 'Train Model',
    description: 'Build a Random Forest classifier from labeled sequences'
  },
  'predict': {
    iconPath: '<circle cx="12" cy="12" r="10"/><path d="M12 16v-4"/><path d="M12 8h.01"/>',
    title: 'Predict',
    description: 'Classify sequences using a trained model'
  },
  'classify': {
    iconPath: '<path d="M9 3H5a2 2 0 0 0-2 2v4m6-6h10a2 2 0 0 1 2 2v4M9 3v18m0 0h10a2 2 0 0 0 2-2v-4M9 21H5a2 2 0 0 1-2-2v-4"/>',
    title: 'Classify',
    description: 'Genotype genomes based on SNP markers'
  },
  'split-fastq': {
    iconPath: '<path d="M14.5 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V7.5L14.5 2z"/><polyline points="14 2 14 8 20 8"/><line x1="16" y1="13" x2="8" y2="13"/><line x1="16" y1="17" x2="8" y2="17"/>',
    title: 'Split FASTQ',
    description: 'Alignment-free genotyping from raw reads'
  },
  'match': {
    iconPath: '<circle cx="11" cy="11" r="8"/><path d="m21 21-4.35-4.35"/>',
    title: 'Match Reference',
    description: 'Find best matching reference for each sample'
  }
};

function setupNavigation() {
  document.querySelectorAll('.nav-item').forEach(item => {
    item.addEventListener('click', () => {
      navigateToPanel(item.dataset.panel);
    });
  });
}

function navigateToPanel(panelId) {
  document.querySelectorAll('.nav-item').forEach(nav => {
    nav.classList.toggle('active', nav.dataset.panel === panelId);
  });

  document.querySelectorAll('.panel').forEach(panel => {
    panel.classList.toggle('active', panel.id === `panel-${panelId}`);
  });

  updateHeader(panelId);
  updateCategoryAccent(panelId);
}

function updateCategoryAccent(panelId) {
  // Map panels to categories
  const categoryMap = {
    'home': null,
    'train': 'ml',
    'predict': 'ml',
    'classify': 'genotyping',
    'split-fastq': 'genotyping',
    'match': 'utils'
  };

  const category = categoryMap[panelId];
  if (category) {
    document.body.setAttribute('data-category', category);
  } else {
    document.body.removeAttribute('data-category');
  }
}

function updateHeader(panelId) {
  const info = panelInfo[panelId];
  if (info) {
    document.getElementById('header-icon').innerHTML = info.iconPath;
    document.getElementById('panel-title-text').textContent = info.title;
    document.getElementById('panel-description').textContent = info.description;
  }
}

// ============================================================================
// Quick Actions
// ============================================================================

function setupQuickActions() {
  document.querySelectorAll('.action-card[data-goto]').forEach(card => {
    card.addEventListener('click', () => {
      navigateToPanel(card.dataset.goto);
    });
  });
}

// ============================================================================
// Dropzones
// ============================================================================

function setupDropzones() {
  document.querySelectorAll('.dropzone').forEach(dropzone => {
    const targetId = dropzone.dataset.target;
    const isMultiple = dropzone.dataset.multiple === 'true';
    const validExtensions = dropzone.dataset.extensions?.split(',') || [];
    const removeBtn = dropzone.querySelector('.file-remove');

    dropzone.addEventListener('mouseenter', () => { currentDropTarget = dropzone; });
    dropzone.addEventListener('mouseleave', () => {
      if (currentDropTarget === dropzone) currentDropTarget = null;
    });

    dropzone.addEventListener('click', async (e) => {
      if (e.target.classList.contains('file-remove') || e.target.closest('.file-chip')) return;

      dropzone.style.transform = 'scale(0.98)';
      setTimeout(() => dropzone.style.transform = '', 100);

      try {
        const selected = await openFileDialog({ multiple: isMultiple });
        if (selected) {
          const paths = Array.isArray(selected) ? selected : [selected];

          // Validate extensions
          const invalidFiles = paths.filter(p => !validateFileExtension(p, validExtensions));
          if (invalidFiles.length > 0) {
            showExtensionError(dropzone, validExtensions);
            return;
          }

          if (isMultiple) {
            setDropzoneFiles(dropzone, targetId, paths);
          } else {
            setDropzoneFile(dropzone, targetId, paths[0]);
          }
        }
      } catch (err) {
        logMessage('Error selecting file: ' + err, 'error');
      }
    });

    if (removeBtn) {
      removeBtn.addEventListener('click', (e) => {
        e.stopPropagation();
        clearDropzone(dropzone, targetId);
        logMessage('File removed', 'info');
      });
    }
  });
}

async function setDropzoneFile(dropzone, targetId, filePath) {
  const input = document.getElementById(targetId);
  const fileDisplay = dropzone.querySelector('.dropzone-file');
  const fileName = fileDisplay.querySelector('.file-name');
  const fileMeta = fileDisplay.querySelector('.file-meta');

  input.value = filePath;
  // Dispatch event for listeners
  input.dispatchEvent(new CustomEvent('change', { detail: { path: filePath } }));

  const name = filePath.split(/[/\\]/).pop();
  fileName.textContent = name;

  // Get file metadata
  const metadata = await getFileMetadata(filePath);
  if (metadata) {
    const sizeStr = formatFileSize(metadata.size || 0);
    fileMeta.innerHTML = `<span class="meta-badge">${sizeStr}</span>`;
  } else {
    fileMeta.innerHTML = '';
  }

  fileDisplay.classList.add('has-file');
  logMessage(`Selected: ${name}`, 'success');
}

async function setDropzoneFiles(dropzone, targetId, filePaths, replace = false) {
  const input = document.getElementById(targetId);
  const fileDisplay = dropzone.querySelector('.dropzone-file');
  const fileName = fileDisplay.querySelector('.file-name');
  const fileMeta = fileDisplay.querySelector('.file-meta');
  const chipsContainer = dropzone.querySelector('.file-chips');

  // Get existing files and merge with new ones (avoiding duplicates)
  let allFiles = filePaths;
  if (!replace && input?.dataset.files) {
    try {
      const existingFiles = JSON.parse(input.dataset.files);
      // Merge: add new files that aren't already in the list
      const newFiles = filePaths.filter(p => !existingFiles.includes(p));
      allFiles = [...existingFiles, ...newFiles];
    } catch (e) {
      // If parsing fails, just use the new files
      allFiles = filePaths;
    }
  }

  input.value = allFiles.join(';');
  input.dataset.files = JSON.stringify(allFiles);

  fileName.textContent = `${allFiles.length} files selected`;
  fileMeta.innerHTML = '';
  fileDisplay.classList.add('has-file');

  // Create file chips
  if (chipsContainer) {
    chipsContainer.innerHTML = '';
    for (const path of allFiles) {
      const name = path.split(/[/\\]/).pop();
      const metadata = await getFileMetadata(path);
      const sizeStr = metadata ? formatFileSize(metadata.size || 0) : '';

      const chip = document.createElement('div');
      chip.className = 'file-chip';
      chip.innerHTML = `
        <span class="chip-name" title="${name}">${name}</span>
        ${sizeStr ? `<span class="chip-size">${sizeStr}</span>` : ''}
        <button type="button" class="chip-remove" data-path="${path}">&times;</button>
      `;
      chipsContainer.appendChild(chip);
    }

    // Setup chip remove handlers
    chipsContainer.querySelectorAll('.chip-remove').forEach(btn => {
      btn.addEventListener('click', (e) => {
        e.stopPropagation();
        const pathToRemove = btn.dataset.path;
        const remaining = allFiles.filter(p => p !== pathToRemove);
        if (remaining.length > 0) {
          setDropzoneFiles(dropzone, targetId, remaining, true); // Replace mode for removal
        } else {
          clearDropzone(dropzone, targetId);
        }
      });
    });
  }

  const addedCount = allFiles.length - (replace ? 0 : (allFiles.length - filePaths.length));
  if (addedCount > 0) {
    logMessage(`Added ${filePaths.length} file(s). Total: ${allFiles.length}`, 'success');
  }
}

function clearDropzone(dropzone, targetId) {
  const input = document.getElementById(targetId);
  const fileDisplay = dropzone.querySelector('.dropzone-file');
  const chipsContainer = dropzone.querySelector('.file-chips');

  input.value = '';
  delete input.dataset.files;
  // Dispatch event for listeners
  input.dispatchEvent(new CustomEvent('change', { detail: { path: null } }));
  fileDisplay.classList.remove('has-file');

  if (chipsContainer) {
    chipsContainer.innerHTML = '';
  }
}

// ============================================================================
// Save Buttons
// ============================================================================

function setupSaveButtons() {
  document.querySelectorAll('.btn-save').forEach(btn => {
    btn.addEventListener('click', async () => {
      const targetId = btn.dataset.target;
      const defaultName = btn.dataset.default || 'output';

      try {
        const selected = await saveFileDialog({ defaultPath: defaultName });
        if (selected) {
          document.getElementById(targetId).value = selected;
          logMessage(`Output: ${selected.split(/[/\\]/).pop()}`, 'success');
        }
      } catch (err) {
        logMessage('Error selecting save location: ' + err, 'error');
      }
    });
  });
}

// ============================================================================
// Radio Button Handlers
// ============================================================================

function setupRadioHandlers() {
  document.querySelectorAll('input[name="classify-input-method"]').forEach(radio => {
    radio.addEventListener('change', () => {
      const isFile = radio.value === 'file';
      document.getElementById('classify-input-file-group').classList.toggle('hidden', !isFile);
      document.getElementById('classify-input-list-group').classList.toggle('hidden', isFile);
    });
  });

  document.querySelectorAll('input[name="splitfq-input-method"]').forEach(radio => {
    radio.addEventListener('change', () => {
      const isFiles = radio.value === 'files';
      document.getElementById('splitfq-input-files-group').classList.toggle('hidden', !isFiles);
      document.getElementById('splitfq-input-list-group').classList.toggle('hidden', isFiles);
    });
  });

  document.querySelectorAll('input[name="match-input-method"]').forEach(radio => {
    radio.addEventListener('change', () => {
      const isFiles = radio.value === 'files';
      document.getElementById('match-input-files-group').classList.toggle('hidden', !isFiles);
      document.getElementById('match-input-list-group').classList.toggle('hidden', isFiles);
    });
  });

  // Match Reference method selector (FASTA vs Pre-built Index)
  document.querySelectorAll('input[name="match-ref-method"]').forEach(radio => {
    radio.addEventListener('change', () => {
      const isFasta = radio.value === 'fasta';
      document.getElementById('match-fasta-group').classList.toggle('hidden', !isFasta);
      document.getElementById('match-index-group').classList.toggle('hidden', isFasta);
    });
  });

  // Build Index button
  const buildIndexBtn = document.getElementById('match-build-index-btn');
  if (buildIndexBtn) {
    buildIndexBtn.addEventListener('click', handleBuildIndex);
  }

  // Enable/disable build index button based on reference file selection
  setupBuildIndexButton();
}

// ============================================================================
// Form Validation
// ============================================================================

function validateForm(fields, formName) {
  const missing = [];

  for (const [fieldId, label] of Object.entries(fields)) {
    const el = document.getElementById(fieldId);
    if (!el || !el.value || el.value.trim() === '') {
      missing.push(label);
      highlightMissingField(fieldId);
    }
  }

  if (missing.length > 0) {
    logMessage(`Cannot run ${formName}: Missing required fields`, 'error');
    missing.forEach(field => {
      logMessage(`  - ${field}`, 'error');
    });
    openConsoleDrawer();
    return false;
  }
  return true;
}

function highlightMissingField(fieldId) {
  const el = document.getElementById(fieldId);
  if (!el) return;

  const dropzone = el.closest('.form-group')?.querySelector('.dropzone');
  const inputWrapper = el.closest('.file-input-wrapper');
  const target = dropzone || inputWrapper || el;

  target.classList.add('field-missing');
  setTimeout(() => target.classList.remove('field-missing'), 3000);
}

// ============================================================================
// Form Submissions
// ============================================================================

function setupForms() {
  document.getElementById('form-train').addEventListener('submit', handleTrainSubmit);
  document.getElementById('form-predict').addEventListener('submit', handlePredictSubmit);
  document.getElementById('form-classify').addEventListener('submit', handleClassifySubmit);
  document.getElementById('form-split-fastq').addEventListener('submit', handleSplitFastqSubmit);
  document.getElementById('form-match').addEventListener('submit', handleMatchSubmit);
}

async function handleTrainSubmit(e) {
  e.preventDefault();

  if (!validateForm({
    'train-input': 'Training FASTA file',
    'train-output': 'Output model file location'
  }, 'Train')) return;

  const input = document.getElementById('train-input').value;
  const output = document.getElementById('train-output').value;
  const kmerSize = parseInt(document.getElementById('train-kmer').value) || 21;
  const testSplit = parseFloat(document.getElementById('train-split').value) || 0.2;
  const threads = parseInt(document.getElementById('train-threads').value) || null;

  startProgress('train');
  logMessage('Starting model training...', 'info');
  logMessage(`Input: ${makeClickablePath(input)}`, 'info');

  try {
    const result = await tauriInvoke('run_train', {
      params: { input, output, kmer_size: kmerSize, test_split: testSplit, threads }
    });
    handleResult(result);
  } catch (err) {
    logMessage('Error: ' + err, 'error');
  } finally {
    stopProgress();
  }
}

async function handlePredictSubmit(e) {
  e.preventDefault();

  if (!validateForm({
    'predict-input': 'FASTA file to classify',
    'predict-model': 'Trained model file (.zst)',
    'predict-output': 'Output results file location'
  }, 'Predict')) return;

  const input = document.getElementById('predict-input').value;
  const model = document.getElementById('predict-model').value;
  const output = document.getElementById('predict-output').value;
  const threads = parseInt(document.getElementById('predict-threads').value) || null;

  startProgress('predict');
  logMessage('Starting lineage prediction...', 'info');

  try {
    const result = await tauriInvoke('run_predict', {
      params: { input, model, output, threads }
    });
    await handleResultWithModal(result, output, 'predict');
  } catch (err) {
    logMessage('Error: ' + err, 'error');
  } finally {
    stopProgress();
  }
}

async function handleClassifySubmit(e) {
  e.preventDefault();

  const inputMethod = document.querySelector('input[name="classify-input-method"]:checked').value;
  const requiredFields = {
    'classify-markers': 'Markers TSV file',
    'classify-reference': 'Reference FASTA file',
    'classify-output': 'Output prefix location'
  };

  if (inputMethod === 'file') {
    requiredFields['classify-input'] = 'Genome FASTA to classify';
  } else {
    requiredFields['classify-input-list'] = 'Sample list file';
  }

  if (!validateForm(requiredFields, 'Classify')) return;

  const markers = document.getElementById('classify-markers').value;
  const reference = document.getElementById('classify-reference').value;
  const outputPrefix = document.getElementById('classify-output').value;
  const kmerSize = parseInt(document.getElementById('classify-kmer').value) || 21;
  const threads = parseInt(document.getElementById('classify-threads').value) || null;
  const nestedClassification = document.getElementById('classify-nested').checked;
  const gff = document.getElementById('classify-gff').value || null;
  const input = inputMethod === 'file' ? document.getElementById('classify-input').value : null;
  const inputList = inputMethod === 'list' ? document.getElementById('classify-input-list').value : null;

  startProgress('classify');
  logMessage('Starting genome classification...', 'info');

  try {
    const result = await tauriInvoke('run_classify', {
      params: {
        markers, reference, input, input_list: inputList, gff,
        output_prefix: outputPrefix, kmer_size: kmerSize, threads,
        nested_classification: nestedClassification
      }
    });
    await handleResultWithModal(result, outputPrefix + '_summary.tsv', 'classify');
  } catch (err) {
    logMessage('Error: ' + err, 'error');
  } finally {
    stopProgress();
  }
}

async function handleSplitFastqSubmit(e) {
  e.preventDefault();

  const inputMethod = document.querySelector('input[name="splitfq-input-method"]:checked').value;
  const requiredFields = {
    'splitfq-reference': 'Reference FASTA file',
    'splitfq-markers': 'Markers TSV file'
  };

  if (inputMethod === 'files') {
    requiredFields['splitfq-input'] = 'FASTQ input files';
  } else {
    requiredFields['splitfq-input-list'] = 'Sample list file';
  }

  if (!validateForm(requiredFields, 'Split-FASTQ')) return;

  const reference = document.getElementById('splitfq-reference').value;
  const markers = document.getElementById('splitfq-markers').value;
  const outputPrefix = document.getElementById('splitfq-output').value || 'split';
  const minDepth = parseInt(document.getElementById('splitfq-min-depth').value) || 10;
  const minAltPercent = parseInt(document.getElementById('splitfq-min-alt').value) || 95;
  const threads = parseInt(document.getElementById('splitfq-threads').value) || null;
  const paired = document.getElementById('splitfq-paired').checked;
  const nestedClassification = document.getElementById('splitfq-nested').checked;

  let input = null;
  let inputList = null;

  if (inputMethod === 'files') {
    const inputEl = document.getElementById('splitfq-input');
    input = inputEl.dataset.files ? JSON.parse(inputEl.dataset.files) : inputEl.value.split(';');
  } else {
    inputList = document.getElementById('splitfq-input-list').value;
  }

  startProgress('split-fastq');
  logMessage('Starting alignment-free genotyping...', 'info');

  try {
    const result = await tauriInvoke('run_split_fastq', {
      params: {
        input, input_list: inputList, paired, reference, markers, threads,
        output_prefix: outputPrefix, min_depth: minDepth, min_alt_percent: minAltPercent,
        nested_classification: nestedClassification
      }
    });
    await handleResultWithModal(result, outputPrefix + '_summary.tsv', 'splitfq');
  } catch (err) {
    logMessage('Error: ' + err, 'error');
  } finally {
    stopProgress();
  }
}

async function handleMatchSubmit(e) {
  e.preventDefault();

  const inputMethod = document.querySelector('input[name="match-input-method"]:checked').value;
  const refMethod = document.querySelector('input[name="match-ref-method"]:checked').value;

  const requiredFields = {
    'match-output': 'Output TSV file location'
  };

  // Require either references or index based on selected method
  if (refMethod === 'fasta') {
    requiredFields['match-references'] = 'References Multi-FASTA file';
  } else {
    requiredFields['match-index'] = 'Pre-built reference index file';
  }

  if (inputMethod === 'files') {
    requiredFields['match-input'] = 'FASTQ input files';
  } else {
    requiredFields['match-input-list'] = 'Sample list file';
  }

  if (!validateForm(requiredFields, 'Match')) return;

  const output = document.getElementById('match-output').value;
  const kmerSize = parseInt(document.getElementById('match-kmer').value) || 21;
  const threads = parseInt(document.getElementById('match-threads').value) || null;

  // Get references or index based on selected method
  let references = null;
  let index = null;
  if (refMethod === 'fasta') {
    references = document.getElementById('match-references').value;
  } else {
    index = document.getElementById('match-index').value;
  }

  let fastqs = null;
  let inputList = null;

  if (inputMethod === 'files') {
    const inputEl = document.getElementById('match-input');
    fastqs = inputEl.dataset.files ? JSON.parse(inputEl.dataset.files) : inputEl.value.split(';');
  } else {
    inputList = document.getElementById('match-input-list').value;
  }

  startProgress('match');
  const usingIndex = refMethod === 'index';
  logMessage(`Starting reference matching${usingIndex ? ' (using pre-built index)' : ''}...`, 'info');

  try {
    const result = await tauriInvoke('run_match', {
      params: { fastqs, input_list: inputList, references, index, output, kmer_size: kmerSize, threads }
    });
    await handleResultWithModal(result, output, 'match');
  } catch (err) {
    logMessage('Error: ' + err, 'error');
  } finally {
    stopProgress();
  }
}

// ============================================================================
// Reference Index Building
// ============================================================================

function setupBuildIndexButton() {
  const buildIndexBtn = document.getElementById('match-build-index-btn');
  const matchReferencesInput = document.getElementById('match-references');
  const matchIndexInput = document.getElementById('match-index');

  if (!buildIndexBtn || !matchReferencesInput) return;

  // Function to update button state
  const updateButtonState = () => {
    const hasReferences = matchReferencesInput.value && matchReferencesInput.value.trim() !== '';
    buildIndexBtn.disabled = !hasReferences;
  };

  // Listen for changes on the references input
  matchReferencesInput.addEventListener('change', updateButtonState);

  // Listen for changes on the index input to load index info
  if (matchIndexInput) {
    matchIndexInput.addEventListener('change', async (e) => {
      const indexPath = e.detail?.path || matchIndexInput.value;
      if (indexPath && indexPath.trim() !== '') {
        await loadIndexInfo(indexPath);
      } else {
        document.getElementById('match-index-info').classList.add('hidden');
      }
    });
  }

  // Initial state
  updateButtonState();
}

async function handleBuildIndex() {
  const buildIndexBtn = document.getElementById('match-build-index-btn');
  const referencesPath = document.getElementById('match-references').value;

  if (!referencesPath) {
    logMessage('Please select a references FASTA file first.', 'warning');
    return;
  }

  // Ask user where to save the index
  let outputPath;
  try {
    const { save } = window.__TAURI__.dialog;
    outputPath = await save({
      defaultPath: referencesPath.replace(/\.(fasta|fa|fna|fas)$/i, '.ptidx'),
      filters: [{ name: 'Pathotypr Index', extensions: ['ptidx'] }],
      title: 'Save Reference Index'
    });

    if (!outputPath) {
      logMessage('Index building cancelled.', 'info');
      return;
    }
  } catch (err) {
    logMessage('Error opening save dialog: ' + err, 'error');
    return;
  }

  // Get k-mer size from options
  const kmerSize = parseInt(document.getElementById('match-kmer').value) || 21;

  // Update button state
  buildIndexBtn.classList.add('building');
  buildIndexBtn.disabled = true;
  const originalText = buildIndexBtn.innerHTML;
  buildIndexBtn.innerHTML = `
    <svg class="spinner" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
      <circle cx="12" cy="12" r="10" stroke-dasharray="60" stroke-dashoffset="20"/>
    </svg>
    Building Index...
  `;

  logMessage(`Building reference index (k=${kmerSize})...`, 'info');
  openConsoleDrawer();

  try {
    const result = await tauriInvoke('build_reference_index', {
      params: {
        references: referencesPath,
        output: outputPath,
        kmer_size: kmerSize
      }
    });

    if (result.success) {
      logMessage(result.message, 'success');
      logMessage(`You can now use the "Pre-built Index" option with: ${makeClickablePath(outputPath)}`, 'info');
    } else {
      logMessage(result.message, 'error');
    }
  } catch (err) {
    logMessage('Error building index: ' + err, 'error');
  } finally {
    // Restore button state
    buildIndexBtn.classList.remove('building');
    buildIndexBtn.disabled = false;
    buildIndexBtn.innerHTML = originalText;
  }
}

async function loadIndexInfo(indexPath) {
  const indexInfoDiv = document.getElementById('match-index-info');
  const kmerSpan = document.getElementById('match-index-kmer');
  const refsSpan = document.getElementById('match-index-refs');

  if (!indexInfoDiv) return;

  try {
    const info = await tauriInvoke('get_index_info', { path: indexPath });

    if (info.valid) {
      kmerSpan.textContent = info.kmer_size;
      refsSpan.textContent = info.num_references;
      indexInfoDiv.classList.remove('hidden');
      logMessage(`Loaded index: ${info.num_references} references, k=${info.kmer_size}`, 'info');
    } else {
      indexInfoDiv.classList.add('hidden');
      logMessage(`Invalid index file: ${info.error}`, 'error');
    }
  } catch (err) {
    indexInfoDiv.classList.add('hidden');
    logMessage('Error reading index info: ' + err, 'error');
  }
}

// ============================================================================
// Utility Functions
// ============================================================================

function makeClickablePath(path) {
  const name = path.split(/[/\\]/).pop();
  return `<span class="clickable-path" data-path="${path}">${name}</span>`;
}

function logMessage(message, type = 'info') {
  const outputContent = document.getElementById('output-content');
  const entry = document.createElement('div');
  entry.className = `log-entry log-${type}`;

  const timestamp = new Date().toLocaleTimeString();
  const icon = type === 'success' ? '✓' : type === 'error' ? '✗' : type === 'warning' ? '⚠' : '→';

  // Process clickable paths in message
  entry.innerHTML = `[${timestamp}] ${icon} ${message}`;

  // Make paths clickable
  entry.querySelectorAll('.clickable-path').forEach(pathEl => {
    pathEl.addEventListener('click', () => {
      navigator.clipboard.writeText(pathEl.dataset.path);
      logMessage('Path copied to clipboard', 'info');
    });
  });

  const placeholder = outputContent.querySelector('.placeholder');
  if (placeholder) placeholder.remove();

  outputContent.appendChild(entry);
  outputContent.scrollTop = outputContent.scrollHeight;
  updateConsoleBadge();
}

function handleResult(result) {
  if (result.success) {
    logMessage(result.message, 'success');
    if (result.output_path) {
      logMessage(`Output: ${makeClickablePath(result.output_path)}`, 'info');
      showFileGeneratedAnimation(result.output_path);
    }
  } else {
    logMessage(result.message, 'error');
    openConsoleDrawer();
  }
}

async function handleResultWithModal(result, outputPath, panelId = null) {
  handleResult(result);

  if (result.success && outputPath) {
    // Determine panel ID if not provided
    if (!panelId) {
      const activePanel = document.querySelector('.panel.active');
      panelId = activePanel ? activePanel.id.replace('panel-', '') : null;
    }

    console.log('[Results] Attempting to show results for panel:', panelId, 'path:', outputPath);

    // Try to read and display results in inline viewer
    // Add a small delay to ensure file is written
    await new Promise(resolve => setTimeout(resolve, 200));

    try {
      const content = await readTextFile(outputPath);
      console.log('[Results] File content length:', content ? content.length : 'null');
      if (content && panelId) {
        showInlineResults(panelId, outputPath, content);
      } else if (!content) {
        console.warn('[Results] Could not read file or file is empty:', outputPath);
        logMessage('Results file could not be read. Check the output location.', 'warning');
      }
    } catch (err) {
      console.error('[Results] Error reading output file:', err);
      logMessage('Could not display results: ' + err.message, 'warning');
    }
  }
}

// ============================================================================
// Stepped Progress
// ============================================================================

function startProgress(command) {
  // Clear any existing interval first (safety measure)
  if (progressInterval) {
    clearInterval(progressInterval);
    progressInterval = null;
  }

  // Force reset state to ensure clean start
  isProcessing = true;
  console.log('[Progress] Starting progress for:', command);

  const progressBar = document.getElementById('progress-bar');
  const progressSteps = document.getElementById('progress-steps');
  const progressText = document.getElementById('progress-text');
  const runBtn = document.getElementById('btn-run');

  // Set up steps for this command
  const steps = processSteps[command] || ['Processing'];
  progressSteps.innerHTML = steps.map((label, i) => `
    <div class="step ${i === 0 ? 'active' : ''}" data-step="${i + 1}">
      <div class="step-dot"></div>
      <span class="step-label">${label}</span>
    </div>
  `).join('');

  progressText.textContent = steps[0] + '...';
  progressBar.classList.add('active');

  if (runBtn) {
    runBtn.classList.add('running');
    runBtn.dataset.originalClick = 'run';
    runBtn.innerHTML = `
      <svg class="btn-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
        <rect x="6" y="6" width="12" height="12" rx="2"/>
      </svg>
      <span>Cancel</span>
    `;
  }

  // Simulate progress through steps - store in global variable for cleanup
  let currentStep = 0;
  progressInterval = setInterval(() => {
    if (!isProcessing || currentStep >= steps.length - 1) {
      clearInterval(progressInterval);
      progressInterval = null;
      return;
    }

    currentStep++;
    updateProgressStep(currentStep, steps[currentStep]);
  }, 2000);
}

function updateProgressStep(stepIndex, label) {
  const steps = document.querySelectorAll('#progress-steps .step');
  const progressText = document.getElementById('progress-text');

  steps.forEach((step, i) => {
    step.classList.remove('active', 'completed');
    if (i < stepIndex) {
      step.classList.add('completed');
    } else if (i === stepIndex) {
      step.classList.add('active');
    }
  });

  progressText.textContent = label + '...';
}

function stopProgress() {
  console.log('[Progress] Stopping progress - isProcessing:', isProcessing);

  // Clear the interval immediately
  if (progressInterval) {
    clearInterval(progressInterval);
    progressInterval = null;
    console.log('[Progress] Cleared progress interval');
  }

  // Reset processing state
  isProcessing = false;

  const progressBar = document.getElementById('progress-bar');
  const runBtn = document.getElementById('btn-run');

  // Mark all steps as completed briefly
  document.querySelectorAll('#progress-steps .step').forEach(step => {
    step.classList.remove('active');
    step.classList.add('completed');
  });

  setTimeout(() => {
    progressBar.classList.remove('active');
  }, 500);

  if (runBtn) {
    runBtn.classList.remove('running');
    runBtn.innerHTML = `
      <svg class="btn-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
        <polygon points="5 3 19 12 5 21 5 3"/>
      </svg>
      <span>Run</span>
    `;
    console.log('[Progress] Run button reset - isProcessing now:', isProcessing);
  }
}

// ============================================================================
// File Generated Animation
// ============================================================================

function showFileGeneratedAnimation(outputPath) {
  document.querySelectorAll('.output-file').forEach(field => {
    if (field.value && outputPath.includes(field.value.split(/[/\\]/).pop())) {
      field.classList.add('file-generated');
      setTimeout(() => field.classList.remove('file-generated'), 2000);
    }
  });

  showSuccessToast(outputPath);
  showConfetti();
}

function showSuccessToast(outputPath) {
  const fileName = outputPath.split(/[/\\]/).pop();
  const toast = document.createElement('div');
  toast.className = 'success-toast';
  toast.innerHTML = `
    <span class="toast-icon">✓</span>
    <div class="toast-content">
      <span class="toast-title">File generated successfully!</span>
      <span class="toast-path">${fileName}</span>
    </div>
  `;
  document.body.appendChild(toast);
  setTimeout(() => toast.remove(), 4000);
}

function showConfetti() {
  const container = document.createElement('div');
  container.className = 'confetti-container';

  for (let i = 0; i < 50; i++) {
    const confetti = document.createElement('div');
    confetti.className = 'confetti';
    confetti.style.left = `${Math.random() * 100}%`;
    confetti.style.animationDelay = `${Math.random() * 0.5}s`;
    confetti.style.animationDuration = `${2 + Math.random() * 2}s`;

    const shapes = ['circle', 'square', 'rectangle'];
    const shape = shapes[Math.floor(Math.random() * shapes.length)];
    if (shape === 'circle') confetti.style.borderRadius = '50%';
    else if (shape === 'rectangle') {
      confetti.style.width = '6px';
      confetti.style.height = '14px';
    }

    container.appendChild(confetti);
  }

  document.body.appendChild(container);
  setTimeout(() => container.remove(), 4000);
}

function clearOutput() {
  document.getElementById('output-content').innerHTML = '<p class="placeholder">Output will appear here...</p>';
  resetConsoleBadge();
}

async function getVersion() {
  try {
    const version = await tauriInvoke('get_version');
    document.getElementById('app-version').textContent = `v${version}`;
  } catch (err) {
    console.error('Failed to get version:', err);
  }
}
