// ============================================================================
// Forms Module - Form Handling and Submissions
// ============================================================================

import { defaultParameters, demoGenotypingData, lineageMarkerData, drMarkerData, mtbcModel, TIMING } from './config.js';
import {
  clearRunCancellationRequest,
  clearRunContext,
  finishRunContext,
  getCurrentRun,
  isProcessing,
  isRunCurrent,
  requestRunCancellation,
  setProcessing,
  startRunContext,
  wasCancellationRequested
} from './state.js';
import { tauriInvoke, saveFileDialog, readTextFile, confirmDialog } from './tauri.js';
import {
  clearFormFieldErrors,
  formatError,
  getFileName,
  makeClickablePath,
  showFieldError,
  validateForm
} from './utils.js';
import { setDropzoneFile, setDropzoneFiles } from './dropzone.js';
import { logMessage, openConsoleDrawer, openConsoleDrawerIfClosed } from './console.js';
import {
  setCancellingState,
  setRunningState,
  showFileGeneratedAnimation,
  startProgress,
  stopProgress
} from './progress.js';
import { showInlineResults, resetInlineResults } from './results.js';
import { refreshVisualizationForTool, resetVisualizationForTool, renderTrainSummaryCard } from './visualization.js';

/**
 * Setup all form handlers
 */
export function setupForms() {
  document.getElementById('form-train')?.addEventListener('submit', handleTrainSubmit);
  document.getElementById('form-predict')?.addEventListener('submit', handlePredictSubmit);
  document.getElementById('form-classify')?.addEventListener('submit', handleClassifySubmit);
  document.getElementById('form-split-fastq')?.addEventListener('submit', handleSplitFastqSubmit);
  document.getElementById('form-match')?.addEventListener('submit', handleMatchSubmit);
}

/**
 * Setup save buttons
 */
export function setupSaveButtons() {
  document.querySelectorAll('.btn-save').forEach(btn => {
    btn.addEventListener('click', async () => {
      const targetId = btn.dataset.target;
      const defaultName = btn.dataset.default || 'output';

      try {
        const selected = await saveFileDialog({ defaultPath: defaultName });
        if (selected) {
          const input = document.getElementById(targetId);
          if (input) input.value = selected;
          logMessage(`Output: ${getFileName(selected)}`, 'success');
        }
      } catch (err) {
        logMessage('Error selecting save location: ' + formatError(err), 'error');
      }
    });
  });
}

/**
 * Setup a radio group that toggles between two input groups (files vs list)
 */
function setupInputMethodRadio(prefix) {
  document.querySelectorAll(`input[name="${prefix}-input-method"]`).forEach(radio => {
    radio.addEventListener('change', () => {
      const isFiles = radio.value === 'files';
      document.getElementById(`${prefix}-input-files-group`)?.classList.toggle('hidden', !isFiles);
      document.getElementById(`${prefix}-input-list-group`)?.classList.toggle('hidden', isFiles);
    });
  });
}

/**
 * Setup radio button handlers
 */
export function setupRadioHandlers() {
  setupInputMethodRadio('classify');
  setupInputMethodRadio('splitfq');
  setupInputMethodRadio('match');

  // Match reference method (FASTA vs Index)
  document.querySelectorAll('input[name="match-ref-method"]').forEach(radio => {
    radio.addEventListener('change', () => {
      const isFasta = radio.value === 'fasta';
      document.getElementById('match-fasta-group')?.classList.toggle('hidden', !isFasta);
      document.getElementById('match-index-group')?.classList.toggle('hidden', isFasta);

      // Lock/unlock k-mer size based on reference method
      const kmerInput = document.getElementById('match-kmer');
      const kmerGroup = kmerInput?.closest('.form-group');
      if (isFasta) {
        // Switching to FASTA - unlock k-mer size
        if (kmerInput) {
          kmerInput.disabled = false;
          kmerInput.title = '';
        }
        kmerGroup?.classList.remove('kmer-locked');
        document.getElementById('match-index-info')?.classList.add('hidden');
      }
      // Note: k-mer will be locked when index is loaded via loadIndexInfo
    });
  });

  // Build index button
  const buildIndexBtn = document.getElementById('match-build-index-btn');
  if (buildIndexBtn) {
    buildIndexBtn.addEventListener('click', handleBuildIndex);
  }

  setupBuildIndexButton();
}

/**
 * Setup default buttons
 */
export function setupDefaultButtons() {
  document.querySelectorAll('.btn-defaults').forEach(btn => {
    btn.addEventListener('click', () => {
      const tool = btn.dataset.tool;
      setDefaultParameters(tool);
    });
  });
}

/**
 * Set default parameters for a tool
 */
function setDefaultParameters(tool) {
  const toolDefaults = defaultParameters[tool];
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

/**
 * Setup demo button
 */
export function setupDemoButton() {
  const btn = document.getElementById('btn-load-demo');
  if (btn) btn.addEventListener('click', loadDemoData);

  // Quick-load buttons in Genotyping panels
  const lineageBtn = document.getElementById('btn-load-lineage-markers');
  if (lineageBtn) lineageBtn.addEventListener('click', () => loadQuickMarkers(lineageMarkerData, lineageBtn, 'lineage'));

  const drBtn = document.getElementById('btn-load-dr-markers');
  if (drBtn) drBtn.addEventListener('click', () => loadQuickMarkers(drMarkerData, drBtn, 'DR'));

  // Download MTBC model button in Predict panel
  const modelBtn = document.getElementById('btn-download-mtbc-model');
  if (modelBtn) modelBtn.addEventListener('click', handleDownloadMtbcModel);
}

/**
 * Download pre-trained MTBC model and load into Predict model dropzone
 */
async function handleDownloadMtbcModel(e) {
  e.preventDefault();

  const btn = document.getElementById('btn-download-mtbc-model');
  const originalContent = btn ? btn.innerHTML : '';

  if (btn) {
    btn.disabled = true;
    btn.innerHTML = `
      <svg class="spinner" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
        <circle cx="12" cy="12" r="10" stroke-dasharray="60" stroke-dashoffset="20"/>
      </svg>
      Downloading...
    `;
  }

  logMessage('Downloading MTBC pre-trained model...', 'info');
  openConsoleDrawerIfClosed(true);

  try {
    const result = await tauriInvoke('download_file', {
      params: { url: mtbcModel.url, filename: mtbcModel.filename }
    });

    if (result.success && result.path) {
      const modelDropzone = document.querySelector('[data-target="predict-model"]');
      if (modelDropzone) {
        await setDropzoneFile(modelDropzone, 'predict-model', result.path);
      }
      logMessage(`MTBC model loaded: ${mtbcModel.filename}`, 'success');
    } else {
      logMessage(`Failed to download MTBC model: ${result.error || 'Unknown error'}`, 'error');
    }
  } catch (err) {
    logMessage(`Error downloading MTBC model: ${formatError(err)}`, 'error');
  } finally {
    if (btn) {
      btn.innerHTML = originalContent;
      btn.disabled = false;
    }
  }
}

/**
 * Download and load marker + reference files into Classify dropzones
 */
async function loadQuickMarkers(dataConfig, btn, label) {
  const { markers, reference } = dataConfig;

  if (markers.url.startsWith('YOUR_')) {
    logMessage(`${label} markers URL not configured yet. Download manually and drag the file.`, 'warning');
    openConsoleDrawerIfClosed(true);
    return;
  }

  const originalHTML = btn.innerHTML;
  btn.disabled = true;
  btn.innerHTML = `
    <svg class="spinner" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
      <circle cx="12" cy="12" r="10" stroke-dasharray="60" stroke-dashoffset="20"/>
    </svg>
    Downloading...
  `;

  logMessage(`Downloading ${label} markers and reference...`, 'info');
  openConsoleDrawerIfClosed(true);

  try {
    const [markersResult, referenceResult] = await Promise.all([
      tauriInvoke('download_file', { params: { url: markers.url, filename: markers.filename } }),
      tauriInvoke('download_file', { params: { url: reference.url, filename: reference.filename } })
    ]);

    // Add markers to classify dropzone (accumulate with existing)
    if (markersResult.success && markersResult.path) {
      const classifyMarkersDropzone = document.querySelector('[data-target="classify-markers"]');
      if (classifyMarkersDropzone) {
        await setDropzoneFiles(classifyMarkersDropzone, 'classify-markers', [markersResult.path]);
      }
      // Also add to split-fastq (accumulate with existing)
      const splitfqMarkersDropzone = document.querySelector('[data-target="splitfq-markers"]');
      if (splitfqMarkersDropzone) {
        await setDropzoneFiles(splitfqMarkersDropzone, 'splitfq-markers', [markersResult.path]);
      }
      logMessage(`${label} markers loaded: ${markers.filename}`, 'success');
    } else {
      logMessage(`Failed to download ${label} markers: ${markersResult.error || 'Unknown error'}`, 'error');
    }

    // Set reference (only if not already set)
    if (referenceResult.success && referenceResult.path) {
      const classifyRefDropzone = document.querySelector('[data-target="classify-reference"]');
      const refInput = document.getElementById('classify-reference');
      if (classifyRefDropzone && !refInput?.value) {
        await setDropzoneFile(classifyRefDropzone, 'classify-reference', referenceResult.path);
      }
      const splitfqRefDropzone = document.querySelector('[data-target="splitfq-reference"]');
      const splitfqRefInput = document.getElementById('splitfq-reference');
      if (splitfqRefDropzone && !splitfqRefInput?.value) {
        await setDropzoneFile(splitfqRefDropzone, 'splitfq-reference', referenceResult.path);
      }
      logMessage(`Reference loaded: ${reference.filename}`, 'success');
    } else {
      logMessage(`Failed to download reference: ${referenceResult.error || 'Unknown error'}`, 'error');
    }

    saveFormState();
  } catch (err) {
    logMessage(`Error downloading ${label} data: ${err}`, 'error');
  } finally {
    btn.disabled = false;
    btn.innerHTML = originalHTML;
  }
}

/**
 * Load demo data
 */
async function loadDemoData() {
  document.getElementById('train-kmer').value = 21;
  document.getElementById('train-split').value = 0.2;
  document.getElementById('train-threads').value = 4;

  document.getElementById('classify-kmer').value = 31;
  document.getElementById('classify-threads').value = 4;
  document.getElementById('classify-nested').checked = true;
  document.getElementById('classify-min-flank').value = 10;
  document.getElementById('classify-masked-fasta').checked = false;

  document.getElementById('splitfq-kmer').value = 21;
  document.getElementById('splitfq-min-depth').value = 10;
  document.getElementById('splitfq-min-alt').value = 95;
  document.getElementById('splitfq-threads').value = 4;

  document.getElementById('match-kmer').value = 31;
  document.getElementById('match-threads').value = 4;
  document.getElementById('match-max-ref-freq').value = 0;
  // These two controls are optional in the DOM; guard the assignments so a
  // missing element does not throw a TypeError and abort the demo loader.
  const matchCoarseTopN = document.getElementById('match-coarse-top-n');
  if (matchCoarseTopN) matchCoarseTopN.value = 128;
  const matchCoarseStride = document.getElementById('match-coarse-stride');
  if (matchCoarseStride) matchCoarseStride.value = 8;
  document.getElementById('match-early-stop-confidence').value = 0;
  document.getElementById('match-early-stop-min-kmers').value = 1000000;
  document.getElementById('match-strict-percentages').checked = true;
  document.getElementById('match-auto-tune').checked = true;

  saveFormState();

  // Download MTB genotyping data (markers + reference) AND pre-trained model
  await Promise.all([
    loadGenotypingDemoData(),
    downloadMtbcModelSilent()
  ]);
}

/**
 * Download MTBC model silently (used by Load MTB Data button)
 */
async function downloadMtbcModelSilent() {
  if (mtbcModel.url.startsWith('YOUR_')) return;

  try {
    const result = await tauriInvoke('download_file', {
      params: { url: mtbcModel.url, filename: mtbcModel.filename }
    });

    if (result.success && result.path) {
      const modelDropzone = document.querySelector('[data-target="predict-model"]');
      if (modelDropzone) {
        await setDropzoneFile(modelDropzone, 'predict-model', result.path);
      }
      logMessage(`MTBC model loaded: ${mtbcModel.filename}`, 'success');
    } else {
      logMessage(`Failed to download MTBC model: ${result.error || 'Unknown error'}`, 'error');
    }
  } catch (err) {
    logMessage(`Error downloading MTBC model: ${formatError(err)}`, 'error');
  }
}

/**
 * Download and load MTB genotyping data (markers TSV + reference FASTA)
 * for both Classify and Split-FASTQ panels
 */
async function loadGenotypingDemoData() {
  const { markers, reference } = demoGenotypingData;

  // Validate that URLs have been configured
  if (markers.url.startsWith('YOUR_') || reference.url.startsWith('YOUR_')) {
    logMessage('MTB data URLs not configured yet.', 'warning');
    return;
  }

  const btn = document.getElementById('btn-load-demo');
  const originalContent = btn ? btn.innerHTML : '';
  if (btn) {
    btn.disabled = true;
    btn.innerHTML = `
      <svg class="spinner" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
        <circle cx="12" cy="12" r="10" stroke-dasharray="60" stroke-dashoffset="20"/>
      </svg>
      Downloading MTB data...
    `;
  }

  logMessage('Downloading M. tuberculosis markers and reference genome...', 'info');
  openConsoleDrawerIfClosed(true);

  try {
    // Download markers and reference in parallel
    const [markersResult, referenceResult] = await Promise.all([
      tauriInvoke('download_file', { params: { url: markers.url, filename: markers.filename } }),
      tauriInvoke('download_file', { params: { url: reference.url, filename: reference.filename } })
    ]);

    let success = true;

    // Set markers on Classify and Split-FASTQ dropzones
    if (markersResult.success && markersResult.path) {
      const classifyMarkersDropzone = document.querySelector('[data-target="classify-markers"]');
      if (classifyMarkersDropzone) {
        await setDropzoneFile(classifyMarkersDropzone, 'classify-markers', markersResult.path);
      }
      const splitfqMarkersDropzone = document.querySelector('[data-target="splitfq-markers"]');
      if (splitfqMarkersDropzone) {
        await setDropzoneFile(splitfqMarkersDropzone, 'splitfq-markers', markersResult.path);
      }
      logMessage(`Markers loaded: ${markers.filename}`, 'success');
    } else {
      logMessage(`Failed to download markers: ${markersResult.error || 'Unknown error'}`, 'error');
      success = false;
    }

    // Set reference on Classify and Split-FASTQ dropzones
    if (referenceResult.success && referenceResult.path) {
      const classifyRefDropzone = document.querySelector('[data-target="classify-reference"]');
      if (classifyRefDropzone) {
        await setDropzoneFile(classifyRefDropzone, 'classify-reference', referenceResult.path);
      }
      const splitfqRefDropzone = document.querySelector('[data-target="splitfq-reference"]');
      if (splitfqRefDropzone) {
        await setDropzoneFile(splitfqRefDropzone, 'splitfq-reference', referenceResult.path);
      }
      logMessage(`Reference loaded: ${reference.filename}`, 'success');
    } else {
      logMessage(`Failed to download reference: ${referenceResult.error || 'Unknown error'}`, 'error');
      success = false;
    }

    if (success) {
      logMessage('MTB data ready! Markers and reference genome loaded in Classify and Split-FASTQ. Select your input genomes/reads to run.', 'success');
    }
  } catch (err) {
    logMessage('Error downloading MTB data: ' + formatError(err), 'error');
  } finally {
    if (btn) {
      btn.disabled = false;
      btn.innerHTML = originalContent;
    }
  }
}

/**
 * Save form state to localStorage
 */
export function saveFormState() {
  const state = {
    'train-kmer': document.getElementById('train-kmer')?.value,
    'train-split': document.getElementById('train-split')?.value,
    'train-threads': document.getElementById('train-threads')?.value,
    'classify-kmer': document.getElementById('classify-kmer')?.value,
    'classify-threads': document.getElementById('classify-threads')?.value,
    'classify-nested': document.getElementById('classify-nested')?.checked,
    'classify-min-flank': document.getElementById('classify-min-flank')?.value,
    'classify-masked-fasta': document.getElementById('classify-masked-fasta')?.checked,
    'splitfq-kmer': document.getElementById('splitfq-kmer')?.value,
    'splitfq-min-depth': document.getElementById('splitfq-min-depth')?.value,
    'splitfq-min-alt': document.getElementById('splitfq-min-alt')?.value,
    'splitfq-threads': document.getElementById('splitfq-threads')?.value,
    'splitfq-paired': document.getElementById('splitfq-paired')?.checked,
    'splitfq-nested': document.getElementById('splitfq-nested')?.checked,
    'match-kmer': document.getElementById('match-kmer')?.value,
    'match-threads': document.getElementById('match-threads')?.value,
    'match-max-ref-freq': document.getElementById('match-max-ref-freq')?.value,
    'match-coarse-top-n': document.getElementById('match-coarse-top-n')?.value,
    'match-coarse-stride': document.getElementById('match-coarse-stride')?.value,
    'match-early-stop-confidence': document.getElementById('match-early-stop-confidence')?.value,
    'match-early-stop-min-kmers': document.getElementById('match-early-stop-min-kmers')?.value,
    'match-strict-percentages': document.getElementById('match-strict-percentages')?.checked,
    'match-auto-tune': document.getElementById('match-auto-tune')?.checked
  };
  localStorage.setItem('pathotypr-form-state', JSON.stringify(state));
}

/**
 * Load saved form state from localStorage
 */
export function loadSavedFormState() {
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

/**
 * Setup form state persistence (debounced to avoid localStorage thrashing)
 */
let saveFormTimeout = null;
export function setupFormStatePersistence() {
  document.addEventListener('change', (e) => {
    if (e.target.closest('form')) {
      if (saveFormTimeout) clearTimeout(saveFormTimeout);
      saveFormTimeout = setTimeout(() => {
        saveFormState();
        saveFormTimeout = null;
      }, TIMING.FORM_SAVE_DEBOUNCE);
    }
  });
}

function isCancellationMessage(text) {
  return /cancel(?:led|ed|ado|ada|ación|acion)/i.test(String(text || ''));
}

function isAbsolutePath(path) {
  return /^([A-Za-z]:[\\/]|\/|\\\\[^\\]+\\[^\\]+)/.test(String(path || '').trim());
}

function warnIfRelativePath(path, label) {
  if (!path || isAbsolutePath(path)) return;
  logMessage(`${label} is relative. The run works, but inline preview may be unavailable.`, 'warning');
}

function parseNumberField(fieldId, options) {
  const {
    label,
    fallback = null,
    min = null,
    max = null,
    integer = false
  } = options;

  const element = document.getElementById(fieldId);
  const raw = String(element?.value ?? '').trim().replace(',', '.');
  if (raw === '') {
    return { ok: true, value: fallback };
  }

  const parsed = Number(raw);
  if (!Number.isFinite(parsed)) {
    showFieldError(fieldId, `${label} must be a valid number.`);
    logMessage(`${label} must be a valid number.`, 'error');
    return { ok: false, value: fallback };
  }

  if (integer && !Number.isInteger(parsed)) {
    showFieldError(fieldId, `${label} must be an integer.`);
    logMessage(`${label} must be an integer.`, 'error');
    return { ok: false, value: fallback };
  }

  if (min !== null && parsed < min) {
    showFieldError(fieldId, `${label} must be at least ${min}.`);
    logMessage(`${label} must be at least ${min}.`, 'error');
    return { ok: false, value: fallback };
  }

  if (max !== null && parsed > max) {
    showFieldError(fieldId, `${label} must be at most ${max}.`);
    logMessage(`${label} must be at most ${max}.`, 'error');
    return { ok: false, value: fallback };
  }

  return { ok: true, value: parsed };
}

function parseInputFilesField(fieldId, label) {
  const inputEl = document.getElementById(fieldId);
  if (!inputEl) return [];

  if (inputEl.dataset.files) {
    try {
      const parsed = JSON.parse(inputEl.dataset.files);
      if (!Array.isArray(parsed)) throw new Error('Invalid file list');
      return parsed.map(path => String(path || '').trim()).filter(Boolean);
    } catch {
      const message = `${label} could not be parsed. Re-select the files.`;
      showFieldError(fieldId, message);
      logMessage(message, 'error');
      return null;
    }
  }

  return String(inputEl.value || '')
    .split(';')
    .map(path => path.trim())
    .filter(Boolean);
}

function getParentDir(path) {
  const value = String(path || '').trim();
  if (!value) return '';
  const slash = Math.max(value.lastIndexOf('/'), value.lastIndexOf('\\'));
  return slash >= 0 ? value.slice(0, slash) : '';
}

function joinPath(baseDir, childPath) {
  const base = String(baseDir || '').replace(/[\\/]+$/, '');
  const child = String(childPath || '').replace(/^[\\/]+/, '');
  if (!base) return child;
  const useBackslash = base.includes('\\') && !base.includes('/');
  const sep = useBackslash ? '\\' : '/';
  return `${base}${sep}${child}`;
}

function resolveRelativePath(baseDir, candidatePath) {
  let value = String(candidatePath || '').trim();
  if (!value) return '';
  if (!isAbsolutePath(value) && baseDir) {
    value = joinPath(baseDir, value);
  }
  return value;
}

function normalizeListHeaderToken(value) {
  return String(value || '')
    .trim()
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '_');
}

function detectInputListHasHeader(columns) {
  if (!Array.isArray(columns) || columns.length < 2) return false;
  // A second column that looks like a filesystem path (contains a separator or
  // ends in a file extension) means this is a data row, not a header. Without
  // this, a legitimate first sample such as "genome_A\t/data/genome_A.fasta"
  // is misread as a header (col0 matches 'genome', col1 matches 'fasta') and
  // silently dropped from the classify track context.
  const secondCol = String(columns[1] || '').trim();
  if (/[\\/]/.test(secondCol) || /\.[A-Za-z0-9]+$/.test(secondCol)) return false;
  const col0 = normalizeListHeaderToken(columns[0]);
  const col1 = normalizeListHeaderToken(columns[1]);
  const hasSample = col0.includes('sample') || col0.includes('genome') || col0.includes('name');
  const hasFasta = col1.includes('fasta') || col1.includes('path') || col1.includes('genome');
  return hasSample && hasFasta;
}

async function parseClassifyInputListMap(inputListPath) {
  const parsed = {
    sampleToFasta: {},
    sampleToGff: {}
  };
  if (!inputListPath || !isAbsolutePath(inputListPath)) return parsed;

  const content = await readTextFile(inputListPath);
  if (!content) return parsed;

  const lines = String(content)
    .split(/\r?\n/)
    .map(line => line.trim())
    .filter(Boolean);
  if (lines.length === 0) return parsed;

  const columns0 = lines[0].split('\t').map(col => col.trim());
  const hasHeader = detectInputListHasHeader(columns0);
  const startIdx = hasHeader ? 1 : 0;
  const inputListDir = getParentDir(inputListPath);

  for (let i = startIdx; i < lines.length; i += 1) {
    const cols = lines[i].split('\t');
    if (cols.length < 2) continue;
    const sample = String(cols[0] || '').trim();
    const fastaPath = resolveRelativePath(inputListDir, cols[1]);
    if (!sample || !fastaPath) continue;
    parsed.sampleToFasta[sample] = fastaPath;

    if (cols.length >= 3) {
      const gffPath = resolveRelativePath(inputListDir, cols[2]);
      if (gffPath) {
        parsed.sampleToGff[sample] = gffPath;
      }
    }
  }
  return parsed;
}

async function buildClassifyTrackContext(reference, inputMethod, inputFiles, inputList, gff, gffFiles) {
  const normalizedInputFiles = Array.isArray(inputFiles) ? [...inputFiles] : [];
  const normalizedGffFiles = Array.isArray(gffFiles)
    ? gffFiles.map(path => String(path || '').trim()).filter(Boolean)
    : [];
  const context = {
    trackContext: {
      referencePath: reference || '',
      inputMethod,
      fastaFiles: normalizedInputFiles,
      sampleToFasta: {},
      gffPath: String(gff || '').trim(),
      gffFiles: normalizedGffFiles,
      sampleToGff: {}
    }
  };

  if (inputMethod === 'list' && inputList) {
    const parsed = await parseClassifyInputListMap(inputList);
    context.trackContext.sampleToFasta = parsed.sampleToFasta;
    context.trackContext.sampleToGff = parsed.sampleToGff;
  } else if (inputMethod === 'files' && normalizedGffFiles.length === 1) {
    context.trackContext.gffPath = normalizedGffFiles[0];
  }

  return context;
}

async function executeRun(options) {
  const {
    command,
    panelId,
    startMessage,
    invokeCommand,
    params,
    onResult
  } = options;

  if (panelId) {
    resetInlineResults(panelId);
    resetVisualizationForTool(panelId);
  }

  const runId = startRunContext(command, panelId);
  startProgress(command, { runId, message: startMessage || 'Starting...' });
  if (startMessage) logMessage(startMessage, 'info');

  let outcome = 'error';
  try {
    const result = await tauriInvoke(invokeCommand, { params });
    if (!isRunCurrent(runId)) return;

    if (wasCancellationRequested(runId)) {
      outcome = 'cancelled';
      logMessage('Cancellation requested. Ignoring late backend result.', 'info');
      if (panelId) {
        resetInlineResults(panelId);
        resetVisualizationForTool(panelId);
      }
      return;
    }

    outcome = await onResult(result, runId);
  } catch (err) {
    const friendly = formatError(err);
    if (wasCancellationRequested(runId) || isCancellationMessage(friendly)) {
      outcome = 'cancelled';
      logMessage('Task cancelled.', 'info');
      if (panelId) {
        resetInlineResults(panelId);
        resetVisualizationForTool(panelId);
      }
    } else {
      outcome = 'error';
      logMessage(friendly, 'error');
      openConsoleDrawer();
    }
  } finally {
    stopProgress({ runId, outcome });
    finishRunContext(runId, outcome);
    clearRunContext(runId);
  }
}

/**
 * Setup run button
 */
export function setupRunButton() {
  const runBtn = document.getElementById('btn-run');
  if (runBtn) runBtn.addEventListener('click', handleRunClick);
}

/**
 * Handle run button click
 */
async function handleRunClick() {
  if (isProcessing()) {
    const run = getCurrentRun();
    if (!run?.id) return;
    if (run.status === 'cancelling') {
      logMessage('Cancellation already requested. Waiting for backend to stop...', 'info');
      return;
    }

    const cancelRequested = requestRunCancellation(run.id);
    if (!cancelRequested) return;
    setCancellingState('Cancelling task...');
    logMessage('Cancelling task...', 'warning');
    try {
      const cancelled = await tauriInvoke('cancel_task');
      if (cancelled) {
        logMessage('Cancel signal sent. Waiting for task to stop...', 'info');
      } else {
        clearRunCancellationRequest(run.id);
        setRunningState('Cancel request failed. Task is still running...');
        logMessage('Cancel request could not be delivered. Waiting for task completion...', 'warning');
      }
    } catch (err) {
      clearRunCancellationRequest(run.id);
      setRunningState('Cancel request failed. Task is still running...');
      logMessage('Error requesting cancellation: ' + formatError(err), 'error');
    }
    return;
  }

  // Guard against double-click: mark processing immediately
  setProcessing(true);

  const activePanel = document.querySelector('.panel.active');
  if (!activePanel) { setProcessing(false); return; }

  if (activePanel.id === 'panel-home') {
    logMessage('Select a tool from the sidebar or quick actions to run', 'warning');
    setProcessing(false);
    return;
  }

  const form = activePanel.querySelector('form');
  if (form) {
    form.dispatchEvent(new Event('submit', { cancelable: true, bubbles: true }));
  } else {
    setProcessing(false);
  }
}

// Form submit handlers
async function handleTrainSubmit(e) {
  e.preventDefault();
  clearFormFieldErrors(e.target);
  let startedRun = false;
  try {
    if (!validateForm({
      'train-input': 'Training FASTA file',
      'train-output': 'Output model file location'
    }, 'Train', logMessage, openConsoleDrawer)) return;

    const confirmed = await confirmDialog(
      'Model training can take a long time depending on the input size. Continue?',
      'Start Training'
    );
    if (!confirmed) return;

    const input = document.getElementById('train-input').value;
    const output = document.getElementById('train-output').value;
    const kmerSize = parseNumberField('train-kmer', {
      label: 'Train k-mer size',
      fallback: 21,
      min: 1,
      max: 31,
      integer: true
    });
    const testSplit = parseNumberField('train-split', {
      label: 'Train test split',
      fallback: 0.2,
      min: 0.01,
      max: 0.99
    });
    const threads = parseNumberField('train-threads', {
      label: 'Train threads',
      fallback: null,
      min: 1,
      integer: true
    });
    const cvFolds = parseNumberField('train-cv-folds', {
      label: 'CV folds',
      fallback: null,
      min: 2,
      max: 50,
      integer: true
    });
    const maxDepth = parseNumberField('train-max-depth', {
      label: 'Max depth',
      fallback: 20,
      min: 1,
      max: 50,
      integer: true
    });
    const minLeaf = parseNumberField('train-min-leaf', {
      label: 'Min leaf samples',
      fallback: 5,
      min: 1,
      max: 100,
      integer: true
    });
    if (!kmerSize.ok || !testSplit.ok || !threads.ok || !cvFolds.ok || !maxDepth.ok || !minLeaf.ok) return;

    warnIfRelativePath(output, 'Train output path');
    logMessage(`Input: ${makeClickablePath(input)}`, 'info');

    startedRun = true;
    await executeRun({
      command: 'train',
      panelId: 'train',
      startMessage: 'Starting model training...',
      invokeCommand: 'run_train',
      params: {
        input,
        output,
        kmer_size: kmerSize.value,
        test_split: testSplit.value,
        threads: threads.value,
        cv_folds: cvFolds.value,
        max_depth: maxDepth.value,
        min_samples_leaf: minLeaf.value
      },
      onResult: async (result, runId) => {
        const outcome = handleResult(result, runId);
        if (outcome === 'success' && result?.extra) {
          renderTrainSummaryCard('train', result.extra);
        }
        return outcome;
      }
    });
  } finally {
    if (!startedRun) setProcessing(false);
  }
}

async function handlePredictSubmit(e) {
  e.preventDefault();
  clearFormFieldErrors(e.target);
  let startedRun = false;
  try {
    if (!validateForm({
      'predict-input': 'FASTA file to classify',
      'predict-model': 'Trained model file (.zst)',
      'predict-output': 'Output results file location'
    }, 'Predict', logMessage, openConsoleDrawer)) return;

    const input = document.getElementById('predict-input').value;
    const model = document.getElementById('predict-model').value;
    const output = document.getElementById('predict-output').value;
    const threads = parseNumberField('predict-threads', {
      label: 'Predict threads',
      fallback: null,
      min: 1,
      integer: true
    });
    if (!threads.ok) return;

    warnIfRelativePath(output, 'Predict output path');

    startedRun = true;
    await executeRun({
      command: 'predict',
      panelId: 'predict',
      startMessage: 'Starting lineage prediction...',
      invokeCommand: 'run_predict',
      params: { input, model, output, threads: threads.value },
      onResult: async (result, runId) =>
        handleResultWithModal(result, output, 'predict', input, runId)
    });
  } finally {
    if (!startedRun) setProcessing(false);
  }
}

async function handleClassifySubmit(e) {
  e.preventDefault();
  clearFormFieldErrors(e.target);
  let startedRun = false;
  try {
    const inputMethod = document.querySelector('input[name="classify-input-method"]:checked').value;
    const requiredFields = {
      'classify-markers': 'Markers TSV file',
      'classify-reference': 'Reference FASTA file',
      'classify-output': 'Output prefix location'
    };

    if (inputMethod === 'files') {
      requiredFields['classify-input'] = 'FASTA files to classify';
    } else {
      requiredFields['classify-input-list'] = 'Sample list file';
    }

    if (!validateForm(requiredFields, 'Classify', logMessage, openConsoleDrawer)) return;

    const markersInput = document.getElementById('classify-markers');
    let markerFiles = null;
    let markers = null;
    if (markersInput.dataset.files) {
      try {
        markerFiles = JSON.parse(markersInput.dataset.files);
      } catch {
        markers = markersInput.value;
      }
    } else {
      markers = markersInput.value;
    }
    const reference = document.getElementById('classify-reference').value;
    const outputPrefix = document.getElementById('classify-output').value;
    const kmerSize = parseNumberField('classify-kmer', {
      label: 'Classify k-mer size',
      fallback: 31,
      min: 1,
      max: 31,
      integer: true
    });
    const threads = parseNumberField('classify-threads', {
      label: 'Classify threads',
      fallback: null,
      min: 1,
      integer: true
    });
    if (!kmerSize.ok || !threads.ok) return;

    const nestedClassification = document.getElementById('classify-nested').checked;
    const outputMaskedFasta = document.getElementById('classify-masked-fasta').checked;
    const minFlankBases = parseNumberField('classify-min-flank', {
      label: 'Min flank bases',
      fallback: 10,
      min: 0,
      max: 15,
      integer: true
    });
    if (!minFlankBases.ok) return;

    let inputFiles = null;
    let inputList = null;
    let gff = null;
    let gffFiles = null;

    if (inputMethod === 'files') {
      inputFiles = parseInputFilesField('classify-input', 'Classify input files');
      if (inputFiles === null) return;

      // Read GFF files as array for batch mode (matched by filename)
      const gffEl = document.getElementById('classify-gff');
      if (gffEl && gffEl.dataset.files) {
        try {
          const parsedGffFiles = JSON.parse(gffEl.dataset.files);
          gffFiles = Array.isArray(parsedGffFiles) ? parsedGffFiles : null;
        } catch {
          // Fallback: single GFF value
          if (gffEl.value) gffFiles = gffEl.value.split(';').filter(Boolean);
        }
      } else if (gffEl && gffEl.value) {
        gffFiles = [gffEl.value];
      }
      if (gffFiles && gffFiles.length === 0) gffFiles = null;
    } else {
      inputList = document.getElementById('classify-input-list').value;
      // For input-list mode, GFF is embedded in the TSV or use single GFF
      gff = document.getElementById('classify-gff').value || null;
    }

    warnIfRelativePath(outputPrefix, 'Classify output prefix');

    const classifyTrackContext = await buildClassifyTrackContext(
      reference,
      inputMethod,
      inputFiles,
      inputList,
      gff,
      gffFiles
    );

    startedRun = true;
    await executeRun({
      command: 'classify',
      panelId: 'classify',
      startMessage: 'Starting genome classification...',
      invokeCommand: 'run_classify',
      params: {
        markers: markerFiles ? markerFiles[0] : markers,
        marker_files: markerFiles,
        reference,
        input: null,
        input_files: inputFiles,
        input_list: inputList,
        gff,
        gff_files: gffFiles,
        output_prefix: outputPrefix,
        kmer_size: kmerSize.value,
        threads: threads.value,
        nested_classification: nestedClassification,
        output_masked_fasta: outputMaskedFasta,
        min_flank_bases: minFlankBases.value
      },
      onResult: async (result, runId) =>
        handleClassifyResult(result, runId, classifyTrackContext)
    });
  } finally {
    if (!startedRun) setProcessing(false);
  }
}

async function handleSplitFastqSubmit(e) {
  e.preventDefault();
  clearFormFieldErrors(e.target);
  let startedRun = false;
  try {
    const inputMethod = document.querySelector('input[name="splitfq-input-method"]:checked').value;
    const requiredFields = {
      'splitfq-reference': 'Reference FASTA file',
      'splitfq-markers': 'Markers TSV file',
      'splitfq-output': 'Output prefix location'
    };

    if (inputMethod === 'files') {
      requiredFields['splitfq-input'] = 'FASTQ input files';
    } else {
      requiredFields['splitfq-input-list'] = 'Sample list file';
    }

    if (!validateForm(requiredFields, 'Split-FASTQ', logMessage, openConsoleDrawer)) return;

    const reference = document.getElementById('splitfq-reference').value;
    const markersInput = document.getElementById('splitfq-markers');
    let markerFiles = null;
    let markers = null;
    if (markersInput.dataset.files) {
      try {
        markerFiles = JSON.parse(markersInput.dataset.files);
      } catch {
        markers = markersInput.value;
      }
    } else {
      markers = markersInput.value;
    }
    const outputPrefix = document.getElementById('splitfq-output').value;
    const minDepth = parseNumberField('splitfq-min-depth', {
      label: 'Split-FASTQ min depth',
      fallback: 10,
      min: 1,
      integer: true
    });
    const minAltPercent = parseNumberField('splitfq-min-alt', {
      label: 'Split-FASTQ min alt %',
      fallback: 95,
      min: 1,
      max: 100,
      integer: true
    });
    const kmerSize = parseNumberField('splitfq-kmer', {
      label: 'Split-FASTQ k-mer size',
      fallback: 21,
      min: 11,
      max: 31,
      integer: true
    });
    const threads = parseNumberField('splitfq-threads', {
      label: 'Split-FASTQ threads',
      fallback: null,
      min: 1,
      integer: true
    });
    if (!kmerSize.ok || !minDepth.ok || !minAltPercent.ok || !threads.ok) return;

    const paired = document.getElementById('splitfq-paired').checked;
    const nestedClassification = document.getElementById('splitfq-nested').checked;

    let input = null;
    let inputList = null;

    if (inputMethod === 'files') {
      input = parseInputFilesField('splitfq-input', 'Split-FASTQ input files');
      if (input === null) return;
    } else {
      inputList = document.getElementById('splitfq-input-list').value;
    }

    warnIfRelativePath(outputPrefix, 'Split-FASTQ output prefix');

    startedRun = true;
    await executeRun({
      command: 'split-fastq',
      panelId: 'splitfq',
      startMessage: 'Starting alignment-free genotyping...',
      invokeCommand: 'run_split_fastq',
      params: {
        input,
        input_list: inputList,
        paired,
        reference,
        markers: markerFiles ? markerFiles[0] : markers,
        marker_files: markerFiles,
        kmer_size: kmerSize.value,
        threads: threads.value,
        output_prefix: outputPrefix,
        min_depth: minDepth.value,
        min_alt_percent: minAltPercent.value,
        nested_classification: nestedClassification
      },
      onResult: async (result, runId) =>
        handleSplitFastqResult(result, runId)
    });
  } finally {
    if (!startedRun) setProcessing(false);
  }
}

async function handleMatchSubmit(e) {
  e.preventDefault();
  clearFormFieldErrors(e.target);
  let startedRun = false;
  try {
    const inputMethod = document.querySelector('input[name="match-input-method"]:checked').value;
    const refMethod = document.querySelector('input[name="match-ref-method"]:checked').value;

    const requiredFields = {
      'match-output': 'Output TSV file location'
    };

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

    if (!validateForm(requiredFields, 'Match', logMessage, openConsoleDrawer)) return;

    const output = document.getElementById('match-output').value;
    const kmerSize = parseNumberField('match-kmer', {
      label: 'Match k-mer size',
      fallback: 21,
      min: 11,
      max: 51,
      integer: true
    });
    const threads = parseNumberField('match-threads', {
      label: 'Match threads',
      fallback: null,
      min: 1,
      integer: true
    });
    const maxRefFreq = parseNumberField('match-max-ref-freq', {
      label: 'Max reference k-mer frequency',
      fallback: 0,
      min: 0,
      max: 100
    });
    const coarseTopN = parseNumberField('match-coarse-top-n', {
      label: 'Coarse candidate references (top N)',
      fallback: 128,
      min: 0,
      integer: true
    });
    const coarseStride = parseNumberField('match-coarse-stride', {
      label: 'Coarse scan stride',
      fallback: 8,
      min: 1,
      integer: true
    });
    const earlyStopConfidence = parseNumberField('match-early-stop-confidence', {
      label: 'Early-stop confidence',
      fallback: 0,
      min: 0,
      max: 1
    });
    const earlyStopMinKmers = parseNumberField('match-early-stop-min-kmers', {
      label: 'Early-stop minimum k-mers',
      fallback: 1000000,
      min: 1,
      integer: true
    });
    const autoTune = document.getElementById('match-auto-tune')?.checked ?? true;
    const strictPercentages = document.getElementById('match-strict-percentages')?.checked ?? true;
    if (!kmerSize.ok || !threads.ok || !maxRefFreq.ok || !coarseTopN.ok || !coarseStride.ok || !earlyStopConfidence.ok || !earlyStopMinKmers.ok) return;

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
      fastqs = parseInputFilesField('match-input', 'Match FASTQ input files');
      if (fastqs === null) return;
    } else {
      inputList = document.getElementById('match-input-list').value;
    }

    const usingIndex = refMethod === 'index';
    warnIfRelativePath(output, 'Match output path');

    startedRun = true;
    await executeRun({
      command: 'match',
      panelId: 'match',
      startMessage: `Starting reference matching${usingIndex ? ' (using pre-built index)' : ''}...`,
      invokeCommand: 'run_match',
      params: {
        fastqs,
        input_list: inputList,
        references,
        index,
        output,
        kmer_size: kmerSize.value,
        threads: threads.value,
        max_ref_freq: maxRefFreq.value,
        coarse_top_n: coarseTopN.value,
        coarse_stride: coarseStride.value,
        early_stop_confidence: earlyStopConfidence.value,
        early_stop_min_kmers: earlyStopMinKmers.value,
        strict_percentages: strictPercentages,
        auto_tune: autoTune
      },
      onResult: async (result, runId) =>
        handleResultWithModal(result, output, 'match', fastqs, runId)
    });
  } finally {
    if (!startedRun) setProcessing(false);
  }
}

// Result handlers
function handleResult(result, runId = 0) {
  if (runId && !isRunCurrent(runId)) return 'cancelled';

  const message = result?.message || '';
  if (isCancellationMessage(message)) {
    logMessage('Task cancelled.', 'info');
    return 'cancelled';
  }

  if (result.success) {
    if (message) logMessage(message, 'success');
    if (result.output_path) {
      logMessage(`Output: ${makeClickablePath(result.output_path)}`, 'info');
      showFileGeneratedAnimation(result.output_path);
    }
    return 'success';
  }

  logMessage(message || 'Task failed.', 'error');
  openConsoleDrawer();
  return 'error';
}

/**
 * Handle classify result with multiple marker set outputs.
 * Shows tabs when there are multiple output sets.
 */
async function handleClassifyResult(result, runId, context) {
  if (runId && !isRunCurrent(runId)) return 'cancelled';

  const outputSets = result?.extra?.output_sets;
  if (outputSets && outputSets.length > 1) {
    // Multiple marker sets — show with result set selector
    const baseOutcome = handleResult(result, runId);
    if (baseOutcome !== 'success') {
      resetInlineResults('classify');
      resetVisualizationForTool('classify');
      return baseOutcome;
    }

    // Wait for filesystem flush
    await new Promise(resolve => setTimeout(resolve, 200));

    // Set up tab bars for results and visualization
    const resultsArea = document.getElementById('classify-results');
    const vizPanel = document.getElementById('classify-visualization');

    let tabBar = resultsArea?.querySelector('.result-tabs');
    if (resultsArea && !tabBar) {
      tabBar = document.createElement('div');
      tabBar.className = 'result-tabs';
      resultsArea.prepend(tabBar);
    }
    if (tabBar) tabBar.innerHTML = '';

    let vizTabBar = vizPanel?.querySelector('.result-tabs');
    if (vizPanel && !vizTabBar) {
      vizTabBar = document.createElement('div');
      vizTabBar.className = 'result-tabs';
      vizPanel.prepend(vizTabBar);
    }
    if (vizTabBar) vizTabBar.innerHTML = '';

    // Sync function: update all tab bars + load content
    const switchToSet = async (index) => {
      [tabBar, vizTabBar].forEach(bar => {
        if (!bar) return;
        bar.querySelectorAll('.result-tab').forEach(t => {
          t.classList.toggle('active', Number(t.dataset.index) === index);
        });
      });
      try {
        const content = await readTextFile(outputSets[index].path);
        if (content) {
          showInlineResults('classify', outputSets[index].path, content, null, context);
          refreshVisualizationForTool('classify');
        }
      } catch (err) {
        logMessage('Could not load results: ' + err, 'warning');
      }
    };

    // Create tabs
    for (let i = 0; i < outputSets.length; i++) {
      if (tabBar) {
        const tab = document.createElement('button');
        tab.className = 'result-tab' + (i === 0 ? ' active' : '');
        tab.textContent = outputSets[i].name;
        tab.dataset.index = i;
        tab.addEventListener('click', () => switchToSet(i));
        tabBar.appendChild(tab);
      }

      if (vizTabBar) {
        const vizTab = document.createElement('button');
        vizTab.className = 'result-tab' + (i === 0 ? ' active' : '');
        vizTab.textContent = outputSets[i].name;
        vizTab.dataset.index = i;
        vizTab.addEventListener('click', () => switchToSet(i));
        vizTabBar.appendChild(vizTab);
      }
    }

    // Load first set
    await switchToSet(0);

    return 'success';
  }

  // Single marker set — use standard handler
  return handleResultWithModal(result, null, 'classify', null, runId, context);
}

/**
 * Handle split-fastq results with multi-marker tab support
 */
async function handleSplitFastqResult(result, runId) {
  if (runId && !isRunCurrent(runId)) return 'cancelled';

  const outputSets = result?.extra?.output_sets;
  if (outputSets && outputSets.length > 1) {
    const baseOutcome = handleResult(result, runId);
    if (baseOutcome !== 'success') {
      resetInlineResults('splitfq');
      resetVisualizationForTool('splitfq');
      return baseOutcome;
    }

    await new Promise(resolve => setTimeout(resolve, 200));

    const resultsArea = document.getElementById('splitfq-results');
    const vizPanel = document.getElementById('splitfq-visualization');

    let tabBar = resultsArea?.querySelector('.result-tabs');
    if (resultsArea && !tabBar) {
      tabBar = document.createElement('div');
      tabBar.className = 'result-tabs';
      resultsArea.prepend(tabBar);
    }
    if (tabBar) tabBar.innerHTML = '';

    let vizTabBar = vizPanel?.querySelector('.result-tabs');
    if (vizPanel && !vizTabBar) {
      vizTabBar = document.createElement('div');
      vizTabBar.className = 'result-tabs';
      vizPanel.prepend(vizTabBar);
    }
    if (vizTabBar) vizTabBar.innerHTML = '';

    const switchToSet = async (index) => {
      [tabBar, vizTabBar].forEach(bar => {
        if (!bar) return;
        bar.querySelectorAll('.result-tab').forEach(t => {
          t.classList.toggle('active', Number(t.dataset.index) === index);
        });
      });
      try {
        const content = await readTextFile(outputSets[index].path);
        if (content) {
          showInlineResults('splitfq', outputSets[index].path, content);
          refreshVisualizationForTool('splitfq');
        }
      } catch (err) {
        logMessage('Could not load results: ' + err, 'warning');
      }
    };

    for (let i = 0; i < outputSets.length; i++) {
      if (tabBar) {
        const tab = document.createElement('button');
        tab.className = 'result-tab' + (i === 0 ? ' active' : '');
        tab.textContent = outputSets[i].name;
        tab.dataset.index = i;
        tab.addEventListener('click', () => switchToSet(i));
        tabBar.appendChild(tab);
      }
      if (vizTabBar) {
        const vizTab = document.createElement('button');
        vizTab.className = 'result-tab' + (i === 0 ? ' active' : '');
        vizTab.textContent = outputSets[i].name;
        vizTab.dataset.index = i;
        vizTab.addEventListener('click', () => switchToSet(i));
        vizTabBar.appendChild(vizTab);
      }
    }

    await switchToSet(0);
    return 'success';
  }

  return handleResultWithModal(result, null, 'splitfq', null, runId);
}

async function handleResultWithModal(
  result,
  outputPath,
  panelId = null,
  inputFiles = null,
  runId = 0,
  context = null
) {
  if (runId && !isRunCurrent(runId)) return 'cancelled';
  const baseOutcome = handleResult(result, runId);

  if (!panelId) {
    const activePanel = document.querySelector('.panel.active');
    panelId = activePanel ? activePanel.id.replace('panel-', '') : null;
  }

  if (baseOutcome !== 'success') {
    if (panelId) {
      resetInlineResults(panelId);
      resetVisualizationForTool(panelId);
    }
    return baseOutcome;
  }

  if (runId && wasCancellationRequested(runId)) {
    if (panelId) {
      resetInlineResults(panelId);
      resetVisualizationForTool(panelId);
    }
    logMessage('Run completed after cancellation request. Results were ignored.', 'warning');
    return 'cancelled';
  }

  const resolvedOutputPath = result?.output_path || outputPath;

  if (resolvedOutputPath) {
    // Backend preview command only accepts absolute paths.
    if (!isAbsolutePath(resolvedOutputPath)) {
      logMessage('Results preview requires an absolute output path. Use "Save As" for output files.', 'warning');
      if (panelId) {
        resetInlineResults(panelId);
        resetVisualizationForTool(panelId);
      }
      return 'success';
    }

    // Wait for file system flush before reading results
    await new Promise(resolve => setTimeout(resolve, 200));

    try {
      if (runId && (!isRunCurrent(runId) || wasCancellationRequested(runId))) return 'cancelled';
      const content = await readTextFile(resolvedOutputPath);
      if (content && panelId) {
        showInlineResults(panelId, resolvedOutputPath, content, inputFiles, context);
        refreshVisualizationForTool(panelId);
      } else if (!content) {
        logMessage('Results file could not be read. Check the output location.', 'warning');
        if (panelId) {
          resetInlineResults(panelId);
          resetVisualizationForTool(panelId);
        }
      }
    } catch (err) {
      logMessage('Could not display results: ' + formatError(err), 'warning');
      if (panelId) {
        resetInlineResults(panelId);
        resetVisualizationForTool(panelId);
      }
    }
  } else if (panelId) {
    resetInlineResults(panelId);
    resetVisualizationForTool(panelId);
  }

  return 'success';
}

// Build index handlers
function setupBuildIndexButton() {
  const buildIndexBtn = document.getElementById('match-build-index-btn');
  const matchReferencesInput = document.getElementById('match-references');
  const matchIndexInput = document.getElementById('match-index');

  if (!buildIndexBtn || !matchReferencesInput) return;

  const updateButtonState = () => {
    const hasReferences = matchReferencesInput.value && matchReferencesInput.value.trim() !== '';
    buildIndexBtn.disabled = !hasReferences;
  };

  matchReferencesInput.addEventListener('change', updateButtonState);

  if (matchIndexInput) {
    matchIndexInput.addEventListener('change', async (e) => {
      const indexPath = e.detail?.path || matchIndexInput.value;
      if (indexPath && indexPath.trim() !== '') {
        await loadIndexInfo(indexPath);
      } else {
        document.getElementById('match-index-info')?.classList.add('hidden');
      }
    });
  }

  updateButtonState();
}

async function handleBuildIndex() {
  const buildIndexBtn = document.getElementById('match-build-index-btn');
  const referencesPath = document.getElementById('match-references').value;

  if (!referencesPath) {
    logMessage('Please select a references FASTA file first.', 'warning');
    return;
  }

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
    logMessage('Error opening save dialog: ' + formatError(err), 'error');
    return;
  }

  const kmerSize = parseInt(document.getElementById('match-kmer').value) || 21;

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
  openConsoleDrawerIfClosed(true);

  try {
    const result = await tauriInvoke('build_reference_index', {
      params: { references: referencesPath, output: outputPath, kmer_size: kmerSize }
    });

    if (result.success) {
      logMessage(result.message, 'success');
      logMessage(`You can now use the "Pre-built Index" option with: ${makeClickablePath(outputPath)}`, 'info');
    } else {
      logMessage(result.message, 'error');
    }
  } catch (err) {
    logMessage('Error building index: ' + formatError(err), 'error');
  } finally {
    buildIndexBtn.classList.remove('building');
    buildIndexBtn.disabled = false;
    buildIndexBtn.innerHTML = originalText;
  }
}

async function loadIndexInfo(indexPath) {
  const indexInfoDiv = document.getElementById('match-index-info');
  const kmerSpan = document.getElementById('match-index-kmer');
  const refsSpan = document.getElementById('match-index-refs');
  const kmerInput = document.getElementById('match-kmer');
  const kmerGroup = kmerInput?.closest('.form-group');

  if (!indexInfoDiv) return;

  try {
    const info = await tauriInvoke('get_index_info', { path: indexPath });

    if (info.valid) {
      kmerSpan.textContent = info.kmer_size;
      refsSpan.textContent = info.num_references;
      indexInfoDiv.classList.remove('hidden');

      // Lock k-mer size to index value
      if (kmerInput) {
        kmerInput.value = info.kmer_size;
        kmerInput.disabled = true;
        kmerInput.title = `K-mer size is locked to ${info.kmer_size} (from index)`;
      }
      kmerGroup?.classList.add('kmer-locked');

      logMessage(`Loaded index: ${info.num_references} references, k=${info.kmer_size} (k-mer size locked)`, 'info');
    } else {
      indexInfoDiv.classList.add('hidden');
      // Unlock k-mer on invalid index
      if (kmerInput) {
        kmerInput.disabled = false;
        kmerInput.title = '';
      }
      kmerGroup?.classList.remove('kmer-locked');
      logMessage(`Invalid index file: ${info.error}`, 'error');
    }
  } catch (err) {
    indexInfoDiv.classList.add('hidden');
    // Unlock k-mer on error
    if (kmerInput) {
      kmerInput.disabled = false;
      kmerInput.title = '';
    }
    kmerGroup?.classList.remove('kmer-locked');
    logMessage('Error reading index info: ' + formatError(err), 'error');
  }
}
