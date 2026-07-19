// ============================================================================
// Progress Module - Real-time Progress from Backend Events
// ============================================================================

import { processSteps, TIMING } from './config.js';
import {
  getCurrentRun,
  isRunCurrent,
  setProcessing
} from './state.js';
import { getSystemUsage } from './tauri.js';
import { getFileName } from './utils.js';
import { updateConsoleBadge } from './console.js';

// Store the event unlisten function
let progressUnlisten = null;
let completeUnlisten = null;
let logUnlisten = null;
let systemUsageInterval = null;
let systemUsageRequestInFlight = false;

function setRunButtonIdle() {
  const runBtn = document.getElementById('btn-run');
  if (!runBtn) return;
  runBtn.disabled = false;
  runBtn.classList.remove('running', 'cancelling');
  runBtn.innerHTML = `
    <svg class="btn-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
      <polygon points="5 3 19 12 5 21 5 3"/>
    </svg>
    <span>Run</span>
  `;
}

function setRunButtonRunning() {
  const runBtn = document.getElementById('btn-run');
  if (!runBtn) return;
  runBtn.disabled = false;
  runBtn.classList.add('running');
  runBtn.classList.remove('cancelling');
  runBtn.innerHTML = `
    <svg class="btn-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5">
      <circle cx="12" cy="12" r="10"/>
      <path d="M15 9l-6 6M9 9l6 6"/>
    </svg>
    <span>Cancel</span>
  `;
}

function setRunButtonCancelling() {
  const runBtn = document.getElementById('btn-run');
  if (!runBtn) return;
  runBtn.disabled = true;
  runBtn.classList.add('running', 'cancelling');
  runBtn.innerHTML = `
    <svg class="btn-icon spinner" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
      <circle cx="12" cy="12" r="10" stroke-dasharray="60" stroke-dashoffset="20"/>
    </svg>
    <span>Cancelling...</span>
  `;
}

function updateProgressMessage(message, stateClass = '') {
  const progressText = document.getElementById('progress-text');
  if (!progressText) return;
  progressText.classList.remove('state-running', 'state-success', 'state-error', 'state-cancelled', 'state-cancelling');
  if (stateClass) progressText.classList.add(stateClass);
  progressText.textContent = message;
}

function setProgressContainerState(stateClass = '') {
  const progressBar = document.getElementById('progress-bar');
  if (!progressBar) return;
  progressBar.classList.remove('state-running', 'state-success', 'state-error', 'state-cancelled', 'state-cancelling');
  if (stateClass) progressBar.classList.add(stateClass);
}

function formatPercent(value, maxDecimals = 1) {
  if (!Number.isFinite(value)) return '--%';
  const safe = Math.max(0, value);
  return `${safe.toFixed(maxDecimals)}%`;
}

function formatMemoryMB(value) {
  if (!Number.isFinite(value)) return '--';
  const safeMb = Math.max(0, value);
  if (safeMb >= 1024) return `${(safeMb / 1024).toFixed(2)} GB`;
  if (safeMb >= 100) return `${safeMb.toFixed(0)} MB`;
  return `${safeMb.toFixed(1)} MB`;
}

function updateSystemUsageDisplay(usage) {
  const panelEl = document.getElementById('system-usage-panel');
  const cpuEl = document.getElementById('system-cpu-value');
  const ramEl = document.getElementById('system-ram-value');
  if (!cpuEl || !ramEl) return;

  const cpuPercent = Number(usage?.cpu_percent);
  const memoryPercent = Number(usage?.memory_percent);
  const memoryMb = Number(usage?.memory_mb);
  const memoryBytes = Number(usage?.memory_rss_bytes);

  cpuEl.textContent = formatPercent(cpuPercent, 1);
  ramEl.textContent = formatMemoryMB(memoryMb);

  if (panelEl) {
    const cpuText = formatPercent(cpuPercent, 1);
    const ramText = formatMemoryMB(memoryMb);
    const memPctText = formatPercent(memoryPercent, 2);
    const bytesText = Number.isFinite(memoryBytes) ? `${Math.max(0, Math.round(memoryBytes)).toLocaleString()} B` : '--';
    panelEl.title = `Pathotypr usage — CPU ${cpuText} · RAM ${ramText} (${memPctText}, RSS ${bytesText})`;
  }
}

function startSystemUsageMonitoring() {
  if (systemUsageInterval) return;

  const pollUsage = async () => {
    if (systemUsageRequestInFlight) return;
    systemUsageRequestInFlight = true;
    try {
      const usage = await getSystemUsage();
      updateSystemUsageDisplay(usage);
    } finally {
      systemUsageRequestInFlight = false;
    }
  };

  pollUsage();
  systemUsageInterval = setInterval(pollUsage, TIMING.SYSTEM_USAGE_POLL);
}

function stopSystemUsageMonitoring() {
  if (systemUsageInterval) {
    clearInterval(systemUsageInterval);
    systemUsageInterval = null;
  }
  systemUsageRequestInFlight = false;
  updateSystemUsageDisplay(null);
}

/**
 * Setup progress event listeners (call once on app init)
 */
export async function setupProgressEvents() {
  try {
    const { listen } = window.__TAURI__.event;

    // Listen for progress events from backend
    progressUnlisten = await listen('progress', (event) => {
      const { step, total_steps, message } = event.payload;
      updateProgressFromEvent(step, total_steps, message);
      // Also log to console drawer
      logToConsole(message);
    });

    // Listen for completion event. Do NOT force every step to 'completed' here:
    // the backend emits 'progress-complete' on all terminal paths (including
    // errors), so marking steps green before the invoke promise resolves would
    // make a failed run appear fully successful. stopProgress() finalizes the
    // step visuals based on the real outcome once the command returns.
    completeUnlisten = await listen('progress-complete', () => {
      const run = getCurrentRun();
      if (!run?.id || run.status !== 'running') return;
    });

    // Listen for log events from backend
    logUnlisten = await listen('log', (event) => {
      const message = event.payload;
      console.log('[Backend Log]', message);
      logToConsole(message);
    });

    console.log('[Progress] Event listeners registered');
  } catch (err) {
    console.warn('[Progress] Failed to setup event listeners:', err);
  }

  startSystemUsageMonitoring();
}

/**
 * Log message to the console drawer
 */
function logToConsole(message) {
  const outputContent = document.getElementById('output-content');
  if (!outputContent) return;

  const placeholder = outputContent.querySelector('.placeholder');
  if (placeholder) placeholder.remove();

  const timestamp = new Date().toLocaleTimeString();
  const entry = document.createElement('div');
  entry.className = 'log-entry log-info';
  entry.textContent = `[${timestamp}] → ${message}`;
  outputContent.appendChild(entry);
  outputContent.scrollTop = outputContent.scrollHeight;
  updateConsoleBadge();
}

/**
 * Start progress display for a command
 */
export function startProgress(command, options = {}) {
  const { runId = 0, message = 'Starting...' } = options;

  const progressStepsEl = document.getElementById('progress-steps');
  const progressBar = document.getElementById('progress-bar');

  setProcessing(true);
  setProgressContainerState('state-running');

  // Set up steps for this command (initial state)
  const steps = processSteps[command] || ['Processing'];
  if (progressStepsEl) {
    progressStepsEl.innerHTML = steps.map((label, i) => `
      <div class="step" data-step="${i + 1}">
        <div class="step-dot"></div>
        <span class="step-label">${label}</span>
      </div>
    `).join('');
  }

  if (progressBar) {
    progressBar.classList.add('active');
    progressBar.dataset.runId = String(runId || '');
    progressBar.dataset.command = command || '';
  }
  updateProgressMessage(message, 'state-running');
  startSystemUsageMonitoring();
  setRunButtonRunning();
}

/**
 * Move UI into cancelling state immediately after user clicks cancel.
 */
export function setCancellingState(message = 'Cancelling task...') {
  const run = getCurrentRun();
  if (!run?.id) return;
  setProgressContainerState('state-cancelling');
  updateProgressMessage(message, 'state-cancelling');
  setRunButtonCancelling();

  document.querySelectorAll('#progress-steps .step.active').forEach(stepEl => {
    stepEl.classList.remove('active');
    stepEl.classList.add('cancelled');
  });
}

/**
 * Restore running state when cancellation could not be delivered.
 */
export function setRunningState(message = 'Task is still running...') {
  const run = getCurrentRun();
  if (!run?.id) return;
  setProgressContainerState('state-running');
  updateProgressMessage(message, 'state-running');
  setRunButtonRunning();

  document.querySelectorAll('#progress-steps .step.cancelled').forEach(stepEl => {
    stepEl.classList.remove('cancelled');
  });
}

/**
 * Update progress based on backend event
 */
function updateProgressFromEvent(stepIndex, totalSteps, message) {
  const run = getCurrentRun();
  if (!run?.id) return;

  const progressStepsEl = document.getElementById('progress-steps');
  const steps = progressStepsEl?.querySelectorAll('.step');

  if (steps && run.status === 'running') {
    const maxIndex = Math.max(0, steps.length - 1);
    const normalizedStep = Math.min(maxIndex, Math.max(0, Number(stepIndex) || 0));
    steps.forEach((stepEl, i) => {
      stepEl.classList.remove('active', 'completed', 'error', 'cancelled');
      if (i < normalizedStep) {
        stepEl.classList.add('completed');
      } else if (i === normalizedStep) {
        stepEl.classList.add('active');
      }
    });
  }

  if (run.status === 'cancelling') {
    updateProgressMessage('Cancelling task... waiting for backend to stop', 'state-cancelling');
  } else if (message) {
    updateProgressMessage(message, 'state-running');
  }
}

/**
 * Stop progress animation
 */
export function stopProgress(options = {}) {
  const { outcome = 'success', runId = 0, message = '' } = options;
  if (runId && !isRunCurrent(runId)) return;

  const progressBar = document.getElementById('progress-bar');
  const steps = document.querySelectorAll('#progress-steps .step');
  setProgressContainerState(
    outcome === 'success'
      ? 'state-success'
      : outcome === 'cancelled'
        ? 'state-cancelled'
        : 'state-error'
  );

  if (outcome === 'success') {
    steps.forEach(stepEl => {
      stepEl.classList.remove('active', 'error', 'cancelled');
      stepEl.classList.add('completed');
    });
  } else if (outcome === 'cancelled') {
    steps.forEach(stepEl => {
      stepEl.classList.remove('active', 'error');
      if (!stepEl.classList.contains('completed')) stepEl.classList.add('cancelled');
    });
  } else {
    const activeStep = document.querySelector('#progress-steps .step.active');
    if (activeStep) {
      activeStep.classList.remove('active');
      activeStep.classList.add('error');
    }
  }

  const finalMessage = message || (
    outcome === 'success'
      ? 'Completed successfully'
      : outcome === 'cancelled'
        ? 'Task cancelled'
        : 'Task failed'
  );
  updateProgressMessage(
    finalMessage,
    outcome === 'success'
      ? 'state-success'
      : outcome === 'cancelled'
        ? 'state-cancelled'
        : 'state-error'
  );

  const hideDelay = outcome === 'success' ? TIMING.PROGRESS_HIDE : TIMING.PROGRESS_HIDE + 600;
  const expectedRunId = runId ? String(runId) : '';
  setTimeout(() => {
    if (progressBar) {
      const activeRunId = progressBar.dataset.runId || '';
      if (expectedRunId && activeRunId && activeRunId !== expectedRunId) {
        return;
      }
    }
    if (progressBar) {
      progressBar.classList.remove('active');
      setProgressContainerState('');
      delete progressBar.dataset.runId;
      delete progressBar.dataset.command;
    }
    updateProgressMessage('Processing...', '');
    document.querySelectorAll('#progress-steps .step').forEach(stepEl => {
      stepEl.classList.remove('active', 'completed', 'error', 'cancelled');
    });
  }, hideDelay);

  setRunButtonIdle();
  setProcessing(false);
}

/**
 * Show file generated animation
 */
export function showFileGeneratedAnimation(outputPath) {
  document.querySelectorAll('.output-file').forEach(field => {
    if (field.value && outputPath.includes(getFileName(field.value))) {
      field.classList.add('file-generated');
      setTimeout(() => field.classList.remove('file-generated'), TIMING.ANIMATION_FEEDBACK);
    }
  });

  showSuccessToast(outputPath);
  showConfetti();
}

/**
 * Show success toast notification
 */
function showSuccessToast(outputPath) {
  const fileName = getFileName(outputPath);
  const toast = document.createElement('div');
  toast.className = 'success-toast';

  const icon = document.createElement('span');
  icon.className = 'toast-icon';
  icon.textContent = '✓';
  toast.appendChild(icon);

  const content = document.createElement('div');
  content.className = 'toast-content';

  const title = document.createElement('span');
  title.className = 'toast-title';
  title.textContent = 'File generated successfully!';
  content.appendChild(title);

  const path = document.createElement('span');
  path.className = 'toast-path';
  path.textContent = fileName;
  content.appendChild(path);

  toast.appendChild(content);
  document.body.appendChild(toast);
  setTimeout(() => toast.remove(), TIMING.TOAST_DURATION);
}

/**
 * Show confetti animation
 */
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
  setTimeout(() => container.remove(), TIMING.TOAST_DURATION);
}

/**
 * Cleanup event listeners (call on app teardown)
 */
export function cleanupProgressEvents() {
  if (progressUnlisten) {
    progressUnlisten();
    progressUnlisten = null;
  }
  if (completeUnlisten) {
    completeUnlisten();
    completeUnlisten = null;
  }
  if (logUnlisten) {
    logUnlisten();
    logUnlisten = null;
  }
  stopSystemUsageMonitoring();
}
