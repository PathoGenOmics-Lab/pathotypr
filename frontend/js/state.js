// ============================================================================
// Global State Management
// ============================================================================

// Global state object
const state = {
  currentDropTarget: null,
  isProcessing: false,
  consoleMessageCount: 0,
  lastOutputPath: null,
  progressInterval: null,
  panelResultsData: {},
  toolCharts: {},
  runSequence: 0,
  currentRun: {
    id: 0,
    command: null,
    panelId: null,
    status: 'idle', // idle | running | cancelling
    cancelRequested: false,
    startedAt: 0,
    finishedAt: 0,
    outcome: null // success | error | cancelled
  }
};

// State getters
export function getCurrentDropTarget() {
  return state.currentDropTarget;
}

export function isProcessing() {
  return state.isProcessing;
}

export function getConsoleMessageCount() {
  return state.consoleMessageCount;
}

export function getLastOutputPath() {
  return state.lastOutputPath;
}

export function getProgressInterval() {
  return state.progressInterval;
}

export function getPanelResultsData(viewerId) {
  return state.panelResultsData[viewerId];
}

export function getAllPanelResultsData() {
  return state.panelResultsData;
}

export function getToolCharts(toolId) {
  return state.toolCharts[toolId];
}

export function getCurrentRun() {
  return state.currentRun;
}

// State setters
export function setCurrentDropTarget(target) {
  state.currentDropTarget = target;
}

export function setProcessing(value) {
  state.isProcessing = value;
}

export function startRunContext(command, panelId) {
  state.runSequence += 1;
  state.currentRun = {
    id: state.runSequence,
    command: command || null,
    panelId: panelId || null,
    status: 'running',
    cancelRequested: false,
    startedAt: Date.now(),
    finishedAt: 0,
    outcome: null
  };
  state.isProcessing = true;
  return state.currentRun.id;
}

export function isRunCurrent(runId) {
  return state.currentRun.id !== 0 && state.currentRun.id === runId;
}

export function requestRunCancellation(runId) {
  if (!isRunCurrent(runId)) return false;
  state.currentRun.cancelRequested = true;
  state.currentRun.status = 'cancelling';
  return true;
}

export function clearRunCancellationRequest(runId) {
  if (!isRunCurrent(runId)) return false;
  state.currentRun.cancelRequested = false;
  state.currentRun.status = 'running';
  return true;
}

export function wasCancellationRequested(runId) {
  return isRunCurrent(runId) && state.currentRun.cancelRequested;
}

export function finishRunContext(runId, outcome = 'error') {
  if (!isRunCurrent(runId)) return;
  state.currentRun.finishedAt = Date.now();
  state.currentRun.outcome = outcome;
  state.currentRun.status = 'idle';
  state.isProcessing = false;
}

export function clearRunContext(runId) {
  if (!isRunCurrent(runId)) return;
  state.currentRun = {
    id: 0,
    command: null,
    panelId: null,
    status: 'idle',
    cancelRequested: false,
    startedAt: 0,
    finishedAt: 0,
    outcome: null
  };
}

export function incrementConsoleMessageCount() {
  state.consoleMessageCount++;
  return state.consoleMessageCount;
}

export function resetConsoleMessageCount() {
  state.consoleMessageCount = 0;
}

export function setLastOutputPath(path) {
  state.lastOutputPath = path;
}

export function setProgressInterval(interval) {
  state.progressInterval = interval;
}

export function clearProgressInterval() {
  if (state.progressInterval) {
    clearInterval(state.progressInterval);
    state.progressInterval = null;
  }
}

export function setPanelResultsData(viewerId, data) {
  state.panelResultsData[viewerId] = data;
}

export function clearPanelResultsData(viewerId) {
  delete state.panelResultsData[viewerId];
}

export function initToolCharts(toolId) {
  state.toolCharts[toolId] = { donut: null, bar: null };
}

export function setToolChart(toolId, chartType, chart) {
  if (!state.toolCharts[toolId]) {
    state.toolCharts[toolId] = { donut: null, bar: null };
  }
  state.toolCharts[toolId][chartType] = chart;
}

export function destroyToolCharts(toolId) {
  if (state.toolCharts[toolId]?.donut) {
    state.toolCharts[toolId].donut.destroy();
    state.toolCharts[toolId].donut = null;
  }
  if (state.toolCharts[toolId]?.bar) {
    state.toolCharts[toolId].bar.destroy();
    state.toolCharts[toolId].bar = null;
  }
}

export function destroyAllCharts() {
  for (const toolId of Object.keys(state.toolCharts)) {
    destroyToolCharts(toolId);
  }
}
