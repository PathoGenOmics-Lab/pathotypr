// ============================================================================
// Navigation Module
// ============================================================================

import { panelInfo, categoryMap } from './config.js';
import { destroyToolCharts } from './state.js';
import { refreshVisualizationForTool } from './visualization.js';

const PANEL_TOOL_ID_MAP = {
  'split-fastq': 'splitfq'
};

function resolveToolId(panelId) {
  if (!panelId) return null;
  return PANEL_TOOL_ID_MAP[panelId] || panelId;
}

/**
 * Setup navigation handlers
 */
export function setupNavigation() {
  document.querySelectorAll('.nav-item').forEach(item => {
    item.addEventListener('click', () => {
      navigateToPanel(item.dataset.panel);
      closeMobileSidebar();
    });
  });

  // Mobile hamburger menu
  const hamburger = document.getElementById('btn-hamburger');
  const overlay = document.getElementById('sidebar-overlay');
  if (hamburger) {
    hamburger.addEventListener('click', toggleMobileSidebar);
  }
  if (overlay) {
    overlay.addEventListener('click', closeMobileSidebar);
  }
}

function toggleMobileSidebar() {
  const sidebar = document.getElementById('sidebar');
  const overlay = document.getElementById('sidebar-overlay');
  sidebar?.classList.toggle('mobile-open');
  overlay?.classList.toggle('visible');
}

function closeMobileSidebar() {
  const sidebar = document.getElementById('sidebar');
  const overlay = document.getElementById('sidebar-overlay');
  sidebar?.classList.remove('mobile-open');
  overlay?.classList.remove('visible');
}

/**
 * Setup quick action cards
 */
export function setupQuickActions() {
  document.querySelectorAll('.action-card[data-goto]').forEach(card => {
    card.addEventListener('click', () => {
      navigateToPanel(card.dataset.goto);
    });
  });
}

/**
 * Navigate to a specific panel
 */
export function navigateToPanel(panelId) {
  // Destroy charts from the panel we're leaving to free memory
  const prevPanel = getActivePanel();
  const prevToolId = resolveToolId(prevPanel);
  const currentToolId = resolveToolId(panelId);
  if (prevPanel && prevPanel !== panelId) {
    destroyToolCharts(prevToolId);
  }

  // Update nav item active state
  document.querySelectorAll('.nav-item').forEach(nav => {
    nav.classList.toggle('active', nav.dataset.panel === panelId);
  });

  // Update panel visibility
  document.querySelectorAll('.panel').forEach(panel => {
    panel.classList.toggle('active', panel.id === `panel-${panelId}`);
  });

  updateHeader(panelId);
  updateCategoryAccent(panelId);

  // Rebuild charts if this panel has an open visualization.
  requestAnimationFrame(() => {
    if (!currentToolId) return;
    refreshVisualizationForTool(currentToolId);
  });
}

/**
 * Update header with panel information
 */
function updateHeader(panelId) {
  const info = panelInfo[panelId];
  if (info) {
    const headerIcon = document.getElementById('header-icon');
    const titleText = document.getElementById('panel-title-text');
    const description = document.getElementById('panel-description');

    if (headerIcon) headerIcon.innerHTML = info.iconPath;
    if (titleText) titleText.textContent = info.title;
    if (description) description.textContent = info.description;
  }

  // Show/hide quick-load buttons for genotyping and predict panels
  const quickLoad = document.getElementById('header-quick-load');
  if (quickLoad) {
    const isGenotyping = panelId === 'classify' || panelId === 'split-fastq';
    const isPredict = panelId === 'predict';
    quickLoad.classList.toggle('visible', isGenotyping || isPredict);

    // Show only relevant buttons per panel
    const markerBtns = quickLoad.querySelectorAll('.btn-glass-lineage, .btn-glass-dr');
    const modelBtn = quickLoad.querySelector('.btn-glass-model');
    markerBtns.forEach(b => b.style.display = isGenotyping ? '' : 'none');
    if (modelBtn) modelBtn.style.display = isPredict ? '' : 'none';
  }
}

/**
 * Update body category accent for styling
 */
function updateCategoryAccent(panelId) {
  const category = categoryMap[panelId];
  if (category) {
    document.body.setAttribute('data-category', category);
  } else {
    document.body.removeAttribute('data-category');
  }
}

/**
 * Get current active panel ID
 */
export function getActivePanel() {
  const activePanel = document.querySelector('.panel.active');
  return activePanel ? activePanel.id.replace('panel-', '') : null;
}
