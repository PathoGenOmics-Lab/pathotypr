// ============================================================================
// Pathotypr GUI - Main Entry Point
// Modular Frontend Application (Tauri v2)
// ============================================================================

// Import modules
import { isTauriAvailable, getVersion, openExternalUrl } from './tauri.js';
import { setupNavigation, setupQuickActions } from './navigation.js';
import { setupThemeToggle, loadSavedTheme } from './theme.js';
import { setupConsoleDrawer, logMessage, clearOutput } from './console.js';
import { setupDropzones, setupTauriDragDrop, clearActivePanelDropzones } from './dropzone.js';
import { startProgress, stopProgress, setupProgressEvents } from './progress.js';
import { setupResultsModal } from './results.js';
import { setupAllVisualizations } from './visualization.js';
import {
  setupForms,
  setupSaveButtons,
  setupRadioHandlers,
  setupDefaultButtons,
  setupDemoButton,
  setupRunButton,
  loadSavedFormState,
  setupFormStatePersistence
} from './forms.js';

/**
 * Initialize the application
 */
async function initApp() {
  console.log('Initializing Pathotypr GUI...');
  console.log('Tauri available:', isTauriAvailable());

  if (!isTauriAvailable()) {
    console.error('Tauri API not available!');
    logMessage('Error: Tauri API not loaded. Please restart the application.', 'error');
    return;
  }

  // Setup all modules
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
  setupAllVisualizations();
  setupDemoButton();
  setupDefaultButtons();
  setupFormStatePersistence();

  // Setup real-time progress events from backend
  await setupProgressEvents();

  // Setup clear output button
  const clearBtn = document.getElementById('btn-clear-output');
  if (clearBtn) {
    clearBtn.addEventListener('click', clearOutput);
  }

  const clearInputsBtn = document.getElementById('btn-clear-inputs');
  if (clearInputsBtn) {
    clearInputsBtn.addEventListener('click', () => {
      clearActivePanelDropzones();
    });
  }

  // Load saved state
  loadVersion();
  loadSavedTheme();
  loadSavedFormState();

  setupImageFallbacks();
  setupSchemeLightbox();
  setupExternalLinkHandling();

  logMessage('Pathotypr GUI ready', 'success');
}

/**
 * Load and display version
 */
async function loadVersion() {
  const version = await getVersion();
  const versionEl = document.getElementById('app-version');
  if (version && versionEl) {
    versionEl.textContent = `v${version}`;
  }
}

/**
 * Setup scheme image lightbox
 */
function setupSchemeLightbox() {
  const schemeImg = document.getElementById('scheme-img');
  const lightbox = document.getElementById('scheme-lightbox');
  if (!schemeImg || !lightbox) return;

  const backdrop = lightbox.querySelector('.scheme-lightbox-backdrop');
  const closeBtn = lightbox.querySelector('.scheme-lightbox-close');

  const openLightbox = () => lightbox.classList.remove('hidden');
  const closeLightbox = () => lightbox.classList.add('hidden');

  schemeImg.addEventListener('click', openLightbox);
  if (backdrop) backdrop.addEventListener('click', closeLightbox);
  if (closeBtn) closeBtn.addEventListener('click', closeLightbox);
  document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape' && !lightbox.classList.contains('hidden')) {
      closeLightbox();
    }
  });
}

/**
 * Handle image loading failures without inline event handlers.
 */
function setupImageFallbacks() {
  const appLogo = document.getElementById('app-logo');
  if (appLogo) {
    appLogo.addEventListener('error', () => {
      appLogo.style.display = 'none';
    }, { once: true });
  }

  const schemeImg = document.getElementById('scheme-img');
  if (schemeImg) {
    schemeImg.addEventListener('error', () => {
      const schemeContainer = schemeImg.closest('.home-scheme');
      if (schemeContainer) schemeContainer.style.display = 'none';
    }, { once: true });
  }
}

/**
 * Route external links through backend so they open in the system browser.
 */
function setupExternalLinkHandling() {
  document.addEventListener('click', async (event) => {
    if (!(event.target instanceof Element)) return;
    const anchor = event.target.closest('a[href]');
    if (!anchor) return;

    const href = String(anchor.getAttribute('href') || '').trim();
    if (!href || !/^(https?:|mailto:)/i.test(href)) return;

    event.preventDefault();
    const opened = await openExternalUrl(href);
    if (!opened) {
      window.open(href, '_blank', 'noopener,noreferrer');
    }
  });
}

// Wait for DOM to be ready, then initialize
document.addEventListener('DOMContentLoaded', async () => {
  // Small delay to ensure Tauri is ready
  await new Promise(resolve => setTimeout(resolve, 100));
  await initApp();
});
