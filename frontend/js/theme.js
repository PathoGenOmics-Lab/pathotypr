// ============================================================================
// Theme Module
// ============================================================================

import { refreshChartsForTheme } from './visualization.js';

const SCHEME_LIGHT_SRC = 'pathotypr_scheme.webp';
const SCHEME_DARK_SRC = 'pathotypr_scheme-dark.webp';

/**
 * Setup theme toggle button
 */
export function setupThemeToggle() {
  const themeToggle = document.getElementById('theme-toggle');
  if (themeToggle) {
    themeToggle.addEventListener('click', toggleTheme);
  }
}

/**
 * Toggle between light and dark theme
 */
export function toggleTheme() {
  const html = document.documentElement;
  const currentTheme = html.getAttribute('data-theme');
  const newTheme = currentTheme === 'dark' ? 'light' : 'dark';
  html.setAttribute('data-theme', newTheme);
  updateThemeUI(newTheme);
  localStorage.setItem('pathotypr-theme', newTheme);

  // Refresh charts to update label colors for new theme
  setTimeout(() => {
    refreshChartsForTheme();
  }, 50);
}

/**
 * Update theme UI elements
 */
function updateThemeUI(theme) {
  const themeLabel = document.querySelector('.theme-label');
  if (themeLabel) {
    themeLabel.innerHTML = theme === 'dark' ? 'DARK<br>MODE' : 'LIGHT<br>MODE';
  }
  const themeToggle = document.getElementById('theme-toggle');
  if (themeToggle) {
    themeToggle.setAttribute('aria-checked', theme === 'dark' ? 'true' : 'false');
  }

  updateSchemeImageForTheme(theme);
}

function updateSchemeImageForTheme(theme) {
  const schemeSrc = theme === 'dark' ? SCHEME_DARK_SRC : SCHEME_LIGHT_SRC;

  const schemeImg = document.getElementById('scheme-img');
  if (schemeImg) {
    schemeImg.setAttribute('src', schemeSrc);
  }

  const lightboxImg = document.querySelector('.scheme-lightbox-img');
  if (lightboxImg) {
    lightboxImg.setAttribute('src', schemeSrc);
  }
}

/**
 * Load saved theme from localStorage
 */
export function loadSavedTheme() {
  const savedTheme = localStorage.getItem('pathotypr-theme') || 'light';
  document.documentElement.setAttribute('data-theme', savedTheme);
  updateThemeUI(savedTheme);
}
