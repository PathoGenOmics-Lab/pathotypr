// ============================================================================
// Console Drawer Module
// ============================================================================

import { incrementConsoleMessageCount, resetConsoleMessageCount } from './state.js';
import { getFileName, splitClickablePathTokens } from './utils.js';
import { openFileLocation } from './tauri.js';

let isSetup = false;
const AUTO_OPEN_DRAWER_HEIGHT = 220;
const MIN_DRAWER_HEIGHT = 150;

function clampDrawerHeight(height, maxHeight = window.innerHeight * 0.7) {
  const parsed = Number.parseInt(String(height || ''), 10);
  if (!Number.isFinite(parsed)) return null;
  const clamped = Math.min(Math.max(parsed, MIN_DRAWER_HEIGHT), Math.floor(maxHeight));
  return `${clamped}px`;
}

/**
 * Setup console drawer functionality
 */
export function setupConsoleDrawer() {
  if (isSetup) return;
  isSetup = true;

  const drawer = document.getElementById('console-drawer');
  const backdrop = document.getElementById('drawer-backdrop');
  const toggleBtn = document.getElementById('btn-console-toggle');
  const closeBtn = document.getElementById('btn-close-drawer');
  const handle = document.getElementById('drawer-handle');

  if (toggleBtn) toggleBtn.addEventListener('click', toggleConsoleDrawer);
  if (closeBtn) closeBtn.addEventListener('click', closeConsoleDrawer);
  if (backdrop) backdrop.addEventListener('click', closeConsoleDrawer);

  // Drawer resize functionality
  if (handle && drawer) {
    let startY, startHeight;

    const onMouseMove = (e) => {
      const diff = startY - e.clientY;
      const newHeight = Math.min(Math.max(startHeight + diff, MIN_DRAWER_HEIGHT), window.innerHeight * 0.7);
      drawer.style.height = `${newHeight}px`;
    };

    const onMouseUp = () => {
      document.removeEventListener('mousemove', onMouseMove);
      document.removeEventListener('mouseup', onMouseUp);
      document.body.style.cursor = '';
      localStorage.setItem('pathotypr-drawer-height', drawer.style.height);
    };

    handle.addEventListener('mousedown', (e) => {
      startY = e.clientY;
      startHeight = drawer.offsetHeight;
      document.body.style.cursor = 'row-resize';
      document.addEventListener('mousemove', onMouseMove);
      document.addEventListener('mouseup', onMouseUp);
      e.preventDefault();
    });

    // Restore saved height
    const savedHeight = localStorage.getItem('pathotypr-drawer-height');
    const restoredHeight = clampDrawerHeight(savedHeight);
    if (restoredHeight) drawer.style.height = restoredHeight;
  }

  // Escape key to close
  document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape' && drawer?.classList.contains('open')) {
      closeConsoleDrawer();
    }
  });
}

/**
 * Toggle console drawer visibility
 */
export function toggleConsoleDrawer() {
  const drawer = document.getElementById('console-drawer');
  const backdrop = document.getElementById('drawer-backdrop');
  if (drawer.classList.contains('open')) {
    closeConsoleDrawer();
  } else {
    drawer.classList.add('open');
    backdrop.classList.add('visible');
    resetBadge();
  }
}

/**
 * Close console drawer
 */
export function closeConsoleDrawer() {
  document.getElementById('console-drawer')?.classList.remove('open');
  document.getElementById('drawer-backdrop')?.classList.remove('visible');
}

/**
 * Open console drawer
 */
export function openConsoleDrawer() {
  document.getElementById('console-drawer')?.classList.add('open');
  document.getElementById('drawer-backdrop')?.classList.add('visible');
  resetBadge();
}

/**
 * Open console only if currently closed (avoids repeated forced-open behavior).
 */
export function openConsoleDrawerIfClosed(compact = false) {
  const drawer = document.getElementById('console-drawer');
  if (!drawer || drawer.classList.contains('open')) return false;
  if (compact) {
    const compactHeight = clampDrawerHeight(AUTO_OPEN_DRAWER_HEIGHT, window.innerHeight * 0.45);
    if (compactHeight) drawer.style.height = compactHeight;
  }
  openConsoleDrawer();
  return true;
}

/**
 * Update console badge with message count
 */
export function updateConsoleBadge() {
  const badge = document.getElementById('console-badge');
  const drawer = document.getElementById('console-drawer');
  if (badge && drawer && !drawer.classList.contains('open')) {
    const count = incrementConsoleMessageCount();
    badge.textContent = count > 99 ? '99+' : count;
    badge.classList.add('has-messages');
  }
}

/**
 * Reset console badge
 */
function resetBadge() {
  const badge = document.getElementById('console-badge');
  resetConsoleMessageCount();
  if (badge) {
    badge.classList.remove('has-messages');
    badge.textContent = '';
  }
}

/**
 * Log a message to the console
 */
export function logMessage(message, type = 'info') {
  const outputContent = document.getElementById('output-content');
  if (!outputContent) return;

  const entry = document.createElement('div');
  entry.className = `log-entry log-${type}`;

  const timestamp = new Date().toLocaleTimeString();
  const icon = type === 'success' ? '✓' : type === 'error' ? '✗' : type === 'warning' ? '⚠' : '→';

  entry.appendChild(document.createTextNode(`[${timestamp}] ${icon} `));

  // Render message safely while preserving clickable path tokens.
  for (const part of splitClickablePathTokens(message)) {
    if (part.type === 'path') {
      const pathEl = document.createElement('span');
      pathEl.className = 'clickable-path';
      pathEl.dataset.path = part.value;
      pathEl.textContent = getFileName(part.value) || part.value;
      pathEl.title = 'Click to open in folder. Shift+Click to copy path.';
      pathEl.addEventListener('click', async (event) => {
        const targetPath = pathEl.dataset.path || '';
        if (!targetPath) return;
        if (event.shiftKey) {
          try {
            await navigator.clipboard.writeText(targetPath);
            logMessage('Path copied to clipboard', 'info');
          } catch {
            logMessage('Failed to copy path to clipboard', 'warning');
          }
          return;
        }

        try {
          const opened = await openFileLocation(targetPath);
          if (opened) {
            logMessage(`Opened location: ${getFileName(targetPath) || targetPath}`, 'info');
          } else {
            await navigator.clipboard.writeText(targetPath);
            logMessage('Path copied to clipboard', 'info');
          }
        } catch {
          try {
            await navigator.clipboard.writeText(targetPath);
            logMessage('Path copied to clipboard', 'info');
          } catch {
            logMessage('Failed to copy path to clipboard', 'warning');
          }
        }
      });
      entry.appendChild(pathEl);
    } else {
      entry.appendChild(document.createTextNode(part.value));
    }
  }

  const placeholder = outputContent.querySelector('.placeholder');
  if (placeholder) placeholder.remove();

  outputContent.appendChild(entry);
  outputContent.scrollTop = outputContent.scrollHeight;
  updateConsoleBadge();
}

/**
 * Clear console output
 */
export function clearOutput() {
  const outputContent = document.getElementById('output-content');
  if (outputContent) {
    outputContent.innerHTML = '<p class="placeholder">Output will appear here...</p>';
  }
  resetBadge();
}
