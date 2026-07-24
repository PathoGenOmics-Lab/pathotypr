// ============================================================================
// Dropzone Module - File Drag & Drop Handling
// ============================================================================

import { TIMING } from './config.js';
import { getCurrentDropTarget, setCurrentDropTarget } from './state.js';
import { openFileDialog, getFileMetadata } from './tauri.js';
import { validateFileExtension, formatFileSize, showExtensionError, getFileName } from './utils.js';
import { logMessage } from './console.js';

/**
 * Setup all dropzone elements
 */
export function setupDropzones() {
  document.querySelectorAll('.dropzone').forEach(dropzone => {
    const targetId = dropzone.dataset.target;
    const isMultiple = dropzone.dataset.multiple === 'true';
    const validExtensions = dropzone.dataset.extensions?.split(',') || [];
    const removeBtn = dropzone.querySelector('.file-remove');

    // Track current drop target on hover
    dropzone.addEventListener('mouseenter', () => setCurrentDropTarget(dropzone));
    dropzone.addEventListener('mouseleave', () => {
      if (getCurrentDropTarget() === dropzone) setCurrentDropTarget(null);
    });

    // Click to open file dialog
    dropzone.addEventListener('click', async (e) => {
      if (e.target.classList.contains('file-remove') || e.target.closest('.file-chip')) return;

      dropzone.style.transform = 'scale(0.98)';
      setTimeout(() => dropzone.style.transform = '', TIMING.SCROLL_SETTLE);

      try {
        const selected = await openFileDialog({ multiple: isMultiple });
        if (selected) {
          const paths = Array.isArray(selected) ? selected : [selected];

          // Validate extensions
          const invalidFiles = paths.filter(p => !validateFileExtension(p, validExtensions));
          if (invalidFiles.length > 0) {
            showExtensionError(dropzone, validExtensions);
            logMessage(`Invalid file type. Expected: ${validExtensions.join(', ')}`, 'error');
            return;
          }

          if (isMultiple) {
            await setDropzoneFiles(dropzone, targetId, paths);
          } else {
            await setDropzoneFile(dropzone, targetId, paths[0]);
          }
        }
      } catch (err) {
        logMessage('Error selecting file: ' + err, 'error');
      }
    });

    // Remove button handler
    if (removeBtn) {
      removeBtn.addEventListener('click', (e) => {
        e.stopPropagation();
        clearDropzone(dropzone, targetId);
        logMessage('File removed', 'info');
      });
    }
  });
}

/**
 * Setup Tauri native drag & drop
 */
export async function setupTauriDragDrop() {
  try {
    const { getCurrentWindow } = window.__TAURI__.window;
    const appWindow = getCurrentWindow();

    await appWindow.onDragDropEvent((event) => {
      if (event.payload.type === 'over') {
        showDragOverlay();
      } else if (event.payload.type === 'drop') {
        hideDragOverlay();
        handleSmartDrop(event.payload);
      } else if (event.payload.type === 'leave' || event.payload.type === 'cancel') {
        hideDragOverlay();
      }
    });

    logMessage('Drag & drop enabled', 'info');
  } catch (err) {
    console.warn('Tauri drag-drop setup failed:', err);
  }
}

/**
 * Show the drag overlay on the active panel
 */
function showDragOverlay() {
  if (document.getElementById('drag-overlay')) return;

  const overlay = document.createElement('div');
  overlay.id = 'drag-overlay';
  overlay.innerHTML = `
    <div class="drag-overlay-content">
      <svg class="drag-overlay-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
        <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/>
        <polyline points="7 10 12 15 17 10"/>
        <line x1="12" y1="15" x2="12" y2="3"/>
      </svg>
      <span class="drag-overlay-title">Drop your files</span>
      <span class="drag-overlay-hint">Files will be auto-assigned by type</span>
    </div>
  `;
  document.body.appendChild(overlay);
  // Trigger animation
  requestAnimationFrame(() => overlay.classList.add('visible'));
}

/**
 * Hide the drag overlay
 */
function hideDragOverlay() {
  const overlay = document.getElementById('drag-overlay');
  if (overlay) {
    overlay.classList.remove('visible');
    overlay.addEventListener('transitionend', () => overlay.remove(), { once: true });
    // Fallback removal
    setTimeout(() => overlay.remove(), 300);
  }
}

/**
 * Handle drop — always smart-routes by extension
 */
async function handleSmartDrop(payload) {
  const paths = payload.paths;
  if (!paths || paths.length === 0) return;

  const routed = await smartRouteFiles(paths);
  if (!routed) {
    logMessage(`Could not match ${paths.length} file(s) to any input field. Check file extensions.`, 'warning');
  }
}

/**
 * Set single file on dropzone
 */
export async function setDropzoneFile(dropzone, targetId, filePath) {
  const input = document.getElementById(targetId);
  const fileDisplay = dropzone.querySelector('.dropzone-file');
  const fileName = fileDisplay?.querySelector('.file-name');
  const fileMeta = fileDisplay?.querySelector('.file-meta');

  if (input) {
    input.value = filePath;
    // Drop any multi-file selection left over from a previous run: the form
    // builders prefer dataset.files over input.value, so a stale list would
    // silently override the single file the user just picked.
    delete input.dataset.files;
    input.dispatchEvent(new CustomEvent('change', { detail: { path: filePath } }));
  }

  const name = getFileName(filePath);
  if (fileName) fileName.textContent = name;
  if (fileDisplay) fileDisplay.classList.add('has-file');
  if (fileMeta) fileMeta.replaceChildren();

  // Load metadata in background — don't block UI
  getFileMetadata(filePath).then(metadata => {
    if (fileMeta && metadata) {
      const sizeStr = formatFileSize(metadata.size || 0);
      fileMeta.replaceChildren();
      const badge = document.createElement('span');
      badge.className = 'meta-badge';
      badge.textContent = sizeStr;
      fileMeta.appendChild(badge);
    }
  }).catch(() => {});

  logMessage(`Selected: ${name}`, 'success');
}

/**
 * Set multiple files on dropzone (accumulates with existing files)
 * @param {boolean} replace - If true, replace existing files instead of accumulating
 */
export async function setDropzoneFiles(dropzone, targetId, filePaths, replace = false) {
  const input = document.getElementById(targetId);
  const fileDisplay = dropzone.querySelector('.dropzone-file');
  const fileName = fileDisplay?.querySelector('.file-name');
  const fileMeta = fileDisplay?.querySelector('.file-meta');
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

  if (input) {
    input.value = allFiles.join(';');
    input.dataset.files = JSON.stringify(allFiles);
  }

  if (fileName) fileName.textContent = `${allFiles.length} files selected`;
  if (fileMeta) fileMeta.replaceChildren();
  if (fileDisplay) fileDisplay.classList.add('has-file');

  // Create file chips instantly, load metadata in background
  if (chipsContainer) {
    chipsContainer.innerHTML = '';
    for (const path of allFiles) {
      const name = getFileName(path);

      const chip = document.createElement('div');
      chip.className = 'file-chip';
      const chipName = document.createElement('span');
      chipName.className = 'chip-name';
      chipName.title = name;
      chipName.textContent = name;
      chip.appendChild(chipName);

      const chipSize = document.createElement('span');
      chipSize.className = 'chip-size';
      chip.appendChild(chipSize);

      // Load size in background
      getFileMetadata(path).then(metadata => {
        if (metadata) chipSize.textContent = formatFileSize(metadata.size || 0);
      }).catch(() => {});

      const removeBtn = document.createElement('button');
      removeBtn.type = 'button';
      removeBtn.className = 'chip-remove';
      removeBtn.dataset.path = path;
      removeBtn.textContent = '×';
      chip.appendChild(removeBtn);

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

/**
 * Clear dropzone
 */
export function clearDropzone(dropzone, targetId) {
  const input = document.getElementById(targetId);
  const fileDisplay = dropzone.querySelector('.dropzone-file');
  const chipsContainer = dropzone.querySelector('.file-chips');

  if (input) {
    input.value = '';
    delete input.dataset.files;
    input.dispatchEvent(new CustomEvent('change', { detail: { path: null } }));
  }

  if (fileDisplay) fileDisplay.classList.remove('has-file');
  if (chipsContainer) chipsContainer.innerHTML = '';
}

/**
 * Clear all loaded file inputs in the currently active panel.
 * Returns the number of cleared input fields.
 */
export function clearActivePanelDropzones() {
  const activePanel = document.querySelector('.panel.active');
  if (!activePanel || activePanel.id === 'panel-home') {
    logMessage('Select a tool panel first to clear its input files.', 'info');
    return 0;
  }

  const dropzones = Array.from(activePanel.querySelectorAll('.dropzone[data-target]'));
  if (dropzones.length === 0) {
    logMessage('No input file fields found in the active panel.', 'info');
    return 0;
  }

  let clearedCount = 0;
  for (const dropzone of dropzones) {
    const targetId = dropzone.dataset.target;
    if (!targetId) continue;
    const input = document.getElementById(targetId);
    const hasFiles = Boolean(
      String(input?.value || '').trim() ||
      (input?.dataset?.files && input.dataset.files !== '[]')
    );
    if (!hasFiles) continue;
    clearDropzone(dropzone, targetId);
    clearedCount += 1;
  }

  if (clearedCount > 0) {
    logMessage(`Cleared ${clearedCount} input field${clearedCount === 1 ? '' : 's'}.`, 'success');
  } else {
    logMessage('No loaded input files to clear in this panel.', 'info');
  }

  return clearedCount;
}

/**
 * Smart-route dropped files to the correct dropzone based on file extension.
 * Groups files by matching extension → dropzone, then assigns them.
 * Returns true if at least one file was routed successfully.
 */
async function smartRouteFiles(paths) {
  const activePanel = document.querySelector('.panel.active');
  if (!activePanel) return false;

  // Collect all visible dropzones in the active panel
  const dropzones = Array.from(activePanel.querySelectorAll('.dropzone[data-target]'))
    .filter(dz => !dz.closest('.hidden'));

  if (dropzones.length === 0) return false;

  // Group files by which dropzone they match
  const routeMap = new Map(); // dropzone → [paths]
  const unrouted = [];

  // A dropzone the user has already filled before this drop.
  const isPreFilled = (dz) => {
    const input = document.getElementById(dz.dataset.target);
    return Boolean(
      String(input?.value || '').trim() ||
      (input?.dataset?.files && input.dataset.files !== '[]')
    );
  };
  // A single-file dropzone can only take one file from this drop.
  const canAccept = (dz) => dz.dataset.multiple === 'true' || !routeMap.has(dz);

  for (const filePath of paths) {
    let matched = false;
    // Two passes: prefer a dropzone that is still empty, so dropping a second
    // file of the same type does not overwrite an input the user already
    // filled (e.g. landing a sample FASTA on top of the reference). Only if no
    // empty dropzone matches do we fall back to any matching one.
    for (const preferEmpty of [true, false]) {
      for (const dz of dropzones) {
        const extensions = dz.dataset.extensions?.split(',') || [];
        if (extensions.length === 0 || !validateFileExtension(filePath, extensions)) continue;
        if (!canAccept(dz)) continue;
        if (preferEmpty && isPreFilled(dz)) continue;
        if (!routeMap.has(dz)) routeMap.set(dz, []);
        routeMap.get(dz).push(filePath);
        matched = true;
        break;
      }
      if (matched) break;
    }
    if (!matched) unrouted.push(filePath);
  }

  if (routeMap.size === 0) return false;

  // Assign files to their matched dropzones in parallel
  const assignments = [...routeMap.entries()].map(([dz, filePaths]) => {
    const targetId = dz.dataset.target;
    const isMultiple = dz.dataset.multiple === 'true';

    const promise = isMultiple
      ? setDropzoneFiles(dz, targetId, filePaths)
      : setDropzoneFile(dz, targetId, filePaths[0]).then(() => {
          if (filePaths.length > 1) {
            logMessage(`Only first file assigned to ${targetId} (single file input).`, 'warning');
          }
        });

    dz.classList.add('drop-received');
    setTimeout(() => dz.classList.remove('drop-received'), TIMING.DROP_FEEDBACK);
    return promise;
  });
  await Promise.all(assignments);

  // Report results
  const routedCount = paths.length - unrouted.length;
  const dzNames = [...routeMap.keys()].map(dz => {
    const label = dz.closest('.form-group')?.querySelector('label')?.textContent?.trim() || dz.dataset.target;
    return label.replace(/\s*\?.*$/, ''); // Remove tooltip text
  });
  logMessage(`Auto-routed ${routedCount} file(s) to: ${dzNames.join(', ')}`, 'success');

  if (unrouted.length > 0) {
    logMessage(`${unrouted.length} file(s) could not be matched: ${unrouted.map(getFileName).join(', ')}`, 'warning');
  }

  return true;
}
