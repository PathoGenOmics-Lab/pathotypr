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
        highlightDropzone(dropzoneFromPosition(event.payload.position));
      } else if (event.payload.type === 'drop') {
        const aimed = dropzoneFromPosition(event.payload.position);
        hideDragOverlay();
        highlightDropzone(null);
        handleDrop(event.payload, aimed);
      } else if (event.payload.type === 'leave' || event.payload.type === 'cancel') {
        hideDragOverlay();
        highlightDropzone(null);
      }
    });

    logMessage('Drag & drop enabled', 'info');
  } catch (err) {
    console.warn('Tauri drag-drop setup failed:', err);
  }
}

// Which dropzone the pointer is over during a drag. Several fields in a panel accept the
// same extension (classify takes .fasta for both the reference and the samples, and .tsv
// for both the markers and the input list), so extension alone cannot say where a file
// belongs. Aiming does: drop on the field you mean and it goes there.
let hoveredDropzone = null;

/**
 * The dropzone under a drag pointer, or null.
 *
 * The position is typed as physical in Tauri, but on this platform it arrives already in
 * CSS pixels, so it is hit-tested as-is and only scaled by the device pixel ratio as a
 * fallback, for a platform that does report physical pixels. The overlay does not
 * interfere with the hit test: it is pointer-events: none.
 */
function dropzoneFromPosition(position) {
  if (!position || typeof position.x !== 'number' || typeof position.y !== 'number') return null;

  const ratio = window.devicePixelRatio || 1;
  const candidates = ratio === 1
    ? [[position.x, position.y]]
    : [[position.x, position.y], [position.x / ratio, position.y / ratio]];

  for (const [x, y] of candidates) {
    if (x < 0 || y < 0 || x > window.innerWidth || y > window.innerHeight) continue;
    const dropzone = document.elementFromPoint(x, y)?.closest('.dropzone[data-target]');
    // Only fields on the panel in view can take a drop.
    if (dropzone?.closest('.panel.active')) return dropzone;
  }
  return null;
}

/** Human-readable name of a dropzone, taken from its form label. */
function dropzoneLabel(dropzone) {
  const label = dropzone.closest('.form-group')?.querySelector('label')?.textContent?.trim();
  return (label || dropzone.dataset.target).replace(/\s*\?.*$/, '');   // drop tooltip text
}

/** Mark the dropzone being aimed at, so it is clear where the file will land. */
function highlightDropzone(dropzone) {
  if (hoveredDropzone === dropzone) return;   // don't touch the DOM on every drag event
  hoveredDropzone?.classList.remove('drag-over');
  dropzone?.classList.add('drag-over');
  hoveredDropzone = dropzone;
}

// The overlay is tracked here rather than looked up in the DOM, because it outlives a
// single drag: Tauri emits `over` continuously while the cursor moves, and the element
// lingers during its fade-out. Both cases have to be told apart from a fresh drag.
let dragOverlay = null;
let dragOverlayRemoval = null;

/**
 * Show the drag overlay on the active panel
 */
function showDragOverlay() {
  // A fade-out may be in flight from a drag that just left; keep the element.
  if (dragOverlayRemoval) {
    clearTimeout(dragOverlayRemoval);
    dragOverlayRemoval = null;
  }

  if (dragOverlay) {
    // Repeat `over` events, or a drag that left and came straight back. Re-adding a
    // class already present is a no-op, so the animation is not restarted on every
    // mouse move.
    dragOverlay.classList.add('visible');
    return;
  }

  // Just a tint. The card that used to sit here covered the fields the drag is aimed at,
  // and its promise of automatic assignment is no longer what happens.
  const overlay = document.createElement('div');
  overlay.id = 'drag-overlay';
  document.body.appendChild(overlay);
  dragOverlay = overlay;
  // Trigger animation
  requestAnimationFrame(() => overlay.classList.add('visible'));
}

/**
 * Hide the drag overlay
 */
function hideDragOverlay() {
  if (!dragOverlay) return;
  const overlay = dragOverlay;
  overlay.classList.remove('visible');
  // Removal is deferred so the fade-out plays, and is cancelled by showDragOverlay if
  // the drag returns first. No transitionend listener: fading back in would fire it and
  // tear down an overlay that is on its way back.
  dragOverlayRemoval = setTimeout(() => {
    overlay.remove();
    if (dragOverlay === overlay) dragOverlay = null;
    dragOverlayRemoval = null;
  }, 300);
}

/**
 * Handle a drop: files land on the field they were dropped on, and anything dropped
 * outside a field — or that the aimed field cannot take — falls back to routing by
 * extension.
 */
async function handleDrop(payload, aimedDropzone) {
  const paths = payload.paths;
  if (!paths || paths.length === 0) return;

  let remaining = paths;

  if (aimedDropzone) {
    const targetId = aimedDropzone.dataset.target;
    const extensions = aimedDropzone.dataset.extensions?.split(',') || [];
    const accepted = paths.filter(p => validateFileExtension(p, extensions));
    const rejected = paths.filter(p => !accepted.includes(p));

    if (accepted.length > 0) {
      const isMultiple = aimedDropzone.dataset.multiple === 'true';
      if (isMultiple) {
        await setDropzoneFiles(aimedDropzone, targetId, accepted);
        remaining = rejected;
      } else {
        await setDropzoneFile(aimedDropzone, targetId, accepted[0]);
        if (accepted.length > 1) {
          logMessage(`${dropzoneLabel(aimedDropzone)} takes one file; routing the rest by type`, 'info');
        }
        // Anything the field could not hold is routed rather than dropped on the floor.
        remaining = [...accepted.slice(1), ...rejected];
      }
      aimedDropzone.classList.add('drop-received');
      setTimeout(() => aimedDropzone.classList.remove('drop-received'), TIMING.DROP_FEEDBACK);
    } else if (rejected.length > 0) {
      // Dropped on a field that cannot take these files. Say so, then let the fallback
      // try, rather than silently sending them somewhere the user did not point at.
      logMessage(
        `${dropzoneLabel(aimedDropzone)} only accepts ${extensions.join(', ')}`,
        'warning'
      );
    }
  }

  if (remaining.length === 0) return;

  const routed = await smartRouteFiles(remaining);
  if (!routed) {
    logMessage(`Could not match ${remaining.length} file(s) to any input field. Check file extensions.`, 'warning');
  }
}

/**
 * Set single file on dropzone
 */
export async function setDropzoneFile(dropzone, targetId, filePath, quiet = false) {
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
  // A dropzone that accepts several files may still be showing chips from an earlier
  // multi-file selection. The list was just dropped from the input, so the chips would
  // claim files that are no longer selected.
  const staleChips = dropzone.querySelector('.file-chips');
  if (staleChips) staleChips.innerHTML = '';

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

  if (!quiet) logMessage(`Selected: ${name}`, 'success');
}

/**
 * Set multiple files on dropzone (accumulates with existing files)
 * @param {boolean} replace - If true, replace existing files instead of accumulating
 */
export async function setDropzoneFiles(dropzone, targetId, filePaths, replace = false, quiet = false) {
  const input = document.getElementById(targetId);
  const fileDisplay = dropzone.querySelector('.dropzone-file');
  const fileName = fileDisplay?.querySelector('.file-name');
  const fileMeta = fileDisplay?.querySelector('.file-meta');
  const chipsContainer = dropzone.querySelector('.file-chips');

  // Get existing files and merge with new ones (avoiding duplicates)
  let allFiles = filePaths;
  let existingCount = 0;
  if (!replace && input?.dataset.files) {
    try {
      const existingFiles = JSON.parse(input.dataset.files);
      existingCount = existingFiles.length;
      // Merge: add new files that aren't already in the list
      const newFiles = filePaths.filter(p => !existingFiles.includes(p));
      allFiles = [...existingFiles, ...newFiles];
    } catch (e) {
      // If parsing fails, just use the new files
      allFiles = filePaths;
      existingCount = 0;
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
          // Replace mode, and quiet: re-rendering the shorter list is a removal, and
          // letting it announce itself would report the leftovers as newly added.
          setDropzoneFiles(dropzone, targetId, remaining, true, true);
        } else {
          clearDropzone(dropzone, targetId);
        }
        logMessage(`Removed: ${getFileName(pathToRemove)}`, 'info');
      });
    });
  }

  // Files already in the list are dropped by the merge above, so report what actually
  // landed rather than how many were handed in.
  const addedCount = replace ? allFiles.length : allFiles.length - existingCount;
  if (addedCount > 0 && !quiet) {
    logMessage(`Added ${addedCount} file(s). Total: ${allFiles.length}`, 'success');
  } else if (addedCount === 0 && !quiet) {
    logMessage(`Already selected: ${filePaths.map(getFileName).join(', ')}`, 'info');
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
 * Route dropped files by extension, but only where the extension is unambiguous.
 *
 * A panel can offer several fields for the same extension — classify takes .fasta for
 * both the reference and the samples, and .tsv for both the markers and the input list.
 * Guessing between them is worse than not routing: a reference genome silently loaded as
 * a sample is a wrong run, not a small annoyance. So a file that matches more than one
 * field is left alone and reported, and the user drops it on the field they mean.
 *
 * Returns true if at least one file was routed.
 */
async function smartRouteFiles(paths) {
  const activePanel = document.querySelector('.panel.active');
  if (!activePanel) return false;

  const dropzones = Array.from(activePanel.querySelectorAll('.dropzone[data-target]'))
    .filter(dz => !dz.closest('.hidden'));
  if (dropzones.length === 0) return false;

  const routeMap = new Map();      // dropzone -> [paths]
  const ambiguous = [];            // [{ path, candidates }]
  const unmatched = [];
  const overflow = [];             // matched a single-file field already taken by this drop

  for (const filePath of paths) {
    const candidates = dropzones.filter(dz => {
      const extensions = dz.dataset.extensions?.split(',') || [];
      return extensions.length > 0 && validateFileExtension(filePath, extensions);
    });

    if (candidates.length === 0) {
      unmatched.push(filePath);
    } else if (candidates.length > 1) {
      ambiguous.push({ path: filePath, candidates });
    } else {
      const dz = candidates[0];
      const isMultiple = dz.dataset.multiple === 'true';
      if (!isMultiple && routeMap.has(dz)) {
        overflow.push({ path: filePath, dropzone: dz });
      } else {
        if (!routeMap.has(dz)) routeMap.set(dz, []);
        routeMap.get(dz).push(filePath);
      }
    }
  }

  if (routeMap.size > 0) {
    const assignments = [...routeMap.entries()].map(([dz, filePaths]) => {
      const targetId = dz.dataset.target;
      const promise = dz.dataset.multiple === 'true'
        ? setDropzoneFiles(dz, targetId, filePaths)
        : setDropzoneFile(dz, targetId, filePaths[0]);
      dz.classList.add('drop-received');
      setTimeout(() => dz.classList.remove('drop-received'), TIMING.DROP_FEEDBACK);
      return promise;
    });
    await Promise.all(assignments);
    logMessage(
      `Auto-routed ${paths.length - ambiguous.length - unmatched.length - overflow.length} ` +
      `file(s) to: ${[...routeMap.keys()].map(dropzoneLabel).join(', ')}`,
      'success'
    );
  }

  for (const { path, candidates } of ambiguous) {
    logMessage(
      `${getFileName(path)} fits ${candidates.map(dropzoneLabel).join(' or ')}. ` +
      `Drop it on the field you want.`,
      'warning'
    );
  }
  for (const { path, dropzone } of overflow) {
    logMessage(`${dropzoneLabel(dropzone)} takes one file; ${getFileName(path)} was not loaded.`, 'warning');
  }
  if (unmatched.length > 0) {
    logMessage(
      `${unmatched.length} file(s) could not be matched: ${unmatched.map(getFileName).join(', ')}`,
      'warning'
    );
  }

  return routeMap.size > 0;
}
