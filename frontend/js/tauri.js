// ============================================================================
// Tauri API Helpers
// ============================================================================

/**
 * Check if Tauri is available
 */
export function isTauriAvailable() {
  return typeof window.__TAURI__ !== 'undefined';
}

/**
 * Invoke a Tauri command
 */
export async function tauriInvoke(cmd, args = {}) {
  try {
    return await window.__TAURI__.core.invoke(cmd, args);
  } catch (err) {
    console.error(`Tauri invoke error (${cmd}):`, err);
    throw err;
  }
}

/**
 * Open file dialog
 */
export async function openFileDialog(options = {}) {
  try {
    return await window.__TAURI__.dialog.open(options);
  } catch (err) {
    console.error('Open dialog error:', err);
    throw err;
  }
}

/**
 * Save file dialog
 */
export async function saveFileDialog(options = {}) {
  try {
    return await window.__TAURI__.dialog.save(options);
  } catch (err) {
    console.error('Save dialog error:', err);
    throw err;
  }
}

/**
 * Show confirmation dialog
 */
export async function confirmDialog(message, title = 'Confirm') {
  try {
    return await window.__TAURI__.dialog.ask(message, { title, kind: 'warning' });
  } catch (err) {
    console.error('Confirm dialog error:', err);
    return false;
  }
}

/**
 * Read text file content
 */
export async function readTextFile(path) {
  try {
    const { invoke } = window.__TAURI__.core;
    return await invoke('read_text_file', { path });
  } catch (err) {
    console.error('[readTextFile] Error reading file:', path, err);
    return null;
  }
}

/**
 * Read a FASTA record subsequence by genomic range.
 * Range is interpreted as 1-based inclusive [start, end].
 */
export async function readFastaRange(path, start, end, recordName = null) {
  try {
    const { invoke } = window.__TAURI__.core;
    return await invoke('read_fasta_range', {
      path,
      start: Math.max(1, Number.parseInt(start, 10) || 1),
      end: Math.max(1, Number.parseInt(end, 10) || 1),
      record_name: recordName || null
    });
  } catch (err) {
    console.warn('[readFastaRange] Failed:', path, recordName, start, end, err);
    return null;
  }
}

/**
 * Get file metadata (size, modified time, etc.)
 */
export async function getFileMetadata(path) {
  try {
    const { stat } = window.__TAURI__.fs;
    return await stat(path);
  } catch (err) {
    return null;
  }
}

/**
 * Open file location in system file manager
 */
export async function openFileLocation(path) {
  try {
    const { invoke } = window.__TAURI__.core;
    await invoke('open_file_location', { path });
    return true;
  } catch (err) {
    return false;
  }
}

/**
 * Open an external URL in the system browser/mail app.
 */
export async function openExternalUrl(url) {
  try {
    const { invoke } = window.__TAURI__.core;
    await invoke('open_external_url', { url: String(url || '') });
    return true;
  } catch (err) {
    console.warn('Failed to open external URL:', url, err);
    return false;
  }
}

/**
 * Get Pathotypr process usage snapshot from backend.
 */
export async function getSystemUsage() {
  try {
    const { invoke } = window.__TAURI__.core;
    return await invoke('get_system_usage');
  } catch (err) {
    return null;
  }
}

/**
 * Get application version
 */
export async function getVersion() {
  try {
    return await tauriInvoke('get_version');
  } catch (err) {
    console.error('Failed to get version:', err);
    return null;
  }
}
