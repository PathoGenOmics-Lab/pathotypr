// ============================================================================
// Utility Functions
// ============================================================================

import { TIMING } from './config.js';

/**
 * Escape HTML special characters to prevent XSS.
 *
 * Escapes quotes as well as angle brackets so the result is safe inside
 * double- or single-quoted HTML attribute values, not just text nodes.
 */
export function escapeHtml(text) {
  return String(text ?? '')
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
}

/**
 * Format file size in human-readable format
 */
export function formatFileSize(bytes) {
  if (bytes === 0) return '0 B';
  const k = 1024;
  const sizes = ['B', 'KB', 'MB', 'GB'];
  const i = Math.floor(Math.log(bytes) / Math.log(k));
  return parseFloat((bytes / Math.pow(k, i)).toFixed(1)) + ' ' + sizes[i];
}

/**
 * Validate file extension against allowed extensions
 */
export function validateFileExtension(filePath, validExtensions) {
  if (!validExtensions || validExtensions.length === 0) return true;
  const fileName = filePath.toLowerCase();
  return validExtensions.some(ext => fileName.endsWith(ext.toLowerCase().trim()));
}

/**
 * Create a clickable path element string
 */
export function makeClickablePath(path) {
  const safePath = String(path ?? '');
  // Return a token that logMessage can render safely as a clickable path.
  return `[[PATHOTYPR_PATH:${encodeURIComponent(safePath)}]]`;
}

/**
 * Split message text into text and clickable-path token parts.
 */
export function splitClickablePathTokens(message) {
  const msg = String(message ?? '');
  const tokenRegex = /\[\[PATHOTYPR_PATH:([^\]]+)]]/g;
  const parts = [];
  let lastIndex = 0;
  let match;

  while ((match = tokenRegex.exec(msg)) !== null) {
    if (match.index > lastIndex) {
      parts.push({ type: 'text', value: msg.slice(lastIndex, match.index) });
    }
    let decodedPath = '';
    try {
      decodedPath = decodeURIComponent(match[1]);
    } catch {
      decodedPath = match[1];
    }
    parts.push({ type: 'path', value: decodedPath });
    lastIndex = tokenRegex.lastIndex;
  }

  if (lastIndex < msg.length) {
    parts.push({ type: 'text', value: msg.slice(lastIndex) });
  }

  return parts;
}

/**
 * Get file name from path
 */
export function getFileName(path) {
  return path.split(/[/\\]/).pop();
}

/**
 * Show extension error animation on dropzone
 */
export function showExtensionError(dropzone, validExtensions) {
  dropzone.classList.add('extension-error');
  setTimeout(() => {
    dropzone.classList.remove('extension-error');
  }, TIMING.ANIMATION_FEEDBACK);
}

/**
 * Highlight a form field as missing
 */
export function highlightMissingField(fieldId) {
  const el = document.getElementById(fieldId);
  if (!el) return;

  const dropzone = el.closest('.form-group')?.querySelector('.dropzone');
  const inputWrapper = el.closest('.file-input-wrapper');
  const target = dropzone || inputWrapper || el;

  target.classList.add('field-missing');
  setTimeout(() => target.classList.remove('field-missing'), TIMING.FIELD_MISSING);
}

function getFieldErrorTarget(fieldId) {
  const el = document.getElementById(fieldId);
  if (!el) return null;
  const dropzone = el.closest('.form-group')?.querySelector('.dropzone');
  const inputWrapper = el.closest('.file-input-wrapper');
  return dropzone || inputWrapper || el;
}

function getFieldErrorHost(fieldId) {
  const el = document.getElementById(fieldId);
  if (!el) return null;
  return el.closest('.form-group') || el.parentElement;
}

export function clearFieldError(fieldId) {
  const host = getFieldErrorHost(fieldId);
  const target = getFieldErrorTarget(fieldId);
  if (!host) return;

  host.querySelectorAll(`.field-error-msg[data-for="${fieldId}"]`).forEach(node => node.remove());
  if (target) target.classList.remove('field-invalid');
}

export function showFieldError(fieldId, message) {
  const host = getFieldErrorHost(fieldId);
  const target = getFieldErrorTarget(fieldId);
  if (!host) return;

  clearFieldError(fieldId);
  if (target) target.classList.add('field-invalid');

  const error = document.createElement('div');
  error.className = 'field-error-msg';
  error.dataset.for = fieldId;
  error.textContent = message;
  host.appendChild(error);
}

export function clearFormFieldErrors(formOrElement) {
  if (!formOrElement) return;
  const scope = formOrElement.matches?.('form')
    ? formOrElement
    : formOrElement.closest?.('form');
  if (!scope) return;

  scope.querySelectorAll('.field-error-msg').forEach(node => node.remove());
  scope.querySelectorAll('.field-invalid').forEach(node => node.classList.remove('field-invalid'));
}

/**
 * Format backend error messages into user-friendly text
 */
export function formatError(err) {
  const msg = typeof err === 'string' ? err : (err?.message || String(err));

  // Map common backend error patterns to friendly messages
  const patterns = [
    [/Failed to open.*?:\s*No such file/i, 'File not found. Check that the file path is correct.'],
    [/Failed to open.*?:\s*Permission denied/i, 'Cannot access file. Check file permissions.'],
    [/Invalid position in TSV/i, 'Markers file has invalid format. The first column must be numeric positions.'],
    [/No FASTQ files found/i, 'No FASTQ files found in the provided list. Check the input file.'],
    [/not enough data/i, 'Not enough data to process. Check that input files are not empty.'],
    [/connection.*?reset|connection.*?refused|Download failed/i, 'Network error. Check your internet connection and try again.'],
    [/HTTP error: 4\d\d/i, 'Server rejected the request. The URL may be invalid or require authentication.'],
    [/HTTP error: 5\d\d/i, 'Server error. Try again later.'],
    [/cancelled|canceled/i, 'Task was cancelled.'],
  ];

  for (const [pattern, friendly] of patterns) {
    if (pattern.test(msg)) return friendly;
  }

  // Trim long Rust error chains: keep the last meaningful part
  if (msg.includes(': ') && msg.length > 120) {
    const parts = msg.split(': ');
    return parts[parts.length - 1];
  }

  return msg;
}

/**
 * Validate form fields
 */
export function validateForm(fields, formName, logMessageFn, openConsoleFn = null, options = {}) {
  const { openConsole = false } = options;
  const missing = [];

  for (const [fieldId, label] of Object.entries(fields)) {
    clearFieldError(fieldId);
    const el = document.getElementById(fieldId);
    if (!el || !el.value || el.value.trim() === '') {
      missing.push(label);
      highlightMissingField(fieldId);
      showFieldError(fieldId, `${label} is required.`);
    }
  }

  if (missing.length > 0) {
    logMessageFn(`Cannot run ${formName}: Missing required fields`, 'error');
    missing.forEach(field => {
      logMessageFn(`  - ${field}`, 'error');
    });
    if (openConsole && typeof openConsoleFn === 'function') {
      openConsoleFn();
    }
    return false;
  }
  return true;
}
