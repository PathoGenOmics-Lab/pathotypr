//! Managed state and shared types for the Tauri GUI backend.

use parking_lot::Mutex;
use pathotypr_core::CancellationToken;
use serde::{Deserialize, Serialize};
use std::process::Child;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use tauri::{AppHandle, Emitter};

/// Managed state for task tracking.
pub struct TaskState {
    pub is_running: AtomicBool,
    pub cancel_token: Mutex<CancellationToken>,
    pub match_worker_child: Mutex<Option<Child>>,
}

impl Default for TaskState {
    fn default() -> Self {
        Self {
            is_running: AtomicBool::new(false),
            cancel_token: Mutex::new(CancellationToken::new()),
            match_worker_child: Mutex::new(None),
        }
    }
}

/// Result type for GUI operations.
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct TaskResult {
    pub success: bool,
    pub message: String,
    #[serde(default)]
    pub output_path: Option<String>,
    /// Optional structured data (e.g., TrainReport) for rich GUI display.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub extra: Option<serde_json::Value>,
}

/// Progress event sent to frontend.
#[derive(Debug, Clone, Serialize)]
pub struct ProgressEvent {
    /// Current step index (0-based)
    pub step: usize,
    /// Total number of steps
    pub total_steps: usize,
    /// Step label/message
    pub message: String,
    /// Progress within current step (0.0 - 1.0), -1 for indeterminate
    pub progress: f32,
}

/// Emit a progress event to the frontend.
pub fn emit_progress(app: &AppHandle, step: usize, total_steps: usize, message: &str) {
    let _ = app.emit(
        "progress",
        ProgressEvent {
            step,
            total_steps,
            message: message.to_string(),
            progress: -1.0,
        },
    );
}

/// Emit task completion event.
pub fn emit_complete(app: &AppHandle) {
    let _ = app.emit("progress-complete", ());
}

/// Get a clone of the cancellation token for use in tasks.
pub fn get_cancel_token(state: &Arc<TaskState>) -> CancellationToken {
    state.cancel_token.lock().clone()
}

/// Reset cancellation flag and mark task as running, failing if another task is active.
pub fn try_start_task(state: &Arc<TaskState>) -> Result<(), String> {
    if state
        .is_running
        .compare_exchange(false, true, Ordering::SeqCst, Ordering::SeqCst)
        .is_err()
    {
        return Err(
            "Another task is already running. Please wait or cancel the current task.".to_string(),
        );
    }
    state.cancel_token.lock().reset();
    Ok(())
}

/// Mark task as finished.
pub fn finish_task(state: &Arc<TaskState>) {
    if let Some(mut child) = state.match_worker_child.lock().take() {
        let _ = child.kill();
        let _ = child.wait();
    }
    state.is_running.store(false, Ordering::SeqCst);
    state.cancel_token.lock().reset();
}

/// Resolve a possibly-relative output path to an absolute path.
pub fn resolve_output_path(path: &str) -> String {
    let candidate = std::path::PathBuf::from(path);
    if candidate.is_absolute() {
        return candidate.to_string_lossy().to_string();
    }
    match std::env::current_dir() {
        Ok(cwd) => cwd.join(candidate).to_string_lossy().to_string(),
        Err(_) => path.to_string(),
    }
}
