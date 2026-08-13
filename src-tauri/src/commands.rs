//! Tauri command handlers for all GUI-invokable operations.

use pathotypr_core::{
    classify, classify_split_fastq, defaults, predict, r#match, train, ClassifyArgs,
    MatchArgs, PredictArgs, SplitFastqArgs, TrainArgs,
};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::sync::atomic::Ordering;
use std::sync::Arc;
use std::time::{Duration, SystemTime, UNIX_EPOCH};
use tauri::{AppHandle, Emitter, Manager, State};

use crate::state::*;
use crate::util;

// ============================================================================
// TRAIN
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct TrainParams {
    pub input: String,
    pub output: String,
    pub kmer_size: Option<usize>,
    pub test_split: Option<f64>,
    pub threads: Option<usize>,
    /// Number of cross-validation folds (e.g. 5 or 10). None = single split.
    pub cv_folds: Option<usize>,
    /// Maximum tree depth (default: 20).
    pub max_depth: Option<usize>,
    /// Minimum samples per leaf node (default: 5).
    pub min_samples_leaf: Option<usize>,
}

#[tauri::command]
pub async fn run_train(
    app: AppHandle,
    params: TrainParams,
    state: State<'_, Arc<TaskState>>,
) -> Result<TaskResult, String> {
    try_start_task(&state)?;
    let cancel_token = get_cancel_token(&state);

    emit_progress(&app, 0, 3, "Loading training data...");

    let args = TrainArgs {
        input: params.input,
        output: params.output.clone(),
        kmer_size: params.kmer_size.unwrap_or(21),
        test_split: params.test_split.unwrap_or(0.2),
        threads: params.threads,
        cv_folds: params.cv_folds,
        max_depth: params.max_depth.unwrap_or(20),
        min_samples_leaf: params.min_samples_leaf.unwrap_or(5),
        cancel_token: Some(cancel_token),
    };

    let cv_active = args.cv_folds.is_some();
    emit_progress(&app, 1, 3, if cv_active {
        "Running cross-validation & training model..."
    } else {
        "Training model..."
    });

    let train_result = tauri::async_runtime::spawn_blocking(move || train::run(args)).await;
    let result = match train_result {
        Err(e) => {
            emit_complete(&app);
            finish_task(&state);
            return Err(format!("Task failed: {}", e));
        }
        Ok(Ok(report)) => {
            emit_progress(&app, 2, 3, "Saving model...");
            let extra = serde_json::to_value(&report).ok();
            let accuracy_msg = if let (Some(mean), Some(std)) = (report.cv_mean_accuracy_pct, report.cv_std_accuracy_pct) {
                format!("CV Accuracy: {:.2}% ± {:.2}%", mean, std)
            } else {
                format!("Accuracy: {:.2}%", report.accuracy_pct)
            };
            let oob_msg = report.oob_accuracy_pct
                .map(|oob| format!(" | OOB: {:.2}%", oob))
                .unwrap_or_default();
            TaskResult {
                success: true,
                message: format!("Model trained successfully! {}{}", accuracy_msg, oob_msg),
                output_path: Some(resolve_output_path(&params.output)),
                extra,
            }
        }
        Ok(Err(pathotypr_core::AppError::Cancelled)) => TaskResult {
            success: false,
            message: "Training cancelled by user.".to_string(),
            output_path: None,
            extra: None,
        },
        Ok(Err(e)) => TaskResult {
            success: false,
            message: format!("Training failed: {}", e),
            output_path: None,
            extra: None,
        },
    };

    emit_complete(&app);
    finish_task(&state);
    Ok(result)
}

// ============================================================================
// PREDICT
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct PredictParams {
    pub input: String,
    pub model: String,
    pub output: String,
    pub threads: Option<usize>,
    pub excel: Option<bool>,
}

#[tauri::command]
pub async fn run_predict(
    app: AppHandle,
    params: PredictParams,
    state: State<'_, Arc<TaskState>>,
) -> Result<TaskResult, String> {
    try_start_task(&state)?;
    let cancel_token = get_cancel_token(&state);

    emit_progress(&app, 0, 3, "Loading model...");

    let args = PredictArgs {
        input: params.input,
        model: params.model,
        output: params.output.clone(),
        threads: params.threads,
        excel: params.excel.unwrap_or(true),
        cancel_token: Some(cancel_token),
    };

    emit_progress(&app, 1, 3, "Processing sequences...");

    let predict_result = tauri::async_runtime::spawn_blocking(move || predict::run(args)).await;
    let result = match predict_result {
        Err(e) => {
            emit_complete(&app);
            finish_task(&state);
            return Err(format!("Task failed: {}", e));
        }
        Ok(Ok(_)) => {
            emit_progress(&app, 2, 3, "Classifying...");
            TaskResult {
                success: true,
                message: "Prediction completed successfully!".to_string(),
                output_path: Some(resolve_output_path(&params.output)),
                extra: None,
            }
        }
        Ok(Err(pathotypr_core::AppError::Cancelled)) => TaskResult {
            success: false,
            message: "Prediction cancelled by user.".to_string(),
            output_path: None,
            extra: None,
        },
        Ok(Err(e)) => TaskResult {
            success: false,
            message: format!("Prediction failed: {}", e),
            output_path: None,
            extra: None,
        },
    };

    emit_complete(&app);
    finish_task(&state);
    Ok(result)
}

// ============================================================================
// CLASSIFY
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct ClassifyParams {
    /// Single marker file (legacy) or first marker file
    pub markers: Option<String>,
    /// Multiple marker files (new)
    pub marker_files: Option<Vec<String>>,
    pub reference: String,
    pub input_list: Option<String>,
    pub input: Option<String>,
    /// Multiple FASTA files (GUI batch mode)
    pub input_files: Option<Vec<String>>,
    pub gff: Option<String>,
    /// Multiple GFF files (GUI batch mode, matched by filename to FASTA files)
    pub gff_files: Option<Vec<String>>,
    pub output_prefix: String,
    pub kmer_size: Option<usize>,
    pub threads: Option<usize>,
    pub nested_classification: Option<bool>,
    pub output_masked_fasta: Option<bool>,
    pub min_flank_bases: Option<usize>,
    pub excel: Option<bool>,
}

#[tauri::command]
pub async fn run_classify(
    app: AppHandle,
    params: ClassifyParams,
    state: State<'_, Arc<TaskState>>,
) -> Result<TaskResult, String> {
    try_start_task(&state)?;
    let cancel_token = get_cancel_token(&state);

    // Resolve marker files: support both single and multiple
    let marker_files: Vec<String> = if let Some(files) = &params.marker_files {
        files.clone()
    } else if let Some(single) = &params.markers {
        vec![single.clone()]
    } else {
        emit_complete(&app);
        finish_task(&state);
        return Err("No marker files provided.".to_string());
    };

    let total_steps = marker_files.len() * 2 + 1;
    let mut output_paths: Vec<serde_json::Value> = Vec::new();
    let mut first_output: Option<String> = None;

    for (idx, marker_file) in marker_files.iter().enumerate() {
        let marker_name = std::path::Path::new(marker_file)
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        let step_base = idx * 2;

        emit_progress(
            &app,
            step_base,
            total_steps,
            &format!("Loading markers: {}...", marker_name),
        );

        // Generate unique output prefix per marker file. Include the index so
        // two marker files with the same basename (from different folders) do
        // not collide and silently overwrite each other's output.
        let ofile = if marker_files.len() > 1 {
            format!("{}_{}_{}", params.output_prefix, idx + 1, marker_name)
        } else {
            params.output_prefix.clone()
        };

        let args = ClassifyArgs {
            tsv_pos: marker_file.clone(),
            ref_fasta: params.reference.clone(),
            tsv_genomes: params.input_list.clone(),
            fasta_genomes: params.input.clone(),
            fasta_files: params.input_files.clone(),
            gff_file: params.gff.clone(),
            gff_files: params.gff_files.clone(),
            ofile: ofile.clone(),
            kmer_size: params.kmer_size.unwrap_or(31),
            num_cpu: params.threads,
            nested_classification: params.nested_classification.unwrap_or(false),
            output_masked_fasta: params.output_masked_fasta.unwrap_or(false),
            min_flank_bases: params.min_flank_bases.unwrap_or(10),
            excel: params.excel.unwrap_or(true),
            cancel_token: Some(cancel_token.clone()),
        };

        emit_progress(
            &app,
            step_base + 1,
            total_steps,
            &format!("Classifying with {}...", marker_name),
        );

        let classify_result =
            tauri::async_runtime::spawn_blocking(move || classify::run(args)).await;

        match classify_result {
            Err(e) => {
                emit_complete(&app);
                finish_task(&state);
                return Err(format!("Task failed: {}", e));
            }
            Ok(Err(pathotypr_core::AppError::Cancelled)) => {
                emit_complete(&app);
                finish_task(&state);
                return Ok(TaskResult {
                    success: false,
                    message: "Task cancelled by user.".to_string(),
                    output_path: None,
                    extra: None,
                });
            }
            Ok(Err(e)) => {
                emit_complete(&app);
                finish_task(&state);
                return Ok(TaskResult {
                    success: false,
                    message: format!("Classification failed ({}): {}", marker_name, e),
                    output_path: None,
                    extra: None,
                });
            }
            Ok(Ok(_)) => {
                let detailed_path = if ofile.ends_with(".tsv") {
                    ofile.clone()
                } else {
                    format!("{}.tsv", ofile)
                };
                let resolved = resolve_output_path(&detailed_path);
                if first_output.is_none() {
                    first_output = Some(resolved.clone());
                }
                output_paths.push(serde_json::json!({
                    "name": marker_name,
                    "path": resolved,
                }));
            }
        }
    }

    emit_progress(&app, total_steps, total_steps, "Done!");

    let result = TaskResult {
        success: true,
        message: format!(
            "Classification completed successfully! ({} marker set{})",
            marker_files.len(),
            if marker_files.len() > 1 { "s" } else { "" }
        ),
        output_path: first_output,
        extra: if output_paths.len() > 1 {
            Some(serde_json::json!({ "output_sets": output_paths }))
        } else {
            None
        },
    };

    emit_complete(&app);
    finish_task(&state);
    Ok(result)
}

// ============================================================================
// SPLIT-FASTQ
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct SplitFastqParams {
    pub input: Option<Vec<String>>,
    pub input_list: Option<String>,
    pub paired: Option<bool>,
    pub reference: String,
    /// Single marker file (legacy)
    pub markers: Option<String>,
    /// Multiple marker files (new)
    pub marker_files: Option<Vec<String>>,
    pub threads: Option<usize>,
    pub output_prefix: Option<String>,
    pub min_depth: Option<u32>,
    pub min_alt_percent: Option<u32>,
    pub nested_classification: Option<bool>,
    pub kmer_size: Option<usize>,
    pub excel: Option<bool>,
}

#[tauri::command]
pub async fn run_split_fastq(
    app: AppHandle,
    params: SplitFastqParams,
    state: State<'_, Arc<TaskState>>,
) -> Result<TaskResult, String> {
    try_start_task(&state)?;
    let cancel_token = get_cancel_token(&state);

    // Resolve marker files: support both single and multiple
    let marker_files: Vec<String> = if let Some(files) = &params.marker_files {
        files.clone()
    } else if let Some(single) = &params.markers {
        vec![single.clone()]
    } else {
        emit_complete(&app);
        finish_task(&state);
        return Err("No marker files provided.".to_string());
    };

    let output_prefix = params.output_prefix.clone().unwrap_or("split".to_string());
    let total_steps = marker_files.len() * 2 + 1;
    let mut output_paths: Vec<serde_json::Value> = Vec::new();
    let mut first_output: Option<String> = None;

    let input_files = params.input.unwrap_or_default();
    let input_list = params.input_list.clone();

    for (idx, marker_file) in marker_files.iter().enumerate() {
        let marker_name = std::path::Path::new(marker_file)
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        let step_base = idx * 2;

        emit_progress(
            &app,
            step_base,
            total_steps,
            &format!("Indexing reference for {}...", marker_name),
        );

        // Generate unique output prefix per marker file. Include the index so
        // two marker files with the same basename (from different folders) do
        // not collide and silently overwrite each other's output.
        let oprefix = if marker_files.len() > 1 {
            format!("{}_{}_{}", output_prefix, idx + 1, marker_name)
        } else {
            output_prefix.clone()
        };

        let args = SplitFastqArgs {
            input: input_files.clone(),
            input_list: input_list.clone(),
            paired: params.paired.unwrap_or(false),
            no_auto_paired: false,
            reference: params.reference.clone(),
            markers: marker_file.clone(),
            threads: params.threads,
            output_prefix: oprefix.clone(),
            min_depth: params.min_depth.unwrap_or(10),
            min_alt_percent: params.min_alt_percent.unwrap_or(95),
            nested_classification: params.nested_classification.unwrap_or(false),
            kmer_size: params.kmer_size.unwrap_or(pathotypr_core::split_kmer::DEFAULT_MARKER_KMER_LEN),
            excel: params.excel.unwrap_or(true),
            cancel_token: Some(cancel_token.clone()),
        };

        emit_progress(
            &app,
            step_base + 1,
            total_steps,
            &format!("Processing reads with {}...", marker_name),
        );

        let split_result =
            tauri::async_runtime::spawn_blocking(move || classify_split_fastq::run(args)).await;

        match split_result {
            Err(e) => {
                emit_complete(&app);
                finish_task(&state);
                return Err(format!("Task failed: {}", e));
            }
            Ok(Err(pathotypr_core::AppError::Cancelled)) => {
                emit_complete(&app);
                finish_task(&state);
                return Ok(TaskResult {
                    success: false,
                    message: "Task cancelled by user.".to_string(),
                    output_path: None,
                    extra: None,
                });
            }
            Ok(Err(e)) => {
                emit_complete(&app);
                finish_task(&state);
                return Ok(TaskResult {
                    success: false,
                    message: format!("Split-FASTQ failed ({}): {}", marker_name, e),
                    output_path: None,
                    extra: None,
                });
            }
            Ok(Ok(_)) => {
                let summary_path = format!("{}_summary.tsv", oprefix);
                let resolved = resolve_output_path(&summary_path);
                if first_output.is_none() {
                    first_output = Some(resolved.clone());
                }
                output_paths.push(serde_json::json!({
                    "name": marker_name,
                    "path": resolved,
                }));
            }
        }
    }

    emit_progress(&app, total_steps, total_steps, "Done!");

    let result = TaskResult {
        success: true,
        message: format!(
            "Split-FASTQ analysis completed successfully! ({} marker set{})",
            marker_files.len(),
            if marker_files.len() > 1 { "s" } else { "" }
        ),
        output_path: first_output,
        extra: if output_paths.len() > 1 {
            Some(serde_json::json!({ "output_sets": output_paths }))
        } else {
            None
        },
    };

    emit_complete(&app);
    finish_task(&state);
    Ok(result)
}

// ============================================================================
// MATCH (with worker subprocess)
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct MatchParams {
    pub fastqs: Option<Vec<String>>,
    pub input_list: Option<String>,
    /// Path to multi-FASTA reference file
    pub references: Option<String>,
    pub output: Option<String>,
    pub kmer_size: Option<u8>,
    pub threads: Option<usize>,
    pub early_stop_confidence: Option<f64>,
    pub early_stop_min_kmers: Option<u64>,
    pub strict_percentages: Option<bool>,
    pub min_kmer_count: Option<u32>,
    pub excel: Option<bool>,
}

pub const MATCH_WORKER_FLAG: &str = "--pathotypr-match-worker";

#[derive(Debug, Serialize, Deserialize)]
pub struct MatchWorkerRequest {
    pub fastqs: Vec<String>,
    pub input_list: Option<String>,
    pub paired: bool,
    pub no_auto_paired: bool,
    pub references: Option<String>,
    pub output: Option<String>,
    pub kmer_size: u8,
    pub threads: Option<usize>,
    pub early_stop_confidence: f64,
    pub early_stop_min_kmers: u64,
    pub strict_percentages: bool,
    pub min_kmer_count: u32,
    pub excel: bool,
}

impl MatchWorkerRequest {
    fn into_match_args(self) -> MatchArgs {
        MatchArgs {
            fastqs: self.fastqs,
            input_list: self.input_list,
            paired: self.paired,
            no_auto_paired: self.no_auto_paired,
            references: self.references,
            output: self.output,
            kmer_size: self.kmer_size,
            threads: self.threads,
            early_stop_confidence: self.early_stop_confidence,
            early_stop_min_kmers: self.early_stop_min_kmers,
            strict_percentages: self.strict_percentages,
            min_kmer_count: self.min_kmer_count,
            excel: self.excel,
            cancel_token: None,
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MatchWorkerResponse {
    pub success: bool,
    pub cancelled: bool,
    pub message: String,
}

pub fn run_match_worker_mode(request_path: &str, response_path: &str) -> Result<(), String> {
    let payload_text = std::fs::read_to_string(request_path)
        .map_err(|e| format!("Failed to read worker request: {}", e))?;
    let payload: MatchWorkerRequest = serde_json::from_str(&payload_text)
        .map_err(|e| format!("Failed to parse worker request JSON: {}", e))?;
    let args = payload.into_match_args();

    let response = match r#match::run(args) {
        Ok(()) => MatchWorkerResponse {
            success: true,
            cancelled: false,
            message: "Reference matching completed successfully.".to_string(),
        },
        Err(pathotypr_core::AppError::Cancelled) => MatchWorkerResponse {
            success: false,
            cancelled: true,
            message: "Task cancelled by user.".to_string(),
        },
        Err(e) => MatchWorkerResponse {
            success: false,
            cancelled: false,
            message: format!("Reference matching failed: {}", e),
        },
    };

    let response_json = serde_json::to_string(&response)
        .map_err(|e| format!("Failed to serialize worker response: {}", e))?;
    std::fs::write(response_path, response_json)
        .map_err(|e| format!("Failed to write worker response: {}", e))?;
    Ok(())
}

pub fn maybe_run_match_worker_from_cli() -> Option<i32> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() >= 4 && args.get(1).map(String::as_str) == Some(MATCH_WORKER_FLAG) {
        let code = match run_match_worker_mode(&args[2], &args[3]) {
            Ok(()) => 0,
            Err(e) => {
                eprintln!("{}", e);
                1
            }
        };
        return Some(code);
    }
    None
}

fn build_match_worker_temp_paths() -> (PathBuf, PathBuf) {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_else(|_| Duration::from_secs(0))
        .as_nanos();
    let pid = std::process::id();
    let base = std::env::temp_dir();
    let req = base.join(format!("pathotypr_match_req_{}_{}.json", pid, ts));
    let resp = base.join(format!("pathotypr_match_resp_{}_{}.json", pid, ts));
    (req, resp)
}

fn run_match_via_worker(
    state: &Arc<TaskState>,
    request: MatchWorkerRequest,
) -> Result<MatchWorkerResponse, String> {
    let current_exe = std::env::current_exe()
        .map_err(|e| format!("Failed to resolve current executable path: {}", e))?;
    let (request_path, response_path) = build_match_worker_temp_paths();
    let request_json = serde_json::to_string(&request)
        .map_err(|e| format!("Failed to serialize match worker request: {}", e))?;
    std::fs::write(&request_path, request_json)
        .map_err(|e| format!("Failed to write match worker request: {}", e))?;

    let spawn_result = Command::new(current_exe)
        .arg(MATCH_WORKER_FLAG)
        .arg(&request_path)
        .arg(&response_path)
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn();

    let child = match spawn_result {
        Ok(child) => child,
        Err(e) => {
            let _ = std::fs::remove_file(&request_path);
            return Err(format!(
                "Failed to start isolated match worker process: {}",
                e
            ));
        }
    };

    *state.match_worker_child.lock() = Some(child);

    let status = loop {
        let maybe_status = {
            let mut guard = state.match_worker_child.lock();
            match guard.as_mut() {
                Some(child) => child
                    .try_wait()
                    .map_err(|e| format!("Failed while waiting for match worker: {}", e))?,
                None => {
                    return Err("Match worker handle disappeared unexpectedly.".to_string());
                }
            }
        };

        if let Some(status) = maybe_status {
            break status;
        }
        std::thread::sleep(Duration::from_millis(80));
    };

    let _ = state.match_worker_child.lock().take();
    let cancelled_by_user = state.cancel_token.lock().is_cancelled();

    let mut response = if response_path.exists() {
        let text = std::fs::read_to_string(&response_path)
            .map_err(|e| format!("Failed to read match worker response: {}", e))?;
        serde_json::from_str::<MatchWorkerResponse>(&text)
            .map_err(|e| format!("Failed to parse match worker response: {}", e))?
    } else if cancelled_by_user {
        MatchWorkerResponse {
            success: false,
            cancelled: true,
            message: "Task cancelled by user.".to_string(),
        }
    } else if status.success() {
        MatchWorkerResponse {
            success: true,
            cancelled: false,
            message: "Reference matching completed successfully.".to_string(),
        }
    } else {
        MatchWorkerResponse {
            success: false,
            cancelled: false,
            message: format!("Match worker exited with status {}", status),
        }
    };

    if cancelled_by_user && !response.success {
        response.cancelled = true;
        response.message = "Task cancelled by user.".to_string();
    }

    let _ = std::fs::remove_file(&request_path);
    let _ = std::fs::remove_file(&response_path);

    Ok(response)
}

#[tauri::command]
pub async fn run_match(
    app: AppHandle,
    params: MatchParams,
    state: State<'_, Arc<TaskState>>,
) -> Result<TaskResult, String> {
    if params.references.is_none() {
        return Ok(TaskResult {
            success: false,
            message: "References FASTA must be provided.".to_string(),
            output_path: None,
            extra: None,
        });
    }

    try_start_task(&state)?;

    emit_progress(
        &app,
        0,
        5,
        "Preparing reference matching...",
    );

    let request = MatchWorkerRequest {
        fastqs: params.fastqs.unwrap_or_default(),
        input_list: params.input_list,
        paired: false,
        no_auto_paired: false,
        references: params.references,
        output: params.output.clone(),
        kmer_size: params.kmer_size.unwrap_or(31),
        threads: params.threads,
        early_stop_confidence: params.early_stop_confidence.unwrap_or(0.0),
        early_stop_min_kmers: params.early_stop_min_kmers.unwrap_or(1_000_000),
        strict_percentages: params.strict_percentages.unwrap_or(true),
        min_kmer_count: params.min_kmer_count.unwrap_or(2),
        excel: params.excel.unwrap_or(true),
    };

    let output_path = request.output.as_deref().map(resolve_output_path);
    let app_clone = app.clone();
    let state_arc = state.inner().clone();

    let join_result = tauri::async_runtime::spawn_blocking(move || {
        let _ = app_clone.emit("log", "Loading reference k-mers...");
        run_match_via_worker(&state_arc, request)
    })
    .await;

    let result = match join_result {
        Ok(Ok(worker_response)) => {
            if worker_response.success {
                let _ = app.emit("log", "Matching complete!");
                TaskResult {
                    success: true,
                    message: "Reference matching completed successfully!".to_string(),
                    output_path,
                    extra: None,
                }
            } else if worker_response.cancelled {
                let _ = app.emit("log", "Task cancelled by user.");
                TaskResult {
                    success: false,
                    message: "Task cancelled by user.".to_string(),
                    output_path: None,
                    extra: None,
                }
            } else {
                let _ = app.emit("log", format!("Error: {}", worker_response.message));
                TaskResult {
                    success: false,
                    message: worker_response.message,
                    output_path: None,
                    extra: None,
                }
            }
        }
        Ok(Err(e)) => {
            let _ = app.emit("log", format!("Error: {}", e));
            TaskResult {
                success: false,
                message: format!("Reference matching failed: {}", e),
                output_path: None,
                extra: None,
            }
        }
        Err(e) => {
            emit_complete(&app);
            finish_task(&state);
            return Err(format!("Task failed: {}", e));
        }
    };

    emit_complete(&app);
    finish_task(&state);
    Ok(result)
}

// ============================================================================
// TASK CONTROL
// ============================================================================

#[tauri::command]
pub fn cancel_task(state: State<'_, Arc<TaskState>>) -> bool {
    if state.is_running.load(Ordering::SeqCst) {
        state.cancel_token.lock().cancel();
        if let Some(child) = state.match_worker_child.lock().as_mut() {
            let _ = child.kill();
        }
        true
    } else {
        false
    }
}

#[tauri::command]
pub fn is_task_running(state: State<'_, Arc<TaskState>>) -> bool {
    state.is_running.load(Ordering::SeqCst)
}

// ============================================================================
// DOWNLOAD FILE
// ============================================================================

#[derive(Debug, Deserialize)]
pub struct DownloadParams {
    pub url: String,
    pub filename: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct DownloadResult {
    pub success: bool,
    pub path: Option<String>,
    pub error: Option<String>,
}

#[tauri::command]
pub async fn download_file(
    app: AppHandle,
    params: DownloadParams,
) -> Result<DownloadResult, String> {
    use futures_util::StreamExt;
    use std::io::Write;
    use std::path::Path;
    use std::time::Duration;

    // --- Validate filename: prevent path traversal ---
    let safe_filename = Path::new(&params.filename)
        .file_name()
        .and_then(|n| n.to_str())
        .ok_or_else(|| "Invalid filename".to_string())?;
    if safe_filename.contains("..") {
        return Err("Invalid filename: path traversal not allowed".to_string());
    }

    // --- Validate URL: only allow HTTPS and block internal/reserved targets ---
    let parsed_url = reqwest::Url::parse(&params.url).map_err(|_| "Invalid URL".to_string())?;
    util::validate_download_url(&parsed_url)?;

    // Save to app data dir / demo_data /
    let app_data_dir = app
        .path()
        .app_data_dir()
        .map_err(|e| format!("Failed to get app data dir: {}", e))?;
    let demo_dir = app_data_dir.join("demo_data");
    std::fs::create_dir_all(&demo_dir)
        .map_err(|e| format!("Failed to create demo_data dir: {}", e))?;

    let dest_path = demo_dir.join(safe_filename);

    // If file already exists, return its path directly
    if dest_path.exists() {
        return Ok(DownloadResult {
            success: true,
            path: Some(dest_path.to_string_lossy().to_string()),
            error: None,
        });
    }

    let _ = app.emit("log", format!("Downloading {}...", safe_filename));

    // Build HTTP client with User-Agent, timeout, and redirect policy.
    // The custom DNS resolver validates the resolved IP for every connection
    // (initial and redirects), so a low-TTL DNS record cannot rebind to a
    // local/reserved address between validation and connect.
    let client = reqwest::Client::builder()
        .user_agent(format!("Pathotypr/{}", env!("CARGO_PKG_VERSION")))
        .dns_resolver(Arc::new(util::SsrfGuardResolver))
        .redirect(reqwest::redirect::Policy::custom(|attempt| {
            if attempt.previous().len() >= 10 {
                return attempt.stop();
            }
            if let Err(e) = util::validate_download_url(attempt.url()) {
                return attempt.error(e);
            }
            attempt.follow()
        }))
        .connect_timeout(Duration::from_secs(30))
        .timeout(Duration::from_secs(600))
        .build()
        .map_err(|e| format!("Failed to create HTTP client: {}", e))?;

    let response = client
        .get(&params.url)
        .send()
        .await
        .map_err(|e| format!("Download failed: {}", e))?;

    if !response.status().is_success() {
        return Ok(DownloadResult {
            success: false,
            path: None,
            error: Some(format!("HTTP error: {}", response.status())),
        });
    }

    // --- Check file size limit (max 2 GB) ---
    const MAX_DOWNLOAD_SIZE: u64 = 2_000_000_000;
    let total_size = response.content_length();
    if let Some(total) = total_size {
        if total > MAX_DOWNLOAD_SIZE {
            return Ok(DownloadResult {
                success: false,
                path: None,
                error: Some(format!(
                    "File too large: {:.0} MB (max {} MB)",
                    total as f64 / 1_000_000.0,
                    MAX_DOWNLOAD_SIZE / 1_000_000
                )),
            });
        }
    }

    let mut stream = response.bytes_stream();
    let tmp_path = demo_dir.join(format!("{}.tmp", safe_filename));
    let mut file =
        std::fs::File::create(&tmp_path).map_err(|e| format!("Failed to create file: {}", e))?;

    let mut downloaded: u64 = 0;

    while let Some(chunk) = stream.next().await {
        let chunk = chunk.map_err(|e| {
            let _ = std::fs::remove_file(&tmp_path);
            format!("Download error: {}", e)
        })?;
        file.write_all(&chunk).map_err(|e| {
            let _ = std::fs::remove_file(&tmp_path);
            format!("Write error: {}", e)
        })?;
        downloaded += chunk.len() as u64;

        // Enforce size limit even if Content-Length was missing/wrong
        if downloaded > MAX_DOWNLOAD_SIZE {
            let _ = std::fs::remove_file(&tmp_path);
            return Ok(DownloadResult {
                success: false,
                path: None,
                error: Some(format!(
                    "Download exceeded max size ({} MB)",
                    MAX_DOWNLOAD_SIZE / 1_000_000
                )),
            });
        }

        if let Some(total) = total_size {
            let pct = (downloaded as f64 / total as f64 * 100.0) as u32;
            let _ = app.emit("log", format!("Downloading {}... {}%", safe_filename, pct));
        }
    }

    // Rename temp file to final path
    std::fs::rename(&tmp_path, &dest_path).map_err(|e| {
        let _ = std::fs::remove_file(&tmp_path);
        format!("Failed to finalize download: {}", e)
    })?;

    let _ = app.emit("log", format!("Downloaded: {}", safe_filename));

    Ok(DownloadResult {
        success: true,
        path: Some(dest_path.to_string_lossy().to_string()),
        error: None,
    })
}

// ============================================================================
// FILE UTILITY COMMANDS
// ============================================================================

#[tauri::command]
pub async fn read_text_file(app: AppHandle, path: String) -> Result<String, String> {
    const MAX_TEXT_FILE_SIZE: u64 = 16 * 1024 * 1024;
    let canonical = util::resolve_allowed_file_path(&app, &path)?;

    let metadata = std::fs::metadata(&canonical)
        .map_err(|e| format!("Failed to get metadata '{}': {}", canonical.display(), e))?;
    if metadata.len() > MAX_TEXT_FILE_SIZE {
        return Err(format!(
            "File is too large to preview (max {} MB): {}",
            MAX_TEXT_FILE_SIZE / (1024 * 1024),
            canonical.display()
        ));
    }

    std::fs::read_to_string(&canonical)
        .map_err(|e| format!("Failed to read file '{}': {}", canonical.display(), e))
}

#[derive(Debug, Clone, Serialize)]
pub struct FastaRangeResult {
    pub path: String,
    pub record_name: String,
    pub start: usize,
    pub end: usize,
    pub total_length: usize,
    pub sequence: String,
}

#[tauri::command]
pub async fn read_fasta_range(
    app: AppHandle,
    path: String,
    start: usize,
    end: usize,
    record_name: Option<String>,
) -> Result<FastaRangeResult, String> {
    use std::io::{BufRead, BufReader};
    const MAX_FASTA_SLICE_BP: usize = 50_000;

    let canonical = util::resolve_allowed_file_path(&app, &path)?;
    let file = std::fs::File::open(&canonical)
        .map_err(|e| format!("Failed to open file '{}': {}", canonical.display(), e))?;
    let reader = BufReader::new(file);

    let target_record = record_name
        .as_ref()
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty());
    let mut selected_record = String::new();
    let mut collecting = false;
    let mut sequence = String::new();

    for line_result in reader.lines() {
        let line = line_result
            .map_err(|e| format!("Failed to read file '{}': {}", canonical.display(), e))?;
        if let Some(header) = line.strip_prefix('>') {
            if collecting && !sequence.is_empty() {
                break;
            }

            let full_header = header.trim();
            let record_id = full_header.split_whitespace().next().unwrap_or("").trim();
            let matches_target = match &target_record {
                Some(target) => {
                    record_id.eq_ignore_ascii_case(target)
                        || full_header.eq_ignore_ascii_case(target)
                }
                None => selected_record.is_empty(),
            };

            if matches_target {
                collecting = true;
                selected_record = if record_id.is_empty() {
                    full_header.to_string()
                } else {
                    record_id.to_string()
                };
            } else {
                collecting = false;
            }
            continue;
        }

        if collecting {
            sequence.extend(
                line.chars()
                    .filter(|ch| !ch.is_whitespace())
                    .map(|ch| ch.to_ascii_uppercase()),
            );
        }
    }

    if selected_record.is_empty() || sequence.is_empty() {
        if let Some(target) = target_record {
            return Err(format!(
                "Could not find FASTA record '{}' in '{}'",
                target,
                canonical.display()
            ));
        }
        return Err(format!(
            "No FASTA records found in '{}'",
            canonical.display()
        ));
    }

    let total_length = sequence.len();
    let mut bounded_start = start.max(1).min(total_length);
    let mut bounded_end = end.max(1).min(total_length);
    if bounded_end < bounded_start {
        std::mem::swap(&mut bounded_start, &mut bounded_end);
    }
    if bounded_end.saturating_sub(bounded_start) + 1 > MAX_FASTA_SLICE_BP {
        bounded_end = (bounded_start + MAX_FASTA_SLICE_BP - 1).min(total_length);
    }

    let slice_start = bounded_start.saturating_sub(1);
    let slice_end = bounded_end;
    let slice = sequence
        .get(slice_start..slice_end)
        .ok_or_else(|| "Requested FASTA range is invalid".to_string())?
        .to_string();

    Ok(FastaRangeResult {
        path: canonical.to_string_lossy().to_string(),
        record_name: selected_record,
        start: bounded_start,
        end: bounded_end,
        total_length,
        sequence: slice,
    })
}

#[tauri::command]
pub fn open_file_location(path: String) -> Result<(), String> {
    let target = PathBuf::from(&path);
    if !target.exists() {
        return Err(format!("Path does not exist: {}", path));
    }

    #[cfg(target_os = "macos")]
    {
        let status = std::process::Command::new("open")
            .arg("-R")
            .arg(&target)
            .status()
            .map_err(|e| format!("Failed to open Finder: {}", e))?;
        if !status.success() {
            return Err("Failed to reveal file in Finder".to_string());
        }
        return Ok(());
    }

    #[cfg(target_os = "windows")]
    {
        let mut cmd_arg = String::from("/select,");
        cmd_arg.push_str(&target.to_string_lossy());
        let status = std::process::Command::new("explorer")
            .arg(cmd_arg)
            .status()
            .map_err(|e| format!("Failed to open Explorer: {}", e))?;
        if !status.success() {
            return Err("Failed to reveal file in Explorer".to_string());
        }
        return Ok(());
    }

    #[cfg(target_os = "linux")]
    {
        let open_path = if target.is_dir() {
            target.clone()
        } else {
            target
                .parent()
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from("."))
        };
        let status = std::process::Command::new("xdg-open")
            .arg(&open_path)
            .status()
            .map_err(|e| format!("Failed to open file manager: {}", e))?;
        if !status.success() {
            return Err("Failed to open file manager".to_string());
        }
        return Ok(());
    }

    #[allow(unreachable_code)]
    Err("Open file location is not supported on this platform".to_string())
}

#[tauri::command]
pub fn open_external_url(url: String) -> Result<(), String> {
    let target = url.trim();
    if target.is_empty() {
        return Err("URL is empty".to_string());
    }

    let lower = target.to_ascii_lowercase();
    let allowed = lower.starts_with("https://")
        || lower.starts_with("http://")
        || lower.starts_with("mailto:");
    if !allowed {
        return Err("Only http(s) and mailto URLs are allowed".to_string());
    }

    #[cfg(target_os = "macos")]
    {
        let status = std::process::Command::new("open")
            .arg(target)
            .status()
            .map_err(|e| format!("Failed to open URL in browser: {}", e))?;
        if !status.success() {
            return Err("Failed to open URL in browser".to_string());
        }
        return Ok(());
    }

    #[cfg(target_os = "windows")]
    {
        let status = std::process::Command::new("explorer")
            .arg(target)
            .status()
            .map_err(|e| format!("Failed to open URL in browser: {}", e))?;
        if !status.success() {
            return Err("Failed to open URL in browser".to_string());
        }
        return Ok(());
    }

    #[cfg(target_os = "linux")]
    {
        let status = std::process::Command::new("xdg-open")
            .arg(target)
            .status()
            .map_err(|e| format!("Failed to open URL in browser: {}", e))?;
        if !status.success() {
            return Err("Failed to open URL in browser".to_string());
        }
        return Ok(());
    }

    #[allow(unreachable_code)]
    Err("Open external URL is not supported on this platform".to_string())
}

// ============================================================================
// SYSTEM INFO
// ============================================================================

#[tauri::command]
pub fn get_system_usage() -> Result<util::SystemUsage, String> {
    util::read_system_usage()
}

#[tauri::command]
pub fn get_version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

#[tauri::command]
pub fn get_app_info() -> serde_json::Value {
    serde_json::json!({
        "name": "Pathotypr",
        "version": env!("CARGO_PKG_VERSION"),
        "description": "A high-performance tool for genome classification and variant genotyping",
        "authors": "Paula Ruiz Rodriguez & Mireia Coscolla (PathoGenOmics Lab)",
        "license": "AGPL-3.0"
    })
}

// ============================================================================
// MARKER DEPOSIT RESOLUTION
// ============================================================================

#[derive(Debug, Serialize)]
pub struct ZenodoAsset {
    pub kind: String,
    pub filename: String,
    pub url: String,
    /// True when the deposit does not publish this asset and a pinned URL is served instead.
    pub fallback: bool,
}

#[derive(Debug, Serialize)]
pub struct ZenodoResolution {
    pub record_id: String,
    pub version: Option<String>,
    pub assets: Vec<ZenodoAsset>,
    /// True when Zenodo could not be reached and the pinned URLs are being served.
    pub fallback: bool,
    /// Why the newest version could not be resolved, so the caller can report it.
    pub reason: Option<String>,
}

fn pinned_resolution(reason: &str) -> ZenodoResolution {
    log::warn!("Could not resolve the newest marker deposit ({reason}); using pinned URLs");
    let assets = ["lineage_markers", "dr_markers", "rf_model"]
        .iter()
        .filter_map(|kind| {
            defaults::fallback_asset(kind).map(|(url, filename)| ZenodoAsset {
                kind: kind.to_string(),
                filename: filename.to_string(),
                url: url.to_string(),
                fallback: true,
            })
        })
        .collect();
    ZenodoResolution {
        record_id: String::new(),
        version: None,
        assets,
        fallback: true,
        reason: Some(reason.to_string()),
    }
}

/// Resolve the newest published version of the marker deposit.
///
/// Assets are matched by filename prefix rather than by a version written into the code,
/// so publishing a new catalogue on Zenodo is enough to make the app offer it. If the
/// deposit cannot be reached the pinned URLs are returned instead, flagged as a fallback
/// so the caller can say so.
#[tauri::command]
pub async fn resolve_marker_assets() -> Result<ZenodoResolution, String> {
    // Zenodo's API answers 403 to requests without a User-Agent.
    let client = match reqwest::Client::builder()
        .user_agent(format!("Pathotypr/{}", env!("CARGO_PKG_VERSION")))
        .timeout(std::time::Duration::from_secs(15))
        .build()
    {
        Ok(c) => c,
        Err(e) => return Ok(pinned_resolution(&e.to_string())),
    };

    let response = match client.get(defaults::ZENODO_LATEST_VERSION_API).send().await {
        Ok(r) if r.status().is_success() => r,
        Ok(r) => return Ok(pinned_resolution(&format!("HTTP {}", r.status()))),
        Err(e) => return Ok(pinned_resolution(&e.to_string())),
    };

    // reqwest is built without its json feature here, so the body is parsed explicitly.
    let body = match response.text().await {
        Ok(b) => b,
        Err(e) => return Ok(pinned_resolution(&e.to_string())),
    };
    let record: serde_json::Value = match serde_json::from_str(&body) {
        Ok(v) => v,
        Err(e) => return Ok(pinned_resolution(&e.to_string())),
    };

    let files = record["files"].as_array().cloned().unwrap_or_default();
    let filenames: Vec<String> = files
        .iter()
        .filter_map(|f| f["key"].as_str().map(str::to_string))
        .collect();
    let borrowed: Vec<&str> = filenames.iter().map(String::as_str).collect();

    let record_id = record["id"]
        .as_u64()
        .map(|n| n.to_string())
        .or_else(|| record["id"].as_str().map(str::to_string))
        .unwrap_or_default();

    let mut assets = Vec::new();
    for kind in ["lineage_markers", "dr_markers", "rf_model"] {
        // A deposit can be published missing one of the assets, as happened when a new
        // version carried the resistance catalogue but not the lineage panel or the model.
        // That is no reason to fail the other downloads, so each asset falls back on its own.
        let Some(filename) = defaults::select_asset(kind, borrowed.iter().copied()) else {
            if let Some((url, filename)) = defaults::fallback_asset(kind) {
                log::warn!("The current marker deposit does not publish a {kind} file; using the pinned URL");
                assets.push(ZenodoAsset {
                    kind: kind.to_string(),
                    filename: filename.to_string(),
                    url: url.to_string(),
                    fallback: true,
                });
            }
            continue;
        };
        let url = files
            .iter()
            .find(|f| f["key"].as_str() == Some(filename))
            .and_then(|f| f["links"]["self"].as_str())
            .map(str::to_string)
            .unwrap_or_else(|| {
                format!("https://zenodo.org/records/{record_id}/files/{filename}?download=1")
            });
        assets.push(ZenodoAsset {
            kind: kind.to_string(),
            filename: filename.to_string(),
            url,
            fallback: false,
        });
    }

    if assets.iter().all(|a| a.fallback) {
        return Ok(pinned_resolution("the deposit listed no recognisable assets"));
    }

    Ok(ZenodoResolution {
        record_id,
        version: record["metadata"]["version"]
            .as_str()
            .map(str::to_string),
        assets,
        fallback: false,
        reason: None,
    })
}
