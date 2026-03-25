// Prevents additional console window on Windows in release
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod commands;
mod state;
mod util;

use std::sync::Arc;

fn main() {
    env_logger::init();
    if let Some(exit_code) = commands::maybe_run_match_worker_from_cli() {
        std::process::exit(exit_code);
    }

    let task_state = Arc::new(state::TaskState::default());

    tauri::Builder::default()
        .manage(task_state)
        .plugin(tauri_plugin_dialog::init())
        .invoke_handler(tauri::generate_handler![
            commands::run_train,
            commands::run_predict,
            commands::run_classify,
            commands::run_split_fastq,
            commands::run_match,
            commands::download_file,
            commands::cancel_task,
            commands::is_task_running,
            commands::read_text_file,
            commands::read_fasta_range,
            commands::open_file_location,
            commands::open_external_url,
            commands::get_system_usage,
            commands::get_version,
            commands::get_app_info,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
