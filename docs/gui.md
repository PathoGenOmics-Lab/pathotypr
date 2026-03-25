# Desktop GUI

pathotypr includes a native desktop application built with [Tauri](https://tauri.app/).

## Features

- All five workflows: train, predict, classify, split-fastq, match
- Drag-and-drop file selection
- Interactive result tables with sorting and filtering
- Real-time progress bars with cancellation support
- Training summary card (accuracy, OOB, CV metrics, model size)
- Live CPU/RAM usage indicator
- Light and dark themes
- Excel export toggle per workflow
- Configurable parameters with sensible defaults and reset buttons

## Building

### Prerequisites

1. **Rust** stable toolchain
2. **tauri-cli**:
   ```bash
   cargo install tauri-cli --version "^2"
   ```
3. **System dependencies** (varies by platform):

   **Linux (Debian/Ubuntu):**
   ```bash
   sudo apt install libgtk-3-dev libwebkit2gtk-4.1-dev libjavascriptcoregtk-4.1-dev
   ```

   **macOS:** Xcode command line tools (usually pre-installed)

   **Windows:** WebView2 runtime (usually pre-installed on Windows 10+)

### Development

```bash
cargo tauri dev
```

Opens the app in development mode with hot-reload for the frontend.

### Production build

```bash
cargo tauri build
```

Produces a distributable binary:
- **macOS**: `.dmg` in `src-tauri/target/release/bundle/dmg/`
- **Linux**: `.deb` and `.AppImage` in `src-tauri/target/release/bundle/`
- **Windows**: `.msi` in `src-tauri/target/release/bundle/msi/`

## Architecture

```
src-tauri/
├── src/
│   ├── main.rs       # Tauri app setup
│   ├── commands.rs   # Tauri command handlers (bridge GUI → core)
│   ├── state.rs      # Task state management + cancellation
│   └── util.rs       # Progress events + helpers
├── Cargo.toml
└── tauri.conf.json

frontend/
├── index.html        # Main page layout
├── styles.css        # Styling
├── app.js            # App bootstrap
└── js/
    ├── main.js       # Initialization
    ├── forms.js      # Form submission handlers
    ├── results.js    # Result table rendering
    ├── visualization.js  # Charts + summary cards
    ├── progress.js   # Progress bars + cancellation
    ├── config.js     # Default parameters
    ├── state.js      # Run state management
    ├── navigation.js # Tab navigation
    ├── dropzone.js   # Drag-and-drop file handling
    ├── console.js    # Log console drawer
    ├── tauri.js      # Tauri API wrappers
    ├── theme.js      # Light/dark theme toggle
    └── utils.js      # Shared utilities
```

## How it works

1. Frontend collects form parameters and sends them to Tauri via `invoke()`
2. `commands.rs` deserializes parameters and calls `pathotypr_core` functions
3. Core library runs the workflow, emitting progress events back to the frontend
4. Results are displayed in interactive tables and summary cards
5. Cancellation token allows stopping any running workflow via the UI
