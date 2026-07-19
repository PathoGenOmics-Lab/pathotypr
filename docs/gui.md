# Desktop GUI

pathotypr includes a native desktop application built with [Tauri](https://tauri.app/). It wraps the same `pathotypr-core` engine as the CLI, so every workflow produces identical results — just with drag-and-drop, live progress, and interactive tables.

!!! tip "You probably don't need to build this"
    Signed-format installers for macOS, Linux, and Windows are published for
    every release. Grab one from the [Installation](installation.md) page and
    skip straight to using the app — the steps below are only for building from
    source.

## Features

- All five workflows: [train](train.md), [predict](predict.md), [classify](classify.md), [split-fastq](split-fastq.md), [match](match.md)
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

1. **Rust** — a recent stable toolchain (install via [rustup](https://rustup.rs/)).
2. **Tauri CLI v2**:

    ```bash
    cargo install tauri-cli --version "^2"
    ```

    !!! warning "Must be the 2.x CLI"
        pathotypr targets **Tauri 2**. A `tauri-cli` from the 1.x series already
        on your `PATH` will fail to build this project. Confirm your version
        with `cargo tauri --version` before continuing.

3. **System dependencies** — these vary by platform:

    === ":material-linux: Linux (Debian/Ubuntu)"

        ```bash
        sudo apt install libgtk-3-dev libwebkit2gtk-4.1-dev libjavascriptcoregtk-4.1-dev
        ```

        !!! note "These are the 4.1 packages"
            Tauri 2 links against the **4.1** WebKitGTK / JavaScriptCore
            packages. The older `4.0` variants (all that some distributions,
            e.g. Ubuntu 20.04, ship by default) will **not** work — use a newer
            release or a backport that provides the `4.1` packages.

    === ":material-apple: macOS"

        The Xcode Command Line Tools are usually pre-installed. If they are
        missing, install them with:

        ```bash
        xcode-select --install
        ```

    === ":material-microsoft-windows: Windows"

        The **WebView2 runtime** is required. It is pre-installed on Windows 11
        and most up-to-date Windows 10 machines. If it is missing, install the
        Evergreen runtime from
        [Microsoft](https://developer.microsoft.com/en-us/microsoft-edge/webview2/).

### Development

```bash
cargo tauri dev
```

Opens the app in development mode with hot-reload for the frontend.

!!! note "First build is slow"
    The initial `cargo tauri dev` compiles the full Rust dependency tree and can
    take several minutes. Subsequent runs are incremental and start quickly.

### Production build

```bash
cargo tauri build
```

With `targets: "all"` (set in `tauri.conf.json`), this produces distributable
installers under `src-tauri/target/release/bundle/`:

- **macOS**: `.dmg` in `bundle/dmg/` (and the `.app` in `bundle/macos/`)
- **Linux**: `.deb` in `bundle/deb/`, `.rpm` in `bundle/rpm/`, and `.AppImage` in `bundle/appimage/`
- **Windows**: `.msi` in `bundle/msi/` and the NSIS `.exe` in `bundle/nsis/`

!!! warning "Self-built apps are unsigned"
    Binaries you build yourself are not signed with a paid developer
    certificate, so macOS Gatekeeper and Windows SmartScreen will warn you the
    first time you open them. See the
    [first-launch notes](installation.md#desktop-gui-pre-built) on the
    Installation page for how to get past the warning.

## Architecture

```text
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

## See also

- [Installation](installation.md) — pre-built desktop installers for macOS, Linux, and Windows
- [Input formats](input-formats.md) — accepted assemblies, reads, and marker files
- Command references: [train](train.md) · [predict](predict.md) · [classify](classify.md) · [split-fastq](split-fastq.md) · [match](match.md)
