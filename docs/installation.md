# Installation

pathotypr ships as a **native desktop app**, a **command-line binary**, and as
**source** you can build yourself.

## Desktop GUI (pre-built)

[:material-download: **Download the latest release**](https://github.com/PathoGenOmics-Lab/pathotypr/releases/latest){ .md-button .md-button--primary }

Then pick the file for your platform:

| Platform | File | Notes |
|---|---|---|
| :material-apple: macOS (Apple Silicon) | `Pathotypr_<version>_aarch64.dmg` | M1 / M2 / M3 / M4 Macs |
| :material-apple: macOS (Intel) | `Pathotypr_<version>_x64.dmg` | Pre-2020 Macs |
| :material-linux: Linux (.deb) | `Pathotypr_<version>_amd64.deb` | Debian / Ubuntu |
| :material-linux: Linux (.rpm) | `Pathotypr-<version>-1.x86_64.rpm` | Fedora / RHEL |
| :material-linux: Linux (AppImage) | `Pathotypr_<version>_amd64.AppImage` | Any distro, no install needed |
| :material-microsoft-windows: Windows (installer) | `Pathotypr_<version>_x64-setup.exe` | Windows 10+ |
| :material-microsoft-windows: Windows (.msi) | `Pathotypr_<version>_x64_en-US.msi` | Windows 10+ (MSI) |

Older versions are on the [releases page](https://github.com/PathoGenOmics-Lab/pathotypr/releases).

!!! warning "First launch on macOS and Windows"
    The app is **not** signed with a paid developer certificate, so the OS may
    warn you the first time you open it.

    - **macOS**: right-click the app → **Open** → click **Open** in the dialog.
      See [Apple support](https://support.apple.com/en-us/HT202491) for details.
    - **Windows**: if SmartScreen appears, click **More info** → **Run anyway**.

## CLI (Bioconda)

```bash
conda create -n pathotypr -c bioconda pathotypr
conda activate pathotypr
pathotypr --help
```

## From source

=== "CLI"

    ```bash
    git clone https://github.com/PathoGenOmics-Lab/pathotypr.git
    cd pathotypr
    cargo build --release -p pathotypr-core --bin pathotypr
    ./target/release/pathotypr --help
    ```

    !!! note "Requirements"
        A recent [Rust toolchain](https://rustup.rs/) (stable). No other system
        dependencies are needed for the CLI.

=== "Desktop GUI"

    !!! tip "Most people should download the app instead"
        Ready-built installers are published for every release — see
        [Desktop GUI (pre-built)](#desktop-gui-pre-built) above. Build from
        source only if you are modifying pathotypr itself.

    ```bash
    # Development build with hot reload
    cargo tauri dev

    # Production build (installers/bundles)
    cargo tauri build
    ```

    The GUI has additional system dependencies (WebView, etc.). See the
    [Desktop GUI](gui.md) guide for the full per-platform setup.

## MTBC marker files & pre-trained model

Ready-to-use marker panels and a pre-trained Random Forest model for
*Mycobacterium tuberculosis* complex (MTBC) are published on Zenodo:

| File | Description |
|---|---|
| `pathotypr_lineage_markers_*.tsv` | 3,707 lineage SNPs (L1–L10, A1–A4) |
| `pathotypr_dr_markers_ancestor_*.tsv` | DR mutations from the WHO catalogue (2nd edition, 2023), in ancestor coordinates |
| `pathotypr_dr_markers_H37Rv_*.tsv` | The same catalogue in H37Rv coordinates |
| `pathotypr_rf_model_*.pathotypr` | Pre-trained RF model (k=31, 100 trees) |

!!! info "DOI"
    [10.5281/zenodo.19210043](https://doi.org/10.5281/zenodo.19210043)

Once installed, head to the [command guides](train.md) or the
[input formats](input-formats.md) reference to get started.
