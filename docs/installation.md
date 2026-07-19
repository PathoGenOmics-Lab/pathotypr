# Installation

pathotypr ships as a **native desktop app**, a **command-line binary**, and as
**source** you can build yourself.

## Desktop GUI (pre-built)

Download the latest release for your platform:

| Platform | Download | Notes |
|---|---|---|
| :material-apple: macOS (Apple Silicon) | [**Pathotypr_1.0.0_aarch64.dmg**](https://github.com/PathoGenOmics-Lab/pathotypr/releases/latest/download/Pathotypr_1.0.0_aarch64.dmg) | M1 / M2 / M3 / M4 Macs |
| :material-apple: macOS (Intel) | [**Pathotypr_1.0.0_x64.dmg**](https://github.com/PathoGenOmics-Lab/pathotypr/releases/latest/download/Pathotypr_1.0.0_x64.dmg) | Pre-2020 Macs |
| :material-linux: Linux (.deb) | [**Pathotypr_1.0.0_amd64.deb**](https://github.com/PathoGenOmics-Lab/pathotypr/releases/latest/download/Pathotypr_1.0.0_amd64.deb) | Debian / Ubuntu |
| :material-linux: Linux (.rpm) | [**Pathotypr-1.0.0-1.x86_64.rpm**](https://github.com/PathoGenOmics-Lab/pathotypr/releases/latest/download/Pathotypr-1.0.0-1.x86_64.rpm) | Fedora / RHEL |
| :material-linux: Linux (AppImage) | [**Pathotypr_1.0.0_amd64.AppImage**](https://github.com/PathoGenOmics-Lab/pathotypr/releases/latest/download/Pathotypr_1.0.0_amd64.AppImage) | Any distro, no install needed |
| :material-microsoft-windows: Windows (installer) | [**Pathotypr_1.0.0_x64-setup.exe**](https://github.com/PathoGenOmics-Lab/pathotypr/releases/latest/download/Pathotypr_1.0.0_x64-setup.exe) | Windows 10+ |
| :material-microsoft-windows: Windows (.msi) | [**Pathotypr_1.0.0_x64_en-US.msi**](https://github.com/PathoGenOmics-Lab/pathotypr/releases/latest/download/Pathotypr_1.0.0_x64_en-US.msi) | Windows 10+ (MSI) |

All builds are on the [**Releases page**](https://github.com/PathoGenOmics-Lab/pathotypr/releases).

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

| File | Description | Download |
|---|---|---|
| `pathotypr_lineage_markers_v1.0.0.tsv` | 3,707 lineage SNPs (L1–L10, A1–A4) | [:material-download: Download](https://zenodo.org/records/19210044/files/pathotypr_lineage_markers_v1.0.0.tsv?download=1) |
| `pathotypr_dr_markers_v1.0.0.tsv` | 102,213 DR mutations (WHO catalogue 2021) | [:material-download: Download](https://zenodo.org/records/19210044/files/pathotypr_dr_markers_v1.0.0.tsv?download=1) |
| `pathotypr_rf_model_v1.0.0.pathotypr` | Pre-trained RF model (k=31, 100 trees) | [:material-download: Download](https://zenodo.org/records/19210044/files/pathotypr_rf_model_v1.0.0.pathotypr?download=1) |

!!! info "DOI"
    [10.5281/zenodo.19210044](https://zenodo.org/records/19210044)

Once installed, head to the [command guides](train.md) or the
[input formats](input-formats.md) reference to get started.
