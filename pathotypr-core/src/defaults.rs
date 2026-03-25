//! Default URLs and filenames for Pathotypr data resources.
//!
//! These constants provide a single source of truth for download URLs,
//! shared between the CLI and GUI frontends.

/// Zenodo record for v1.0.0 markers and model
pub const ZENODO_RECORD: &str = "https://zenodo.org/records/19210044/files";

/// Lineage markers TSV
pub const LINEAGE_MARKERS_URL: &str =
    "https://zenodo.org/records/19210044/files/pathotypr_lineage_markers_v1.0.0.tsv?download=1";
pub const LINEAGE_MARKERS_FILENAME: &str = "pathotypr_lineage_markers_v1.0.0.tsv";

/// Drug resistance markers TSV
pub const DR_MARKERS_URL: &str =
    "https://zenodo.org/records/19210044/files/pathotypr_dr_markers_v1.0.0.tsv?download=1";
pub const DR_MARKERS_FILENAME: &str = "pathotypr_dr_markers_v1.0.0.tsv";

/// Pre-trained Random Forest model
pub const RF_MODEL_URL: &str =
    "https://zenodo.org/records/19210044/files/pathotypr_rf_model_v1.0.0.pathotypr?download=1";
pub const RF_MODEL_FILENAME: &str = "pathotypr_rf_model_v1.0.0.pathotypr";

/// MTB ancestor reference genome
pub const REFERENCE_URL: &str =
    "https://zenodo.org/records/3497110/files/MTB_ancestor_reference.fasta?download=1";
pub const REFERENCE_FILENAME: &str = "MTB_ancestor_reference.fasta";
