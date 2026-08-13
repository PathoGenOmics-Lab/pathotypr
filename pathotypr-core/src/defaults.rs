//! Default locations for Pathotypr data resources.
//!
//! Markers and the pre-trained model live in a versioned Zenodo deposit. The constants
//! here name the deposit and describe how to recognise each asset inside it, rather than
//! naming one particular version, so that publishing a new catalogue does not require a
//! code change. Callers resolve the newest version at download time; the pinned URLs are
//! the fallback for when Zenodo cannot be reached.

/// Concept record of the marker deposit. Unlike a version record, this identifier is
/// stable across releases and always has a newest version.
pub const ZENODO_CONCEPT_RECORD: &str = "19210043";

/// API endpoint that redirects to the newest published version of the deposit.
pub const ZENODO_LATEST_VERSION_API: &str =
    "https://zenodo.org/api/records/19210043/versions/latest";

/// Landing page of the newest version, for documentation and links.
pub const ZENODO_CONCEPT_DOI: &str = "https://doi.org/10.5281/zenodo.19210043";

/// Assets are recognised by the prefix of their filename, because filenames carry the
/// catalogue version (`pathotypr_dr_markers_ancestor_v1.0.2.tsv`).
pub const LINEAGE_MARKERS_PREFIX: &str = "pathotypr_lineage_markers";
pub const DR_MARKERS_PREFIX: &str = "pathotypr_dr_markers";
pub const RF_MODEL_PREFIX: &str = "pathotypr_rf_model";

/// The resistance catalogue is published in two coordinate frames. The desktop app pairs
/// markers with the ancestor reference below, so it wants the ancestor frame; deposits
/// before v1.0.2 carry a single unsuffixed file, which is also the ancestor frame.
pub const DR_MARKERS_ANCESTOR_PREFIX: &str = "pathotypr_dr_markers_ancestor";
pub const DR_MARKERS_H37RV_MARKER: &str = "H37Rv";

/// Pinned fallbacks, used only when the newest version cannot be resolved.
pub const LINEAGE_MARKERS_FALLBACK_URL: &str =
    "https://zenodo.org/records/19210044/files/pathotypr_lineage_markers_v1.0.0.tsv?download=1";
pub const LINEAGE_MARKERS_FALLBACK_FILENAME: &str = "pathotypr_lineage_markers_v1.0.0.tsv";

pub const DR_MARKERS_FALLBACK_URL: &str =
    "https://zenodo.org/records/19210044/files/pathotypr_dr_markers_v1.0.0.tsv?download=1";
pub const DR_MARKERS_FALLBACK_FILENAME: &str = "pathotypr_dr_markers_v1.0.0.tsv";

pub const RF_MODEL_FALLBACK_URL: &str =
    "https://zenodo.org/records/19210044/files/pathotypr_rf_model_v1.0.0.pathotypr?download=1";
pub const RF_MODEL_FALLBACK_FILENAME: &str = "pathotypr_rf_model_v1.0.0.pathotypr";

/// MTBC ancestor reference genome. Its deposit is external and not versioned by us.
pub const REFERENCE_URL: &str =
    "https://zenodo.org/records/3497110/files/MTB_ancestor_reference.fasta?download=1";
pub const REFERENCE_FILENAME: &str = "MTB_ancestor_reference.fasta";

/// Pick the asset of a given kind from the filenames published in a deposit.
///
/// Returns the filename, which the caller turns into a download URL. `dr_markers`
/// prefers the ancestor coordinate frame and never returns the H37Rv file.
pub fn select_asset<'a, I>(kind: &str, filenames: I) -> Option<&'a str>
where
    I: IntoIterator<Item = &'a str> + Clone,
{
    match kind {
        "lineage_markers" => filenames
            .into_iter()
            .find(|f| f.starts_with(LINEAGE_MARKERS_PREFIX)),
        "rf_model" => filenames
            .into_iter()
            .find(|f| f.starts_with(RF_MODEL_PREFIX)),
        "dr_markers" => filenames
            .clone()
            .into_iter()
            .find(|f| f.starts_with(DR_MARKERS_ANCESTOR_PREFIX))
            .or_else(|| {
                filenames.into_iter().find(|f| {
                    f.starts_with(DR_MARKERS_PREFIX) && !f.contains(DR_MARKERS_H37RV_MARKER)
                })
            }),
        _ => None,
    }
}

/// Pinned fallback for a kind, as (url, filename).
pub fn fallback_asset(kind: &str) -> Option<(&'static str, &'static str)> {
    match kind {
        "lineage_markers" => Some((LINEAGE_MARKERS_FALLBACK_URL, LINEAGE_MARKERS_FALLBACK_FILENAME)),
        "dr_markers" => Some((DR_MARKERS_FALLBACK_URL, DR_MARKERS_FALLBACK_FILENAME)),
        "rf_model" => Some((RF_MODEL_FALLBACK_URL, RF_MODEL_FALLBACK_FILENAME)),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dr_markers_prefers_the_ancestor_frame() {
        let files = vec![
            "pathotypr_dr_markers_H37Rv_v1.0.2.tsv",
            "pathotypr_dr_markers_ancestor_v1.0.2.tsv",
            "pathotypr_lineage_markers_v1.0.0.tsv",
            "README.md",
        ];
        assert_eq!(
            select_asset("dr_markers", files.iter().copied()),
            Some("pathotypr_dr_markers_ancestor_v1.0.2.tsv")
        );
    }

    #[test]
    fn dr_markers_accepts_the_unsuffixed_name_of_older_deposits() {
        let files = vec!["pathotypr_dr_markers_v1.0.0.tsv", "README.md"];
        assert_eq!(
            select_asset("dr_markers", files.iter().copied()),
            Some("pathotypr_dr_markers_v1.0.0.tsv")
        );
    }

    #[test]
    fn dr_markers_never_returns_the_h37rv_frame() {
        let files = vec!["pathotypr_dr_markers_H37Rv_v1.0.2.tsv", "README.md"];
        assert_eq!(select_asset("dr_markers", files.iter().copied()), None);
    }

    #[test]
    fn lineage_markers_and_model_are_matched_by_prefix() {
        let files = vec![
            "pathotypr_lineage_markers_v1.0.0.tsv",
            "pathotypr_rf_model_v1.0.0.pathotypr",
        ];
        assert_eq!(
            select_asset("lineage_markers", files.iter().copied()),
            Some("pathotypr_lineage_markers_v1.0.0.tsv")
        );
        assert_eq!(
            select_asset("rf_model", files.iter().copied()),
            Some("pathotypr_rf_model_v1.0.0.pathotypr")
        );
    }
}
