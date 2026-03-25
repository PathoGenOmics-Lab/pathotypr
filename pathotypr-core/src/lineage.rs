//! Shared lineage classification logic.
//!
//! Used by both `classify` (FASTA-based) and `classify_split_fastq` (FASTQ-based)
//! workflows to determine the final lineage call from marker counts.

use std::collections::{HashMap, HashSet};

/// Determines the best lineage call using nested hierarchical classification.
///
/// Selects the deepest lineage with a fully supported ancestor path.
/// Ties are broken by SNP count, then lexicographic order.
/// Falls back to the most abundant lineage when no valid nested path exists.
pub fn get_final_lineage_call(lineage_counts: &HashMap<String, usize>) -> String {
    if lineage_counts.is_empty() {
        return "Unclassified".to_string();
    }

    let supported_lineages: HashSet<&str> = lineage_counts.keys().map(|s| s.as_str()).collect();
    let mut candidate_lineages: Vec<&String> = lineage_counts.keys().collect();
    candidate_lineages.sort();
    let mut valid_candidates: Vec<&String> = Vec::new();

    // 1. Identify all lineages with a valid, fully supported path.
    for candidate in candidate_lineages {
        let mut is_path_valid = true;
        if candidate.contains(';') {
            let components: Vec<&str> = candidate.split(';').collect();
            for i in 1..components.len() {
                let parent_path = components[0..i].join(";");
                if !supported_lineages.contains(parent_path.as_str()) {
                    is_path_valid = false;
                    break;
                }
            }
        }
        if is_path_valid {
            valid_candidates.push(candidate);
        }
    }

    // If no valid paths exist, fallback to the most abundant lineage.
    if valid_candidates.is_empty() {
        let mut sorted_lineages: Vec<(&String, &usize)> = lineage_counts.iter().collect();
        sorted_lineages.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
        return sorted_lineages
            .first()
            .map(|(lineage, _)| (*lineage).clone())
            .unwrap_or_else(|| "Unclassified".to_string());
    }

    // 2. From the valid candidates, find the deepest one.
    // If there's a tie in depth, the higher SNP count for that specific deep lineage wins.
    valid_candidates.sort_by(|a, b| {
        let depth_a = a.split(';').count();
        let depth_b = b.split(';').count();
        depth_b
            .cmp(&depth_a)
            .then_with(|| lineage_counts[*b].cmp(&lineage_counts[*a]))
            .then_with(|| a.cmp(b))
    });
    let best_deep_lineage = valid_candidates[0];

    // 3. Find the most abundant lineage overall.
    let mut sorted_lineages: Vec<(&String, &usize)> = lineage_counts.iter().collect();
    sorted_lineages.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
    let most_abundant_lineage = sorted_lineages[0].0;

    // 4. The final decision logic.
    let best_deep_count = lineage_counts[best_deep_lineage];
    let most_abundant_count = lineage_counts[most_abundant_lineage];

    // Check if the most abundant lineage is from a different branch than the deepest one.
    // A shared branch means one is a prefix of the other, delimited by ';'.
    let shares_branch = best_deep_lineage == most_abundant_lineage
        || best_deep_lineage.starts_with(&format!("{};", most_abundant_lineage))
        || most_abundant_lineage.starts_with(&format!("{};", best_deep_lineage));

    if !shares_branch && most_abundant_count > best_deep_count {
        most_abundant_lineage.clone()
    } else {
        best_deep_lineage.clone()
    }
}

/// Determines the major lineage from lineage counts.
///
/// When `nested_classification` is true, uses hierarchical path validation.
/// Otherwise, returns the most abundant lineage (lexicographic tiebreak).
pub fn determine_major_lineage(
    lineage_counts: &HashMap<String, usize>,
    nested_classification: bool,
) -> String {
    if nested_classification {
        get_final_lineage_call(lineage_counts)
    } else {
        let mut sorted_lineages: Vec<(&String, &usize)> = lineage_counts.iter().collect();
        sorted_lineages.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
        if sorted_lineages.is_empty() {
            "Unclassified".to_string()
        } else {
            sorted_lineages[0].0.clone()
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nested_call_is_deterministic_on_ties() {
        let mut counts = HashMap::new();
        counts.insert("L1".to_string(), 5);
        counts.insert("L2".to_string(), 5);
        assert_eq!(get_final_lineage_call(&counts), "L1");
    }

    #[test]
    fn nested_call_uses_lexical_tiebreak_for_same_depth() {
        let mut counts = HashMap::new();
        counts.insert("L1;A".to_string(), 3);
        counts.insert("L1;B".to_string(), 3);
        counts.insert("L1".to_string(), 1);
        assert_eq!(get_final_lineage_call(&counts), "L1;A");
    }

    #[test]
    fn empty_counts_returns_unclassified() {
        let counts = HashMap::new();
        assert_eq!(get_final_lineage_call(&counts), "Unclassified");
    }

    #[test]
    fn determine_major_without_nesting() {
        let mut counts = HashMap::new();
        counts.insert("L1".to_string(), 3);
        counts.insert("L2".to_string(), 5);
        assert_eq!(determine_major_lineage(&counts, false), "L2");
    }

    #[test]
    fn lineage_prefix_does_not_match_different_branch() {
        // L1 should NOT be considered the same branch as L12
        let mut counts = HashMap::new();
        counts.insert("L1".to_string(), 3);
        counts.insert("L12".to_string(), 10);
        // L12 has more counts and is a different branch, so it should win
        assert_eq!(get_final_lineage_call(&counts), "L12");
    }

    #[test]
    fn lineage_prefix_matches_child_with_delimiter() {
        // L1 IS the same branch as L1;L1.1
        let mut counts = HashMap::new();
        counts.insert("L1".to_string(), 3);
        counts.insert("L1;L1.1".to_string(), 10);
        // L1;L1.1 is deeper and same branch, should win
        assert_eq!(get_final_lineage_call(&counts), "L1;L1.1");
    }
}
