//! Reading a WHO drug-resistance panel for what it is.
//!
//! A marker panel is a position, two alleles and a path. For a lineage panel
//! that path is a lineage, and tallying it gives a lineage call. For the WHO
//! resistance catalogue the same columns hold
//!
//! ```text
//! drug ; resistance ; marker ; WHO grade ; gene ; mutation
//! ```
//!
//! which is not a hierarchy at all. Counting it the lineage way answers "which
//! drug has the most markers", and most markers of the catalogue are graded 4
//! or 5, meaning *not* associated with resistance. A genome carrying
//! rpoB D435Y and katG S315T, both grade 1, came out as `OTHER` because an
//! ungraded marker of another drug happened to be more numerous.
//!
//! So a resistance panel is recognised here and answered on its own terms: the
//! drugs with a grade 1 or 2 marker, which is the call the WHO catalogue is
//! built to give.

use std::collections::HashMap;

/// What one marker says about its drug. Ordered: a stronger verdict wins.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Verdict {
    /// WHO grades 4 and 5: not associated with resistance.
    NotAssociated,
    /// No grade that can be read.
    Ungraded,
    /// WHO grade 3: uncertain significance.
    Uncertain,
    /// WHO grades 1 and 2: associated with resistance.
    Resistant,
}

/// Reads a WHO grade cell, as the catalogue writes it: "1) Assoc w R",
/// "4) Not assoc w R - Interim", and the free text some entries carry.
pub fn classify_grade(grade: &str) -> Verdict {
    let grade = grade.trim();
    if grade.is_empty() {
        return Verdict::Ungraded;
    }
    if let Some(digit) = grade
        .split(')')
        .next()
        .map(str::trim)
        .and_then(|d| d.parse::<u8>().ok())
    {
        return match digit {
            1 | 2 => Verdict::Resistant,
            3 => Verdict::Uncertain,
            4 | 5 => Verdict::NotAssociated,
            _ => Verdict::Ungraded,
        };
    }
    let lower = grade.to_ascii_lowercase();
    if lower.contains("not assoc") {
        Verdict::NotAssociated
    } else if lower.contains("uncertain") {
        Verdict::Uncertain
    } else if lower.contains("assoc") {
        Verdict::Resistant
    } else {
        Verdict::Ungraded
    }
}

/// Whether a cell looks like a WHO grade. No lineage level ever does, which is
/// what makes it a reliable way to tell the two kinds of panel apart.
fn looks_like_grade(cell: &str) -> bool {
    let cell = cell.trim();
    if cell.is_empty() {
        return false;
    }
    if cell
        .split(')')
        .next()
        .map(str::trim)
        .and_then(|d| d.parse::<u8>().ok())
        .is_some_and(|d| (1..=5).contains(&d))
        && cell.contains(')')
    {
        return true;
    }
    let lower = cell.to_ascii_lowercase();
    lower.contains("assoc") || lower.contains("uncertain")
}

/// The drug and verdict of a resistance marker path, or `None` when the path
/// is not one.
pub fn parse_path(path: &str) -> Option<(String, Verdict)> {
    let parts: Vec<&str> = path.split(';').map(str::trim).collect();
    if parts.len() < 4 {
        return None;
    }
    if parts[0].is_empty() || !looks_like_grade(parts[3]) {
        return None;
    }
    Some((parts[0].to_ascii_uppercase(), classify_grade(parts[3])))
}

/// Whether a panel is the WHO resistance catalogue rather than a lineage one.
///
/// A clear majority has to parse as a resistance path, so a lineage panel with
/// one oddly annotated marker is never mistaken for a resistance panel.
pub fn is_resistance_panel<'a>(paths: impl Iterator<Item = &'a str>) -> bool {
    let (mut total, mut dr_like) = (0usize, 0usize);
    for path in paths {
        total += 1;
        if parse_path(path).is_some() {
            dr_like += 1;
        }
    }
    total > 0 && dr_like * 2 >= total
}

/// The strongest verdict found per drug, from the markers of one genome.
pub fn fold_marker(per_drug: &mut HashMap<String, Verdict>, path: &str) {
    if let Some((drug, verdict)) = parse_path(path) {
        per_drug
            .entry(drug)
            .and_modify(|current| {
                if verdict > *current {
                    *current = verdict;
                }
            })
            .or_insert(verdict);
    }
}

/// The `lineage:count` tally of a resistance run: markers found per drug.
///
/// A lineage path's prefixes are its ancestors, so counting every prefix
/// tallies each level of the hierarchy. A catalogue path is not a hierarchy,
/// and its prefixes run into the grade, which has spaces in it: the pairs of a
/// space-separated column came apart as `...;1)`, `Assoc`, `w` and `R:1`, and
/// the column could not be read back. The count at the drug level is already
/// there, since every marker counts towards each of its prefixes, so it is the
/// only one kept.
pub fn per_drug_counts(prefix_counts: &HashMap<String, usize>) -> HashMap<String, usize> {
    prefix_counts
        .iter()
        .filter(|(path, _)| !path.contains(';'))
        .map(|(drug, count)| (drug.clone(), *count))
        .collect()
}

/// The call for one genome: the drugs it carries a grade 1 or 2 marker for.
///
/// Says "found" rather than "none", because a marker that could not be read is
/// not a marker that is absent; the `undetermined` columns beside this one say
/// how much was not read.
pub fn call(per_drug: &HashMap<String, Verdict>) -> String {
    let mut resistant: Vec<&str> = per_drug
        .iter()
        .filter(|(_, v)| **v == Verdict::Resistant)
        .map(|(d, _)| d.as_str())
        .collect();
    if !resistant.is_empty() {
        resistant.sort_unstable();
        return format!("Resistant: {}", resistant.join(","));
    }
    let mut uncertain: Vec<&str> = per_drug
        .iter()
        .filter(|(_, v)| **v == Verdict::Uncertain)
        .map(|(d, _)| d.as_str())
        .collect();
    if !uncertain.is_empty() {
        uncertain.sort_unstable();
        return format!("No resistance markers found; uncertain: {}", uncertain.join(","));
    }
    "No resistance markers found".to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn who_grades_are_read_the_way_the_catalogue_writes_them() {
        assert_eq!(classify_grade("1) Assoc w R"), Verdict::Resistant);
        assert_eq!(classify_grade("2) Assoc w R - Interim"), Verdict::Resistant);
        assert_eq!(classify_grade("3) Uncertain significance"), Verdict::Uncertain);
        assert_eq!(classify_grade("4) Not assoc w R - Interim"), Verdict::NotAssociated);
        assert_eq!(classify_grade("5) Not assoc w R"), Verdict::NotAssociated);
        assert_eq!(classify_grade("TBProfiler_literature"), Verdict::Ungraded);
        assert_eq!(classify_grade(""), Verdict::Ungraded);
    }

    #[test]
    fn a_lineage_panel_is_never_taken_for_a_resistance_one() {
        let lineage = ["L4", "L4;L4.9", "A4;ancestorM_caprae", "L2;L2.2;L2.2.1"];
        assert!(!is_resistance_panel(lineage.into_iter()));
        for path in lineage {
            assert_eq!(parse_path(path), None, "{path}");
        }
    }

    #[test]
    fn a_resistance_panel_is_recognised_and_parsed() {
        let dr = [
            "RIF;RIF-R;rpoB_S450L;1) Assoc w R;rpoB;S450L",
            "INH;INH-R;ndh_c.168C>T;4) Not assoc w R - Interim;ndh;c.168C>T",
            "EMB;EMB-R;embB_M306V;3) Uncertain significance;embB;M306V",
        ];
        assert!(is_resistance_panel(dr.into_iter()));
        assert_eq!(parse_path(dr[0]), Some(("RIF".to_string(), Verdict::Resistant)));
        assert_eq!(parse_path(dr[1]), Some(("INH".to_string(), Verdict::NotAssociated)));
    }

    #[test]
    fn one_odd_marker_does_not_flip_a_lineage_panel() {
        let mostly_lineage = [
            "L4", "L4;L4.9", "L2", "L2;L2.2", "L1",
            "RIF;RIF-R;rpoB_S450L;1) Assoc w R;rpoB;S450L",
        ];
        assert!(!is_resistance_panel(mostly_lineage.into_iter()));
    }

    /// The case that prompted this module: a genome with two grade 1 markers
    /// used to be summarised as whichever drug had the most markers, most of
    /// them graded "not associated".
    #[test]
    fn the_call_names_the_drugs_with_grade_1_or_2_markers() {
        let mut per_drug = HashMap::new();
        for path in [
            "RIF;RIF-R;rpoB_D435Y;1) Assoc w R;rpoB;D435Y",
            "INH;INH-R;katG_S315T;1) Assoc w R;katG;S315T",
            "INH;INH-R;ndh_c.168C>T;4) Not assoc w R - Interim;ndh;c.168C>T",
            "OTHER;OTHER-R;mtrB_c.589C>T;4) Not assoc w R - Interim;mtrB;c.589C>T",
            "OTHER;OTHER-R;PPE35_I417V;3) Uncertain significance;PPE35;I417V",
            "EMB;EMB-R;embB_S565G;3) Uncertain significance;embB;S565G",
        ] {
            fold_marker(&mut per_drug, path);
        }
        assert_eq!(call(&per_drug), "Resistant: INH,RIF");
        // A grade 4 marker must not weaken a grade 1 one for the same drug.
        assert_eq!(per_drug["INH"], Verdict::Resistant);
    }

    #[test]
    fn uncertain_markers_are_reported_without_claiming_resistance() {
        let mut per_drug = HashMap::new();
        fold_marker(&mut per_drug, "EMB;EMB-R;embB_S565G;3) Uncertain significance;embB;S565G");
        fold_marker(&mut per_drug, "INH;INH-R;ndh_c.168C>T;5) Not assoc w R;ndh;c.168C>T");
        assert_eq!(call(&per_drug), "No resistance markers found; uncertain: EMB");
    }

    /// The column must split back into the pairs it was made of, so it keeps
    /// the drugs and drops the prefixes that run into a grade with spaces.
    #[test]
    fn a_resistance_tally_counts_markers_per_drug() {
        let mut prefixes = HashMap::new();
        for path in [
            "RIF;RIF-R;rpoB_D435F;1) Assoc w R;rpoB;D435F",
            "RIF;RIF-R;rpoB_D435V;1) Assoc w R;rpoB;D435V",
            "INH;INH-R;katG_S315T;1) Assoc w R;katG;S315T",
        ] {
            let parts: Vec<&str> = path.split(';').collect();
            for depth in 1..=parts.len() {
                *prefixes.entry(parts[..depth].join(";")).or_insert(0usize) += 1;
            }
        }
        let per_drug = per_drug_counts(&prefixes);
        assert_eq!(per_drug.len(), 2);
        assert_eq!(per_drug["RIF"], 2);
        assert_eq!(per_drug["INH"], 1);
        assert!(per_drug.keys().all(|drug| !drug.contains(' ')));
    }

    #[test]
    fn a_genome_with_only_ungraded_or_unassociated_markers_claims_nothing() {
        let mut per_drug = HashMap::new();
        fold_marker(&mut per_drug, "INH;INH-R;ndh_c.168C>T;4) Not assoc w R - Interim;ndh;c.168C>T");
        assert_eq!(call(&per_drug), "No resistance markers found");
        assert_eq!(call(&HashMap::new()), "No resistance markers found");
    }
}
