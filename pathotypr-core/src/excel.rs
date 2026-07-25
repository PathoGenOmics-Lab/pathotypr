//! Excel export utilities.
//!
//! Provides streaming and batch Excel writing alongside TSV output.

use rust_xlsxwriter::{
    Color, ConditionalFormatCell, ConditionalFormatCellRule, Format, Workbook,
};

// ---------------------------------------------------------------------------
// Path helpers
// ---------------------------------------------------------------------------

/// Derives the `.xlsx` path that the Excel writers use for a given TSV path.
///
/// Cleanup routines must use this too: deriving the path differently (e.g. with
/// `Path::with_extension`) makes them delete the wrong file when the output does
/// not end in `.tsv`, leaving the real partial file behind.
pub(crate) fn excel_path_from_tsv(tsv_path: &str) -> String {
    // Replace only the trailing ".tsv" extension, not every occurrence in the
    // path (a directory name or sample name may legitimately contain ".tsv").
    if let Some(stem) = tsv_path.strip_suffix(".tsv") {
        format!("{}.xlsx", stem)
    } else {
        format!("{}.xlsx", tsv_path)
    }
}

// ---------------------------------------------------------------------------
// Streaming writer
// ---------------------------------------------------------------------------

/// Streaming Excel writer that avoids holding all rows in memory.
pub struct ExcelStreamWriter {
    workbook: Workbook,
    worksheet_index: usize,
    xlsx_path: String,
    data_format: Format,
    max_widths: Vec<usize>,
    next_row: u32,
    /// Deferred conditional formatting rules: (col_idx, low_threshold, high_threshold).
    conditional_formats: Vec<(u16, f64, f64)>,
}

impl ExcelStreamWriter {
    /// Creates a streaming Excel writer and writes the header row.
    pub fn new(
        tsv_path: &str,
        headers: &[&str],
    ) -> Result<Self, Box<dyn std::error::Error + Send + Sync>> {
        let xlsx_path = excel_path_from_tsv(tsv_path);
        let mut workbook = Workbook::new();
        let worksheet_index = 0usize;
        let _ = workbook.add_worksheet();

        let header_format = Format::new()
            .set_bold()
            .set_background_color(Color::RGB(0xE4E4E7))
            .set_border(rust_xlsxwriter::FormatBorder::Thin);
        let data_format = Format::new().set_border(rust_xlsxwriter::FormatBorder::Thin);

        {
            let worksheet = workbook.worksheet_from_index(worksheet_index)?;
            for (col, header) in headers.iter().enumerate() {
                worksheet.write_string_with_format(0, col as u16, *header, &header_format)?;
            }
            worksheet.set_freeze_panes(1, 0)?;
        }

        let max_widths = headers.iter().map(|h| h.len()).collect();

        Ok(Self {
            workbook,
            worksheet_index,
            xlsx_path,
            data_format,
            max_widths,
            next_row: 1,
            conditional_formats: Vec::new(),
        })
    }

    /// Writes a single data row.
    pub fn write_row(
        &mut self,
        row: &[String],
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        let worksheet = self.workbook.worksheet_from_index(self.worksheet_index)?;
        for (col_idx, cell) in row.iter().enumerate() {
            if col_idx >= self.max_widths.len() {
                self.max_widths.push(0);
            }
            self.max_widths[col_idx] = self.max_widths[col_idx].max(cell.len());

            if let Ok(num) = cell.parse::<f64>() {
                worksheet.write_number_with_format(
                    self.next_row,
                    col_idx as u16,
                    num,
                    &self.data_format,
                )?;
            } else {
                worksheet.write_string_with_format(
                    self.next_row,
                    col_idx as u16,
                    cell,
                    &self.data_format,
                )?;
            }
        }
        self.next_row += 1;
        Ok(())
    }

    /// Registers conditional formatting for a numeric column.
    ///
    /// Values below `low` get red, above `high` get green, in between get yellow.
    /// Applied when `finish()` is called.
    pub fn set_conditional_formatting(&mut self, col_idx: u16, low: f64, high: f64) {
        self.conditional_formats.push((col_idx, low, high));
    }

    /// Finalizes and saves the workbook.
    pub fn finish(mut self) -> Result<String, Box<dyn std::error::Error + Send + Sync>> {
        {
            let worksheet = self.workbook.worksheet_from_index(self.worksheet_index)?;

            for (col, width) in self.max_widths.iter().enumerate() {
                worksheet.set_column_width(col as u16, (width + 2).min(50) as f64)?;
            }

            if self.next_row > 1 && !self.max_widths.is_empty() {
                worksheet.autofilter(
                    0,
                    0,
                    self.next_row - 1,
                    (self.max_widths.len() - 1) as u16,
                )?;
            }

            let last_data_row = self.next_row.saturating_sub(1);
            let fmt_red = Format::new().set_background_color(Color::RGB(0xFECDD3));
            let fmt_yellow = Format::new().set_background_color(Color::RGB(0xFEF3C7));
            let fmt_green = Format::new().set_background_color(Color::RGB(0xBBF7D0));

            for &(col, low, high) in &self.conditional_formats {
                if last_data_row < 1 {
                    continue;
                }
                let rule_green = ConditionalFormatCell::new()
                    .set_rule(ConditionalFormatCellRule::GreaterThanOrEqualTo(high))
                    .set_format(&fmt_green);
                let rule_yellow = ConditionalFormatCell::new()
                    .set_rule(ConditionalFormatCellRule::Between(low, high))
                    .set_format(&fmt_yellow);
                let rule_red = ConditionalFormatCell::new()
                    .set_rule(ConditionalFormatCellRule::LessThan(low))
                    .set_format(&fmt_red);

                worksheet.add_conditional_format(1, col, last_data_row, col, &rule_green)?;
                worksheet.add_conditional_format(1, col, last_data_row, col, &rule_yellow)?;
                worksheet.add_conditional_format(1, col, last_data_row, col, &rule_red)?;
            }
        }

        // A failed save can leave a partially written workbook behind. Callers
        // only warn on failure, so without this the run would be reported as
        // successful next to a corrupt .xlsx.
        if let Err(e) = self.workbook.save(&self.xlsx_path) {
            if let Err(rm) = std::fs::remove_file(&self.xlsx_path) {
                if rm.kind() != std::io::ErrorKind::NotFound {
                    log::warn!("Failed to remove partial Excel file {}: {}", self.xlsx_path, rm);
                }
            }
            return Err(e.into());
        }
        Ok(self.xlsx_path)
    }
}

// ---------------------------------------------------------------------------
// Batch writer (convenience)
// ---------------------------------------------------------------------------

/// Writes data to an Excel file alongside the TSV output.
pub fn write_excel_file(
    tsv_path: &str,
    headers: &[&str],
    rows: &[Vec<String>],
) -> Result<String, Box<dyn std::error::Error + Send + Sync>> {
    let mut writer = ExcelStreamWriter::new(tsv_path, headers)?;
    for row in rows {
        writer.write_row(row)?;
    }
    writer.finish()
}

/// Parses TSV content into headers and rows for Excel export.
pub fn parse_tsv_for_excel(tsv_content: &str) -> (Vec<&str>, Vec<Vec<String>>) {
    let mut lines = tsv_content.lines();

    let headers: Vec<&str> = lines
        .next()
        .map(|h| h.split('\t').collect())
        .unwrap_or_default();

    let rows: Vec<Vec<String>> = lines
        .filter(|line| !line.is_empty())
        .map(|line| line.split('\t').map(String::from).collect())
        .collect();

    (headers, rows)
}
