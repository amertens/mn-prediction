gc()
d <- readRDS("data/IPD/Malawi/Malawi_merged_dataset.rds")
cat("Rows:", nrow(d), "\n")
cat("Cols:", ncol(d), "\n\n")

# Search for date/interview columns
all_cols <- colnames(d)
date_cols <- grep("date|month|year|interview|hv00[5-8]|v00[5-8]|doi|dob|intd|intm|inty|cm[cy]",
                  all_cols, ignore.case = TRUE, value = TRUE)
cat("Date-related columns:\n")
cat(paste(date_cols, collapse = "\n"), "\n\n")

for (dc in date_cols[1:min(length(date_cols), 15)]) {
  vals <- d[[dc]]
  vals <- vals[!is.na(vals)]
  if (length(vals) > 0) {
    cat(sprintf("--- %s (class: %s) ---\n", dc, class(vals)[1]))
    if (is.numeric(vals)) {
      cat(sprintf("  Range: %s to %s\n", min(vals), max(vals)))
    }
    tt <- head(sort(table(vals), decreasing = TRUE), 10)
    print(tt)
    cat("\n")
  }
}

# Also check first 20 column names for DHS-standard interview date fields
cat("First 50 column names:\n")
cat(paste(head(all_cols, 50), collapse = "\n"), "\n")
