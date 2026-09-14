#!/usr/bin/env Rscript
# Render the country briefs — run from the repo root:
#   Rscript dashboard/report/render_reports.R            # all countries
#   Rscript dashboard/report/render_reports.R ghana      # one
#
# Emits, per country, into dashboard/report/out/:
#   <country>_brief.pdf    for external reviewers, and readable page-by-page
#   <country>_brief.html   self-contained; keeps nothing interactive but travels
#   <country>_brief.md     the ~4k-token text version (keep-md), cheap to read
#   <country>_districts.csv, <country>_national.csv
#
# The CSVs exist so a number can be checked against the source rather than read
# off a chart, and because a 2 MB self-contained HTML is ~99% machinery: of
# docs/admin2_predictions_report.html's 2.1M characters, 14k were prose.

suppressPackageStartupMessages({library(here)})

# quarto's R wrapper reports every failure as "Error running quarto CLI from R"
# and then breaks while formatting it ("Could not evaluate cli {} expression:
# `captions`"), so the actual cause never reaches the console. The real message
# is in the xelatex .log. Print the interesting lines of it.
show_latex_log <- function(stem) {
  lg <- here::here("dashboard", "report", paste0(stem, ".log"))
  if (!file.exists(lg)) { cat("    (no .log at", lg, ")\n"); return(invisible()) }
  ln <- readLines(lg, warn = FALSE)
  hit <- grep("^!|LaTeX Error|Undefined control|Emergency stop|File ended while", ln)
  if (!length(hit)) { cat("    (.log has no error lines)\n"); return(invisible()) }
  cat("    --- from", basename(lg), "---\n")
  for (i in head(hit, 3)) cat(paste0("    ", ln[seq(i, min(i + 4, length(ln)))]), sep = "\n")
  cat("    ---\n")
}

# Quarto ships inside RStudio and is not on PATH on this machine.
if (!nzchar(Sys.getenv("QUARTO_PATH"))) {
  cand <- c(Sys.which("quarto"),
            "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe",
            file.path(Sys.getenv("LOCALAPPDATA"), "Programs/Quarto/bin/quarto.exe"),
            "C:/Program Files/Quarto/bin/quarto.exe")
  cand <- cand[nzchar(cand) & file.exists(cand)]
  if (!length(cand)) stop("Quarto CLI not found. Install it, or set QUARTO_PATH.")
  Sys.setenv(QUARTO_PATH = cand[1])
}
cat("Quarto:", Sys.getenv("QUARTO_PATH"), "\n")

QMD <- here::here("dashboard", "report", "country_brief.qmd")
OUT <- here::here("dashboard", "report", "out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
countries <- if (length(args)) args else c("gambia", "ghana", "sierraleone", "malawi")

# ── Companion CSVs, straight from the same tables the brief prints ──────────
source(here::here("dashboard", "report", "write_brief_csvs.R"))
write_csvs <- function(ck) write_brief_csvs(ck, OUT)

failed <- character(0)
for (ck in countries) {
  cat(sprintf("\n== %s ==\n", ck))
  # Every country renders through the same country_brief.tex / .pdf /
  # country_brief_files paths, so a slow cleanup from the previous one can leave
  # a stale intermediate and xelatex fails on the next. Seen once in four.
  # Worse, quarto's own error formatter breaks reporting it ("object 'captions'
  # not found"), so the real cause never surfaces. Clear the intermediates and
  # retry once before giving up.
  render_once <- function() {
    for (f in list.files(dirname(QMD), pattern = "^country_brief\\.(tex|pdf|html|log|aux|md)$",
                         full.names = TRUE)) unlink(f)
    unlink(file.path(dirname(QMD), "country_brief_files"), recursive = TRUE)
    quarto::quarto_render(input = QMD, output_format = "all",
                          execute_params = list(country = ck), quiet = FALSE)
  }
  ok <- tryCatch({ render_once(); TRUE },
                 error = function(e) {
                   cat("  render failed:", conditionMessage(e), "\n")
                   show_latex_log("country_brief")
                   cat("  retrying once...\n")
                   tryCatch({ Sys.sleep(2); render_once(); TRUE },
                            error = function(e2) {
                              cat("  retry failed:", conditionMessage(e2), "\n")
                              show_latex_log("country_brief"); FALSE })
                 })
  if (!ok) { failed <- c(failed, ck); next }

  # quarto writes beside the .qmd; move the artefacts under out/ with a
  # country-stamped name so four renders do not overwrite each other.
  # keep-md writes country_brief.pdf.md / .html.md, not country_brief.md. Keep
  # the pdf-flavoured one: it is the plain-text rendering, ~13 KB against the
  # 2.2 MB self-contained HTML, and it is what makes the brief cheap to read
  # back later without parsing base64.
  moves <- c("country_brief.pdf" = "pdf", "country_brief.html" = "html",
             "country_brief.pdf.md" = "md")
  for (nm in names(moves)) {
    src <- here::here("dashboard", "report", nm)
    if (file.exists(src)) {
      dst <- file.path(OUT, sprintf("%s_brief.%s", ck, moves[[nm]]))
      file.rename(src, dst)
      cat(sprintf("    %-4s %s (%.0f KB)\n", moves[[nm]], basename(dst),
                  file.size(dst) / 1024))
    }
  }
  unlink(here::here("dashboard", "report", "country_brief.html.md"))
  unlink(here::here("dashboard", "report", "country_brief_files"), recursive = TRUE)
  tryCatch(write_csvs(ck), error = function(e) cat("    csv failed:", conditionMessage(e), "\n"))
}

# ── The document that is not per-country ───────────────────────────────────
# overview.qmd — all four countries side by side (the portfolio view). The
# technical annex was retired on 2026-09-13 with the tabs it printed.
if (!length(args)) {
  for (doc in c("overview")) {
    cat(sprintf("\n== %s ==\n", doc))
    src_qmd <- here::here("dashboard", "report", paste0(doc, ".qmd"))
    if (!file.exists(src_qmd)) { cat("  missing:", src_qmd, "\n"); next }
    ok <- tryCatch({
      for (f in list.files(dirname(src_qmd),
                           pattern = sprintf("^%s\\.(tex|log|aux)$", doc),
                           full.names = TRUE)) unlink(f)
      quarto::quarto_render(input = src_qmd, output_format = "all", quiet = FALSE)
      TRUE
    }, error = function(e) {
      cat("  render failed:", conditionMessage(e), "\n"); show_latex_log(doc)
      cat("  retrying once...\n")
      tryCatch({ Sys.sleep(2)
                 quarto::quarto_render(input = src_qmd, output_format = "all",
                                       quiet = FALSE)
                 TRUE },
               error = function(e2) {
                 cat("  retry failed:", conditionMessage(e2), "\n")
                 show_latex_log(doc); FALSE })
    })
    if (!ok) { failed <- c(failed, doc); next }
    for (nm in c("pdf", "html", "pdf.md")) {
      src <- here::here("dashboard", "report", sprintf("%s.%s", doc, nm))
      if (file.exists(src)) {
        dst <- file.path(OUT, sprintf("%s.%s", doc, sub("^pdf\\.", "", nm)))
        file.rename(src, dst)
        cat(sprintf("    %-4s %s (%.0f KB)\n", sub("^pdf\\.", "", nm),
                    basename(dst), file.size(dst) / 1024))
      }
    }
    unlink(here::here("dashboard", "report", sprintf("%s.html.md", doc)))
    unlink(here::here("dashboard", "report", sprintf("%s_files", doc)), recursive = TRUE)
  }
}

if (length(failed)) {
  cat("\nFAILED:", paste(failed, collapse = ", "), "\n"); quit(status = 1)
}
cat(sprintf("\nDone. Output in %s\n", OUT))
