# =============================================================================
# R/lint_admin2_joins.R
#
# A LINT FOR THE DEFECT CLASS THAT HAS NOW RECURRED TEN TIMES.
#
# WHY THIS EXISTS
# ---------------
# Malawi has six Admin-2 names that occur in more than one Admin-1 region
# (`Lake Chilwa`, `Lake Malawi`, `TA Lundu`, `TA Malemia`, `TA Ngabu`,
# `TA Pemba`): GADM level 2 holds 256 polygons under 243 distinct names. Joining
# two Admin-2 tables on the district NAME therefore fans rows, and the fan is
# silent. Nine consumers were migrated to the pair key (Admin1, Admin2) earlier
# in the project. A tenth, scripts/covariates/08_admin1_arms.R, was found only
# by comparing a deployed dashboard against the application's own prediction
# table, 90 districts against 87.
#
# Section 8 of docs/SESSION_FINDINGS_FOR_REVIEW.md records the conclusion: the
# codebase should not depend on a developer remembering to call
# admin2_join_by(). This file is the structural answer.
#
# WHY IT IS A RATCHET AND NOT A BAN
# ---------------------------------
# Measured 2026-08-31 on this branch: there are 65 name-only `by = "Admin2"`
# joins across R/, scripts/ and dashboard/. A test that fails on all of them
# fails on first contact and gets switched off within a day, which leaves the
# codebase less protected than before. The test therefore fails on any site NOT
# in a recorded baseline. Existing sites are grandfathered and visible; a new
# one cannot be added without either fixing it or amending the baseline
# deliberately, which is a reviewable act rather than an oversight.
#
# The baseline is keyed on (file, trimmed source line), not on a line number, so
# inserting or deleting lines elsewhere in a file does not invalidate it.
#
# NOT EVERY NAME-ONLY JOIN IS A DEFECT. A join is safe when neither side can
# carry a duplicate name, for example a table already restricted to one country
# whose districts are uniquely named. The baseline records that these sites were
# seen, not that they were cleared. `audit` in the baseline carries the
# assessment where one has been made.
# =============================================================================

#' Directories the lint scans.
#' @keywords internal
ADMIN2_LINT_DIRS <- c("R", "scripts", "dashboard")

#' Files exempt from the lint.
#'
#' R/admin2_key_hygiene.R defines admin2_join_by(), whose documented behaviour is
#' to RETURN "Admin2" when a caller's table has no Admin1 column. That fallback
#' is the sanctioned name-only path and cannot be written any other way.
#' R/lint_admin2_joins.R contains the patterns themselves.
#' @keywords internal
ADMIN2_LINT_EXEMPT <- c(
  file.path("R", "admin2_key_hygiene.R"),
  file.path("R", "lint_admin2_joins.R")
)

#' Patterns that identify a name-only Admin-2 join.
#'
#' Comment lines are stripped before matching, so the many explanatory comments
#' that quote `by = "Admin2"` while describing the defect are not themselves
#' reported as defects.
#' @keywords internal
.admin2_lint_patterns <- function() {
  list(
    # by = "Admin2" / by="Admin2" / by = c("Admin2")
    list(kind = "by_name_only",
         rx   = "\\bby(\\.[xy])?\\s*=\\s*(c\\s*\\(\\s*)?[\"']Admin2[\"']\\s*\\)?"),
    # dplyr positional join: left_join(x, y, "Admin2")
    list(kind = "positional_join",
         rx   = "\\b(inner|left|right|full|semi|anti)_join\\s*\\([^)]*,\\s*[\"']Admin2[\"']"),
    # base positional merge: merge(x, y, "Admin2")
    list(kind = "positional_merge",
         rx   = "\\bmerge\\s*\\([^)]*,\\s*[\"']Admin2[\"']\\s*[,)]")
  )
}

#' Does this line join on a composite key whose other element is not Admin1?
#'
#' `by = c("country", "Admin2")` and `by = c("country_label", "Admin2")` are
#' composite, so they are not bare name joins, and the first version of this
#' lint mislabelled them "positional". They are also not safe: `country` and
#' `country_label` disambiguate Malawi from Ghana, and the Malawi fan is WITHIN
#' Malawi, where six district names each occur in more than one region. They are
#' therefore reported under their own kind rather than exempted, so the audit can
#' record what each one joins against.
#' @keywords internal
.admin2_composite_no_admin1 <- function(code) {
  has_multi <- grepl("c\\s*\\(\\s*[\"'][^\"']+[\"']\\s*,\\s*[\"'][^\"']+[\"']", code)
  has_a1    <- grepl("[\"']Admin1[\"']", code)
  has_multi & !has_a1
}

#' Strip whole-line and trailing comments without breaking quoted strings.
#'
#' A `#` inside a quoted string is not a comment. Walking the line character by
#' character is slower than a regex and correct, which matters because a false
#' strip would hide a real join.
#' @keywords internal
.strip_r_comment <- function(line) {
  chars <- strsplit(line, "", fixed = TRUE)[[1]]
  if (!length(chars)) return(line)
  in_s <- FALSE; in_d <- FALSE; esc <- FALSE
  for (i in seq_along(chars)) {
    ch <- chars[i]
    if (esc) { esc <- FALSE; next }
    if (ch == "\\" && (in_s || in_d)) { esc <- TRUE; next }
    if (ch == "'"  && !in_d) { in_s <- !in_s; next }
    if (ch == '"'  && !in_s) { in_d <- !in_d; next }
    if (ch == "#"  && !in_s && !in_d) {
      return(if (i == 1L) "" else paste(chars[seq_len(i - 1L)], collapse = ""))
    }
  }
  line
}

#' Scan the codebase for name-only Admin-2 joins.
#'
#' @param root repository root; defaults to here::here()
#' @param dirs directories to scan
#' @return data.frame(file, line, kind, code) with `file` relative to root and
#'   forward-slash separated, so the baseline is portable across platforms.
scan_admin2_joins <- function(root = here::here(), dirs = ADMIN2_LINT_DIRS) {
  pats <- .admin2_lint_patterns()
  out <- list()
  for (d in dirs) {
    dd <- file.path(root, d)
    if (!dir.exists(dd)) next
    files <- list.files(dd, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
    for (f in files) {
      rel <- sub("\\\\", "/", sub(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", root), "/?"),
                                 "", f))
      rel <- gsub("\\\\", "/", rel)
      if (rel %in% gsub("\\\\", "/", ADMIN2_LINT_EXEMPT)) next
      txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
      if (!length(txt)) next
      code <- vapply(txt, .strip_r_comment, character(1), USE.NAMES = FALSE)
      # A line that names Admin1 as well as Admin2 is joining on the pair key,
      # which is the correct form. The positional patterns would otherwise match
      # inside `by = c("Admin1", "Admin2")`, because the character class before
      # the token happily spans the Admin1 argument. Caught by this file's own
      # self-test in tests/testthat/test-admin2-join-lint.R.
      is_pair <- grepl("[\"']Admin1[\"']", code)
      is_comp <- .admin2_composite_no_admin1(code)
      for (p in pats) {
        hit <- which(grepl(p$rx, code, perl = TRUE) & !is_pair)
        for (i in hit) {
          out[[length(out) + 1L]] <- data.frame(
            file = rel, line = i,
            kind = if (is_comp[i]) "composite_without_admin1" else p$kind,
            code = trimws(code[i]), stringsAsFactors = FALSE)
        }
      }
    }
  }
  res <- if (length(out)) do.call(rbind, out) else
    data.frame(file = character(0), line = integer(0), kind = character(0),
               code = character(0), stringsAsFactors = FALSE)
  # One row per site: a line matching two patterns is one defect, not two.
  res <- res[!duplicated(res[, c("file", "line")]), , drop = FALSE]
  res[order(res$file, res$line), , drop = FALSE]
}

#' Compare a scan against the recorded baseline.
#'
#' @param scan output of scan_admin2_joins()
#' @param baseline data.frame with at least `file` and `code`
#' @return list(new = <sites absent from the baseline>,
#'              fixed = <baseline sites no longer present>)
diff_admin2_baseline <- function(scan, baseline) {
  key <- function(d) paste(d$file, trimws(d$code), sep = " :: ")
  ks <- key(scan); kb <- key(baseline)
  list(new   = scan[!ks %in% kb, , drop = FALSE],
       fixed = baseline[!kb %in% ks, , drop = FALSE])
}
