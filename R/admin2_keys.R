# =============================================================================
# R/admin2_keys.R   [JK-01, 2026-09-15]
#
# CHECKED ADMIN-2 JOINS
#
# The project's Admin-2 tables are keyed by names. R/admin2_key_hygiene.R made
# the (Admin1, Admin2) pair the canonical key after name-only joins fanned
# rows in Malawi; this file adds the three things a join still needed to be
# safe:
#
#   admin2_spine()          the canonical unit list with GADM codes
#                           (metadata/admin2_spine.csv, built by
#                           scripts/covariates/build_admin2_spine.R)
#   join_admin2_v2()        a left/inner join on the pair key that REFUSES to
#                           fan rows, reports unmatched keys, and validates
#                           both sides against the spine when asked
#   admin2_population_v2()  the dashboard population table for one country on
#                           the pair key, country label normalised (the file
#                           spells "Sierra Leone" with a space; three scripts
#                           carried the same one-line fix)
#   admin2_match_v2()       the one fuzzy matcher: exact on a normalised key,
#                           then explicit aliases from a CSV, then Jaro-Winkler
#                           within a threshold, and it always writes a review
#                           file of what it decided
#
# Sourced by targets::tar_source("R"); scripts that only source
# R/protocol_v2.R call source("R/admin2_keys.R") themselves.
# =============================================================================

#' Find a project file upward from the working directory
.a2_project_file <- function(...) {
  d <- getwd()
  for (i in 1:6) { f <- file.path(d, ...); if (file.exists(f)) return(f); d <- dirname(d) }
  NULL
}

#' Normalised key for name matching: lower case, letters and digits only
admin2_kk <- function(x) tolower(gsub("[^a-z0-9]", "", tolower(as.character(x))))

#' The canonical Admin-2 spine (country, Admin1, Admin2, gid_1, gid_2, ...)
admin2_spine <- function() {
  f <- .a2_project_file("metadata", "admin2_spine.csv")
  if (is.null(f)) stop("metadata/admin2_spine.csv not found; run scripts/covariates/build_admin2_spine.R")
  sp <- read.csv(f, stringsAsFactors = FALSE)
  stopifnot(all(c("country", "Admin1", "Admin2", "gid_2") %in% names(sp)), !anyDuplicated(paste(sp$country, sp$Admin1, sp$Admin2)))
  sp
}

#' Normalise the country labels different files use to the spine's spelling
admin2_country_label <- function(x) {
  x <- as.character(x)
  k <- admin2_kk(x)
  map <- c(gambia = "Gambia", thegambia = "Gambia", gmb = "Gambia", ghana = "Ghana", gha = "Ghana",
           malawi = "Malawi", mwi = "Malawi", sierraleone = "SierraLeone", sle = "SierraLeone")
  out <- unname(map[k]); out[is.na(out)] <- x[is.na(out)]
  out
}

#' Join two Admin-2 tables on the pair key without fanning rows.
#'
#' @param x,y data.frames carrying Admin1 and Admin2 (and optionally country)
#' @param how "left" keeps every row of x; "inner" keeps the matched rows
#' @param what label for the log line
#' @param check_spine if TRUE, every key of x and y must be a spine unit of
#'   `country` (x's country column, or the `country` argument)
#' @param country the country when x has no country column
#' @return the joined data.frame; stops if y's keys are not unique (the fan),
#'   if the pair key is missing on either side, or if a key is off the spine
join_admin2_v2 <- function(x, y, how = c("left", "inner"), what = "join", check_spine = FALSE, country = NULL, quiet = FALSE) {
  how <- match.arg(how)
  key <- c("Admin1", "Admin2")
  if (!all(key %in% names(x)) || !all(key %in% names(y)))
    stop(sprintf("[join_admin2_v2] %s: both sides need Admin1 and Admin2 (x: %s | y: %s)", what,
                 paste(intersect(key, names(x)), collapse = ","), paste(intersect(key, names(y)), collapse = ",")))
  by <- key; if ("country" %in% names(x) && "country" %in% names(y)) by <- c("country", key)
  ky <- do.call(paste, c(y[by], sep = "||"))
  if (anyDuplicated(ky)) stop(sprintf("[join_admin2_v2] %s: y has %d duplicated key(s), the join would fan rows: %s", what,
                                      sum(duplicated(ky)), paste(head(unique(ky[duplicated(ky)]), 5), collapse = "; ")))
  kx <- do.call(paste, c(x[by], sep = "||"))
  if (check_spine) {
    sp <- admin2_spine()
    cn <- if ("country" %in% names(x)) admin2_country_label(x$country) else if (!is.null(country)) rep(admin2_country_label(country), nrow(x)) else stop("[join_admin2_v2] check_spine needs a country")
    ks <- paste(sp$country, sp$Admin1, sp$Admin2, sep = "||")
    off <- unique(paste(cn, x$Admin1, x$Admin2, sep = "||")[!paste(cn, x$Admin1, x$Admin2, sep = "||") %in% ks])
    if (length(off)) stop(sprintf("[join_admin2_v2] %s: %d key(s) of x are not spine units: %s", what, length(off), paste(head(off, 5), collapse = "; ")))
  }
  out <- if (how == "left") dplyr::left_join(x, y, by = by) else dplyr::inner_join(x, y, by = by)
  n_un <- sum(!kx %in% ky)
  if (how == "left" && nrow(out) != nrow(x)) stop(sprintf("[join_admin2_v2] %s: row count changed %d -> %d", what, nrow(x), nrow(out)))
  if (!quiet) cat(sprintf("[join_admin2_v2] %s: %d x %d rows on (%s) -> %d rows; %d of x unmatched%s\n", what, nrow(x), nrow(y), paste(by, collapse = ","), nrow(out), n_un,
                          if (n_un) paste0(" (", paste(head(unique(kx[!kx %in% ky]), 3), collapse = "; "), ")") else ""))
  out
}

#' The dashboard population table for one country on the pair key.
#'
#' @param POP the table read from dashboard/data/admin2_population.rds
#' @param cn country as the spine spells it ("SierraLeone")
#' @param col the population column to return as `pop`
#' @return data.frame(Admin1, Admin2, pop) with unique keys, spine-validated
admin2_population_v2 <- function(POP, cn, col) {
  stopifnot(col %in% names(POP), all(c("country", "Admin1", "Admin2") %in% names(POP)))
  p <- POP[admin2_country_label(POP$country) == admin2_country_label(cn), c("Admin1", "Admin2", col)]
  names(p)[3] <- "pop"; p$Admin1 <- as.character(p$Admin1); p$Admin2 <- as.character(p$Admin2)
  if (!nrow(p)) stop("[admin2_population_v2] no population rows for ", cn)
  k <- paste(p$Admin1, p$Admin2)
  if (anyDuplicated(k)) stop("[admin2_population_v2] duplicated units for ", cn, ": ", paste(head(unique(k[duplicated(k)]), 5), collapse = "; "))
  sp <- admin2_spine(); sp <- sp[sp$country == admin2_country_label(cn), ]
  off <- k[!k %in% paste(sp$Admin1, sp$Admin2)]
  if (length(off)) stop("[admin2_population_v2] units not on the spine for ", cn, ": ", paste(head(off, 5), collapse = "; "))
  p
}

#' Match source names to a target vocabulary, with aliases and a review file.
#'
#' Exact match on the normalised key first, then the alias table (source key ->
#' target name; metadata/crosswalks/<aliases>.csv with columns source, target),
#' then Jaro-Winkler within `max_jw` AND ahead of the runner-up by `min_margin`.
#' Every decision is written to `review_csv` (source, target, method, jw, jw2,
#' margin, runner_up) so it can be checked line by line.
#'
#' WHY THE MARGIN (JOIN_REVIEW_2026-09-29). `max_jw` alone is far looser than the
#' separation between genuinely different districts: Ghana has 147 distinct
#' district pairs within 0.15 of each other, Malawi 322. `Ahafo Ano South East`
#' and `Ahafo Ano South West` are 0.024 apart, so a source spelling `Ahafo Ano
#' South Est` sits 0.012 from BOTH and the old code silently took whichever
#' which.min() returned first. Measured margins: the three fuzzy matches this
#' project has ever made scored 0.318, 0.053 and 0.041, while every East/West
#' typo scores exactly 0.000. A threshold of 0.02 accepts all the former and
#' rejects all the latter.
#'
#' An ambiguous source is left UNMATCHED (`method = "ambiguous"`) rather than
#' guessed. The review file names the runner-up, so the fix is one line in the
#' alias CSV. `min_margin = 0` restores the old behaviour.
#'
#' @param max_jw largest accepted Jaro-Winkler distance to the best target
#' @param min_margin how far the runner-up must sit beyond the best match
#' @return character vector of target names (NA where unmatched or ambiguous)
admin2_match_v2 <- function(src, tgt, aliases_csv = NULL, max_jw = 0.15, min_margin = 0.02, review_csv = NULL, label = "match") {
  src <- as.character(src); tgt <- unique(as.character(tgt))
  ks <- admin2_kk(src); kt <- admin2_kk(tgt)
  out <- rep(NA_character_, length(src)); method <- rep(NA_character_, length(src)); jw <- rep(NA_real_, length(src))
  hit <- ks %in% kt; out[hit] <- tgt[match(ks[hit], kt)]; method[hit] <- "exact"
  if (!is.null(aliases_csv) && file.exists(aliases_csv)) {
    al <- read.csv(aliases_csv, stringsAsFactors = FALSE); stopifnot(all(c("source", "target") %in% names(al)))
    a <- stats::setNames(al$target, admin2_kk(al$source))
    h <- is.na(out) & ks %in% names(a); out[h] <- unname(a[ks[h]]); method[h] <- "alias"
    bad <- unique(out[h][!out[h] %in% tgt]); if (length(bad)) stop("[admin2_match_v2] alias targets not in the target vocabulary: ", paste(bad, collapse = "; "))
  }
  jw2 <- rep(NA_real_, length(src)); runner <- rep(NA_character_, length(src))
  bestc <- rep(NA_character_, length(src))   # best candidate, recorded even when rejected
  i <- which(is.na(out) & !is.na(ks) & nzchar(ks))   # NA / empty sources stay unmatched
  if (length(i) && length(tgt)) {
    dm <- stringdist::stringdistmatrix(ks[i], kt, method = "jw", p = 0.1)
    dm <- matrix(dm, nrow = length(i))
    # best and runner-up per source; with a single target there is no runner-up,
    # so the margin is infinite and only max_jw binds.
    best <- vapply(seq_len(nrow(dm)), function(r) { z <- dm[r, ]; if (all(is.na(z))) NA_integer_ else which.min(z) }, 1L)
    d    <- vapply(seq_len(nrow(dm)), function(r) { z <- dm[r, ]; if (all(is.na(z))) NA_real_ else min(z, na.rm = TRUE) }, 0)
    sec  <- vapply(seq_len(nrow(dm)), function(r) { z <- sort(dm[r, ], na.last = NA); if (length(z) < 2) Inf else z[2] }, 0)
    sec2 <- vapply(seq_len(nrow(dm)), function(r) { z <- dm[r, ]; if (all(is.na(z)) || length(z) < 2) NA_integer_ else order(z, na.last = NA)[2] }, 1L)
    near  <- !is.na(d) & d <= max_jw
    clear <- near & (sec - d) >= min_margin
    amb   <- near & !clear
    out[i[clear]] <- tgt[best[clear]]; method[i[clear]] <- "fuzzy"
    method[i[amb]] <- "ambiguous"      # left unmatched on purpose; see the docs above
    jw[i] <- round(d, 3); jw2[i] <- round(sec, 3)
    bestc[i[!is.na(best)]] <- tgt[best[!is.na(best)]]
    runner[i[!is.na(sec2)]] <- tgt[sec2[!is.na(sec2)]]
  }
  rev <- data.frame(source = src, target = out, method = method, jw = jw,
                    best_candidate = bestc, jw2 = jw2, margin = round(jw2 - jw, 3),
                    runner_up = runner, stringsAsFactors = FALSE)
  rev <- rev[!duplicated(rev$source), ]
  if (!is.null(review_csv)) { dir.create(dirname(review_csv), showWarnings = FALSE, recursive = TRUE); write.csv(rev[order(is.na(rev$target), rev$method, rev$source), ], review_csv, row.names = FALSE) }
  n_amb <- sum(rev$method %in% "ambiguous")
  cat(sprintf("[admin2_match_v2] %-24s %3d exact, %3d alias, %3d fuzzy, %3d ambiguous, %3d unmatched of %3d%s\n", label, sum(rev$method %in% "exact"), sum(rev$method %in% "alias"),
              sum(rev$method %in% "fuzzy"), n_amb, sum(is.na(rev$target)), nrow(rev), if (!is.null(review_csv)) paste0(" -> ", review_csv) else ""))
  if (n_amb) {
    a <- rev[rev$method %in% "ambiguous", ]
    warning(sprintf("[admin2_match_v2] %s: %d source name(s) too close to call and left unmatched; add them to the alias CSV. %s",
                    label, n_amb, paste(sprintf("%s (%s %.3f vs %s %.3f)", a$source, a$best_candidate, a$jw, a$runner_up, a$jw2)[seq_len(min(3, nrow(a)))], collapse = "; ")), call. = FALSE)
  }
  out
}
