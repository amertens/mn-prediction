# =============================================================================
# R/mrp.R
#
# Multilevel regression and poststratification (MRP) for admin-2 prevalence.
#
# WHY THIS COMPARATOR EXISTS
# --------------------------
# MRP is the reference small-area method a survey-statistics reviewer expects to
# see, and its absence from a 21-method benchmark table is conspicuous. It is
# also the individual-level analogue of the area-level models already compared
# here (Fay-Herriot, BYM2): rather than modelling the district's direct estimate
# and its sampling variance, it models INDIVIDUALS and rebuilds the district
# from a poststratification table.
#
# It additionally addresses a bias the manuscript already concedes: aggregating
# shrunken person-level predicted probabilities up to a district mean is not an
# unbiased district estimator. Poststratifying to a frame is the principled fix.
#
# HOW IT WORKS HERE
#   1. fit  Pr(deficient) ~ individual poststratification variables
#                         + a small screened set of area covariates
#                         + (1 | Admin1) + (1 | Admin1:Admin2)
#      unweighted, which is standard MRP practice: the design information enters
#      through the poststratification table, not through the likelihood.
#   2. build the poststratification table: survey-weighted cell counts per
#      district over the same individual variables.
#   3. predict every cell in every district and take the cell-count-weighted
#      mean.
#
# THE FRAME, AND ITS LIMITATION -- READ THIS BEFORE QUOTING THE NUMBERS
# --------------------------------------------------------------------
# Textbook MRP poststratifies to a CENSUS frame. No census microdata is
# available for these countries in this project, so the frame is built from the
# survey's own design weights: within each district, the weighted distribution
# of the poststratification cells estimates the population distribution.
#
# That is the standard fallback when no census frame exists, and it is not
# vacuous -- the biomarker subsample is not compositionally representative of
# the district that the full survey describes, and this corrects for that. But
# it cannot correct for coverage error in the survey frame itself, and it makes
# MRP here a model-assisted smoother of the design-weighted estimate rather than
# a fully independent estimator. Reported as such.
#
# For a district with no survey data, the cell distribution falls back to its
# Admin-1 parent (then to the national distribution), and the unobserved
# district random effect is set to its prior mean of zero -- so unsurveyed
# districts are separated from one another only by their area covariates. That
# is the honest behaviour of MRP out of sample, not a defect of this
# implementation.
# =============================================================================

suppressPackageStartupMessages({library(dplyr)})

# Poststratification variables are discovered by PATTERN, one per conceptual
# dimension, because the individual-level survey vocabulary drifts between
# countries exactly as the covariate vocabulary does:
#
#   Gambia       gw_urban              gw_hWealthquintile   gw_cSex
#   Ghana        "gw_urban/ rural"     gw_hWealthquintile   gw_cSex
#   Sierra Leone -                     gw_hWealthQuintile   gw_cSex
#   Tanzania     -                     -                    gw_sex
#
# A hardcoded name list silently reduced Ghana and Sierra Leone to no usable
# variables at all. Patterns fix that; the per-country resolution is logged so
# it is never silent.
MRP_PS_PATTERNS <- c(
  residence = "^gw_(urban|rural|urban/ ?rural)$",
  wealth    = "^gw_h?[Ww]ealth[Qq]uintile$",
  sex       = "^gw_(cSex|sex|wSex|cSex01)$",
  child_age = "^gw_cAgeMonthsCat(IYCF)?$",
  woman_age = "^gw_wAgeYearsCat$",
  education = "^gw_(wEducCat|wEducLevel|hEducLevel)$"
)

#' Resolve one usable variable per poststratification dimension.
#'
#' @param d individual data, already restricted to rows with a usable outcome
#' @param max_levels reject anything with more levels than this -- a
#'   high-cardinality variable shatters the poststratification table into cells
#'   with one observation each
#' @param max_missing reject anything missing in more than this fraction of rows
#' @return named character vector: dimension -> column name
.mrp_usable_ps_vars <- function(d, patterns = MRP_PS_PATTERNS,
                                max_levels = 8L, max_missing = 0.4) {
  keep <- character()
  for (dim in names(patterns)) {
    hits <- grep(patterns[[dim]], names(d), value = TRUE, perl = TRUE)
    best <- NULL; best_na <- Inf
    for (v in hits) {
      x <- d[[v]]
      n_na <- sum(is.na(x))
      if (n_na / max(length(x), 1) > max_missing) next
      lv <- unique(x[!is.na(x)])
      if (length(lv) < 2L || length(lv) > max_levels) next
      if (n_na < best_na) { best <- v; best_na <- n_na }
    }
    if (!is.null(best)) keep[dim] <- best
  }
  keep
}

#' Screen area covariates down to a handful by |correlation| with the district
#' survey prevalence. With 14-87 areas, anything more than a few area terms in
#' the fixed part overfits immediately.
.mrp_screen_area_covs <- function(area_cov, svy_admin2, k = 3L) {
  covs <- setdiff(names(area_cov), c("Admin1", "Admin2"))
  if (!length(covs) || is.null(svy_admin2) || !nrow(svy_admin2)) return(character())
  by <- if (all(c("Admin1", "Admin2") %in% names(area_cov)) &&
            all(c("Admin1", "Admin2") %in% names(svy_admin2))) c("Admin1", "Admin2") else "Admin2"
  m <- dplyr::inner_join(svy_admin2[, c(by, "svy_prev")], area_cov, by = by)
  if (nrow(m) < 8L) return(character())
  r <- vapply(covs, function(v) {
    x <- suppressWarnings(as.numeric(m[[v]]))
    if (sum(is.finite(x)) < 8L || stats::sd(x, na.rm = TRUE) == 0) return(0)
    abs(suppressWarnings(stats::cor(x, m$svy_prev, use = "complete.obs")))
  }, numeric(1))
  r[!is.finite(r)] <- 0
  names(sort(r, decreasing = TRUE))[seq_len(min(k, sum(r > 0)))]
}

#' Fit MRP and return admin-2 prevalence for every area in `area_cov`.
#'
#' @param outcome_data output of build_outcome_dataset() (list with $data)
#' @param cc,oc country and outcome configs
#' @param area_cov admin-2 covariate table (Admin1/Admin2 + covariates); supplies
#'   both the set of districts to predict and the area-level fixed effects
#' @param svy_admin2 district survey prevalences, used only to screen area
#'   covariates
#' @param k_area_cov number of area covariates admitted to the fixed part
#' @return data.frame(Admin1, Admin2, mrp_prev, n_cells) or NULL
fit_mrp_admin2 <- function(outcome_data, cc, oc, area_cov, svy_admin2 = NULL,
                           k_area_cov = 3L) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    cat("[mrp] lme4 not installed - skipping\n"); return(NULL)
  }
  d <- outcome_data$data
  ycol <- oc$binary
  if (is.null(ycol) || !ycol %in% names(d)) {
    cat(sprintf("[mrp] outcome column '%s' absent - skipping\n", ycol)); return(NULL)
  }
  # Resolve the SAME uniform outcome compute_svy_admin2() uses. Without this,
  # MRP models the survey's native binary while the district target is the
  # WHO-cutoff-derived one, and the two are different quantities.
  derived <- resolve_uniform_outcome(d, cc, oc, label = "[mrp]")
  if (!is.null(derived)) d[[ycol]] <- derived
  a2 <- cc$admin2_col %||% "Admin2"
  if (!a2 %in% names(d)) { cat("[mrp] no Admin2 column - skipping\n"); return(NULL) }

  d$.y <- suppressWarnings(as.numeric(d[[ycol]]))
  d$.y[!d$.y %in% c(0, 1)] <- NA_real_
  d$.a2 <- trimws(as.character(d[[a2]]))
  d$.a1 <- if ("Admin1" %in% names(d)) trimws(as.character(d$Admin1)) else "ALL"

  wt <- cc$weight_col
  d$.wt <- if (!is.null(wt) && wt %in% names(d)) suppressWarnings(as.numeric(d[[wt]])) else 1
  d$.wt[!is.finite(d$.wt) | d$.wt <= 0] <- 1

  # Resolve on the rows that can actually be fitted, so a variable that is only
  # populated outside the outcome's subpopulation is not counted as usable.
  usable <- !is.na(d$.y) & !is.na(d$.a2)
  ps_map <- .mrp_usable_ps_vars(d[usable, , drop = FALSE])
  if (!length(ps_map)) {
    cat(sprintf(paste0("[mrp] %s - %s: no usable poststratification variables ",
                       "(the analytic table carries %d gw_ columns) - skipping\n"),
                cc$country, oc$tag, sum(grepl("^gw_", names(d)))))
    return(NULL)
  }
  # Copy into sanitised .ps* columns: source names include a slash and a space
  # ("gw_urban/ rural"), which no formula should have to quote.
  ps <- paste0(".ps_", names(ps_map))
  for (i in seq_along(ps_map))
    d[[ps[i]]] <- factor(ifelse(is.na(d[[ps_map[i]]]), "(missing)",
                                as.character(d[[ps_map[i]]])))

  fit_rows <- d[usable, , drop = FALSE]
  if (nrow(fit_rows) < 50L || length(unique(fit_rows$.a2)) < 5L) {
    cat("[mrp] too few observations or districts - skipping\n"); return(NULL)
  }

  # ── area covariates: screen on the districts we can see ──────────────────
  acov <- .mrp_screen_area_covs(area_cov, svy_admin2, k = k_area_cov)
  area_std <- NULL
  if (length(acov)) {
    area_std <- area_cov[, c(intersect(c("Admin1", "Admin2"), names(area_cov)), acov),
                         drop = FALSE]
    for (v in acov) {
      x <- suppressWarnings(as.numeric(area_std[[v]]))
      mu <- mean(x, na.rm = TRUE); sdv <- stats::sd(x, na.rm = TRUE)
      area_std[[v]] <- if (is.finite(sdv) && sdv > 0) (x - mu) / sdv else 0
      area_std[[v]][!is.finite(area_std[[v]])] <- 0
    }
    key_a <- if ("Admin1" %in% names(area_std)) paste(area_std$Admin1, area_std$Admin2, sep = "\u001f") else area_std$Admin2
    key_f <- if ("Admin1" %in% names(area_std)) paste(fit_rows$.a1, fit_rows$.a2, sep = "\u001f") else fit_rows$.a2
    idx <- match(key_f, key_a)
    for (v in acov) fit_rows[[v]] <- area_std[[v]][idx]
    acov <- acov[vapply(acov, function(v) sum(is.finite(fit_rows[[v]])) > 0, logical(1))]
  }

  rhs <- c(ps, acov)
  re  <- if (length(unique(fit_rows$.a1)) > 2L) "(1 | .a1) + (1 | .a1:.a2)" else "(1 | .a2)"
  form <- stats::as.formula(paste(".y ~", paste(c(rhs, re), collapse = " + ")))

  fit <- tryCatch(
    lme4::glmer(form, data = fit_rows, family = stats::binomial(),
                control = lme4::glmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 2e5))),
    error = function(e) { cat("[mrp] glmer failed: ", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(NULL)

  # ── poststratification table: survey-weighted cell counts per district ────
  cells <- d[!is.na(d$.a2), , drop = FALSE] %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(".a1", ".a2", ps)))) %>%
    dplyr::summarise(n_cell = sum(.wt), .groups = "drop")

  # National and Admin-1 fallback distributions for districts with no survey
  # data. Without these, an unsurveyed district has no cells and no estimate.
  nat <- cells %>% dplyr::group_by(dplyr::across(dplyr::all_of(ps))) %>%
    dplyr::summarise(n_cell = sum(n_cell), .groups = "drop")
  a1_cells <- cells %>% dplyr::group_by(dplyr::across(dplyr::all_of(c(".a1", ps)))) %>%
    dplyr::summarise(n_cell = sum(n_cell), .groups = "drop")

  target <- area_cov[, intersect(c("Admin1", "Admin2"), names(area_cov)), drop = FALSE]
  if (!"Admin1" %in% names(target)) target$Admin1 <- "ALL"
  target <- target[!duplicated(paste(target$Admin1, target$Admin2, sep = "\u001f")), , drop = FALSE]

  KEY_SEP <- "\u001f"   # see R/covariates/area_vocabulary.R
  have <- unique(paste(cells$.a1, cells$.a2, sep = KEY_SEP))
  out <- vector("list", nrow(target))

  for (i in seq_len(nrow(target))) {
    a1 <- trimws(as.character(target$Admin1[i])); a2 <- trimws(as.character(target$Admin2[i]))
    cl <- if (paste(a1, a2, sep = KEY_SEP) %in% have) {
      cells[cells$.a1 == a1 & cells$.a2 == a2, , drop = FALSE]
    } else if (a1 %in% a1_cells$.a1) {
      x <- a1_cells[a1_cells$.a1 == a1, , drop = FALSE]; x$.a2 <- a2; x
    } else {
      x <- nat; x$.a1 <- a1; x$.a2 <- a2; x
    }
    if (!nrow(cl)) { out[[i]] <- NULL; next }
    cl$.a1 <- a1; cl$.a2 <- a2
    if (length(acov)) {
      j <- match(paste(a1, a2, sep = KEY_SEP),
                 paste(area_std$Admin1 %||% "ALL", area_std$Admin2, sep = KEY_SEP))
      for (v in acov) cl[[v]] <- if (is.na(j)) 0 else area_std[[v]][j]
    }
    p <- tryCatch(
      stats::predict(fit, newdata = as.data.frame(cl), type = "response",
                     allow.new.levels = TRUE),
      error = function(e) rep(NA_real_, nrow(cl)))
    ok <- is.finite(p) & is.finite(cl$n_cell)
    out[[i]] <- data.frame(
      Admin1 = a1, Admin2 = a2,
      mrp_prev = if (any(ok)) sum(p[ok] * cl$n_cell[ok]) / sum(cl$n_cell[ok]) else NA_real_,
      n_cells = sum(ok), stringsAsFactors = FALSE)
  }

  res <- dplyr::bind_rows(out)
  cat(sprintf("[mrp] %s - %s: %d PS vars (%s), %d area covs (%s), %d districts predicted\n",
              cc$country, oc$tag, length(ps_map),
              paste(sprintf("%s=%s", names(ps_map), ps_map), collapse = ", "),
              length(acov), if (length(acov)) paste(acov, collapse = "/") else "none",
              sum(is.finite(res$mrp_prev))))
  attr(res, "ps_vars") <- ps_map
  attr(res, "area_covs") <- acov
  res
}
