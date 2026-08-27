# =============================================================================
# sandbox_parsimony/R/00_assemble.R
#
# Assemble a pooled Admin-2 dataset straight from the cached `_targets_full`
# store (no pipeline re-run). One row per country x Admin-2, per outcome:
#   svy_prev, svy_prev_se, n_svy  +  all proxy covariates (gee_/ihme_/MAP_/fsec_)
#
# Cached to sandbox_parsimony/out/pooled_<outcome>.rds
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})

STORE <- "_targets_full/objects"
COUNTRIES <- c("gambia", "ghana", "sierraleone", "malawi", "tanzania")

.rd <- function(nm) {
  p <- file.path(STORE, nm)
  if (!file.exists(p)) return(NULL)
  tryCatch(readRDS(p), error = function(e) NULL)
}
.num <- function(v) {
  if (inherits(v, "haven_labelled")) return(as.double(unclass(v)))
  if (is.factor(v) || is.character(v)) return(suppressWarnings(as.numeric(as.character(v))))
  suppressWarnings(as.numeric(v))
}

# Same year-stripping harmonization the production code uses, so multi-year
# layers collapse to one stem and the cross-country common set is not ~13 vars.
harmonize_names <- function(x) {
  x <- gsub("_(19|20)[0-9]{2}([0-9]{4})?", "", x)
  x <- gsub("__+", "_", x)
  sub("_$", "", x)
}

#' Country-level Admin-2 covariate table (GEE polygons + merged cluster means)
build_cov <- function(ctry) {
  ac <- .rd(paste0("area_covariates_", ctry))
  gee <- if (!is.null(ac)) ac$gee_admin2 else .rd(paste0("gee_admin2_", ctry))
  if (is.null(gee)) return(NULL)

  gee_cols <- grep("^gee_", names(gee), value = TRUE)
  stems <- harmonize_names(gee_cols)
  gmat <- vapply(gee_cols, function(c) .num(gee[[c]]), numeric(nrow(gee)))
  ustem <- unique(stems)
  gee_df <- as.data.frame(vapply(ustem, function(s) {
    cols <- which(stems == s)
    if (length(cols) == 1) gmat[, cols] else rowMeans(gmat[, cols, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(gee))))
  names(gee_df) <- ustem
  gee_df$Admin2 <- as.character(gee$Admin2)

  merged <- .rd(paste0("merged_", ctry))
  if (!is.null(merged) && "Admin2" %in% names(merged)) {
    dom <- grep("^(ihme_|MAP_|fsec_)", names(merged), value = TRUE)
    dom <- dom[!is.na(dom)]
    agg <- data.frame(Admin2 = as.character(merged$Admin2), stringsAsFactors = FALSE)
    for (c in dom) agg[[c]] <- .num(merged[[c]])
    agg <- agg |> group_by(Admin2) |>
      summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
      as.data.frame()
    gee_df <- dplyr::full_join(gee_df, agg, by = "Admin2")
  }

  # Admin-2 centroids for spatial terms / neighbour structures
  ctr <- tryCatch({
    cen <- .rd("admin2_centroids")
    NULL
  }, error = function(e) NULL)
  if (!is.null(ac) && !is.null(ac$polygons)) {
    pg <- ac$polygons
    xy <- tryCatch(sf::st_coordinates(sf::st_centroid(sf::st_geometry(pg))),
                   error = function(e) NULL)
    if (!is.null(xy)) {
      nm2 <- if ("Admin2" %in% names(pg)) as.character(pg$Admin2) else
             if ("NAME_2" %in% names(pg)) as.character(pg$NAME_2) else NULL
      nm1 <- if ("Admin1" %in% names(pg)) as.character(pg$Admin1) else
             if ("NAME_1" %in% names(pg)) as.character(pg$NAME_1) else NA_character_
      if (!is.null(nm2))
        gee_df <- dplyr::left_join(
          gee_df,
          data.frame(Admin2 = nm2, Admin1 = nm1, lon = xy[, 1], lat = xy[, 2],
                     stringsAsFactors = FALSE) |> distinct(Admin2, .keep_all = TRUE),
          by = "Admin2")
    }
  }
  gee_df
}

#' Pool one outcome across every country that has it
assemble_outcome <- function(outcome, cov_list) {
  frames <- list()
  for (ctry in names(cov_list)) {
    s <- .rd(paste0("svy_admin2_", ctry, "_", outcome))
    if (is.null(s) || !any(is.finite(s$svy_prev))) next
    sv <- data.frame(
      Admin2      = as.character(s$Admin2),
      svy_prev    = as.numeric(s$svy_prev),
      svy_prev_se = if ("svy_prev_se" %in% names(s)) as.numeric(s$svy_prev_se) else NA_real_,
      n_svy       = if ("n_svy" %in% names(s)) as.numeric(s$n_svy) else NA_real_,
      stringsAsFactors = FALSE)
    sv <- sv[is.finite(sv$svy_prev), , drop = FALSE]
    d <- dplyr::inner_join(sv, cov_list[[ctry]], by = "Admin2")
    if (nrow(d) < 4) next
    d$country <- ctry
    frames[[ctry]] <- d
  }
  if (length(frames) < 1) return(NULL)

  covnames <- function(df) {
    cv <- setdiff(names(df), c("Admin2", "Admin1", "svy_prev", "svy_prev_se",
                               "n_svy", "country"))
    cv[vapply(cv, function(c) any(is.finite(df[[c]])), logical(1))]
  }
  common <- Reduce(intersect, lapply(frames, covnames))
  keep <- c("country", "Admin2", "Admin1", "svy_prev", "svy_prev_se", "n_svy", common)
  pooled <- bind_rows(lapply(frames, function(df) {
    for (k in setdiff(keep, names(df))) df[[k]] <- NA
    df[, keep, drop = FALSE]
  }))
  list(data = as.data.frame(pooled), predictors = setdiff(common, c("lon", "lat")),
       countries = names(frames), outcome = outcome)
}

# Set .ASSEMBLE_FNS_ONLY <- TRUE before sourcing to load the helpers without
# rebuilding the cached pooled dataset.
if (!isTRUE(get0(".ASSEMBLE_FNS_ONLY", ifnotfound = FALSE))) {
  message("Building per-country covariate tables ...")
  cov_list <- list()
  for (ctry in COUNTRIES) {
    cv <- build_cov(ctry)
    if (!is.null(cv)) {
      cov_list[[ctry]] <- cv
      message(sprintf("  %-12s %d areas x %d covariates", ctry, nrow(cv), ncol(cv) - 1))
    }
  }
  saveRDS(cov_list, "sandbox_parsimony/out/cov_list.rds")

  outcomes <- c("child_vitA", "women_vitA", "child_iron", "women_iron",
                "women_folate", "women_b12", "child_zinc", "women_zinc")
  pooled_all <- list()
  for (oc in outcomes) {
    p <- assemble_outcome(oc, cov_list)
    if (is.null(p)) { message(sprintf("  %-14s (skipped)", oc)); next }
    pooled_all[[oc]] <- p
    message(sprintf("  %-14s n=%3d areas, %d countries, %d common predictors",
                    oc, nrow(p$data), length(p$countries), length(p$predictors)))
  }
  saveRDS(pooled_all, "sandbox_parsimony/out/pooled_all.rds")
  message("Saved -> sandbox_parsimony/out/pooled_all.rds")
}
