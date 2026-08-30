# =============================================================================
# sandbox_parsimony/R/17_admin1_layer_fns.R
# Level-aggregation helpers shared by 17_admin1_layer.R and
# 18_mixed_level_check.R.
# =============================================================================

WATER <- "^Lake |^Water|Reservoir|^Lac |^Sea$"
STORE <- "_targets_full/objects"

covariates_admin1 <- function(d) {
  if (!"Admin1" %in% names(d)) return(NULL)
  wcol <- intersect(c("gee_ghspop", "gee_gpw_demographic_annual_mean",
                      "gee_popdensity"), names(d))
  w <- if (length(wcol)) suppressWarnings(as.numeric(d[[wcol[1]]])) else rep(1, nrow(d))
  w[!is.finite(w) | w < 0] <- stats::median(w[is.finite(w) & w > 0], na.rm = TRUE)
  if (!any(is.finite(w)) || all(w == 0)) w <- rep(1, nrow(d))
  d$.w <- w
  num <- names(d)[vapply(d, is.numeric, logical(1))]
  num <- setdiff(num, ".w")
  d |> filter(!is.na(Admin1)) |> group_by(Admin1) |>
    summarise(across(all_of(num),
                     ~ if (all(!is.finite(.x))) NA_real_ else
                       stats::weighted.mean(.x, .w * is.finite(.x), na.rm = TRUE)),
              n_admin2 = n(), .groups = "drop") |>
    rename(Admin2 = Admin1) |>            # reuse the Admin2 key downstream
    as.data.frame()
}

# ---------------------------------------------------------------------------
# Survey outcome at a given level, from the same outcome_data_* targets the
# production pipeline uses.
# ---------------------------------------------------------------------------
outcome_at <- function(ctry, oc, cc, level = c("Admin1", "Admin2")) {
  level <- match.arg(level)
  f <- file.path(STORE, paste0("outcome_data_", ctry, "_", oc$tag))
  if (!file.exists(f)) return(NULL)
  od <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(od)) return(NULL)
  d <- od$data
  gcol <- if (level == "Admin1") cc$admin1_col else cc$admin2_col
  if (!all(c(gcol, cc$weight_col, oc$binary) %in% names(d))) return(NULL)
  g <- as.character(d[[gcol]])
  y <- suppressWarnings(as.numeric(d[[oc$binary]]))
  w <- suppressWarnings(as.numeric(d[[cc$weight_col]])); w[!is.finite(w) | w <= 0] <- 1
  ok <- !is.na(g) & is.finite(y)
  if (!sum(ok)) return(NULL)
  out <- data.frame(Admin2 = g[ok], y = y[ok], w = w[ok]) |>
    group_by(Admin2) |>
    summarise(svy_prev = stats::weighted.mean(y, w), n_svy = n(), .groups = "drop") |>
    as.data.frame()
  # never model a lake
  out[!grepl(WATER, out$Admin2, ignore.case = TRUE), , drop = FALSE]
}

