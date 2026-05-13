# =============================================================================
# src/GEE/seasonality_features.R
#
# Derive seasonality summary features from monthly GEE buffer columns.
# Replaces 12 monthly columns per (family × scale) with 5 dense features:
#   *_annual_mean   — mean across months
#   *_seasonal_sd   — sd across months (raw variability)
#   *_seasonal_cv   — coefficient of variation (sd / abs(mean))
#   *_seasonal_range— max - min across months
#   *_peak_month    — index 1..12 of max month (numeric, treated as ordinal)
#
# Existing monthly columns are kept (the model may still use them) but the
# summary columns capture the seasonal *shape*, which the existing learners
# can pick up more directly than via 12 tree-splits.
# =============================================================================

#' Add seasonality features to a data frame.
#'
#' Detects column families of the form `gee_<family>_<Month>_<scale>` (e.g.
#' `gee_trmm_Jan_10km`) by scanning column names. For each detected
#' (family, scale) pair, computes 5 row-wise summary columns.
#'
#' @param df data frame containing monthly buffer columns
#' @param verbose if TRUE, print which families were detected
#' @return df with added `gee_<family>_<scale>_<stat>` columns
add_seasonality_features <- function(df, verbose = TRUE) {

  months <- c("Jan","Feb","Mar","Apr","May","Jun",
              "Jul","Aug","Sep","Oct","Nov","Dec")
  month_pat <- paste(months, collapse = "|")
  # Match e.g. "gee_trmm_Jan_10km" — capture family + scale
  rgx <- sprintf("^gee_([a-z][a-z0-9]*)_(%s)_([0-9]+km)$", month_pat)

  hits <- regmatches(colnames(df), regexec(rgx, colnames(df)))
  hits <- hits[lengths(hits) == 4L]
  if (length(hits) == 0) {
    if (verbose) message("[seasonality] no monthly columns detected")
    return(df)
  }
  fam_scale <- vapply(hits, function(m) paste(m[2], m[4], sep = "__"),
                      character(1))

  for (key in unique(fam_scale)) {
    parts  <- strsplit(key, "__", fixed = TRUE)[[1]]
    family <- parts[1]; scale <- parts[2]

    cols <- sprintf("gee_%s_%s_%s", family, months, scale)
    cols <- cols[cols %in% colnames(df)]
    if (length(cols) < 6) next  # need at least 6 months for a seasonality signal

    M <- as.matrix(df[, cols, drop = FALSE])
    storage.mode(M) <- "numeric"
    n_complete <- rowSums(!is.na(M))

    row_mean   <- rowMeans(M, na.rm = TRUE)
    # Row sd via formula (avoids the apply overhead)
    row_var    <- rowSums((M - row_mean)^2 * !is.na(M), na.rm = TRUE) /
                  pmax(n_complete - 1L, 1L)
    row_sd     <- sqrt(pmax(row_var, 0))
    row_min    <- suppressWarnings(apply(M, 1L, min, na.rm = TRUE))
    row_max    <- suppressWarnings(apply(M, 1L, max, na.rm = TRUE))
    row_range  <- ifelse(is.finite(row_min) & is.finite(row_max),
                          row_max - row_min, NA_real_)
    row_peak   <- suppressWarnings(apply(M, 1L,
                          function(r) if (all(is.na(r))) NA_integer_
                                       else which.max(r)))
    row_cv     <- ifelse(abs(row_mean) > .Machine$double.eps,
                          row_sd / abs(row_mean), NA_real_)

    # Replace -Inf / Inf from all-NA rows with NA
    row_min[!is.finite(row_min)] <- NA
    row_max[!is.finite(row_max)] <- NA

    pfx <- sprintf("gee_%s_%s", family, scale)
    df[[paste0(pfx, "_annual_mean")]]    <- row_mean
    df[[paste0(pfx, "_seasonal_sd")]]    <- row_sd
    df[[paste0(pfx, "_seasonal_cv")]]    <- row_cv
    df[[paste0(pfx, "_seasonal_range")]] <- row_range
    df[[paste0(pfx, "_peak_month")]]     <- row_peak
  }
  if (verbose) {
    new_cols <- grep("_(annual_mean|seasonal_sd|seasonal_cv|seasonal_range|peak_month)$",
                      colnames(df), value = TRUE)
    message(sprintf("[seasonality] added %d derived features across %d (family, scale) pairs",
                     length(new_cols), length(unique(fam_scale))))
  }
  df
}
