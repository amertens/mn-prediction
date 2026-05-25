# =============================================================================
# R/feature_engineering.R  (PROTOTYPE)
#
# Derive transportability-friendly engineered features from the EXISTING raw
# Admin-2 GEE extraction (no new raster extraction required):
#
#   1. Anomalies / trends   — from year-suffixed columns (gee_<var>_YYYY):
#        z-score of the survey-year value vs the area's own multi-year
#        climatology, plus the multi-year mean and a simple trend. Anomalies
#        are unit-free and comparable ACROSS countries, directly targeting the
#        "absolute level differs by country" failure mode in LOCO transfer.
#
#   2. Seasonality          — from month-suffixed columns (..._monthlyMM_...):
#        intra-annual amplitude, coefficient of variation, peak month, and the
#        first Fourier harmonic (amplitude + phase) of the annual cycle.
#
#   3. Soil x crop          — products of soil-mineral content and cropland
#        intensity, a mechanistic proxy for dietary mineral availability.
#
# All engineered columns are prefixed "fe_" and keyed on Admin2.
# =============================================================================

# --- helpers ----------------------------------------------------------------
.fe_num <- function(v) {
  if (inherits(v, "haven_labelled")) return(as.double(unclass(v)))
  suppressWarnings(as.numeric(v))
}

# Group year-suffixed columns: "gee_ndvi_2017" / "gee_fldas_2016_annual_mean_x"
# -> stem (year token removed) + the year. Returns list(stem = data.frame(col, year)).
.fe_year_groups <- function(cols) {
  m <- regmatches(cols, regexpr("(?<![0-9])(19|20)[0-9]{2}(?![0-9])", cols, perl = TRUE))
  has_year <- vapply(cols, function(c) grepl("(?<![0-9])(19|20)[0-9]{2}(?![0-9])", c, perl = TRUE),
                     logical(1))
  cols2 <- cols[has_year]; if (!length(cols2)) return(list())
  years <- as.integer(regmatches(cols2, regexpr("(?<![0-9])(19|20)[0-9]{2}(?![0-9])", cols2, perl = TRUE)))
  stems <- gsub("_(19|20)[0-9]{2}", "", cols2)
  stems <- gsub("__+", "_", sub("_$", "", stems))
  split(data.frame(col = cols2, year = years, stringsAsFactors = FALSE), stems)
}

# Group month-suffixed columns: "gee_fldas_monthly01_lwnet_tavg" -> stem + month.
.fe_month_groups <- function(cols) {
  has_m <- grepl("monthly[0-9]{2}", cols)
  cols2 <- cols[has_m]; if (!length(cols2)) return(list())
  months <- as.integer(regmatches(cols2, regexpr("(?<=monthly)[0-9]{2}", cols2, perl = TRUE)))
  stems <- gsub("monthly[0-9]{2}_?", "", cols2)
  stems <- gsub("__+", "_", sub("_$", "", stems))
  split(data.frame(col = cols2, month = months, stringsAsFactors = FALSE), stems)
}

#' Engineer Admin-2 features from a raw gee_admin2 data.frame
#'
#' @param gee_admin2 raw Admin-2 GEE table (Admin2 + per-year / per-month cols)
#' @param survey_year integer survey year (for anomaly alignment)
#' @param max_groups cap on number of year/month groups processed (speed)
#' @return data.frame keyed on Admin2 with fe_* columns
engineer_admin2_features <- function(gee_admin2, survey_year = NA_integer_,
                                     max_groups = 400L) {
  gcols <- grep("^gee_", names(gee_admin2), value = TRUE)
  out <- data.frame(Admin2 = as.character(gee_admin2$Admin2), stringsAsFactors = FALSE)

  # ---- 1. anomalies / trend from year groups ----
  yg <- .fe_year_groups(gcols)
  yg <- Filter(function(g) nrow(g) >= 3, yg)               # need >=3 years
  if (length(yg) > max_groups) yg <- yg[seq_len(max_groups)]
  for (stem in names(yg)) {
    g <- yg[[stem]]; g <- g[order(g$year), ]
    M <- sapply(g$col, function(c) .fe_num(gee_admin2[[c]]))
    if (is.null(dim(M))) next
    mu <- rowMeans(M, na.rm = TRUE)
    sdv <- apply(M, 1, sd, na.rm = TRUE)
    # survey-year aligned value (nearest available year)
    yj <- if (!is.na(survey_year)) which.min(abs(g$year - survey_year)) else ncol(M)
    cur <- M[, yj]
    anom <- (cur - mu) / ifelse(is.finite(sdv) & sdv > 1e-8, sdv, NA)
    trend <- (M[, ncol(M)] - M[, 1]) / (g$year[nrow(g)] - g$year[1])
    base <- gsub("[^A-Za-z0-9]+", "_", stem)
    out[[paste0("fe_", base, "_anom")]]  <- anom
    out[[paste0("fe_", base, "_mean")]]  <- mu
    out[[paste0("fe_", base, "_trend")]] <- trend
  }

  # ---- 2. seasonality from month groups ----
  mg <- .fe_month_groups(gcols)
  mg <- Filter(function(g) nrow(g) >= 6, mg)               # need a decent cycle
  if (length(mg) > max_groups) mg <- mg[seq_len(max_groups)]
  ang <- function(mo) 2 * pi * mo / 12
  for (stem in names(mg)) {
    g <- mg[[stem]]; g <- g[order(g$month), ]
    M <- sapply(g$col, function(c) .fe_num(gee_admin2[[c]]))
    if (is.null(dim(M))) next
    mu <- rowMeans(M, na.rm = TRUE)
    amp <- apply(M, 1, function(r) { r <- r[is.finite(r)]; if (length(r) < 2) NA else max(r) - min(r) })
    cv  <- apply(M, 1, function(r) { r <- r[is.finite(r)]; if (length(r) < 2 || mean(r) == 0) NA else sd(r) / abs(mean(r)) })
    peak <- apply(M, 1, function(r) if (all(is.na(r))) NA else g$month[which.max(r)])
    # first Fourier harmonic
    cosw <- cos(ang(g$month)); sinw <- sin(ang(g$month))
    a <- as.numeric(M %*% cosw) / length(g$month) * 2
    b <- as.numeric(M %*% sinw) / length(g$month) * 2
    base <- gsub("[^A-Za-z0-9]+", "_", stem)
    out[[paste0("fe_", base, "_amp")]]   <- amp
    out[[paste0("fe_", base, "_cv")]]    <- cv
    out[[paste0("fe_", base, "_h1amp")]] <- sqrt(a^2 + b^2)
    out[[paste0("fe_", base, "_phase")]] <- atan2(b, a)
    out[[paste0("fe_", base, "_peakmo")]] <- peak
  }

  # ---- 3. soil x crop interactions ----
  soil_cols <- grep("^gee_soil.*mean", gcols, value = TRUE)
  crop_cols <- grep("worldcereal.*(annual_mean|annual_max)", gcols, value = TRUE)
  if (length(soil_cols) && length(crop_cols)) {
    crop <- .fe_num(gee_admin2[[crop_cols[1]]])
    cropz <- (crop - mean(crop, na.rm = TRUE)) / sd(crop, na.rm = TRUE)
    for (sc in head(soil_cols, 25L)) {
      s <- .fe_num(gee_admin2[[sc]])
      sz <- (s - mean(s, na.rm = TRUE)) / sd(s, na.rm = TRUE)
      base <- gsub("[^A-Za-z0-9]+", "_", sub("^gee_", "", sc))
      out[[paste0("fe_", base, "_x_crop")]] <- sz * cropz
    }
  }

  out
}
