# =============================================================================
# fe_features.R — engineered-feature constructors operating on a numeric matrix
# (column names from the cache). Each returns a NEW numeric matrix of features.
# These are deterministic per-row transforms (no train stats) => leak-free to
# compute up-front and then feed through the normal prep pipeline.
# =============================================================================
suppressMessages(library(matrixStats))

MONTHS <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

# Identify monthly-climatology groups: gee_<base>_<Month>_<suffix>
.month_groups <- function(cols) {
  pat <- paste0("_(", paste(MONTHS, collapse="|"), ")_")
  has <- grepl(pat, cols)
  cm <- cols[has]; if (!length(cm)) return(list())
  mo <- regmatches(cm, regexpr(paste0("(", paste(MONTHS,collapse="|"), ")"), cm))
  stem <- sub(pat, "_MO_", cm)               # collapse month token
  split(data.frame(col=cm, month=match(mo,MONTHS), stringsAsFactors=FALSE), stem)
}

#' Build seasonality summary features from monthly climatology columns.
#' For each base variable with >=8 months present: mean, amplitude, CV,
#' peak month (cyclic via sin/cos), 1st Fourier harmonic amplitude.
build_seasonality <- function(M) {
  cols <- colnames(M); grps <- .month_groups(cols)
  grps <- Filter(function(g) nrow(g) >= 8, grps)
  if (!length(grps)) return(matrix(numeric(0), nrow=nrow(M), ncol=0))
  ang <- function(mo) 2*pi*mo/12
  out <- list()
  for (stem in names(grps)) {
    g <- grps[[stem]]; g <- g[order(g$month),]
    X <- M[, g$col, drop=FALSE]
    mu  <- rowMeans(X, na.rm=TRUE)
    amp <- matrixStats::rowMaxs(X, na.rm=TRUE) - matrixStats::rowMins(X, na.rm=TRUE)
    sdv <- matrixStats::rowSds(X, na.rm=TRUE)
    cv  <- sdv / ifelse(abs(mu) < 1e-8, NA, abs(mu))
    cosw <- cos(ang(g$month)); sinw <- sin(ang(g$month))
    a <- as.numeric(X %*% cosw) * 2/length(g$month)
    b <- as.numeric(X %*% sinw) * 2/length(g$month)
    h1 <- sqrt(a^2 + b^2)
    phase <- atan2(b, a)
    base <- gsub("[^A-Za-z0-9]+","_", sub("^gee_","", stem))
    base <- sub("_MO_.*$","", base)
    out[[paste0("fe_",base,"_mean")]]  <- mu
    out[[paste0("fe_",base,"_amp")]]   <- amp
    out[[paste0("fe_",base,"_cv")]]    <- cv
    out[[paste0("fe_",base,"_h1")]]    <- h1
    out[[paste0("fe_",base,"_phsin")]] <- sin(phase)
    out[[paste0("fe_",base,"_phcos")]] <- cos(phase)
  }
  M2 <- do.call(cbind, out); colnames(M2) <- names(out); M2
}

# Columns that are population-scaled COUNTS (vs rates/percent/prevalence).
is_count_col <- function(cols) {
  grepl("count|_counts$|_Count$|_number_", cols, ignore.case=TRUE) &
  !grepl("percent|rate|prevalence|_pr_|proportion|prop_", cols, ignore.case=TRUE)
}

# columns that are monthly climatology (the raw 12-per-var)
is_monthly_col <- function(cols) grepl(paste0("_(", paste(MONTHS,collapse="|"), ")_"), cols)
