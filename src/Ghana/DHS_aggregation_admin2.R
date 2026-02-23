#-------------------------------------------------------------------------------
# Robust indicator building + Admin-2 estimation for Ghana 2016 MIS
#-------------------------------------------------------------------------------

rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(labelled)
  library(surveyPrev)  # surveyPrevIndicators, getDHSindicator(), directEST()
  # library(expss)     # optional; we stub if not installed
})

# --- Load your MIS objects ----------------------------------------------------
d <- readRDS(here("data/DHS/dhs_Ghana_2016.RDS"))

PRdata <- d$Ghana_2016$PRdata
HRdata <- d$Ghana_2016$HRdata
KRdata <- d$Ghana_2016$KRdata
IRdata <- d$Ghana_2016$IRdata
MRdata <- NULL   # set if available
ARdata <- NULL   # set if available

DHS_REPO_ZIP <- if (file.exists(here("src/DHS/DHS-Indicators-R-main.zip"))) {
  here("src/DHS/DHS-Indicators-R-main.zip")
} else {
  "/mnt/data/DHS-Indicators-R-main.zip"
}

#-------------------------------------------------------------------------------
# Helpers
#-------------------------------------------------------------------------------
.make_dhs_list <- function(PRdata=NULL, KRdata=NULL, HRdata=NULL, IRdata=NULL, MRdata=NULL, ARdata=NULL, BRdata=NULL, CRdata=NULL) {
  lst <- list()
  if (!is.null(PRdata)) lst$PRdata <- PRdata
  if (!is.null(HRdata)) lst$HRdata <- HRdata
  if (!is.null(KRdata)) lst$KRdata <- KRdata
  if (!is.null(IRdata)) lst$IRdata <- IRdata
  if (!is.null(MRdata)) lst$MRdata <- MRdata
  if (!is.null(ARdata)) lst$ARdata <- ARdata
  if (!is.null(BRdata)) lst$BRdata <- BRdata
  if (!is.null(CRdata)) lst$CRdata <- CRdata
  lst
}

# Standardize one binary variable for surveyPrev/directEST()
.standardize_one <- function(df, var, which) {
  stopifnot(var %in% names(df))
  if (which %in% c("PRdata","HRdata")) {
    cluster <- df[["hv001"]]; householdID <- df[["hv002"]]
    region  <- df[["hv024"]]; weight <- df[["hv005"]] / 1e6
    strata  <- if ("hv023" %in% names(df)) df[["hv023"]] else df[["hv025"]]
  } else {
    cluster <- df[["v001"]];  householdID <- df[["v002"]]
    region  <- df[["v024"]];  weight <- df[["v005"]] / 1e6
    strata  <- if ("v023" %in% names(df)) df[["v023"]] else df[["v025"]]
  }
  out <- tibble(
    cluster     = cluster,
    householdID = householdID,
    v024        = region,
    weight      = weight,
    strata      = strata,
    value       = df[[var]]
  ) %>% filter(!is.na(value))
  attr(out, "varname") <- var
  attr(out, "recode")  <- which
  out
}

#-------------------------------------------------------------------------------
# Lenient DHS repo sourcing with no-op expss stubs
#-------------------------------------------------------------------------------
.attach_expss_stubs <- function(E, verbose = TRUE) {
  if (!requireNamespace("expss", quietly = TRUE)) {
    if (verbose) message("[note] 'expss' not installed; creating no-op stubs for common functions.")
    no_op <- function(x, ...) x
    for (fn in c("apply_labels","add_val_labels","drop_val_labels","drop_var_labs",
                 "add_labels","remove_labels","val_lab","var_lab","to_factor",
                 "set_value_labels","set_variable_labels")) {
      if (!exists(fn, envir = E, inherits = FALSE)) assign(fn, no_op, envir = E)
    }
  }
}

.source_lenient <- function(file, envir, verbose = TRUE) {
  exprs <- tryCatch(parse(file = file, keep.source = TRUE),
                    error = function(e) { if (verbose) message("[parse fail] ", basename(file), " :: ", e$message)
                      return(expression()) })
  logs <- tibble(file = basename(file), idx = integer(), ok = logical(), message = character(), head = character())

  for (i in seq_along(exprs)) {
    expr <- exprs[[i]]
    head_txt <- tryCatch(paste(deparse(expr, width.cutoff = 80)[1]), error = function(e) "")
    warn_msgs <- character(0)

    res <- tryCatch({
      withCallingHandlers(
        {
          eval(expr, envir = envir)
          NULL
        },
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          if (!is.null(base::findRestart("muffleWarning"))) invokeRestart("muffleWarning")
        }
      )
      list(ok = TRUE, msg = if (length(warn_msgs)) paste("warn:", paste(unique(warn_msgs), collapse=" | ")) else "")
    }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))

    if (verbose && !res$ok) message(sprintf("[skip expr] %s [%d]  %s", basename(file), i, res$msg))
    logs <- bind_rows(logs, tibble(file = basename(file), idx = i, ok = res$ok, message = res$msg, head = head_txt))
  }
  logs
}

# Binary detector: numeric 0/1, logical, or factor with specific 2-level encodings
.is_binary_column <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.logical(x)) return(any(!is.na(x)))
  if (is.numeric(x) || is.integer(x)) return(all(x %in% c(0,1,NA), na.rm = TRUE) && any(!is.na(x)))
  if (is.factor(x)) {
    lev <- levels(x)
    # Safe cases only: levels exactly "0"/"1" or "no"/"yes" (case-insensitive)
    if (length(lev) == 2) {
      if (all(tolower(lev) %in% c("0","1"))) return(TRUE)
      if (all(tolower(lev) %in% c("no","yes"))) return(TRUE)
    }
  }
  FALSE
}

# Convert factor/logical safely to 0/1 (only for the safe cases above)
.to_binary01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))       # FALSE=0, TRUE=1
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  if (is.factor(x)) {
    lev <- levels(x)
    if (length(lev) == 2) {
      if (all(tolower(lev) %in% c("0","1"))) {
        # map "0"->0, "1"->1
        return(as.integer(as.character(x)))
      }
      if (all(tolower(lev) %in% c("no","yes"))) {
        # map "no"->0, "yes"->1
        return(as.integer(tolower(as.character(x)) == "yes"))
      }
    }
  }
  x
}

.discover_and_standardize <- function(E, dhs_list_in) {
  out <- list()
  disc_log <- tibble(source="repo", indicator = character(), ok = logical(), message = character())

  for (nm in names(dhs_list_in)) {
    if (!exists(nm, envir = E)) next
    df <- get(nm, envir = E)
    old <- names(dhs_list_in[[nm]])
    cand <- setdiff(names(df), old)

    bin_vars <- cand[vapply(cand, function(v) .is_binary_column(df[[v]]), logical(1))]

    # Coerce safe non-numeric binaries to 0/1
    if (length(bin_vars)) {
      for (v in bin_vars) {
        if (!is.numeric(df[[v]])) {
          df[[v]] <- .to_binary01(df[[v]])
        }
      }
      assign(nm, df, envir = E)
    }

    for (v in bin_vars) {
      res <- tryCatch({ .standardize_one(df, v, nm) },
                      error = function(e) structure(NULL, error = conditionMessage(e)))
      ok  <- !is.null(res)
      msg <- if (ok) "" else attr(res, "error", exact = TRUE)
      key <- paste0("CUSTOM::", nm, "::", v)
      if (ok) out[[key]] <- res
      disc_log <- bind_rows(disc_log, tibble(source="repo", indicator = key, ok = ok, message = msg))
    }
  }

  list(data = out, log = disc_log)
}

.compute_from_repo_lenient <- function(dhs_list, zip_path, verbose = TRUE, file_patterns = NULL) {
  if (is.null(zip_path) || !file.exists(zip_path)) {
    return(list(data = list(), log = tibble(source="repo", indicator=character(), ok=logical(), message=character())))
  }
  E <- rlang::env()
  for (nm in names(dhs_list)) assign(nm, dhs_list[[nm]], envir = E)
  for (nm in c("PRdata","HRdata","IRdata","KRdata","MRdata","ARdata","BRdata","CRdata")) {
    if (!exists(nm, envir = E, inherits = FALSE)) assign(nm, tibble(), envir = E)
  }
  .attach_expss_stubs(E, verbose = verbose)

  tmpd <- file.path(tempdir(), paste0("DHS-Indicators-R_", as.integer(runif(1,1e6,1e7))))
  dir.create(tmpd, showWarnings = FALSE, recursive = TRUE)
  utils::unzip(zipfile = zip_path, exdir = tmpd)

  r_files <- list.files(tmpd, pattern = "\\.[rR]$", recursive = TRUE, full.names = TRUE)
  if (!is.null(file_patterns)) {
    keep <- vapply(r_files, function(p) any(str_detect(basename(p), file_patterns)), logical(1))
    r_files <- r_files[keep]
  }

  src_logs <- map_dfr(r_files, ~.source_lenient(.x, envir = E, verbose = verbose))
  harvested <- .discover_and_standardize(E, dhs_list)

  list(
    data = harvested$data,
    log  = bind_rows(
      tibble(source="repo-source",
             indicator = paste0(basename(src_logs$file), ":", src_logs$idx),
             ok = src_logs$ok, message = src_logs$message),
      harvested$log
    )
  )
}

#-------------------------------------------------------------------------------
# Built-ins (22 indicators): try sensible recode sets; keep first success
#-------------------------------------------------------------------------------
.safe_get_builtin <- function(id, dhs_list, verbose = TRUE) {
  combos <- list(
    c("HRdata"), c("PRdata"),              # MIS-friendly first
    c("IRdata"), c("KRdata"), c("BRdata"),
    c("MRdata"), c("ARdata"),
    c("IRdata","MRdata","ARdata")         # HIV
  )
  msgs <- character(0)
  for (need in combos) {
    if (!all(need %in% names(dhs_list))) next
    dd <- dhs_list[need]
    res <- tryCatch({
      surveyPrev::getDHSindicator(dd, indicator = id)
    }, error = function(e) { msgs <<- c(msgs, sprintf("[%s] %s", paste(need, collapse="+"), e$message)); NULL })
    if (is.data.frame(res) && "value" %in% names(res)) {
      if (verbose) message(sprintf("[surveyPrev] %s -> OK using %s", id, paste(need, collapse="+")))
      return(list(ok=TRUE, data=res, msg=""))
    }
  }
  if (verbose) message(sprintf("[surveyPrev] %s -> SKIP", id))
  list(ok=FALSE, data=NULL, msg=paste(unique(msgs), collapse=" | "))
}

.compute_builtin <- function(dhs_list, ids = surveyPrevIndicators$ID, verbose = TRUE) {
  logs <- tibble(source="built-in", indicator = character(), ok = logical(), message = character())
  out  <- vector("list", length(ids)); names(out) <- ids
  for (id in ids) {
    got <- .safe_get_builtin(id, dhs_list, verbose)
    logs <- bind_rows(logs, tibble(source="built-in", indicator=id, ok=got$ok, message=got$msg))
    out[[id]] <- got$data
  }
  list(data = out, log = logs)
}

#-------------------------------------------------------------------------------
# High-level wrapper
#-------------------------------------------------------------------------------
compute_all_indicators <- function(
    PRdata=NULL, KRdata=NULL, HRdata=NULL, IRdata=NULL, MRdata=NULL, ARdata=NULL, BRdata=NULL, CRdata=NULL,
    dhs_repo_zip = NULL,
    include_builtin = TRUE,
    include_repo   = TRUE,
    verbose = TRUE,
    repo_file_patterns = c("Ch10","Ch11","Ch15","ML","WE","WS")
) {
  dhs_list <- .make_dhs_list(PRdata, KRdata, HRdata, IRdata, MRdata, ARdata, BRdata, CRdata)

  built <- if (include_builtin) .compute_builtin(dhs_list, ids = surveyPrevIndicators$ID, verbose = verbose)
  else list(data=list(), log=tibble())
  repo  <- if (include_repo && !is.null(dhs_repo_zip))
    .compute_from_repo_lenient(dhs_list, zip_path = dhs_repo_zip, verbose = verbose,
                               file_patterns = repo_file_patterns)
  else list(data=list(), log=tibble())

  data_list <- c(built$data, repo$data)
  # keep only those that actually exist
  data_list <- data_list[!vapply(data_list, is.null, logical(1))]

  list(
    data = data_list,
    log  = bind_rows(built$log, repo$log)
  )
}

#-------------------------------------------------------------------------------
# Admin-2 direct estimates for all successful indicators
#-------------------------------------------------------------------------------
.safe_direct_admin2 <- function(dat, cluster.info) {
  if (is.null(dat) || !is.data.frame(dat)) return(list(error="NULL input"))
  need <- c("cluster","householdID","v024","weight","strata","value")
  if (!all(need %in% names(dat))) return(list(error="Missing required columns"))
  if (all(is.na(dat$value))) return(list(error="All NA 'value'"))
  tryCatch({
    res <- surveyPrev::directEST(data = dat, cluster.info = cluster.info, admin = 2)
    list(result = res$res.admin2)
  }, error=function(e) list(error = conditionMessage(e)))
}

estimate_all_admin2 <- function(indicators, cluster.info, prefer_MIS = TRUE, only = NULL, exclude = NULL, drop_empty = TRUE) {
  processed <- if (is.list(indicators) && !is.null(indicators$data)) indicators$data else indicators
  if (!is.list(processed) || length(processed) == 0) stop("No indicators supplied.")

  if (prefer_MIS) {
    processed <- processed[vapply(processed, function(dat) {
      rec <- attr(dat, "recode", exact=TRUE)
      isTRUE(rec %in% c("PRdata","HRdata"))
    }, logical(1))]
  }
  if (!is.null(only))    processed <- processed[names(processed) %in% only]
  if (!is.null(exclude)) processed <- processed[!names(processed) %in% exclude]

  out_rows <- list(); fails <- tibble(indicator=character(), reason=character())
  for (nm in names(processed)) {
    dat   <- processed[[nm]]
    n_raw <- sum(!is.na(dat$value))
    rec   <- attr(dat,"recode",  exact=TRUE)
    var   <- attr(dat,"varname", exact=TRUE)
    ans   <- .safe_direct_admin2(dat, cluster.info)
    if (!is.null(ans$error)) { fails <- bind_rows(fails, tibble(indicator=nm, reason=ans$error)); next }
    df <- ans$result
    df$indicator <- nm; df$recode <- rec; df$varname <- var; df$n_raw <- n_raw
    if (drop_empty && "direct.est" %in% names(df)) df <- df[!is.na(df$direct.est), , drop=FALSE]
    out_rows[[nm]] <- df
  }
  list(
    admin2 = if (length(out_rows)) bind_rows(out_rows) else tibble(),
    failed = fails,
    kept   = names(processed)
  )
}

#-------------------------------------------------------------------------------
# RUN
#-------------------------------------------------------------------------------

all_ind <- compute_all_indicators(
  PRdata = PRdata, HRdata = HRdata,
  # Include the others only if you truly have them:
  # IRdata = IRdata, KRdata = KRdata, MRdata = MRdata, ARdata = ARdata,
  dhs_repo_zip       = DHS_REPO_ZIP,
  include_builtin    = TRUE,
  include_repo       = TRUE,
  verbose            = TRUE,
  repo_file_patterns = c("Ch10","Ch11","Ch15","ML","WE","WS")
)

# Inspect success/failure and what we actually harvested
print(dplyr::count(all_ind$log, source, ok))
cat("Number of indicators kept:", length(all_ind$data), "\n")
cat("First few indicator names:\n"); print(utils::head(names(all_ind$data), 30))

# Example Admin-2 direct estimates (once you have 'cluster.info')
# adm2 <- estimate_all_admin2(all_ind, cluster.info, prefer_MIS = TRUE)
# head(adm2$admin2); adm2$failed

all_ind$data
