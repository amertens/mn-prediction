# =============================================================================
# scripts/protocol_v2/71_person_level_honest_lib.R   [IL-02e]
# Sourced by 71_person_level_honest.R and each of its workers.
#
# THE JANUARY PERSON-LEVEL FIGURE, UPDATED TO THE CURRENT DATA AND MODELS
#
# Current data: the targets store's outcome_data_* for the country, with the
# pipeline's own outcome definitions (resolve_uniform_outcome for the flags,
# VITA_RULE = rbp070 by default; the 01_build_targets_v2.R level definitions for
# concentrations, on the log scale). Current models: IL-01's person-level
# SuperLearner (fit_area_superlearner: mean / elastic net / ranger, survey
# weights, district-blocked inner folds, discrete pick, family binomial for
# flags and gaussian for concentrations), IL-01's district domain PCs as the
# proxies, and the protocol-v2 PCA domain index. Scored out of fold on IL-01's
# district-blocked 5-fold assignments (group_folds seeds), REPS draws.
#   Survey only        the respondent's own gw_ columns, after the project guard
#                      allowed_under_arm("questionnaire") (blood draw and Hb),
#                      IL-01's regex with ages exempted from its date patterns,
#                      and the identifier / sampling-design / region / team /
#                      date / RDT / genotype columns IL-02 found in IL-01's sets.
#                      No 80-column cap. 70 percent coverage, median-imputed.
#   Proxies            IL-01's px_ domain PCs of the district
#   Survey + proxies   both
#   Index              domain PCs of the surveyed districts weighted by their
#                      training-district Spearman with the district outcome
#                      (survey-weighted prevalence on the logit scale, or mean log
#                      concentration), then a person-level logistic / linear
#                      calibration fitted on the training respondents
#   Ceiling            not a model: each respondent's district value among the
#                      OTHER respondents, shrunk to the national value
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(targets)})
invisible(capture.output(suppressMessages(targets::tar_source("R"))))
STORE71 <- "_targets_full"
COUNTRY71 <- Sys.getenv("IL_HONEST_COUNTRY", "Ghana")
cfg71 <- get_country_configs()[[COUNTRY71]]
OUTS71 <- intersect(c("child_vitA", "women_vitA", "child_iron", "women_iron", "women_folate", "women_b12"), names(cfg71$outcomes))

S71  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD71 <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS71 <- drop_near_outcome_v2(intersect(MD71$column, names(S71)), MD71)
dom71 <- stats::setNames(MD71$domain, MD71$column)

# IL-01's regex (46_individual_level_models.R), split so that ages survive its date patterns
IL01_LEAK <- paste0("RBP|rbp|VAD|vad|LogFer|logfer|Ferr|ferr|IDA|ida$|Brinda|BRINDA|Thurn|Folate|folate|Fol|B12|b12|Zinc|zinc|zn_|Hb|hb$|Hgb|hgb|Haem|Hem|anaem|anem|Anem|CRP|crp|AGP|agp|infl|Infl|UIC|uic|Iod|iod|Salt|salt|MUAC|muac|weight|Weight|Wt$|wt$|_id$|ID$|Id$|cnum|clust|Clust|psu|PSU|strat|Strat|Admin|admin|lat|Lat|lon|Lon|GPS|gps|Team|team|line|Line|hhid|caseid|Retinol|retinol|Vit|vit|Anemia|Malaria|malaria|RDT|rdt|Plasmod|Sickle|G6PD|Transferrin|sTfR|stfr|ZPP|zpp|Ret$|Def$|def$|Adj|adj",
                    "|^sf_|ferritin|_nmol|^fol|^zn|^map2_|^rbp|^vit|_def$")
IL01_DATE <- "Date|date|month|Month|year|Year"
# found by IL-02 in IL-01's published sets (il02_il01_survey_columns.csv), plus the regex's misses
IL02_DROP <- paste0("cnum|^gw_cn$|bccn|b_cn$|hhln|bchnc|bclnr|bcgln|b_hhn|wlnr|^gw_in1$|^gw_wp1$|b_wp1|indivID|momID|childid|",
                    "pcn$|^gw_mcn$|Region|region|Number_of|probab|Check_|sWeight|Total.population|Households|cttp|datasource|",
                    "INT_m|a_mon|gcmst|MalariaYN|wrmal|malref|sickle|thal|Thal|mrdr|pctopendef|Phleb|phleb|Urine|urine|Specs|",
                    "IntNumb|InterNumb|HHNumb|WomanNumb|RandomSelect|^gw_fer$|^gw_bis$|TFR|^gw_HH1$|^gw_HH7$|^gw_lga$|",
                    # fieldwork timing (encodes the team's route), a line number, and the whole biomarker-collection
                    # modules: child gc (RDT, Hb), hc (malaria referral), cp (phlebotomy); women wm (RDT, blood volume), wo
                    "week$|intmon|fieldwork|cbamln|^gw_(gc|hc|cp)|^gw_w[mo]_")

lg71 <- function(v) { v <- .v2_num(v); v[!is.finite(v) | v <= 0] <- NA; log(v) }
.cache71 <- new.env()
load_cell71 <- function(outcome) {
  if (!is.null(.cache71[[outcome]])) return(.cache71[[outcome]])
  oc <- cfg71$outcomes[[outcome]]
  od <- tar_read_raw(paste0("outcome_data_", tolower(COUNTRY71), "_", outcome), store = STORE71); d <- od$data
  ybin <- tryCatch(suppressMessages(resolve_uniform_outcome(d, cfg71, oc, label = "[il02e]")), error = function(e) NULL)
  ybin <- if (!is.null(ybin)) .v2_num(ybin) else .v2_num(d[[oc$binary]])
  # concentration, as 01_build_targets_v2.R builds the level target (not negated here)
  ycont <- rep(NA_real_, nrow(d)); cont_src <- "none"
  adj <- if (grepl("vitA", outcome)) tryCatch(suppressMessages(brinda_vad_adjusted(d, cfg71, oc, label = "[il02e level]")), error = function(e) NULL) else NULL
  if (!is.null(adj)) { ycont <- lg71(adj); cont_src <- "BRINDA-adjusted RBP"
  } else if (!is.null(oc$continuous) && oc$continuous %in% names(d)) {
    v <- .v2_num(d[[oc$continuous]]); ycont <- if (identical(oc$cutoff_scale, "log")) v else lg71(v); cont_src <- oc$continuous }
  w <- .v2_num(d[[cfg71$weight_col]]); w[!is.finite(w) | w <= 0] <- NA
  dist <- paste(as.character(d$Admin1), as.character(d$Admin2), sep = "||")
  psu <- as.character(d[[cfg71$psu_col]])
  # survey columns
  outcome_vars <- unique(unlist(lapply(cfg71$outcomes, function(o) c(o$binary, o$continuous))))
  gw <- names(d)[grepl("^gw_", names(d))]
  gw <- gw[allowed_under_arm(gw, "questionnaire") & !gw %in% outcome_vars & !grepl(IL01_LEAK, gw) &
           !(grepl(IL01_DATE, gw) & !grepl("age|Age", gw)) & !grepl(IL02_DROP, gw)]
  gw <- gw[vapply(gw, function(k) { v <- .v2_num(d[[k]]); mean(is.finite(v)) >= 0.7 && isTRUE(stats::sd(v, na.rm = TRUE) > 0) && length(unique(v[is.finite(v)])) >= 2 }, TRUE)]
  # Malawi's store data carries no gw_ questionnaire columns: an empty survey set, and the survey arms are skipped
  Xs <- if (length(gw)) as.matrix(as.data.frame(lapply(d[gw], .v2_num))) else matrix(numeric(0), nrow = nrow(d), ncol = 0)
  if (ncol(Xs)) { Xs <- apply(Xs, 2, function(v) { v[!is.finite(v)] <- stats::median(v, na.rm = TRUE); v })
    if (is.null(dim(Xs))) Xs <- matrix(Xs, ncol = 1, dimnames = list(NULL, gw)) }
  # IL-01's proxies: domain PCs over the country's rows of the shared table
  Sx <- S71[S71$country == COUNTRY71, ]; Dm <- domain_representation_v2(prep_predictors_v2(as.matrix(Sx[, PREDS71, drop = FALSE])), dom71)
  colnames(Dm) <- paste0("px_", colnames(Dm)); Xp <- Dm[match(dist, paste(Sx$Admin1, Sx$Admin2, sep = "||")), , drop = FALSE]
  # the index basis: domain PCs of the surveyed districts only (build_cell's rule)
  keys <- sort(unique(dist[is.finite(Xp[, 1])])); sc <- Sx[match(keys, paste(Sx$Admin1, Sx$Admin2, sep = "||")), PREDS71, drop = FALSE]
  Dix <- domain_representation_v2(prep_predictors_v2(as.matrix(sc)), dom71)
  cell <- list(d_n = nrow(d), ybin = ybin, ycont = ycont, cont_src = cont_src, w = w, dist = dist, psu = psu, Xs = Xs, Xp = Xp,
               survey_cols = gw, keys = keys, Dix = Dix)
  assign(outcome, cell, envir = .cache71); cell
}
rows71 <- function(cell, type) { y <- if (type == "bin") cell$ybin else cell$ycont
  which(is.finite(y) & is.finite(cell$w) & !is.na(cell$dist) & is.finite(cell$Xp[, 1]) & cell$dist %in% cell$keys) }
group_folds71 <- function(groups, k, rep_id) { set.seed(20260951L + rep_id); g <- unique(groups); f <- sample(rep(seq_len(min(k, length(g))), length.out = length(g))); f[match(groups, g)] }

# one SL arm, one draw (IL-01's learner call, unchanged)
run_sl71 <- function(type, outcome, set, rep) {
  cell <- load_cell71(outcome); r <- rows71(cell, type); y <- (if (type == "bin") cell$ybin else cell$ycont)[r]
  w <- cell$w[r]; dist <- cell$dist[r]
  X <- switch(set, survey = cell$Xs[r, , drop = FALSE], proxies = cell$Xp[r, , drop = FALSE],
              both = cbind(cell$Xs[r, , drop = FALSE], cell$Xp[r, , drop = FALSE]))
  folds <- group_folds71(dist, 5, rep); p <- p_null <- rep(NA_real_, length(y)); picks <- character(0)
  fam <- if (type == "bin") stats::binomial() else stats::gaussian()
  for (f in sort(unique(folds))) { te <- which(folds == f); tr <- which(folds != f)
    fit <- tryCatch(fit_area_superlearner(y[tr], X[tr, , drop = FALSE], newX = X[te, , drop = FALSE], weights = w[tr], block = dist[tr],
                                          library = c("mean", "enet", "ranger"), V = 3L, discrete = TRUE, family = fam, meta = "mse"),
                    error = function(e) NULL)
    p[te] <- if (!is.null(fit) && length(fit$pred_new) == length(te)) fit$pred_new else mean(y[tr])
    if (!is.null(fit)) picks <- c(picks, fit$pick)
    p_null[te] <- mean(y[tr]) }
  if (type == "bin") p <- pmin(pmax(p, 1e-4), 1 - 1e-4)
  data.frame(type = type, outcome = outcome, set = set, rep = rep, row = r, y = y, pred = p, null = p_null, fold = folds,
             picks = paste(names(table(picks)), table(picks), collapse = ";"))
}

# the index, one draw, on the same folds
run_index71 <- function(type, outcome, rep) {
  cell <- load_cell71(outcome); r <- rows71(cell, type); y <- (if (type == "bin") cell$ybin else cell$ycont)[r]
  w <- cell$w[r]; pd <- match(cell$dist[r], cell$keys); D <- cell$Dix
  folds <- group_folds71(cell$dist[r], 5, rep); pred <- p_null <- rep(NA_real_, length(y))
  for (f in sort(unique(folds))) {
    trp <- which(folds != f); tep <- which(folds == f); trd <- sort(unique(pd[trp]))
    yd <- vapply(trd, function(j) { i <- trp[pd[trp] == j]; stats::weighted.mean(y[i], w[i]) }, 0)
    if (type == "bin") yd <- .v2_logit(yd)
    z <- .index_weights_v2(D[trd, , drop = FALSE], yd); s <- as.numeric(D %*% z); s <- (s - mean(s[trd])) / max(stats::sd(s[trd]), 1e-9)
    dtr <- data.frame(y = y[trp], s = s[pd[trp]], w = w[trp] / mean(w[trp])); dte <- data.frame(s = s[pd[tep]])
    fit <- if (type == "bin") suppressWarnings(stats::glm(y ~ s, family = stats::binomial(), data = dtr, weights = w)) else stats::lm(y ~ s, data = dtr, weights = w)
    pred[tep] <- as.numeric(stats::predict(fit, newdata = dte, type = "response")); p_null[tep] <- mean(y[trp])
  }
  data.frame(type = type, outcome = outcome, set = "index", rep = rep, row = r, y = y, pred = pred, null = p_null, fold = folds, picks = "")
}

# the ceiling: other respondents in the same district, shrunk to the national value
run_ceiling71 <- function(type, outcome) {
  cell <- load_cell71(outcome); r <- rows71(cell, type); y <- (if (type == "bin") cell$ybin else cell$ycont)[r]; g <- cell$dist[r]
  s <- ave(y, g, FUN = sum); n <- ave(y, g, FUN = length); mu <- mean(y); gm <- tapply(y, g, mean); gn <- tapply(y, g, length)
  if (type == "bin") { tau2 <- max(stats::var(gm) - mean(mu * (1 - mu) / gn), 1e-6); k <- max(mu * (1 - mu) / tau2 - 1, 1)
  } else { sw2 <- sum(tapply(y, g, function(v) sum((v - mean(v))^2))) / (length(y) - length(gm)); tau2 <- max(stats::var(gm) - mean(sw2 / gn), 1e-6); k <- sw2 / tau2 }
  data.frame(type = type, outcome = outcome, set = "ceiling", rep = 0, row = r, y = y, pred = (s - y + k * mu) / (n - 1 + k),
             null = mu, fold = NA_integer_, picks = "", tau2 = tau2)
}
psu71 <- function(type, outcome) { cell <- load_cell71(outcome); cell$psu[rows71(cell, type)] }
dist71 <- function(type, outcome) { cell <- load_cell71(outcome); cell$dist[rows71(cell, type)] }
# The ceiling for ANY predictor that is constant within district: its skill cannot exceed the
# between-district share of person-level variance, tau^2 / total (method of moments).
mom_share71 <- function(y, g, type) {
  mu <- mean(y); gm <- tapply(y, g, mean); gn <- tapply(y, g, length)
  if (type == "bin") { vt <- mu * (1 - mu); tau2 <- stats::var(gm) - mean(mu * (1 - mu) / gn)
  } else { sw2 <- sum(tapply(y, g, function(v) sum((v - mean(v))^2))) / (length(y) - length(gm)); vt <- stats::var(y); tau2 <- stats::var(gm) - mean(sw2 / gn) }
  max(tau2, 0) / vt
}
