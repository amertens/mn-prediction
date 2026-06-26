# =============================================================================
# sandbox/iron_reharmonization_loco.R
#
# PROTOTYPE — test whether re-harmonizing Gambia's child-iron biomarker onto the
# other countries' adjustment/scale improves leave-one-country-out (LOCO)
# transport. Hold out Gambia; train an area-level model on Ghana+Malawi+SL;
# predict Gambia's admin-2 iron-deficiency prevalence; compare predictions
# against Gambia's prevalence computed three ways:
#   1. CURRENT     — gw_LogFerAdj (log-BRINDA), the pipeline's current choice.
#   2. THURNHAM    — gw_cFerrAdjThurn < 12 (same adjustment family as Ghana).
#   3. LEVEL-MATCH — Gambia raw ferritin rescaled so its median equals the pooled
#                    training raw-ferritin median (illustrative: removes the
#                    cross-survey LEVEL offset while preserving within-Gambia
#                    ranking). NOT a real fix — may erase true biology — but it
#                    upper-bounds how much a pure level correction could help.
#
# Diagnostic question: is the LOCO failure an ADJUSTMENT-method artifact (then #2
# helps) or a raw cross-survey LEVEL offset (then only #3 helps, and it is not
# fixable by adjustment harmonization)?
# =============================================================================
suppressWarnings(suppressMessages({
  library(here); library(dplyr); library(targets)
  source(here::here("R","config.R")); source(here::here("R","data_prep.R"))
  source(here::here("R","benchmark_models.R"))
}))
STORE <- "_targets_full"
cfgs  <- get_country_configs()

# ── Gambia admin-2 prevalence under the three definitions ────────────────────
g  <- load_merged_data(cfgs$Gambia$data_path)
gc <- g[which(g$gw_child_flag == 1), , drop = FALSE]
gc$.w  <- suppressWarnings(as.numeric(gc$gw_svy_weight)); gc$.w[!is.finite(gc$.w)] <- 1
gc$.a2 <- gc$Admin2

# pooled TRAIN raw-ferritin median (children): Ghana gw_fer, Malawi fer, SL gw_cFerritin
raw_tr <- c(
  suppressWarnings(as.numeric(load_merged_data(cfgs$Ghana$data_path) |>
      (\(d) d[which(d$gw_child_flag==1),"gw_fer"])() )),
  { m <- load_merged_data(cfgs$Malawi$data_path);
    suppressWarnings(as.numeric(m[which(m$population=="preschool children"),"fer"])) },
  { s <- load_merged_data(cfgs$SierraLeone$data_path);
    suppressWarnings(as.numeric(s[which(s$gw_child_flag==1),"gw_cFerritin"])) }
)
raw_tr <- raw_tr[is.finite(raw_tr)]
train_med <- median(raw_tr)
gam_med   <- median(suppressWarnings(as.numeric(gc$gw_fer)), na.rm = TRUE)
cat(sprintf("pooled train raw-ferritin median = %.1f ug/L; Gambia = %.1f; rescale factor = %.2f\n",
            train_med, gam_med, train_med/gam_med))

def_thurn <- as.numeric(suppressWarnings(as.numeric(gc$gw_cFerrAdjThurn)) < 12)
fer_std   <- suppressWarnings(as.numeric(gc$gw_fer)) * (train_med / gam_med)
def_std   <- as.numeric(fer_std < 12)
def_curr  <- as.numeric(exp(suppressWarnings(as.numeric(gc$gw_LogFerAdj))) < 12)

agg <- function(def) {
  ok <- is.finite(def) & is.finite(gc$.w)
  d <- data.frame(Admin2 = gc$.a2[ok], def = def[ok], w = gc$.w[ok])
  d |> group_by(Admin2) |>
    summarise(prev = weighted.mean(def, w), n = n(), .groups = "drop")
}
gam_curr  <- agg(def_curr)  |> rename(prev_current  = prev)
gam_thurn <- agg(def_thurn) |> rename(prev_thurnham = prev) |> select(Admin2, prev_thurnham)
gam_std   <- agg(def_std)   |> rename(prev_levelmatch = prev) |> select(Admin2, prev_levelmatch)
gam <- gam_curr |> left_join(gam_thurn, "Admin2") |> left_join(gam_std, "Admin2")
cat(sprintf("\nGambia admin-2 mean prevalence: current=%.1f%%  thurnham=%.1f%%  level-match=%.1f%%\n",
            mean(gam$prev_current)*100, mean(gam$prev_thurnham)*100, mean(gam$prev_levelmatch)*100))

# ── Training set: Ghana + Malawi + SL admin-2 (existing svy targets + gee) ────
tr <- c(ghana="Ghana", malawi="Malawi", sierraleone="Sierra Leone")
gee <- lapply(c("gambia", names(tr)), function(k) tar_read_raw(paste0("gee_admin2_", k), store=STORE))
names(gee) <- c("gambia", names(tr))
common <- Reduce(intersect, lapply(gee, function(d) grep("^gee_", names(d), value=TRUE)))
cat(sprintf("common GEE covariates across 4 countries: %d\n", length(common)))

train_rows <- lapply(names(tr), function(k){
  svy <- tar_read_raw(paste0("svy_admin2_", k, "_child_iron"), store=STORE)
  dplyr::inner_join(svy[,c("Admin2","svy_prev","n_svy")], gee[[k]][,c("Admin2",common)], "Admin2")
})
train <- dplyr::bind_rows(train_rows)
test  <- dplyr::inner_join(gam, gee[["gambia"]][,c("Admin2",common)], "Admin2")
cat(sprintf("train areas=%d (3 countries), Gambia test areas=%d\n", nrow(train), nrow(test)))

# ── Fit penalized area model on training countries, predict Gambia ───────────
suppressMessages(library(glmnet))
Xtr <- as.matrix(train[, common]); Xte <- as.matrix(test[, common])
med <- apply(Xtr, 2, function(z) { m <- median(z[is.finite(z)], na.rm=TRUE); if (is.finite(m)) m else 0 })
for (j in seq_along(common)) {
  zt <- Xtr[,j]; zt[!is.finite(zt)] <- med[j]; Xtr[,j] <- zt
  ze <- Xte[,j]; ze[!is.finite(ze)] <- med[j]; Xte[,j] <- ze
}
keep <- which(apply(Xtr, 2, function(z) sd(z) > 0))
cvfit <- cv.glmnet(Xtr[,keep,drop=FALSE], train$svy_prev, alpha = 0.5,
                   weights = pmax(train$n_svy, 1))
test$pred <- pmin(pmax(as.numeric(predict(cvfit, Xte[,keep,drop=FALSE], s = "lambda.min")), 0), 1)

metr <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  data.frame(
    r        = round(cor(truth[ok], pred[ok]), 3),
    mae_pp   = round(mean(abs(truth[ok]-pred[ok]))*100, 1),
    bias_pp  = round(mean(pred[ok]-truth[ok])*100, 1),
    truth_mean_pp = round(mean(truth[ok])*100,1),
    pred_mean_pp  = round(mean(pred[ok])*100,1))
}
out <- rbind(
  cbind(target="1. CURRENT (log-BRINDA)",  metr(test$prev_current,    test$pred)),
  cbind(target="2. THURNHAM (unify adj.)", metr(test$prev_thurnham,   test$pred)),
  cbind(target="3. LEVEL-MATCH (illustr.)",metr(test$prev_levelmatch, test$pred)))
cat("\n=== Gambia held-out LOCO: predictions vs target under 3 harmonizations ===\n")
print(out, row.names = FALSE)
