# Rebuild dashboard/data/admin2_predictions.rds from the targets store, using
# the DISTRICT-level conformal intervals (conformal_ci$admin2_ci). Mirrors
# section 1 of 01_prepare_dashboard_data.R but avoids the heavy GADM/population
# steps. Run after the conformal-CI fix to restore SuperLearner intervals.
suppressWarnings(suppressMessages({ library(dplyr) }))
root <- "C:/Users/andre/OneDrive/Documents/mn-prediction"
if (!dir.exists(file.path(root, "_targets_full"))) root <- getwd()
S    <- file.path(root, "_targets_full", "objects")
DASH <- file.path(root, "dashboard", "data")
safe_read <- function(nm){ p<-file.path(S,nm); if(file.exists(p)) tryCatch(readRDS(p),error=function(e)NULL) else NULL }

countries <- c(gambia="Gambia", ghana="Ghana", sierraleone="Sierra Leone", malawi="Malawi")
outcome_labels <- c(child_vitA=NA, women_vitA=NA, child_iron=NA, women_iron=NA,
                    women_folate=NA, women_b12=NA, child_zinc=NA, women_zinc=NA)
who_thresholds <- list(
  child_vitA=c(none=0.02,mild=0.10,moderate=0.20), women_vitA=c(none=0.02,mild=0.10,moderate=0.20),
  child_iron=c(none=0.05,mild=0.20,moderate=0.40), women_iron=c(none=0.05,mild=0.20,moderate=0.40),
  women_folate=c(none=0.05,mild=0.20,moderate=0.40), women_b12=c(none=0.05,mild=0.20,moderate=0.40),
  child_zinc=c(none=0.10,mild=0.20,moderate=0.30), women_zinc=c(none=0.10,mild=0.20,moderate=0.30))
classify_who <- function(prev, oc){ th<-who_thresholds[[oc]]; if(is.null(th)||is.na(prev)) return(NA_character_)
  if(prev<th["none"])return("Low"); if(prev<th["mild"])return("Mild"); if(prev<th["moderate"])return("Moderate"); "Severe" }

rows <- list(); n_ci <- 0
for (low in names(countries)) {
  for (oc in names(outcome_labels)) {
    suf <- paste0(low, "_", oc)
    sl <- safe_read(paste0("admin2_sl_", suf)); if (is.null(sl) || nrow(sl)==0) next
    if (!"Admin2" %in% colnames(sl)) next
    pred_col <- intersect(c("sl_prev","pred_prev","predicted_prev","yhat","fit"), colnames(sl))[1]
    if (is.na(pred_col)) next
    df <- data.frame(country=unname(countries[low]), outcome=oc,
                     Admin1 = if ("Admin1" %in% colnames(sl)) sl$Admin1 else NA_character_,
                     Admin2 = sl$Admin2, pred_prev = sl[[pred_col]], stringsAsFactors=FALSE)
    svy <- safe_read(paste0("svy_admin2_", suf))
    if (!is.null(svy) && nrow(svy)>0) {
      pc <- intersect(c("svy_prev","obs_prev","prev","prevalence"), colnames(svy))[1]
      nc <- intersect(c("n_svy","n","n_obs","n_indiv"), colnames(svy))[1]
      df <- left_join(df, data.frame(Admin2=svy$Admin2,
        obs_prev = if(!is.na(pc)) svy[[pc]] else NA_real_,
        n_survey = if(!is.na(nc)) svy[[nc]] else NA_integer_), by="Admin2")
    } else { df$obs_prev <- NA_real_; df$n_survey <- NA_integer_ }
    df$ci_lo <- NA_real_; df$ci_hi <- NA_real_; df$ci_width <- NA_real_
    cc <- safe_read(paste0("conformal_ci_", suf))
    if (!is.null(cc) && !is.null(cc$admin2_ci) && nrow(cc$admin2_ci)>0 &&
        "Admin2" %in% colnames(cc$admin2_ci)) {
      a2 <- cc$admin2_ci
      k <- data.frame(Admin2=a2$Admin2, ci_lo=a2$ci_lo, ci_hi=a2$ci_hi, stringsAsFactors=FALSE)
      k$ci_width <- k$ci_hi - k$ci_lo
      df <- df |> select(-ci_lo,-ci_hi,-ci_width) |> left_join(k, by="Admin2")
      n_ci <- n_ci + 1
    }
    df$who_class <- vapply(df$pred_prev, classify_who, character(1), oc=oc)
    rows[[suf]] <- df
  }
}
out <- bind_rows(rows)
saveRDS(out, file.path(DASH, "admin2_predictions.rds"))
cat(sprintf("admin2_predictions.rds: %d rows; combos with district CIs: %d\n", nrow(out), n_ci))
print(tapply(out$ci_lo, out$country, function(x) sum(!is.na(x))))
