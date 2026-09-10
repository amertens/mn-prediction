# =============================================================================
# scripts/policy_deck/07_worst_fifth_probability.R
#
# HOW SURE IS THE MODEL THAT A DISTRICT IS IN THE WORST FIFTH?
#
# A ranking without uncertainty is hard to act on. The protocol already fits the
# zero-tuning index under replicated 5-fold cross-validation by district, so the
# uncertainty a ranking admits is available for free: refit the index on R
# different fold draws, rank each district's held-out prediction within its
# country on every draw, and record how often it lands in the worst fifth. A
# district that is in the worst fifth on 90 percent of draws is a firm priority;
# one that is there on 30 percent is a coin toss and should be treated as such.
#
# Every district's prediction is out-of-fold on every draw (the district is
# never in the training set that predicts it), so the probability is honest
# about in-fill uncertainty. It does not cover the survey's own noise in the
# district value, which the reliability-ceiling slide reports separately.
#
# Outputs
#   results/tables/policy_deck/worst_fifth_probability.csv   every country x outcome x district
#   results/tables/policy_deck/worst_fifth_calibration.csv    is p calibrated against the survey's own worst fifth?
#   results/figures/policy_deck/fig11_worst_fifth_probability.png   Ghana, child vitamin A and women's iron
#
#   Rscript scripts/policy_deck/07_worst_fifth_probability.R
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
P2 <- "results/tables/protocol_v2"; OUTT <- "results/tables/policy_deck"; OUTF <- "results/figures/policy_deck"
dir.create(OUTT, recursive = TRUE, showWarnings = FALSE); dir.create(OUTF, recursive = TRUE, showWarnings = FALSE)
R_DRAWS <- as.integer(Sys.getenv("WF_DRAWS", "40"))

TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv", stringsAsFactors = FALSE)
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")

rows <- list()
for (cn in COUNTRIES) for (on in unique(TG$outcome[TG$country == cn])) {
  t <- TG[TG$country == cn & TG$outcome == on & is.finite(TG$y_prev) & is.finite(TG$n_eff), ]
  if (nrow(t) < 15) next                              # in-fill needs enough districts for 5 folds of 3 or more
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS])); if (ncol(Xr) < 20) next
  Y <- .v2_logit(m$y_prev); aux <- list(Admin1 = m$Admin1, y_nat = Y)
  D <- domain_representation_v2(Xr, domain_of)      # the in-country basis, as in the protocol
  k_worst <- max(1L, ceiling(0.2 * nrow(m)))
  inworst <- matrix(NA, nrow(m), R_DRAWS); pred_all <- matrix(NA_real_, nrow(m), R_DRAWS)
  for (r in seq_len(R_DRAWS)) {
    folds <- make_folds_v2("kfold_district", nrow(m), k = 5, rep_id = r); pred <- rep(NA_real_, nrow(m))
    for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 12) next
      pred[te] <- ARMS_V2[["domain_index"]](tr, te, Y, NULL, D, aux) }
    pred_all[, r] <- pred
    ok <- is.finite(pred); worst <- order(pred, decreasing = TRUE)[seq_len(k_worst)]
    inworst[, r] <- FALSE; inworst[worst, r] <- TRUE; inworst[!ok, r] <- NA
  }
  p <- rowMeans(inworst, na.rm = TRUE)
  obs_worst <- rank(-m$y_prev, ties.method = "first") <= k_worst
  rows[[paste(cn, on)]] <- data.frame(country = cn, outcome = on, Admin1 = m$Admin1, Admin2 = m$Admin2, n_districts = nrow(m), k_worst = k_worst,
                                       obs_prev = m$y_prev, obs_worst_fifth = obs_worst, pred_mean = .v2_expit(rowMeans(pred_all, na.rm = TRUE)),
                                       p_worst_fifth = p, draws = R_DRAWS, stringsAsFactors = FALSE)
  cat(sprintf("%-12s %-13s districts %3d | p >= 0.8: %2d | p <= 0.2: %2d | of the survey's worst fifth, share with p >= 0.5: %.2f\n",
              cn, on, nrow(m), sum(p >= 0.8), sum(p <= 0.2), mean(p[obs_worst] >= 0.5)))
}
WF <- bind_rows(rows)
write.csv(WF, file.path(OUTT, "worst_fifth_probability.csv"), row.names = FALSE)

# calibration: within probability bands, how often is the district in the survey's own worst fifth?
CAL <- WF |> mutate(band = cut(p_worst_fifth, c(-0.01, 0.2, 0.5, 0.8, 1.0), labels = c("0 to 20%", "20 to 50%", "50 to 80%", "80 to 100%"))) |>
  group_by(band) |> summarise(districts = n(), share_in_survey_worst_fifth = mean(obs_worst_fifth), .groups = "drop")
write.csv(CAL, file.path(OUTT, "worst_fifth_calibration.csv"), row.names = FALSE)
cat("\ncalibration (all countries and outcomes):\n"); print(as.data.frame(CAL), row.names = FALSE)

# figure: Ghana, child vitamin A and women's iron; fill = probability, survey's worst fifth outlined
B <- readRDS("dashboard/data/admin2_boundaries.rds")[["ghana"]]
lab <- c(child_vitA = "Children: vitamin A deficiency", women_iron = "Women: iron deficiency")
g <- WF |> filter(country == "Ghana", outcome %in% names(lab)) |> mutate(panel = factor(lab[outcome], levels = lab))
gs <- dplyr::left_join(B, g, by = c("Admin1", "Admin2")) |> filter(!is.na(panel))
outline <- gs |> filter(obs_worst_fifth)
p <- ggplot(sf::st_as_sf(gs)) +
  geom_sf(aes(fill = p_worst_fifth), colour = "white", linewidth = 0.18) +
  geom_sf(data = sf::st_as_sf(outline), fill = NA, colour = "#1A1A1A", linewidth = 0.9) +
  facet_wrap(~ panel) +
  scale_fill_gradientn(colours = c("#F1F5F7", "#9CC8D3", "#0F7B8A", "#083D45"), limits = c(0, 1), labels = scales::percent, name = "Chance the district\nis in the worst fifth") +
  labs(caption = "Each district predicted with itself left out, 40 fold draws. Black outline: the survey's own worst fifth. Grey: no survey clusters.") +
  theme_void(base_size = 18) + theme(legend.position = "right", strip.text = element_text(face = "bold", size = 18), plot.caption = element_text(size = 12, colour = "grey30", hjust = 0))
ggsave(file.path(OUTF, "fig11_worst_fifth_probability.png"), p, width = 12.2, height = 6.6, dpi = 200, bg = "white")
cat("wrote fig11_worst_fifth_probability.png\n")
