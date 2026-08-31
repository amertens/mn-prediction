# =============================================================================
# scripts/accuracy_impact/wsb2_tournament.R
#
# WS-B2 and WS-B3. The shipping decision, made against simulated truth.
#
#   PROFILE=smoke   Ghana child_iron, 30 replicates
#   Rscript scripts/accuracy_impact/wsb2_tournament.R
# -> results/tables/estimator_tournament_truth.csv
# -> results/figures/estimator_tournament.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
PROFILE <- Sys.getenv("PROFILE", "full")
STORE <- here("_targets_full"); SUF <- if (PROFILE == "smoke") "_SMOKE" else ""
R_REP <- as.integer(Sys.getenv("WSB2_R", if (PROFILE == "smoke") "30" else "100"))
# The tournament costs about half an hour per cell: every replicate refits the
# ridge and the tilt once per region. Twenty-four cells is not affordable, so the
# full run is a named eight-cell subset spanning all four countries and both
# nutrient families, and the scope is stated rather than implied.
SUBSET <- list(c("Ghana","child_iron"), c("Ghana","women_vitA"),
               c("Gambia","child_iron"), c("Gambia","women_iron"),
               c("Malawi","child_vitA"), c("Malawi","women_folate"),
               c("SierraLeone","child_iron"), c("SierraLeone","women_iron"))
kk <- function(x) tolower(gsub("[^a-z]", "", tolower(x)))
RHO <- c(0, 0.2, 0.35, 0.6); SEED <- 20260922L
TDIR <- here("results","tables"); FDIR <- here("results","figures")
dir.create(FDIR, showWarnings = FALSE, recursive = TRUE)

H <- suppressMessages(readr::read_csv(
  here("data","covariates","harmonized","predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS <- setdiff(names(H), c("country","Admin1","Admin2"))
cfgs <- get_country_configs()
if (PROFILE == "smoke") cfgs <- cfgs["Ghana"]

rows <- list()
for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; hc <- H[H$country == cn, , drop = FALSE]
  if (!nrow(hc)) next
  ocs <- names(cc$outcomes)
  if (PROFILE == "smoke") ocs <- intersect(ocs, "child_iron") else
    ocs <- ocs[vapply(ocs, function(o) any(vapply(SUBSET,
      function(z) kk(z[1]) == kk(cn) && z[2] == o, logical(1))), logical(1))]
  for (ocn in ocs) {
    oc <- cc$outcomes[[ocn]]
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", tolower(cn), "_", ocn),
                                         store = STORE), error = function(e) NULL)
    if (is.null(od) || is.null(oc$binary)) next
    d <- od$data
    need <- c(cc$admin1_col, cc$admin2_col, cc$cluster_id, cc$weight_col, oc$binary)
    if (!all(need %in% names(d))) next
    # district covariates and the region map, aligned on the pair key that the
    # individual data uses for its Admin2 column
    key <- unique(data.frame(Admin1 = trimws(as.character(d[[cc$admin1_col]])),
                             Admin2 = trimws(as.character(d[[cc$admin2_col]])),
                             stringsAsFactors = FALSE))
    key <- key[!is.na(key$Admin1) & !is.na(key$Admin2), , drop = FALSE]
    m <- dplyr::inner_join(key, hc, by = admin2_join_by(key, hc))
    if (nrow(m) < 12) next
    X <- as.matrix(m[, COVS, drop = FALSE])
    keepc <- apply(X, 2, function(z) all(is.finite(z)) && stats::sd(z) > 0)
    X <- X[, keepc, drop = FALSE]
    if (ncol(X) < 5) next
    rownames(X) <- m$Admin2
    reg <- stats::setNames(m$Admin1, m$Admin2)
    r <- tryCatch(run_tournament(d, cc$admin2_col, cc$cluster_id, cc$weight_col,
                                 oc$binary, X, reg, rho = RHO, R = R_REP,
                                 seed = SEED),
                  error = function(e) { message("  ", cn, " ", ocn, ": ",
                                                conditionMessage(e)); NULL })
    if (is.null(r)) next
    r$country <- cc$country; r$outcome <- ocn
    rows[[length(rows)+1L]] <- r
    cat(sprintf("  [ok] %-13s %-13s districts=%3d covariates=%3d\n",
                cn, ocn, nrow(m), ncol(X)))
  }
}
res <- dplyr::bind_rows(rows)
if (!nrow(res)) stop("No tournament rows.")
front <- c("country","outcome","estimator","rho")
res <- res[, c(front, setdiff(names(res), front))]
readr::write_csv(res, file.path(TDIR, sprintf("estimator_tournament_truth%s.csv", SUF)))

cat("\n=== WS-B2: scored against TRUTH, pooled over cells ===\n")
s <- res |> group_by(rho, estimator) |>
  summarise(cells = dplyr::n(),
            r_truth = round(mean(r_truth, na.rm=TRUE), 3),
            mae_truth = round(mean(mae_truth, na.rm=TRUE), 2),
            bias_truth = round(mean(bias_truth, na.rm=TRUE), 2),
            r_obs = round(mean(r_obs, na.rm=TRUE), 3), .groups="drop") |>
  arrange(rho, desc(r_truth))
print(as.data.frame(s), row.names = FALSE)

cat("\n=== winner by rho, on r against truth ===\n")
w <- s |> group_by(rho) |> slice_max(r_truth, n = 1) |> ungroup()
print(as.data.frame(w[, c("rho","estimator","r_truth","mae_truth")]), row.names=FALSE)
cat("\n=== winner by rho, on MAE against truth ===\n")
w2 <- s |> group_by(rho) |> slice_min(mae_truth, n = 1) |> ungroup()
print(as.data.frame(w2[, c("rho","estimator","mae_truth","r_truth")]), row.names=FALSE)

cat("\n=== does the EB blend beat the flat regional mean against truth? ===\n")
cmp <- res |> filter(estimator %in% c("eb_blend","flat_region","direct")) |>
  select(country, outcome, rho, estimator, r_truth, mae_truth) |>
  tidyr::pivot_wider(names_from = estimator, values_from = c(r_truth, mae_truth))
if (all(c("r_truth_eb_blend","r_truth_flat_region") %in% names(cmp)))
  print(as.data.frame(cmp |> group_by(rho) |> summarise(
    cells = dplyr::n(),
    eb_beats_flat_r = sprintf("%d/%d", sum(r_truth_eb_blend > r_truth_flat_region, na.rm=TRUE), dplyr::n()),
    eb_beats_direct_r = sprintf("%d/%d", sum(r_truth_eb_blend > r_truth_direct, na.rm=TRUE), dplyr::n()),
    eb_beats_flat_mae = sprintf("%d/%d", sum(mae_truth_eb_blend < mae_truth_flat_region, na.rm=TRUE), dplyr::n()),
    .groups="drop")), row.names = FALSE)

png(file.path(FDIR, sprintf("estimator_tournament%s.png", SUF)),
    width = 1200, height = 800, res = 130)
op <- par(mfrow = c(1,2), mar = c(4.5,4.5,3,1))
ests <- unique(s$estimator); cols <- grDevices::hcl.colors(length(ests), "Dark 3")
for (metric in c("r_truth","mae_truth")) {
  plot(range(RHO), range(s[[metric]], na.rm=TRUE), type="n", xlab="rho (covariate signal share)",
       ylab = metric, main = paste("against truth:", metric))
  for (i in seq_along(ests)) { z <- s[s$estimator==ests[i],]; z <- z[order(z$rho),]
    lines(z$rho, z[[metric]], col=cols[i], lwd=2, type="b", pch=19) }
  if (metric == "r_truth") legend("topleft", bty="n", lwd=2, col=cols, legend=ests, cex=0.8)
}
par(op); dev.off()
cat(sprintf("\n-> %s\n", file.path("results","figures", sprintf("estimator_tournament%s.png", SUF))))
