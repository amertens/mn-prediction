# =============================================================================
# shap_kernel/run_shap_kernel.R  — Option B driver (checkpointed, resumable)
#
# Computes model-agnostic SHAP per country x outcome, one checkpoint RDS per
# slice in results/shap_kernel/. Re-running SKIPS completed slices, so it is
# resume-on-restart safe. After all slices, combines into
# results/tables/shap_district_factors.csv + shap_global_importance.csv and
# rebuilds dashboard/data/importance.rds (keeping varimp + ablation).
#
# Env vars:
#   SHAP_SMOKE=1            tiny settings (1 slice quick sanity)
#   SHAP_ONLY=ghana_women_iron   run a single slice only
#   SHAP_NSIM / SHAP_TOPK / SHAP_NPD / SHAP_NMAX / SHAP_BG  override defaults
# =============================================================================
suppressWarnings(suppressMessages({ library(recipes) }))
source("shap_kernel/shap_kernel_lib.R")
source("R/config.R")

OUT <- file.path("results", "shap_kernel")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
smoke <- Sys.getenv("SHAP_SMOKE", "") == "1"
only  <- Sys.getenv("SHAP_ONLY", "")
geti  <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) as.integer(v) else d }
P <- list(nsim = geti("SHAP_NSIM", 20L), top_k = geti("SHAP_TOPK", 20L),
          n_per_district = geti("SHAP_NPD", 8L), n_explain_max = geti("SHAP_NMAX", 400L),
          bg_size = geti("SHAP_BG", 100L))

cfgs <- get_country_configs()
slices <- list()
for (clab in names(cfgs)) {
  low <- tolower(clab)
  for (oc_name in names(cfgs[[clab]]$outcomes)) {
    oc <- cfgs[[clab]]$outcomes[[oc_name]]
    tag <- oc$tag %||% oc_name
    if (!file.exists(file.path(PROD_STORE, "objects", paste0("sl_fit_", low, "_", tag)))) next
    suf <- paste0(low, "_", tag)
    if (nzchar(only) && suf != only) next
    slices[[suf]] <- list(low = low, cc = cfgs[[clab]], oc = oc)
  }
}
cat(sprintf("[shap] %d slices to consider | smoke=%s | nsim=%d top_k=%d npd=%d nmax=%d bg=%d\n",
            length(slices), smoke, P$nsim, P$top_k, P$n_per_district, P$n_explain_max, P$bg_size))

# ── per-slice compute (skip completed) ──────────────────────────────────────
for (suf in names(slices)) {
  ck <- file.path(OUT, paste0("shap_", suf, ".rds"))
  if (file.exists(ck)) { cat(sprintf("  [skip] %s (checkpoint exists)\n", suf)); next }
  s <- slices[[suf]]
  t0 <- proc.time()
  res <- tryCatch(
    compute_shap_slice(s$low, s$cc, s$oc, top_k = P$top_k, n_per_district = P$n_per_district,
                       n_explain_max = P$n_explain_max, bg_size = P$bg_size,
                       nsim = P$nsim, smoke = smoke),
    error = function(e) list(error = conditionMessage(e)))
  saveRDS(res, ck)
  el <- round((proc.time() - t0)["elapsed"] / 60, 1)
  if (!is.null(res$error)) cat(sprintf("  [ERR ] %-26s %s (%.1f min)\n", suf, res$error, el))
  else cat(sprintf("  [ok  ] %-26s districts=%d feats=%d explained=%d (%.1f min)\n",
                   suf, res$n_districts %||% 0, res$n_features %||% 0,
                   res$n_explained %||% 0, el))
}

# ── combine + rebuild dashboard bundle ──────────────────────────────────────
cat("[shap] combining checkpoints...\n")
cks <- list.files(OUT, pattern = "^shap_.*\\.rds$", full.names = TRUE)
df_rows <- list(); gi_rows <- list()
for (f in cks) {
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r) || !is.null(r$error)) next
  if (!is.null(r$district_factors)) df_rows[[f]] <- r$district_factors
  if (!is.null(r$global_importance)) gi_rows[[f]] <- r$global_importance
}
shap_factors <- if (length(df_rows)) do.call(rbind, df_rows) else NULL
shap_global  <- if (length(gi_rows)) do.call(rbind, gi_rows) else NULL

if (!is.null(shap_factors)) {
  utils::write.csv(shap_factors, file.path("results", "tables", "shap_district_factors.csv"),
                   row.names = FALSE)
  utils::write.csv(shap_global, file.path("results", "tables", "shap_global_importance.csv"),
                   row.names = FALSE)
  cat(sprintf("[shap] wrote shap_district_factors.csv (%d rows, %d slices)\n",
              nrow(shap_factors), length(df_rows)))

  # rebuild dashboard/data/importance.rds, preserving varimp + ablation
  rd <- function(f) { p <- file.path("results", "tables", f)
                      if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL }
  saveRDS(list(shap = shap_factors, varimp = rd("single_var_importance.csv"),
               ablation = rd("domain_ablation_all.csv")),
          file.path("dashboard", "data", "importance.rds"))
  cat("[shap] rebuilt dashboard/data/importance.rds — District Factors will populate.\n")
} else {
  cat("[shap] no successful slices yet — nothing to combine.\n")
}
cat("[shap] done.\n")
