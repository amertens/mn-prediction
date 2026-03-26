# =============================================================================
# scripts/test_conceptual_ablation.R
#
# Test the conceptual domain permutation importance on Ghana child_iron
# (the best-performing country x outcome combination).
#
# Run interactively in RStudio -- NOT inside the targets pipeline.
# =============================================================================

library(targets)
library(ggplot2)
source(here::here("R", "conceptual_ablation.R"))

# ── Settings ────────────────────────────────────────────────────────────────
N_PERM      <- 50        # number of permutations per domain
RUN_LEVEL2  <- FALSE     # set TRUE to also run level-2 analysis

# ── 1. Load targets artifacts ────────────────────────────────────────────────

cat("Loading Ghana child_iron model...\n")
sl_fit       <- tar_read(sl_fit_ghana_child_iron)
outcome_data <- tar_read(outcome_data_ghana_child_iron)

# Outcome config for Ghana child_iron
oc <- list(
  tag          = "child_iron",
  label        = "Iron Deficiency (Children)",
  population   = "children",
  continuous   = "gw_cLogFer",
  binary       = "gw_cIDA",
  cutoff       = 12,
  cutoff_dir   = "less",
  cutoff_scale = "ug/L"
)

# ── 2. Check variable classification ────────────────────────────────────────

mapping <- load_conceptual_domains()
xvars   <- sl_fit$cont_fit$Xvars
cat("\nTotal model Xvars:", length(xvars), "\n")

# Level-1 classification
l1 <- classify_variables(xvars, mapping)
cat("\nLevel-1 domain distribution:\n")
print(sort(table(l1), decreasing = TRUE))
cat("\nUnknown variables:", sum(l1 == "Unknown"), "\n")
if (sum(l1 == "Unknown") > 0) {
  cat("Unknown vars:\n")
  cat(paste("  ", names(l1[l1 == "Unknown"]), collapse = "\n"), "\n")
}

# Level-2 classification (always show distribution even if not running analysis)
l2 <- classify_variables_l2(xvars, mapping)
cat("\nLevel-2 domain distribution:\n")
print(sort(table(l2), decreasing = TRUE))

# ── 3. Run Level-1 permutation importance (permute-one) ─────────────────────

cat("\n=== Running Level-1 conceptual permutation (permute-one) ===\n")
perm_l1 <- run_conceptual_permutation(
  sl_fit       = sl_fit,
  outcome_data = outcome_data,
  oc           = oc,
  level        = "level1",
  n_perm       = N_PERM,
  seed         = 12345
)

cat("\nLevel-1 permute-one results:\n")
print(perm_l1[order(-perm_l1[[grep("drop", colnames(perm_l1), value = TRUE)[1]]]), ])

# ── 4. Run Level-2 permutation importance (optional) ─────────────────────────

perm_l2 <- NULL
if (RUN_LEVEL2) {
  cat("\n=== Running Level-2 conceptual permutation ===\n")
  perm_l2 <- run_conceptual_permutation(
    sl_fit       = sl_fit,
    outcome_data = outcome_data,
    oc           = oc,
    level        = "level2",
    n_perm       = N_PERM,
    seed         = 12345
  )
  cat("\nLevel-2 results:\n")
  print(perm_l2[order(-perm_l2[[grep("drop", colnames(perm_l2), value = TRUE)[1]]]), ])
}

# ── 5. Run Level-1 "leave-one-in" importance ─────────────────────────────────

cat("\n=== Running Level-1 leave-one-in ===\n")
loi_l1 <- run_leave_one_in(
  sl_fit       = sl_fit,
  outcome_data = outcome_data,
  oc           = oc,
  level        = "level1",
  n_perm       = N_PERM,
  seed         = 12345
)

cat("\nLeave-one-in Level-1 results:\n")
sort_col <- if ("recovery_auc" %in% colnames(loi_l1)) "recovery_auc" else "recovery_r2"
print(loi_l1[order(-loi_l1[[sort_col]]), ])

# ── 6. Forest plots ──────────────────────────────────────────────────────────

p1 <- plot_conceptual_forest(
  perm_l1,
  title = "Ghana Child Iron — Level-1 Domain Importance (Permute-One)"
)
print(p1)

p3 <- plot_leave_one_in_forest(
  loi_l1,
  title = "Ghana Child Iron — Level-1 Standalone Contribution (Leave-One-In)"
)
print(p3)

if (RUN_LEVEL2 && !is.null(perm_l2)) {
  p2 <- plot_conceptual_forest(
    perm_l2,
    title = "Ghana Child Iron — Level-2 Domain Importance (Permute-One)"
  )
  print(p2)
}

# ── 7. Save results ──────────────────────────────────────────────────────────

results_dir <- here::here("results", "conceptual_ablation")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(perm_l1, file.path(results_dir, "ghana_child_iron_level1.csv"),
          row.names = FALSE)
write.csv(loi_l1, file.path(results_dir, "ghana_child_iron_leave_one_in_l1.csv"),
          row.names = FALSE)

ggsave(file.path(results_dir, "ghana_child_iron_forest_l1.png"), p1,
       width = 10, height = 8, dpi = 150)
ggsave(file.path(results_dir, "ghana_child_iron_forest_loi_l1.png"), p3,
       width = 10, height = 8, dpi = 150)

if (RUN_LEVEL2 && !is.null(perm_l2)) {
  write.csv(perm_l2, file.path(results_dir, "ghana_child_iron_level2.csv"),
            row.names = FALSE)
  ggsave(file.path(results_dir, "ghana_child_iron_forest_l2.png"), p2,
         width = 10, height = 10, dpi = 150)
}

cat("\nResults saved to:", results_dir, "\n")
cat("Done!\n")
