# =============================================================================
# sandbox/09_h6_with_dag.R
#
# Experimental: causal-DAG prior before the H6 filter. Use the outcome-
# specific DAG (.dag_iron, .dag_vita, .dag_zinc, .dag_b12, .dag_folate
# in R/benchmark_models.R) to restrict the candidate pool to DAG-implied
# causal predictor groups, THEN apply the H6 composite filter (within-
# country invariance + outcome-shared score).
#
# Rationale: DAG selection encodes biological knowledge; H6 filter
# encodes statistical transportability. Their intersection should be
# both interpretable AND transportable.
#
# Not promoted to the main pipeline: this stacks two filters and
# requires the per-outcome DAG which is currently a starter template.
# =============================================================================

source(here::here("sandbox", "00_setup.R"))
source(here::here("R", "benchmark_models.R"))

pooled_all <- load_all_pooled()
cross_pooled <- setNames(lapply(pooled_all, function(x) x$pooled_data),
                          names(pooled_all))

make_dag_h6_predict <- function(this_outcome) {
  function(train, test, vars) {
    # Step 1: DAG ancestor set → covariates.
    adj_nodes <- .dag_adjustment_set(this_outcome)
    if (length(adj_nodes) == 0) {
      # No DAG defined for this outcome — fall through to H6 alone.
      restricted <- vars
    } else {
      dag_covars <- unique(unlist(lapply(adj_nodes, .match_dag_node, vars)))
      if (length(dag_covars) < 4) {
        cat(sprintf("    [%s] DAG mapped to %d vars (<4), falling back to all\n",
                    this_outcome, length(dag_covars)))
        restricted <- vars
      } else {
        cat(sprintf("    [%s] DAG restricted: %d/%d candidates\n",
                    this_outcome, length(dag_covars), length(vars)))
        restricted <- dag_covars
      }
    }
    # Step 2: H6 combined filter on the restricted pool.
    out <- fit_predict_combined_filter(
      train, test, restricted,
      cross_outcome_pooled = cross_pooled,
      this_outcome = this_outcome,
      min_within_ratio = 0.50,   # slightly looser because DAG already restricts
      top_k = 8L)
    if (is.null(out)) NULL else out$pred
  }
}

res <- list()
for (oc in OUTCOMES) {
  cat(sprintf("\n[H6+DAG] LOCO on %s ...\n", oc))
  r <- loco_evaluate(pooled_all[[oc]], make_dag_h6_predict(oc),
                      outcome = oc, method_name = "h6_plus_dag")
  res[[oc]] <- r
}
out <- dplyr::bind_rows(res)
print_loco(out)
readr::write_csv(out, file.path(SANDBOX_RESULTS, "09_h6_with_dag.csv"))
cat(sprintf("\nwritten: %s\n", file.path(SANDBOX_RESULTS, "09_h6_with_dag.csv")))
