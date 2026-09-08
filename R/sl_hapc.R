# =============================================================================
# R/sl_hapc.R
#
# SuperLearner wrappers around the hapc Python package (Wang, Schuler,
# van der Laan, Garcia Meixide, arXiv 2602.10613: principal-component Highly
# Adaptive Ridge / Lasso). The HAL kernel is eigendecomposed and the outcome
# regressed on the leading scores; lambda is chosen by hapc's own 5-fold CV
# inside the training rows SuperLearner hands the wrapper, so the outer
# cross-validation stays honest.
#
#   SL.hapc        PCHAR (norm = "2", closed-form ridge in the PC basis)
#   SL.hapc_lasso  PCHAL (norm = "1", soft-thresholding)
#
# Requires the reticulate venv with hapc 2.6.0 installed
# (C:/Users/andre/OneDrive/Documents/.virtualenvs/r-reticulate). Survey
# weights are NOT used: hapc has no observation-weight argument, so the
# wrapper ignores obsWeights and says so once. Smoke test: HP-01 in
# docs/findings/SANDBOX_LOG_2026-09.md (scripts/protocol_v2/50, 51).
# =============================================================================
.hapc_python <- "C:/Users/andre/OneDrive/Documents/.virtualenvs/r-reticulate/Scripts/python.exe"
if (!nzchar(Sys.getenv("RETICULATE_PYTHON")) && file.exists(.hapc_python)) Sys.setenv(RETICULATE_PYTHON = .hapc_python)
.hapc_env <- new.env()
.hapc_mod <- function(what = "cv") {
  key <- paste0("mod_", what)
  if (!is.null(.hapc_env[[key]])) return(.hapc_env[[key]])
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("SL.hapc needs the reticulate package")
  m <- reticulate::import(paste0("hapc.", what), delay_load = FALSE)
  .hapc_env[[key]] <- m; m
}
.hapc_np <- function() { if (is.null(.hapc_env$np)) .hapc_env$np <- reticulate::import("numpy", delay_load = FALSE); .hapc_env$np }

SL.hapc <- function(Y, X, newX, family, obsWeights, id, norm = "2", max_degree = 1L,
                    log_lambda_min = -8, log_lambda_max = 2, grid_length = 12L, nfolds = 5L, ...) {
  if (!identical(family$family, "gaussian")) stop("SL.hapc: gaussian family only")
  if (!isTRUE(.hapc_env$warned) && !is.null(obsWeights) && length(unique(obsWeights)) > 1) {
    message("[SL.hapc] observation weights are ignored (hapc has no weights argument)"); .hapc_env$warned <- TRUE }
  cv <- .hapc_mod("cv"); np <- .hapc_np()
  Xtr <- as.matrix(X); Xte <- as.matrix(newX); storage.mode(Xtr) <- "double"; storage.mode(Xte) <- "double"
  res <- cv$cv_hapc(np$array(Xtr), np$array(as.numeric(Y)), family = "gaussian", max_degree = as.integer(max_degree),
                    npcs = NULL, norm = norm, nfolds = as.integer(nfolds), predict = np$array(Xte),
                    log_lambda_min = log_lambda_min, log_lambda_max = log_lambda_max, grid_length = as.integer(grid_length))
  pred <- as.numeric(reticulate::py_to_r(res$predictions))
  if (length(pred) != nrow(Xte) || !all(is.finite(pred))) stop("SL.hapc: hapc returned non-finite or mis-sized predictions")
  fit <- list(X = Xtr, Y = as.numeric(Y), norm = norm, max_degree = as.integer(max_degree),
              lambda = tryCatch(as.numeric(reticulate::py_to_r(res$best_lambda)), error = function(e) NA_real_))
  class(fit) <- "SL.hapc"
  list(pred = pred, fit = fit)
}
predict.SL.hapc <- function(object, newdata, ...) {
  single <- .hapc_mod("single"); np <- .hapc_np()
  Xte <- as.matrix(newdata); storage.mode(Xte) <- "double"
  res <- tryCatch(single$single_lambda_fit(np$array(object$X), np$array(object$Y), max_degree = object$max_degree,
                                           npcs = as.integer(nrow(object$X)), lambda_ = object$lambda, predict = np$array(Xte),
                                           l1 = identical(object$norm, "1")), error = function(e) NULL)
  if (is.null(res)) stop("SL.hapc: refit for prediction failed; predict via newX at fit time instead")
  as.numeric(reticulate::py_to_r(res$predictions))
}
SL.hapc_lasso <- function(...) SL.hapc(..., norm = "1")
predict.SL.hapc_lasso <- predict.SL.hapc
