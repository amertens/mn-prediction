


rm(list=ls())
library(dplyr)
library(tidyverse)
library(haven)
library(here)
library(purrr)
library(labelled)
library(sl3)
library(origami)
library(tlverse)
library(caret)
library(data.table)
library(ck37r)
library(rdhs)
library(maps)
library(SuperLearner)
library(washb)
library(future)
options(future.globals.maxSize = 5 * 1024^3)


source(paste0(here::here(),"/src/0-functions.R"))
source(paste0(here::here(),"/src/0-SL-setup.R"))
source(paste0(here::here(),"/src/DHS/DHS_functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_variable_recode.R"))

res = readRDS(here("results/models/res_GW_Ghana_SL_child_vitA_full.rds"))
res_gw = readRDS(here("results/models/res_GW_Ghana_SL_child_vitA_gwPreds.rds"))
res_nogw = readRDS(here("results/models/res_GW_Ghana_SL_child_vitA.rds"))
metadata= readRDS(here("metadata/variable_categories.rds"))
names(metadata$variable_sources)

sl_fit=res$sl_fit
eval_fun = loss_loglik_binomial
importance.metric="ratio"
n_vars=20

fit = sl_fit
type = "permute"
fold_number = "validation"
#importance_metric = c("difference", "ratio"),
importance_metric = "difference"

task <- fit$training_task
d <- task$data
X <- task$nodes$covariates
Y <- task$Y
weights <- task$weights


covariate_groups=list("MN survey" = metadata$variable_sources$gw_vars[metadata$variable_sources$gw_vars %in% X],
                      "DHS" = metadata$variable_sources$dhs_vars[metadata$variable_sources$dhs_vars %in% X],
                      "MICS" = metadata$variable_sources$mics_vars[metadata$variable_sources$mics_vars %in% X],
                      "IHME" = metadata$variable_sources$ihme_vars[metadata$variable_sources$ihme_vars %in% X],
                      "LSMS" = metadata$variable_sources$lsms_vars[metadata$variable_sources$lsms_vars %in% X],
                      "WFP" = metadata$variable_sources$wfp_vars[metadata$variable_sources$wfp_vars %in% X],
                      "FluNet" = metadata$variable_sources$flunet_vars[metadata$variable_sources$flunet_vars %in% X])


X <- res_nogw$sl_fit$training_task$nodes$covariates
covariate_groups_noGW=list("DHS" = metadata$variable_sources$dhs_vars[metadata$variable_sources$dhs_vars %in% X],
                           "MICS" = metadata$variable_sources$mics_vars[metadata$variable_sources$mics_vars %in% X],
                           "IHME" = metadata$variable_sources$ihme_vars[metadata$variable_sources$ihme_vars %in% X],
                           "LSMS" = metadata$variable_sources$lsms_vars[metadata$variable_sources$lsms_vars %in% X],
                           "WFP" = metadata$variable_sources$wfp_vars[metadata$variable_sources$wfp_vars %in% X],
                           "FluNet" = metadata$variable_sources$flunet_vars[metadata$variable_sources$flunet_vars %in% X])




calc_mn_importance <- function(sl_fit, eval_fun = loss_loglik_binomial,
                               importance.metric="ratio", n_vars=20,
                               covariate_groups=NULL){

  set.seed(983)

  fit = sl_fit
  type = "permute"
  fold_number = "validation"
  #importance_metric = c("difference", "ratio"),
  importance_metric = "difference"

  task <- fit$training_task
  d <- task$data
  X <- task$nodes$covariates
  Y <- task$Y
  weights <- task$weights
  if (!is.null(covariate_groups)) {
    if (!is.list(covariate_groups)) {
      stop("Covariate groups must be a list.")
    }
    if (!all(unlist(covariate_groups) %in% X)) {
      stop("Groups contain covariates that are not in the task's covariates.")
    }
    # if (!all(X %in% unlist(covariate_groups))) {
    #   missingX <- as.list(X[which(!X %in% unlist(covariate_groups))])
    #   names(missingX) <- X[which(!X %in% unlist(covariate_groups))]
    #   covariate_groups <- c(covariate_groups, missingX)
    # }
    if(any(is.null(names(covariate_groups))) | any(names(covariate_groups) =="") | any(is.na(names(covariate_groups)))){
      no_name <- unique(which(is.null(names(covariate_groups)) |names(covariate_groups) == "" | is.na(names(covariate_groups))))
      if (any(sapply(covariate_groups[unique(no_name)],length)) != 1) {
        stop("Covariate groups with more than one covariate must be named.")
      }else{
        if(all(sapply(covariate_groups[unique(no_name)],length)) == 1){
          names(covariate_groups[no_name]) <- unlist(covariate_groups[no_name])
        }
      }
    }
    X <- covariate_groups
  }else{
    names(X) <- X
  }

  original_pred <- fit$predict_fold(task, fold_number = fold_number)
  original_eval <- eval_fun(original_pred, Y)
  if (!is.null(attr(original_eval, "loss")) && !attr(original_eval,"loss")) {
    original_risk <- original_eval
  }else {
    original_losses <- original_eval
    original_risk <- weighted.mean(original_losses, weights)
  }

  res_list <- lapply(X, function(x){
    if (type == "permute"){
      perm <- sample(1:nrow(d), nrow(d))
      x_perm <- d[perm, x, with = FALSE]
      data.table::setnames(x_perm, x)
      x_perm_name <- task$add_columns(x_perm)
      task_x_perm <- task$next_in_chain(column_names = x_perm_name)
      x_perm_pred <- fit$predict_fold(task_x_perm, fold_number = fold_number)
      x_perm_eval <- eval_fun(x_perm_pred, Y)

      if (!is.null(attr(original_eval, "loss")) && !attr(original_eval,"loss")) {
        no_x_risk <- x_perm_eval
      }else{
        no_x_losses <- x_perm_eval
        no_x_risk <- weighted.mean(no_x_losses, weights)
      }
    }

    if (importance_metric == "ratio") {
      result <- no_x_risk/original_risk
    }
    else if (importance_metric == "difference") {
      result <- no_x_risk - original_risk
    }
    return(result)
  })

  result <- data.table::data.table(covariate = names(X), metric = unlist(res_list))
  if (!is.null(covariate_groups)) {
    colnames(result)[1] <- "covariate_group"
  }
  data.table::setorderv(result, cols = "metric", order = -1L)
  metric_name <- paste0("risk_", importance_metric)
  if (!is.null(attr(original_eval, "name"))) {
    metric_name <- gsub("risk", attr(original_eval, "name"),
                        metric_name)
  }
  colnames(result)[2] <- metric_name


  varimp20=result %>% arrange(.[,2]) %>% tail(n=n_vars)
  p=importance_plot(x = varimp20)
  return(list(varimp=result, varimp20=varimp20, p=p))
}










calc_mn_importance_par <- function(sl_fit, eval_fun = loss_loglik_binomial,
                                   importance.metric="ratio", n_vars=20,
                                   covariate_groups=NULL, n_cores=NULL){

  # Load required packages
  require(parallel)

  set.seed(983)
  fit = sl_fit
  type = "permute"
  fold_number = "validation"
  importance_metric = "difference"

  task <- fit$training_task
  d <- task$data
  X <- task$nodes$covariates
  Y <- task$Y
  weights <- task$weights

  # Set up number of cores (defaults to all available cores - 1)
  if(is.null(n_cores)) {
    n_cores <- max(1, detectCores() - 1)
  }

  if (!is.null(covariate_groups)) {
    if (!is.list(covariate_groups)) {
      stop("Covariate groups must be a list.")
    }
    if (!all(unlist(covariate_groups) %in% X)) {
      stop("Groups contain covariates that are not in the task's covariates.")
    }
    # if (!all(X %in% unlist(covariate_groups))) {
    #   missingX <- as.list(X[which(!X %in% unlist(covariate_groups))])
    #   names(missingX) <- X[which(!X %in% unlist(covariate_groups))]
    #   covariate_groups <- c(covariate_groups, missingX)
    # }
    if(any(is.null(names(covariate_groups))) | any(names(covariate_groups) ==
                                                   "") | any(is.na(names(covariate_groups)))) {
      no_name <- unique(which(is.null(names(covariate_groups)) |
                                names(covariate_groups) == "" | is.na(names(covariate_groups))))
      if (any(sapply(covariate_groups[unique(no_name)],
                     length)) != 1) {
        stop("Covariate groups with more than one covariate must be named.")
      }
      else if (all(sapply(covariate_groups[unique(no_name)],length)) == 1) {
        names(covariate_groups[no_name]) <- unlist(covariate_groups[no_name])
      }
    }
    X <- covariate_groups
  }else{
    names(X) <- X
  }

  original_pred <- fit$predict_fold(task, fold_number = fold_number)
  original_eval <- eval_fun(original_pred, Y)

  if (!is.null(attr(original_eval, "loss")) && !attr(original_eval,"loss")) {
    original_risk <- original_eval
  }else {
    original_losses <- original_eval
    original_risk <- weighted.mean(original_losses, weights)
  }

  # Create cluster for Windows
  cl <- makeCluster(n_cores, type = "PSOCK")

  # Export necessary objects to cluster workers
  clusterExport(cl, c("fit", "task", "d", "Y", "weights", "eval_fun",
                      "original_eval", "original_risk", "type",
                      "fold_number", "importance_metric"),
                envir = environment())

  # Load required packages on each worker
  clusterEvalQ(cl, {
    library(data.table)
    # Add any other required libraries here
  })

  # Define the function to be applied in parallel
  process_covariate <- function(x) {
    if (type == "permute") {
      perm <- sample(1:nrow(d), nrow(d))
      x_perm <- d[perm, x, with = FALSE]
      data.table::setnames(x_perm, x)
      x_perm_name <- task$add_columns(x_perm)
      task_x_perm <- task$next_in_chain(column_names = x_perm_name)
      x_perm_pred <- fit$predict_fold(task_x_perm, fold_number = fold_number)
      x_perm_eval <- eval_fun(x_perm_pred, Y)

      if (!is.null(attr(original_eval, "loss")) && !attr(original_eval, "loss")) {
        no_x_risk <- x_perm_eval
      } else {
        no_x_losses <- x_perm_eval
        no_x_risk <- weighted.mean(no_x_losses, weights)
      }
    }

    if (importance_metric == "ratio") {
      result <- no_x_risk/original_risk
    } else if (importance_metric == "difference") {
      result <- no_x_risk - original_risk
    }

    return(result)
  }

  # Run parallel computation
  tryCatch({
    res_list <- parLapply(cl, X, process_covariate)
  }, finally = {
    # Always stop the cluster
    stopCluster(cl)
  })

  result <- data.table::data.table(covariate = names(X), metric = unlist(res_list))

  if (!is.null(covariate_groups)) {
    colnames(result)[1] <- "covariate_group"
  }

  data.table::setorderv(result, cols = "metric", order = -1L)
  metric_name <- paste0("risk_", importance_metric)

  if (!is.null(attr(original_eval, "name"))) {
    metric_name <- gsub("risk", attr(original_eval, "name"), metric_name)
  }

  colnames(result)[2] <- metric_name


  varimp <- result  # placeholder - replace with actual importance() call
  varimp20 <- result %>% arrange(.[,2]) %>% tail(n=n_vars)
  p <- NULL  # placeholder - replace with actual importance_plot() call

  return(list(varimp=result, varimp20=varimp20, p=p))
}




res_vim_diff_data_sources=try(calc_mn_importance(res$sl_fit,  eval_fun = loss_squared_error, importance.metric="difference",
                                                 covariate_groups=covariate_groups))
res_vim_diff_data_sources
res_vim_diff_data_sources$p

# res_vim_diff_data_sources2=try(calc_mn_importance_par(res$sl_fit,  eval_fun = loss_squared_error, importance.metric="difference",
#                                                  covariate_groups=covariate_groups))
# res_vim_diff_data_sources2


start_time <- Sys.time()
res_vim_diff=try(calc_mn_importance(res_gw$sl_fit,  eval_fun = loss_squared_error, importance.metric="difference", covariate_groups=NULL))
res_vim_diff
end_time <- Sys.time()
execution_time <- end_time - start_time
execution_time

# res_vim_diff2=try(calc_mn_importance_par(res$sl_fit,  eval_fun = loss_squared_error, importance.metric="difference",covariate_groups=NULL))
# res_vim_diff2


start_time <- Sys.time()
res_vim_diff_noGW=try(calc_mn_importance(res_nogw$sl_fit,  eval_fun = loss_squared_error, importance.metric="difference", covariate_groups=NULL))
res_vim_diff_noGW
end_time <- Sys.time()
execution_time <- end_time - start_time
execution_time

res_vim_diff_data_sources_noGW=try(calc_mn_importance(res_nogw$sl_fit,  eval_fun = loss_squared_error, importance.metric="difference",
                                                      covariate_groups=covariate_groups_noGW))
res_vim_diff_data_sources_noGW
res_vim_diff_data_sources_noGW$p
