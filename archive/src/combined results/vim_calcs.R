


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

res = readRDS(here("results/models/res_GW_Ghana_SL_mom_iron_v2.rds"))


metadata= readRDS(here("metadata/ghana_variable_categories.rds"))
dhs_indicators = read.csv( here("data/DHS/clean/dhs_indicators_metadata.csv"))
dhs_var_groups = unique(dhs_indicators$Level1)

names(metadata)

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

X = res$task$nodes$covariates

covariate_groups=list(#"MN survey" = metadata$variable_sources$gw_vars[metadata$variable_sources$gw_vars %in% X],
  "Geography" =c("Admin1" ),
  "DHS" = metadata$dhs_vars[metadata$dhs_vars %in% X],
  "MICS" = metadata$mics_vars[metadata$mics_vars %in% X],
  "IHME" = metadata$ihme_vars[metadata$ihme_vars %in% X],
  "LSMS" = metadata$lsms_vars[metadata$lsms_vars %in% X],
  "MAP" = metadata$map_vars[metadata$map_vars %in% X],
  "WFP" = metadata$wfp_vars[metadata$wfp_vars %in% X],
  "GEE" = metadata$gee_vars[metadata$gee_vars %in% X],
  "FluNet" = metadata$flunet_vars[metadata$flunet_vars %in% X])

#get DHS variables
dhs_var_list = list()
for(i in 1:length(dhs_var_groups)){
  dhs_var_group=dhs_var_groups[i]
  vars = dhs_indicators$IndicatorId[dhs_indicators$Level1==dhs_var_group]

  temp = NULL

  for(j in vars){

    temp =c(temp, X[grepl(j, X)])

  }
  dhs_var_list[[i]]=temp
  names(dhs_var_list)[i] = paste0("DHS: ", dhs_var_group)
}

dhs_var_list


variable_domains_nogw <- c(dhs_var_list, list(

  # ============================================================
  # 2. MICS
  # ============================================================


  # ---- Household composition & assets
  mics_household = grep("^mics_hh", metadata$mics_vars, value = TRUE),

  # ---- Child health & illness
  mics_child_health = grep("^mics_hc", metadata$mics_vars, value = TRUE),

  # ---- Water & sanitation
  mics_wash = grep("^mics_ws", metadata$mics_vars, value = TRUE),

  # ---- Environmental contamination (e.g. E. coli)
  mics_environment = grep("^mics_ec", metadata$mics_vars, value = TRUE),

  # ---- Wealth indices
  mics_wealth = grep("windex|wscore", metadata$mics_vars, value = TRUE),

  # ---- MICS other
  mics_other = setdiff(
    metadata$mics_vars,
    grep("^(mics_hh|mics_hc|mics_ws|mics_ec|windex|wscore)",
         metadata$mics_vars, value = TRUE)
  ),


  # ============================================================
  # 3. LSMS
  # ============================================================

  # ---- Food consumption & diet
  lsms_food_consumption = grep(
    "fd_|food|hhds",
    metadata$lsms_vars, value = TRUE
  ),

  # ---- Household income & labor
  lsms_income = grep(
    "inc|wage|aginc|rent|remit",
    metadata$lsms_vars, value = TRUE
  ),

  # ---- Poverty & welfare
  lsms_poverty = grep(
    "pov|pl_|welfare|quintile|decile",
    metadata$lsms_vars, value = TRUE
  ),

  # ---- Assets & housing
  lsms_assets_housing = grep(
    "house|rent|fuel|water|sewer|elec",
    metadata$lsms_vars, value = TRUE
  ),

  # ---- Education expenditure
  lsms_education = grep(
    "educ",
    metadata$lsms_vars, value = TRUE
  ),

  # ---- LSMS other
  lsms_other = setdiff(
    metadata$lsms_vars,
    grep("(fd_|food|hhds|inc|wage|aginc|rent|remit|pov|pl_|welfare|educ|house|fuel|water|sewer|elec)",
         metadata$lsms_vars, value = TRUE)
  ),


  # ============================================================
  # 4. IHME
  # ============================================================



  # ---- HIV prevalence & PLHIV
  ihme_hiv = grep("hiv", metadata$ihme_vars, value = TRUE),

  # ---- Child mortality
  ihme_child_mortality = grep("death|mortality", metadata$ihme_vars, value = TRUE),

  # ---- Diarrhea burden
  ihme_diarrhea = grep("diarrhea", metadata$ihme_vars, value = TRUE),

  # ---- Malaria prevalence & incidence
  ihme_malaria = grep("malaria", metadata$ihme_vars, value = TRUE),

  # ---- Undernutrition (stunting, wasting, anemia)
  ihme_nutrition = grep(
    "stunting|wasting|underweight|anemia",
    metadata$ihme_vars, value = TRUE
  ),

  # ---- Education attainment
  ihme_education = grep("education|school", metadata$ihme_vars, value = TRUE),

  # ---- WASH access
  ihme_wash = grep("_w_", metadata$ihme_vars, value = TRUE),

  # ---- IHME other
  ihme_other = setdiff(
    metadata$ihme_vars,
    grep("(hiv|death|mortality|diarrhea|malaria|stunting|wasting|underweight|anemia|education|_w_)",
         metadata$ihme_vars, value = TRUE)
  ),

  # ============================================================
  # 5. Malaria Atlas Project (MAP)
  # ============================================================



  # ---- Malaria burden
  map_burden = grep("Incidence|Mortality|Parasite", metadata$map_vars, value = TRUE),

  # ---- Vector control & treatment
  map_interventions = grep("Net|IRS|Treatment", metadata$map_vars, value = TRUE),

  # ---- Transmission intensity
  map_transmission = grep("Reproductive", metadata$map_vars, value = TRUE),

  # ============================================================
  # 6. WFP food availability
  # ============================================================


  # ---- Staples
  wfp_staples = grep("rice|maize|millet|sorghum|cassava",
                     metadata$wfp_vars, value = TRUE),

  # ---- Animal-source foods
  wfp_animal_source = grep("meat|fish|eggs",
                           metadata$wfp_vars, value = TRUE),

  # ---- Legumes & oilseeds
  wfp_legumes_oilseeds = grep("soy|cowpea",
                              metadata$wfp_vars, value = TRUE),

  # ---- Fruits & vegetables
  wfp_fruits_vegetables = grep("tomatoes|onions|peppers|eggplants",
                               metadata$wfp_vars, value = TRUE),

  # ---- WFP other
  wfp_other = setdiff(
    metadata$wfp_vars,
    grep("(rice|maize|millet|sorghum|cassava|meat|fish|eggs|soy|cowpea|tomatoes|onions|peppers|eggplants)",
         metadata$wfp_vars, value = TRUE)
  ),

  # ============================================================
  # 7. FluNet
  # ============================================================

  flunet = list(
    flunet_influenza = metadata$flunet_vars
  ),

  # ============================================================
  # 8. Google Earth Engine (GEE)
  # ============================================================


  # ---- Rainfall
  gee_rainfall = grep("trmm", metadata$gee_vars, value = TRUE),

  # ---- Vegetation & greenness
  gee_vegetation = grep("ndvi", metadata$gee_vars, value = TRUE),

  # ---- Temperature
  gee_temperature = grep("temp", metadata$gee_vars, value = TRUE),

  # ---- Primary productivity
  gee_productivity = grep("npp|tpp|fpp", metadata$gee_vars, value = TRUE),

  # ---- Built environment & urbanicity
  gee_built_env = grep("built|smod", metadata$gee_vars, value = TRUE),

  # ---- Land cover
  gee_landcover = grep("grassland|cropland|ghm", metadata$gee_vars, value = TRUE)

))


prune_domains <- function(
    domains,
    X,
    drop_empty = TRUE,
    return_removed = FALSE
) {

  removed <- list()

  prune_one <- function(obj, path = character()) {

    if (is.list(obj)) {
      out <- list()

      for (nm in names(obj)) {
        res <- prune_one(obj[[nm]], c(path, nm))
        if (!is.null(res) && (!drop_empty || length(res) > 0)) {
          out[[nm]] <- res
        }
      }

      return(out)
    }

    if (is.character(obj)) {
      keep <- intersect(obj, X)
      drop <- setdiff(obj, X)

      if (return_removed && length(drop) > 0) {
        removed[[paste(path, collapse = "::")]] <<- drop
      }

      return(keep)
    }

    stop("All leaves must be character vectors.")
  }

  pruned <- prune_one(domains)

  if (return_removed) {
    return(list(
      pruned = pruned,
      removed = removed
    ))
  }

  pruned
}

variable_domains <- prune_domains(
  domains = variable_domains_nogw,
  X = X
)



sl_fit=res$sl_fit
eval_fun = loss_squared_error
importance.metric="difference"
covariate_groups=covariate_groups
importance.metric="difference"
n_vars=10
covariate_groups=covariate_groups


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


res_vim_domains=try(calc_mn_importance(res$sl_fit,  eval_fun = loss_squared_error, importance.metric="difference",
                                       covariate_groups=variable_domains))
res_vim_domains
res_vim_domains$p

saveRDS(res_vim_domains, here("results/vim/res_GW_Ghana_SL_mom_iron_vim_variable_domains.rds"))





res_vim_diff_data_sources=try(calc_mn_importance(res$sl_fit,  eval_fun = loss_squared_error, importance.metric="difference",
                                                 covariate_groups=covariate_groups))
res_vim_diff_data_sources
res_vim_diff_data_sources$p

saveRDS(res_vim_diff_data_sources, here("results/vim/res_GW_Ghana_SL_mom_iron_vim_data_sources.rds"))


res_vim=try(calc_mn_importance(res$sl_fit,  eval_fun = loss_squared_error, importance.metric="difference",
                               covariate_groups=NULL, n_vars=20))
res_vim
res_vim$p
