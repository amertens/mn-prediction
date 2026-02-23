
# --- setup ----
rm(list = ls())
library(tidymodels)
library(tidyverse)
library(haven)
library(here)
library(janitor)
library(skimr)
library(vip)

set.seed(123)

# --- load data (your paths from the prompt) ----
df_women    <- read_dta(here("data","IPD","Ghana","GMS_WomenAllData_DavisBerkeley.dta")) #%>% clean_names()
df_children <- read_dta(here("data","IPD","Ghana","GMS_ChildrenAllData_DavisBerkeley.dta")) #%>% clean_names()

# --- handy helpers ----
fix_missing_sentinels <- function(df) {
  # drop the common '999' sentinel in doubles (seen for Hb etc)
  df %>% mutate(across(where(is.double), ~na_if(.x, 999)))
}

has_col <- function(df, nm) nm %in% names(df)

pick_existing <- function(df, candidates) {
  # keep only the names that actually exist in df
  candidates[candidates %in% names(df)]
}

logify <- function(df, vars) {
  for (v in vars) {
    if (v %in% names(df)) {
      df[[paste0("log1p_", v)]] <- scale(log(df[[v]]))
    }
  }
  df
}

# # --- clean ----
# df_women    <- fix_missing_sentinels(df_women)
# df_children <- fix_missing_sentinels(df_children)

# --- define outcomes present in each file (names from codebooks) ----
outcomes_women <- c(
  hb      = "wm_whbc",
  ferritin= "fer",
  stfr    = "stfr",
  bis     = "bis",
  rbp     = "rbp",
  crp     = "crp",
  agp     = "agp",
  b12     = "wb12",
  folate  = "wfolate"
)

outcomes_child <- c(
  hb      = "gchb",
  ferritin= "fer",
  stfr    = "stfr",
  bis     = "bis",
  rbp     = "rbp",
  crp     = "crp",
  agp     = "agp"
)

summary(df_women$wFolate)


# Women: use these columns
outcomes_women <- c(
  hb_for_anemia = "wm_whbc",        # apply WHO altitude adjustment at classification time
  ferritin_adj  = "wFerrAdjThurn",  # adjusted ferritin available; use this for ID prevalence
  stfr_raw      = "wSTFR",          # used as-is
  rbp_raw       = "wRBP",           # women typically unadjusted for VAD; if you prefer adjusted: "wRBPAdjThurn"
  bis_raw       = "bis",            # recompute BIS from adjusted ferritin if you need BIS_adj
  b12_raw       = "wB12",
  folate_raw    = "wFolate"
)


summary(df_children$bis)

# Children: use these columns
outcomes_child <- c(
  hb_for_anemia = "gchb",              # apply WHO altitude adjustment at classification time
  ferritin_adj  = "cFerrAdjThurn",     # adjusted ferritin available; use for ID prevalence
  stfr_raw      = "cSTFR",             # used as-is
  rbp_adj       = "cRBPAdjBrinda",     # adjusted RBP available; use for VAD prevalence
  bis_raw       = "bis"                # recompute BIS from adjusted ferritin if you need BIS_adj
)

# restrict to outcomes that actually exist in your files
outcomes_women  <- pick_existing(df_women,  outcomes_women)
outcomes_child  <- pick_existing(df_children, outcomes_child)


for(i in outcomes_women){
  print(summary(df_women[[i]]))
}

for(i in outcomes_child){
  print(summary(df_children[[i]]))
}


# create  scaled and log10 versions for skewed outcomes
df_women    <- logify(df_women,    names(outcomes_women))
df_children <- logify(df_children, names(outcomes_child))

for(i in outcomes_child){
  print(summary(df_children[[i]]))
}

for(i in outcomes_women){
  print(summary(df_women[[i]]))
}

# --- choose predictor sets (flexible fallbacks) ----
# preds_women_cands <- c(
#   "age_years", "wbmi",                         # demography/anthro
#   "region", "stratum", "belt", "urban", "residence", # geography
#   "wealth_quintile", "wealthq", "hwealthquintile",   # SES
#   "wmalariardtyn", "wsickleyn", "walphathalyn",      # infection/hemoglobinopathies
#   "crp", "agp"                                     # inflammation
# )
# preds_child_cands <- c(
#   "age_months","agemonths","agem","sex","child_sex",
#   "region","stratum","belt","urban","residence",
#   "wealth_quintile","wealthq","hwealthquintile",
#   "cmalariardtyn","cmalariayn","csickleyn","calphathalyn",
#   "crp","agp"
# )

# ---- WOMEN: survey-only predictors including diet variables ----
preds_women_cands <- c(
  # Demography & anthropometry
  "wAgeYears","wAgeYearsCat","wBMI","wBMICat","wWeightKg","wHeightCm","wHeightM",
  # Geography / SES / WASH
  "region","hWealthquintile","hDrinkWaterSafeYN","hSanitation","hHandwashingH2OYN",
  # Infection / hemoglobinopathies
  "wMalariaRDTYN","wMalariaRDTCat","wMalariaYN",
  "sickle_genotype","thal_genotype","wSickle","wSickleYN","wSickleAS",
  "wAlphaThal","wAlphaThalYN",
  # Diet & supplementation
  "wFoodGrp5plus","wSuppIron","wSuppFolic","wSuppMulti","wSupplastPregFeFA","wSuppVitA",
  # Education & literacy / marital
  "wEducCat","wRead","wMarital",
  # Inflammation
  "crp","agp",
  # Survey weight
  "sWeight"
)

# ---- CHILDREN: survey-only predictors including diet variables ----
preds_child_cands <- c(
  # Demography & anthropometry
  "age_days","age_months","cAgeMonths","cAgeMonthsCat","cSex","cSex01",
  "cBirthweight","cBirthweightCat",
  "haz06","waz06","whz06",
  "cStuntCat","cStuntYN","cStuntSev",
  "cWasteCat","cWasteYN","cWasteSev","cWasteOWCat","cOverweightYN",
  # Geography / SES / WASH
  "region","hWealthquintile","hDrinkWaterSafeYN","hSanitation","hHandwashingH2OYN",
  # Infection / hemoglobinopathies
  "cMalariaRDTYN","cMalariaRDTCat","cMalariaYN",
  "sickle_genotype","thal_genotype","cSickle","cSickleYN","cSickle01",
  "cAlphaThal","cAlphaThalYN","cAlphaThal01",
  # Morbidity
  "cFever","cCough","cDiarrhea","cDiarrheaBlood","cFeverMalaria","cFeverMalariaResult",
  # Diet & supplementation
  "cFoodGrp5plus","cVitASupp","cSuppIron","cSuppFolic","cSuppMulti",
  # Caregiver education
  "wEducCat",
  # Survey weight
  "sWeight"
)


preds_women <- pick_existing(df_women, preds_women_cands)
preds_child <- pick_existing(df_children, preds_child_cands)

# optional: include survey weights if present
w_col_w <- if (has_col(df_women, "sweight")) "sweight" else NULL
w_col_c <- if (has_col(df_children, "sweight")) "sweight" else NULL

# --- quick EDA summaries and plots ----
summarize_outcomes <- function(df, outcomes_named) {
  cols <- unname(outcomes_named[outcomes_named %in% names(df)])
  df %>%
    select(all_of(cols)) %>%
    skim()
}

plot_outcomes <- function(df, outcomes_named, title = "Outcome distributions") {
  cols <- unname(outcomes_named[outcomes_named %in% names(df)])
  df %>%
    select(all_of(cols)) %>%
    pivot_longer(everything(), names_to = "outcome", values_to = "value") %>%
    ggplot(aes(value)) +
    geom_histogram(bins = 40) +
    facet_wrap(~ outcome, scales = "free") +
    labs(title = title, x = NULL, y = "Count")
}

fit_benchmarks <- function(df, outcome, predictors, weight_col = NULL) {

  # 1) subset & basic checks
  model_df <- df %>%
    select(all_of(c(outcome, predictors, weight_col))) %>%
    drop_na(all_of(outcome))
  if (nrow(model_df) < 50) return(NULL)

  # 2) recipe
  rec <- recipe(as.formula(paste(outcome, "~ .")), data = model_df) %>%
    update_role(all_of(predictors), new_role = "predictor") %>%
    update_role(all_of(outcome),    new_role = "outcome") %>%
    step_rm(any_of(weight_col)) %>%                 # keep weights out of X
    step_zv(all_predictors()) %>%
    step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_normalize(all_numeric_predictors())

  folds <- vfold_cv(model_df, v = 5)

  # -------- Elastic net --------
  en_spec <- linear_reg(penalty = tune(), mixture = tune()) %>% set_engine("glmnet")
  en_wf   <- workflow() %>% add_model(en_spec) %>% add_recipe(rec)

  en_res <- tune_grid(
    en_wf,
    resamples = folds,
    grid = grid_space_filling(penalty(), mixture(), size = 20),
    metrics = metric_set(rmse, rsq)
  )
  best_en  <- select_best(en_res, metric = "rmse")
  en_final <- finalize_workflow(en_wf, best_en) %>% fit(model_df)

  en_rmse <- collect_metrics(en_res) %>%
    filter(.metric == "rmse",
           penalty == best_en$penalty,
           mixture == best_en$mixture) %>%
    pull(mean)
  en_mse  <- en_rmse^2

  en_imp <- broom::tidy(extract_fit_parsnip(en_final)) %>%
    filter(term != "(Intercept)") %>%
    transmute(term, importance = abs(estimate)) %>%
    arrange(desc(importance))

  # -------- Random forest --------
  rf_spec <- rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>%
    set_engine(
      "ranger",
      importance   = "permutation",
      case.weights = if (!is.null(weight_col) && weight_col %in% names(model_df))
        model_df[[weight_col]] else NULL
    ) %>%
    set_mode("regression")

  rf_wf <- workflow() %>% add_model(rf_spec) %>% add_recipe(rec)

  rf_res <- tune_grid(
    rf_wf,
    resamples = folds,
    grid = grid_space_filling(
      mtry(range = c(1L, max(1L, length(predictors) - 1L))),
      min_n(),
      size = 20
    ),
    metrics = metric_set(rmse, rsq)
  )
  best_rf  <- select_best(rf_res, metric = "rmse")
  rf_final <- finalize_workflow(rf_wf, best_rf) %>% fit(model_df)

  rf_rmse <- collect_metrics(rf_res) %>%
    filter(.metric == "rmse",
           mtry == best_rf$mtry,
           min_n == best_rf$min_n) %>%
    pull(mean)
  rf_mse  <- rf_rmse^2

  rf_imp <- vip::vi(extract_fit_engine(rf_final)) %>%
    arrange(desc(Importance)) %>%
    transmute(term = Variable, importance = Importance)

  # 3) return results
  tibble(
    outcome    = outcome,
    model      = c("elastic_net", "random_forest"),
    mse        = c(en_mse, rf_mse),
    importance = list(en_imp, rf_imp)
  )
}

run_benchmarks_for <- function(df, outcomes, predictors, weight_col = NULL) {
  keep_outcomes <- unname(outcomes[outcomes %in% names(df)])
  map_dfr(keep_outcomes, ~fit_benchmarks(df, .x, predictors, weight_col))
}

# --- RUN: EDA ----
cat("\n==== WOMEN: outcome summary ====\n")
print(summarize_outcomes(df_women, outcomes_women))
print(plot_outcomes(df_women, outcomes_women, "Women's outcomes"))

cat("\n==== CHILDREN: outcome summary ====\n")
print(summarize_outcomes(df_children, outcomes_child))
print(plot_outcomes(df_children, outcomes_child, "Children's outcomes"))

# --- RUN: Benchmarks ----
bench_women   <- run_benchmarks_for(df_women,   outcomes_women, preds_women, w_col_w)
bench_children<- run_benchmarks_for(df_children, outcomes_child, preds_child, w_col_c)

for(i in outcomes_women){
  cat(i,"\n")
  print(sd(df_women[[i]], na.rm=T))
}

# Show MSE table per outcome (both models)
cat("\n==== WOMEN: cross-validated MSE ====\n")
bench_women %>% mutate(rmse=sqrt(mse)) %>% select(outcome, model, rmse) %>% arrange(outcome, rmse) %>% print(n = Inf)

cat("\n==== CHILDREN: cross-validated MSE ====\n")
bench_children %>% mutate(rmse=sqrt(mse)) %>% select(outcome, model, mse) r%>% arrange(outcome, rmse) %>% print(n = Inf)

# Example: top 10 predictors for the best model per outcome
top_imp <- function(bench_df) {
  bench_df %>%
    group_by(outcome) %>%
    slice_min(order_by = mse, n = 1, with_ties = FALSE) %>%
    transmute(outcome, model, importance) %>%
    rowwise() %>%
    mutate(top10 = list(importance %>% slice_max(importance, n = 10))) %>%
    ungroup()
}

cat("\n==== WOMEN: top 10 predictors for best model per outcome ====\n")
print(top_imp(bench_women))

cat("\n==== CHILDREN: top 10 predictors for best model per outcome ====\n")
print(top_imp(bench_children))








library(tidyverse)
library(glue)

# Normalize importance data frames coming from different models
.norm_imp <- function(x) {
  nm <- names(x)
  if (all(c("Variable", "Importance") %in% nm)) {
    x %>% transmute(term = Variable, importance = as.numeric(Importance))
  } else if (all(c("term", "importance") %in% nm)) {
    x %>% transmute(term, importance = as.numeric(importance))
  } else {
    # Fallback: use first 2 columns
    x %>% setNames(c("term","importance", names(x)[-(1:2)])) %>%
      transmute(term, importance = as.numeric(importance))
  }
}

# Build a clean long table of top-N importances for each outcome/model
top_importance_table <- function(bench_df, n = 10) {
  bench_df %>%
    mutate(importance = purrr::map(importance, .norm_imp)) %>%
    mutate(top = purrr::map(importance, ~ .x %>% slice_max(importance, n = n))) %>%
    select(outcome, model, top) %>%
    unnest(top) %>%
    arrange(outcome, model, desc(importance)) %>%
    mutate(importance = round(importance, 4))
}

# Pretty print to console (grouped by outcome)
print_importance <- function(bench_df, n = 10) {
  tt <- top_importance_table(bench_df, n)
  purrr::walk(unique(tt$outcome), function(yy) {
    cat(glue("\n==== Top {n} predictors — {yy} ====\n"))
    tt %>%
      filter(outcome == yy) %>%
      select(model, term, importance) %>%
      arrange(model, desc(importance)) %>%
      print(n = Inf)
  })
}

# OPTIONAL: quick ggplot for each outcome (saved to files)
plot_importance <- function(bench_df, n = 10, outdir = "outputs/importance_plots") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  tt <- top_importance_table(bench_df, n)
  purrr::walk(unique(tt$outcome), function(yy) {
    p <- tt %>%
      filter(outcome == yy) %>%
      mutate(term = forcats::fct_reorder(term, importance)) %>%
      ggplot(aes(term, importance, fill = model)) +
      geom_col(position = position_dodge(width = 0.8)) +
      coord_flip() +
      labs(title = glue("Top {n} predictors — {yy}"), x = NULL, y = "Importance") +
      theme_minimal(base_size = 12)
    ggsave(file.path(outdir, glue("{yy}_top{n}_importance.png")), p, width = 7, height = 5, dpi = 300)
  })
}


# Console-friendly tables
print_importance(bench_women, n = 10)
print_importance(bench_children, n = 10)

# If you want a single combined data frame for export:
imp_women    <- top_importance_table(bench_women, n = 10)
imp_children <- top_importance_table(bench_children, n = 10)
# readr::write_csv(imp_women,    "outputs/top10_importance_women.csv")
# readr::write_csv(imp_children, "outputs/top10_importance_children.csv")

# Optional plots saved as PNGs
plot_importance(bench_women, n = 10)
plot_importance(bench_children, n = 10)





# ---- slide-friendly tables with gt ----
library(gt)
library(scales)

slide_friendly_tables <- function(bench_df, n = 10, subtitle = NULL) {
  tt <- top_importance_table(bench_df, n) %>%
    group_by(outcome, model) %>%
    arrange(desc(importance), .by_group = TRUE) %>%
    ungroup() %>%
    mutate(rank = ave(-importance, outcome, model, FUN = rank),
           rank = as.integer(rank)) %>%
    relocate(rank, .after = model)

  # Split into a gt table per outcome for easy per‑slide copy/paste
  outcomes <- unique(tt$outcome)

  tabs <- purrr::map(outcomes, function(yy) {
    dat <- tt %>%
      filter(outcome == yy) %>%
      select(Model = model, Rank = rank, Predictor = term, Importance = importance)

    gt(dat) |>
      tab_header(
        title = md(glue::glue("**Top {n} predictors — {yy}**")),
        subtitle = subtitle
      ) |>
      fmt_number(columns = "Importance", decimals = 3) |>
      data_color(
        columns = "Importance",
        colors = col_numeric(c("#e8f5e9", "#fff3e0", "#ffebee"),  # light green → light orange → light red
                             domain = range(dat$Importance, na.rm = TRUE))
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels(everything())
      ) |>
      opt_table_font(
        font = list(google_font("Source Sans Pro"), default_fonts())
      ) |>
      tab_options(
        table.font.size = px(16),
        data_row.padding = px(6),
        table.width = px(900),
        column_labels.background.color = "#F7F7F7"
      ) |>
      opt_row_striping()
  })

  # Return a named list of gt tables so you can print or gtsave() as needed
  names(tabs) <- outcomes
  tabs
}

# Example usage:
gt_women_tabs    <- slide_friendly_tables(bench_women, n = 10)
gt_children_tabs <- slide_friendly_tables(bench_children, n = 10)
gt_women_tabs$rbp          # prints the RBP table for women in the Viewer
#gtsave(gt_women_tabs$rbp, "rbp_top10_women.png", vwidth = 1200, vheight = 700, expand = 2)


# ---- output plots instead of saving ----
library(ggplot2)
library(forcats)

plot_importance <- function(bench_df, n = 10) {
  tt <- top_importance_table(bench_df, n)

  plots <- purrr::map(unique(tt$outcome), function(yy) {
    p <- tt %>%
      filter(outcome == yy) %>%
      mutate(term = fct_reorder(term, importance)) %>%
      ggplot(aes(term, importance, fill = model)) +
      geom_col(position = position_dodge(width = 0.85), width = 0.8) +
      coord_flip() +
      labs(
        title = glue::glue("Top {n} predictors — {yy}"),
        x = NULL, y = "Importance"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "top",
        plot.title = element_text(face = "bold")
      )
    print(p)  # output to the graphics device
    p
  })

  names(plots) <- unique(tt$outcome)
  invisible(plots)  # return plots (named list) without spamming the console
}

# Example usage:
p_women    <- plot_importance(bench_women, n = 10)       # prints each plot, returns list
p_children <- plot_importance(bench_children, n = 10)
p_women$rbp                                            # reprint a specific plot
p_children$rbp

