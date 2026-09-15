




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
library(cowplot)


source(paste0(here::here(),"/src/0-functions.R"))
source(paste0(here::here(),"/src/0-SL-setup.R"))
source(paste0(here::here(),"/src/DHS/DHS_functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_variable_recode.R"))


res_vim = readRDS(here("results/vim/res_GW_Ghana_SL_mom_b12_vim_individual.rds"))
res_vim_diff_plot_df <- res_vim$varimp
res_vim_diff_plot_df <- res_vim_diff_plot_df %>% filter(covariate!="Dataid")
res_vim_diff_plot_df$covariate = rename_vars_for_plotting(res_vim_diff_plot_df$covariate, extract_labels(c(ir_rename_spec, pr_rename_spec, hr_rename_spec)) )
res_vim_diff_plot_df$covariate = gsub("Global ", "", res_vim_diff_plot_df$covariate)
res_vim_diff_plot_df$covariate = gsub("Africa ", "", res_vim_diff_plot_df$covariate)
res_vim_diff_plot_df$covariate = gsub("Pf", "Malaria", res_vim_diff_plot_df$covariate)
p=ggplot(res_vim_diff_plot_df %>%
            arrange(MSE_difference) %>% tail(n=20), aes(x = reorder(covariate, -MSE_difference), y = MSE_difference)) +
  geom_point() +
  coord_flip() +
  labs(title = "Top 20 Variables by effect on MSE",
       subtitle = "Women's vitamin B12 concentration in Ghana",
       x = "Covariate",
       y = "MSE_difference") +
  theme_minimal()
p




res_vim_domains = readRDS(here("results/vim/res_GW_Ghana_SL_mom_b12_vim_variable_domains.rds"))
res_vim_diff_data_sources = readRDS(here("results/vim/res_GW_Ghana_SL_mom_b12_vim_data_sources.rds"))
#res_vim_diff_individual = readRDS(here("results/vim/res_GW_Ghana_SL_mom_b12_vim_individual.rds"))

res_vim_diff_plot_df <- res_vim_domains$varimp
res_vim_diff_plot_df$covariate = rename_vars_for_plotting(res_vim_diff_plot_df$covariate, extract_labels(c(ir_rename_spec, pr_rename_spec, hr_rename_spec)) )
res_vim_diff_plot_df$covariate = gsub("Global ", "", res_vim_diff_plot_df$covariate)
res_vim_diff_plot_df$covariate = gsub("Africa ", "", res_vim_diff_plot_df$covariate)
res_vim_diff_plot_df$covariate = gsub("Pf", "Malaria", res_vim_diff_plot_df$covariate)
p1=ggplot(res_vim_diff_plot_df %>%
            arrange(MSE_difference) %>% tail(n=20), aes(x = reorder(covariate, -MSE_difference), y = MSE_difference)) +
  geom_point() +
  coord_flip() +
  labs(title = "Top 20 Variable Domains by effect on MSE",
       subtitle = "Women's vitamin B12 concentration in Ghana",
       x = "Covariate",
       y = "MSE_difference") +
  theme_minimal()
p1



res_vim_diff_plot_df <- res_vim_diff_data_sources$varimp
res_vim_diff_plot_df$covariate = rename_vars_for_plotting(res_vim_diff_plot_df$covariate, extract_labels(c(ir_rename_spec, pr_rename_spec, hr_rename_spec)) )
res_vim_diff_plot_df$covariate = gsub("Global ", "", res_vim_diff_plot_df$covariate)
res_vim_diff_plot_df$covariate = gsub("Africa ", "", res_vim_diff_plot_df$covariate)
res_vim_diff_plot_df$covariate = gsub("Pf", "Malaria", res_vim_diff_plot_df$covariate)
p2=ggplot(res_vim_diff_plot_df %>%
            arrange(MSE_difference) %>% tail(n=20), aes(x = reorder(covariate, -MSE_difference), y = MSE_difference)) +
  geom_point() +
  coord_flip() +
  labs(title = "Data sources by effect on MSE",
       subtitle = "Women's vitamin B12 concentration in Ghana",
       x = "Covariate",
       y = "MSE_difference") +
  theme_minimal()
p2


res_vim_domains = readRDS(here("results/vim/res_GW_Ghana_SL_mom_iron_vim_data_sources.rds"))
res_vim_diff_data_sources = readRDS(here("results/vim/res_GW_Ghana_SL_mom_iron_vim_variable_domains.rds"))

res_vim_diff_plot_df <- res_vim_domains$varimp
res_vim_diff_plot_df$covariate = rename_vars_for_plotting(res_vim_diff_plot_df$covariate, extract_labels(c(ir_rename_spec, pr_rename_spec, hr_rename_spec)) )
res_vim_diff_plot_df$covariate = gsub("Global ", "", res_vim_diff_plot_df$covariate)
res_vim_diff_plot_df$covariate = gsub("Africa ", "", res_vim_diff_plot_df$covariate)
res_vim_diff_plot_df$covariate = gsub("Pf", "Malaria", res_vim_diff_plot_df$covariate)
p3=ggplot(res_vim_diff_plot_df %>%
            arrange(MSE_difference) %>% tail(n=20), aes(x = reorder(covariate, -MSE_difference), y = MSE_difference)) +
  geom_point() +
  coord_flip() +
  labs(title = "",
       subtitle = "Women's iron concentration in Ghana",
       x = "",
       y = "MSE_difference") +
  theme_minimal()
p3



res_vim_diff_plot_df <- res_vim_diff_data_sources$varimp
res_vim_diff_plot_df$covariate = rename_vars_for_plotting(res_vim_diff_plot_df$covariate, extract_labels(c(ir_rename_spec, pr_rename_spec, hr_rename_spec)) )
res_vim_diff_plot_df$covariate = gsub("Global ", "", res_vim_diff_plot_df$covariate)
res_vim_diff_plot_df$covariate = gsub("Africa ", "", res_vim_diff_plot_df$covariate)
res_vim_diff_plot_df$covariate = gsub("Pf", "Malaria", res_vim_diff_plot_df$covariate)
p4=ggplot(res_vim_diff_plot_df %>%
            arrange(MSE_difference) %>% tail(n=20), aes(x = reorder(covariate, -MSE_difference), y = MSE_difference)) +
  geom_point() +
  coord_flip() +
  labs(title = "",
       subtitle = "Women's iron concentration in Ghana",
       x = "",
       y = "MSE_difference") +
  theme_minimal()
p4

cowplot::plot_grid(p1, p4)
cowplot::plot_grid(p2, p3)


#
# res_vim_domains = readRDS(here("results/vim/res_GW_Malawi_SL_mom_iron_vim_data_sources.rds"))
 res_vim_diff_data_sources = readRDS(here("results/vim/res_GW_Malawi_SL_mom_iron_vim_variable_domains.rds"))
#
# res_vim_diff_plot_df <- res_vim_domains$varimp
# res_vim_diff_plot_df$covariate = rename_vars_for_plotting(res_vim_diff_plot_df$covariate, extract_labels(c(ir_rename_spec, pr_rename_spec, hr_rename_spec)) )
# res_vim_diff_plot_df$covariate = gsub("Global ", "", res_vim_diff_plot_df$covariate)
# res_vim_diff_plot_df$covariate = gsub("Africa ", "", res_vim_diff_plot_df$covariate)
# res_vim_diff_plot_df$covariate = gsub("Pf", "Malaria", res_vim_diff_plot_df$covariate)
# p=ggplot(res_vim_diff_plot_df %>%
#            arrange(MSE_difference) %>% tail(n=20), aes(x = reorder(covariate, -MSE_difference), y = MSE_difference)) +
#   geom_point() +
#   coord_flip() +
#   labs(title = "Top 20 Variable Domains by effect on MSE",
#        subtitle = "Women's iron concentration in Malawi",
#        x = "Covariate",
#        y = "MSE_difference") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p
#
#
#
# res_vim_diff_plot_df <- res_vim_diff_data_sources$varimp
# res_vim_diff_plot_df$covariate = rename_vars_for_plotting(res_vim_diff_plot_df$covariate, extract_labels(c(ir_rename_spec, pr_rename_spec, hr_rename_spec)) )
# res_vim_diff_plot_df$covariate = gsub("Global ", "", res_vim_diff_plot_df$covariate)
# res_vim_diff_plot_df$covariate = gsub("Africa ", "", res_vim_diff_plot_df$covariate)
# res_vim_diff_plot_df$covariate = gsub("Pf", "Malaria", res_vim_diff_plot_df$covariate)
# p=ggplot(res_vim_diff_plot_df %>%
#            arrange(MSE_difference) %>% tail(n=20), aes(x = reorder(covariate, -MSE_difference), y = MSE_difference)) +
#   geom_point() +
#   coord_flip() +
#   labs(title = "Data sources by effect on MSE",
#        subtitle = "Women's iron concentration in Malawi",
#        x = "Covariate",
#        y = "MSE_difference") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p

 res_vim_diff_data_sources = readRDS(here("results/vim/res_GW_Ghana_SL_mom_b12_vim_variable_domains.rds"))
 res_vim_diff_data_sources$varimp
 res_vim_diff_data_sources = readRDS(here("results/vim/res_GW_Ghana_SL_mom_iron_vim_variable_domains.rds"))
 res_vim_diff_data_sources$varimp
 res_vim_diff_data_sources = readRDS(here("results/vim/res_GW_Malawi_SL_mom_iron_vim_variable_domains.rds"))
 res_vim_diff_data_sources$varimp
