


rm(list=ls())
library(tidyverse)
library(here)
library(pROC)

source(paste0(here::here(),"/src/0-functions.R"))


res1 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_child_vitA_full.rds"))
res2 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_child_vitA_gwPreds.rds"))
res3 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_child_vitA.rds"))
res4 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_women_vitA_full.rds"))
res5 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_women_vitA_gwPreds.rds"))
res6 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_women_vitA.rds"))
res7 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_women_b12_full.rds"))
res8 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_women_b12_gwPreds.rds"))
res9 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_women_b12.rds"))
res10 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_women_folate_full.rds"))
res11 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_women_folate_gwPreds.rds"))
res12 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_women_folate.rds"))
res13 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_child_iron_full.rds"))
res14 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_child_iron_gwPreds.rds"))
res15 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_child_iron.rds"))
res16 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_mom_iron_full.rds"))
res17 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_mom_iron_gwPreds.rds"))
res18 <- readRDS(here("results/models/res_bin_GW_Ghana_SL_mom_iron.rds"))






res=res1
pred_source = c("auto")
weights = NULL
ci = TRUE
ci_method = c("delong")
quiet = F


auc_pr_brier_res1 = auc_pr_brier_from_sl3(res1, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-A deficiency", Xvars="All", population="children")
auc_pr_brier_res2 = auc_pr_brier_from_sl3(res2, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-A deficiency", Xvars="Survey", population="children")
auc_pr_brier_res3 = auc_pr_brier_from_sl3(res3, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-A deficiency", Xvars="Proxy", population="children")
auc_pr_brier_res4 = auc_pr_brier_from_sl3(res4, pred_source = c("auto"),ci = TRUE,quiet = TRUE)  %>% mutate(Y="Vit-A deficiency", Xvars="All", population="women")
auc_pr_brier_res5 = auc_pr_brier_from_sl3(res5, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-A deficiency", Xvars="Survey", population="women")
auc_pr_brier_res6 = auc_pr_brier_from_sl3(res6, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-A deficiency", Xvars="Proxy", population="women")
auc_pr_brier_res7 = auc_pr_brier_from_sl3(res7, pred_source = c("auto"),ci = TRUE,quiet = TRUE)  %>% mutate(Y="Vit-B12 deficiency", Xvars="All", population="women")
auc_pr_brier_res8 = auc_pr_brier_from_sl3(res8, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-B12 deficiency", Xvars="Survey", population="women")
auc_pr_brier_res9 = auc_pr_brier_from_sl3(res9, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-B12 deficiency", Xvars="Proxy", population="women")
auc_pr_brier_res10 = auc_pr_brier_from_sl3(res10, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Folate deficiency", Xvars="All", population="women")
auc_pr_brier_res11 = auc_pr_brier_from_sl3(res11, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Folate deficiency", Xvars="Survey", population="women")
auc_pr_brier_res12 = auc_pr_brier_from_sl3(res12, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Folate deficiency", Xvars="Proxy", population="women")
auc_pr_brier_res13 = auc_pr_brier_from_sl3(res13, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Ferratin deficiency", Xvars="All", population="children")
auc_pr_brier_res14 = auc_pr_brier_from_sl3(res14, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Ferratin deficiency", Xvars="Survey", population="children")
auc_pr_brier_res15 = auc_pr_brier_from_sl3(res15, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Ferratin deficiency", Xvars="Proxy", population="children")
auc_pr_brier_res16 = auc_pr_brier_from_sl3(res16, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Ferratin deficiency", Xvars="All", population="women")
auc_pr_brier_res17 = auc_pr_brier_from_sl3(res17, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Ferratin deficiency", Xvars="Survey", population="women")
auc_pr_brier_res18 = auc_pr_brier_from_sl3(res18, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Ferratin deficiency", Xvars="Proxy", population="women")

res <- bind_rows(auc_pr_brier_res1, auc_pr_brier_res2, auc_pr_brier_res3, auc_pr_brier_res4, auc_pr_brier_res5, auc_pr_brier_res6,
                 auc_pr_brier_res7, auc_pr_brier_res8, auc_pr_brier_res9, auc_pr_brier_res10,
                 auc_pr_brier_res11, auc_pr_brier_res12, auc_pr_brier_res13, auc_pr_brier_res14,
                 auc_pr_brier_res15, auc_pr_brier_res16, auc_pr_brier_res17, auc_pr_brier_res18)
head(res)

saveRDS(res, file = here("results/res_bin_GW_Ghana_SL_all.rds"))









res1 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_vitA_full.rds"))
res2 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_vitA_gwPreds.rds"))
res3 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_vitA.rds"))
res4 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_vitA_full.rds"))
res5 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_vitA_gwPreds.rds"))
res6 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_vitA.rds"))
res7 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_b12_full.rds"))
res8 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_b12_gwPreds.rds"))
res9 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_b12.rds"))
res10 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_folate_full.rds"))
res11 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_folate_gwPreds.rds"))
res12 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_folate.rds"))
res13 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_iron_full.rds"))
res14 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_iron_gwPreds.rds"))
res15 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_child_iron.rds"))
res16 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_mom_iron_full.rds"))
res17 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_mom_iron_gwPreds.rds"))
res18 <- readRDS(here("results/models/res_full_bin_GW_Ghana_SL_mom_iron.rds"))


auc_pr_brier_res1 = auc_pr_brier_from_sl3(res1, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-A deficiency", Xvars="All", population="children")
auc_pr_brier_res2 = auc_pr_brier_from_sl3(res2, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-A deficiency", Xvars="Survey", population="children")
auc_pr_brier_res3 = auc_pr_brier_from_sl3(res3, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-A deficiency", Xvars="Proxy", population="children")
auc_pr_brier_res4 = auc_pr_brier_from_sl3(res4, pred_source = c("auto"),ci = TRUE,quiet = TRUE)  %>% mutate(Y="Vit-A deficiency", Xvars="All", population="women")
auc_pr_brier_res5 = auc_pr_brier_from_sl3(res5, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-A deficiency", Xvars="Survey", population="women")
auc_pr_brier_res6 = auc_pr_brier_from_sl3(res6, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-A deficiency", Xvars="Proxy", population="women")
auc_pr_brier_res7 = auc_pr_brier_from_sl3(res7, pred_source = c("auto"),ci = TRUE,quiet = TRUE)  %>% mutate(Y="Vit-B12 deficiency", Xvars="All", population="women")
auc_pr_brier_res8 = auc_pr_brier_from_sl3(res8, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-B12 deficiency", Xvars="Survey", population="women")
auc_pr_brier_res9 = auc_pr_brier_from_sl3(res9, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Vit-B12 deficiency", Xvars="Proxy", population="women")
auc_pr_brier_res10 = auc_pr_brier_from_sl3(res10, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Folate deficiency", Xvars="All", population="women")
auc_pr_brier_res11 = auc_pr_brier_from_sl3(res11, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Folate deficiency", Xvars="Survey", population="women")
auc_pr_brier_res12 = auc_pr_brier_from_sl3(res12, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Folate deficiency", Xvars="Proxy", population="women")
auc_pr_brier_res13 = auc_pr_brier_from_sl3(res13, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Ferratin deficiency", Xvars="All", population="children")
auc_pr_brier_res14 = auc_pr_brier_from_sl3(res14, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Ferratin deficiency", Xvars="Survey", population="children")
auc_pr_brier_res15 = auc_pr_brier_from_sl3(res15, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Ferratin deficiency", Xvars="Proxy", population="children")
auc_pr_brier_res16 = auc_pr_brier_from_sl3(res16, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Ferratin deficiency", Xvars="All", population="women")
auc_pr_brier_res17 = auc_pr_brier_from_sl3(res17, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Ferratin deficiency", Xvars="Survey", population="women")
auc_pr_brier_res18 = auc_pr_brier_from_sl3(res18, pred_source = c("auto"),ci = TRUE,quiet = TRUE) %>% mutate(Y="Ferratin deficiency", Xvars="Proxy", population="women")

res <- bind_rows(auc_pr_brier_res1, auc_pr_brier_res2, auc_pr_brier_res3, auc_pr_brier_res4, auc_pr_brier_res5, auc_pr_brier_res6,
                 auc_pr_brier_res7, auc_pr_brier_res8, auc_pr_brier_res9, auc_pr_brier_res10,
                 auc_pr_brier_res11, auc_pr_brier_res12, auc_pr_brier_res13, auc_pr_brier_res14,
                 auc_pr_brier_res15, auc_pr_brier_res16, auc_pr_brier_res17, auc_pr_brier_res18)
head(res)

saveRDS(res, file = here("results/res_full_bin_GW_Ghana_SL_all.rds"))
