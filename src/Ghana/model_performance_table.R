


rm(list=ls())
library(tidyverse)
library(here)


res1 <- readRDS(here("results/models/res_GW_Ghana_SL_child_vitA_full.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Vit-A", Xvars="All", population="children")
res2 <- readRDS(here("results/models/res_GW_Ghana_SL_child_vitA_gwPreds.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Vit-A", Xvars="Survey", population="children")
res3 <- readRDS(here("results/models/res_GW_Ghana_SL_child_vitA.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Vit-A", Xvars="Proxy", population="children")
res4 <- readRDS(here("results/models/res_GW_Ghana_SL_women_vitA_full.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Vit-A", Xvars="All", population="women")
res5 <- readRDS(here("results/models/res_GW_Ghana_SL_women_vitA_gwPreds.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Vit-A", Xvars="Survey", population="women")
res6 <- readRDS(here("results/models/res_GW_Ghana_SL_women_vitA.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Vit-A", Xvars="Proxy", population="women")
res7 <- readRDS(here("results/models/res_GW_Ghana_SL_women_b12_full.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Vit-B12", Xvars="All", population="women")
res8 <- readRDS(here("results/models/res_GW_Ghana_SL_women_b12_gwPreds.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Vit-B12", Xvars="Survey", population="women")
res9 <- readRDS(here("results/models/res_GW_Ghana_SL_women_b12.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Vit-B12", Xvars="Proxy", population="women")
res10 <- readRDS(here("results/models/res_GW_Ghana_SL_women_folate_full.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Folate", Xvars="All", population="women")
res11 <- readRDS(here("results/models/res_GW_Ghana_SL_women_folate_gwPreds.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Folate", Xvars="Survey", population="women")
res12 <- readRDS(here("results/models/res_GW_Ghana_SL_women_folate.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Folate", Xvars="Proxy", population="women")
res13 <- readRDS(here("results/models/res_GW_Ghana_SL_child_iron_full.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Ferratin", Xvars="All", population="children")
res14 <- readRDS(here("results/models/res_GW_Ghana_SL_child_iron_gwPreds.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Ferratin", Xvars="Survey", population="children")
res15 <- readRDS(here("results/models/res_GW_Ghana_SL_child_iron.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Ferratin", Xvars="Proxy", population="children")
res16 <- readRDS(here("results/models/res_GW_Ghana_SL_mom_iron_full.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Ferratin", Xvars="All", population="women")
res17 <- readRDS(here("results/models/res_GW_Ghana_SL_mom_iron_gwPreds.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Ferratin", Xvars="Survey", population="women")
res18 <- readRDS(here("results/models/res_GW_Ghana_SL_mom_iron.rds"))$cv_risk_w_sl_revere %>% mutate(Y="Ferratin", Xvars="Proxy", population="women")


res <- bind_rows(res1, res2, res3, res4, res5, res6, res7, res8, res9, res10, res11, res12, res13, res14, res15, res16, res17, res18)
head(res)

unique(res$learner)
unique(res$Y)

plotdf <- res %>% filter(learner=="SuperLearner")
ggplot(plotdf, aes(x=sqrt(MSE), y=paste0(population,"-", Xvars))) + geom_point() + facet_wrap(~Y, scales="free")


plotdf <- res %>% filter(Xvars=="All", population=="women", Y=="Vit-A")
ggplot(plotdf, aes(x=sqrt(MSE), y=learner)) + geom_point() + facet_wrap(~Y, scales="free")

plotdf <- res %>% filter(Xvars=="Proxy", population=="women", Y=="Vit-A")
ggplot(plotdf, aes(x=sqrt(MSE), y=learner)) + geom_point() + facet_wrap(~Y, scales="free") + coord_cartesian(xlim=c(0.575, 0.595))



plotdf <- res %>% filter(Xvars=="Proxy", population=="children", Y=="Vit-A")
ggplot(plotdf, aes(x=sqrt(MSE), y=learner)) + geom_point() + facet_wrap(~Y, scales="free") + coord_cartesian(xlim=c(0.31, 0.35))

(plotdf$MSE[1] -plotdf$MSE[length(plotdf$MSE)])/plotdf$MSE[1] *100
