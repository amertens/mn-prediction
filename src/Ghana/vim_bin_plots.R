


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

metadata= readRDS(here("metadata/variable_categories.rds"))

for(i in 1:length(metadata)){
  cat(names(metadata)[i],": ", length(metadata[[i]]),"\n")
}

res = readRDS(here("results/models/res_full_bin_GW_Ghana_SL_women_folate_full.rds"))
dhs_indicators = read.csv( here("data/DHS/clean/dhs_indicators_metadata.csv"))
unique(dhs_indicators$Level1)
unique(dhs_indicators$Level2)
unique(dhs_indicators$Level2)[grepl("nutr",unique(dhs_indicators$Level2))]

dhs_indicators %>% filter(Level2 %in% c("Micronutrient intake among children", "Micronutrient intake among mothers" )) %>% select(Label)

dhs_indicators %>% filter(IndicatorId=="ML_AMLD_C_CHL")  %>% select(Label)
#Children who took Chloroquine


#NOTE! Need to update the variable importance function to have survey variables groups,
#otherwise Missingness is capturing the GW vars

#child_vitA_full_vim=readRDS(file = here("results/child_vitA_def_varimp_group.rds"))
child_vitA_vim=readRDS(file = here("results/child_vitA_def_proxy_varimp_group.rds"))
#child_iron_full_vim=readRDS(file = here("results/child_iron_def_varimp_group.rds"))
child_iron_vim=readRDS(file = here("results/child_iron_def_proxy_varimp_group.rds"))
#women_vitA_full_vim= readRDS(file = here("results/women_vitA_def_varimp_group.rds"))
women_vitA_vim= readRDS(file = here("results/women_vitA_def_proxy_varimp_group.rds"))
#women_b12_full_vim= readRDS(file = here("results/women_b12_def_varimp_group.rds"))
women_b12_vim= readRDS(file = here("results/women_b12_def_proxy_varimp_group.rds"))
#women_iron_full_vim= readRDS(file = here("results/women_iron_def_varimp_group.rds"))
women_iron_vim= readRDS(file = here("results/women_iron_def_proxy_varimp_group.rds"))
#women_folate_full_vim= readRDS(file = here("results/women_folate_def_varimp_group.rds"))
women_folate_vim= readRDS(file = here("results/women_folate_def_proxy_varimp_group.rds"))

plot_vim <- function(vim, title="", top_n=5){

  plotdf=vim$varimp %>% filter(NLL_difference>0) %>%
    arrange(NLL_difference) %>%
    mutate(covariate_group=factor(covariate_group , levels=unique(covariate_group )))

  if(!is.null(top_n)){
    plotdf=plotdf %>% arrange(desc(NLL_difference)) %>% slice(1:top_n)
  }
  p = ggplot(plotdf,
             aes(x=covariate_group, y=NLL_difference)) +
    geom_bar(stat="identity", fill="steelblue")+
    coord_flip()+
    labs(x="", y="", title=title)+
    theme_minimal() +
    theme(plot.title = element_text(size=16, face='bold', hjust = 0.5),
          axis.text.y = element_text(size=12))
  return(p)
}


p1=plot_vim(child_vitA_vim, title="Child Vitamin A Deficiency")
p2=plot_vim(women_vitA_vim, title="Womens' Vitamin A Deficiency")
p3=plot_vim(child_iron_vim, title="Child Iron Deficiency")
p4=plot_vim(women_iron_vim, title="Womens' Iron Deficiency")
p5=plot_vim(women_b12_vim, title="Womens' B12 Deficiency")
p6=plot_vim(women_folate_vim, title="Womens' Folate Deficiency")


library(patchwork)
# (p1 | p2) / (p3 | p4) / (p5 | p6)

#combine plots into a single grid with one X-axis "Change in Negative Log Likelihood
p <- p1 + theme(axis.title.x=element_blank()) + p2 + theme(axis.title.x=element_blank()) +
  p3 + theme(axis.title.x=element_blank()) + p4 + theme(axis.title.x=element_blank()) +
  p5 + theme(axis.title.x=element_blank()) + p6 + theme(axis.title.x=element_blank()) +
  plot_layout(ncol=2, guides = 'collect') & theme(legend.position='bottom') &
  plot_annotation(#title = 'Top 5 Most Important Proxy Variable Domains by Outcome',
    caption = 'Variable importance measured by change in negative log likelihood ') &
  theme(plot.title = element_text(size=16, face='bold', hjust = 0.5),
        plot.caption = element_text(size=14, face='bold', hjust = 0.5)) &
  ylab("Change in Negative Log Likelihood")
p

ggsave(p,
       file= here("figures/varimp_bin_proxy_def_all_outcomes.png"),
       width = 14.03,
       height = 5.16,
       units = "in")

