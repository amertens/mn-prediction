


rm(list=ls())
library(tidyverse)
library(here)
library(pROC)

source(paste0(here::here(),"/src/0-functions.R"))

res_simple = readRDS( here("results/res_bin_GW_Ghana_SL_all.rds"))
res = readRDS( here("results/res_full_bin_GW_Ghana_SL_all.rds"))

unique(res$Y)

head(res)


#format plot data
plotdf <- res

unique(plotdf$Xvars)
unique(plotdf$population)


plotdf <- plotdf %>%
  mutate(Y = case_when(
    Y == "Vit-A deficiency" ~ "Vitamin A deficiency",
    Y == "Vit-B12 deficiency" ~ "Vitamin B12 deficiency",
    Y == "Folate deficiency" ~ "Folate deficiency",
    Y == "Ferratin deficiency" ~ "Iron deficiency",
    TRUE ~ Y
  ),
  Y = factor(Y, levels=c("Vitamin A deficiency", "Iron deficiency", "Folate deficiency", "Vitamin B12 deficiency")),
  Xvars = case_when(
    Xvars == "Survey" ~ "Survey only",
    Xvars == "Proxy" ~ "Proxies",
    TRUE ~ Xvars
  ),
  population = str_to_title(population),
  facet_var=paste0(population, " - ", Xvars),
  facet_var=factor(facet_var, levels=rev(c("Children - Survey only","Children - All", "Children - Proxies","Women - Survey only","Women - All", "Women - Proxies" )))
  )

unique(plotdf$facet_var)


table(plotdf$Xvars)
plotdf1 <- plotdf %>%
  mutate(brier_skill=ifelse(Xvars=="Survey only", brier_skill, NA),
         brier_skill_ci_low=ifelse(Xvars=="Survey only", brier_skill_ci_low, NA),
         brier_skill_ci_high=ifelse(Xvars=="Survey only", brier_skill_ci_high, NA))

plotdf2 <- plotdf %>%
  mutate(brier_skill=ifelse(Xvars!="Proxies", brier_skill, NA),
         brier_skill_ci_low=ifelse(Xvars!="Proxies", brier_skill_ci_low, NA),
         brier_skill_ci_high=ifelse(Xvars!="Proxies", brier_skill_ci_high, NA))



p1 <- ggplot(plotdf1, aes(x = brier_skill * 100,
                         y = facet_var)) +
  geom_point(size = 3, color = "#2C3E50") +
  geom_errorbar(aes(xmin = brier_skill_ci_low * 100,
                    xmax = brier_skill_ci_high * 100),
                width = 0.2, color = "#2C3E50") +
  facet_wrap(~Y, scales = "free") +
  labs(
    x = "Percent improvement over null model (Brier Skill Score)",
    y = "Population and predictor set",
    title = "Model results by outcome and predictor set"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),   # remove horizontal grid lines
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),     # remove grey facet strip box
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")  # ref line at null
p1


p2 <- ggplot(plotdf2, aes(x = brier_skill * 100,
                          y = facet_var, color=facet_var)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = brier_skill_ci_low * 100,
                    xmax = brier_skill_ci_high * 100),
                width = 0.2) +
  scale_color_manual(values=rep(c("#2C3E50", "#E74C3C","#2C3E50"),2)) +
  facet_wrap(~Y, scales = "free") +
  labs(
    x = "Percent improvement over null model (Brier Skill Score)",
    y = "Population and predictor set",
    title = "Model results by outcome and predictor set"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),   # remove horizontal grid lines
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),     # remove grey facet strip box
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")  # ref line at null
p2


p3 <- ggplot(plotdf, aes(x = brier_skill * 100,
                         y = facet_var, color=facet_var)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = brier_skill_ci_low * 100,
                    xmax = brier_skill_ci_high * 100),
                width = 0.2) +
  scale_color_manual(values=rep(c("#35b779", "#E74C3C","#2C3E50"),2)) +

  facet_wrap(~Y, scales = "free") +
  labs(
    x = "Percent improvement over null model (Brier Skill Score)",
    y = "Population and predictor set",
    title = "Model results by outcome and predictor set"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),   # remove horizontal grid lines
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),     # remove grey facet strip box
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")  # ref line at null
p3
ggsave(here("figures/model_performance_bin_brier_skill_all1.png"), p1, width=8, height=4)
ggsave(here("figures/model_performance_bin_brier_skill_all2.png"), p2, width=8, height=4)
ggsave(here("figures/model_performance_bin_brier_skill_all3.png"), p3, width=8, height=4)

ggplot(plotdf, aes(x=roc_auc , y=paste0(population,"-", Xvars))) +
  geom_point() +
  geom_errorbar(aes(xmin=roc_ci_low , xmax=roc_ci_high), width=0.1) +
  geom_point(aes(x=pr_baseline ), color="blue") +
  geom_point(aes(x=pr_auc ), color="red") +
  geom_errorbar(aes(xmin=pr_ci_low  , xmax=pr_ci_high), color="red", width=0.1) +
  geom_vline(xintercept = 0.5, linetype="dashed") +
  facet_wrap(~Y, scales="free") + theme_bw()

ggplot(plotdf%>% filter(Xvars=="Proxies"), aes(x=roc_auc, y=paste0(population,"-", Xvars))) +
  geom_point() +
  geom_errorbar(aes(xmin=roc_ci_low , xmax=roc_ci_high), width=0.1) +
  geom_point(aes(x=pr_baseline ), color="blue") +
  geom_point(aes(x=pr_auc ), color="red") +
  geom_errorbar(aes(xmin=pr_ci_low  , xmax=pr_ci_high), color="red", width=0.1) +
  geom_vline(xintercept = 0.5, linetype="dashed") +
  facet_wrap(~Y, scales="free") + theme_bw()


#
# ggplot(plotdf, aes(x=brier_skill*100 , y=paste0(population,"-", Xvars))) +
#   geom_point() +
#   geom_errorbar(aes(xmin=brier_skill_ci_low*100, xmax=brier_skill_ci_high*100), width=0.1) +
#   facet_wrap(~Y, scales="free") + theme_bw()
#
#
#
#
# plotdf2 <- res %>% filter(Xvars=="Proxy")
# ggplot(plotdf2, aes(x=roc_auc, y=Y)) +
#   geom_point() +
#   geom_point(aes(x=pr_auc ), color="red") +
#   geom_vline(xintercept = 0.5, linetype="dashed") +
#   geom_errorbar(aes(xmin=roc_ci_low , xmax=roc_ci_high), width=0.1) +
#   geom_errorbar(aes(xmin=pr_ci_low  , xmax=pr_ci_high), color="red", width=0.1) +
#   facet_wrap(~paste0(population,"-", Xvars), scales="free") + theme_bw()
#
# ggplot(plotdf2, aes(x=brier_skill*100 , y=paste0(population,"-", Xvars))) +
#   geom_point() +
#   geom_errorbar(aes(xmin=brier_skill_ci_low*100, xmax=brier_skill_ci_high*100), width=0.1) +
#   facet_wrap(~Y, scales="free") + theme_bw()
#
#
# plotdf3 <- res %>% filter(Xvars!="Proxy")
#
# ggplot(plotdf3, aes(x=brier_skill*100 , y=paste0(population,"-", Xvars))) +
#   geom_point() +
#   geom_errorbar(aes(xmin=brier_skill_ci_low*100, xmax=brier_skill_ci_high*100), width=0.1) +
#   facet_wrap(~Y, scales="free") + theme_bw()
#
#
#
#
# #gain curves
#
# gain_curve_df <- function(y, p, model, grid = seq(0.01, 1, by = 0.01)) {
#   o <- order(p, decreasing = TRUE)
#   y <- y[o]
#   n <- length(y); n_pos <- sum(y == 1)
#   stopifnot(n_pos > 0)
#   k <- ceiling(grid * n)
#   cum_tp <- cumsum(y)[k]
#   data.frame(model, target_frac = grid, recall = cum_tp / n_pos)
# }
#
# Y <- unclass(res1$Y)
# P <- res1$yhat_full
#
#
# dfA <- gain_curve_df(Y, P, model = "SL (OOF)")
# ggplot(dfA, aes(target_frac, recall, color = model)) +
#   geom_line(size = 1) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   scale_x_continuous(labels = scales::percent) +
#   scale_y_continuous(labels = scales::percent) +
#   labs(x = "% of population targeted", y = "Recall (cases captured)",
#        title = "Capture curve: Vitamin A deficiency") +
#   theme_minimal()
#
# P_modelA <- res1$yhat_full
# P_modelB <- res2$yhat_full
# P_modelC <- res3$yhat_full
#
# dfs <- rbind(
#   gain_curve_df(Y, P_modelA, "Model A"),
#   gain_curve_df(Y, P_modelB, "Model B"),
#   gain_curve_df(Y, P_modelC, "Model C")
# )
# ggplot(dfs, aes(target_frac, recall, color = model)) + geom_line() + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
