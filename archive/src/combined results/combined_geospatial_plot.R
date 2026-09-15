resfull =readRDS(file = here("results/models/res_GW_transportability_SL_women_vitA_geospatial_only.rds"))
res_ghana =readRDS(file = here("results/models/res_GW_Ghana_SL_women_vitA_geospatial_only.rds"))
res_gambia =readRDS(file = here("results/models/res_GW_Gambia_SL_women_vitA_geospatial_only.rds"))
res_malawi =readRDS(file = here("results/models/res_GW_Malawi_SL_women_vitA_geospatial_only.rds"))
res_sierraleone =readRDS(file = here("results/models/res_GW_SierraLeone_SL_women_vitA_geospatial_only.rds"))


res_ghana_full =readRDS(file = here("results/models/res_GW_Ghana_SL_women_vitA_V2.rds"))
res_gambia_full =readRDS(file = here("results/models/res_GW_Gambia_SL_women_vitA.rds"))
res_malawi_full =readRDS(file = here("results/models/res_Malawi_SL_women_vitA.rds"))
res_sierraleone_full =readRDS(file = here("results/models/res_GW_Sierra_Leone_SL_women_vitA.rds"))

plotdf= bind_rows(
  data.frame(country="Ghana", Covariates="Geospatial_only", MSE=res_ghana$cv_risk_w_sl_revere$MSE[nrow(res_ghana$cv_risk_w_sl_revere)], SE=res_ghana$cv_risk_w_sl_revere$se[nrow(res_ghana$cv_risk_w_sl_revere)]),
  data.frame(country="Gambia", Covariates="Geospatial_only", MSE=res_gambia$cv_risk_w_sl_revere$MSE[nrow(res_gambia$cv_risk_w_sl_revere)], SE=res_gambia$cv_risk_w_sl_revere$se[nrow(res_gambia$cv_risk_w_sl_revere)]),
  data.frame(country="Malawi", Covariates="Geospatial_only", MSE=res_malawi$cv_risk_w_sl_revere$MSE[nrow(res_malawi$cv_risk_w_sl_revere)], SE=res_malawi$cv_risk_w_sl_revere$se[nrow(res_malawi$cv_risk_w_sl_revere)]),
  data.frame(country="Sierra Leone", Covariates="Geospatial_only", MSE=res_sierraleone$cv_risk_w_sl_revere$MSE[nrow(res_sierraleone$cv_risk_w_sl_revere)], SE=res_sierraleone$cv_risk_w_sl_revere$se[nrow(res_sierraleone$cv_risk_w_sl_revere)]),
  data.frame(country="Ghana", Covariates="Full covariates", MSE=res_ghana_full$cv_risk_w_sl_revere$MSE[nrow(res_ghana_full$cv_risk_w_sl_revere)], SE=res_ghana_full$cv_risk_w_sl_revere$se[nrow(res_ghana_full$cv_risk_w_sl_revere)]),
  data.frame(country="Gambia", Covariates="Full covariates", MSE=res_gambia_full$cv_risk_w_sl_revere$MSE[nrow(res_gambia_full$cv_risk_w_sl_revere)], SE=res_gambia_full$cv_risk_w_sl_revere$se[nrow(res_gambia_full$cv_risk_w_sl_revere)]),
  data.frame(country="Malawi", Covariates="Full covariates", MSE=res_malawi_full$cv_risk_w_sl_revere$MSE[nrow(res_malawi_full$cv_risk_w_sl_revere)], SE=res_malawi_full$cv_risk_w_sl_revere$se[nrow(res_malawi_full$cv_risk_w_sl_revere)]),
  data.frame(country="Sierra Leone", Covariates="Full covariates", MSE=res_sierraleone_full$cv_risk_w_sl_revere$MSE[nrow(res_sierraleone_full$cv_risk_w_sl_revere)], SE=res_sierraleone_full$cv_risk_w_sl_revere$se[nrow(res_sierraleone_full$cv_risk_w_sl_revere)]),
  data.frame(country="All", Covariates="Geospatial_only", MSE=resfull$cv_risk_w_sl_revere$MSE[nrow(resfull$cv_risk_w_sl_revere)], SE=resfull$cv_risk_w_sl_revere$se[nrow(resfull$cv_risk_w_sl_revere)])
)


plotdf <- plotdf %>%
  mutate(
    ci_lwr = MSE - 1.96 * SE,
    ci_upr = MSE + 1.96 * SE
  )

plotdf

ggplot(plotdf, aes(x=country, y=MSE, fill=Covariates, color=Covariates)) +
  geom_point( position=position_dodge(.5)) +
  geom_errorbar(aes(ymin=ci_lwr, ymax=ci_upr), width=.2,
                position=position_dodge(.5)) +
  labs(title="Model Performance by Country and Covariates",
       x="Country",
       y="Mean Squared Error (MSE)") +
  theme_minimal()

#make plot with and without the pooled performance

#need to get the national prevalence from the external models
