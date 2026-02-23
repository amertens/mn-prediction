
rm(list=ls())
# load  packages
library(rdhs)
library(ggplot2)


ghana_dhs <- readRDS(here("data/DHS/dhs_Ghana_2014.RDS"))
ghana_dhs$Ghana_2014$BRdata



# Function to extract names + labels
get_var_labels <- function(df) {
  data.frame(
    name  = names(df),
    label = sapply(df, function(x) {
      lab <- attr(x, "label", exact = TRUE)
      if (is.null(lab)) "" else lab
    }),
    stringsAsFactors = FALSE
  )
}

# look at diet data and make sure it's being aggregated
labels_df <- get_var_labels(ghana_dhs$Ghana_2014$BRdata)
head(labels_df)


#Test
d2014 <- dhs_data(countryIds = "GH",
              indicatorIds = "FE_FRTR_W_A15",
              surveyYearStart = 2014,
              breakdown = "subnational")
d2016 <- dhs_data(countryIds = "GH",
                  indicatorIds = "FE_FRTR_W_A15",
                  surveyYearStart = 2016,
                  breakdown = "subnational")
d2017 <- dhs_data(countryIds = "GH",
                  indicatorIds = "FE_FRTR_W_A15",
                  surveyYearStart = 2017,
                  breakdown = "subnational")

# get our related spatial data frame object
sp2014 <- download_boundaries(surveyId = d2014$SurveyId[1])
sp2016 <- download_boundaries(surveyId = d2016$SurveyId[1])
sp2017 <- download_boundaries(surveyId = d2017$SurveyId[1])

m2014 <- d2014$Value[match(sp2014$sdr_subnational_boundaries$REG_ID, d2014$RegionId)]
m2016 <- d2016$Value[match(sp2016$sdr_subnational_boundaries$REG_ID, d2016$RegionId)]
m2017 <- d2017$Value[match(sp2017$sdr_subnational_boundaries$REG_ID, d2017$RegionId)]



dhs_indicators()
dhs_indicators_df <- dhs_indicators(returnFields=c("IndicatorId","ShortName","Label","Definition"," Level1", "Level2", "Level3", "MeasurementType", "Denominator"))
write.csv(dhs_indicators_df, file = here("data/DHS/clean/dhs_indicators_metadata.csv"), row.names = FALSE)

unique(dhs_indicators_df$Level1)
unique(dhs_indicators_df$Level2)
unique(dhs_indicators_df$Level3)

for(i in 1:nrow(dhs_indicators_df)) {
  m2014<- m2016 <- m2017 <- NA
  cat(paste0(dhs_indicators_df$IndicatorId[i], ": ", dhs_indicators_df$Label[i]), "\n")
  # make request
  d2014 <-d2016 <-d2017 <- NULL
  try(
  d2014 <- dhs_data(countryIds = "GH",
                indicatorIds = dhs_indicators_df$IndicatorId[i],
                surveyYearStart = 2014,
                breakdown = "subnational"))
  try(
    d2016 <- dhs_data(countryIds = "GH",
                      indicatorIds = dhs_indicators_df$IndicatorId[i],
                      surveyYearStart = 2016,
                      breakdown = "subnational"))
  try(
    d2017 <- dhs_data(countryIds = "GH",
                      indicatorIds = dhs_indicators_df$IndicatorId[i],
                      surveyYearStart = 2017,
                      breakdown = "subnational"))
  if(!is.null(d2014) & length(unique(d2014$Value)) > 1) {
    # match our values to the regions

    m2014 <- d2014$Value[match(sp2014$sdr_subnational_boundaries$REG_ID, d2014$RegionId)]
    sp2014$sdr_subnational_boundaries$Value <- m2014
    colnames(sp2014$sdr_subnational_boundaries)[colnames(sp2014$sdr_subnational_boundaries) == "Value"] <- dhs_indicators_df$IndicatorId[i]

  }

  if(!is.null(d2016) & length(unique(d2016$Value)) > 1) {
    # match our values to the regions

    m2016 <- d2016$Value[match(sp2016$sdr_subnational_boundaries$REG_ID, d2016$RegionId)]
    sp2016$sdr_subnational_boundaries$Value <- m2016
    colnames(sp2016$sdr_subnational_boundaries)[colnames(sp2016$sdr_subnational_boundaries) == "Value"] <- dhs_indicators_df$IndicatorId[i]

  }

  if(!is.null(d2017) & length(unique(d2017$Value)) > 1) {
    # match our values to the regions

    m2017 <- d2017$Value[match(sp2017$sdr_subnational_boundaries$REG_ID, d2017$RegionId)]
    sp2017$sdr_subnational_boundaries$Value <- m2017
    colnames(sp2017$sdr_subnational_boundaries)[colnames(sp2017$sdr_subnational_boundaries) == "Value"] <- dhs_indicators_df$IndicatorId[i]

  }


}




# # Use ggplot and geom_sf to plot
# head(sp$sdr_subnational_boundaries)[,1:30]
# sp$sdr_subnational_boundaries$SV_BACK_M_UNW
#
# indicator="SV_BACK_M_UNW"
# ggplot(sp$sdr_subnational_boundaries) +
#   #geom_sf(aes(fill = !!(indicator)   )) +
#   geom_sf(aes(fill =  SV_BACK_M_UNW  )) +
#   ggtitle(dhs_indicators_df$Definition[dhs_indicators_df$IndicatorId == indicator])


#save the spatial data frames
saveRDS(sp2014$sdr_subnational_boundaries, file = here("data/DHS/clean/Ghana_2014_dhs_aggregation.rds"))
saveRDS(sp2016$sdr_subnational_boundaries, file = here("data/DHS/clean/Ghana_2016_dhs_aggregation.rds"))
saveRDS(sp2017$sdr_subnational_boundaries, file = here("data/DHS/clean/Ghana_2017_dhs_aggregation.rds"))

