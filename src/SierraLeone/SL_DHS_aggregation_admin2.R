
rm(list=ls())
# load  packages
library(rdhs)
library(ggplot2)


SL_dhs <- readRDS(here("data/DHS/dhs_Sierra Leone_2013.RDS"))
SL_dhs$`Sierra Leone_2013`$BRdata



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
labels_df <- get_var_labels(SL_dhs$`Sierra Leone_2013`$BRdata)
head(labels_df)


#Test
d2013 <- dhs_data(countryIds = "SL",
              indicatorIds = "FE_FRTR_W_A15",
              surveyYearStart = 2013,
              breakdown = "subnational")

# get our related spatial data frame object
sp2013 <- download_boundaries(surveyId = d2013$SurveyId[1])

m2013 <- d2013$Value[match(sp2013$sdr_subnational_boundaries$REG_ID, d2013$RegionId)]



dhs_indicators()
dhs_indicators_df <- dhs_indicators(returnFields=c("IndicatorId","ShortName","Label","Definition"," Level1", "Level2", "Level3", "MeasurementType", "Denominator"))
write.csv(dhs_indicators_df, file = here("data/DHS/clean/SL_dhs_indicators_metadata.csv"), row.names = FALSE)

unique(dhs_indicators_df$Level1)
unique(dhs_indicators_df$Level2)
unique(dhs_indicators_df$Level3)

for(i in 1:nrow(dhs_indicators_df)) {
  m2013<-  NA
  cat(paste0(dhs_indicators_df$IndicatorId[i], ": ", dhs_indicators_df$Label[i]), "\n")
  # make request
  d2013 <- NULL
  try(
  d2013 <- dhs_data(countryIds = "SL",
                indicatorIds = dhs_indicators_df$IndicatorId[i],
                surveyYearStart = 2013,
                breakdown = "subnational"))

  if(!is.null(d2013) & length(unique(d2013$Value)) > 1) {
    # match our values to the regions

    m2013 <- d2013$Value[match(sp2013$sdr_subnational_boundaries$REG_ID, d2013$RegionId)]
    sp2013$sdr_subnational_boundaries$Value <- m2013
    colnames(sp2013$sdr_subnational_boundaries)[colnames(sp2013$sdr_subnational_boundaries) == "Value"] <- dhs_indicators_df$IndicatorId[i]

  }


}


#save the spatial data frames
saveRDS(sp2013$sdr_subnational_boundaries, file = here("data/DHS/clean/SL_2013_dhs_aggregation_admin2.rds"))

