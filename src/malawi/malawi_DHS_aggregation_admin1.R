
rm(list=ls())
# load  packages
library(rdhs)
library(ggplot2)


malawi_dhs <- readRDS(here("data/DHS/dhs_Malawi_2015.RDS"))
#malawi_dhs$Malawi_2015$MRdata



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
labels_df <- get_var_labels(malawi_dhs$malawi_2015$BRdata)
head(labels_df)


#Test
d2015 <- dhs_data(countryIds = "MW",
              indicatorIds = "FE_FRTR_W_A15",
              surveyYearStart = 2015,
              breakdown = "subnational")

# get our related spatial data frame object
sp2015 <- download_boundaries(surveyId = d2015$SurveyId[1])

m2015 <- d2015$Value[match(sp2015$sdr_subnational_boundaries$REG_ID, d2015$RegionId)]



dhs_indicators()
dhs_indicators_df <- dhs_indicators(returnFields=c("IndicatorId","ShortName","Label","Definition"," Level1", "Level2", "Level3", "MeasurementType", "Denominator"))
write.csv(dhs_indicators_df, file = here("data/DHS/clean/malawi_dhs_indicators_metadata.csv"), row.names = FALSE)

unique(dhs_indicators_df$Level1)
unique(dhs_indicators_df$Level2)
unique(dhs_indicators_df$Level3)

for(i in 1:nrow(dhs_indicators_df)) {
  m2015<-  NA
  cat(paste0(dhs_indicators_df$IndicatorId[i], ": ", dhs_indicators_df$Label[i]), "\n")
  # make request
  d2015 <- NULL
  try(
  d2015 <- dhs_data(countryIds = "MW",
                indicatorIds = dhs_indicators_df$IndicatorId[i],
                surveyYearStart = 2015,
                breakdown = "subnational"))

  if(!is.null(d2015) & length(unique(d2015$Value)) > 1) {
    # match our values to the regions

    m2015 <- d2015$Value[match(sp2015$sdr_subnational_boundaries$REG_ID, d2015$RegionId)]
    sp2015$sdr_subnational_boundaries$Value <- m2015
    colnames(sp2015$sdr_subnational_boundaries)[colnames(sp2015$sdr_subnational_boundaries) == "Value"] <- dhs_indicators_df$IndicatorId[i]

  }


}


#save the spatial data frames
saveRDS(sp2015$sdr_subnational_boundaries, file = here("data/DHS/clean/malawi_2015_dhs_aggregation.rds"))

