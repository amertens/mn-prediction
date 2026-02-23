
rm(list=ls())
# load  packages
library(rdhs)
library(ggplot2)


gambia_dhs <- readRDS(here("data/DHS/dhs_Gambia_2019.RDS"))
gambia_dhs$Gambia_2019$BRdata



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
labels_df <- get_var_labels(gambia_dhs$Gambia_2019$BRdata)
head(labels_df)


#Test
d2019 <- dhs_data(countryIds = "GM",
              indicatorIds = "FE_FRTR_W_A15",
              surveyYearStart = 2019,
              breakdown = "subnational")

# get our related spatial data frame object
sp2019 <- download_boundaries(surveyId = d2019$SurveyId[1])

m2019 <- d2019$Value[match(sp2019$sdr_subnational_boundaries$REG_ID, d2019$RegionId)]



dhs_indicators()
dhs_indicators_df <- dhs_indicators(returnFields=c("IndicatorId","ShortName","Label","Definition"," Level1", "Level2", "Level3", "MeasurementType", "Denominator"))
write.csv(dhs_indicators_df, file = here("data/DHS/clean/gambia_dhs_indicators_metadata.csv"), row.names = FALSE)

unique(dhs_indicators_df$Level1)
unique(dhs_indicators_df$Level2)
unique(dhs_indicators_df$Level3)

for(i in 1:nrow(dhs_indicators_df)) {
  m2019<-  NA
  cat(paste0(dhs_indicators_df$IndicatorId[i], ": ", dhs_indicators_df$Label[i]), "\n")
  # make request
  d2019 <- NULL
  try(
  d2019 <- dhs_data(countryIds = "GM",
                indicatorIds = dhs_indicators_df$IndicatorId[i],
                surveyYearStart = 2019,
                breakdown = "subnational"))

  if(!is.null(d2019) & length(unique(d2019$Value)) > 1) {
    # match our values to the regions

    m2019 <- d2019$Value[match(sp2019$sdr_subnational_boundaries$REG_ID, d2019$RegionId)]
    sp2019$sdr_subnational_boundaries$Value <- m2019
    colnames(sp2019$sdr_subnational_boundaries)[colnames(sp2019$sdr_subnational_boundaries) == "Value"] <- dhs_indicators_df$IndicatorId[i]

  }


}


#save the spatial data frames
saveRDS(sp2019$sdr_subnational_boundaries, file = here("data/DHS/clean/Gambia_2019_dhs_aggregation.rds"))

