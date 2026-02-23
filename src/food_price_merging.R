
rm(list=ls())
library(dplyr)
library(haven)
library(here)
library(purrr)
library(labelled)
library(sl3)
library(origami)
library(tlverse)
library(caret)
library(data.table)
library(lubridate)
library(ck37r)
library(maps)

dhs <- paste0("dhs_Ghana_2022.RDS")


ghana_fp <- read.csv("C:/Users/andre/OneDrive/Documents/mn-proxies/data/food_price/wfp_food_prices_gha.csv")
ghana_fp$year <- year(ymd(ghana_fp$date))
ghana_fp <- ghana_fp %>% filter(year==2022, pricetype=="Retail")
#transform to wide
ghana_fp_wide <- ghana_fp %>% spread( category)
ghana_fp_wide
head(ghana_fp)
dim(ghana_fp)
table(ghana_fp$priceflag )



ghana_wide <- ghana_fp %>%
  # Select only the columns we need
  select(date, admin1, admin2, market, latitude, longitude, commodity, usdprice) %>%
  # Spread the commodity column into multiple columns with usdprice values
  pivot_wider(
    names_from = commodity,
    values_from = usdprice,
    values_fill = NA
  )
write.csv(ghana_wide, "C:/Users/andre/OneDrive/Documents/mn-proxies/data/food_price/ghana_wide_example.csv")

# Assuming ghana_fp is already loaded in your R environment
# Extract unique GPS coordinates
unique_locations <- ghana_wide %>%
  select(market, latitude, longitude) %>%
  mutate(latitude=as.numeric(latitude),
         longitude=as.numeric(longitude)) %>%
  distinct()

# Get Ghana map data
ghana_map <- map_data("world") %>%
  filter(region == "Ghana")

# Create the plot
ggplot() +
  # Add Ghana map as background
  geom_polygon(data = ghana_map, aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "darkgray") +
  # Add the unique GPS points
  geom_point(data = unique_locations,
             aes(x = longitude, y = latitude, color = market),
             size = 3, alpha = 0.8) +
  # Add labels for markets
  geom_text(data = unique_locations,
            aes(x = longitude, y = latitude, label = market),
            size = 3, hjust = 0, vjust = -1) +
  # Set appropriate theme and labels
  theme_minimal() +
  labs(title = "Unique Market Locations in Ghana",
       x = "Longitude", y = "Latitude",
       color = "Market") +
  # Set appropriate map bounds for Ghana
  coord_map(xlim = c(-3.5, 1.5), ylim = c(4.5, 11.5))

