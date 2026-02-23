
# Load required libraries
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(haven)
library(here)
library(patchwork)

d <- read_dta(here("data/VMNIS/VMNISIndicator_long_format.dta"))
head(d)

unique(d$Country)


#filter to LMIC's in Africa, Asia, and South America only

LMICs = c("Malawi", "Mozambique", "Zambia", "Zimbabwe", "Angola", "Benin",
          "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
          "Chad", "Comoros", "Congo", "Côte d'Ivoire", "Democratic Republic of the Congo",
          "Djibouti", "Equatorial Guinea", "Eritrea", "Ethiopia", "Gabon",
          "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Kenya", "Lesotho",
          "Liberia", "Madagascar", "Mali", "Mauritania", "Mauritius", "Namibia",
          "Niger", "Nigeria", "Rwanda", "Sao Tome and Principe", "Senegal",
          "Seychelles", "Sierra Leone", "Somalia", "South Africa", "South Sudan",
          "Sudan", "Swaziland", "Tanzania", "Togo", "Uganda", "Afghanistan",
          "Bangladesh", "Bhutan", "Cambodia", "China", "Democratic People's Republic of Korea",
          "India", "Indonesia", "Iran (Islamic Republic of)", "Iraq", "Jordan",
          "Kiribati", "Kyrgyzstan", "Lao People's Democratic Republic",
          "Lebanon", "Malaysia", "Maldives", "Marshall Islands", "Micronesia (Federated States of)",
          "Mongolia", "Myanmar", "Nauru", "Nepal", "Pakistan", "Palau",
          "Papua New Guinea", "Philippines", "Samoa", "Solomon Islands",
          "Sri Lanka", "Tajikistan", "Thailand", "Timor-Leste", "Tonga",
          "Turkmenistan", "Tuvalu", "Uzbekistan", "Vanuatu", "Vietmam",
          "Yemen", "Albania", "Armenia", "Azerbaijan","Belize","Colombia",
          "Costa Rica","Ecuador", "Guatemala" , "Mexico", "Panama", "West Bank and Gaza Strip",
          "Viet Nam", "Peru", "United Republic of Tanzania","Bolivia (Plurinational State of)", "Dominican Republic",
          "Guyana", "Haiti", "Jamaica",  "Argentina", "Brazil", "Chile", "Paraguay",  "Morocco",
          "Nicaragua", "Egypt",  "Brazil", "El Salvador", "Tunisia",  "Botswana" , "Kazakhstan" , "Algeria","Paraguay",  "Venezuela (Bolivarian Republic of)",
          "Syrian Arab Republic" )


unique(d$Country)[!(unique(d$Country) %in% LMICs)]

unique(d$Population)

plotdf <- d %>% filter(
  Representativeness=="national",Areacovered=="both urban and rural",
  Country %in% LMICs) %>%
  mutate(Population=case_when(Population %in% c("Infants (1-11 months)","Children (12-23 months)","Neonates (0-28 days)",
                                                "Young infants (0-5 months)","Children (24-59 months)" ) ~"Children <5",
         Population==Population~Population))
head(plotdf)

#to do:
#collapse MN categories and populations, and drop LMIC countries

survey_summary <- plotdf %>%
  group_by(Indicator, Country, Population) %>%
  filter(Date == max(Date)) %>%
  summarize(n_surveys = n(),
            mean=mean(Mean,na.rm=T)) %>%
  filter(!is.na(mean)) %>%
  ungroup() %>% droplevels()
head(survey_summary)

#plotdf = survey_summary %>% filter(Population=="Pregnant women")
write.csv(plotdf, here("data/VMNIS/vmnis_heatmap_data.csv"), row.names = F)


  # ------------------------------------------------------------------
  # 1.  Recode Indicator → mn_group  (broad micronutrient categories)
  # ------------------------------------------------------------------
mn_map <- function(x) {
  case_when(
    x %in% c("Ferritin", "Serum transferrin receptor",
             "Iron (serum or plasma)", "Total iron binding capacity (TIBC)",
             "Transferrin saturation", "Mean corpuscular volume (MCV)")        ~ "Iron",
    x %in% c("Retinol (plasma or serum)", "Beta-carotene")                     ~ "Vitamin A",
    x == "Zinc (plasma or serum)"                                              ~ "Zinc",
    x %in% c("Goitre by ultrasound", "Iodine")                                 ~ "Iodine",
    x %in% c("Folate (plasma or serum)", "Folate (red blood cell)")            ~ "Folate",
    x == "Vitamin B12"                                                         ~ "Vitamin B12",
    x == "25-Hydroxyvitamin D"                                                 ~ "Vitamin D",
    x == "Vitamin E (alpha-tocopherol)"                                        ~ "Vitamin E",
    x == "Vitamin C"                                                           ~ "Vitamin C",
    x == "Haemoglobin"                                                         ~ "Haemoglobin",
    TRUE                                                                       ~ NA_character_
  )
}

plotdf <- d %>%
  filter(
    Representativeness == "national",
    Areacovered        == "both urban and rural",
    Country            %in% LMICs
  ) %>%
  mutate(
    Population = case_when(
      Population %in% c("Infants (1-11 months)", "Children (12-23 months)",
                        "Neonates (0-28 days)",  "Young infants (0-5 months)",
                        "Children (24-59 months)") ~ "Children <5",
      TRUE ~ Population),
    mn_group  = mn_map(Indicator)
  ) %>%
  drop_na(mn_group)                                            # toss very rare / misc. assays

# ------------------------------------------------------------------
# 2.  Keep the most‑recent survey per Country × Population × micronutrient
# ------------------------------------------------------------------
survey_summary <- plotdf %>%
  group_by(Indicator  ) %>%
  mutate(Mean=scale(Mean)) %>%
  group_by(Country, Population, mn_group) %>%
  filter(Date == max(Date)) %>%            # latest survey for that triad
  summarise(
    n_surveys = n(),
    mean      = mean(Mean, na.rm = TRUE),  # or whatever statistic you prefer
    .groups   = "drop"
  ) %>%
  filter(!is.na(mean))

# ------------------------------------------------------------------
# 3.  OPTIONAL: focus on one population
# ------------------------------------------------------------------
plotdf_pw <- survey_summary %>% filter(Population == "Pregnant women")

# ------------------------------------------------------------------
# 4.  Heat‑map style dot‑plot
# ------------------------------------------------------------------
p <- ggplot(plotdf_pw, aes(y = mn_group, x = Country, fill = mean)) +
  geom_point(aes(size = n_surveys), shape = 21, colour = "black") +
  scale_size_continuous(range = c(1, 6), name = "Number of surveys") +
  coord_flip() +
  scale_fill_gradient(low = "navy", high = "skyblue",
                      na.value = "grey90", name = "Scaled mean (latest)") +
  theme_minimal(base_size = 12) +
  labs(title = "VMNIS data availability:\nPregnant women",
       subtitle = "Most recent national survey\ncovering both urban & rural areas",
       x = "Country", y = "Micronutrient (collapsed)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

p_alt <- ggplot(plotdf_pw %>% filter(mn_group!="Haemoglobin") %>% droplevels(),
       aes(y = mn_group, x = Country, fill = mean)) +
  geom_point(aes(size = n_surveys), shape = 21, colour = "black") +
  scale_size_continuous(range = c(1, 6), name = "Number of surveys") +
  coord_flip() +
  scale_fill_gradient(low = "navy", high = "skyblue",
                      na.value = "grey90", name = "Scaled mean (latest survey)") +
  theme_minimal(base_size = 12) +
  labs(title = "VMNIS data availability:\nPregnant women",
       subtitle = "Most recent national survey\ncovering both urban & rural areas",
       x = "Country", y = "Micronutrient (collapsed)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

#save plot objects
ggsave(here("figures/VMNIS_heatmap_pw.png"), p, width = 10, height = 8)
ggsave(here("figures/VMNIS_heatmap_pw_alt.png"), p_alt, width = 10, height = 8)

#all populations
ggplot(survey_summary, aes(y = mn_group, x = Country, fill = mean)) +
  geom_point(aes(size = n_surveys), shape = 21, colour = "black") +
  scale_size_continuous(range = c(1, 6), name = "Number of surveys") +
  facet_wrap(~Population) +
  coord_flip() +
  scale_fill_gradient(low = "navy", high = "skyblue",
                      na.value = "grey90", name = "Scaled mean (latest)") +
  theme_minimal(base_size = 12) +
  labs(title = "VMNIS data availability – Pregnant women",
       subtitle = "Most recent national survey covering both urban & rural areas",
       x = "Country", y = "Micronutrient (collapsed)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")


