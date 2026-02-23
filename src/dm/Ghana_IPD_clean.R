

library(tidyverse)
library(haven)
library(here)
source(paste0(here::here(),"/src/0-functions.R"))


# Load the datasets

dchild <- read_dta(file=here("data/IPD/Ghana/GMS_ChildrenAllData_DavisBerkeley.dta"))
dwomen <- read_dta(file=here("data/IPD/Ghana/GMS_WomenAllData_DavisBerkeley.dta"))

head(dchild)
head(dwomen)


labels_child <- makeVlist(dchild)
labels_women <- makeVlist(dwomen)
write.csv(labels_child, file=here("data/IPD/Ghana/labels_child.csv"), row.names = FALSE)
write.csv(labels_women, file=here("data/IPD/Ghana/labels_women.csv"), row.names = FALSE)
