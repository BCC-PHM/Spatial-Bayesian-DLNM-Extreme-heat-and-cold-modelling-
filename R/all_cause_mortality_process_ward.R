
library(tidyverse)
library(readxl)
library(sf)

all_cause_mortality_2014_2024 = read_excel("data/raw/all_cause_mortality_2014_2024.xlsx")


ward_map = read_sf("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp")

all_cause_mortality_2014_2024 %>% 
  distinct(WARD_OF_RESIDENCE_CODE) %>% 
  filter(WARD_OF_RESIDENCE_CODE %in% ward_map$Ward_Code)


#####################################################
#death data that not yet transformed to lsoa11
all_cause_mortality_2014_2024 = all_cause_mortality_2014_2024 %>%
  filter(WARD_OF_RESIDENCE_CODE %in% ward_map$Ward_Code) %>% 
  #remove covid deaths
  filter(! S_UNDERLYING_COD_ICD10 %in% c("U071", "U072")) %>% 
  rename(ward22cd =WARD_OF_RESIDENCE_CODE) %>% 
  group_by(REG_DATE, ward22cd) %>% 
  summarise( deaths = n(),
             .groups = "drop")


######################################################

#aggregate air temp 
air_temp =  readRDS("data/processed/tmean_ward_daily.rds")



######################################################
#join daily by age group to air temp

dates = air_temp %>%
  distinct(date)


full_panel_no_agegroup = dates %>% 
  crossing(
    ward22cd = ward_map$Ward_Code
  )


#join the deaths and air temp to the full panel.

model_data_no_agegroup = full_panel_no_agegroup %>% 
  left_join(air_temp, by = c("date", "ward22cd" = "Ward_Code")) %>% 
  left_join(all_cause_mortality_2014_2024, by = c("date" = "REG_DATE", "ward22cd")) %>% 
  mutate(deaths = ifelse(is.na(deaths), 0, deaths))





write_rds(model_data_no_agegroup , "data/processed/model_data_ward22.rds")








