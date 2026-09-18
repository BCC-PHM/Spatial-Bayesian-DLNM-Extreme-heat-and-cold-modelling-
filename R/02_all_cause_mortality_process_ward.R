setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Spatial Bayesian DNLM Extreme heat and cold")

library(sf)
library(dplyr)
library(readr)
library(tidyverse)

#=======================================================
#load the death extract and ward boundaries
#=======================================================
daily_death_20052026 = read_csv("data/raw/daily_death_20052026.csv")

ward_map = read_sf("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp")

#---completeness of the geography fields---------------
daily_death_20052026 %>% 
  count(is.na(`USUAL ADDRESS OF DECEASED`))

daily_death_20052026 %>% 
  count(is.na(`ELECTORAL WARD CODE OF USUAL RESIDENCE OF DECEASED`))

#---how many distinct wards per year-------------------
wardnum = daily_death_20052026 %>% 
  filter(`COUNTY DISTRICT CODE OF USUAL RESIDENCE OF DECEASED` == "E08000025") %>% 
  distinct(Year, `ELECTORAL WARD CODE OF USUAL RESIDENCE OF DECEASED`) %>% 
  count(Year)

#=======================================================
#check death ward codes against the boundary file
#=======================================================
#i guess Mo has remapped it to 69 wards
codes_death = daily_death_20052026 %>% 
  filter(`COUNTY DISTRICT CODE OF USUAL RESIDENCE OF DECEASED` == "E08000025") %>% 
  distinct(ward_code = `ELECTORAL WARD CODE OF USUAL RESIDENCE OF DECEASED`) %>% 
  pull(ward_code)

codes_map = ward_map %>% 
  st_drop_geometry() %>% 
  pull(Ward_Code)

setdiff(codes_death, codes_map)
setdiff(codes_map, codes_death)

#=======================================================
#daily ward counts, covid underlying cause removed
#=======================================================
#---all ages-------------------------------------------
all_cause_mortality2005_2025 = daily_death_20052026 %>% 
  filter(`ELECTORAL WARD CODE OF USUAL RESIDENCE OF DECEASED` %in% ward_map$Ward_Code) %>% 
  filter(Year %in% c(2005:2025)) %>% 
  filter(!`UNDERLYING CAUSE OF DEATH CODE` %in% c("U071", "U072")) %>% 
  mutate(dod = as.Date(sprintf("%08d", `DATE OF DEATH OF DECEASED`), format = "%d%m%Y")) %>% 
  rename(ward22cd = `ELECTORAL WARD CODE OF USUAL RESIDENCE OF DECEASED`) %>% 
  group_by(dod, ward22cd) %>% 
  summarise(deaths = n(), .groups = "drop")

#---by age group---------------------------------------
all_cause_mortality2005_2025_age_group = daily_death_20052026 %>% 
  filter(`ELECTORAL WARD CODE OF USUAL RESIDENCE OF DECEASED` %in% ward_map$Ward_Code) %>% 
  filter(Year %in% c(2005:2025)) %>% 
  filter(!`UNDERLYING CAUSE OF DEATH CODE` %in% c("U071", "U072")) %>% 
  mutate(dod       = as.Date(sprintf("%08d", `DATE OF DEATH OF DECEASED`), format = "%d%m%Y"),
         age_group = case_when(AGE %in% c(0:64)         ~ "0-64",
                               AGE %in% c(65:74)        ~ "65-74",
                               AGE %in% c(75:84)        ~ "75-84",
                               AGE %in% c(85:max(AGE))  ~ "85+",
                               TRUE                     ~ NA_character_)) %>% 
  rename(ward22cd = `ELECTORAL WARD CODE OF USUAL RESIDENCE OF DECEASED`) %>% 
  group_by(dod, age_group, ward22cd) %>% 
  summarise(deaths = n(), .groups = "drop")

#=======================================================
#ward daily mean air temperature
#=======================================================
air_temp = readRDS("data/processed/01_tmean_ward_daily.rds")

#=======================================================
#build the full ward x date panel
#=======================================================
dates = air_temp %>% 
  distinct(date)

full_panel_no_agegroup = dates %>% 
  crossing(ward22cd = ward_map$Ward_Code)

#---join deaths and air temp onto the panel------------
model_data_no_agegroup = full_panel_no_agegroup %>% 
  left_join(air_temp, by = c("date", "ward22cd" = "Ward_Code")) %>% 
  left_join(all_cause_mortality2005_2025, by = c("date" = "dod", "ward22cd")) %>% 
  mutate(deaths = ifelse(is.na(deaths), 0, deaths))

#---one panel per age group----------------------------
model_data_0_64 = full_panel_no_agegroup %>% 
  left_join(air_temp, by = c("date", "ward22cd" = "Ward_Code")) %>% 
  left_join(all_cause_mortality2005_2025_age_group %>% filter(age_group == "0-64"), by = c("date" = "dod", "ward22cd")) %>% 
  mutate(deaths    = ifelse(is.na(deaths), 0, deaths),
         age_group = ifelse(is.na(age_group), "0-64", age_group))

model_data_65_74 = full_panel_no_agegroup %>% 
  left_join(air_temp, by = c("date", "ward22cd" = "Ward_Code")) %>% 
  left_join(all_cause_mortality2005_2025_age_group %>% filter(age_group == "65-74"), by = c("date" = "dod", "ward22cd")) %>% 
  mutate(deaths    = ifelse(is.na(deaths), 0, deaths),
         age_group = ifelse(is.na(age_group), "65-74", age_group))

model_data_75_84 = full_panel_no_agegroup %>% 
  left_join(air_temp, by = c("date", "ward22cd" = "Ward_Code")) %>% 
  left_join(all_cause_mortality2005_2025_age_group %>% filter(age_group == "75-84"), by = c("date" = "dod", "ward22cd")) %>% 
  mutate(deaths    = ifelse(is.na(deaths), 0, deaths),
         age_group = ifelse(is.na(age_group), "75-84", age_group))

model_data_85 = full_panel_no_agegroup %>% 
  left_join(air_temp, by = c("date", "ward22cd" = "Ward_Code")) %>% 
  left_join(all_cause_mortality2005_2025_age_group %>% filter(age_group == "85+"), by = c("date" = "dod", "ward22cd")) %>% 
  mutate(deaths    = ifelse(is.na(deaths), 0, deaths),
         age_group = ifelse(is.na(age_group), "85+", age_group))

#---bundle and save------------------------------------
model_data_no_agegroup = list(model_data_no_agegroup,
                              model_data_0_64,
                              model_data_65_74,
                              model_data_75_84,
                              model_data_85)

write_rds(model_data_no_agegroup, "data/processed/02_model_data_ward22.rds")






