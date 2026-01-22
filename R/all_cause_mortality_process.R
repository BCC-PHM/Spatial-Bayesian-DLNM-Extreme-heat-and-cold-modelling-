setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Extreme heat and cold")


library(tidyverse)
library(readxl)
library(sf)

all_cause_mortality_2014_2024 = read_excel("data/raw/all_cause_mortality_2014_2024.xlsx")


#check duplicates
# duplicated_ID = all_cause_mortality_2014_2024 %>%
#   count(PatientId) %>%
#   filter(n>1) %>%
#   pull(PatientId)


#####################################################
LSOA_21_map = read_sf(
  "data/external/boundaries/boundaries-lsoa-2021-birmingham/boundaries-lsoa-2021-birmingham.shp")

LSOA_11_map = read_sf(
  "data/external/boundaries/LSOA_2011_map/Birmingham LSOA shape.shp")

#extra the extra 20 codes that lsoa21 has and lsoa11 does not
extra_code = LSOA_21_map$LSOA21CD[
  !LSOA_21_map$LSOA21CD %in% LSOA_11_map$LSOA11CD
]

#Confirm whether any extra codes appear
any(all_cause_mortality_2014_2024$LSOA_OF_RESIDENCE_CODE %in% extra_code)

#it reutns false, it seems like the death register never used lsoa21 till now 
#it means that we will have to apply the transformation of weigtings to all the lsoa11 to lsoa21 
#####################################################
#weightings data
lsoa11_to_lsoa21_weights = read_csv("data/processed/lsoa11_to_lsoa21_weights.csv")

#####################################################
#death data that not yet transformed to lsoa11
death_data_lsoa11 = all_cause_mortality_2014_2024 %>%
  filter(LSOA_OF_RESIDENCE_CODE %in% LSOA_11_map$LSOA11CD) %>% 
  #remove covid deaths
  filter(! S_UNDERLYING_COD_ICD10 %in% c("U071", "U072")) %>% 
  rename(lsoa11 =LSOA_OF_RESIDENCE_CODE ) %>% 
  mutate(age_group = case_when(
    DEC_AGEC >=0 & DEC_AGEC <=64 ~ "0-64",
    DEC_AGEC >=65 & DEC_AGEC <=74 ~"65-74",
    DEC_AGEC >=75 & DEC_AGEC <=84 ~"75-84",
    DEC_AGEC >=85 ~"85+",
    TRUE ~ NA_character_)) %>% 
  group_by(REG_DATE, age_group, lsoa11) %>% 
  summarise( deaths = n(),
             .groups = "drop")
  
  
############################################### 


options(scipen = 999)

death_data = all_cause_mortality_2014_2024 %>% 
  filter(LSOA_OF_RESIDENCE_CODE %in% LSOA_11_map$LSOA11CD) %>% 
  mutate(REG_DATE = as.Date(REG_DATE)) %>% 
  left_join(lsoa11_to_lsoa21_weights, by = c("LSOA_OF_RESIDENCE_CODE" = "lsoa11"),  relationship = "many-to-many") %>% 
  mutate(age_group = case_when(
    DEC_AGEC >=0 & DEC_AGEC <=64 ~ "0-64",
    DEC_AGEC >=65 & DEC_AGEC <=74 ~"65-74",
    DEC_AGEC >=75 & DEC_AGEC <=84 ~"75-84",
    DEC_AGEC >=85 ~"85+",
    TRUE ~ NA_character_)) %>% 
  group_by(REG_DATE, age_group, lsoa21) %>% 
  summarise( deaths = sum(weight),
             .groups = "drop")

#####################################################
#make one without age group
deathdata_no_agegroup = all_cause_mortality_2014_2024 %>% 
  filter(LSOA_OF_RESIDENCE_CODE %in% LSOA_11_map$LSOA11CD) %>% 
  filter(! S_UNDERLYING_COD_ICD10 %in% c("U071", "U072")) %>% 
  rename(lsoa11 =LSOA_OF_RESIDENCE_CODE ) %>% 
  mutate(REG_DATE = as.Date(REG_DATE)) %>% 
  group_by(REG_DATE, lsoa11) %>% 
  summarise( deaths = n(),
             .groups = "drop")
  
######################################################

#aggregate air temp 
air_temp =  readRDS("data/processed/tmean_lsoa11_daily.rds")

#####################################################

#join daily by age group to air temp

age_levels = c("0-64", "65-74", "75-84", "85+")

dates = air_temp %>%
  distinct(date)

lsoa_list = LSOA_11_map %>%
  st_drop_geometry() %>%
  select(lsoa11 = LSOA11CD)


full_panel = dates %>%
  crossing(
    age_group = age_levels,
    lsoa11 = lsoa_list$lsoa11
  )


#join the deaths and air temp to the full panel.

model_data = full_panel %>% 
  left_join(air_temp, by = c("date", "lsoa11" = "LSOA11CD")) %>% 
  left_join(death_data_lsoa11, by = c("date" = "REG_DATE", "age_group", "lsoa11")) %>% 
  mutate(
    deaths = if_else(is.na(deaths), 0, deaths)
  )
#####################################################
#make one without age group

full_panel_no_agegroup = dates %>% 
  crossing(
    lsoa11 = lsoa_list$lsoa11
  )


model_data_no_agegroup = full_panel_no_agegroup %>% 
  left_join(air_temp, by = c("date", "lsoa11" = "LSOA11CD")) %>% 
  left_join(deathdata_no_agegroup, by = c("date" = "REG_DATE", "lsoa11")) %>% 
  mutate(deaths = ifelse(is.na(deaths), 0, deaths))



write_rds(model_data, "data/processed/model_data_lsoa11.rds")

# write_rds(model_data, "data/processed/model_data.rds")

write_rds(model_data_no_agegroup, "data/processed/model_data_no_agegroup.rds")



##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################





# #fitting with a GAM cubic spline 
# ggplot(combined_daily, aes(x=date, y= count))+
#   geom_point(colour = "grey40")+
#   geom_smooth(
#     method = "gam",
#     method.args = list(family = quasipoisson),   #quasi poission allow variance to be larger than mean 
#     formula = y ~ s(x, bs = "cr", k = 21),
#     linetype = 2,
#     colour = "#962558",
#     fill = "#f391bd"
#   )+
#   scale_x_date(date_labels = "%e %b %y", 
#                limits = as.Date(c("2022-01-01", "2024-12-31")),
#                breaks = seq(
#                  as.Date("2022-01-01"),
#                  as.Date("2025-01-01"),
#                  by = "6 months"
#                ))+
#   ggtitle("Cubic spline model")+
#   labs(
#     x="Date",
#     y= " N deaths"
#   )+
#   theme_bw()+
#   theme(
#     plot.title = element_text(hjust = 0.5)
#   )


# A GLM method by deciding the number of knots yourself and fitting the exact number of knots   
# spl = bs(death_count_daily$REG_DATE, degree = 3,df=10)
#   
#   
# model =  glm(count~spl, death_count_daily, family = quasipoisson())
#   
# pred = predict(model, type="response")
#   
# plot(death_count_daily$REG_DATE,death_count_daily$count,pch=19,cex=0.2,col=grey(0.6),
#        main="Flexible cubic spline model",ylab="Daily number of deaths",
#        xlab="Date")
# 
# lines(death_count_daily$REG_DATE,pred,lwd=2,col=2)



