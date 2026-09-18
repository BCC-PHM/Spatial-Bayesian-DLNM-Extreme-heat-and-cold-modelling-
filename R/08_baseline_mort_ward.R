library(tidyverse)
library(readr)

#------------------------------------------------------
#load the data
#------------------------------------------------------
df_complete = read_rds("output/03_df_complete_for_inla.rds")
inla_model = read_rds("output/03_inla_model.rds")

#------------------------------------------------------
#median fitted value, averaged across ward,
#to obtain baseline deaths
#------------------------------------------------------
stopifnot(nrow(inla_model$summary.fitted.values) == nrow(df_complete))

baseline_mort = inla_model$summary.fitted.values[, "0.5quant", drop = FALSE] %>% 
  as.data.frame() %>% 
  rename(expected_deaths = `0.5quant`) %>% 
  cbind(df_complete) %>% 
  drop_na(expected_deaths, ward22cd, new_id) %>% 
  group_by(ward22cd, new_id) %>% 
  summarise(mean_baseline = mean(expected_deaths), .groups = "drop")


write_rds(baseline_mort, "output/08_baseline_mort_ward.rds")