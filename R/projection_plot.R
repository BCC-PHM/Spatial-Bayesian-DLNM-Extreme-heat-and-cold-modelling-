setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Spatial Bayesian DNLM Extreme heat and cold")
# setwd("~/R_project/Extreme heat and cold")
library(readr)
library(tidyverse)
library(sf)
library(tmap)
library(readxl)

RCP26 = read_csv("output/RCP26_heat_cold_final.csv")
RCP45 = read_csv("output/RCP45_heat_cold_final.csv")
RCP60 = read_csv("output/RCP60_heat_cold_final.csv")
RCP85 = read_csv("output/RCP85_heat_cold_final.csv")


ward_pop = read_excel("data/external/population/sapewardstablefinal.xlsx", 
                      sheet = "Mid-2022 Ward 2022", skip = 3)


ward_map = read_sf("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp")


all_combined = rbind(RCP26,RCP45,RCP60,RCP85)
all_combined = ward_pop %>% 
  filter(`LAD 2022 Name` == "Birmingham") %>% 
  pivot_longer(cols = c(-`LAD 2022 Name`,-`LAD 2022 Name`,-`Ward 2022 Code`,-`Ward 2022 Name`,-Total,
                        -`LAD 2022 Code`),
               names_to = "age",
               values_to = "count") %>% 
  group_by(`Ward 2022 Code`,`Ward 2022 Name`) %>% 
  summarise(count = sum(count),
            .groups = "drop") %>% 
  rename(Ward_code = `Ward 2022 Code`,
         Ward_name = `Ward 2022 Name`) %>% 
  left_join(all_combined, by =c("Ward_code" = "Ward_Code")) %>% 
  mutate(heat_med = as.numeric(heat_med),
         heat_LL = as.numeric(heat_LL),
         heat_UL = as.numeric(heat_UL),
         cold_med = as.numeric(cold_med),
         cold_LL = as.numeric(cold_LL),
         cold_UL = as.numeric(cold_UL)) %>% 
  group_by(Ward_code, Ward_name,Decade,emission) %>%
  mutate(heat_med = heat_med/count*100000,
         heat_LL = heat_LL/count*100000,
         heat_UL = heat_UL/count*100000,
         cold_med = cold_med/count*100000,
         cold_LL = cold_LL/count*100000,
         cold_UL = cold_UL/count*100000,
         Decade = paste0(Decade,"s")
  ) %>% 
  ungroup()
  
  

#-------------------------------------------------------------
#baseline heat and cold from observed deaths

ward_EM_summary = read_rds("output/heat_and_cold_related_EM_ward.rds")


obs = ward_pop %>% 
  filter(`LAD 2022 Name` == "Birmingham") %>% 
  pivot_longer(cols = c(-`LAD 2022 Name`,-`LAD 2022 Name`,-`Ward 2022 Code`,-`Ward 2022 Name`,-Total,
                        -`LAD 2022 Code`),
               names_to = "age",
               values_to = "count") %>% 
  group_by(`Ward 2022 Code`,`Ward 2022 Name`) %>% 
  summarise(count = sum(count)) %>% 
  rename(Ward_code = `Ward 2022 Code`,
         Ward_name = `Ward 2022 Name`) %>% 
  left_join(ward_EM_summary, by ="Ward_code") %>%  #join em summary with pop estimate
  group_by(Ward_code, Ward_name, Ward_id) %>% 
  #standardisation
  mutate(heat_med = heat_med/count*100000,
         heat_LL = heat_LL/count*100000,
         heat_UL = heat_UL/count*100000,
         cold_med = cold_med/count*100000,
         cold_LL = cold_LL/count*100000,
         cold_UL = cold_UL/count*100000,
         Decade = "baseline"
  ) %>% 
  ungroup()


#grid to duplicate baseline for each scenario 
scenario_grid = expand.grid(emission = c("RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"),
            Ward_id = 1:69)


obs = obs %>% 
  left_join(scenario_grid, by = "Ward_id", relationship = "many-to-many") %>% 
  mutate(Decade = "Baseline")

#separate cold
obs_cold =obs %>% 
  select(Decade, Ward_code,Ward_name,cold_med,emission)
  

#separate heat
obs_heat =obs %>% 
  select(Decade, Ward_code,Ward_name,heat_med,emission)

#-------------------------------------------------------------
#cold

all_combined_cold = all_combined %>% 
  select(Decade, Ward_code,Ward_name,cold_med,emission) %>% 
  rbind(obs_cold) %>% 
  mutate(
    Decade = factor(Decade, levels= c("Baseline", "2030s", "2040s", "2050s", "2060s","2070s")),
    emission = factor(emission, levels = c("RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"))
  )

all_combined_cold = ward_map %>% 
  left_join(all_combined_cold, by = c("Ward_Code"="Ward_code"))

#-------------------------------------------------------------
#heat

all_combined_heat= all_combined %>% 
  select(Decade, Ward_code,Ward_name,heat_med,emission)%>% 
  rbind(obs_heat) %>% 
  mutate(
    Decade = factor(Decade, levels= c("Baseline", "2030s", "2040s", "2050s", "2060s","2070s")),
    emission = factor(emission, levels = c("RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"))
  )

all_combined_heat = ward_map %>% 
  left_join(all_combined_heat, by = c("Ward_Code"="Ward_code")) 

#-------------------------------------------------------------
tmap_mode("plot")

proj_cold = tm_shape(all_combined_cold)+
  tm_polygons(
    fill= "cold_med",
    fill.scale = tm_scale_continuous(values = "blues"),
    fill.legend = tm_legend("per 100,000", group_id = "bottom", frame = TRUE,bg.alpha=0)
  )+
  tm_facets_grid(rows = "emission", columns = "Decade")+
  tm_layout(
    # facet strip labels (row/column labels)
    panel.label.size = 1.4,
    panel.label.bg.color = "white",
    panel.label.bg.alpha = 0.85,
    
    # legend text
    legend.title.size = 1.2,
    legend.text.size  = 1.05,
    
    # overall font scaling (applies broadly)
    fontfamily = "sans",
    
    legend.outside = TRUE,
    legend.outside.position = "bottom",
    frame = FALSE
  )+
  tm_title("Posterior Median Projections of Annual Cold-Related Excess Mortality by Scenario and Decade") 

tmap_save(
  proj_cold,
  filename = "figs/cold_projection_posterior_median.png",
  width = 12, height = 7, units = "in",
  dpi = 600
)

proj_heat = tm_shape(all_combined_heat)+
  tm_polygons(
    fill= "heat_med",
    fill.scale = tm_scale_continuous(values = "reds"),
    fill.legend = tm_legend("per 100,000", group_id = "top", frame = FALSE,bg.alpha=0)
  )+
  tm_facets_grid(rows = "emission", columns = "Decade")+
  tm_facets_grid(rows = "emission", columns = "Decade")+
  tm_layout(
    # facet strip labels (row/column labels)
    panel.label.size = 1.4,
    panel.label.bg.color = "white",
    panel.label.bg.alpha = 0.85,
    
    # legend text
    legend.title.size = 1.2,
    legend.text.size  = 1.05,
    
    # overall font scaling (applies broadly)
    fontfamily = "sans",
    
    legend.outside = TRUE,
    legend.outside.position = "bottom",
    frame = FALSE
  )+
  tm_title("Posterior Median Projections of Annual Heat-Related Excess Mortality by Scenario and Decade") 

tmap_save(
  proj_heat,
  filename = "figs/heat_projection_posterior_median.png",
  width = 12, height = 7, units = "in",
  dpi = 600
)
