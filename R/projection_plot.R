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
BRCP26 = read.csv("output/Baseline_RCP26_heat_cold_final.csv")
BRCP45 = read.csv("output/Baseline_RCP45_heat_cold_final.csv")
BRCP60 = read.csv("output/Baseline_RCP60_heat_cold_final.csv")
BRCP85 = read.csv("output/Baseline_RCP85_heat_cold_final.csv")

ward_pop = read_excel("data/external/population/sapewardstablefinal.xlsx", 
                      sheet = "Mid-2022 Ward 2022", skip = 3)


ward_map = read_sf("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp")


all_combined = rbind(RCP26[,-1],BRCP26[,-1],
                     RCP45[,-1],BRCP45[,-1],
                     RCP60[,-1],BRCP60[,-1],
                     RCP85[,-1],BRCP85[,-1]
)
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
         Decade = ifelse(Decade != "Baseline",paste0(Decade,"s"), Decade)
  ) %>% 
  ungroup()




#-------------------------------------------------------------
#baseline heat and cold from observed deaths

ward_EM_summary = read_rds("output/heat_and_cold_related_EM_ward.rds")



#-------------------------------------------------------------
#cold

all_combined_cold = all_combined %>% 
  select(Decade, Ward_code,Ward_name,cold_med,emission) %>% 
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
  mutate(
    Decade = factor(Decade, levels= c("Baseline", "2030s", "2040s", "2050s", "2060s","2070s")),
    emission = factor(emission, levels = c("RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"))
  )

all_combined_heat = ward_map %>% 
  left_join(all_combined_heat, by = c("Ward_Code"="Ward_code")) 

#-------------------------------------------------------------
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
    legend.title.size = 1,
    legend.text.size  = 0.8,
    
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
    legend.title.size = 1,
    legend.text.size  = 0.8,
    
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


################################################################################

RCP26AF = read_csv("output/AF_RCP26_heat_cold_final.csv")
RCP45AF = read_csv("output/AF_RCP45_heat_cold_final.csv")
RCP60AF = read_csv("output/AF_RCP60_heat_cold_final.csv")
RCP85AF = read_csv("output/AF_RCP85_heat_cold_final.csv")
BRCP26AF = read.csv("output/Baseline_AF_RCP26_heat_cold_final.csv")
BRCP45AF = read.csv("output/Baseline_AF_RCP45_heat_cold_final.csv")
BRCP60AF = read.csv("output/Baseline_AF_RCP60_heat_cold_final.csv")
BRCP85AF = read.csv("output/Baseline_AF_RCP85_heat_cold_final.csv")






all_combined_AF = rbind(RCP26AF[,-1],BRCP26AF[,-1],
                     RCP45AF[,-1],BRCP45AF[,-1],
                     RCP60AF[,-1],BRCP60AF[,-1],
                     RCP85AF[,-1],BRCP85AF[,-1])%>% 
  mutate(
    Decade = ifelse(Decade != "Baseline",paste0(Decade,"s"), Decade),
    Decade = factor(Decade, levels= c("Baseline", "2030s", "2040s", "2050s", "2060s","2070s")),
    emission = factor(emission, levels = c("RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"))
    )
  




all_combined_AF = all_combined_AF %>%
  mutate(Ward_Name = coalesce(Ward_Name.x, Ward_Name.y)) %>% 
  select(-Ward_Name.x, -Ward_Name.y) %>%
  rename_with(~ gsub("coldAF_", "cold_", .x)) %>%
  rename_with(~ gsub("heatAF_", "heat_", .x)) %>%
  pivot_longer(
    cols = c(cold_med, cold_LL, cold_UL, heat_med, heat_LL, heat_UL),
    names_to = c("type", "stat"),
    names_sep = "_",
    values_to = "AF"
  ) %>%
  pivot_wider(names_from = stat, values_from = AF, names_prefix = "AF_")
  
  
lvl <- c("Baseline","2030s","2040s","2050s","2060s","2070s")

plot_df <- all_combined_AF %>%
  filter(ID == 69, emission == "RCP4.5") %>%
  mutate(
    Decade = factor(Decade, levels = lvl),
    x_dec  = as.numeric(Decade)   # 1..6 (Baseline=1, 2030s=2, ...)
  )





ggplot(plot_df, aes(x = x_dec, y = AF_med, group = type)) +
  geom_ribbon(aes(ymin = AF_LL, ymax = AF_UL, fill = type),
              alpha = 0.2, colour = NA) +
  geom_line(aes(colour = type), linewidth = 1) +
  geom_point(aes(colour = type), size = 2) +
  scale_x_continuous(breaks = seq_along(lvl), labels = lvl) +
  labs(x = "Decade", y = "AF", colour = NULL, fill = NULL)+
  coord_cartesian(ylim = c(0, 20)) 




















