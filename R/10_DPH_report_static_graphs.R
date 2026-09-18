#==================================================
#packages
#==================================================
library(here)
library(readr)
library(readxl)
library(tidyverse)
library(sf)
library(tmap)

#==================================================
#load rr / mmt outputs and ward boundaries
#==================================================
real_plot_df = readRDS(here("output", "04_RR_MMT_plot_data.rds"))
RR_prob_df = readRDS(here("output", "04_RR_exceedance_prob.rds"))
ward_map = read_sf(here("data", "external", "boundaries", "boundaries-wards-2022-birmingham", "boundaries-wards-2022-birmingham.shp"))

#==================================================
#full-distribution excess mortality
#==================================================
heat_annual_post = read_rds(here("output", "05_heat_annual_post_EM_ward.rds"))
cold_annual_post = read_rds(here("output", "05_cold_annual_post_EM_ward.rds"))
ward_EM_summary = read_rds(here("output", "05_heat_and_cold_related_EM_ward.rds"))
ward_pop = read_excel(here("data", "external", "population", "sapewardstablefinal.xlsx"), sheet = "Mid-2022 Ward 2022", skip = 3)

#---------------------------------------------
#birmingham cold draws
#---------------------------------------------
bham_excess_cold_list = vector("list", 1000)

for (i in 1:1000){
  #2. use a temporary name for the math (e.g., current_sum)
  current_sum = sum(sapply(cold_annual_post, `[`, i))
  #3. save it into the storage list
  bham_excess_cold_list [[i]] = data.frame(draws = i, value = current_sum)
}

bham_excess_cold_df = bind_rows(bham_excess_cold_list)

#---------------------------------------------
#birmingham heat draws
#---------------------------------------------
bham_excess_heat_list = vector("list", 1000)

for (i in 1:1000){
  #2. use a temporary name for the math (e.g., current_sum)
  current_sum = sum(sapply(heat_annual_post, `[`, i))
  #3. save it into the storage list
  bham_excess_heat_list[[i]] = data.frame(draws = i, value = current_sum)
}

bham_excess_heat_df = bind_rows(bham_excess_heat_list)

#---------------------------------------------
#birmingham total population
#---------------------------------------------
bham_pop = ward_pop %>% 
  filter(`LAD 2022 Name` == "Birmingham") %>% 
  pivot_longer(
    cols      = c(-`LAD 2022 Name`, -`LAD 2022 Name`, -`Ward 2022 Code`, -`Ward 2022 Name`, -Total, -`LAD 2022 Code`),
    names_to  = "age",
    values_to = "count"
  ) %>% 
  group_by(`Ward 2022 Code`, `Ward 2022 Name`) %>% 
  summarise(count = sum(count)) %>% 
  rename(Ward_code = `Ward 2022 Code`, Ward_name = `Ward 2022 Name`) %>% 
  pull(count) %>% 
  sum()

#---------------------------------------------
#birmingham median rates and cris
#---------------------------------------------
bham_mean_cold_EM = median(bham_excess_cold_df$value)/bham_pop*100000
bham_mean_heat_EM = median(bham_excess_heat_df$value)/bham_pop*100000
bham_mean_cold_EM
bham_mean_heat_EM

quantile(bham_excess_cold_df$value, 0.975)
median(bham_excess_cold_df$value)
quantile(bham_excess_cold_df$value, 0.025)

quantile(bham_excess_heat_df$value, 0.975)
median(bham_excess_heat_df$value)
quantile(bham_excess_heat_df$value, 0.025)

#==================================================
#extreme temperature excess mortality
#==================================================
X_heat_annual_post = read_rds(here("output", "05_X_heat_annual_post_EM_ward.rds"))
X_cold_annual_post = read_rds(here("output", "05_X_cold_annual_post_EM_ward.rds"))

#---------------------------------------------
#birmingham extreme cold draws
#---------------------------------------------
bham_excess_Xcold_list = vector("list", 1000)

for (i in 1:1000){
  #2. use a temporary name for the math (e.g., current_sum)
  current_sum = sum(sapply(X_cold_annual_post, `[`, i))
  #3. save it into the storage list
  bham_excess_Xcold_list [[i]] = data.frame(draws = i, value = current_sum)
}

bham_excess_Xcold_df = bind_rows(bham_excess_Xcold_list)

#---------------------------------------------
#birmingham extreme heat draws
#---------------------------------------------
bham_excess_Xheat_list = vector("list", 1000)

for (i in 1:1000){
  #2. use a temporary name for the math (e.g., current_sum)
  current_sum = sum(sapply(X_heat_annual_post, `[`, i))
  #3. save it into the storage list
  bham_excess_Xheat_list [[i]] = data.frame(draws = i, value = current_sum)
}

bham_excess_Xheat_df = bind_rows(bham_excess_Xheat_list)

#---------------------------------------------
#birmingham total population
#---------------------------------------------
bham_pop = ward_pop %>% 
  filter(`LAD 2022 Name` == "Birmingham") %>% 
  pivot_longer(
    cols      = c(-`LAD 2022 Name`, -`LAD 2022 Name`, -`Ward 2022 Code`, -`Ward 2022 Name`, -Total, -`LAD 2022 Code`),
    names_to  = "age",
    values_to = "count"
  ) %>% 
  group_by(`Ward 2022 Code`, `Ward 2022 Name`) %>% 
  summarise(count = sum(count)) %>% 
  rename(Ward_code = `Ward 2022 Code`, Ward_name = `Ward 2022 Name`) %>% 
  pull(count) %>% 
  sum()

#---------------------------------------------
#birmingham extreme median rates and cris
#---------------------------------------------
bham_mean_Xcold_EM = median(bham_excess_Xcold_df$value)/bham_pop*100000
bham_mean_Xheat_EM = median(bham_excess_Xheat_df$value)/bham_pop*100000
bham_mean_Xcold_EM
bham_mean_Xheat_EM

quantile(bham_excess_Xcold_df$value, 0.975)
median(bham_excess_Xcold_df$value)
quantile(bham_excess_Xcold_df$value, 0.025)

quantile(bham_excess_Xheat_df$value, 0.975)
median(bham_excess_Xheat_df$value)
quantile(bham_excess_Xheat_df$value, 0.025)

#==================================================
#figure 6
#==================================================
ward_X_EM_summary = read_rds(here("output", "05_X_heat_and_cold_related_EM_ward.rds"))
ward_pop = read_excel("data/external/population/sapewardstablefinal.xlsx", sheet = "Mid-2022 Ward 2022", skip = 3)
exceed_prob_EM_gt_bham_X_heatcold = read_rds(here("output", "05_exceed_prob_EM_gt_bham_X_heatcold.rds"))

#---------------------------------------------
#evidence bands for the exceedance maps
#---------------------------------------------
evidence_levels = c("Very strong evidence (≥0.95)", "Strong evidence (0.90–0.95)", "Some evidence (0.80–0.90)", "No evidence (<0.80)")

em_exceedance_Xplot_df = ward_map %>% 
  left_join(exceed_prob_EM_gt_bham_X_heatcold, by = c("Ward_Code" = "ward22cd"))

em_exceedance_Xplot_df = em_exceedance_Xplot_df %>% 
  mutate(evidence_Xcold = case_when(p_Xcold_gt_bham >= 0.95 ~ "Very strong evidence (≥0.95)",
                                    p_Xcold_gt_bham >= 0.90 ~ "Strong evidence (0.90–0.95)",
                                    p_Xcold_gt_bham >= 0.80 ~ "Some evidence (0.80–0.90)",
                                    TRUE                    ~ "No evidence (<0.80)"),
         evidence_Xheat = case_when(p_Xheat_gt_bham >= 0.95 ~ "Very strong evidence (≥0.95)",
                                    p_Xheat_gt_bham >= 0.90 ~ "Strong evidence (0.90–0.95)",
                                    p_Xheat_gt_bham >= 0.80 ~ "Some evidence (0.80–0.90)",
                                    TRUE                    ~ "No evidence (<0.80)"),
         evidence_Xcold = factor(evidence_Xcold, levels = evidence_levels),
         evidence_Xheat = factor(evidence_Xheat, levels = evidence_levels))

#---------------------------------------------
#ward rates per 100,000 and tooltips
#---------------------------------------------
ward_pop = ward_pop %>% 
  filter(`LAD 2022 Name` == "Birmingham") %>% 
  pivot_longer(
    cols      = c(-`LAD 2022 Name`, -`LAD 2022 Name`, -`Ward 2022 Code`, -`Ward 2022 Name`, -Total, -`LAD 2022 Code`),
    names_to  = "age",
    values_to = "count"
  ) %>% 
  group_by(`Ward 2022 Code`, `Ward 2022 Name`) %>% 
  summarise(count = sum(count), .groups = "drop") %>% 
  rename(Ward_code = `Ward 2022 Code`, Ward_name = `Ward 2022 Name`) %>% 
  left_join(ward_X_EM_summary, by = "Ward_code") %>%  #join em summary with pop estimate
  group_by(Ward_code, Ward_name, Ward_id) %>% 
  #standardisation
  mutate(heat_med = heat_med/count*100000,
         heat_LL  = heat_LL/count*100000,
         heat_UL  = heat_UL/count*100000,
         cold_med = cold_med/count*100000,
         cold_LL  = cold_LL/count*100000,
         cold_UL  = cold_UL/count*100000,
         #cold tootip
         tooltip_cold = paste0("<B>Ward Name:</B> ", Ward_name, "\n",
                               "<B>Ward Code:</B> ", Ward_code, "\n",
                               "<B>Median EM:</B> ", round(cold_med, 2), "\n",
                               "<B>Upper95CrI:</B> ", round(cold_UL, 2), "\n",
                               "<B>Lower95CrI:</B> ", round(cold_LL, 2)),
         #heat tooltip
         tooltip_heat = paste0("<B>Ward Name:</B> ", Ward_name, "\n",
                               "<B>Ward Code:</B> ", Ward_code, "\n",
                               "<B>Median EM:</B> ", round(heat_med, 2), "\n",
                               "<B>Upper95CrI:</B> ", round(heat_UL, 2), "\n",
                               "<B>Lower95CrI:</B> ", round(heat_LL, 2))
  )

tmap_options(component.autoscale = FALSE)

#---------------------------------------------
#extreme cold rate map
#---------------------------------------------
figure_EM_Xcold = tm_shape(ward_map %>% 
                             left_join(ward_pop, by = c("Ward_Code" = "Ward_code")))+
  tm_polygons(
    fill        = "cold_med",
    col         = "white",
    fill.scale  = tm_scale_continuous(values = "brewer.blues"),
    fill.legend = tm_legend("per 100,000", group_id = "top", frame = FALSE)
  )+
  tm_layout(
    frame             = FALSE,
    legend.position   = c("left", "top"),
    legend.frame      = FALSE,
    legend.text.size  = 0.6,
    legend.title.size = 0.8,
    inner.margins     = c(0.07, 0, 0.01, 0)
  )+
  tm_title("Estimated non-age-standardised annual median \nexcess mortality rate from extreme cold by ward", size = 1)+
  tm_compass(type = "8star", size = 4, position = c("RIGHT", "bottom"), color.light = "white")+
  tm_credits(
    text     = paste0("Birmingham median: ", round(bham_mean_Xcold_EM, 2)),
    position = c(0, 1.02),
    col      = "grey40",
    size     = 0.9
  )+
  tm_credits(
    text     = paste("Contains OS data \u00A9 Crown copyright and database right", format(Sys.Date(), "%Y"), ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."),
    position = c("LEFT", "BOTTOM")
  )

#---------------------------------------------
#extreme heat rate map
#---------------------------------------------
figure_EM_Xheat = tm_shape(ward_map %>% 
                             left_join(ward_pop, by = c("Ward_Code" = "Ward_code")))+
  tm_polygons(
    fill        = "heat_med",
    col         = "white",
    fill.scale  = tm_scale_continuous(values = "brewer.reds", ticks = c(2, 4, 6)),
    fill.legend = tm_legend("per 100,000", group_id = "top", frame = FALSE)
  )+
  tm_layout(
    frame             = FALSE,
    legend.position   = c("left", "top"),
    legend.frame      = FALSE,
    egend.text.size   = 0.6,
    legend.title.size = 0.8,
    inner.margins     = c(0.07, 0, 0.01, 0)
  )+
  tm_title("Estimated non-age-standardised annual median \nexcess mortality rate from extreme heat by ward", size = 1)+
  tm_compass(type = "8star", size = 4, position = c("RIGHT", "bottom"), color.light = "white")+
  tm_credits(
    text     = paste0("Birmingham median: ", round(bham_mean_Xheat_EM, 2)),
    position = c(0, 1.02),
    col      = "grey40",
    size     = 0.9
  )+
  tm_credits(
    text     = paste("Contains OS data \u00A9 Crown copyright and database right", format(Sys.Date(), "%Y"), ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."),
    position = c("LEFT", "BOTTOM")
  )

#---------------------------------------------
#extreme cold evidence map
#---------------------------------------------
figure_EM_prob_Xcold = tm_shape(em_exceedance_Xplot_df)+
  tm_polygons(
    fill        = "evidence_Xcold",
    col         = "white",
    fill.scale  = tm_scale_categorical(values = c("#08519c", "#4292c6", "#9ecae1", "grey85")),
    fill.legend = tm_legend("Probability > Birmingham median", group_id = "top", frame = FALSE)
  )+
  tm_layout(
    frame             = FALSE,
    legend.position   = c("left", "top"),
    legend.frame      = FALSE,
    legend.text.size  = 0.6,
    legend.title.size = 0.8,
    inner.margins     = c(0.07, 0, 0.01, 0)
  )+
  tm_title("Evidence that excess mortality rate during extreme cold \nis above the Birmingham median, by ward", size = 1)+
  tm_compass(type = "8star", size = 4, position = c("RIGHT", "bottom"), color.light = "white")+
  tm_credits(
    text     = paste("Contains OS data \u00A9 Crown copyright and database right", format(Sys.Date(), "%Y"), ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."),
    position = c("LEFT", "BOTTOM")
  )

#---------------------------------------------
#extreme heat evidence map
#---------------------------------------------
figure_EM_prob_Xheat = tm_shape(em_exceedance_Xplot_df)+
  tm_polygons(
    fill        = "evidence_Xheat",
    col         = "white",
    fill.scale  = tm_scale_categorical(values = c("#a50f15", "#ef3b2c", "#fc9272", "grey85")),
    fill.legend = tm_legend("Probability > Birmingham median", group_id = "top", frame = FALSE)
  )+
  tm_layout(
    frame             = FALSE,
    legend.position   = c("left", "top"),
    legend.frame      = FALSE,
    legend.text.size  = 0.6,
    legend.title.size = 0.8,
    inner.margins     = c(0.07, 0, 0.01, 0)
  )+
  tm_title("Evidence that excess mortality rate during extreme heat \nis above the Birmingham median, by ward", size = 1)+
  tm_compass(type = "8star", size = 4, position = c("RIGHT", "bottom"), color.light = "white")+
  tm_credits(
    text     = paste("Contains OS data \u00A9 Crown copyright and database right", format(Sys.Date(), "%Y"), ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."),
    position = c("LEFT", "BOTTOM")
  )

#---------------------------------------------
#combine and save
#---------------------------------------------
tmap_arrange(figure_EM_Xcold, figure_EM_prob_Xcold)

combined_cold = tmap_arrange(figure_EM_Xcold, figure_EM_prob_Xcold, ncol = 2)
combined_heat = tmap_arrange(figure_EM_Xheat, figure_EM_prob_Xheat, ncol = 2)

tmap_save(combined_cold, filename = "figs/DPH_2026/10_DPH_figure_EM_Xcold_combined.png", width = 10, height = 6, units = "in", dpi = 600)
tmap_save(combined_heat, filename = "figs/DPH_2026/10_DPH_figure_EM_Xheat_combined.png", width = 10, height = 6, units = "in", dpi = 600)

