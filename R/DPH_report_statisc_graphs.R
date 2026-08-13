setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Spatial Bayesian DNLM Extreme heat and cold")

library(here)
library(readr)
library(tmap)

real_plot_df = readRDS(here("output", "RR_MMT_plot_data.rds"))
RR_prob_df = readRDS(here("output", "RR_exceedance_prob.rds"))
ward_map = read_sf(here(
  "data", "external", "boundaries",
  "boundaries-wards-2022-birmingham",
  "boundaries-wards-2022-birmingham.shp"
))

###########################################################################
#figure 6
ward_X_EM_summary = read_rds(here("output", "X_heat_and_cold_related_EM_ward.rds"))



ward_pop = read_excel(here("data", "external", "population", "sapewardstablefinal.xlsx"),
                      sheet = "Mid-2022 Ward 2022",
                      skip  = 3
)

exceed_prob_EM_gt_bham_X_heatcold = read_rds(here("output","exceed_prob_EM_gt_bham_X_heatcold.rds"))





em_exceedance_Xplot_df = ward_map %>% 
  left_join(exceed_prob_EM_gt_bham_X_heatcold, by = c("Ward_Code" = "ward22cd"))




ward_pop =  ward_pop %>% 
  filter(`LAD 2022 Name` == "Birmingham") %>% 
  pivot_longer(cols = c(-`LAD 2022 Name`,-`LAD 2022 Name`,-`Ward 2022 Code`,-`Ward 2022 Name`,-Total,
                        -`LAD 2022 Code`),
               names_to = "age",
               values_to = "count") %>% 
  group_by(`Ward 2022 Code`,`Ward 2022 Name`) %>% 
  summarise(count = sum(count), .groups = "drop") %>% 
  rename(Ward_code = `Ward 2022 Code`,
         Ward_name = `Ward 2022 Name`) %>% 
  left_join(ward_X_EM_summary, by ="Ward_code") %>%  #join em summary with pop estimate
  group_by(Ward_code, Ward_name, Ward_id) %>% 
  #standardisation
  mutate(heat_med = heat_med/count*100000,
         heat_LL = heat_LL/count*100000,
         heat_UL = heat_UL/count*100000,
         cold_med = cold_med/count*100000,
         cold_LL = cold_LL/count*100000,
         cold_UL = cold_UL/count*100000,
         #cold tootip
         tooltip_cold  = paste0("<B>Ward Name:</B> ", Ward_name, "\n",
                                "<B>Ward Code:</B> ", Ward_code, "\n",
                                "<B>Median EM:</B> ", round(cold_med,2), "\n",
                                "<B>Upper95CrI:</B> ", round(cold_UL,2), "\n",
                                "<B>Lower95CrI:</B> ", round(cold_LL,2)),
         #heat tooltip
         tooltip_heat  = paste0("<B>Ward Name:</B> ", Ward_name, "\n",
                                "<B>Ward Code:</B> ", Ward_code, "\n",
                                "<B>Median EM:</B> ", round(heat_med,2), "\n",
                                "<B>Upper95CrI:</B> ", round(heat_UL,2), "\n",
                                "<B>Lower95CrI:</B> ", round(heat_LL,2))
  )



tmap_options(component.autoscale = FALSE)


figure_EM_Xcold =
  tm_shape(ward_map %>% 
             left_join(ward_pop, by = c("Ward_Code"="Ward_code")))+
  tm_polygons(
    fill = "cold_med",
    fill.scale = tm_scale_continuous(values = "brewer.blues"),
    fill.legend = tm_legend("per 100,000", group_id = "top", frame = FALSE)
  )+
  tm_layout(
    frame = FALSE,
    legend.position = c("left", "top"),
    legend.frame = FALSE,
    inner.margins = c(0.07, 0, 0.01, 0)
  )+
  tm_title("Estimated non-age-standardised annual median \nexcess mortality rate from extreme cold by ward", size = 1)+
  tm_compass(
    type = "8star",
    size = 4,
    position = c("RIGHT", "bottom"),
    color.light = "white"
  ) +
  tm_credits(
    text = paste(
      "Contains OS data \u00A9 Crown copyright and database right",
      format(Sys.Date(), "%Y"),
      ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."
    ),
    position = c("LEFT", "BOTTOM")
  )



figure_EM_Xheat =
  tm_shape(ward_map %>% 
             left_join(ward_pop, by = c("Ward_Code"="Ward_code")))+
  tm_polygons(
    fill = "heat_med",
    fill.scale = tm_scale_continuous(values = "brewer.reds",ticks = c(2, 4, 6)),
    fill.legend = tm_legend("per 100,000", group_id = "top", frame = FALSE)
  )+
  tm_layout(
    frame = FALSE,
    legend.position = c("left", "top"),
    legend.frame = FALSE,
    inner.margins = c(0.07, 0, 0.01, 0)
  )+
  tm_title("Estimated non-age-standardised annual median \nexcess mortality rate from extreme heat by ward", size = 1)+
  tm_compass(
    type = "8star",
    size = 4,
    position = c("RIGHT", "bottom"),
    color.light = "white"
  ) +
  tm_credits(
    text = paste(
      "Contains OS data \u00A9 Crown copyright and database right",
      format(Sys.Date(), "%Y"),
      ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."
    ),
    position = c("LEFT", "BOTTOM")
  )


figure_EM_prob_Xcold =tm_shape(em_exceedance_Xplot_df)+
  tm_polygons(fill = "p_Xcold_gt_bham",
              fill.scale = tm_scale_continuous(values = "brewer.blues", ticks = c(0.00,0.25,0.50,0.75,1)),
              fill.legend = tm_legend("Pr(EM>mean EM)", group_id = "top", frame = FALSE))+
  tm_layout(
    frame = FALSE,
    legend.position = c("left", "top"),
    legend.frame = FALSE,
    inner.margins = c(0.07, 0, 0.01, 0)
  )+
  tm_title("Probability that excess mortality rate during extreme cold \nis above Birmingham average by ward", size = 1)+
  tm_compass(
    type = "8star",
    size = 4,
    position = c("RIGHT", "bottom"),
    color.light = "white"
  ) +
  tm_credits(
    text = paste(
      "Contains OS data \u00A9 Crown copyright and database right",
      format(Sys.Date(), "%Y"),
      ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."
    ),
    position = c("LEFT", "BOTTOM")
  )


figure_EM_prob_Xheat =tm_shape(em_exceedance_Xplot_df)+
  tm_polygons(fill = "p_Xheat_gt_bham",
              fill.scale = tm_scale_continuous(values = "brewer.reds", ticks = c(0.00,0.25,0.50,0.75,1)),
              fill.legend = tm_legend("Pr(EM>mean EM)", group_id = "top", frame = FALSE))+
  tm_layout(
    frame = FALSE,
    legend.position = c("left", "top"),
    legend.frame = FALSE,
    inner.margins = c(0.07, 0, 0.01, 0)
  )+
  tm_title("Probability that excess mortality rate during extreme heat \nis above Birmingham average by ward", size = 1)+
  tm_compass(
    type = "8star",
    size = 4,
    position = c("RIGHT", "bottom"),
    color.light = "white"
  ) +
  tm_credits(
    text = paste(
      "Contains OS data \u00A9 Crown copyright and database right",
      format(Sys.Date(), "%Y"),
      ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."
    ),
    position = c("LEFT", "BOTTOM")
  )

  


tmap_arrange(figure_EM_Xcold,figure_EM_prob_Xcold)

combined_cold = tmap_arrange(
  figure_EM_Xcold,
  figure_EM_prob_Xcold,
  ncol = 2
)

combined_heat = tmap_arrange(
  figure_EM_Xheat,
  figure_EM_prob_Xheat,
  ncol = 2
)

tmap_save(
  combined_cold,
  filename = "figs/DPH_figure_EM_Xcold_combined.png",
  width = 10,
  height = 6,
  units = "in",
  dpi = 600
)


tmap_save(
  combined_heat,
  filename = "figs/DPH_figure_EM_Xheat_combined.png",
  width = 10,
  height = 6,
  units = "in",
  dpi = 600
)


###########################################################################



tm_shape(ward_map)+
  tm_