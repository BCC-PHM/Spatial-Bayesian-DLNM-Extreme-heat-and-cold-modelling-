#==============================================================================
# CALCULATE HEAT-ATTRIBUTABLE DEATHS
#==============================================================================
library(tidyverse)
library(data.table)
library(bcctheme)

#===========================================================
# just handy function adjusting ggsave
pixel_2_in = function(width,
                      height){
  width_in = width/96
  height_in = height/96
  
  result = c(width_in,height_in)
  return(result)
}



#============================================================
#------------------------------------------------------
#load heat-only cumlogRR files
#------------------------------------------------------
parts_dir = "output/cumlogRR_parts_2026_heatonly_newMMT/"
files = list.files(parts_dir, pattern = "^ward_[0-9]+_cumlogRR_heat\\.rds$", full.names = TRUE)
ward_id = as.integer(stringr::str_extract(basename(files), "(?<=ward_)\\d+"))
file_index = tibble(
  Ward_id = ward_id,
  list_id = as.character(ward_id),
  file    = files
) %>% 
  arrange(Ward_id)
stopifnot(
  nrow(file_index) == 69,
  !anyDuplicated(file_index$Ward_id),
  !anyNA(file_index$Ward_id)
)
#------------------------------------------------------
#load baseline mortality
#------------------------------------------------------
baseline_mort = read_rds("output/baseline_mort_ward.rds") %>% 
  transmute(
    Ward_id       = as.integer(new_id),
    list_id       = as.character(new_id),
    mean_baseline = as.numeric(mean_baseline)
  ) %>% 
  distinct()
stopifnot(
  nrow(baseline_mort) == 69,
  !anyDuplicated(baseline_mort$Ward_id)
)
#------------------------------------------------------
#calculate daily AN by ward
#------------------------------------------------------
heat_list = vector("list", 69)
names(heat_list) = as.character(1:69)

for (k in seq_len(nrow(file_index))) {
  
  i = file_index$Ward_id[k]
  
  cat("\n-----------------------------------------\n")
  cat(sprintf("[%s] %d / %d Processing ward %d\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              k, nrow(file_index), i))
  
  result_i = readRDS(file_index$file[k])
  cumlogRR_i = as.matrix(result_i$cumlogRR)
  baseline_i = baseline_mort %>% 
    filter(Ward_id == i) %>% 
    pull(mean_baseline)
  
  stopifnot(
    length(baseline_i) == 1,
    nrow(cumlogRR_i) == length(result_i$Date)
  )
  #--------------------------------------------
  #backward attributable fraction
  #--------------------------------------------
  AF_i = 1 - exp(-cumlogRR_i)
  AN_i = AF_i * baseline_i
  
  # Optional positive-burden sensitivity analysis:
  #
  # cumlogRR_positive = pmax(cumlogRR_i, 0)
  # AF_i = 1 - exp(-cumlogRR_positive)
  # AN_i = AF_i * baseline_i
  
  colnames(AN_i) = paste0("V", seq_len(ncol(AN_i)))
  ward_name_i = NA_character_
  
  heat_list[[as.character(i)]] = bind_cols(
    tibble(
      ID        = i,
      Date      = as.Date(result_i$Date),
      tasmean   = result_i$tasmean,
      Ward_Code = result_i$Ward_Code,
      list_id   = as.character(i)
    ),
    as.data.frame(AN_i)
  )
  
  cat(sprintf("DONE ward %d | baseline=%.4f | days=%d | draws=%d\n",
              i, baseline_i, nrow(AN_i), ncol(AN_i)))
}
#------------------------------------------------------
#combine all wards
#------------------------------------------------------
heat_daily = data.table::rbindlist(heat_list, fill = TRUE) %>% 
  as_tibble()
#------------------------------------------------------
#define heatwave outcome-date windows
#------------------------------------------------------
heat_windows = heat_daily %>% 
  mutate(period = case_when(
    Date >= as.Date("2026-05-21") & Date <= as.Date("2026-05-29") ~ "May heatwave",
    Date >= as.Date("2026-06-18") & Date <= as.Date("2026-06-30") ~ "June heatwave",
    Date >= as.Date("2026-07-06") & Date <= as.Date("2026-07-16") ~ "July heatwave",
    TRUE                                                          ~ NA_character_
  )) %>% 
  filter(!is.na(period))
#==========================================================
# Birmingham totals by heatwave period
# sum first within each posterior draw
#==========================================================
bham_total_draws = heat_windows %>% 
  group_by(period) %>% 
  summarise(across(starts_with("V"), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>% 
  pivot_longer(
    cols      = starts_with("V"),
    names_to  = "draw",
    values_to = "AN"
  ) %>% 
  mutate(draw = as.integer(sub("V", "", draw)))

bham_total = bham_total_draws %>% 
  group_by(period) %>% 
  summarise(
    AN_LL       = quantile(AN, 0.025, na.rm = TRUE),
    AN_med      = median(AN, na.rm = TRUE),
    AN_UL       = quantile(AN, 0.975, na.rm = TRUE),
    Pr_positive = mean(AN > 0, na.rm = TRUE),
    .groups     = "drop"
  )
print(bham_total)
#------------------------------------------------------
#Birmingham daily time series
#------------------------------------------------------
bham_daily_draws = heat_daily %>% 
  group_by(Date) %>% 
  summarise(across(starts_with("V"), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>% 
  pivot_longer(
    cols      = starts_with("V"),
    names_to  = "draw",
    values_to = "AN"
  ) %>% 
  mutate(draw = as.integer(sub("V", "", draw)))

daily_ts = bham_daily_draws %>% 
  group_by(Date) %>% 
  summarise(
    AN_LL       = quantile(AN, 0.025, na.rm = TRUE),
    AN_med      = median(AN, na.rm = TRUE),
    AN_UL       = quantile(AN, 0.975, na.rm = TRUE),
    Pr_positive = mean(AN > 0, na.rm = TRUE),
    .groups     = "drop"
  )
print(daily_ts)

june_heat_plot = daily_ts %>% 
  mutate(
    period = case_when(
      Date >= as.Date("2026-06-18") & Date <= as.Date("2026-06-30") ~ "June heatwave",
      TRUE                                                          ~ NA_character_
    ),
    sufficient_evidence = if_else(Pr_positive >= 0.95, "Sufficient evidence", "Insufficient evidence")
  ) %>% 
  filter(!is.na(period)) %>% 
  filter(!is.na(period)) %>% 
  ggplot(aes(x = Date, y = AN_med)) +
  geom_ribbon(aes(ymin = AN_LL, ymax = AN_UL), alpha = 0.2,fill = "#E9997E") +
  geom_line(colour="#DC582A",linewidth=1) +
  geom_point(aes(colour = sufficient_evidence), size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c(
    "Sufficient evidence"   = "#DC582A",
    "Insufficient evidence" = "grey30"
  )) +
  scale_x_date(date_breaks = "1 day", date_labels = "%d %b") +
  labs(
    title    = "Estimated number of heat-attributable deaths during the \nJune heatwave in Birmingham",
    x        = "Date",
    y        = "Heat-attributable deaths",
    colour   = NULL
  ) +
  theme_bcc(base_size = 11,
    gridline_x = FALSE) +
  theme(
    axis.text.x     = element_text(angle = 0, hjust = 1),
    legend.position = "bottom"
  )

ggsave("figs/DPH_2026/june_heat_plot.png",
       june_heat_plot  ,
       width = pixel_2_in(707,488)[1], height = pixel_2_in(707,488)[2],
       units = "in",
       dpi=600)


bcc_pal("orange")(10)

july_heat_plot = daily_ts %>% 
  mutate(
    period = case_when(
      Date >= as.Date("2026-07-06") & Date <= as.Date("2026-07-17") ~ "July heatwave",
      TRUE                                                          ~ NA_character_
    ),
    sufficient_evidence = if_else(Pr_positive >= 0.95, "Sufficient evidence", "Insufficient evidence")
  ) %>% 
  filter(!is.na(period)) %>% 
  filter(!is.na(period)) %>% 
  ggplot(aes(x = Date, y = AN_med)) +
  geom_ribbon(aes(ymin = AN_LL, ymax = AN_UL), alpha = 0.2,fill = "#E9997E") +
  geom_line(colour="#DC582A",linewidth=1) +
  geom_point(aes(colour = sufficient_evidence), size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c(
    "Sufficient evidence"   = "#DC582A",
    "Insufficient evidence" = "grey30"
  )) +
  scale_x_date(date_breaks = "1 day", date_labels = "%d %b") +
  labs(
    title    = "Estimated number of heat-attributable deaths during the \nJuly heatwave in Birmingham",
    x        = "Date",
    y        = "Heat-attributable deaths",
    colour   = NULL
  ) +
  theme_bcc(base_size = 11,
            gridline_x = FALSE) +
  theme(
    axis.text.x     = element_text(angle = 0, hjust = 1),
    legend.position = "bottom"
  )
  


ggsave("figs/DPH_2026/july_heat_plot.png",
       july_heat_plot  ,
       width = pixel_2_in(707,488)[1], height = pixel_2_in(707,488)[2],
       units = "in",
       dpi=600)



#------------------------------------------------------
#ward-level totals
#------------------------------------------------------
ward_total_draws = heat_windows %>% 
  group_by(ID, Ward_Code, period) %>% 
  summarise(across(starts_with("V"), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>% 
  pivot_longer(
    cols      = starts_with("V"),
    names_to  = "draw",
    values_to = "AN"
  ) %>% 
  mutate(draw = as.integer(sub("V", "", draw)))

ward_total = ward_total_draws %>% 
  group_by(ID, Ward_Code, period) %>% 
  summarise(
    AN_LL       = quantile(AN, 0.025, na.rm = TRUE),
    AN_med      = median(AN, na.rm = TRUE),
    AN_UL       = quantile(AN, 0.975, na.rm = TRUE),
    Pr_positive = mean(AN > 0, na.rm = TRUE),
    .groups     = "drop"
  )
#==========================================================
# save results
#==========================================================
write_rds(heat_daily, "output/heat_daily_2026_newMMT.rds")
write_rds(daily_ts, "output/Birmingham_daily_heat_AN_2026_newMMT.rds")
write_rds(bham_total, "output/Birmingham_heatwave_total_AN_2026_newMMT.rds")
write_rds(ward_total, "output/Ward_heatwave_total_AN_2026_newMMT.rds")
write_csv(bham_total, "output/Birmingham_heatwave_total_AN_2026_newMMT.csv")
write_csv(ward_total, "output/Ward_heatwave_total_AN_2026_newMMT.csv")







