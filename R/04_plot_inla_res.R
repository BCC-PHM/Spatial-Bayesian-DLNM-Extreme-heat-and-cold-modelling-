setwd("D:/R project D/Extreme heat and cold")

library(tidyverse)
library(dlnm)
library("sf")
library(qpdf)
library(tmap)
library(qpdf)

#------------------------------------------------------
#load the data
#------------------------------------------------------
df_complete = read_rds("output/03_df_complete_for_inla.rds")
cb_res = read_rds("output/03_predicted_inla_spatial_casecrossover.rds")

#------------------------------------------------------
#load ward shp file
#------------------------------------------------------
ward_map = read_sf("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp")

#------------------------------------------------------
#percentile grid, high resolution at the tails
#------------------------------------------------------
percentiles = c(
  seq(0, 1, by = 0.1),    #extreme lows (0.0% to 1.0%)
  seq(2, 98, by = 1),     #core distribution (2.0% to 98.0%)
  seq(99, 100, by = 0.1)  #extreme highs (99.0% to 100.0%)
) / 100

#------------------------------------------------------
#ward-specific temperature distributions
#------------------------------------------------------
x_temp = vector("list", 69)

for (i in 1:69) {
  #subset data for the current ward id
  df_w = df_complete[df_complete$new_id == i, ]
  
  #temperature values at the custom percentiles
  x_temp[[i]] = quantile(df_w$tasmean, percentiles, na.rm = TRUE)
}

#------------------------------------------------------
#cross-basis matrices for prediction
#------------------------------------------------------
cb_pred = vector("list", 69)

for (i in 1:69) {
  #subset data again to access ward-specific knots
  df_w = df_complete[df_complete$new_id == i, ]
  
  #interaction between exposure and lag space
  #21 + 1 accounts for lags 0 through 21
  cb_pred[[i]] = crossbasis(
    matrix(rep(x_temp[[i]], 21 + 1), ncol = 21 + 1),
    argvar = list(
      fun   = "bs",
      knots = quantile(df_w$tasmean, probs = c(0.1, 0.75, 0.9), na.rm = TRUE)
    ),
    arglag = list(fun = "ns", knots = logknots(21, 3), intercept = TRUE)
  )
}

#======================================================
#posterior draws for ward id 1
#======================================================
cb_res[[1]]

#------------------------------------------------------
#posterior samples of the coefficients, 1000 x 30
#------------------------------------------------------
beta_reg = cb_res[[1]]

#------------------------------------------------------
#cross-basis prediction matrix, 119 x 30
#------------------------------------------------------
cb_i = cb_pred[[1]]

rr = cb_i %*% t(beta_reg)

#------------------------------------------------------
#reference temperature at the 90th percentile
#------------------------------------------------------
i_cen = which(percentiles == 0.9)

#------------------------------------------------------
#centre each draw so log-rr is 0 at the reference
#------------------------------------------------------
rr_cen = apply(rr, 2, function(x) x - x[i_cen])

full_1000_rr = exp(rr_cen)

#------------------------------------------------------
#median and 95% credible interval
#------------------------------------------------------
RR_med = apply(exp(rr_cen), 1, median)
RR_lo  = apply(exp(rr_cen), 1, quantile, probs = 0.025, na.rm = TRUE)
RR_hi  = apply(exp(rr_cen), 1, quantile, probs = 0.975, na.rm = TRUE)

#------------------------------------------------------
#tidy the draws for plotting
#------------------------------------------------------
full_ward_1_plot_df = as.data.frame(full_1000_rr)
colnames(full_ward_1_plot_df) = gsub("V", "RR_", colnames(full_ward_1_plot_df))

individuals_line_plot = data.frame(
  Percentage = names(x_temp[[i]]),
  Temp       = as.numeric(x_temp[[i]]),
  RR_median  = RR_med,
  RR_LL      = RR_lo,
  RR_UL      = RR_hi
)

ward_map$Ward_Name[ward_map$Ward_Code == unique(df_complete$ward22cd[df_complete$new_id == 1])]

#------------------------------------------------------
#spaghetti plot of the 1000 draws
#------------------------------------------------------
posterior_draws_for_acocks_green = cbind(
  data.frame(Percentage = names(x_temp[[i]]), Temp = as.numeric(x_temp[[i]])),
  full_ward_1_plot_df
) %>% 
  pivot_longer(cols = c(-Temp, -Percentage), names_to = "Sim", values_to = "RR") %>% 
  ggplot(aes(x = Temp, y = RR, group = Sim)) +
  #the 1000 grey lines
  geom_line(colour = "grey90", size = 0.6) +
  #reference lines
  geom_hline(aes(yintercept = 1), colour = "black") +
  geom_vline(aes(xintercept = as.numeric(x_temp[[1]][[100]])), linetype = "dotted", size = 0.2, colour = "black") +
  #median and cis
  geom_line(data = individuals_line_plot, aes(x = Temp, y = RR_median), inherit.aes = FALSE, size = 0.8) +
  geom_line(data = individuals_line_plot, aes(x = Temp, y = RR_LL), inherit.aes = FALSE, size = 0.8, linetype = "dashed") +
  geom_line(data = individuals_line_plot, aes(x = Temp, y = RR_UL), inherit.aes = FALSE, size = 0.8, linetype = "dashed") +
  #annotate() for the text label
  annotate(
    "text",
    x      = as.numeric(x_temp[[1]][[100]]),
    y      = 3,
    label  = paste0("Centering         \ntemperature = ", round(as.numeric(x_temp[[1]][[100]]), 2)),
    angle  = 0,
    vjust  = 0,
    hjust  = 1.2,
    size   = 3,
    colour = "black"
  ) +
  #scales and theme
  scale_y_continuous(breaks = 0:6) +
  scale_x_continuous(breaks = seq(-5, 30, 5)) +
  coord_cartesian(ylim = c(0, 6), xlim = c(-5, 30)) +
  theme_classic(base_size = 11) +
  labs(x = "Temperature °C", y = "Relative Risk") +
  ggtitle(paste0(
    "Posterior draws of exposure-response curves for ",
    ward_map$Ward_Name[ward_map$Ward_Code == unique(df_complete$ward22cd[df_complete$new_id == 1])]
  )) +
  theme(
    plot.margin     = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5)
  )

ggsave("figs/04_posterior_draws_for_acocks_green.jpg", posterior_draws_for_acocks_green, width = 6, height = 4, dpi = 600, units = "in")

#======================================================
#mmt-centred curves for every ward
#======================================================
mmt_draws_by_ward = vector("list", 69)
RR_MMT_plot_list = list()
rr_mmt_centered = vector("list", 69)
log_rr_mmt_centered = vector("list", 69)

for (w in 1:69) {
  
  #----------------------------------------------------
  #select location index
  #----------------------------------------------------
  i = w
  current_ward_code = unique(df_complete$ward22cd[df_complete$new_id == i])
  current_ward_name = ward_map$Ward_Name[ward_map$Ward_Code == current_ward_code]
  
  #----------------------------------------------------
  #posterior samples from inla results, 1000 x 30
  #----------------------------------------------------
  beta_reg = cb_res[[i]]
  
  #----------------------------------------------------
  #cross-basis prediction matrix, 119 x 30
  #----------------------------------------------------
  cb_i = cb_pred[[i]]
  
  #transpose beta_reg to 30 x 1000 so it aligns with cb_i,
  #giving a 119 x 1000 matrix
  rr = cb_i %*% t(beta_reg)
  
  #----------------------------------------------------
  #reference temperature at the 90th percentile
  #----------------------------------------------------
  i_cen = which(percentiles == 0.9)
  
  #----------------------------------------------------
  #centre each draw at the reference
  #----------------------------------------------------
  rr_cen = apply(rr, 2, function(x) x - x[i_cen])
  
  #----------------------------------------------------
  #median and 95% credible interval on the rr scale
  #----------------------------------------------------
  RR_med = apply(exp(rr_cen), 1, median)
  RR_lo  = apply(exp(rr_cen), 1, quantile, probs = 0.025, na.rm = TRUE)
  RR_hi  = apply(exp(rr_cen), 1, quantile, probs = 0.975, na.rm = TRUE)
  
  #----------------------------------------------------
  #find mmt from posterior median log-risk
  #----------------------------------------------------
  logRR_med = apply(rr, 1, median, na.rm = TRUE)
  
  valid_idx = which(percentiles >= 0.01 & percentiles <= 0.99)
  
  i_mmt     = valid_idx[which.min(logRR_med[valid_idx])]
  mmt_value = x_temp[[i]][i_mmt]
  min_risk  = exp(logRR_med[i_mmt])
  
  #----------------------------------------------------
  #centre all draws at the ward mmt
  #----------------------------------------------------
  mmt_position_draws = apply(rr[valid_idx, , drop = FALSE], 2, which.min)
  mmt_position_draws = valid_idx[mmt_position_draws]
  
  mmt_draws_by_ward[[i]] = tibble(
    nsim      = seq_len(ncol(rr)),
    MMT       = as.numeric(x_temp[[i]][mmt_position_draws]),
    Ward_code = current_ward_code,
    Ward_id   = i
  )
  
  rr_cen_mmt = sapply(seq_len(ncol(rr)), function(j) rr[, j] - rr[mmt_position_draws[j], j])
  
  #keep it as logrr
  log_rr_mmt_centered[[i]] = rr_cen_mmt
  
  #exponentiate it and store it
  rr_mmt_centered[[i]] = exp(rr_cen_mmt)
  
  #----------------------------------------------------
  #recalculate median and cis on the new centre
  #----------------------------------------------------
  RR_med_mmt = apply(exp(rr_cen_mmt), 1, median)
  RR_lo_mmt  = apply(exp(rr_cen_mmt), 1, quantile, 0.025, na.rm = TRUE)
  RR_hi_mmt  = apply(exp(rr_cen_mmt), 1, quantile, 0.975, na.rm = TRUE)
  
  #----------------------------------------------------
  #assemble the ward plot frame
  #----------------------------------------------------
  plot_df = data.frame(
    Percentage = names(x_temp[[i]]),
    Temp       = as.numeric(x_temp[[i]]),
    RR_median  = RR_med_mmt,
    RR_LL      = RR_lo_mmt,
    RR_UL      = RR_hi_mmt,
    MMT        = as.numeric(mmt_value),
    Ward_id    = i,
    Ward_code  = current_ward_code,
    Ward_name  = current_ward_name
  ) %>% 
    mutate(
      legend_label = ifelse(
        Percentage %in% c("1.0%", "99.0%"),
        paste0(
          ifelse(Percentage == "1.0%", " P1", " P99"),
          ": RR ",
          sprintf("%.2f", RR_median),
          " (",
          sprintf("%.2f", RR_LL),
          "–",
          sprintf("%.2f", RR_UL),
          ")"
        ),
        NA
      )
    )
  
  RR_MMT_plot_list[[w]] = plot_df
}

#------------------------------------------------------
#bind the ward frames and save
#------------------------------------------------------
real_plot_df = data.table::rbindlist(RR_MMT_plot_list)

write_rds(real_plot_df, file = "output/04_RR_MMT_plot_data.rds")
write_rds(mmt_draws_by_ward, "output/04_mmt_draws_by_ward_new.rds")
#------------------------------------------------------
#batch starts, every group of 15: 1, 16, 31, 46, 61
#------------------------------------------------------
starts = seq(1, 69, by = 15)

for (start_id in starts) {
  
  #----------------------------------------------------
  #range for this batch
  #----------------------------------------------------
  end_id = min(start_id + 14, 69)
  current_ids = start_id:end_id
  
  batch_data = subset(real_plot_df, Ward_id %in% current_ids)
  
  #----------------------------------------------------
  #dynamic height, rows needed rounding up by 3
  #----------------------------------------------------
  n_wards = length(unique(batch_data$Ward_id))
  n_rows = ceiling(n_wards / 3)
  
  #----------------------------------------------------
  #create the plot
  #----------------------------------------------------
  p = ggplot(batch_data, aes(Temp, y = RR_median)) +
    geom_ribbon(aes(ymin = RR_LL, ymax = RR_UL), alpha = 0.2, fill = "#A092AA", colour = NA) +
    geom_line(size = 0.8, colour = "#A092AA") +
    geom_hline(yintercept = 1) +
    geom_vline(data = batch_data, aes(xintercept = MMT), linetype = "dashed", size = 0.5, alpha = 0.5) +
    geom_point(
      data         = batch_data,
      aes(x = -4.8, y = ifelse(Percentage == "1.0%", 5.78, 5.18), shape = Percentage),
      size         = 2,
      fill         = "#7F6588",
      colour       = "black",
      stroke       = 0.8,
      inherit.aes  = FALSE
    ) +
    geom_point(data = batch_data, aes(shape = Percentage), size = 3, fill = "#7F6588", colour = "black", stroke = 0.8) +
    geom_text(
      data        = batch_data,
      aes(x = -4.0, y = ifelse(Percentage == "1.0%", 5.8, 5.2), label = legend_label),
      hjust       = 0,
      size        = 3,
      inherit.aes = FALSE,
      fontface    = "bold"
    ) +
    geom_text(
      data        = batch_data %>% group_by(Ward_name) %>% summarise(MMT = mean(MMT)),
      aes(x = MMT, y = 3, label = paste0("MMT = ", sprintf("%.1f", MMT))),
      angle       = 0,
      vjust       = 0,
      hjust       = 1.2,
      size        = 3,
      colour      = "black",
      inherit.aes = FALSE
    ) +
    #let ggplot decide the layout now that we control the file size
    facet_wrap(~ Ward_name, ncol = 3, scales = "free") +
    scale_y_continuous(breaks = 0:6) +
    scale_x_continuous(breaks = seq(-5, 30, 5)) +
    coord_cartesian(ylim = c(0, 6), xlim = c(-5, 30)) +
    scale_shape_manual(values = c("1.0%" = 21, "99.0%" = 24)) +
    theme_classic(base_size = 11) +
    labs(x = "Temperature °C", y = "Relative Risk") +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), legend.position = "none")
  
  #----------------------------------------------------
  #save with dynamic height
  #a full page (5 rows) is 297mm, so one row is approx 59.4mm
  #----------------------------------------------------
  dynamic_height = 297 * (n_rows / 5)
  
  filename = paste0("figs/", "04_Overall_cumul_RR_ward_plots_", start_id, "_to_", end_id, ".pdf")
  
  ggsave(filename, plot = p, width = 210, height = dynamic_height, units = "mm")
  
  message(paste("Saved:", filename, "(Height:", round(dynamic_height, 1), "mm)"))
}

#------------------------------------------------------
#merge the files (optional)
#------------------------------------------------------
pdf_combine(
  input = c(
    "figs/04_Overall_cumul_RR_ward_plots_1_to_15.pdf",
    "figs/04_Overall_cumul_RR_ward_plots_16_to_30.pdf",
    "figs/04_Overall_cumul_RR_ward_plots_31_to_45.pdf",
    "figs/04_Overall_cumul_RR_ward_plots_46_to_60.pdf",
    "figs/04_Overall_cumul_RR_ward_plots_61_to_69.pdf"
  ),
  output = "figs/04_Combined_overall_cumul_RR_wards.pdf"
)

real_plot_df %>% 
  filter(Ward_id == 1)

#======================================================
#map the mmt of the median ensemble response curve
#======================================================
mmt_plot = real_plot_df %>% 
  group_by(Ward_name, Ward_code) %>% 
  summarise(MMT = mean(MMT, .group = "drop"))

tmap::tmap_mode("plot")

mmt_plot = tm_shape(ward_map %>% left_join(mmt_plot, by = c("Ward_Code" = "Ward_code"))) +
  tm_polygons(
    fill        = "MMT",
    fill.scale  = tm_scale_continuous(values = "brewer.RdPu"),
    fill.legend = tm_legend("MMT", group_id = "top", frame = FALSE, bg.alpha = 0)
  ) +
  tm_layout(
    title.size      = 1.2,
    legend.position = c(0.02, 0.92),
    frame           = FALSE,
    #increased bottom margin (3rd number) to make room for caption
    inner.margins   = c(0.07, 0, 0.15, 0)
  ) +
  tm_title("Minimum Mortality Temperature") +
  tm_credits(
    text     = "Minimum Mortality Temperature of the ward temperature distribution \nat the median exposure response curve.",
    position = c("LEFT", "TOP")
  ) +
  tm_compass(type = "8star", size = 4, position = c("RIGHT", "bottom"), color.light = "white") +
  tm_credits(
    text = paste(
      "Contains OS data \u00A9 Crown copyright and database right",
      #get current year
      format(Sys.Date(), "%Y"),
      ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."
    ),
    position = c("LEFT", "BOTTOM")
  )

tmap_save(mmt_plot, filename = "figs/04_mmt_plot.png", height = 7, width = 6, unit = "in", dpi = 600)

#======================================================
#map rr of the median curve at the 99th percentile
#======================================================
RR_99th_map = real_plot_df %>% 
  filter(Percentage == "99.0%")

RR_99th_map = tm_shape(ward_map %>% left_join(RR_99th_map, by = c("Ward_Code" = "Ward_code"))) +
  tm_polygons(
    fill        = "RR_median",
    fill.scale  = tm_scale_continuous(values = "reds"),
    fill.legend = tm_legend("Relative Risk", group_id = "top", frame = FALSE, bg.alpha = 0)
  ) +
  tm_layout(
    title.size      = 1.2,
    legend.position = c(0.02, 0.92),
    frame           = FALSE,
    #increased bottom margin (3rd number) to make room for caption
    inner.margins   = c(0.07, 0, 0.15, 0)
  ) +
  tm_title("Relative Risk at 99th Percentile Temperature") +
  tm_credits(
    text     = "Median Relative Risk of death at the 99th percentile of the ward temperature distribution \ncompared to the risk at the MMT of the median exposure response curve.",
    position = c("LEFT", "TOP")
  ) +
  tm_compass(type = "8star", size = 4, position = c("RIGHT", "bottom"), color.light = "white") +
  tm_credits(
    text = paste(
      "Contains OS data \u00A9 Crown copyright and database right",
      #get current year
      format(Sys.Date(), "%Y"),
      ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."
    ),
    position = c("LEFT", "BOTTOM")
  )

tmap_save(RR_99th_map, filename = "figs/04_RR_99th_map.png", height = 7, width = 6, unit = "in", dpi = 600)

#======================================================
#map rr of the median curve at the 1st percentile
#======================================================
RR_1st_map = real_plot_df %>% 
  filter(Percentage == "1.0%")

tmap::tmap_mode("plot")

RR_1st_map = tm_shape(ward_map %>% left_join(RR_1st_map, by = c("Ward_Code" = "Ward_code"))) +
  tm_polygons(
    fill        = "RR_median",
    fill.scale  = tm_scale_continuous(values = "blues"),
    fill.legend = tm_legend("Relative Risk", group_id = "top", frame = FALSE, bg.alpha = 0)
  ) +
  tm_layout(
    title.size      = 1.2,
    legend.position = c(0.02, 0.92),
    frame           = FALSE,
    #increased bottom margin (3rd number) to make room for caption
    inner.margins   = c(0.07, 0, 0.15, 0)
  ) +
  tm_title("Relative Risk at 1st Percentile Temperature") +
  tm_credits(
    text     = "Median Relative Risk of death at the 1st percentile of the ward temperature distribution \ncompared to the risk at the MMT of the median exposure response curve.",
    position = c("LEFT", "TOP")
  ) +
  tm_compass(type = "8star", size = 4, position = c("RIGHT", "bottom"), color.light = "white") +
  tm_credits(
    text = paste(
      "Contains OS data \u00A9 Crown copyright and database right",
      #get current year
      format(Sys.Date(), "%Y"),
      ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."
    ),
    position = c("LEFT", "BOTTOM")
  )

tmap_save(RR_1st_map, filename = "figs/04_RR_1st_map.png", height = 7, width = 6, unit = "in", dpi = 600)

#------------------------------------------------------
#merge the plots
#------------------------------------------------------
merged_1st_99th_RR_plot = tmap::tmap_arrange(RR_1st_map, RR_99th_map)

tmap_save(merged_1st_99th_RR_plot, filename = "figs/04_merged_1st_99th_RR_plot.png", height = 7, width = 12, unit = "in", dpi = 600)

#======================================================
#mmt posterior uncertainty by ward
#rr is [n_temp x n_draws], rows = temp, cols = draws
#======================================================
mmt_uncertainty_plot_df = data.table::rbindlist(mmt_draws_by_ward) %>% 
  group_by(Ward_code, Ward_id) %>% 
  mutate(
    prob_gt_18   = round(mean(MMT > 18, na.rm = TRUE) * 100, 1),
    upperCrI     = quantile(MMT, 0.975),
    lowerCrI     = quantile(MMT, 0.025),
    label_pr18   = paste0("Pr(MMT>18) ", "\n= ", prob_gt_18, "%"),
    label_crI    = paste0("MMT CrI ", "95%: ", "(", round(lowerCrI, 1), ", ", round(upperCrI, 1), ")"),
    bin_boundary = floor(min(MMT, na.rm = TRUE)),
    ymax         = max(
      c(
        hist(MMT, breaks = seq(floor(min(MMT)), ceiling(max(MMT)), by = 1), plot = FALSE)$density,
        #kernel density heights
        density(MMT, adjust = 1)$y
      ),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>% 
  mutate(
    label_pr18 = ifelse(nsim == 1, label_pr18, NA),
    label_crI  = ifelse(nsim == 1, label_crI, NA)
  ) %>% 
  left_join(ward_map, by = c("Ward_code" = "Ward_Code"))

#------------------------------------------------------
#batch the mmt uncertainty panels
#------------------------------------------------------
starts = seq(1, 69, by = 15)

for (start_id in starts) {
  
  #----------------------------------------------------
  #range for this batch
  #----------------------------------------------------
  end_id = min(start_id + 14, 69)
  current_ids = start_id:end_id
  
  batch_data = subset(mmt_uncertainty_plot_df, Ward_id %in% current_ids)
  
  #----------------------------------------------------
  #dynamic height, rows needed rounding up by 3
  #----------------------------------------------------
  n_wards = length(unique(batch_data$Ward_id))
  n_rows = ceiling(n_wards / 3)
  
  #----------------------------------------------------
  #create the plot
  #----------------------------------------------------
  p = batch_data %>% 
    ggplot(aes(x = MMT)) +
    geom_histogram(aes(y = after_stat(density)), colour = "grey40", fill = "grey90", binwidth = 1) +
    geom_density(lwd = 1, colour = "#578f9f", fill = "#83b3c0", alpha = 0.5, adjust = 1) +
    geom_vline(xintercept = 18) +
    geom_segment(
      aes(x = 18, y = ymax * 1.18, xend = 28, yend = ymax * 1.18),
      arrow     = arrow(type = "open", length = unit(0.1, "cm")),
      linewidth = 0.15
    ) +
    geom_text(aes(x = 19, y = ymax * 1.3, label = label_pr18), hjust = 0, size = 2.5, inherit.aes = FALSE, fontface = "bold") +
    geom_text(aes(x = -5, y = ymax * 1.3, label = label_crI), hjust = 0, size = 2.5, inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq(-5, 30, 5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    coord_cartesian(xlim = c(-5, 30)) +
    labs(x = "Temperature °C", y = "Density posterior distribution MMT") +
    theme_classic(base_size = 11) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), legend.position = "none") +
    facet_wrap(~ Ward_Name, ncol = 3, scales = "free")
  
  #----------------------------------------------------
  #save with dynamic height
  #a full page (5 rows) is 297mm, so one row is approx 59.4mm
  #----------------------------------------------------
  dynamic_height = 297 * (n_rows / 5)
  
  filename = paste0("figs/", "04_MMT_uncertainty_ward_plots_", start_id, "_to_", end_id, ".pdf")
  
  ggsave(filename, plot = p, width = 210, height = dynamic_height, units = "mm")
  
  message(paste("Saved:", filename, "(Height:", round(dynamic_height, 1), "mm)"))
}

#------------------------------------------------------
#merge the files (optional)
#------------------------------------------------------
pdf_combine(
  input = c(
    "figs/04_MMT_uncertainty_ward_plots_1_to_15.pdf",
    "figs/04_MMT_uncertainty_ward_plots_16_to_30.pdf",
    "figs/04_MMT_uncertainty_ward_plots_31_to_45.pdf",
    "figs/04_MMT_uncertainty_ward_plots_46_to_60.pdf",
    "figs/04_MMT_uncertainty_ward_plots_61_to_69.pdf"
  ),
  output = "figs/04_Combined_MMT_uncertainty_wards.pdf"
)

#======================================================
#exceedance prob for rr at p99/p1 > 1
#======================================================
RR_Prob_plot_list = list()

for (i in 1:69) {
  
  current_ward_code = unique(df_complete$ward22cd[df_complete$new_id == i])
  current_ward_name = ward_map$Ward_Name[ward_map$Ward_Code == current_ward_code]
  
  idx_99 = which(names(x_temp[[i]]) == "99.0%")
  idx_1 = which(names(x_temp[[i]]) == "1.0%")
  
  RR_draws_P1 = rr_mmt_centered[[i]][idx_1, ]
  RR_draws_P99 = rr_mmt_centered[[i]][idx_99, ]
  
  plot_df = data.frame(
    RR_draws_P1  = mean(RR_draws_P1 > 1),
    RR_draws_P99 = mean(RR_draws_P99 > 1),
    Ward_id      = i,
    Ward_code    = current_ward_code,
    Ward_name    = current_ward_name
  )
  
  RR_Prob_plot_list[[i]] = plot_df
}

RR_prob_df = dplyr::bind_rows(RR_Prob_plot_list)

write_rds(RR_prob_df, "output/04_RR_exceedance_prob.rds")


#------------------------------------------------------
#join to the boundaries and build tooltip fields
#------------------------------------------------------
plot_sf = ward_map %>%
  left_join(RR_prob_df, by = c("Ward_Code" = "Ward_code")) %>%
  mutate(
    RR_draws_P99        = as.numeric(RR_draws_P99),
    RR_draws_P99_grey90 = if_else(RR_draws_P99 < 0.90, NA_real_, RR_draws_P99),
    evidence_99th       = case_when(
      is.na(RR_draws_P99)  ~ "Missing",
      RR_draws_P99 >= 0.95 ~ "Strong evidence (≥0.95)",
      RR_draws_P99 >= 0.90 ~ "Some evidence (0.90–0.95)",
      TRUE                 ~ "No evidence (<0.90)"
    ),
    tooltip_99th        = paste0(
      "<B>Pr(RR at P99)>1:</B> ", round(RR_draws_P99, 2),
      "\n<B>Inference:</B> ", evidence_99th
    ),
    RR_draws_P1         = as.numeric(RR_draws_P1),
    RR_draws_P1_grey90  = if_else(RR_draws_P1 < 0.90, NA_real_, RR_draws_P1),
    evidence_P1         = case_when(
      is.na(RR_draws_P1)  ~ "Missing",
      RR_draws_P1 >= 0.95 ~ "Strong evidence (≥0.95)",
      RR_draws_P1 >= 0.90 ~ "Some evidence (0.90–0.95)",
      TRUE                ~ "No evidence (<0.90)"
    ),
    tooltip_P1          = paste0(
      "<B>Pr(RR at P1)>1:</B> ", round(RR_draws_P1, 3),
      "\n<B>Inference:</B> ", evidence_P1
    )
  )

#------------------------------------------------------
#custom css
#------------------------------------------------------
tooltip_css = "
  background: rgba(255, 255, 255, 0.97);
  color: #1f2937;
  padding: 10px 12px;
  border-radius: 12px;
  border: 1px solid rgba(17, 24, 39, 0.12);
  box-shadow: 0 10px 24px rgba(0, 0, 0, 0.18);
  font-size: 16px;
  line-height: 1.35;
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Arial, sans-serif;
"

hover_css = "
  cursor: pointer;
  stroke: #111827 ;    /* darker border */
  stroke-width: 1.5px ;  /* thicker border */
  opacity: 1 ;         /* keep fill the same */
  transition: all 0.15s ease-out;
"

#------------------------------------------------------
#interactive exceedance maps
#------------------------------------------------------
figure_1_Prob = ggplot(plot_sf) +
  geom_sf_interactive(aes(fill = RR_draws_P99_grey90, data_id = Ward_Code, tooltip = tooltip_99th), color = "white") +
  scale_fill_distiller(
    palette   = "Reds",
    direction = 1,
    na.value  = "grey85",
    #keep the legend focused on 0.90-1.00
    limits    = c(0.90, 1.00),
    breaks    = c(0.90, 0.95, 1.00)
  ) +
  labs(fill = "Probability") +
  ggtitle("Exceedance probability that relative risk exceeds 1 \nat the 99th percentile temperature") +
  theme_void(base_size = 16) +
  theme()

figure_2_Prob = ggplot(plot_sf) +
  geom_sf_interactive(aes(fill = RR_draws_P1_grey90, data_id = Ward_Code, tooltip = tooltip_P1), color = "white") +
  scale_fill_distiller(
    palette   = "Blues",
    direction = 1,
    na.value  = "grey85",
    #keep the legend focused on 0.90-1.00
    limits    = c(0.90, 1.00),
    breaks    = c(0.90, 0.95, 1.00)
  ) +
  labs(fill = "Probability") +
  ggtitle("Exceedance probability that relative risk exceeds 1 \nat the 99th percentile temperature") +
  theme_void(base_size = 16) +
  theme()

ggiraph::girafe(
  ggobj      = figure_1_Prob,
  width_svg  = 14,
  height_svg = 10,
  options    = list(opts_hover(css = hover_css), opts_tooltip(css = tooltip_css), opts_hover_inv(css = ""))
)

ggiraph::girafe(
  ggobj      = figure_2_Prob,
  width_svg  = 14,
  height_svg = 10,
  options    = list(opts_hover(css = hover_css), opts_tooltip(css = tooltip_css), opts_hover_inv(css = ""))
)

#======================================================
#combine wards into one array
#ward x temperature percentile x posterior draw
#======================================================
n_wards = length(cb_res)
n_temp  = length(x_temp[[1]])
n_sim   = nrow(cb_res[[1]])

ward_logRR_array = array(NA_real_, dim = c(n_wards, n_temp, n_sim))

for (i in seq_len(n_wards)) {
  stopifnot(
    nrow(log_rr_mmt_centered[[i]]) == n_temp,
    ncol(log_rr_mmt_centered[[i]]) == n_sim
  )
  
  ward_logRR_array[i, , ] = log_rr_mmt_centered[[i]]
}

stopifnot(!anyNA(ward_logRR_array))

centering_check = lapply(log_rr_mmt_centered, function(x) apply(x, 2, min, na.rm = TRUE))

summary(unlist(centering_check))

#------------------------------------------------------
#birmingham overall cumulative log-rr
#equal weighting across wards, each draw aggregated separately
#------------------------------------------------------
ward_weights = rep(1 / n_wards, n_wards)

bham_logRR_draws = apply(ward_logRR_array, MARGIN = c(2, 3), FUN = function(x) sum(x * ward_weights))

#dimensions should be temperature percentiles x posterior draws
stopifnot(
  nrow(bham_logRR_draws) == n_temp,
  ncol(bham_logRR_draws) == n_sim
)

bham_RR_draws = exp(bham_logRR_draws)

#------------------------------------------------------
#birmingham posterior summary at every percentile
#------------------------------------------------------
bham_RR_by_pct = tibble(
  Percentile = names(x_temp[[1]]),
  RR_median  = apply(bham_RR_draws, 1, median, na.rm = TRUE),
  RR_lo      = apply(bham_RR_draws, 1, quantile, probs = 0.025, na.rm = TRUE),
  RR_hi      = apply(bham_RR_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
)

write_rds(bham_RR_by_pct, "output/04_Birmingham_overall_RR_by_percentile.rds")

#------------------------------------------------------
#ward exceedance probability relative to birmingham
#for each draw: ward log-rr > birmingham overall log-rr
#------------------------------------------------------
exceedance_prob_list = vector("list", n_wards)

for (i in seq_len(n_wards)) {
  
  ward_logRR = log_rr_mmt_centered[[i]]
  
  stopifnot(identical(dim(ward_logRR), dim(bham_logRR_draws)))
  
  exceedance_vec = rowMeans(ward_logRR > bham_logRR_draws, na.rm = TRUE)
  
  current_ward_code = unique(df_complete$ward22cd[df_complete$new_id == i])
  current_ward_name = ward_map$Ward_Name[ward_map$Ward_Code == current_ward_code]
  
  exceedance_prob_list[[i]] = tibble(
    Ward_id       = i,
    Ward_code     = current_ward_code,
    Ward_name     = current_ward_name,
    Percentile    = names(x_temp[[i]]),
    p_exceed_bham = exceedance_vec
  )
}

#------------------------------------------------------
#combine ward results
#------------------------------------------------------
exceedance_prob_df = bind_rows(exceedance_prob_list) %>%
  mutate(
    Percentile = factor(Percentile, levels = names(x_temp[[1]])),
    Ward_name  = factor(Ward_name, levels = rev(sort(unique(Ward_name))))
  )

write_rds(exceedance_prob_df, "output/04_Exceedance_prob_RR_greater_Birmingham.rds")

#------------------------------------------------------
#heatmap
#------------------------------------------------------
exceedance_heatmap = ggplot(exceedance_prob_df, aes(x = Percentile, y = Ward_name, fill = p_exceed_bham)) +
  geom_tile() +
  scale_fill_distiller(palette = "Greens", direction = 1, limits = c(0, 1), oob = scales::squish) +
  scale_x_discrete(
    breaks = c("0.0%", "1.0%", "10.0%", "25.0%", "50.0%", "75.0%", "90.0%", "99.0%", "100.0%")
  ) +
  labs(
    title    = paste0(
      "Posterior probability that ward relative risk exceeds ",
      "the Birmingham overall relative risk"
    ),
    subtitle = "Comparison at each ward-specific temperature percentile",
    x        = "Temperature percentile",
    y        = NULL,
    fill     = "Probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid  = element_blank()
  )

exceedance_heatmap

ggsave(
  filename = "figs/04_Exceedance_probability_ward_vs_Birmingham.png",
  plot     = exceedance_heatmap,
  width    = 12,
  height   = 14,
  units    = "in",
  dpi      = 600
)

#======================================================
#calculate excess mortality
#======================================================
# ward_id =1
# 
# calc_excess_bayes = function(
    #     df_complete,
#     cb_res,
#     ward_id,
#     lag_max = 21,
#     temp_bounds = c(0.01, 0.99)
# ) {
#   
#   # -----------------------------
#   # 1. Subset data for ward
#   # -----------------------------
#   df_w = df_complete[df_complete$new_id == ward_id, ]
#   
#   deaths = df_w$deaths
#   temp   = df_w$tasmean
#   year   = df_w$year
#   
#   X_day = as.matrix(df_w[, paste0("cb", 1:30)])   # daily CB matrix
#   
#   # -----------------------------
#   # 2. Temperature grid for MMT
#   # -----------------------------
#   t_lo = quantile(temp, temp_bounds[1], na.rm = TRUE)
#   t_hi = quantile(temp, temp_bounds[2], na.rm = TRUE)
#   
#   temp_grid = seq(t_lo, t_hi, by = 0.1)
#   
#   # build CB on grid 
#   cb_grid = crossbasis(
#     matrix(rep(temp_grid, 21 + 1), ncol = 21 + 1),
#     argvar = list(
#       fun   = "bs",
#       knots = quantile(temp, probs = c(0.1, 0.75, 0.9), na.rm = TRUE)
#     ),
#     arglag = list(
#       fun = "ns",
#       knots = logknots(lag_max, 3),
#       intercept = TRUE
#     )
#   )
#   
#   X_grid = as.matrix(cb_grid)
#   
#   # -----------------------------
#   # 3. Posterior draws
#   # -----------------------------
#   beta_mat = cb_res[[ward_id]]
#   n_draws  = nrow(beta_mat)
#   
#   years_u  = sort(unique(year))
#   
#   # storage (small!)
#   heat_year = matrix(NA, nrow = length(years_u), ncol = n_draws,
#                      dimnames = list(years_u, NULL))
#   cold_year = heat_year
#   
#   heat_total = numeric(n_draws)
#   cold_total = numeric(n_draws)
#   
#   # -----------------------------
#   # 4. Loop over posterior draws
#   # -----------------------------
#   for (b in seq_len(n_draws)) {
#     
#     beta_b = beta_mat[b, ]
#     
#     # daily log-risk
#     eta_day = drop(X_day %*% beta_b)
#     
#     # grid log-risk → MMT
#     eta_grid = drop(X_grid %*% beta_b)
#     
#     i_mmt    = which.min(eta_grid)  #reference mmt
#     eta_ref  = eta_grid[i_mmt]
#     mmt_b    = temp_grid[i_mmt]
#     
#     # daily RR
#     RR_day = exp(eta_day - eta_ref)
#     
#     # daily attributable deaths
#     attr_day = deaths * (RR_day - 1) / RR_day
#     
#     # heat / cold masks
#     is_heat = temp > mmt_b
#     is_cold = temp < mmt_b
#     
#     # yearly aggregation
#     heat_year[, b] = tapply(
#       attr_day[is_heat],
#       factor(year[is_heat], levels = years_u),
#       sum,
#       na.rm = TRUE
#     )
#     
#     cold_year[, b] = tapply(
#       attr_day[is_cold],
#       factor(year[is_cold], levels = years_u),
#       sum,
#       na.rm = TRUE
#     )
#     
#     
#     
#     # whole-period totals
#     heat_total[b] = sum(attr_day[is_heat], na.rm = TRUE)
#     cold_total[b] = sum(attr_day[is_cold], na.rm = TRUE)
#   }
#   
#   # -----------------------------
#   # 5. Summaries
#   # -----------------------------
#   summarise_draws = function(x) {
#     c(
#       median = median(x, na.rm = TRUE),
#       lo     = quantile(x, 0.025, na.rm = TRUE),
#       hi     = quantile(x, 0.975, na.rm = TRUE)
#     )
#   }
#   
#   summary = list(
#     yearly = list(
#       heat = t(apply(heat_year, 1, summarise_draws)),
#       cold = t(apply(cold_year, 1, summarise_draws))
#     ),
#     total = list(
#       heat = summarise_draws(heat_total),
#       cold = summarise_draws(cold_total)
#     )
#   )
#   
#   draws = list(
#     yearly = list(
#       heat = heat_year,
#       cold = cold_year
#     ),
#     total = list(
#       heat = heat_total,
#       cold = cold_total
#     )
#   )
#   
#   return(list(
#     draws   = draws,
#     summary = summary
#   ))
# }
# 
# 
# excess_check = calc_excess_bayes(df_complete = df_complete,
#                                  cb_res      = cb_res,
#                                  ward_id     = 20)
# 
# 
# 
# 
# 
# 
# 



# # Plot with MMT as the reference
# plot(
#   x_temp[[i]], RR_med_mmt, type = "l", log = "y",
#   ylim = range(RR_lo_mmt, RR_hi_mmt),
#   xlab = "Temperature", ylab = "Relative Risk (Ref = MMT)",
#   main = paste("RR centered at MMT:", round(mmt_value, 1))
# )
# lines(x_temp[[i]], RR_lo_mmt, lty = 2)
# lines(x_temp[[i]], RR_hi_mmt, lty = 2)
# abline(h = 1)
# abline(v = mmt_value, col = "red", lty = 3) # Mark the MMT




