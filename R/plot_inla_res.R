library(tidyverse)
library(dlnm)
library("sf")
library(qpdf)
library(tmap)

#load the data
df_complete = read_rds("output/df_complete_for_inla.rds")
cb_res = read_rds("output/predicted_inla_spatial_casecrossover.rds")



# Load Ward shp file
ward_map = read_sf("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp")


# Create a sequence of probabilities to capture the distribution of temperature.
# This specific sequence has high resolution at the tails (0-1% and 99-100%)
percentiles = c(
  seq(0, 1, by = 0.1),    # Extreme lows (0.0% to 1.0%)
  seq(2, 98, by = 1),     # Core distribution (2.0% to 98.0%)
  seq(99, 100, by = 0.1)  # Extreme highs (99.0% to 100.0%)
) / 100

# Initialize a list to store the temperature distribution for each of the 69 wards.
x_temp = vector("list", 69)

for (i in 1:69) {
  # Subset data for the current ward ID
  df_w = df_complete[df_complete$new_id == i, ]
  
  # Calculate the actual temperature values corresponding to the custom percentiles
  # defined above for this specific location.
  x_temp[[i]] = quantile(df_w$tasmean, percentiles, na.rm = TRUE)
}

# Generate Cross-Basis Matrices for Prediction
# Initialize a list to store the DLNM cross-basis objects.
cb_pred = vector("list", 69)

for (i in 1:69) {
  
  # Subset data again to access ward-specific knots
  df_w = df_complete[df_complete$new_id == i, ]
  
  # Create the Cross-Basis matrix (interaction between exposure and lag space)
  cb_pred[[i]] = crossbasis(
    
    # PREDICTION MATRIX:
    # Repeat the temperature quantiles (x_temp) across all lag periods.
    # We use '21 + 1' to account for lags 0 through 21.
    matrix(
      rep(x_temp[[i]], 21 + 1),
      ncol = 21 + 1
    ),
    argvar = list(
      fun = "bs",
      knots = quantile(df_w$tasmean, probs = c(0.1, 0.75, 0.9), na.rm = TRUE)
    ),
    arglag = list(
      fun = "ns",
      knots = logknots(21, 3),
      intercept = TRUE
    )
  )
}


######################################################################
#plot all thje posterior draws for ward id 1 

cb_res[[1]]


# Extract the posterior samples of the coefficients (e.g., from INLA or MCMC).
# Dimensions: 1000 rows x 30 columns.
beta_reg = cb_res[[1]]

# Extract the cross-basis prediction matrix for this location.
# Dimensions: 119 rows x 30 columns 
cb_i = cb_pred[[1]]

#

rr = cb_i %*% t(beta_reg)

      
# 3. Define Reference Temperature (Centering)
# Find the index in the 'percentiles' vector corresponding to the 90th percentile.
# This will be our reference point (RR = 1).
i_cen = which(percentiles == 0.9)

# Center the predictions.
# For each simulation (column of rr), subtract the value at the reference index.
# This ensures that at the 30th percentile, the log-RR is always 0.
rr_cen = apply(rr, 2, function(x) x - x[i_cen])

full_1000_rr = exp(rr_cen)

# Median Estimate
RR_med = apply(exp(rr_cen), 1, median)

# Lower Bound (2.5%)
RR_lo  = apply(exp(rr_cen), 1, quantile, probs = 0.025, na.rm = TRUE)

# Upper Bound (97.5%)
RR_hi  = apply(exp(rr_cen), 1, quantile, probs = 0.975, na.rm = TRUE)

##########
full_ward_1_plot_df = as.data.frame(full_1000_rr)
colnames(full_ward_1_plot_df) = gsub("V","RR_", colnames(full_ward_1_plot_df))

##########
individuals_line_plot = data.frame(
  Percentage = names(x_temp[[i]]),
  Temp = as.numeric(x_temp[[i]]),
  RR_median = RR_med,
  RR_LL = RR_lo,
  RR_UL = RR_hi )
  
ward_map$Ward_Name[ ward_map$Ward_Code ==  unique(df_complete$ward22cd[df_complete$new_id==1])]


posterior_draws_for_acocks_green = cbind(data.frame(Percentage = names(x_temp[[i]]),
                 Temp = as.numeric(x_temp[[i]])),
      full_ward_1_plot_df) %>% 
  pivot_longer(
    cols = c(-Temp,-Percentage),        
    names_to = "Sim",    
    values_to = "RR"     
  ) %>% 
  ggplot(aes(x = Temp, y = RR, group = Sim)) +
  
  # 1. The 1000 Grey Lines
  geom_line(colour = "grey90", size=0.6) +
  
  # 2. Reference Lines
  geom_hline(aes(yintercept = 1), colour = "black") +
  geom_vline(aes(xintercept = as.numeric(x_temp[[1]][[100]])), 
             linetype = "dotted", size = 0.2, colour = "black") +
  
  # 3. Median and CIs
  geom_line(data = individuals_line_plot, aes(x = Temp, y = RR_median), 
            inherit.aes = FALSE, size = 0.8) +
  geom_line(data = individuals_line_plot, aes(x = Temp, y = RR_LL), 
            inherit.aes = FALSE, size = 0.8, linetype = "dashed") +
  geom_line(data = individuals_line_plot, aes(x = Temp, y = RR_UL), 
            inherit.aes = FALSE, size = 0.8, linetype = "dashed") +
  
  # 4. THE FIX: Use annotate() for the text label
  annotate("text", 
           x = as.numeric(x_temp[[1]][[100]]), 
           y = 3, 
           label = paste0("Centering         \ntemperature = ", round(as.numeric(x_temp[[1]][[100]]), 2)),
           angle = 0, vjust = 0, hjust = 1.2, size = 3, colour = "black") +
  
  # Scales and Theme
  scale_y_continuous(breaks = 0:6) +
  scale_x_continuous(breaks = seq(-5, 30, 5)) +
  coord_cartesian(ylim = c(0, 6), xlim = c(-5, 30)) +
  theme_classic(base_size = 11) +
  labs(x = "Temperature °C", y = "Relative Risk") +
  ggtitle(paste0("Posterior draws of exposure-response curves for ", 
                 ward_map$Ward_Name[ward_map$Ward_Code == unique(df_complete$ward22cd[df_complete$new_id == 1])])) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))


ggsave("figs/posterior_draws_for_acocks_green.jpg",posterior_draws_for_acocks_green, width=6, height=4, dpi =600,units = "in")



######################################################################
RR_MMT_plot_list = list()


for (w in 1:69){


# 1. Select location index.
i = w
current_ward_code = unique(df_complete$ward22cd[df_complete$new_id == i])
current_ward_name = ward_map$Ward_Name[ward_map$Ward_Code == current_ward_code]


# Extract the posterior samples of the coefficients fron inla results.
# Dimensions: 1000 rows x 30 columns.
beta_reg = cb_res[[i]]

# Extract the cross-basis prediction matrix for this location.
# Dimensions: 119 rows x 30 columns 
cb_i = cb_pred[[i]]

#we want to multiply the above matrices to obtain 119 rows x 1000 columns
#we transpose beta_reg to 30 rows * 1000 colums, such that align perfectly with the cb_i
#resulting a 119 rows x 1000 matrix 
rr = cb_i %*% t(beta_reg)


# 3. Define Reference Temperature (Centering)
# Find the index in the 'percentiles' vector corresponding to the 90th percentile.
# This will be our reference point (RR = 1).
i_cen = which(percentiles == 0.9)

# Center the predictions.
# For each simulation (column of rr), subtract the value at the reference index.
# This ensures that at the 30th percentile, the log-RR is always 0.
rr_cen = apply(rr, 2, function(x) x - x[i_cen])


# 4. Calculate Final Summary Statistics
# Convert log-RR back to the normal scale using exp() and calculate
# the median and 95% Credible Intervals (2.5% and 97.5%).

# Median Estimate
RR_med = apply(exp(rr_cen), 1, median)

# Lower Bound (2.5%)
RR_lo  = apply(exp(rr_cen), 1, quantile, probs = 0.025, na.rm = TRUE)

# Upper Bound (97.5%)
RR_hi  = apply(exp(rr_cen), 1, quantile, probs = 0.975, na.rm = TRUE)

#############################################################
# Find the index (position) where the Relative Risk is lowest
i_mmt = which.min(RR_med)

# Extract the actual MMT value from your temperature vector
mmt_value = x_temp[[i]][i_mmt]

# (Optional) Check what the Risk is at that point (should be the minimum)
min_risk = RR_med[i_mmt]


#Use the MMT index we just found as the new centering point
rr_cen_mmt = apply(rr, 2, function(x) x - x[i_mmt])


# Recalculate Median and CIs based on this new center
RR_med_mmt = apply(exp(rr_cen_mmt), 1, median)
RR_lo_mmt  = apply(exp(rr_cen_mmt), 1, quantile, 0.025, na.rm = TRUE)
RR_hi_mmt  = apply(exp(rr_cen_mmt), 1, quantile, 0.975, na.rm = TRUE)


#########################################

plot_df = data.frame(
  Percentage = names(x_temp[[i]]),
  Temp = as.numeric(x_temp[[i]]),
  RR_median = RR_med_mmt,
  RR_LL = RR_lo_mmt,
  RR_UL = RR_hi_mmt,
  MMT = as.numeric(mmt_value),
  Ward_id = i,
  Ward_code = current_ward_code,
  Ward_name = current_ward_name) %>% 
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



#########################################
#  Ensure your data is ready
real_plot_df = data.table::rbindlist(RR_MMT_plot_list)

write_rds(real_plot_df, file ="output/RR_MMT_plot_data.rds")
#  Define your batches (Start index for every group of 15)
# This sequence will be: 1, 16, 31, 46, 61
starts = seq(1, 69, by = 15)

# Loop through each batch
for (start_id in starts) {
  
  # Calculate the range for this batch
  end_id = min(start_id + 14, 69)
  current_ids = start_id:end_id
  
  # Filter data
  batch_data = subset(real_plot_df, Ward_id %in% current_ids)
  
  # Calculate dynamic height 
  # Count how many unique wards are in this batch
  n_wards = length(unique(batch_data$Ward_id))
  # Calculate rows needed (round up division by 3)
  n_rows = ceiling(n_wards / 3)
  
  # Create the plot
  p = ggplot(batch_data, aes(Temp, y = RR_median)) +
    geom_ribbon(aes(ymin = RR_LL, ymax = RR_UL), alpha = 0.2, fill = "#A092AA", colour = NA) +
    geom_line(size = 0.8, colour = "#A092AA") +
    geom_hline(yintercept = 1) +
    geom_vline(data = batch_data, aes(xintercept = MMT), linetype = "dashed", size = 0.5, alpha = 0.5) +
    
    geom_point(data = batch_data, aes(x = -4.8, y = ifelse(Percentage == "1.0%", 5.78, 5.18), shape = Percentage),
               size = 2, fill = "#7F6588", colour = "black", stroke = 0.8, inherit.aes = FALSE) +
    
    geom_point(data = batch_data, aes(shape = Percentage),
               size = 3, fill = "#7F6588", colour = "black", stroke = 0.8) +
    
    geom_text(data = batch_data, aes(x = -4.0, y = ifelse(Percentage == "1.0%", 5.8, 5.2), label = legend_label),
              hjust = 0, size = 3, inherit.aes = FALSE, fontface = "bold") +
    
    geom_text(data = batch_data %>% group_by(Ward_name) %>% summarise(MMT = mean(MMT)),
              aes(x = MMT, y = 3, label = paste0("MMT = ", sprintf("%.1f", MMT))),
              angle = 0, vjust = 0, hjust = 1.2, size = 3, colour = "black", inherit.aes = FALSE) +
    
    # We let ggplot decide the layout naturally now that we control the file size
    facet_wrap(~ Ward_name, ncol = 3, scales = 'free') +
    
    scale_y_continuous(breaks = 0:6) +
    scale_x_continuous(breaks = seq(-5, 30, 5)) +
    coord_cartesian(ylim = c(0, 6), xlim = c(-5, 30)) +
    scale_shape_manual(values = c("1.0%" = 21, "99.0%" = 24)) +
    
    theme_classic(base_size = 11) +
    labs(x = "Temperature °C", y = "Relative Risk") +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          legend.position = "none") 
  
  # --- STEP 3: Save with Dynamic Height ---
  # If a full page (5 rows) is 297mm, then one row is approx 59.4mm.
  # We calculate the height required for the current number of rows.
  dynamic_height = 297 * (n_rows / 5)
  
  filename = paste0("figs/","Overall_cumul_RR_ward_plots_", start_id, "_to_", end_id, ".pdf")
  
  ggsave(filename, plot = p, width = 210, height = dynamic_height, units = "mm")
  
  message(paste("Saved:", filename, "(Height:", round(dynamic_height, 1), "mm)"))
}

# merged the files (optional)
pdf_combine(
  input = c("figs/Overall_cumul_RR_ward_plots_1_to_15.pdf",
            "figs/Overall_cumul_RR_ward_plots_16_to_30.pdf",
            "figs/Overall_cumul_RR_ward_plots_31_to_45.pdf",
            "figs/Overall_cumul_RR_ward_plots_46_to_60.pdf",
            "figs/Overall_cumul_RR_ward_plots_61_to_69.pdf"),
  output = "figs/Combined_overall_cumul_RR_wards.pdf"
)



#################################################################################


real_plot_df %>% 
  filter(Ward_id == 1)










#################################################################################
#plot the mmt on the median of the ensemble response curve 

mmt_plot = real_plot_df%>% 
  group_by(Ward_name, Ward_code) %>% 
  summarise(MMT = mean(MMT,
                       .group = "drop"))


tmap::tmap_mode("plot")



mmt_plot = tm_shape(ward_map %>% left_join(mmt_plot, by = c("Ward_Code" = "Ward_code" )))+
  tm_polygons(
    fill= "MMT",
    fill.scale = tm_scale_continuous(values = "brewer.RdPu"),
    fill.legend = tm_legend("MMT", group_id = "top", frame = FALSE,bg.alpha=0)
  )+
  tm_layout(
    title.size = 1.2,
    legend.position = c(0.02, 0.92),
    frame = FALSE,
    inner.margins = c(0.07, 0, 0.15, 0) # Increased bottom margin (3rd number) to make room for caption
  ) +
  tm_title("Minimum Mortality Temperature")+
  tm_credits(text = "Minimum Mortality Temperature of the ward temperature distribution \nat the median exposure response curve.", position = c("LEFT","TOP"))+
  tm_compass(type = "8star",
             size = 4,
             position = c("RIGHT", "bottom"),
             color.light = "white")+
  tm_credits(
    text = paste("Contains OS data \u00A9 Crown copyright and database right",
                 # Get current year
                 format(Sys.Date(), "%Y"),
                 ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."),
    position = c("LEFT", "BOTTOM")
  )




tmap_save(mmt_plot ,filename = "figs/mmt_plot.png", height = 7,width =6, unit="in",dpi = 600 )











#################################################################################
#plot the RR of the posterior median exposure response curve at 99th percentile of temp
#################################################################################
RR_99th_map = real_plot_df %>% 
  filter(Percentage == "99.0%")
  


RR_99th_map = tm_shape(ward_map %>% left_join(RR_99th_map, by = c("Ward_Code" = "Ward_code" )))+
  tm_polygons(
    fill= "RR_median",
    fill.scale = tm_scale_continuous(values = "reds"),
    fill.legend = tm_legend("Relative Risk", group_id = "top", frame = FALSE,bg.alpha=0)
  )+
  tm_layout(
    title.size = 1.2,
    legend.position = c(0.02, 0.92),
    frame = FALSE,
    inner.margins = c(0.07, 0, 0.15, 0) # Increased bottom margin (3rd number) to make room for caption
  ) +
  tm_title("Relative Risk at 99th Percentile Temperature")+
  tm_credits(text = "Median Relative Risk of death at the 99th percentile of the ward temperature distribution \ncompared to the risk at the MMT of the median exposure response curve.", position = c("LEFT","TOP"))+
  tm_compass(type = "8star",
             size = 4,
             position = c("RIGHT", "bottom"),
             color.light = "white")+
  tm_credits(
    text = paste("Contains OS data \u00A9 Crown copyright and database right",
                 # Get current year
                 format(Sys.Date(), "%Y"),
                 ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."),
    position = c("LEFT", "BOTTOM")
  )
  


tmap_save( RR_99th_map,filename = "figs/RR_99th_map.png", height = 7,width =6, unit="in",dpi = 600 )


#################################################################################
#plot the RR of the posterior median exposure response curve at 1th percentile of temp
#################################################################################

RR_1st_map = real_plot_df %>% 
  filter(Percentage == "1.0%")


tmap::tmap_mode("plot")

RR_1st_map = tm_shape(ward_map %>% left_join(RR_1st_map, by = c("Ward_Code" = "Ward_code" )))+
  tm_polygons(
    fill= "RR_median",
    fill.scale = tm_scale_continuous(values = "blues"),
    fill.legend = tm_legend("Relative Risk", group_id = "top", frame = FALSE,bg.alpha=0)
  )+
  tm_layout(
    title.size = 1.2,
    legend.position = c(0.02, 0.92),
    frame = FALSE,
    inner.margins = c(0.07, 0, 0.15, 0) # Increased bottom margin (3rd number) to make room for caption
  ) +
  tm_title("Relative Risk at 1st Percentile Temperature")+
  tm_credits(text = "Median Relative Risk of death at the 1st percentile of the ward temperature distribution \ncompared to the risk at the MMT of the median exposure response curve.", position = c("LEFT","TOP"))+
  tm_compass(type = "8star",
             size = 4,
             position = c("RIGHT", "bottom"),
             color.light = "white")+
  tm_credits(
    text = paste("Contains OS data \u00A9 Crown copyright and database right",
                 # Get current year
                 format(Sys.Date(), "%Y"),
                 ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."),
    position = c("LEFT", "BOTTOM")
  )


tmap_save(RR_1st_map,filename = "figs/RR_1st_map.png", height = 7,width =6, unit="in",dpi = 600 )

#merge the plots

merged_1st_99th_RR_plot = tmap::tmap_arrange(RR_1st_map ,RR_99th_map )

tmap_save(merged_1st_99th_RR_plot, filename ="figs/merged_1st_99th_RR_plot.png", height = 7,width =12, unit="in",dpi = 600)

#################################################################################
#calculate excess mortality \
ward_id =1

calc_excess_bayes = function(
    df_complete,
    cb_res,
    ward_id,
    temp_bounds = c(0.01, 0.99)
) {
  
  # -----------------------------
  # 1. Subset data for ward
  # -----------------------------
  df_w = df_complete[df_complete$new_id == ward_id, ]
  
  deaths = df_w$deaths
  temp   = df_w$tasmean
  year   = df_w$year
  
  X_day = as.matrix(df_w[, paste0("cb", 1:30)])   # daily CB matrix
  
  # -----------------------------
  # 2. Temperature grid for MMT
  # -----------------------------
  t_lo = quantile(temp, temp_bounds[1], na.rm = TRUE)
  t_hi = quantile(temp, temp_bounds[2], na.rm = TRUE)
  
  temp_grid = seq(t_lo, t_hi, by = 0.1)
  
  # build CB on grid 
  cb_grid = crossbasis(
    matrix(rep(temp_grid, 21 + 1), ncol = 21 + 1),
    argvar = list(
      fun   = "bs",
      knots = quantile(temp, probs = c(0.1, 0.75, 0.9), na.rm = TRUE)
    ),
    arglag = list(
      fun = "ns",
      knots = logknots(lag_max, 3),
      intercept = TRUE
    )
  )
  
  X_grid = as.matrix(cb_grid)
  
  # -----------------------------
  # 3. Posterior draws
  # -----------------------------
  beta_mat = cb_res[[ward_id]]
  n_draws  = nrow(beta_mat)
  
  years_u  = sort(unique(year))
  
  # storage (small!)
  heat_year = matrix(NA, nrow = length(years_u), ncol = n_draws,
                      dimnames = list(years_u, NULL))
  cold_year = heat_year
  
  heat_total = numeric(n_draws)
  cold_total = numeric(n_draws)
  
  # -----------------------------
  # 4. Loop over posterior draws
  # -----------------------------
  for (b in seq_len(n_draws)) {
    
    beta_b = beta_mat[b, ]
    
    # daily log-risk
    eta_day = drop(X_day %*% beta_b)
    
    # grid log-risk → MMT
    eta_grid = drop(X_grid %*% beta_b)
    
    i_mmt    = which.min(eta_grid)  #reference mmt
    eta_ref  = eta_grid[i_mmt]
    mmt_b    = temp_grid[i_mmt]
    
    # daily RR
    RR_day = exp(eta_day - eta_ref)
    
    # daily attributable deaths
    attr_day = deaths * (RR_day - 1) / RR_day
    
    # heat / cold masks
    is_heat = temp > mmt_b
    is_cold = temp < mmt_b
    
    # yearly aggregation
    heat_year[, b] = tapply(
      attr_day[is_heat],
      factor(year[is_heat], levels = years_u),
      sum,
      na.rm = TRUE
    )
    
    cold_year[, b] = tapply(
      attr_day[is_cold],
      factor(year[is_cold], levels = years_u),
      sum,
      na.rm = TRUE
    )
    

    
    # whole-period totals
    heat_total[b] = sum(attr_day[is_heat], na.rm = TRUE)
    cold_total[b] = sum(attr_day[is_cold], na.rm = TRUE)
  }
  
  # -----------------------------
  # 5. Summaries
  # -----------------------------
  summarise_draws = function(x) {
    c(
      median = median(x, na.rm = TRUE),
      lo     = quantile(x, 0.025, na.rm = TRUE),
      hi     = quantile(x, 0.975, na.rm = TRUE)
    )
  }
  
  summary = list(
    yearly = list(
      heat = t(apply(heat_year, 1, summarise_draws)),
      cold = t(apply(cold_year, 1, summarise_draws))
    ),
    total = list(
      heat = summarise_draws(heat_total),
      cold = summarise_draws(cold_total)
    )
  )
  
  draws = list(
    yearly = list(
      heat = heat_year,
      cold = cold_year
    ),
    total = list(
      heat = heat_total,
      cold = cold_total
    )
  )
  
  return(list(
    draws   = draws,
    summary = summary
  ))
}


excess_check = calc_excess_bayes(df_complete = df_complete,
                  cb_res      = cb_res,
                  ward_id     = 20)



























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



###################################################################
# rr is [n_temp × n_draws]
# transpose so each column = one posterior draw (if needed)
# Here rr rows = temp, cols = draws → good

i=24

temp_bounds = quantile(x_temp[[i]], probs = c(0, 1), na.rm = TRUE)

# Iterate through every simulation (column of rr)
mmt_draws_constrained = apply(rr, 2, function(rr_draw) {
  
  # A. Create a "fair" uniform grid based on the current ward's range
  # Spaced every 0.1°C (or 0.05°C for higher precision)
  fair_grid = seq(min(x_temp[[i]]), max(x_temp[[i]]), by = 0.1)
  
  # B. Interpolate the biased 'rr_draw' onto this new fair grid
  # approx() connects the dots linearly; spline() is smoother but approx is usually sufficient
  rr_fair = approx(x = x_temp[[i]], y = rr_draw, xout = fair_grid)$y
  
  # C. Apply the valid boundary constraints to the FAIR grid
  # We only look for the minimum inside the safe 1-99% zone
  ok = fair_grid >= temp_bounds[1] & fair_grid <= temp_bounds[2]
  
  # D. Find the minimum index within the valid, fair subset
  valid_fair_grid = fair_grid[ok]
  valid_rr_fair   = rr_fair[ok]
  
  i_min = which.min(valid_rr_fair)
  
  # E. Return the temperature at that minimum
  return(valid_fair_grid[i_min])
})

# 3. Calculate Stats for Plotting
mmt_median = median(mmt_draws_constrained)
mmt_cri    = coda::HPDinterval(coda::as.mcmc(mmt_draws_constrained), prob = 0.90)



data.frame(mmt_draws_constrained ) %>% 
  ggplot(aes(x=mmt_draws_constrained ))+
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "grey40", 
                 fill = "grey90")+
  geom_density(lwd = 1,colour = "#578f9f", 
               fill = "#83b3c0",, alpha = 0.5)+
  geom_vline(xintercept=mmt_median,size=1,colour = "black",lwd = 1.2)+
  scale_x_continuous(breaks = seq(-5, 30, 5)) +
  coord_cartesian(xlim = c(-5, 30)) +
  
  theme_classic()


