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


###############################################################################
#===============================================================================
#plot all thje posterior draws for ward id 1 
#===============================================================================


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

################################################################################
#===============================================================================
#calculate ward sepcific mmt doe each posterior simulation for each ward
#===============================================================================
mmt_draws_by_ward = vector("list", 69)
rr_mmt_centered  = vector("list", 69)


for (i in 1:69){


# Extract the posterior samples of the coefficients (e.g., from INLA or MCMC).
# Dimensions: 1000 rows x 30 columns.
beta_reg = cb_res[[i]]

# Extract the cross-basis prediction matrix for this location.
# Dimensions: 119 rows x 30 columns 
cb_i = cb_pred[[i]]

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


# temperature grid used to compute RR
temp_grid <- x_temp[[i]]

# tocontraint temp bounds (change if needed)
t_bounds = quantile(x_temp[[i]], probs = c(0.01, 0.99), na.rm = TRUE)

# indices of grid points inside bounds
ok_idx = x_temp[[i]] >= t_bounds[1] &  x_temp[[i]] <= t_bounds[2]


#Get global indices once
allowed_idx = which(ok_idx)



#find the lowest RR for each simulated posterior from full_1000_rr

min_RR_position =  apply(full_1000_rr[allowed_idx, , drop = FALSE],2, function(x) allowed_idx[ which.min(x) ])


#initialise an empty dataframe to store each specific mmt below
store_specific_mmt = matrix(nrow=length(x_temp[[i]]), ncol=1000)
             

for(n in 1:1000){ #for every simulated posterior 

#apply the simulation-specific min rr to each col in the rr matrix to re-centre, after that, write in into the matrix
store_specific_mmt[,n]= rr[,n]-rr[,n][min_RR_position[n]]




}

mmt_draws_by_ward[[i]] = data.frame(nsim = seq_len(ncol(rr)),
                       MMT = as.numeric(x_temp[[i]][min_RR_position]),
                       Ward_code = ward_map$Ward_Code[ward_map$Ward_Code == unique(df_complete$ward22cd[df_complete$new_id == i])],
                       Ward_id = i
)

rr_mmt_centered[[i]] = exp(store_specific_mmt)

}


write_rds(mmt_draws_by_ward, "output/mmt_draws_by_ward.rds")


###############################################################################
###############################################################################
###############################################################################
#===============================================================================
#plot all the ward median exposure response curves
#===============================================================================

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
#################################################################################
#################################################################################
# MMT = mmt_draws_by_ward


mmt_uncertainty_plot_df = data.table::rbindlist(MMT) %>% 
  group_by(Ward_code,Ward_id) %>% 
  mutate(
    prob_gt_18 = round(mean(MMT > 18, na.rm = TRUE)*100,1),
    upperCrI = quantile(MMT, 0.975),
    lowerCrI = quantile(MMT, 0.025),
    label_pr18 = paste0("Pr(MMT>18) ", "\n= ", prob_gt_18, "%"),
    label_crI = paste0("MMT CrI ", "95%: ", "(",round(lowerCrI,1), ", ", round(upperCrI,1),")"),
    bin_boundary = floor(min(MMT, na.rm = TRUE)),
    ymax = max(
      c(hist(MMT,
            breaks = seq(floor(min(MMT)), ceiling(max(MMT)), by = 1),
            plot   = FALSE)$density,
        # kernel density heights
        density(MMT, adjust = 1)$y), na.rm = TRUE),
    .groups = "drop"
  )%>% 
  mutate(label_pr18 = ifelse(nsim ==1, label_pr18, NA),
         label_crI = ifelse(nsim ==1, label_crI, NA)) %>% 
  left_join(ward_map, by = c("Ward_code" = "Ward_Code"))

#--------------------------------------------------------

starts = seq(1, 69, by = 15)

# Loop through each batch
for (start_id in starts) {
  
  # Calculate the range for this batch
  end_id = min(start_id + 14, 69)
  current_ids = start_id:end_id
  
  # Filter data
  batch_data = subset(mmt_uncertainty_plot_df, Ward_id %in% current_ids)
  
  # Calculate dynamic height 
  # Count how many unique wards are in this batch
  n_wards = length(unique(batch_data$Ward_id))
  # Calculate rows needed (round up division by 3)
  n_rows = ceiling(n_wards / 3)
  

p =batch_data %>% 
  ggplot(aes(x=MMT))+
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "grey40", 
                 fill = "grey90",
                 binwidth = 1)+
  geom_density(lwd = 1,colour = "#578f9f", 
               fill = "#83b3c0",alpha = 0.5,adjust = 1)+
  geom_vline(xintercept = 18)+
  geom_segment(aes(x = 18, y = ymax*1.18, xend = 28, yend = ymax*1.18),
               arrow = arrow(type = "open", length = unit(0.1, "cm")),
               linewidth = 0.15)+
  geom_text(aes(x = 19, y = ymax*1.3, label = label_pr18),
            hjust = 0, size = 2.5, inherit.aes = FALSE,fontface = "bold")  +
  geom_text(aes(x = -5, y = ymax*1.3, label = label_crI),
            hjust = 0, size = 2.5, inherit.aes = FALSE)  +
  scale_x_continuous(breaks = seq(-5, 30, 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
  coord_cartesian(xlim = c(-5, 30))+
  labs(x="Temperature °C",
       y="Density posterior distribution MMT")+
  theme_classic(base_size = 11)+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "none")+
  facet_wrap(~ Ward_Name, ncol = 3, scales = 'free') 

# --- STEP 3: Save with Dynamic Height ---
# If a full page (5 rows) is 297mm, then one row is approx 59.4mm.
# We calculate the height required for the current number of rows.
dynamic_height = 297 * (n_rows / 5)

filename = paste0("figs/","MMT_uncertainty_ward_plots_", start_id, "_to_", end_id, ".pdf")

ggsave(filename, plot = p, width = 210, height = dynamic_height, units = "mm")

message(paste("Saved:", filename, "(Height:", round(dynamic_height, 1), "mm)"))
}


# merged the files (optional)
pdf_combine(
  input = c("figs/MMT_uncertainty_ward_plots_1_to_15.pdf",
            "figs/MMT_uncertainty_ward_plots_16_to_30.pdf",
            "figs/MMT_uncertainty_ward_plots_31_to_45.pdf",
            "figs/MMT_uncertainty_ward_plots_46_to_60.pdf",
            "figs/MMT_uncertainty_ward_plots_61_to_69.pdf"),
  output = "figs/Combined_MMT_uncertainty_wards.pdf"
)


#--------------------------------------------------------
#========================================================
#exceedance prob for RR at P99/P1 >1
#--------------------------------------------------------




RR_Prob_plot_list = list()
for(i in 1:69){

current_ward_code = unique(df_complete$ward22cd[df_complete$new_id == i])
current_ward_name = ward_map$Ward_Name[ward_map$Ward_Code == current_ward_code]
  
  
  
idx_99 =which(names(x_temp[[i]]) == "99.0%")
idx_1 = which(names(x_temp[[i]]) == "1.0%")



RR_draws_P1  = rr_mmt_centered[[i]][idx_1, ] 
RR_draws_P99  = rr_mmt_centered[[i]][idx_99, ]



plot_df = data.frame(
  RR_draws_P1 = mean(RR_draws_P1>1),
  RR_draws_P99 = mean(RR_draws_P99 >1),
  Ward_id = i,
  Ward_code = current_ward_code,
  Ward_name = current_ward_name) 


RR_Prob_plot_list[[i]] = plot_df

}

RR_prob_df <- dplyr::bind_rows(RR_Prob_plot_list)



write_rds(RR_prob_df, "output/RR_exceedance_prob.rds")




plot_sf = ward_map %>%
  left_join(RR_prob_df, by = c("Ward_Code" = "Ward_code")) %>%
  mutate(
    RR_draws_P99 = as.numeric(RR_draws_P99),
    RR_draws_P99_grey90 = if_else(RR_draws_P99 < 0.90, NA_real_, RR_draws_P99),
    
    evidence_99th = case_when(
      is.na(RR_draws_P99)  ~ "Missing",
      RR_draws_P99 >= 0.95 ~ "Strong evidence (≥0.95)",
      RR_draws_P99 >= 0.90 ~ "Some evidence (0.90–0.95)",
      TRUE                 ~ "No evidence (<0.90)"
    ),
    
    tooltip_99th = paste0(
      "<B>Pr(RR at P99)>1:</B> ", round(RR_draws_P99, 2),
      "\n<B>Inference:</B> ", evidence_99th
    ),
    RR_draws_P1 = as.numeric(RR_draws_P1),
    RR_draws_P1_grey90 = if_else(RR_draws_P1 < 0.90, NA_real_, RR_draws_P1),
    evidence_P1 = case_when(
      is.na(RR_draws_P1)  ~ "Missing",
      RR_draws_P1 >= 0.95 ~ "Strong evidence (≥0.95)",
      RR_draws_P1 >= 0.90 ~ "Some evidence (0.90–0.95)",
      TRUE                 ~ "No evidence (<0.90)"
    ),
    tooltip_P1 = paste0(
      "<B>Pr(RR at P1)>1:</B> ", round(RR_draws_P1, 3),
      "\n<B>Inference:</B> ", evidence_P1
    ))
    
    



#custom css
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




figure_1_Prob = ggplot(plot_sf) +
  geom_sf_interactive(aes(fill = RR_draws_P99_grey90, data_id = Ward_Code, tooltip = tooltip_99th),color = "white") +
  scale_fill_distiller(
    palette  = "Reds",
    direction = 1,
    na.value = "grey85",
    limits   = c(0.90, 1.00),                 # keep the legend focused on 0.90–1.00
    breaks   = c(0.90, 0.95, 1.00)
  ) +
  labs(fill = "Probability") +
  ggtitle("Exceedance probability that relative risk exceeds 1 \nat the 99th percentile temperature")+
  theme_void(base_size = 16) +
  theme()




figure_2_Prob = ggplot(plot_sf) +
  geom_sf_interactive(aes(fill = RR_draws_P1_grey90, data_id = Ward_Code, tooltip = tooltip_P1),color = "white") +
  scale_fill_distiller(
    palette  = "Blues",
    direction = 1,
    na.value = "grey85",
    limits   = c(0.90, 1.00),                 # keep the legend focused on 0.90–1.00
    breaks   = c(0.90, 0.95, 1.00)
  ) +
  labs(fill = "Probability") +
  ggtitle("Exceedance probability that relative risk exceeds 1 \nat the 99th percentile temperature")+
  theme_void(base_size = 16) +
  theme()




ggiraph::girafe(
  ggobj = figure_1_Prob,
  width_svg = 14,
  height_svg = 10,
  options = list(
    opts_hover(css = hover_css),
    opts_tooltip(css = tooltip_css),
    opts_hover_inv(css = "")))




ggiraph::girafe(
  ggobj = figure_2_Prob ,
  width_svg = 14,
  height_svg = 10,
  options = list(
    opts_hover(css = hover_css),
    opts_tooltip(css = tooltip_css),
    opts_hover_inv(css = "")))

#========================================================
#to obtain the heatmap for RR and get the birmingham overall RR heat and cold  
mmt_draws_by_ward = vector("list", 69)
log_rr_mmt_centered  = vector("list", 69)


for (i in 1:69){
  
  
  # Extract the posterior samples of the coefficients (e.g., from INLA or MCMC).
  # Dimensions: 1000 rows x 30 columns.
  beta_reg = cb_res[[i]]
  
  # Extract the cross-basis prediction matrix for this location.
  # Dimensions: 119 rows x 30 columns 
  cb_i = cb_pred[[i]]
  
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
  
  
  # temperature grid used to compute RR
  temp_grid <- x_temp[[i]]
  
  # tocontraint temp bounds (change if needed)
  t_bounds = quantile(x_temp[[i]], probs = c(0.01, 0.99), na.rm = TRUE)
  
  # indices of grid points inside bounds
  ok_idx = x_temp[[i]] >= t_bounds[1] &  x_temp[[i]] <= t_bounds[2]
  
  
  #Get global indices once
  allowed_idx = which(ok_idx)
  
  
  
  #find the lowest RR for each simulated posterior from full_1000_rr
  
  min_RR_position =  apply(full_1000_rr[allowed_idx, , drop = FALSE],2, function(x) allowed_idx[ which.min(x) ])
  
  
  #initialise an empty dataframe to store each specific mmt below
  store_specific_mmt = matrix(nrow=length(x_temp[[i]]), ncol=1000)
  
  
  for(n in 1:1000){ #for every simulated posterior 
    
    #apply the simulation-specific min rr to each col in the rr matrix to re-centre, after that, write in into the matrix
    store_specific_mmt[,n]= rr[,n]-rr[,n][min_RR_position[n]]
    
    
    
    
  }
  
  mmt_draws_by_ward[[i]] = data.frame(nsim = seq_len(ncol(rr)),
                                      MMT = as.numeric(x_temp[[i]][min_RR_position]),
                                      Ward_code = ward_map$Ward_Code[ward_map$Ward_Code == unique(df_complete$ward22cd[df_complete$new_id == i])],
                                      Ward_id = i
  )
  
  log_rr_mmt_centered[[i]] = store_specific_mmt
  
}

#-------------------------------------------------------------

n_wards = 69

# Use log scale for aggregation (stable), then exp back to RR
# Step 1: ward-level posterior medians at each percentile (logRR)
ward_logRR_med_mat = sapply(1:n_wards, function(i){
  apply(log_rr_mmt_centered[[i]], 1,  function (x) mean(x))  # length nTemp
})
ward_logRR_med_mat = t(ward_logRR_med_mat)  # wards x nTemp

# Step 2: Birmingham "overall" at each percentile across wards (median + CrI)
bham_logRR_mean = apply(ward_logRR_med_mat, 2, mean, na.rm = TRUE)
bham_logRR_lo  = apply(ward_logRR_med_mat, 2, quantile, 0.025, na.rm = TRUE)
bham_logRR_hi  = apply(ward_logRR_med_mat, 2, quantile, 0.975, na.rm = TRUE)

bham_RR_by_pct = data.frame(
  Percentile = names(x_temp[[1]]),          
  RR_mean = exp(bham_logRR_mean),
  RR_lo  = exp(bham_logRR_lo),
  RR_hi  = exp(bham_logRR_hi)
)

#--------------------------------------------
exceedance_prob_list = vector("list", 69)

nrow(bham_RR_by_pct)

for(i in 1:69){

ward_rr_mmt_centered = rr_mmt_centered[[i]]
  

exceedance_vec = sapply(1:119, function(n){
  mean(ward_rr_mmt_centered[n,] >bham_RR_by_pct$RR_mean[[n]])})  

exceedance_prob_list[[i]] = data.frame(
  Ward_id = i,
  Ward_code = unique(df_complete$ward22cd[df_complete$new_id == i]),
  Ward_name = ward_map$Ward_Name[ward_map$Ward_Code == unique(df_complete$ward22cd[df_complete$new_id == i])],
  Percentile = bham_RR_by_pct$Percentile,
  p_exceed_bham = exceedance_vec
  
)

}


#-------------------------------------------------------------------------

exceedance_prob_df <- dplyr::bind_rows(exceedance_prob_list) %>%
  mutate(
    Percentile = factor(Percentile, levels = unique(Percentile)),
    Ward_name  = factor(Ward_name, levels = rev(sort(unique(Ward_name))))
  )




write_rds(exceedance_prob_df, "output/Exceedance_prob_RR_greater_mean.rds")






ggplot(exceedance_prob_df, aes(x = Percentile, y = Ward_name, fill = p_exceed_bham)) +
  geom_tile() +
  scale_fill_distiller(
    palette = "Greens",
    direction = 1,
    limits = c(0, 1),
    oob = squish
  ) +
  scale_x_discrete(breaks = c("0.0%", "1.0%","10.0%","25.0%","50.0%","75.0%","90.0%","99.0%" ,"100.0%")) +
  labs(
    title = "Exceedance probability of ward RR is higher than the mean RR of Brimingham at \nevery temperature percentile",
    x = "Temperature percentile",
    y = NULL,
    fill = "Probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )


#==============================================================================





































