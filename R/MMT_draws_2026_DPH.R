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
  
  # fitted 0-21 cross-basis
  cb_i = cb_pred[[i]]
  
  # temperature grid for this ward
  temp_grid_i = as.numeric(x_temp[[i]])
  
  # lag range used to define the short-term MMT
  mmt_lags = 0:3
  
  # evaluate the original fitted 0-21 basis
  # only at lags 0, 1, 2 and 3
  Xpred_03 = dlnm:::mkXpred(
    type    = "cb",
    basis   = cb_i,
    at      = temp_grid_i,
    predvar = temp_grid_i,
    predlag = mmt_lags,
    cen     = NULL
  )
  
  # sum lag 0-3 contributions for each temperature
  g = rep(
    seq_along(temp_grid_i),
    times = length(mmt_lags)
  )
  
  Xsum_03 = rowsum(
    Xpred_03,
    group = g,
    reorder = FALSE
  )
  
  # check compatibility with fitted coefficients
  stopifnot(
    ncol(Xsum_03) == ncol(beta_reg)
  )
  
  # cumulative log-RR over lags 0-3
  # temperature grid x posterior draws
  rr = Xsum_03 %*% t(beta_reg)
  
  # Search between P1 and P99
  allowed_idx = which(
    percentiles >= 0.01 &
      percentiles <= 0.99
  )
  
  # Find draw-specific MMT
  min_RR_position = apply(
    rr[allowed_idx, , drop = FALSE],
    2,
    function(x) {
      
      minimum_logRR = min(
        x,
        na.rm = TRUE )
      
      # All temperatures whose RR is within 1%
      # of the minimum RR
      near_minimum = which(
        x <= minimum_logRR + log(1.01)
      )
      # Select warmest point in minimum-risk plateau
      allowed_idx[max(near_minimum)]
    }
  )
  
  # Log-RR at each draw-specific MMT
  logRR_at_MMT = rr[
    cbind(
      min_RR_position,
      seq_len(ncol(rr))
    )
  ]
  
  
  # Centre each posterior draw at its own MMT
  store_specific_mmt = sweep(
    rr,
    MARGIN = 2,
    STATS = logRR_at_MMT,
    FUN = "-"
  )
  
  current_ward_code = unique(
    df_complete$ward22cd[
      df_complete$new_id == i
    ]
  )
  
  mmt_draws_by_ward[[i]] = tibble(
    nsim = seq_len(ncol(rr)),
    MMT = as.numeric(
      x_temp[[i]][min_RR_position]
    ),
    Ward_code = current_ward_code,
    Ward_id = i
  )
  
  rr_mmt_centered[[i]] = exp(
    store_specific_mmt
  )
  
}
mmt_diagnostic = map_dfr(
  mmt_draws_by_ward,
  function(x) {
    tibble(
      Ward_id = unique(x$Ward_id),
      MMT_min = min(x$MMT),
      MMT_p025 = quantile(x$MMT, 0.025),
      MMT_median = median(x$MMT),
      MMT_p975 = quantile(x$MMT, 0.975),
      MMT_max = max(x$MMT)
    )
  }
)

mmt_diagnostic
summary(mmt_diagnostic$MMT_median)

write_rds(mmt_draws_by_ward, "output/mmt_draws_by_ward_new_0_3.rds")

