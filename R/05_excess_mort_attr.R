setwd("D:/R project D/Extreme heat and cold")

library(readr)
library(dlnm)
library(tidyverse)
library(sf)
library(tmap)
library(readxl)
library(doParallel)
library(foreach)

#======================================================
# load data
#======================================================
df_complete = read_rds("output/03_df_complete_for_inla.rds")
cb_res = read_rds("output/03_predicted_inla_spatial_casecrossover.rds")

#------------------------------------------------------
# load ward boundaries
#------------------------------------------------------
ward_map = read_sf("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp")

#------------------------------------------------------
# load posterior mmt draws
#------------------------------------------------------
MMT = read_rds("output/04_mmt_draws_by_ward_new.rds")

#======================================================
# calculate ward-specific cumulative rr
#======================================================
ward_specific_cumul_RR = vector("list", 69)

cat("Starting cumulative RR calculation for 69 wards...\n")

for (i in 1:69) {
  
  cat("\n--------------------------------------------------\n")
  cat(sprintf("Ward %d / 69 started at %s\n", i, Sys.time()))
  cat("--------------------------------------------------\n")
  
  #------------------------------------------------------
  # ward-specific posterior inputs
  #------------------------------------------------------
  mmt_i = MMT[[i]] %>% 
    arrange(nsim) %>% 
    pull(MMT)
  
  # posterior coefficient draws
  # dimensions: simulations x basis coefficients
  beta_reg = cb_res[[i]]
  
  stopifnot(length(mmt_i) == nrow(beta_reg))
  
  #------------------------------------------------------
  # subset ward-specific data
  #------------------------------------------------------
  df_w = df_complete[df_complete$new_id == i, ]
  
  obs_temp   = df_w$tasmean
  obs_deaths = df_w$deaths
  
  #------------------------------------------------------
  # recreate original cross-basis
  #------------------------------------------------------
  cb_daily = crossbasis(
    obs_temp,
    lag = 21,
    argvar = list(
      fun   = "bs",
      knots = quantile(df_complete$tasmean, probs = c(0.1, 0.75, 0.9), na.rm = TRUE)
    ),
    arglag = list(
      fun   = "ns",
      knots = logknots(21, 3)
    )
  )
  
  #------------------------------------------------------
  # observed lagged temperature matrix
  #------------------------------------------------------
  at = df_w %>% 
    select(tasmean, contains("lag")) %>% 
    rename(lag0 = tasmean)
  
  stopifnot(ncol(cb_daily) == ncol(beta_reg))
  
  at_mat = as.matrix(at)
  
  #------------------------------------------------------
  # prediction setup
  #------------------------------------------------------
  N = nrow(at_mat)
  
  predvar = seq_len(N)
  predlag = dlnm:::seqlag(c(0, 21))
  
  # first 21 days do not have complete lag history
  ok = seq_len(N) > 21 & !is.na(obs_deaths)
  
  n_sim = nrow(beta_reg)
  
  cumlogRR = matrix(
    NA_real_,
    nrow = sum(ok),
    ncol = n_sim
  )
  
  #------------------------------------------------------
  # identify unique posterior mmt values
  #------------------------------------------------------
  unique_mmt = unique(mmt_i)
  
  cat(sprintf("Ward %d: %d valid days after lag exclusion\n", i, sum(ok)))
  cat(sprintf(
    "Ward %d: %d posterior draws, %d unique MMT values\n",
    i,
    n_sim,
    length(unique_mmt)
  ))
  
  #======================================================
  # calculate by unique mmt
  #======================================================
  for (m in seq_along(unique_mmt)) {
    
    cen_mmt = unique_mmt[m]
    
    #------------------------------------------------------
    # draws sharing this mmt
    #------------------------------------------------------
    draw_idx = which(mmt_i == cen_mmt)
    
    #------------------------------------------------------
    # construct prediction matrix once for this mmt
    #------------------------------------------------------
    Xpred_m = dlnm:::mkXpred(
      type    = "cb",
      basis   = cb_daily,
      at      = at_mat,
      predvar = predvar,
      predlag = predlag,
      cen     = cen_mmt
    )
    
    #------------------------------------------------------
    # sum effects across lags
    #------------------------------------------------------
    Xsum_m = matrix(
      0,
      nrow = N,
      ncol = ncol(Xpred_m)
    )
    
    for (l in seq_along(predlag)) {
      ind = predvar + N * (l - 1)
      Xsum_m = Xsum_m + Xpred_m[ind, , drop = FALSE]
    }
    
    #------------------------------------------------------
    # remove incomplete lag days
    #------------------------------------------------------
    Xsum_m = Xsum_m[ok, , drop = FALSE]
    
    #------------------------------------------------------
    # apply all coefficient draws sharing this mmt
    #------------------------------------------------------
    cumlogRR[, draw_idx] = Xsum_m %*% t(beta_reg[draw_idx, , drop = FALSE])
    
    if (m %% 10 == 0 || m == length(unique_mmt)) {
      cat(sprintf(
        "Ward %d: MMT %d / %d completed\n",
        i,
        m,
        length(unique_mmt)
      ))
    }
  }
  
  #------------------------------------------------------
  # store ward results
  #------------------------------------------------------
  stopifnot(!anyNA(cumlogRR))
  
  ward_specific_cumul_RR[[i]] = cumlogRR
  
  cat(sprintf("Ward %d completed at %s\n", i, Sys.time()))
  
  #------------------------------------------------------
  # save progress after each ward
  #------------------------------------------------------
  write_rds(
    ward_specific_cumul_RR,
    "output/05_ward_specific_cumul_RR_progress.rds"
  )
}

#======================================================
# save final result
#======================================================
# write_rds( ward_specific_cumul_RR,"output/05_ward_specific_cumul_RR.rds")

cat("\n==============================================\n")
cat("ALL WARDS COMPLETED SUCCESSFULLY\n")
cat(sprintf("Finished at %s\n", Sys.time()))
cat("==============================================\n")


# ward_specific_cumul_RR = read_rds("output/05_ward_specific_cumul_RR.rds")

##########################################################################
# ======================================================================
# Compute daily temperature-attributable excess mortality by ward
# ----------------------------------------------------------------------

excess_mortality_daily_by_ward = list()

for (i in 1:69){
  
  df_w = df_complete[df_complete$new_id == i, ]
  obs_deaths = df_w$deaths
  at = df_w %>% 
    select(tasmean, contains("lag")) %>% 
    rename(lag0 = tasmean)
  
  ok = seq_len(nrow(at)) > 21 & !is.na(obs_deaths)  
  
  # ------------------------------------------------------------------
  # AF is computed as:
  #   AF = 1 - exp(-logRR)
  # Backward attributable fraction (AF)
  # ------------------------------------------------------------------
  AF_daily = 1 - exp(-ward_specific_cumul_RR[[i]])
  
  # ------------------------------------------------------------------
  # Multiplying AF by observed deaths gives the attributable
  # number of deaths for each day and posterior simulation.
  #Backward attributable number (AN)
  # ------------------------------------------------------------------
  deaths_ok = obs_deaths[ok]
  AN_daily = AF_daily * deaths_ok
  
  excess_mortality_daily_by_ward[[i]] = AN_daily 
  
}

write_rds(excess_mortality_daily_by_ward, "output/05_excess_mortality_daily_by_ward.rds")

excess_mortality_daily_by_ward = readRDS("output/05_excess_mortality_daily_by_ward.rds")
##########################################################
#==================================================================
#total excess deaths for ward 
# posterior distributions (lists of length 69)
heat_annual_post = vector("list", 69)
cold_annual_post = vector("list", 69)

# summary table
ward_EM_summary = data.frame(
  Ward_id   = integer(69),
  Ward_code = character(69),
  
  heat_med  = numeric(69),
  heat_LL   = numeric(69),
  heat_UL   = numeric(69),
  
  cold_med  = numeric(69),
  cold_LL   = numeric(69),
  cold_UL   = numeric(69)
)



for (i in 1:69) {
  
  # --------------------------------------------------
  # Ward-specific data
  # --------------------------------------------------
  df_w = df_complete[df_complete$new_id == i, ]
  obs_deaths = df_w$deaths
  
  at = df_w %>%
    select(tasmean, contains("lag")) %>%
    rename(lag0 = tasmean)
  
  ok = seq_len(nrow(at)) > 21 & !is.na(obs_deaths)
  
  testing = data.frame(
    Date = as.Date(df_w$date[ok]),
    year = year(as.Date(df_w$date[ok])),
    daily_tasmean = df_w$tasmean[ok],
    Ward_code = unique(df_w$ward22cd)
  )
  
  n_years = length(unique(testing$year))
  
  # --------------------------------------------------
  # Ward-specific MMT (posterior)
  # --------------------------------------------------
  mmt_i = MMT[[i]] %>%
    arrange(nsim) %>%
    pull(MMT)
  
  nsim = length(mmt_i)
  nday = nrow(testing)
  
  # temperature classification per posterior draw
  temp_class_mat = matrix(NA_character_, nrow = nday, ncol = nsim)
  
  for (j in 1:nsim) {
    temp_class_mat[, j] = ifelse(
      testing$daily_tasmean > mmt_i[j], "heat",
      ifelse(testing$daily_tasmean < mmt_i[j], "cold", "atMMT")
    )
  }
  
  # --------------------------------------------------
  # Daily excess mortality (already computed)
  # --------------------------------------------------
  EM_for_ward_i = excess_mortality_daily_by_ward[[i]]
  
  heat_annual_i = numeric(nsim)
  cold_annual_i = numeric(nsim)
  
  for (j in 1:nsim) {
    
    AN_j   = EM_for_ward_i[, j]
    temp_j = temp_class_mat[, j]
    
    heat_total_j = sum(AN_j[temp_j == "heat"], na.rm = TRUE)
    cold_total_j = sum(AN_j[temp_j == "cold"], na.rm = TRUE)
    
    heat_annual_i[j] = heat_total_j / n_years
    cold_annual_i[j] = cold_total_j / n_years
  }
  
  # --------------------------------------------------
  # Store posterior distributions
  # --------------------------------------------------
  heat_annual_post[[i]] = heat_annual_i
  cold_annual_post[[i]] = cold_annual_i
  
  # --------------------------------------------------
  # Store summary statistics (median + 95% CrI)
  # --------------------------------------------------
  ward_EM_summary$Ward_id[i]   = i
  ward_EM_summary$Ward_code[i] = unique(testing$Ward_code)
  
  ward_EM_summary$heat_med[i] = median(heat_annual_i)
  ward_EM_summary$heat_LL[i]  = quantile(heat_annual_i, 0.025)
  ward_EM_summary$heat_UL[i]  = quantile(heat_annual_i, 0.975)
  
  ward_EM_summary$cold_med[i] = median(cold_annual_i)
  ward_EM_summary$cold_LL[i]  = quantile(cold_annual_i, 0.025)
  ward_EM_summary$cold_UL[i]  = quantile(cold_annual_i, 0.975)
}


write_rds(heat_annual_post,"output/05_heat_annual_post_EM_ward.rds")
write_rds(cold_annual_post,"output/05_cold_annual_post_EM_ward.rds")
write_rds(ward_EM_summary, "output/05_heat_and_cold_related_EM_ward.rds")


# heat_annual_post = read_rds("output/05_heat_annual_post_EM_ward.rds")
# cold_annual_post = read_rds("output/05_cold_annual_post_EM_ward.rds")
# ward_EM_summary = read_rds("output/05_heat_and_cold_related_EM_ward.rds")
##########################################################################
X_heat_annual_post = vector("list", 69)
X_cold_annual_post = vector("list", 69)
# summary table
ward_X_EM_summary <- data.frame(
  Ward_id   = integer(69),
  Ward_code = character(69),
  
  heat_med  = numeric(69),
  heat_LL   = numeric(69),
  heat_UL   = numeric(69),
  
  cold_med  = numeric(69),
  cold_LL   = numeric(69),
  cold_UL   = numeric(69)
)

nsim <- 1000  # set once

for (i in 1:69) {
  
  # --------------------------------------------------
  # Ward-specific data
  # --------------------------------------------------
  df_w <- df_complete[df_complete$new_id == i, ]
  obs_deaths <- df_w$deaths
  
  at <- df_w %>%
    dplyr::select(tasmean, dplyr::contains("lag")) %>%
    dplyr::rename(lag0 = tasmean)
  
  ok <- seq_len(nrow(at)) > 21 & !is.na(obs_deaths)
  
  testing <- data.frame(
    Date         = as.Date(df_w$date[ok]),
    year         = lubridate::year(as.Date(df_w$date[ok])),
    daily_tasmean= df_w$tasmean[ok],
    Ward_code    = unique(df_w$ward22cd)
  )
  
  n_years <- length(unique(testing$year))
  
  # --------------------------------------------------
  # Extreme temperature classification (ward-specific, NOT draw-specific)
  # --------------------------------------------------
  T <- testing$daily_tasmean
  
  q_hi <- as.numeric(stats::quantile(T, 0.975, na.rm = TRUE))
  q_lo <- as.numeric(stats::quantile(T, 0.025, na.rm = TRUE))
  
  is_heat <- T > q_hi
  is_cold <- T < q_lo
  
  # Optional: if you want to ensure extremes exist
  # stopifnot(any(is_heat), any(is_cold))
  
  # --------------------------------------------------
  # Daily excess mortality (already computed)
  # --------------------------------------------------
  EM_for_ward_i <- excess_mortality_daily_by_ward[[i]]
  
  # Safety check: rows must align with testing days
  if (nrow(EM_for_ward_i) != nrow(testing)) {
    stop(sprintf(
      "Row mismatch in ward %d: EM has %d rows, testing has %d rows",
      i, nrow(EM_for_ward_i), nrow(testing)
    ))
  }
  
  # Safety check: columns (draws)
  if (ncol(EM_for_ward_i) < nsim) {
    stop(sprintf(
      "Draw mismatch in ward %d: EM has %d draws/cols but nsim=%d",
      i, ncol(EM_for_ward_i), nsim
    ))
  }
  
  heat_annual_i <- numeric(nsim)
  cold_annual_i <- numeric(nsim)
  
  for (j in 1:nsim) {
    
    AN_j <- EM_for_ward_i[, j]
    
    heat_total_j <- sum(AN_j[is_heat], na.rm = TRUE)
    cold_total_j <- sum(AN_j[is_cold], na.rm = TRUE)
    
    heat_annual_i[j] <- heat_total_j / n_years
    cold_annual_i[j] <- cold_total_j / n_years
  }
  
  # --------------------------------------------------
  # Store posterior distributions
  # --------------------------------------------------
  X_heat_annual_post[[i]] <- heat_annual_i
  X_cold_annual_post[[i]] <- cold_annual_i
  
  # --------------------------------------------------
  # Store summary statistics (median + 95% CrI)
  # --------------------------------------------------
  ward_X_EM_summary$Ward_id[i]   <- i
  ward_X_EM_summary$Ward_code[i] <- unique(testing$Ward_code)
  
  ward_X_EM_summary$heat_med[i] <- stats::median(heat_annual_i, na.rm = TRUE)
  ward_X_EM_summary$heat_LL[i]  <- stats::quantile(heat_annual_i, 0.025, na.rm = TRUE)
  ward_X_EM_summary$heat_UL[i]  <- stats::quantile(heat_annual_i, 0.975, na.rm = TRUE)
  
  ward_X_EM_summary$cold_med[i] <- stats::median(cold_annual_i, na.rm = TRUE)
  ward_X_EM_summary$cold_LL[i]  <- stats::quantile(cold_annual_i, 0.025, na.rm = TRUE)
  ward_X_EM_summary$cold_UL[i]  <- stats::quantile(cold_annual_i, 0.975, na.rm = TRUE)
}

write_rds(X_heat_annual_post,"output/05_X_heat_annual_post_EM_ward.rds")
write_rds(X_cold_annual_post,"output/05_X_cold_annual_post_EM_ward.rds")
write_rds(ward_X_EM_summary, "output/05_X_heat_and_cold_related_EM_ward.rds")

# X_heat_annual_post = read_rds("output/05_X_heat_annual_post_EM_ward.rds")
# X_cold_annual_post = read_rds("output/05_X_cold_annual_post_EM_ward.rds")
# ward_X_EM_summary = read_rds("output/05_X_heat_and_cold_related_EM_ward.rds")

#====================================================================
#calculate per 100,000 rate
#--------------------------------------------------------------------

ward_EM_summary = read_rds("output/05_heat_and_cold_related_EM_ward.rds")



ward_pop = read_excel("data/external/population/sapewardstablefinal.xlsx", 
                      sheet = "Mid-2022 Ward 2022", skip = 3)


ward_pop = ward_pop %>% 
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
         cold_UL = cold_UL/count*100000
  )

print(ward_pop, n=69)


#plot the map



tmap_mode("plot")

cold_related_em = tm_shape(ward_map %>% 
                             left_join(ward_pop, by = c("Ward_Code"="Ward_code")))+
  tm_polygons(
    fill= "cold_med",
    fill.scale = tm_scale_continuous(values = "blues"),
    fill.legend = tm_legend("per 100,000", group_id = "top", frame = FALSE,bg.alpha=0)
  )+
  tm_layout(
    title.size = 1.2,
    legend.position = c(0.02, 0.92),
    frame = FALSE,
    inner.margins = c(0.07, 0, 0.15, 0) # Increased bottom margin (3rd number) to make room for caption
  ) +
  tm_title("Posterior Median Estimates of Annual Cold-Related Excess Mortality")+
  tm_credits(text = "Median cold-related excess mortality rate, calculated by comparing observed ward temperatures to \nward-specific MMTs (used as counterfactuals across posterior simulations)", 
             position = c("LEFT","TOP"))+
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


tmap_save(cold_related_em,filename = "figs/05_cold_related_em.png", height = 7,width =6, unit="in",dpi = 600 )



heat_related_em = tm_shape(ward_map %>% 
                             left_join(ward_pop, by = c("Ward_Code"="Ward_code")))+
  tm_polygons(
    fill= "heat_med",
    fill.scale = tm_scale_continuous(values = "reds"),
    fill.legend = tm_legend("per 100,000", group_id = "top", frame = FALSE,bg.alpha=0)
  )+
  tm_layout(
    title.size = 1.2,
    legend.position = c(0.02, 0.92),
    frame = FALSE,
    inner.margins = c(0.07, 0, 0.15, 0) # Increased bottom margin (3rd number) to make room for caption
  ) +
  tm_title("Posterior Median Estimates of Annual Heat-Related Excess Mortality")+
  tm_credits(text = "Median heat-related excess mortality rate, calculated by comparing observed ward temperatures to \nward-specific MMTs (used as counterfactuals across posterior simulations)", 
             position = c("LEFT","TOP"))+
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




tmap_save(heat_related_em,filename = "figs/05_heat_related_em.png", height = 7,width =6, unit="in",dpi = 600 )


#merge the plots
merged_cold_heat_em = tmap_arrange(cold_related_em,heat_related_em )
tmap_save(merged_cold_heat_em, filename ="figs/05_merged_cold_heat_em.png",height = 7,width =12, unit="in",dpi = 600)

##############################################################################################################
#calculate extreme heat and cold
#only use to >=97.5% temp for extreme heat
#only use <=2.5% temp for extreme cold

ward_X_EM_summary = read_rds("output/05_X_heat_and_cold_related_EM_ward.rds")

ward_pop = ward_pop %>% 
  filter(`LAD 2022 Name` == "Birmingham") %>% 
  pivot_longer(cols = c(-`LAD 2022 Name`,-`LAD 2022 Name`,-`Ward 2022 Code`,-`Ward 2022 Name`,-Total,
                        -`LAD 2022 Code`),
               names_to = "age",
               values_to = "count") %>% 
  group_by(`Ward 2022 Code`,`Ward 2022 Name`) %>% 
  summarise(count = sum(count)) %>% 
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
         cold_UL = cold_UL/count*100000)





tmap_mode("plot")

Xcold_related_em = tm_shape(ward_map %>% 
                              left_join(ward_pop, by = c("Ward_Code"="Ward_code")))+
  tm_polygons(
    fill= "cold_med",
    fill.scale = tm_scale_continuous(values = "blues"),
    fill.legend = tm_legend("per 100,000", group_id = "top", frame = FALSE,bg.alpha=0)
  )+
  tm_layout(
    legend.position = c(0.02, 0.92),
    frame = FALSE,
    inner.margins = c(0.07, 0, 0.15, 0), # Increased bottom margin (3rd number) to make room for caption
    component.autoscale = TRUE
  ) +
  tm_title("Posterior Median Estimates of Annual Excess Mortality due to Extreme Cold")+
  tm_credits(
    text = "Posterior median annual excess deaths on extreme-cold days, defined as days with ward-specific daily mean \ntemperature at or below the 2.5th percentile",
    position = c("LEFT","TOP")    # = make this bigger (try 0.9–1.1)
    
  )+
  tm_compass(type = "8star",
             size = 3,
             position = c("RIGHT", "bottom"),
             color.light = "white")+
  tm_credits(
    text = paste("Contains OS data \u00A9 Crown copyright and database right",
                 # Get current year
                 format(Sys.Date(), "%Y"),
                 ". Source:\nOffice for National Statistics licensed under the Open Government Licence v.3.0."),
    position = c("LEFT", "BOTTOM")
  )


tmap_save(Xcold_related_em ,filename = "figs/05_Xcold_related_em.png", height = 7,width =6, unit="in",dpi = 600 )




Xheat_related_em = tm_shape(ward_map %>% 
                              left_join(ward_pop, by = c("Ward_Code"="Ward_code")))+
  tm_polygons(
    fill= "heat_med",
    fill.scale = tm_scale_continuous(values = "reds"),
    fill.legend = tm_legend("per 100,000", group_id = "top", frame = FALSE,bg.alpha=0)
  )+
  tm_layout(
    legend.position = c(0.02, 0.92),
    frame = FALSE,
    inner.margins = c(0.07, 0, 0.15, 0) # Increased bottom margin (3rd number) to make room for caption
  ) +
  tm_title("Posterior Median Estimates of Annual Excess Mortality due to Extreme Heat")+
  tm_credits(text = "Posterior median annual excess deaths on extreme-cold days, defined as days with ward-specific daily mean \ntemperature at or above the 97.5th percentile", 
             position = c("LEFT","TOP"))+
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



tmap_save(Xheat_related_em,filename = "figs/05_Xheat_related_em.png", height = 7,width =6, unit="in",dpi = 600 )


#merge the plots
merged_Xcold_Xheat_em = tmap_arrange(Xcold_related_em,Xheat_related_em )
tmap_save(merged_Xcold_Xheat_em, filename ="figs/05_merged_Xcold_Xheat_em.png",height = 7,width =12, unit="in",dpi = 600)

##############################################################################################################
#to calculate the exceedance probability > the bimrigham mean excess mortality 


ward_pop = read_excel("data/external/population/sapewardstablefinal.xlsx", 
                      sheet = "Mid-2022 Ward 2022", skip = 3)



bham_excess_cold_list = vector("list", 1000)

for (i in 1:1000){
  # 2. Use a temporary name for the math (e.g., current_sum)
  current_sum = sum(sapply(cold_annual_post, `[`, i))
  
  # 3. Save it into the storage list
  bham_excess_cold_list [[i]] <- data.frame(
    draws = i,
    value = current_sum
  )
}

bham_excess_cold_df = bind_rows(bham_excess_cold_list)



bham_excess_heat_list = vector("list", 1000)

for (i in 1:1000){
  # 2. Use a temporary name for the math (e.g., current_sum)
  current_sum = sum(sapply(heat_annual_post, `[`, i))
  
  # 3. Save it into the storage list
  bham_excess_heat_list[[i]] <- data.frame(
    draws = i,
    value = current_sum
  )
}

bham_excess_heat_df = bind_rows(bham_excess_heat_list)


bham_pop = ward_pop %>% 
  filter(`LAD 2022 Name` == "Birmingham") %>% 
  pivot_longer(cols = c(-`LAD 2022 Name`,-`LAD 2022 Name`,-`Ward 2022 Code`,-`Ward 2022 Name`,-Total,
                        -`LAD 2022 Code`),
               names_to = "age",
               values_to = "count") %>% 
  group_by(`Ward 2022 Code`,`Ward 2022 Name`) %>% 
  summarise(count = sum(count)) %>% 
  rename(Ward_code = `Ward 2022 Code`,
         Ward_name = `Ward 2022 Name`) %>% 
  pull(count) %>% 
  sum()




bham_mean_cold_EM = mean(bham_excess_cold_df$value)/bham_pop*100000


bham_mean_heat_EM = mean(bham_excess_heat_df$value)/bham_pop*100000

quantile(bham_excess_cold_df$value, 0.975)

quantile(bham_excess_cold_df$value, 0.025)


quantile(bham_excess_heat_df$value, 0.975)

quantile(bham_excess_heat_df$value, 0.025)





bham_ward_pop = ward_pop %>% 
  filter(`LAD 2022 Name` == "Birmingham") %>% 
  pivot_longer(cols = c(-`LAD 2022 Name`,-`LAD 2022 Name`,-`Ward 2022 Code`,-`Ward 2022 Name`,-Total,
                        -`LAD 2022 Code`),
               names_to = "age",
               values_to = "count") %>% 
  group_by(`Ward 2022 Code`,`Ward 2022 Name`) %>% 
  summarise(count = sum(count)) %>% 
  rename(Ward_code = `Ward 2022 Code`,
         Ward_name = `Ward 2022 Name`) %>% 
  left_join(df_complete %>% 
              distinct(new_id,ward22cd), by = c("Ward_code" = "ward22cd"))

# apply(cold_annual_post[[1]],1,function(x) x/(bham_ward_pop$count[[1]])*100000)





exceedance_EM_cold = vector("list",69)

for (i in 1:69){
  
  exceedance_EM_cold[[i]] = mean(cold_annual_post[[i]]/(bham_ward_pop$count[[i]])*100000>bham_mean_cold_EM)
  
  
}


exceedance_EM_heat = vector("list",69)

for (i in 1:69){
  
  exceedance_EM_heat[[i]] = mean(heat_annual_post[[i]]/(bham_ward_pop$count[[i]])*100000>bham_mean_heat_EM)
  
  
}



exceed_heat = unlist(exceedance_EM_heat, use.names = FALSE)
exceed_cold = unlist(exceedance_EM_cold, use.names = FALSE)

exceed_df = tibble::tibble(
  Ward_id = 1:69,
  p_heat_gt_bham = exceed_heat,
  p_cold_gt_bham = exceed_cold
)



exceed_prob_EM_gt_bham_heatcold = exceed_df %>% 
  left_join(df_complete %>% 
              distinct(new_id,ward22cd), by = c("Ward_id" = "new_id"))


write_rds(exceed_prob_EM_gt_bham_heatcold, "output/05_exceed_prob_EM_gt_bham_heatcold.rds")

#------------------------------------------------------------------

em_exceedance_plot_df = ward_map %>% 
  left_join(exceed_prob_EM_gt_bham_heatcold, by = c("Ward_Code" = "ward22cd")) %>% 
  mutate(evidence_cold = case_when(is.na(p_cold_gt_bham)  ~ "Missing",
                                   p_cold_gt_bham  >= 0.95 ~ "Strong evidence (≥0.95)",
                                   p_cold_gt_bham  >= 0.90 ~ "Some evidence (0.90–0.95)",
                                   TRUE                 ~ "No evidence (<0.90)"),
         evidence_heat = case_when(is.na(p_heat_gt_bham)  ~ "Missing",
                                   p_heat_gt_bham  >= 0.95 ~ "Strong evidence (≥0.95)",
                                   p_heat_gt_bham  >= 0.90 ~ "Some evidence (0.90–0.95)",
                                   TRUE                 ~ "No evidence (<0.90)"),
         
         tooltip_cold = paste0(
           "<B>Ward Name:</B> ",Ward_Name, "\n",
           "<B>Ward Code:</B> ", Ward_Code, "\n",
           "<B>Pr(EM > mean EM):</B> ", p_cold_gt_bham, "\n",
           "<B>Inference:</B> ",evidence_cold),
         tooltip_heat = paste0(
           "<B>Ward Name:</B> ",Ward_Name, "\n",
           "<B>Ward Code:</B> ", Ward_Code, "\n",
           "<B>Pr(EM > mean EM):</B> ", p_heat_gt_bham, "\n",
           "<B>Inference:</B> ",evidence_heat) )


ggplot(em_exceedance_plot_df )+
  geom_sf_interactive(aes(fill = p_heat_gt_bham, tooltip = tooltip_heat, data_id = Ward_Code))+
  scale_fill_distiller(
    palette  = "Reds",
    direction = 1,
    na.value = "grey85",
    limits   = c(0.0,1.00)
  )+
  labs(fill = "Pr(EM>mean EM)") +
  ggtitle("Exceedance probability that heat-related excess mortality rate \nin the ward is higher than the overall Birmingham mean \nheat-realted excess mortality")+
  theme_void(base_size = 16) +
  theme()



ggplot(em_exceedance_plot_df )+
  geom_sf_interactive(aes(fill = p_cold_gt_bham,tooltip = tooltip_cold, data_id = Ward_Code))+
  scale_fill_distiller(
    palette  = "Blues",
    direction = 1,
    na.value = "grey85",
    limits   = c(0.0,1.00)
  )+
  labs(fill = "Pr(EM>mean EM)") +
  ggtitle("Exceedance probability that cold-related excess mortality rate \nin the ward is higher than the overall Birmingham mean \ncold-realted excess mortality")+
  theme_void(base_size = 16) +
  theme()


##############################################################################################################
#to calculate the exceedance probability > the bimrigham mean excess mortality at 2.5th and 97.5th percentile 

ward_pop = read_excel("data/external/population/sapewardstablefinal.xlsx", 
                      sheet = "Mid-2022 Ward 2022", skip = 3)


X_heat_annual_post = read_rds("output/05_X_heat_annual_post_EM_ward.rds")
X_cold_annual_post = read_rds("output/05_X_cold_annual_post_EM_ward.rds")




bham_excess_Xcold_list = vector("list", 1000)

for (i in 1:1000){
  # 2. Use a temporary name for the math (e.g., current_sum)
  current_sum = sum(sapply(X_cold_annual_post, `[`, i))
  
  # 3. Save it into the storage list
  bham_excess_Xcold_list [[i]] <- data.frame(
    draws = i,
    value = current_sum
  )
}

bham_excess_Xcold_df = bind_rows(bham_excess_Xcold_list)





bham_excess_Xheat_list = vector("list", 1000)

for (i in 1:1000){
  # 2. Use a temporary name for the math (e.g., current_sum)
  current_sum = sum(sapply(X_heat_annual_post, `[`, i))
  
  # 3. Save it into the storage list
  bham_excess_Xheat_list [[i]] <- data.frame(
    draws = i,
    value = current_sum
  )
}

bham_excess_Xheat_df = bind_rows(bham_excess_Xheat_list)




bham_pop = ward_pop %>% 
  filter(`LAD 2022 Name` == "Birmingham") %>% 
  pivot_longer(cols = c(-`LAD 2022 Name`,-`LAD 2022 Name`,-`Ward 2022 Code`,-`Ward 2022 Name`,-Total,
                        -`LAD 2022 Code`),
               names_to = "age",
               values_to = "count") %>% 
  group_by(`Ward 2022 Code`,`Ward 2022 Name`) %>% 
  summarise(count = sum(count)) %>% 
  rename(Ward_code = `Ward 2022 Code`,
         Ward_name = `Ward 2022 Name`) %>% 
  pull(count) %>% 
  sum()




bham_mean_Xcold_EM = mean(bham_excess_Xcold_df$value)/bham_pop*100000


bham_mean_Xheat_EM = mean(bham_excess_Xheat_df$value)/bham_pop*100000

quantile(bham_excess_Xcold_df$value, 0.975)

quantile(bham_excess_Xcold_df$value, 0.025)


quantile(bham_excess_Xheat_df$value, 0.975)

quantile(bham_excess_Xheat_df$value, 0.025)

#------------------------------------------------------
# birmingham median extreme mortality rates
#------------------------------------------------------
bham_median_Xcold_EM = median(
  bham_excess_Xcold_df$value / bham_pop * 100000,
  na.rm = TRUE
)

bham_median_Xheat_EM = median(
  bham_excess_Xheat_df$value / bham_pop * 100000,
  na.rm = TRUE
)


#======================================================
# ward population indexed by new_id
#======================================================
bham_ward_pop = ward_pop %>% 
  filter(`LAD 2022 Name` == "Birmingham") %>% 
  pivot_longer(
    cols = c(
      -`LAD 2022 Name`,
      -`Ward 2022 Code`,
      -`Ward 2022 Name`,
      -Total,
      -`LAD 2022 Code`
    ),
    names_to = "age",
    values_to = "count"
  ) %>% 
  group_by(`Ward 2022 Code`, `Ward 2022 Name`) %>% 
  summarise(count = sum(count), .groups = "drop") %>% 
  rename(
    Ward_code = `Ward 2022 Code`,
    Ward_name = `Ward 2022 Name`
  ) %>% 
  left_join(
    df_complete %>% 
      distinct(new_id, ward22cd),
    by = c("Ward_code" = "ward22cd")
  ) %>% 
  arrange(new_id)

stopifnot(
  nrow(bham_ward_pop) == 69,
  !anyNA(bham_ward_pop$new_id),
  !anyDuplicated(bham_ward_pop$new_id),
  identical(bham_ward_pop$new_id, 1:69)
)

#======================================================
# extreme exceedance probabilities
#======================================================
exceedance_EM_Xcold = vector("list", 69)
exceedance_EM_Xheat = vector("list", 69)

for (i in 1:69) {
  
  ward_count = bham_ward_pop %>% 
    filter(new_id == i) %>% 
    pull(count)
  
  stopifnot(length(ward_count) == 1)
  
  cold_rate_draws = X_cold_annual_post[[i]] / ward_count * 100000
  heat_rate_draws = X_heat_annual_post[[i]] / ward_count * 100000
  
  exceedance_EM_Xcold[[i]] = mean(
    cold_rate_draws > bham_mean_Xcold_EM,
    na.rm = TRUE
  )
  
  exceedance_EM_Xheat[[i]] = mean(
    heat_rate_draws > bham_mean_Xheat_EM,
    na.rm = TRUE
  )
}

exceed_Xheat = unlist(exceedance_EM_Xheat, use.names = FALSE)
exceed_Xcold = unlist(exceedance_EM_Xcold, use.names = FALSE)

exceed_Xdf = tibble(
  Ward_id = 1:69,
  p_Xheat_gt_bham = exceed_Xheat,
  p_Xcold_gt_bham = exceed_Xcold
)

exceed_prob_EM_gt_bham_X_heatcold = exceed_Xdf %>% 
  left_join(
    df_complete %>% 
      distinct(new_id, ward22cd),
    by = c("Ward_id" = "new_id")
  )



write_rds(exceed_prob_EM_gt_bham_X_heatcold, "output/05_exceed_prob_EM_gt_bham_X_heatcold.rds")



em_exceedance_Xplot_df = ward_map %>% 
  left_join(exceed_prob_EM_gt_bham_X_heatcold, by = c("Ward_Code" = "ward22cd")) %>% 
  mutate(evidence_cold = case_when(is.na(p_Xcold_gt_bham)  ~ "Missing",
                                   p_Xcold_gt_bham  >= 0.95 ~ "Strong evidence (≥0.95)",
                                   p_Xcold_gt_bham  >= 0.90 ~ "Some evidence (0.90–0.95)",
                                   TRUE                 ~ "No evidence (<0.90)"),
         evidence_heat = case_when(is.na(p_Xheat_gt_bham)  ~ "Missing",
                                   p_Xheat_gt_bham  >= 0.95 ~ "Strong evidence (≥0.95)",
                                   p_Xheat_gt_bham  >= 0.90 ~ "Some evidence (0.90–0.95)",
                                   TRUE                 ~ "No evidence (<0.90)"),
         
         tooltip_cold = paste0(
           "<B>Ward Name:</B> ",Ward_Name, "\n",
           "<B>Ward Code:</B> ", Ward_Code, "\n",
           "<B>Pr(EM > mean EM):</B> ", p_Xcold_gt_bham, "\n",
           "<B>Inference:</B> ",evidence_cold),
         tooltip_heat = paste0(
           "<B>Ward Name:</B> ",Ward_Name, "\n",
           "<B>Ward Code:</B> ", Ward_Code, "\n",
           "<B>Pr(EM > mean EM):</B> ", p_Xheat_gt_bham, "\n",
           "<B>Inference:</B> ",evidence_heat) )




ggplot(em_exceedance_Xplot_df )+
  geom_sf_interactive(aes(fill = p_Xheat_gt_bham, tooltip = tooltip_heat, data_id = Ward_Code))+
  scale_fill_distiller(
    palette  = "Reds",
    direction = 1,
    na.value = "grey85",
    limits   = c(0.0,1.00)
  )+
  labs(fill = "Pr(EM>mean EM)") +
  ggtitle("Exceedance probability that excess mortality rate above 97.5th \npercentile in the ward is higher than the overall Birmingham mean \nheat-realted excess mortality")+
  theme_void(base_size = 16) +
  theme()




ggplot(em_exceedance_Xplot_df )+
  geom_sf_interactive(aes(fill = p_Xcold_gt_bham, tooltip = tooltip_heat, data_id = Ward_Code))+
  scale_fill_distiller(
    palette  = "Blues",
    direction = 1,
    na.value = "grey85",
    limits   = c(0.0,1.00)
  )+
  labs(fill = "Pr(EM>mean EM)") +
  ggtitle("Exceedance probability that excess mortality rate below 2.5th \npercentile in the ward is higher than the overall Birmingham mean \nheat-realted excess mortality")+
  theme_void(base_size = 16) +
  theme()




#======================================================
# birmingham annual extreme excess mortality
#======================================================
bham_extreme_EM_summary = tibble(
  exposure = c("Extreme cold", "Extreme heat"),
  median = c(
    median(bham_excess_Xcold_df$value, na.rm = TRUE),
    median(bham_excess_Xheat_df$value, na.rm = TRUE)
  ),
  LL = c(
    quantile(bham_excess_Xcold_df$value, 0.025, na.rm = TRUE),
    quantile(bham_excess_Xheat_df$value, 0.025, na.rm = TRUE)
  ),
  UL = c(
    quantile(bham_excess_Xcold_df$value, 0.975, na.rm = TRUE),
    quantile(bham_excess_Xheat_df$value, 0.975, na.rm = TRUE)
  )
)

bham_extreme_EM_summary


#======================================================
# updated ward extreme mortality results
#======================================================
ward_population = read_excel(
  "data/external/population/sapewardstablefinal.xlsx",
  sheet = "Mid-2022 Ward 2022",
  skip = 3
) %>% 
  filter(`LAD 2022 Name` == "Birmingham") %>% 
  pivot_longer(
    cols = c(
      -`LAD 2022 Name`,
      -`Ward 2022 Code`,
      -`Ward 2022 Name`,
      -Total,
      -`LAD 2022 Code`
    ),
    names_to = "age",
    values_to = "count"
  ) %>% 
  group_by(`Ward 2022 Code`, `Ward 2022 Name`) %>% 
  summarise(count = sum(count), .groups = "drop") %>% 
  rename(
    Ward_code = `Ward 2022 Code`,
    Ward_name = `Ward 2022 Name`
  )

extreme_ward_results = ward_X_EM_summary %>% 
  left_join(ward_population, by = "Ward_code") %>% 
  mutate(
    heat_med = heat_med / count * 100000,
    heat_LL  = heat_LL / count * 100000,
    heat_UL  = heat_UL / count * 100000,
    cold_med = cold_med / count * 100000,
    cold_LL  = cold_LL / count * 100000,
    cold_UL  = cold_UL / count * 100000
  ) %>% 
  left_join(
    exceed_prob_EM_gt_bham_X_heatcold,
    by = c("Ward_id", "Ward_code" = "ward22cd")
  )



#------------------------------------------------------
# highest extreme cold mortality
#------------------------------------------------------
extreme_ward_results %>% 
  arrange(desc(cold_med)) %>% 
  select(
    Ward_name,
    cold_med,
    cold_LL,
    cold_UL,
    p_Xcold_gt_bham
  ) %>% 
  slice_head(n = 10)


#------------------------------------------------------
# highest extreme heat mortality
#------------------------------------------------------
extreme_ward_results %>% 
  arrange(desc(heat_med)) %>% 
  select(
    Ward_name,
    heat_med,
    heat_LL,
    heat_UL,
    p_Xheat_gt_bham
  ) %>% 
  slice_head(n = 10)


















