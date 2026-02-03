setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Spatial Bayesian DNLM Extreme heat and cold")

library(readr)
library(dlnm)
library(tidyverse)
library(sf)
library(tmap)
library(readxl)        
library(doParallel)
library(foreach)
###################################################################################

#load the data
df_complete = read_rds("output/df_complete_for_inla.rds")
cb_res = read_rds("output/predicted_inla_spatial_casecrossover.rds")


# Load Ward shp file
ward_map = read_sf("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp")


#load mmt resutls 
MMT = read_rds("output/mmt_draws_by_ward.rds")
####################################################################################
ward_specific_cumul_RR = list()

cat("Starting cumulative RR calculation for 69 wards...\n")

for (i in 1:69){
  
  cat("\n--------------------------------------------------\n")
  cat(sprintf("Ward %d / 69 started at %s\n", i, Sys.time()))
  cat("--------------------------------------------------\n")

# ==============================================================
# Ward-specific posterior inputs
# ==============================================================

cb_res[[i]]

        
#specify which ward's MMT to apply to the centering temp in Xpred
#it contains length of 1000
mmt_i = MMT[[i]] %>%
  arrange(nsim) %>%
  pull(MMT)


# Extract the posterior samples of the coefficients 
# This generates a matrix: [Number of posterior draws (1000)  x (30) basis coeficients]
beta_reg = cb_res[[i]]

stopifnot(length(mmt_i) == nrow(beta_reg))

# ==============================================================
# Subset ward-specific data
# ==============================================================

df_w = df_complete[df_complete$new_id == i, ]


# Extract the observed daily temperature and deaths for the analysis period
# Note: Ensure these vectors align perfectly with your model's lagged timeframe
obs_temp   = df_w$tasmean
obs_deaths = df_w$deaths


# Re-create the cross-basis for the ACTUAL daily time series
# (Use the exact same knots/lag settings as your main model)
# This generates a matrix: [Number of days (2853)  x (30) basis coeficients]
cb_daily = crossbasis(
  obs_temp,
  lag = 21,
  argvar = list(fun = "bs", knots = quantile(df_complete$tasmean, probs = c(0.1, 0.75, 0.9), na.rm=T)), 
  arglag = list(fun = "ns", knots = logknots(21, 3)) # Match your original model!
)


# use the lag exposure matrix the have pre-lagged before running inla 
#subset those columns
#if not use the following 
# at = tsModel:::Lag(obs_temp, seq(0, 21))  
#inhave already specified the lags before running inla
at = df_w %>% 
  select(tasmean, contains("lag")) %>% 
  rename(lag0 = tasmean)


#quick checks to make sure beta_reg has [nsim x K] and cb_daily should have K columns
stopifnot(ncol(cb_daily) == ncol(beta_reg))



# ==============================================================
# DLNM prediction setup
# Warning:: long running time!
# ==============================================================
#covert at to numeric matrix [number of days(2853) x number of lags(22)]
at_mat = as.matrix(at)

# #select a centering temp 
# cen = unique(MMT$MMT[MMT$Ward_id==i])

# PREPARE THE ARGUMENTS FOR TH BASIS TRANSFORMATION
predvar = seq(nrow(at_mat))
predlag = dlnm:::seqlag(c(0,21))
# ======================================================================
# lag-specific prediction blocks stacked vertically
#For today’s mortality, the model looks at today’s temperature AND the temperatures from the previous 21 days, and applies the cross-basis to all of them.
#for each posterior draw (1000), reconstruct the daily cumulative temp-mortality association across lags, 
#centred at the specific MMT that correspond to each draw, and project them to the draw's regression coefficient 
# ======================================================================
#all rows of the observed temp
N = nrow(at_mat)

#drop first 21 days because they do not have full lag history and any missing data
ok = seq_len(N) > 21 & !is.na(obs_deaths)    

#Initialise storage for posterior predictions
# Rows    = number of valid days (after lag exclusion)
# Columns = number of posterior simulations (1000)
#
# Column j will correspond to posterior draw j, conditional on:
#   - regression coefficients beta_reg[j, ]
#   - minimum mortality temperature mmt_i[j]
cumlogRR = matrix(NA, nrow = sum(ok), ncol = 1000)
# ------------------------------------------------------------------

cat(sprintf("Ward %d: %d valid days after lag exclusion\n", i, sum(ok)))
cat("Running 1000 posterior simulations...\n")

for (j in 1:1000) {
# ------------------------------------------------------------------
# mkXpred() rebuilds the DLNM design matrix using the observed
# temperature series, but re-centres the exposure–response function
# so that log-relative risk = 0 at mmt_i[j].
#
# The output Xpred_j is stacked by lag:
#   rows 1..N       -> lag 0
#   rows N+1..2N    -> lag 1
#   ...
#   rows 21N+1..22N -> lag 21
#
# Columns correspond to the DLNM basis coefficients.
# ------------------------------------------------------------------
  # --- simulation-specific centring ---
  Xpred_j = dlnm:::mkXpred(
    type     = "cb",
    basis   = cb_daily,
    at      = at_mat,
    predvar = predvar,
    predlag = predlag,
    cen     = mmt_i[j]
  )
  
  # ------------------------------------------------------------------
  # For each lag, extract the rows corresponding to that lag and
  # add them together. This yields the cumulative lag effect for
  # each day, already centred at mmt_i[j].
  # ------------------------------------------------------------------
  Xsum_j = 0
  for (l in seq_along(predlag)) {
    ind = predvar + N * (l - 1)
    Xsum_j = Xsum_j + Xpred_j[ind, , drop = FALSE]
  }
  
  # --- drop incomplete lag days ---
  Xsum_j = Xsum_j[ok, , drop = FALSE]
  
  # ------------------------------------------------------------------
  # --- apply posterior coefficients (row = simulation) ---
  # This generates a matrix: [Number of valid days  x (30) basis coeficients]
  # beta_reg[j, ] is what determines that only j-th posterior draw is used and only
  # j-th column of cumulative log-relative risk is produced.
  # ------------------------------------------------------------------
  cumlogRR[, j] = Xsum_j %*% beta_reg[j, ]
}
  
ward_specific_cumul_RR[[i]] = cumlogRR

cat(sprintf("Ward %d completed at %s\n", i, Sys.time()))

}

cat("\n==============================================\n")
cat("ALL WARDS COMPLETED SUCCESSFULLY \n")
cat(sprintf("Finished at %s\n", Sys.time()))
cat("==============================================\n")


write_rds(ward_specific_cumul_RR, "output/ward_specific_cumul_RR.rds")


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

write_rds(excess_mortality_daily_by_ward, "output/excess_mortality_daily_by_ward.rds")

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


write_rds(heat_annual_post,"output/heat_annual_post_EM_ward.rds")
write_rds(cold_annual_post,"output/cold_annual_post_EM_ward.rds")
write_rds(ward_EM_summary, "output/heat_and_cold_related_EM_ward.rds")

##########################################################################
#====================================================================
#calculate per 100,000 rate
#--------------------------------------------------------------------

ward_EM_summary = read_rds("output/heat_and_cold_related_EM_ward.rds")



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


tmap_save(cold_related_em,filename = "figs/cold_related_em.png", height = 7,width =6, unit="in",dpi = 600 )



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




tmap_save(heat_related_em,filename = "figs/heat_related_em.png", height = 7,width =6, unit="in",dpi = 600 )


#merge the plots
merged_cold_heat_em = tmap_arrange(cold_related_em,heat_related_em )
tmap_save(merged_cold_heat_em, filename ="figs/merged_cold_heat_em.png",height = 7,width =12, unit="in",dpi = 600)

##############################################################################################################
#calculate extreme heat and cold
#only use to >=95% temp for extreme heat
#only use <=5 % temp for extreme cold



















