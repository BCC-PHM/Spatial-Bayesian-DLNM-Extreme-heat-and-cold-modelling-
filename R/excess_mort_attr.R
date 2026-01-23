setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Spatial Bayesian DNLM Extreme heat and cold")

library(readr)
library(dlnm)
library(tidyverse)
library(sf)
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

for (i in 1:69){
  
  

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

}

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


##########################################################
#total excess deaths for ward i
testing = data.frame(
  Date = df_w$date[ok],
  Ward_code = unique(df_w$ward22cd[df_w$new_id ==1]),
  daily_tasmean = df_w$tasmean[ok]
  
) %>% 
  mutate(temp_class= case_when(
    daily_tasmean > cen ~ "heat",
    daily_tasmean < cen ~ "cold",
    TRUE                    ~ "atMMT"
  ),
  day_id = row_number(),
  year = year(Date))

####################################################

is_cold = testing$temp_class == "cold"

cold_EM = testing %>% 
  filter(temp_class == "cold") 


AN_cold_daily  =  AN_daily
AN_cold_daily[!is_cold, ] = 0


cold_yearly_sim =testing %>%
  group_by(year) %>%
  summarise(
    cold_sim = list(colSums(AN_cold_daily[day_id, , drop = FALSE])),
    .groups = "drop"
  )


cold_EM = cold_yearly_sim %>%
  mutate(
    EM_median = sapply(cold_sim, median),
    EM_LL     = sapply(cold_sim, quantile, 0.025),
    EM_UL     = sapply(cold_sim, quantile, 0.975)
  ) %>%
  select(year, EM_median, EM_LL, EM_UL)




