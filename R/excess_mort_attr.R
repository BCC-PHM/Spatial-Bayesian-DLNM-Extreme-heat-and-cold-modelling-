

#load the data
df_complete = read_rds("output/df_complete_for_inla.rds")
cb_res = read_rds("output/predicted_inla_spatial_casecrossover.rds")



# Load Ward shp file
ward_map = read_sf("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp")


#load mmt resutls 
MMT = read_rds("output/RR_MMT_plot_data.rds")


i=1
#plot all thje posterior draws for ward id 1 

cb_res[[i]]


# Extract the posterior samples of the coefficients (e.g., from INLA or MCMC).
# Dimensions: 1000 rows x 30 columns.
beta_reg = cb_res[[i]]

# Extract the cross-basis prediction matrix for this location.
# Dimensions: 119 rows x 30 columns 

attr_mortality_res <- data.frame(
  sim = 1:nrow(beta_reg),
  total_an = NA,     # Total Excess Deaths
  heat_an  = NA,     # All Heat
  cold_an  = NA,     # All Cold
  heat_ext_an = NA,  # Extreme Heat Only
  cold_ext_an = NA   # Extreme Cold Only
)


# 1. Subset data for ward
# -----------------------------
df_w = df_complete[df_complete$new_id == i, ]


# Extract the observed daily temperature and deaths for the analysis period
# Note: Ensure these vectors align perfectly with your model's lagged timeframe
obs_temp   <- df_w$tasmean
obs_deaths <- df_w$deaths


# Re-create the cross-basis for the ACTUAL daily time series
# (Use the exact same knots/lag settings as your main model)
# This generates a matrix: [N_days x 30 basis_vars]
cb_daily <- crossbasis(
  obs_temp,
  lag = 21,
  argvar = list(fun = "bs", knots = quantile(df_complete$tasmean, probs = c(0.1, 0.75, 0.9), na.rm=T)), 
  arglag = list(fun = "ns", knots = logknots(21, 3)) # Match your original model!
)


# use the lag exposure matrix the have pre-lagged before runing inla 
#subset those columns
#if not use the following 
# at = tsModel:::Lag(obs_temp, seq(0, 21))  

at = df_w %>% 
  select(tasmean, contains("lag")) %>% 
  rename(lag0 = tasmean)


#quick checks to make sure beta_reg has [nsim x K] and cb_daily should have K columns
stopifnot(ncol(cb_daily) == ncol(beta_reg))



##############################################################
#covert at to numeric matrix 
at_mat = as.matrix(at)

#select a centering temp 
cen = unique(MMT$MMT[MMT$Ward_id==i])

# PREPARE THE ARGUMENTS FOR TH BASIS TRANSFORMATION
predvar = seq(nrow(at_mat))
predlag = dlnm:::seqlag(c(0,21))
##############################################################
# lag-specific prediction blocks stacked vertically
#For today’s mortality, the model looks at today’s temperature AND the temperatures from the previous 21 days, and applies the cross-basis to all of them.

Xpred = dlnm:::mkXpred(
  type = "cb",
  basis     = cb_daily,
  at        = at_mat,
  predvar   = predvar,
  predlag   = predlag,
  cen       = cen)

#sum over lags to get one row per day 

K = ncol(beta_reg)          # number of basis coefficients
N = nrow(at_mat)            # number of days


#initialise the matrix as empty 
Xpredall = 0

# CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES 
for (l in seq_along(predlag)) {
  ind = predvar + N * (l - 1)
  Xpredall = Xpredall + Xpred[ind, , drop = FALSE]
}


#remove first 21 rows since they do not have full 21 lags 

ok = complete.cases(Xpredall) & !is.na(obs_deaths)

Xpredall_ok = Xpredall[ok, , drop = FALSE]
deaths_ok   = obs_deaths[ok]


#apply posterior draws
# cumulative log-relative risk per day per draw
cumlogRR = Xpredall_ok %*% t(beta_reg)

##########################################################

# Backward attributable fraction (AF)

AF_daily = 1 - exp(-cumlogRR)


##########################################################
#Backward attributable number (AN)
AN_daily = AF_daily * deaths_ok

##########################################################
#total excess deaths for ward i
attr_mortality_res$total_an = colSums(AN_daily, na.rm = TRUE)

attr_mortality_res %>% 
  pull(total_an) %>% 
  median()




