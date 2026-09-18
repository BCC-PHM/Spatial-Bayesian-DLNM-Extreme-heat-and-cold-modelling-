# ============================================================
# Parallel cumulative RR for UKCP18 
# ============================================================

library(tidyverse)
library(dlnm)
library(tsModel)
library(doParallel)
library(doSNOW)
library(foreach)
library(iterators)

#load the data

df_complete = read_rds("output/df_complete_for_inla.rds")

cb_res = read_rds("output/predicted_inla_spatial_casecrossover.rds")



# Load Ward shp file

ward_map = read_sf("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp")


#load mmt resutls

MMT = read_rds("output/mmt_draws_by_ward.rds")


#load the fitted model

load("output/R_inla_eb_gaussian_results.RData")

rm(cb_df,cb_list,df,first_stage_LAD_age_sepcific,g,lag_df,lag_mat,mcnty_nb)


#extracted the median of fitted value and average across ward to obtain baseline deaths

baseline_mort = inla_model$summary.fitted.values[, "0.5quant", drop = FALSE] %>%
  as.data.frame() %>%
  rename(expected_deaths = `0.5quant`) %>%
  cbind(df_complete) %>%
  drop_na() %>%
  group_by(ward22cd, new_id) %>%
  summarise(mean_baseline = mean(expected_deaths))



#load the daily future time series for each ward

UKCP18CM_ward = read_rds("data/processed/chessscape85_daily_ward.rds")


#----------------------------------------------------

# Create a combined grouping factor based on both 'ward22cd' and 'Member'

group_factor = interaction(UKCP18CM_ward$Ward_Code , UKCP18CM_ward$Member)



# create lag matrix (includes lag 0)

UKCP18CM_ward_lags = tsModel::Lag(
   UKCP18CM_ward$tasmean,
   group = group_factor,
    k     = 0:21)


# Name columns for clarity

UKCP18CM_ward_lags = as.data.frame(UKCP18CM_ward_lags)

colnames(UKCP18CM_ward_lags) = paste0("lag", 0:21)


UKCP18CM_ward = bind_cols(UKCP18CM_ward, UKCP18CM_ward_lags)

UKCP18CM_ward = UKCP18CM_ward %>% 
  arrange(ID,Member)


#remove lag0
# UKCP18CM_ward = UKCP18CM_ward %>%
#   select(-lag0)


# ------------------------------------------------------------
# 0) Create task table (ward x member)
# ------------------------------------------------------------
ward_codes = unique(df_complete$ward22cd)  # should be 69

task_df = map_dfr(seq_along(ward_codes), function(i) {
  w = ward_codes[i]
  members = unique(UKCP18CM_ward$Member[UKCP18CM_ward$Ward_Code == w])
  tibble(i = i, Ward_Code = w, Member = members)
})

cat("Total tasks (ward-member combos):", nrow(task_df), "\n")

# ------------------------------------------------------------
# 1) Start cluster + register doSNOW for progress
# ------------------------------------------------------------
n_cores = max(1, parallel::detectCores() - 4)
cl = parallel::makeCluster(n_cores)

doSNOW::registerDoSNOW(cl)  # = enables progress callback
cat("Registered", n_cores, "cores\n")

# Progress bar (ticks once per completed task)
pb = txtProgressBar(min = 0, max = nrow(task_df), style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

# ------------------------------------------------------------
# 2) Run parallel foreach
# ------------------------------------------------------------
start = Sys.time()

ward_specific_cumul_RR_UKCP18 = foreach(
  t = iterators::iter(task_df, by = "row"),
  .packages = c("dlnm", "tsModel", "dplyr"),
  .export   = c("cb_res", "MMT", "baseline_mort", "df_complete", "UKCP18CM_ward", "ward_codes"),
  .options.snow = opts,
  .errorhandling = "stop"
) %dopar% {
  
  # ---- unpack task FIRST ----
  i = t$i
  w = t$Ward_Code
  m = t$Member
  
  # ---- logging setup ----
  dir.create("output/logs", showWarnings = FALSE, recursive = TRUE)
  log_file = file.path("output/logs", paste0("ward_", i, "_member_", m, ".log"))
  cat(sprintf("[%s] START (ward=%s member=%s)\n", Sys.time(), w, m),
      file = log_file, append = TRUE)
  
  # ---- Ward-specific posterior inputs ----
  beta_reg = cb_res[[i]]  # [nsim x K]
  
  mmt_i = MMT[[i]] %>%
    dplyr::arrange(nsim) %>%
    dplyr::pull(MMT)
  
  stopifnot(length(mmt_i) == nrow(beta_reg))
  
  # ---- Baseline for ward (not used in RR calc; kept for later AN calc) ----
  baseline_i = baseline_mort %>%
    dplyr::filter(ward22cd == w) %>%
    dplyr::pull(mean_baseline)
  
  stopifnot(length(baseline_i) == 1)
  
  # ---- Future temps for ward-member ----
  df_future_w = UKCP18CM_ward %>%
    dplyr::filter(Ward_Code == w, Member == m) %>%
    dplyr::arrange(Date)
  
  future_temp = df_future_w$tasmean
  
  # ---- Lag matrix ----
  at_mat = tsModel:::Lag(future_temp, 0:21)
  colnames(at_mat) = paste0("lag", 0:21)
  
  N = nrow(at_mat)
  ok = seq_len(N) > 21 & complete.cases(at_mat)
  
  # ---- Crossbasis ----
  cb_daily = crossbasis(
    future_temp,
    lag    = 21,
    argvar = list(
      fun   = "bs",
      knots = quantile(df_complete$tasmean, probs = c(0.1, 0.75, 0.9), na.rm = TRUE)
    ),
    arglag = list(fun = "ns", knots = logknots(21, 3))
  )
  
  stopifnot(ncol(cb_daily) == ncol(beta_reg))
  
  predvar = seq(N)
  predlag = dlnm:::seqlag(c(0, 21))
  
  # ---- Posterior simulations ----
  nsim = nrow(beta_reg)
  cumlogRR = matrix(NA_real_, nrow = sum(ok), ncol = nsim)
  
  for (j in seq_len(nsim)) {
    
    Xpred_j = dlnm:::mkXpred(
      type     = "cb",
      basis    = cb_daily,
      at       = at_mat,
      predvar  = predvar,
      predlag  = predlag,
      cen      = mmt_i[j]
    )
    
    Xsum_j = 0
    for (l in seq_along(predlag)) {
      ind = predvar + N * (l - 1)
      Xsum_j = Xsum_j + Xpred_j[ind, , drop = FALSE]
    }
    
    Xsum_j = Xsum_j[ok, , drop = FALSE]
    cumlogRR[, j] = Xsum_j %*% beta_reg[j, ]
  }
  
  # ---- log DONE BEFORE returning ----
  cat(sprintf("[%s] DONE (n_days=%d, nsim=%d)\n", Sys.time(), nrow(cumlogRR), ncol(cumlogRR)),
      file = log_file, append = TRUE)
  
  # ---- return named element ----
  setNames(list(cumlogRR), paste0(i, "_", m))
}

# Close progress bar
close(pb)

# Combine into one list
ward_specific_cumul_RR_UKCP18 = unlist(ward_specific_cumul_RR_UKCP18, recursive = FALSE)

end = Sys.time()
cat("\nFinished. Runtime:\n")
print(end - start)

# ------------------------------------------------------------
# 3) Stop cluster
# ------------------------------------------------------------
parallel::stopCluster(cl)
cat("Cluster stopped.\n")



write_rds(ward_specific_cumul_RR_UKCP18,"output/ward_specific_cumul_RR_UKCP18.rds")



