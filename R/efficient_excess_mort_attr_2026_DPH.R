#==============================================================================
# GENERATE HEAT-ONLY CUMULATIVE LOG-RR FOR 2026
# Uses the new draw-specific MMT distribution
#==============================================================================
library(tidyverse)
library(dlnm)
library(doParallel)
library(doSNOW)
library(foreach)
library(iterators)
library(RhpcBLASctl)
library(sf)

# Set one working directory only
setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Spatial Bayesian DNLM Extreme heat and cold")
#------------------------------------------------------
#load model inputs
#------------------------------------------------------
df_complete = read_rds("output/df_complete_for_inla.rds")
cb_res = read_rds("output/predicted_inla_spatial_casecrossover.rds")
MMT = read_rds("output/mmt_draws_by_ward_new_0_3.rds")
uk2026_heatwave_data = read_rds("data/processed/tmean_ward_daily_2026_DPH.rds") %>% 
  mutate(date = as.Date(date)) %>% 
  arrange(ID, date)
ward_map = read_sf("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp")
#------------------------------------------------------
#basic checks
#------------------------------------------------------
stopifnot(length(cb_res) == 69, length(MMT) == 69)

for (i in seq_along(cb_res)) {
  
  MMT[[i]] = MMT[[i]] %>% 
    mutate(nsim = as.integer(nsim)) %>% 
    arrange(nsim)
  
  stopifnot(
    all(c("nsim", "MMT") %in% names(MMT[[i]])),
    nrow(MMT[[i]]) == nrow(cb_res[[i]]),
    !anyDuplicated(MMT[[i]]$nsim),
    identical(MMT[[i]]$nsim, seq_len(nrow(cb_res[[i]])))
  )
}
#------------------------------------------------------
#ward lookup
#ensure list position i corresponds to new_id i
#------------------------------------------------------
ward_lookup = df_complete %>% 
  st_drop_geometry() %>% 
  distinct(new_id, Ward_Code = ward22cd) %>% 
  arrange(new_id)
stopifnot(nrow(ward_lookup) == 69, identical(ward_lookup$new_id, 1:69))
#------------------------------------------------------
#check temperature IDs align with model IDs
#------------------------------------------------------
temperature_lookup = uk2026_heatwave_data %>% 
  distinct(ID, Ward_Code) %>% 
  arrange(ID)
alignment_check = ward_lookup %>% 
  left_join(temperature_lookup, by = c("new_id" = "ID", "Ward_Code" = "Ward_Code"))
stopifnot(nrow(alignment_check) == 69, !anyNA(alignment_check$new_id))
#------------------------------------------------------
#MMT diagnostic
#------------------------------------------------------
mmt_diagnostic = map_dfr(seq_along(MMT), function(i) {
  tibble(
    Ward_id             = i,
    MMT_LL              = quantile(MMT[[i]]$MMT, 0.025, na.rm = TRUE),
    MMT_median          = median(MMT[[i]]$MMT, na.rm = TRUE),
    MMT_UL              = quantile(MMT[[i]]$MMT, 0.975, na.rm = TRUE),
    proportion_below_12 = mean(MMT[[i]]$MMT < 12, na.rm = TRUE)
  )
})
print(mmt_diagnostic %>% arrange(MMT_median))
print(summary(mmt_diagnostic$MMT_median))
#------------------------------------------------------
#use one BLAS thread per worker
#------------------------------------------------------
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
Sys.setenv(
  OMP_NUM_THREADS        = "1",
  OPENBLAS_NUM_THREADS   = "1",
  MKL_NUM_THREADS        = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS    = "1"
)
#==========================================================
# precompute the ward-specific basis settings
#
# These match your MMT calculation, where each ward used its own historical
# temperature distribution for the spline knots.
#==========================================================
var_knots_by_ward = vector("list", 69)
var_bounds_by_ward = vector("list", 69)

for (i in 1:69) {
  
  historical_temperature_i = df_complete %>% 
    st_drop_geometry() %>% 
    filter(new_id == i) %>% 
    pull(tasmean)
  
  var_knots_by_ward[[i]] = quantile(historical_temperature_i, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
  var_bounds_by_ward[[i]] = range(historical_temperature_i, na.rm = TRUE)
}

lag_knots = logknots(21, 3)
# Maximum lag used for heat-attributable mortality
attr_lags = 0:3
#------------------------------------------------------
#output folders
#use a new folder so the old total-temperature results are preserved
#------------------------------------------------------
out_dir = "output/cumlogRR_parts_2026_heatonly_newMMT/"
log_dir = "output/logs_2026_heatonly_newMMT/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# Remove previous files from this new output folder
old_files = list.files(out_dir, pattern = "^ward_[0-9]+_cumlogRR_heat\\.rds$", full.names = TRUE)
if (length(old_files) > 0) {
  file.remove(old_files)
}
#------------------------------------------------------
#task table
#------------------------------------------------------
task_df = ward_lookup %>% 
  transmute(i = new_id, Ward_Code)
cat("Total wards:", nrow(task_df), "\n")
#------------------------------------------------------
#start cluster
#------------------------------------------------------
n_cores = max(1, parallel::detectCores() - 14)
cl = parallel::makeCluster(n_cores)
doSNOW::registerDoSNOW(cl)
cat("Registered", n_cores, "cores\n")

pb = txtProgressBar(min = 0, max = nrow(task_df), style = 3)
progress = function(n) {
  setTxtProgressBar(pb, n)
}
opts = list(progress = progress)
start_time = Sys.time()
#==========================================================
# generate ward-specific heat-only cumulative log-RR
#==========================================================
task_files = foreach(
  t              = iterators::iter(task_df, by = "row"),
  .combine       = dplyr::bind_rows,
  .packages      = c("dlnm", "dplyr", "tibble"),
  .export        = c("cb_res", "MMT", "uk2026_heatwave_data",
                     "var_knots_by_ward", "var_bounds_by_ward",
                     "lag_knots", "attr_lags", "out_dir", "log_dir"),
  .options.snow  = opts,
  .errorhandling = "stop"
) %dopar% {
  
  i = as.integer(t$i)
  current_ward_code = as.character(t$Ward_Code)
  
  log_file = file.path(log_dir, sprintf("ward_%02d.log", i))
  cat(sprintf("[%s] START ward=%s\n", Sys.time(), current_ward_code),
      file = log_file, append = TRUE)
  #--------------------------------------------
  #posterior coefficient draws and corresponding MMT draws
  #--------------------------------------------
  beta_reg = cb_res[[i]]
  nsim = nrow(beta_reg)
  mmt_i = MMT[[i]] %>% 
    dplyr::arrange(nsim) %>% 
    dplyr::pull(MMT)
  stopifnot(length(mmt_i) == nsim)
  #--------------------------------------------
  #future observed temperatures
  #--------------------------------------------
  df_future_i = uk2026_heatwave_data %>% 
    dplyr::filter(Ward_Code == current_ward_code) %>% 
    dplyr::arrange(date)
  future_temp = df_future_i$tasmean
  N = length(future_temp)
  stopifnot(N > 21)
  #--------------------------------------------
  #backward lag-history matrix
  #
  # Row t:
  # lag 0  = temperature on day t
  # lag 1  = temperature on day t - 1
  # ...
  # lag 21 = temperature on day t - 21
  #--------------------------------------------
  at_mat = do.call(cbind, lapply(attr_lags, function(l) dplyr::lag(future_temp, n = l)))
  colnames(at_mat) = paste0("lag", attr_lags)
  ok = complete.cases(at_mat)
  #--------------------------------------------
  #crossbasis template
  #uses the ward's historical knots and historical boundary knots
  #--------------------------------------------
  cb_daily = crossbasis(
    future_temp,
    lag    = 21,
    argvar = list(fun = "bs", knots = var_knots_by_ward[[i]], Boundary.knots = var_bounds_by_ward[[i]]),
    arglag = list(fun = "ns", knots = lag_knots, intercept = TRUE)
  )
  stopifnot(ncol(cb_daily) == ncol(beta_reg))
  
  predvar = seq_len(N)
  predlag = attr_lags
  g = rep(seq_len(N), times = length(predlag))
  #--------------------------------------------
  #output matrix: valid outcome days x posterior draws
  #--------------------------------------------
  cumlogRR_heat = matrix(NA_real_, nrow = sum(ok), ncol = nsim)
  #--------------------------------------------
  #draw-specific heat-only prediction
  #--------------------------------------------
  for (j in seq_len(nsim)) {
    
    # Heat-only exposure history:
    # all lagged temperatures below the draw-specific MMT
    # are replaced by that MMT.
    #
    # Therefore, below-MMT exposure makes no heat contribution.
    at_heat_j = pmax(at_mat, mmt_i[j])
    
    Xpred_j = dlnm:::mkXpred(
      type    = "cb",
      basis   = cb_daily,
      at      = at_heat_j,
      predvar = predvar,
      predlag = predlag,
      cen     = mmt_i[j]
    )
    
    Xsum_j = rowsum(Xpred_j, group = g, reorder = FALSE)
    Xsum_j = Xsum_j[ok, , drop = FALSE]
    cumlogRR_heat[, j] = drop(Xsum_j %*% beta_reg[j, ])
  }
  #--------------------------------------------
  #diagnostics
  #--------------------------------------------
  min_value = min(cumlogRR_heat, na.rm = TRUE)
  proportion_negative = mean(cumlogRR_heat < 0, na.rm = TRUE)
  #--------------------------------------------
  #save matrix plus identifying information
  #--------------------------------------------
  output_object = list(
    Ward_id   = i,
    Ward_Code = current_ward_code,
    Date      = as.Date(df_future_i$date[ok]),
    tasmean   = future_temp[ok],
    cumlogRR  = cumlogRR_heat
  )
  
  output_file = file.path(out_dir, sprintf("ward_%02d_cumlogRR_heat.rds", i))
  saveRDS(output_object, output_file, compress = "gzip")
  
  cat(sprintf(paste0("[%s] DONE days=%d, draws=%d, ",
                     "minimum=%.6f, proportion_negative=%.4f\n"),
              Sys.time(), nrow(cumlogRR_heat), ncol(cumlogRR_heat),
              min_value, proportion_negative),
      file = log_file, append = TRUE)
  
  tibble::tibble(
    Ward_id             = i,
    Ward_Code           = current_ward_code,
    file                = output_file,
    min_cumlogRR        = min_value,
    proportion_negative = proportion_negative
  )
}
#------------------------------------------------------
#close cluster
#------------------------------------------------------
close(pb)
parallel::stopCluster(cl)
end_time = Sys.time()
cat("\nFinished. Runtime:\n")
print(end_time - start_time)
#------------------------------------------------------
#save diagnostics
#------------------------------------------------------
write_rds(task_files, "output/cumlogRR_heatonly_newMMT_diagnostic.rds")
print(
  task_files %>% 
    arrange(min_cumlogRR) %>% 
    slice_head(n = 15)
)
cat("\nNumber of completed ward files:",
    length(list.files(out_dir, pattern = "^ward_[0-9]+_cumlogRR_heat\\.rds$")),
    "\n")



