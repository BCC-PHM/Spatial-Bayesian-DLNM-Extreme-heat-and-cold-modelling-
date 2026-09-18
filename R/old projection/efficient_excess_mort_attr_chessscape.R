library(tidyverse)
library(dlnm)
library(tsModel)
library(doParallel)
library(doSNOW)
library(foreach)
library(iterators)
library(fst)
library(RhpcBLASctl)
library(sf)

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

UKCP18CM_ward = read_rds("data/processed/chessscape65_daily_ward.rds")  #change your corresponding file for different rcp scenrio here 

#----------------------------------------------------

# Create a combined grouping factor based on both 'ward22cd' and 'Member'

UKCP18CM_ward = UKCP18CM_ward %>%
  arrange(Ward_Code, Member, Date)

group_factor = interaction(UKCP18CM_ward$Ward_Code, UKCP18CM_ward$Member, drop = TRUE)

UKCP18CM_ward_lags = tsModel::Lag(
  UKCP18CM_ward$tasmean,
  group = group_factor,
  k = 0:21
)

UKCP18CM_ward = bind_cols(
  UKCP18CM_ward,
  setNames(as.data.frame(UKCP18CM_ward_lags), paste0("lag", 0:21))
)


# ---------------------------
# BLAS threads = 1 (IMPORTANT)
# ---------------------------
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

# ---------------------------
# Precompute knots ONCE
# ---------------------------
var_knots = quantile(df_complete$tasmean, probs = c(0.1, 0.75, 0.9), na.rm = TRUE)
lag_knots = logknots(21, 3)

# ---------------------------
# Output folders
# ---------------------------
out_dir = "output/cumulRR_parts_chesscape65"
log_dir = "output/logs_ukcp18"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

# (Optional) drop sf geometry if df_complete ever becomes sf
if (inherits(df_complete, "sf")) df_complete = sf::st_drop_geometry(df_complete)

#--------------------------------------------------------------------------
#build task table 

ward_codes = unique(df_complete$ward22cd)  # 69

task_df = map_dfr(seq_along(ward_codes), function(i) {
  w = ward_codes[i]
  members = unique(UKCP18CM_ward$Member[UKCP18CM_ward$Ward_Code == w])
  tibble(i = i, Ward_Code = w, Member = members)
})

cat("Total tasks (ward-member combos):", nrow(task_df), "\n")

#--------------------------------------------------------------------------
n_cores = max(1, parallel::detectCores() - 14)
cl = parallel::makeCluster(n_cores)
doSNOW::registerDoSNOW(cl)
cat("Registered", n_cores, "cores\n")

pb = txtProgressBar(min = 0, max = nrow(task_df), style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

start = Sys.time()

task_files = foreach(
  t = iterators::iter(task_df, by = "row"),
  .packages = c("dlnm", "dplyr"),
  .export   = c("cb_res", "MMT", "UKCP18CM_ward",
                "var_knots", "lag_knots", "out_dir", "log_dir"),
  .options.snow = opts,
  .errorhandling = "stop"
) %dopar% {
  
  i = t$i
  w = t$Ward_Code
  m = t$Member
  
  log_file = file.path(log_dir, sprintf("ward_%02d_member_%s.log", i, m))
  cat(sprintf("[%s] START (ward=%s member=%s)\n", Sys.time(), w, m),
      file = log_file, append = TRUE)
  
  # --- posterior coefs + MMT ---
  beta_reg = cb_res[[i]]               # [nsim x K]
  nsim = nrow(beta_reg)
  
  mmt_i = MMT[[i]] %>%
    dplyr::arrange(nsim) %>%
    dplyr::pull(MMT)
  
  stopifnot(length(mmt_i) == nsim)
  
  # --- future temps for ward-member ---
  df_future_w = UKCP18CM_ward %>%
    dplyr::filter(Ward_Code == w, Member == m) %>%
    dplyr::arrange(Date)
  
  future_temp = df_future_w$tasmean
  
  # use precomputed lag columns (no tsModel inside)
  lag_cols = paste0("lag", 0:21)
  if (!all(lag_cols %in% names(df_future_w))) {
    stop(sprintf(
      "Missing lag0..lag21 for ward=%s member=%s. Precompute lags BEFORE cluster.",
      w, m
    ))
  }
  at_mat = as.matrix(df_future_w[, lag_cols])
  
  
  
  N = nrow(at_mat)
  ok = seq_len(N) > 21 & complete.cases(at_mat)
  
  # --- crossbasis (same spec, knots precomputed) ---
  cb_daily = crossbasis(
    future_temp,
    lag    = 21,
    argvar = list(fun = "bs", knots = var_knots),
    arglag = list(fun = "ns", knots = lag_knots)
  )
  
  stopifnot(ncol(cb_daily) == ncol(beta_reg))
  
  predvar = seq_len(N)
  predlag = dlnm:::seqlag(c(0, 21))
  
  # allocate
  cumlogRR = matrix(NA_real_, nrow = sum(ok), ncol = nsim)
  mode(cumlogRR) = "single"
  
  # grouping for rowsum to collapse stacked lags back to N rows
  g = rep(seq_len(N), times = length(predlag))
  
  for (j in seq_len(nsim)) {
    
    Xpred_j = dlnm:::mkXpred(
      type     = "cb",
      basis    = cb_daily,
      at       = at_mat,
      predvar  = predvar,
      predlag  = predlag,
      cen      = mmt_i[j]
    )
    
    Xsum_j = rowsum(Xpred_j, group = g, reorder = FALSE)  # N x K
    Xsum_j = Xsum_j[ok, , drop = FALSE]
    
    cumlogRR[, j] = drop(Xsum_j %*% beta_reg[j, ])
  }
  
  # save per task
  f = file.path(out_dir, sprintf("ward_%02d_member_%s_cumlogRR.rds", i, m))
  saveRDS(cumlogRR, f, compress = "gzip")
  
  
  cat(sprintf("[%s] DONE (n_days=%d, nsim=%d) -> %s\n",
              Sys.time(), nrow(cumlogRR), ncol(cumlogRR), basename(f)),
      file = log_file, append = TRUE)
  
  f
}

close(pb)
parallel::stopCluster(cl)

end = Sys.time()
cat("\nFinished. Runtime:\n")
print(end - start)



#------------------------------------------------------------------------------
#failed at 224

# task_df_resume = task_df[224:nrow(task_df), , drop = FALSE]

#keeping failing
#just to filter out the remaining task to run on
# 
# files =list.files(path = "output/cumulRR_parts_chesscape65/", pattern = ".rds", full.names = TRUE)
# 
# ward_id = as.numeric(str_extract(files, "(?<=ward_)\\d+"))
# 
# member_id = str_extract(files, "(?<=member_)\\d+")
# 
# keys = paste0(ward_id, "_", member_id)
# 
# finished_datset = data.frame(finished = keys)
# finished_datset = finished_datset %>% 
#   mutate(Status = "finished")
# 
# task_df_resume  = task_df %>% 
#   mutate(list_id = paste0(i,"_",Member)) %>% 
#   left_join(finished_datset, by = c("list_id" = "finished")) %>% 
#   filter(is.na(Status))
# 
# 
# # --- (Re)make cluster ---
# n_cores = min(6, max(1, parallel::detectCores() - 14))
# cl = parallel::makeCluster(n_cores)
# doSNOW::registerDoSNOW(cl)
# cat("Registered", n_cores, "cores\n")
# pb = txtProgressBar(min = 0, max = nrow(task_df_resume), style = 3)
# progress = function(n) setTxtProgressBar(pb, n)
# opts = list(progress = progress)
# 
# task_files_resume = foreach(
#   t = iterators::iter(task_df_resume, by = "row"),
#   .packages = c("dlnm", "dplyr"),
#   .export   = c("cb_res", "MMT", "UKCP18CM_ward", "var_knots", "lag_knots", "out_dir", "log_dir"),
#   .options.snow = opts,
#   .errorhandling = "stop"
# ) %dopar% {
#   
#   i = t$i
#   w = t$Ward_Code
#   m = t$Member
#   
#   log_file = file.path(log_dir, sprintf("ward_%02d_member_%s.log", i, m))
#   cat(sprintf("[%s] START (ward=%s member=%s)\n", Sys.time(), w, m),
#       file = log_file, append = TRUE)
#   
#   beta_reg = cb_res[[i]]
#   nsim = nrow(beta_reg)
#   
#   mmt_i = MMT[[i]] %>% arrange(nsim) %>% pull(MMT)
#   stopifnot(length(mmt_i) == nsim)
#   
#   df_future_w = UKCP18CM_ward %>%
#     filter(Ward_Code == w, Member == m) %>%
#     arrange(Date)
#   
#   future_temp = df_future_w$tasmean
#   
#   lag_cols = paste0("lag", 0:21)
#   if (!all(lag_cols %in% names(df_future_w))) {
#     stop(sprintf("Missing lag0..lag21 for ward=%s member=%s.", w, m))
#   }
#   at_mat = as.matrix(df_future_w[, lag_cols])
#   
#   N = nrow(at_mat)
#   ok = seq_len(N) > 21 & complete.cases(at_mat)
#   
#   cb_daily = crossbasis(
#     future_temp,
#     lag    = 21,
#     argvar = list(fun = "bs", knots = var_knots),
#     arglag = list(fun = "ns", knots = lag_knots)
#   )
#   stopifnot(ncol(cb_daily) == ncol(beta_reg))
#   
#   predvar = seq_len(N)
#   predlag = dlnm:::seqlag(c(0, 21))
#   g = rep(seq_len(N), times = length(predlag))
#   
#   cumlogRR = matrix(NA_real_, nrow = sum(ok), ncol = nsim)
#   mode(cumlogRR) = "single"
#   
#   for (j in seq_len(nsim)) {
#     Xpred_j = dlnm:::mkXpred(
#       type     = "cb",
#       basis    = cb_daily,
#       at       = at_mat,
#       predvar  = predvar,
#       predlag  = predlag,
#       cen      = mmt_i[j]
#     )
#     
#     Xsum_j = rowsum(Xpred_j, group = g, reorder = FALSE)
#     Xsum_j = Xsum_j[ok, , drop = FALSE]
#     cumlogRR[, j] = drop(Xsum_j %*% beta_reg[j, ])
#   }
#   
#   f   = file.path(out_dir, sprintf("ward_%02d_member_%s_cumlogRR.rds", i, m))
#   tmp = paste0(f, ".tmp")
#   
#   if (file.exists(tmp)) file.remove(tmp)
#   saveRDS(cumlogRR, tmp, compress = "gzip")
#   if (!file.rename(tmp, f)) stop(sprintf("Rename failed for %s", basename(f)))
#   
#   cat(sprintf("[%s] DONE (n_days=%d, nsim=%d) -> %s\n",
#               Sys.time(), nrow(cumlogRR), ncol(cumlogRR), basename(f)),
#       file = log_file, append = TRUE)
#   
#   NULL  # =- key: don't return anything meaningful
# }
# 
# close(pb)
# parallel::stopCluster(cl)
# 
# 
# end = Sys.time()
# cat("\nFinished. Runtime:\n")
# print(end - start)
# 
# 


#-----------------------------------------------------------------------------
# files = unlist(task_files)
# 
# # sort nicely by ward then member
# ward_id   = as.integer(sub(".*ward_(\\d+)_.*", "\\1", basename(files)))
# member_id = sub(".*member_([^_]+)_.*", "\\1", basename(files))
# o = order(ward_id, member_id)
# files = files[o]
# 
# ward_specific_cumul_RR_UKCP18 = lapply(files, readRDS)
# names(ward_specific_cumul_RR_UKCP18) = paste0(ward_id[o], "_", member_id[o])
# 
# saveRDS(ward_specific_cumul_RR_UKCP18,
#         "output/ward_specific_cumul_CHESSSCAPE_RCP85.rds",
#         compress = "gzip")

# or RDS:
# saveRDS(ward_specific_cumul_RR_UKCP18, "output/ward_specific_cumul_CHESSSCAPE_RCP85.rds", compress="gzip")












