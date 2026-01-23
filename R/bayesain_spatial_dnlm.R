setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Spatial Bayesian DNLM Extreme heat and cold")

library(tidyverse)
library(readr)
library(dlnm)
library(mgcv)
library(gnm)
library(splines)

library("INLA")
library("tmap")
library("sf")
library("spdep")
library("readxl")

library(tsModel)

# Load data

first_stage_LAD_age_sepcific = readRDS("data/processed/model_data_ward22.rds")

################################################################################
#simply the data and arrange it by ward and date 
df = first_stage_LAD_age_sepcific %>% 
  select(ward22cd, date, deaths, tasmean) %>% 
  arrange(ward22cd, date)


# create lag matrix (includes lag 0)
lag_mat = tsModel::Lag(
  df$tasmean,
  group = df$ward22cd,
  k     = 0:21
)

lag_df = as.data.frame(lag_mat)
colnames(lag_df) <- paste0("lag", 0:21)

df = bind_cols(df, lag_df)

#remove lag0 
df = df %>% 
  select(-lag0)

###############################################################################
# #DROP incomplete lag rows before crossbasis
# df_complete = df %>%
#   group_by(ward22cd) %>%
#   filter(row_number() > 21) %>%   # removes first 21 days per ward
#   ungroup()


################################################################################
# Load Ward shp file
ward_map = read_sf("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp")

# 1) Build adjacency graph with poly2nb
ward_map = st_transform(ward_map, 27700)
ward_map= st_make_valid(ward_map)
which(!st_is_valid(ward_map))

# Create new IDs in order of data
ward_map$new_id = 1:nrow(ward_map)


# Export to INLA graph
mcnty_nb = poly2nb(ward_map, row.names = ward_map$new_id, queen = T)

nb2INLA("ward_map.adj", mcnty_nb )
g = inla.read.graph("ward_map.adj")


################################################################################

df  = df %>% left_join(ward_map %>% select(new_id, Ward_Code), by = c("ward22cd" = "Ward_Code")) %>% 
  arrange(new_id,date)

stopifnot(!is.unsorted(df$new_id))


cb_list = lapply(seq_len(max(df$new_id)), function(i) {
  
  # 1. Subset the data for the current ward
  d = df[df$new_id == i, ]
  
  # 2. Calculate local knots (This part was correct)
  temp_knots = quantile(
    d$tasmean,
    probs = c(0.1, 0.75, 0.9),
    na.rm = TRUE
  )
  
  temp_bounds = range(d$tasmean, na.rm = TRUE)
  
  # 3. Create the Cross-Basis using ONLY the subset data 'd'
  # [CRITICAL FIX]: Changed 'df$tasmean' to 'd$tasmean'
  cb = crossbasis(
    d$tasmean,     # <--- THIS WAS THE ERROR
    lag = 21,
    argvar = list(
      fun = "bs",
      knots = temp_knots,
      Boundary.knots = temp_bounds
    ),
    arglag = list(
      fun = "ns",
      knots = logknots(21, 3),
      intercept = TRUE # [NOTE]: Marco uses this. If you want to match, add it.
    )
  )
  
  as.data.frame(cb)
})

# Now this will bind correctly
cb_df = do.call(rbind, cb_list)



####################################################
colnames(cb_df) <- paste0("cb", seq_len(ncol(cb_df)))


stopifnot(nrow(cb_df) == nrow(df))

df = bind_cols(df, cb_df)



###################################################
#DROP incomplete lag rows before crossbasis
# df_complete = df %>%
#   group_by(ward22cd) %>%
#   filter(row_number() > 21) %>%   # removes first 21 days per ward
#   ungroup()
# 



###################################################
#create strata 

df_complete = df %>% 
  mutate(
    year = year(date),
    #this is the part to be conditioned out of the regression
    strata = factor(paste(ward22cd, year, formatC(month(date), width = 2, flag = "0"),
                          wday(date, week_start = 1),sep = "-"))
  )

#remove strata with 0 deaths 

df_complete = df_complete %>% 
  group_by(strata) %>% 
  filter(sum(deaths, na.rm = TRUE) > 0) %>%
  ungroup()


################################################################################
#Strata index (for conditional intercepts)
df_complete = df_complete %>%
  mutate(strata_id = as.integer(strata))

#Create spatial slope indices (Marcos β[a,c])
#One estimate per region for cb coefficients

n_cb = length(grep("^cb", names(df)))

for (j in seq_len(n_cb)) {
  df_complete[[paste0("ward_cb", j)]] <- df_complete$new_id
}


################################################################################
#...............................................................................
### MODEL 3 (SB-DLNM – case-crossover design) ####
#...............................................................................
### INLA - MODEL 3 SB-DLNM – case-crossover design

inla_formula = deaths ~ -1 + 
  cb1 + cb2 + cb3 + cb4 + cb5 + cb6 + 
  cb7 + cb8 + cb9 + cb10 + cb11 + cb12 +
  cb13 + cb14 +cb15 + cb16 + cb17 + cb18 + cb19 + cb20 +
  cb21 + cb22 + cb23 + cb24 + cb25 + cb26 + cb27 + cb28 + cb29 + cb30 +

  f(strata_id, model = "iid",
    hyper = list(prec = list(initial = log(1e-04), fixed = TRUE)))+

  f(ward_cb1, cb1, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb2, cb2, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb3, cb3, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb4, cb4, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb5, cb5, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb6, cb6, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb7, cb7, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb8, cb8, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb9, cb9, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb10, cb10, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb11, cb11, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb12, cb12, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb13, cb13, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb14, cb14, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb15, cb15, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb16, cb16, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb17, cb17, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb18, cb18, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb19, cb19, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb20, cb20, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb21, cb21, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb22, cb22, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb23, cb23, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb24, cb24, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb25, cb25, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb26, cb26, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb27, cb27, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb28, cb28, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb29, cb29, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) +
  
  f(ward_cb30, cb30, model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi  = list(prior = "pc", param = c(0.5, 2/3)))) 
  


# In order to sample from the posterior distribution model need to be fitted 
# with config = TRUE

start = Sys.time()

inla_model = inla(inla_formula,
                   data = df_complete, family = "poisson",
                   control.compute = list(config = TRUE,dic = TRUE, waic = TRUE),
                   control.inla = list(strategy = "gaussian",
                                       int.strategy = "eb"),
                  control.fixed = list(correlation.matrix=TRUE),
                  num.threads = round(parallel::detectCores()*0.4),
                  verbose = TRUE)

end = Sys.time()

end-start 


#################################################################################
# Extract the ensemble of sample coefficients of the crossbasis
inla_res = inla.posterior.sample(1000, inla_model,
                                  selection = list(cb1 = 1,
                                                   cb2 = 1,
                                                   cb3 = 1,
                                                   cb4 = 1,
                                                   cb5 = 1,
                                                   cb6 = 1,
                                                   cb7 = 1,
                                                   cb8 = 1,
                                                   cb9 = 1,
                                                   cb10 = 1,
                                                   cb11 = 1,
                                                   cb12 = 1,
                                                   cb13 = 1,
                                                   cb14 = 1,
                                                   cb15= 1,
                                                   cb16 =1,
                                                   cb17 =1,
                                                   cb18 =1,
                                                   cb19 =1,
                                                   cb20 =1,
                                                   cb21 =1,
                                                   cb22 =1,
                                                   cb23 =1,
                                                   cb24 =1,
                                                   cb25 =1,
                                                   cb26 =1,
                                                   cb27 =1,
                                                   cb28 =1,
                                                   cb29 =1,
                                                   cb30 =1,
                                                   "ward_cb1" = 1:69,
                                                   "ward_cb2" = 1:69,
                                                   "ward_cb3" = 1:69,
                                                   "ward_cb4" = 1:69,
                                                   "ward_cb5" = 1:69,
                                                   "ward_cb6" = 1:69,
                                                   "ward_cb7" = 1:69,
                                                   "ward_cb8" = 1:69,
                                                   "ward_cb9" = 1:69,
                                                   "ward_cb10" = 1:69,
                                                   "ward_cb11" = 1:69,
                                                   "ward_cb12" = 1:69,
                                                   "ward_cb13" = 1:69,
                                                   "ward_cb14" = 1:69,
                                                   "ward_cb15" = 1:69,
                                                   "ward_cb16" = 1:69,
                                                   "ward_cb17" = 1:69,
                                                   "ward_cb18" = 1:69,
                                                   "ward_cb19" = 1:69,
                                                   "ward_cb20" = 1:69,
                                                   "ward_cb21" = 1:69,
                                                   "ward_cb22" = 1:69,
                                                   "ward_cb23" = 1:69,
                                                   "ward_cb24" = 1:69,
                                                   "ward_cb25" = 1:69,
                                                   "ward_cb26" = 1:69,
                                                   "ward_cb27" = 1:69,
                                                   "ward_cb28" =1:69,
                                                   "ward_cb29" = 1:69,
                                                   "ward_cb30" = 1:69))
###################################################################
#look how it structured
s1 = inla_res[[1]]

str(s1$latent)


head(rownames(s1$latent), 50)

tail(rownames(s1$latent), 50)
###################################################################
#export the results
cb_res = vector("list", 69)

for (i in 1:69){
  
  beta_mat = matrix(NA, nrow = 1000, ncol = 30)
  
for (s in 1:1000){
  latent_s = inla_res[[s]]$latent
  
  for (k in 1:30){
    
    # global coefficient
    beta_global = latent_s[paste0("cb", k, ":1"), 1]
    
    # ward-specific deviation
    beta_ward = latent_s[paste0("ward_cb", k,":",i),1]
    
    # total coefficient
    beta_mat[s, k] = beta_global + beta_ward
  }
}
  
  cb_res[[i]] = beta_mat
  
}

cb_res[[1]]

write_rds(cb_res, file = "output/predicted_inla_spatial_casecrossover.rds")

write_rds(df_complete, file = "output/df_complete_for_inla.rds")



