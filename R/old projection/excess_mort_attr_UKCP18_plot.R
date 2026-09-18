setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Spatial Bayesian DNLM Extreme heat and cold")
setwd("~/R_project/Extreme heat and cold")
library(readr)
library(tidyverse)
library(readxl)

#
ward_specific_cumul_RR_UKCP18 = read_rds("output/ward_specific_cumul_RR_UKCP18.rds")

#-------------------------------------------------------------------------
#load mmt resutls 
MMT = read_rds("output/mmt_draws_by_ward.rds")




#--------------------------------------------------------------------------
#load baseline mortality by ward
baseline_mort = read_rds("output/baseline_mort_ward.rds")

#prepare the baseline mortality for later loop
mort_grid = expand.grid(ID=seq(1,69,by=1),Member=c(1,4,7,8))

baseline_mort = baseline_mort %>% 
  left_join(mort_grid, by = c("new_id" = "ID")) %>% 
  mutate(list_id = paste0(new_id, "_", Member))

#--------------------------------------------------------------------------
#load the daily future time series for each ward
UKCP18CM_ward = read_rds("data/processed/UKCP18CPM_tmean_ward_daily.rds")

# create lag matrix (includes lag 0)

UKCP18CM_ward_lags = tsModel::Lag(
  UKCP18CM_ward$tasmean,
  group = interaction(UKCP18CM_ward$Ward_Code , UKCP18CM_ward$Member),
  k     = 0:21)

# Name columns for clarity

UKCP18CM_ward_lags = as.data.frame(UKCP18CM_ward_lags)
colnames(UKCP18CM_ward_lags) = paste0("lag", 0:21)


UKCP18CM_ward = bind_cols(UKCP18CM_ward, UKCP18CM_ward_lags)

#remove lag0 
UKCP18CM_ward = UKCP18CM_ward %>% 
  select(-lag0)


# Ensure we remove the first 21 days due to incomplete lag history
UKCP18CM_ward = UKCP18CM_ward[seq_len(nrow(UKCP18CM_ward)) > 21, ]

#sine in the loop i called ts model again so i lagged the first set 1-1 twice
#but it is not a problem here because i dont need the 2020 data anyways 
#therefore below is to make sure it macthes the output to cbind


UKCP18CM_ward_adjusted = UKCP18CM_ward %>% 
  mutate(
    Date  = as.Date(Date),
    group = interaction(Ward_Code,Member)) %>% 
  group_by(group) %>% 
  mutate(
    idx_in_group = row_number()) %>%
  filter(idx_in_group > 21) %>% 
  ungroup() %>% 
  select(ID, Date,tasmean,Ward_Code,Ward_Name,Member) %>% 
  mutate(list_id = paste0(ID,"_",Member))



# Backward attributable fraction (AF) and AN
 
keys = unique(UKCP18CM_ward_adjusted$list_id)

AN_daily = vector("list", length(keys))
names(AN_daily) = keys   # iso AF_daily[[i]] works with character keys

start <- Sys.time()
cat("Starting AN calculation for", length(keys), "ward-member groups at", format(start), "\n")


  for (k in seq_along(keys)){
    
    
    i = keys[k]
    
    cat("\n-----------------------------------------\n")
    cat(sprintf("[%s] %d / %d  Processing %s\n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                k, length(keys), i))
    
    AF_calculation = 1-exp(-ward_specific_cumul_RR_UKCP18[[i]])
    
    baseline_mort_value = baseline_mort$mean_baseline[baseline_mort$list_id==i]
    
    AN_calculation = AF_calculation*baseline_mort_value
    
    # 
    # df_i =UKCP18CM_ward_adjusted %>% 
    #             filter(list_id == i)
    # 
    # combined = cbind(df_i, AN_calculation)
    
    AN_daily[[i]] = AN_calculation
    
    cat(sprintf("DONE %s | baseline=%.4f | rows=%d | sims=%d\n",
                i, baseline_mort_value, nrow(AN_calculation), ncol(AN_calculation)))
  }

end <- Sys.time()
cat("\nAll done at", format(end), " | Runtime:", as.character(end - start), "\n")

#---------------------------------------------------------------------------  


#===========================================================================
#create matrix to define heat or cold days base on ward specific MMTsfor each 
#posterior distribution for each future day 
#===========================================================================
MMT_table = data.table::rbindlist(MMT) %>% 
  left_join(mort_grid, by = c("Ward_id" = "ID"),relationship = "many-to-many") %>% 
  mutate(list_id = paste0(Ward_id, "_", Member))



# Backward attributable fraction (AF) and AN

keys = unique(UKCP18CM_ward_adjusted$list_id)

temp_daily = vector("list", length(keys))
names(temp_daily ) = keys   # so temp_daily[[i]] works with character keys


for (i in seq_along(keys)){
m = keys[i]

cat("\n-----------------------------------------\n")
cat(sprintf("[%s] %d / %d  Processing %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            i, length(keys), m))

df_im = UKCP18CM_ward_adjusted[UKCP18CM_ward_adjusted$list_id == m, ]

mmt_im = MMT_table[MMT_table$list_id ==m,] %>%
  pull(MMT)


# temperature classification per posterior draw
temp_class_mat = matrix(NA_character_, nrow = nrow(df_im), ncol = 1000)

for(n in 1:1000){

temp_class_mat[,n] = ifelse(df_im$tasmean > mmt_im[n],
                            "heat",
                            ifelse(df_im$tasmean < mmt_im[n], "cold", "atMMT"))
}

temp_daily[[m]] = temp_class_mat

cat(sprintf("DONE %s | rows=%d | sims=%d\n",
            m, nrow(df_im), ncol(temp_class_mat)))

}

#------------------------------------------------------------------
#to filter out and separate for heat and cold data set
#------------------------------------------------------------------
#specifiy the key
keys = unique(UKCP18CM_ward_adjusted$list_id)

#initialise the empty dataset fot heat
heat_list = vector("list", length(keys))
names(heat_list) = keys

#initialise the empty dataset fot cold
cold_list = vector("list", length(keys))
names(cold_list) = keys


for (i in seq_along(keys)){
  m = keys[i]
  
  cat("\n-----------------------------------------\n")
  cat(sprintf("[%s] %d / %d  Processing %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              i, length(keys), m))

#just to add decade 
df_im = UKCP18CM_ward_adjusted %>% 
  filter(list_id == m) %>% 
  mutate(Decade = case_when(
    Date >= as.Date("2030-12-01") & Date <= as.Date("2040-11-30") ~ "2030s",
    Date >= as.Date("2060-12-01") & Date <= as.Date("2070-11-30") ~ "2060s",
    Date >= as.Date("2070-12-01") & Date <= as.Date("2080-11-30") ~ "2070s",
    TRUE ~ NA_character_
  ))


# rows to drop (2060-12-01 to 2060-12-21)
drop_idx = which(df_im$Date >= as.Date("2060-12-01") &
                   df_im$Date <= as.Date("2060-12-21"))

#define the lkeep index
keep_idx = setdiff(seq_len(nrow(df_im)), drop_idx)


  
#apply the keep index to the corresponding AN daily and temp class
AN_mat =as.matrix(AN_daily[[m]])        

cls_mat = temp_daily[[m]]   

#make sure matching row
stopifnot(nrow(df_im) == nrow(AN_mat), nrow(df_im) == nrow(cls_mat))

#now use the keep index to filter the dataset, AN and temp class
df_im2  = df_im[keep_idx, , drop = FALSE]
AN2    = AN_mat[keep_idx, , drop = FALSE]
cls2   = cls_mat[keep_idx, , drop = FALSE]


#filter out heat and cold 
AN_heat = AN2
AN_heat[cls2 != "heat"] = 0

AN_cold = AN2
AN_cold[cls2 != "cold"] = 0


#add date back
heat = cbind(df_im2, as.data.frame(AN_heat))

cold = cbind(df_im2, as.data.frame(AN_cold))


#aggreate and calculate annual median for heat and cold 


check2 = heat %>% 
  filter(!is.na(Decade)) %>% 
  mutate(year = year(Date)) %>% 
  group_by(ID, Ward_Code, Ward_Name, Member, list_id, Decade, year) %>% 
  summarise(across(where(is.numeric),~ sum(.x, na.rm = TRUE))) %>% 
  group_by(ID, Ward_Code, Ward_Name, Member, list_id, Decade) %>% 
  summarise(across(where(is.numeric),~ mean(.x, na.rm = TRUE)))


check3 = cold %>% 
  filter(!is.na(Decade)) %>% 
  mutate(year = year(Date)) %>% 
  group_by(ID, Ward_Code, Ward_Name, Member, list_id, Decade, year) %>% 
  summarise(across(where(is.numeric),~ sum(.x, na.rm = TRUE))) %>% 
  group_by(ID, Ward_Code, Ward_Name, Member, list_id, Decade) %>% 
  summarise(across(where(is.numeric),~ mean(.x, na.rm = TRUE)))

#store them 
heat_list[[m]] = check2

cold_list[[m]] = check3

cat(sprintf("DONE %s | decades=%d | posterior draws=%d\n",
            m,
            nrow(check3),
            sum(grepl("^V", names(check3)))))

}


#------------------------------------------------------------------
#process the output from previos step
#-----------------------------------------------------------------

cold_long = data.table::rbindlist(cold_list) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "draw",
    values_to = "AN"
  ) %>%
  group_by(ID, Ward_Code, Ward_Name, Decade) %>% 
  mutate(draw = as.integer(sub("V", "", draw)))%>%  # pool all the ensemble 
  summarise(
    cold_LL = quantile(AN, 0.025),
    cold_med  = quantile(AN, 0.5),
    cold_UL = quantile(AN, 0.975),
    .groups = "drop"
  )




heat_long = data.table::rbindlist(heat_list) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "draw",
    values_to = "AN"
  ) %>%
  group_by(ID, Ward_Code, Ward_Name, Decade) %>% 
  mutate(draw = as.integer(sub("V", "", draw)))%>% # pool all the ensemble 
  summarise(
    heat_LL = quantile(AN, 0.025),
    heat_med  = quantile(AN, 0.5),
    heat_UL = quantile(AN, 0.975),
    .groups = "drop"
  )


#merge heat and cold 

heath_and_cold_future_long = cold_long %>% 
  left_join(heat_long, by = c("ID", "Ward_Code","Decade"))


write_rds(heath_and_cold_future_long, "output/heat_and_cold_related_EM_ward_UKCP18.rds")

##########################################################################
#====================================================================
#calculate per 100,000 rate
#--------------------------------------------------------------------

heath_and_cold_future_long 



ward_pop = read_excel("data/external/population/sapewardstablefinal.xlsx", 
                      sheet = "Mid-2022 Ward 2022", skip = 3)














