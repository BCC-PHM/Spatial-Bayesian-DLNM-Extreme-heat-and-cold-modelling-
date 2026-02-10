setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Spatial Bayesian DNLM Extreme heat and cold")
setwd("~/R_project/Extreme heat and cold")
library(readr)
library(tidyverse)
library(readxl)
library(parallel)
library(foreach)
#the file that store the individual outputs 
parts_dir = "output/cumulRR_parts_chessscape_85/"

#Finds all files ending in .rds
files = list.files(parts_dir, pattern = ".rds", full.names = TRUE)



# parse ward + member from filename
ward_id   = as.numeric(str_extract(files, "(?<=ward_)\\d+"))
member_id = str_extract(files, "(?<=member_)\\d+")

# create keys like wardid_memberid e.g. ("1_1")
keys = paste0(ward_id, "_", member_id)

idx = tibble(key = keys, file = files) %>%
  arrange(ward_id, member_id)

# saveRDS(idx, "output/cumlogRR_file_index.rds", compress = "gzip")
#---------------------------------------------------------------------
#Create the named vector
ward_specific_cumul_RR_chessscape = stats::setNames(idx$file, idx$key)



#Define the automatic reader
`[[.filebacked_list` = function(x, i) {
  readRDS(NextMethod())
}

#Apply the label
class(ward_specific_cumul_RR_chessscape) = "filebacked_list"

ward_specific_cumul_RR_chessscape[["1_6"]]
#==========================================================================

#load mmt resutls 
MMT = read_rds("output/mmt_draws_by_ward.rds")

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#load baseline mortality by ward
baseline_mort = read_rds("output/baseline_mort_ward.rds")

#prepare the baseline mortality for later loop
mort_grid = expand.grid(ID=seq(1,69,by=1),Member=c(1,15,4,6))

baseline_mort = baseline_mort %>% 
  left_join(mort_grid, by = c("new_id" = "ID")) %>% 
  mutate(list_id = paste0(new_id, "_", Member))

#--------------------------------------------------------------------------
#load the daily future time series for each ward
chessscapeCM_ward = read_rds("data/processed/chessscape85_daily_ward.rds")

chessscapeCM_ward = chessscapeCM_ward %>% 
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




#=====================================================================================
# Backward attributable fraction (AF) and AN

keys = unique(chessscapeCM_ward$list_id)

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
  
  AF_calculation = 1-exp(-ward_specific_cumul_RR_chessscape[[i]])
  
  baseline_mort_value = baseline_mort$mean_baseline[baseline_mort$list_id==i]
  
  AN_calculation = AF_calculation*baseline_mort_value
  
  # 
  # df_i =chessscapeCM_ward %>% 
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
rm(AF_calculation)
gc()
#===========================================================================
#create matrix to define heat or cold days base on ward specific MMTsfor each 
#posterior distribution for each future day 
#===========================================================================
MMT_table = data.table::rbindlist(MMT) %>% 
  left_join(mort_grid, by = c("Ward_id" = "ID"),relationship = "many-to-many") %>% 
  mutate(list_id = paste0(Ward_id, "_", Member))



# Backward attributable fraction (AF) and AN

keys = unique(chessscapeCM_ward$list_id)

temp_daily = vector("list", length(keys))
names(temp_daily ) = keys   # so temp_daily[[i]] works with character keys


for (i in seq_along(keys)){
  m = keys[i]
  
  cat("\n-----------------------------------------\n")
  cat(sprintf("[%s] %d / %d  Processing %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              i, length(keys), m))
  
  df_im = chessscapeCM_ward[chessscapeCM_ward$list_id == m, ]
  
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
keys = unique(chessscapeCM_ward$list_id)

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
  df_im = chessscapeCM_ward %>% 
    filter(list_id == m) %>% 
    mutate(Decade = case_when(
      Date >= as.Date("2030-01-01") & Date <= as.Date("2039-12-31") ~ "2030s",
      Date >= as.Date("2040-01-01") & Date <= as.Date("2049-12-31") ~ "2040s",
      Date >= as.Date("2050-01-01") & Date <= as.Date("2059-12-31") ~ "2050s",
      Date >= as.Date("2060-01-01") & Date <= as.Date("2069-12-31") ~ "2060s",
      Date >= as.Date("2070-01-01") & Date <= as.Date("2079-12-31") ~ "2070s",
      TRUE ~ NA_character_
    ))
  
  
  #apply the keep index to the corresponding AN daily and temp class
  AN_mat =as.matrix(AN_daily[[m]])        
  
  cls_mat = temp_daily[[m]]   
  
  #make sure matching row
  stopifnot(nrow(df_im) == nrow(AN_mat), nrow(df_im) == nrow(cls_mat))
  
  
  #filter out heat and cold 
  AN_heat = AN_mat
  AN_heat[cls_mat != "heat"] = 0
  
  AN_cold = AN_mat
  AN_cold[cls_mat != "cold"] = 0
  
  
  #add date back
  heat = cbind(df_im, as.data.frame(AN_heat))
  
  cold = cbind(df_im, as.data.frame(AN_cold))
  
  
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
  left_join(heat_long, by = c("ID", "Ward_Code","Decade")) %>% 
  mutate(emission = "RCP8.5")

write.csv(heath_and_cold_future_long, "output/RCP85_heat_cold_final.csv")

















