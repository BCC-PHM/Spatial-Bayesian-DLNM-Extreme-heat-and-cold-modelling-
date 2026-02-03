setwd("~/R_project/Extreme heat and cold")
library(terra)
library(tidyterra)
library(dplyr)
library(stringr)
library(sf)
library(data.table)
library(readr)
###########################################
# --- PROJ setup (once) ---
#point R at the folder that contains proj.db

list.files(.libPaths(), pattern = "proj\\.db$", recursive = TRUE, full.names = TRUE)
# Sys.setenv(PROJ_LIB = "C:/Users/TMPACGAG/AppData/Local/Programs/R/R-4.4.2/library/sf/proj")
Sys.setenv(PROJ_LIB = "C:/Users/chung/AppData/Local/R/win-library/4.4/sf/proj/proj.db")
Sys.getenv("PROJ_LIB")

# --- Read Birmingham Ward shapefile ---

bham_sf = st_read("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp") # or .shp
#convert to terra for raster ops

bham =  terra::vect(bham_sf)

#before running the extraction 

bham_id_ward = bham_sf %>%
  st_drop_geometry() %>%
  mutate(ID = row_number()) %>%             
  select(ID, Ward_Code, Ward_Name) 

#Project Birmingham to British National Grid
bham_r = terra::project(bham, "EPSG:27700")

#============================================================
#data from https://zenodo.org/records/16213003
#set a list of path of the .nc files
# --- List NetCDF files (make sure ordering matches) ---

tas_file = sort(list.files(path = "data/raw/UKCP18_CPM_bias_corrected/", full.names = TRUE, pattern = ".nc"))

#initiate empty list 

output = list()

#loop over all files
for(i in seq_along(tas_file)){
  
  #--------------------------------------------------------
  #initiate the nc file
  
  tasfirst = rast(tas_file[i])
  
  # Assign British National Grid
  crs(tasfirst) = "EPSG:27700"
  
  #Pre-crop the entire 7,200-layer stack to Birmingham to save memory
  # This significantly speeds up the extraction process
  tas_bham_crop = crop(tasfirst, bham_r)
  
  # We use 'fun = mean' to get the average temperature for each ward polygon
  ext_raw = terra::extract(tas_bham_crop, bham_r, fun = mean, na.rm = TRUE)
  
  #--------------------------------------------------------
  #process the output
  
  #Get the dates from the raster metadata
  raster_dates = time(tasfirst)
  
  # We format them as "2020-12-01" etc.
  colnames(ext_raw)[2:ncol(ext_raw)] = as.character(raster_dates)
  
  #turn the data into long format and add the ward name id and member info back
  ext_raw =  ext_raw %>% 
    pivot_longer(cols = c(-ID),
                 names_to = "Date",
                 values_to = "tasmean") %>% 
    left_join(bham_id_ward, by = "ID") %>% 
    mutate(Member = as.character(str_extract(names(tasfirst)[1], "(?<=member=)\\d+(?=_)")))
  
  
  
  output[[i]] = ext_raw 
  
  
  
}
#============================================================

UKCPCPM_ward = rbindlist(output)

write_rds(UKCPCPM_ward, "data/processed/UKCP18CPM_tmean_ward_daily.rds")
