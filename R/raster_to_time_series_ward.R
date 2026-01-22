setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Extreme heat and cold")

library(terra)
library(tidyterra)
library(dplyr)
library(stringr)
library(doParallel)
library(foreach)
library(doSNOW)

###########################################
# --- PROJ setup (once) ---
#point R at the folder that contains proj.db

list.files(.libPaths(), pattern = "proj\\.db$", recursive = TRUE, full.names = TRUE)
# Sys.setenv(PROJ_LIB = "C:/Users/TMPACGAG/AppData/Local/Programs/R/R-4.4.2/library/sf/proj")
Sys.setenv(PROJ_LIB = "C:/Users/chung/AppData/Local/R/win-library/4.4/sf/proj/proj.db")
Sys.getenv("PROJ_LIB")

# --- Read Birmingham LSOA shapefile (once) ---
#Read Birmingham boundary (shapefile / geopackage)
# --- Read Birmingham LSOA shapefile (once) ---
#Read Birmingham boundary (shapefile / geopackage)

bham_sf = st_read("data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp") # or .shp
#convert to terra for raster ops

bham =  terra::vect(bham_sf)
###########################################
###########################################
# Stable lookup table (row order at time of extract)
#get the associated LSOA21CD for each row id
#before running the extraction 

bham_id_ward = bham_sf %>%
  st_drop_geometry() %>%
  mutate(ID = row_number()) %>%             
  select(ID, Ward_Code, Ward_Name) 
###########################################
#set a list of path of the .nc files
# --- List NetCDF files (make sure ordering matches) ---

tmin_file = sort(list.files(path = "data/raw/tasmin/", full.names = TRUE, pattern = ".nc"))
tmax_file = sort(list.files(path = "data/raw/tasmax/", full.names = TRUE, pattern = ".nc"))

# --- Use first raster to set CRS + bham reprojection ONCE (assuming all files same CRS) ---

tmin0 = rast(tmin_file[1])

#Reproject LSOAs to match raster CRS

bham_r = terra::project(bham, crs(tmin0))

# Dissolve to single boundary for fast crop/mask

bham_boundary = aggregate(bham_r)


###########################################
# Initialize empty lists to hold the output 
output = vector("list", length(tmin_file))


###########################################
n_cores = min(6, parallel::detectCores() - 2)
cl = makeCluster(n_cores)
registerDoParallel(cl)

#export only the shapefile path
bham_path = "data/external/boundaries/boundaries-wards-2022-birmingham/boundaries-wards-2022-birmingham.shp"


message(format(Sys.time(), "%H:%M:%S"), " | Cluster initiated with ", n_cores, " cores.")

output = foreach(
  i = seq_along(tmin_file),
  .packages = c("terra", "sf", "dplyr", "tidyr", "stringr"),
  .export = c("bham_path", "bham_id_ward")
) %dopar% {
  
  # --- Read rasters ---
  tmin = rast(tmin_file[i])
  tmax = rast(tmax_file[i])
  
  # --- Read boundary INSIDE worker ---
  bham_sf  = sf::st_read(bham_path, quiet = TRUE)
  bham     = terra::vect(bham_sf)
  bham_r   = terra::project(bham, crs(tmin))
  bham_boundary = terra::aggregate(bham_r)
  
  # --- Dates ---
  date_data = data.frame(
    layer_id = seq_along(time(tmin)),
    date     = as.Date(time(tmin))
  )
  
  # --- Crop + mask ---
  tmin_bham = mask(crop(tmin, bham_boundary), bham_boundary)
  tmax_bham = mask(crop(tmax, bham_boundary), bham_boundary)
  
  # --- Extract ---
  tmin_ex = terra::extract(
    tmin_bham, bham_r,
    fun = mean, na.rm = TRUE, exact = TRUE
  )
  
  tmax_ex = terra::extract(
    tmax_bham, bham_r,
    fun = mean, na.rm = TRUE, exact = TRUE
  )
  
  # --- Long tasmin ---
  tasmin_long = as.data.frame(tmin_ex) %>%
    left_join(bham_id_ward, by = "ID") %>%
    pivot_longer(
      cols = -c(ID, Ward_Code, Ward_Name),
      names_to  = "layers",
      values_to = "tasmin"
    ) %>%
    mutate(layer_id = as.integer(str_extract(layers, "\\d+$"))) %>%
    select(-layers)
  
  # --- Long tasmax ---
  tasmax_long = as.data.frame(tmax_ex) %>%
    left_join(bham_id_ward, by = "ID") %>%
    pivot_longer(
      cols = -c(ID, Ward_Code, Ward_Name),
      names_to  = "layers",
      values_to = "tasmax"
    ) %>%
    mutate(layer_id = as.integer(str_extract(layers, "\\d+$"))) %>%
    select(-layers)
  
  # --- Combine + mean ---
  out = tasmin_long %>%
    left_join(
      tasmax_long,
      by = c("ID", "Ward_Code", "Ward_Name", "layer_id")
    ) %>%
    left_join(date_data, by = "layer_id") %>%
    mutate(tasmean = (tasmin + tasmax) / 2)
  
  out
}


final_output = dplyr::bind_rows(output)

message(
  format(Sys.time(), "%H:%M:%S"),
  " | Processing complete. Total rows: ",
  nrow(final_output)
)

stopCluster(cl)
registerDoSEQ()




write_rds(final_output, "data/processed/tmean_ward_daily.rds")

























