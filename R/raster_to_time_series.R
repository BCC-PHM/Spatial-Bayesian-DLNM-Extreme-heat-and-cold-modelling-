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
Sys.setenv(PROJ_LIB = "C:/Users/TMPACGAG/AppData/Local/Programs/R/R-4.4.2/library/sf/proj")
Sys.getenv("PROJ_LIB")

# --- Read Birmingham LSOA shapefile (once) ---
#Read Birmingham boundary (shapefile / geopackage)

bham_sf = st_read("data/external/boundaries/boundaries-lsoa-2021-birmingham/boundaries-lsoa-2021-birmingham.shp") # or .shp
#convert to terra for raster ops

bham =  terra::vect(bham_sf)
###########################################
# Stable lookup table (row order at time of extract)
#get the associated LSOA21CD for each row id
#before running the extraction 

bham_id_lsoa21 = bham_sf %>%
  st_drop_geometry() %>%
  mutate(ID = row_number()) %>%             
  select(ID, LSOA21CD, LSOA21NM) 

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


 for (i in seq_along(tmin_file)){
  
  message(
    format(Sys.time(), "%H:%M:%S"),
    " | Processing ", i, "/", length(tmin_file),
    " | ", basename(tmin_file[i]), " and ", basename(tmax_file[i])
  )
  
  
  
  tmin = rast(tmin_file[i])
  tmax = rast(tmax_file[i])

#extracting the date 
  
date_data = data.frame( date = as.Date(time(tmin)))
date_data = date_data %>% 
  mutate(layers_date_id = as.character(as.numeric(str_extract(date, "\\d\\d$"))))
  


# Crop + mask rasters to Birmingham
tmin_bham = mask(crop(tmin, bham_boundary), bham_boundary)
tmax_bham = mask(crop(tmax, bham_boundary), bham_boundary)


# Extract daily mean per LSOA (area-weighted)
# In terra, exact=TRUE means: for each polygon, it computes the coverage fraction of every raster cell 
# that intersects it (partial cells on the boundary get fractional weights), and when you supply fun = mean, 
# terra returns the area-weighted mean for each layer (each date).

tmin_ex = extract(tmin_bham, bham_r, fun = mean, na.rm = TRUE, exact = TRUE)
tmax_ex = extract(tmax_bham, bham_r, fun = mean, na.rm = TRUE, exact = TRUE)


# Compute tmean at LSOA/day level (wide)
#copies the whole tmin_ex to keep the same id 
tmean_ex = tmin_ex
#turn into matrix and calucalte the mean
tmean_ex[,-1] = (as.matrix(tmin_ex[,-1]) + as.matrix(tmax_ex[,-1])) / 2



##################################################################################

#rename the colnames
colnames(tmean_ex) = gsub("tasmin", "tasmean", colnames(tmean_ex))

#covert back to a dataframe and pviot it to long format 
tmean_long = as.data.frame(tmean_ex) %>% 
  left_join(bham_id_lsoa21, by = ("ID")) %>% 
  pivot_longer(cols = c(-ID, -LSOA21CD,-LSOA21NM),
               names_to = "layers",
               values_to = "tasmean"
  ) %>% 
  mutate(layers_date_id = str_extract(layers, "\\d+$")) %>% 
  left_join(date_data, by = ("layers_date_id"))

output[[i]] = tmean_long

message("  ✓ Done ", i, " (rows: ", nrow(tmean_long), ")")
  
}


#combine all the sublists inside output
all_combined = bind_rows(output)


#write the data to csv
write.csv(all_combined, file = "data/processed/Raster_to_time_Series.csv")





####################################################################################
n_cores = min(6, parallel::detectCores() - 2)
cl = makeCluster(n_cores)
registerDoParallel(cl)

#export only the shapefile path
bham_path = "data/external/boundaries/boundaries-lsoa-2021-birmingham/boundaries-lsoa-2021-birmingham.shp"


message(format(Sys.time(), "%H:%M:%S"), " | Cluster initiated with ", n_cores, " cores.")

output = foreach(
  i = seq_along(tmin_file),
  .packages = c("terra", "sf", "dplyr", "tidyr", "stringr"),
  .export   = c("bham_path", "bham_id_lsoa21")
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
  date_data = data.frame(date = as.Date(time(tmin))) %>%
    mutate(layers_date_id = as.character(as.numeric(str_extract(date, "\\d\\d$"))))
  
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
  
  # --- Mean ---
  tmean_ex = tmin_ex
  tmean_ex[, -1] = (as.matrix(tmin_ex[, -1]) + as.matrix(tmax_ex[, -1])) / 2
  
  colnames(tmean_ex) = gsub("tasmin", "tasmean", colnames(tmean_ex))
  
  # --- Long ---
  tmean_long = as.data.frame(tmean_ex) %>%
    left_join(bham_id_lsoa21, by = "ID") %>%
    pivot_longer(
      cols = c(-ID, -LSOA21CD, -LSOA21NM),
      names_to = "layers",
      values_to = "tasmean"
    ) %>%
    mutate(layers_date_id = str_extract(layers, "\\d+$")) %>%
    left_join(date_data, by = "layers_date_id")
  
  tmean_long
}



final_output = dplyr::bind_rows(output)

message(
  format(Sys.time(), "%H:%M:%S"),
  " | Processing complete. Total rows: ",
  nrow(final_output)
)

stopCluster(cl)
registerDoSEQ()










