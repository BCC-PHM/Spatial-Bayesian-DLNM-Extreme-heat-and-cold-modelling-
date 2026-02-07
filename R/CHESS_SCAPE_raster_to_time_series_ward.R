# setwd("~/R_project/Extreme heat and cold")
setwd("C:/Users/TMPACGAG/OneDrive - Birmingham City Council/Documents/R projects/PHM/Spatial Bayesian DNLM Extreme heat and cold")

library(terra)
library(tidyterra)
library(dplyr)
library(stringr)
library(sf)
library(data.table)
library(readr)
library(foreach)
library(doParallel)
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


#-----------------------------------------------------------
#for the use of parallel computing 
terra::writeVector(
  bham_r,
  filename = "data/external/boundaries/bham_wards_27700.gpkg",
  filetype = "GPKG",
  overwrite = TRUE
)

spatvector_path_abs <- normalizePath(
  "data/external/boundaries/bham_wards_27700.gpkg",
  winslash = "/",
  mustWork = TRUE
)
#-----------------------------------------------------------
#============================================================
#data from https://data.ceda.ac.uk/badc/deposited2021/chess-scape/data/rcp85_bias-corrected/01/daily/tas
# tas_file = sort(list.files(path = "data/raw/chess_scape/chess85baiscorrect/01/", full.names = TRUE, pattern = ".nc"))

spatvector_path_abs = normalizePath("data/external/boundaries/bham_wards_2022.gpkg", winslash = "/", mustWork = TRUE)

source("R/chess_scape_function_raster2ts.R")
#====================================================================================================
#----------------------------------------------------------------------------------------------------
#RCP 8.5 bias corrected member 1 
#----------------------------------------------------------------------------------------------------
start = Sys.time()

chessscape85_01 = raster_to_ts_chessscape(
  tas_file_path = sort(list.files(path = "data/raw/chess_scape/chess85baiscorrect/01/", full.names = TRUE, pattern = ".nc")),
  member_id = "1",
  spatvector = bham_r,
  spatvector_path = spatvector_path_abs,
  id_match_file = bham_id_ward,
  parallel = TRUE,
  workers = 10,
  show_progress = TRUE
)

end = Sys.time()

end-start

#----------------------------------------------------------------------------------------------------
#RCP 8.5 bias corrected member 4 
#----------------------------------------------------------------------------------------------------
start = Sys.time()

chessscape85_04 = raster_to_ts_chessscape(
  tas_file_path = sort(list.files(path = "data/raw/chess_scape/chess85baiscorrect/04/", full.names = TRUE, pattern = ".nc")),
  member_id = "4",
  spatvector = bham_r,
  spatvector_path = spatvector_path_abs,
  id_match_file = bham_id_ward,
  parallel = TRUE,
  workers = 10,
  show_progress = TRUE
)

end = Sys.time()

end-start

#----------------------------------------------------------------------------------------------------
#RCP 8.5 bias corrected member 6 
#----------------------------------------------------------------------------------------------------
start = Sys.time()

chessscape85_06 = raster_to_ts_chessscape(
  tas_file_path = sort(list.files(path = "data/raw/chess_scape/chess85baiscorrect/06/", full.names = TRUE, pattern = ".nc")),
  member_id = "6",
  spatvector = bham_r,
  spatvector_path = spatvector_path_abs,
  id_match_file = bham_id_ward,
  parallel = TRUE,
  workers = 12,
  show_progress = TRUE
)

end = Sys.time()

end-start

#----------------------------------------------------------------------------------------------------
#RCP 8.5 bias corrected member 15 
#----------------------------------------------------------------------------------------------------
start = Sys.time()

chessscape85_15 = raster_to_ts_chessscape(
  tas_file_path = sort(list.files(path = "data/raw/chess_scape/chess85baiscorrect/15/", full.names = TRUE, pattern = ".nc")),
  member_id = "15",
  spatvector = bham_r,
  spatvector_path = spatvector_path_abs,
  id_match_file = bham_id_ward,
  parallel = TRUE,
  workers = 12,
  show_progress = TRUE
)

end = Sys.time()

end-start

#----------------------------------------------------------------------------------------------------
#bind the chess scape rcp85 data
chessscape85_daily_ward = rbind(chessscape_01,
      chessscape_04,
      chessscape_06,
      chessscape_15)

write_rds(chessscape85_daily_ward, "data/processed/chessscape85_daily_ward.rds")
#====================================================================================================











