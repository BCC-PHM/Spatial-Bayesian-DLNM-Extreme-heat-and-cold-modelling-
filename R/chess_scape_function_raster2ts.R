
raster_to_ts_chessscape =function(
    tas_file_path = sort(list.files(
      path = "data/raw/chess_scape/chess85baiscorrect/01/",
      full.names = TRUE,
      pattern = "\\.nc$"
    )),
    member_id,
    spatvector        = NULL,          # used in non parallel for loop
    spatvector_path   = NULL,          # REQUIRED for parallel (recommended)
    id_match_file     = bham_id_ward,
    parallel          = FALSE,
    workers           = max(1, future::availableCores() - 6),
    show_progress     = TRUE
){
  
  if (missing(member_id) || is.null(member_id) ||
      length(member_id) != 1 || is.na(member_id) || trimws(member_id) == "") {
    warning("member_id was not provided (or is blank/NA). Supply e.g. member_id = '01'.")
    stop("Stopping: member_id is required to label this extraction.", call. = FALSE)
  }
  
  if (parallel == TRUE && (is.null(spatvector_path) || trimws(spatvector_path) == "")) {
    stop("Parallel mode requires spatvector_path (path to .shp/.gpkg) so each worker can load it.", call. = FALSE)
  }
  
  library(terra); library(dplyr); library(tidyr)
  
  # ----- helper: process the non parallel for loop -----
  process_one =function(tasfirst, member_id, id_match_file, spatvector, spatvector_path){
    
    tas_stack =terra::rast(tasfirst)
    
    # only set CRS if missing
    if (is.na(terra::crs(tas_stack))) terra::crs(tas_stack) ="EPSG:27700"
    
    # get polygons (robust for parallel)
    if (!is.null(spatvector_path)) {
      sv =terra::vect(spatvector_path)
    } else {
      sv =spatvector
      if (!inherits(sv, "SpatVector")) sv =terra::vect(sv)
    }
    
    # ensure CRS matches raster (use same.crs for robustness)
    if (!is.na(terra::crs(sv)) && !is.na(terra::crs(tas_stack)) &&
        !terra::same.crs(sv, tas_stack)) {
      sv =terra::project(sv, terra::crs(tas_stack))
    }
    
    tas_bham_crop =terra::crop(tas_stack, sv)
    ext_raw =terra::extract(tas_bham_crop, sv, fun = mean, na.rm = TRUE)
    
    # time vector (fallback if missing)
    raster_time =terra::time(tas_bham_crop)
    if (is.null(raster_time)) {
      raster_dates = as.Date(seq_len(nlyr(tas_bham_crop)) - 1, origin = "1970-01-01")
    } else {
      raster_dates = as.Date(raster_time)
    }
    
    #name the column with their corresponding date
    colnames(ext_raw)[2:ncol(ext_raw)] =as.character(raster_dates)
    
    out =ext_raw %>%
      tidyr::pivot_longer(cols = -ID, names_to = "Date", values_to = "tasmean") %>%
      dplyr::mutate(
        Date    = as.Date(Date),
        tasmean = tasmean - 273.15,
        Member  = as.character(member_id)
      ) %>%
      dplyr::left_join(id_match_file, by = "ID")
    
    return(out)
  }
  
  # ----- PARALLEL -----
  if (parallel == TRUE) {
    
    library(foreach)
    library(doSNOW)
    
    cl =parallel::makeCluster(workers)
    doSNOW::registerDoSNOW(cl)
    
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)
    
    pb =NULL
    opts =list()
    if (show_progress) {
      pb =txtProgressBar(min = 0, max = length(tas_file_path), style = 3)
      on.exit(try(close(pb), silent = TRUE), add = TRUE)
      opts =list(progress = function(n) setTxtProgressBar(pb, n))
    }
    
    output =foreach(
      i = seq_along(tas_file_path),
      .packages = c("terra", "dplyr", "tidyr"),
      .export   = c("process_one", "member_id", "id_match_file", "spatvector_path"),
      .options.snow = opts,
      .errorhandling = "stop"
    ) %dopar% {
      
      process_one(
        tasfirst = tas_file_path[i],
        member_id = member_id,
        id_match_file = id_match_file,
        spatvector = NULL,
        spatvector_path = spatvector_path
      )
    }
    
    return(dplyr::bind_rows(output))
  }
  
  # ----- SERIAL -----
  if (is.null(spatvector)) stop("Serial mode requires spatvector (your bham_r).", call. = FALSE)
  
  pb =NULL
  if (show_progress) {
    pb =txtProgressBar(min = 0, max = length(tas_file_path), style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  
  output =vector("list", length(tas_file_path))
  for (i in seq_along(tas_file_path)) {
    if (show_progress) setTxtProgressBar(pb, i)
    
    output[[i]] =process_one(
      tasfirst = tas_file_path[i],
      member_id = member_id,
      id_match_file = id_match_file,
      spatvector = spatvector,
      spatvector_path = NULL
    )
  }
  
  dplyr::bind_rows(output)
  
}


