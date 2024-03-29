
library(terra)
library(ncdf4)
library(parallel)

### roi_file can be a path to a shapefile or a manually created polygon using vect()

terraOptions(memfrac = 0.8) # Fraction of memory to allow terra

tmpdir          <- "/mnt/c/Rwork"
out_dir         <- "/mnt/g/GOSAT_SIF/extracted/africa"
out_name        <- "/Africa_gosat_LtSIF_"
f_list          <- list.files("/mnt/g/GOSAT_SIF/original", pattern = "*.nc", full.names = TRUE, recursive = TRUE)
notes           <- ""

# Africa project
roi_file         <- "/mnt/g/Africa/Tropical_Africa_Ecoregions/Tropical_Africa_Ecoregions1.shp"
forest_mask      <- vect("/mnt/g/Africa/Forest_Masks/dissolved/Africa_merged_2019_2.5km_Buffer.shp")

tmp_create <- function(tmpdir) {
  
  p_tmp_dir <- paste0(tmpdir, "/", as.character(Sys.getpid())) # Process ID
  
  if (!dir.exists(p_tmp_dir)) {
    dir.create(p_tmp_dir, recursive = TRUE)
  }
  
  terraOptions(tempdir = p_tmp_dir)
}
tmp_remove <- function(tmpdir) {
  
  p_tmp_dir <- paste0(tmpdir, "/", as.character(Sys.getpid())) # Process ID
  unlink(p_tmp_dir, recursive = TRUE)
}

clip_nc <- function(input_file, roi_file, forest_mask, out_dir, out_name, tmpdir) {
  
  tmp_create(tmpdir)
  
  if (!dir.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Get ecoregions
  roi_all   <- vect(roi_file)
  eco_names <- unique(roi_all$ECO_NAME)
  
  time_s <- Sys.time()
  
  t_data <- nc_open(input_file)
  
  # Get spatial and time
  coords           <- cbind(ncvar_get(t_data, "Longitude"), ncvar_get(t_data, "Latitude"))
  colnames(coords) <- c("lon", "lat")
  t                <- basename(input_file)
  t                <- substr(t, 13, 18)
  t                <- paste0("20", t)
  t                <- gsub("(\\d{4})(\\d{2})(\\d{2})$","\\1-\\2-\\3",t) # add dashes
  
  # Get variables and transform to vect for clipping to ROI
  df_var                 <- data.frame(lon = ncvar_get(t_data, "Longitude"))
  df_var$lat             <- ncvar_get(t_data, "Latitude")
  
  df_var$Daily_SIF_740nm <- colMeans(ncvar_get(t_data, "Daily_SIF_740nm"))
  df_var$Daily_SIF_757nm <- colMeans(ncvar_get(t_data, "Daily_SIF_757nm"))
  df_var$Daily_SIF_771nm <- colMeans(ncvar_get(t_data, "Daily_SIF_771nm"))
  
  df_var$phase_angle      <- ncvar_get(t_data, "PA")
  df_var$SZA              <- ncvar_get(t_data, "SZA")
  df_var$SAz              <- ncvar_get(t_data, "SAz")
  df_var$VZA              <- ncvar_get(t_data, "VZA")
  df_var$VAz              <- ncvar_get(t_data, "VAz")
  
  df_var$mode             <- ncvar_get(t_data, "Metadata/MeasurementMode")
  df_var$qc               <- ncvar_get(t_data, "Quality_Flag")[1,]
  df_var$cloud            <- ncvar_get(t_data, "Cloud/cloud_flag_abp")
  
  nc_close(t_data)
  
  # Put coords in their own
  coords <- cbind(df_var$lon, df_var$lat)
  df_var <- subset(df_var, select = -c(lon,lat))
  
  # Convert to vector
  t_vec      <- vect(coords, atts = df_var, crs = "+proj=longlat +datum=WGS84")
  
  for (e in 1:length(eco_names)) {
    
    time_s <- Sys.time()
    
    # Create dirs
    eco_region <- gsub(" ", "_", eco_names[e])
    eco_dir    <- paste0(out_dir, "/", eco_region)
    if (!dir.exists(eco_dir)){
      dir.create(eco_dir, recursive = TRUE)
    }
    
    roi <- vect(roi_file, query = paste0("SELECT * FROM Tropical_Africa_Ecoregions1 WHERE ECO_NAME = ", "'", eco_names[e], "'"))
    roi <- aggregate(roi)
    roi <- crop(forest_mask, roi)
    
    var_roi  <- intersect(t_vec, roi)
    
    # If number of soundings > 0, then proceed
    if (nrow(crds(var_roi, df = TRUE)) == 0) {
      message(paste0("File for this date is being skipped as it has 0 soundings for the region: ", t))
      
    } else {
      # Build data frame for writing to nc file
      df <- crds(var_roi, df = TRUE)
      
      for (i in 1:length(names(var_roi))) {
        df <- cbind(df, var_roi[[i]])
      }
      
      #### Create NC file ####
      ### Note: When creating point files, use number of points as a dim
      ### rather than lon and lat, and make lon and lat variables.
      ###
      
      # Create dimensions nc file
      elemdim <- ncdim_def("obs", "", seq(1, nrow(df)))
      
      t_num   <- as.numeric(julian(as.Date(t), origin = as.Date("1970-01-01")))
      # timedim <- ncdim_def("time", "days since 1970-01-01", t_num)
      # londim  <- ncdim_def("lon", "degrees_east", as.double(df$x)) 
      # latdim  <- ncdim_def("lat", "degrees_north", as.double(df$y))
      
      # define variables
      fillvalue     <- -9999
      dlname        <- "time"
      time_def      <- ncvar_def("time", "days since 1970-01-01", elemdim, fillvalue, dlname, prec = "float")
      
      dlname        <- "Longitude"
      lon_def       <- ncvar_def("lon", "degrees_east", elemdim, fillvalue, dlname, prec = "float")
      
      dlname        <- "Latitude"
      lat_def       <- ncvar_def("lat", "degrees_north", elemdim, fillvalue, dlname, prec = "float")
      
      dlname        <- "Daily_SIF_740nm"
      sif740_def    <- ncvar_def("Daily_SIF_740nm", "mW/m2/sr/nm", elemdim, fillvalue, dlname, prec = "float")
      
      dlname        <- "Daily_SIF_757nm"
      sif757_def    <- ncvar_def("Daily_SIF_757nm", "mW/m2/sr/nm", elemdim, fillvalue, dlname, prec = "float")
      
      dlname        <- "Daily_SIF_771nm"
      sif771_def    <- ncvar_def("Daily_SIF_771nm", "mW/m2/sr/nm", elemdim, fillvalue, dlname, prec = "float")
      
      dlname        <- "phase angle"
      pa_def        <- ncvar_def("PA", "degrees", elemdim, fillvalue, dlname, prec = "float")
      
      dlname        <- "solar zenith angle"
      sza_def       <- ncvar_def("SZA", "degrees", elemdim, fillvalue, dlname, prec = "float")
      
      dlname        <- "solar azimuth angle"
      saz_def       <- ncvar_def("SAz", "degrees", elemdim, fillvalue, dlname, prec = "float")
      
      dlname        <- "viewing zenith angle"
      vza_def       <- ncvar_def("VZA", "degrees", elemdim, fillvalue, dlname, prec = "float")
      
      dlname        <- "viewing azimuth angle"
      vaz_def       <- ncvar_def("VAz", "degrees", elemdim, fillvalue, dlname, prec = "float")
      
      dlname        <- "Instrument Measurement Mode, 0=Nadir, 1=Glint, 2=Target, 3=AreaMap, 4=Transition; users might consider to separate these for analysis"
      mm_def        <- ncvar_def("Mode", "", elemdim, fillvalue, dlname, prec = "float")
      
      dlname        <- "SIF Lite Quality Flag: 0 = best (passes quality control + cloud fraction = 0.0); 1 = good (passes quality control); 2 = bad (failed quality control); -1 = not investigated"
      qc_def        <- ncvar_def("qc", "", elemdim, fillvalue, dlname, prec = "float")
      
      dlname        <- "Indicator of whether the sounding contained clouds: 0 - Classified clear, 1 - Classified cloudy, 2 - Not classified, all other values undefined; not used in SIF processing"
      cloud_def     <- ncvar_def("cloud_flag_abp", "flag", elemdim, fillvalue, dlname, prec = "float")
      
      # create netCDF file and put arrays
      out_f <- paste0(eco_dir, out_name, eco_region, "_", t, ".nc")
      
      ncout <- nc_create(out_f,
                         list(time_def, lon_def, lat_def, sif740_def, sif757_def, sif771_def,
                              pa_def, sza_def, saz_def, vza_def, vaz_def,
                              mm_def, qc_def, cloud_def), 
                         force_v4 = TRUE)
      
      # put variables
      ncvar_put(ncout, time_def, rep(t_num, times = nrow(df)))
      ncvar_put(ncout, lon_def, df$x)
      ncvar_put(ncout, lat_def, df$y)
      ncvar_put(ncout, sif740_def, df$Daily_SIF_740nm)
      ncvar_put(ncout, sif757_def, df$Daily_SIF_757nm)
      ncvar_put(ncout, sif771_def, df$Daily_SIF_771nm)
      ncvar_put(ncout, pa_def, df$phase_angle)
      ncvar_put(ncout, sza_def, df$SZA)
      ncvar_put(ncout, saz_def, df$SAz)
      ncvar_put(ncout, vza_def, df$VZA)
      ncvar_put(ncout, vaz_def, df$VAz)
      ncvar_put(ncout, mm_def, df$mode)
      ncvar_put(ncout, qc_def, df$qc)
      ncvar_put(ncout, cloud_def, df$cloud)
      
      # put additional attributes into dimension and data variables
      ncatt_put(ncout,"lon","axis","X")
      ncatt_put(ncout,"lat","axis","Y")
      ncatt_put(ncout,"time","axis","T")
      
      # add global attributes
      ncatt_put(ncout,0,"title", "GOSAT SIF")
      ncatt_put(ncout,0,"institution", "University of Oklahoma")
      ncatt_put(ncout,0,"source", "Russell Doughty, PhD")
      ncatt_put(ncout,0,"date_created", date())
      ncatt_put(ncout,0,"notes", notes)
      
      # Close input file
      nc_close(ncout)
      
      time_e   <- Sys.time()
      time_dif <- difftime(time_e, time_s)
      
      message(paste0("Saved ", out_f, ". Time elapsed: ", time_dif))
    }
  }
  tmp_remove(tmpdir)
}

mclapply(f_list, clip_nc, mc.cores = 10, mc.preschedule = FALSE, roi_file = roi_file,
         forest_mask = forest_mask, out_dir = out_dir, out_name = out_name, tmpdir = tmpdir)