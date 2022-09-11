library(terra)
library(ncdf4)
library(parallel)

### roi_file can be a path to a shapefile or a manually created polygon using vect()

terraOptions(memfrac = 0.8) # Fraction of memory to allow terra

tmpdir          <- "/mnt/c/Rwork"
out_dir         <- "/mnt/g/GOME2/extracted/amazon"
out_name        <- "/Amazon_NSIFv2.6.2.GOME-2A_"
f_list          <- list.files("/mnt/g/GOME2", pattern = "*.nc", full.names = TRUE, recursive = TRUE)
land_cover      <- 2    # Set to NULL if not filtering land cover class
land_cover_var  <- "LC_MASK_2020" # Can be default or one we added
land_cover_perc <- "LC_PERC_2020"
notes           <- "This data has been filtered to include only soundings that fall within EBF"

### Polygons for clipping
# roi_file       <- vect("POLYGON ((-18 -11, 52 -11, 52 14, -18 14, -18 -11))", crs="+proj=longlat +datum=WGS84") # Africa
roi_file       <- "/mnt/g/Amazon_shp/Amazon_poly.shp" # Amazon
# roi_file       <- vect("POLYGON ((-180 -23.5, 180 -23.5, 180 23.5, -180 23.5, -180 -23.5))", crs="+proj=longlat +datum=WGS84") # Tropics
# roi_file       <- "/mnt/g/SIF_comps/figs/EC_Sites/K67/K67_ebf.shp" # K67
# roi_file       <- "/mnt/g/SIF_comps/figs/EC_Sites/K34/K34_ebf.shp" # RJA

# Asia: Many soundings appear to be in the sea, so we need to also clip by coastlines
# roi_seasia <- vect("POLYGON ((72 -23.5, 180 -23.5, 180 23.5, 72 23.5, 72 -23.5))", crs="+proj=longlat +datum=WGS84") # Asia & Pacific
# coastlines <- vect("/mnt/c/Russell/R_Scripts/TROPOMI_2/mapping/GSHHS_shp/c/GSHHS_c_L1.shp")
# roi_file   <- intersect(roi_seasia, coastlines)

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

clip_nc <- function(input_file, roi_file, out_dir, out_name, land_cover,
                          land_cover_var, land_cover_perc, cloud_fraction, tmpdir) {
  
  tmp_create(tmpdir)
  
  if (!dir.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }
  
  time_s <- Sys.time()
  
  # if needed, vectorize roi shp file for clipping
  if (typeof(roi_file) != "S4"){
    roi <- vect(roi_file)
    roi <- aggregate(roi)
  } else {
    roi <- aggregate(roi_file)
  }
  
  t_data <- nc_open(input_file)
  
  # Get spatial and time
  coords           <- cbind(ncvar_get(t_data, "Longitude"), ncvar_get(t_data, "Latitude"))
  colnames(coords) <- c("lon", "lat")
  t                <- basename(input_file)
  t                <- substr(t, 20, 27)
  t                <- gsub("(\\d{4})(\\d{2})(\\d{2})$","\\1-\\2-\\3",t) # add dashes
  
  # Get variables and transform to vect for clipping to ROI
  df_var                 <- data.frame(lon = ncvar_get(t_data, "Longitude"))
  df_var$lat             <- ncvar_get(t_data, "Latitude")
  
  df_var$SIF_740            <- ncvar_get(t_data, "SIF_740")
  df_var$Daily_Averaged_SIF <- ncvar_get(t_data, "Daily_Averaged_SIF")

  df_var$phase_angle      <- ncvar_get(t_data, "PA")
  df_var$SZA              <- ncvar_get(t_data, "SZA")
  df_var$SAz              <- ncvar_get(t_data, "SAz")
  df_var$VZA              <- ncvar_get(t_data, "VZA")
  df_var$VAz              <- ncvar_get(t_data, "VAz")
  df_var$cloud            <- ncvar_get(t_data, "Cloud_Fraction")
  df_var$qc               <- ncvar_get(t_data, "Quality_Flag")
  df_var$LC_MASK          <- ncvar_get(t_data, land_cover_var)
  
  if (!is.null(land_cover_perc)) {
    df_var$LC_PERC <- ncvar_get(t_data, land_cover_perc)
  }
  nc_close(t_data)
  
  if (!is.null(land_cover)) {
    df_var <- df_var[df_var$LC_MASK == land_cover, ]
  }
  
  # Put coords in their own
  coords <- cbind(df_var$lon, df_var$lat)
  df_var <- subset(df_var, select = -c(lon,lat))
  
  # Clip data
  vec      <- vect(coords, atts = df_var, crs = "+proj=longlat +datum=WGS84")
  var_roi  <- intersect(vec, roi)
  
  # If number of soundings > 0, then proceed
  if (nrow(crds(var_roi, df = TRUE)) == 0) {
    message(paste0("File for this date is being skipped as it has 0 soundings for the region: ", t))
    
  } else {
    # Build data frame for writing to nc file
    df <- crds(var_roi, df = TRUE)
    
    for (i in 1:length(names(var_roi))) {
      df <- cbind(df, var_roi[[i]])
    }
    
    # kick out
    rm(df_var, vec, var_roi)
    
    invisible(gc())
    
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
    
    dlname        <- "SIF_740"
    sif740_def    <- ncvar_def("SIF_740", "mW/m2/sr/nm", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "Daily_Averaged_SIF"
    sifd_def      <- ncvar_def("Daily_Averaged_SIF", "mW/m2/sr/nm", elemdim, fillvalue, dlname, prec = "float")
    
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
    
    dlname        <- "0 - bad; 1 - good_passed_all_QC_checks; 2 - good_and_passed_cloud_check"
    qc_def        <- ncvar_def("Quality_Flag", "flag", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "Effective_cloud_fraction_MLER_model"
    cloud_def     <- ncvar_def("Cloud_Fraction", "fraction", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- basename(land_cover_var)
    lc_def        <- ncvar_def(basename(land_cover_var), "Majority IGBP Land Cover Class",
                               elemdim, fillvalue, dlname, prec = "float")
    
    if (!is.null(land_cover_perc)) {
      dlname        <- basename(land_cover_perc)
      lc_perc_def   <- ncvar_def(basename(land_cover_perc), "% of Majority IGBP Land Cover Class",
                                 elemdim, fillvalue, dlname, prec = "float")
    }
    
    # create netCDF file and put arrays
    out_f <- paste0(out_dir, out_name, t, ".nc")
    
    if (!is.null(land_cover_perc)) {
      ncout <- nc_create(out_f,
                         list(time_def, lon_def, lat_def, sif740_def, sifd_def,
                              pa_def, sza_def, saz_def, vza_def, vaz_def,
                              qc_def, cloud_def, lc_def, lc_perc_def), 
                         force_v4 = TRUE)
    } else {
      ncout <- nc_create(out_f,
                         list(time_def, lon_def, lat_def, sif740_def, sifd_def,
                              pa_def, sza_def, saz_def, vza_def, vaz_def,
                              qc_def, cloud_def, lc_def), 
                         force_v4 = TRUE)
    }
    
    # put variables
    ncvar_put(ncout, time_def, rep(t_num, times = nrow(df)))
    ncvar_put(ncout, lon_def, df$x)
    ncvar_put(ncout, lat_def, df$y)
    ncvar_put(ncout, sif740_def, df$SIF_740)
    ncvar_put(ncout, sifd_def, df$Daily_Averaged_SIF)
    ncvar_put(ncout, pa_def, df$phase_angle)
    ncvar_put(ncout, sza_def, df$SZA)
    ncvar_put(ncout, saz_def, df$SAz)
    ncvar_put(ncout, vza_def, df$VZA)
    ncvar_put(ncout, vaz_def, df$VAz)
    ncvar_put(ncout, qc_def, df$qc)
    ncvar_put(ncout, cloud_def, df$cloud)
    ncvar_put(ncout, lc_def, df$LC_MASK)
    if (!is.null(land_cover_perc)) {
      ncvar_put(ncout, lc_perc_def, df$LC_PERC)
    }
    
    # put additional attributes into dimension and data variables
    ncatt_put(ncout,"lon","axis","X")
    ncatt_put(ncout,"lat","axis","Y")
    ncatt_put(ncout,"time","axis","T")
    
    # add global attributes
    ncatt_put(ncout,0,"title", "GOME-2A NSIF v2.6.2")
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
  
  tmp_remove(tmpdir)
}

mclapply(f_list, clip_nc, mc.cores = 10, mc.preschedule = FALSE, roi_file = roi_file,
         out_dir = out_dir, out_name = out_name, land_cover = land_cover, land_cover_var = land_cover_var,
         land_cover_perc = land_cover_perc, cloud_fraction = cloud_fraction,  tmpdir = tmpdir)