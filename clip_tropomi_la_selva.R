library(terra)
library(ncdf4)
library(parallel)

### roi_file can be a path to a shapefile or a manually created polygon using vect()

terraOptions(memfrac = 0.8) # Fraction of memory to allow terra

tmpdir          <- "/mnt/c/Rwork"
out_dir         <- "/mnt/g/TROPOMI/esa/extracted/ebf/La_Selva"
out_name        <- "/La_Selva_TROPOSIF_L2B_"
f_list          <- list.files("/mnt/g/TROPOMI/esa/original/v2.1/l2b", pattern = "*.nc", full.names = TRUE, recursive = TRUE)
land_cover      <- 2    # Set to NULL if not filtering land cover class
land_cover_var  <- "PRODUCT/LC_MASK_2020" # Can be default or one we added
land_cover_perc <- "PRODUCT/LC_PERC_2020"
cloud_fraction  <- NULL # Set to NULL if not filtering cloud fraction
notes           <- "This data has been filtered to include only soundings in area of La Selva and the Braulio Carrillo National Park"

### Polygons for clipping
roi_file       <- "/mnt/g/Costa_Rica/La_Selva/La_Selva_Area.shp"

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

clip_TROPOSIF <- function(input_file, roi_file, out_dir, out_name, land_cover,
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
  coords           <- cbind(ncvar_get(t_data, "PRODUCT/longitude"), ncvar_get(t_data, "PRODUCT/latitude"))
  colnames(coords) <- c("lon", "lat")
  t                <- basename(input_file)
  t                <- substr(t, 14, 23)
  
  # Get variables and transform to vect for clipping to ROI
  df_var                  <- data.frame(lon = ncvar_get(t_data, "PRODUCT/longitude"))
  df_var$lat              <- ncvar_get(t_data, "PRODUCT/latitude")
  df_var$SIF_743          <- ncvar_get(t_data, "PRODUCT/SIF_743")
  df_var$SIF_Corr_743     <- ncvar_get(t_data, "PRODUCT/SIF_Corr_743")
  df_var$SIF_ERROR_743    <- ncvar_get(t_data, "PRODUCT/SIF_ERROR_743")
  df_var$phase_angle      <- ncvar_get(t_data, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/phase_angle")
  df_var$cloud_fraction   <- ncvar_get(t_data, "PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_fraction_L2")
  df_var$LC_MASK          <- ncvar_get(t_data, land_cover_var)
  if (!is.null(land_cover_perc)) {
    df_var$LC_PERC <- ncvar_get(t_data, land_cover_perc)
  }
  nc_close(t_data)
  
  if (!is.null(land_cover)) {
    df_var <- df_var[df_var$LC_MASK == land_cover, ]
  }
  
  if (!is.null(cloud_fraction)) {
    df_var <- df_var[df_var$cloud_fraction < cloud_fraction, ]
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
    elemdim <- ncdim_def("n_elem", "", seq(1, nrow(df)))
    
    t_num   <- as.numeric(julian(as.Date(t), origin = as.Date("1970-01-01")))
    # timedim <- ncdim_def("time", "days since 1970-01-01", t_num)
    # londim  <- ncdim_def("lon", "degrees_east", as.double(df$x)) 
    # latdim  <- ncdim_def("lat", "degrees_north", as.double(df$y))
    
    # define variables
    fillvalue     <- -9999
    dlname        <- "time"
    time_def      <- ncvar_def("time", "days since 1970-01-01", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "longitude"
    lon_def       <- ncvar_def("lon", "degrees_east", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "latitude"
    lat_def       <- ncvar_def("lat", "degrees_north", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "SIF@740 (743-758nm)"
    sif_def       <- ncvar_def("SIF_743", "mW/m2/sr/nm", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "daylength-corr SIF@740 (743-758nm)"
    sif_d_def     <- ncvar_def("SIF_Corr_743", "mW/m2/sr/nm", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "1-sigma SIF retrieval error (743-758nm)"
    sif_e_def     <- ncvar_def("SIF_ERROR_743", "mW/m2/sr/nm", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "phase angle"
    pa_def        <- ncvar_def("phase_angle", "degrees", elemdim, fillvalue, dlname, prec = "float")
    
    dlname        <- "cloud fraction"
    cf_def        <- ncvar_def("cloud_fraction_L2", "fraction", elemdim, fillvalue, dlname, prec = "float")
    
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
    message(out_f)
    
    if (!is.null(land_cover_perc)) {
      ncout <- nc_create(out_f,
                         list(time_def, lon_def, lat_def, sif_def,
                              sif_d_def, sif_e_def, pa_def, cf_def, lc_def, lc_perc_def), 
                         force_v4 = TRUE)
    } else {
      ncout <- nc_create(out_f,
                         list(time_def, lon_def, lat_def, sif_def,
                              sif_d_def, sif_e_def, pa_def, cf_def, lc_def), 
                         force_v4 = TRUE)
    }
    
    # put variables
    ncvar_put(ncout, time_def, rep(t_num, times = nrow(df)))
    ncvar_put(ncout, lon_def, df$x)
    ncvar_put(ncout, lat_def, df$y)
    ncvar_put(ncout, sif_def, df$SIF_743)
    ncvar_put(ncout, sif_d_def, df$SIF_Corr_743)
    ncvar_put(ncout, sif_e_def, df$SIF_ERROR_743)
    ncvar_put(ncout, pa_def, df$phase_angle)
    ncvar_put(ncout, cf_def, df$cloud_fraction)
    ncvar_put(ncout, lc_def, df$LC_MASK)
    if (!is.null(land_cover_perc)) {
      ncvar_put(ncout, lc_perc_def, df$LC_PERC)
    }
    
    # put additional attributes into dimension and data variables
    ncatt_put(ncout,"lon","axis","X")
    ncatt_put(ncout,"lat","axis","Y")
    ncatt_put(ncout,"time","axis","T")
    
    # add global attributes
    ncatt_put(ncout,0,"title", "TROPOSIF_L2B")
    ncatt_put(ncout,0,"institution", "University of Oklahoma")
    ncatt_put(ncout,0,"source", "Russell Doughty, PhD")
    ncatt_put(ncout,0,"references", "http://ftp.sron.nl/open-access-data-2/TROPOMI/tropomi/sif/v2.1/l2b/)")
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

mclapply(f_list, clip_TROPOSIF, mc.cores = 10, mc.preschedule = FALSE, roi_file = roi_file,
         out_dir = out_dir, out_name = out_name, land_cover = land_cover, land_cover_var = land_cover_var,
         land_cover_perc = land_cover_perc, cloud_fraction = cloud_fraction,  tmpdir = tmpdir)