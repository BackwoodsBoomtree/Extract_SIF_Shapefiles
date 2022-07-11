library(terra)
library(ncdf4)

in_dir    <- "G:/TROPOMI/esa/original/v2.1/l2b/2019"
roi_file  <- "G:/Africa/Ghana/Ghana_Disturbance_Data/Ghana_Protected_Reserves.shp"
out_dir   <- "G:/TROPOMI/esa/extracted/Ghana/protected_reserves"

clip_TROPOSIF <- function (in_dir, roi_file, out_dir){
  
  if (!dir.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }
  
  # List tropo files
  f_list <- list.files(in_dir, pattern = "*.nc", full.names = TRUE, recursive = TRUE)
  
  # vectorize roi shp file for clipping
  roi <- vect(roi_file)
  
  for (i in 1:length(f_list)) {
    
    time_s <- Sys.time()
    
    print()
    t_data <- nc_open(f_list[i])
    
    # Get spatial and time
    lon <- ncvar_get(t_data, "PRODUCT/longitude")
    lat <- ncvar_get(t_data, "PRODUCT/latitude")
    t   <- basename(f_list[i])
    t   <- substr(t, 14, 23)
    
    # Get variables
    ndvi   <- data.frame(NDVI = ncvar_get(t_data, "PRODUCT/NDVI"))
    nirv   <- data.frame(NIRv = ncvar_get(t_data, "PRODUCT/NIRv"))
    nirv_r <- data.frame(NIRv_RAD = ncvar_get(t_data, "PRODUCT/NIRv_RAD"))
    sif    <- data.frame(SIF_743 = ncvar_get(t_data, "PRODUCT/SIF_743"))
    sif_d  <- data.frame(SIF_Corr_743 = ncvar_get(t_data, "PRODUCT/SIF_Corr_743"))
    sif_e  <- data.frame(SIF_ERROR_743 = ncvar_get(t_data, "PRODUCT/SIF_ERROR_743"))
    rad    <- data.frame(Mean_TOA_RAD_743 = ncvar_get(t_data, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/Mean_TOA_RAD_743"))
    pa     <- data.frame(phase_angle = ncvar_get(t_data, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/phase_angle"))
    cf     <- data.frame(cloud_fraction_L2 = ncvar_get(t_data, "PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_fraction_L2"))
    
    # Transform to vect for clipping to ROI
    ndvi_v     <- vect(cbind(lon, lat), atts = ndvi, crs = "+proj=longlat +datum=WGS84")
    ndvi_roi   <- intersect(ndvi_v, roi)
    nirv_v     <- vect(cbind(lon, lat), atts = nirv, crs = "+proj=longlat +datum=WGS84")
    nirv_roi   <- intersect(nirv_v, roi)
    nirv_r_v   <- vect(cbind(lon, lat), atts = nirv_r, crs = "+proj=longlat +datum=WGS84")
    nirv_r_roi <- intersect(nirv_r_v, roi)
    sif_v      <- vect(cbind(lon, lat), atts = sif, crs = "+proj=longlat +datum=WGS84")
    sif_roi    <- intersect(sif_v, roi)
    sif_d_v    <- vect(cbind(lon, lat), atts = sif_d, crs = "+proj=longlat +datum=WGS84")
    sif_d_roi  <- intersect(sif_d_v, roi)
    sif_e_v    <- vect(cbind(lon, lat), atts = sif_e, crs = "+proj=longlat +datum=WGS84")
    sif_e_roi  <- intersect(sif_e_v, roi)
    rad_v      <- vect(cbind(lon, lat), atts = rad, crs = "+proj=longlat +datum=WGS84")
    rad_roi    <- intersect(rad_v, roi)
    rad_v      <- vect(cbind(lon, lat), atts = rad, crs = "+proj=longlat +datum=WGS84")
    rad_roi    <- intersect(rad_v, roi)
    pa_v       <- vect(cbind(lon, lat), atts = pa, crs = "+proj=longlat +datum=WGS84")
    pa_roi     <- intersect(pa_v, roi)
    cf_v       <- vect(cbind(lon, lat), atts = cf, crs = "+proj=longlat +datum=WGS84")
    cf_roi     <- intersect(cf_v, roi)
    
    # Construct df of all values
    df <- crds(ndvi_roi, df = TRUE)
    df <- cbind(df, ndvi_roi[[1]])
    df <- cbind(df, nirv_roi[[1]])
    df <- cbind(df, nirv_r_roi[[1]])
    df <- cbind(df, sif_roi[[1]])
    df <- cbind(df, sif_d_roi[[1]])
    df <- cbind(df, sif_e_roi[[1]])
    df <- cbind(df, rad_roi[[1]])
    df <- cbind(df, pa_roi[[1]])
    df <- cbind(df, cf_roi[[1]])
    
    #### Create NC file ####
    ### Note: When creating point files, use number of points as a dim
    ### rather than lon and lat, and make lon and lat variables.
    ###
    
    # Create dimensions nc file
    elemdim <- ncdim_def("n_elem", "", seq(1,length(df$NDVI)))
    
    t_num   <- as.numeric(julian(as.Date(t), origin = as.Date("1970-01-01")))
    # timedim <- ncdim_def("time", "days since 1970-01-01", t_num)
    # londim  <- ncdim_def("lon", "degrees_east", as.double(df$x)) 
    # latdim  <- ncdim_def("lat", "degrees_north", as.double(df$y))
    
    # define variables
    fillvalue  <- -9999
    dlname     <- "time"
    time_def   <- ncvar_def("time", "days since 1970-01-01", elemdim, fillvalue, dlname, prec = "float")
    dlname     <- "longitude"
    lon_def    <- ncvar_def("lon", "degrees_east", elemdim, fillvalue, dlname, prec = "float")
    dlname     <- "latitude"
    lat_def    <- ncvar_def("lat", "degrees_north", elemdim, fillvalue, dlname, prec = "float")
    dlname     <- "Normalized Difference Vegetation Index"
    ndvi_def   <- ncvar_def("NDVI", "Index", elemdim, fillvalue, dlname, prec = "float")
    dlname     <- "NIR Reflectance of Vegetation"
    nirv_def   <- ncvar_def("NIRv", "Index", elemdim, fillvalue, dlname, prec = "float")
    dlname     <- "NIRv Radiance"
    nirv_r_def <- ncvar_def("NIRv", "Index", elemdim, fillvalue, dlname, prec = "float")
    dlname     <- "retrieved SIF@740 (743-758nm)"
    sif_def    <- ncvar_def("SIF_743", "mW/m2/sr/nm", elemdim, fillvalue, dlname, prec = "float")
    dlname     <- "daylength-corr SIF@740 (743-758nm)"
    sif_d_def  <- ncvar_def("SIF_Corr_743", "mW/m2/sr/nm", elemdim, fillvalue, dlname, prec = "float")
    dlname     <- "1-sigma SIF retrieval error (743-758nm)"
    sif_e_def  <- ncvar_def("SIF_ERROR_743", "mW/m2/sr/nm", elemdim, fillvalue, dlname, prec = "float")
    dlname     <- "Mean TOA Radiance in 743-758 nm fitting window"
    rad_def    <- ncvar_def("Mean_TOA_RAD_743", "mW/m2/sr/nm", elemdim, fillvalue, dlname, prec = "float")
    dlname     <- "Phase Angle"
    pa_def     <- ncvar_def("phase_angle", "degrees", elemdim, fillvalue, dlname, prec = "float")
    dlname     <- "cloud fraction"
    cf_def     <- ncvar_def("cloud_fraction_L2", "fraction", elemdim, fillvalue, dlname, prec = "float")
    
    # create netCDF file and put arrays
    out_f <- paste0(out_dir, "/Ghana_Reserves_TROPOSIF_L2B_", t, ".nc")
    ncout <- nc_create(out_f,
                       list(time_def, lon_def, lat_def, ndvi_def, nirv_def, sif_def,
                            sif_d_def, sif_e_def, rad_def, pa_def, cf_def), 
                       force_v4 = TRUE)
    
    # put variables
    ncvar_put(ncout, time_def, rep(t_num, times = length(df$NDVI)))
    ncvar_put(ncout, lon_def, df$x)
    ncvar_put(ncout, lat_def, df$y)
    ncvar_put(ncout, ndvi_def, df$NDVI)
    ncvar_put(ncout, nirv_def, df$NIRv)
    ncvar_put(ncout, nirv_r_def, df$NIRv_RAD)
    ncvar_put(ncout, sif_def, df$SIF_743)
    ncvar_put(ncout, sif_d_def, df$SIF_Corr_743)
    ncvar_put(ncout, sif_e_def, df$SIF_ERROR_743)
    ncvar_put(ncout, rad_def, df$Mean_TOA_RAD_743)
    ncvar_put(ncout, pa_def, df$phase_angle)
    ncvar_put(ncout, cf_def, df$cloud_fraction_L2)
    
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
    
    # Get a summary of the created file:
    # ncout
    
    # Close input file
    nc_close(ncout)
    nc_close(t_data)
    
    time_e <- Sys.time()
    time_dif <- difftime(time_e, time_s)
    print(paste0("Saved ", out_f, ". Time elapsed: ", time_dif))
    print(paste0("Estimated time remaining: ", format(time_dif * (length(f_list) - i))), units = "hours")
    
  }
}

clip_TROPOSIF(in_dir, roi_file, out_dir)