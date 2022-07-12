library(terra)
library(ncdf4)
library(parallel)

terraOptions(memfrac = 0.8) # Fraction of memory to allow terra

tmpdir    <- "/mnt/c/Rwork"
roi_file  <- "/mnt/g/Africa/Ghana/Ghana_Disturbance_Data/Ghana_Protected_Reserves.shp"
out_dir   <- "/mnt/g/TROPOMI/esa/extracted/Ghana/protected_reserves/2021"
# f_list    <- list.files("/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021", pattern = "*.nc", full.names = TRUE, recursive = TRUE)
f_list <- c("/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/01/TROPOSIF_L2B_2021-01-01.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/01/TROPOSIF_L2B_2021-01-06.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/01/TROPOSIF_L2B_2021-01-17.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/01/TROPOSIF_L2B_2021-01-19.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/01/TROPOSIF_L2B_2021-01-25.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/01/TROPOSIF_L2B_2021-01-27.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/02/TROPOSIF_L2B_2021-02-02.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/02/TROPOSIF_L2B_2021-02-07.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/02/TROPOSIF_L2B_2021-02-18.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/02/TROPOSIF_L2B_2021-02-19.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/02/TROPOSIF_L2B_2021-02-27.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/03/TROPOSIF_L2B_2021-03-06.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/03/TROPOSIF_L2B_2021-03-14.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/03/TROPOSIF_L2B_2021-03-16.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-07.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-09.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-10.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-11.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-12.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-13.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-14.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-15.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-16.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-17.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-19.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-20.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-21.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-22.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-23.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-28.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-29.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/04/TROPOSIF_L2B_2021-04-30.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/05/TROPOSIF_L2B_2021-05-23.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/07/TROPOSIF_L2B_2021-07-30.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/08/TROPOSIF_L2B_2021-08-22.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/10/TROPOSIF_L2B_2021-10-09.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/10/TROPOSIF_L2B_2021-10-21.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/10/TROPOSIF_L2B_2021-10-27.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/11/TROPOSIF_L2B_2021-11-01.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/11/TROPOSIF_L2B_2021-11-06.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/11/TROPOSIF_L2B_2021-11-19.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/11/TROPOSIF_L2B_2021-11-24.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/12/TROPOSIF_L2B_2021-12-06.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/12/TROPOSIF_L2B_2021-12-08.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/12/TROPOSIF_L2B_2021-12-12.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/12/TROPOSIF_L2B_2021-12-19.nc", 
            "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2021/12/TROPOSIF_L2B_2021-12-24.nc"
)

# in_dir    <- "/mnt/g/TROPOMI/esa/original/v2.1/l2b/2019"
# roi_file  <- "/mnt/g/Africa/Ghana/Ghana_Disturbance_Data/Ghana_Protected_Reserves.shp"
# out_dir   <- "/mnt/g/TROPOMI/esa/extracted/Ghana/protected_reserves"

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

clip_TROPOSIF <- function(input_file, roi_file, out_dir, tmpdir) {
  
  tmp_create(tmpdir)
  
  if (!dir.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }
  
  roi <- vect(roi_file) # vectorize roi shp file for clipping
  
  time_s <- Sys.time()

  t_data <- nc_open(input_file)
  
  # Get spatial and time
  lon <- ncvar_get(t_data, "PRODUCT/longitude")
  lat <- ncvar_get(t_data, "PRODUCT/latitude")
  t   <- basename(input_file)
  t   <- substr(t, 14, 23)

  # Get variables and transform to vect for clipping to ROI
  var      <- data.frame(NDVI = ncvar_get(t_data, "PRODUCT/NDVI"))
  vec      <- vect(cbind(lon, lat), atts = var, crs = "+proj=longlat +datum=WGS84")
  var_roi  <- intersect(vec, roi)
  
  # If number of soundings > 0, then proceed
  if (nrow(crds(var_roi, df = TRUE)) == 0) {
    print(paste0("File for this date is being skipped as it has 0 soundings for the region: ", t))
    # Close input file
    nc_close(t_data)
    tmp_remove(tmpdir)
    
  } else {
    # Build data frame for writing to nc file
    df <- crds(var_roi, df = TRUE)
    df <- cbind(df, var_roi[[1]])
    
    # Repeat for each variable
    var      <- data.frame(NIRv = ncvar_get(t_data, "PRODUCT/NIRv"))
    vec      <- vect(cbind(lon, lat), atts = var, crs = "+proj=longlat +datum=WGS84")
    var_roi  <- intersect(vec, roi)
    df       <- cbind(df, var_roi[[1]])
    
    var      <- data.frame(NIRv_RAD = ncvar_get(t_data, "PRODUCT/NIRv_RAD"))
    vec      <- vect(cbind(lon, lat), atts = var, crs = "+proj=longlat +datum=WGS84")
    var_roi  <- intersect(vec, roi)
    df       <- cbind(df, var_roi[[1]])
    
    var      <- data.frame(SIF_743 = ncvar_get(t_data, "PRODUCT/SIF_743"))
    vec      <- vect(cbind(lon, lat), atts = var, crs = "+proj=longlat +datum=WGS84")
    var_roi  <- intersect(vec, roi)
    df       <- cbind(df, var_roi[[1]])
    
    var      <- data.frame(SIF_Corr_743 = ncvar_get(t_data, "PRODUCT/SIF_Corr_743"))
    vec      <- vect(cbind(lon, lat), atts = var, crs = "+proj=longlat +datum=WGS84")
    var_roi  <- intersect(vec, roi)
    df       <- cbind(df, var_roi[[1]])
    
    var      <- data.frame(SIF_ERROR_743 = ncvar_get(t_data, "PRODUCT/SIF_ERROR_743"))
    vec      <- vect(cbind(lon, lat), atts = var, crs = "+proj=longlat +datum=WGS84")
    var_roi  <- intersect(vec, roi)
    df       <- cbind(df, var_roi[[1]])
    
    var      <- data.frame(Mean_TOA_RAD_743 = ncvar_get(t_data, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/Mean_TOA_RAD_743"))
    vec      <- vect(cbind(lon, lat), atts = var, crs = "+proj=longlat +datum=WGS84")
    var_roi  <- intersect(vec, roi)
    df       <- cbind(df, var_roi[[1]])
    
    var      <- data.frame(phase_angle = ncvar_get(t_data, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/phase_angle"))
    vec      <- vect(cbind(lon, lat), atts = var, crs = "+proj=longlat +datum=WGS84")
    var_roi  <- intersect(vec, roi)
    df       <- cbind(df, var_roi[[1]])
    
    var      <- data.frame(cloud_fraction_L2 = ncvar_get(t_data, "PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_fraction_L2"))
    vec      <- vect(cbind(lon, lat), atts = var, crs = "+proj=longlat +datum=WGS84")
    var_roi  <- intersect(vec, roi)
    df       <- cbind(df, var_roi[[1]])

    var      <- c() # kick out of memory
    vec      <- c()
    var_roi  <- c()
    
    #### Create NC file ####
    ### Note: When creating point files, use number of points as a dim
    ### rather than lon and lat, and make lon and lat variables.
    ###
    
    # Create dimensions nc file
    elemdim <- ncdim_def("n_elem", "", seq(1, length(df$NDVI)))
    
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
    
    # Close input file
    nc_close(ncout)
    nc_close(t_data)
    
    time_e <- Sys.time()
    time_dif <- difftime(time_e, time_s)
    print(paste0("Saved ", out_f, ". Time elapsed: ", time_dif))
  }
  
  tmp_remove(tmpdir)
}

mclapply(f_list, clip_TROPOSIF, mc.cores = 3, mc.preschedule = FALSE, roi_file = roi_file, out_dir = out_dir, tmpdir = tmpdir)