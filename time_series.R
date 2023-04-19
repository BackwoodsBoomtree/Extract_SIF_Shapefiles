
library(ncdf4)

file_df <- function(input_dir, year, time) {
  file_list <- list.files(input_dir, pattern = "*.nc", full.names = TRUE, recursive = TRUE)
  
  if (time == "8-day") {
    dates <- seq(as.Date(paste0(year,"-01-01")), as.Date(paste0((year),"-12-31")), by="days")
    
    # Create data frame with column for each 8-day file list
    for (i in 1:46) {
      
      sub_dates <- dates[(i * 8 - 7):(i * 8)]
      
      sub_files <- c()
      
      for (j in 1:length(sub_dates)) {
        
        check_file <- file_list[grepl(sub_dates[j], file_list)]
        
        if (length(check_file) != 0) {
          sub_files <- c(sub_files, check_file)
        } else {
          sub_files <- c(sub_files, NA)
        }
      }
      if (i == 1) {
        df <- cbind(sub_files)
      } else {
        df <- cbind(df, sub_files)
      }
    }
  }
  
  if (time == "16-day") {
    dates <- seq(as.Date(paste0(year,"-01-01")), as.Date(paste0((year + 1),"-12-31")), by="days")
    
    # Create data frame with column for each 16-day file list
    for (i in 1:23) {
      
      sub_dates <- dates[(i * 16 - 15):(i * 16)]
      
      sub_files <- c()
      
      for (j in 1:length(sub_dates)) {
        
        check_file <- file_list[grepl(sub_dates[j], file_list)]
        
        if (length(check_file) != 0) {
          sub_files <- c(sub_files, check_file)
        } else {
          sub_files <- c(sub_files, NA)
        }
      }
      if (i == 1) {
        df <- cbind(sub_files)
      } else {
        df <- cbind(df, sub_files)
      }
    }
  }
  
  if (time == "month") {
    dates <- seq(as.Date(paste0(year,"-01-01")), as.Date(paste0(year,"-12-31")), by="days")
    
    df <- data.frame(matrix(ncol = 12, nrow = 31))
    # Create data frame with column for each month
    for (i in 1:12){
      if (i < 10) {
        m <- paste0("0", i)
      } else {
        m <- as.character(i)
      }
      
      sub_dates <- subset(dates, format.Date(dates, "%m") == m)
      sub_files <- c()
      
      for (j in 1:length(sub_dates)) {
        
        check_file <- file_list[grepl(sub_dates[j], file_list)]
        
        if (length(check_file) != 0) {
          sub_files <- c(sub_files, check_file)
        } else {
          sub_files <- c(sub_files, NA)
        }
      }
      
      # Force length to 31
      if (length(sub_files) < 31) {
        sub_files <- sub_files[1:31]
      }
      
      if (i == 1) {
        df <- cbind(sub_files)
      } else {
        df <- cbind(df, sub_files)
      }
    }
  }
  
  return(df)
}

files_2019_8day <- file_df("G:/TROPOMI/esa/extracted/Ghana/protected_reserves/2019", 2019, "8-day")
files_2020_8day <- file_df("G:/TROPOMI/esa/extracted/Ghana/protected_reserves/2020", 2020, "8-day")
files_2021_8day <- file_df("G:/TROPOMI/esa/extracted/Ghana/protected_reserves/2021", 2021, "8-day")

files_2019_16day <- file_df("G:/TROPOMI/esa/extracted/Ghana/protected_reserves/2019", 2019, "16-day")
files_2020_16day <- file_df("G:/TROPOMI/esa/extracted/Ghana/protected_reserves/2020", 2020, "16-day")
files_2021_16day <- file_df("G:/TROPOMI/esa/extracted/Ghana/protected_reserves/2021", 2021, "16-day")

get_ts <- function(df, variable, time) {
  
  annual_df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(annual_df) <- c("Mean", "SD", "SEM")
  
  if (time == "8-day") {
    for (i in 1:46) {
      
      df_t <- (df[,i])
      df_t <- df_t[!is.na(df_t)]
      
      if (length(df_t) != 0) {
        for (j in 1:length(df_t)) {
          nc <- nc_open(df_t[j])
          
          if (j == 1){
            data <- ncvar_get(nc, variable)
          } else {
            data <- c(data, ncvar_get(nc, variable))
          }
        }
        
        nc_close(nc)
        
        annual_df[nrow(annual_df) + 1,] <- c(mean(data, rm.na = TRUE), sd(data), sd(data) / (sqrt(length(data))))
      
      } else {
      annual_df[nrow(annual_df) + 1,] <- c(NA, NA, NA)
      }
    }
  }
  
  if (time == "16-day") {
    for (i in 1:23) {
      
      df_t <- (df[,i])
      df_t <- df_t[!is.na(df_t)]
      
      if (length(df_t) != 0) {
        for (j in 1:length(df_t)) {
          nc <- nc_open(df_t[j])
          
          if (j == 1){
            data <- ncvar_get(nc, variable)
          } else {
            data <- c(data, ncvar_get(nc, variable))
          }
        }
        
        nc_close(nc)
        
        annual_df[nrow(annual_df) + 1,] <- c(mean(data, rm.na = TRUE), sd(data), sd(data) / (sqrt(length(data))))
        
      } else {
        annual_df[nrow(annual_df) + 1,] <- c(NA, NA, NA)
      }
    }
  }
  
  return(annual_df)
}

ts_2019_8day <- get_ts(files_2019_8day, "SIF_Corr_743", "8-day")
ts_2020_8day <- get_ts(files_2020_8day, "SIF_Corr_743", "8-day")
ts_2021_8day <- get_ts(files_2021_8day, "SIF_Corr_743", "8-day")
ts_sif_8day  <- c(ts_2019_8day$Mean, ts_2020_8day$Mean, ts_2021_8day$Mean)
ts_sem_8day  <- c(ts_2019_8day$SEM, ts_2020_8day$SEM, ts_2021_8day$SEM)

ts_2019_16day <- get_ts(files_2019_16day, "SIF_Corr_743", "16-day")
ts_2020_16day <- get_ts(files_2020_16day, "SIF_Corr_743", "16-day")
ts_2021_16day <- get_ts(files_2021_16day, "SIF_Corr_743", "16-day")
ts_sif_16day  <- c(ts_2019_16day$Mean, ts_2020_16day$Mean, ts_2021_16day$Mean)
ts_sem_16day  <- c(ts_2019_16day$SEM, ts_2020_16day$SEM, ts_2021_16day$SEM)

#### 8-day plot ####
cairo_pdf("G:/TROPOMI/esa/extracted/Ghana/figs/protected_reserves_sif_e_8day.pdf", width = 7.5, height = 4.25)

op <- par(mar = c(3,3,2,1))
x      <- seq(1, 138)
x_year <- c(1, 47, 93, 138)
x_half <- c(24, 70, 116)
plot(x, ts_sif_8day, type = "l", axes = FALSE, lwd = 1.5, xlab = NA, ylab = NA, col = "red")
arrows(x0 = x, y0 = ts_sif_8day - ts_sem_8day, x1 = x, y1 = ts_sif_8day + ts_sem_8day, code = 3, angle = 90, length = 0.1)
axis(1, tck = 0.04, labels = TRUE, at = x_year, mgp=c(3, 0.2, 0))
axis(1, tck = 0.02, labels = FALSE, at = x_half)
axis(2, tck = 0.02, mgp=c(3, 0.2, 0), las = 2)

mtext(1, text = "2019 - 2021", line = 2)
mtext(2, text = "Daily Corrected SIF (mW/m2/sr/nm)", line = 2)
mtext(3, text = "8-day Mean SIFd Protected Reserves Ghana", line = 0.5)

box()

dev.off()


#### 16-day plot ####
cairo_pdf("G:/TROPOMI/esa/extracted/Ghana/figs/protected_reserves_sif_e_16day.pdf", width = 7.5, height = 4.25)

op <- par(mar = c(3,3,2,1))
x      <- seq(1, 69)
x_year <- c(1, 24, 47, 69)
x_half <- c(12, 35, 58)
plot(x, ts_sif_16day, type = "l", axes = FALSE, lwd = 1.5, xlab = NA, ylab = NA, col = "red")
arrows(x0 = x, y0 = ts_sif_16day - ts_sem_16day, x1 = x, y1 = ts_sif_16day + ts_sem_16day, code = 3, angle = 90, length = 0.1)
axis(1, tck = 0.04, labels = TRUE, at = x_year, mgp=c(3, 0.2, 0))
axis(1, tck = 0.02, labels = FALSE, at = x_half)
axis(2, tck = 0.02, mgp=c(3, 0.2, 0), las = 2)

mtext(1, text = "2019 - 2021", line = 2)
mtext(2, text = "Daily Corrected SIF (mW/m2/sr/nm)", line = 2)
mtext(3, text = "16-day Mean SIFd Protected Reserves Ghana", line = 0.5)

box()

dev.off()