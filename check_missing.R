
origin_list <- list.files("G:/TROPOMI/esa/original/v2.1/l2b/2021", pattern = "*.nc", full.names = TRUE, recursive = TRUE)
out_list    <- list.files("G:/TROPOMI/esa/extracted/Ghana/protected_reserves/2021", pattern = "*.nc", full.names = TRUE, recursive = TRUE)

for (i in 1:length(origin_list)) {
  t   <- basename(origin_list[i])
  t   <- substr(t, 14, 23)
  if (i == 1) {
    origin_dates <- t
  } else {
    origin_dates <- c(origin_dates, t)
  }
}

for (i in 1:length(out_list)) {
  t   <- basename(out_list[i])
  t   <- substr(t, 29, 38)
  if (i == 1) {
    out_dates <- t
  } else {
    out_dates <- c(out_dates, t)
  }
}

# compare output files with input and get list of missing dates
missing_dates  <- origin_dates[!(origin_dates %in% out_dates)]

# Now do partial matching to get complete file names for input files
missing_origin <- origin_list[sapply(missing_dates, function(x) { grep(x, origin_list) })]

# print so we can copy/paste
dput(missing_origin)
