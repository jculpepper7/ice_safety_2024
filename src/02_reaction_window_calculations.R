
# 1. Load libraries -------------------------------------------------------

library(ncdf4)
library(tidyverse)
library(lubridate)
library(here)
library(janitor)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(fields)
library(reshape)
library(introdataviz)
library(FSA)

# 2a. Read data-----------------------------------------------------------------

#read in the lake ice time series file
seasonal_ice_depth <- nc_open(here('data/cesm2_le_netcdf_files/icethick_daily_targetyears_ensmean.nc'))

# 2b. Extract variables of interest---------------------------------------------

#Extract ice duration data
ice_thickness <- ncvar_get(seasonal_ice_depth, 'LAKEICETHICK')
ilat_id <- ncvar_get(seasonal_ice_depth, 'lat')
ilon_id <- ncvar_get(seasonal_ice_depth, 'lon')

# 2c. Reorient time variable----------------------------------------------------
tt_id <- ncvar_get(seasonal_ice_depth, 'time')

#Convert time from days since 1850-01-01 to a vector of dates
itime_id <- as.Date(tt_id, origin = "1850-01-01")

#Extract the year
iyear_id <- year(itime_id)

#Extract seasonal cycle for specific years
idx_year1 <-  iyear_id >= 1851 & iyear_id <= 1851
idx_year2 <-  iyear_id >= 1899 & iyear_id <= 1901
idx_year3 <-  iyear_id >= 2009 & iyear_id <= 2011
idx_year4 <-  iyear_id >= 2042 & iyear_id <= 2044
idx_year5 <-  iyear_id >= 2086 & iyear_id <= 2088

#Get the day of year for each seasonal cycle

#index the time variable by the year of interest and extract the day of year
#These variables will be the x-axis of the time series plots
iyear1 <- itime_id[idx_year1] 
iyear2 <- itime_id[idx_year2]
iyear3 <- itime_id[idx_year3]
iyear4 <- itime_id[idx_year4]
iyear5 <- itime_id[idx_year5]

#Must revise the day of year so that the days are sequential and do not repeat.
idoy1 <- iyear1 - iyear1[1] + 1  
idoy2 <- iyear2 - iyear2[1] + 1  
idoy3 <- iyear3 - iyear3[1] + 1  
idoy4 <- iyear4 - iyear4[1] + 1  
idoy5 <- iyear5 - iyear5[1] + 1  

# 2d. Extract lat and long of interest from global data-------------------------

#Index locations (gives T/F)
idx_lat_id2 = ilat_id >= 40 & ilat_id <= 90
idx_lon_id2 = ilon_id >= 0 & ilon_id <= 360

# Filter ice thickness time series to lat and lon of interest

#North America
seasonal_ice_thick = ice_thickness[idx_lon_id2,idx_lat_id2,]

# 2e. Create time series based on temperature threshold-------------------------

seasonal_ice_thick2 = seasonal_ice_thick[,,idx_year2]
seasonal_ice_thick3 = seasonal_ice_thick[,,idx_year3]
seasonal_ice_thick4 = seasonal_ice_thick[,,idx_year4]
seasonal_ice_thick5 = seasonal_ice_thick[,,idx_year5]

seasonal_ice_thick2 <- seasonal_ice_thick2[,,225:589]
seasonal_ice_thick3 <- seasonal_ice_thick3[,,225:589]
seasonal_ice_thick4 <- seasonal_ice_thick4[,,225:589]
seasonal_ice_thick5 <- seasonal_ice_thick5[,,225:589]


# 2f. Ice phenology functions -------------------------------------------------

ice_on_func <- function(y){
  
  y2 <- NA
  if(all(!is.na(y))){
    y2 <- -1 #Eliminate lakes where not all time series has a value
    if(any(y>0.02)){
      y2 <- first(which(y>0.02))
    }
    if(y2==1){
      y2 <- -1 #eliminate pixels where lake is 2cm on day 1
    }
  }
  return(y2)
}

ice_on_safe_func <- function(y){
  
  y2 <- NA
  if(all(!is.na(y))){
    y2 <- -1 #Eliminate lakes where not all time series has a value
    if(any(y>0.1)){
      y2 <- first(which(y>0.1))
    }
    if(y2==1){
      y2 <- -1 #eliminate pixels where lake is 10cm on day 1
    }
  }
  return(y2)
}

ice_on_safe_15cm_func <- function(y){
  
  y2 <- NA
  if(all(!is.na(y))){
    y2 <- -1 #Eliminate lakes where not all time series has a value
    if(any(y>0.15)){
      y2 <- first(which(y>0.15))
    }
    if(y2==1){
      y2 <- -1 #eliminate pixels where lake is 15cm on day 1
    }
  }
  return(y2)
}

ice_on_safe_20cm_func <- function(y){
  
  y2 <- NA
  if(all(!is.na(y))){
    y2 <- -1 #Eliminate lakes where not all time series has a value
    if(any(y>0.2)){
      y2 <- first(which(y>0.2))
    }
    if(y2==1){
      y2 <- -1 #eliminate pixels where lake is 20cm on day 1
    }
  }
  return(y2)
}

ice_off_func <- function(y){
  
  y2 <- NA
  if(all(!is.na(y))){
    y2 <- -1 #Eliminate lakes where not all time series has a value
    if(any(y>0.02)){
      y2 <- last(which(y>0.02))
    }
    if(y2==1){
      y2 <- -1 #eliminate pixels where lake is 2cm on day 1
    }
  }
  return(y2)
}

ice_off_safe_func <- function(y){
  
  y2 <- NA
  if(all(!is.na(y))){
    y2 <- -1 #Eliminate lakes where not all time series has a value
    if(any(y>0.1)){
      y2 <- last(which(y>0.1))
    }
    if(y2==1){
      y2 <- -1 #eliminate pixels where lake is 10cm on day 1
    }
  }
  return(y2)
}

ice_off_safe_15cm_func <- function(y){
  
  y2 <- NA
  if(all(!is.na(y))){
    y2 <- -1 #Eliminate lakes where not all time series has a value
    if(any(y>0.15)){
      y2 <- last(which(y>0.15))
    }
    if(y2==1){
      y2 <- -1 #eliminate pixels where lake is 15cm on day 1
    }
  }
  return(y2)
}

ice_off_safe_20cm_func <- function(y){
  
  y2 <- NA
  if(all(!is.na(y))){
    y2 <- -1 #Eliminate lakes where not all time series has a value
    if(any(y>0.2)){
      y2 <- last(which(y>0.2))
    }
    if(y2==1){
      y2 <- -1 #eliminate pixels where lake is is 20cm on day 1
    }
  }
  return(y2)
}

ice_max_func <- function(y){
  
  y2 <- NA
  if(all(!is.na(y))){
    y2 <- -1 #Eliminate lakes where not all time series has a value
    if(any(y>0.1)){
      y2 <- which.max(y)
    }
    if(y2==1){
      y2 <- -1 
    }
  }
  return(y2)
}


# 2g. Reaction window data ------------------------------------------------

ice_on_2cm_0c <- apply(seasonal_ice_thick2, c(1,2), ice_on_func)
ice_on_2cm_1c <- apply(seasonal_ice_thick3, c(1,2), ice_on_func)
ice_on_2cm_2c <- apply(seasonal_ice_thick4, c(1,2), ice_on_func)
ice_on_2cm_4c <- apply(seasonal_ice_thick5, c(1,2), ice_on_func)

ice_off_2cm_0c <- apply(seasonal_ice_thick2, c(1,2), ice_off_func)
ice_off_2cm_1c <- apply(seasonal_ice_thick3, c(1,2), ice_off_func)
ice_off_2cm_2c <- apply(seasonal_ice_thick4, c(1,2), ice_off_func)
ice_off_2cm_4c <- apply(seasonal_ice_thick5, c(1,2), ice_off_func)

ice_on_10cm_0c <- apply(seasonal_ice_thick2, c(1,2), ice_on_safe_func)
ice_on_10cm_1c <- apply(seasonal_ice_thick3, c(1,2), ice_on_safe_func)
ice_on_10cm_2c <- apply(seasonal_ice_thick4, c(1,2), ice_on_safe_func)
ice_on_10cm_4c <- apply(seasonal_ice_thick5, c(1,2), ice_on_safe_func)

ice_off_10cm_0c <- apply(seasonal_ice_thick2, c(1,2), ice_off_safe_func)
ice_off_10cm_1c <- apply(seasonal_ice_thick3, c(1,2), ice_off_safe_func)
ice_off_10cm_2c <- apply(seasonal_ice_thick4, c(1,2), ice_off_safe_func)
ice_off_10cm_4c <- apply(seasonal_ice_thick5, c(1,2), ice_off_safe_func)

ice_on_15cm_0c <- apply(seasonal_ice_thick2, c(1,2), ice_on_safe_15cm_func)
ice_on_15cm_1c <- apply(seasonal_ice_thick3, c(1,2), ice_on_safe_15cm_func)
ice_on_15cm_2c <- apply(seasonal_ice_thick4, c(1,2), ice_on_safe_15cm_func)
ice_on_15cm_4c <- apply(seasonal_ice_thick5, c(1,2), ice_on_safe_15cm_func)

ice_off_15cm_0c <- apply(seasonal_ice_thick2, c(1,2), ice_off_safe_15cm_func)
ice_off_15cm_1c <- apply(seasonal_ice_thick3, c(1,2), ice_off_safe_15cm_func)
ice_off_15cm_2c <- apply(seasonal_ice_thick4, c(1,2), ice_off_safe_15cm_func)
ice_off_15cm_4c <- apply(seasonal_ice_thick5, c(1,2), ice_off_safe_15cm_func)

ice_on_20cm_0c <- apply(seasonal_ice_thick2, c(1,2), ice_on_safe_20cm_func)
ice_on_20cm_1c <- apply(seasonal_ice_thick3, c(1,2), ice_on_safe_20cm_func)
ice_on_20cm_2c <- apply(seasonal_ice_thick4, c(1,2), ice_on_safe_20cm_func)
ice_on_20cm_4c <- apply(seasonal_ice_thick5, c(1,2), ice_on_safe_20cm_func)

ice_off_20cm_0c <- apply(seasonal_ice_thick2, c(1,2), ice_off_safe_20cm_func)
ice_off_20cm_1c <- apply(seasonal_ice_thick3, c(1,2), ice_off_safe_20cm_func)
ice_off_20cm_2c <- apply(seasonal_ice_thick4, c(1,2), ice_off_safe_20cm_func)
ice_off_20cm_4c <- apply(seasonal_ice_thick5, c(1,2), ice_off_safe_20cm_func)

ice_max_0c <- apply(seasonal_ice_thick2, c(1,2), ice_max_func)
ice_max_1c <- apply(seasonal_ice_thick3, c(1,2), ice_max_func)
ice_max_2c <- apply(seasonal_ice_thick4, c(1,2), ice_max_func)
ice_max_4c <- apply(seasonal_ice_thick5, c(1,2), ice_max_func)


# 2h. Remove -1 values --------------------------------------------------------


ice_on_2cm_0c[ice_on_2cm_0c==-1] <- NA
ice_on_2cm_1c[ice_on_2cm_1c==-1] <- NA
ice_on_2cm_2c[ice_on_2cm_2c==-1] <- NA
ice_on_2cm_4c[ice_on_2cm_4c==-1] <- NA

ice_off_2cm_0c[ice_off_2cm_0c==-1] <- NA
ice_off_2cm_1c[ice_off_2cm_1c==-1] <- NA
ice_off_2cm_2c[ice_off_2cm_2c==-1] <- NA
ice_off_2cm_4c[ice_off_2cm_4c==-1] <- NA

ice_on_10cm_0c[ice_on_10cm_0c==-1] <- NA
ice_on_10cm_1c[ice_on_10cm_1c==-1] <- NA
ice_on_10cm_2c[ice_on_10cm_2c==-1] <- NA
ice_on_10cm_4c[ice_on_10cm_4c==-1] <- NA

ice_off_10cm_0c[ice_off_10cm_0c==-1] <- NA
ice_off_10cm_1c[ice_off_10cm_1c==-1] <- NA
ice_off_10cm_2c[ice_off_10cm_2c==-1] <- NA
ice_off_10cm_4c[ice_off_10cm_4c==-1] <- NA

ice_on_15cm_0c[ice_on_15cm_0c==-1] <- NA
ice_on_15cm_1c[ice_on_15cm_1c==-1] <- NA
ice_on_15cm_2c[ice_on_15cm_2c==-1] <- NA
ice_on_15cm_4c[ice_on_15cm_4c==-1] <- NA

ice_off_15cm_0c[ice_off_15cm_0c==-1] <- NA
ice_off_15cm_1c[ice_off_15cm_1c==-1] <- NA
ice_off_15cm_2c[ice_off_15cm_2c==-1] <- NA
ice_off_15cm_4c[ice_off_15cm_4c==-1] <- NA

ice_on_20cm_0c[ice_on_20cm_0c==-1] <- NA
ice_on_20cm_1c[ice_on_20cm_1c==-1] <- NA
ice_on_20cm_2c[ice_on_20cm_2c==-1] <- NA
ice_on_20cm_4c[ice_on_20cm_4c==-1] <- NA

ice_off_20cm_0c[ice_off_20cm_0c==-1] <- NA
ice_off_20cm_1c[ice_off_20cm_1c==-1] <- NA
ice_off_20cm_2c[ice_off_20cm_2c==-1] <- NA
ice_off_20cm_4c[ice_off_20cm_4c==-1] <- NA

ice_max_0c[ice_max_0c==-1] <- NA
ice_max_1c[ice_max_1c==-1] <- NA
ice_max_2c[ice_max_2c==-1] <- NA
ice_max_4c[ice_max_4c==-1] <- NA


# 2i. Calculate reaction windows ------------------------------------------

rw_ice_on_0c <- ice_on_10cm_0c - ice_on_2cm_0c
rw_ice_on_1c <- ice_on_10cm_1c - ice_on_2cm_1c
rw_ice_on_2c <- ice_on_10cm_2c - ice_on_2cm_2c
rw_ice_on_4c <- ice_on_10cm_4c - ice_on_2cm_4c

rw_ice_off_0c <- ice_off_2cm_0c - ice_off_10cm_0c
rw_ice_off_1c <- ice_off_2cm_1c - ice_off_10cm_1c
rw_ice_off_2c <- ice_off_2cm_2c - ice_off_10cm_2c
rw_ice_off_4c <- ice_off_2cm_4c - ice_off_10cm_4c

rw_ice_on_0.5wht_0c <- ice_on_15cm_0c - ice_on_2cm_0c
rw_ice_on_0.5wht_1c <- ice_on_15cm_1c - ice_on_2cm_1c
rw_ice_on_0.5wht_2c <- ice_on_15cm_2c - ice_on_2cm_2c
rw_ice_on_0.5wht_4c <- ice_on_15cm_4c - ice_on_2cm_4c

rw_ice_on_wht_0c <- ice_on_20cm_0c - ice_on_2cm_0c
rw_ice_on_wht_1c <- ice_on_20cm_1c - ice_on_2cm_1c
rw_ice_on_wht_2c <- ice_on_20cm_2c - ice_on_2cm_2c
rw_ice_on_wht_4c <- ice_on_20cm_4c - ice_on_2cm_4c

rw_ice_off_0.5wht_0c <- ice_off_2cm_0c - ice_off_15cm_0c
rw_ice_off_0.5wht_1c <- ice_off_2cm_1c - ice_off_15cm_1c
rw_ice_off_0.5wht_2c <- ice_off_2cm_2c - ice_off_15cm_2c
rw_ice_off_0.5wht_4c <- ice_off_2cm_4c - ice_off_15cm_4c

rw_ice_off_wht_0c <- ice_off_2cm_0c - ice_off_20cm_0c
rw_ice_off_wht_1c <- ice_off_2cm_1c - ice_off_20cm_1c
rw_ice_off_wht_2c <- ice_off_2cm_2c - ice_off_20cm_2c
rw_ice_off_wht_4c <- ice_off_2cm_4c - ice_off_20cm_4c

rw_ice_on_max_0c <- ice_max_0c - ice_on_10cm_0c
rw_ice_on_max_1c <- ice_max_1c - ice_on_10cm_1c
rw_ice_on_max_2c <- ice_max_2c - ice_on_10cm_2c
rw_ice_on_max_4c <- ice_max_4c - ice_on_10cm_4c

rw_ice_off_max_0c <- ice_off_10cm_0c - ice_max_0c
rw_ice_off_max_1c <- ice_off_10cm_1c - ice_max_1c
rw_ice_off_max_2c <- ice_off_10cm_2c - ice_max_2c
rw_ice_off_max_4c <- ice_off_10cm_4c - ice_max_4c

ice_duration_0c <- ice_off_2cm_0c - ice_on_2cm_0c
ice_duration_1c <- ice_off_2cm_1c - ice_on_2cm_1c
ice_duration_2c <- ice_off_2cm_2c - ice_on_2cm_2c
ice_duration_4c <- ice_off_2cm_4c - ice_on_2cm_4c

ice_duration_safe_0c <- ice_off_10cm_0c - ice_on_10cm_0c
ice_duration_safe_1c <- ice_off_10cm_1c - ice_on_10cm_1c
ice_duration_safe_2c <- ice_off_10cm_2c - ice_on_10cm_2c
ice_duration_safe_4c <- ice_off_10cm_4c - ice_on_10cm_4c

# 2j. Calculate reaction window anomalies --------------------------------

#Compare reaction windows to historical condition with black ice

#Ice on anomaly (black ice)
#NOTE: negative value = shorter shoulder season; pos val = longer shoulder season
rw_ice_on_1c_anom <- rw_ice_on_1c - rw_ice_on_0c
rw_ice_on_2c_anom <- rw_ice_on_2c - rw_ice_on_0c
rw_ice_on_4c_anom <- rw_ice_on_4c - rw_ice_on_0c

#Ice off anom (black ice)
#NOTE: neg val = shorter shoulder season; pos val = longer shoulder season
rw_ice_off_1c_anom <- rw_ice_off_1c - rw_ice_off_0c
rw_ice_off_2c_anom <- rw_ice_off_2c - rw_ice_off_0c
rw_ice_off_4c_anom <- rw_ice_off_4c - rw_ice_off_0c

#Ice on anomaly (50% white ice)
#NOTE: negative value = shorter shoulder season; pos val = longer shoulder season
rw_ice_on_1c_0.5wht_anom <- rw_ice_on_0.5wht_1c - rw_ice_on_0c
rw_ice_on_2c_0.5wht_anom <- rw_ice_on_0.5wht_2c - rw_ice_on_0c
rw_ice_on_4c_0.5wht_anom <- rw_ice_on_0.5wht_4c - rw_ice_on_0c

#Ice off anom (50% white ice)
#NOTE: neg val = shorter shoulder season; pos val = longer shoulder season
rw_ice_off_1c_0.5wht_anom <- rw_ice_off_0.5wht_1c - rw_ice_off_0c
rw_ice_off_2c_0.5wht_anom <- rw_ice_off_0.5wht_2c - rw_ice_off_0c
rw_ice_off_4c_0.5wht_anom <- rw_ice_off_0.5wht_4c - rw_ice_off_0c

#Ice on anomaly (white ice)
#NOTE: negative value = shorter shoulder season; pos val = longer shoulder season
rw_ice_on_1c_wht_anom <- rw_ice_on_wht_1c - rw_ice_on_0c
rw_ice_on_2c_wht_anom <- rw_ice_on_wht_2c - rw_ice_on_0c
rw_ice_on_4c_wht_anom <- rw_ice_on_wht_4c - rw_ice_on_0c

#Ice off anom (white ice)
#NOTE: neg val = shorter shoulder season; pos val = longer shoulder season
rw_ice_off_1c_wht_anom <- rw_ice_off_wht_1c - rw_ice_off_0c
rw_ice_off_2c_wht_anom <- rw_ice_off_wht_2c - rw_ice_off_0c
rw_ice_off_4c_wht_anom <- rw_ice_off_wht_4c - rw_ice_off_0c


# 3a. Convert reaction window to tibble -----------------------------------

rw_abs_val <- tibble(
  warming = c(
    as.factor(c('0\u00B0C', '1\u00B0C', '2\u00B0C', '4\u00B0C')) #This code (\u00B0) adds a degree symbol
  ),
  ice_on_doy = c(
    mean(ice_on_2cm_0c, na.rm = T), 
    mean(ice_on_2cm_1c, na.rm = T), 
    mean(ice_on_2cm_2c, na.rm = T), 
    mean(ice_on_2cm_4c, na.rm = T)
  ),
  ice_off_doy = c(
    mean(ice_off_2cm_0c, na.rm = T), 
    mean(ice_off_2cm_1c, na.rm = T), 
    mean(ice_off_2cm_2c, na.rm = T), 
    mean(ice_off_2cm_4c, na.rm = T)
  ),
  ice_max_doy =  c(
    mean(ice_max_0c, na.rm = T), 
    mean(ice_max_1c, na.rm = T), 
    mean(ice_max_2c, na.rm = T), 
    mean(ice_max_4c, na.rm = T)
  ),
  rw_on_days = c(
    mean(rw_ice_on_0c, na.rm = T),
    mean(rw_ice_on_1c, na.rm = T),
    mean(rw_ice_on_2c, na.rm = T),
    mean(rw_ice_on_4c, na.rm = T)
  ),
  rw_on_sd = c(
    mean(rw_ice_on_0c, na.rm = T),
    mean(rw_ice_on_1c, na.rm = T),
    mean(rw_ice_on_2c, na.rm = T),
    mean(rw_ice_on_4c, na.rm = T)
  ),
  rw_off_days = c(
    mean(rw_ice_off_0c, na.rm = T),
    mean(rw_ice_off_1c, na.rm = T),
    mean(rw_ice_off_2c, na.rm = T),
    mean(rw_ice_off_4c, na.rm = T)
  ),
  rw_on_max_days = c(
    mean(rw_ice_on_max_0c, na.rm = T),
    mean(rw_ice_on_max_1c, na.rm = T),
    mean(rw_ice_on_max_2c, na.rm = T),
    mean(rw_ice_on_max_4c, na.rm = T)
  ),
  rw_off_max_days = c(
    mean(rw_ice_off_max_0c, na.rm = T),
    mean(rw_ice_off_max_1c, na.rm = T),
    mean(rw_ice_off_max_2c, na.rm = T),
    mean(rw_ice_off_max_4c, na.rm = T)
  ),
  ice_duration_days = c(
    mean(ice_duration_0c, na.rm = T),
    mean(ice_duration_1c, na.rm = T),
    mean(ice_duration_2c, na.rm = T),
    mean(ice_duration_4c, na.rm = T)
  ),
  ice_duration_safe_days = c(
    mean(ice_duration_safe_0c, na.rm = T),
    mean(ice_duration_safe_1c, na.rm = T),
    mean(ice_duration_safe_2c, na.rm = T),
    mean(ice_duration_safe_4c, na.rm = T)
  )
)
rw_abs_val

rw_abs_val_long <- rw_abs_val %>% 
  pivot_longer(
    !warming, 
    names_to = 'phenology',
    values_to = 'value'
  ) %>% 
  arrange(
    phenology
  ) %>% 
  mutate(
    phenology = as.factor(phenology)
  )

# 3b. Convert RW anomaly to tibble ----------------------------------------

rw_anom <- tibble(
  warming = c(
    as.factor(c('1\u00B0C', '2\u00B0C', '4\u00B0C')) #This code ("\u00B0") adds a degree symbol
  ),
  rw_on_blk_anom_days = c(
    mean(rw_ice_on_1c_anom, na.rm = T),
    mean(rw_ice_on_2c_anom, na.rm = T),
    mean(rw_ice_on_4c_anom, na.rm = T)
  ),
  rw_on_blk_nval = c(
    sum(!is.na(rw_ice_on_1c_anom)),
    sum(!is.na(rw_ice_on_2c_anom)),
    sum(!is.na(rw_ice_on_4c_anom))
  ),
  rw_off_blk_anom_days = c(
    mean(rw_ice_off_1c_anom, na.rm = T),
    mean(rw_ice_off_2c_anom, na.rm = T),
    mean(rw_ice_off_4c_anom, na.rm = T)
  ),
  rw_off_blk_nval = c(
    sum(!is.na(rw_ice_off_1c_anom)),
    sum(!is.na(rw_ice_off_2c_anom)),
    sum(!is.na(rw_ice_off_4c_anom))
  ),
  rw_on_0.5wht_anom_days = c(
    mean(rw_ice_on_1c_0.5wht_anom, na.rm = T),
    mean(rw_ice_on_2c_0.5wht_anom, na.rm = T),
    mean(rw_ice_on_4c_0.5wht_anom, na.rm = T)
  ),
  rw_on_0.5wht_nval = c(
    sum(!is.na(rw_ice_on_1c_0.5wht_anom)),
    sum(!is.na(rw_ice_on_2c_0.5wht_anom)),
    sum(!is.na(rw_ice_on_4c_0.5wht_anom))
  ),
  rw_off_0.5wht_anom_days = c(
    mean(rw_ice_off_1c_0.5wht_anom, na.rm = T),
    mean(rw_ice_off_2c_0.5wht_anom, na.rm = T),
    mean(rw_ice_off_4c_0.5wht_anom, na.rm = T)
  ),
  rw_off_0.5wht_nval = c(
    sum(!is.na(rw_ice_off_1c_0.5wht_anom)),
    sum(!is.na(rw_ice_off_2c_0.5wht_anom)),
    sum(!is.na(rw_ice_off_4c_0.5wht_anom))
  ),
  rw_on_wht_anom_days = c(
    mean(rw_ice_on_1c_wht_anom, na.rm = T),
    mean(rw_ice_on_2c_wht_anom, na.rm = T),
    mean(rw_ice_on_4c_wht_anom, na.rm = T)
  ),
  rw_on_wht_nval = c(
    sum(!is.na(rw_ice_on_1c_wht_anom)),
    sum(!is.na(rw_ice_on_2c_wht_anom)),
    sum(!is.na(rw_ice_on_4c_wht_anom))
  ),
  rw_off_wht_anom_days = c(
    mean(rw_ice_off_1c_wht_anom, na.rm = T),
    mean(rw_ice_off_2c_wht_anom, na.rm = T),
    mean(rw_ice_off_4c_wht_anom, na.rm = T)
  ),
  rw_off_wht_nval = c(
    sum(!is.na(rw_ice_off_1c_wht_anom)),
    sum(!is.na(rw_ice_off_2c_wht_anom)),
    sum(!is.na(rw_ice_off_4c_wht_anom))
  )
)

head(rw_anom)

rw_anom_long <- rw_anom %>% 
  pivot_longer(
    !warming, 
    names_to = 'phenology',
    values_to = 'value'
  ) %>% 
  arrange(
    phenology
  ) %>% 
  mutate(
    phenology = as.factor(phenology)
  )


# 3c. Standard deviations for reaction window anomalies -------------------

rw_anom_sd <- tibble(
  warming = c(
    as.factor(c('1\u00B0C', '2\u00B0C', '4\u00B0C')) #This code ("\u00B0") adds a degree symbol
  ),
  rw_on_blk_anom_days_sd = c(
    sd(rw_ice_on_1c_anom, na.rm = T),
    sd(rw_ice_on_2c_anom, na.rm = T),
    sd(rw_ice_on_4c_anom, na.rm = T)
  ),
  rw_off_blk_anom_days_sd = c(
    sd(rw_ice_off_1c_anom, na.rm = T),
    sd(rw_ice_off_2c_anom, na.rm = T),
    sd(rw_ice_off_4c_anom, na.rm = T)
  ),
  rw_on_0.5wht_anom_days_sd = c(
    sd(rw_ice_on_1c_0.5wht_anom, na.rm = T),
    sd(rw_ice_on_2c_0.5wht_anom, na.rm = T),
    sd(rw_ice_on_4c_0.5wht_anom, na.rm = T)
  ),
  rw_off_0.5wht_anom_days_sd = c(
    sd(rw_ice_off_1c_0.5wht_anom, na.rm = T),
    sd(rw_ice_off_2c_0.5wht_anom, na.rm = T),
    sd(rw_ice_off_4c_0.5wht_anom, na.rm = T)
  ),
  rw_on_wht_anom_days_sd = c(
    sd(rw_ice_on_1c_wht_anom, na.rm = T),
    sd(rw_ice_on_2c_wht_anom, na.rm = T),
    sd(rw_ice_on_4c_wht_anom, na.rm = T)
  ),
  rw_off_wht_anom_days_sd = c(
    sd(rw_ice_off_1c_wht_anom, na.rm = T),
    sd(rw_ice_off_2c_wht_anom, na.rm = T),
    sd(rw_ice_off_4c_wht_anom, na.rm = T)
  )
)

head(rw_anom_sd)

rw_anom_sd_long <- rw_anom_sd %>% 
  pivot_longer(
    !warming, 
    names_to = 'phenology',
    values_to = 'value'
  ) %>% 
  arrange(
    phenology
  ) %>% 
  mutate(
    phenology = as.factor(phenology)
  )


# 5. RW boxplots ----------------------------------------------------------

# 5a. Convert grids into a dataframe --------------------------------------

#**100% black ice----
rw_ice_on_1c_blk_anom_df <- melt(rw_ice_on_1c_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('1\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('black')
    ),
    event = as.factor('on')
  )

rw_ice_on_2c_blk_anom_df <- melt(rw_ice_on_2c_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('2\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('black')
    ),
    event = as.factor('on')
  )

rw_ice_on_4c_blk_anom_df <- melt(rw_ice_on_4c_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('4\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('black')
    ),
    event = as.factor('on')
  )

rw_ice_off_1c_blk_anom_df <- melt(rw_ice_off_1c_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('1\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('black')
    ),
    event = as.factor('off')
  )

rw_ice_off_2c_blk_anom_df <- melt(rw_ice_off_2c_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('2\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('black')
    ),
    event = as.factor('off')
  )

rw_ice_off_4c_blk_anom_df <- melt(rw_ice_off_4c_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('4\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('black')
    ),
    event = as.factor('off')
  )

# **50% white ice----
rw_ice_on_1c_0.5wht_anom_df <- melt(rw_ice_on_1c_0.5wht_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('1\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('mix')
    ),
    event = as.factor('on')
  )

rw_ice_on_2c_0.5wht_anom_df <- melt(rw_ice_on_2c_0.5wht_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('2\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('mix')
    ),
    event = as.factor('on')
  )

rw_ice_on_4c_0.5wht_anom_df <- melt(rw_ice_on_4c_0.5wht_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('4\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('mix')
    ),
    event = as.factor('on')
  )

rw_ice_off_1c_0.5wht_anom_df <- melt(rw_ice_off_1c_0.5wht_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('1\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('mix')
    ),
    event = as.factor('off')
  )

rw_ice_off_2c_0.5wht_anom_df <- melt(rw_ice_off_2c_0.5wht_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('2\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('mix')
    ),
    event = as.factor('off')
  )

rw_ice_off_4c_0.5wht_anom_df <- melt(rw_ice_off_4c_0.5wht_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('4\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('mix')
    ),
    event = as.factor('off')
  )

# **100% white ice----
rw_ice_on_1c_wht_anom_df <- melt(rw_ice_on_1c_wht_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('1\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('white')
    ),
    event = as.factor('on')
  )

rw_ice_on_2c_wht_anom_df <- melt(rw_ice_on_2c_wht_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('2\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('white')
    ),
    event = as.factor('on')
  )

rw_ice_on_4c_wht_anom_df <- melt(rw_ice_on_4c_wht_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('4\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('white')
    ),
    event = as.factor('on')
  )

rw_ice_off_1c_wht_anom_df <- melt(rw_ice_off_1c_wht_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('1\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('white')
    ),
    event = as.factor('off')
  )

rw_ice_off_2c_wht_anom_df <- melt(rw_ice_off_2c_wht_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('2\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = c(
      as.factor('white')
    ),
    event = as.factor('off')
  )

rw_ice_off_4c_wht_anom_df <- melt(rw_ice_off_4c_wht_anom) %>% 
  select(
    lon = 1,
    lat = 2, 
    anomaly_days = 3
  ) %>% 
  mutate(
    warming = c(
      as.factor('4\u00B0C') #This code ("\u00B0") adds a degree symbol
    ),
    ice_type = as.factor('white'),
    event = as.factor('off')
  )


# **Clip outliers (>2 std dev) --------------------------------------------

#Function to get standard deviations for clipping the data
#Removes points greater than 2 std. devs. 

sd_func <- function(df_c, df){
  #df is the column we're intersted in: df$anomaly_days
  
  #First take the mean of the anomaly_days
  mean = mean(df_c, na.rm = T);
  
  #Then take the standard deviation of the anomaly_days
  sd = sd(df_c, na.rm = T);
  
  #Get the upper end of the values
  upper = mean + (2*sd)
  
  #Get the lower value
  lower = mean - (2*sd)
  
  #create a simple list
  x <- c(lower,upper)
  
  #Filter the dataframe by the upper and lower boundaries
  df %>% filter(
    df_c > x[1] & df_c < x[2]
    )

}

# Ice on - Black ice clipping
rw_ice_on_1c_blk_anom_df_cut <- sd_func(df_c = rw_ice_on_1c_blk_anom_df$anomaly_days, 
                                        df = rw_ice_on_1c_blk_anom_df)

rw_ice_on_2c_blk_anom_df_cut <- sd_func(df_c = rw_ice_on_2c_blk_anom_df$anomaly_days, 
                                        df = rw_ice_on_2c_blk_anom_df)

rw_ice_on_4c_blk_anom_df_cut <- sd_func(df_c = rw_ice_on_4c_blk_anom_df$anomaly_days, 
                                        df = rw_ice_on_4c_blk_anom_df)

# Ice on - 50% White ice clipping
rw_ice_on_1c_0.5wht_anom_df_cut <- sd_func(df_c = rw_ice_on_1c_0.5wht_anom_df$anomaly_days, 
                                        df = rw_ice_on_1c_0.5wht_anom_df)

rw_ice_on_2c_0.5wht_anom_df_cut <- sd_func(df_c = rw_ice_on_2c_0.5wht_anom_df$anomaly_days, 
                                        df = rw_ice_on_2c_0.5wht_anom_df)

rw_ice_on_4c_0.5wht_anom_df_cut <- sd_func(df_c = rw_ice_on_4c_0.5wht_anom_df$anomaly_days, 
                                        df = rw_ice_on_4c_0.5wht_anom_df)

#Ice on - 100% White ice clipping
rw_ice_on_1c_wht_anom_df_cut <- sd_func(df_c = rw_ice_on_1c_wht_anom_df$anomaly_days, 
                                           df = rw_ice_on_1c_wht_anom_df)

rw_ice_on_2c_wht_anom_df_cut <- sd_func(df_c = rw_ice_on_2c_wht_anom_df$anomaly_days, 
                                           df = rw_ice_on_2c_wht_anom_df)

rw_ice_on_4c_wht_anom_df_cut <- sd_func(df_c = rw_ice_on_4c_wht_anom_df$anomaly_days, 
                                           df = rw_ice_on_4c_wht_anom_df)

# Ice off - Black ice clipping
rw_ice_off_1c_blk_anom_df_cut <- sd_func(df_c = rw_ice_off_1c_blk_anom_df$anomaly_days, 
                                        df = rw_ice_off_1c_blk_anom_df)

rw_ice_off_2c_blk_anom_df_cut <- sd_func(df_c = rw_ice_off_2c_blk_anom_df$anomaly_days, 
                                        df = rw_ice_off_2c_blk_anom_df)

rw_ice_off_4c_blk_anom_df_cut <- sd_func(df_c = rw_ice_off_4c_blk_anom_df$anomaly_days, 
                                        df = rw_ice_off_4c_blk_anom_df)

# Ice off - 50% White ice clipping
rw_ice_off_1c_0.5wht_anom_df_cut <- sd_func(df_c = rw_ice_off_1c_0.5wht_anom_df$anomaly_days, 
                                           df = rw_ice_off_1c_0.5wht_anom_df)

rw_ice_off_2c_0.5wht_anom_df_cut <- sd_func(df_c = rw_ice_off_2c_0.5wht_anom_df$anomaly_days, 
                                           df = rw_ice_off_2c_0.5wht_anom_df)

rw_ice_off_4c_0.5wht_anom_df_cut <- sd_func(df_c = rw_ice_off_4c_0.5wht_anom_df$anomaly_days, 
                                           df = rw_ice_off_4c_0.5wht_anom_df)

#Ice off - 100% White ice clipping
rw_ice_off_1c_wht_anom_df_cut <- sd_func(df_c = rw_ice_off_1c_wht_anom_df$anomaly_days, 
                                        df = rw_ice_off_1c_wht_anom_df)

rw_ice_off_2c_wht_anom_df_cut <- sd_func(df_c = rw_ice_off_2c_wht_anom_df$anomaly_days, 
                                        df = rw_ice_off_2c_wht_anom_df)

rw_ice_off_4c_wht_anom_df_cut <- sd_func(df_c = rw_ice_off_4c_wht_anom_df$anomaly_days, 
                                        df = rw_ice_off_4c_wht_anom_df)
# **Combine DFs -----------------------------------------------------------

ice_safety_df <- rw_ice_on_1c_blk_anom_df_cut %>% 
  bind_rows(
    #rw ice on
    rw_ice_on_2c_blk_anom_df_cut,
    rw_ice_on_4c_blk_anom_df_cut,
    rw_ice_on_1c_0.5wht_anom_df_cut,
    rw_ice_on_2c_0.5wht_anom_df_cut,
    rw_ice_on_4c_0.5wht_anom_df_cut,
    rw_ice_on_1c_wht_anom_df_cut,
    rw_ice_on_2c_wht_anom_df_cut,
    rw_ice_on_4c_wht_anom_df_cut,
    #rw ice off
    rw_ice_off_1c_blk_anom_df_cut,
    rw_ice_off_2c_blk_anom_df_cut,
    rw_ice_off_4c_blk_anom_df_cut,
    rw_ice_off_1c_0.5wht_anom_df_cut,
    rw_ice_off_2c_0.5wht_anom_df_cut,
    rw_ice_off_4c_0.5wht_anom_df_cut,
    rw_ice_off_1c_wht_anom_df_cut,
    rw_ice_off_2c_wht_anom_df_cut,
    rw_ice_off_4c_wht_anom_df_cut
  )

# **5b. Split violin plot -------------------------------------------------

ggplot(data = ice_safety_df, aes(x = warming, y = anomaly_days, fill = event))+ 
  geom_split_violin(
    width = 1,
    scale = 'count',
    adjust = 4,
    alpha = 0.4)+
  geom_boxplot(width = 0.2,
               alpha = 0.6,
               fatten = NULL,
               show.legend = NULL,
               outlier.alpha = 0)+
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(0.2)) +
  facet_wrap(~ice_type,
             #scales = 'free',
             labeller = as_labeller(c(
               `black` = "100% Black Ice", 
               `mix` = "50/50 Ice Ratio", 
               `white` = "100% White Ice")))+
  theme_classic()+
  theme(
    legend.title = element_blank(),
    legend.position = 'bottom',
    legend.text = element_text(size = 15),
    axis.text = element_text(size = 20),
    axis.title.y = element_text(size = 22),
    strip.text = element_text(size = 20)
  )+
  xlab('')+
  ylab('Transition Period Anomaly (days)')+
  scale_fill_manual(values=c( "#F5A773FF", "#B32E5FFF"),labels = c('Formation', 'Melt'))

#ggsave(here('results/transition_period_split_violin_filter_2024.03.19_w_box.png'), dpi = 300, width = 5.2, height = 3.467, units = 'in')

# 6. Means tests ----------------------------------------------------------

#**Ice on normality test----

#Recall: if p-val > 0.05, then distribution is normal

#***black ice----
shapiro.test(rw_ice_on_1c_blk_anom_df_cut$anomaly_days) #p.val<0.05, n = 2734
shapiro.test(rw_ice_on_2c_blk_anom_df_cut$anomaly_days) #p.val<0.05, n = 2641
shapiro.test(rw_ice_on_4c_blk_anom_df_cut$anomaly_days) #p.val<0.05, n = 2529

#Get n value for each test
#length(rw_ice_on_4c_blk_anom_df_cut$anomaly_days)

#***50/50----
shapiro.test(rw_ice_on_1c_0.5wht_anom_df_cut$anomaly_days) #p.val<0.05, n = 2616
shapiro.test(rw_ice_on_2c_0.5wht_anom_df_cut$anomaly_days) #p.val<0.05, n = 2559
shapiro.test(rw_ice_on_4c_0.5wht_anom_df_cut$anomaly_days) #p.val<0.05, n = 2404

#Get n value for each test
#length(rw_ice_on_4c_0.5wht_anom_df_cut$anomaly_days)

#***white ice----
shapiro.test(rw_ice_on_1c_wht_anom_df_cut$anomaly_days) #p.val<0.05, n = 2537
shapiro.test(rw_ice_on_2c_wht_anom_df_cut$anomaly_days) #p.val<0.05, n = 2464
shapiro.test(rw_ice_on_4c_wht_anom_df_cut$anomaly_days) #p.val<0.05, n = 2258

#NOTE: None of the ice-on distributions are normally distributed

#**Ice off normality test----

#***black ice----
shapiro.test(rw_ice_off_1c_blk_anom_df_cut$anomaly_days) #p.val<0.05
shapiro.test(rw_ice_off_2c_blk_anom_df_cut$anomaly_days) #p.val<0.05
shapiro.test(rw_ice_off_4c_blk_anom_df_cut$anomaly_days) #p.val<0.05

#Get n value for each test
#length(rw_ice_off_1c_blk_anom_df_cut$anomaly_days)

#***50/50----
shapiro.test(rw_ice_off_1c_0.5wht_anom_df_cut$anomaly_days) #p.val<0.05
shapiro.test(rw_ice_off_2c_0.5wht_anom_df_cut$anomaly_days) #p.val<0.05
shapiro.test(rw_ice_off_4c_0.5wht_anom_df_cut$anomaly_days) #p.val<0.05

#Get n value for each test
length(rw_ice_off_4c_0.5wht_anom_df_cut$anomaly_days)

#***white ice----
shapiro.test(rw_ice_off_1c_wht_anom_df_cut$anomaly_days) #p.val<0.05
shapiro.test(rw_ice_off_2c_wht_anom_df_cut$anomaly_days) #p.val<0.05
shapiro.test(rw_ice_off_4c_wht_anom_df_cut$anomaly_days) #p.val<0.05

#Get n value for each test
length(rw_ice_off_4c_wht_anom_df_cut$anomaly_days)

#NOTE: None of the ice-off distributions are normally distributed.

#**Ice on means test within ice quality----

#Recall: if p-val < 0.05, then at least one warming scenario is different

#NOTE: To test means across warming or ice scenarios, create data frames
#      that have data combined either across warming scenarios or ice quality (as below)

#**Ice on within black ice category

ice_on_blk <- rw_ice_on_1c_blk_anom_df_cut %>% 
  bind_rows(
    rw_ice_on_2c_blk_anom_df_cut,
    rw_ice_on_4c_blk_anom_df_cut
  )

ice_on_50 <- rw_ice_on_1c_0.5wht_anom_df_cut %>% 
  bind_rows(
    rw_ice_on_2c_0.5wht_anom_df_cut,
    rw_ice_on_4c_0.5wht_anom_df_cut
  )

ice_on_wht <- rw_ice_on_1c_wht_anom_df_cut %>% 
  bind_rows(
    rw_ice_on_2c_wht_anom_df_cut,
    rw_ice_on_4c_wht_anom_df_cut
  )

ice_on_1C_diff_qual <- rw_ice_on_1c_blk_anom_df_cut %>% 
  bind_rows(
    rw_ice_on_1c_0.5wht_anom_df_cut,
    rw_ice_on_1c_wht_anom_df_cut
  )

ice_on_2C_diff_qual <- rw_ice_on_2c_blk_anom_df_cut %>% 
  bind_rows(
    rw_ice_on_2c_0.5wht_anom_df_cut,
    rw_ice_on_2c_wht_anom_df_cut
  )

ice_on_4C_diff_qual <- rw_ice_on_4c_blk_anom_df_cut %>% 
  bind_rows(
    rw_ice_on_4c_0.5wht_anom_df_cut,
    rw_ice_on_4c_wht_anom_df_cut
  )

#Now bind ice on across quality within the same temp increase
#Then do the same for ice off

#***black ice----
kruskal.test(
  data = ice_on_blk,
  anomaly_days ~ warming,
) #p.val<0.05, therefore at least one mean is different

#posthoc Dunn test to determine which warming scenario is different

dunnTest(
  data = ice_on_blk,
  anomaly_days ~ warming,
  method = 'holm'
  #method = 'hochberg'
) #All mean anomalies among warming scenarios within black ice quality are significantly different from one another

#
#***50/50----
kruskal.test(
  data = ice_on_50,
  anomaly_days ~ warming,
) #p.val<0.05, therefore at least one mean is different

#posthoc Dunn test to determine which warming scenario is different

dunnTest(
  data = ice_on_50,
  anomaly_days ~ warming,
  method = 'holm'
) #All mean anomalies among warming scenarios within black ice quality are significantly different from one another

#***white ice----
kruskal.test(
  data = ice_on_wht,
  anomaly_days ~ warming,
) #p.val<0.05, therefore at least one mean is different

#posthoc Dunn test to determine which warming scenario is different

dunnTest(
  data = ice_on_wht,
  anomaly_days ~ warming,
  method = 'holm'
) #All mean anomalies among warming scenarios within black ice quality are significantly different from one another

#**Ice off means test within ice quality----

#NOTE: To test means across warming or ice scenarios,create data frames
#      that have data combines either across warming scenarios or ice quality (as below)

#**Ice off within black ice category

ice_off_blk <- rw_ice_off_1c_blk_anom_df_cut %>% 
  bind_rows(
    rw_ice_off_2c_blk_anom_df_cut,
    rw_ice_off_4c_blk_anom_df_cut
  )

ice_off_50 <- rw_ice_off_1c_0.5wht_anom_df_cut %>% 
  bind_rows(
    rw_ice_off_2c_0.5wht_anom_df_cut,
    rw_ice_off_4c_0.5wht_anom_df_cut
  )

ice_off_wht <- rw_ice_off_1c_wht_anom_df_cut %>% 
  bind_rows(
    rw_ice_off_2c_wht_anom_df_cut,
    rw_ice_off_4c_wht_anom_df_cut
  )

ice_off_1C_diff_qual <- rw_ice_off_1c_blk_anom_df_cut %>% 
  bind_rows(
    rw_ice_off_1c_0.5wht_anom_df_cut,
    rw_ice_off_1c_wht_anom_df_cut
  )

ice_off_2C_diff_qual <- rw_ice_off_2c_blk_anom_df_cut %>% 
  bind_rows(
    rw_ice_off_2c_0.5wht_anom_df_cut,
    rw_ice_off_2c_wht_anom_df_cut
  )

ice_off_4C_diff_qual <- rw_ice_off_4c_blk_anom_df_cut %>% 
  bind_rows(
    rw_ice_off_4c_0.5wht_anom_df_cut,
    rw_ice_off_4c_wht_anom_df_cut
  )

#***black ice----
kruskal.test(
  data = ice_off_blk,
  anomaly_days ~ warming,
) #p.val<0.05, therefore at least one mean is different

#posthoc Dunn test to determine which warming scenario is different

dunnTest(
  data = ice_off_blk,
  anomaly_days ~ warming,
  method = 'holm'
) #All mean anomalies among warming scenarios within black ice quality are significantly different from one another

#
#***50/50----
kruskal.test(
  data = ice_off_50,
  anomaly_days ~ warming,
) #p.val<0.05, therefore at least one mean is different

#posthoc Dunn test to determine which warming scenario is different

dunnTest(
  data = ice_off_50,
  anomaly_days ~ warming,
  method = 'holm'
) #All mean anomalies among warming scenarios within black ice quality are significantly different from one another

#***white ice----
kruskal.test(
  data = ice_off_wht,
  anomaly_days ~ warming,
) #p.val<0.05, therefore at least one mean is different

#posthoc Dunn test to determine which warming scenario is different

dunnTest(
  data = ice_off_wht,
  anomaly_days ~ warming,
  method = 'holm'
) #All mean anomalies among warming scenarios within black ice quality are significantly different from one another

#**Ice on means test between ice quality----

#***1C----
kruskal.test(
  data = ice_on_1C_diff_qual,
  anomaly_days ~ ice_type,
) #p.val<0.05, therefore at least one mean is different

#posthoc Dunn test to determine which warming scenario is different

dunnTest(
  data = ice_on_1C_diff_qual,
  anomaly_days ~ ice_type,
  method = 'holm'
) #All mean anomalies among warming scenarios within black ice quality are significantly different from one another

#
#***2C----
kruskal.test(
  data = ice_on_2C_diff_qual,
  anomaly_days ~ ice_type,
) #p.val<0.05, therefore at least one mean is different

#posthoc Dunn test to determine which warming scenario is different

dunnTest(
  data = ice_on_2C_diff_qual,
  anomaly_days ~ ice_type,
  method = 'holm'
) #All mean anomalies among warming scenarios within black ice quality are significantly different from one another

#***4C----
kruskal.test(
  data = ice_on_4C_diff_qual,
  anomaly_days ~ ice_type,
) #p.val<0.05, therefore at least one mean is different

#posthoc Dunn test to determine which warming scenario is different

dunnTest(
  data = ice_on_4C_diff_qual,
  anomaly_days ~ ice_type,
  method = 'holm'
) #All mean anomalies among warming scenarios within black ice quality are significantly different from one another

#**Ice off means test between ice quality----

#***1C----
kruskal.test(
  data = ice_off_1C_diff_qual,
  anomaly_days ~ ice_type,
) #p.val<0.05, therefore at least one mean is different

#posthoc Dunn test to determine which warming scenario is different

dunnTest(
  data = ice_off_1C_diff_qual,
  anomaly_days ~ ice_type,
  method = 'holm'
) #All mean anomalies among warming scenarios within black ice quality are significantly different from one another

#
#***2C----
kruskal.test(
  data = ice_off_2C_diff_qual,
  anomaly_days ~ ice_type,
) #p.val<0.05, therefore at least one mean is different

#posthoc Dunn test to determine which warming scenario is different

dunnTest(
  data = ice_off_2C_diff_qual,
  anomaly_days ~ ice_type,
  method = 'holm'
) #All mean anomalies among warming scenarios within black ice quality are significantly different from one another

#***4C----
kruskal.test(
  data = ice_off_4C_diff_qual,
  anomaly_days ~ ice_type,
) #p.val<0.05, therefore at least one mean is different

#posthoc Dunn test to determine which warming scenario is different

dunnTest(
  data = ice_off_4C_diff_qual,
  anomaly_days ~ ice_type,
  method = 'holm'
) #All mean anomalies among warming scenarios within black ice quality are significantly different from one another


# 7. Kolmogorov-Smirnov test for distribution -----------------------------

#**7a. Ice on across temps, within quality----

#NOTE: p.val < 0.05 implies that the sample data does not come from the same distribution

#Ice on within 100% black ice
ks.test(
  rw_ice_on_1c_blk_anom_df_cut$anomaly_days, 
  rw_ice_on_2c_blk_anom_df_cut$anomaly_days
  ) #p.val < 0.05
ks.test(
  rw_ice_on_1c_blk_anom_df_cut$anomaly_days, 
  rw_ice_on_4c_blk_anom_df_cut$anomaly_days
  ) #p.val < 0.05
ks.test(
  rw_ice_on_2c_blk_anom_df_cut$anomaly_days, 
  rw_ice_on_4c_blk_anom_df_cut$anomaly_days
) #p.val < 0.05

#Ice on within 50% white ice
ks.test(
  rw_ice_on_1c_0.5wht_anom_df_cut$anomaly_days, 
  rw_ice_on_2c_0.5wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_on_1c_0.5wht_anom_df_cut$anomaly_days, 
  rw_ice_on_4c_0.5wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_on_2c_0.5wht_anom_df_cut$anomaly_days, 
  rw_ice_on_4c_0.5wht_anom_df_cut$anomaly_days
) #p.val < 0.05

#Ice on within 100% white ice
ks.test(
  rw_ice_on_1c_wht_anom_df_cut$anomaly_days, 
  rw_ice_on_2c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_on_1c_wht_anom_df_cut$anomaly_days, 
  rw_ice_on_4c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_on_2c_wht_anom_df_cut$anomaly_days, 
  rw_ice_on_4c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05

#**7b. Ice on across quality, within same warming----

#NOTE: p.val < 0.05 means that the sample data does not come from the same distribtion

#Ice on within 1C warming
ks.test(
  rw_ice_on_1c_blk_anom_df_cut$anomaly_days, 
  rw_ice_on_1c_0.5wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_on_1c_blk_anom_df_cut$anomaly_days, 
  rw_ice_on_1c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_on_1c_0.5wht_anom_df_cut$anomaly_days,
  rw_ice_on_1c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05

#Ice on within 2C warming
ks.test(
  rw_ice_on_2c_blk_anom_df_cut$anomaly_days, 
  rw_ice_on_2c_0.5wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_on_2c_blk_anom_df_cut$anomaly_days, 
  rw_ice_on_2c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_on_2c_0.5wht_anom_df_cut$anomaly_days,
  rw_ice_on_2c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05


#Ice on within 4C warming
ks.test(
  rw_ice_on_4c_blk_anom_df_cut$anomaly_days, 
  rw_ice_on_4c_0.5wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_on_4c_blk_anom_df_cut$anomaly_days, 
  rw_ice_on_4c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_on_4c_0.5wht_anom_df_cut$anomaly_days,
  rw_ice_on_4c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05

#**7c. Ice off across temps, within quality----

#NOTE: p.val < 0.05 means that the sample data does not come from the same distribtion

#Ice off within 100% black ice
ks.test(
  rw_ice_off_1c_blk_anom_df_cut$anomaly_days, 
  rw_ice_off_2c_blk_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_off_1c_blk_anom_df_cut$anomaly_days, 
  rw_ice_off_4c_blk_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_off_2c_blk_anom_df_cut$anomaly_days, 
  rw_ice_off_4c_blk_anom_df_cut$anomaly_days
) #p.val < 0.05

#Ice off within 50% white ice
ks.test(
  rw_ice_off_1c_0.5wht_anom_df_cut$anomaly_days, 
  rw_ice_off_2c_0.5wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_off_1c_0.5wht_anom_df_cut$anomaly_days, 
  rw_ice_off_4c_0.5wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_off_2c_0.5wht_anom_df_cut$anomaly_days, 
  rw_ice_off_4c_0.5wht_anom_df_cut$anomaly_days
) #p.val < 0.05

#Ice off within 100% white ice
ks.test(
  rw_ice_off_1c_wht_anom_df_cut$anomaly_days, 
  rw_ice_off_2c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_off_1c_wht_anom_df_cut$anomaly_days, 
  rw_ice_off_4c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_off_2c_wht_anom_df_cut$anomaly_days, 
  rw_ice_off_4c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05

#**7d. Ice off across quality, within same warming----

#NOTE: p.val < 0.05 means that the sample data does not come from the same distribtion

#Ice off within 1C warming
ks.test(
  rw_ice_off_1c_blk_anom_df_cut$anomaly_days, 
  rw_ice_off_1c_0.5wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_off_1c_blk_anom_df_cut$anomaly_days, 
  rw_ice_off_1c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_off_1c_0.5wht_anom_df_cut$anomaly_days,
  rw_ice_off_1c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05

#Ice on within 2C warming
ks.test(
  rw_ice_off_2c_blk_anom_df_cut$anomaly_days, 
  rw_ice_off_2c_0.5wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_off_2c_blk_anom_df_cut$anomaly_days, 
  rw_ice_off_2c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_off_2c_0.5wht_anom_df_cut$anomaly_days,
  rw_ice_off_2c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05


#Ice off within 4C warming
ks.test(
  rw_ice_off_4c_blk_anom_df_cut$anomaly_days, 
  rw_ice_off_4c_0.5wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_off_4c_blk_anom_df_cut$anomaly_days, 
  rw_ice_off_4c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05
ks.test(
  rw_ice_off_4c_0.5wht_anom_df_cut$anomaly_days,
  rw_ice_off_4c_wht_anom_df_cut$anomaly_days
) #p.val < 0.05


# plot(ecdf(rw_ice_off_1c_blk_anom_df_cut$anomaly_days), verticals = TRUE, do.points = FALSE, col.01line = NULL, xlab = "", main = "ECDFs")
# lines(ecdf(rw_ice_off_4c_wht_anom_df_cut$anomaly_days), verticals = TRUE, do.points = FALSE, col.01line = NULL, col = 4)



# 8. ECDF Plots --------------------------------------------------------------


test_diff_qual <- ice_on_1C_diff_qual %>% 
  bind_rows(ice_on_2C_diff_qual, ice_on_4C_diff_qual,
            ice_off_1C_diff_qual, ice_off_2C_diff_qual, ice_off_4C_diff_qual) %>% 
  mutate(
    ice_type = as.factor(if_else(ice_type == 'mix', 'ice_mix', ice_type))
  )

test_diff_warm <- ice_on_blk %>% 
  bind_rows(ice_on_50, ice_on_wht,
            ice_off_blk, ice_off_50, ice_off_wht) %>% 
  mutate(
    ice_type = as.factor(if_else(ice_type == 'mix', 'ice_mix', ice_type))
  )

p0 <- ggplot(data = test_diff_qual, aes(x = anomaly_days, color = ice_type))
p0 + stat_ecdf() +
  theme_classic() +
  scale_color_viridis_d(begin = 0.3, end = 0.8)+
  theme(
    legend.title = element_blank()
  )+
  labs(x = 'Transition Period Anomaly (days)',
       y = 'ECDF (%)')+
  facet_wrap(~warming + event, ncol = 2)

#ggsave(here('results/ecdf/ecdf_warming.png'), dpi = 300, width = 8, height = 10, units = 'in')

p1 <- ggplot(data = test_diff_warm, aes(x = anomaly_days, color = warming))
p1 + stat_ecdf() +
  theme_classic() +
  scale_color_viridis_d(begin = 0.3, end = 0.8)+
  theme(
    legend.title = element_blank()
  )+
  labs(x = 'Transition Period Anomaly (days)',
       y = 'ECDF (%)')+
  facet_wrap(~ice_type + event, ncol = 2)

#ggsave(here('results/ecdf/ecdf_ice_type.png'), dpi = 300, width = 8, height = 10, units = 'in')

# 9. Std. Dev. of transition periods --------------------------------------

#Calculate the std. dev. for the transition period times, not the anomalies
#as recommended by the reviewer to identify changes in variability of the
#transition periods.

rw_sd <- tibble(
  warming = c(
    as.factor(c('0\u00B0C', '1\u00B0C', '2\u00B0C', '4\u00B0C')) #This code (\u00B0) adds a degree symbol
  ),
  ice_on_sd_blk = c(
    mean(rw_ice_on_0c, na.rm = T), 
    mean(rw_ice_on_1c, na.rm = T), 
    mean(rw_ice_on_2c, na.rm = T), 
    mean(rw_ice_on_4c, na.rm = T)
  ),
  ice_off_sd_blk = c(
    mean(rw_ice_off_0c, na.rm = T), 
    mean(rw_ice_off_1c, na.rm = T), 
    mean(rw_ice_off_2c, na.rm = T), 
    mean(rw_ice_off_4c, na.rm = T)
  ),
  ice_on_sd_50 = c(
    mean(rw_ice_on_0.5wht_0c, na.rm = T), 
    mean(rw_ice_on_0.5wht_1c, na.rm = T), 
    mean(rw_ice_on_0.5wht_2c, na.rm = T), 
    mean(rw_ice_on_0.5wht_4c, na.rm = T)
  ),
  ice_off_sd_50 = c(
    mean(rw_ice_off_0.5wht_0c, na.rm = T), 
    mean(rw_ice_off_0.5wht_1c, na.rm = T), 
    mean(rw_ice_off_0.5wht_2c, na.rm = T), 
    mean(rw_ice_off_0.5wht_4c, na.rm = T)
  ),
  ice_on_sd_wht = c(
    mean(rw_ice_on_wht_0c, na.rm = T), 
    mean(rw_ice_on_wht_1c, na.rm = T), 
    mean(rw_ice_on_wht_2c, na.rm = T), 
    mean(rw_ice_on_wht_4c, na.rm = T)
  ),
  ice_off_sd_wht = c(
    mean(rw_ice_off_wht_0c, na.rm = T), 
    mean(rw_ice_off_wht_1c, na.rm = T), 
    mean(rw_ice_off_wht_2c, na.rm = T), 
    mean(rw_ice_off_wht_4c, na.rm = T)
  )
)
#View(rw_sd)