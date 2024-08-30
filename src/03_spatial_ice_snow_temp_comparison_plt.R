
# 1. Load libraries -------------------------------------------------------

library(ncdf4)
library(tidyverse)
library(here)
library(janitor)
library(fields)

ice_ann_mean <- nc_open(
  here(
    'data/cesm2_le_netcdf_files/iceANNmean_mod_ensmean.nc'
  )
)

ice_ann_max <- nc_open(
  here(
    'data/cesm2_le_netcdf_files/iceANNmax_mod_ensmean.nc'
  )
)

snow_ann_mean <- nc_open(
  here(
    'data/cesm2_le_netcdf_files/snowdpANNmean_mod_ensmean.nc'
  )
)

snow_ann_max <- nc_open(
  here(
    'data/cesm2_le_netcdf_files/snowdpANNmax_mod_ensmean.nc'
  )
)

tsa_ann_1c <- nc_open(
  here(
    'data/cesm2_le_netcdf_files/tsa_ann_ano_1C.nc'
  )
)

tsa_ann_2c <- nc_open(
  here(
    'data/cesm2_le_netcdf_files/tsa_ann_ano_2C.nc'
  )
)

tsa_ann_4c <- nc_open(
  here(
    'data/cesm2_le_netcdf_files/tsa_ann_ano_4C.nc'
  )
)

# 3. Extract lat, lon, year -----------------------------------------------

#Extract lat, long, and year data 
#NOTE: Values will be the same for all, so extraction is from ice duration only
#NOTE: Except tsa_ann_Xc, b/c the 'time' dimension for the other variables
#      is actually an ensemble number.
#      The dimensions for tsa_ann_Xc are 288x192x90 for lon x lat x ensemble member

lat <- ncvar_get(ice_ann_max,'lat')
lon <- ncvar_get(ice_ann_max,'lon')
year <- ncvar_get(ice_ann_max,'year')
ensemble <- ncvar_get(tsa_ann_1c, )


# 4. Extract vars of interest ---------------------------------------------

li_mean <- ncvar_get(ice_ann_mean, 'LAKEICETHICK')
li_max <- ncvar_get(ice_ann_max, 'LAKEICETHICK')

snw_mean <- ncvar_get(snow_ann_mean, 'SNOWDP')
snw_max <- ncvar_get(snow_ann_max, 'SNOWDP')

t_1c <- ncvar_get(tsa_ann_1c, 'TREFHT')
t_2c <- ncvar_get(tsa_ann_2c, 'TREFHT')
t_4c <- ncvar_get(tsa_ann_4c, 'TREFHT')

# **4b. Apply a lake mask -------------------------------------------------

#Obtain grids where there are lakes

#Create a dataframe for lake pixels
lake_mask <- li_mean[,,1]

#When a value does not == NA, assign it 1 (lake present)
lake_mask[!is.na(lake_mask)] <- TRUE

#When a value == NA, assign it 0 (no lake present)
lake_mask[is.na(lake_mask)] <- FALSE

#Apply lake mask to mean snow data
snw_mean_lm <- snw_mean

snw_mean_lm[lake_mask == FALSE] <- NA

#Apply lake mask to max snow data
snw_max_lm <- snw_max

snw_max_lm[lake_mask == FALSE] <- NA

#Apply lake mask to temperature data
t_1c_lm <- t_1c
t_2c_lm <- t_2c
t_4c_lm <- t_4c

t_1c_lm[lake_mask == FALSE] <- NA
t_2c_lm[lake_mask == FALSE] <- NA
t_4c_lm[lake_mask == FALSE] <- NA


# 5. Confine lon/lat to Northern Hemisphere -------------------------------

lon_nh = lon >= 0 & lon <= 360
lat_nh = lat >= 40 & lat <= 90


# 6. Determine anomalies for ice and snow ---------------------------------

#Years against historical anomaly
hist_years <- year>=1851 & year <= 1880

# **6a. Lake ice historical condition----

li_hist <- li_mean[lon_nh, lat_nh, hist_years]
li_hist_max <- li_max[lon_nh, lat_nh, hist_years]
#dim(li_hist) #check dimensions to ensure correct matrix slices

#Get the mean of the historical period
li_hist_mean <- apply(li_hist, c(1,2), mean, na.rm = T)
li_max_hist_mean <- apply(li_hist_max, c(1,2), mean, na.rm = T) #This is the hist. condition of max ice thickness


# **6b. Snow historical condition----
snw_hist <- snw_mean_lm[lon_nh, lat_nh, hist_years]
snw_hist_max <- snw_max_lm[lon_nh, lat_nh, hist_years]

#Get the mean of the historical period for snow depth
snw_hist_mean <- apply(snw_hist, c(1,2), mean, na.rm = T)
snw_max_hist_mean <- apply(snw_hist_max, c(1,2), mean, na.rm = T)

#Repeat historical means for length of 'year' so that an anomaly can be obtained

li_mean_hist_rep <- replicate(length(year), li_hist_mean)
li_max_hist_rep <- replicate(length(year), li_max_hist_mean)
snw_mean_hist_rep <- replicate(length(year), snw_hist_mean)
snw_max_hist_rep <- replicate(length(year), snw_max_hist_mean)


# **6c. Determine anomalies----

li_mean_anom <- li_mean[lon_nh, lat_nh,] - li_mean_hist_rep
li_max_anom <- li_max[lon_nh, lat_nh,] - li_max_hist_rep

snw_mean_anom <- snw_mean_lm[lon_nh, lat_nh,] - snw_mean_hist_rep
snw_max_anom <- snw_max_lm[lon_nh, lat_nh,] - snw_max_hist_rep

# 7. Isolate warming scenarios for snow and ice ---------------------------


# **7a. Lake ice ----------------------------------------------------------

#Anomaly of mean ice thickness
li_anom_1c <- li_mean_anom[,,year == 2011] 
li_anom_2c <- li_mean_anom[,,year == 2044]
li_anom_4c <- li_mean_anom[,,year == 2088]

li_max_1c <- li_max_anom[,,year == 2011]
li_max_2c <- li_max_anom[,,year == 2044]
li_max_4c <- li_max_anom[,,year == 2088]


# **7b. Snow depth --------------------------------------------------------

snw_mean_1c <- snw_mean_anom[,,year == 2011]
snw_mean_2c <- snw_mean_anom[,,year == 2044]
snw_mean_4c <- snw_mean_anom[,,year == 2088]

snw_max_1c <- snw_max_anom[,,year == 2011]
snw_max_2c <- snw_max_anom[,,year == 2044]
snw_max_4c <- snw_max_anom[,,year == 2088]

# **7c. Temperature -------------------------------------------------------

t1_1c <- t_1c_lm[lon_nh, lat_nh, ]
t1_2c <- t_2c_lm[lon_nh, lat_nh, ]
t1_4c <- t_4c_lm[lon_nh, lat_nh, ]


# 8. Get mean for ensemble members ----------------------------------------

t_mean_1c <- apply(t1_1c, c(1,2), mean, na.rm = T)
t_mean_2c <- apply(t1_2c, c(1,2), mean, na.rm = T)
t_mean_4c <- apply(t1_4c, c(1,2), mean, na.rm = T)


# 9. Convert matrices to dataframes ---------------------------------------


# **9a. Create mean and max lake ice dataframes ---------------------------

li_1c_df <- reshape2::melt(li_anom_1c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    li_mean_1c = value
  ) %>% 
  select(
    -value
  )

li_2c_df <- reshape2::melt(li_anom_2c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    li_mean_2c = value
  ) %>% 
  select(
    -value
  )

li_4c_df <- reshape2::melt(li_anom_4c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    li_mean_4c = value
  ) %>% 
  select(
    -value
  )

li_max_1c_df <- reshape2::melt(li_max_1c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    li_max_1c = value
  ) %>% 
  select(
    -value
  )

li_max_2c_df <- reshape2::melt(li_max_2c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    li_max_2c = value
  ) %>% 
  select(
    -value
  )

li_max_4c_df <- reshape2::melt(li_max_4c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    li_max_4c = value
  ) %>% 
  select(
    -value
  )


# **9b. Snow depth dataframes ---------------------------------------------

snw_1c_df <- reshape2::melt(snw_mean_1c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    snw_mean_1c = value
  ) %>% 
  select(
    -value
  )

snw_2c_df <- reshape2::melt(snw_mean_2c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    snw_mean_2c = value
  ) %>% 
  select(
    -value
  )

snw_4c_df <- reshape2::melt(snw_mean_4c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    snw_mean_4c = value
  ) %>% 
  select(
    -value
  )

snw_max_1c_df <- reshape2::melt(snw_max_1c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    snw_max_1c = value
  ) %>% 
  select(
    -value
  )

snw_max_2c_df <- reshape2::melt(snw_max_2c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    snw_max_2c = value
  ) %>% 
  select(
    -value
  )

snw_max_4c_df <- reshape2::melt(snw_max_4c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    snw_max_4c = value
  ) %>% 
  select(
    -value
  )

# **9c. Temperature dataframes --------------------------------------------

temp_1c_df <- reshape2::melt(t_mean_1c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    temp_1c = value
  ) %>% 
  select(
    -value
  )
#dim(test_temp)

temp_2c_df <- reshape2::melt(t_mean_2c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    temp_2c = value
  ) %>% 
  select(
    -value
  )

temp_4c_df <- reshape2::melt(t_mean_4c, varnames = c('lon', 'lat')) %>% 
  mutate(
    lon = as.factor(lon),
    lat = as.factor(lat),
    temp_4c = value
  ) %>% 
  select(
    -value
  )


# **9d. Combine dataframes ------------------------------------------------

anomalies_combined <- li_1c_df %>%  #by = c('lon','lat')
  #mean lake ice dfs
  left_join(li_2c_df,) %>% left_join(li_4c_df) %>% 
  #max lake ice dfs
  left_join(li_max_1c_df) %>% left_join(li_max_2c_df) %>% left_join(li_max_4c_df) %>% 
  #mean snow dfs
  left_join(snw_1c_df) %>% left_join(snw_2c_df) %>% left_join(snw_4c_df) %>% 
  #max_snow dfs
  left_join(snw_max_1c_df) %>% left_join(snw_max_2c_df) %>% left_join(snw_max_4c_df) %>% 
  #temperature dfs
  left_join(temp_1c_df) %>% left_join(temp_2c_df) %>% left_join(temp_4c_df)


# 10. Plot anomalies ------------------------------------------------------


ggplot(data = anomalies_combined, aes(x = temp_4c, y = snw_mean_4c, color = li_mean_4c))+
  geom_point(alpha  = 0.6)+
  theme_classic()+
  xlab('Air Temperature Anomaly (\u00B0C)')+
  ylab('Snow Depth Anomaly (m)')+
  scale_color_continuous(name = "Ice Thickness\nAnomaly (m)")+
  guides(color=guide_colorbar(ticks.colour = NA))




