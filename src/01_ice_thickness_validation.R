
# 1. Load libraries -------------------------------------------------------

library(ncdf4)
library(tidyverse)
library(here)
library(janitor)
library(fields)
library(Metrics)
library(ggpmisc)

# 3. Import lake ice thickness files --------------------------------------

# 2. Import data ----------------------------------------------------------

#I have commented out the individual lake observations, because the individual 
#data are not provided in the figshare account.
#This is because the data derive from other sources, not all of which are openly
#available. Some authors have not published their data and have restricted access
#To be available upon request. 
#To access those data to complete the validation section, please go to the individual
#datasets.

#Below are the sources for the lake ice observations, also cited in the manuscript.

#Finnish ice data:  	
#Finnish Environment Institute (Syke). 
#Ice thickness observation network. 
#Open Environmental Information Systems: 
#Hertta; 2024. Available: https://wwwp2.ymparisto.fi/scripts/kirjaudu.asp

#North Temperate Lakes:
#Magnuson JJ, Carpenter SR, Stanley EH. 
#North Temperate Lakes LTER: Snow and Ice Depth 1982 - current. 
#[object Object]; 2023. doi:10.6073/PASTA/038D117F5244B520F7F87EF41DBA7A31

#MacDonald and Clear (Haliburton area):
#Ariano SS, Brown LC. Ice processes on medium‐sized north‐temperate lakes.
#Hydrological Processes. 2019;33:2434–2448. https://doi.org/10.1002/ hyp.13481
#Data available upon request from the authors.

#Lake Vendyurskoe:
#Zdorovennova, G.; Palshin, N.; Golosov, S.; Efremova, T.; Belashev, B.; 
#Bogdanov, S.; Fedorova, I.; Zverev, I.; Zdorovennov, R.; Terzhevik, A. 
#Dissolved Oxygen in a Shallow Ice-Covered Lake in Winter: 
#Effect of Changes in Light, Thermal and Ice Regimes. 
#Water 2021, 13, 2435. https://doi.org/10.3390/w13172435
#Data available upon request from authors

#Canadian lakes:
#Canadian Ice Service. Ice Thickness Program Collection, 1947–2002 
#(Environment and Climate Change Canada, 2015). 
#available at https://www.canada.ca/en/environment-climate-change/services/ice-forecasts-observations/latest-conditions/archive-overview/thickness-data.html

#Additional note: The data path (e.g., 'data/validation_data') is relative.
#The 'here()' function chooses your working directory. You can save the data from
#Figshare wherever you like and simply update the data path within the 'here()' function. 

# fin_total_ice <- read_csv(here('data/validation_data/finnish_total_ice.csv'))
# 
# fin_geo_loc <- read_csv(here('data/validation_data/finnish_standard.csv'))
# 
# ntl_ice <- read_csv(here('data/validation_data/ntl_ice_quality.csv'))
# 
# hal <- read_csv(here('data/validation_data/haliburton_field_data_2016_17.csv'))
# 
# rus <- read_csv(here('data/validation_data/ice_quality_zdorovennova_2021.csv'))
# 
# can1 <- read_csv(here('data/validation_data/citp_to_2002.csv'),
#                  locale = locale(encoding = 'latin1'))
# 
# can2 <- read_csv(here('data/validation_data/citp_2002_2024_ver2.csv'),
#                  locale = locale(encoding = 'latin1'))
# 
# can_meta <- read_csv(here('data/validation_data/citp_metadata.csv')) %>% 
#   select(
#     site = station_ID,
#     lat = latitude,
#     lon = longitude
#   ) %>% 
#   mutate(
#     site = as.factor(site)
#   )

# 3. Clean data -----------------------------------------------------------

#**Clean Finnish data----

#Clean total ice data
fin_tot_ice <- fin_total_ice %>% 
  separate(lake_id, c('site', 'I')) %>%  #removing leading 0 and trailing "_I" to aggregate with tibble "finnish_standard.csv" obtained from Aman Basu
  mutate(
    site = as.numeric(site),
    site = as.factor(site)
  ) %>% 
  select(
    -I
  ) 

#Retain only geo location data for Finnish lakes
fin_std <- fin_geo_loc %>% 
  arrange(site) %>% 
  select(site, lat, lon) %>% 
  mutate(
    site = as.factor(site)
  ) %>% 
  group_by(site) %>% 
  slice(1) 

#aggregate Finnish data

fin_clean <- fin_tot_ice %>% 
  full_join(fin_std) %>% 
  select(
    site,
    lat,
    lon,
    date,
    year,
    total_thickness
  ) %>% 
  mutate(
    year = if_else(month(date)>=10, year, year-1)
  )


#**Clean NTL data----

ntl_clean <- ntl_ice %>% 
  dplyr::select(
    lakeid,
    year4, daynum, sampledate,
    lat, lon,
    totice
  ) %>% 
  mutate(
    site = as.factor(lakeid),
    date = mdy(sampledate),
    year = if_else(month(date)>=10, year4, year4-1), #water year, but calling year to fit with Finnish dataframe
    total_thickness = totice
  ) %>% 
  select(
    site,
    lat, lon,
    date, year,
    total_thickness
  )

#**Clean Haliburton----

hal_clean <- hal %>%
  clean_names() %>% 
  dplyr::select(
    1:3, 7, 10, 11
  ) %>% 
  # mutate(
  #   site = if_else(site_id == 'm1' | site_id == 'm4' | site_id == 'm5' | site_id == 'm6', 'm', 'c')
  # ) %>% 
  mutate(
    site = as.factor(site_id),
    site = if_else(site_id == 'm1' | site_id == 'm4' | site_id == 'm5' | site_id == 'm6', 'm', 'c'),
    lon = lon*-1,
    total_thickness = total,
    date = mdy(date)
  ) %>% 
  dplyr::select(
    site,
    lat, lon,
    date, year,
    total_thickness
  )

#**Clean Lake Vendyurskoe data----

rus_clean <- rus %>% 
  select(
    1, 4, 8, 9
  ) %>% 
  mutate(
    site = as.factor('vendyurskoe'),
    date = mdy(date),
    year = year(date),
    total_thickness = total_ice_avg_cm
  ) %>% 
  select(
    site, 
    lat, lon,
    date, year,
    total_thickness
  )

#**Clean Canadian Ice Thickness Project lakes----

can1_clean <- can1 %>% 
  select(
    site = 1,
    date = 3, 
    total_thickness = 4
  ) %>% 
  mutate(
    site = as.factor(site),
    date = mdy(date),
    year = if_else(month(date)>=10, year(date), year(date)-1)
  ) %>% 
  inner_join(
    can_meta
  ) %>% 
  select(
    site,
    lat, lon,
    date, year,
    total_thickness
  )

can2_clean <- can2 %>% 
  select(
    site = 1,
    date = 3,
    total_thickness = 4
  ) %>% 
  mutate(
    site = as.factor(site),
    year = if_else(month(date)>=10, year(date), year(date)-1)
  ) %>% 
  inner_join(
    can_meta
  ) %>% 
  select(
    site,
    lat, lon,
    date, year,
    total_thickness
  )

# 4. Combine the dataframes -----------------------------------------------

ice_thickness_clean <- can1_clean %>% 
  bind_rows(
    can2_clean,
    fin_clean,
    hal_clean,
    ntl_clean,
    rus_clean
  ) %>% 
  filter(
    year>=1980
  ) %>% 
  mutate(
    site = as.factor(site)
  ) %>% 
  arrange(site, date)

loc_data <- ice_thickness_clean %>% 
  select(
    site, lat, lon
  ) %>%
  group_by(site) %>% 
  slice(1)

#You can write the 
#write_csv(ice_thickness_clean, here('data/validation_data/ice_thickness_combined.csv'))
#write_csv(loc_data, here('data/validation_data/loc_data.csv'))


# 5. CESM2-LE Validataion files -------------------------------------------

val_df <- nc_open(here('data/validation_data/model_lake_ice_thick_validation.nc'))

li_thick <- ncvar_get(val_df, 'LAKEICETHICK')
time <- ncvar_get(val_df, 'time') #in days since 1980-01-01 without leap years
lon <- ncvar_get(val_df, 'lon')
lat <- ncvar_get(val_df, 'lat')
lake_number <- ncvar_get(val_df, 'lake_number')

#dim(li_thick)

mean <- apply(li_thick, c(1,2), mean, na.rm = TRUE)
#dim(mean)

median <- apply(li_thick, c(1,2), median, na.rm = TRUE)
#dim(median)

max <- apply(li_thick, c(1,2), max, na.rm = TRUE)
#dim(max)

min <- apply(li_thick, c(1,2), min, na.rm = TRUE)
#dim(min)

sd <- apply(li_thick, c(1,2), sd, na.rm = TRUE)
#dim(sd)

mean_long <- reshape2::melt(mean, id = 1) %>% 
  arrange(Var1) %>% 
  select(
    lake = Var1,
    day = Var2,
    ice_mean_m = value
  ) %>% 
  mutate(
    lake = as.factor(lake),
    day = as.factor(day)
  )

median_long <- reshape2::melt(median, id = 1) %>% 
  arrange(Var1) %>% 
  select(
    lake = Var1,
    day = Var2,
    ice_median_m = value
  ) %>% 
  mutate(
    lake = as.factor(lake),
    day = as.factor(day)
  )

max_long <- reshape2::melt(max, id = 1) %>% 
  arrange(Var1) %>% 
  select(
    lake = Var1,
    day = Var2,
    ice_max_m = value
  ) %>% 
  mutate(
    lake = as.factor(lake),
    day = as.factor(day)
  )

min_long <- reshape2::melt(min, id = 1) %>% 
  arrange(Var1) %>% 
  select(
    lake = Var1,
    day = Var2,
    ice_min_m = value
  ) %>% 
  mutate(
    lake = as.factor(lake),
    day = as.factor(day)
  )

sd_long <- reshape2::melt(sd, id = 1) %>% 
  arrange(Var1) %>% 
  select(
    lake = Var1,
    day = Var2, 
    ice_sd_m = value
  ) %>% 
  mutate(
    lake = as.factor(lake),
    day = as.factor(day)
  )

#Need to create a date sequence for the data. Code from: https://stackoverflow.com/questions/55576747/remove-leap-day-from-date-sequence
myDates=data.frame(seq(as.Date("1980-01-01"), to=as.Date("2024-12-30"),by="days"))
names(myDates)= "Dates"
nrow(myDates) #16436 -- This value means that leap days have been included, but I need to exclude them.

#Identify and remove leap days
myDates <- myDates[!(format(myDates$Dates,"%m") == "02" & format(myDates$Dates, "%d") == "29"), ,drop = FALSE]
nrow(myDates) #16424 -- Date sequence with leap days removed.

#Now repeat the dates 72 times (once for each lake) to make the dfs the correct sizes
myDates2 <- as.data.frame(rep(c(myDates), times = 72)) %>% 
  mutate(
    values = 1:16424
  ) %>% 
  pivot_longer(
    !values,
    names_to = 'rep_names',
    values_to = 'date'
  ) %>% 
  arrange(
    rep_names
  ) %>% 
  select(
    date
  )

cesm_model_val <- mean_long %>% 
  bind_cols(median_long, max_long, min_long, sd_long, myDates2) %>% 
  select(
    lake = 1,
    16,3,6,9,12,15
  ) 

#write_csv(cesm_model_val, here('data/validation_data/cesm_val_data.csv'))

# Comparison between in situ and modeled data -----------------------------

#Again, file paths are relative. Simply correct the below read_csv() functions
#to your desired file path.

obs <- read_csv(here('data/validation_data/ice_thickness_combined.csv'))

mod <- read_csv(here('data/validation_data/cesm_val_data.csv'))

#####

obs1 <- obs %>% 
  mutate(
    site = if_else(
      site == 'c', as.factor('Hal_C'), site
    ),
    site = if_else(
      site == 'm', as.factor('Hal_M'), site
    ),
    site = if_else(
      site == 'vendyurskoe', as.factor('Vendyurskoe'), site
    )
  ) %>% 
  filter(
    site != '67111'
  )

#unique(obs1$site)

obs2 <- obs1 %>% 
  group_by(site) %>% 
  mutate(
    lake = cur_group_id()
  )

# unique(obs2$lake)
# length(unique(obs2$lake))
# unique(mod$lake)
# length(unique(mod$lake))

#Combine observed and modelled data

val_data <- obs2 %>% 
  inner_join(mod) %>% 
  #need to change mean, median, and sd to cm from m
  mutate(
    ice_mean_cm = ice_mean_m*100,
    ice_median_cm = ice_median_m*100,
    ice_max_cm = ice_max_m*100,
    ice_min_cm = ice_min_m*100,
    ice_sd_cm = ice_sd_m*100
  ) %>% 
  na.omit()

rmse(val_data$total_thickness, val_data$ice_mean_cm) #28.5 cm
rmse(val_data$total_thickness, val_data$ice_median_cm) #29.3 cm

val_daily_rmse <- val_data %>% 
  group_by(lake) %>% 
  summarise(
    rmse_mean = rmse(total_thickness, ice_mean_cm),
    rmse_med = rmse(total_thickness, ice_median_cm),
    obs = n()
  )

max(val_daily_rmse$rmse_mean) #80.6 cm
min(val_daily_rmse$rmse_mean) #9.5 cm
max(val_daily_rmse$rmse_med) #81.1 cm
min(val_daily_rmse$rmse_med) #8.6 cm
sum(val_daily_rmse$obs) #28,772

val_plt <- ggplot(data = val_data, aes(x = total_thickness, y = ice_mean_cm))+
  geom_point(, alpha = 0.2)+
  theme_classic()+
  geom_smooth( se = FALSE, method = 'lm', color = 'grey50', alpha = 0.5)+
  stat_poly_eq(label.y = 0.98, label.x = 0.05)+
  stat_poly_eq(use_label('eq'), label.y = 0.93, label.x = 0.05)+
  xlab('Observed (cm)')+
  ylab('Modeled (cm)')+
  theme(
    text = element_text(size = 25)
  )
#val_plt


#Look at data by year
val_yearly <- val_data %>% 
  group_by(year, lake) %>% 
  summarise(
    mean_obs = mean(total_thickness),
    med_obs = median(total_thickness),
    mean_mod = mean(ice_mean_cm),
    med_mod = median(ice_median_cm),
    sd = mean(ice_sd_cm)
  ) %>% 
  arrange(lake)

rmse(val_yearly$mean_obs, val_yearly$mean_mod) #20.0 cm
rmse(val_yearly$med_obs, val_yearly$med_mod) #21.4 cm

val_yearly_rmse <- val_yearly %>% 
  group_by(lake) %>% 
  summarise(
    rmse_mean = rmse(mean_obs, mean_mod),
    rmse_med = rmse(med_obs, med_mod),
    obs = n()
  )

max(val_yearly_rmse$rmse_mean) #68.7 cm
min(val_yearly_rmse$rmse_mean) #7.5 cm
max(val_yearly_rmse$rmse_med) #76.9 cm
min(val_yearly_rmse$rmse_med) #6.5 cm
sum(val_yearly_rmse$obs) #2,296 observations

val_yearly_plt <- ggplot(data = val_yearly, aes(x = mean_obs, y = mean_mod))+
  geom_errorbar(aes(ymin = mean_mod-sd, ymax = mean_mod+sd), width = 0.2, alpha = 0.2)+
  geom_point(alpha = 0.2)+
  theme_classic()+
  geom_smooth( se = FALSE, method = 'lm', color = 'grey50', alpha = 0.5)+
  stat_poly_eq(label.y = 0.98, label.x = 0.05)+
  stat_poly_eq(use_label('eq'), label.y = 0.93, label.x = 0.05)+
  xlab('Observed (cm)')+
  ylab('Modeled (cm)')+
  theme(
    text = element_text(size = 25)
  )
#val_yearly_plt

#Look at data by lake
val_all <- val_data %>% 
  group_by(lake) %>% 
  summarise(
    mean_obs = mean(total_thickness),
    med_obs = median(total_thickness),
    mean_mod = mean(ice_mean_cm),
    med_mod = median(ice_median_cm),
    max_obs = max(total_thickness),
    max_mod = mean(ice_max_cm),
    min_obs = min(total_thickness),
    min_mod = mean(ice_min_cm),
    sd = mean(ice_sd_cm),
    sd2 = mean(ice_sd_cm)*2
  ) %>% 
  arrange(lake)

rmse(val_all$mean_obs, val_all$mean_mod) #18.5 cm
rmse(val_all$med_obs, val_all$med_mod) #18.5 cm
rmse(val_all$max_obs, val_all$max_mod) #39.1 cm
rmse(val_all$min_obs, val_all$min_mod) #11.0 cm

val_lake_rmse <- val_all %>% 
  group_by(lake) %>% 
  summarise(
    rmse_mean = rmse(mean_obs, mean_mod),
    rmse_med = rmse(med_obs, med_mod),
    rmse_max = rmse(max_obs, max_mod),
    obs = n()
  )

max(val_lake_rmse$rmse_mean) #69.8 cm
min(val_lake_rmse$rmse_mean) #3.6 cm
max(val_lake_rmse$rmse_med) #76.3 cm
min(val_lake_rmse$rmse_med) #0 cm
sum(val_lake_rmse$obs) #71 observations

val_all_plt <- ggplot(data = val_all, aes(x = mean_obs, y = mean_mod))+
#val_all_plt <- ggplot(data = val_all, aes(x = max_obs, y = max_mod))+
  geom_errorbar(aes(ymin = mean_mod-sd2, ymax = mean_mod+sd2), width = 0.2, alpha = 0.2)+
  #geom_errorbar(aes(ymin = min_mod, ymax = max_mod), width = 0.2, alpha = 0.2)+
  geom_point(alpha = 0.6)+
  theme_classic()+
  geom_smooth( se = FALSE, method = 'lm', color = 'grey50', alpha = 0.5)+
  stat_poly_eq(label.y = 0.98, label.x = 0.05)+
  stat_poly_eq(use_label('eq'), label.y = 0.93, label.x = 0.05)+
  xlab('Observed (cm)')+
  ylab('Modeled (cm)')+
  theme(
    text = element_text(size = 20)
  )+
  xlim(c(0,155))
val_all_plt

#ggsave(here('results/validation_plts/val_lake_max_plt.png'), dpi = 300, width = 6, height = 4, units = 'in')

# Example lake timeseries section -----------------------------------------

#First filter data to just show lake 7
lake7_obs <- obs2 %>% 
  filter(lake == '7')

lake7_mod <- mod %>% 
  filter(lake == '7') %>% 
  mutate(
    mean_cm = ice_mean_m*100,
    max_cm = ice_max_m*100,
    min_cm = ice_min_m*100,
    sd_cm = ice_sd_m*100,
    sd_min = mean_cm-sd_cm,
    sd_max = mean_cm+sd_cm,
    sd_min = if_else(sd_min<0, 0, sd_min),
    sd2_min = mean_cm-(sd_cm*2),
    sd2_max = mean_cm+(sd_cm*2),
    sd2_min = if_else(sd2_min<0, 0, sd2_min)
  )

#Now plot the mod as a continuous line and obs as points on top of the modeled ts
val_ts_plt <- ggplot()+
  geom_ribbon(
    data = lake7_mod,
    aes(x = date, ymin = sd2_min, ymax = sd2_max),
    alpha = 0.5, fill = '#67B7D1'
  )+
  geom_ribbon(
    data = lake7_mod,
    aes(x = date, ymin = sd_min, ymax = sd_max),
    alpha = 0.5, fill = 'grey50'
  )+
  geom_line(data = lake7_mod, aes(x = date, y = mean_cm))+
  geom_point(data = lake7_obs, aes(x = date, y = total_thickness), alpha = 0.6)+
  scale_x_date(limits = as.Date(c("2009-10-01","2015-12-30")))+
  theme_classic()+
  ylim(c(0,65))+
  xlab('')+
  ylab('Ice Thickness (cm)')+
  theme(
    text = element_text(size = 20)
  )
  
#val_ts_plt  

#ggsave(here('results/validation_plts/validation_time_series_minmax.png'), dpi = 300, width = 6, height = 4, units = 'in')
