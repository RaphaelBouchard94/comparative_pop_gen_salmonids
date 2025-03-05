rm(list=ls())

setwd(dir = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/")

#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)

#---------------#
#----- Data ----#
#---------------#

weather <- read.csv("data/env/all_stations_daily_weather.csv")

colnames(weather)

weather %>% 
  group_by(station_name) %>% 
  summarise(avg_el = mean(elev))

weather %>% filter(!year == 2024 & year > 1960) %>% 
  group_by(station_name, year) %>% 
  filter(mean_temp > 5.6) %>% 
  summarise(n_days = n()) %>% View()
  ungroup %>% 
  group_by(station_name) %>% 
  summarise(avg_dd = mean(n_days))
