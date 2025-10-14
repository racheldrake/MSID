## ---------------------------
##
## Script name: Macaulay library processing
##
## Purpose of script: use extracted numbers of photos and audio recordings from
##            the macaulay online library to construct an audio index per species
##
## Author: ****
##
## Date Created: 08/05/2025
##
## Email: ****
##
## ---------------------------
##
## Notes: 
##
##  
##
##
##
##
##
##
## ---------------------------

## set working directory for Mac and PC

#setwd("~/Google Drive/")  		# working directory (mac)
setwd("G:/R Coding/MSID")  	# working directory (windows)

## ---------------------------

## load up the packages we will need:

library(tidyverse)
library(terra)
library(sf)

# resolve namespace conflicts
select <- dplyr::select

## ---------------------------

# data tag
data_string <- "bcr23_MSID"

# datasets path
data_path <- 'data/'

# where I want to save results
results_path <- 'proc_data/'

# load in study region
bcr <- st_read(paste0(data_path, 'modis/BCR_Terrestrial_master.shp')) %>% 
  filter(BCR == 23) %>% st_union() %>% st_transform(crs = 4326)

bcr_bbox <- st_bbox(bcr)

# load in the raw data
data <- read_csv(paste0(data_path, 'macaulay_library_data.csv'))

universal_AVI <- data %>% 
  group_by(species_code, asset_format_code) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = asset_format_code, values_from = n, values_fill = 0)

write_csv(universal_AVI, paste0(results_path, 'universal_AVI.csv'))

BCR_23_AVI <- data %>% 
  drop_na(latitude, longitude) %>%
  # filter breeding season
  filter(obs_month %in% 3:8) %>%
  # filter lat long by bbox first for processing speed
  filter(longitude > bcr_bbox[1],
         latitude > bcr_bbox[2],
         longitude < bcr_bbox[3],
         latitude < bcr_bbox[4]) %>%
  st_as_sf(coords = c(longitude, latitude), crs = 4326) %>%
  st_intersection(bcr) %>% st_drop_geometry() %>%
  # calculate statistics
  group_by(species_code, asset_format_code) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = asset_format_code, values_from = n, values_fill = 0)


write_csv(BCR_23_AVI, paste0(results_path, 'BCR23_AVI.csv'))

universal_AVI <- read_csv(paste0(results_path, 'universal_AVI.csv'))
BCR_23_AVI <- read_csv(paste0(results_path, 'BCR23_AVI.csv'))

BCR_23_AVI %>% rename(audio_23 = audio, photo_23 = photo, video_23 = video) %>%
  left_join(universal_AVI, by = 'species_code') %>%
  left_join(taxon, by = 'species_code') %>%
  filter(common_name %in% model_data$common_name) %>%
  mutate(AVI_23 = audio_23/(audio_23 + photo_23),
         AVI = audio/(audio+photo)) %>% 
  ggplot() + geom_point(aes(AVI, AVI_23)) + theme_bw() +
  xlab('Universal Audio Visual Index') + 
  ylab('BCR23 Audio Visual Index')

BCR_23_AVI %>% 
  left_join(taxon, by = 'species_code') %>%
  # species with no photo records are always lost/extinct species not
  # relevant to the metric
  #calculate index
  mutate(audio_index = (audio)/(audio + photo)) %>%
  # remove unneccessary columns
  select(-c(photo, video, audio, species_code)) %>%
  # write to csv
  write_csv(paste0('analysis_data/audio_index.csv'))
