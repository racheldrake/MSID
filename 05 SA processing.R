## ---------------------------
##
## Script name: Species availability processing
##
## Purpose of script: Calculate each species occurrence per observer across
##                    all checklists
##
## Author: ****
##
## Date Created: 10/07/2025
##
## Email: ****
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac and PC

#setwd("~/Google Drive/")  		# working directory (mac)
setwd("G:/R Coding/MSID")  	# working directory (windows)

## ---------------------------

## load up the packages we will need:

library(tidyverse)
library(ebirdst)
library(sf)
library(terra)
library(exactextractr)
#set_ebirdst_access_key("xxx")

# resolve namespace conflicts
select <- dplyr::select

## ---------------------------

# data tag
data_string <- "bcr23_MSID"

# eBird path
ebd_path <- 'proc_data/'

# where I want to save results
results_path <- 'proc_data/'

# load in all observers for 2023-24
ebd <- read_csv(paste0(ebd_path, data_string, '_2023_24_all.csv'))

# load in status CRS
crs <- readRDS('status_crs.rds')

# find all unique locations in ebd
ebd_pts <- ebd %>%
  distinct(locality_id, latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% st_transform(crs = crs)

# filter to give necessary species
species <- unique(ebd$common_name)
status_species <- ebirdst_runs %>% filter(!is.na(breeding_quality))
species <- species[species %in% status_species$common_name]

# some species that don't process
species <- species[species != 'Whooper Swan']
species <- species[species != 'Whooping Crane']
species <- species[species != 'Great Black-backed Gull']
species <- species[species != 'Northern Saw-whet Owl']
species <- species[species != 'Red Crossbill']
species <- species[species != 'Fish Crow']
species <- species[species != 'African Collared-Dove']
species <- species[species != 'Budgerigar']
# yebsap needs code because it has the example files I believe
species <- species[species != 'Yellow-bellied Sapsucker']
species <- append(species, 'yebsap')

name_code <- data.frame()

# find speceis occurrence for defined species and locations
for (target_species in species){
  # minimise amount of data to download
  ebirdst_download_status(target_species, 
                          path = 'data/', 
                          pattern = 'occurrence_seasonal_mean_3km', 
                          download_abundance = FALSE, 
                          download_occurrence = TRUE)
  
  # load status as raster
  status <- load_raster(target_species, product = "occurrence",
                        period = "seasonal", metric = "mean",
                        resolution = "3km", path = 'data/')
  status <- status[['breeding']]
  
  # extract occurrence values
  occurrence_vals <- terra::extract(status, ebd_pts, method = 'bilinear')[,2]
  
  # name the new column using the value in target_species
  species_colname <- make.names(target_species)  # ensures it's a valid column name
  ebd_pts[[species_colname]] <- occurrence_vals
  
  name_code <- rbind(name_code, c(target_species, species_colname))
  
}

# rename title row (doesn't affect data it's just for cleaning)
name_code <- name_code %>% rename(names = X.Red.tailed.Hawk..1, common_name = X.Red.tailed.Hawk.)

# join with locations
ebd_join <- ebd_pts %>% st_drop_geometry() %>%
  pivot_longer(cols = -1, names_to = 'names', values_to = 'occur') %>%
  left_join(name_code, by = 'names') %>% select(-names) %>% write_csv(paste0(results_path, 'species_occurrence_breeding.csv'))

# repeat for the pre-breeding season
name_code <- data.frame()

for (target_species in species){
  # ebirdst_download_status(target_species, path = 'data/', pattern = 'occurrence_seasonal_mean_3km', download_abundance = FALSE, download_occurrence = TRUE)
  status <- load_raster(target_species, product = "occurrence",
                        period = "seasonal", metric = "mean",
                        resolution = "3km", path = 'data/')
  status <- status[['prebreeding_migration']]
  
  # extract occurrence values
  occurrence_vals <- terra::extract(status, ebd_pts, method = 'bilinear')[,2]
  
  # name the new column using the value in target_species
  species_colname <- make.names(target_species)  # ensures it's a valid column name
  ebd_pts[[species_colname]] <- occurrence_vals
  
  name_code <- rbind(name_code, c(target_species, species_colname))
  
}

name_code <- name_code %>% rename(names = X.Red.tailed.Hawk..1, common_name = X.Red.tailed.Hawk.)

ebd_join <- ebd_pts %>% st_drop_geometry() %>%
  pivot_longer(cols = -1, names_to = 'names', values_to = 'occur') %>%
  left_join(name_code, by = 'names') %>% select(-names) %>% write_csv(paste0(results_path, 'species_occurrence_prebreeding.csv'))


# call in both datasets
ebd_join_breeding <- read_csv(paste0(results_path, 'species_occurrence_breeding.csv'))
ebd_join_prebreeding <- read_csv(paste0(results_path, 'species_occurrence_prebreeding.csv'))

# join and average
ebd_join <- ebd_join_breeding %>% rename(breeding_occur = occur) %>%
  left_join(ebd_join_prebreeding, by = c('locality_id', 'common_name')) %>%
  mutate(occur = (occur + breeding_occur)/2) %>% select(-breeding_occur) %>% drop_na()

# then aggregate to observer level
# note this takes ages to run to calculate realised occurrence don't do it 
# unless you absolutely must
ebd_checklist <- ebd %>% select(observer_id, checklist_id, locality_id) %>%
  distinct(observer_id, checklist_id, locality_id) %>%
  crossing(common_name = species) %>% 
  left_join(ebd_join, by = c('locality_id', 'common_name')) %>%
  select(-checklist_id, -locality_id) %>%
  group_by(observer_id, common_name) %>%
  summarise(realised_occur = 1 - prod(1 - occur),
            occur = mean(occur)) %>% 
  write_csv(paste0(results_path, 'species_occurrence.csv'))
 