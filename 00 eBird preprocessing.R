## ---------------------------
##
## Script name: eBird pre-processing
##
## Purpose of script: extract data for the citizen science schemes used,
##     currently eBird and something else. Format data to CSVs ready for use.
##
## Author: ****
##
## Date Created: 05/06/2025
##
## Email: ****
##
## ---------------------------
##
## Notes: eBird from raw file to dataset processed in line with analytical
##        guidelines:
## - Complete checklists, 
## - Stationary and travelling protocols only, 
## - Starting times between 5 a.m. and 8 p.m., 
## - Durations of continuous observation of 0–3 hr, 
## - Observers travelling distances during the checklist of up to 5 miles (8 km).
##
## ---------------------------

## set working directory for Mac and PC

#setwd("~/Google Drive/")  		# working directory (mac)
setwd("G:/")  	# working directory (windows)

## ---------------------------

## load up the packages we will need:
library(tidyverse)
library(devtools)
library(auk)
library(sf)

# resolve namespace conflicts
select <- dplyr::select


## ---------------------------

##- eBird processing ####

# ebird file location
data_file <- 'ebd_US_relMar-2025.txt'

# ebird sampling file
data_file_s <- 'ebd_US_relMar-2025_sampling.txt'

# where I want processed data to end up
data_path <- 'R Coding/MSID/data/' 

# data string (used to call same set of files throughout all scripts)
data_string <- 'bcr23_2025'


# big data process
ebd_filtered <- auk_ebd(file = data_file, file_sampling = data_file_s) %>% # load data
  auk_bcr(23) %>% # spatial filter to BCR 23
  auk_complete() %>% # filter complete checklists
  auk_distance(distance = c(0, 8)) %>% # effort distance of max 8km
  auk_duration(c(0, 120)) %>% # max duration of 3 hours
  auk_protocol(c('Stationary', 'Traveling')) %>% # only stationary or traveling protocols
  auk_time(c('05:00', '20:00')) %>% # checklist started between 5am and 8pm
  # auk_date(date = c("2016-01-01", "2024-12-31")) %>% # filter to years
  # auk_date(date = c('*-03-01', '*-08-31')) %>% # filter to breeding season
  auk_filter(file = paste0(data_path, data_string, '.txt'),
             file_sampling = paste0(data_path, data_string, '_s.txt'),
             overwrite = TRUE) # write to output files

# reading in an ebird file takes longer so I'll read it once and save the dataframe to a csv in future
ebird_og <- read_ebd(paste0(data_path, data_string, '.txt'), unique = FALSE) %>% 
  write.csv(paste0(data_path, data_string, '.csv'))





