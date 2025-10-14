## ---------------------------
##
## Script name: eBird processing
##
## Purpose of script: extract subset of eBird data needed for analysis
##
## Author: ****
##
## Date Created: original 10/01/25, updated 01/04/2025, finalised 07/10/2025
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

# resolve namespace conflicts
select <- dplyr::select

## ---------------------------

# where to find eBird csv from last script
raw_data_path <- 'data/'

# where I want processed datasets to end up
data_path <- 'proc_data/' 

# name of data extraction
data_tag <- "bcr23_2025"

# ebird data
ebd <- read_csv(paste0(raw_data_path, data_tag, ".csv"))

# MSID data
MSID <- read_csv(paste0(raw_data_path, 'MSID.csv')) %>%
  rename(checklist_id = CHECKLIST_ID, MSID = WasMSIDrunning)

# tag to label resulting dataset
data_tag <- 'bcr23_MSID'

# filter years and months
ebd %>% 
  rename(checklist_id = sampling_event_identifier) %>%
  filter(year(ymd(observation_date)) %in% c(2023, 2024)) %>%
  filter(month(ymd(observation_date)) %in% 3:8) %>% 
  write_csv(paste0(data_path, data_tag, '_2023_24_all.csv'))

# save data summary

# set up column titles
dataset_summary <- c('No. Observers', 'No. Checklists', 'No. Observations')

# load in eBird data
ebd <- read_csv(paste0(data_path, data_tag, '_2023_24_all.csv'))
dataset_summary <- rbind(dataset_summary, c(length(unique(ebd$observer_id)), length(unique(ebd$checklist_id)), length(ebd$checklist_id)))

# save to CSV
write_csv(as.data.frame(dataset_summary), file = paste0(data_path, data_tag, '_model_summary.csv'))



