## ---------------------------
##
## Script name: MSID dataset processing
##
## Purpose of script: aggregate MSID status to observer level
##
## Author: ****
##
## Date Created: 07/07/2025
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
## ---------------------------

## set working directory for Mac and PC

#setwd("~/Google Drive/")  		# working directory (mac)
setwd("G:/R Coding/MSID")  	# working directory (windows)

## ---------------------------

## load up the packages we will need:

# library(geodata)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select

## ---------------------------

# data tag
data_string <- "bcr23_MSID"

# SAC path
SAC_path <- 'SAC/'

# MSID path
MSID_path <- 'data/'

# where the processed eBird datasets are
ebird_path <- 'proc_data/'

# where I want the outputs to go
output_path <- 'proc_data/'


# To make a Merlin usage total I need to separate each observer in
# the 2023-24 dataset and then use total their MSID column

MSID <- read_csv(paste0(MSID_path, 'MSID.csv')) %>%
  select(checklist_id = CHECKLIST_ID, MSID_bool = WasMSIDrunning)
data <- read_csv(paste0(ebird_path, data_string, '_2023_24_all.csv'))

# make dataframe of observer MSID levels
data_proc <- data %>% 
  select(observer_id, checklist_id) %>%
  # keep one row per checklist (multiples bc many bird recorded)
  group_by(checklist_id) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  # join the MSID status to checklists
  left_join(MSID, by = 'checklist_id') %>%
  # remove records without MSID info
  drop_na() %>%
  # summarise for each observer
  group_by(observer_id) %>%
  mutate(MSID_total = sum(MSID_bool),
         MSID_prop = sum(MSID_bool) / n()) %>%
  slice_sample(n = 1) %>%
  select(-checklist_id, -MSID_bool) %>% 
  # save file
  write_csv(paste0(output_path, 'observer_MSID.csv'))





