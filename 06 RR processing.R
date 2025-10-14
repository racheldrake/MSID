## ---------------------------
##
## Script name: Reporting rate processing
##
## Purpose of script: Calculate each individuals reporting rate for analysis 2b MSID
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

# load in observers MSID score
MSID <- read_csv(paste0(ebd_path, 'observer_MSID.csv'))

# load in observer histories
ebd_previous <- read_csv(paste0(ebd_path, 'observer_histories_global_all.csv')) %>% 
  select(observer_id, previous_checklists = n_complete)

# load in species availability
ebd_checklist <- read_csv(paste0(ebd_path, 'species_occurrence.csv'))


#- Observer-species reporting rates ####

# filter the dataset to give us the reporting rate for each species
ebd_filter <- ebd %>% select(observer_id, checklist_id, locality_id, common_name) %>%
  filter(observer_id %in% MSID$observer_id) %>%
  group_by(observer_id) %>%
  # calculate the denominator of the RR
  mutate(no.checklists = n_distinct(checklist_id)) %>%
  ungroup() %>%
  filter(no.checklists >= 10) %>%
  group_by(observer_id, common_name) %>%
  summarise(no.checklists = mean(no.checklists),
            successes = n()) %>% ungroup()

# fill in zeros where observers haven't seen a species
ebd_filter <- ebd_filter %>%
  pivot_wider(
    names_from = common_name,
    values_from = successes,
    values_fill = 0) %>%
  pivot_longer(
    cols = -c(observer_id, no.checklists),
    names_to = "common_name",
    values_to = "successes"
  ) %>% mutate(failures = no.checklists - successes,
               RR = successes/no.checklists)

#- model data processing ####
 
# add species availability to reoprting rates
ebd_filter <- ebd_filter %>% left_join(ebd_checklist, by = c('observer_id', 'common_name'))

# add observer history to reporting rates and species availability
ebd_filter <- ebd_filter %>% left_join(ebd_previous, by = 'observer_id') 

# zero fill previous checklists
ebd_filter <- ebd_filter %>% mutate(previous_checklists = ifelse(is.na(previous_checklists), 0, previous_checklists))

# add observer MSID score to SA and OH
ebd_filter <- ebd_filter %>% left_join(MSID, by = 'observer_id')

# save data
ebd_filter %>% write_csv(paste0(results_path, 'observer_reporting_rates_occur.csv'))


