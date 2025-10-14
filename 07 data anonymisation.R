## ---------------------------
##
## Script name: Reporting rate processing
##
## Purpose of script: Anonymise all MSID datasets used in analyses
##
## Author: ****
##
## Date Created: 15/07/2025
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

## Path locations

# data tag
data_string <- "bcr23_MSID"


# datasets path
data_path <- 'proc_data/'

# where I want to save results
output_path <- 'analysis_data/'

## ---------------------------

# Firstly, we need to find the observers which should be anonymised. I had
# already aggregated this information here

# Overall file with observer MSID
MSID <- read_csv(paste0(data_path, 'observer_MSID.csv'))
MSID_key <- data.frame(observer_id = sample(MSID$observer_id), pseudo_id = as.character(seq(from = 1, to = length(MSID$observer_id))))

# Save outside the repository
setwd("G:/")
MSID_key %>% write_csv('MSID_key.csv')
setwd("G:/R Coding/MSID")

# Save anonymised MSID
MSID %>% left_join(MSID_key, by = 'observer_id') %>% 
  select(-observer_id) %>% 
  rename(observer_id = pseudo_id) %>%
  write_csv(paste0(output_path, 'observer_MSID.csv'))

# Model data
data_path = 'proc_data/'

reporting_rates <- read_csv(paste0(data_path, 'observer_reporting_rates_occur.csv'))

reporting_rates %>% left_join(MSID_key, by = 'observer_id') %>% 
  select(-observer_id) %>% 
  rename(observer_id = pseudo_id) %>%
  write_csv(paste0(output_path, 'observer_reporting_rates_occur.csv'))


