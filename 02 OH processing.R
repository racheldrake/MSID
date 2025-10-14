## ---------------------------
##
## Script name: Observer history processing
##
## Purpose of script: extract global checklist counts for each observer
##
## Author: ***
##
## Date Created: 15/08/2025
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac and PC

#setwd("~/Google Drive/")  		# working directory (mac)
setwd("G:/")  	# working directory (windows)

## ---------------------------

## load up the packages we will need:

library(tidyverse)
library(auk)

# resolve namespace conflicts
select <- dplyr::select


## ---------------------------

# ebird sampling path
ebd_samp_loc <- "ebd_US_relMar-2025_sampling.txt"

# datasets path
results_path <- '/R Coding/MSID/proc_data/'


# -------------------------------------------------------
# filter sampling data 

# Read sampling dataset
samp <- auk_sampling(ebd_samp_loc)

# set filters
filters <- samp %>%
  auk_year(year = c(1900:2022))

# apply filters and keep only observer_id and complete flag
f_out <- paste0(results_path, "ebd_samp_obs_alltime.txt")
samp_filtered <- auk_filter(filters, 
                            keep = c("group_identifier", "sampling_event_identifier", "observer_id", "all_species_reported", "observation_date"),
                            file = f_out, overwrite = TRUE)


# -------------------------------------------------------
# summarise number of observations

tab <- read_tsv(f_out, n_max = 5) %>%
  dplyr::select("OBSERVER ID", "ALL SPECIES REPORTED") %>%
  rename(observer_id = "OBSERVER ID") %>%
  rename(all_species_reported = "ALL SPECIES REPORTED") %>%
  mutate(all_species_reported = as.numeric(all_species_reported)) %>%
  filter(!is.na(all_species_reported)) %>%
  group_by(observer_id) %>%
  summarise(n_all = n(), n_complete = sum(all_species_reported))

write_loc <- paste0(results_path, "observer_histories_global_all.csv")
write_csv(tab, write_loc)





