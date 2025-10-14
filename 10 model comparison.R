## ---------------------------
##
## Script name: model comparison
##
## Purpose of script: comparison of individual and community model results
##
## Author: ****
##
## Date Created: 23/07/2025
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

# data tag
data_string <- "bcr23_MSID"


# datasets path
data_path <- 'analysis_data/'

# where I want to save results
results_path <- 'comparison/'


#- PROP MODEL ####
#- High experience ####

prop_species_full <- read_csv('community/prop/prop_species.csv') %>% filter(checklist_label == 'high') %>% select(species = common_name, index_full = new_diff)
prop_species_ind <- read_csv('individual/index_prop.csv') %>% filter(checklist_labels == 'high') %>% select(species = common_name, index_ind = new_diff)

joint <- prop_species_full %>% left_join(prop_species_ind, by = 'species')
lims <- range(c(joint$index_ind, joint$index_full), na.rm = TRUE)
plot <- ggplot(data = joint) +  
  geom_abline(intercept = 0, slope = 1, colour = 'grey', linetype = 'dashed') + theme_bw() + 
  geom_point(aes(index_ind, index_full)) + 
  geom_smooth(aes(index_ind, index_full), method = 'lm', colour = 'darkred') + 
  xlab('Change in RR Individual Model') + ylab('Change in RR Community Model') +
  coord_equal(xlim = lims, ylim = lims) + 
  ggtitle('Proportional Model - High Experience')

ggsave(paste0(results_path, 'prop_high.png'), plot, height = 6.5, width = 6)

#- Low experience ####

prop_species_full <- read_csv('community/prop/prop_species.csv') %>% filter(checklist_label == 'low') %>% select(species = common_name, index_full = new_diff)
prop_species_ind <- read_csv('individual/index_prop.csv') %>% filter(checklist_labels == 'low') %>% select(species = common_name, index_ind = new_diff)

joint <- prop_species_full %>% left_join(prop_species_ind, by = 'species')
lims <- range(c(joint$index_ind, joint$index_full), na.rm = TRUE)
plot <- ggplot(data = joint) + 
  geom_abline(intercept = 0, slope = 1, colour = 'grey', linetype = 'dashed') + theme_bw() + 
  geom_point(aes(index_ind, index_full)) + 
  geom_smooth(aes(index_ind, index_full), method = 'lm', colour = 'darkred') + 
  xlab('Change in RR Individual Model') + ylab('Change in RR Community Model') +
  coord_equal(xlim = lims, ylim = lims) + 
  ggtitle('Proportional Model - Low Experience')

ggsave(paste0(results_path, 'prop_low.png'), plot, height = 6.5, width = 6)


#- TOTAL MODEL ####
#- High experience ####

total_species_full <-read_csv('community/total/prop_species.csv') %>% filter(checklist_label == 'high') %>% select(species = common_name, index_full = new_diff)
total_species_ind <- read_csv('individual/index_total.csv') %>% filter(checklist_labels == 'high') %>% select(species = common_name, index_ind = new_diff)

joint <- total_species_full %>% left_join(total_species_ind, by = 'species')
lims <- range(c(joint$index_ind, joint$index_full), na.rm = TRUE)
plot <- ggplot(data = joint) +  
  geom_abline(intercept = 0, slope = 1, colour = 'grey', linetype = 'dashed') + theme_bw() + 
  geom_point(aes(index_ind, index_full)) + 
  geom_smooth(aes(index_ind, index_full), method = 'lm', colour = 'darkred') + 
  xlab('Change in RR Individual Model') + ylab('Change in RR Community Model') +
  coord_equal(xlim = lims, ylim = lims) + 
  ggtitle('Total Model - High Experience')

ggsave(paste0(results_path, 'total_high.png'), plot, height = 6.5, width = 6)

#- Low experience ####

total_species_full <- read_csv('community/total/prop_species.csv') %>% filter(checklist_label == 'low') %>% select(species = common_name, index_full = new_diff)
total_species_ind <- read_csv('individual/index_total.csv') %>% filter(checklist_labels == 'low') %>% select(species = common_name, index_ind = new_diff)

joint <- total_species_full %>% left_join(total_species_ind, by = 'species')
lims <- range(c(joint$index_ind, joint$index_full), na.rm = TRUE)
plot <- ggplot(data = joint) + 
  geom_abline(intercept = 0, slope = 1, colour = 'grey', linetype = 'dashed') + theme_bw() + 
  geom_point(aes(index_ind, index_full)) + 
  geom_smooth(aes(index_ind, index_full), method = 'lm', colour = 'darkred') + 
  xlab('Change in RR Individual Model') + ylab('Change in RR Community Model') +
  coord_equal(xlim = lims, ylim = lims) + 
  ggtitle('Total Model - Low Experience')

ggsave(paste0(results_path, 'total_low.png'), plot, height = 6.5, width = 6)