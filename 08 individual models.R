## ---------------------------
##
## Script name: individual species models
##
## Purpose of script: fit individual models for all species
##
## Author: ****
##
## Date Created: 26/06/2025
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
library(glmmTMB)

# resolve namespace conflicts
select <- dplyr::select

# helper funtion for scaling model variables
apply_scalers <- function(df, scalers, vars) {
  for (var in vars) {
    m <- scalers[[paste0(var, "_mean")]]
    s <- scalers[[paste0(var, "_sd")]]
    df[[var]] <- (df[[var]] - m) / s
  }
  df
}

options(glmmTMB.cores = 8)

## ---------------------------

# data tag
data_string <- "bcr23_MSID"


# datasets path
data_path <- 'analysis_data/'

# where I want to save results
results_path <- 'individual/'


# load in datasets
reporting_rates <- read_csv(paste0(data_path, 'observer_reporting_rates_occur.csv')) %>% rename(n = no.checklists)

# create model data frame
model_data <- reporting_rates %>% 
  group_by(common_name) %>% 
  # need to have at least 1% reporting rate to be a stable model
  filter(mean(RR) >= 0.01) %>% 
  ungroup() %>%
  drop_na() %>%
  mutate(common_name = factor(common_name),
         # reduce where RR == 1
         RR = ifelse(RR == 1, 0.9999, RR))

# inititalise empty dataframes to store species results
cols <- c('common_name', 'checklist_labels', 'diff', 'ratio', 'new_diff')
write_csv(as.data.frame(matrix(ncol = length(cols), nrow = 0, dimnames = list(NULL, cols))), paste0(results_path, 'index_total.csv'), col_names =  TRUE)
write_csv(as.data.frame(matrix(ncol = length(cols), nrow = 0, dimnames = list(NULL, cols))), paste0(results_path, 'index_prop.csv'),  col_names =  TRUE)

# find species to iterate over from model data
species <- unique(as.character(model_data$common_name))
species_results <- data.frame()

for (target_species in species){
  #- MSID PROP ####
  set.seed(1)
  
  # filter model data for individual species data
  species_data <- model_data %>% 
    filter(common_name == target_species) %>%
    mutate(RR = ifelse(RR == 1, 0.99999, RR))
  
  scale_vars <- c("realised_occur")
  
  scalers <- species_data %>%
    summarise(across(all_of(scale_vars),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd   = ~sd(.x, na.rm = TRUE))))
  
  species_data_scaled <- apply_scalers(species_data, scalers, scale_vars)
  
  species_model <- glmmTMB(formula = RR ~ log10(previous_checklists + 1)*MSID_prop*occur,
                           data = species_data_scaled,
                           family = beta_family(),
                           ziformula = ~ log10(previous_checklists + 1)*MSID_prop*realised_occur)
  
  high_mean_occur <- mean(species_data$occur[species_data$occur >= quantile(species_data$occur, 0.8)])
  low_mean_occur <- mean(species_data$occur[species_data$occur <= quantile(species_data$occur, 0.2)])
  
  high_mean_realised_occur <- mean(species_data$realised_occur[species_data$realised_occur >= quantile(species_data$realised_occur, 0.8)])
  low_mean_realised_occur <- mean(species_data$realised_occur[species_data$realised_occur <= quantile(species_data$realised_occur, 0.2)])
  
  
  species_pred <- data.frame(MSID_prop = seq(min(model_data$MSID_prop), max(model_data$MSID_prop), 
                                             length.out = 100)) %>%
    mutate(low = quantile(model_data$previous_checklists, 0.1),
           high = quantile(model_data$previous_checklists, 0.9), 
           observer_id = 'new') %>%
    pivot_longer(cols = c('low', 'high'), 
                 names_to = 'checklist_labels', 
                 values_to = 'previous_checklists') %>%
    mutate(occur = ifelse(checklist_labels == 'high', high_mean_occur, low_mean_occur),
           realised_occur = ifelse(checklist_labels == 'high', high_mean_realised_occur, low_mean_realised_occur))
  
  
  pred <- predict(species_model, newdata = apply_scalers(species_pred, scalers, scale_vars), type = 'response', allow.new.levels = TRUE, se.fit = TRUE)
  
  species_pred$prediction <- pred$fit
  species_pred$se <- pred$se.fit
  
  species_index_pred <- data.frame(n = c(1)) %>%
    mutate(low = quantile(model_data$previous_checklists, 0.1),
           high = quantile(model_data$previous_checklists, 0.9),
           min = quantile(model_data$MSID_prop, 0.1),
           max = quantile(model_data$MSID_prop, 0.9),
           observer_id = 'new') %>%
    pivot_longer(cols = c('low', 'high'), names_to = 'checklist_labels', values_to = 'previous_checklists') %>%
    pivot_longer(cols = c('min', 'max'), 
                 names_to = 'MSID_level', values_to = 'MSID_prop')  %>%
    mutate(occur = ifelse(checklist_labels == 'high', high_mean_occur, low_mean_occur),
           realised_occur = ifelse(checklist_labels == 'high', high_mean_realised_occur, low_mean_realised_occur))
  
  
  species_index_pred$prediction <- predict(species_model, newdata = apply_scalers(species_index_pred, scalers, scale_vars), type = 'response')
  
  species_index_pred %>% select(-MSID_prop) %>%
    pivot_wider(names_from = MSID_level, values_from = prediction) %>%
    mutate(diff = max - min,
           ratio = max / min,
           new_diff = ratio - 1,
           common_name = target_species) %>% select(common_name, checklist_labels, diff, ratio, new_diff) %>%
    write_csv(paste0(results_path, 'index_prop.csv'), append = TRUE)
  
  
  plot_data <- species_data %>%
    mutate(
      # Calculate cut-off's
      low  = quantile(previous_checklists, 0.2, na.rm = TRUE),
      high = quantile(previous_checklists, 0.8, na.rm = TRUE),
      # Label based on cut-off's
      category = case_when(
        previous_checklists <= low  ~ "low",
        previous_checklists >= high ~ "high"
      )
    ) %>%
    # Keep only the top and bottom 20%
    filter(!is.na(category)) %>%
    select(-low, -high) %>% # remove helper columns
    group_by(category) %>%
    mutate(
      bin = cut(MSID_prop, breaks = 10)
    ) %>%
    group_by(category, bin) %>%
    summarise(
      mean_prop = mean(MSID_prop),
      mean_RR = mean(RR),
      n = n(),
      .groups = "drop"
    ) %>% select(-bin)
  
  species_graph <- species_pred %>% 
    mutate(low_se = prediction - 1.96 * se,
           high_se = prediction + 1.96 * se,
           prediction = prediction,
           low_se = low_se,
           high_se = high_se,
           low_se = ifelse(low_se < 0, 0, low_se)) %>%
    ggplot() + geom_ribbon(aes(x = MSID_prop, ymin = low_se, ymax = high_se, 
                               fill = checklist_labels), colour = NA, alpha = 0.1) + 
    geom_point(aes(mean_prop, mean_RR, size = n, colour = category), alpha = 0.5, data = plot_data) + 
    geom_line(aes(MSID_prop, prediction, colour = checklist_labels)) + theme_bw() +
    xlab('Proportion of checklists using Merlin') + ylab('Predicted Reporting Rate') +
    labs(colour = 'Previously submitted\nchecklists (expertise)', 
         fill = 'Previously submitted\nchecklists (expertise)',
         size = 'Number of observers') + 
    ggtitle(target_species) + 
    coord_cartesian(ylim = c(0, NA)) +
    scale_colour_discrete(labels = c('low' = 'Low', 'high' = 'High')) + 
    scale_fill_discrete(labels = c('low' = 'Low', 'high' = 'High'))
  
  ggsave(paste0(results_path, 'prop/', target_species, '.png'), height = 4, width = 6)
  
  pred <- predict(species_model, newdata = apply_scalers(species_pred, scalers, scale_vars), type = 'conditional', allow.new.levels = TRUE, se.fit = TRUE)
  
  species_pred$prediction <- pred$fit
  species_pred$se <- pred$se.fit
  
  species_graph <- species_pred %>% 
    mutate(low_se = prediction - 1.96 * se,
           high_se = prediction + 1.96 * se,
           prediction = prediction,
           low_se = low_se,
           high_se = high_se,
           low_se = ifelse(low_se < 0, 0, low_se)) %>%
    ggplot() + geom_ribbon(aes(x = MSID_prop, ymin = low_se, ymax = high_se, 
                               fill = checklist_labels), colour = NA, alpha = 0.1) + 
    geom_line(aes(MSID_prop, prediction, colour = checklist_labels)) + theme_bw() +
    xlab('Proportion of checklists using Merlin') + ylab('Predicted Reporting Rate') +
    labs(colour = 'Previously submitted\nchecklists (expertise)', 
         fill = 'Previously submitted\nchecklists (expertise)',
         size = 'Number of observers') + 
    ggtitle(target_species) + 
    coord_cartesian(ylim = c(0, NA)) +
    scale_colour_discrete(labels = c('low' = 'Low', 'high' = 'High')) + 
    scale_fill_discrete(labels = c('low' = 'Low', 'high' = 'High'))
  
  ggsave(paste0(results_path, 'prop/cond/', target_species, '.png'), height = 4, width = 6)
  
  pred <- predict(species_model, newdata = apply_scalers(species_pred, scalers, scale_vars), type = 'zlink', allow.new.levels = TRUE, se.fit = TRUE)
  
  species_pred$prediction <- pred$fit
  species_pred$se <- pred$se.fit
  
  species_graph <- species_pred %>% 
    mutate(low_se = prediction - 1.96 * se,
           high_se = prediction + 1.96 * se,
           prediction = 1 - plogis(prediction),
           low_se = 1 - plogis(low_se),
           high_se = 1 - plogis(high_se),
           low_se = ifelse(low_se < 0, 0, low_se)) %>%
    ggplot() + geom_ribbon(aes(x = MSID_prop, ymin = low_se, ymax = high_se, 
                               fill = checklist_labels), colour = NA, alpha = 0.1) + 
    geom_line(aes(MSID_prop, prediction, colour = checklist_labels)) + theme_bw() +
    xlab('Proportion of checklists using Merlin') + ylab('Predicted Reporting Rate') +
    labs(colour = 'Previously submitted\nchecklists (expertise)', 
         fill = 'Previously submitted\nchecklists (expertise)',
         size = 'Number of observers') + 
    ggtitle(target_species) + 
    coord_cartesian(ylim = c(0, NA)) +
    scale_colour_discrete(labels = c('low' = 'Low', 'high' = 'High')) + 
    scale_fill_discrete(labels = c('low' = 'Low', 'high' = 'High'))
  
  ggsave(paste0(results_path, 'prop/zi/', target_species, '.png'), height = 4, width = 6)
  
  species_summary <- data.frame(summary(species_model)$coefficient$cond, row.names = rownames(summary(species_model)$coefficient$cond)) %>%
    mutate(species = target_species)
  
  species_results <- rbind(species_results, rownames_to_column(species_summary))
  
  species_results %>% write_csv(paste0(results_path, 'coef_prop.csv'))
  
}

