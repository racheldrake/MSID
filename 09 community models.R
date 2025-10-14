## ---------------------------
##
## Script name: Community models
##
## Purpose of script: Model using audio index to compare across >1% RR species
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
library(DHARMa)
library(MuMIn)

# resolve namespace conflicts
select <- dplyr::select

options(glmmTMB.cores = 8)

## ---------------------------

# data tag
data_string <- "bcr23_MSID"

# datasets path
data_path <- 'analysis_data/'

# where I want to save results
results_path <- 'community/'


# load in datasets
reporting_rates <- read_csv(paste0(data_path, 'observer_reporting_rates_occur.csv')) %>% 
  rename(n = no.checklists)

audio_index <- read_csv(paste0(data_path, 'audio_index.csv')) %>%
  mutate(audio_index = audio/(audio + photo)) %>%
  left_join(read_csv(paste0(data_path, 'eBird_taxonomy.csv')) %>% 
              select(species_code = SPECIES_CODE, common_name = PRIMARY_COM_NAME)) %>%
  select(common_name, audio_index)

# create model data frame
model_data <- reporting_rates %>% 
  group_by(common_name) %>% 
  # need to have at least 1% reporting rate to be a stable model
  filter(mean(RR) >= 0.01) %>% 
  ungroup() %>%
  drop_na() %>%
  left_join(audio_index, by = 'common_name') %>%
  filter(audio_index != 0) %>%
  mutate(common_name = factor(common_name),
         # scale audio index
         audio_index = log(audio_index) - min(log(audio_index)),
         # reduce where RR == 1
         RR = ifelse(RR == 1, 0.9999, RR))

model_prop <- glmmTMB(formula = RR ~ log10(previous_checklists+1)*MSID_prop*occur +
                        audio_index*log10(previous_checklists+1) + audio_index*MSID_prop + (1 + MSID_prop | common_name),
                      data = model_data,
                      family = beta_family(),
                      ziformula = ~ log10(previous_checklists+1)*MSID_prop*realised_occur +
                        audio_index*log10(previous_checklists+1) + audio_index*MSID_prop)

saveRDS(model_prop, paste0(results_path, 'model_prop.rds'))
model_prop <- readRDS(paste0(results_path, 'model_prop.rds'))


# model diagnostics
# need to fit model without random effects because otherwise there isn't
# enough memory
DHARMa_model <- glmmTMB(formula = RR ~ log10(previous_checklists+1)*MSID_prop*occur +
                        audio_index*log10(previous_checklists+1) + audio_index*MSID_prop,
                      data = model_data,
                      family = beta_family(),
                      ziformula = ~ log10(previous_checklists+1)*MSID_prop*realised_occur +
                        audio_index*log10(previous_checklists+1) + audio_index*MSID_prop)

saveRDS(DHARMa_model, paste0(results_path, 'model_DHARMa.rds'))
DHARMa_model <- readRDS(paste0(results_path, 'model_DHARMa.rds'))

simres <- simulateResiduals(DHARMa_model, n = 250)
png(filename = paste0(results_path, 'DHARMa_resid.png'), width = 2*480)
plot(simres)
dev.off()

# model summaries
species_summary <- data.frame(summary(fit2)$coefficient$cond, row.names = rownames(summary(fit2)$coefficient$cond))
zi_summary <- data.frame(summary(fit2)$coefficient$zi, row.names = rownames(summary(fit2)$coefficient$zi))

#- MODEL PROP ####

# Assume your model data is called `model_data`
n <- 100  # continuous audio index samples
m <- 6  # discrete levels of merlin use

# Separate zero and non-zero values
nonzero_vals <- model_data$MSID_prop[model_data$MSID_prop > 0]

# Get quantiles for non-zero values only
quantile_vals <- quantile(nonzero_vals, probs = seq(0, 1, length.out = m - 1), na.rm = TRUE)

# Combine 0 with quantile values
prop_seq <- c(0, as.numeric(quantile_vals))

# Build data frame
prop_df <- data.frame(
  MSID_prop = prop_seq,
  MSID_label = as.character(round(prop_seq, 2))
)

audio_index_seq <- seq(min(model_data$audio_index, na.rm = TRUE),
                       max(model_data$audio_index, na.rm = TRUE),
                       length.out = n)

# Get quantiles for no.checklists
no_checklists_vals <- quantile(model_data$previous_checklists, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
no_checklists_df <- data.frame(
  previous_checklists = as.numeric(no_checklists_vals),
  checklist_label = c("low", "medium", "high")
)

# Create prediction dataframe
pred_data <- expand.grid(
  MSID_prop = prop_seq,
  audio_index = audio_index_seq,
  previous_checklists = no_checklists_df$previous_checklists
)

# Add labels for previous_checklists
pred_data$checklist_label <- rep(no_checklists_df$checklist_label, each = m * n)
pred_data$MSID_label <- rep(prop_df$MSID_label)

# Add static columns
pred_data$common_name <- "new"
pred_data$observer_id <- "new"
pred_data$occur <- mean(model_data$occur)
pred_data$realised_occur <- mean(model_data$realised_occur)

pred_data$prediction <- predict(model_prop, newdata = pred_data, type = 'response', allow.new.levels = TRUE)

# Colours for MSID levels
prop_seq <- as.character(round(prop_seq, 2))
gradient_cols <- colorRampPalette(c("lightblue", "darkblue"))(length(prop_seq) - 1)
custom_cols <- c("0" = "red3", setNames(gradient_cols, prop_seq[prop_seq != 0]))


prop_graph <- pred_data %>% mutate(checklist_label = factor(checklist_label, levels = c('low', 'medium', 'high')),
                                   MSID_label = factor(MSID_label),
                                   MSID_label = fct_reorder(MSID_label, MSID_prop)) %>% ggplot() +
  geom_line(aes(audio_index, prediction, colour = MSID_label)) +
  facet_wrap(~checklist_label, labeller = as_labeller(c('low' = 'Low', 'medium' = 'Medium', 'high' = 'High'))) + theme_bw() + xlab('Audio Index') + 
  ylab('Estimated Reporting Rate') +
  scale_colour_manual(
    values = custom_cols,
    name = "Proportion of\nMerlin use") +
  ggtitle('Previously submitted checklists (expertise)') + 
  coord_cartesian(ylim = c(0, NA)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11), 
        legend.title = element_text(size = 11)) + 
  labs(colour ='Proportion of\nMerlin use') 

ggsave(paste0(results_path, 'prop/prop_graph.png'), prop_graph, width = 12, height = 4)

# conditional only
pred_data$prediction <- predict(model_prop, newdata = pred_data, type = 'conditional', allow.new.levels = TRUE)

prop_graph <- pred_data %>% mutate(checklist_label = factor(checklist_label, levels = c('low', 'medium', 'high')),
                                   MSID_label = factor(MSID_label),
                                   MSID_label = fct_reorder(MSID_label, MSID_prop)) %>% ggplot() +
  geom_line(aes(audio_index, prediction, colour = MSID_label)) +
  facet_wrap(~checklist_label, labeller = as_labeller(c('low' = 'Low', 'medium' = 'Medium', 'high' = 'High'))) + theme_bw() + xlab('Audio Index') + 
  ylab('Estimated Reporting Rate') +
  scale_colour_manual(
    values = custom_cols,
    name = "Proportion of\nMerlin use") +
  ggtitle('Previously submitted checklists (expertise)') + 
  coord_cartesian(ylim = c(0, NA)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11), 
        legend.title = element_text(size = 11)) + 
  labs(colour ='Proportion of\nMerlin use') 

ggsave(paste0(results_path, 'prop/prop_graph_cond.png'), prop_graph, width = 12, height = 4)

# and then for zero inflated portion
pred_data$prediction <- predict(model_prop, newdata = pred_data, type = 'zprob', allow.new.levels = TRUE)

prop_graph <- pred_data %>% mutate(checklist_label = factor(checklist_label, levels = c('low', 'medium', 'high')),
                                   MSID_label = factor(MSID_label),
                                   MSID_label = fct_reorder(MSID_label, MSID_prop),
                                   prediction = 1 - prediction) %>% ggplot() +
  geom_line(aes(audio_index, prediction, colour = MSID_label)) +
  facet_wrap(~checklist_label, labeller = as_labeller(c('low' = 'Low', 'medium' = 'Medium', 'high' = 'High'))) + theme_bw() + xlab('Audio Index') + 
  ylab('Estimated Reporting Rate') +
  scale_colour_manual(
    values = custom_cols,
    name = "Proportion of\nMerlin use") +
  ggtitle('Previously submitted checklists (expertise)') +
  coord_cartesian(ylim = c(0, NA)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11), 
        legend.title = element_text(size = 11)) + 
  labs(colour ='Proportion of\nMerlin use') 

ggsave(paste0(results_path, 'prop/prop_graph_zi.png'), prop_graph, width = 12, height = 4)

# prediction for each species
pred_data_spec <- data.frame(
  common_name = unique(model_data$common_name),
  audio_index = unique(model_data$audio_index)
) %>%
  left_join(model_data %>% group_by(common_name) %>% 
              summarise(occur = mean(occur),
                        realised_occur = mean(realised_occur)))%>%
  mutate(observer_id = 'new',
         n = mean(model_data$n),
         min = as.numeric(quantile(model_data$MSID_prop, 0.1)),
         max = as.numeric(quantile(model_data$MSID_prop, 0.9))) %>% 
  pivot_longer(cols = c('min', 'max'), names_to = 'MSID_vals', values_to = 'MSID_prop') %>%
  mutate(low = as.numeric(quantile(model_data$previous_checklists, 0.1)),
         medium = as.numeric(quantile(model_data$previous_checklists, 0.5)),
         high = as.numeric(quantile(model_data$previous_checklists, 0.9))) %>%
  pivot_longer(cols = c('low', 'medium', 'high'), names_to = 'checklist_label', values_to = 'previous_checklists')

pred_data_spec$common_name <- factor(
  pred_data_spec$common_name,
  levels = levels(model_prop$frame$common_name))

pred_data_spec$prediction <- predict(model_prop, newdata = pred_data_spec, type = 'response')

pred_data_spec %>% select(-MSID_prop) %>% 
  pivot_wider(names_from = MSID_vals, values_from = prediction) %>%
  mutate(diff = max - min,
         ratio = max/min,
         new_diff = ratio - 1) %>% write_csv(paste0(results_path, 'prop/prop_species.csv'))

species_prop <- pred_data_spec %>% select(-MSID_prop) %>% 
  mutate(checklist_label = factor(checklist_label, levels = c('low', 'medium', 'high'))) %>%
  pivot_wider(names_from = MSID_vals, values_from = prediction) %>%
  mutate(diff = max - min,
         ratio = max/min,
         new_diff = ratio - 1) %>% 
  ggplot() + geom_point(aes(audio_index, ratio)) + theme_bw() +
  facet_wrap(~checklist_label, labeller = as_labeller(c('low' = 'Low', 'medium' = 'Medium', 'high' = 'High'))) +
  geom_hline(yintercept = 1, colour = 'grey', linetype = 'dashed') + 
  geom_smooth(aes(audio_index, ratio), method = 'lm', colour = 'darkred') + 
  xlab('Audio Index') + ylab('Ratio of predicted reporting rate for high/low proportion of Merlin use') +
  theme(plot.title = element_text(hjust = 0.5, size = 11), 
        legend.title = element_text(size = 11)) + 
  ggtitle('Previously submitted checklists (expertise)') + xlim(c(0, 10)) + ylim(0, 4.5)

ggsave(paste0(results_path, 'prop/species_prop.png'), species_prop, width = 11, height = 4)
