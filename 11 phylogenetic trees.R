## ---------------------------
##
## Script name: phylogenetic tree plots
##
## Purpose of script: Phylogenetic tree plots of results
##
## Author: ****
##
## Date Created: 23/07/2025
##
## Email: ****
##
## ---------------------------
##
## Notes: This code is a bit repetitive because it was difficult to get one
##        tree to work and I didn't want anything to break!
##   
##
## ---------------------------

## set working directory for Mac and PC

#setwd("~/Google Drive/")  		# working directory (mac)
setwd("G:/R Coding/MSID")  	# working directory (windows)

## ---------------------------

## load up the packages we will need:

library(tidyverse)
library(devtools)
# install_github("eliotmiller/clootl")
library(clootl)
library(ggtree)
library(taxize)
library(treeio)
library(scales)

# resolve namespace conflicts
select <- dplyr::select

## ---------------------------

# data tag
data_string <- "bcr23_MSID"


# datasets path
data_path <- 'analysis_data/'

# where I want to save results
results_path <- 'phylo_trees/'

taxon <- read_csv(paste0(data_path, 'eBird_taxonomy.csv')) %>%
  select(species_code = SPECIES_CODE, common_name = PRIMARY_COM_NAME, scientific_name = SCI_NAME)

species <- read_csv(paste0('community/prop/prop_species.csv')) %>%
  distinct(common_name) %>% left_join(taxon, by = 'common_name')

orders <- read_csv(paste0(data_path, 'Clements2024.csv')) %>%
  select(species = `scientific name`, order)


#-Individual Model ####
species <- read_csv(paste0('individual/index_prop.csv')) %>%
  distinct(common_name) %>% left_join(taxon, by = 'common_name')

tree <- extractTree(species = species$scientific_name)
#- High experience ####

prop_species <- read_csv('individual/index_prop.csv') %>% 
  left_join(species, by = 'common_name') %>% filter(checklist_labels == 'high') %>% select(species = scientific_name, index = new_diff)

group_info <- data.frame(
  id = tree$tip.label,
  species = gsub('_', ' ', tree$tip.label)
) %>%
  left_join(prop_species) %>%
  left_join(orders) %>%
  mutate(order = factor(order))

# Count species per order
order_counts <- group_info %>%
  count(order, name = "n_species")

# Split into multi- and single-species orders
multi_species <- order_counts %>% filter(n_species > 1) 
single_species <- order_counts %>% filter(n_species == 1)

label_data <- group_info %>% filter(order %in% single_species$order)

p <- ggtree(tree, layout = 'circular', linewidth = 0.3)

tree_data <- p$data

# Filter tips only and calculate hjust based on angle
tip_labels <- tree_data %>%
  filter(isTip) %>%
  mutate(
    angle = angle,
    hjust = ifelse(angle < 90 | angle > 270, 0, 1)  # Right-align left side, left-align right side
  ) %>% select(id = label, hjust, angle)

strip_data <- group_info %>% filter(order %in% multi_species$order) %>%
  left_join(tip_labels) %>%
  group_by(order) %>%
  arrange(angle) %>%
  summarise(
    taxa1 = first(id),
    taxa2 = last(id),
    hjust = max(hjust),
    index = mean(index),
    .groups = "drop"
  )
strip_data$order_factor <- factor(strip_data$order, levels = unique(strip_data$order))

p <- p %<+% group_info + 
  geom_tippoint(aes(x = x + 1.5,
                    color = index, 
                    fill = index), 
                size = 0.65) +
  scale_color_gradient2("Change in Relative\nReporting Rate",
                        low = "darkblue",   # Low end of the scale
                        mid = "white",       # Midpoint color
                        high = "darkred",    # High end of the scale
                        limits = c(-1, 1),
                        midpoint = 0
  ) +
  scale_fill_gradient2("Change in Relative\nReporting Rate",
                       low = "darkblue",   # Low end of the scale
                       mid = "white",       # Midpoint color
                       high = "darkred",    # High end of the scale
                       limits = c(-1, 1),
                       midpoint = 0
  ) +
  theme(
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank()  #remove minor gridlines
  ) + ggtitle('High experience Observers')

# 2. Create the color scale function with fixed endpoints
col_fun <- scales::gradient_n_pal(c("darkblue", "white", "darkred"))

# 3. Normalize strip_data$index to [0, 1] using fixed limits
norm_audio <- scales::rescale(strip_data$index, to = c(0, 1), from = c(-1, 1))

# 4. Apply the color function
strip_colors <- col_fun(norm_audio)


for (i in seq_len(nrow(strip_data))) {
  p <- p + geom_strip(
    taxa1 = strip_data$taxa1[i],
    taxa2 = strip_data$taxa2[i],
    label = strip_data$order[i],
    color = strip_colors[i],
    fill = strip_colors[i],
    inherit.aes = FALSE,
    barsize = 2,
    offset = 4,
    offset.text = 5,
    fontsize = 2,
    hjust = strip_data$hjust[i]
  )
}

p
ggsave(paste0(results_path, 'individual/phylo_tree_high.png'), p, width = 6, height = 4)


#-- Low experience  ####

prop_species <- read_csv('individual/index_prop.csv') %>% 
  left_join(species, by = 'common_name') %>% filter(checklist_labels == 'low') %>% select(species = scientific_name, index = new_diff)

group_info <- data.frame(
  id = tree$tip.label,
  species = gsub('_', ' ', tree$tip.label)
) %>%
  left_join(prop_species) %>%
  left_join(orders) %>%
  mutate(order = factor(order))

# Count species per order
order_counts <- group_info %>%
  count(order, name = "n_species")

# Split into multi- and single-species orders
multi_species <- order_counts %>% filter(n_species > 1) 
single_species <- order_counts %>% filter(n_species == 1)

label_data <- group_info %>% filter(order %in% single_species$order)

p <- ggtree(tree, layout = 'circular', linewidth = 0.3)

tree_data <- p$data

# Filter tips only and calculate hjust based on angle
tip_labels <- tree_data %>%
  filter(isTip) %>%
  mutate(
    angle = angle,
    hjust = ifelse(angle < 90 | angle > 270, 0, 1)  # Right-align left side, left-align right side
  ) %>% select(id = label, hjust, angle)

strip_data <- group_info %>% filter(order %in% multi_species$order) %>%
  left_join(tip_labels) %>%
  group_by(order) %>%
  arrange(angle) %>%
  summarise(
    taxa1 = first(id),
    taxa2 = last(id),
    hjust = max(hjust),
    index = mean(index),
    .groups = "drop"
  )
strip_data$order_factor <- factor(strip_data$order, levels = unique(strip_data$order))

p <- p %<+% group_info + 
  geom_tippoint(aes(x = x + 1.5,
                    color = index, 
                    fill = index), 
                size = 0.65) +
  scale_color_gradient2("Change in Relative\nReporting Rate",
                        low = "darkblue",   # Low end of the scale
                        mid = "white",       # Midpoint color
                        high = "darkred",    # High end of the scale
                        limits = c(-1, 1),
                        midpoint = 0
  ) +
  scale_fill_gradient2("Change in Relative\nReporting Rate",
                       low = "darkblue",   # Low end of the scale
                       mid = "white",       # Midpoint color
                       high = "darkred",    # High end of the scale
                       limits = c(-1, 1),
                       midpoint = 0
  ) +
  theme(
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank()  #remove minor gridlines
  ) + ggtitle('Low experience Observers')

# 2. Create the color scale function with fixed endpoints
col_fun <- scales::gradient_n_pal(c("darkblue", "white", "darkred"))

# 3. Normalize strip_data$index to [0, 1] using fixed limits
norm_audio <- scales::rescale(strip_data$index, to = c(0, 1), from = c(-1, 1))

# 4. Apply the color function
strip_colors <- col_fun(norm_audio)


for (i in seq_len(nrow(strip_data))) {
  p <- p + geom_strip(
    taxa1 = strip_data$taxa1[i],
    taxa2 = strip_data$taxa2[i],
    label = strip_data$order[i],
    color = strip_colors[i],
    fill = strip_colors[i],
    inherit.aes = FALSE,
    barsize = 2,
    offset = 4,
    offset.text = 5,
    fontsize = 2,
    hjust = strip_data$hjust[i]
  )
}
p
ggsave(paste0(results_path, 'individual/phylo_tree_low.png'), p, width = 6, height = 4)


#- Community model ####
#- High experience ####

prop_species <- read_csv('community/prop/prop_species.csv') %>% 
  left_join(species, by = 'common_name') %>% filter(checklist_label == 'high') %>% select(species = scientific_name, index = new_diff)

tree <- extractTree(species = species$scientific_name)

group_info <- data.frame(
  id = tree$tip.label,
  species = gsub('_', ' ', tree$tip.label)
) %>%
  left_join(prop_species) %>%
  left_join(orders) %>%
  mutate(order = factor(order))

# Count species per order
order_counts <- group_info %>%
  count(order, name = "n_species")

# Split into multi- and single-species orders
multi_species <- order_counts %>% filter(n_species > 1) 
single_species <- order_counts %>% filter(n_species == 1)

label_data <- group_info %>% filter(order %in% single_species$order)

p <- ggtree(tree, layout = 'circular', linewidth = 0.3)

tree_data <- p$data

# Filter tips only and calculate hjust based on angle
tip_labels <- tree_data %>%
  filter(isTip) %>%
  mutate(
    angle = angle,
    hjust = ifelse(angle < 90 | angle > 270, 0, 1)  # Right-align left side, left-align right side
  ) %>% select(id = label, hjust, angle)

strip_data <- group_info %>% filter(order %in% multi_species$order) %>%
  left_join(tip_labels) %>%
  group_by(order) %>%
  arrange(angle) %>%
  summarise(
    taxa1 = first(id),
    taxa2 = last(id),
    hjust = max(hjust),
    index = mean(index),
    .groups = "drop"
  )
strip_data$order_factor <- factor(strip_data$order, levels = unique(strip_data$order))

p <- p %<+% group_info + 
  geom_tippoint(aes(x = x + 1.5,
                    color = index, 
                    fill = index), 
                size = 0.65) +
  scale_color_gradient2("Change in Relative\nReporting Rate",
                        low = "darkblue",   # Low end of the scale
                        mid = "white",       # Midpoint color
                        high = "darkred",    # High end of the scale
                        limits = c(-1, 1),
                        midpoint = 0
  ) +
  scale_fill_gradient2("Change in Relative\nReporting Rate",
                       low = "darkblue",   # Low end of the scale
                       mid = "white",       # Midpoint color
                       high = "darkred",    # High end of the scale
                       limits = c(-1, 1),
                       midpoint = 0
  ) +
  theme(
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank()  #remove minor gridlines
  ) + ggtitle('High experience Observers')

# 2. Create the color scale function with fixed endpoints
col_fun <- scales::gradient_n_pal(c("darkblue", "white", "darkred"))

# 3. Normalize strip_data$index to [0, 1] using fixed limits
norm_audio <- scales::rescale(strip_data$index, to = c(0, 1), from = c(-1, 1))

# 4. Apply the color function
strip_colors <- col_fun(norm_audio)


for (i in seq_len(nrow(strip_data))) {
  p <- p + geom_strip(
    taxa1 = strip_data$taxa1[i],
    taxa2 = strip_data$taxa2[i],
    label = strip_data$order[i],
    color = strip_colors[i],
    fill = strip_colors[i],
    inherit.aes = FALSE,
    barsize = 2,
    offset = 4,
    offset.text = 5,
    fontsize = 2,
    hjust = strip_data$hjust[i]
  )
}

p
ggsave(paste0(results_path, '/community/phylo_tree_high.png'), p, width = 6, height = 4)


#-- Low experience  ####

prop_species <- read_csv('community/prop/prop_species.csv') %>% 
  left_join(species, by = 'common_name') %>% filter(checklist_label == 'low') %>% select(species = scientific_name, index = new_diff)

tree <- extractTree(species = species$scientific_name)

group_info <- data.frame(
  id = tree$tip.label,
  species = gsub('_', ' ', tree$tip.label)
) %>%
  left_join(prop_species) %>%
  left_join(orders) %>%
  mutate(order = factor(order))

# Count species per order
order_counts <- group_info %>%
  count(order, name = "n_species")

# Split into multi- and single-species orders
multi_species <- order_counts %>% filter(n_species > 1) 
single_species <- order_counts %>% filter(n_species == 1)

label_data <- group_info %>% filter(order %in% single_species$order)

p <- ggtree(tree, layout = 'circular', linewidth = 0.3)

tree_data <- p$data

# Filter tips only and calculate hjust based on angle
tip_labels <- tree_data %>%
  filter(isTip) %>%
  mutate(
    angle = angle,
    hjust = ifelse(angle < 90 | angle > 270, 0, 1)  # Right-align left side, left-align right side
  ) %>% select(id = label, hjust, angle)

strip_data <- group_info %>% filter(order %in% multi_species$order) %>%
  left_join(tip_labels) %>%
  group_by(order) %>%
  arrange(angle) %>%
  summarise(
    taxa1 = first(id),
    taxa2 = last(id),
    hjust = max(hjust),
    index = mean(index),
    .groups = "drop"
  )
strip_data$order_factor <- factor(strip_data$order, levels = unique(strip_data$order))


p <- p %<+% group_info + 
  geom_tippoint(aes(x = x + 1.5,
                    color = index, 
                    fill = index), 
                size = 0.65) +
  scale_color_gradient2("Change in Relative\nReporting Rate",
                        low = "darkblue",   # Low end of the scale
                        mid = "white",       # Midpoint color
                        high = "darkred",    # High end of the scale
                        limits = c(-1, 1),
                        midpoint = 0
  ) +
  scale_fill_gradient2("Change in Relative\nReporting Rate",
                       low = "darkblue",   # Low end of the scale
                       mid = "white",       # Midpoint color
                       high = "darkred",    # High end of the scale
                       limits = c(-1, 1),
                       midpoint = 0
  ) +
  theme(
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank()  #remove minor gridlines
  ) + ggtitle('Low experience Observers')

# 2. Create the color scale function with fixed endpoints
col_fun <- scales::gradient_n_pal(c("darkblue", "white", "darkred"))

# 3. Normalize strip_data$index to [0, 1] using fixed limits
norm_audio <- scales::rescale(strip_data$index, to = c(0, 1), from = c(-1, 1))

# 4. Apply the color function
strip_colors <- col_fun(norm_audio)


for (i in seq_len(nrow(strip_data))) {
  p <- p + geom_strip(
    taxa1 = strip_data$taxa1[i],
    taxa2 = strip_data$taxa2[i],
    label = strip_data$order[i],
    color = strip_colors[i],
    fill = strip_colors[i],
    inherit.aes = FALSE,
    barsize = 2,
    offset = 4,
    offset.text = 5,
    fontsize = 2,
    hjust = strip_data$hjust[i]
  )
}

p
ggsave(paste0(results_path, '/community/phylo_tree_low.png'), p, width = 6, height = 4)


#- Community model max ####
#- High experience ####

prop_species <- read_csv('community_max/prop/prop_species.csv') %>% 
  left_join(species, by = 'common_name') %>% filter(checklist_label == 'high') %>% select(species = scientific_name, index = new_diff)

tree <- extractTree(species = species$scientific_name)

group_info <- data.frame(
  id = tree$tip.label,
  species = gsub('_', ' ', tree$tip.label)
) %>%
  left_join(prop_species) %>%
  left_join(orders) %>%
  mutate(order = factor(order))

# Count species per order
order_counts <- group_info %>%
  count(order, name = "n_species")

# Split into multi- and single-species orders
multi_species <- order_counts %>% filter(n_species > 1) 
single_species <- order_counts %>% filter(n_species == 1)

label_data <- group_info %>% filter(order %in% single_species$order)

p <- ggtree(tree, layout = 'circular', linewidth = 0.3)

tree_data <- p$data

# Filter tips only and calculate hjust based on angle
tip_labels <- tree_data %>%
  filter(isTip) %>%
  mutate(
    angle = angle,
    hjust = ifelse(angle < 90 | angle > 270, 0, 1)  # Right-align left side, left-align right side
  ) %>% select(id = label, hjust, angle)

strip_data <- group_info %>% filter(order %in% multi_species$order) %>%
  left_join(tip_labels) %>%
  group_by(order) %>%
  arrange(angle) %>%
  summarise(
    taxa1 = first(id),
    taxa2 = last(id),
    hjust = max(hjust),
    index = mean(index),
    .groups = "drop"
  )
strip_data$order_factor <- factor(strip_data$order, levels = unique(strip_data$order))

p <- p %<+% group_info + 
  geom_tippoint(aes(x = x + 1.5,
                    color = index, 
                    fill = index), 
                size = 0.65) +
  scale_color_gradient2("Change in Relative\nReporting Rate",
                        low = "darkblue",   # Low end of the scale
                        mid = "white",       # Midpoint color
                        high = "darkred",    # High end of the scale
                        limits = c(-1, 1),
                        midpoint = 0
  ) +
  scale_fill_gradient2("Change in Relative\nReporting Rate",
                       low = "darkblue",   # Low end of the scale
                       mid = "white",       # Midpoint color
                       high = "darkred",    # High end of the scale
                       limits = c(-1, 1),
                       midpoint = 0
  ) +
  theme(
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank()  #remove minor gridlines
  ) + ggtitle('High experience Observers')

# 2. Create the color scale function with fixed endpoints
col_fun <- scales::gradient_n_pal(c("darkblue", "white", "darkred"))

# 3. Normalize strip_data$index to [0, 1] using fixed limits
norm_audio <- scales::rescale(strip_data$index, to = c(0, 1), from = c(-1, 1))

# 4. Apply the color function
strip_colors <- col_fun(norm_audio)


for (i in seq_len(nrow(strip_data))) {
  p <- p + geom_strip(
    taxa1 = strip_data$taxa1[i],
    taxa2 = strip_data$taxa2[i],
    label = strip_data$order[i],
    color = strip_colors[i],
    fill = strip_colors[i],
    inherit.aes = FALSE,
    barsize = 2,
    offset = 4,
    offset.text = 5,
    fontsize = 2,
    hjust = strip_data$hjust[i]
  )
}
p
ggsave(paste0(results_path, '/community_max/phylo_tree_high.png'), p, width = 6, height = 4)


#-- Low experience  ####

prop_species <- read_csv('community_max/prop/prop_species.csv') %>% 
  left_join(species, by = 'common_name') %>% filter(checklist_label == 'low') %>% select(species = scientific_name, index = new_diff)

tree <- extractTree(species = species$scientific_name)

group_info <- data.frame(
  id = tree$tip.label,
  species = gsub('_', ' ', tree$tip.label)
) %>%
  left_join(prop_species) %>%
  left_join(orders) %>%
  mutate(order = factor(order))

# Count species per order
order_counts <- group_info %>%
  count(order, name = "n_species")

# Split into multi- and single-species orders
multi_species <- order_counts %>% filter(n_species > 1) 
single_species <- order_counts %>% filter(n_species == 1)

label_data <- group_info %>% filter(order %in% single_species$order)

p <- ggtree(tree, layout = 'circular', linewidth = 0.3)

tree_data <- p$data

# Filter tips only and calculate hjust based on angle
tip_labels <- tree_data %>%
  filter(isTip) %>%
  mutate(
    angle = angle,
    hjust = ifelse(angle < 90 | angle > 270, 0, 1)  # Right-align left side, left-align right side
  ) %>% select(id = label, hjust, angle)

strip_data <- group_info %>% filter(order %in% multi_species$order) %>%
  left_join(tip_labels) %>%
  group_by(order) %>%
  arrange(angle) %>%
  summarise(
    taxa1 = first(id),
    taxa2 = last(id),
    hjust = max(hjust),
    index = mean(index),
    .groups = "drop"
  )
strip_data$order_factor <- factor(strip_data$order, levels = unique(strip_data$order))


p <- p %<+% group_info + 
  geom_tippoint(aes(x = x + 1.5,
                    color = index, 
                    fill = index), 
                size = 0.65) +
  scale_color_gradient2("Change in Relative\nReporting Rate",
                        low = "darkblue",   # Low end of the scale
                        mid = "white",       # Midpoint color
                        high = "darkred",    # High end of the scale
                        limits = c(-1, 1),
                        midpoint = 0
  ) +
  scale_fill_gradient2("Change in Relative\nReporting Rate",
                       low = "darkblue",   # Low end of the scale
                       mid = "white",       # Midpoint color
                       high = "darkred",    # High end of the scale
                       limits = c(-1, 1),
                       midpoint = 0
  ) +
  theme(
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank()  #remove minor gridlines
  ) + ggtitle('Low experience Observers')

# 2. Create the color scale function with fixed endpoints
col_fun <- scales::gradient_n_pal(c("darkblue", "white", "darkred"))

# 3. Normalize strip_data$index to [0, 1] using fixed limits
norm_audio <- scales::rescale(strip_data$index, to = c(0, 1), from = c(-1, 1))

# 4. Apply the color function
strip_colors <- col_fun(norm_audio)


for (i in seq_len(nrow(strip_data))) {
  p <- p + geom_strip(
    taxa1 = strip_data$taxa1[i],
    taxa2 = strip_data$taxa2[i],
    label = strip_data$order[i],
    color = strip_colors[i],
    fill = strip_colors[i],
    inherit.aes = FALSE,
    barsize = 2,
    offset = 4,
    offset.text = 5,
    fontsize = 2,
    hjust = strip_data$hjust[i]
  )
}

p
ggsave(paste0(results_path, '/community_max/phylo_tree_low.png'), p, width = 6, height = 4)

#- Audio Index  ####

audio_index <- read_csv(paste0(data_path, 'audio_index.csv')) %>%
  mutate(audio_index = audio/(audio + photo)) %>%
  left_join(read_csv(paste0(data_path, 'eBird_taxonomy.csv')) %>% 
              select(species_code = SPECIES_CODE, common_name = PRIMARY_COM_NAME)) %>%
  select(common_name, audio_index)

prop_species <- read_csv(paste0('community_max/prop/prop_species.csv')) %>%
  distinct(common_name) %>% left_join(taxon, by = 'common_name') %>% 
  left_join(audio_index, by = 'common_name') %>% 
  mutate(index = log(audio_index) - min(log(audio_index))) %>%
  rename(species = scientific_name)

tree <- extractTree(species = species$scientific_name)

group_info <- data.frame(
  id = tree$tip.label,
  species = gsub('_', ' ', tree$tip.label)
) %>%
  left_join(prop_species) %>%
  left_join(orders) %>%
  mutate(order = factor(order))

# Count species per order
order_counts <- group_info %>%
  count(order, name = "n_species")

# Split into multi- and single-species orders
multi_species <- order_counts %>% filter(n_species > 1) 
single_species <- order_counts %>% filter(n_species == 1)

label_data <- group_info %>% filter(order %in% single_species$order)


p <- ggtree(tree, layout = 'circular', linewidth = 0.3)

tree_data <- p$data

# Filter tips only and calculate hjust based on angle
tip_labels <- tree_data %>%
  filter(isTip) %>%
  mutate(
    angle = angle,
    hjust = ifelse(angle < 90 | angle > 270, 0, 1)  # Right-align left side, left-align right side
  ) %>% select(id = label, hjust, angle)

strip_data <- group_info %>% filter(order %in% multi_species$order) %>%
  left_join(tip_labels) %>%
  group_by(order) %>%
  arrange(angle) %>%
  summarise(
    taxa1 = first(id),
    taxa2 = last(id),
    hjust = max(hjust),
    index = mean(index),
    .groups = "drop"
  )
strip_data$order_factor <- factor(strip_data$order, levels = unique(strip_data$order))


p <- p %<+% group_info + 
  geom_tippoint(aes(x = x + 1.5,
                    color = index, 
                    fill = index), 
                size = 0.65) +
  scale_color_viridis_c("Audio Index") +
  scale_fill_viridis_c("Audio Index") +
  theme(
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank()  #remove minor gridlines
  )

# Normalize audio_index for color mapping
norm_audio <- scales::rescale(strip_data$index, to = c(0, 1), from = c(min(group_info$index), max(group_info$index)))

# Create viridis palette function
col_fun <- viridis_pal(option = "D")

# Generate colors for each strip based on normalized audio_index
strip_colors <- col_fun(100)[floor(norm_audio * 99) + 1]  # 100 colors interpolated


for (i in seq_len(nrow(strip_data))) {
  p <- p + geom_strip(
    taxa1 = strip_data$taxa1[i],
    taxa2 = strip_data$taxa2[i],
    label = strip_data$order[i],
    color = strip_colors[i],
    fill = strip_colors[i],
    inherit.aes = FALSE,
    barsize = 2,
    offset = 4,
    offset.text = 5,
    fontsize = 2,
    hjust = strip_data$hjust[i]
  )
}

p
ggsave(paste0(results_path, 'AVI.png'), p, width = 6, height = 4)



