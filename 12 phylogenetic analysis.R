## ---------------------------
##
## Script name: Analysis 2b phylogenetic
##
## Purpose of script: Phylogenetic analysis on the Merlin effect index
##
## Author: ****
##
## Date Created: 09/08/2025
##
## Email: ****
##
## ---------------------------
##
## Notes: Merlin effect index is the predicted change in species reporting rate
##          given otherwise standardised covariates. We use trait data from
##          Schuetz Johnston 2015.
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
#install_github("eliotmiller/clootl")
library(clootl)
library(ggtree)
library(taxize)
library(treeio)
library(phytools)
library(nlme)

# resolve namespace conflicts
select <- dplyr::select

## ---------------------------

# data tag
data_string <- "bcr23_MSID"

# datasets path
data_path <- 'analysis_data/'

# where I want to save results
results_path <- 'phylo_analysis/'

# load in datasets
audio_index <- read_csv(paste0(data_path, 'audio_index.csv'))

taxon <- read_csv(paste0(data_path, 'eBird_taxonomy.csv')) %>%
  select(species_code = SPECIES_CODE, common_name = PRIMARY_COM_NAME, scientific_name = SCI_NAME)

# list of all species in region with corrected names
species <- read_csv(paste0('individual/index_prop.csv')) %>%
  distinct(common_name) %>% left_join(taxon, by = 'common_name')

# load in trait data for species
trait_data <- read_csv(paste0(data_path, 'traits.csv')) %>% 
  select(common_name = common.name, mass = log.10.mass, 
         scientific_name = scientific.name,
         crest, contrast = max.color.contrast) %>%
  mutate(common_name = str_to_title(common_name),
         common_name = str_replace_all(common_name, "-(\\p{L})", function(m) {          # Match hyphen + letter
           paste0("-", tolower(str_sub(m, 2, 2)))            # Hyphen + lowercase letter
         }),
         common_name = recode(common_name, 'Herring Gull' = 'American Herring Gull','Eastern Wood-pewee' = 'Eastern Wood-Pewee'))

# load in the merlin index data
prop_species <- read_csv('individual/index_prop.csv') %>% 
  left_join(species, by = 'common_name') %>% 
  left_join(audio_index, by = 'common_name') %>%
  filter(checklist_labels == 'high') %>% 
  select(species = scientific_name, index = new_diff, audio_index)

traits <- species %>% left_join(trait_data %>% select(-scientific_name), 
                                by = 'common_name')

# remove anything as a character
traits <- traits %>% mutate(crest = factor(crest, levels = c('none', 'short', 'long'))) %>%
  left_join(prop_species, by = join_by(scientific_name == species)) %>%
  drop_na() %>%
  mutate(label = gsub(' ', '_', scientific_name)) %>%
  select(-c(common_name, scientific_name))

# load in phylogenetic tree for relevant species
tree <- extractTree(species = gsub('_', ' ', traits$label))

# To perform PGLS in R, we must first estimate the phylogenetic covariance matrix V (C in matrix form):
spc <- tree$tip.label
V<-corBrownian(phy = tree,form = ~spc)

fit <- gls(index ~ mass + crest + contrast + audio_index, correlation = V, data=traits)
summary <- summary(fit)
write_csv(data.frame(cbind(rownames(summary$tTable),summary$tTable)), paste0(results_path, 'ind_model_prop_high.csv'))

# load in the merlin index data
prop_species <- read_csv('individual/index_prop.csv')%>% 
  left_join(species, by = 'common_name') %>% 
  left_join(audio_index, by = 'common_name') %>%
  filter(checklist_labels == 'low') %>% 
  select(species = scientific_name, index = new_diff, audio_index)

traits <- species %>% left_join(trait_data %>% select(-scientific_name), 
                                by = 'common_name')

# remove anything as a character
traits <- traits %>% mutate(crest = factor(crest, levels = c('none', 'short', 'long'))) %>%
  left_join(prop_species, by = join_by(scientific_name == species)) %>%
  drop_na() %>%
  mutate(label = gsub(' ', '_', scientific_name)) %>%
  select(-c(common_name, scientific_name))

# load in phylogenetic tree for relevant species
tree <- extractTree(species = gsub('_', ' ', traits$label))

# To perform PGLS in R, we must first estimate the phylogenetic covariance matrix V (C in matrix form):
spc <- tree$tip.label
V<-corBrownian(phy = tree,form = ~spc)

fit <- gls(index ~ mass + crest + contrast + audio_index, correlation = V, data=traits)
summary <- summary(fit)
write_csv(data.frame(cbind(rownames(summary$tTable),summary$tTable)), paste0(results_path, 'ind_model_prop_low.csv'))
