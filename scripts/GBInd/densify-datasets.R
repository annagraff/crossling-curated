#### This script densifies the GBInd logical and statistical datasets
rm(list=ls())

# load packages
library(densify)
library(tidyverse)
library(gmt)
library(ggplot2)

# read functions
source("../functions.R")

# read in logical GBInd data
logical <- read.csv("../../curated_data/GBInd/logicalGBI/logicalGBI.csv") %>% select(-X)

# read in statistical GBInd data
statistical <- read.csv("../../curated_data/GBInd/statisticalGBI/statisticalGBI.csv") %>% select(-X)

# generate taxonomy matrix
taxonomy_matrix <- as_flat_taxonomy_matrix(glottolog_languoids)

# for densification, ensure all blanks, ? and "NA" are coded as NA
logical_for_pruning <- na_convert(logical)
statistical_for_pruning <- na_convert(statistical)

# specify parameters for densification 
min_variability <- 3 # each variable must have at least 3 languages in its second-largest state
density_mean <- "log_odds"

# comment on weights: GB is quite dense, and part of the "NA"s on the column/variable side in both curations is explicitly wanted
# densification should thus be biased towards the taxonomic diversity criterion, expressed in a higher weight

# run densify, set seed for reproducibility
set.seed(1111)
logical_log <-
  densify(data = logical_for_pruning,
                min_variability = min_variability,
                density_mean = density_mean,
                taxonomy = glottolog_languoids,
                taxon_id = "glottocode",
                density_mean_weights = list(coding = 0.999, taxonomy = 1))

statistical_log <-
  densify(data = statistical_for_pruning,
          min_variability = min_variability,
          density_mean = density_mean,
          taxonomy = glottolog_languoids,
          taxon_id = "glottocode",
          density_mean_weights = list(coding = 0.999, taxonomy = 1))

# prune to optima
# we include minimum row coding density, since NAs on language end should largely be random
# we include taxonomic index since densification here explicitly seeks to increase taxonomic diversity
logical_pruned <- prune(logical_log, 
                         scoring_function = n_data_points*coding_density*row_coding_density_min*taxonomic_index^3)

statistical_pruned <- prune(statistical_log, 
                            scoring_function = n_data_points*coding_density*row_coding_density_min*taxonomic_index^3)

# retrieve corresponding data from input (to re-establish differences between ? and NA)
logical_pruned <- logical[which(logical$glottocode%in%logical_pruned$glottocode), which(names(logical)%in%names(logical_pruned))]
statistical_pruned <- statistical[which(statistical$glottocode%in%statistical_pruned$glottocode), which(names(statistical)%in%names(statistical_pruned))]

# save pruned matrices
write.csv(logical_pruned,"../../curated_data/GBInd/logicalGBI/logicalGBI_pruned.csv")
write.csv(statistical_pruned,"../../curated_data/GBInd/statisticalGBI/statisticalGBI_pruned.csv")
