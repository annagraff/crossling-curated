#### This script densifies the GBInd logical and statistical datasets

rm(list=ls())

# load packages
library(densify)
library(tidyverse)
library(gmt)
library(ggplot2)

# read functions
source("../../functions.R")

# read in logical GBInd data
logical <- read.csv("output/logicalGBI/logicalGBI.csv") %>% select(-X)
lgs_lo <- logical$glottocode
logical <- logical %>% select(-glottocode)
rownames(logical) <- lgs_lo

# read in statistical GBInd data
statistical <- read.csv("output/statisticalGBI/statisticalGBI.csv") %>% select(-X)
lgs_st <- statistical$glottocode
statistical <- statistical %>% select(-glottocode)
rownames(statistical) <- lgs_st

# generate taxonomy matrix
taxonomy_matrix <- as_flat_taxonomy_matrix(glottolog_languoids)

# subset logical and statistical GBInd data to those languages with language coordinates (for plotting later)
macroareas <- read.csv("../../crossling-curated/input/glottolog_v.4.8/languages_and_dialects_geo.csv")
names(macroareas)[1]<-"id"
tax_with_coords <- merge(taxonomy_matrix,select(macroareas, c("id","longitude","latitude")),by="id")
tax_with_coords <- tax_with_coords %>% filter(!is.na(latitude))
logical_for_pruning <- logical[which(rownames(logical)%in%tax_with_coords$id),]
statistical_for_pruning <- statistical[which(rownames(statistical)%in%tax_with_coords$id),]

# for densification, ensure all blanks, ? and "NA" are coded as NA
logical_for_pruning <- na_convert(logical_for_pruning)
statistical_for_pruning <- na_convert(statistical_for_pruning)

# specify parameters for densification 
variability_threshold <- 3 # each variable must have at least 3 languages in its second-largest state
mean_type <- "log_odds"

# comment on weights: GB is quite dense, and part of the "NA"s on the column/variable side in both curations is explicitly wanted
# densification should thus be biased towards the taxonomic diversity criterion, expressed in a higher weight
taxonomy_weight <- 0.999 
coding_weight <- 0.995

# run matrix optimisation (note this takes a while), set seed for reproducibility
set.seed(8765)
logfile_gblogical_logodds_taxtrue_09990995_seed8765 <-
  densify_steps(original_data = logical_for_pruning,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight,
                coding_weight = coding_weight)
write.csv(logfile_gblogical_logodds_taxtrue_09990995_seed8765,"output/logicalGBI/logfile_gbilogical_logodds_taxtrue_09990995_seed8765.csv")

logfile_gbstatistical_logodds_taxtrue_09990995_seed8765 <-
  densify_steps(original_data = statistical_for_pruning,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight,
                coding_weight = coding_weight)
write.csv(logfile_gbstatistical_logodds_taxtrue_09990995_seed8765,"output/statisticalGBI/logfile_gbistatistical_logodds_taxtrue_09990995_seed8765.csv")

# set exponents for identifying optimal number of iterations
exponent_prop_coded_data = 1 # default
exponent_available_data_points = 1 # default
exponent_lowest_taxon_coding_score = 1 # since NAs on language end should largely be random
exponent_lowest_variable_coding_score = 0 # since part of the NAs on the variable end are explicitly wanted
exponent_taxonomic_diversity = 3 # since densification here explicitly seeks to increase taxonomic diversity

# identify optima
optimum_logical <- densify_score(logfile_gblogical_logodds_taxtrue_09990995_seed8765, 
                                exponent_prop_coded_data = exponent_prop_coded_data, 
                                exponent_available_data_points = exponent_available_data_points,
                                exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score, 
                                exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score, 
                                exponent_taxonomic_diversity = exponent_taxonomic_diversity, 
                                plot = TRUE)

optimum_statistical <- densify_score(logfile_gbstatistical_logodds_taxtrue_09990995_seed8765, 
                                    exponent_prop_coded_data = exponent_prop_coded_data, 
                                    exponent_available_data_points = exponent_available_data_points, 
                                    exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score, 
                                    exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score, 
                                    exponent_taxonomic_diversity = exponent_taxonomic_diversity,
                                    plot = TRUE)

# prune matrices to optima
pruned_logical <- densify_prune(logical_for_pruning, logfile_gblogical_logodds_taxtrue_09990995_seed8765, optimum_logical)
pruned_statistical <- densify_prune(statistical_for_pruning, logfile_gbstatistical_logodds_taxtrue_09990995_seed8765, optimum_statistical)

# save pruned matrices
write.csv(pruned_logical,"output/logicalGBI/logicalGBI_pruned.csv")
write.csv(pruned_statistical,"output/statisticalGBI/statisticalGBI_pruned.csv")

