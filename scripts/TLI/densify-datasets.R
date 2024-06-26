#### This script densifies the full TLI logical and statistical datasets, as well as their domain-wise subsets (grammar, phonology, lexicon)
rm(list=ls())

# load packages
library(densify)
library(tidyverse)
library(gmt)
library(ggplot2)

# read functions
source("scripts/functions.R")

# read in logical TLI data
logical <- read.csv("curated_data/TLI/logicalTLI/full/logicalTLI_full.csv") %>% select(-X)

# read in statistical TLI data
statistical <- read.csv("curated_data/TLI/statisticalTLI/full/statisticalTLI_full.csv") %>% select(-X)

# generate taxonomy matrix
taxonomy_matrix <- as_flat_taxonomy_matrix(glottolog_languoids)

# for densification, ensure all blanks, ? and "NA" are coded as NA
logical_for_pruning <- na_convert(logical)
statistical_for_pruning <- na_convert(statistical)

# group both logical and statistical data in 4 sets: (1) full, (2) grammar only (Grammar+Grammatical categories), (3) phonology only, (4) lexicon only (Lexical semantics + Lexical classes)
logical_parameters <- read.csv("curated_data/TLI/logicalTLI/cldf/parameters.csv")
logical_full_for_pruning <- logical_for_pruning
logical_grammar_for_pruning <- logical_for_pruning[which(apply(logical_for_pruning[which(names(logical_for_pruning)%in%filter(logical_parameters,grouping %in% c("Grammar_other", "Grammar_linear_order", "Grammatical_categories"))$short.name)],1,function(x)length(na.omit(x)))>0), c(1,which(names(logical_for_pruning)%in%filter(logical_parameters,grouping %in% c("Grammar_other", "Grammar_linear_order", "Grammatical_categories"))$short.name))]
logical_phonology_for_pruning <- logical_for_pruning[which(apply(logical_for_pruning[which(names(logical_for_pruning)%in%filter(logical_parameters,grouping %in% c("Phonology", "Phonology_prosodic", "Phonology_segmental_c", "Phonology_segmental_v"))$short.name)],1,function(x)length(na.omit(x)))>0),c(1,which(names(logical_for_pruning)%in%filter(logical_parameters,grouping %in% c("Phonology", "Phonology_prosodic", "Phonology_segmental_c", "Phonology_segmental_v"))$short.name))]
logical_lexicon_for_pruning <- logical_for_pruning[which(apply(logical_for_pruning[which(names(logical_for_pruning)%in%filter(logical_parameters,grouping %in% c("Lexical_semantics", "Lexical_categories"))$short.name)],1,function(x)length(na.omit(x)))>0),c(1,which(names(logical_for_pruning)%in%filter(logical_parameters,grouping %in% c("Lexical_semantics", "Lexical_categories"))$short.name))]

statistical_parameters <- read.csv("curated_data/TLI/statisticalTLI/cldf/parameters.csv")
statistical_full_for_pruning <- statistical_for_pruning
statistical_grammar_for_pruning <- statistical_for_pruning[which(apply(statistical_for_pruning[which(names(statistical_for_pruning)%in%filter(statistical_parameters,grouping %in% c("Grammar_other", "Grammar_linear_order", "Grammatical_categories"))$short.name)],1,function(x)length(na.omit(x)))>0),c(1,which(names(statistical_for_pruning)%in%filter(statistical_parameters,grouping %in% c("Grammar_other", "Grammar_linear_order", "Grammatical_categories"))$short.name))]
statistical_phonology_for_pruning <- statistical_for_pruning[which(apply(statistical_for_pruning[which(names(statistical_for_pruning)%in%filter(statistical_parameters,grouping %in% c("Phonology", "Phonology_prosodic", "Phonology_segmental_c", "Phonology_segmental_v"))$short.name)],1,function(x)length(na.omit(x)))>0),c(1,which(names(statistical_for_pruning)%in%filter(statistical_parameters,grouping %in% c("Phonology", "Phonology_prosodic", "Phonology_segmental_c", "Phonology_segmental_v"))$short.name))]
statistical_lexicon_for_pruning <- statistical_for_pruning[which(apply(statistical_for_pruning[which(names(statistical_for_pruning)%in%filter(statistical_parameters,grouping %in% c("Lexical_semantics", "Lexical_categories"))$short.name)],1,function(x)length(na.omit(x)))>0),c(1,which(names(statistical_for_pruning)%in%filter(statistical_parameters,grouping %in% c("Lexical_semantics", "Lexical_categories"))$short.name))]

summarize_matrix(logical_full_for_pruning)
summarize_matrix(logical_grammar_for_pruning)
summarize_matrix(logical_phonology_for_pruning)
summarize_matrix(logical_lexicon_for_pruning)

summarize_matrix(statistical_full_for_pruning)
summarize_matrix(statistical_grammar_for_pruning)
summarize_matrix(statistical_phonology_for_pruning)
summarize_matrix(statistical_lexicon_for_pruning)

# specify parameters for densification
min_variability <- 3 # each variable must have at least 3 languages in its second-largest state
density_mean <- "log_odds"

##### densify full curation
# comment on weights: TLI is very sparse in terms of coding density and in terms of taxonomic diversity
# we do not specify that taxonomy should be more important in trimming and provide equal weights for densify()

# run matrix optimisation, set a seed for reproducibility
set.seed(1111)
logical_full_log <-
  densify(data = logical_full_for_pruning,
          taxon_id = "glottocode",
          min_variability = min_variability,
          density_mean = density_mean,
          taxonomy = glottolog_languoids,
          density_mean_weights = list(coding = 0.999, taxonomy = 0.999))

set.seed(1111)
statistical_full_log <-
  densify(data = statistical_full_for_pruning,
          taxon_id = "glottocode",
          min_variability = min_variability,
          density_mean = density_mean,
          taxonomy = glottolog_languoids,
          density_mean_weights = list(coding = 0.999, taxonomy = 0.999))

# identify optima
densified_logical_full_large <- prune(logical_full_log,
                                   scoring_function = n_data_points * coding_density * taxonomic_index)

densified_logical_full_small <- prune(logical_full_log,
                                   scoring_function = n_data_points^3 * coding_density^3 * row_coding_density_min * taxonomic_index)

densified_statistical_full_large <- prune(statistical_full_log,
                                       scoring_function = n_data_points * coding_density * taxonomic_index)

densified_statistical_full_small <- prune(statistical_full_log,
                                       scoring_function = n_data_points^3 * coding_density^3 * row_coding_density_min * taxonomic_index)

# retrieve corresponding data from input (to re-establish differences between ? and NA)
densified_logical_full_large <- logical[which(logical$glottocode%in%densified_logical_full_large$glottocode), which(names(logical)%in%names(densified_logical_full_large))]
densified_logical_full_small <- logical[which(logical$glottocode%in%densified_logical_full_small$glottocode), which(names(logical)%in%names(densified_logical_full_small))]
densified_statistical_full_large <- statistical[which(statistical$glottocode%in%densified_statistical_full_large$glottocode), which(names(statistical)%in%names(densified_statistical_full_large))]
densified_statistical_full_small <- statistical[which(statistical$glottocode%in%densified_statistical_full_small$glottocode), which(names(statistical)%in%names(densified_statistical_full_small))]

# save densified matrices
write.csv(densified_logical_full_large,"curated_data/TLI/logicalTLI/full/logicalTLI_full_densified_large.csv")
write.csv(densified_logical_full_small,"curated_data/TLI/logicalTLI/full/logicalTLI_full_densified_small.csv")
write.csv(densified_statistical_full_large,"curated_data/TLI/statisticalTLI/full/statisticalTLI_full_densified_large.csv")
write.csv(densified_statistical_full_small,"curated_data/TLI/statisticalTLI/full/statisticalTLI_full_densified_small.csv")

# describe densified matrices and their relation to the full ones
densified_logical_full_large <- na_convert(densified_logical_full_large)
densified_logical_full_small <- na_convert(densified_logical_full_small)
summarize_matrix(densified_logical_full_large)
summarize_matrix(densified_logical_full_small)
summarize_matrix(densified_logical_full_large)/summarize_matrix(logical_full_for_pruning)
summarize_matrix(densified_logical_full_small)/summarize_matrix(logical_full_for_pruning)

densified_statistical_full_large <- na_convert(densified_statistical_full_large)
densified_statistical_full_small <- na_convert(densified_statistical_full_small)
summarize_matrix(densified_statistical_full_large)
summarize_matrix(densified_statistical_full_small)
summarize_matrix(densified_statistical_full_large)/summarize_matrix(statistical_full_for_pruning)
summarize_matrix(densified_statistical_full_small)/summarize_matrix(statistical_full_for_pruning)

##### densify grammar curation
# comment on weights: TLI grammar is very sparse in terms of coding density (much like the full dataset) and also uneven in terms of taxonomic diversity
# we do not specify that taxonomy should be more important in trimming and provide equal weights for densify()

# run matrix optimisation, set seed for reproducibility
set.seed(2222)
logical_grammar_log <-
  densify(data = logical_grammar_for_pruning,
          taxon_id = "glottocode",
          min_variability = min_variability,
          density_mean = density_mean,
          taxonomy = glottolog_languoids,
          density_mean_weights = list(coding = 0.999, taxonomy = 0.999))

set.seed(2222)
statistical_grammar_log <-
  densify(data = statistical_grammar_for_pruning,
          taxon_id = "glottocode",
          min_variability = min_variability,
          density_mean = density_mean,
          taxonomy = glottolog_languoids,
          density_mean_weights = list(coding = 0.999, taxonomy = 0.999))

# set scoring function for identifying optimal number of iterations
# we want to densify the matrix whilst maintaining taxonomic diversity and create two densifications, one being larger but sparser and one being smaller but denser
# --> n_data_points * coding_density * taxonomic_index
# --> n_data_points^3 * coding_density^3 * row_coding_density_min * taxonomic_index

# identify optima
densified_logical_grammar_large <- prune(logical_grammar_log,
                                            scoring_function = n_data_points * coding_density * taxonomic_index)

densified_logical_grammar_small <- prune(logical_grammar_log,
                                           scoring_function = n_data_points^3 * coding_density^3 * row_coding_density_min * taxonomic_index)

densified_statistical_grammar_large <- prune(statistical_grammar_log,
                                               scoring_function = n_data_points * coding_density * taxonomic_index)

densified_statistical_grammar_small <- prune(statistical_grammar_log,
                                               scoring_function = n_data_points^3 * coding_density^3 * row_coding_density_min * taxonomic_index)

# retrieve corresponding data from input (to re-establish differences between ? and NA)
densified_logical_grammar_large <- logical[which(logical$glottocode%in%densified_logical_grammar_large$glottocode), which(names(logical)%in%names(densified_logical_grammar_large))]
densified_logical_grammar_small <- logical[which(logical$glottocode%in%densified_logical_grammar_small$glottocode), which(names(logical)%in%names(densified_logical_grammar_small))]
densified_statistical_grammar_large <- statistical[which(statistical$glottocode%in%densified_statistical_grammar_large$glottocode), which(names(statistical)%in%names(densified_statistical_grammar_large))]
densified_statistical_grammar_small <- statistical[which(statistical$glottocode%in%densified_statistical_grammar_small$glottocode), which(names(statistical)%in%names(densified_statistical_grammar_small))]

# save densified matrices
write.csv(densified_logical_grammar_large,"curated_data/TLI/logicalTLI/grammar/logicalTLI_grammar_densified_large.csv")
write.csv(densified_logical_grammar_small,"curated_data/TLI/logicalTLI/grammar/logicalTLI_grammar_densified_small.csv")
write.csv(densified_statistical_grammar_large,"curated_data/TLI/statisticalTLI/grammar/statisticalTLI_grammar_densified_large.csv")
write.csv(densified_statistical_grammar_small,"curated_data/TLI/statisticalTLI/grammar/statisticalTLI_grammar_densified_small.csv")

# describe densified matrices and their relation to the full ones
densified_logical_grammar_large <- na_convert(densified_logical_grammar_large)
densified_logical_grammar_small <- na_convert(densified_logical_grammar_small)
summarize_matrix(densified_logical_grammar_large)
summarize_matrix(densified_logical_grammar_small)
summarize_matrix(densified_logical_grammar_large)/summarize_matrix(logical_grammar_for_pruning)
summarize_matrix(densified_logical_grammar_small)/summarize_matrix(logical_grammar_for_pruning)

densified_statistical_grammar_large <- na_convert(densified_statistical_grammar_large)
densified_statistical_grammar_small <- na_convert(densified_statistical_grammar_small)
summarize_matrix(densified_statistical_grammar_large)
summarize_matrix(densified_statistical_grammar_small)
summarize_matrix(densified_statistical_grammar_large)/summarize_matrix(statistical_grammar_for_pruning)
summarize_matrix(densified_statistical_grammar_small)/summarize_matrix(statistical_grammar_for_pruning)


##### densify lexicon curation
# comment on weights: TLI lexicon is quite sparse in terms of coding density (a bit less than the datasets above) and also uneven in terms of taxonomic diversity (but less strongly so than the full dataset).
# There are not many variables (just 30 (logical) or 28 (statistical))
# we do not specify that taxonomy should be more important for trimming and provide equal weights for densify()

# run matrix optimisation, set seed for reproducibility
set.seed(3333)
logical_lexicon_log <-
  densify(data = logical_lexicon_for_pruning,
          taxon_id = "glottocode",
          min_variability = min_variability,
          density_mean = density_mean,
          taxonomy = glottolog_languoids,
          density_mean_weights = list(coding = 0.999, taxonomy = 0.999))

set.seed(3333)
statistical_lexicon_log <-
  densify(data = statistical_lexicon_for_pruning,
          taxon_id = "glottocode",
          min_variability = min_variability,
          density_mean = density_mean,
          taxonomy = glottolog_languoids,
          density_mean_weights = list(coding = 0.999, taxonomy = 0.999))

# set scoring function for identifying optimal number of iterations
# we want to densify the matrix whilst maintaining taxonomic diversity; not many badly coded languages
# --> n_data_points * coding_density * taxonomic_index

# prune matrices to optima
densified_logical_lexicon <- prune(logical_lexicon_log, scoring_function = n_data_points * coding_density * taxonomic_index)
densified_statistical_lexicon <- prune(statistical_lexicon_log, scoring_function = n_data_points * coding_density * taxonomic_index)

# retrieve corresponding data from input (to re-establish differences between ? and NA)
densified_logical_lexicon <- logical[which(logical$glottocode%in%densified_logical_lexicon$glottocode), which(names(logical)%in%names(densified_logical_lexicon))]
densified_statistical_lexicon <- statistical[which(statistical$glottocode%in%densified_statistical_lexicon$glottocode), which(names(statistical)%in%names(densified_statistical_lexicon))]

# save densified matrices
write.csv(densified_logical_lexicon,"curated_data/TLI/logicalTLI/lexicon/logicalTLI_lexicon_densified.csv")
write.csv(densified_statistical_lexicon,"curated_data/TLI/statisticalTLI/lexicon/statisticalTLI_lexicon_densified.csv")

# describe densified matrices and their relation to the full ones
densified_logical_lexicon <- na_convert(densified_logical_lexicon)
summarize_matrix(densified_logical_lexicon)
summarize_matrix(densified_logical_lexicon)/summarize_matrix(logical_lexicon_for_pruning)

densified_statistical_lexicon <- na_convert(densified_statistical_lexicon)
summarize_matrix(densified_statistical_lexicon)
summarize_matrix(densified_statistical_lexicon)/summarize_matrix(statistical_lexicon_for_pruning)

##### densify phonology curation
# comment on weights: TLI phonology is rather dense (not particularly sparse in terms of coding density), but rather uneven in terms of taxonomic diversity
# still, we do not specify that taxonomy should be more important for trimming and provide equal weights for densify()

# run matrix optimisation, set seed for reproducibility
set.seed(4444)
logical_phonology_log <-
  densify(data = logical_phonology_for_pruning,
          taxon_id = "glottocode",
          min_variability = min_variability,
          density_mean = density_mean,
          taxonomy = glottolog_languoids,
          density_mean_weights = list(coding = 0.99, taxonomy = 0.99))

set.seed(4444)
statistical_phonology_log <-
  densify(data = statistical_phonology_for_pruning,
          taxon_id = "glottocode",
          min_variability = min_variability,
          density_mean = density_mean,
          taxonomy = glottolog_languoids,
          density_mean_weights = list(coding = 0.99, taxonomy = 0.99))

# set scoring function for identifying optimal number of iterations
# we want to densify the matrix whilst maintaining taxonomic diversity; the focus should strongly lie on trimming away languages from well-represented families
# --> n_data_points * coding_density * row_coding_density_min * taxonomic_index

densified_logical_phonology <- prune(logical_phonology_log, scoring_function = n_data_points * coding_density * row_coding_density_min * taxonomic_index)
densified_statistical_phonology <- prune(statistical_phonology_log, scoring_function = n_data_points * coding_density * row_coding_density_min * taxonomic_index)

# retrieve corresponding data from input (to re-establish differences between ? and NA)
densified_logical_phonology <- logical[which(logical$glottocode%in%densified_logical_phonology$glottocode), which(names(logical)%in%names(densified_logical_phonology))]
densified_statistical_phonology <- statistical[which(statistical$glottocode%in%densified_statistical_phonology$glottocode), which(names(statistical)%in%names(densified_statistical_phonology))]

# save densified matrices
write.csv(densified_logical_phonology,"curated_data/TLI/logicalTLI/phonology/logicalTLI_phonology_densified.csv")
write.csv(densified_statistical_phonology,"curated_data/TLI/statisticalTLI/phonology/statisticalTLI_phonology_densified.csv")

# describe densified matrices and their relation to the full ones
densified_logical_phonology <- na_convert(densified_logical_phonology)
summarize_matrix(densified_logical_phonology)
summarize_matrix(densified_logical_phonology)/summarize_matrix(logical_phonology_for_pruning)

densified_statistical_phonology <- na_convert(densified_statistical_phonology)
summarize_matrix(densified_statistical_phonology)
summarize_matrix(densified_statistical_phonology)/summarize_matrix(statistical_phonology_for_pruning)

