#### This script densifies the full TypLinkInd logical and statistical datasets, as well as their domain-wise subsets (morphosyntax, phonology, lexicon)

rm(list=ls())

# load packages
library(densify)
library(tidyverse)
library(gmt)
library(ggplot2)

# read functions
source("../functions.R")

# read in logical TLI data
logical <- read.csv("output/logicalTLI/full/logicalTLI_full.csv") %>% select(-X)
lgs_lo <- logical$glottocode
logical <- logical %>% select(-glottocode)
rownames(logical) <- lgs_lo
names(logical) <- str_replace_all(names(logical),"\\.","\\+")

# read in statistical TLI data
statistical <- read.csv("output/statisticalTLI/full/statisticalTLI_full.csv") %>% select(-X)
lgs_st <- statistical$glottocode
statistical <- statistical %>% select(-glottocode)
rownames(statistical) <- lgs_st
names(statistical) <- str_replace_all(names(statistical),"\\.","\\+")

# generate taxonomy matrix
taxonomy_matrix <- as_flat_taxonomy_matrix(glottolog_languoids)

# subset logical and statistical TLI data to those languages with language coordinates (for plotting later)
macroareas <- read.csv("../../raw/glottolog_v.4.8/languages_and_dialects_geo.csv")
names(macroareas)[1]<-"id"
tax_with_coords <- merge(taxonomy_matrix,select(macroareas, c("id","longitude","latitude")),by="id")
tax_with_coords <- tax_with_coords %>% filter(!is.na(latitude))

# for densification, ensure all blanks, ? and "NA" are coded as NA
logical_for_pruning <- na_convert(logical)
statistical_for_pruning <- na_convert(statistical)

# group both logical and statistical data in 4 sets: (1) full (morphosyntax+phonology+lexicon), (2) morphosyntax only, (3) phonology only, (4) lexicon only
logical_parameters <- read.csv("output/logicalTLI/cldf/parameters.csv")
logical_full_for_pruning <- logical_for_pruning
logical_morphosyntax_for_pruning <- logical_for_pruning[which(apply(logical_for_pruning[,which(names(logical_for_pruning)%in%filter(logical_parameters,domain == "Morphosyntax")$new.name)],1,function(x)length(na.omit(x)))>0), which(names(logical_for_pruning)%in%filter(logical_parameters,domain == "Morphosyntax")$new.name)]
logical_phonology_for_pruning <- logical_for_pruning[which(apply(logical_for_pruning[,which(names(logical_for_pruning)%in%filter(logical_parameters,domain == "Phonology")$new.name)],1,function(x)length(na.omit(x)))>0),which(names(logical_for_pruning)%in%filter(logical_parameters,domain == "Phonology")$new.name)]
logical_lexicon_for_pruning <- logical_for_pruning[which(apply(logical_for_pruning[,which(names(logical_for_pruning)%in%filter(logical_parameters,domain == "Lexicon")$new.name)],1,function(x)length(na.omit(x)))>0),which(names(logical_for_pruning)%in%filter(logical_parameters,domain == "Lexicon")$new.name)]

statistical_parameters <- read.csv("output/statisticalTLI/cldf/parameters.csv")
statistical_full_for_pruning <- statistical_for_pruning
statistical_morphosyntax_for_pruning <- statistical_for_pruning[which(apply(statistical_for_pruning[,which(names(statistical_for_pruning)%in%filter(statistical_parameters,domain == "Morphosyntax")$new.name)],1,function(x)length(na.omit(x)))>0),which(names(statistical_for_pruning)%in%filter(statistical_parameters,domain == "Morphosyntax")$new.name)]
statistical_phonology_for_pruning <- statistical_for_pruning[which(apply(statistical_for_pruning[,which(names(statistical_for_pruning)%in%filter(statistical_parameters,domain == "Phonology")$new.name)],1,function(x)length(na.omit(x)))>0),which(names(statistical_for_pruning)%in%filter(statistical_parameters,domain == "Phonology")$new.name)]
statistical_lexicon_for_pruning <- statistical_for_pruning[which(apply(statistical_for_pruning[,which(names(statistical_for_pruning)%in%filter(statistical_parameters,domain == "Lexicon")$new.name)],1,function(x)length(na.omit(x)))>0),which(names(statistical_for_pruning)%in%filter(statistical_parameters,domain == "Lexicon")$new.name)]

# specify parameters for densification 
variability_threshold <- 3 # each variable must have at least 3 languages in its second-largest state, for all densifications
mean_type <- "log_odds"

##### densify full curation
# comment on weights: TLI is very sparse in terms of coding density and in terms of taxonomic diversity
# we do not specify that taxonomy should be more important in trimming and provide equal weights for taxonomy_steps
taxonomy_weight_full <- 0.999
coding_weight_full <- 0.999

# run matrix optimisation (note this takes a while), set a seed for reproducibility
set.seed(2222)
logfile_full_tlilogical_logodds_taxtrue_09990999_seed2222 <-
  densify_steps(original_data = logical_full_for_pruning,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_full,
                coding_weight = coding_weight_full)
write.csv(logfile_full_tlilogical_logodds_taxtrue_09990999_seed2222,"output/logicalTLI/full/logfile_full_tlilogical_logodds_taxtrue_09990999_seed2222.csv")

set.seed(2222)
logfile_full_tlistatistical_logodds_taxtrue_09990999_seed2222 <-
  densify_steps(original_data = statistical_full_for_pruning,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_full,
                coding_weight = coding_weight_full)
write.csv(logfile_full_tlistatistical_logodds_taxtrue_09990999_seed2222,"output/statisticalTLI/full/logfile_full_tlistatistical_logodds_taxtrue_09990999_seed2222.csv")

# set exponents for identifying optimal number of iterations
# we want to densify the matrix whilst maintaining taxonomic diversity
# A=B	This trade-off is not touched 
# E=C or E>C	Taxonomic diversity is key here too, so E can be > C. But because we have many badly coded languages, we can also set E=C
# D=0	Differential codedness in variables is partly induced by design, so D=0

## first set of exponents, resulting in rather large matrix
exponent_prop_coded_data_large = 1 # default
exponent_available_data_points_large = 1 # default
exponent_lowest_taxon_coding_score_large = 0
exponent_lowest_variable_coding_score_large = 0
exponent_taxonomic_diversity_large = 1 # included but not increased relative to first two exponents

## second set of exponents, resulting in rather small matrix
exponent_prop_coded_data_small = 3 # increased for denser matrix
exponent_available_data_points_small = 3 # increased for denser matrix
exponent_lowest_taxon_coding_score_small = 1 # included for less badly-coded languages, denser matrix
exponent_lowest_variable_coding_score_small = 0
exponent_taxonomic_diversity_small = 1 # included and decreased relative to first two exponents 

# identify optima
optimum_logical_large <- densify_score(logfile_full_tlilogical_logodds_taxtrue_09990999_seed0000, 
                                       exponent_prop_coded_data = exponent_prop_coded_data_large, 
                                       exponent_available_data_points = exponent_available_data_points_large,
                                       exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_large, 
                                       exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_large, 
                                       exponent_taxonomic_diversity = exponent_taxonomic_diversity_large, 
                                       plot = TRUE)

optimum_logical_small <- densify_score(logfile_full_tlilogical_logodds_taxtrue_09990999_seed0000, 
                                       exponent_prop_coded_data = exponent_prop_coded_data_small, 
                                       exponent_available_data_points = exponent_available_data_points_small,
                                       exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_small, 
                                       exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_small, 
                                       exponent_taxonomic_diversity = exponent_taxonomic_diversity_small, 
                                       plot = TRUE)


optimum_statistical_large <- densify_score(logfile_full_tlistatistical_logodds_taxtrue_09990999_seed0000, 
                                           exponent_prop_coded_data = exponent_prop_coded_data_large, 
                                           exponent_available_data_points = exponent_available_data_points_large, 
                                           exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_large, 
                                           exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_large, 
                                           exponent_taxonomic_diversity = exponent_taxonomic_diversity_large,
                                           plot = TRUE)

optimum_statistical_small <- densify_score(logfile_full_tlistatistical_logodds_taxtrue_09990999_seed0000, 
                                           exponent_prop_coded_data = exponent_prop_coded_data_small, 
                                           exponent_available_data_points = exponent_available_data_points_small, 
                                           exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_small, 
                                           exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_small, 
                                           exponent_taxonomic_diversity = exponent_taxonomic_diversity_small,
                                           plot = TRUE)

# prune matrices to optima
pruned_logical_large <- densify_prune(logical_full_for_pruning, logfile_full_tlilogical_logodds_taxtrue_09990999_seed0000, optimum_logical_large)
pruned_logical_small <- densify_prune(logical_full_for_pruning, logfile_full_tlilogical_logodds_taxtrue_09990999_seed0000, optimum_logical_small)
pruned_statistical_large <- densify_prune(statistical_full_for_pruning, logfile_full_tlistatistical_logodds_taxtrue_09990999_seed0000, optimum_statistical_large)
pruned_statistical_small <- densify_prune(statistical_full_for_pruning, logfile_full_tlistatistical_logodds_taxtrue_09990999_seed0000, optimum_statistical_small)

# save pruned matrices
write.csv(pruned_logical_large,"output/logicalTLI/full/logicalTLI_full_pruned_large.csv")
write.csv(pruned_logical_small,"output/logicalTLI/full/logicalTLI_full_pruned_small.csv")
write.csv(pruned_statistical_large,"output/statisticalTLI/full/statisticalTLI_full_pruned_large.csv")
write.csv(pruned_statistical_small,"output/statisticalTLI/full/statisticalTLI_full_pruned_small.csv")

##### densify morphosyntax curation
# comment on weights: TLI morphosyntax is very sparse in terms of coding density (much like the full dataset) and also uneven in terms of taxonomic diversity
# we do not specify that taxonomy should be more important in trimming and provide equal weights for taxonomy_steps
taxonomy_weight_morphosyntax <- 0.9999
coding_weight_morphosyntax <- 0.9999

# run matrix optimisation (note this takes a while), set seed for reproducibility
set.seed(1111)
logfile_morphosyntax_tlilogical_logodds_taxtrue_0999909999_seed1111 <-
  densify_steps(original_data = logical_morphosyntax_for_pruning,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_morphosyntax,
                coding_weight = coding_weight_morphosyntax)
write.csv(logfile_morphosyntax_tlilogical_logodds_taxtrue_0999909999_seed1111,"output/logicalTLI/morphosyntax/logfile_morphosyntax_tlilogical_logodds_taxtrue_0999909999_seed1111.csv")

set.seed(1111)
logfile_morphosyntax_tlistatistical_logodds_taxtrue_0999909999_seed1111 <-
  densify_steps(original_data = statistical_morphosyntax_for_pruning,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_morphosyntax,
                coding_weight = coding_weight_morphosyntax)
write.csv(logfile_morphosyntax_tlistatistical_logodds_taxtrue_0999909999_seed1111,"output/statisticalTLI/morphosyntax/logfile_morphosyntax_tlistatistical_logodds_taxtrue_0999909999_seed1111.csv")

# set exponents for identifying optimal number of iterations
# we want to densify the matrix whilst maintaining taxonomic diversity
# A=B	This trade-off is not touched 
# E=C or E>C	Taxonomic diversity is key here too, so E can be > C. But because we have many badly coded languages, we can also set E=C
# D=0	Differential codedness in variables is partly induced by design, so D=0

## exponents for morphosyntax, large
exponent_prop_coded_data_morphosyntax_large = 1 # default
exponent_available_data_points_morphosyntax_large = 1 # default
exponent_lowest_taxon_coding_score_morphosyntax_large = 0
exponent_lowest_variable_coding_score_morphosyntax_large = 0
exponent_taxonomic_diversity_morphosyntax_large = 1 # included but not increased relative to first two exponents

## exponents for morphosyntax, small
exponent_prop_coded_data_morphosyntax_small = 3 # 
exponent_available_data_points_morphosyntax_small = 3 # 
exponent_lowest_taxon_coding_score_morphosyntax_small = 1
exponent_lowest_variable_coding_score_morphosyntax_small = 0
exponent_taxonomic_diversity_morphosyntax_small = 1 # 

# identify optima
optimum_logical_morphosyntax_large <- densify_score(logfile_morphosyntax_tlilogical_logodds_taxtrue_0999909999_seed1111, 
                                       exponent_prop_coded_data = exponent_prop_coded_data_morphosyntax_large, 
                                       exponent_available_data_points = exponent_available_data_points_morphosyntax_large,
                                       exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_morphosyntax_large, 
                                       exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_morphosyntax_large, 
                                       exponent_taxonomic_diversity = exponent_taxonomic_diversity_morphosyntax_large, 
                                       plot = TRUE)

optimum_logical_morphosyntax_small <- densify_score(logfile_morphosyntax_tlistatistical_logodds_taxtrue_0999909999_seed1111, 
                                               exponent_prop_coded_data = exponent_prop_coded_data_morphosyntax_small, 
                                               exponent_available_data_points = exponent_available_data_points_morphosyntax_small,
                                               exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_morphosyntax_small, 
                                               exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_morphosyntax_small, 
                                               exponent_taxonomic_diversity = exponent_taxonomic_diversity_morphosyntax_small, 
                                               plot = TRUE)

optimum_statistical_morphosyntax_large <- densify_score(logfile_morphosyntax_tlilogical_logodds_taxtrue_0999909999_seed1111, 
                                           exponent_prop_coded_data = exponent_prop_coded_data_morphosyntax_large, 
                                           exponent_available_data_points = exponent_available_data_points_morphosyntax_large, 
                                           exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_morphosyntax_large, 
                                           exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_morphosyntax_large, 
                                           exponent_taxonomic_diversity = exponent_taxonomic_diversity_morphosyntax_large,
                                           plot = TRUE)

optimum_statistical_morphosyntax_small <- densify_score(logfile_morphosyntax_tlistatistical_logodds_taxtrue_0999909999_seed1111, 
                                                   exponent_prop_coded_data = exponent_prop_coded_data_morphosyntax_small, 
                                                   exponent_available_data_points = exponent_available_data_points_morphosyntax_small, 
                                                   exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_morphosyntax_small, 
                                                   exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_morphosyntax_small, 
                                                   exponent_taxonomic_diversity = exponent_taxonomic_diversity_morphosyntax_small,
                                                   plot = TRUE)


# prune matrices to optima
pruned_logical_morphosyntax_large <- densify_prune(logical_morphosyntax_for_pruning, logfile_morphosyntax_tlilogical_logodds_taxtrue_0999909999_seed1111, optimum_logical_morphosyntax_large)
pruned_logical_morphosyntax_small <- densify_prune(logical_morphosyntax_for_pruning, logfile_morphosyntax_tlilogical_logodds_taxtrue_0999909999_seed1111, optimum_logical_morphosyntax_small)
pruned_statistical_morphosyntax_large <- densify_prune(statistical_morphosyntax_for_pruning, logfile_morphosyntax_tlistatistical_logodds_taxtrue_0999909999_seed1111, optimum_statistical_morphosyntax_large)
pruned_statistical_morphosyntax_small <- densify_prune(statistical_morphosyntax_for_pruning, logfile_morphosyntax_tlistatistical_logodds_taxtrue_0999909999_seed1111, optimum_statistical_morphosyntax_small)


# save pruned matrices
write.csv(pruned_logical_morphosyntax_large,"output/logicalTLI/morphosyntax/logicalTLI_morphosyntax_pruned_large.csv")
write.csv(pruned_logical_morphosyntax_small,"output/logicalTLI/morphosyntax/logicalTLI_morphosyntax_pruned_small.csv")
write.csv(pruned_statistical_morphosyntax_large,"output/statisticalTLI/morphosyntax/statisticalTLI_morphosyntax_pruned_large.csv")
write.csv(pruned_statistical_morphosyntax_small,"output/statisticalTLI/morphosyntax/statisticalTLI_morphosyntax_pruned_small.csv")

##### densify lexicon curation
# comment on weights: TLI lexicon is quite sparse in terms of coding density (a bit less than the datasets above) and also uneven in terms of taxonomic diversity (but less strongly so than the full dataset). 
# There are not many variables (just 27 (logical) or 25 (statistical))
# we do not specify that taxonomy should be more important for trimming and provide equal weights for taxonomy_steps
taxonomy_weight_lexicon <- 0.9999
coding_weight_lexicon <- 0.9999

# run matrix optimisation (note this takes a while), set seed for reproducibility
set.seed(2222)
logfile_lexicon_tlilogical_logodds_taxtrue_0999909999_seed2222 <-
  densify_steps(original_data = logical_lexicon_for_pruning,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_lexicon,
                coding_weight = coding_weight_lexicon)
write.csv(logfile_lexicon_tlilogical_logodds_taxtrue_0999909999_seed2222,"output/logicalTLI/lexicon/logfile_lexicon_tlilogical_logodds_taxtrue_0999909999_seed2222.csv")

set.seed(2222)
logfile_lexicon_tlistatistical_logodds_taxtrue_0999909999_seed2222 <-
  densify_steps(original_data = statistical_lexicon_for_pruning,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_lexicon,
                coding_weight = coding_weight_lexicon)
write.csv(logfile_lexicon_tlistatistical_logodds_taxtrue_0999909999_seed2222,"output/statisticalTLI/lexicon/logfile_lexicon_tlistatistical_logodds_taxtrue_0999909999_seed2222.csv")

# set exponents for identifying optimal number of iterations
# we want to densify the matrix whilst maintaining taxonomic diversity
# A=B	This trade-off is not touched 
# E>C	Taxonomic diversity is key here too, and there are not many badly coded languages, so E > C.
# D=0	Differential codedness in variables is partly induced by design, so D=0

## exponents for lexicon
exponent_prop_coded_data_lexicon = 1 # default
exponent_available_data_points_lexicon = 1 # default
exponent_lowest_taxon_coding_score_lexicon = 0
exponent_lowest_variable_coding_score_lexicon = 0
exponent_taxonomic_diversity_lexicon = 1 # included but not increased relative to first two exponents

# identify optima
optimum_logical_lexicon <- densify_score(logfile_lexicon_tlilogical_logodds_taxtrue_0999909999_seed2222, 
                                           exponent_prop_coded_data = exponent_prop_coded_data_lexicon, 
                                           exponent_available_data_points = exponent_available_data_points_lexicon,
                                           exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_lexicon, 
                                           exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_lexicon, 
                                           exponent_taxonomic_diversity = exponent_taxonomic_diversity_lexicon, 
                                           plot = TRUE)

optimum_statistical_lexicon <- densify_score(logfile_lexicon_tlistatistical_logodds_taxtrue_0999909999_seed2222, 
                                               exponent_prop_coded_data = exponent_prop_coded_data_lexicon, 
                                               exponent_available_data_points = exponent_available_data_points_lexicon, 
                                               exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_lexicon, 
                                               exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_lexicon, 
                                               exponent_taxonomic_diversity = exponent_taxonomic_diversity_lexicon,
                                               plot = TRUE)


# prune matrices to optima
pruned_logical_lexicon <- densify_prune(logical_lexicon_for_pruning, logfile_lexicon_tlilogical_logodds_taxtrue_0999909999_seed2222, optimum_logical_lexicon)
pruned_statistical_lexicon <- densify_prune(statistical_lexicon_for_pruning, logfile_lexicon_tlistatistical_logodds_taxtrue_0999909999_seed2222, optimum_statistical_lexicon)

# save pruned matrices
write.csv(pruned_logical_lexicon,"output/logicalTLI/lexicon/logicalTLI_lexicon_pruned.csv")
write.csv(pruned_statistical_lexicon,"output/statisticalTLI/lexicon/statisticalTLI_lexicon_pruned.csv")



##### densify phonology curation
# comment on weights: TLI phonology is rather dense (not particularly sparse in terms of coding density), but rather uneven in terms of taxonomic diversity
# still, we do not specify that taxonomy should be more important for trimming and provide equal weights for taxonomy_steps, without penalising 
taxonomy_weight_phonology <- 0.9
coding_weight_phonology <- 0.9

# run matrix optimisation (note this takes a while), set seed for reproducibility
set.seed(3333)
logfile_phonology_tlilogical_logodds_taxtrue_0909_seed3333 <-
  densify_steps(original_data = logical_phonology_for_pruning,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_phonology,
                coding_weight = coding_weight_phonology)
write.csv(logfile_phonology_tlilogical_logodds_taxtrue_0909_seed3333,"output/logicalTLI/phonology/logfile_phonology_tlilogical_logodds_taxtrue_0909_seed3333.csv")

set.seed(3333)
logfile_phonology_tlistatistical_logodds_taxtrue_0909_seed3333 <-
  densify_steps(original_data = statistical_phonology_for_pruning,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_phonology,
                coding_weight = coding_weight_phonology)
write.csv(logfile_phonology_tlistatistical_logodds_taxtrue_0909_seed3333,"output/statisticalTLI/phonology/logfile_phonology_tlistatistical_logodds_taxtrue_0909_seed3333.csv")

# set exponents for identifying optimal number of iterations
# we want to densify the matrix whilst maintaining taxonomic diversity
# A=B	This trade-off is not touched 
# E>C	Taxonomic diversity is key here too, and there are not many badly coded languages, so E > C.
# D=0	Differential codedness in variables is partly induced by design, so D=0

## exponents for phonology
exponent_prop_coded_data_phonology = 1 # default
exponent_available_data_points_phonology = 1 # default
exponent_lowest_taxon_coding_score_phonology = 0
exponent_lowest_variable_coding_score_phonology = 0
exponent_taxonomic_diversity_phonology = 1 # included but not increased relative to first two exponents

# identify optima
optimum_logical_phonology <- densify_score(logfile_phonology_tlilogical_logodds_taxtrue_0909_seed3333, 
                                           exponent_prop_coded_data = exponent_prop_coded_data_phonology, 
                                           exponent_available_data_points = exponent_available_data_points_phonology,
                                           exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_phonology, 
                                           exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_phonology, 
                                           exponent_taxonomic_diversity = exponent_taxonomic_diversity_phonology, 
                                           plot = TRUE)

optimum_statistical_phonology <- densify_score(logfile_phonology_tlistatistical_logodds_taxtrue_0909_seed3333, 
                                               exponent_prop_coded_data = exponent_prop_coded_data_phonology, 
                                               exponent_available_data_points = exponent_available_data_points_phonology, 
                                               exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_phonology, 
                                               exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_phonology, 
                                               exponent_taxonomic_diversity = exponent_taxonomic_diversity_phonology,
                                               plot = TRUE)

# prune matrices to optima
pruned_logical_phonology <- densify_prune(logical_phonology_for_pruning, logfile_phonology_tlilogical_logodds_taxtrue_0909_seed3333, optimum_logical_phonology)
pruned_statistical_phonology <- densify_prune(statistical_phonology_for_pruning, logfile_phonology_tlistatistical_logodds_taxtrue_0909_seed3333, optimum_statistical_phonology)

# save pruned matrices
write.csv(pruned_logical_phonology,"output/logicalTLI/phonology/logicalTLI_phonology_pruned.csv")
write.csv(pruned_statistical_phonology,"output/statisticalTLI/phonology/statisticalTLI_phonology_pruned.csv")
