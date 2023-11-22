rm(list=ls())

# load packages
library(densify)
library(tidyverse)
library(gmt)
library(ggplot2)

# read functions
source("../functions.R")

# read in logical MM data
logical <- read.csv("output/logicalMM/full/logicalMM_full.csv") %>% select(-X)
lgs_lo <- logical$glottocode
logical <- logical %>% select(-glottocode)
rownames(logical) <- lgs_lo

# read in statistical MM data
statistical <- read.csv("output/statisticalMM/full/statisticalMM_full.csv") %>% select(-X)
lgs_st <- statistical$glottocode
statistical <- statistical %>% select(-glottocode)
rownames(statistical) <- lgs_st

# generate taxonomy matrix
taxonomy_matrix <- as_flat_taxonomy_matrix(glottolog_languoids)

# subset logical and statistical MM data to those languages with language coordinates (for plotting later)
macroareas <- read.csv("../../crossling-curated/input/glottolog_v.4.8/languages_and_dialects_geo.csv")
names(macroareas)[1]<-"id"
tax_with_coords <- merge(taxonomy_matrix,select(macroareas, c("id","longitude","latitude")),by="id")
tax_with_coords <- tax_with_coords %>% filter(!is.na(latitude))
logical_for_pruning <- logical[which(rownames(logical)%in%tax_with_coords$id),]
statistical_for_pruning <- statistical[which(rownames(statistical)%in%tax_with_coords$id),]

# for densification, ensure all blanks, ? and "NA" are coded as NA
logical_for_pruning <- na_convert(logical_for_pruning)
statistical_for_pruning <- na_convert(statistical_for_pruning)

# group both logical and statistical data in 4 sets: (1) full (grammar+phonology+lexicon), (2) grammar only, (3) phonology only, (4) lexicon only
logical_parameters <- read.csv("output/logicalMM/cldf/parameters.csv")
logical_full_for_pruning <- logical_for_pruning
logical_grammar_for_pruning <- logical_for_pruning[which(apply(logical_for_pruning[,which(names(logical_for_pruning)%in%filter(logical_parameters,domain == "Grammar")$new.name)],1,function(x)length(na.omit(x)))>0), which(names(logical_for_pruning)%in%filter(logical_parameters,domain == "Grammar")$new.name)]
logical_phonology_for_pruning <- logical_for_pruning[which(apply(logical_for_pruning[,which(names(logical_for_pruning)%in%filter(logical_parameters,domain == "Phonology")$new.name)],1,function(x)length(na.omit(x)))>0),which(names(logical_for_pruning)%in%filter(logical_parameters,domain == "Phonology")$new.name)]
logical_lexicon_for_pruning <- logical_for_pruning[which(apply(logical_for_pruning[,which(names(logical_for_pruning)%in%filter(logical_parameters,domain == "Lexicon")$new.name)],1,function(x)length(na.omit(x)))>0),which(names(logical_for_pruning)%in%filter(logical_parameters,domain == "Lexicon")$new.name)]

statistical_parameters <- read.csv("output/statisticalMM/cldf/parameters.csv")
statistical_full_for_pruning <- statistical_for_pruning
statistical_grammar_for_pruning <- statistical_for_pruning[which(apply(statistical_for_pruning[,which(names(statistical_for_pruning)%in%filter(statistical_parameters,domain == "Grammar")$new.name)],1,function(x)length(na.omit(x)))>0),which(names(statistical_for_pruning)%in%filter(statistical_parameters,domain == "Grammar")$new.name)]
statistical_phonology_for_pruning <- statistical_for_pruning[which(apply(statistical_for_pruning[,which(names(statistical_for_pruning)%in%filter(statistical_parameters,domain == "Phonology")$new.name)],1,function(x)length(na.omit(x)))>0),which(names(statistical_for_pruning)%in%filter(statistical_parameters,domain == "Phonology")$new.name)]
statistical_lexicon_for_pruning <- statistical_for_pruning[which(apply(statistical_for_pruning[,which(names(statistical_for_pruning)%in%filter(statistical_parameters,domain == "Lexicon")$new.name)],1,function(x)length(na.omit(x)))>0),which(names(statistical_for_pruning)%in%filter(statistical_parameters,domain == "Lexicon")$new.name)]

# specify parameters for densification 
variability_threshold <- 3 # each variable must have at least 3 languages in its second-largest state, for all densifications
mean_type <- "log_odds"

# maximal logically possible number of iterations as max_steps for each set
max_steps_logical_full <- nrow(logical_full_for_pruning)+ncol(logical_full_for_pruning)-2 
max_steps_statistical_full <- nrow(statistical_full_for_pruning)+ncol(statistical_full_for_pruning)-2

max_steps_logical_grammar <- nrow(logical_grammar_for_pruning)+ncol(logical_grammar_for_pruning)-2 
max_steps_statistical_grammar <- nrow(statistical_grammar_for_pruning)+ncol(statistical_grammar_for_pruning)-2 

max_steps_logical_phonology <- nrow(logical_phonology_for_pruning)+ncol(logical_phonology_for_pruning)-2 
max_steps_statistical_phonology <- nrow(statistical_phonology_for_pruning)+ncol(statistical_phonology_for_pruning)-2 

max_steps_logical_lexicon <- nrow(logical_lexicon_for_pruning)+ncol(logical_lexicon_for_pruning)-2 
max_steps_statistical_lexicon <- nrow(statistical_lexicon_for_pruning)+ncol(statistical_lexicon_for_pruning)-2 

##### densify full curation
# comment on weights: MM is very sparse in terms of coding density and in terms of taxonomic diversity
# we do not specify that taxonomy should be more important in trimming and provide equal weights for taxonomy_steps
taxonomy_weight_full <- 0.999
coding_weight_full <- 0.999

# run matrix optimisation (note this takes a while), set seed for reproducibility
set.seed(2222)
logfile_full_mmlogical_logodds_taxtrue_09990999_seed2222 <-
  densify_steps(original_data = logical_full_for_pruning,
                max_steps = max_steps_logical_full,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_full,
                coding_weight = coding_weight_full)
write.csv(logfile_full_mmlogical_logodds_taxtrue_09990999_seed2222,"output/logicalMM/full/logfile_full_mmlogical_logodds_taxtrue_09990999_seed2222.csv")

set.seed(2222)
logfile_full_mmstatistical_logodds_taxtrue_09990999_seed2222 <-
  densify_steps(original_data = statistical_full_for_pruning,
                max_steps = max_steps_statistical_full,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_full,
                coding_weight = coding_weight_full)
write.csv(logfile_full_mmstatistical_logodds_taxtrue_09990999_seed2222,"output/statisticalMM/full/logfile_full_mmstatistical_logodds_taxtrue_09990999_seed2222.csv")

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
optimum_logical_large <- densify_score(logfile_full_mmlogical_logodds_taxtrue_09990999_seed2222, 
                                       exponent_prop_coded_data = exponent_prop_coded_data_large, 
                                       exponent_available_data_points = exponent_available_data_points_large,
                                       exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_large, 
                                       exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_large, 
                                       exponent_taxonomic_diversity = exponent_taxonomic_diversity_large, 
                                       plot = TRUE)

optimum_logical_small <- densify_score(logfile_full_mmlogical_logodds_taxtrue_09990999_seed2222, 
                                       exponent_prop_coded_data = exponent_prop_coded_data_small, 
                                       exponent_available_data_points = exponent_available_data_points_small,
                                       exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_small, 
                                       exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_small, 
                                       exponent_taxonomic_diversity = exponent_taxonomic_diversity_small, 
                                       plot = TRUE)


optimum_statistical_large <- densify_score(logfile_full_mmstatistical_logodds_taxtrue_09990999_seed2222, 
                                           exponent_prop_coded_data = exponent_prop_coded_data_large, 
                                           exponent_available_data_points = exponent_available_data_points_large, 
                                           exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_large, 
                                           exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_large, 
                                           exponent_taxonomic_diversity = exponent_taxonomic_diversity_large,
                                           plot = TRUE)

optimum_statistical_small <- densify_score(logfile_full_mmstatistical_logodds_taxtrue_09990999_seed2222, 
                                           exponent_prop_coded_data = exponent_prop_coded_data_small, 
                                           exponent_available_data_points = exponent_available_data_points_small, 
                                           exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_small, 
                                           exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_small, 
                                           exponent_taxonomic_diversity = exponent_taxonomic_diversity_small,
                                           plot = TRUE)

# prune matrices to optima
pruned_logical_large <- densify_prune(logical_full_for_pruning, logfile_full_mmlogical_logodds_taxtrue_09990999_seed2222, optimum_logical_large)
pruned_logical_small <- densify_prune(logical_full_for_pruning, logfile_full_mmlogical_logodds_taxtrue_09990999_seed2222, optimum_logical_small)
pruned_statistical_large <- densify_prune(statistical_full_for_pruning, logfile_full_mmstatistical_logodds_taxtrue_09990999_seed2222, optimum_statistical_large)
pruned_statistical_small <- densify_prune(statistical_full_for_pruning, logfile_full_mmstatistical_logodds_taxtrue_09990999_seed2222, optimum_statistical_small)

# save pruned matrices
write.csv(pruned_logical_large,"output/logicalMM/full/logicalMM_full_pruned_large.csv")
write.csv(pruned_logical_small,"output/logicalMM/full/logicalMM_full_pruned_small.csv")
write.csv(pruned_statistical_large,"output/statisticalMM/full/statisticalMM_full_pruned_large.csv")
write.csv(pruned_statistical_small,"output/statisticalMM/full/statisticalMM_full_pruned_small.csv")

# compare datasets
comp_logical <- compare_full_to_densified(full=logical,
                                          densified=pruned_logical_large,
                                          densified2=pruned_logical_small,
                                          taxonomy_matrix=taxonomy_matrix,
                                          directory="output/logicalMM/full/")

comp_statistical <- compare_full_to_densified(full=statistical,
                                              densified=pruned_statistical_large,
                                              densified2=pruned_statistical_small,
                                              taxonomy_matrix=taxonomy_matrix,
                                              directory="output/statisticalMM/full/")

##### densify grammar curation
# comment on weights: MM grammar is very sparse in terms of coding density (much like the full dataset) and also uneven in terms of taxonomic diversity (but less strongly so than the full dataset)
# we do not specify that taxonomy should be more important in trimming and provide equal weights for taxonomy_steps
taxonomy_weight_grammar <- 0.999
coding_weight_grammar <- 0.999

# run matrix optimisation (note this takes a while), set seed for reproducibility
set.seed(2222)
logfile_grammar_mmlogical_logodds_taxtrue_09990999_seed2222 <-
  densify_steps(original_data = logical_grammar_for_pruning,
                max_steps = max_steps_logical_grammar,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_grammar,
                coding_weight = coding_weight_grammar)
write.csv(logfile_grammar_mmlogical_logodds_taxtrue_09990999_seed2222,"output/logicalMM/grammar/logfile_grammar_mmlogical_logodds_taxtrue_09990999_seed2222.csv")

set.seed(2222)
logfile_grammar_mmstatistical_logodds_taxtrue_09990999_seed2222 <-
  densify_steps(original_data = statistical_grammar_for_pruning,
                max_steps = max_steps_statistical_grammar,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_grammar,
                coding_weight = coding_weight_grammar)
write.csv(logfile_grammar_mmstatistical_logodds_taxtrue_09990999_seed2222,"output/statisticalMM/grammar/logfile_grammar_mmstatistical_logodds_taxtrue_09990999_seed2222.csv")

# set exponents for identifying optimal number of iterations
# we want to densify the matrix whilst maintaining taxonomic diversity
# A=B	This trade-off is not touched 
# E=C or E>C	Taxonomic diversity is key here too, so E can be > C. But because we have many badly coded languages, we can also set E=C
# D=0	Differential codedness in variables is partly induced by design, so D=0

## exponents for grammar, large
exponent_prop_coded_data_grammar_large = 1 # default
exponent_available_data_points_grammar_large = 1 # default
exponent_lowest_taxon_coding_score_grammar_large = 0
exponent_lowest_variable_coding_score_grammar_large = 0
exponent_taxonomic_diversity_grammar_large = 1 # included but not increased relative to first two exponents

## exponents for grammar, small
exponent_prop_coded_data_grammar_small = 3 # 
exponent_available_data_points_grammar_small = 3 # 
exponent_lowest_taxon_coding_score_grammar_small = 1
exponent_lowest_variable_coding_score_grammar_small = 0
exponent_taxonomic_diversity_grammar_small = 1 # 


# identify optima
optimum_logical_grammar_large <- densify_score(logfile_grammar_mmlogical_logodds_taxtrue_09990999_seed2222, 
                                       exponent_prop_coded_data = exponent_prop_coded_data_grammar_large, 
                                       exponent_available_data_points = exponent_available_data_points_grammar_large,
                                       exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_grammar_large, 
                                       exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_grammar_large, 
                                       exponent_taxonomic_diversity = exponent_taxonomic_diversity_grammar_large, 
                                       plot = TRUE)

optimum_logical_grammar_small <- densify_score(logfile_grammar_mmlogical_logodds_taxtrue_09990999_seed2222, 
                                               exponent_prop_coded_data = exponent_prop_coded_data_grammar_small, 
                                               exponent_available_data_points = exponent_available_data_points_grammar_small,
                                               exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_grammar_small, 
                                               exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_grammar_small, 
                                               exponent_taxonomic_diversity = exponent_taxonomic_diversity_grammar_small, 
                                               plot = TRUE)

optimum_statistical_grammar_large <- densify_score(logfile_grammar_mmstatistical_logodds_taxtrue_09990999_seed2222, 
                                           exponent_prop_coded_data = exponent_prop_coded_data_grammar_large, 
                                           exponent_available_data_points = exponent_available_data_points_grammar_large, 
                                           exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_grammar_large, 
                                           exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_grammar_large, 
                                           exponent_taxonomic_diversity = exponent_taxonomic_diversity_grammar_large,
                                           plot = TRUE)

optimum_statistical_grammar_small <- densify_score(logfile_grammar_mmstatistical_logodds_taxtrue_09990999_seed2222, 
                                                   exponent_prop_coded_data = exponent_prop_coded_data_grammar_small, 
                                                   exponent_available_data_points = exponent_available_data_points_grammar_small, 
                                                   exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_grammar_small, 
                                                   exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_grammar_small, 
                                                   exponent_taxonomic_diversity = exponent_taxonomic_diversity_grammar_small,
                                                   plot = TRUE)


# prune matrices to optima
pruned_logical_grammar_large <- densify_prune(logical_grammar_for_pruning, logfile_grammar_mmlogical_logodds_taxtrue_09990999_seed2222, optimum_logical_grammar_large)
pruned_logical_grammar_small <- densify_prune(logical_grammar_for_pruning, logfile_grammar_mmlogical_logodds_taxtrue_09990999_seed2222, optimum_logical_grammar_small)
pruned_statistical_grammar_large <- densify_prune(statistical_grammar_for_pruning, logfile_grammar_mmstatistical_logodds_taxtrue_09990999_seed2222, optimum_statistical_grammar_large)
pruned_statistical_grammar_small <- densify_prune(statistical_grammar_for_pruning, logfile_grammar_mmstatistical_logodds_taxtrue_09990999_seed2222, optimum_statistical_grammar_small)


# save pruned matrices
write.csv(pruned_logical_grammar_large,"output/logicalMM/grammar/logicalMM_grammar_pruned_large.csv")
write.csv(pruned_logical_grammar_small,"output/logicalMM/grammar/logicalMM_grammar_pruned_small.csv")
write.csv(pruned_statistical_grammar_large,"output/statisticalMM/grammar/statisticalMM_grammar_pruned_large.csv")
write.csv(pruned_statistical_grammar_small,"output/statisticalMM/grammar/statisticalMM_grammar_pruned_small.csv")

# compare datasets
comp_logical <- compare_full_to_densified(full=logical_grammar_for_pruning,
                                          densified=pruned_logical_grammar_large,
                                          densified2=pruned_logical_grammar_small,
                                          taxonomy_matrix=taxonomy_matrix,
                                          directory="output/logicalMM/grammar/")

comp_statistical <- compare_full_to_densified(full=statistical_grammar_for_pruning,
                                              densified=pruned_statistical_grammar_large,
                                              densified2=pruned_statistical_grammar_small,
                                              taxonomy_matrix=taxonomy_matrix,
                                              directory="output/statisticalMM/grammar/")

##### densify phonology curation
# comment on weights: MM phonology is rather dense (not particularly sparse in terms of coding density), but rather uneven in terms of taxonomic diversity
# still, we do not specify that taxonomy should be more important for trimming and provide equal weights for taxonomy_steps, without penalising 
taxonomy_weight_phonology <- 0.9
coding_weight_phonology <- 0.9

# run matrix optimisation (note this takes a while), set seed for reproducibility
set.seed(4444)
logfile_phonology_mmlogical_logodds_taxtrue_0909_seed4444 <-
  densify_steps(original_data = logical_phonology_for_pruning,
                max_steps = max_steps_logical_phonology,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_phonology,
                coding_weight = coding_weight_phonology)
write.csv(logfile_phonology_mmlogical_logodds_taxtrue_0909_seed4444,"output/logicalMM/phonology/logfile_phonology_mmlogical_logodds_taxtrue_0909_seed4444.csv")

set.seed(4444)
logfile_phonology_mmstatistical_logodds_taxtrue_0909_seed4444 <-
  densify_steps(original_data = statistical_phonology_for_pruning,
                max_steps = max_steps_statistical_phonology,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_phonology,
                coding_weight = coding_weight_phonology)
write.csv(logfile_phonology_mmstatistical_logodds_taxtrue_0909_seed4444,"output/statisticalMM/phonology/logfile_phonology_mmstatistical_logodds_taxtrue_0909_seed4444.csv")

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
optimum_logical_phonology <- densify_score(logfile_phonology_mmlogical_logodds_taxtrue_0909_seed4444, 
                                         exponent_prop_coded_data = exponent_prop_coded_data_phonology, 
                                         exponent_available_data_points = exponent_available_data_points_phonology,
                                         exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_phonology, 
                                         exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_phonology, 
                                         exponent_taxonomic_diversity = exponent_taxonomic_diversity_phonology, 
                                         plot = TRUE)

optimum_statistical_phonology <- densify_score(logfile_phonology_mmstatistical_logodds_taxtrue_0909_seed4444, 
                                             exponent_prop_coded_data = exponent_prop_coded_data_phonology, 
                                             exponent_available_data_points = exponent_available_data_points_phonology, 
                                             exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_phonology, 
                                             exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_phonology, 
                                             exponent_taxonomic_diversity = exponent_taxonomic_diversity_phonology,
                                             plot = TRUE)

# prune matrices to optima
pruned_logical_phonology <- densify_prune(logical_phonology_for_pruning, logfile_phonology_mmlogical_logodds_taxtrue_0909_seed4444, optimum_logical_phonology)
pruned_statistical_phonology <- densify_prune(statistical_phonology_for_pruning, logfile_phonology_mmstatistical_logodds_taxtrue_0909_seed4444, optimum_statistical_phonology)

# save pruned matrices
write.csv(pruned_logical_phonology,"output/logicalMM/phonology/logicalMM_phonology_pruned.csv")
write.csv(pruned_statistical_phonology,"output/statisticalMM/phonology/statisticalMM_phonology_pruned.csv")

# compare datasets
comp_logical <- compare_full_to_densified(full=logical_phonology_for_pruning,
                                          densified=pruned_logical_phonology,
                                          taxonomy_matrix=taxonomy_matrix,
                                          directory="output/logicalMM/phonology/")

comp_statistical <- compare_full_to_densified(full=statistical_phonology_for_pruning,
                                              densified=pruned_statistical_phonology,
                                              taxonomy_matrix=taxonomy_matrix,
                                              directory="output/statisticalMM/phonology/")


##### densify lexicon curation
# comment on weights: MM lexicon is quite sparse in terms of coding density (a bit less than the full dataset) and also uneven in terms of taxonomic diversity (but less strongly so than the full dataset). There are hardly any variables (just 21)
# we do not specify that taxonomy should be more important for trimming and provide equal weights for taxonomy_steps
taxonomy_weight_lexicon <- 0.99
coding_weight_lexicon <- 0.99

# run matrix optimisation (note this takes a while), set seed for reproducibility
set.seed(3333)
logfile_lexicon_mmlogical_logodds_taxtrue_099099_seed3333 <-
  densify_steps(original_data = logical_lexicon_for_pruning,
                max_steps = max_steps_logical_lexicon,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_lexicon,
                coding_weight = coding_weight_lexicon)
write.csv(logfile_lexicon_mmlogical_logodds_taxtrue_099099_seed3333,"output/logicalMM/lexicon/logfile_lexicon_mmlogical_logodds_taxtrue_099099_seed3333.csv")

set.seed(3333)
logfile_lexicon_mmstatistical_logodds_taxtrue_099099_seed3333 <-
  densify_steps(original_data = statistical_lexicon_for_pruning,
                max_steps = max_steps_statistical_lexicon,
                variability_threshold = variability_threshold,
                mean_type = mean_type,
                use_taxonomy = T,
                taxonomy = glottolog_languoids,
                taxonomy_weight = taxonomy_weight_lexicon,
                coding_weight = coding_weight_lexicon)
write.csv(logfile_lexicon_mmstatistical_logodds_taxtrue_099099_seed3333,"output/statisticalMM/lexicon/logfile_lexicon_mmstatistical_logodds_taxtrue_099099_seed3333.csv")

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
optimum_logical_lexicon <- densify_score(logfile_lexicon_mmlogical_logodds_taxtrue_099099_seed3333, 
                                           exponent_prop_coded_data = exponent_prop_coded_data_lexicon, 
                                           exponent_available_data_points = exponent_available_data_points_lexicon,
                                           exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_lexicon, 
                                           exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_lexicon, 
                                           exponent_taxonomic_diversity = exponent_taxonomic_diversity_lexicon, 
                                           plot = TRUE)

optimum_statistical_lexicon <- densify_score(logfile_lexicon_mmstatistical_logodds_taxtrue_099099_seed3333, 
                                               exponent_prop_coded_data = exponent_prop_coded_data_lexicon, 
                                               exponent_available_data_points = exponent_available_data_points_lexicon, 
                                               exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score_lexicon, 
                                               exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score_lexicon, 
                                               exponent_taxonomic_diversity = exponent_taxonomic_diversity_lexicon,
                                               plot = TRUE)


# prune matrices to optima
pruned_logical_lexicon <- densify_prune(logical_lexicon_for_pruning, logfile_lexicon_mmlogical_logodds_taxtrue_099099_seed3333, optimum_logical_lexicon)
pruned_statistical_lexicon <- densify_prune(statistical_lexicon_for_pruning, logfile_lexicon_mmstatistical_logodds_taxtrue_099099_seed3333, optimum_statistical_lexicon)

# save pruned matrices
write.csv(pruned_logical_lexicon,"output/logicalMM/lexicon/logicalMM_lexicon_pruned.csv")
write.csv(pruned_statistical_lexicon,"output/statisticalMM/lexicon/statisticalMM_lexicon_pruned.csv")

# compare datasets
comp_logical <- compare_full_to_densified(full=logical_lexicon_for_pruning,
                                          densified=pruned_logical_lexicon,
                                          taxonomy_matrix=taxonomy_matrix,
                                          directory="output/logicalMM/lexicon/")

comp_statistical <- compare_full_to_densified(full=statistical_lexicon_for_pruning,
                                              densified=pruned_statistical_lexicon,
                                              taxonomy_matrix=taxonomy_matrix,
                                              directory="output/statisticalMM/lexicon/")


