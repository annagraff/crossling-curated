# compare raw to densified datasets
library(tidyverse)
source("../functions.R")

# GBInd
gbi_lf <- read.csv("../data/GBInd/output/logicalGBI/logicalGBI.csv", row.names = "glottocode") %>% select(-X)
gbi_lp <- read.csv("../data/GBInd/output/logicalGBI/logicalGBI_pruned.csv", row.names = "glottocode") %>% select(-X)
gbi_sf <- read.csv("../data/GBInd/output/statisticalGBI/statisticalGBI.csv", row.names = "glottocode") %>% select(-X)
gbi_sp <- read.csv("../data/GBInd/output/statisticalGBI/statisticalGBI_pruned.csv", row.names = "glottocode") %>% select(-X)


gbi_comp_logical <- compare_full_to_densified(full = logical_for_pruning,
                                          densified = pruned_logical,
                                          taxonomy_matrix = taxonomy_matrix)

gbi_comp_statistical <- compare_full_to_densified(full = statistical_for_pruning,
                                              densified = pruned_statistical,
                                              taxonomy_matrix = taxonomy_matrix)

