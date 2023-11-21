rm(list=ls())
library(densify)
library(tidyverse)
library(cluster)
library(pcaMethods)
library(gridExtra)
library(rnaturalearth)
library(lattice)

source("../functions.R")

# generate taxonomy
taxonomy_matrix <- as_flat_taxonomy_matrix(glottolog_languoids)
# link coordinates to lgs
macroareas <- read.csv("../input/glottolog_v.4.8/languages_and_dialects_geo.csv")
names(macroareas)[1] <- "id"
taxonomy_matrix_for_plot <- merge(taxonomy_matrix,select(macroareas, c("id","latitude","longitude")),by="id")
names(taxonomy_matrix_for_plot)[29:30] <- c("lat","lon") # change names to lat and lon for downstream functions

# read in matrices, convert all NA and ? to NA, factorize all variables
# full
full_logical_large <- read.csv("output/logicalMM/full/logicalMM_full_pruned_large.csv",row.names = "X")
full_logical_large <- na_convert(full_logical_large)
full_logical_large <- factorise(full_logical_large)

full_logical_small <- read.csv("output/logicalMM/full/logicalMM_full_pruned_small.csv", row.names = "X")
full_logical_small <- na_convert(full_logical_small)
full_logical_small <- factorise(full_logical_small)

full_statistical_large <- read.csv("output/statisticalMM/full/statisticalMM_full_pruned_large.csv",row.names = "X")
full_statistical_large <- na_convert(full_statistical_large)
full_statistical_large <- factorise(full_statistical_large)

full_statistical_small <- read.csv("output/statisticalMM/full/statisticalMM_full_pruned_small.csv", row.names = "X")
full_statistical_small <- na_convert(full_statistical_small)
full_statistical_small <- factorise(full_statistical_small)

# grammar
grammar_logical_large <- read.csv("output/logicalMM/grammar/logicalMM_grammar_pruned_large.csv",row.names = "X")
grammar_logical_large <- na_convert(grammar_logical_large)
grammar_logical_large <- factorise(grammar_logical_large)

grammar_logical_small <- read.csv("output/logicalMM/grammar/logicalMM_grammar_pruned_small.csv", row.names = "X")
grammar_logical_small <- na_convert(grammar_logical_small)
grammar_logical_small <- factorise(grammar_logical_small)

grammar_statistical_large <- read.csv("output/statisticalMM/grammar/statisticalMM_grammar_pruned_large.csv",row.names = "X")
grammar_statistical_large <- na_convert(grammar_statistical_large)
grammar_statistical_large <- factorise(grammar_statistical_large)

grammar_statistical_small <- read.csv("output/statisticalMM/grammar/statisticalMM_grammar_pruned_small.csv", row.names = "X")
grammar_statistical_small <- na_convert(grammar_statistical_small)
grammar_statistical_small <- factorise(grammar_statistical_small)

# phonology
phonology_logical <- read.csv("output/logicalMM/phonology/logicalMM_phonology_pruned.csv", row.names = "X")
phonology_logical <- na_convert(phonology_logical)
phonology_logical <- factorise(phonology_logical)

phonology_statistical <- read.csv("output/statisticalMM/phonology/statisticalMM_phonology_pruned.csv",row.names = "X")
phonology_statistical <- na_convert(phonology_statistical)
phonology_statistical <- factorise(phonology_statistical)

# lexicon
lexicon_logical <- read.csv("output/logicalMM/lexicon/logicalMM_lexicon_pruned.csv", row.names = "X")
lexicon_logical <- na_convert(lexicon_logical)
lexicon_logical <- factorise(lexicon_logical)

lexicon_statistical <- read.csv("output/statisticalMM/lexicon/statisticalMM_lexicon_pruned.csv",row.names = "X")
lexicon_statistical <- na_convert(lexicon_statistical)
lexicon_statistical <- factorise(lexicon_statistical)

#### Gower distances
# compute gower distances, using daisy()
full_logical_large_gower <- daisy(full_logical_large, metric = "gower") # %>% as.matrix() %>% na.omit() %>% as.dist()
full_logical_small_gower <- daisy(full_logical_small, metric = "gower")
full_statistical_large_gower <- daisy(full_statistical_large, metric = "gower")
full_statistical_small_gower <- daisy(full_statistical_small, metric = "gower")

##################
## PCA, aggregated to family level
# generate per family proportion table
proportions_full_logical_large <- generate_per_family_prop_table(data = full_logical_large, taxonomy_matrix = taxonomy_matrix)
proportions_full_logical_small <- generate_per_family_prop_table(data = full_logical_small, taxonomy_matrix = taxonomy_matrix)
proportions_full_statistical_large <- generate_per_family_prop_table(data = full_statistical_large, taxonomy_matrix = taxonomy_matrix)
proportions_full_statistical_small <- generate_per_family_prop_table(data = full_statistical_small, taxonomy_matrix = taxonomy_matrix)
proportions_grammar_logical_large <- generate_per_family_prop_table(data = grammar_logical_large, taxonomy_matrix = taxonomy_matrix)
proportions_grammar_logical_small <- generate_per_family_prop_table(data = grammar_logical_small, taxonomy_matrix = taxonomy_matrix)
proportions_grammar_statistical_large <- generate_per_family_prop_table(data = grammar_statistical_large, taxonomy_matrix = taxonomy_matrix)
proportions_grammar_statistical_small <- generate_per_family_prop_table(data = grammar_statistical_small, taxonomy_matrix = taxonomy_matrix)
proportions_phonology_logical <- generate_per_family_prop_table(data = phonology_logical, taxonomy_matrix = taxonomy_matrix)
proportions_phonology_statistical <- generate_per_family_prop_table(data = phonology_statistical, taxonomy_matrix = taxonomy_matrix)
proportions_lexicon_logical <- generate_per_family_prop_table(data = lexicon_logical, taxonomy_matrix = taxonomy_matrix)
proportions_lexicon_statistical <- generate_per_family_prop_table(data = lexicon_statistical, taxonomy_matrix = taxonomy_matrix)

# # proportions of NA per family
# boxplot(apply(proportions_full_logical_large,1,function(x)sum(is.na(x))/ncol(proportions_full_logical_large)))
# title("Proportion of NA per family, large densified logical dataset")
# 
# boxplot(apply(proportions_full_logical_small,1,function(x)sum(is.na(x))/ncol(proportions_full_logical_small)))
# title("Proportion of NA per family, small densified logical dataset")
# 
# boxplot(apply(proportions_full_statistical_large,1,function(x)sum(is.na(x))/ncol(proportions_full_statistical_large)))
# title("Proportion of NA per family, large densified statistical dataset")
# 
# boxplot(apply(proportions_full_statistical_small,1,function(x)sum(is.na(x))/ncol(proportions_full_statistical_small)))
# title("Proportion of NA per family, small densified statistical dataset")
# 
# # proportion of NA per variable
# boxplot(apply(proportions_full_logical_large,2,function(x)sum(is.na(x))/nrow(proportions_full_logical_large)))
# title("Proportion of NA per variable, large densified logical dataset")
# 
# boxplot(apply(proportions_full_logical_small,2,function(x)sum(is.na(x))/nrow(proportions_full_logical_small)))
# title("Proportion of NA per variable, small densified logical dataset")
# 
# boxplot(apply(proportions_full_statistical_large,2,function(x)sum(is.na(x))/nrow(proportions_full_statistical_large)))
# title("Proportion of NA per variable, large densified statistical dataset")
# 
# boxplot(apply(proportions_full_statistical_small,2,function(x)sum(is.na(x))/nrow(proportions_full_statistical_small)))
# title("Proprotion of NA per variable, small statistical dataset")

# perform PCA
proportions_full_logical_large[proportions_full_logical_large=="NaN"]<-NA
proportions_full_logical_small[proportions_full_logical_small=="NaN"]<-NA
proportions_full_statistical_large[proportions_full_statistical_large=="NaN"]<-NA
proportions_full_statistical_small[proportions_full_statistical_small=="NaN"]<-NA
proportions_grammar_logical_large[proportions_grammar_logical_large=="NaN"]<-NA
proportions_grammar_logical_small[proportions_grammar_logical_small=="NaN"]<-NA
proportions_grammar_statistical_large[proportions_grammar_statistical_large=="NaN"]<-NA
proportions_grammar_statistical_small[proportions_grammar_statistical_small=="NaN"]<-NA
proportions_phonology_logical[proportions_phonology_logical=="NaN"]<-NA
proportions_phonology_statistical[proportions_phonology_statistical=="NaN"]<-NA
proportions_lexicon_logical[proportions_lexicon_logical=="NaN"]<-NA
proportions_lexicon_statistical[proportions_lexicon_statistical=="NaN"]<-NA

per_family_props_PCA_full_logical_large <- pca(proportions_full_logical_large, method='rnipals', nPcs=10)
per_family_props_PCA_full_logical_small <- pca(proportions_full_logical_small, method='rnipals', nPcs=10)
per_family_props_PCA_full_statistical_large <- pca(proportions_full_statistical_large, method='rnipals', nPcs=10)
per_family_props_PCA_full_statistical_small <- pca(proportions_full_statistical_small, method='rnipals', nPcs=10)

per_family_props_PCA_grammar_logical_large <- pca(proportions_grammar_logical_large, method='rnipals', nPcs=10)
per_family_props_PCA_grammar_logical_small <- pca(proportions_grammar_logical_small, method='rnipals', nPcs=10)
per_family_props_PCA_grammar_statistical_large <- pca(proportions_grammar_statistical_large, method='rnipals', nPcs=10)
per_family_props_PCA_grammar_statistical_small <- pca(proportions_grammar_statistical_small, method='rnipals', nPcs=10)

per_family_props_PCA_phonology_logical <- pca(proportions_phonology_logical, method='rnipals', nPcs=10)
per_family_props_PCA_phonology_statistical <- pca(proportions_phonology_statistical, method='rnipals', nPcs=10)
per_family_props_PCA_lexicon_logical <- pca(proportions_lexicon_logical, method='rnipals', nPcs=10)
per_family_props_PCA_lexicon_statistical <- pca(proportions_lexicon_statistical, method='rnipals', nPcs=10)

# explained variances
explained_variance(per_family_props_PCA_full_logical_large)[1:9]
explained_variance(per_family_props_PCA_full_logical_small)[1:9]
explained_variance(per_family_props_PCA_full_statistical_large)[1:9]
explained_variance(per_family_props_PCA_full_statistical_small)[1:9]

explained_variance(per_family_props_PCA_grammar_logical_large)[1:9]
explained_variance(per_family_props_PCA_grammar_logical_small)[1:9]
explained_variance(per_family_props_PCA_grammar_statistical_large)[1:9]
explained_variance(per_family_props_PCA_grammar_statistical_small)[1:9]

explained_variance(per_family_props_PCA_phonology_logical)[1:9]
explained_variance(per_family_props_PCA_phonology_statistical)[1:9]
explained_variance(per_family_props_PCA_lexicon_logical)[1:9]
explained_variance(per_family_props_PCA_lexicon_statistical)[1:9]

# most relevant loadings
report_loadings(per_family_props_PCA_full_logical_large, directory = "logicalMM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_full_logical_large, directory = "logicalMM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_full_logical_large, directory = "logicalMM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_full_logical_small, directory = "logicalMM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_full_logical_small, directory = "logicalMM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_full_logical_small, directory = "logicalMM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_full_statistical_large, directory = "statisticalMM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_full_statistical_large, directory = "statisticalMM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_full_statistical_large, directory = "statisticalMM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_full_statistical_small, directory = "statisticalMM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_full_statistical_small, directory = "statisticalMM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_full_statistical_small, directory = "statisticalMM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_grammar_logical_large, directory = "logicalMM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_grammar_logical_large, directory = "logicalMM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_grammar_logical_large, directory = "logicalMM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_grammar_logical_small, directory = "logicalMM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_grammar_logical_small, directory = "logicalMM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_grammar_logical_small, directory = "logicalMM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_grammar_statistical_large, directory = "statisticalMM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_grammar_statistical_large, directory = "statisticalMM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_grammar_statistical_large, directory = "statisticalMM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_grammar_statistical_small, directory = "statisticalMM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_grammar_statistical_small, directory = "statisticalMM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_grammar_statistical_small, directory = "statisticalMM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_phonology_logical, directory = "logicalMM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_phonology_logical, directory = "logicalMM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_phonology_logical, directory = "logicalMM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_phonology_statistical, directory = "statisticalMM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_phonology_statistical, directory = "statisticalMM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_phonology_statistical, directory = "statisticalMM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_lexicon_logical, directory = "logicalMM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_lexicon_logical, directory = "logicalMM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_lexicon_logical, directory = "logicalMM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_lexicon_statistical, directory = "statisticalMM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_lexicon_statistical, directory = "statisticalMM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_lexicon_statistical, directory = "statisticalMM", pc = 3, nr_variables_from_top_and_bottom = 5)

# illustration of rgb cube
steps=.02
mygrid <- as.data.frame(expand.grid(x=seq(0,1,steps), y=seq(0,1,steps), z=seq(0,1,steps)))
m <- mygrid %>%  mutate(rgbcolors=rgb(red=x, green=y, blue=z, alpha=1))
cex = 1
rgb_cube <- cloud(z~x*y, m,
        xlab=list(label='PC1: XXXX', rot=30,  cex=cex),
        ylab=list(label="PC2: XXXX", rot=-39, cex=cex), # note: the label flips the true order because the arrow goes from right to left!
        zlab=list(label="PC3: XXXX", rot=95,  cex=cex),
        pch=15,
        cex=cex+1/3*cex,
        col=m$rgbcolors,
        par.box=list(lty=0))

# map languages to colours
rgb_mapping_full_logical_large <- RGB_mapping(per_family_props_PCA_full_logical_large,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(full_logical_large)) %>% na.omit()
rgb_mapping_full_logical_small <- RGB_mapping(per_family_props_PCA_full_logical_small,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(full_logical_small)) %>% na.omit()
rgb_mapping_full_statistical_large <- RGB_mapping(per_family_props_PCA_full_statistical_large,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(full_statistical_large)) %>% na.omit()
rgb_mapping_full_statistical_small <- RGB_mapping(per_family_props_PCA_full_statistical_small,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(full_statistical_small)) %>% na.omit()

rgb_mapping_grammar_logical_large <- RGB_mapping(per_family_props_PCA_grammar_logical_large,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(grammar_logical_large)) %>% na.omit()
rgb_mapping_grammar_logical_small <- RGB_mapping(per_family_props_PCA_grammar_logical_small,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(grammar_logical_small)) %>% na.omit()
rgb_mapping_grammar_statistical_large <- RGB_mapping(per_family_props_PCA_grammar_statistical_large,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(grammar_statistical_large)) %>% na.omit()
rgb_mapping_grammar_statistical_small <- RGB_mapping(per_family_props_PCA_grammar_statistical_small,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(grammar_statistical_small)) %>% na.omit()

rgb_mapping_phonology_logical <- RGB_mapping(per_family_props_PCA_phonology_logical,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(phonology_logical)) %>% na.omit()
rgb_mapping_phonology_statistical <- RGB_mapping(per_family_props_PCA_phonology_statistical,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(phonology_statistical)) %>% na.omit()
rgb_mapping_lexicon_logical <- RGB_mapping(per_family_props_PCA_lexicon_logical,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(lexicon_logical)) %>% na.omit()
rgb_mapping_lexicon_statistical <- RGB_mapping(per_family_props_PCA_lexicon_statistical,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(lexicon_statistical)) %>% na.omit()

# global PC maps for PCs 1-3

### download Natural Earth data for sf, including islands
if(!exists("world_map_initial")){
  world_map <- ne_download(scale = 110, category = "physical", type = "land", returnclass = "sf")
  minor_islands <- ne_download(scale = 10, category = "physical", type = "minor_islands", returnclass = "sf")
  world_map_initial <- rbind(world_map, minor_islands)
}

rgb_map_full_logical_large <- rgb_map(rgb_mapping_full_logical_large, world_map_initial, main="Logical curation, large densified dataset", plot=F)
rgb_map_full_logical_small <- rgb_map(rgb_mapping_full_logical_small, world_map_initial, main="Logical curation, small densified dataset", plot=F)
rgb_map_full_statistical_large <- rgb_map(rgb_mapping_full_statistical_large, world_map_initial, main="Statistical curation, large densified dataset", plot=F)
rgb_map_full_statistical_small <- rgb_map(rgb_mapping_full_statistical_small, world_map_initial, main="Statistical curation, large densified dataset", plot=F)

rgb_map_grammar_logical_large <- rgb_map(rgb_mapping_grammar_logical_large, world_map_initial, main="Logical curation, large densified dataset", plot=F)
rgb_map_grammar_logical_small <- rgb_map(rgb_mapping_grammar_logical_small, world_map_initial, main="Logical curation, small densified dataset", plot=F)
rgb_map_grammar_statistical_large <- rgb_map(rgb_mapping_grammar_statistical_large, world_map_initial, main="Statistical curation, large densified dataset", plot=F)
rgb_map_grammar_statistical_small <- rgb_map(rgb_mapping_grammar_statistical_small, world_map_initial, main="Statistical curation, large densified dataset", plot=F)

rgb_map_phonology_logical <- rgb_map(rgb_mapping_phonology_logical, world_map_initial, main="Logical curation, large densified dataset", plot=F)
rgb_map_phonology_statistical <- rgb_map(rgb_mapping_phonology_statistical, world_map_initial, main="Statistical curation, large densified dataset", plot=F)
rgb_map_lexicon_logical <- rgb_map(rgb_mapping_lexicon_logical, world_map_initial, main="Logical curation, small densified dataset", plot=F)
rgb_map_lexicon_statistical <- rgb_map(rgb_mapping_lexicon_statistical, world_map_initial, main="Statistical curation, large densified dataset", plot=F)


# arrange the plots side by side, print and save
combined_plot_full <- grid.arrange(rgb_map_full_logical_large, rgb_map_full_statistical_large,
                              rgb_map_full_logical_small, rgb_map_full_statistical_small,
                              ncol = 2)
ggsave("output/family_collapsed_rgb_pca_maps_full.png", plot = combined_plot_full, width = 14, height = 8, dpi = 300)

combined_plot_grammar <- grid.arrange(rgb_map_grammar_logical_large, rgb_map_grammar_statistical_large,
                                   rgb_map_grammar_logical_small, rgb_map_grammar_statistical_small,
                                   ncol = 2)
ggsave("output/family_collapsed_rgb_pca_maps_grammar.png", plot = combined_plot_grammar, width = 14, height = 8, dpi = 300)

combined_plot_phonology <- grid.arrange(rgb_map_phonology_logical, rgb_map_phonology_statistical,
                                   ncol = 2)
ggsave("output/family_collapsed_rgb_pca_maps_phonology.png", plot = combined_plot_phonology, width = 14, height = 4, dpi = 300)

combined_plot_lexicon <- grid.arrange(rgb_map_lexicon_logical, rgb_map_lexicon_statistical,
                                   ncol = 2)
ggsave("output/family_collapsed_rgb_pca_maps_lexicon.png", plot = combined_plot_lexicon, width = 14, height = 4, dpi = 300)

##################
## Grid points and entropies at grid points

# coordinate_scaling 20 yields 715 grid points if buffer_distance is 150'000, changing coordinate_scaling gives more or less grid points
coordinate_scaling <- 20  # 715 grid points if coordinate_scaling=20 and buffer_distance=150000
buffer_distance <- 150000  # defining buffer distance for map for some Pacific islands to be kept
lon <- 155

# trim taxonomy to appropriate languages
full_logical_large_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(full_logical_large)) %>% filter(!is.na(lon))
full_logical_small_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(full_logical_small)) %>% filter(!is.na(lon))
full_statistical_large_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(full_statistical_large)) %>% filter(!is.na(lon))
full_statistical_small_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(full_statistical_small)) %>% filter(!is.na(lon))
grammar_logical_large_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(grammar_logical_large)) %>% filter(!is.na(lon))
grammar_logical_small_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(grammar_logical_small)) %>% filter(!is.na(lon))
grammar_statistical_large_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(grammar_statistical_large)) %>% filter(!is.na(lon))
grammar_statistical_small_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(grammar_statistical_small)) %>% filter(!is.na(lon))
phonology_logical_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(phonology_logical)) %>% filter(!is.na(lon))
phonology_statistical_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(phonology_statistical)) %>% filter(!is.na(lon))
lexicon_logical_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(lexicon_logical)) %>% filter(!is.na(lon))
lexicon_statistical_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(lexicon_statistical)) %>% filter(!is.na(lon))

# trim data to appropriate languages
full_logical_large_with_lg_coords <- full_logical_large[which(rownames(full_logical_large)%in%full_logical_large_taxonomy$id),]
full_logical_small_with_lg_coords <- full_logical_small[which(rownames(full_logical_small)%in%full_logical_small_taxonomy$id),]
full_statistical_large_with_lg_coords <- full_statistical_large[which(rownames(full_statistical_large)%in%full_statistical_large_taxonomy$id),]
full_statistical_small_with_lg_coords <- full_statistical_small[which(rownames(full_statistical_small)%in%full_statistical_small_taxonomy$id),]
grammar_logical_large_with_lg_coords <- grammar_logical_large[which(rownames(grammar_logical_large)%in%grammar_logical_large_taxonomy$id),]
grammar_logical_small_with_lg_coords <- grammar_logical_small[which(rownames(grammar_logical_small)%in%grammar_logical_small_taxonomy$id),]
grammar_statistical_large_with_lg_coords <- grammar_statistical_large[which(rownames(grammar_statistical_large)%in%grammar_statistical_large_taxonomy$id),]
grammar_statistical_small_with_lg_coords <- grammar_statistical_small[which(rownames(grammar_statistical_small)%in%grammar_statistical_small_taxonomy$id),]
phonology_logical_with_lg_coords <- phonology_logical[which(rownames(phonology_logical)%in%phonology_logical_taxonomy$id),]
phonology_statistical_with_lg_coords <- phonology_statistical[which(rownames(phonology_statistical)%in%phonology_statistical_taxonomy$id),]
lexicon_logical_with_lg_coords <- lexicon_logical[which(rownames(lexicon_logical)%in%lexicon_logical_taxonomy$id),]
lexicon_statistical_with_lg_coords <- lexicon_statistical[which(rownames(lexicon_statistical)%in%lexicon_statistical_taxonomy$id),]

full_logical_large_gps <- gridpointwise_entropies(taxonomy_matrix = full_logical_large_taxonomy,
                                                   data_trimmed_to_lgs_with_coords = full_logical_large_with_lg_coords, 
                                                   world_map = world_map_including_islands, 
                                                   lon = 155, 
                                                   coordinate_scaling = coordinate_scaling, 
                                                   buffer_distance = buffer_distance, 
                                                   verbose = T, 
                                                   directory = "logicalMM/full",
                                                   data_type = "densified_large",
                                                  robust = TRUE, robust_errors = "with_robust_errors")

full_logical_small_gps <- gridpointwise_entropies(taxonomy_matrix = full_logical_small_taxonomy,
                                                 data_trimmed_to_lgs_with_coords = full_logical_small_with_lg_coords, 
                                                 world_map = world_map_including_islands, 
                                                 lon = 155, 
                                                 coordinate_scaling = coordinate_scaling, 
                                                 buffer_distance = buffer_distance, 
                                                 verbose = T, 
                                                 directory = "logicalMM/full",
                                                 data_type = "densified_small",
                                                 robust = TRUE, robust_errors = "with_robust_errors")

full_statistical_large_gps <- gridpointwise_entropies(taxonomy_matrix = full_statistical_large_taxonomy,
                                                data_trimmed_to_lgs_with_coords = full_statistical_large_with_lg_coords, 
                                                world_map = world_map_including_islands, 
                                                lon = 155, 
                                                coordinate_scaling = coordinate_scaling, 
                                                buffer_distance = buffer_distance, 
                                                verbose = T, 
                                                directory = "statisticalMM/full",
                                                data_type = "densified_large",
                                                robust = TRUE, robust_errors = "with_robust_errors")

full_statistical_small_gps <-gridpointwise_entropies(taxonomy_matrix = full_statistical_small_taxonomy,
                                                    data_trimmed_to_lgs_with_coords = full_statistical_small_with_lg_coords, 
                                                    world_map = world_map_including_islands, 
                                                    lon = 155, 
                                                    coordinate_scaling = coordinate_scaling, 
                                                    buffer_distance = buffer_distance, 
                                                    verbose = T, 
                                                    directory = "statisticalMM/full",
                                                    data_type = "densified_small",
                                                    robust = TRUE, robust_errors = "with_robust_errors")

grammar_logical_large_gps <- gridpointwise_entropies(taxonomy_matrix = grammar_logical_large_taxonomy,
                                                  data_trimmed_to_lgs_with_coords = grammar_logical_large_with_lg_coords, 
                                                  world_map = world_map_including_islands, 
                                                  lon = 155, 
                                                  coordinate_scaling = coordinate_scaling, 
                                                  buffer_distance = buffer_distance, 
                                                  verbose = T, 
                                                  directory = "logicalMM/grammar",
                                                  data_type = "densified_large",
                                                  robust = TRUE, robust_errors = "with_robust_errors")

grammar_logical_small_gps <- gridpointwise_entropies(taxonomy_matrix = grammar_logical_small_taxonomy,
                                                  data_trimmed_to_lgs_with_coords = grammar_logical_small_with_lg_coords, 
                                                  world_map = world_map_including_islands, 
                                                  lon = 155, 
                                                  coordinate_scaling = coordinate_scaling, 
                                                  buffer_distance = buffer_distance, 
                                                  verbose = T, 
                                                  directory = "logicalMM/grammar",
                                                  data_type = "densified_small",
                                                  robust = TRUE, robust_errors = "with_robust_errors")

grammar_statistical_large_gps <- gridpointwise_entropies(taxonomy_matrix = grammar_statistical_large_taxonomy,
                                                      data_trimmed_to_lgs_with_coords = grammar_statistical_large_with_lg_coords, 
                                                      world_map = world_map_including_islands, 
                                                      lon = 155, 
                                                      coordinate_scaling = coordinate_scaling, 
                                                      buffer_distance = buffer_distance, 
                                                      verbose = T, 
                                                      directory = "statisticalMM/grammar",
                                                      data_type = "densified_large",
                                                      robust = TRUE, robust_errors = "with_robust_errors")

grammar_statistical_small_gps <-gridpointwise_entropies(taxonomy_matrix = grammar_statistical_small_taxonomy,
                                                     data_trimmed_to_lgs_with_coords = grammar_statistical_small_with_lg_coords, 
                                                     world_map = world_map_including_islands, 
                                                     lon = 155, 
                                                     coordinate_scaling = coordinate_scaling, 
                                                     buffer_distance = buffer_distance, 
                                                     verbose = T, 
                                                     directory = "statisticalMM/grammar",
                                                     data_type = "densified_small",
                                                     robust = TRUE, robust_errors = "with_robust_errors")


phonology_logical_gps <- gridpointwise_entropies(taxonomy_matrix = phonology_logical_taxonomy,
                                                  data_trimmed_to_lgs_with_coords = phonology_logical_with_lg_coords, 
                                                  world_map = world_map_including_islands, 
                                                  lon = 155, 
                                                  coordinate_scaling = coordinate_scaling, 
                                                  buffer_distance = buffer_distance, 
                                                  verbose = T, 
                                                  directory = "logicalMM/phonology",
                                                  data_type = "densified",
                                                  robust = TRUE, robust_errors = "with_robust_errors")

phonology_statistical_gps <- gridpointwise_entropies(taxonomy_matrix = phonology_statistical_taxonomy,
                                                      data_trimmed_to_lgs_with_coords = phonology_statistical_with_lg_coords, 
                                                      world_map = world_map_including_islands, 
                                                      lon = 155, 
                                                      coordinate_scaling = coordinate_scaling, 
                                                      buffer_distance = buffer_distance, 
                                                      verbose = T, 
                                                      directory = "statisticalMM/phonology",
                                                      data_type = "densified",
                                                      robust = TRUE, robust_errors = "with_robust_errors")

lexicon_logical_gps <- gridpointwise_entropies(taxonomy_matrix = lexicon_logical_taxonomy,
                                                 data_trimmed_to_lgs_with_coords = lexicon_logical_with_lg_coords, 
                                                 world_map = world_map_including_islands, 
                                                 lon = 155, 
                                                 coordinate_scaling = coordinate_scaling, 
                                                 buffer_distance = buffer_distance, 
                                                 verbose = T, 
                                                 directory = "logicalMM/lexicon",
                                                 data_type = "densified",
                                                 robust = TRUE, robust_errors = "with_robust_errors")

lexicon_statistical_gps <- gridpointwise_entropies(taxonomy_matrix = lexicon_statistical_taxonomy,
                                                     data_trimmed_to_lgs_with_coords = lexicon_statistical_with_lg_coords, 
                                                     world_map = world_map_including_islands, 
                                                     lon = 155, 
                                                     coordinate_scaling = coordinate_scaling, 
                                                     buffer_distance = buffer_distance, 
                                                     verbose = T, 
                                                     directory = "statisticalMM/lexicon",
                                                     data_type = "densified",
                                                     robust = TRUE, robust_errors = "with_robust_errors")

# compare curations (here: logical vs. statistical): we expect there to be systematic differences
# we do not plot lg/fam comparisons here, because we know nothing changes
grid_point_comparison_horizontal(baseline_gp_frame = full_logical_large_gps, 
                                 name_baseline = "large densified (logical curation)",
                                 comparison_gp_frame = full_statistical_large_gps, 
                                 name_comparison = "large densified (statistical curation)",
                                 lg_and_fam_comparison = TRUE, fig_width = 15, fig_height = 10,
                                 name_for_plot = "_full_large_horizontal_comparison",
                                 world_map = world_map_including_islands, lon = 155)

grid_point_comparison_horizontal(baseline_gp_frame = full_logical_small_gps, 
                                 name_baseline = "small densified (logical curation)",
                                 comparison_gp_frame = full_statistical_small_gps, 
                                 name_comparison = "small densified (statistical curation)",
                                 lg_and_fam_comparison = TRUE, fig_width = 15, fig_height = 10,
                                 name_for_plot = "_full_small_horizontal_comparison",
                                 world_map = world_map_including_islands, lon = 155)

grid_point_comparison_horizontal(baseline_gp_frame = grammar_logical_large_gps, 
                                 name_baseline = "large densified (logical curation)",
                                 comparison_gp_frame = grammar_statistical_large_gps, 
                                 name_comparison = "large densified (statistical curation)",
                                 lg_and_fam_comparison = TRUE, fig_width = 15, fig_height = 10,
                                 name_for_plot = "_grammar_large_horizontal_comparison",
                                 world_map = world_map_including_islands, lon = 155)

grid_point_comparison_horizontal(baseline_gp_frame = grammar_logical_small_gps, 
                                 name_baseline = "small densified (logical curation)",
                                 comparison_gp_frame = grammar_statistical_small_gps, 
                                 name_comparison = "small densified (statistical curation)",
                                 lg_and_fam_comparison = TRUE, fig_width = 15, fig_height = 10,
                                 name_for_plot = "_grammar_small_horizontal_comparison",
                                 world_map = world_map_including_islands, lon = 155)

grid_point_comparison_horizontal(baseline_gp_frame = phonology_logical_gps, 
                                 name_baseline = "densified (logical curation)",
                                 comparison_gp_frame = phonology_statistical_gps, 
                                 name_comparison = "densified (statistical curation)",
                                 lg_and_fam_comparison = TRUE, fig_width = 15, fig_height = 10,
                                 name_for_plot = "_phonology_horizontal_comparison",
                                 world_map = world_map_including_islands, lon = 155)

grid_point_comparison_horizontal(baseline_gp_frame = lexicon_logical_gps, 
                                 name_baseline = "densified (logical curation)",
                                 comparison_gp_frame = lexicon_statistical_gps, 
                                 name_comparison = "densified (statistical curation)",
                                 lg_and_fam_comparison = TRUE, fig_width = 15, fig_height = 10,
                                 name_for_plot = "_lexicon_horizontal_comparison",
                                 world_map = world_map_including_islands, lon = 155)



# compare large densified to small densified: we expect there to be no systematic differences
# full
grid_point_comparison_vertical(baseline_gp_frame = full_logical_large_gps, 
                               name_baseline = "large densified dataset (full, logical curation)", 
                               comparison_gp_frame = full_logical_small_gps,
                               name_comparison = "small densified dataset (full, logical curation)", 
                               world_map = world_map_including_islands, lon = 155,
                               directory = "logicalMM/full/")
grid_point_comparison_vertical(baseline_gp_frame = full_statistical_large_gps, 
                               name_baseline = "large densified dataset (full, statistical curation)", 
                               comparison_gp_frame = full_statistical_small_gps,
                               name_comparison = "small densified dataset (full, statistical curation)", 
                               world_map = world_map_including_islands, lon = 155,
                               directory = "statisticalMM/full/")

# grammar
grid_point_comparison_vertical(baseline_gp_frame = grammar_logical_large_gps, 
                               name_baseline = "large densified dataset (grammar, logical curation)", 
                               comparison_gp_frame = grammar_logical_small_gps,
                               name_comparison = "small densified dataset (grammar, logical curation)", 
                               world_map = world_map_including_islands, lon = 155,
                               directory = "logicalMM/grammar/")
grid_point_comparison_vertical(baseline_gp_frame = grammar_statistical_large_gps, 
                               name_baseline = "large densified dataset (grammar, statistical curation)", 
                               comparison_gp_frame = grammar_statistical_small_gps,
                               name_comparison = "small densified dataset (grammar, statistical curation)", 
                               world_map = world_map_including_islands, lon = 155,
                               directory = "statisticalMM/grammar/")

