rm(list=ls())
library(densify)
library(tidyverse)
library(cluster)
library(pcaMethods)
library(geodist)
library(sp)
library(ape)
library(reshape)
library(ggcorrplot)
library(gridExtra)
library(rnaturalearth)
library(lattice)

source("../functions.R")

# generate taxonomy
taxonomy_matrix <- as_flat_taxonomy_matrix(glottolog_languoids)
# link coordinates to lgs
macroareas <- read.csv("../input/glottolog_v.4.8/languages_and_dialects_geo.csv")
names(macroareas)[1] <- "id"
taxonomy_matrix_for_plot <- merge(taxonomy_matrix,select(macroareas, c("id","latitude","longitude","macroarea")),by="id") %>% filter(!is.na(latitude))
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

# explained variances and R2 among all PCs
explained_variance(per_family_props_PCA_full_logical_large)[,1:10]
explained_variance(per_family_props_PCA_full_logical_small)[,1:10]
explained_variance(per_family_props_PCA_full_statistical_large)[,1:10]
explained_variance(per_family_props_PCA_full_statistical_small)[,1:10]

explained_variance(per_family_props_PCA_grammar_logical_large)[,1:10]
explained_variance(per_family_props_PCA_grammar_logical_small)[,1:10]
explained_variance(per_family_props_PCA_grammar_statistical_large)[,1:10]
explained_variance(per_family_props_PCA_grammar_statistical_small)[,1:10]

explained_variance(per_family_props_PCA_phonology_logical)[,1:10]
explained_variance(per_family_props_PCA_phonology_statistical)[,1:10]
explained_variance(per_family_props_PCA_lexicon_logical)[,1:10]
explained_variance(per_family_props_PCA_lexicon_statistical)[,1:10]

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

trellis.device(png, filename = "../figures-for-paper/main/rgb_cube.png", type="cairo", units="in", width=5, height=5,  pointsize=10, res=300)
print(rgb_cube)
dev.off()

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

rgb_map_full_logical_large <- rgb_map(rgb_mapping_full_logical_large, world_map_initial, main="Mega-Merge, all variables, logical curation, large densified dataset", plot=F)
rgb_map_full_logical_small <- rgb_map(rgb_mapping_full_logical_small, world_map_initial, main="Mega-Merge, all variables, logical curation, small densified dataset", plot=F)
rgb_map_full_statistical_large <- rgb_map(rgb_mapping_full_statistical_large, world_map_initial, main="Mega-Merge, all variables, statistical curation, large densified dataset", plot=F)
rgb_map_full_statistical_small <- rgb_map(rgb_mapping_full_statistical_small, world_map_initial, main="Mega-Merge, all variables, statistical curation, small densified dataset", plot=F)

rgb_map_grammar_logical_large <- rgb_map(rgb_mapping_grammar_logical_large, world_map_initial, main="Mega-Merge, grammar only, logical curation, large densified dataset", plot=F)
rgb_map_grammar_logical_small <- rgb_map(rgb_mapping_grammar_logical_small, world_map_initial, main="Mega-Merge, grammar only, logical curation, small densified dataset", plot=F)
rgb_map_grammar_statistical_large <- rgb_map(rgb_mapping_grammar_statistical_large, world_map_initial, main="Mega-Merge, grammar only, statistical curation, large densified dataset", plot=F)
rgb_map_grammar_statistical_small <- rgb_map(rgb_mapping_grammar_statistical_small, world_map_initial, main="Mega-Merge, grammar only, statistical curation, small densified dataset", plot=F)

rgb_map_phonology_logical <- rgb_map(rgb_mapping_phonology_logical, world_map_initial, main="Mega-Merge, phonology only, logical curation, densified dataset", plot=F)
rgb_map_phonology_statistical <- rgb_map(rgb_mapping_phonology_statistical, world_map_initial, main="Mega-Merge, phonology only, statistical curation, densified dataset", plot=F)
rgb_map_lexicon_logical <- rgb_map(rgb_mapping_lexicon_logical, world_map_initial, main="Mega-Merge, lexicon only, logical curation, densified dataset", plot=F)
rgb_map_lexicon_statistical <- rgb_map(rgb_mapping_lexicon_statistical, world_map_initial, main="Mega-Merge, lexicon only, statistical curation, densified dataset", plot=F)

#### family relationships to one another in PC-space, by macroarea
f2d_full_logical_large <- dim_mapping(pca.object = per_family_props_PCA_full_logical_large, taxonomy_matrix = taxonomy_matrix_for_plot)
f2d_full_logical_small <- dim_mapping(pca.object = per_family_props_PCA_full_logical_small, taxonomy_matrix = taxonomy_matrix_for_plot)
f2d_full_statistical_large <- dim_mapping(pca.object = per_family_props_PCA_full_statistical_large, taxonomy_matrix = taxonomy_matrix_for_plot)
f2d_full_statistical_small <- dim_mapping(pca.object = per_family_props_PCA_full_statistical_small, taxonomy_matrix = taxonomy_matrix_for_plot)

f2d_grammar_logical_large <- dim_mapping(pca.object = per_family_props_PCA_grammar_logical_large, taxonomy_matrix = taxonomy_matrix_for_plot)
f2d_grammar_logical_small <- dim_mapping(pca.object = per_family_props_PCA_grammar_logical_small, taxonomy_matrix = taxonomy_matrix_for_plot)
f2d_grammar_statistical_large <- dim_mapping(pca.object = per_family_props_PCA_grammar_statistical_large, taxonomy_matrix = taxonomy_matrix_for_plot)
f2d_grammar_statistical_small <- dim_mapping(pca.object = per_family_props_PCA_grammar_statistical_small, taxonomy_matrix = taxonomy_matrix_for_plot)

f2d_phonology_logical <- dim_mapping(pca.object = per_family_props_PCA_phonology_logical, taxonomy_matrix = taxonomy_matrix_for_plot)
f2d_phonology_statistical <- dim_mapping(pca.object = per_family_props_PCA_phonology_statistical, taxonomy_matrix = taxonomy_matrix_for_plot)
f2d_lexicon_logical <- dim_mapping(pca.object = per_family_props_PCA_lexicon_logical, taxonomy_matrix = taxonomy_matrix_for_plot)
f2d_lexicon_statistical <- dim_mapping(pca.object = per_family_props_PCA_lexicon_statistical, taxonomy_matrix = taxonomy_matrix_for_plot)

# arrange the plots side by side, print and save
logical_plots <- grid.arrange(
  arrangeGrob(rgb_map_full_logical_large, f2d_full_logical_large, rgb_map_full_logical_small, f2d_full_logical_small, nrow = 2, ncol = 2),
  arrangeGrob(rgb_map_grammar_logical_large, f2d_grammar_logical_large, rgb_map_grammar_logical_small, f2d_grammar_logical_small, nrow = 2, ncol = 2),
  arrangeGrob(rgb_map_phonology_logical, f2d_phonology_logical, rgb_map_lexicon_logical, f2d_lexicon_logical, nrow = 2, ncol = 2),
  ncol = 3
)
ggsave("../figures-for-paper/main/logical_family_collapsed_rgb_pca_maps_mm.png", plot = logical_plots, width = 14*3, height = 8, dpi = 300)

all_plots <- grid.arrange(
  arrangeGrob(rgb_map_full_logical_large, f2d_full_logical_large, rgb_map_full_statistical_large, f2d_full_statistical_large, rgb_map_full_logical_small, f2d_full_logical_small, rgb_map_full_statistical_small, f2d_full_statistical_small, nrow = 2, ncol = 4),
  arrangeGrob(rgb_map_grammar_logical_large, f2d_grammar_logical_large, rgb_map_grammar_statistical_large, f2d_grammar_statistical_large, rgb_map_grammar_logical_small, f2d_grammar_logical_small, rgb_map_grammar_statistical_small, f2d_grammar_statistical_small, nrow = 2, ncol = 4),
  arrangeGrob(rgb_map_phonology_logical, f2d_phonology_logical, rgb_map_phonology_statistical, f2d_phonology_statistical, rgb_map_lexicon_logical, f2d_lexicon_logical, rgb_map_lexicon_statistical, f2d_lexicon_statistical, nrow = 2, ncol = 4),
  ncol = 1,
  nrow = 3
)
ggsave("../figures-for-paper/supplementary/all_family_collapsed_rgb_pca_maps_mm.png", plot = all_plots, width = 14*2, height = 8*3, dpi = 300)


##################
## Grid points and entropies at grid points

# coordinate_scaling 20 yields ~750 grid points if buffer_distance is 150'000, changing coordinate_scaling gives more or less grid points
coordinate_scaling <- 20  # ~750 grid points if coordinate_scaling=20 and buffer_distance=150000
buffer_distance <- 150000  # defining buffer distance for map for some Pacific islands to be kept

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
                                                   world_map = world_map_initial, 
                                                   coordinate_scaling = coordinate_scaling, 
                                                   buffer_distance = buffer_distance, 
                                                   verbose = T, 
                                                   main_fig = F,
                                                   data_type = "Mega-Merge (all variables, logical, large densified)")

full_statistical_large_gps <- gridpointwise_entropies(taxonomy_matrix = full_statistical_large_taxonomy,
                                                data_trimmed_to_lgs_with_coords = full_statistical_large_with_lg_coords, 
                                                world_map = world_map_initial, 
                                                coordinate_scaling = coordinate_scaling, 
                                                buffer_distance = buffer_distance, 
                                                verbose = T, 
                                                main_fig = T,
                                                data_type = "Mega-Merge (all variables, statistical, large densified)")

grammar_logical_large_gps <- gridpointwise_entropies(taxonomy_matrix = grammar_logical_large_taxonomy,
                                                  data_trimmed_to_lgs_with_coords = grammar_logical_large_with_lg_coords, 
                                                  world_map = world_map_initial, 
                                                  coordinate_scaling = coordinate_scaling, 
                                                  buffer_distance = buffer_distance, 
                                                  verbose = T, 
                                                  main_fig = F,
                                                  data_type = "Mega-Merge (grammar only, logical, large densified)")

grammar_statistical_large_gps <- gridpointwise_entropies(taxonomy_matrix = grammar_statistical_large_taxonomy,
                                                      data_trimmed_to_lgs_with_coords = grammar_statistical_large_with_lg_coords, 
                                                      world_map = world_map_initial, 
                                                      coordinate_scaling = coordinate_scaling, 
                                                      buffer_distance = buffer_distance, 
                                                      verbose = T, 
                                                      main_fig = F,
                                                      data_type = "Mega-Merge (grammar only, statistical, large densified)")

phonology_logical_gps <- gridpointwise_entropies(taxonomy_matrix = phonology_logical_taxonomy,
                                                  data_trimmed_to_lgs_with_coords = phonology_logical_with_lg_coords, 
                                                  world_map = world_map_initial, 
                                                  coordinate_scaling = coordinate_scaling, 
                                                  buffer_distance = buffer_distance, 
                                                  verbose = T, 
                                                 main_fig = F,
                                                 data_type = "Mega-Merge (phonology only, logical, densified)")

phonology_statistical_gps <- gridpointwise_entropies(taxonomy_matrix = phonology_statistical_taxonomy,
                                                      data_trimmed_to_lgs_with_coords = phonology_statistical_with_lg_coords, 
                                                      world_map = world_map_initial, 
                                                      coordinate_scaling = coordinate_scaling, 
                                                      buffer_distance = buffer_distance, 
                                                      verbose = T, 
                                                     main_fig = T,
                                                     data_type = "Mega-Merge (phonology only, statistical, densified)")

lexicon_logical_gps <- gridpointwise_entropies(taxonomy_matrix = lexicon_logical_taxonomy,
                                                 data_trimmed_to_lgs_with_coords = lexicon_logical_with_lg_coords, 
                                                 world_map = world_map_initial, 
                                                 coordinate_scaling = coordinate_scaling, 
                                                 buffer_distance = buffer_distance, 
                                                 verbose = T, 
                                               main_fig = F,
                                               data_type = "Mega-Merge (lexicon only, logical, densified)")

lexicon_statistical_gps <- gridpointwise_entropies(taxonomy_matrix = lexicon_statistical_taxonomy,
                                                     data_trimmed_to_lgs_with_coords = lexicon_statistical_with_lg_coords, 
                                                     world_map = world_map_initial, 
                                                     coordinate_scaling = coordinate_scaling, 
                                                     buffer_distance = buffer_distance, 
                                                     verbose = T, 
                                                   main_fig = F,
                                                   data_type = "Mega-Merge (lexicon only, statistical, densified)")

# compare curations (here: logical vs. statistical): we expect there to be systematic differences
# we do not plot lg/fam comparisons here, because we know nothing changes
grid_point_comparison_horizontal(baseline_gp_frame = full_logical_large_gps, 
                                 name_baseline = "large densified (logical curation)",
                                 comparison_gp_frame = full_statistical_large_gps, 
                                 name_comparison = "large densified (statistical curation)",
                                 lg_and_fam_comparison = TRUE, fig_width = 17, fig_height = 10,
                                 name_for_plot = "horizontal_comparison_megamerge_full_large",
                                 world_map = world_map_initial)

grid_point_comparison_horizontal(baseline_gp_frame = phonology_logical_gps, 
                                 name_baseline = "phonology only (logical curation)",
                                 comparison_gp_frame = phonology_statistical_gps, 
                                 name_comparison = "phonology only (statistical curation)",
                                 lg_and_fam_comparison = TRUE, fig_width = 17, fig_height = 10,
                                 name_for_plot = "horizontal_comparison_megamerge_phonology",
                                 world_map = world_map_initial)


# compare large densified to small densified: we expect there to be no systematic differences
# # full
# grid_point_comparison_vertical(baseline_gp_frame = full_statistical_large_gps, 
#                                name_baseline = "large densified dataset (full, statistical curation)", 
#                                comparison_gp_frame = full_statistical_small_gps,
#                                name_comparison = "small densified dataset (full, statistical curation)", 
#                                world_map = world_map_initial)
# 
# # grammar
# grid_point_comparison_vertical(baseline_gp_frame = grammar_statistical_large_gps, 
#                                name_baseline = "large densified dataset (grammar, statistical curation)", 
#                                comparison_gp_frame = grammar_statistical_small_gps,
#                                name_comparison = "small densified dataset (grammar, statistical curation)", 
#                                world_map = world_map_initial)

#### Gower distances, PCoA for small densified datasets only
# compute gower distances, using daisy()
full_statistical_small_gower <- daisy(full_statistical_small_with_lg_coords, metric = "gower") %>% as.matrix()

library(vegan)
set.seed(78)
full_statistical_small_metaMDS <- metaMDS(full_statistical_small_gower, 
                                          k = 3, 
                                          try = 20, trymax = 500,
                                          engine = "monoMDS")


##### genetics vs typology in Eurasia
gelato <- read.csv("../gelato/GELATO_population_glottocode_mapping.csv")

gelato_pops <- filter(gelato,!is.na(MM_statistical_large_proxy)) %>% select(c(PopName,Location,country,lat,lon,glottocodeBase,glottologFamily,comment_MM_proxies,MM_statistical_large_proxy,proxy_mm_statistical_large_density))
gelato_pops <- left_join(gelato_pops,select(full_statistical_large_taxonomy,"id","macroarea"),by=c("MM_statistical_large_proxy" = "id"))
gelato_pops <- gelato_pops %>% filter(!is.na(lat)) %>% filter(macroarea == "Eurasia") 

lgs_in_gelato <- full_statistical_large[which(rownames(full_statistical_large)%in%gelato_pops$MM_statistical_large_proxy),] 

fst_data <- read.delim("../gelato/fst_GelatoHO_mergedSetFeb2022_.txt", header = T)
fst_data <- rbind(c("Abazin","Abazin",NA),fst_data,c("Zoro","Zoro",NA)) # to enable pivot_wider to be 558x558
fst_matrix <- pivot_wider(fst_data, names_from = Pop1, values_from = FST) %>% select(-Pop2) %>% as.matrix()
fst_matrix[fst_matrix < 0] <- 0 # set all negative FST-values to zero (these are artifacts)
fst_matrix<-apply(fst_matrix,1,as.numeric)
# fst_matrix[lower.tri(fst_matrix)] <- t(fst_matrix)[lower.tri(fst_matrix)]
colnames(fst_matrix)<-unique(fst_data$Pop2)
rownames(fst_matrix)<-unique(fst_data$Pop2)
fst_matrix <- fst_matrix[which(rownames(fst_matrix)%in%gelato_pops$PopName),which(colnames(fst_matrix)%in%gelato_pops$PopName)]

# coordinate_scaling 7 yields 30 grid points in Eurasia if buffer_distance is 300'000 
coordinate_scaling <- 7  # ~750 grid points if coordinate_scaling=20 and buffer_distance=150000
buffer_distance <- 300000  # defining buffer distance for map for some Pacific islands to be kept

mm_eurasia_all <- genetics_vs_typology(gelato_pops = gelato_pops,
                                       lgs_in_gelato = lgs_in_gelato,
                                       fst_matrix = fst_matrix,
                                       seed=7,
                                       world_map = world_map_initial, 
                                       coordinate_scaling = coordinate_scaling, 
                                       buffer_distance = buffer_distance,
                                       verbose=T, 
                                       data_type = "mm_full_eurasia")

# repeat procedure without Onge, because they are outliers
gelato_pops_no_onge <- gelato_pops %>% filter(PopName != "Onge")
lgs_in_gelato_no_onge <- lgs_in_gelato[which(rownames(lgs_in_gelato) != "onge1236"),]
fst_matrix_no_onge <- fst_matrix[which(rownames(fst_matrix)!="Onge"),which(colnames(fst_matrix)!="Onge")]

mm_eurasia_all_no_onge <- genetics_vs_typology(gelato_pops = gelato_pops_no_onge,
                                               lgs_in_gelato = lgs_in_gelato_no_onge,
                                               fst_matrix = fst_matrix_no_onge,
                                               seed=7,
                                               world_map = world_map_initial, 
                                               coordinate_scaling = coordinate_scaling, 
                                               buffer_distance = buffer_distance,
                                               verbose=T, 
                                               data_type = "mm_full_eurasia_no_onge")
