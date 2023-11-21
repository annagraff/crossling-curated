rm(list=ls())
library(densify)
library(tidyverse)
library(cluster)
library(pcaMethods)

source("../functions.R")

# generate taxonomy
taxonomy_matrix <- as_flat_taxonomy_matrix(glottolog_languoids)
# link coordinates to lgs
macroareas <- read.csv("../input/glottolog_v.4.8/languages_and_dialects_geo.csv")
names(macroareas)[1] <- "id"
taxonomy_matrix_for_plot <- merge(taxonomy_matrix,select(macroareas, c("id","latitude","longitude")),by="id") %>% filter(!is.na(latitude))
names(taxonomy_matrix_for_plot)[29:30] <- c("lat","lon") # change names to lat and lon for downstream functions

# read in matrices, convert all NA and ? to NA, factorize all variables

# read in original grambank data
original <- read.csv("../input/grambank_v.1.0.3.csv", row.names = "Glottocode")
original <- na_convert(original)
original <- factorise(original)

logical_full <- read.csv("output/logicalGBM/logicalGBM.csv",row.names = "glottocode") %>% select(-X)
logical_full <- na_convert(logical_full)
logical_full <- factorise(logical_full)
logical_full <- logical_full[which(rownames(logical_full)%in%taxonomy_matrix_for_plot$id),]

logical_densified <- read.csv("output/logicalGBM/logicalGBM_pruned.csv", row.names = "X")
logical_densified <- na_convert(logical_densified)
logical_densified <- factorise(logical_densified)

statistical_full <- read.csv("output/statisticalGBM/statisticalGBM.csv",row.names = "glottocode") %>% select(-X)
statistical_full <- na_convert(statistical_full)
statistical_full <- factorise(statistical_full)
statistical_full <- statistical_full[which(rownames(statistical_full)%in%taxonomy_matrix_for_plot$id),]

statistical_densified <- read.csv("output/statisticalGBM/statisticalGBM_pruned.csv", row.names = "X")
statistical_densified <- na_convert(statistical_densified)
statistical_densified <- factorise(statistical_densified)

#### Gower distances
# compute gower distances, using daisy()
logical_full_gower <- daisy(logical_full, metric = "gower") # %>% as.matrix() %>% na.omit() %>% as.dist()
logical_densified_gower <- daisy(logical_densified, metric = "gower")
statistical_full_gower <- daisy(statistical_full, metric = "gower")
statistical_densified_gower <- daisy(statistical_densified, metric = "gower")

##################
## PCA, aggregated to family level
proportions_original <- generate_per_family_prop_table(data = original, taxonomy_matrix = taxonomy_matrix)
proportions_logical_full <- generate_per_family_prop_table(data = logical_full, taxonomy_matrix = taxonomy_matrix)
proportions_logical_densified <- generate_per_family_prop_table(data = logical_densified, taxonomy_matrix = taxonomy_matrix)
proportions_statistical_full <- generate_per_family_prop_table(data = statistical_full, taxonomy_matrix = taxonomy_matrix)
proportions_statistical_densified <- generate_per_family_prop_table(data = statistical_densified, taxonomy_matrix = taxonomy_matrix)

# proportions of NA per family
boxplot(apply(proportions_logical_full,1,function(x)sum(is.na(x))/nrow(proportions_logical_full)))
title("Proportion of NA per family, full logical dataset")

boxplot(apply(proportions_logical_densified,1,function(x)sum(is.na(x))/nrow(proportions_logical_densified)))
title("Proportion of NA per family, densified logical dataset")

boxplot(apply(proportions_statistical_full,1,function(x)sum(is.na(x))/nrow(proportions_statistical_full)))
title("Proportion of NA per family, full statistical dataset")

boxplot(apply(proportions_statistical_densified,1,function(x)sum(is.na(x))/nrow(proportions_statistical_densified)))
title("Proportion of NA per family, densified statistical dataset")

# proportion of NA per variable
boxplot(apply(proportions_logical_full,2,function(x)sum(is.na(x))/nrow(proportions_logical_full)))
title("Proportion of NA per variable, full logical dataset")

boxplot(apply(proportions_logical_densified,2,function(x)sum(is.na(x))/nrow(proportions_logical_densified)))
title("Proportion of NA per variable, densified logical dataset")

boxplot(apply(proportions_statistical_full,2,function(x)sum(is.na(x))/nrow(proportions_statistical_full)))
title("Proportion of NA per variable, full statistical dataset")

boxplot(apply(proportions_statistical_densified,2,function(x)sum(is.na(x))/nrow(proportions_statistical_densified)))
title("Proprotion of NA per variable, densified statistical dataset")

# perform PCA
proportions_original[proportions_original=="NaN"]<-NA
proportions_logical_full[proportions_logical_full=="NaN"]<-NA
proportions_logical_densified[proportions_logical_densified=="NaN"]<-NA
proportions_statistical_full[proportions_statistical_full=="NaN"]<-NA
proportions_statistical_densified[proportions_statistical_densified=="NaN"]<-NA

per_family_props_PCA_original <- pca(proportions_original, method='rnipals', nPcs=10)
per_family_props_PCA_logical_full <- pca(proportions_logical_full, method='rnipals', nPcs=10)
per_family_props_PCA_logical_densified <- pca(proportions_logical_densified, method='rnipals', nPcs=10)
per_family_props_PCA_statistical_full <- pca(proportions_statistical_full, method='rnipals', nPcs=10)
per_family_props_PCA_statistical_densified <- pca(proportions_statistical_densified, method='rnipals', nPcs=10)

# explained variances
explained_variance(per_family_props_PCA_logical_full)[1:9]
explained_variance(per_family_props_PCA_logical_densified)[1:9]
explained_variance(per_family_props_PCA_statistical_full)[1:9]
explained_variance(per_family_props_PCA_statistical_densified)[1:9]

# most relevant loadings 
report_loadings(per_family_props_PCA_logical_full, directory = "logicalGBM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_logical_full, directory = "logicalGBM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_logical_full, directory = "logicalGBM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_logical_densified, directory = "logicalGBM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_logical_densified, directory = "logicalGBM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_logical_densified, directory = "logicalGBM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_statistical_full, directory = "statisticalGBM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_statistical_full, directory = "statisticalGBM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_statistical_full, directory = "statisticalGBM", pc = 3, nr_variables_from_top_and_bottom = 5)

report_loadings(per_family_props_PCA_statistical_densified, directory = "statisticalGBM", pc = 1, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_statistical_densified, directory = "statisticalGBM", pc = 2, nr_variables_from_top_and_bottom = 5)
report_loadings(per_family_props_PCA_statistical_densified, directory = "statisticalGBM", pc = 3, nr_variables_from_top_and_bottom = 5)

# illustration of rgb cube
library(lattice)
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
plot(rgb_cube)

# map languages to colours
rgb_mapping_original <- RGB_mapping(per_family_props_PCA_original,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(original)) %>% na.omit()
rgb_mapping_logical_full <- RGB_mapping(per_family_props_PCA_logical_full,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(logical_full)) %>% na.omit()
rgb_mapping_logical_densified <- RGB_mapping(per_family_props_PCA_logical_densified,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(logical_densified)) %>% na.omit()
rgb_mapping_statistical_full <- RGB_mapping(per_family_props_PCA_statistical_full,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(statistical_full)) %>% na.omit()
rgb_mapping_statistical_densified <- RGB_mapping(per_family_props_PCA_statistical_densified,taxonomy_matrix_for_plot) %>% filter(id%in%rownames(statistical_densified)) %>% na.omit()

# global PC maps for PCs 1-3

# generate world map, centred around pacific (lon = 155)
library(rnaturalearth)
lon=155
world_map <- ne_download(scale = 110, category = "physical", type = "land", returnclass = "sf")
minor_islands <- ne_download(scale = 10, category = "physical", type = "minor_islands", returnclass = "sf")
world_map_including_islands <- rbind(world_map,minor_islands)
pacific_centred_world_map <- shift_map(world_map_including_islands, lon=lon)

rgb_map_original <- rgb_map(rgb_mapping_original, pacific_centred_world_map, lon=lon, main="Original Grambank", plot=F)
rgb_map_logical_full <- rgb_map(rgb_mapping_logical_full, pacific_centred_world_map, lon=lon, main="Logical curation, full dataset with coordinates", plot=F)
rgb_map_logical_densified <- rgb_map(rgb_mapping_logical_densified, pacific_centred_world_map, lon=lon, main="Logical curation, densified dataset", plot=F)
rgb_map_statistical_full <- rgb_map(rgb_mapping_statistical_full, pacific_centred_world_map, lon=lon, main="Statistical curation, full dataset with coordinates", plot=F)
rgb_map_statistical_densified <- rgb_map(rgb_mapping_statistical_densified, pacific_centred_world_map, lon=lon, main="Statistical curation, densified dataset", plot=F)

# arrange the plots side by side
library(gridExtra)

combined_plot <- grid.arrange(
  arrangeGrob(rgb_map_original, ncol = 1),
  arrangeGrob(rgb_map_logical_full, rgb_map_logical_densified, ncol = 1),   
  arrangeGrob(rgb_map_statistical_full, rgb_map_statistical_densified, ncol = 1), 
  ncol = 3
)

# print and save the combined plot
ggsave("output/family_collapsed_rgb_pca_maps.png", plot = combined_plot, width = 14, height = 8, dpi = 300)

##################
## Grid points and entropies at grid points

# coordinate_scaling 20 yields 715 grid points if buffer_distance is 150'000, changing coordinate_scaling gives more or less grid points
coordinate_scaling <- 20  # 715 grid points if coordinate_scaling=20 and buffer_distance=150000
buffer_distance <- 150000  # defining buffer distance for map for some Pacific islands to be kept
lon <- 155

# trim taxonomy to appropriate languages
original_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(original)) %>% filter(!is.na(lon))
logical_full_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(logical_full)) %>% filter(!is.na(lon))
logical_densified_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(logical_densified)) %>% filter(!is.na(lon))
statistical_full_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(statistical_full)) %>% filter(!is.na(lon))
statistical_densified_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(statistical_densified)) %>% filter(!is.na(lon))

# trim data to appropriate languages
original_with_lg_coords <- original[which(rownames(original)%in%logical_full_taxonomy$id),]
logical_full_with_lg_coords <- logical_full[which(rownames(logical_full)%in%logical_full_taxonomy$id),]
logical_densified_with_lg_coords <- logical_densified[which(rownames(logical_densified)%in%logical_densified_taxonomy$id),]
statistical_full_with_lg_coords <- statistical_full[which(rownames(statistical_full)%in%statistical_full_taxonomy$id),]
statistical_densified_with_lg_coords <- statistical_densified[which(rownames(statistical_densified)%in%statistical_densified_taxonomy$id),]

original_gps <- gridpointwise_entropies(taxonomy_matrix = original_taxonomy,
                                            data_trimmed_to_lgs_with_coords = original_with_lg_coords, 
                                            world_map = world_map_including_islands, 
                                            lon = 155, 
                                            coordinate_scaling = coordinate_scaling, 
                                            buffer_distance = buffer_distance, 
                                            verbose = T, 
                                            directory = "",
                                            data_type = "original",
                                            robust = T, robust_errors = "with_robust_errors")

logical_full_gps <- gridpointwise_entropies(taxonomy_matrix = logical_full_taxonomy,
                                                   data_trimmed_to_lgs_with_coords = logical_full_with_lg_coords, 
                                                   world_map = world_map_including_islands, 
                                                   lon = 155, 
                                                   coordinate_scaling = coordinate_scaling, 
                                                   buffer_distance = buffer_distance, 
                                                   verbose = T, 
                                                   directory = "logicalGBM",
                                                   data_type = "full",
                                                   robust = T, robust_errors = "with_robust_errors")

logical_densified_gps <- gridpointwise_entropies(taxonomy_matrix = logical_densified_taxonomy,
                                                 data_trimmed_to_lgs_with_coords = logical_densified_with_lg_coords, 
                                                 world_map = world_map_including_islands, 
                                                 lon = 155, 
                                                 coordinate_scaling = coordinate_scaling, 
                                                 buffer_distance = buffer_distance, 
                                                 verbose = T, 
                                                 directory = "logicalGBM",
                                                 data_type = "densified",
                                                 robust = T, robust_errors = "with_robust_errors")

statistical_full_gps <- gridpointwise_entropies(taxonomy_matrix = statistical_full_taxonomy,
                                                data_trimmed_to_lgs_with_coords = statistical_full_with_lg_coords, 
                                                world_map = world_map_including_islands, 
                                                lon = 155, 
                                                coordinate_scaling = coordinate_scaling, 
                                                buffer_distance = buffer_distance, 
                                                verbose = T, 
                                                directory = "statisticalGBM",
                                                data_type = "full",
                                                robust = T, robust_errors = "with_robust_errors")

statistical_densified_gps <-gridpointwise_entropies(taxonomy_matrix = statistical_densified_taxonomy,
                                                    data_trimmed_to_lgs_with_coords = statistical_densified_with_lg_coords, 
                                                    world_map = world_map_including_islands, 
                                                    lon = 155, 
                                                    coordinate_scaling = coordinate_scaling, 
                                                    buffer_distance = buffer_distance, 
                                                    verbose = T, 
                                                    directory = "statisticalGBM",
                                                    data_type = "densified",
                                                    robust = T, robust_errors = "with_robust_errors")

# compare curations (here: original vs. logical vs. statistical): we expect there to be systematic differences
# we do not plot lg/fam comparisons here, because we know nothing changes
grid_point_comparison_horizontal(baseline_gp_frame = original_gps, 
                                 name_baseline = "original Grambank",
                                 comparison_gp_frame = logical_full_gps, 
                                 name_comparison = "logical curation (full)",
                                 comparison_gp_frame_2 = statistical_full_gps, 
                                 name_comparison_2 = "statistical curation (full)",
                                 lg_and_fam_comparison = FALSE, fig_width = 15, fig_height = 10,
                                 name_for_plot = "_original_logical_statistical",
                                 world_map = world_map_including_islands, lon = 155)


# compare densified to full: we expect there to be no systematic differences
# logical dataset
grid_point_comparison_vertical(baseline_gp_frame = logical_full_gps, 
                               name_baseline = "full dataset (logical curation)", 
                               comparison_gp_frame = logical_densified_gps,
                               name_comparison = "densified dataset (logical curation)", 
                               world_map = world_map_including_islands, lon = 155,
                               directory = "logicalGBM")

# statistical dataset
grid_point_comparison_vertical(baseline_gp_frame = statistical_full_gps, 
                               name_baseline = "full dataset (statistical curation)", 
                               comparison_gp_frame = statistical_densified_gps,
                               name_comparison = "densified dataset (statistical curation)", 
                               world_map = world_map_including_islands, lon = 155,
                               directory = "statisticalGBM")
