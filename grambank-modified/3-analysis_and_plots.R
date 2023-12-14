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

##################
## PCA, aggregated to family level
proportions_original <- generate_per_family_prop_table(data = original, taxonomy_matrix = taxonomy_matrix)
proportions_logical_full <- generate_per_family_prop_table(data = logical_full, taxonomy_matrix = taxonomy_matrix)
proportions_logical_densified <- generate_per_family_prop_table(data = logical_densified, taxonomy_matrix = taxonomy_matrix)
proportions_statistical_full <- generate_per_family_prop_table(data = statistical_full, taxonomy_matrix = taxonomy_matrix)
proportions_statistical_densified <- generate_per_family_prop_table(data = statistical_densified, taxonomy_matrix = taxonomy_matrix)

# # proportions of NA per family
# boxplot(apply(proportions_logical_full,1,function(x)sum(is.na(x))/nrow(proportions_logical_full)))
# title("Proportion of NA per family, full logical dataset")

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

# explained variances and R2 among all PCs
explained_variance(per_family_props_PCA_logical_densified)[,1:10]
explained_variance(per_family_props_PCA_logical_densified)[,1:10]
explained_variance(per_family_props_PCA_statistical_full)[,1:10]
explained_variance(per_family_props_PCA_statistical_densified)[,1:10]

# plot(per_family_props_PCA_original,choix="var",new.plot=FALSE)

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

### download Natural Earth data for sf, including islands
if(!exists("world_map_initial")){
  library(rnaturalearth)
  world_map <- ne_download(scale = 110, category = "physical", type = "land", returnclass = "sf")
  minor_islands <- ne_download(scale = 10, category = "physical", type = "minor_islands", returnclass = "sf")
  world_map_initial <- rbind(world_map, minor_islands)
}

rgb_map_original <- rgb_map(rgb_mapping_original, world_map_initial, main="Original Grambank", plot=F)
rgb_map_logical_full <- rgb_map(rgb_mapping_logical_full, world_map_initial, main="Grambank, logical curation, full dataset", plot=F)
rgb_map_logical_densified <- rgb_map(rgb_mapping_logical_densified, world_map_initial, main="Grambank, logical curation, densified dataset", plot=F)
rgb_map_statistical_full <- rgb_map(rgb_mapping_statistical_full, world_map_initial, main="Grambank, statistical curation, full dataset", plot=F)
rgb_map_statistical_densified <- rgb_map(rgb_mapping_statistical_densified, world_map_initial, main="Grambank, statistical curation, densified dataset", plot=F)

#### family relationships to one another in PC-space, by macroarea
f2d_original <- dim_mapping(pca.object = per_family_props_PCA_original, taxonomy_matrix = taxonomy_matrix_for_plot)
f2d_logical_full <- dim_mapping(pca.object = per_family_props_PCA_logical_full, taxonomy_matrix = taxonomy_matrix_for_plot)
f2d_logical_densified <- dim_mapping(pca.object = per_family_props_PCA_logical_densified, taxonomy_matrix = taxonomy_matrix_for_plot)
f2d_statistical_full <- dim_mapping(pca.object = per_family_props_PCA_statistical_full, taxonomy_matrix = taxonomy_matrix_for_plot)
f2d_statistical_densified <- dim_mapping(pca.object = per_family_props_PCA_statistical_densified, taxonomy_matrix = taxonomy_matrix_for_plot)

# arrange the plots side by side, print and plot
logical_plot <- grid.arrange(arrangeGrob(rgb_map_logical_full, f2d_logical_full, rgb_map_logical_densified, f2d_logical_densified, nrow = 2, ncol = 2))
ggsave("../figures-for-paper/main/logical_family_collapsed_rgb_pca_maps_grambank.png", plot = logical_plot, width = 14, height = 8, dpi = 300)

combined_plot <- grid.arrange(
  arrangeGrob(rgb_map_original, f2d_original, ncol = 2),
  arrangeGrob(rgb_map_logical_full, f2d_logical_full, rgb_map_logical_densified, f2d_logical_densified, nrow = 2, ncol = 2),   
  arrangeGrob(rgb_map_statistical_full, f2d_statistical_full, rgb_map_statistical_densified, f2d_statistical_densified, ncol = 2, nrow = 2), 
  ncol = 3
)
ggsave("../figures-for-paper/supplementary/all_family_collapsed_rgb_pca_maps_grambank.png", plot = combined_plot, width = 14*3, height = 8, dpi = 300)

##################
## Grid points and entropies at grid points

# coordinate_scaling 20 yields ~750 grid points if buffer_distance is 150'000, changing coordinate_scaling gives more or less grid points
coordinate_scaling <- 20  # ~750 grid points if coordinate_scaling=20 and buffer_distance=150000
buffer_distance <- 150000  # defining buffer distance for map for some Pacific islands to be kept

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
                                            world_map = world_map_initial, 
                                            coordinate_scaling = coordinate_scaling, 
                                            buffer_distance = buffer_distance, 
                                            verbose = T, 
                                            main_fig = F,
                                            data_type = "Grambank (original)")

logical_full_gps <- gridpointwise_entropies(taxonomy_matrix = logical_full_taxonomy,
                                                   data_trimmed_to_lgs_with_coords = logical_full_with_lg_coords, 
                                                   world_map = world_map_initial, 
                                                   coordinate_scaling = coordinate_scaling, 
                                                   buffer_distance = buffer_distance, 
                                                   verbose = T, 
                                                  main_fig = F,
                                                   data_type = "Grambank (logical)")

logical_densified_gps <- gridpointwise_entropies(taxonomy_matrix = logical_densified_taxonomy,
                                                 data_trimmed_to_lgs_with_coords = logical_densified_with_lg_coords, 
                                                 world_map = world_map_initial, 
                                                 coordinate_scaling = coordinate_scaling, 
                                                 buffer_distance = buffer_distance, 
                                                 verbose = T, 
                                                 main_fig = F,
                                                 data_type = "Grambank (logical, densified)")

statistical_full_gps <- gridpointwise_entropies(taxonomy_matrix = statistical_full_taxonomy,
                                                data_trimmed_to_lgs_with_coords = statistical_full_with_lg_coords, 
                                                world_map = world_map_initial, 
                                                coordinate_scaling = coordinate_scaling, 
                                                buffer_distance = buffer_distance, 
                                                verbose = T, 
                                                main_fig = T,
                                                data_type = "Grambank (statistical)")

statistical_densified_gps <-gridpointwise_entropies(taxonomy_matrix = statistical_densified_taxonomy,
                                                    data_trimmed_to_lgs_with_coords = statistical_densified_with_lg_coords, 
                                                    world_map = world_map_initial, 
                                                    coordinate_scaling = coordinate_scaling, 
                                                    buffer_distance = buffer_distance, 
                                                    verbose = T, 
                                                    main_fig = F,
                                                    data_type = "Grambank (statistical, densified)")

# compare curations (here: original vs. logical vs. statistical): we expect there to be systematic differences
# we do not plot lg/fam comparisons here, because we know nothing changes
grid_point_comparison_horizontal(baseline_gp_frame = original_gps, 
                                 name_baseline = "original Grambank",
                                 comparison_gp_frame = logical_full_gps, 
                                 name_comparison = "logical curation (full)",
                                 comparison_gp_frame_2 = statistical_full_gps, 
                                 name_comparison_2 = "statistical curation (full)",
                                 lg_and_fam_comparison = FALSE, fig_width = 10, fig_height = 10,
                                 name_for_plot = "grambank_original_logical_statistical",
                                 world_map = world_map_initial)


# compare densified to full: we expect there to be no systematic differences
# logical dataset
grid_point_comparison_vertical(baseline_gp_frame = logical_full_gps, 
                               name_baseline = "full dataset (logical curation)", 
                               comparison_gp_frame = logical_densified_gps,
                               name_comparison = "densified dataset (logical curation)", 
                               world_map = world_map_initial,
                               directory = "logicalGBM")

# statistical dataset
grid_point_comparison_vertical(baseline_gp_frame = statistical_full_gps, 
                               name_baseline = "full dataset (statistical curation)", 
                               comparison_gp_frame = statistical_densified_gps,
                               name_comparison = "densified dataset (statistical curation)", 
                               world_map = world_map_initial,
                               directory = "statisticalGBM")


# compute gower distances, using daisy()
statistical_full_gower <- daisy(statistical_full_with_lg_coords, metric = "gower") %>% as.matrix()
statistical_densified_gower <- daisy(statistical_densified_with_lg_coords, metric = "gower") %>% as.matrix()

library(vegan)
set.seed(78)
# statistical_full_metaMDS <- metaMDS(statistical_full_gower, 
#                                     k = 3, 
#                                     try = 20, trymax = 500,
#                                     engine = "monoMDS")

statistical_densified_metaMDS <- metaMDS(statistical_densified_gower, 
                                    k = 3, 
                                    try = 20, trymax = 500,
                                    engine = "monoMDS")

##### genetics vs typology in Eurasia
gelato <- read.csv("../gelato/GELATO_population_glottocode_mapping.csv")

gelato_pops <- filter(gelato,!is.na(GB_statistical_full_proxy)) %>% select(c(PopName,Location,country,lat,lon,glottocodeBase,glottologFamily,comment_GB_proxy,GB_statistical_full_proxy,proxy_gb_statistical_full_density))
gelato_pops <- left_join(gelato_pops,select(statistical_full_taxonomy,"id","macroarea"),by=c("GB_statistical_full_proxy" = "id"))
gelato_pops <- gelato_pops %>% filter(!is.na(lat)) %>% filter(macroarea == "Eurasia") 

lgs_in_gelato <- statistical_full[which(rownames(statistical_full)%in%gelato_pops$GB_statistical_full_proxy),] 

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

gb_eurasia_all <- genetics_vs_typology(gelato_pops = gelato_pops,
                     lgs_in_gelato = lgs_in_gelato,
                     fst_matrix = fst_matrix,
                     seed=7,
                     world_map = world_map_initial, 
                     coordinate_scaling = coordinate_scaling, 
                     buffer_distance = buffer_distance,
                     verbose=T, 
                     data_type = "grambank_eurasia")

# repeat procedure without Onge, because they are outliers
gelato_pops_no_onge <- gelato_pops %>% filter(PopName != "Onge")
lgs_in_gelato_no_onge <- lgs_in_gelato[which(rownames(lgs_in_gelato) != "onge1236"),]
fst_matrix_no_onge <- fst_matrix[which(rownames(fst_matrix)!="Onge"),which(colnames(fst_matrix)!="Onge")]

gb_eurasia_all_no_onge <- genetics_vs_typology(gelato_pops = gelato_pops_no_onge,
                                       lgs_in_gelato = lgs_in_gelato_no_onge,
                                       fst_matrix = fst_matrix_no_onge,
                                       seed=7,
                                       world_map = world_map_initial, 
                                       coordinate_scaling = coordinate_scaling, 
                                       buffer_distance = buffer_distance,
                                       verbose=T, 
                                       data_type = "grambank_eurasia_no_onge")
