# code for Fig 10
rm(list=ls())

library(densify)
library(tidyverse)
library(rnaturalearth)
library(cowplot)

source("scripts/functions.R")

## Retrieve taxonomy, coordinates and map

# generate taxonomy
taxonomy_matrix <- as_flat_taxonomy_matrix(glottolog_languoids)

# link coordinates to lgs
macroareas <- read.csv("scripts/lang-metadata.csv")
names(macroareas)[1] <- "id"
macroareas <- na_convert(macroareas)
taxonomy_matrix_for_plot <- merge(taxonomy_matrix,select(macroareas, c("id","latitude","longitude","glottolog.macroarea")),by="id") %>% filter(!is.na(latitude))
names(taxonomy_matrix_for_plot)[29:30] <- c("lat","lon") # change names to lat and lon for downstream functions

# download Natural Earth data for sf, including islands
world_map <- ne_download(scale = 110, category = "physical", type = "land", returnclass = "sf")
minor_islands <- ne_download(scale = 10, category = "physical", type = "minor_islands", returnclass = "sf")
world_map_initial <- rbind(world_map, minor_islands)

## Grambank and GBI data: suitable given high coding densities and high number of languages
# read in matrices (original, logical and statistical), convert all NA and ? to NA, factorize all variables
gb_original <- read.csv("raw_data/grambank_v.1.0.3.csv", row.names = "Glottocode")
gb_original <- na_convert(gb_original)
gb_original <- factorise(gb_original)
sum(!is.na(gb_original))/(nrow(gb_original)*ncol(gb_original)) # coding density: 75%

gbi_logical <- read.csv("curated_data/GBI/logicalGBI/logicalGBI.csv",row.names = "glottocode") %>% select(-X)
gbi_logical <- na_convert(gbi_logical)
gbi_logical <- factorise(gbi_logical)
gbi_logical <- gbi_logical[which(rownames(gbi_logical)%in%taxonomy_matrix_for_plot$id),]
sum(!is.na(gbi_logical))/(nrow(gbi_logical)*ncol(gbi_logical)) # coding density: 62%

gbi_statistical <- read.csv("curated_data/GBI/statisticalGBI/statisticalGBI.csv",row.names = "glottocode") %>% select(-X)
gbi_statistical <- na_convert(gbi_statistical)
gbi_statistical <- factorise(gbi_statistical)
gbi_statistical <- gbi_statistical[which(rownames(gbi_statistical)%in%taxonomy_matrix_for_plot$id),]
sum(!is.na(gbi_statistical))/(nrow(gbi_statistical)*ncol(gbi_statistical)) # coding density: 58%

## Grid points and entropies at grid points

# trim taxonomy to appropriate languages
gb_original_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(gb_original)) %>% filter(!is.na(lon))
gbi_logical_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(gbi_logical)) %>% filter(!is.na(lon))
gbi_statistical_taxonomy <- taxonomy_matrix_for_plot %>% filter(id %in% rownames(gbi_statistical)) %>% filter(!is.na(lon))

# trim data to appropriate languages
gb_original_with_lg_coords <- gb_original[which(rownames(gb_original)%in%gbi_logical_taxonomy$id),]
gbi_logical_with_lg_coords <- gbi_logical[which(rownames(gbi_logical)%in%gbi_logical_taxonomy$id),]
gbi_statistical_with_lg_coords <- gbi_statistical[which(rownames(gbi_statistical)%in%gbi_statistical_taxonomy$id),]

# define function that create a network of grid points, associated to a map and a taxonomy matrix of coordinates, and maps mean feature entropies per grid point
gridpointwise_entropies <- function(taxonomy_matrix, 
                                    data_trimmed_to_lgs_with_coords, 
                                    world_map, 
                                    coordinate_scaling, 
                                    buffer_distance, 
                                    verbose=T, 
                                    plot=T, title=""){
  library(geosphere)
  library(sp)
  library(sf)
  library(viridis)
  library(ggplot2)
  library(posterior)
  library(ggfortify)
  library(gridExtra)
  library(RColorBrewer)
  library(lmtest)
  library(sandwich)
  library(jtools)
  library(modelsummary)
  library(lmls)
  
  # generate grid points using the regularCoordinates() function from the geosphere package - the function solves the Thompson Problem, cf. Derungs et al. 2018
  reg_coords <- regularCoordinates(coordinate_scaling) %>% as.data.frame()
  
  # shift all (coordinates, languages, maps)
  shifted <- project_data(df = reg_coords, world_map_initial = world_map_initial)
  shifted2 <- project_data(df = taxonomy_matrix, world_map_initial = world_map_initial)
  shifted$language_locations <- shifted2$data
  
  # intersect regular coordinates with continental polygons, including a buffer to keep islands
  shifted$base_map <- st_buffer(shifted$base_map, dist = buffer_distance)
  overlap <- st_intersection(shifted$data, shifted$base_map)
  
  # there are some duplicates that have to be removed
  coordinates <- st_coordinates(overlap)
  duplicates <- duplicated(coordinates)
  overlap <- overlap[!duplicates, ]
  overlap$grid_point <- 1:nrow(overlap)
  
  # now compute distances between grid points
  gp_dist <- st_distance(overlap) # compute distance between each grid point
  diag(gp_dist) <- NA # diagonal is NA
  gp_dist_min <- apply(gp_dist,2,function(x){min(x,na.rm=T)}) # retrieve distance between each grid point and its nearest neighbour
  
  # compute distance between grid points and all languages
  language_gp_dist <- st_distance(overlap, shifted$language_locations) # compute distance between each language and each grid point 
  colnames(language_gp_dist) <- shifted$language_locations$id
  rownames(language_gp_dist) <- c(1:nrow(language_gp_dist))
  language_gp_dist_min <- apply(language_gp_dist,2,min) # compute distance for each language to its nearest neighbour grid point
  
  # determine closest grid point for each language
  nearest_grid_point_coord <- data.frame(glottocode=shifted$language_locations$id,
                                         grid_point=apply(language_gp_dist,2,function(x)which(x==min(x))))
  
  # iterate through all grid points and count nearest languages and families
  nearest_lgs_and_fams_coord <- data.frame(grid_point=1:nrow(overlap),  mean.feature.entropy=NA)
  
  #iteration over grid points
  for(i in 1:nrow(nearest_lgs_and_fams_coord)){
    if(verbose==T){
      print(i)
    }
    # get nearest languages and families, given the distance
    lgs <- filter(nearest_grid_point_coord,grid_point==i)$glottocode

    # get gridpoint-wise entropies for all features and log the arithmetic mean
    data_at_gp <- data_trimmed_to_lgs_with_coords[rownames(data_trimmed_to_lgs_with_coords)%in%lgs,]
    
    featurewise_entropies <- apply(data_at_gp,2,function(x)
      if(length(na.omit(x))<=0.5*length(x) | length(x)==1) {NA} # if a grid point only is associated to one language, no entropy should be computed. entropies should also not be computed for features coded for under 50% of languages at the grid point 
      else{posterior::entropy(na.omit(x))}) # this is the normalized entropy, computed for those features coded for over 50% of languages at a grid point

    nearest_lgs_and_fams_coord$mean.feature.entropy[i] <- featurewise_entropies %>% na.omit() %>% mean() # compute and store mean across features
  }
  
  # add data to grid points
  grid_points_with_data <- merge(overlap, nearest_lgs_and_fams_coord)
  
  # create plot if required:
  if(plot){
    # this would show just the language locations with respect to the grid points
    # raw_plot <- ggplot() +
    #   geom_sf(data = shifted2$base_map, fill = "white", color = "darkgrey") +
    #   geom_sf(data = grid_points_with_data) +
    #   geom_sf(data = shifted$language_locations, color = "red", shape = 1, size = 0.5) +
    #   ggtitle("") +
    #   theme_minimal() +
    #   theme(panel.grid.major = element_blank())
    
    entropy_at_point_plot <- ggplot() +
      geom_sf(data = shifted2$base_map, fill = "white", color = "darkgrey") +
      # geom_sf(data = shifted$language_locations, color = "darkgreen", shape = 1, size = 0.2, alpha = 0.3) + # these are the language locations
      geom_sf(data = grid_points_with_data, aes(colour = mean.feature.entropy)) +
      scale_color_viridis_c(option = "magma", direction = -1, na.value = "lightgrey", limits = c(0,1)) + 
      ggtitle(title) +
      labs(color = "Mean entropy") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            panel.grid.major = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    
    entropy_at_point_plot
    
    return(list(grid_points_with_data, entropy_at_point_plot))
  } else{
    return(list(grid_points_with_data, NULL))
  }
}

# define function that compares grid point plots across curations
grid_point_comparison_horizontal <- function(baseline_gp_frame, name_baseline, 
                                             comparison_gp_frame, 
                                             comparison_gp_frame_2=NULL, 
                                             world_map){
  
  # shift world map
  shifted <- project_data(df = baseline_gp_frame, world_map_initial = world_map_initial)
  
  # new data frame
  delta <- comparison_gp_frame
  
  # compute delta of mean entropy across features, at each grid point
  delta$delta.1.minus.0.mean.feature.entropy <- comparison_gp_frame$mean.feature.entropy-baseline_gp_frame$mean.feature.entropy
  
  # if there are three curations to compare (grambank):
  if(!is.null(comparison_gp_frame_2)){
    delta$delta.2.minus.1.mean.feature.entropy <- comparison_gp_frame_2$mean.feature.entropy-comparison_gp_frame$mean.feature.entropy
    delta$delta.2.minus.0.mean.feature.entropy <- comparison_gp_frame_2$mean.feature.entropy-baseline_gp_frame$mean.feature.entropy
  }
  
  # min and max values are straightforward to compute if only two-way comparison requested (computed only for plots to work in any case)
  max_delta_entropy <- max(na.omit(delta$delta.1.minus.0.mean.feature.entropy))
  min_delta_entropy <- min(na.omit(delta$delta.1.minus.0.mean.feature.entropy))
  
  # determine minimum and maximum deltas for each type of comparison, if there are three curations to compare
  if(!is.null(comparison_gp_frame_2)){
    max_delta_entropy <- max(max(na.omit(delta$delta.1.minus.0.mean.feature.entropy)), max(na.omit(delta$delta.2.minus.1.mean.feature.entropy)), max(na.omit(delta$delta.2.minus.0.mean.feature.entropy)))
    min_delta_entropy <- min(min(na.omit(delta$delta.1.minus.0.mean.feature.entropy)), min(na.omit(delta$delta.2.minus.1.mean.feature.entropy)), min(na.omit(delta$delta.2.minus.0.mean.feature.entropy)))
  }
  
  boundaries <- sort(c(max_delta_entropy, 0-min_delta_entropy))[2]
  
  col_strip <- brewer.pal(11, "RdBu")

  # delta entropies maps
  delta_entropies_1_minus_0 <- ggplot() +
    geom_sf(data = shifted$base_map, fill = "white", color = "darkgrey") +
    geom_sf(data = delta, aes(colour = delta.1.minus.0.mean.feature.entropy)) +
    scale_color_gradientn(colors = rev(col_strip), na.value = "lightgrey", limits = c(0-boundaries, boundaries)) +
    labs(color = "Difference") +
    ggtitle("b.") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.major = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  if(!is.null(comparison_gp_frame_2)){
    delta_entropies_2_minus_1 <- ggplot() +
      geom_sf(data = shifted$base_map, fill = "white", color = "darkgrey") +
      geom_sf(data = delta, aes(colour = delta.2.minus.1.mean.feature.entropy)) +
      scale_color_gradientn(colors = rev(col_strip), na.value = "lightgrey", limits = c(0-boundaries, boundaries)) +
      labs(color = "Difference") +
      ggtitle("c.") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            panel.grid.major = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    
    delta_entropies_2_minus_0 <- ggplot() +
      geom_sf(data = shifted$base_map, fill = "white", color = "darkgrey") +
      geom_sf(data = delta, aes(colour = delta.2.minus.0.mean.feature.entropy)) +
      scale_color_gradientn(colors = rev(col_strip), na.value = "lightgrey", limits = c(0-boundaries, boundaries)) +
      labs(color = "Difference") +
      ggtitle("d.") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            panel.grid.major = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())  }
  
  # create combined plot according to specifications
  if(!is.null(comparison_gp_frame_2)){
      combined_plot <- plot_grid(delta_entropies_1_minus_0, delta_entropies_2_minus_1,delta_entropies_2_minus_0, 
                                 align = "h", axis = "l", nrow = 3, rel_heights = c(1/3,1/3,1/3))
      return(list(combined_plot, delta_entropies_1_minus_0, delta_entropies_2_minus_1, delta_entropies_2_minus_0))
      
  }else{
      combined_plot <- grid.arrange(
        arrangeGrob(delta_entropies_1_minus_0, ncol = 1),
        ncol = 1
      )
      return(list(combined_plot, delta_entropies_1_minus_0, NULL, NULL))
  }
}


## maps for Grambank and GBI:
# coordinate_scaling 20 yields ~750 grid points if buffer_distance is 150'000, changing coordinate_scaling gives more or less grid points
coordinate_scaling <- 20  # ~750 grid points if coordinate_scaling=20 and buffer_distance=150000
buffer_distance <- 150000  # defining buffer distance for map for some Pacific islands to be kept

gb_original_gps <- gridpointwise_entropies(taxonomy_matrix = gb_original_taxonomy,
                                        data_trimmed_to_lgs_with_coords = gb_original_with_lg_coords, 
                                        world_map = world_map_initial, 
                                        coordinate_scaling = coordinate_scaling, 
                                        buffer_distance = buffer_distance, 
                                        verbose = T, plot = F, title = "")

gbi_logical_gps <- gridpointwise_entropies(taxonomy_matrix = gbi_logical_taxonomy,
                                       data_trimmed_to_lgs_with_coords = gbi_logical_with_lg_coords, 
                                       world_map = world_map_initial, 
                                       coordinate_scaling = coordinate_scaling, 
                                       buffer_distance = buffer_distance, 
                                       verbose = T, plot = F, title = "")

gbi_statistical_gps <- gridpointwise_entropies(taxonomy_matrix = gbi_statistical_taxonomy,
                                           data_trimmed_to_lgs_with_coords = gbi_statistical_with_lg_coords, 
                                           world_map = world_map_initial, 
                                           coordinate_scaling = coordinate_scaling, 
                                           buffer_distance = buffer_distance, 
                                           verbose = T, plot = T, title = "a.")

# compare across curations
gbi_comparison_plot <- grid_point_comparison_horizontal(baseline_gp_frame = gb_original_gps[[1]], 
                                                        comparison_gp_frame = gbi_logical_gps[[1]], 
                                                        comparison_gp_frame_2 = gbi_statistical_gps[[1]], 
                                                        world_map = world_map_initial)
# plot
entropies <- plot_grid(gbi_statistical_gps[[2]],
          gbi_comparison_plot[[2]],
          gbi_comparison_plot[[3]],
          gbi_comparison_plot[[4]], ncol = 1)


ggsave(file="plots/Fig10.pdf", entropies, dpi = 500, width = 8, height = 9, units = "in")


