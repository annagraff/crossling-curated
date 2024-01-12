# this script serves to validate the automated columns (columns 21 to 73) in the file GELATO_population_glottocode_mapping.csv
# note that columns 1 to 20 are asscociated with GeLaTo [perpopMASTER_653pops2022.xlsx] (columns 1 to 8) or manually entered (columns 9 to 20)

rm(list=ls())

library(tidyverse)
library(testthat)
library(GetoptLong)
library(densify)
library(readxl)

source("../functions.R")

# read in glottolog taxonomy
taxonomy <- as_flat_taxonomy_matrix(glottolog_languoids)
names(taxonomy)[1]<-"glottocode"

glottolog_level <- read.csv("../input/glottolog_v.4.8/languages_and_dialects_geo.csv")

# however note that in gelato we still have a location column...
gelato_mapping_decisions <- read.csv("GELATO_population_glottocode_mapping.csv")
gelato_mapping_decisions <- na_convert(gelato_mapping_decisions)

### check all densities are correct
gelato <- gelato_mapping_decisions %>% select(1:20)

# retrieve proxies
proxies <- gelato[,c(11:20)]
proxies <- na_convert(proxies)

# sanity checks
# expect_true(all(perpopMASTER$PopName==gelato$PopName))
# expect_true(all(perpopMASTER$glottocodeBase==gelato$glottocodeBase))
# expect_true(all(perpopMASTER$glottologFamily==gelato$glottologFamily)) 

base_status <- apply(as.data.frame(gelato),1,function(x)if(x[7]%in%glottolog_level$glottocode==F){NA}else{filter(glottolog_level,glottocode==x[7])[1,4]})
upstream_glottocode <- apply(as.data.frame(gelato),1,function(x)if(x[7]%in%taxonomy$glottocode==F){NA}else{filter(taxonomy,glottocode==x[7])[1,sum(!is.na(filter(taxonomy,glottocode==x[7])))-1]})
upstream_glottocode_status <- apply(as.data.frame(upstream_glottocode),1,function(x)if(x%in%glottolog_level$glottocode==F){NA}else{filter(glottolog_level,glottocode==x[1])[1,4]})

gelato <- cbind(gelato,data.frame(status_glottolog = base_status,
                                  glottocodeUpstream = upstream_glottocode,
                                  glottocodeUpstream_status = upstream_glottocode_status))

# read in GB data: logical and statistical, full and densified
gblogfull <- read.csv("../grambank-modified/output/logicalGBM/logicalGBM.csv", row.names = "glottocode") %>% select(-X)
gblogdens <- read.csv("../grambank-modified/output/logicalGBM/logicalGBM_pruned.csv", row.names = "X")

gbstafull <- read.csv("../grambank-modified/output/statisticalGBM/statisticalGBM.csv", row.names = "glottocode") %>% select(-X)
gbstadens <- read.csv("../grambank-modified/output/statisticalGBM/statisticalGBM_pruned.csv", row.names = "X")

# read in MM data: logical and statistical, full, large densified and small densified 
mmlogfull <- read.csv("../mega-merge/output/logicalMM/full/logicalMM_full.csv", row.names = "glottocode") %>% select(-X)
mmlogl <- read.csv("../mega-merge/output/logicalMM/full/logicalMM_full_pruned_large.csv", row.names = "X")
mmlogs <- read.csv("../mega-merge/output/logicalMM/full/logicalMM_full_pruned_small.csv", row.names = "X")

mmstafull <- read.csv("../mega-merge/output/statisticalMM/full/statisticalMM_full.csv", row.names = "glottocode") %>% select(-X)
mmstal <- read.csv("../mega-merge/output/statisticalMM/full/statisticalMM_full_pruned_large.csv", row.names = "X")
mmstas <- read.csv("../mega-merge/output/statisticalMM/full/statisticalMM_full_pruned_small.csv", row.names = "X")

# turn all ? or missing data to NA
gblogfull <- na_convert(gblogfull)
gblogdens <- na_convert(gblogdens)
gbstafull <- na_convert(gbstafull)
gbstadens <- na_convert(gbstadens)

mmlogfull <- na_convert(mmlogfull)
mmlogl <- na_convert(mmlogl)
mmlogs <- na_convert(mmlogs)
mmstafull <- na_convert(mmstafull)
mmstal <- na_convert(mmstal)
mmstas <- na_convert(mmstas)

# which languages are in any of the GB and/or MM curations?
glottocodesInAnyGB <- unique(c(rownames(gblogfull),rownames(gblogdens),rownames(gbstafull),rownames(gbstadens)))
glottocodesInAnyMM <- unique(c(rownames(mmlogfull),rownames(mmlogl),rownames(mmlogs),rownames(mmstafull),rownames(mmstal),rownames(mmstas)))

# variable counts for each dataset
gblogfull.count <- apply(as.data.frame(glottocodesInAnyGB),1,function(x)if(x%in%rownames(gblogfull)==F){NA}else{length(na.omit(t(gblogfull[which(rownames(gblogfull)==x),])))})
gblogdens.count <- apply(as.data.frame(glottocodesInAnyGB),1,function(x)if(x%in%rownames(gblogdens)==F){NA}else{length(na.omit(t(gblogdens[which(rownames(gblogdens)==x),])))})
gbstafull.count <- apply(as.data.frame(glottocodesInAnyGB),1,function(x)if(x%in%rownames(gbstafull)==F){NA}else{length(na.omit(t(gbstafull[which(rownames(gbstafull)==x),])))})
gbstadens.count <- apply(as.data.frame(glottocodesInAnyGB),1,function(x)if(x%in%rownames(gbstadens)==F){NA}else{length(na.omit(t(gbstadens[which(rownames(gbstadens)==x),])))})

mmlogfull.count <- apply(as.data.frame(glottocodesInAnyMM),1,function(x)if(x%in%rownames(mmlogfull)==F){NA}else{length(na.omit(t(mmlogfull[which(rownames(mmlogfull)==x),])))})
mmlogl.count <- apply(as.data.frame(glottocodesInAnyMM),1,function(x)if(x%in%rownames(mmlogl)==F){NA}else{length(na.omit(t(mmlogl[which(rownames(mmlogl)==x),])))})
mmlogs.count <- apply(as.data.frame(glottocodesInAnyMM),1,function(x)if(x%in%rownames(mmlogs)==F){NA}else{length(na.omit(t(mmlogs[which(rownames(mmlogs)==x),])))})
mmstafull.count <- apply(as.data.frame(glottocodesInAnyMM),1,function(x)if(x%in%rownames(mmstafull)==F){NA}else{length(na.omit(t(mmstafull[which(rownames(mmstafull)==x),])))})
mmstal.count <- apply(as.data.frame(glottocodesInAnyMM),1,function(x)if(x%in%rownames(mmstal)==F){NA}else{length(na.omit(t(mmstal[which(rownames(mmstal)==x),])))})
mmstas.count <- apply(as.data.frame(glottocodesInAnyMM),1,function(x)if(x%in%rownames(mmstas)==F){NA}else{length(na.omit(t(mmstas[which(rownames(mmstas)==x),])))})

# densities for all datasets
gblogfull.density <- apply(as.data.frame(glottocodesInAnyGB),1,function(x)if(x%in%rownames(gblogfull)==F){NA}else{length(na.omit(t(gblogfull[which(rownames(gblogfull)==x),])))/ncol(gblogfull)})
gblogdens.density <- apply(as.data.frame(glottocodesInAnyGB),1,function(x)if(x%in%rownames(gblogdens)==F){NA}else{length(na.omit(t(gblogdens[which(rownames(gblogdens)==x),])))/ncol(gblogdens)})
gbstafull.density <- apply(as.data.frame(glottocodesInAnyGB),1,function(x)if(x%in%rownames(gbstafull)==F){NA}else{length(na.omit(t(gbstafull[which(rownames(gbstafull)==x),])))/ncol(gbstafull)})
gbstadens.density <- apply(as.data.frame(glottocodesInAnyGB),1,function(x)if(x%in%rownames(gbstadens)==F){NA}else{length(na.omit(t(gbstadens[which(rownames(gbstadens)==x),])))/ncol(gbstadens)})

mmlogfull.density <- apply(as.data.frame(glottocodesInAnyMM),1,function(x)if(x%in%rownames(mmlogfull)==F){NA}else{length(na.omit(t(mmlogfull[which(rownames(mmlogfull)==x),])))/ncol(mmlogfull)})
mmlogl.density <- apply(as.data.frame(glottocodesInAnyMM),1,function(x)if(x%in%rownames(mmlogl)==F){NA}else{length(na.omit(t(mmlogl[which(rownames(mmlogl)==x),])))/ncol(mmlogl)})
mmlogs.density <- apply(as.data.frame(glottocodesInAnyMM),1,function(x)if(x%in%rownames(mmlogs)==F){NA}else{length(na.omit(t(mmlogs[which(rownames(mmlogs)==x),])))/ncol(mmlogs)})
mmstafull.density <- apply(as.data.frame(glottocodesInAnyMM),1,function(x)if(x%in%rownames(mmstafull)==F){NA}else{length(na.omit(t(mmstafull[which(rownames(mmstafull)==x),])))/ncol(mmstafull)})
mmstal.density <- apply(as.data.frame(glottocodesInAnyMM),1,function(x)if(x%in%rownames(mmstal)==F){NA}else{length(na.omit(t(mmstal[which(rownames(mmstal)==x),])))/ncol(mmstal)})
mmstas.density <- apply(as.data.frame(glottocodesInAnyMM),1,function(x)if(x%in%rownames(mmstas)==F){NA}else{length(na.omit(t(mmstas[which(rownames(mmstas)==x),])))/ncol(mmstas)})

# subset taxonomy to lgs that are in gelato or mm (each curation)
taxonomy_for_densities <- taxonomy
taxonomy_for_densities$GelatoBase <- apply(taxonomy,1,function(x)x[1]%in%gelato$glottocodeBase)
taxonomy_for_densities$GelatoUpstream <- apply(taxonomy,1,function(x)x[1]%in%gelato$glottocodeUpstream)
taxonomy_for_densities$InGB <- apply(taxonomy,1,function(x)x[1]%in%glottocodesInAnyGB)
taxonomy_for_densities$InMM <- apply(taxonomy,1,function(x)x[1]%in%glottocodesInAnyMM)
taxonomy_for_densities <- taxonomy_for_densities %>% filter(GelatoBase == T | GelatoUpstream == T | InGB == T | InMM == T)

lgs_in_any_GB <- data.frame(glottocode = glottocodesInAnyGB,
                            features_count_gb_logical_full = gblogfull.count,
                            features_count_gb_logical_densified = gblogdens.count,
                            features_count_gb_statistical_full = gbstafull.count,
                            features_count_gb_statistical_densified = gbstadens.count, 
                            density_gb_logical_full= gblogfull.density,
                            density_gb_logical_densified= gblogdens.density,
                            density_gb_statistical_full= gbstafull.density,
                            density_gb_statistical_densified = gbstadens.density)

lgs_in_any_MM <- data.frame(glottocode = glottocodesInAnyMM,
                            features_count_mm_logical_full = mmlogfull.count,
                            features_count_mm_logical_large = mmlogl.count,
                            features_count_mm_logical_small = mmlogs.count,
                            features_count_mm_statistical_full = mmstafull.count,
                            features_count_mm_statistical_large= mmstal.count, 
                            features_count_mm_statistical_small = mmstas.count, 
                            density_mm_logical_full= mmlogfull.density,
                            density_mm_logical_large= mmlogl.density,
                            density_mm_logical_small = mmlogs.density,                             
                            density_mm_statistical_full= mmstafull.density,
                            density_mm_statistical_large = mmstal.density,
                            density_mm_statistical_small = mmstas.density)
                            
taxonomy_for_densities <- left_join(taxonomy_for_densities,lgs_in_any_GB)
taxonomy_for_densities <- left_join(taxonomy_for_densities,lgs_in_any_MM)

## helper columns: counts and densities for all datasets for glottocodeBase
gelato <- left_join(gelato,select(filter(taxonomy_for_densities,GelatoBase == T),c(glottocode,
                                                                            features_count_gb_logical_full,
                                                                            features_count_gb_logical_densified,
                                                                            features_count_gb_statistical_full,
                                                                            features_count_gb_statistical_densified,
                                                                            features_count_mm_logical_full,
                                                                            features_count_mm_logical_large,
                                                                            features_count_mm_logical_small,
                                                                            features_count_mm_statistical_full,
                                                                            features_count_mm_statistical_large,
                                                                            features_count_mm_statistical_small,
                                                                            density_gb_logical_full,
                                                                            density_gb_logical_densified,
                                                                            density_gb_statistical_full,
                                                                            density_gb_statistical_densified,
                                                                            density_mm_logical_full,
                                                                            density_mm_logical_large,
                                                                            density_mm_logical_small,
                                                                            density_mm_statistical_full,
                                                                            density_mm_statistical_large,
                                                                            density_mm_statistical_small)), by=c("glottocodeBase"="glottocode"))


## helper columns: counts and densities for all datasets for glottocodeUpstream
gelato <- left_join(gelato,select(filter(taxonomy_for_densities,GelatoUpstream == T),c(glottocode,
                                                                                       features_count_gb_logical_full,
                                                                                       features_count_gb_logical_densified,
                                                                                       features_count_gb_statistical_full,
                                                                                       features_count_gb_statistical_densified,
                                                                                       features_count_mm_logical_full,
                                                                                       features_count_mm_logical_large,
                                                                                       features_count_mm_logical_small,
                                                                                       features_count_mm_statistical_full,
                                                                                       features_count_mm_statistical_large,
                                                                                       features_count_mm_statistical_small,
                                                                                       density_gb_logical_full,
                                                                                       density_gb_logical_densified,
                                                                                       density_gb_statistical_full,
                                                                                       density_gb_statistical_densified,
                                                                                       density_mm_logical_full,
                                                                                       density_mm_logical_large,
                                                                                       density_mm_logical_small,
                                                                                       density_mm_statistical_full,
                                                                                       density_mm_statistical_large,
                                                                                       density_mm_statistical_small)), by=c("glottocodeUpstream"="glottocode"))


names(gelato)[24:63]<-c("base_features_count_gb_logical_full","base_features_count_gb_logical_densified","base_features_count_gb_statistical_full","base_features_count_gb_statistical_full","base_features_count_mm_logical_full","base_features_count_mm_logical_large","base_features_count_mm_logical_small","base_features_count_mm_statistical_full","base_features_count_mm_statistical_large","base_features_count_mm_statistical_small",
                       "base_density_gb_logical_full","base_density_gb_logical_densified","base_density_gb_statistical_full","base_density_gb_statistical_densified","base_density_mm_logical_full","base_density_mm_logical_large","base_density_mm_logical_small","base_density_mm_statistical_full","base_density_mm_statistical_large","base_density_mm_statistical_small",
                       "upstream_features_count_gb_logical_full","upstream_features_count_gb_logical_densified","upstream_features_count_gb_statistical_full","upstream_features_count_gb_statistical_densified","upstream_features_count_mm_logical_full","upstream_features_count_mm_logical_large","upstream_features_count_mm_logical_small","upstream_features_count_mm_statistical_full","upstream_features_count_mm_statistical_large","upstream_features_count_mm_statistical_small",
                       "upstream_density_gb_logical_full","upstream_density_gb_logical_densified","upstream_density_gb_statistical_full","upstream_density_gb_statistical_densified","upstream_density_mm_logical_full","upstream_density_mm_logical_large","upstream_density_mm_logical_small","upstream_density_mm_statistical_full","upstream_density_mm_statistical_large","upstream_density_mm_statistical_small")


proxy_gb_logical_full_density <- unlist(apply(proxies[,1],1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density_gb_logical_full"]}))
proxy_gb_logical_densified_density <- unlist(apply(proxies[,2],1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density_gb_logical_densified"]}))
proxy_gb_statistical_full_density <- unlist(apply(proxies[,3],1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density_gb_statistical_full"]}))
proxy_gb_statistical_densified_density <- unlist(apply(proxies[,4],1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density_gb_statistical_densified"]}))
proxy_mm_logical_full_density <- unlist(apply(proxies[,5],1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density_mm_logical_full"]}))
proxy_mm_logical_large_density <- unlist(apply(proxies[,6],1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density_mm_logical_large"]}))
proxy_mm_logical_small_density <- unlist(apply(proxies[,7],1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density_mm_logical_small"]}))
proxy_mm_statistical_full_density <- unlist(apply(proxies[,8],1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density_mm_statistical_full"]}))
proxy_mm_statistical_large_density <- unlist(apply(proxies[,9],1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density_mm_statistical_large"]}))
proxy_mm_statistical_small_density <- unlist(apply(proxies[,10],1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density_mm_statistical_small"]}))

# add columns specifying coverage for all proxies
gelato <- cbind(gelato,proxy_gb_logical_full_density,proxy_gb_logical_densified_density,proxy_gb_statistical_full_density,proxy_gb_statistical_densified_density,
                proxy_mm_logical_full_density,proxy_mm_logical_large_density,proxy_mm_logical_small_density,
                proxy_mm_statistical_full_density,proxy_mm_statistical_large_density,proxy_mm_statistical_small_density)


gelato_automated_glottologStatus <- gelato[,21:23]
gelato_automated_numeric <- gelato[,24:ncol(gelato)]

# now check everything is correct:
expect_true(all(gelato_automated_glottologStatus==gelato_mapping_decisions[21:23], na.rm = T))
expect_true(all(round(gelato_automated_numeric,digits=5)==round(gelato_mapping_decisions[24:ncol(gelato_mapping_decisions)],digits=5), na.rm = T))

###### density plots
densities <- left_join(taxonomy_for_densities, select(glottolog_level,"glottocode","macroarea","latitude","longitude")) %>% filter(!is.na(latitude))

source("../functions.R")
shifted <- project_data(densities, world_map_initial = world_map_initial)


ggplot() +
  geom_sf(data = shifted$base_map, fill = "white", color = "darkgrey") +
  geom_sf(data = shifted$data, aes(colour = density_mm_statistical_large), size = 1) +
  scale_color_viridis(option = "magma", direction = -1, na.value = "lightgrey") + 
  ggtitle("Mega-Merge, statistical, densified (large)") +
  labs(color = "Coverage") +
  theme_minimal() +
  theme(panel.grid.major = element_blank())                                                                                                                                 

ggplot() +
  geom_sf(data = shifted$base_map, fill = "white", color = "darkgrey") +
  geom_sf(data = shifted$data, aes(colour = density_gb_statistical_full), size = 1) +
  scale_color_viridis(option = "magma", direction = -1, na.value = "lightgrey") + 
  ggtitle("GB, statistical, full") +
  labs(color = "Coverage") +
  theme_minimal() +
  theme(panel.grid.major = element_blank())                                                                                                                                 
