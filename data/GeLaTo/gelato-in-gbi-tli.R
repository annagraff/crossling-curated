# this script serves to validate the automated columns (columns 21 to 73) in the file GELATO-population-glottocode-mapping.csv
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

glottolog_level <- read.csv("../../raw/glottolog_v.4.8/languages_and_dialects_geo.csv")

# however note that in gelato we still have a location column...
gelato_mapping_decisions <- read.csv("GELATO-population-glottocode-mapping.csv")
gelato_mapping_decisions <- na_convert(gelato_mapping_decisions)

### check all densities are correct
gelato <- gelato_mapping_decisions %>% select(1:20)

# retrieve proxies
proxies <- gelato[,c(11:20)]
proxies <- na_convert(proxies)

base_status <- apply(as.data.frame(gelato),1,function(x)if(x[7]%in%glottolog_level$glottocode==F){NA}else{filter(glottolog_level,glottocode==x[7])[1,4]})
upstream_glottocode <- apply(as.data.frame(gelato),1,function(x)if(x[7]%in%taxonomy$glottocode==F){NA}else{filter(taxonomy,glottocode==x[7])[1,sum(!is.na(filter(taxonomy,glottocode==x[7])))-1]})
upstream_glottocode_status <- apply(as.data.frame(upstream_glottocode),1,function(x)if(x%in%glottolog_level$glottocode==F){NA}else{filter(glottolog_level,glottocode==x[1])[1,4]})

gelato <- cbind(gelato,data.frame(status.glottolog = base_status,
                                  glottocodeUpstream = upstream_glottocode,
                                  glottocodeUpstream.status = upstream_glottocode_status))

# read in GBI data: logical and statistical, full and densified
gbilogfull <- read.csv("../GBInd/output/logicalGBI/logicalGBI.csv", row.names = "glottocode") %>% select(-X)
gbilogdens <- read.csv("../GBInd/output/logicalGBI/logicalGBI_pruned.csv", row.names = "glottocode") %>% select(-X)

gbistafull <- read.csv("../GBInd/output/statisticalGBI/statisticalGBI.csv", row.names = "glottocode") %>% select(-X)
gbistadens <- read.csv("../GBInd/output/statisticalGBI/statisticalGBI_pruned.csv", row.names = "glottocode") %>% select(-X)

# read in TLI data: logical and statistical, full, large densified and small densified 
tlilogfull <- read.csv("../TypLinkInd/output/logicalTLI/full/logicalTLI_full.csv", row.names = "glottocode") %>% select(-X)
tlilogl <- read.csv("../TypLinkInd/output/logicalTLI/full/logicalTLI_full_pruned_large.csv", row.names = "glottocode") %>% select(-X)
tlilogs <- read.csv("../TypLinkInd/output/logicalTLI/full/logicalTLI_full_pruned_small.csv", row.names = "glottocode") %>% select(-X)

tlistafull <- read.csv("../TypLinkInd/output/statisticalTLI/full/statisticalTLI_full.csv", row.names = "glottocode") %>% select(-X)
tlistal <- read.csv("../TypLinkInd/output/statisticalTLI/full/statisticalTLI_full_pruned_large.csv", row.names = "glottocode") %>% select(-X)
tlistas <- read.csv("../TypLinkInd/output/statisticalTLI/full/statisticalTLI_full_pruned_small.csv", row.names = "glottocode") %>% select(-X)

# turn all ? or missing data to NA
gbilogfull <- na_convert(gbilogfull)
gbilogdens <- na_convert(gbilogdens)
gbistafull <- na_convert(gbistafull)
gbistadens <- na_convert(gbistadens)

tlilogfull <- na_convert(tlilogfull)
tlilogl <- na_convert(tlilogl)
tlilogs <- na_convert(tlilogs)
tlistafull <- na_convert(tlistafull)
tlistal <- na_convert(tlistal)
tlistas <- na_convert(tlistas)

# which languages are in any of the GBI and/or TLI curations?
glottocodesInAnyGBI <- unique(c(rownames(gbilogfull),rownames(gbilogdens),rownames(gbistafull),rownames(gbistadens)))
glottocodesInAnyTLI <- unique(c(rownames(tlilogfull),rownames(tlilogl),rownames(tlilogs),rownames(tlistafull),rownames(tlistal),rownames(tlistas)))

# variable counts for each dataset
gbilogfull.count <- apply(as.data.frame(glottocodesInAnyGBI),1,function(x)if(x%in%rownames(gbilogfull)==F){NA}else{length(na.omit(t(gbilogfull[which(rownames(gbilogfull)==x),])))})
gbilogdens.count <- apply(as.data.frame(glottocodesInAnyGBI),1,function(x)if(x%in%rownames(gbilogdens)==F){NA}else{length(na.omit(t(gbilogdens[which(rownames(gbilogdens)==x),])))})
gbistafull.count <- apply(as.data.frame(glottocodesInAnyGBI),1,function(x)if(x%in%rownames(gbistafull)==F){NA}else{length(na.omit(t(gbistafull[which(rownames(gbistafull)==x),])))})
gbistadens.count <- apply(as.data.frame(glottocodesInAnyGBI),1,function(x)if(x%in%rownames(gbistadens)==F){NA}else{length(na.omit(t(gbistadens[which(rownames(gbistadens)==x),])))})

tlilogfull.count <- apply(as.data.frame(glottocodesInAnyTLI),1,function(x)if(x%in%rownames(tlilogfull)==F){NA}else{length(na.omit(t(tlilogfull[which(rownames(tlilogfull)==x),])))})
tlilogl.count <- apply(as.data.frame(glottocodesInAnyTLI),1,function(x)if(x%in%rownames(tlilogl)==F){NA}else{length(na.omit(t(tlilogl[which(rownames(tlilogl)==x),])))})
tlilogs.count <- apply(as.data.frame(glottocodesInAnyTLI),1,function(x)if(x%in%rownames(tlilogs)==F){NA}else{length(na.omit(t(tlilogs[which(rownames(tlilogs)==x),])))})
tlistafull.count <- apply(as.data.frame(glottocodesInAnyTLI),1,function(x)if(x%in%rownames(tlistafull)==F){NA}else{length(na.omit(t(tlistafull[which(rownames(tlistafull)==x),])))})
tlistal.count <- apply(as.data.frame(glottocodesInAnyTLI),1,function(x)if(x%in%rownames(tlistal)==F){NA}else{length(na.omit(t(tlistal[which(rownames(tlistal)==x),])))})
tlistas.count <- apply(as.data.frame(glottocodesInAnyTLI),1,function(x)if(x%in%rownames(tlistas)==F){NA}else{length(na.omit(t(tlistas[which(rownames(tlistas)==x),])))})

# densities for all datasets
gbilogfull.density <- apply(as.data.frame(glottocodesInAnyGBI),1,function(x)if(x%in%rownames(gbilogfull)==F){NA}else{length(na.omit(t(gbilogfull[which(rownames(gbilogfull)==x),])))/ncol(gbilogfull)})
gbilogdens.density <- apply(as.data.frame(glottocodesInAnyGBI),1,function(x)if(x%in%rownames(gbilogdens)==F){NA}else{length(na.omit(t(gbilogdens[which(rownames(gbilogdens)==x),])))/ncol(gbilogdens)})
gbistafull.density <- apply(as.data.frame(glottocodesInAnyGBI),1,function(x)if(x%in%rownames(gbistafull)==F){NA}else{length(na.omit(t(gbistafull[which(rownames(gbistafull)==x),])))/ncol(gbistafull)})
gbistadens.density <- apply(as.data.frame(glottocodesInAnyGBI),1,function(x)if(x%in%rownames(gbistadens)==F){NA}else{length(na.omit(t(gbistadens[which(rownames(gbistadens)==x),])))/ncol(gbistadens)})

tlilogfull.density <- apply(as.data.frame(glottocodesInAnyTLI),1,function(x)if(x%in%rownames(tlilogfull)==F){NA}else{length(na.omit(t(tlilogfull[which(rownames(tlilogfull)==x),])))/ncol(tlilogfull)})
tlilogl.density <- apply(as.data.frame(glottocodesInAnyTLI),1,function(x)if(x%in%rownames(tlilogl)==F){NA}else{length(na.omit(t(tlilogl[which(rownames(tlilogl)==x),])))/ncol(tlilogl)})
tlilogs.density <- apply(as.data.frame(glottocodesInAnyTLI),1,function(x)if(x%in%rownames(tlilogs)==F){NA}else{length(na.omit(t(tlilogs[which(rownames(tlilogs)==x),])))/ncol(tlilogs)})
tlistafull.density <- apply(as.data.frame(glottocodesInAnyTLI),1,function(x)if(x%in%rownames(tlistafull)==F){NA}else{length(na.omit(t(tlistafull[which(rownames(tlistafull)==x),])))/ncol(tlistafull)})
tlistal.density <- apply(as.data.frame(glottocodesInAnyTLI),1,function(x)if(x%in%rownames(tlistal)==F){NA}else{length(na.omit(t(tlistal[which(rownames(tlistal)==x),])))/ncol(tlistal)})
tlistas.density <- apply(as.data.frame(glottocodesInAnyTLI),1,function(x)if(x%in%rownames(tlistas)==F){NA}else{length(na.omit(t(tlistas[which(rownames(tlistas)==x),])))/ncol(tlistas)})

# subset taxonomy to lgs that are in gelato or mm (each curation)
taxonomy_for_densities <- taxonomy
taxonomy_for_densities$GelatoBase <- apply(taxonomy,1,function(x)x[1]%in%gelato$glottocodeBase)
taxonomy_for_densities$GelatoUpstream <- apply(taxonomy,1,function(x)x[1]%in%gelato$glottocodeUpstream)
taxonomy_for_densities$InGBI <- apply(taxonomy,1,function(x)x[1]%in%glottocodesInAnyGBI)
taxonomy_for_densities$InTLI <- apply(taxonomy,1,function(x)x[1]%in%glottocodesInAnyTLI)
taxonomy_for_densities <- taxonomy_for_densities %>% filter(GelatoBase == T | GelatoUpstream == T | InGBI == T | InTLI == T)

lgs_in_any_GB <- data.frame(glottocode = glottocodesInAnyGBI,
                            features.count.gbi.logical.full = gbilogfull.count,
                            features.count.gbi.logical.densified = gbilogdens.count,
                            features.count.gbi.statistical.full = gbistafull.count,
                            features.count.gbi.statistical.densified = gbistadens.count, 
                            density.gbi.logical.full= gbilogfull.density,
                            density.gbi.logical.densified= gbilogdens.density,
                            density.gbi.statistical.full= gbistafull.density,
                            density.gbi.statistical.densified = gbistadens.density)

lgs_in_any_TLI <- data.frame(glottocode = glottocodesInAnyTLI,
                            features.count.tli.logical.full = tlilogfull.count,
                            features.count.tli.logical.large = tlilogl.count,
                            features.count.tli.logical.small = tlilogs.count,
                            features.count.tli.statistical.full = tlistafull.count,
                            features.count.tli.statistical.large= tlistal.count, 
                            features.count.tli.statistical.small = tlistas.count, 
                            density.tli.logical.full= tlilogfull.density,
                            density.tli.logical.large= tlilogl.density,
                            density.tli.logical.small = tlilogs.density,                             
                            density.tli.statistical.full= tlistafull.density,
                            density.tli.statistical.large = tlistal.density,
                            density.tli.statistical.small = tlistas.density)
                            
taxonomy_for_densities <- left_join(taxonomy_for_densities,lgs_in_any_GB)
taxonomy_for_densities <- left_join(taxonomy_for_densities,lgs_in_any_TLI)

## helper columns: counts and densities for all datasets for glottocodeBase
gelato <- left_join(gelato,select(filter(taxonomy_for_densities,GelatoBase == T),c(glottocode,
                                                                            features.count.gbi.logical.full,
                                                                            features.count.gbi.logical.densified,
                                                                            features.count.gbi.statistical.full,
                                                                            features.count.gbi.statistical.densified,
                                                                            features.count.tli.logical.full,
                                                                            features.count.tli.logical.large,
                                                                            features.count.tli.logical.small,
                                                                            features.count.tli.statistical.full,
                                                                            features.count.tli.statistical.large,
                                                                            features.count.tli.statistical.small,
                                                                            density.gbi.logical.full,
                                                                            density.gbi.logical.densified,
                                                                            density.gbi.statistical.full,
                                                                            density.gbi.statistical.densified,
                                                                            density.tli.logical.full,
                                                                            density.tli.logical.large,
                                                                            density.tli.logical.small,
                                                                            density.tli.statistical.full,
                                                                            density.tli.statistical.large,
                                                                            density.tli.statistical.small)), by=c("glottocodeBase"="glottocode"))


## helper columns: counts and densities for all datasets for glottocodeUpstream
gelato <- left_join(gelato,select(filter(taxonomy_for_densities,GelatoUpstream == T),c(glottocode,
                                                                                       features.count.gbi.logical.full,
                                                                                       features.count.gbi.logical.densified,
                                                                                       features.count.gbi.statistical.full,
                                                                                       features.count.gbi.statistical.densified,
                                                                                       features.count.tli.logical.full,
                                                                                       features.count.tli.logical.large,
                                                                                       features.count.tli.logical.small,
                                                                                       features.count.tli.statistical.full,
                                                                                       features.count.tli.statistical.large,
                                                                                       features.count.tli.statistical.small,
                                                                                       density.gbi.logical.full,
                                                                                       density.gbi.logical.densified,
                                                                                       density.gbi.statistical.full,
                                                                                       density.gbi.statistical.densified,
                                                                                       density.tli.logical.full,
                                                                                       density.tli.logical.large,
                                                                                       density.tli.logical.small,
                                                                                       density.tli.statistical.full,
                                                                                       density.tli.statistical.large,
                                                                                       density.tli.statistical.small)), by=c("glottocodeUpstream"="glottocode"))


names(gelato)[24:63]<-c("base.features.count.gbi.logical.full","base.features.count.gbi.logical.densified","base.features.count.gbi.statistical.full","base.features.count.gbi.statistical.densified","base.features.count.tli.logical.full","base.features.count.tli.logical.large","base.features.count.tli.logical.small","base.features.count.tli.statistical.full","base.features.count.tli.statistical.large","base.features.count.tli.statistical.small",
                       "base.density.gbi.logical.full","base.density.gbi.logical.densified","base.density.gbi.statistical.full","base.density.gbi.statistical.densified","base.density.tli.logical.full","base.density.tli.logical.large","base.density.tli.logical.small","base.density.tli.statistical.full","base.density.tli.statistical.large","base.density.tli.statistical.small",
                       "upstream.features.count.gbi.logical.full","upstream.features.count.gbi.logical.densified","upstream.features.count.gbi.statistical.full","upstream.features.count.gbi.statistical.densified","upstream.features.count.tli.logical.full","upstream.features.count.tli.logical.large","upstream.features.count.tli.logical.small","upstream.features.count.tli.statistical.full","upstream.features.count.tli.statistical.large","upstream.features.count.tli.statistical.small",
                       "upstream.density.gbi.logical.full","upstream.density.gbi.logical.densified","upstream.density.gbi.statistical.full","upstream.density.gbi.statistical.densified","upstream.density.tli.logical.full","upstream.density.tli.logical.large","upstream.density.tli.logical.small","upstream.density.tli.statistical.full","upstream.density.tli.statistical.large","upstream.density.tli.statistical.small")


proxy.gbi.logical.full.density <- unlist(apply(as.data.frame(proxies[,1]),1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density.gbi.logical.full"]}))
proxy.gbi.logical.densified.density <- unlist(apply(as.data.frame(proxies[,2]),1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density.gbi.logical.densified"]}))
proxy.gbi.statistical.full.density <- unlist(apply(as.data.frame(proxies[,3]),1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density.gbi.statistical.full"]}))
proxy.gbi.statistical.densified.density <- unlist(apply(as.data.frame(proxies[,4]),1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density.gbi.statistical.densified"]}))
proxy.tli.logical.full.density <- unlist(apply(as.data.frame(proxies[,5]),1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density.tli.logical.full"]}))
proxy.tli.logical.large.density <- unlist(apply(as.data.frame(proxies[,6]),1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density.tli.logical.large"]}))
proxy.tli.logical.small.density <- unlist(apply(as.data.frame(proxies[,7]),1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density.tli.logical.small"]}))
proxy.tli.statistical.full.density <- unlist(apply(as.data.frame(proxies[,8]),1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density.tli.statistical.full"]}))
proxy.tli.statistical.large.density <- unlist(apply(as.data.frame(proxies[,9]),1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density.tli.statistical.large"]}))
proxy.tli.statistical.small.density <- unlist(apply(as.data.frame(proxies[,10]),1,function(x)if(is.na(x)){NA}else{taxonomy_for_densities[which(taxonomy_for_densities$glottocode==x),"density.tli.statistical.small"]}))

# add columns specifying coverage for all proxies
gelato <- cbind(gelato,proxy.gbi.logical.full.density,proxy.gbi.logical.densified.density,proxy.gbi.statistical.full.density,proxy.gbi.statistical.densified.density,
                proxy.tli.logical.full.density,proxy.tli.logical.large.density,proxy.tli.logical.small.density,
                proxy.tli.statistical.full.density,proxy.tli.statistical.large.density,proxy.tli.statistical.small.density)


gelato_automated_glottologStatus <- gelato[,21:23]
gelato_automated_numeric <- gelato[,24:ncol(gelato)]

# now check everything is correct:
expect_true(all(gelato_automated_glottologStatus==gelato_mapping_decisions[21:23], na.rm = T))
expect_true(all(round(gelato_automated_numeric,digits=5)==round(gelato_mapping_decisions[24:ncol(gelato_mapping_decisions)],digits=5), na.rm = T))
