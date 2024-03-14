# code for Figs 10, 11 and 12
rm(list=ls())

library(densify)
library(tidyverse)
library(rnaturalearth)
library(cowplot)
library(gridExtra)
library(pcaMethods)

source("../data/functions.R")

## Retrieve taxonomy, coordinates and map

# generate taxonomy
taxonomy_matrix_densify <- as_flat_taxonomy_matrix(glottolog_languoids)

# link coordinates to lgs
macroareas <- read.csv("../data/lang-metadata.csv")
names(macroareas)[1] <- "id"
macroareas <- na_convert(macroareas)
taxonomy_matrix <- merge(taxonomy_matrix_densify,select(macroareas, c("id","latitude","longitude","glottolog.macroarea")),by="id")
taxonomy_matrix_for_plot <- merge(taxonomy_matrix_densify,select(macroareas, c("id","latitude","longitude","glottolog.macroarea")),by="id") %>% filter(!is.na(latitude))
names(taxonomy_matrix)[29:30] <- c("lat","lon") # change names to lat and lon for downstream functions
names(taxonomy_matrix_for_plot)[29:30] <- c("lat","lon") # change names to lat and lon for downstream functions

# download Natural Earth data for sf, including islands
world_map <- ne_download(scale = 110, category = "physical", type = "land", returnclass = "sf")
minor_islands <- ne_download(scale = 10, category = "physical", type = "minor_islands", returnclass = "sf")
world_map_initial <- rbind(world_map, minor_islands)

## GBI logical data, full vs. densified
# read in matrices, convert all NA and ? to NA, factorize all variables
gbi_logical_full <- read.csv("../data/GBInd/output/logicalGBI/logicalGBI.csv",row.names = "glottocode") %>% select(-X)
gbi_logical_full <- na_convert(gbi_logical_full)
gbi_logical_full <- factorise(gbi_logical_full)

gbi_logical_densified <- read.csv("../data/GBInd/output/logicalGBI/logicalGBI_pruned.csv",row.names = "glottocode") %>% select(-X)
gbi_logical_densified <- na_convert(gbi_logical_densified)
gbi_logical_densified <- factorise(gbi_logical_densified)

## TLI logical data: large vs. small densified for full and morphosyntax data, densified for phonology and lexicon
# full set of variables
tli_logical_large <- read.csv("../data/TypLinkInd/output/logicalTLI/full/logicalTLI_full_pruned_large.csv",row.names = "glottocode") %>% select(-X)
tli_logical_large <- na_convert(tli_logical_large)
tli_logical_large <- factorise(tli_logical_large)

tli_logical_small <- read.csv("../data/TypLinkInd/output/logicalTLI/full/logicalTLI_full_pruned_small.csv",row.names = "glottocode") %>% select(-X)
tli_logical_small <- na_convert(tli_logical_small)
tli_logical_small <- factorise(tli_logical_small)

# morphosyntax
morphosyntax_logical_large <- read.csv("../data/TypLinkInd/output/logicalTLI/morphosyntax/logicalTLI_morphosyntax_pruned_large.csv",row.names = "glottocode") %>% select(-X)
morphosyntax_logical_large <- na_convert(morphosyntax_logical_large)
morphosyntax_logical_large <- factorise(morphosyntax_logical_large)

# phonology
phonology_logical_pruned <- read.csv("../data/TypLinkInd/output/logicalTLI/phonology/logicalTLI_phonology_pruned.csv",row.names = "glottocode") %>% select(-X)
phonology_logical_pruned <- na_convert(phonology_logical_pruned)
phonology_logical_pruned <- factorise(phonology_logical_pruned)

# lexicon
lexicon_logical_pruned <- read.csv("../data/TypLinkInd/output/logicalTLI/lexicon/logicalTLI_lexicon_pruned.csv",row.names = "glottocode") %>% select(-X)
lexicon_logical_pruned <- na_convert(lexicon_logical_pruned)
lexicon_logical_pruned <- factorise(lexicon_logical_pruned)

## define functions
# create per family proportion table for pca
generate_per_family_prop_table <- function(data,taxonomy_matrix){
  # For multinomial variables we compute the proportion presence of each level (feature, value) except whichever is the largest level (since one proportion is always redundant), all variables are `factor`.
  lvls <- apply(data, 2, function(x)list(as.factor(names(sort(table(x),decreasing=T))[-1])))
  data$family <- as.factor(unlist(apply(as.data.frame(rownames(data)),1,function(x)filter(taxonomy_matrix,id==x)$level1))) # add family
  data <- select(data,c(ncol(data),1:ncol(data)-1)) # reorder for family to be first column
  props_table <- data.frame(family=sort(unique(data$family)))
  for (vbl in 2:ncol(data)){ # variable-wise
    tbl <- table(data[,"family"], data[,vbl])
    levs <- unname(unlist(lvls[vbl-1]))
    for (lvl in 1:length(levs)){ # level-wise
      props_table<-cbind(props_table,tbl[,levs[lvl]]/rowSums(tbl[,]))
      names(props_table)[ncol(props_table)]<-paste(names(data)[vbl],levs[lvl],sep="_")
    }
  }
  ## remove "families", which actually aren't genealogical units: Artificial Language, Bookkeeping, Mixed Language, Pidgin, Sign Language, Speech Register, Unattested, Unclassifiable
  props_table <- props_table %>% filter(family %in% c("arti1236","book1242","mixe1287","pidg1258","sign1238","spee1234","unat1236","uncl1493") == FALSE)
  
  rownames(props_table)<-props_table$family
  props_table <- select(props_table,-family)
  return(props_table)
}

# explained variances from PCA
explained_variance <- function(pca_object) {
  variances <- pca_object@sDev^2
  r2cum <- pca_object@R2cum
  
  explvar <- rbind(cumsum(variances/sum(variances)),r2cum)
  rownames(explvar) <- c("explained variance", "R2cum")
  return(explvar)
}

# report loadings from PCA 
# For reporting the heaviest contributors per PC, get full names
report_loadings <- function(prop_pca, directory, pc, nr_variables_from_top_and_bottom) {
  prop_pca_loadings <- loadings(prop_pca)
  descriptions <- read.csv(paste("output/", directory,"/cldf/parameters.csv", collapse = "", sep = "")) %>% select(c("new.name","description"))
  
  # but focusing only on major differences:
  loadings <- data.frame(
    PC1.max.name.and.value=names(sort(prop_pca_loadings[,1], decr=T)[1:nr_variables_from_top_and_bottom]),
    PC1.max.variance.explained=(sort(prop_pca_loadings[,1], decr=T)[1:nr_variables_from_top_and_bottom])^2, 
    PC1.max.description=apply(array(names(sort(prop_pca_loadings[,1], decr=T)[1:nr_variables_from_top_and_bottom])),1,function(x)filter(descriptions,new.name%in%sub("_[^_]+$", "", gsub("\\.", "+", x)))$description),
    PC1.min.name.and.value=names(sort(prop_pca_loadings[,1], decr=F)[1:nr_variables_from_top_and_bottom]),
    PC1.min.variance.explained=(sort(prop_pca_loadings[,1], decr=F)[1:nr_variables_from_top_and_bottom])^2,
    PC1.min.description=apply(array(names(sort(prop_pca_loadings[,1], decr=F)[1:nr_variables_from_top_and_bottom])),1,function(x)filter(descriptions,new.name%in%sub("_[^_]+$", "", gsub("\\.", "+", x)))$description),
    
    PC2.max.name.and.value=names(sort(prop_pca_loadings[,2], decr=T)[1:nr_variables_from_top_and_bottom]),
    PC2.max.variance.explained=(sort(prop_pca_loadings[,2], decr=T)[1:nr_variables_from_top_and_bottom])^2, 
    PC2.max.description=apply(array(names(sort(prop_pca_loadings[,2], decr=T)[1:nr_variables_from_top_and_bottom])),1,function(x)filter(descriptions,new.name%in%sub("_[^_]+$", "", gsub("\\.", "+", x)))$description),
    PC2.min.name.and.value=names(sort(prop_pca_loadings[,2], decr=F)[1:nr_variables_from_top_and_bottom]),
    PC2.min.variance.explained=(sort(prop_pca_loadings[,2], decr=F)[1:nr_variables_from_top_and_bottom])^2,
    PC2.min.description=apply(array(names(sort(prop_pca_loadings[,2], decr=F)[1:nr_variables_from_top_and_bottom])),1,function(x)filter(descriptions,new.name%in%sub("_[^_]+$", "", gsub("\\.", "+", x)))$description),
    
    PC3.max.name.and.value=names(sort(prop_pca_loadings[,3], decr=T)[1:nr_variables_from_top_and_bottom]),
    PC3.max.variance.explained=(sort(prop_pca_loadings[,3], decr=T)[1:nr_variables_from_top_and_bottom])^2,     
    PC3.max.description=apply(array(names(sort(prop_pca_loadings[,3], decr=T)[1:nr_variables_from_top_and_bottom])),1,function(x)filter(descriptions,new.name%in%sub("_[^_]+$", "", gsub("\\.", "+", x)))$description),
    PC3.min.name.and.value=names(sort(prop_pca_loadings[,3], decr=F)[1:nr_variables_from_top_and_bottom]),
    PC3.min.variance.explained=(sort(prop_pca_loadings[,3], decr=F)[1:nr_variables_from_top_and_bottom])^2,
    PC3.min.description=apply(array(names(sort(prop_pca_loadings[,3], decr=F)[1:nr_variables_from_top_and_bottom])),1,function(x)filter(descriptions,new.name%in%sub("_[^_]+$", "", gsub("\\.", "+", x)))$description),
    
    PC4.max.name.and.value=names(sort(prop_pca_loadings[,4], decr=T)[1:nr_variables_from_top_and_bottom]),
    PC4.max.variance.explained=(sort(prop_pca_loadings[,4], decr=T)[1:nr_variables_from_top_and_bottom])^2, 
    PC4.max.description=apply(array(names(sort(prop_pca_loadings[,4], decr=T)[1:nr_variables_from_top_and_bottom])),1,function(x)filter(descriptions,new.name%in%sub("_[^_]+$", "", gsub("\\.", "+", x)))$description),
    PC4.min.name.and.value=names(sort(prop_pca_loadings[,4], decr=F)[1:nr_variables_from_top_and_bottom]),
    PC4.min.variance.explained=(sort(prop_pca_loadings[,4], decr=F)[1:nr_variables_from_top_and_bottom])^2,
    PC4.min.description=apply(array(names(sort(prop_pca_loadings[,4], decr=F)[1:nr_variables_from_top_and_bottom])),1,function(x)filter(descriptions,new.name%in%sub("_[^_]+$", "", gsub("\\.", "+", x)))$description),
    
    PC5.max.name.and.value=names(sort(prop_pca_loadings[,5], decr=T)[1:nr_variables_from_top_and_bottom]),
    PC5.max.variance.explained=(sort(prop_pca_loadings[,5], decr=T)[1:nr_variables_from_top_and_bottom])^2, 
    PC5.max.description=apply(array(names(sort(prop_pca_loadings[,5], decr=T)[1:nr_variables_from_top_and_bottom])),1,function(x)filter(descriptions,new.name%in%sub("_[^_]+$", "", gsub("\\.", "+", x)))$description),
    PC5.min.name.and.value=names(sort(prop_pca_loadings[,5], decr=F)[1:nr_variables_from_top_and_bottom]),
    PC5.min.variance.explained=(sort(prop_pca_loadings[,5], decr=F)[1:nr_variables_from_top_and_bottom])^2,
    PC5.min.description=apply(array(names(sort(prop_pca_loadings[,5], decr=F)[1:nr_variables_from_top_and_bottom])),1,function(x)filter(descriptions,new.name%in%sub("_[^_]+$", "", gsub("\\.", "+", x)))$description),
    
    PC6.max.name.and.value=names(sort(prop_pca_loadings[,6], decr=T)[1:nr_variables_from_top_and_bottom]),
    PC6.max.variance.explained=(sort(prop_pca_loadings[,6], decr=T)[1:nr_variables_from_top_and_bottom])^2,
    PC6.max.description=apply(array(names(sort(prop_pca_loadings[,6], decr=T)[1:nr_variables_from_top_and_bottom])),1,function(x)filter(descriptions,new.name%in%sub("_[^_]+$", "", gsub("\\.", "+", x)))$description),
    PC6.min.name.and.value=names(sort(prop_pca_loadings[,6], decr=F)[1:nr_variables_from_top_and_bottom]),
    PC6.min.variance.explained=(sort(prop_pca_loadings[,6], decr=F)[1:nr_variables_from_top_and_bottom])^2,
    PC6.min.description=apply(array(names(sort(prop_pca_loadings[,6], decr=F)[1:nr_variables_from_top_and_bottom])),1,function(x)filter(descriptions,new.name%in%sub("_[^_]+$", "", gsub("\\.", "+", x)))$description)
  )
  rownames(loadings) <- NULL
  
  loadings[((pc-1)*6+1):((pc-1)*6+6)]
}

# groups of three PCs (PC1-PC3, PC4-PC6, PC7-PC9) are mapped to RGB color space and the results are linked back to geographical and genealogical information:
# we allow for manually introducing pc-value-flipping to embellish colours or to make corresponding plots align visually
RGB_mapping <- function(pca.object, taxonomy_matrix, pc1flip=F, pc2flip=F, pc3flip=F) {
  library(scales)
  gg <- as.data.frame(taxonomy_matrix)[ , c('level1', 'id' ,'lat', 'lon')]
  df <- data.frame(scores(pca.object), level1=rownames(scores(pca.object)))
 
  # map PCs to RGB:
  if(pc1flip){
    pcrange <- range(df$PC1)
    pcmid <- mean(pcrange)
    pcflip <- df$PC1-2*(df$PC1-pcmid)
    df$PC1 <- pcflip
  }
  
  if(pc2flip){
    pcrange <- range(df$PC2)
    pcmid <- mean(pcrange)
    pcflip <- df$PC2-2*(df$PC2-pcmid)
    df$PC2 <- pcflip
  }
  
  if(pc3flip){
    pcrange <- range(df$PC3)
    pcmid <- mean(pcrange)
    pcflip <- df$PC3-2*(df$PC3-pcmid)
    df$PC3 <- pcflip
  }

  df$pc.colors1to3 <- with(df, rgb(red=rescale_mid(PC1), 
                                   green=rescale_mid(PC2), 
                                   blue=rescale_mid(PC3), alpha=0.8))
  df$pc.colors4to6 <- with(df, rgb(red=rescale_mid(PC4), 
                                   green=rescale_mid(PC5), 
                                   blue=rescale_mid(PC6), alpha=0.8))
  df$pc.colors7to9 <- with(df, rgb(red=rescale_mid(PC7), 
                                   green=rescale_mid(PC8), 
                                   blue=rescale_mid(PC9), alpha=0.8))
  geo <- subset(gg, gg[,1] %in% df$level1)
  geo.world <- merge(df[,c('level1','pc.colors1to3', 'pc.colors4to6','pc.colors7to9')], 
                     geo, by.x='level1', by.y=names(gg)[1])
  return(geo.world)
}

# groups of PCs are mapped to 2-D space and visualised according to macroarea
dim_mapping <- function(pca.object, taxonomy_matrix, size) {
  library(scales)
  library(lattice)
  
  gg <- as.data.frame(taxonomy_matrix)[ , c('level1','glottolog.macroarea')] %>% filter(level1 %in% rownames(pca.object@scores))
  gg <- data.frame(Macroarea = unlist(apply(table(gg),1,function(x)colnames(table(gg))[which(x==max(x))])))
  gg$family <- rownames(gg)
  
  df <- data.frame(scores(pca.object), family=rownames(scores(pca.object)))
  df <- left_join(df, gg)
  
  # scale PCs:
  df$pc1 <- with(df, rescale_mid(PC1))
  df$pc2 <- with(df, rescale_mid(PC2))
  
  # retrieve percentages
  prop1 <- paste("PC1 (",round(pca.object@R2cum[1],4)*100," %)",sep="")
  prop2 <- paste("PC2 (",round(pca.object@R2cum[2]-pca.object@R2cum[1],4)*100," %)",sep="")
  
  p <- ggplot(df, aes(x=pc1,y=pc2,group=Macroarea))+
    geom_point(size=size, aes(shape=Macroarea))+
    scale_shape_manual(values=c(1,15,2,3,8,4))+
    theme_bw()+
    xlab(prop1)+
    ylab(prop2)+
    labs(subtitle = "", shape = "")

  # # if 3d-plot desired:
  # library(lattice)
  # df$pc3 <- with(df, rescale_mid(PC3))
  # cube <- cloud(pc1~pc2*pc3, data=df,
  #       pch=1,
  #       cex=.4,
  #       groups = macroarea,
  #       auto.key = TRUE,
  #       par.box=list(lty=3),
  #       main = "Families in PC-space, by macroregion")
  
  return(p)
}

# plotting RGB map 
rgb_map <- function(rgb_data_frame, world_map_initial, main, plot, size, pca.object){
  main <- paste(main,", PC1-PC3 (",round(pca.object@R2cum[3],4)*100,"%)",sep="")
  rgb_shifted <- project_data(df = rgb_data_frame, world_map_initial = world_map_initial)
  pacific_centered_map <- rgb_shifted$base_plot +
    geom_sf(data = rgb_shifted$data, aes(fill = pc.colors1to3, color = pc.colors1to3), shape = 21, size = size) +
    labs(title = main) +
    theme_minimal() +
    scale_fill_identity() +  # use color values as specified in the data
    scale_color_identity() +  # use color values as specified in the data
    guides(fill = "none", color = "none")+
    theme(panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5))
  if(plot==T){
    print(pacific_centered_map)
  }
  return(pacific_centered_map)
}

## PCA, aggregated to family level
# generate per family proportion table
proportions_gbi_logical_full <- generate_per_family_prop_table(data = gbi_logical_full, taxonomy_matrix = taxonomy_matrix)
proportions_gbi_logical_densified <- generate_per_family_prop_table(data = gbi_logical_densified, taxonomy_matrix = taxonomy_matrix)
proportions_tli_logical_large <- generate_per_family_prop_table(data = tli_logical_large, taxonomy_matrix = taxonomy_matrix)
proportions_tli_logical_small <- generate_per_family_prop_table(data = tli_logical_small, taxonomy_matrix = taxonomy_matrix)
proportions_morphosyntax_logical_large <- generate_per_family_prop_table(data = morphosyntax_logical_large, taxonomy_matrix = taxonomy_matrix)
proportions_phonology_logical_pruned <- generate_per_family_prop_table(data = phonology_logical_pruned, taxonomy_matrix = taxonomy_matrix)
proportions_lexicon_logical_pruned <- generate_per_family_prop_table(data = lexicon_logical_pruned, taxonomy_matrix = taxonomy_matrix)

# # proportions of NA per family
# boxplot(apply(proportions_tli_logical_large,1,function(x)sum(is.na(x))/ncol(proportions_tli_logical_large)))
# title("Proportion of NA per family, large densified logical dataset")

# perform PCA, first replace "NaN" with NA
proportions_gbi_logical_full[proportions_gbi_logical_full=="NaN"]<-NA
proportions_gbi_logical_densified[proportions_gbi_logical_densified=="NaN"]<-NA
proportions_tli_logical_large[proportions_tli_logical_large=="NaN"]<-NA
proportions_tli_logical_small[proportions_tli_logical_small=="NaN"]<-NA
proportions_morphosyntax_logical_large[proportions_morphosyntax_logical_large=="NaN"]<-NA
proportions_phonology_logical_pruned[proportions_phonology_logical_pruned=="NaN"]<-NA
proportions_lexicon_logical_pruned[proportions_lexicon_logical_pruned=="NaN"]<-NA

per_family_props_PCA_gbi_logical_full <- pca(proportions_gbi_logical_full, method='rnipals', nPcs=10)
per_family_props_PCA_gbi_logical_densified <- pca(proportions_gbi_logical_densified, method='rnipals', nPcs=10)
per_family_props_PCA_tli_logical_large <- pca(proportions_tli_logical_large, method='rnipals', nPcs=10)
per_family_props_PCA_tli_logical_small <- pca(proportions_tli_logical_small, method='rnipals', nPcs=10)
per_family_props_PCA_morphosyntax_logical_large <- pca(proportions_morphosyntax_logical_large, method='rnipals', nPcs=10)
per_family_props_PCA_phonology_logical_pruned <- pca(proportions_phonology_logical_pruned, method='rnipals', nPcs=10)
per_family_props_PCA_lexicon_logical_pruned <- pca(proportions_lexicon_logical_pruned, method='rnipals', nPcs=10)

# explained variances and R2 among all PCs
explained_variance(per_family_props_PCA_gbi_logical_full)[,1:10]
explained_variance(per_family_props_PCA_gbi_logical_densified)[,1:10]
explained_variance(per_family_props_PCA_tli_logical_large)[,1:10]
explained_variance(per_family_props_PCA_tli_logical_small)[,1:10]
explained_variance(per_family_props_PCA_morphosyntax_logical_large)[,1:10]
explained_variance(per_family_props_PCA_phonology_logical_pruned)[,1:10]
explained_variance(per_family_props_PCA_lexicon_logical_pruned)[,1:10]

# most relevant loadings
# report_loadings(per_family_props_PCA_full_logical_large, directory = "logicalMM", pc = 1, nr_variables_from_top_and_bottom = 5)
# report_loadings(per_family_props_PCA_full_logical_large, directory = "logicalMM", pc = 2, nr_variables_from_top_and_bottom = 5)
# report_loadings(per_family_props_PCA_full_logical_large, directory = "logicalMM", pc = 3, nr_variables_from_top_and_bottom = 5)

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

# trellis.device(png, filename = "rgb_cube.png", type="cairo", units="in", width=5, height=5,  pointsize=10, res=300)
print(rgb_cube)
# dev.off()

# map languages to colours
rgb_mapping_gbi_logical_full_kalyan <- RGB_mapping(per_family_props_PCA_gbi_logical_full,taxonomy_matrix_for_plot, pc1flip = T, pc2flip = F, pc3flip = F) %>% filter(id%in%rownames(gbi_logical_full)) %>% na.omit()
rgb_mapping_gbi_logical_densified_kalyan <- RGB_mapping(per_family_props_PCA_gbi_logical_densified,taxonomy_matrix_for_plot, pc1flip = T, pc2flip = F, pc3flip = F) %>% filter(id%in%rownames(gbi_logical_densified)) %>% na.omit()
rgb_mapping_gbi_logical_full_panel10 <- RGB_mapping(per_family_props_PCA_gbi_logical_full,taxonomy_matrix_for_plot, pc1flip = F, pc2flip = F, pc3flip = T) %>% filter(id%in%rownames(gbi_logical_full)) %>% na.omit()
rgb_mapping_gbi_logical_densified_panel10 <- RGB_mapping(per_family_props_PCA_gbi_logical_densified,taxonomy_matrix_for_plot, pc1flip = F, pc2flip = F, pc3flip = T) %>% filter(id%in%rownames(gbi_logical_densified)) %>% na.omit()
rgb_mapping_tli_logical_large <- RGB_mapping(per_family_props_PCA_tli_logical_large,taxonomy_matrix_for_plot, pc1flip = F, pc2flip = T, pc3flip = T) %>% filter(id%in%rownames(tli_logical_large)) %>% na.omit()
rgb_mapping_tli_logical_small <- RGB_mapping(per_family_props_PCA_tli_logical_small,taxonomy_matrix_for_plot, pc1flip = F, pc2flip = T, pc3flip = T) %>% filter(id%in%rownames(tli_logical_small)) %>% na.omit()
rgb_mapping_morphosyntax_logical_large <- RGB_mapping(per_family_props_PCA_morphosyntax_logical_large,taxonomy_matrix_for_plot, pc1flip = F, pc2flip = T, pc3flip = T) %>% filter(id%in%rownames(morphosyntax_logical_large)) %>% na.omit()
rgb_mapping_phonology_logical_pruned <- RGB_mapping(per_family_props_PCA_phonology_logical_pruned,taxonomy_matrix_for_plot, pc1flip = T, pc2flip = T, pc3flip = T) %>% filter(id%in%rownames(phonology_logical_pruned)) %>% na.omit()
rgb_mapping_lexicon_logical_pruned <- RGB_mapping(per_family_props_PCA_lexicon_logical_pruned,taxonomy_matrix_for_plot, pc1flip = T, pc2flip = F, pc3flip = F) %>% filter(id%in%rownames(lexicon_logical_pruned)) %>% na.omit()

# global PC maps for PCs 1-3
rgb_map_gbi_logical_full_kalyan <- rgb_map(rgb_mapping_gbi_logical_full_kalyan, world_map_initial, main="GBInd (PC1 flipped)", plot=F, size = 1.2, pca.object = per_family_props_PCA_gbi_logical_full)
rgb_map_gbi_logical_densified_kalyan <- rgb_map(rgb_mapping_gbi_logical_densified_kalyan, world_map_initial, main="GBInd (densified)", plot=F, size = 1.2, pca.object = per_family_props_PCA_gbi_logical_densified)
rgb_map_gbi_logical_full_panel10 <- rgb_map(rgb_mapping_gbi_logical_full_panel10, world_map_initial, main="GBInd", plot=F, size = 1.2, pca.object = per_family_props_PCA_gbi_logical_full)
rgb_map_gbi_logical_full_panel12 <- rgb_map(rgb_mapping_gbi_logical_full_panel10, world_map_initial, main="GBInd (PC3 flipped)", plot=F, size = 1.2, pca.object = per_family_props_PCA_gbi_logical_full)
rgb_map_gbi_logical_densified_panel10 <- rgb_map(rgb_mapping_gbi_logical_densified_panel10, world_map_initial, main="GBInd (densified)", plot=F, size = 1.2, pca.object = per_family_props_PCA_gbi_logical_full)
rgb_map_tli_logical_large <- rgb_map(rgb_mapping_tli_logical_large, world_map_initial, main="TypLinkInd (large densified)", plot=F, size = 1.2, pca.object = per_family_props_PCA_tli_logical_large)
rgb_map_tli_logical_small <- rgb_map(rgb_mapping_tli_logical_small, world_map_initial, main="TypLinkInd (small densified)", plot=F, size = 1.2, pca.object = per_family_props_PCA_tli_logical_small)
rgb_map_morphosyntax_logical_large <- rgb_map(rgb_mapping_morphosyntax_logical_large, world_map_initial, main="TypLinkInd [morphosyntax]", plot=F, size = 1.2, pca.object = per_family_props_PCA_morphosyntax_logical_large)
rgb_map_phonology_logical_pruned <- rgb_map(rgb_mapping_phonology_logical_pruned, world_map_initial, main="TypLinkInd [phonology]", plot=F, size = 1.2, pca.object = per_family_props_PCA_phonology_logical_pruned)
rgb_map_lexicon_logical_pruned <- rgb_map(rgb_mapping_lexicon_logical_pruned, world_map_initial, main="TypLinkInd [lexicon]", plot=F, size = 1.2, per_family_props_PCA_lexicon_logical_pruned) 


#### family relationships to one another in PC-space, by macroarea
f2d_full_gbi_logical_full <- dim_mapping(pca.object = per_family_props_PCA_gbi_logical_full, taxonomy_matrix = taxonomy_matrix_for_plot, size = 1.2)
f2d_full_gbi_logical_densified <- dim_mapping(pca.object = per_family_props_PCA_gbi_logical_densified, taxonomy_matrix = taxonomy_matrix_for_plot, size = 1.2)
f2d_full_tli_logical_large <- dim_mapping(pca.object = per_family_props_PCA_tli_logical_large, taxonomy_matrix = taxonomy_matrix_for_plot, size = 1.2)
f2d_full_tli_logical_small <- dim_mapping(pca.object = per_family_props_PCA_tli_logical_small, taxonomy_matrix = taxonomy_matrix_for_plot, size = 1.2)
f2d_morphosyntax_logical_large <- dim_mapping(pca.object = per_family_props_PCA_morphosyntax_logical_large, taxonomy_matrix = taxonomy_matrix_for_plot, size = 1.2)
f2d_phonology_logical_pruned <- dim_mapping(pca.object = per_family_props_PCA_phonology_logical_pruned, taxonomy_matrix = taxonomy_matrix_for_plot, size = 1.2)
f2d_lexicon_logical_pruned <- dim_mapping(pca.object = per_family_props_PCA_lexicon_logical_pruned, taxonomy_matrix = taxonomy_matrix_for_plot, size = 1.2)

# arrange the plots side by side, print and save
logical_plots <- grid.arrange(
  arrangeGrob(rgb_map_gbi_logical_full_panel10, f2d_full_gbi_logical_full, nrow = 1, ncol = 2, widths=c(3,2)), 
  arrangeGrob(rgb_map_morphosyntax_logical_large, f2d_morphosyntax_logical_large, nrow = 1, ncol = 2, widths=c(3,2)), 
  arrangeGrob(rgb_map_phonology_logical_pruned, f2d_phonology_logical_pruned, nrow = 1, ncol = 2, widths=c(3,2)), 
  arrangeGrob(rgb_map_lexicon_logical_pruned, f2d_lexicon_logical_pruned, nrow = 1, ncol = 2, widths=c(3,2)), 
  nrow = 4
)
ggsave("Fig10.png", plot = logical_plots, width = 14, height = 8*2, dpi = 500)

logical_plots_densify <- grid.arrange(
  arrangeGrob(rgb_map_gbi_logical_full_panel10, f2d_full_gbi_logical_full, rgb_map_gbi_logical_densified_panel10, f2d_full_gbi_logical_densified, nrow = 2, ncol = 2, widths=c(3,2)),
  arrangeGrob(rgb_map_tli_logical_large, f2d_full_tli_logical_large, rgb_map_tli_logical_small, f2d_full_tli_logical_small, nrow = 2, ncol = 2, widths=c(3,2)),
  nrow = 2
)
ggsave("Fig11.png", plot = logical_plots_densify, width = 14, height = 8*2, dpi = 500)

comp_plots <- grid.arrange(
  arrangeGrob(rgb_map_gbi_logical_full_panel12, rgb_map_gbi_logical_full_kalyan, nrow = 2, ncol = 1)
)
ggsave("Fig12.png", plot = comp_plots, width = 9, height = 8, dpi = 500)


