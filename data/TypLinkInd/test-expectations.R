#### This script generates 1000 diversity samples of 120 languages from the languages in the TypLinkInd logical dataset.
#### It reproduces the expectation-testing performed in the statistical curation step and validates the meta-data on all statistical modification IDs in the file "decisions-log.csv")

rm(list=ls())

# load relevant libraries
library(tidyverse)
library(testthat)
library(readr)
library(densify)

# load functions
source("../functions.R")

### prepare data ###
# load all formulated expectations and stored results
expectations <- read.csv("input/decisions-log.csv") %>% filter(modification.type == "statistical")

# load logical grambank data, glottolog taxonomy and glottolog macroarea data
logical_data <- read.csv("output/logicalTLI/full/logicalTLI_full.csv") %>% select(-X)
taxonomy <- as_flat_taxonomy_matrix(glottolog_languoids)
macroareas <- read.csv("../../raw/glottolog_v.4.8/languages_and_dialects_geo.csv")
names(macroareas)[1]<-"id"
taxonomy <- merge(taxonomy,select(macroareas, c("id","macroarea")),by="id")

# to create the samples, take into consideration all glottocodes for which at least 10 features are coded (not ? or NA) in the design.logical matrix
data_for_sample_generation <- logical_data
data_for_sample_generation[data_for_sample_generation=="?"]<-NA
data_for_sample_generation[data_for_sample_generation=="NA"]<-NA
data_for_sample_generation$sum.non.na <- apply(data_for_sample_generation,1,function(x)length(na.omit(x))-1)
data_for_sample_generation <- data_for_sample_generation %>% filter(sum.non.na > 9)

# select these glottocodes; subset taxonomy to relevant glottocodes
lgs <- as.character(data_for_sample_generation$glottocode)
lgs <- filter(taxonomy,id%in%lgs)

# load data for stats
recoded_data <- read.csv("output/statisticalTLI/data_for_stats.csv") %>% select(-X)

# replace . by + in variable names
names(recoded_data) <- str_replace_all(names(recoded_data), c("\\."="\\+"))

### create 1000 diversity samples; each row of diversity_samples is a diversity sample of 120 lgs; 1 from each family (glottolog-based) 
# list all six glottolog macroareas
regions <- data.frame(regions=c("Africa","Eurasia","Papunesia","Australia","South America","North America"))

# create empty matrix of 120 rows (languages) for 1000 samples (columns)
lg_samples_1000 <- matrix(NA, ncol=1000, nrow=120)

# iteratively generate 1000 samples
set.seed(2023) # this is a random operation; we pick seed 2023 for reproducibility
for (smpl in 1:1000){
  cat("Generating sample ", smpl, "\n", sep="")
  diversity_samples_lgs <- slice(data.frame(lg=NA),0) # store the languages here
  diversity_samples_fams <- slice(data.frame(family=NA),0) # store the corresponding families here
  for (rg in 1:nrow(regions)){ # for each macroarea
    reg_fam <- unique(filter(lgs,macroarea==regions[rg,])$level1) # these are the families in this macroarea
    fam_add <- data.frame(family=reg_fam[sample(length(reg_fam),size=20,replace=F)]) # sample 20 different random families in the macroarea
    diversity_samples_fams<-rbind(diversity_samples_fams,fam_add) # log the 20 families for this sample
    diversity_samples_lgs<-rbind(diversity_samples_lgs, # sample and save one random language per family
                                 data.frame(lg=(apply(fam_add,1,function(x)
                                   filter(lgs,level1==x[1]&macroarea==regions[rg,])[sample(nrow(filter(lgs,level1==x[1]&macroarea==regions[rg,])),size=1,replace=F),]$id))))
  }
  lg_samples_1000[,smpl]<-diversity_samples_lgs[,"lg"] # at the end of each iteration, store the full language sample
}

# store the diversity samples
write.csv(lg_samples_1000,"output/1000_diversity_samples_seed2023.csv")

### rerun all tests to reproduce stored results using the input data, expectations, taxonomy and language samples
# specify test condition, fix threshold
proportion_languages_must_be_in_applicable_state <- 1/3 # for a sample to be considered for evaluation, at least 1/3 of languages need to be coded for the second variable given the relevant state(s) of the first

# run all tests
results <- evaluate_XOR_AND_THEN(recoded_data = recoded_data, 
                                 expectations = expectations, 
                                 taxonomy = taxonomy, 
                                 diversity_samples = lg_samples_1000,
                                 proportion_languages_must_be_in_applicable_state = proportion_languages_must_be_in_applicable_state)

results <- na_convert(results)

# impute maximum and minimum values for mu+sigma and Cohen's D
results$maximum.sum <- apply(results,1,function(x)if(is.na(x[13])){x[10]}else if(x[10]>x[13]){x[10]}else{x[13]})
results$minimum.cohensd <- apply(results,1,function(x)if(is.na(x[17])){x[16]}else if(x[16]>x[17]){x[17]}else{x[16]})
results$minimum.sum <- apply(results,1,function(x)if(is.na(x[13])){x[10]}else if(x[10]>x[13]){x[13]}else{x[10]})
results$maximum.cohensd <- apply(results,1,function(x)if(is.na(x[17])){x[16]}else if(x[16]>x[17]){x[16]}else{x[17]})

# check the results logged in the file are correct
expect_true(all(round(results,digits=7) == round(expectations[,26:46],digits=7), na.rm = T))

# impute summary values ("yes" vs. "no") and verify they are correct
summaries1 <-  data.frame(y=apply(results,1,function(x)if(is.na(x[21])){"no"}else if(x[21]<1.3){"no"}else{"yes"}),
                         b=apply(results,1,function(x)if(x[4]>0.3){"yes"}else{"no"}),
                         c=apply(results,1,function(x)if(x[7]>19){"yes"}else{"no"}),
                         d=apply(results,1,function(x)if(is.na(x[18])){"no"}else if(x[18]>0.2){"no"}else{"yes"}),
                         e=apply(results,1,function(x)if(is.na(x[19])){"no"}else if(x[19]<1.3){"no"}else{"yes"}))

summaries2 <- data.frame(x=apply(summaries1,1,function(x)if(x[2]==x[3]&x[3]==x[4]&x[4]==x[5]&x[5]=="yes"){"yes"}else{"no"}))

expect_true(all(cbind(summaries2,summaries1) == expectations[,8:13]))
