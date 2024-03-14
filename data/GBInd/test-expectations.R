#### This script generates 500 diversity samples of 120 languages from the languages in the GBInd logical dataset.
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
logical_data <- read.csv("output/logicalGBI/logicalGBI.csv") %>% select(-X)
taxonomy <- read.csv("output/logicalGBI/cldf/languages.csv")

# subset taxonomy to languages in logical data with a macroarea-label
lgs <- filter(taxonomy,glottocode%in%logical_data$glottocode)

# load data for stats
recoded_data <- read.csv("output/statisticalGBI/data_for_stats.csv") %>% select(-X)

### create 500 diversity samples; each row of diversity_samples is a diversity sample of 120 lgs; 1 from each family (glottolog-based) ###
# list all six glottolog macroareas
regions <- data.frame(regions=c("Africa","Eurasia","Papunesia","Australia","South America","North America"))

# create empty matrix of 120 rows (languages) for 500 samples (columns)
lg_samples_500<-matrix(NA, ncol=500, nrow=120)

# iteratively generate 500 samples
set.seed(25) # this is a random operation; we pick seed 25 for reproducibility
for (smpl in 1:500){
  cat("Generating sample ", smpl, "\n", sep="")
  diversity_samples_lgs <- slice(data.frame(lg=NA),0) # store the languages here
  diversity_samples_fams <- slice(data.frame(family=NA),0) # store the corresponding families here
  for (rg in 1:nrow(regions)){ # for each macroarea
    reg_fam <- unique(filter(lgs,glottolog.macroarea==regions[rg,])$level1) # these are the families in this macroarea
    fam_add <- data.frame(family=reg_fam[sample(length(reg_fam),size=20,replace=F)]) # sample 20 different random families in the macroarea
    diversity_samples_fams<-rbind(diversity_samples_fams,fam_add) # log the 20 families for this sample
    diversity_samples_lgs<-rbind(diversity_samples_lgs, # sample and save one random language per family
                                 data.frame(lg=(apply(fam_add,1,function(x)
                                   filter(lgs,level1==x[1]&glottolog.macroarea==regions[rg,])[sample(nrow(filter(lgs,level1==x[1]&glottolog.macroarea==regions[rg,])),size=1,replace=F),]$glottocode))))
  }
  lg_samples_500[,smpl]<-diversity_samples_lgs[,"lg"] # at the end of each iteration, store the full language sample
}

# save this diversity sample
write.csv(lg_samples_500,"output/500_diversity_samples25.csv")

### rerun all tests to reproduce stored results using the input data, expectations, taxonomy and language samples
# specify test condition, fix threshold
proportion_languages_must_be_in_applicable_state <- 1/3 # for a sample to be considered for evaluation, at least 1/3 of languages need to be coded for the second variable given the relevant state(s) of the first

# run all tests
results <- evaluate_XOR_AND_THEN(recoded_data = recoded_data, 
                                 expectations = expectations, 
                                 taxonomy = taxonomy, 
                                 diversity_samples = lg_samples_500,
                                 proportion_languages_must_be_in_applicable_state = proportion_languages_must_be_in_applicable_state)

results <- na_convert(results)

# impute maximum and minimum values for mu+sigma and Cohen's D
results$maximum.sum <- apply(results,1,function(x)if(is.na(x[13])){x[10]}else if(x[10]>x[13]){x[10]}else{x[13]})
results$minimum.cohensd <- apply(results,1,function(x)if(is.na(x[17])){x[16]}else if(x[16]>x[17]){x[17]}else{x[16]})
results$minimum.sum <- apply(results,1,function(x)if(is.na(x[13])){x[10]}else if(x[10]>x[13]){x[13]}else{x[10]})
results$maximum.cohensd <- apply(results,1,function(x)if(is.na(x[17])){x[16]}else if(x[16]>x[17]){x[16]}else{x[17]})

## this should have reproduced all results
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