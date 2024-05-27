#### This script reads the Grambank data, parses it via all documented curation steps and generates and saves the 
#### GBInd logical and statistical datasets. Before saving, logical quality checks are performed on the outputs.

rm(list=ls())

# load relevant libraries
library(tidyverse)
library(testthat)
library(readr)
library(densify)
library(reshape2)
library(data.table)

# load functions
source("../functions.R")

########## load and prepare data ########## 
# read in original grambank data
original_feature_matrix <- as.data.frame(read_csv("../../raw/grambank_v.1.0.3.csv", 
                                                   trim_ws = FALSE, col_types = "f"))

names(original_feature_matrix)[1]<-"glottocode"

# replace missing data by ? (--> because these data points are unknown, not "not applicable")
original_feature_matrix[is.na(original_feature_matrix)] <- "?"

# read in manual language meta-data, check all languages are documented
lang_metadata <- read.csv("../lang-metadata.csv")
expect_true(all(original_feature_matrix$glottocode %in% lang_metadata$glottocode))

# read in the file specifying the maintained features, the and the recoded features with their recoding patterns
recode_patterns <- read.csv("input/feature-recode-patterns.csv")

########## parse all recodings in the appropriate order ########## 
## include without modification ##
# extract the features that we don't need to recode because we recode.operation.type them as they are
retained <- filter(recode_patterns, recode.operation.type=="include without modification")
retained_data <- select(original_feature_matrix,c(glottocode, filter(recode_patterns, recode.operation.type=="include without modification")$`original.names`))

# change names to new.names
for(i in 2:ncol(retained_data)){
  names(retained_data)[i]<-filter(retained,`original.names`==names(retained_data)[i])$new.name
}

recode_patterns <- filter(recode_patterns, recode.operation.type!="include without modification")

## recode group 1 (simple recode) ##
# subset to features requiring simple recoding only
first_set <- filter(recode_patterns, recode.operation.type=="recode group 1 (simple recode)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 1 (simple recode)")

# recode all features that require simple recoding
first_set_rec <- rowwise(first_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original feature is present in the original feature matrix
  expect_true(.$`original.names` %in% names(original_feature_matrix)[-1])
  original_data <- na.omit(as.character(original_feature_matrix[[.$`original.names`]]))
  
  # extract relevant attributes for recoding
  expected_states <- unlist(strsplit(unlist(.$original.states), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_states <- unlist(strsplit(.$new.states, ";"))
  
  # recode (single feature)
  new_data <- implement_recode(original_data, expected_states, recoding_groups, recoded_states, 
                               nvar="single", recode_mode="simple")
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = na.omit(original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%.$`original.names`))])$glottocode, 
    value = new_data,  
    stringsAsFactors=FALSE)  
}) %>%
  spread(feature, value)

# save the recoded data from retained features and the simple recodings as recoded_data for further use
recoded_data <- full_join(first_set_rec, retained_data, by=c(glottocode="glottocode"))

## recode group 2 (merge features - recode via logical arguments) ##
# subset to features that have to be merged via logical arguments
second_set <- filter(recode_patterns, recode.operation.type=="recode group 2 (merge features - recode via logical arguments)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 2 (merge features - recode via logical arguments)")

# merge and recode features via logical arguments
second_set_rec <- rowwise(second_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original features are present in the original feature matrix
  original_data <- original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%unlist(strsplit(.$`original.names`,"&"))))]
  
  if (ncol(original_data) != length(unlist(strsplit(.$`original.names`,"&")))+1){
    original_data<-cbind(original_data,recoded_data[,which(names(recoded_data)%in%unlist(strsplit(.$`original.names`,"&")))])
  }
  
  # extract relevant attributes for recoding
  expected_states <- unlist(strsplit(unlist(.$recode.operation.Rcode), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_states <- unlist(strsplit(.$new.states, ";"))
  
  # recode (multiple features)
  new_data <- implement_recode(original_data, expected_states, recoding_groups, recoded_states, 
                               nvar="multiple", recode_mode="logical_arguments")
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = original_data$glottocode, 
    value = new_data,  
    stringsAsFactors=FALSE)  
  
}) %>%
  spread(feature, value)

# merge new features with recoded_data for further use
recoded_data <- full_join(second_set_rec, recoded_data, by=c(glottocode="glottocode"))

## recode group 3 (merge features - recode via logical arguments - simple conditioning) ##
# subset to features that have to be recoded via logical arguments if a condition applies
third_set <- filter(recode_patterns, recode.operation.type=="recode group 3 (merge features - recode via logical arguments - simple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 3 (merge features - recode via logical arguments - simple conditioning)")

# merge and recode via logical arguments if a condition applies
third_set_rec <- rowwise(third_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original features are present in the original feature matrix
  original_data <- original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%unlist(strsplit(.$`original.names`,"&"))))]
  
  if (ncol(original_data) != length(unlist(strsplit(.$`original.names`,"&")))+1){
    original_data<-cbind(original_data,recoded_data[,which(names(recoded_data)%in%unlist(strsplit(.$`original.names`,"&")))])
  }
  
  # extract relevant attributes for recoding
  expected_states <- unlist(strsplit(unlist(.$recode.operation.Rcode), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_states <- unlist(strsplit(.$new.states, ";"))
  
  # recode (multiple features)
  new_data <- implement_recode(original_data, expected_states, recoding_groups, recoded_states, nvar="multiple", recode_mode="logical_arguments")
  
  # create table (unconditioned)
  unconditioned <- data.frame(
    glottocode = original_data$glottocode,
    value = new_data,
    stringsAsFactors=FALSE)
  
  # select condition statement
  condition_statement <-  .$condition.if.feature.conditioned
  
  # extract condition and equator
  condition_and_equator <- extract_condition_and_equator(condition_statement)
  
  # apply condition
  conditioned_data <- implement_conditioning(unconditioned, 
                                             condition=condition_and_equator[[1]], 
                                             equator=condition_and_equator[[2]])
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = conditioned_data$glottocode, 
    value = as.character(conditioned_data[,2]),  
    stringsAsFactors=FALSE)  
  
}) %>%
  spread(feature, value)

# merge new features with recoded_data for further use
recoded_data <- full_join(third_set_rec, recoded_data, by=c(glottocode="glottocode"))

## recode group 4 (merge features - recode via logical arguments - multiple conditioning) ##
# subset to features that have to be recoded via logical arguments if several conditions apply
fourth_set <- filter(recode_patterns, recode.operation.type=="recode group 4 (merge features - recode via logical arguments - multiple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 4 (merge features - recode via logical arguments - multiple conditioning)")

# merge and recode via logical arguments if a condition applies
fourth_set_rec <- rowwise(fourth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original features are present in the original feature matrix
  original_data <- original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%unlist(strsplit(.$`original.names`,"&"))))]
  
  if (ncol(original_data) != length(unlist(strsplit(.$`original.names`,"&")))+1){
    original_data<-cbind(original_data,recoded_data[,which(names(recoded_data)%in%unlist(strsplit(.$`original.names`,"&")))])
  }
  
  # extract relevant attributes for recoding
  expected_states <- unlist(strsplit(unlist(.$recode.operation.Rcode), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_states <- unlist(strsplit(.$new.states, ";"))
  
  # recode (multiple features)
  new_data <- implement_recode(original_data, expected_states, recoding_groups, recoded_states, 
                               nvar="multiple", recode_mode="logical_arguments")
  
  # create table (unconditioned)
  unconditioned <- data.frame(
    glottocode = original_data$glottocode,
    value = new_data,
    stringsAsFactors=FALSE)
  
  # select conditions
  conditions <- unlist(strsplit(.$condition.if.feature.conditioned," & "))
  nr_conditions <- length(conditions)
  
  # determine and extract conditions and equators for condition 1, apply condition 1
  condition_statement_1 <- conditions[1]
  condition_and_equator_1 <- extract_condition_and_equator(condition_statement_1)
  conditioned_data <- implement_conditioning(unconditioned, 
                                             condition=condition_and_equator_1[[1]], 
                                             equator=condition_and_equator_1[[2]])
  
  # determine and extract conditions and equators for condition 2, apply condition 2
  condition_statement_2 <- conditions[2]
  condition_and_equator_2 <- extract_condition_and_equator(condition_statement_2)
  conditioned_data <- implement_conditioning(conditioned_data, 
                                             condition=condition_and_equator_2[[1]], 
                                             equator=condition_and_equator_2[[2]])
  
  # determine and extract conditions and equators for condition 3, apply condition 3, if there are more than 2 conditions
  if (nr_conditions>2){
    condition_statement_3 <- conditions[3]
    condition_and_equator_3 <- extract_condition_and_equator(condition_statement_3)
    conditioned_data <- implement_conditioning(conditioned_data, 
                                               condition=condition_and_equator_3[[1]], 
                                               equator=condition_and_equator_3[[2]])
  }
  
  # determine and extract conditions and equators for condition 4, apply condition 4, if there are more than 3 conditions
  if (nr_conditions>3){
    condition_statement_4 <- conditions[4]
    condition_and_equator_4 <- extract_condition_and_equator(condition_statement_4)
    conditioned_data <- implement_conditioning(conditioned_data, 
                                               condition=condition_and_equator_4[[1]], 
                                               equator=condition_and_equator_4[[2]])
  }
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = conditioned_data$glottocode, 
    value = as.character(conditioned_data[,2]),  
    stringsAsFactors=FALSE)  
}) %>%
  spread(feature, value)

# merge new features with recoded_data for further use
recoded_data <- full_join(fourth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## recode group 5 (simple conditioning) ##
# subset to features that have to be conditioned on another feature
fifth_set <- filter(recode_patterns, recode.operation.type=="recode group 5 (simple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 5 (simple conditioning)")

# condition feature on another feature
fifth_set_rec <- rowwise(fifth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original feature is present in the original data
  expect_true(.$`original.names` %in% names(original_feature_matrix)[-1])
  feature_to_be_conditioned <- original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%.$`original.names`))]
  
  # select condition statement
  condition_statement <-  .$condition.if.feature.conditioned
  
  # extract condition and equator
  condition_and_equator <- extract_condition_and_equator(condition_statement)
  
  # apply condition
  conditioned_data <- implement_conditioning(feature_to_be_conditioned, condition=condition_and_equator[[1]], equator=condition_and_equator[[2]])
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = conditioned_data$glottocode, 
    value = as.character(conditioned_data[,2]),  
    stringsAsFactors=FALSE)  
  
}) %>%
  spread(feature, value)

# merge new features with recoded_data for further use
recoded_data <- full_join(fifth_set_rec, recoded_data, by=c(glottocode="glottocode"))


## recode group 6 (multiple conditioning) ## 
# subset to features that have to be conditioned on several features
sixth_set <- filter(recode_patterns, recode.operation.type=="recode group 6 (multiple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 6 (multiple conditioning)")

# condition on several features
sixth_set_rec <- rowwise(sixth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original feature is present in the original data
  expect_true(.$`original.names` %in% names(original_feature_matrix)[-1])
  feature_to_be_conditioned <- original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%.$`original.names`))]
  
  # select conditions
  conditions <- unlist(strsplit(.$condition.if.feature.conditioned," & "))
  nr_conditions <- length(conditions)
  
  # determine and extract conditions and equators for condition 1, apply condition 1
  condition_statement_1 <- conditions[1]
  condition_and_equator_1 <- extract_condition_and_equator(condition_statement_1)
  conditioned_data <- implement_conditioning(feature_to_be_conditioned, 
                                             condition=condition_and_equator_1[[1]], 
                                             equator=condition_and_equator_1[[2]])
  
  # determine and extract conditions and equators for condition 2, apply condition 2
  condition_statement_2 <- conditions[2]
  condition_and_equator_2 <- extract_condition_and_equator(condition_statement_2)
  conditioned_data <- implement_conditioning(conditioned_data, 
                                             condition=condition_and_equator_2[[1]], 
                                             equator=condition_and_equator_2[[2]])
  
  # determine and extract conditions and equators for condition 3, apply condition 3, if there are more than 2 conditions
  if (nr_conditions>2){
    condition_statement_3 <- conditions[3]
    condition_and_equator_3 <- extract_condition_and_equator(condition_statement_3)
    conditioned_data <- implement_conditioning(conditioned_data, 
                                               condition=condition_and_equator_3[[1]], 
                                               equator=condition_and_equator_3[[2]])
  }
  
  # determine and extract conditions and equators for condition 4, apply condition 4, if there are more than 3 conditions
  if (nr_conditions>3){
    condition_statement_4 <- conditions[4]
    condition_and_equator_4 <- extract_condition_and_equator(condition_statement_4)
    conditioned_data <- implement_conditioning(conditioned_data, 
                                               condition=condition_and_equator_4[[1]], 
                                               equator=condition_and_equator_4[[2]])
  }
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = conditioned_data$glottocode, 
    value = as.character(conditioned_data[,2]),  
    stringsAsFactors=FALSE)  
  
}) %>%
  spread(feature, value)

# merge new features with recoded_data for further use
recoded_data <- full_join(sixth_set_rec, recoded_data, by=c(glottocode="glottocode"))


## recode group 7 (simple conditioning [conditioned feature]) ##
# subset to features that have to be conditioned on another conditioned feature
seventh_set <- filter(recode_patterns, recode.operation.type=="recode group 7 (simple conditioning [conditioned feature])")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 7 (simple conditioning [conditioned feature])")

# condition on conditioned feature
seventh_set_rec <- rowwise(seventh_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original feature is present in the original data
  expect_true(.$`original.names` %in% names(original_feature_matrix)[-1])
  feature_to_be_conditioned <- original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%.$`original.names`))]
  
  # select condition statement
  condition_statement <-  .$condition.if.feature.conditioned
  
  # extract condition and equator
  condition_and_equator <- extract_condition_and_equator(condition_statement)
  
  # apply condition
  conditioned_data <- implement_conditioning(feature_to_be_conditioned, condition=condition_and_equator[[1]], equator=condition_and_equator[[2]])
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = conditioned_data$glottocode, 
    value = as.character(conditioned_data[,2]),  
    stringsAsFactors=FALSE)  
  
}) %>%
  spread(feature, value)

# merge new features with recoded_data for further use
recoded_data <- full_join(seventh_set_rec, recoded_data, by=c(glottocode="glottocode"))

## recode group 8 (multiple conditioning [conditioned feature]) ## 
# subset to features that have to be conditioned on several features
eighth_set <- filter(recode_patterns, recode.operation.type=="recode group 8 (multiple conditioning [conditioned feature])")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 8 (multiple conditioning [conditioned feature])")

# condition on several features
eighth_set_rec <- rowwise(eighth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original feature is present in the original data
  expect_true(.$`original.names` %in% names(original_feature_matrix)[-1])
  feature_to_be_conditioned <- original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%.$`original.names`))]
  
  # select conditions
  conditions <- unlist(strsplit(.$condition.if.feature.conditioned," & "))
  nr_conditions <- length(conditions)
  
  # determine and extract conditions and equators for condition 1, apply condition 1
  condition_statement_1 <- conditions[1]
  condition_and_equator_1 <- extract_condition_and_equator(condition_statement_1)
  conditioned_data <- implement_conditioning(feature_to_be_conditioned, 
                                             condition=condition_and_equator_1[[1]], 
                                             equator=condition_and_equator_1[[2]])
  
  # determine and extract conditions and equators for condition 2, apply condition 2
  condition_statement_2 <- conditions[2]
  condition_and_equator_2 <- extract_condition_and_equator(condition_statement_2)
  conditioned_data <- implement_conditioning(conditioned_data, 
                                             condition=condition_and_equator_2[[1]], 
                                             equator=condition_and_equator_2[[2]])
  
  # determine and extract conditions and equators for condition 3, apply condition 3, if there are more than 2 conditions
  if (nr_conditions>2){
    condition_statement_3 <- conditions[3]
    condition_and_equator_3 <- extract_condition_and_equator(condition_statement_3)
    conditioned_data <- implement_conditioning(conditioned_data, 
                                               condition=condition_and_equator_3[[1]], 
                                               equator=condition_and_equator_3[[2]])
  }
  
  # determine and extract conditions and equators for condition 4, apply condition 4, if there are more than 3 conditions
  if (nr_conditions>3){
    condition_statement_4 <- conditions[4]
    condition_and_equator_4 <- extract_condition_and_equator(condition_statement_4)
    conditioned_data <- implement_conditioning(conditioned_data, 
                                               condition=condition_and_equator_4[[1]], 
                                               equator=condition_and_equator_4[[2]])
  }
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = conditioned_data$glottocode, 
    value = as.character(conditioned_data[,2]),  
    stringsAsFactors=FALSE)  
  
}) %>%
  spread(feature, value)

# merge new features with recoded_data for further use
recoded_data <- full_join(eighth_set_rec, recoded_data, by=c(glottocode="glottocode"))


## recode group 9 (merge features - recode via logical arguments - simple conditioning [conditioned feature]) ##
# subset to features that have to be recoded via logical arguments if a condition applies (conditioned feature)
ninth_set <- filter(recode_patterns, recode.operation.type=="recode group 9 (merge features - recode via logical arguments - simple conditioning [conditioned feature])")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 9 (merge features - recode via logical arguments - simple conditioning [conditioned feature])")

# merge and recode via logical arguments if a condition applies
ninth_set_rec <- rowwise(ninth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original features are present in the original feature matrix
  original_data <- original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%unlist(strsplit(.$`original.names`,"&"))))]
  
  if (ncol(original_data) != length(unlist(strsplit(.$`original.names`,"&")))+1){
    original_data<-cbind(original_data,recoded_data[,which(names(recoded_data)%in%unlist(strsplit(.$`original.names`,"&")))])
  }
  
  # extract relevant attributes for recoding
  expected_states <- unlist(strsplit(unlist(.$recode.operation.Rcode), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_states <- unlist(strsplit(.$new.states, ";"))
  
  # recode (multiple features)
  new_data <- implement_recode(original_data, expected_states, recoding_groups, recoded_states, nvar="multiple", recode_mode="logical_arguments")
  
  # create table (unconditioned)
  unconditioned <- data.frame(
    glottocode = original_data$glottocode,
    value = new_data,
    stringsAsFactors=FALSE)
  
  # select condition statement
  condition_statement <-  .$condition.if.feature.conditioned
  
  # extract condition and equator
  condition_and_equator <- extract_condition_and_equator(condition_statement)
  
  # apply condition
  conditioned_data <- implement_conditioning(unconditioned, 
                                             condition=condition_and_equator[[1]], 
                                             equator=condition_and_equator[[2]])
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = conditioned_data$glottocode, 
    value = as.character(conditioned_data[,2]),  
    stringsAsFactors=FALSE)  
  
}) %>%
  spread(feature, value)

# merge new features with recoded_data for further use
recoded_data <- full_join(ninth_set_rec, recoded_data, by=c(glottocode="glottocode"))


## recode group 10 (merge features - recode via logical arguments - multiple conditioning [conditioned feature]) ##
# subset to features that have to be recoded via logical arguments if a condition applies (conditioned feature)
tenth_set <- filter(recode_patterns, recode.operation.type=="recode group 10 (merge features - recode via logical arguments - multiple conditioning [conditioned feature])")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 10 (merge features - recode via logical arguments - multiple conditioning [conditioned feature])")

# merge and recode via logical arguments if a condition applies
tenth_set_rec <- rowwise(tenth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original features are present in the original feature matrix
  original_data <- original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%unlist(strsplit(.$`original.names`,"&"))))]
  
  if (ncol(original_data) != length(unlist(strsplit(.$`original.names`,"&")))+1){
    original_data<-cbind(original_data,recoded_data[,which(names(recoded_data)%in%unlist(strsplit(.$`original.names`,"&")))])
  }
  
  # extract relevant attributes for recoding
  expected_states <- unlist(strsplit(unlist(.$recode.operation.Rcode), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_states <- unlist(strsplit(.$new.states, ";"))
  
  # recode (multiple features)
  new_data <- implement_recode(original_data, expected_states, recoding_groups, recoded_states, 
                               nvar="multiple", recode_mode="logical_arguments")
  
  # create table (unconditioned)
  unconditioned <- data.frame(
    glottocode = original_data$glottocode,
    value = new_data,
    stringsAsFactors=FALSE)
  
  # select conditions
  conditions <- unlist(strsplit(.$condition.if.feature.conditioned," & "))
  nr_conditions <- length(conditions)
  
  # determine and extract conditions and equators for condition 1, apply condition 1
  condition_statement_1 <- conditions[1]
  condition_and_equator_1 <- extract_condition_and_equator(condition_statement_1)
  conditioned_data <- implement_conditioning(unconditioned, 
                                             condition=condition_and_equator_1[[1]], 
                                             equator=condition_and_equator_1[[2]])
  
  # determine and extract conditions and equators for condition 2, apply condition 2
  condition_statement_2 <- conditions[2]
  condition_and_equator_2 <- extract_condition_and_equator(condition_statement_2)
  conditioned_data <- implement_conditioning(conditioned_data, 
                                             condition=condition_and_equator_2[[1]], 
                                             equator=condition_and_equator_2[[2]])
  
  # determine and extract conditions and equators for condition 3, apply condition 3, if there are more than 2 conditions
  if (nr_conditions>2){
    condition_statement_3 <- conditions[3]
    condition_and_equator_3 <- extract_condition_and_equator(condition_statement_3)
    conditioned_data <- implement_conditioning(conditioned_data, 
                                               condition=condition_and_equator_3[[1]], 
                                               equator=condition_and_equator_3[[2]])
  }
  
  # determine and extract conditions and equators for condition 4, apply condition 4, if there are more than 3 conditions
  if (nr_conditions>3){
    condition_statement_4 <- conditions[4]
    condition_and_equator_4 <- extract_condition_and_equator(condition_statement_4)
    conditioned_data <- implement_conditioning(conditioned_data, 
                                               condition=condition_and_equator_4[[1]], 
                                               equator=condition_and_equator_4[[2]])
  }
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = conditioned_data$glottocode, 
    value = as.character(conditioned_data[,2]),  
    stringsAsFactors=FALSE)  
}) %>%
  spread(feature, value)


# merge new features with recoded_data for further use
recoded_data <- full_join(tenth_set_rec, recoded_data, by=c(glottocode="glottocode"))

# replace all NA as explicit "NA"
recoded_data[is.na(recoded_data)]<-"NA"

# this full set of all input and recoded features needs to be stored to perform statistical tests
write.csv(recoded_data, "output/statisticalGBI/data_for_stats.csv")

# subset full data into original layer; logical layer and statistical layer
recode_patterns <- read.csv("input/feature-recode-patterns.csv")
all_decisions <- read.csv("input/decisions-log.csv")
logical_decisions <- all_decisions %>% filter(modification.type %in% c("logical","design"))
statistical_decisions <- all_decisions %>% filter(modification.type == "statistical")

original_layer <- recode_patterns%>%filter(original.features=="TRUE")
logical_layer <- recode_patterns%>%filter(design.logical=="TRUE")
statistical_layer <- recode_patterns%>%filter(design.logical.statistical=="TRUE")

original_data <- recoded_data %>% select(c("glottocode",original_layer$new.name))
logical_data <- recoded_data %>% select(c("glottocode",logical_layer$new.name))
statistical_data <- recoded_data %>% select(c("glottocode",statistical_layer$new.name))


########## sanity checks ########## 
### check all original features that should be in the original_feature filter are in there and vice versa
original_names_is <- original_layer$new.name
original_names_should <- names(retained_data)[-1]
expect_true(all(original_names_is%in%original_names_should))
expect_true(all(original_names_should%in%original_names_is))

### logical dataset: check all original features that should be in the logical filter are in there and vice versa
design_add <- c(unique(unlist(str_split(filter(logical_decisions,is.added.feature.kept == "yes")$resulting.added.features,", "))),NA)
design_remove <- unique(unlist(str_split(logical_decisions$resulting.removed.features,", ")))
# logical dataset is: a) original features, plus b) all design/logical additions, minus c) all design/logical removals
logical_names_should <- unique(c(original_names_is,design_add))[unique(c(original_names_is,design_add))%in%design_remove==F]
logical_names_is <- names(logical_data)[names(logical_data)!="glottocode"]
# sanity checks
expect_true(all(logical_names_is %in% logical_names_should))
expect_true(all(logical_names_should %in% logical_names_is))

### statistical dataset: check all original features that should be in the statistical filter are in there and vice versa
statistical_add <- unique(statistical_decisions$resulting.added.features) 
statistical_remove <- unique(unlist(str_split(statistical_decisions$resulting.removed.features,", ")))
# statistical dataset is: a) the features from the final logical dataset, plus b) all statistical additions, minus c) all statistical removals
statistical_names_should <- unique(c(logical_names_is,statistical_add))[unique(c(logical_names_is,statistical_add))%in%statistical_remove==F]
statistical_names_is <- names(statistical_data)[names(statistical_data)!="glottocode"]
# sanity checks
expect_true(all(statistical_names_is %in% statistical_names_should))
expect_true(all(statistical_names_should %in% statistical_names_is))

### modification ID match --> ensure all modification IDs in the features sheet are in the modification sheet and vice versa
mod_IDs_is <- na.omit(unique(c(unlist(str_split(recode_patterns$modification.IDs,";")),(unlist(str_split(recode_patterns$associated.modification.IDs.without.resulting.action,";"))))))
mod_IDs_is <- mod_IDs_is[mod_IDs_is!=""]
mod_IDs_should <- na.omit(all_decisions$modification.ID)
expect_true(all(mod_IDs_is %in% mod_IDs_should))
expect_true(all(mod_IDs_should %in% mod_IDs_is))

# specific modification ID match --> ensure that each modification ID in the features sheet is in the modification sheet, associated via the correct columns and features; and vice versa
ids <- na.omit(all_decisions$modification.ID)
for (id in ids){
  type <- filter(all_decisions,modification.ID == id)$modification.type
  if (type == "statistical"){
    # check that each instance of a modification ID in the spreadsheet ("is") is foreseen in the decisions_log ("should") and vice versa
    should_all <- all_decisions %>% filter(modification.ID == id) %>% select(c("feature.1.for.test","feature.2.for.test","resulting.added.features","resulting.removed.features")) %>% as.character() %>% unique()
    should_all <- na.omit(unique(unlist(strsplit(should_all[should_all!="NA"],", "))))
    is_all <- recode_patterns %>% slice(c(which(grepl(id,recode_patterns$modification.IDs)),which(grepl(id,recode_patterns$associated.modification.IDs.without.resulting.action))))
    is_all <- is_all$new.name
    expect_true(all(is_all%in%should_all))
    expect_true(all(should_all%in%is_all))
    
    # check that each instance of a modification ID WITH EFFECT in the spreadsheet ("is") is foreseen in the decisions_log ("should") and vice versa
    should_actedupon <- all_decisions %>% filter(modification.ID == id) %>% select(c("resulting.added.features","resulting.removed.features")) %>% as.character() %>% unique()
    should_actedupon <- na.omit(unique(unlist(strsplit(should_actedupon[should_actedupon!="NA"],", "))))
    is_actedupon <- recode_patterns %>% slice(which(grepl(id,recode_patterns$modification.IDs)))
    is_actedupon <- is_actedupon$new.name
    expect_true(all(is_actedupon%in%should_actedupon))
    expect_true(all(should_actedupon%in%is_actedupon))
    
    # check that each instance of a modification ID WITHOUT EFFECT in the spreadsheet ("is") is foreseen in the decisions_log ("should") and vice versa
    should_associated <- all_decisions %>% filter(modification.ID == id) %>% select(c("feature.1.for.test","feature.2.for.test")) %>% as.character() %>% unique()
    should_associated <- setdiff(na.omit(unique(unlist(strsplit(should_associated[should_associated!="NA"],", ")))),is_actedupon)
    is_associated <- recode_patterns %>% slice(which(grepl(id,recode_patterns$associated.modification.IDs.without.resulting.action)))
    is_associated <- is_associated$new.name
    expect_true(all(is_associated%in%should_associated))
    expect_true(all(should_associated%in%is_associated))
  }
  else if (type %in% c("logical","design-automated","design-manual")){ 
    # check that each instance of a modification ID in the spreadsheet ("is") is foreseen in the decisions_log ("should") and vice versa
    should_all <- all_decisions %>% filter(modification.ID == id) %>% select(c("relevant.features","resulting.added.features","resulting.removed.features")) %>% as.character() %>% unique()
    should_all <- na.omit(unique(unlist(strsplit(should_all[should_all!="NA"],", "))))
    is_all <- recode_patterns %>% slice(c(which(grepl(id,recode_patterns$modification.IDs)),which(grepl(id,recode_patterns$associated.modification.IDs.without.resulting.action))))
    is_all <- is_all$new.name
    expect_true(all(is_all%in%should_all))
    expect_true(all(should_all%in%is_all))
    
    # check that each instance of a modification ID WITH EFFECT in the spreadsheet ("is") is foreseen in the decisions_log ("should") and vice versa
    should_actedupon <- all_decisions %>% filter(modification.ID == id) %>% select(c("resulting.added.features","resulting.removed.features")) %>% as.character() %>% unique()
    should_actedupon <- na.omit(unique(unlist(strsplit(should_actedupon[should_actedupon!="NA"],", "))))
    is_actedupon <- recode_patterns %>% slice(which(grepl(id,recode_patterns$modification.IDs)))
    is_actedupon <- is_actedupon$new.name
    expect_true(all(is_actedupon%in%should_actedupon))
    expect_true(all(should_actedupon%in%is_actedupon))
    
    # check that each instance of a modification ID WITHOUT EFFECT in the spreadsheet ("is") is foreseen in the decisions_log ("should") and vice versa
    should_associated <- all_decisions %>% filter(modification.ID == id) %>% select("relevant.features") %>% as.character() %>% unique()
    should_associated <- setdiff(na.omit(unique(unlist(strsplit(should_associated[should_associated!="NA"],", ")))),is_actedupon)
    is_associated <- recode_patterns %>% slice(which(grepl(id,recode_patterns$associated.modification.IDs.without.resulting.action)))
    is_associated <- is_associated$new.name
    expect_true(all(is_associated%in%should_associated))
    expect_true(all(should_associated%in%is_associated))
  }
}

## known.remaining dependencies match - check remaining dependencies are appropriately logged (one-sided, because there are dependencies not associated with statistical tests!)
rds <- filter(all_decisions, resulting.modification == "tendency logged in known.remaining.dependencies")
for (rd in 1:nrow(rds)){
  id <- rds %>% slice(rd) %>% select(modification.ID)
  should_associated <- rds %>% slice(rd) %>% select(c(feature.1.for.test,feature.2.for.test)) %>% as.character
  is_associated <- recode_patterns %>% slice(which(grepl(id,recode_patterns$known.remaining.dependencies.after.statistical.treatment))) %>% select(new.name) %>% unlist %>% as.character
  expect_true(all(is_associated%in%should_associated))
  expect_true(all(should_associated%in%is_associated))
}
expect_true(all(rds$modification.ID%in%unlist(strsplit(recode_patterns$known.remaining.dependencies.after.statistical.treatment,";"))))
expect_true(all(na.omit(unique(unlist(strsplit(recode_patterns$known.remaining.dependencies.after.statistical.treatment,";"))))%in%rds$modification.ID))

########## save data as language-feature matrices ########## 
# save logical and statistical datasets as language-feature matrices (.csv)
write.csv(logical_data,"output/logicalGBI/logicalGBI.csv")
write.csv(statistical_data,"output/statisticalGBI/statisticalGBI.csv")

########## make, check and save cldf  ########## 
# languages.csv
lang_metadata <- read.csv("../lang-metadata.csv")

taxonomy_logical <- data.frame(glottocode = logical_data$glottocode)
taxonomy_logical <- left_join(taxonomy_logical,lang_metadata)
taxonomy_logical[taxonomy_logical==""] <- NA

taxonomy_statistical <- data.frame(glottocode = statistical_data$glottocode)
taxonomy_statistical <- left_join(taxonomy_statistical,lang_metadata)
taxonomy_statistical[taxonomy_statistical==""] <- NA

# parameters.csv
parameters <- read.csv("input/feature-recode-patterns.csv")

parameters_logical <- parameters %>% filter(design.logical==T)
parameters_statistical <- parameters %>% filter(design.logical.statistical==T)

# values.csv
logical_long <- melt(setDT(logical_data), id.vars = "glottocode", variable.name = "new.name")
logical_long$value_ID <- apply(logical_long,1,function(x) paste(x[2],x[1],sep="-"))
logical_long$code_ID <- apply(logical_long,1,function(x) paste(x[2],x[3],sep="-"))
logical_long <- logical_long %>% select(c("value_ID","glottocode","new.name","value","code_ID"))
logical_long$glottocode <- as.character(logical_long$glottocode)
logical_long$new.name <- as.character(logical_long$new.name)

statistical_long <- melt(setDT(statistical_data), id.vars = "glottocode", variable.name = "new.name")
statistical_long$value_ID <- apply(statistical_long,1,function(x) paste(x[2],x[1],sep="-"))
statistical_long$code_ID <- apply(statistical_long,1,function(x) paste(x[2],x[3],sep="-"))
statistical_long$glottocode <- as.character(statistical_long$glottocode)
statistical_long$new.name <- as.character(statistical_long$new.name)
statistical_long <- statistical_long %>% select(c("value_ID","glottocode","new.name","value","code_ID"))

# codes.csv
logical_codes <- logical_long %>% select(c("code_ID","new.name","value")) %>% unique()

statistical_codes <- statistical_long %>% select(c("code_ID","new.name","value")) %>% unique()

# modifications.csv
modifications <- read.csv("input/decisions-log.csv")

# cldf quality checks:
expect_true(all(unique(logical_long$glottocode) %in% taxonomy_logical$glottocode))
expect_true(all(unique(statistical_long$glottocode) %in% taxonomy_statistical$glottocode))
expect_true(all(taxonomy_logical$glottocode %in% unique(logical_long$glottocode)))
expect_true(all(taxonomy_statistical$glottocode %in% unique(statistical_long$glottocode)))

expect_true(all(unique(logical_long$new.name) %in% parameters$new.name))
expect_true(all(unique(statistical_long$new.name) %in% parameters$new.name))
expect_true(all(filter(parameters,design.logical == "TRUE")$new.name %in% unique(logical_long$new.name)))
expect_true(all(filter(parameters,design.logical.statistical == "TRUE")$new.name %in% unique(statistical_long$new.name)))

expect_true(all(unique(logical_long$code_ID) %in% logical_codes$code_ID))
expect_true(all(unique(statistical_long$code_ID) %in% statistical_codes$code_ID))
expect_true(all(logical_codes$code_ID %in% unique(logical_long$code_ID)))
expect_true(all(statistical_codes$code_ID %in% unique(statistical_long$code_ID)))

expect_true(all(unique(unlist(strsplit(filter(parameters,modification.IDs!="")$modification.IDs,";"))) %in% modifications$modification.ID))
expect_true(all(unique(unlist(strsplit(filter(parameters,associated.modification.IDs.without.resulting.action!="")$associated.modification.IDs.without.resulting.action,";"))) %in% modifications$modification.ID))
expect_true(all(modifications$modification.ID %in% c(unique(unlist(strsplit(filter(parameters,modification.IDs!="")$modification.IDs,";"))),
                                                     unique(unlist(strsplit(filter(parameters,associated.modification.IDs.without.resulting.action!="")$associated.modification.IDs.without.resulting.action,";"))))))

# write all cldf components:
write.csv(taxonomy_logical,"output/logicalGBI/cldf/languages.csv",fileEncoding="UTF-8",row.names = F)
write.csv(taxonomy_statistical,"output/statisticalGBI/cldf/languages.csv",fileEncoding="UTF-8",row.names = F)
write.csv(parameters_logical,"output/logicalGBI/cldf/parameters.csv",fileEncoding="UTF-8",row.names = F)
write.csv(parameters_statistical,"output/statisticalGBI/cldf/parameters.csv",fileEncoding="UTF-8",row.names = F)
write.csv(logical_long,"output/logicalGBI/cldf/values.csv",fileEncoding="UTF-8",row.names = F)
write.csv(statistical_long,"output/statisticalGBI/cldf/values.csv",fileEncoding="UTF-8",row.names = F)
write.csv(logical_codes,"output/logicalGBI/cldf/codes.csv",fileEncoding="UTF-8",row.names = F)
write.csv(statistical_codes,"output/statisticalGBI/cldf/codes.csv",fileEncoding="UTF-8",row.names = F)
write.csv(modifications,"output/logicalGBI/cldf/modifications.csv",fileEncoding="UTF-8",row.names = F)
write.csv(modifications,"output/statisticalGBI/cldf/modifications.csv",fileEncoding="UTF-8",row.names = F)
