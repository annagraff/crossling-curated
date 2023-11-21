rm(list=ls())

# load relevant libraries
library(tidyverse)
library(testthat)
library(readr)
library(densify)

# load functions
source("../functions.R")

# read in original variables and taxonomy
original_variable_matrix <- as.data.frame(read_csv("output/merged_original_variables.csv", 
                                                   trim_ws = FALSE, col_types = "f"))

taxonomy <- as_flat_taxonomy_matrix(glottolog_languoids)

# some glottocodes from WALS not found in glottolog v. 4.8. because they have been assigned new glottocodes - replace them here
original_variable_matrix$glottocode[which(original_variable_matrix$glottocode %in% taxonomy$id == F)]

original_variable_matrix <- update_glottocode_and_data(data=original_variable_matrix, old="woro1255", new="worr1237") # Worora from WALS, manually assigned via name
original_variable_matrix <- update_glottocode_and_data(data=original_variable_matrix, old="jogl1236", new="tase1235") # Jugli from WALS, manually assigned via name, references and iso-code
original_variable_matrix <- update_glottocode_and_data(data=original_variable_matrix, old="sanm1259", new="sanm1295") # Mixtec (Molinos) from WALS, manually assigned via name, references and iso-code
original_variable_matrix <- update_glottocode_and_data(data=original_variable_matrix, old="cypr1245", new="mode1248") # Greek (Cypriot) from WALS, manually assigned via name, references and iso-code
original_variable_matrix <- update_glottocode_and_data(data=original_variable_matrix, old="yidd1255", new="east2295") # Mixtec (Molinos) from WALS, manually assigned via name, references and iso-code

expect_true(all(original_variable_matrix$glottocode %in% taxonomy$id))

# read in the file specifying the maintained variables, the and the recoded variables with their recoding patterns
recode_patterns <- read.csv("input/recode-patterns-variables.csv")

########## include without modification ########## 
# extract the features that we don't need to recode because we recode.operation.type them as they are
retained <- filter(recode_patterns, recode.operation.type=="include without modification")
retained_data <- select(original_variable_matrix,c(glottocode, filter(recode_patterns, recode.operation.type=="include without modification")$`original.names`))

# change names to new.names
for(i in 2:ncol(retained_data)){
  names(retained_data)[i]<-filter(retained,`original.names`==names(retained_data)[i])$new.name
}

recode_patterns <- filter(recode_patterns, recode.operation.type!="include without modification")

########## recode group 1 (simple recode) ########## 
# subset to features requiring simple recoding only
first_set <- filter(recode_patterns, recode.operation.type=="recode group 1 (simple recode)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 1 (simple recode)")

# recode all features that require simple recoding
first_set_rec <- rowwise(first_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variable is present in the original variable matrix
  expect_true(.$`original.names` %in% names(original_variable_matrix)[-1])
  original_data <- na.omit(as.character(original_variable_matrix[[.$`original.names`]]))
  
  # extract relevant attributes for recoding
  expected_levels <- unlist(strsplit(unlist(.$original.levels), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_levels <- unlist(strsplit(.$new.levels, ";"))
  
  # recode (single variable)
  new_data <- implement_recode(original_data, expected_levels, recoding_groups, recoded_levels, 
                               nvar="single", recode_mode="simple")
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = na.omit(original_variable_matrix[,c(1,which(names(original_variable_matrix)%in%.$`original.names`))])$glottocode, 
    value = new_data,  
    stringsAsFactors=FALSE)  
}) %>%
  spread(feature, value)

# save the recoded data from retained variables and the simple recodings as recoded_data for further use
recoded_data <- full_join(first_set_rec, retained_data, by=c(glottocode="glottocode"))

## !! note that variable priorisation is currently implemented for up to 5 variables !!
########## recode group 2 (merge variables - simple recode) ########## 
# subset to features that are merged
second_set <- filter(recode_patterns, recode.operation.type=="recode group 2 (merge variables - simple recode)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 2 (merge variables - simple recode)")

# merge and recode variables
second_set <- rowwise(second_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variables are present in the original variable matrix
  expect_true(unique(unlist(strsplit(.$`original.names`,"&")) %in% names(original_variable_matrix)[-1]))
  original_data <- original_variable_matrix[,c("glottocode", unlist(strsplit(.$`original.names`,"&")))]
  
  # define starting variable states according to merging priorities
  original_data$merged <- ifelse((original_data[,.$first.priority.variable.if.merged] != "?"), 
                                 as.character(original_data[,.$first.priority.variable.if.merged]), 
                                 as.character(original_data[,.$second.priority.variable.if.merged]))
  
  if(ncol(original_data)==5){ # if there is a third variable
    original_data$merged <- ifelse((original_data$merged != "?"), 
                                   as.character(original_data$merged),
                                   as.character(original_data[,.$third.priority.variable.if.merged]))
  }
  
  if(ncol(original_data)==6){ # if there is a fourth variable
    original_data$merged <- ifelse((original_data$merged != "?"), 
                                   as.character(original_data$merged),
                                   as.character(original_data[,.$fourth.priority.variable.if.merged]))
  }
  
  original_data <- na.omit(select(original_data,c(glottocode,merged)))
  
  # extract relevant attributes for recoding
  expected_levels <- unlist(strsplit(unlist(.$original.levels), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_levels <- unlist(strsplit(.$new.levels, ";"))
  
  # recode (multiple variables)
  new_data <- implement_recode(original_data, expected_levels, recoding_groups, recoded_levels, 
                               nvar="multiple", recode_mode="simple")
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = original_data$glottocode, 
    value = new_data,  
    stringsAsFactors=FALSE)  
  
}) %>%
  spread(feature, value)

# merge new variables with recoded_data for further use
recoded_data <- full_join(second_set, recoded_data, by=c(glottocode="glottocode"))

########## recode group 3 (merge variables - recode via logical arguments) ########## 
# subset to features that have to be merged via logical arguments
third_set <- filter(recode_patterns, recode.operation.type=="recode group 3 (merge variables - recode via logical arguments)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 3 (merge variables - recode via logical arguments)")

# merge and recode variables via logical arguments
third_set_rec <- rowwise(third_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variables are present in the original variable matrix
  original_data <- original_variable_matrix[,c(1,which(names(original_variable_matrix)%in%unlist(strsplit(.$`original.names`,"&"))))]
  
  if (ncol(original_data) != length(unlist(strsplit(.$`original.names`,"&")))+1){
    original_data<-cbind(original_data,recoded_data[,which(names(recoded_data)%in%unlist(strsplit(.$`original.names`,"&")))])
  }
  
  # extract relevant attributes for recoding
  expected_levels <- unlist(strsplit(unlist(.$recode.operation.Rcode), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_levels <- unlist(strsplit(.$new.levels, ";"))
  
  # recode (multiple variables)
  new_data <- implement_recode(original_data, expected_levels, recoding_groups, recoded_levels, 
                               nvar="multiple", recode_mode="logical_arguments")
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = original_data$glottocode, 
    value = new_data,  
    stringsAsFactors=FALSE)  
  
}) %>%
  spread(feature, value)

# merge new variables with recoded_data for further use
recoded_data <- full_join(third_set_rec, recoded_data, by=c(glottocode="glottocode"))

########## recode group 4 (recode via logical arguments) ########## 
# subset to features that have to be merged via logical arguments
fourth_set <- filter(recode_patterns, recode.operation.type=="recode group 4 (recode via logical arguments)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 4 (recode via logical arguments)")

# recode variables via logical arguments
fourth_set_rec <- rowwise(fourth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variables are present in the original variable matrix
  original_data<-original_variable_matrix[,c(1,which(names(original_variable_matrix)%in%unlist(strsplit(.$`original.names`,"&"))))]
  
  # extract relevant attributes for recoding
  expected_levels <- unlist(strsplit(unlist(.$recode.operation.Rcode), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_levels <- unlist(strsplit(.$new.levels, ";"))
  
  # recode (multiple variables)
  new_data <- implement_recode(original_data, expected_levels, recoding_groups, recoded_levels, nvar="single", recode_mode="logical_arguments")
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = original_data$glottocode, 
    value = new_data,  
    stringsAsFactors=FALSE)  
  
}) %>%
  spread(feature, value)

# merge new variables with recoded_data for further use
recoded_data <- full_join(fourth_set_rec, recoded_data, by=c(glottocode="glottocode"))

########## recode group 5 (simple conditioning) ########## 
# subset to features that have to be conditioned on another feature
fifth_set <- filter(recode_patterns, recode.operation.type=="recode group 5 (simple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 5 (simple conditioning)")

# condition variable on another variable
fifth_set_rec <- rowwise(fifth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variable is present in the original data
  expect_true(.$`original.names` %in% names(original_variable_matrix)[-1])
  variable_to_be_conditioned <- original_variable_matrix[,c(1,which(names(original_variable_matrix)%in%.$`original.names`))]
  
  # select condition statement
  condition_statement <-  .$condition.if.variable.conditioned
  
  # extract condition and equator
  condition_and_equator <- extract_condition_and_equator(condition_statement)
  
  # apply condition
  conditioned_data <- implement_conditioning(variable_to_be_conditioned, condition=condition_and_equator[[1]], equator=condition_and_equator[[2]])
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = conditioned_data$glottocode, 
    value = as.character(conditioned_data[,2]),  
    stringsAsFactors=FALSE)  
  
}) %>%
  spread(feature, value)

# merge new variables with recoded_data for further use
recoded_data <- full_join(fifth_set_rec, recoded_data, by=c(glottocode="glottocode"))

########## recode group 6 (merge variables - recode via logical arguments - simple conditioning) ########## 
# subset to features that have to be recoded via logical arguments if a condition applies
sixth_set <- filter(recode_patterns, recode.operation.type=="recode group 6 (merge variables - recode via logical arguments - simple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 6 (merge variables - recode via logical arguments - simple conditioning)")

# merge and recode via logical arguments if a condition applies
sixth_set_rec <- rowwise(sixth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variables are present in the original variable matrix
  original_data <- original_variable_matrix[,c(1,which(names(original_variable_matrix)%in%unlist(strsplit(.$`original.names`,"&"))))]
  
  if (ncol(original_data) != length(unlist(strsplit(.$`original.names`,"&")))+1){
    original_data<-cbind(original_data,recoded_data[,which(names(recoded_data)%in%unlist(strsplit(.$`original.names`,"&")))])
  }
  
  # extract relevant attributes for recoding
  expected_levels <- unlist(strsplit(unlist(.$recode.operation.Rcode), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_levels <- unlist(strsplit(.$new.levels, ";"))
  
  # recode (multiple variables)
  new_data <- implement_recode(original_data, expected_levels, recoding_groups, recoded_levels, nvar="multiple", recode_mode="logical_arguments")
  
  # create table (unconditioned)
  unconditioned <- data.frame(
    glottocode = original_data$glottocode,
    value = new_data,
    stringsAsFactors=FALSE)
  
  # select condition statement
  condition_statement <-  .$condition.if.variable.conditioned
  
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

# merge new variables with recoded_data for further use
recoded_data <- full_join(sixth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that multiple conditioning is currently only implemented for up to 4 conditions !!
########## recode group 7 (multiple conditioning) ########## 
# subset to features that have to be conditioned on several features
seventh_set <- filter(recode_patterns, recode.operation.type=="recode group 7 (multiple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 7 (multiple conditioning)")

# condition on several variables
seventh_set_rec <- rowwise(seventh_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variable is present in the original data
  expect_true(.$`original.names` %in% names(original_variable_matrix)[-1])
  variable_to_be_conditioned <- original_variable_matrix[,c(1,which(names(original_variable_matrix)%in%.$`original.names`))]
  
  # select conditions
  conditions <- unlist(strsplit(.$condition.if.variable.conditioned," & "))
  nr_conditions <- length(conditions)
  
  # determine and extract conditions and equators for condition 1, apply condition 1
  condition_statement_1 <- conditions[1]
  condition_and_equator_1 <- extract_condition_and_equator(condition_statement_1)
  conditioned_data <- implement_conditioning(variable_to_be_conditioned, 
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

# merge new variables with recoded_data for further use
recoded_data <- full_join(seventh_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that variable priorisation is currently implemented for up to 5 variables !!
########## recode group 8 (merge variables - simple recode - simple conditioning) ########## 
# subset to features that have to be created by merging and recoding, but only if a condition applies
eighth_set <- filter(recode_patterns, recode.operation.type=="recode group 8 (merge variables - simple recode - simple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 8 (merge variables - simple recode - simple conditioning)")

# merge and recode if condition applies
eighth_set_rec <- rowwise(eighth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variables are present in the original variable matrix
  expect_true(unique(unlist(strsplit(.$`original.names`,"&")) %in% names(original_variable_matrix)[-1]))
  original_data <- original_variable_matrix[,c("glottocode", unlist(strsplit(.$`original.names`,"&")))]
  
  # define starting variable states according to merging priorities
  original_data$merged <- ifelse((original_data[,.$first.priority.variable.if.merged] != "?"), 
                                 as.character(original_data[,.$first.priority.variable.if.merged]), 
                                 as.character(original_data[,.$second.priority.variable.if.merged]))
  
  if(ncol(original_data)==5){ # if there is a third variable
    original_data$merged <- ifelse((original_data$merged != "?"), 
                                   as.character(original_data$merged),
                                   as.character(original_data[,.$third.priority.variable.if.merged]))
  }
  
  if(ncol(original_data)==6){ # if there is a fourth variable
    original_data$merged <- ifelse((original_data$merged != "?"), 
                                   as.character(original_data$merged),
                                   as.character(original_data[,.$fourth.priority.variable.if.merged]))
  }
  
  original_data <- na.omit(select(original_data,c(glottocode,merged)))
  
  # extract relevant attributes for recoding
  expected_levels <- unlist(strsplit(unlist(.$original.levels), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_levels <- unlist(strsplit(.$new.levels, ";"))
  
  # recode (multiple variables)
  new_data <- implement_recode(original_data, expected_levels, recoding_groups, recoded_levels, 
                               nvar="multiple", recode_mode="simple")
  
  # prepare output as data frame (unconditioned)
  unconditioned <- data.frame(
    glottocode = original_data$glottocode,
    value = new_data,
    stringsAsFactors=FALSE)
  
  # select condition statement
  condition_statement <-  .$condition.if.variable.conditioned
  
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

# merge new variables with recoded_data for further use
recoded_data <- full_join(eighth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that variable priorisation is currently implemented for up to 5 variables !!
## !! note that multiple conditioning is currently only implemented for up to 4 conditions !!

########## recode group 9 (merge variables - simple recode - multiple conditioning) ########## 
# subset to features that have to be merged and recoded, but only if several conditions apply
ninth_set <- filter(recode_patterns, recode.operation.type=="recode group 9 (merge variables - simple recode - multiple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 9 (merge variables - simple recode - multiple conditioning)")

# merge and recode if several conditions apply
ninth_set_rec <- rowwise(ninth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variables are present in the original variable matrix
  expect_true(unique(unlist(strsplit(.$`original.names`,"&")) %in% names(original_variable_matrix)[-1]))
  original_data <- original_variable_matrix[,c("glottocode", unlist(strsplit(.$`original.names`,"&")))]
  
  # define starting variable states according to merging priorities
  original_data$merged <- ifelse((original_data[,.$first.priority.variable.if.merged] != "?"), 
                                 as.character(original_data[,.$first.priority.variable.if.merged]), 
                                 as.character(original_data[,.$second.priority.variable.if.merged]))
  
  if(ncol(original_data)==5){ # if there is a third variable
    original_data$merged <- ifelse((original_data$merged != "?"), 
                                   as.character(original_data$merged),
                                   as.character(original_data[,.$third.priority.variable.if.merged]))
  }
  
  if(ncol(original_data)==6){ # if there is a fourth variable
    original_data$merged <- ifelse((original_data$merged != "?"), 
                                   as.character(original_data$merged),
                                   as.character(original_data[,.$fourth.priority.variable.if.merged]))
  }
  
  original_data <- na.omit(select(original_data,c(glottocode,merged)))
  
  # extract relevant attributes for recoding
  expected_levels <- unlist(strsplit(unlist(.$original.levels), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_levels <- unlist(strsplit(.$new.levels, ";"))
  
  # recode (multiple variables)
  new_data <- implement_recode(original_data, expected_levels, recoding_groups, recoded_levels, 
                               nvar="multiple", recode_mode="simple")
  
  # prepare output as data frame (unconditioned)
  unconditioned <- data.frame(
    glottocode = original_data$glottocode,
    value = new_data,
    stringsAsFactors=FALSE)
  
  # select conditions
  conditions <- unlist(strsplit(.$condition.if.variable.conditioned," & "))
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

# merge new variables with recoded_data for further use
recoded_data <- full_join(ninth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that multiple conditioning is currently only implemented for up to 4 conditions !!

########## recode group 10 (merge variables - recode via logical arguments - multiple conditioning) ########## 
# subset to features that have to be recoded via logical arguments if several conditions apply
tenth_set <- filter(recode_patterns, recode.operation.type=="recode group 10 (merge variables - recode via logical arguments - multiple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 10 (merge variables - recode via logical arguments - multiple conditioning)")

# merge and recode via logical arguments if a condition applies
tenth_set_rec <- rowwise(tenth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variables are present in the original variable matrix
  original_data <- original_variable_matrix[,c(1,which(names(original_variable_matrix)%in%unlist(strsplit(.$`original.names`,"&"))))]
  
  if (ncol(original_data) != length(unlist(strsplit(.$`original.names`,"&")))+1){
    original_data<-cbind(original_data,recoded_data[,which(names(recoded_data)%in%unlist(strsplit(.$`original.names`,"&")))])
  }
  
  # extract relevant attributes for recoding
  expected_levels <- unlist(strsplit(unlist(.$recode.operation.Rcode), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_levels <- unlist(strsplit(.$new.levels, ";"))
  
  # recode (multiple variables)
  new_data <- implement_recode(original_data, expected_levels, recoding_groups, recoded_levels, 
                               nvar="multiple", recode_mode="logical_arguments")
  
  # create table (unconditioned)
  unconditioned <- data.frame(
    glottocode = original_data$glottocode,
    value = new_data,
    stringsAsFactors=FALSE)
  
  # select conditions
  conditions <- unlist(strsplit(.$condition.if.variable.conditioned," & "))
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

# merge new variables with recoded_data for further use
recoded_data <- full_join(tenth_set_rec, recoded_data, by=c(glottocode="glottocode"))

########## recode group 11 (simple recode - simple conditioning) ########## 
# subset to features that have to be recoded, but only if a condition applies
eleventh_set <- filter(recode_patterns, recode.operation.type=="recode group 11 (simple recode - simple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 11 (simple recode - simple conditioning)")

# recode if condition applies
eleventh_set_rec <- rowwise(eleventh_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variable is present in the original variable matrix
  expect_true(.$`original.names` %in% names(original_variable_matrix)[-1])
  original_data <- na.omit(as.character(original_variable_matrix[[.$`original.names`]]))
  
  # extract relevant attributes for recoding
  expected_levels <- unlist(strsplit(unlist(.$original.levels), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_levels <- unlist(strsplit(.$new.levels, ";"))
  
  # recode (single variable)
  new_data <- implement_recode(original_data, expected_levels, recoding_groups, recoded_levels, 
                               nvar="single", recode_mode="simple")
  
  # prepare output as data frame (unconditioned)
  unconditioned <- data.frame(
    glottocode = na.omit(original_variable_matrix[,c(1,which(names(original_variable_matrix)%in%.$`original.names`))])$glottocode, 
    value = new_data,  
    stringsAsFactors=FALSE)  
  
  # select condition statement
  condition_statement <-  .$condition.if.variable.conditioned
  
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

# merge new variables with recoded_data for further use
recoded_data <- full_join(eleventh_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that multiple conditioning is currently only implemented for up to 4 conditions !!
########## recode group 12 (simple recode - multiple conditioning) ########## 
# subset to features that have to be recoded, but only if several conditions apply
twelfth_set <- filter(recode_patterns, recode.operation.type=="recode group 12 (simple recode - multiple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 12 (simple recode - multiple conditioning)")

# recode if several conditions apply
twelfth_set_rec <- rowwise(twelfth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variable is present in the original variable matrix
  expect_true(.$`original.names` %in% names(original_variable_matrix)[-1])
  original_data <- na.omit(as.character(original_variable_matrix[[.$`original.names`]]))
  
  # extract relevant attributes for recoding
  expected_levels <- unlist(strsplit(unlist(.$original.levels), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_levels <- unlist(strsplit(.$new.levels, ";"))
  
  # recode (single variable)
  new_data <- implement_recode(original_data, expected_levels, recoding_groups, recoded_levels, 
                               nvar="single", recode_mode="simple")
  
  # prepare output as data frame (unconditioned)
  unconditioned <- data.frame(
    glottocode = na.omit(original_variable_matrix[,c(1,which(names(original_variable_matrix)%in%.$`original.names`))])$glottocode, 
    value = new_data,  
    stringsAsFactors=FALSE)  
  
  # select conditions
  conditions <- unlist(strsplit(.$condition.if.variable.conditioned," & "))
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

# merge new variables with recoded_data for further use
recoded_data <- full_join(twelfth_set_rec, recoded_data, by=c(glottocode="glottocode"))

########## recode group 13 (simple conditioning [conditioned variable]) ########## 
# subset to features that have to be conditioned on another conditioned feature
thirteenth_set <- filter(recode_patterns, recode.operation.type=="recode group 13 (simple conditioning [conditioned variable])")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 13 (simple conditioning [conditioned variable])")

# condition on conditioned variable
thirteenth_set_rec <- rowwise(thirteenth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variable is present in the original data
  expect_true(.$`original.names` %in% names(original_variable_matrix)[-1])
  variable_to_be_conditioned <- original_variable_matrix[,c(1,which(names(original_variable_matrix)%in%.$`original.names`))]
  
  # select condition statement
  condition_statement <-  .$condition.if.variable.conditioned
  
  # extract condition and equator
  condition_and_equator <- extract_condition_and_equator(condition_statement)
  
  # apply condition
  conditioned_data <- implement_conditioning(variable_to_be_conditioned, condition=condition_and_equator[[1]], equator=condition_and_equator[[2]])
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = conditioned_data$glottocode, 
    value = as.character(conditioned_data[,2]),  
    stringsAsFactors=FALSE)  
  
}) %>%
  spread(feature, value)

# merge new variables with recoded_data for further use
recoded_data <- full_join(thirteenth_set_rec, recoded_data, by=c(glottocode="glottocode"))

########## recode group 14 (simple recode - simple conditioning [conditioned variable]) ########## 
# subset to features that have to be recoded, but only if a condition applies involving conditioned variable
fourteenth_set <- filter(recode_patterns, recode.operation.type=="recode group 14 (simple recode - simple conditioning [conditioned variable])")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 14 (simple recode - simple conditioning [conditioned variable])")

# recode if condition on conditioned variable applies
fourteenth_set_rec <- rowwise(fourteenth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variable is present in the original variable matrix
  expect_true(.$`original.names` %in% names(original_variable_matrix)[-1])
  original_data <- na.omit(as.character(original_variable_matrix[[.$`original.names`]]))
  
  # extract relevant attributes for recoding
  expected_levels <- unlist(strsplit(unlist(.$original.levels), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_levels <- unlist(strsplit(.$new.levels, ";"))
  
  # recode (single variable)
  new_data <- implement_recode(original_data, expected_levels, recoding_groups, recoded_levels, 
                               nvar="single", recode_mode="simple")
  
  # prepare output as data frame (unconditioned)
  unconditioned <- data.frame(
    glottocode = na.omit(original_variable_matrix[,c(1,which(names(original_variable_matrix)%in%.$`original.names`))])$glottocode, 
    value = new_data,  
    stringsAsFactors=FALSE)  
  
  # select condition statement
  condition_statement <-  .$condition.if.variable.conditioned
  
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

# merge new variables with recoded_data for further use
recoded_data <- full_join(fourteenth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that multiple conditioning is currently only implemented for up to 4 conditions !!
########## recode group 15 (multiple conditioning [conditioned variable]) ########## 
# subset to features that have to be conditioned on several features, including a condition with a conditioned variable
fifteenth_set <- filter(recode_patterns, recode.operation.type=="recode group 15 (multiple conditioning [conditioned variable])")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 15 (multiple conditioning [conditioned variable])")

# condition on several variables
fifteenth_set_rec <- rowwise(fifteenth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variable is present in the original data
  expect_true(.$`original.names` %in% names(original_variable_matrix)[-1])
  variable_to_be_conditioned <- original_variable_matrix[,c(1,which(names(original_variable_matrix)%in%.$`original.names`))]
  
  # select conditions
  conditions <- unlist(strsplit(.$condition.if.variable.conditioned," & "))
  nr_conditions <- length(conditions)
  
  # determine and extract conditions and equators for condition 1, apply condition 1
  condition_statement_1 <- conditions[1]
  condition_and_equator_1 <- extract_condition_and_equator(condition_statement_1)
  conditioned_data <- implement_conditioning(variable_to_be_conditioned, 
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

# merge new variables with recoded_data for further use
recoded_data <- full_join(fifteenth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that multiple conditioning is currently only implemented for up to 4 conditions !!
########## recode group 16 (simple recode - multiple conditioning [conditioned variable]) ########## 
# subset to features that have to be conditioned on several features, including multiple conditions, including with conditioned variables
sixteenth_set <- filter(recode_patterns, recode.operation.type=="recode group 16 (simple recode - multiple conditioning [conditioned variable])")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 16 (simple recode - multiple conditioning [conditioned variable])")

# condition on several variables
sixteenth_set_rec <- rowwise(sixteenth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original variable is present in the original variable matrix
  expect_true(.$`original.names` %in% names(original_variable_matrix)[-1])
  original_data <- na.omit(as.character(original_variable_matrix[[.$`original.names`]]))
  
  # extract relevant attributes for recoding
  expected_levels <- unlist(strsplit(unlist(.$original.levels), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_levels <- unlist(strsplit(.$new.levels, ";"))
  
  # recode (single variable)
  new_data <- implement_recode(original_data, expected_levels, recoding_groups, recoded_levels, 
                               nvar="single", recode_mode="simple")
  
  # prepare output as data frame (unconditioned)
  unconditioned <- data.frame(
    glottocode = na.omit(original_variable_matrix[,c(1,which(names(original_variable_matrix)%in%.$`original.names`))])$glottocode, 
    value = new_data,  
    stringsAsFactors=FALSE)  
  
  # select conditions
  conditions <- unlist(strsplit(.$condition.if.variable.conditioned," & "))
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

# merge new variables with recoded_data for further use
recoded_data <- full_join(sixteenth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that multiple conditioning is currently only implemented for up to 4 conditions !!
########## save data ########## 
# first replace NA as explicit "NA"
recoded_data[is.na(recoded_data)]<-"NA"

# subset full data into original layer; logical/design layer and logical/design/statistical layer
recode_patterns <- read.csv("input/recode-patterns-variables.csv")
all_decisions <- read.csv("input/recode-patterns-decisions-log.csv")
logical_decisions <- all_decisions %>% filter(modification.type %in% c("logical","design-automated","design-manual"))
statistical_decisions <- all_decisions %>% filter(modification.type == "statistical")

original_layer <- recode_patterns%>%filter(original.variables=="TRUE")
logical_layer <- recode_patterns%>%filter(design.logical=="TRUE")
statistical_layer <- recode_patterns%>%filter(design.logical.statistical=="TRUE")

# save logical/design and logical/design/statistical data
logical_data <- recoded_data %>% select(c("glottocode",logical_layer$new.name))
logical_lg_filter <- logical_data
logical_lg_filter[logical_lg_filter=="?"]<-NA
logical_lg_filter[logical_lg_filter=="NA"]<-NA
logical_lg_filter$sum.non.na <- apply(logical_lg_filter,1,function(x)length(na.omit(x))-1)
logical_lg_filter <- logical_lg_filter %>% filter(sum.non.na > 0)
logical_data <- logical_data %>% filter(glottocode %in% logical_lg_filter$glottocode)
write.csv(logical_data,"output/logicalMM/full/logicalMM_full.csv")

statistical_data <- recoded_data %>% select(c("glottocode",statistical_layer$new.name))
statistical_lg_filter <- statistical_data
statistical_lg_filter[statistical_lg_filter=="?"]<-NA
statistical_lg_filter[statistical_lg_filter=="NA"]<-NA
statistical_lg_filter$sum.non.na <- apply(statistical_lg_filter,1,function(x)length(na.omit(x))-1)
statistical_lg_filter <- statistical_lg_filter %>% filter(sum.non.na > 0)
statistical_data <- statistical_data %>% filter(glottocode %in% statistical_lg_filter$glottocode)

write.csv(statistical_data,"output/statisticalMM/full/statisticalMM_full.csv")

########## sanity checks ########## 
# sanity checks for automated design decisions
input_nas <- as.data.frame(read_csv("output/merged_original_variables.csv", 
                                    trim_ws = FALSE, col_types = "f"))
input_nas[input_nas=="?"] <- NA
nrlevels<-data.frame(variable=names(input_nas),
                     number_of_lgs_coded=apply(input_nas,2,function(x)length(na.omit(x))),
                     number_of_variable_states=apply(input_nas,2,function(x)length(table(as.factor(x)))),
                     count_largest_variable_state=apply(input_nas,2,function(x)sort(table(as.factor(x)),decreasing=T)[1]),
                     count_second_largest_variable_state=apply(input_nas,2,function(x)sort(table(as.factor(x)),decreasing=T)[2]),
                     count_smallest_variable_state=apply(input_nas,2,function(x)sort(table(as.factor(x)),decreasing=F)[1])
                     )

# sanity check for "no variation"
no.variation.is <- filter(nrlevels,number_of_variable_states%in%c("0","1"))$variable
no.variation.should <- filter(logical_decisions,modification.ID == "A-IV0") %>% select(resulting.removed.variables) %>% str_split(pattern=", ") %>% unlist %>% str_split("o_") %>% unlist %>% unique
no.variation.should <- no.variation.should[no.variation.should!=""]
expect_true(all(no.variation.should%in%no.variation.is))
expect_true(all(no.variation.is%in%no.variation.should))

# sanity check too little variation // 1 
tl1.variation.is <- filter(nrlevels,count_second_largest_variable_state=="1")$variable
tl1.variation.is <- tl1.variation.is[tl1.variation.is!="glottocode"] # glottocode is of course fine
# --> WALS variables need to be able to be crossreferenced
tl1.variation.is[which(lapply(str_split(tl1.variation.is," "),length)>1)] <- unlist(lapply(strsplit(tl1.variation.is[which(lapply(str_split(tl1.variation.is," "),length)>1)]," "),function(x)(x[1])))
tl1.variation.is <- tl1.variation.is[tl1.variation.is!="glottocode"]
tl1.variation.should <- filter(logical_decisions,modification.ID == "A-IV1") %>% select(resulting.removed.variables) %>% str_split(pattern=", ") %>% unlist %>% str_split("o_") %>% unlist %>% unique
tl1.variation.should <- tl1.variation.should[tl1.variation.should!=""]
expect_true(all(tl1.variation.should%in%tl1.variation.is))
expect_true(all(tl1.variation.is%in%tl1.variation.should))

# sanity check too little variation // 2
tl2.variation.is <- filter(nrlevels,count_second_largest_variable_state=="2")$variable
# --> WALS variables need to be able to be crossreferenced
tl2.variation.is[which(lapply(str_split(tl2.variation.is," "),length)>1)] <- unlist(lapply(strsplit(tl2.variation.is[which(lapply(str_split(tl2.variation.is," "),length)>1)]," "),function(x)(x[1])))
tl2.variation.should <- filter(logical_decisions,modification.ID == "A-IV2") %>% select(resulting.removed.variables) %>% str_split(pattern=", ") %>% unlist %>% str_split("o_") %>% unlist %>% unique
tl2.variation.should <- tl2.variation.should[tl2.variation.should!=""]
expect_true(all(tl2.variation.should%in%tl2.variation.is))
expect_true(all(tl2.variation.is%in%tl2.variation.should))

# sanity check <100 lgs
remaining_lgs_nrlevels <- nrlevels %>% filter(count_second_largest_variable_state%in%c("0","1","2")==F) %>% filter(number_of_variable_states!="1")
too.few.lgs.is <- filter(remaining_lgs_nrlevels,number_of_lgs_coded < 100)$variable
# --> WALS variables need to be able to be crossreferenced
too.few.lgs.is[which(lapply(str_split(too.few.lgs.is," "),length)>1)] <- unlist(lapply(strsplit(too.few.lgs.is[which(lapply(str_split(too.few.lgs.is," "),length)>1)]," "),function(x)(x[1])))
too.few.lgs.should <- filter(logical_decisions,modification.ID == "A-U100LG") %>% select(resulting.removed.variables) %>% str_split(pattern=", ") %>% unlist %>% str_split("o_") %>% unlist %>% unique
too.few.lgs.should <- too.few.lgs.should[too.few.lgs.should!=""]
expect_true(all(too.few.lgs.should%in%too.few.lgs.is))
expect_true(all(too.few.lgs.is%in%too.few.lgs.should))

### check all original variables that should be in the original_variable filter are in there and vice versa
original_names_is <- original_layer$new.name
original_names_should <- names(retained_data)[-1]
expect_true(all(original_names_is%in%original_names_should))
expect_true(all(original_names_should%in%original_names_is))

### logical/design: check all original variables that should be in the logical/design filter are in there and vice versa
design_add <- c(unique(unlist(str_split(filter(logical_decisions,is.added.variable.kept == "yes")$resulting.added.variables,", "))),NA)
design_remove <- unique(unlist(str_split(logical_decisions$resulting.removed.variables,", ")))
# logical/design is: a) original variables, plus b) all design/logical additions, minus c) all design/logical removals
logical_names_should <- unique(c(original_names_is,design_add))[unique(c(original_names_is,design_add))%in%design_remove==F]
logical_names_is <- names(logical_data)[names(logical_data)!="glottocode"]
# sanity checks
expect_true(all(logical_names_is %in% logical_names_should))
expect_true(all(logical_names_should %in% logical_names_is))

### statistical+logical/design: check all original variables that should be in the logical/design/statistical filter are in there and vice versa
statistical_add <- unique(statistical_decisions$resulting.added.variables) 
statistical_remove <- unique(unlist(str_split(statistical_decisions$resulting.removed.variables,", ")))
# statistical is: a) the variables from the final design/logical, plus b) all statistical additions, minus c) all statistical removals
statistical_names_should <- unique(c(logical_names_is,statistical_add))[unique(c(logical_names_is,statistical_add))%in%statistical_remove==F]
statistical_names_is <- names(statistical_data)[names(statistical_data)!="glottocode"]
# sanity checks
expect_true(all(statistical_names_is %in% statistical_names_should))
expect_true(all(statistical_names_should %in% statistical_names_is))

### modification ID match --> ensure all modification IDs in the variables sheet are in the modification sheet and vice versa
mod_IDs_is <- na.omit(unique(c(unlist(str_split(recode_patterns$modification.IDs,";")),(unlist(str_split(recode_patterns$associated.modification.IDs.without.resulting.action,";"))))))
mod_IDs_is <- mod_IDs_is[mod_IDs_is!=""]
mod_IDs_should <- na.omit(all_decisions$modification.ID)
expect_true(all(mod_IDs_is %in% mod_IDs_should))
expect_true(all(mod_IDs_should %in% mod_IDs_is))

# specific modification ID match --> ensure that each modification ID in the variables sheet is in the modification sheet, associated via the correct columns and variables; and vice versa
ids <- na.omit(all_decisions$modification.ID)
for (id in ids){
  type <- filter(all_decisions,modification.ID == id)$modification.type
  if (type == "statistical"){
    # check that each instance of a modification ID in the spreadsheet ("is") is foreseen in the decisions_log ("should") and vice versa
    should_all <- all_decisions %>% filter(modification.ID == id) %>% select(c("variable.1.for.test","variable.2.for.test","resulting.added.variables","resulting.removed.variables")) %>% as.character() %>% unique()
    should_all <- na.omit(unique(unlist(strsplit(should_all[should_all!="NA"],", "))))
    is_all <- recode_patterns %>% slice(c(which(grepl(id,recode_patterns$modification.IDs)),which(grepl(id,recode_patterns$associated.modification.IDs.without.resulting.action))))
    is_all <- is_all$new.name
    expect_true(all(is_all%in%should_all))
    expect_true(all(should_all%in%is_all))
    
    # check that each instance of a modification ID WITH EFFECT in the spreadsheet ("is") is foreseen in the decisions_log ("should") and vice versa
    should_actedupon <- all_decisions %>% filter(modification.ID == id) %>% select(c("resulting.added.variables","resulting.removed.variables")) %>% as.character() %>% unique()
    should_actedupon <- na.omit(unique(unlist(strsplit(should_actedupon[should_actedupon!="NA"],", "))))
    is_actedupon <- recode_patterns %>% slice(which(grepl(id,recode_patterns$modification.IDs)))
    is_actedupon <- is_actedupon$new.name
    expect_true(all(is_actedupon%in%should_actedupon))
    expect_true(all(should_actedupon%in%is_actedupon))
    
    # check that each instance of a modification ID WITHOUT EFFECT in the spreadsheet ("is") is foreseen in the decisions_log ("should") and vice versa
    should_associated <- all_decisions %>% filter(modification.ID == id) %>% select(c("variable.1.for.test","variable.2.for.test")) %>% as.character() %>% unique()
    should_associated <- setdiff(na.omit(unique(unlist(strsplit(should_associated[should_associated!="NA"],", ")))),is_actedupon)
    is_associated <- recode_patterns %>% slice(which(grepl(id,recode_patterns$associated.modification.IDs.without.resulting.action)))
    is_associated <- is_associated$new.name
    expect_true(all(is_associated%in%should_associated))
    expect_true(all(should_associated%in%is_associated))
    print(c(id,"ok (3/3)"))
    
  }
  else if (type %in% c("logical","design-automated","design-manual")){ 
    # check that each instance of a modification ID in the spreadsheet ("is") is foreseen in the decisions_log ("should") and vice versa
    should_all <- all_decisions %>% filter(modification.ID == id) %>% select(c("relevant.variables","resulting.added.variables","resulting.removed.variables")) %>% as.character() %>% unique()
    should_all <- na.omit(unique(unlist(strsplit(should_all[should_all!="NA"],", "))))
    is_all <- recode_patterns %>% slice(c(which(grepl(id,recode_patterns$modification.IDs)),which(grepl(id,recode_patterns$associated.modification.IDs.without.resulting.action))))
    is_all <- is_all$new.name
    expect_true(all(is_all%in%should_all))
    expect_true(all(should_all%in%is_all))
    
    # check that each instance of a modification ID WITH EFFECT in the spreadsheet ("is") is foreseen in the decisions_log ("should") and vice versa
    should_actedupon <- all_decisions %>% filter(modification.ID == id) %>% select(c("resulting.added.variables","resulting.removed.variables")) %>% as.character() %>% unique()
    should_actedupon <- na.omit(unique(unlist(strsplit(should_actedupon[should_actedupon!="NA"],", "))))
    is_actedupon <- recode_patterns %>% slice(which(grepl(id,recode_patterns$modification.IDs)))
    is_actedupon <- is_actedupon$new.name
    expect_true(all(is_actedupon%in%should_actedupon))
    expect_true(all(should_actedupon%in%is_actedupon))
    
    # check that each instance of a modification ID WITHOUT EFFECT in the spreadsheet ("is") is foreseen in the decisions_log ("should") and vice versa
    should_associated <- all_decisions %>% filter(modification.ID == id) %>% select("relevant.variables") %>% as.character() %>% unique()
    should_associated <- setdiff(na.omit(unique(unlist(strsplit(should_associated[should_associated!="NA"],", ")))),is_actedupon)
    is_associated <- recode_patterns %>% slice(which(grepl(id,recode_patterns$associated.modification.IDs.without.resulting.action)))
    is_associated <- is_associated$new.name
    expect_true(all(is_associated%in%should_associated))
    expect_true(all(should_associated%in%is_associated))
    print(c(id,"ok (3/3)"))
  }
}

## known.remaining dependencies match - check remaining dependencies are appropriately logged (one-sided, because there are dependencies not associated with statistical tests!)
rds <- filter(all_decisions, resulting.modification == "tendency logged in known.remaining.dependencies")
for (rd in 1:nrow(rds)){
  id <- rds %>% slice(rd) %>% select(modification.ID)
  should_associated <- rds %>% slice(rd) %>% select(c(variable.1.for.test,variable.2.for.test)) %>% as.character
  is_associated <- recode_patterns %>% slice(which(grepl(id,recode_patterns$known.remaining.dependencies.after.statistical.treatment))) %>% select(new.name) %>% unlist %>% as.character
  expect_true(all(is_associated%in%should_associated))
  expect_true(all(should_associated%in%is_associated))
  print(c(unlist(id),"ok"))
}
expect_true(all(rds$modification.ID%in%unlist(strsplit(recode_patterns$known.remaining.dependencies.after.statistical.treatment,";"))))
expect_true(all(na.omit(unique(unlist(strsplit(recode_patterns$known.remaining.dependencies.after.statistical.treatment,";"))))%in%rds$modification.ID))


## untestable dependencies match - check remaining dependencies are appropriately logged
ntrd <- filter(all_decisions, resulting.modification == "modification.ID logged in not.testable.remaining.dependencies")
for (rd in 1:nrow(ntrd)){
  id <- ntrd %>% slice(rd) %>% select(modification.ID)
  should_associated <- ntrd %>% slice(rd) %>% select(c(variable.1.for.test,variable.2.for.test)) %>% as.character
  is_associated <- recode_patterns %>% slice(which(grepl(id,recode_patterns$suspected.but.untestable.dependencies))) %>% select(new.name) %>% unlist %>% as.character
  expect_true(all(is_associated%in%should_associated))
  expect_true(all(should_associated%in%is_associated))
  print(c(unlist(id),"ok"))
}
expect_true(all(ntrd$modification.ID%in%unlist(strsplit(recode_patterns$suspected.but.untestable.dependencies,";"))))
expect_true(all(na.omit(unique(unlist(strsplit(recode_patterns$suspected.but.untestable.dependencies,";"))))%in%ntrd$modification.ID))

########## statistical dependencies ######### 
# create 1000 diversity samples; each row of diversity_samples is a diversity sample of 120 lgs; 1 from each family (glottolog-based) 
macroareas <- read.csv("../input/glottolog_v.4.8/languages_and_dialects_geo.csv")
names(macroareas)[1]<-"id"
taxonomy <- merge(taxonomy,select(macroareas, c("id","macroarea")),by="id")

# to create the samples, take into consideration all glottocodes for which at least 10 features are coded (not ? or NA) in the design.logical matrix
data_for_sample_generation <- logical_data
data_for_sample_generation[data_for_sample_generation=="?"]<-NA
data_for_sample_generation[data_for_sample_generation=="NA"]<-NA
data_for_sample_generation$sum.non.na <- apply(data_for_sample_generation,1,function(x)length(na.omit(x))-1)
data_for_sample_generation <- data_for_sample_generation %>% filter(sum.non.na > 9)

# select these glottocodes; subset taxonomy to relevant glottocodes
recoded_lgs <- as.character(data_for_sample_generation$glottocode)
recoded_lgs <- filter(taxonomy,id%in%recoded_lgs)

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
    reg_fam <- unique(filter(recoded_lgs,macroarea==regions[rg,])$level1) # these are the families in this macroarea
    fam_add <- data.frame(family=reg_fam[sample(length(reg_fam),size=20,replace=F)]) # sample 20 different random families in the macroarea
    diversity_samples_fams<-rbind(diversity_samples_fams,fam_add) # log the 20 families for this sample
    diversity_samples_lgs<-rbind(diversity_samples_lgs, # sample and save one random language per family
                                 data.frame(lg=(apply(fam_add,1,function(x)
                                   filter(recoded_lgs,level1==x[1]&macroarea==regions[rg,])[sample(nrow(filter(recoded_lgs,level1==x[1]&macroarea==regions[rg,])),size=1,replace=F),]$id))))
  }
  lg_samples_1000[,smpl]<-diversity_samples_lgs[,"lg"] # at the end of each iteration, store the full language sample
}

# store the diversity samples
write.csv(lg_samples_1000,"output/1000_diversity_samples_seed_2023.csv")

# specify test condition, fix threshold
proportion_languages_must_be_in_applicable_state <- 1/3 # for a sample to be considered for evaluation, at least 1/3 of languages need to be coded for the second variable given the relevant state(s) of the first

# load expectations
expectations <- read.csv("input/recode-patterns-decisions-log.csv") %>%
  filter(modification.type == "statistical")

# run all tests
results <- evaluate_XOR_AND_THEN(recoded_data = recoded_data, 
                                 expectations = expectations, 
                                 taxonomy = taxonomy, 
                                 diversity_samples = lg_samples_1000,
                                 proportion_languages_must_be_in_applicable_state = proportion_languages_must_be_in_applicable_state)

results <- na_convert(results)
results$maximum.sum <- apply(results,1,function(x)if(is.na(x[13])){x[10]}else if(x[10]>x[13]){x[10]}else{x[13]})
results$minimum.cohensd <- apply(results,1,function(x)if(is.na(x[17])){x[16]}else if(x[16]>x[17]){x[17]}else{x[16]})
results$minimum.sum <- apply(results,1,function(x)if(is.na(x[13])){x[10]}else if(x[10]>x[13]){x[13]}else{x[10]})
results$maximum.cohensd <- apply(results,1,function(x)if(is.na(x[17])){x[16]}else if(x[16]>x[17]){x[16]}else{x[17]})

# check the results logged in the file are correct
expect_true(all(round(resultstbc,digits=7) == round(expectations[,26:46],digits=7), na.rm = T))

########## make cldf ########## 
# languages.csv
taxonomy_for_csv <- as_flat_taxonomy_matrix(glottolog_languoids)
names(taxonomy_for_csv)[1] <- "glottocode"
macroareas <- read.csv("../input/glottolog_v.4.8/languages_and_dialects_geo.csv")

taxonomy_logical <- data.frame(glottocode = logical_data$glottocode)
taxonomy_logical <- left_join(taxonomy_logical,macroareas)
taxonomy_logical <- left_join(taxonomy_logical,taxonomy_for_csv)
write.csv(taxonomy_logical,"output/logicalMM/cldf/languages.csv",fileEncoding="UTF-8",row.names = F)

taxonomy_statistical <- data.frame(glottocode = statistical_data$glottocode)
taxonomy_statistical <- left_join(taxonomy_statistical,macroareas)
taxonomy_statistical <- left_join(taxonomy_statistical,taxonomy_for_csv)
write.csv(taxonomy_statistical,"output/statisticalMM/cldf/languages.csv",fileEncoding="UTF-8",row.names = F)

# parameters.csv is directly from google sheets
parameters <- read.csv("input/recode-patterns-variables.csv") %>% select(1:23)

parameters_logical <- parameters %>% filter(design.logical==T)
parameters_statistical <- parameters %>% filter(design.logical.statistical==T)

write.csv(parameters_logical,"output/logicalMM/cldf/parameters.csv",fileEncoding="UTF-8",row.names = F)
write.csv(parameters_statistical,"output/statisticalMM/cldf/parameters.csv",fileEncoding="UTF-8",row.names = F)

# values.csv
library(reshape2)
library(data.table)
logical_long <- melt(setDT(logical_data), id.vars = "glottocode", variable.name = "new.name")
logical_long$value_ID <- apply(logical_long,1,function(x) paste(x[2],x[1],sep="-"))
logical_long$code_ID <- apply(logical_long,1,function(x) paste(x[2],x[3],sep="-"))
logical_long <- logical_long %>% select(c("value_ID","glottocode","new.name","value","code_ID"))
logical_long$glottocode <- as.character(logical_long$glottocode)
logical_long$new.name <- as.character(logical_long$new.name)
write.csv(logical_long,"output/logicalMM/cldf/values.csv",fileEncoding="UTF-8",row.names = F)

statistical_long <- melt(setDT(statistical_data), id.vars = "glottocode", variable.name = "new.name")
statistical_long$value_ID <- apply(statistical_long,1,function(x) paste(x[2],x[1],sep="-"))
statistical_long$code_ID <- apply(statistical_long,1,function(x) paste(x[2],x[3],sep="-"))
statistical_long$glottocode <- as.character(statistical_long$glottocode)
statistical_long$new.name <- as.character(statistical_long$new.name)
statistical_long <- statistical_long %>% select(c("value_ID","glottocode","new.name","value","code_ID"))
write.csv(statistical_long,"output/statisticalMM/cldf/values.csv",fileEncoding="UTF-8",row.names = F)

# codes.csv
logical_codes <- logical_long %>% select(c("code_ID","new.name","value")) %>% unique()
write.csv(logical_codes,"output/logicalMM/cldf/codes.csv",fileEncoding="UTF-8",row.names = F)

statistical_codes <- statistical_long %>% select(c("code_ID","new.name","value")) %>% unique()
write.csv(statistical_codes,"output/statisticalMM/cldf/codes.csv",fileEncoding="UTF-8",row.names = F)

# modifications.csv is directly from google sheets (columns B-AI)
modifications <- read.csv("input/recode-patterns-decisions-log.csv")

write.csv(modifications,"output/logicalMM/cldf/modifications.csv",fileEncoding="UTF-8",row.names = F)
write.csv(modifications,"output/statisticalMM/cldf/modifications.csv",fileEncoding="UTF-8",row.names = F)

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
