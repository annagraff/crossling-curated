#### This script reads the external input data, parses it via all documented curation steps and generates and saves the 
#### TypLinkInd logical and statistical datasets. Before saving, logical quality checks are performed on the outputs.

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
# read in original features and taxonomy
original_feature_matrix <- as.data.frame(read_csv("../../curated_data/TypLinkInd/compiled_external_input_features.csv", 
                                                   trim_ws = FALSE, col_types = "f"))

taxonomy <- as_flat_taxonomy_matrix(glottolog_languoids)

# some glottocodes are not found in glottolog v. 5.0. because they have been assigned new glottocodes - replace them here (manual assignment via ISO 639-3 code)
original_feature_matrix$glottocode[which(original_feature_matrix$glottocode %in% taxonomy$id == F)]

original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="woro1255", new="worr1237") # Worora from WALS
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="jogl1236", new="tase1235") # Jugli from WALS
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="sanm1259", new="sanm1295") # Mixtec (Molinos) from WALS
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="cypr1245", new="mode1248") # Greek (Cypriot) from WALS
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="yidd1255", new="east2295") # Yiddish from WALS
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="lenc1244", new="lenc1239") # Lenca-Salvador (Bookkeeping) from Lexibank

expect_true(all(original_feature_matrix$glottocode %in% taxonomy$id))

# some glottocodes are of the kind "bookkeeping" and have been retired / merged with other glottocodes. we update these glottocodes manually
taxonomy %>% filter(id %in% original_feature_matrix$glottocode) %>% filter(level1=="book1242") %>% select(id)
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="band1337", new="darl1243") # Bandjigali from PHOIBLE
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="bubi1249", new="bube1242") # Bubia from PHOIBLE
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="chua1256", new="firs1234") # Chuanqiandian Cluster Miao from PHOIBLE
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="kukn1238", new="dhod1238") # Kukna from WALS 
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="naxi1246", new="yong1270") # Naxi (Yongning) from Lexibank --> Narua is best match (Yongning)
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="nucl1668", new="kana1291") # Katukina from PHOIBLE
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="samr1245", new="cent2314") # Samre from PHOIBLE --> reassigned glottocode via reference
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="sana1281", new="sana1298") # Sanapana from Lexibank
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="sout3125", new="sout2920") # Betsimisaraka from Lexibank
original_feature_matrix <- update_glottocode_and_data(data=original_feature_matrix, old="wela1234", new="ngal1291") # Rawngtu from Lexibank

expect_false(any(select(filter(taxonomy, id %in% original_feature_matrix$glottocode), level1) == "book1242"))

# read in the file specifying the maintained features, the and the recoded features with their recoding patterns
recode_patterns <- read.csv("feature-recode-patterns.csv")

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

## !! note that feature priorisation is currently implemented for up to 5 features !!
## recode group 2 (merge features - simple recode) ##
# subset to features that are merged
second_set <- filter(recode_patterns, recode.operation.type=="recode group 2 (merge features - simple recode)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 2 (merge features - simple recode)")

# merge and recode features
second_set <- rowwise(second_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original features are present in the original feature matrix
  expect_true(unique(unlist(strsplit(.$`original.names`,"&")) %in% names(original_feature_matrix)[-1]))
  original_data <- original_feature_matrix[,c("glottocode", unlist(strsplit(.$`original.names`,"&")))]
  
  # define starting feature states according to merging priorities
  original_data$merged <- ifelse((original_data[,.$first.priority.feature.if.merged] != "?"), 
                                 as.character(original_data[,.$first.priority.feature.if.merged]), 
                                 as.character(original_data[,.$second.priority.feature.if.merged]))
  
  if(ncol(original_data)==5){ # if there is a third feature
    original_data$merged <- ifelse((original_data$merged != "?"), 
                                   as.character(original_data$merged),
                                   as.character(original_data[,.$third.priority.feature.if.merged]))
  }
  
  if(ncol(original_data)==6){ # if there is a fourth feature
    original_data$merged <- ifelse((original_data$merged != "?"), 
                                   as.character(original_data$merged),
                                   as.character(original_data[,.$fourth.priority.feature.if.merged]))
  }
  
  original_data <- na.omit(select(original_data,c(glottocode,merged)))
  
  # extract relevant attributes for recoding
  expected_states <- unlist(strsplit(unlist(.$original.states), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_states <- unlist(strsplit(.$new.states, ";"))
  
  # recode (multiple features)
  new_data <- implement_recode(original_data, expected_states, recoding_groups, recoded_states, 
                               nvar="multiple", recode_mode="simple")
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = original_data$glottocode, 
    value = new_data,  
    stringsAsFactors=FALSE)  
  
}) %>%
  spread(feature, value)

# merge new features with recoded_data for further use
recoded_data <- full_join(second_set, recoded_data, by=c(glottocode="glottocode"))

## recode group 3 (merge features - recode via logical arguments) ##
# subset to features that have to be merged via logical arguments
third_set <- filter(recode_patterns, recode.operation.type=="recode group 3 (merge features - recode via logical arguments)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 3 (merge features - recode via logical arguments)")

# merge and recode features via logical arguments
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
recoded_data <- full_join(third_set_rec, recoded_data, by=c(glottocode="glottocode"))

## recode group 4 (recode via logical arguments) ##
# subset to features that have to be merged via logical arguments
fourth_set <- filter(recode_patterns, recode.operation.type=="recode group 4 (recode via logical arguments)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 4 (recode via logical arguments)")

# recode features via logical arguments
fourth_set_rec <- rowwise(fourth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original features are present in the original feature matrix
  original_data<-original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%unlist(strsplit(.$`original.names`,"&"))))]
  
  # extract relevant attributes for recoding
  expected_states <- unlist(strsplit(unlist(.$recode.operation.Rcode), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_states <- unlist(strsplit(.$new.states, ";"))
  
  # recode (multiple features)
  new_data <- implement_recode(original_data, expected_states, recoding_groups, recoded_states, nvar="single", recode_mode="logical_arguments")
  
  # prepare output as data frame
  data.frame(
    feature = .$new.name, 
    glottocode = original_data$glottocode, 
    value = new_data,  
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

## recode group 6 (merge features - recode via logical arguments - simple conditioning) ## 
# subset to features that have to be recoded via logical arguments if a condition applies
sixth_set <- filter(recode_patterns, recode.operation.type=="recode group 6 (merge features - recode via logical arguments - simple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 6 (merge features - recode via logical arguments - simple conditioning)")

# merge and recode via logical arguments if a condition applies
sixth_set_rec <- rowwise(sixth_set) %>% do({
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
recoded_data <- full_join(sixth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that multiple conditioning is currently only implemented for up to 4 conditions !!
## recode group 7 (multiple conditioning) ##
# subset to features that have to be conditioned on several features
seventh_set <- filter(recode_patterns, recode.operation.type=="recode group 7 (multiple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 7 (multiple conditioning)")

# condition on several features
seventh_set_rec <- rowwise(seventh_set) %>% do({
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
recoded_data <- full_join(seventh_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that feature priorisation is currently implemented for up to 5 features !!
## recode group 8 (merge features - simple recode - simple conditioning) ##
# subset to features that have to be created by merging and recoding, but only if a condition applies
eighth_set <- filter(recode_patterns, recode.operation.type=="recode group 8 (merge features - simple recode - simple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 8 (merge features - simple recode - simple conditioning)")

# merge and recode if condition applies
eighth_set_rec <- rowwise(eighth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original features are present in the original feature matrix
  expect_true(unique(unlist(strsplit(.$`original.names`,"&")) %in% names(original_feature_matrix)[-1]))
  original_data <- original_feature_matrix[,c("glottocode", unlist(strsplit(.$`original.names`,"&")))]
  
  # define starting feature states according to merging priorities
  original_data$merged <- ifelse((original_data[,.$first.priority.feature.if.merged] != "?"), 
                                 as.character(original_data[,.$first.priority.feature.if.merged]), 
                                 as.character(original_data[,.$second.priority.feature.if.merged]))
  
  if(ncol(original_data)==5){ # if there is a third feature
    original_data$merged <- ifelse((original_data$merged != "?"), 
                                   as.character(original_data$merged),
                                   as.character(original_data[,.$third.priority.feature.if.merged]))
  }
  
  if(ncol(original_data)==6){ # if there is a fourth feature
    original_data$merged <- ifelse((original_data$merged != "?"), 
                                   as.character(original_data$merged),
                                   as.character(original_data[,.$fourth.priority.feature.if.merged]))
  }
  
  original_data <- na.omit(select(original_data,c(glottocode,merged)))
  
  # extract relevant attributes for recoding
  expected_states <- unlist(strsplit(unlist(.$original.states), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_states <- unlist(strsplit(.$new.states, ";"))
  
  # recode (multiple features)
  new_data <- implement_recode(original_data, expected_states, recoding_groups, recoded_states, 
                               nvar="multiple", recode_mode="simple")
  
  # prepare output as data frame (unconditioned)
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
recoded_data <- full_join(eighth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that feature priorisation is currently implemented for up to 5 features !!
## !! note that multiple conditioning is currently only implemented for up to 4 conditions !!

## recode group 9 (merge features - simple recode - multiple conditioning) ##
# subset to features that have to be merged and recoded, but only if several conditions apply
ninth_set <- filter(recode_patterns, recode.operation.type=="recode group 9 (merge features - simple recode - multiple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 9 (merge features - simple recode - multiple conditioning)")

# merge and recode if several conditions apply
ninth_set_rec <- rowwise(ninth_set) %>% do({
  cat("Processing ", .$new.name, "\n", sep="")
  
  # check that the original features are present in the original feature matrix
  expect_true(unique(unlist(strsplit(.$`original.names`,"&")) %in% names(original_feature_matrix)[-1]))
  original_data <- original_feature_matrix[,c("glottocode", unlist(strsplit(.$`original.names`,"&")))]
  
  # define starting feature states according to merging priorities
  original_data$merged <- ifelse((original_data[,.$first.priority.feature.if.merged] != "?"), 
                                 as.character(original_data[,.$first.priority.feature.if.merged]), 
                                 as.character(original_data[,.$second.priority.feature.if.merged]))
  
  if(ncol(original_data)==5){ # if there is a third feature
    original_data$merged <- ifelse((original_data$merged != "?"), 
                                   as.character(original_data$merged),
                                   as.character(original_data[,.$third.priority.feature.if.merged]))
  }
  
  if(ncol(original_data)==6){ # if there is a fourth feature
    original_data$merged <- ifelse((original_data$merged != "?"), 
                                   as.character(original_data$merged),
                                   as.character(original_data[,.$fourth.priority.feature.if.merged]))
  }
  
  original_data <- na.omit(select(original_data,c(glottocode,merged)))
  
  # extract relevant attributes for recoding
  expected_states <- unlist(strsplit(unlist(.$original.states), ";"))
  recoding_groups <- unlist(strsplit(.$recode.pattern, "-"))
  recoded_states <- unlist(strsplit(.$new.states, ";"))
  
  # recode (multiple features)
  new_data <- implement_recode(original_data, expected_states, recoding_groups, recoded_states, 
                               nvar="multiple", recode_mode="simple")
  
  # prepare output as data frame (unconditioned)
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
recoded_data <- full_join(ninth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that multiple conditioning is currently only implemented for up to 4 conditions !!

## recode group 10 (merge features - recode via logical arguments - multiple conditioning) ##
# subset to features that have to be recoded via logical arguments if several conditions apply
tenth_set <- filter(recode_patterns, recode.operation.type=="recode group 10 (merge features - recode via logical arguments - multiple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 10 (merge features - recode via logical arguments - multiple conditioning)")

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

## recode group 11 (simple recode - simple conditioning) ##
# subset to features that have to be recoded, but only if a condition applies
eleventh_set <- filter(recode_patterns, recode.operation.type=="recode group 11 (simple recode - simple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 11 (simple recode - simple conditioning)")

# recode if condition applies
eleventh_set_rec <- rowwise(eleventh_set) %>% do({
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
  
  # prepare output as data frame (unconditioned)
  unconditioned <- data.frame(
    glottocode = na.omit(original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%.$`original.names`))])$glottocode, 
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
recoded_data <- full_join(eleventh_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that multiple conditioning is currently only implemented for up to 4 conditions !!
## recode group 12 (simple recode - multiple conditioning) ##
# subset to features that have to be recoded, but only if several conditions apply
twelfth_set <- filter(recode_patterns, recode.operation.type=="recode group 12 (simple recode - multiple conditioning)")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 12 (simple recode - multiple conditioning)")

# recode if several conditions apply
twelfth_set_rec <- rowwise(twelfth_set) %>% do({
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
  
  # prepare output as data frame (unconditioned)
  unconditioned <- data.frame(
    glottocode = na.omit(original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%.$`original.names`))])$glottocode, 
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
recoded_data <- full_join(twelfth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## recode group 13 (simple conditioning [conditioned feature]) ##
# subset to features that have to be conditioned on another conditioned feature
thirteenth_set <- filter(recode_patterns, recode.operation.type=="recode group 13 (simple conditioning [conditioned feature])")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 13 (simple conditioning [conditioned feature])")

# condition on conditioned feature
thirteenth_set_rec <- rowwise(thirteenth_set) %>% do({
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
recoded_data <- full_join(thirteenth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## recode group 14 (simple recode - simple conditioning [conditioned feature]) ##
# subset to features that have to be recoded, but only if a condition applies involving conditioned feature
fourteenth_set <- filter(recode_patterns, recode.operation.type=="recode group 14 (simple recode - simple conditioning [conditioned feature])")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 14 (simple recode - simple conditioning [conditioned feature])")

# recode if condition on conditioned feature applies
fourteenth_set_rec <- rowwise(fourteenth_set) %>% do({
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
  
  # prepare output as data frame (unconditioned)
  unconditioned <- data.frame(
    glottocode = na.omit(original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%.$`original.names`))])$glottocode, 
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
recoded_data <- full_join(fourteenth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that multiple conditioning is currently only implemented for up to 4 conditions !!
## recode group 15 (multiple conditioning [conditioned feature]) ##
# subset to features that have to be conditioned on several features, including a condition with a conditioned feature
fifteenth_set <- filter(recode_patterns, recode.operation.type=="recode group 15 (multiple conditioning [conditioned feature])")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 15 (multiple conditioning [conditioned feature])")

# condition on several features
fifteenth_set_rec <- rowwise(fifteenth_set) %>% do({
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
recoded_data <- full_join(fifteenth_set_rec, recoded_data, by=c(glottocode="glottocode"))

## !! note that multiple conditioning is currently only implemented for up to 4 conditions !!
## recode group 16 (simple recode - multiple conditioning [conditioned feature]) ##
# subset to features that have to be conditioned on several features, including multiple conditions, including with conditioned features
sixteenth_set <- filter(recode_patterns, recode.operation.type=="recode group 16 (simple recode - multiple conditioning [conditioned feature])")
recode_patterns <- filter(recode_patterns, recode.operation.type!="recode group 16 (simple recode - multiple conditioning [conditioned feature])")

# condition on several features
sixteenth_set_rec <- rowwise(sixteenth_set) %>% do({
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
  
  # prepare output as data frame (unconditioned)
  unconditioned <- data.frame(
    glottocode = na.omit(original_feature_matrix[,c(1,which(names(original_feature_matrix)%in%.$`original.names`))])$glottocode, 
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
recoded_data <- full_join(sixteenth_set_rec, recoded_data, by=c(glottocode="glottocode"))

# replace NA as explicit "NA"
recoded_data[is.na(recoded_data)]<-"NA"

# this full set of all input and recoded features needs to be stored to perform statistical tests
write.csv(recoded_data, "../../curated_data/TypLinkInd/statisticalTLI/data_for_stats.csv")

# subset full data into original layer; logical layer and statistical layer
recode_patterns <- read.csv("feature-recode-patterns.csv")
all_decisions <- read.csv("decisions-log.csv")
logical_decisions <- all_decisions %>% filter(modification.type %in% c("logical","design-automated","design-manual"))
statistical_decisions <- all_decisions %>% filter(modification.type == "statistical")

original_layer <- recode_patterns%>%filter(original.features=="TRUE")
logical_layer <- recode_patterns%>%filter(design.logical=="TRUE")
statistical_layer <- recode_patterns%>%filter(design.logical.statistical=="TRUE")

logical_data <- recoded_data %>% select(c("glottocode",logical_layer$new.name))
logical_lg_filter <- logical_data
logical_lg_filter[logical_lg_filter=="?"]<-NA
logical_lg_filter[logical_lg_filter=="NA"]<-NA
logical_lg_filter$sum.non.na <- apply(logical_lg_filter,1,function(x)length(na.omit(x))-1)
logical_lg_filter <- logical_lg_filter %>% filter(sum.non.na > 0)
logical_data <- logical_data %>% filter(glottocode %in% logical_lg_filter$glottocode)

statistical_data <- recoded_data %>% select(c("glottocode",statistical_layer$new.name))
statistical_lg_filter <- statistical_data
statistical_lg_filter[statistical_lg_filter=="?"]<-NA
statistical_lg_filter[statistical_lg_filter=="NA"]<-NA
statistical_lg_filter$sum.non.na <- apply(statistical_lg_filter,1,function(x)length(na.omit(x))-1)
statistical_lg_filter <- statistical_lg_filter %>% filter(sum.non.na > 0)
statistical_data <- statistical_data %>% filter(glottocode %in% statistical_lg_filter$glottocode)


########## sanity checks ########## 
# sanity checks for automated design decisions
input_nas <- as.data.frame(read_csv("../../curated_data/TypLinkInd/compiled_external_input_features.csv", 
                                    trim_ws = FALSE, col_types = "f"))
input_nas[input_nas=="?"] <- NA
nrstates<-data.frame(feature=names(input_nas),
                     number_of_lgs_coded=apply(input_nas,2,function(x)length(na.omit(x))),
                     number_of_feature_states=apply(input_nas,2,function(x)length(table(as.factor(x)))),
                     count_largest_feature_state=apply(input_nas,2,function(x)sort(table(as.factor(x)),decreasing=T)[1]),
                     count_second_largest_feature_state=apply(input_nas,2,function(x)sort(table(as.factor(x)),decreasing=T)[2]),
                     count_smallest_feature_state=apply(input_nas,2,function(x)sort(table(as.factor(x)),decreasing=F)[1])
)

# sanity check for "no variation"
no.variation.is <- filter(nrstates,number_of_feature_states%in%c("0","1"))$feature
no.variation.should <- filter(logical_decisions,modification.ID == "A-IV0") %>% select(resulting.removed.features) %>% str_split(pattern=", ") %>% unlist %>% str_split("o_") %>% unlist %>% unique
no.variation.should <- no.variation.should[no.variation.should!=""]
expect_true(all(no.variation.should%in%no.variation.is))
expect_true(all(no.variation.is%in%no.variation.should))

# sanity check too little variation // 1 
tl1.variation.is <- filter(nrstates,count_second_largest_feature_state=="1")$feature
tl1.variation.is <- tl1.variation.is[tl1.variation.is!="glottocode"] # glottocode is of course fine
# --> WALS features need to be able to be crossreferenced
tl1.variation.is[which(lapply(str_split(tl1.variation.is," "),length)>1)] <- unlist(lapply(strsplit(tl1.variation.is[which(lapply(str_split(tl1.variation.is," "),length)>1)]," "),function(x)(x[1])))
tl1.variation.is <- tl1.variation.is[tl1.variation.is!="glottocode"]
tl1.variation.should <- filter(logical_decisions,modification.ID == "A-IV1") %>% select(resulting.removed.features) %>% str_split(pattern=", ") %>% unlist %>% str_split("o_") %>% unlist %>% unique
tl1.variation.should <- tl1.variation.should[tl1.variation.should!=""]
expect_true(all(tl1.variation.should%in%tl1.variation.is))
expect_true(all(tl1.variation.is%in%tl1.variation.should))

# sanity check too little variation // 2
tl2.variation.is <- filter(nrstates,count_second_largest_feature_state=="2")$feature
# --> WALS features need to be able to be crossreferenced
tl2.variation.is[which(lapply(str_split(tl2.variation.is," "),length)>1)] <- unlist(lapply(strsplit(tl2.variation.is[which(lapply(str_split(tl2.variation.is," "),length)>1)]," "),function(x)(x[1])))
tl2.variation.should <- filter(logical_decisions,modification.ID == "A-IV2") %>% select(resulting.removed.features) %>% str_split(pattern=", ") %>% unlist %>% str_split("o_") %>% unlist %>% unique
tl2.variation.should <- tl2.variation.should[tl2.variation.should!=""]
expect_true(all(tl2.variation.should%in%tl2.variation.is))
expect_true(all(tl2.variation.is%in%tl2.variation.should))

# sanity check <100 lgs
remaining_lgs_nrstates <- nrstates %>% filter(count_second_largest_feature_state%in%c("0","1","2")==F) %>% filter(number_of_feature_states!="1")
too.few.lgs.is <- filter(remaining_lgs_nrstates,number_of_lgs_coded < 100)$feature
# --> WALS features need to be able to be crossreferenced
too.few.lgs.is[which(lapply(str_split(too.few.lgs.is," "),length)>1)] <- unlist(lapply(strsplit(too.few.lgs.is[which(lapply(str_split(too.few.lgs.is," "),length)>1)]," "),function(x)(x[1])))
too.few.lgs.should <- filter(logical_decisions,modification.ID == "A-U100LG") %>% select(resulting.removed.features) %>% str_split(pattern=", ") %>% unlist %>% str_split("o_") %>% unlist %>% unique
too.few.lgs.should <- too.few.lgs.should[too.few.lgs.should!=""]
expect_true(all(too.few.lgs.should%in%too.few.lgs.is))
expect_true(all(too.few.lgs.is%in%too.few.lgs.should))

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


## untestable dependencies match - check remaining dependencies are appropriately logged
ntrd <- filter(all_decisions, resulting.modification == "modification.ID logged in not.testable.remaining.dependencies")
for (rd in 1:nrow(ntrd)){
  id <- ntrd %>% slice(rd) %>% select(modification.ID)
  should_associated <- ntrd %>% slice(rd) %>% select(c(feature.1.for.test,feature.2.for.test)) %>% as.character
  is_associated <- recode_patterns %>% slice(which(grepl(id,recode_patterns$suspected.but.untestable.dependencies))) %>% select(new.name) %>% unlist %>% as.character
  expect_true(all(is_associated%in%should_associated))
  expect_true(all(should_associated%in%is_associated))
}
expect_true(all(ntrd$modification.ID%in%unlist(strsplit(recode_patterns$suspected.but.untestable.dependencies,";"))))
expect_true(all(na.omit(unique(unlist(strsplit(recode_patterns$suspected.but.untestable.dependencies,";"))))%in%ntrd$modification.ID))

########## save data as language-feature matrices ##########
# for language-feature matrices, rename features to short name
names(logical_data)[2:ncol(logical_data)] <- as.character(sapply(names(logical_data)[2:ncol(logical_data)],function(x)filter(recode_patterns,new.name==x)$short.name))
names(statistical_data)[2:ncol(statistical_data)] <- as.character(sapply(names(statistical_data)[2:ncol(statistical_data)],function(x)filter(recode_patterns,new.name==x)$short.name))

# save logical and statistical datasets as language-feature matrices (.csv)
write.csv(logical_data,"../../curated_data/TypLinkInd/logicalTLI/full/logicalTLI_full.csv")
write.csv(statistical_data,"../../curated_data/TypLinkInd/statisticalTLI/full/statisticalTLI_full.csv")

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
parameters <- read.csv("feature-recode-patterns.csv") %>% select(1:24)

parameters_logical <- parameters %>% filter(design.logical==T)
parameters_statistical <- parameters %>% filter(design.logical.statistical==T)

# values.csv
logical_long <- melt(setDT(logical_data), id.vars = "glottocode", variable.name = "short.name")
logical_long$value_ID <- apply(logical_long,1,function(x) paste(x[2],x[1],sep="-"))
logical_long$code_ID <- apply(logical_long,1,function(x) paste(x[2],x[3],sep="-"))
logical_long <- logical_long %>% select(c("value_ID","glottocode","short.name","value","code_ID"))
logical_long$glottocode <- as.character(logical_long$glottocode)
logical_long$short.name <- as.character(logical_long$short.name)

statistical_long <- melt(setDT(statistical_data), id.vars = "glottocode", variable.name = "short.name")
statistical_long$value_ID <- apply(statistical_long,1,function(x) paste(x[2],x[1],sep="-"))
statistical_long$code_ID <- apply(statistical_long,1,function(x) paste(x[2],x[3],sep="-"))
statistical_long$glottocode <- as.character(statistical_long$glottocode)
statistical_long$short.name <- as.character(statistical_long$short.name)
statistical_long <- statistical_long %>% select(c("value_ID","glottocode","short.name","value","code_ID"))

# codes.csv
logical_codes <- logical_long %>% select(c("code_ID","short.name","value")) %>% unique()

statistical_codes <- statistical_long %>% select(c("code_ID","short.name","value")) %>% unique()

# modifications.csv
modifications <- read.csv("decisions-log.csv")

# cldf quality checks:
expect_true(all(unique(logical_long$glottocode) %in% taxonomy_logical$glottocode))
expect_true(all(unique(statistical_long$glottocode) %in% taxonomy_statistical$glottocode))
expect_true(all(taxonomy_logical$glottocode %in% unique(logical_long$glottocode)))
expect_true(all(taxonomy_statistical$glottocode %in% unique(statistical_long$glottocode)))

expect_true(all(unique(logical_long$short.name) %in% parameters$short.name))
expect_true(all(unique(statistical_long$short.name) %in% parameters$short.name))
expect_true(all(filter(parameters,design.logical == "TRUE")$short.name %in% unique(logical_long$short.name)))
expect_true(all(filter(parameters,design.logical.statistical == "TRUE")$short.name %in% unique(statistical_long$short.name)))

expect_true(all(unique(logical_long$code_ID) %in% logical_codes$code_ID))
expect_true(all(unique(statistical_long$code_ID) %in% statistical_codes$code_ID))
expect_true(all(logical_codes$code_ID %in% unique(logical_long$code_ID)))
expect_true(all(statistical_codes$code_ID %in% unique(statistical_long$code_ID)))

expect_true(all(unique(unlist(strsplit(filter(parameters,modification.IDs!="")$modification.IDs,";"))) %in% modifications$modification.ID))
expect_true(all(unique(unlist(strsplit(filter(parameters,associated.modification.IDs.without.resulting.action!="")$associated.modification.IDs.without.resulting.action,";"))) %in% modifications$modification.ID))
expect_true(all(modifications$modification.ID %in% c(unique(unlist(strsplit(filter(parameters,modification.IDs!="")$modification.IDs,";"))),
                                                     unique(unlist(strsplit(filter(parameters,associated.modification.IDs.without.resulting.action!="")$associated.modification.IDs.without.resulting.action,";"))))))

# write all cldf components:
write.csv(taxonomy_logical,"../../curated_data/TypLinkInd/logicalTLI/cldf/languages.csv",fileEncoding="UTF-8",row.names = F)
write.csv(taxonomy_statistical,"../../curated_data/TypLinkInd/statisticalTLI/cldf/languages.csv",fileEncoding="UTF-8",row.names = F)
write.csv(parameters_logical,"../../curated_data/TypLinkInd/logicalTLI/cldf/parameters.csv",fileEncoding="UTF-8",row.names = F)
write.csv(parameters_statistical,"../../curated_data/TypLinkInd/statisticalTLI/cldf/parameters.csv",fileEncoding="UTF-8",row.names = F)
write.csv(logical_long,file="../../curated_data/TypLinkInd/logicalTLI/cldf/values.csv",row.names = F)
write.csv(statistical_long,file="../../curated_data/TypLinkInd/statisticalTLI/cldf/values.csv",row.names = F)
write.csv(logical_codes,"../../curated_data/TypLinkInd/logicalTLI/cldf/codes.csv",fileEncoding="UTF-8",row.names = F)
write.csv(statistical_codes,"../../curated_data/TypLinkInd/statisticalTLI/cldf/codes.csv",fileEncoding="UTF-8",row.names = F)
write.csv(modifications,"../../curated_data/TypLinkInd/logicalTLI/cldf/modifications.csv",fileEncoding="UTF-8",row.names = F)
write.csv(modifications,"../../curated_data/TypLinkInd/statisticalTLI/cldf/modifications.csv",fileEncoding="UTF-8",row.names = F)
