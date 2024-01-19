## this file stores all custom-written functions relevant for this repository.

# NA conversions (? to NA, "NA" to NA, blank to NA)
na_convert <- function(data){
  data[data=="?"]<-NA
  data[data=="NA"]<-NA
  data[data==""]<-NA
  data[is.na(data)]<-NA
  return(data)
}

# make all columns become factors
factorise <- function(data){
  for(i in 1:ncol(data)){
    data[,i]<-as.factor(data[,i])
  }
  return(data)
}

# this function serves to replace redundant glottocodes with new ones and merging data from both glottocodes (prioritising the data for the new glottocode from the original dataframe, if available)
update_glottocode_and_data <- function(data, old, new){
  data$glottocode <- as.character(data$glottocode)
  if(new%in%data$glottocode){ # if the new glottocode is already in the matrix, procede as follows:
    expect_true(new%in%data$glottocode)
    expect_true(old%in%data$glottocode)
    original <- data %>% filter(glottocode%in%c(old,new)) # subset original data to relevant glottocodes
    data <- data %>% filter(glottocode%in%c(old,new)==F)
    original <- original[c(which(original$glottocode==new), which(original$glottocode==old)),] # ensure the new glottocode is on top
    original[1,] <- apply(original,2,function(x){if(x[1]=="?"|is.na(x[1])){x[2]}else{x[1]}}) %>% as.character() # replace top value by bottom value if top value is ? or NA
    data <- rbind(data,original[1,])
    expect_true(new%in%data$glottocode)
    expect_false(old%in%data$glottocode)
  }
  if(new%in%data$glottocode==F){ # if the new glottocode is not in the matrix, just update the name
    expect_false(new%in%data$glottocode)
    expect_true(old%in%data$glottocode)
    data[data$glottocode==old,"glottocode"] <- new # replace glottocode
    expect_true(new%in%data$glottocode)
    expect_false(old%in%data$glottocode)
  }
  return(data)
}

# this function serves to recode states from one or several original variables into other states, as specified in modifications.csv
implement_recode <- function(original_data, expected_levels, recoding_groups, recoded_levels, nvar, recode_mode){
  
  # make a table of original values
  expected_levels <- data.frame(i = as.integer(gsub("^([0-9])*([0-9])+.+$", "\\1\\2", expected_levels)),
                                level = gsub("^[0-9]+\\.? +", "", expected_levels),
                                stringsAsFactors=FALSE)
  
  # sanity checks
  expect_true(all(!is.na(expected_levels$i)))
  expect_true(all(!is.na(expected_levels$level)))
  
  # make sure that the expected values match the original values found (applies only to simple recode)
  if (recode_mode=="simple"){
    if(nvar=="single"){
      expect_true(setequal(expected_levels$level, na.omit(original_data)), info=
                    paste0("Expected:\n", paste0("  ", (expected_levels$level), collapse="\n"), "\n",
                           "Got:\n",  paste0("  ", (unique(original_data)), collapse="\n")))
    }
    if(nvar=="multiple"){
      expect_true(all(unique(na.omit(original_data$merged)) %in% expected_levels$level), info=
                    paste0("Expected:\n", paste0("  ", (expected_levels$level), collapse="\n"), "\n",
                           "Got:\n",  paste0("  ", (unique(original_data)), collapse="\n")))
    }
  }
  
  # parse the recoding pattern
  recoding_groups <- strsplit(recoding_groups, "/") %>% lapply(as.integer)
  
  # sanity checks
  expect_true(length(recoding_groups)>1) # must have at least 2 recoding groups
  expect_true(all(!is.na(unlist(recoding_groups)))) # can't have NAs
  expect_true(all(unlist(recoding_groups) %in% expected_levels$i)) # must correspond to original values
  expect_false(any(duplicated(unlist(recoding_groups)))) # can't have any duplicates
  
  # build the recoding table
  expect_true(length(recoding_groups)==length(recoded_levels)) # must have at least 2 recoding groups
  recoded_levels <- bind_rows(mapply(recoded_levels, recoding_groups, FUN=function(value, ii) {
    data.frame(i = ii, new_level=as.character(value), stringsAsFactors=FALSE)
  }, SIMPLIFY=FALSE))
  
  level_table <- full_join(expected_levels, recoded_levels, by="i")
  
  # sanity checks
  expect_true(all(!is.na(level_table$level)))
  
  # recode the data
  if (recode_mode == "simple"){
    if (nvar=="single"){
      new_data <- level_table$new_level[match(original_data, level_table$level)]
    }
    if (nvar=="multiple"){
      new_data <- level_table$new_level[match(original_data$merged, level_table$level)]
    }
  }
  
  if (recode_mode == "logical_arguments"){
    # recode the data according to prioritised variable
    for (i in 1:nrow(level_table)){
      original_data[original_data$glottocode%in%filter(original_data,eval(parse(text=level_table$level[i])))$glottocode,"merged"]<-level_table$new_level[i]
    }
    # make "other" state become "?" if applicable
    original_data[is.na(original_data)]<-"?"
    new_data <- original_data$merged
  }
  return(new_data)
}

# this function serves to extract "condition" and "equator" from a condition statement for further use
extract_condition_and_equator <- function(condition_statement){
  condition <- unlist(strsplit(condition_statement," == "))
  equator <- " == "
  # if condition has not been split, it is not positive; check whether it is negative
  if (length(condition)==1){ 
    condition <- unlist(strsplit(condition_statement," != "))
    equator <- " != "
  }
  # if condition has not been split, it is not positive or negative, check whether it is positive multiple
  if (length(condition)==1){
    condition <- unlist(strsplit(condition_statement," %in% "))
    equator <- " %in% "
  }
  # if condition has not been split, it is not positive or negative or positive mulitple, so it is negative multiple
  if (length(condition)==1){
    condition <- unlist(strsplit(condition_statement," !%in% "))
    equator <- " !%in% "
  } 
  
  return(list(condition,equator))
}

# this function serves to condition a variable on another -- note that the currently implemented function works for up to 5 desired states in the %in% case
implement_conditioning <- function(variable_to_be_conditioned, condition, equator){
  
  # select conditioned upon variable
  conditioned_upon_variable <- recoded_data[,c(1,which(names(recoded_data)%in%condition[1]))]
  
  # select glottocodes for which condition applies and turn data into "?" where applicable
  if (equator == " == "){ 
    # if the condition in question is positive (" == "), we want to keep languages that have the desired state of conditioned_upon_variable OR which are "?" to both conditioned_upon_variable and variable_to_be_conditioned to not become NA
    # select languages with desired state or "?" in conditioned_upon_variable
    condition_applies_strict <- filter(conditioned_upon_variable,conditioned_upon_variable[,2]==condition[2])$glottocode
    condition_applies_q <- filter(conditioned_upon_variable,conditioned_upon_variable[,2]=="?")$glottocode
    
    # select languages in variable_to_be_conditioned to which condition applies (strict and q)
    conditioned_data <- filter(variable_to_be_conditioned, glottocode %in% as.character(c(condition_applies_strict,condition_applies_q)))
    names(conditioned_data)[2] <- "conditioned_upon_variable"
    
    # the languages, which are "?" to conditioned_upon_variable but specified for variable_to_be_conditioned are recoded into "?"
    conditioned_data$conditioned_upon_variable[conditioned_data$glottocode%in%condition_applies_q] <- rep("?")
    
  } else if (equator == " != "){ ## this applies if the condition in question is negative (" != ")
    # if the condition in question is negative (" != "), we want to keep all languages that do not have the specified state of conditioned_upon_variable
    # select languages which do not have the specified state in conditioned_upon_variable
    condition_applies <- setdiff(conditioned_upon_variable$glottocode,filter(conditioned_upon_variable,conditioned_upon_variable[,2]==condition[2])$glottocode)
    
    # select languages in variable_to_be_conditioned to which condition applies
    conditioned_data <- filter(variable_to_be_conditioned, glottocode %in% as.character(condition_applies))
    
  } else if (equator == " %in% "){ ## this applies if the condition in the question is multiple --> conservative ("%in%")
    
    # if the condition in question is multiple (" %in% "), we want to keep languages that have any of the desired state of conditioned_upon_variable OR which are "?" to both conditioned_upon_variable and variable_to_be_conditioned to not become NA
    desired_states <- unlist(strsplit(condition[2],", "))
    nr_desired_states <- length(desired_states)
    
    # select languages with desired states or "?" in conditioned_upon_variable
    condition_applies_q <- filter(conditioned_upon_variable,conditioned_upon_variable[,2]=="?")$glottocode
    condition_applies_desired_states <- filter(conditioned_upon_variable,conditioned_upon_variable[,2]==desired_states[1]|conditioned_upon_variable[,2]==desired_states[2])$glottocode
    # if there are more than 2 desired states, add the third
    if(nr_desired_states>2){condition_applies_desired_states <- c(condition_applies_desired_states,filter(conditioned_upon_variable,conditioned_upon_variable[,2]==desired_states[3])$glottocode)}
    # if there are more than 3 desired states, add the fourth
    if(nr_desired_states>3){condition_applies_desired_states <- c(condition_applies_desired_states,filter(conditioned_upon_variable,conditioned_upon_variable[,2]==desired_states[4])$glottocode)}
    # if there are more than 4 desired states, add the fifth
    if(nr_desired_states>4){condition_applies_desired_states <- c(condition_applies_desired_states,filter(conditioned_upon_variable,conditioned_upon_variable[,2]==desired_states[5])$glottocode)}
    
    # select languages in variable_to_be_conditioned to which condition applies (strict and q)
    conditioned_data <- filter(variable_to_be_conditioned, glottocode %in% as.character(c(condition_applies_q, condition_applies_desired_states)))
    names(conditioned_data)[2] <- "conditioned_upon_variable"
    
    # the languages, which are "?" to conditioned_upon_variable but specified for variable_to_be_conditioned are recoded into "?"
    conditioned_data$conditioned_upon_variable[conditioned_data$glottocode%in%condition_applies_q] <- rep("?")
    
  } else if (equator == " !%in% "){ ## this applies if the condition in the question is multiple (but negative) --> liberal ("!%in%")
    # if the condition in question is negative multiple (" !%in% "), we want to keep all languages that do not have the specified states of conditioned_upon_variable
    # select languages which do not have the specified state in conditioned_upon_variable
    undesired_states <- unlist(strsplit(condition[2],", "))
    nr_undesired_states <- length(undesired_states)
    condition_applies_undesired_states <- filter(conditioned_upon_variable,conditioned_upon_variable[,2]==undesired_states[1]|conditioned_upon_variable[,2]==undesired_states[2])$glottocode
    # if there are more than 2 undesired states, add the third
    if(nr_undesired_states>2){condition_applies_undesired_states <- c(condition_applies_undesired_states,filter(conditioned_upon_variable,conditioned_upon_variable[,2]==undesired_states[3])$glottocode)}
    # if there are more than 3 undesired states, add the fourth
    if(nr_undesired_states>3){condition_applies_undesired_states <- c(condition_applies_undesired_states,filter(conditioned_upon_variable,conditioned_upon_variable[,2]==undesired_states[4])$glottocode)}
    # if there are more than 4 undesired states, add the fifth
    if(nr_undesired_states>4){condition_applies_undesired_states <- c(condition_applies_undesired_states,filter(conditioned_upon_variable,conditioned_upon_variable[,2]==undesired_states[5])$glottocode)}
    condition_applies <- setdiff(conditioned_upon_variable$glottocode,condition_applies_undesired_states)
    # select languages in variable_to_be_conditioned to which condition applies
    conditioned_data <- filter(variable_to_be_conditioned, glottocode %in% as.character(condition_applies))
  }
  
  return(conditioned_data)
}

# this function runs all tests
evaluate_XOR_AND_THEN <- function(recoded_data, expectations, taxonomy, diversity_samples, proportion_languages_must_be_in_applicable_state){
  library(lsr)
  rslt <- slice(data.frame(nr.lgs.v1=NA,
                           nr.lgs.v2=NA,
                           absolute.overlap.v1.v2=NA,
                           relative.overlap.vsmall.in.vlarge=NA,
                           relative.overlap.vlarge.in.vsmall=NA,
                           samples.disregarded.insufficient.lgs.applicable=NA,
                           samples.assessed=NA,
                           mean.d1.across.tested.samples=NA,
                           sd.d1.across.tested.samples=NA,
                           mean.plus.sd.d1.across.tested.samples=NA,
                           mean.d2.across.tested.samples=NA,
                           sd.d2.across.tested.samples=NA,
                           mean.plus.sd.d2.across.tested.samples=NA,
                           mean.baseline.1=NA,
                           mean.baseline.2=NA,
                           d1.baseline1.cohensd=NA,
                           d2.baseline2.cohensd=NA),0)
  
  for (i in 1:nrow(expectations)){ # we loop through all expectations
    cat("Processing expectation ", i, "\n", sep="")
    expectation <- expectations[i,] # select expectation
    raw_data <- recoded_data %>% select(c(expectation$variable.1.for.test,expectation$variable.2.for.test,glottocode)) # subset data to relevant variables, including ? and NA
    
    full_tbl_withqs <- table(unlist(raw_data[,1]),unlist(raw_data[,2])) # cross-tabulate variables using full data (including ? and NA)
    
    # log characteristics about language overlap: for certain
    # for certain considerations,  we do not care about "?" and "NA" --> convert them to NA
    non_NAq_data <- raw_data
    non_NAq_data[non_NAq_data =="?"] <- NA
    non_NAq_data[non_NAq_data =="NA"] <- NA
    
    # log how many languages each variable is coded for (nrlgs1 for variable 1 and nrlgs2 for variable 2; these values are logged)
    nrlgs1 <- nrow(na.omit(as.data.frame(non_NAq_data[,1])))
    nrlgs2 <- nrow(na.omit(as.data.frame(non_NAq_data[,2])))
    
    # determine the number of languages coded for both variables (this value is logged)
    absoverlap <- non_NAq_data %>% filter(!is.na(get(expectation$variable.1.for.test))) %>% filter(!is.na(get(expectation$variable.2.for.test))) %>% nrow()
    
    # log the number of languages in the smaller and the number of languages in the larger variable; then determine the relative overlap of languages for the smaller and larger variable (these values are logged)
    smaller <- sort(c(nrlgs1,nrlgs2))[1]
    larger <- sort(c(nrlgs1,nrlgs2))[2]
    relative.overlap.vsmall.in.vlarge <- absoverlap/smaller
    relative.overlap.vlarge.in.vsmall <- absoverlap/larger
    
    # convert the feature values into factor format
    non_NAq_data[,expectation$variable.1.for.test] <- as.factor(unlist(non_NAq_data[,expectation$variable.1.for.test]))
    non_NAq_data[,expectation$variable.2.for.test] <- as.factor(unlist(non_NAq_data[,expectation$variable.2.for.test]))
    
    raw_data[,expectation$variable.1.for.test] <- as.factor(unlist(raw_data[,expectation$variable.1.for.test]))
    raw_data[,expectation$variable.2.for.test] <- as.factor(unlist(raw_data[,expectation$variable.2.for.test]))
    
    # prepare a table to log relevant characteristics of the crosstable for each of the diversity samples
    expectation_assessment<-slice(data.frame(proportion_overlap_applicable_vs_non_applicable=NA, # this logs for each sample the proportion of languages coded for v2 in the relevant state of v1
                                             overlap_sufficient_test_power_positive=NA, # this denotes for each sample whether the proportion and number of languages coded for v2 in the relevant state of v1 is sufficient, the proportion is defined by the user and denoted proportion_languages_must_be_in_applicable_state
                                             result1_XOR_AND_THEN=NA, # this denotes the proportion of languages behaving against the expectation in the THEN condition and in the first dimension of the XOR and AND conditions
                                             result2_XOR_AND=NA, # this denotes the proportion of languages behaving against the expectation in the second dimension of the XOR and AND conditions
                                             baseline1_XOR_AND_THEN=NA, # this denotes the proportion of languages in the "unequals" state for variable 2 in the THEN condition overall or in the first dimension of the XOR condition
                                             baseline2_XOR_AND=NA),0) # this denotes whether the expected relationship is confirmed in the subsample in accordance with the defined thresholds
    
    for (ds in 1:ncol(diversity_samples)){ # loop through all diversity samples
      # select the subset of languages in the diversity sample for the relevant variables and tabulate for both the full data (including ? and NA) and the non-Q/NA data.
      expectation_ds_sample_withqs <- raw_data %>% filter(glottocode%in%diversity_samples[,ds])
      sample_table_withqs <- table(unlist(expectation_ds_sample_withqs[,1]),unlist(expectation_ds_sample_withqs[,2])) # this is the table
      
      expectation_ds_sample_noqs <- non_NAq_data %>% filter(glottocode%in%diversity_samples[,ds])
      sample_table_noqs <- table(unlist(expectation_ds_sample_noqs[,1]),unlist(expectation_ds_sample_noqs[,2])) # this is the table
      
      # sanity checks
      expect_true(all(levels(as.data.frame(sample_table_noqs)$Var1) %in% c(unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")), unlist(strsplit(as.character(unlist(expectation$v1.unequals.for.test)),", ")))))
      expect_true(all(c(unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")), unlist(strsplit(as.character(unlist(expectation$v1.unequals.for.test)),", "))) %in% levels(as.data.frame(sample_table_noqs)$Var1)))
      expect_true(all(c(unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", ")), unlist(strsplit(as.character(unlist(expectation$v2.unequals.for.test)),", "))) %in% levels(as.data.frame(sample_table_noqs)$Var2)))
      expect_true(all(levels(as.data.frame(sample_table_noqs)$Var2) %in% c(unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", ")), unlist(strsplit(as.character(unlist(expectation$v2.unequals.for.test)),", ")))))
      
      if(expectation$test == "XOR"){ # for XOR-case
        # in the XOR-case, we evaluate how many languages are coded for both variables in the first dimension (applicable1) and how many languages are coded for just the first variable in the first dimension (all1)
        applicable1 <- sum(sample_table_withqs[unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),c(unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$v2.unequals.for.test)),", ")))])
        all1 <- sum(sample_table_withqs[unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),])
        # this yields the relevant proportion in the first dimension
        relevant.proportion.1 <- applicable1/all1
        
        # in the XOR-case, we also evaluate how many languages are coded for both variables in the second dimension (applicable2) and how many languages are coded for just the second variable in the second dimension (all2)
        applicable2 <- sum(sample_table_withqs[c(unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$v1.unequals.for.test)),", "))),unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", "))])
        all2 <- sum(sample_table_withqs[,unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", "))])
        # this yields the relevant proportion in the second dimension
        relevant.proportion.2 <- applicable2/all2
        
        # the lower of both proportions is relevant and logged
        relevant.proportion <- if(length(na.omit(c(relevant.proportion.1,relevant.proportion.2)))==1){0}else{sort(c(relevant.proportion.1,relevant.proportion.2))[1]}
        expectation_assessment[ds,"proportion_overlap_applicable_vs_non_applicable"] <- relevant.proportion
        expectation_assessment[ds,"overlap_sufficient_test_power_positive"] <- relevant.proportion >= proportion_languages_must_be_in_applicable_state
        
        # in the XOR-condition, the languages that behave against the expectation are those with the "equals" state of both variable 1 and variable 2
        hyp.NOT.applies <- sample_table_noqs[unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", "))]
        
        # the languages that behave according to the expectation are those with the "equals" state of variable 1 and the "unequals" state of variable 2 (hyp.applies.1) and those with the "unequals" state of variable 1 and the "equals" state of variable 2 (hyp.applies.2)
        hyp.applies.1 <- sample_table_noqs[unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$v2.unequals.for.test)),", "))]
        hyp.applies.2 <- sample_table_noqs[unlist(strsplit(as.character(unlist(expectation$v1.unequals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", "))]
        
        # track the proportion of cases, in which the expectation is violated with respect to variable 1 (d1) and with respect to variable 2 (d2)
        d1 <- sum(hyp.NOT.applies)/(sum(hyp.NOT.applies)+sum(hyp.applies.1)) 
        d2 <- sum(hyp.NOT.applies)/(sum(hyp.NOT.applies)+sum(hyp.applies.2))
        expectation_assessment[ds,"result1_XOR_AND_THEN"] <- d1
        expectation_assessment[ds,"result2_XOR_AND"] <- d2
        
        # the proportions must be compared to the baseline expectations: these are the proportions of languages in the "unequals" state of variable 2 and variable 1 overall (independently of variable 1 or 2 being coded or not, and which state it would be in)
        baseline1_v2_equals <- expectation_ds_sample_withqs %>% 
          filter(get(expectation$variable.2.for.test)%in%
                   unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", "))) %>% nrow()
        baseline1_v2_unequals <- expectation_ds_sample_withqs %>% 
          filter(get(expectation$variable.2.for.test)%in%
                   unlist(strsplit(as.character(unlist(expectation$v2.unequals.for.test)),", "))) %>% nrow()
        baseline1_v2 <- baseline1_v2_unequals/(baseline1_v2_unequals+baseline1_v2_equals)
        
        baseline2_v1_equals <- expectation_ds_sample_withqs %>% 
          filter(get(expectation$variable.1.for.test)%in%
                   unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", "))) %>% nrow()
        baseline2_v1_unequals <- expectation_ds_sample_withqs %>% 
          filter(get(expectation$variable.1.for.test)%in%
                   unlist(strsplit(as.character(unlist(expectation$v1.unequals.for.test)),", "))) %>% nrow()
        baseline2_v1 <- baseline2_v1_unequals/(baseline2_v1_unequals+baseline2_v1_equals)
        
        # log the baselines
        expectation_assessment[ds,"baseline1_XOR_AND_THEN"] <- baseline1_v2
        expectation_assessment[ds,"baseline2_XOR_AND"] <- baseline2_v1
        
      }
      else if(expectation$test == "AND"){
        # in the AND-case, we evaluate how many languages are coded for the second variable, given the relevant state of the first variable (applicable), and how many languages are coded for the relevant state of the first variable overall (all), and vice versa.
        applicable.1 <- sum(sample_table_withqs[unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),c(unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$v2.unequals.for.test)),", ")))])
        all.1 <- sum(sample_table_withqs[unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),])
        
        applicable.2 <- sum(sample_table_withqs[c(unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$v1.unequals.for.test)),", "))), unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", "))])
        all.2 <- sum(sample_table_withqs[,unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", "))])
        
        # this yields the relevant proportion in both dimensions
        relevant.proportion.1 <- applicable.1/all.1
        relevant.proportion.2 <- applicable.2/all.2
        
        # the lower of both proportions is relevant and logged
        relevant.proportion <- if(length(na.omit(c(relevant.proportion.1,relevant.proportion.2)))==1){0}else{sort(c(relevant.proportion.1,relevant.proportion.2))[1]}
        expectation_assessment[ds,"proportion_overlap_applicable_vs_non_applicable"] <- relevant.proportion
        expectation_assessment[ds,"overlap_sufficient_test_power_positive"] <- relevant.proportion >= proportion_languages_must_be_in_applicable_state
        
        # in the AND-condition, the languages that behave according to the expectation are those with the "equals" state of both variable 1 and variable 2
        hyp.applies <- sample_table_noqs[unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", "))]
        
        # the languages that behave against the expectation are those with the "equals" state of variable 1 and the "unequals" state of variable 2, and vice versa
        hyp.NOT.applies.1 <- sample_table_noqs[unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$v2.unequals.for.test)),", "))]
        hyp.NOT.applies.2 <- sample_table_noqs[unlist(strsplit(as.character(unlist(expectation$v1.unequals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", "))]
        
        # track the proportion of cases, in which the expectation is violated (d1 and d2)
        d1 <- sum(hyp.NOT.applies.1)/(sum(hyp.applies)+sum(hyp.NOT.applies.1)) # proportion of cases, in which gut expectation is violated with respect to variable 1
        d2 <- sum(hyp.NOT.applies.2)/(sum(hyp.applies)+sum(hyp.NOT.applies.2)) # proportion of cases, in which gut expectation is violated with respect to variable 1
        expectation_assessment[ds,"result1_XOR_AND_THEN"]<-d1
        expectation_assessment[ds,"result2_XOR_AND"]<-d2
        
        # the proportion must be compared to the baseline expectation: this is the proportion of languages in the "unequals" state of variable 2 overall (independently of variable 1 being coded or not, and which state it would be in)
        baseline_v2_equals <- expectation_ds_sample_withqs %>% filter(get(expectation$variable.2.for.test)%in%unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", "))) %>% nrow()
        baseline_v2_unequals <- expectation_ds_sample_withqs %>% filter(get(expectation$variable.2.for.test)%in%unlist(strsplit(as.character(unlist(expectation$v2.unequals.for.test)),", "))) %>% nrow()
        baseline_v2 <- baseline_v2_unequals/(baseline_v2_unequals+baseline_v2_equals)
        
        baseline_v1_equals <- expectation_ds_sample_withqs %>% filter(get(expectation$variable.1.for.test)%in%unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", "))) %>% nrow()
        baseline_v1_unequals <- expectation_ds_sample_withqs %>% filter(get(expectation$variable.1.for.test)%in%unlist(strsplit(as.character(unlist(expectation$v1.unequals.for.test)),", "))) %>% nrow()
        baseline_v1 <- baseline_v1_unequals/(baseline_v1_unequals+baseline_v1_equals)
        
        
        # log the baseline, and whether the baseline is higher than d1
        expectation_assessment[ds,"baseline1_XOR_AND_THEN"] <- baseline_v2
        expectation_assessment[ds,"baseline2_XOR_AND"] <- baseline_v1
      }
      else if(expectation$test == "THEN"){
        # in the THEN-case, we evaluate how many languages are coded for the second variable, given the relevant state of the first variable (applicable), and how many languages are coded for the relevant state of the first variable overall (all)
        applicable <- sum(sample_table_withqs[unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),c(unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$v2.unequals.for.test)),", ")))])
        all <- sum(sample_table_withqs[unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),])
        
        # this yields the relevant proportion, which is logged
        relevant.proportion <- applicable/all
        expectation_assessment[ds,"proportion_overlap_applicable_vs_non_applicable"] <- relevant.proportion
        expectation_assessment[ds,"overlap_sufficient_test_power_positive"] <- relevant.proportion >= proportion_languages_must_be_in_applicable_state
        
        # in the THEN-condition, the languages that behave according to the expectation are those with the "equals" state of both variable 1 and variable 2
        hyp.applies <- sample_table_noqs[unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", "))]
        
        # the languages that behave against the expectation are those with the "equals" state of variable 1 and the "unequals" state of variable 2
        hyp.NOT.applies <- sample_table_noqs[unlist(strsplit(as.character(unlist(expectation$v1.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$v2.unequals.for.test)),", "))]
        
        # track the proportion of cases, in which the expectation is violated (d1)
        d1 <- sum(hyp.NOT.applies)/(sum(hyp.applies)+sum(hyp.NOT.applies)) # proportion of cases, in which gut expectation is violated with respect to variable 1
        expectation_assessment[ds,"result1_XOR_AND_THEN"]<-d1
        
        # the proportion must be compared to the baseline expectation: this is the proportion of languages in the "unequals" state of variable 2 overall (independently of variable 1 being coded or not, and which state it would be in)
        baseline_v2_equals <- expectation_ds_sample_withqs %>% filter(get(expectation$variable.2.for.test)%in%unlist(strsplit(as.character(unlist(expectation$v2.equals.for.test)),", "))) %>% nrow()
        baseline_v2_unequals <- expectation_ds_sample_withqs %>% filter(get(expectation$variable.2.for.test)%in%unlist(strsplit(as.character(unlist(expectation$v2.unequals.for.test)),", "))) %>% nrow()
        baseline_v2 <- baseline_v2_unequals/(baseline_v2_unequals+baseline_v2_equals)
        
        # log the baseline, and whether the baseline is higher than d1
        expectation_assessment[ds,"baseline1_XOR_AND_THEN"] <- baseline_v2
      }
    }
    
    # log relevant characteristics for the expectation:  aggregated for each diversity sample
    rslt=rbind(rslt,data.frame(nr.lgs.v1=nrlgs1, # raw
                               nr.lgs.v2=nrlgs2, # raw
                               absolute.overlap.v1.v2=absoverlap, # raw
                               relative.overlap.vsmall.in.vlarge=relative.overlap.vsmall.in.vlarge, # raw
                               relative.overlap.vlarge.in.vsmall=relative.overlap.vlarge.in.vsmall, # raw
                               samples.disregarded.insufficient.lgs.applicable=nrow(expectation_assessment)-sum(na.omit(expectation_assessment$overlap_sufficient_test_power_positive)), # aggregation: how many subsamples were discarded?
                               samples.assessed=sum(na.omit(expectation_assessment$overlap_sufficient_test_power_positive)), # aggregation: how many subsamples could be computed?
                               mean.d1.across.tested.samples=round(mean(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==T)$result1_XOR_AND_THEN)),digits=2), # aggregation: what was the mean proportion of cases, in which the expectation was violated, across all samples? (with respect to v1 if XOR or AND)
                               sd.d1.across.tested.samples=round(sd(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==T)$result1_XOR_AND_THEN)),digits=2), # aggregation: what was the standard deviation of the proportion of cases, in which the expectation was violated, across all samples? (with respect to v1 if XOR or AND)
                               mean.plus.sd.d1.across.tested.samples=round(mean(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==T)$result1_XOR_AND_THEN))+sd(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==T)$result1_XOR_AND_THEN)),digits=2), # aggregation: what was the sum of the mean plus one standard deviation of the proportion of cases, in which the expectation was violated, across all samples? (with respect to v1 if XOR or AND)
                               mean.d2.across.tested.samples=round(mean(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==T)$result2_XOR_AND)),digits=2), # aggregation: what was the mean proportion of cases, in which the expectation was violated, across all samples, with respect to v2? (only for XOR or AND)
                               sd.d2.across.tested.samples=round(sd(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==T)$result2_XOR_AND)),digits=2), # aggregation: what was the standard deviation of the proportion of cases, in which the expectation was violated, across all samples, with respect to v2? (only for XOR or AND)
                               mean.plus.sd.d2.across.tested.samples=round(mean(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==T)$result2_XOR_AND))+sd(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==T)$result2_XOR_AND)),digits=2), # aggregation: what was sum of the mean plus one standard deviation of the proportion of cases, in which the expectation was violated, across all samples, with respect to v2? (only for XOR or AND)
                               mean.baseline.1=round(mean(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==T)$baseline1_XOR_AND_THEN)),digits=2), # aggregation: what was the mean baseline probability for v2 being in the "unequals" state? (with respect to v1 for XOR and AND)
                               mean.baseline.2=round(mean(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==T)$baseline2_XOR_AND)),digits=2), # aggregation: what was the mean baseline probability for v1 being in the "unequals" state? (only for XOR and AND)
                               d1.baseline1.cohensd=if(length(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==T)$baseline1_XOR_AND_THEN))!=0)
                               {round(cohensD(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==TRUE)$result1_XOR_AND_THEN),na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==TRUE)$baseline1_XOR_AND_THEN), method="paired"), digits = 1)}
                               else{NA}, # what is effect size (with respect to v1 if XOR or AND), measured by paired cohen's D
                               d2.baseline2.cohensd=if(length(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==T)$baseline2_XOR_AND))!=0)
                               {round(cohensD(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==TRUE)$result2_XOR_AND),na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==TRUE)$baseline2_XOR_AND), method="paired"), digits = 1)}
                               else{NA})) # what is effect size (with respect to v2 if XOR or AND), measured by paired cohen's D)
  }
  return(rslt)
}

# comparisons between full and densified datasets
compare_full_to_densified <- function(full, densified, densified2 = NULL, taxonomy_matrix){
  library(gridExtra)
  library(ggplot2)
  
  # set all ?, "NA" and blanks to NA
  full <- na_convert(full)
  densified <- na_convert(densified)
  if(!is.null(densified2)){
    densified2 <- na_convert(densified2)
  }
  
  # coding proportions
  coding_proportion_full <- sum(!is.na(full))/(ncol(full)*nrow(full))
  coding_proportion_densified <- sum(!is.na(densified))/(ncol(densified)*nrow(densified))
  if(!is.null(densified2)){
    coding_proportion_densified2 <- sum(!is.na(densified2))/(ncol(densified2)*nrow(densified2))
  }
  
  # number of languages
  nlg_full <- nrow(full)
  nlg_densified <- nrow(densified)
  if(!is.null(densified2)){
    nlg_densified2 <- nrow(densified2)
  }
  
  # number of families
  nfam_full <- taxonomy_matrix %>% filter(id %in% rownames(full)) %>% select(level1) %>% unique() %>% nrow()
  nfam_densified <- taxonomy_matrix %>% filter(id %in% rownames(densified)) %>% select(level1) %>% unique() %>% nrow()
  if(!is.null(densified2)){
    nfam_densified2 <- taxonomy_matrix %>% filter(id %in% rownames(densified2)) %>% select(level1) %>% unique() %>% nrow()
  }
  
  # number of variables
  nvar_full <- ncol(full)
  nvar_densified <- ncol(densified)
  if(!is.null(densified2)){
    nvar_densified2 <- ncol(densified2)
  }
  
  # generate a table comparing the values
  comparison <- data.frame(data=c("full","densified"),
                           coding.density=c(coding_proportion_full,coding_proportion_densified),
                           nr.languages=c(nlg_full,nlg_densified),
                           nr.families=c(nfam_full,nfam_densified),
                           nr.variables=c(nvar_full,nvar_densified),
                           coding.density.improvement=c(NA,round(coding_proportion_densified/coding_proportion_full,digits=3)),
                           proportion.languages=c(NA,round(nlg_densified/nlg_full,digits=3)),
                           proportion.families=c(NA,round(nfam_densified/nfam_full,digits=3)),
                           proportion.variables=c(NA,round(nvar_densified/nvar_full,digits=3)))
  
  if(!is.null(densified2)){
    comparison <- data.frame(data=c("full","densified_large","densified_small"),
                             coding.density=c(coding_proportion_full,coding_proportion_densified,coding_proportion_densified2),
                             nr.languages=c(nlg_full,nlg_densified,nlg_densified2),
                             nr.families=c(nfam_full,nfam_densified,nfam_densified2),
                             nr.variables=c(nvar_full,nvar_densified,nvar_densified2),
                             coding.density.improvement=c(NA,round(coding_proportion_densified/coding_proportion_full,digits=3),round(coding_proportion_densified2/coding_proportion_full,digits=3)),
                             proportion.languages=c(NA,round(nlg_densified/nlg_full,digits=3),round(nlg_densified2/nlg_full,digits=3)),
                             proportion.families=c(NA,round(nfam_densified/nfam_full,digits=3),round(nfam_densified2/nfam_full,digits=3)),
                             proportion.variables=c(NA,round(nvar_densified/nvar_full,digits=3),round(nvar_densified2/nvar_full,digits=3)))
  }
  
  # generate coding density and family distribution comparisons
  lg_plot_full <- ggplot(data.frame(density=apply(full, 1, function(x) (length(na.omit(x))))/ncol(full)), aes(x=density))+
    geom_histogram(color="lightblue", fill="lightblue", bins=20)+
    xlim(c(0,1.05))+
    theme_minimal()+
    labs(title="Coding density per language, full dataset", 
         x="coding density per language",
         y="number of languages")
  
  var_plot_full <- ggplot(data.frame(density=apply(full, 2, function(x) (length(na.omit(x))))/nrow(full)), aes(x=density))+
    geom_histogram(color="lightblue", fill="lightblue", bins=20)+
    theme_minimal()+
    xlim(c(0,1.05))+
    labs(title="Coding density per variable, full dataset",
         x="coding density per variable",
         y="number of variables")
  
  fams_full <- as.data.frame(table(select(filter(taxonomy_matrix, id%in%rownames(full)),level1)))
  family_plot_full <-  ggplot(fams_full,aes(x=level1,y=Freq))+
    theme_minimal()+
    geom_bar(stat="identity",color="lightblue",fill="lightblue")+
    theme(panel.grid.major = element_blank(), 
          panel.background = element_blank(), axis.text.x = element_blank(),  # Remove x-axis labels
          axis.line.x = element_blank())+    # Remove x-axis line
    labs(title="Language counts per family, full dataset",y="language count per family",x="families")+
    ylim(c(0,max(fams_full$Freq)+20))
  
  if(!is.null(densified2)){
    lg_plot_densified <- ggplot(data.frame(density=apply(densified, 1, function(x) (length(na.omit(x))))/ncol(densified)), aes(x=density))+
      geom_histogram(color="steelblue", fill="steelblue", bins=20)+
      theme_minimal()+
      xlim(c(0,1.05))+
      labs(title="Coding density per language, densified dataset (large)", 
           x="coding density per language",
           y="number of languages")
    
    var_plot_densified <- ggplot(data.frame(density=apply(densified, 2, function(x) (length(na.omit(x))))/nrow(densified)), aes(x=density))+
      geom_histogram(color="steelblue", fill="steelblue", bins=20)+
      theme_minimal()+
      theme(axis.line.x = element_blank())+    # Remove x-axis line
      xlim(c(0,1.05))+
      labs(title="Coding density per variable, densified dataset (large)",
           x="coding density per variable",
           y="number of variables")
    
    fams_densified <- as.data.frame(table(select(filter(taxonomy_matrix, id%in%rownames(densified)),level1)))
    family_plot_densified <-  ggplot(fams_densified,aes(x=level1,y=Freq))+
      geom_bar(stat="identity",color="steelblue",fill="steelblue")+
      theme_minimal()+
      theme(panel.grid.major = element_blank(), 
            panel.background = element_blank(), axis.text.x = element_blank(),  # Remove x-axis labels
            axis.line.x = element_blank())+    # Remove x-axis line
      labs(title="Language counts per family, densified dataset (large)",y="language count per family",x="families")+
      ylim(c(0,max(fams_full$Freq)+20))
    
    
    lg_plot_densified2 <- ggplot(data.frame(density=apply(densified2, 1, function(x) (length(na.omit(x))))/ncol(densified2)), aes(x=density))+
      geom_histogram(color="darkblue", fill="darkblue", bins=20)+
      theme_minimal()+
      xlim(c(0,1.05))+
      labs(title="Coding density per language, densified dataset (small)", 
           x="coding density per language",
           y="number of languages")
    
    var_plot_densified2 <- ggplot(data.frame(density=apply(densified2, 2, function(x) (length(na.omit(x))))/nrow(densified2)), aes(x=density))+
      geom_histogram(color="darkblue", fill="darkblue", bins=20)+
      theme_minimal()+
      theme(axis.line.x = element_blank())+    # Remove x-axis line
      xlim(c(0,1.05))+
      labs(title="Coding density per variable, densified dataset (small)",
           x="coding density per variable",
           y="number of variables")
    
    fams_densified2 <- as.data.frame(table(select(filter(taxonomy_matrix, id%in%rownames(densified2)),level1)))
    family_plot_densified2 <-  ggplot(fams_densified2,aes(x=level1,y=Freq))+
      geom_bar(stat="identity",color="darkblue",fill="darkblue")+
      theme_minimal()+
      theme(panel.grid.major = element_blank(), 
            panel.background = element_blank(), axis.text.x = element_blank(),  # Remove x-axis labels
            axis.line.x = element_blank())+    # Remove x-axis line
      labs(title="Language counts per family, densified dataset (small)",y="language count per family",x="families")+
      ylim(c(0,max(fams_full$Freq)+20))
    
  }else{
    lg_plot_densified <- ggplot(data.frame(density=apply(densified, 1, function(x) (length(na.omit(x))))/ncol(densified)), aes(x=density))+
      geom_histogram(color="darkblue", fill="darkblue", bins=20)+
      theme_minimal()+
      xlim(c(0,1.05))+
      labs(title="Coding density per language, densified dataset", 
           x="coding density per language",
           y="number of languages")
    
    var_plot_densified <- ggplot(data.frame(density=apply(densified, 2, function(x) (length(na.omit(x))))/nrow(densified)), aes(x=density))+
      geom_histogram(color="darkblue", fill="darkblue", bins=20)+
      theme_minimal()+
      theme(axis.line.x = element_blank())+    # Remove x-axis line
      xlim(c(0,1.05))+
      labs(title="Coding density per variable, densified dataset",
           x="coding density per variable",
           y="number of variables")
    
    fams_densified <- as.data.frame(table(select(filter(taxonomy_matrix, id%in%rownames(densified)),level1)))
    family_plot_densified <-  ggplot(fams_densified,aes(x=level1,y=Freq))+
      geom_bar(stat="identity",color="darkblue",fill="darkblue")+
      theme_minimal()+
      theme(panel.grid.major = element_blank(), 
            panel.background = element_blank(), axis.text.x = element_blank(),  # Remove x-axis labels
            axis.line.x = element_blank())+    # Remove x-axis line
      labs(title="Language counts per family, densified dataset",y="language count per family",x="families")+
      ylim(c(0,max(fams_full$Freq)+20))
  }
  
  # arrange the plots side by side
  if(!is.null(densified2)){
    combined_plot <- grid.arrange(lg_plot_full, var_plot_full, family_plot_full,
                                  lg_plot_densified, var_plot_densified, family_plot_densified,
                                  lg_plot_densified2, var_plot_densified2, family_plot_densified2,
                                  ncol = 3)
  }else{
    combined_plot <- grid.arrange(lg_plot_full, var_plot_full, family_plot_full,
                                  lg_plot_densified, var_plot_densified, family_plot_densified,
                                  ncol = 3)
  }
  return(list(comparison,combined_plot))
}

project_data <-  function(
    df,  # a data frame with Longitude, Latitude, and data
    world_map_initial = world_map_initial, # a natural earth map or overlays of several maps
    projection = 8859, # 8859 is the crs code for equal earth WGS84 with 150° as center meridian
    graticules_gap = 20, # gap between single graticule lines in ° (1, 5, 10, 15, 20, 30 are possible) (default 20)
    xmin = -180, # minimum longitude of desired extent (in ° between -180 and 180) (default -180)
    xmax = 180, # maximum longitude of desired extent (in ° between -180 and 180) (default 180)
    ymin = -60, # minimum latitude of desired extent (in ° between -90, 90) (default -60)
    ymax = 85, # maximum latitude of desired extent (in ° between -90, 90) (default 85)
    labels_lat = -175, # position of latitude labels on longitude axis (in ° long) (default -175)
    labels_long = -52.5 # position of latitude labels on latitude axis (in ° lat) (default -52.5)
) {
  library(tidyverse)
  library(dplyr)
  library(sf)
  library(rnaturalearth)
  library(ggplot2)
  
  #### Settings
  # target projection
  df_projection <- 4326 # projection code of the data frame (e.g. 4326 for WGS84) (default 4326)
  # meridian where world map is split up
  split_meridian <- -30
  # deactivating s2 spherical geometry to make following map crops possible
  sf_use_s2(FALSE)
  
  # duplicate base_map for further processing
  base_map <- world_map_initial
  
  # create "split line" to split polygons/linestrings that cross the splitting meridian
  split_line <- st_linestring(x = cbind(split_meridian,c(-90,90)), dim = "XY")
  split_line <- st_geometry(split_line) # makes it possible to assign crs
  st_crs(split_line) <- st_crs(base_map) # assign crs from base map to line
  
  # intersect line with continent polygons/linestrings to identify the ones that cross splitting meridian
  base_map$intersects <- suppressMessages(st_intersects(base_map, split_line, sparse = F))
  base_map_intersects <- filter(base_map, intersects == T) # map with intersecting polygons
  base_map_cleaned <- filter(base_map, intersects == F) # map without intersecting polygons
  
  # crop polygons/linestrings on both sides of splitting meridian separately
  bbox_left <- c(xmin = -180, xmax = split_meridian-0.000001, ymin = -90, ymax = 90)
  bbox_right <- c(xmin = split_meridian+0.000001, xmax = 180, ymin = -90, ymax = 90)
  # 0.000001 ensures that edges of right and left side do not have the exact same coordinates
  base_map_intersects_left <- suppressMessages(suppressWarnings(st_crop(base_map_intersects, bbox_left)))
  base_map_intersects_right <- suppressMessages(suppressWarnings(st_crop(base_map_intersects, bbox_right)))
  
  # combine all three maps
  base_map <- bind_rows(base_map_cleaned, base_map_intersects_left)
  base_map <- bind_rows(base_map, base_map_intersects_right)
  
  # crop map to desired plotting extent
  bbox_map <- c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  base_map <- suppressMessages(suppressWarnings(st_crop(base_map, bbox_map)))
  
  # reproject map
  base_map <- st_transform(base_map, crs = projection)
  
  #### Graticules
  # download graticules from naturalearthdata
  if (!exists(paste0("graticules_", graticules_gap))){
    assign(paste0("graticules_", graticules_gap),
           ne_download(scale = 110,
                       type = paste0("graticules_", graticules_gap),
                       category = "physical"),
           pos = 1)
  }
  
  # load the graticules file from global environment
  graticules <- get(paste0("graticules_", graticules_gap), envir = .GlobalEnv)
  
  # crop graticules to desired plotting extent
  bbox_graticule <- c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  graticules <- suppressMessages(suppressWarnings(st_crop(st_as_sf(graticules), bbox_graticule)))
  
  # reproject graticules
  graticules <- st_transform(graticules, crs = projection)
  
  #### Graticule Labels
  # create data frames with grid labels
  x_labs <- c("-25°", "-20°", "-15°", "-10°", "-5°", "0°", "5°", "10°", "15°", "20°",
              "25°", "30°", "35°", "40°", "45°", "50°", "55°", "60°", "65°", "70°",
              "75°", "80°", "85°", "90°", "95°", "100°", "105°", "110°", "115°", "120°",
              "125°", "130°", "135°", "140°", "145°", "150°", "155°", "160°", "165°",
              "170°", "175°", "180°", "-175°", "-170°", "-165°", "-160°", "-155°", "-150°",
              "-145°", "-140°", "-135°", "-130°", "-125°", "-120°", "-115°", "-110°",
              "-105°", "-100°", "-95°", "-90°", "-85°", "-80°", "-75°", "-70°", "-65°",
              "-60°", "-55°", "-50°", "-45°", "-40°", "-35°")
  y_labs <- c("-85°", "-80°", "-75°", "-70°", "-65°", "-60°", "-55°", "-50°", "-45°",
              "-40°", "-35°", "-30°", "-25°", "-20°", "-15°", "-10°", "-5°", "0°",
              "5°", "10°", "15°", "20°", "25°", "30°", "35°", "40°", "45°", "50°",
              "55°", "60°", "65°", "70°", "75°", "80°", "85°")
  x <- c(seq(-25, 180, by = 5), seq(-175, -35, by = 5))
  long_labels <- data.frame(labs=x_labs, x, y=labels_long)
  y <- c(seq(-85, 0, by = 5), seq(5, 85, by = 5))
  lat_labels <- data.frame(labs=y_labs, x=labels_lat, y=y)
  
  # subset data frames according to graticules gap input variable
  long_labels <- long_labels %>%
    slice(which(x %% graticules_gap == 0))
  lat_labels <- lat_labels %>%
    slice(which(y %% graticules_gap == 0))
  
  # combine data frames into one
  grid_labels <- bind_rows(long_labels, lat_labels)
  
  # create sf object
  grid_labels <- st_as_sf(grid_labels, coords = c("x", "y"), crs = 4326)
  
  # crop graticule labels to desired plotting extent
  bbox_grid <- c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  grid_labels <- suppressMessages(suppressWarnings(st_crop(st_as_sf(grid_labels), bbox_grid)))
  
  # reproject graticule lables
  grid_labels <- st_transform(grid_labels, crs = projection)
  
  #### data
  # change name of columns to "Longitude" and "Latitude"
  if(any(names(df) %in% 'longitude')) {
    df$Longitude <- df$longitude
    df$Latitude <- df$latitude
  }
  
  if(any(names(df) %in% 'lon')) {
    df$Longitude <- df$lon
    df$Latitude <- df$lat
  }
  
  if(any(names(df) %in% 'long')) {
    df$Longitude <- df$long
    df$Latitude <- df$lat
  }
  
  # create an sf object
  df = st_as_sf(df, coords = c("Longitude", "Latitude"), crs = df_projection)
  
  # crop data frame to desired plotting extent
  bbox_data <- c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  df <- suppressMessages(suppressWarnings(st_crop(df, bbox_data)))
  
  # reproject data frame to desired projection
  df <- st_transform(df, crs = projection)
  
  #### ggplot as base plot with base map
  base_plot <- ggplot() +
    geom_sf(data = base_map, colour = "darkgrey", fill = "transparent", size = .25) +
    scale_size(range = c(10, 10)) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.box.background = element_rect(fill = "white", color = "white"),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 6)))
  
  #### objects returned by function
  return(list(base_map = base_map,
              graticules = graticules,
              grid_labels = grid_labels,
              data = df,
              base_plot = base_plot))
}


