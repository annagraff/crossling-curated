# code for Fig 7
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(densify)
library(testthat)
library(lsr)
library(gridExtra)

# this function generates expectation_assessment for each expectation
for_viz <- function(ID, modifications, diversity_samples,recoded_data,taxonomy,proportion_languages_must_be_in_applicable_state, title, textpos, conditional, marginal, xlab){
  
  expectation <- filter(modifications, modification.ID == ID)
  
  raw_data <- recoded_data %>% select(c(expectation$feature.1.for.test,expectation$feature.2.for.test,glottocode)) # subset data to relevant variables, including ? and NA
  
  full_tbl_withqs <- table(unlist(raw_data[,1]),unlist(raw_data[,2])) # cross-tabulate variables using full data (including ? and NA)
  
  # log characteristics about language overlap: for certain
  # for certain considerations,  we do not care about "?" and "NA" --> convert them to NA
  non_NAq_data <- raw_data
  non_NAq_data[non_NAq_data =="?"] <- NA
  non_NAq_data[non_NAq_data =="NA"] <- NA
  
  # log how many languages each variable is coded for (nrlgs1 for variable 1 and nrlgs2 for variable 2; these values are logged)
  nrlgs1 <- length(na.omit(non_NAq_data[,1]))
  nrlgs2 <- length(na.omit(non_NAq_data[,2]))
  
  # determine the number of languages coded for both variables (this value is logged)
  absoverlap <- non_NAq_data %>% filter(!is.na(get(expectation$feature.1.for.test))) %>% filter(!is.na(get(expectation$feature.2.for.test))) %>% nrow()
  
  # log the number of languages in the smaller and the number of languages in the larger variable; then determine the relative overlap of languages for the smaller and larger variable (these values are logged)
  smaller <- sort(c(nrlgs1,nrlgs2))[1]
  larger <- sort(c(nrlgs1,nrlgs2))[2]
  relative.overlap.fsmall.in.flarge <- absoverlap/smaller
  relative.overlap.flarge.in.fsmall <- absoverlap/larger
  
  # convert the feature values into factor format
  non_NAq_data[,expectation$feature.1.for.test] <- as.factor(unlist(non_NAq_data[,expectation$feature.1.for.test]))
  non_NAq_data[,expectation$feature.2.for.test] <- as.factor(unlist(non_NAq_data[,expectation$feature.2.for.test]))
  
  raw_data[,expectation$feature.1.for.test] <- as.factor(unlist(raw_data[,expectation$feature.1.for.test]))
  raw_data[,expectation$feature.2.for.test] <- as.factor(unlist(raw_data[,expectation$feature.2.for.test]))
  
  # prepare a table to log relevant characteristics of the crosstable for each of the diversity samples
  expectation_assessment<-slice(data.frame(proportion_overlap_applicable_vs_non_applicable=NA, # this logs for each sample the proportion of languages coded for v2 in the relevant state of v1
                                           overlap_sufficient_test_power_positive=NA, # this denotes for each sample whether the proportion and number of languages coded for v2 in the relevant state of v1 is sufficient, the proportion is defined by the user and denoted proportion_languages_must_be_in_applicable_state
                                           result1_OR_AND_THEN=NA, # this denotes the proportion of languages behaving against the expectation in the THEN condition and in the first dimension of the OR and AND conditions
                                           result2_OR_AND=NA, # this denotes the proportion of languages behaving against the expectation in the second dimension of the OR and AND conditions
                                           baseline1_OR_AND_THEN=NA, # this denotes the proportion of languages in the "unequals" state for variable 2 in the THEN condition overall or in the first dimension of the OR condition
                                           baseline2_OR_AND=NA),0) # this denotes whether the expected relationship is confirmed in the subsample in accordance with the defined thresholds
  
  for (ds in 1:ncol(diversity_samples)){ # loop through all diversity samples
    # select the subset of languages in the diversity sample for the relevant variables and tabulate for both the full data (including ? and NA) and the non-Q/NA data.
    expectation_ds_sample_withqs <- raw_data %>% filter(glottocode%in%diversity_samples[,ds])
    sample_table_withqs <- table(unlist(expectation_ds_sample_withqs[,1]),unlist(expectation_ds_sample_withqs[,2])) # this is the table
    
    expectation_ds_sample_noqs <- non_NAq_data %>% filter(glottocode%in%diversity_samples[,ds])
    sample_table_noqs <- table(unlist(expectation_ds_sample_noqs[,1]),unlist(expectation_ds_sample_noqs[,2])) # this is the table
    
    # sanity checks
    expect_true(all(levels(as.data.frame(sample_table_noqs)$Var1) %in% c(unlist(strsplit(as.character(unlist(expectation$f1.equals.for.test)),", ")), unlist(strsplit(as.character(unlist(expectation$f1.unequals.for.test)),", ")))))
    expect_true(all(c(unlist(strsplit(as.character(unlist(expectation$f1.equals.for.test)),", ")), unlist(strsplit(as.character(unlist(expectation$f1.unequals.for.test)),", "))) %in% levels(as.data.frame(sample_table_noqs)$Var1)))
    expect_true(all(c(unlist(strsplit(as.character(unlist(expectation$f2.equals.for.test)),", ")), unlist(strsplit(as.character(unlist(expectation$f2.unequals.for.test)),", "))) %in% levels(as.data.frame(sample_table_noqs)$Var2)))
    expect_true(all(levels(as.data.frame(sample_table_noqs)$Var2) %in% c(unlist(strsplit(as.character(unlist(expectation$f2.equals.for.test)),", ")), unlist(strsplit(as.character(unlist(expectation$f2.unequals.for.test)),", ")))))
    
    # in the THEN-case, we evaluate how many languages are coded for the second variable, given the relevant state of the first variable (applicable), and how many languages are coded for the relevant state of the first variable overall (all)
    applicable <- sum(sample_table_withqs[unlist(strsplit(as.character(unlist(expectation$f1.equals.for.test)),", ")),c(unlist(strsplit(as.character(unlist(expectation$f2.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$f2.unequals.for.test)),", ")))])
    all <- sum(sample_table_withqs[unlist(strsplit(as.character(unlist(expectation$f1.equals.for.test)),", ")),])
    
    # this yields the relevant proportion, which is logged
    relevant.proportion <- applicable/all
    expectation_assessment[ds,"proportion_overlap_applicable_vs_non_applicable"] <- relevant.proportion
    expectation_assessment[ds,"overlap_sufficient_test_power_positive"] <- relevant.proportion >= proportion_languages_must_be_in_applicable_state
    
    # in the THEN-condition, the languages that behave according to the expectation are those with the "equals" state of both variable 1 and variable 2
    hyp.applies <- sample_table_noqs[unlist(strsplit(as.character(unlist(expectation$f1.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$f2.equals.for.test)),", "))]
    
    # the languages that behave against the expectation are those with the "equals" state of variable 1 and the "unequals" state of variable 2
    hyp.NOT.applies <- sample_table_noqs[unlist(strsplit(as.character(unlist(expectation$f1.equals.for.test)),", ")),unlist(strsplit(as.character(unlist(expectation$f2.unequals.for.test)),", "))]
    
    # track the proportion of cases, in which the expectation is violated (d1)
    d1 <- sum(hyp.NOT.applies)/(sum(hyp.applies)+sum(hyp.NOT.applies)) # proportion of cases, in which gut expectation is violated with respect to variable 1
    expectation_assessment[ds,"result1_OR_AND_THEN"]<-d1
    
    # the proportion must be compared to the baseline expectation: this is the proportion of languages in the "unequals" state of variable 2 overall (independently of variable 1 being coded or not, and which state it would be in)
    baseline_v2_equals <- expectation_ds_sample_withqs %>% filter(get(expectation$feature.2.for.test)%in%unlist(strsplit(as.character(unlist(expectation$f2.equals.for.test)),", "))) %>% nrow()
    baseline_v2_unequals <- expectation_ds_sample_withqs %>% filter(get(expectation$feature.2.for.test)%in%unlist(strsplit(as.character(unlist(expectation$f2.unequals.for.test)),", "))) %>% nrow()
    baseline_v2 <- baseline_v2_unequals/(baseline_v2_unequals+baseline_v2_equals)
    
    # log the baseline
    expectation_assessment[ds,"baseline1_OR_AND_THEN"] <- baseline_v2
  }
  
  dd <- data.frame(Proportion=c(expectation_assessment$result1_OR_AND_THEN,expectation_assessment$baseline1_OR_AND_THEN), Probability = c(rep(conditional,nrow(expectation_assessment)),rep(marginal,nrow(expectation_assessment)))) %>% na.omit()
  
  p <- ggplot(dd, aes(x = Proportion, fill = Probability))+
    geom_density(alpha = 0.8, na.rm = TRUE,  aes(x=Proportion,fill=Probability))+
    scale_fill_manual( values = c("#FFDB6D","#4E84C4"))+
    geom_vline(xintercept = mean(expectation_assessment$result1_OR_AND_THEN,na.rm=T)+sd(expectation_assessment$result1_OR_AND_THEN,na.rm=T), linetype="dotted") +
    geom_vline(xintercept = 0.2, linetype="dashed") +
    annotate("text", x=0.2, y=textpos, label="P = 0.2", 
             size=3, color="black") +
    annotate("text", x=mean(expectation_assessment$result1_OR_AND_THEN,na.rm=T)+sd(expectation_assessment$result1_OR_AND_THEN,na.rm=T), y=textpos, label="µ+σ", 
             size=3, color="black") +
    
    xlim(c(0,1)) +
    labs(x=xlab,
         y="Density",
         title = title, 
         subtitle = paste("µ + σ = ", 
                          round(mean(expectation_assessment$result1_OR_AND_THEN,na.rm=T)+sd(expectation_assessment$result1_OR_AND_THEN,na.rm=T),2),
                          ", Cohen's D = ",round(cohensD(na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==TRUE)$result1_OR_AND_THEN),
                                                         na.omit(filter(expectation_assessment,overlap_sufficient_test_power_positive==TRUE)$baseline1_OR_AND_THEN), 
                                                         method="paired"), digits = 1),collapse = "", sep = ""))+
    theme_bw()
  return(p)
}


modifications <- read.csv("curated_data/TLI/statisticalTLI/cldf/modifications.csv")
diversity_samples <- read.csv("curated_data/TLI/1000_diversity_samples.csv", row.names = 1)
recoded_data <- read.csv("curated_data/TLI/statisticalTLI/data_for_stats.csv", row.names = "X")
names(recoded_data) <- gsub("\\.","+",names(recoded_data))
taxonomy <- as_flat_taxonomy_matrix(glottolog_languoids)
proportion_languages_must_be_in_applicable_state <- 1/3

swoo23a <- for_viz(ID = "S-WOO-23a", 
                   modifications = modifications,
                   diversity_samples = diversity_samples,
                   recoded_data = recoded_data,
                   taxonomy = taxonomy,
                   proportion_languages_must_be_in_applicable_state = proportion_languages_must_be_in_applicable_state,
                   title = "'if prepositions, then word order noun-genitive'",
                   conditional = "conditional: P(prepositions | genitive-noun)",
                   marginal = "marginal: P(genitive-noun)",
                   xlab = "",
                   textpos = 8.25)

plot(swoo23a)

swoo23b <- for_viz(ID = "S-WOO-23b", 
                   modifications = modifications,
                   diversity_samples = diversity_samples,
                   recoded_data = recoded_data,
                   taxonomy = taxonomy,
                   proportion_languages_must_be_in_applicable_state = proportion_languages_must_be_in_applicable_state,
                   title = "'if postpositions, then word order genitive-noun'",
                   conditional = "conditional: P(postpositions | noun-genitive)",
                   marginal = "marginal: P(noun-genitive)",
                   xlab = "",
                   textpos = 32)

plot(swoo23b)

smor10 <- for_viz(ID = "S-MOR-10", 
                  modifications = modifications,
                  diversity_samples = diversity_samples,
                  recoded_data = recoded_data,
                  taxonomy = taxonomy,
                  proportion_languages_must_be_in_applicable_state = proportion_languages_must_be_in_applicable_state,
                  title = "'if verb agreement, then tense-aspect inflection'",
                  conditional = "conditional: P(agreement | no TA-inflection)",
                  marginal = "marginal: P(no TA-inflection)",
                  xlab = "Probability distributions across 1000 samples",
                  textpos = 14)

plot(smor10)

# arrange the plots side by side, print and save
combined_plot_full <- grid.arrange(swoo23b,swoo23a,smor10, ncol = 1)
ggsave("plots/Fig7 associations.png", plot = combined_plot_full, width = 14, height = 8, dpi = 700)
