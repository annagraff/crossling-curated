#### This script compiles external data from WALS, AUTOTYP, Lexibank and PHOIBLE to generate TypLinkInd

rm(list=ls())

#load relevant libraries
library(tidyverse)
library(testthat)
library(GetoptLong)
library(readr)
library(densify)

# create glottolog taxonomy
taxonomy <- as_flat_taxonomy_matrix(glottolog_languoids)

######### collect and compile typological data from WALS ######### 
temp <- tempfile()
download.file("https://zenodo.org/records/7385533/files/cldf-datasets/wals-v2020.3.zip?download=1>cldf-datasets/wals-v2020.3.zip",temp)
wals_languages <- read.csv(unz(temp, "cldf-datasets-wals-878ea47/cldf/languages.csv"))
wals_values <- read.csv(unz(temp, "cldf-datasets-wals-878ea47/cldf/values.csv"))
wals_parameters <- read.csv(unz(temp, "cldf-datasets-wals-878ea47/cldf/parameters.csv"))
wals_codes <- read.csv(unz(temp, "cldf-datasets-wals-878ea47/cldf/codes.csv"))
unlink(temp)
rm(temp)

wals_parameters$Feature <- apply(wals_parameters, 1, function(x) paste("WALS_",x[1]," ",x[2],sep=""))
wals_parameters <- wals_parameters %>% select(c("ID","Feature"))
wals_codes <- wals_codes %>% select(c("ID","Name"))

names(wals_languages)[1]<-"Language_ID"
names(wals_parameters)[1]<-"Parameter_ID"
names(wals_codes)<-c("Code_ID","Code_Name")

wals <- wals_languages %>% merge(wals_values, by = "Language_ID", all = TRUE) %>%
  merge(wals_parameters, by = "Parameter_ID", all = TRUE) %>%
  merge(wals_codes, by = "Code_ID", all = TRUE) %>%
  dplyr::filter(Glottocode!="") %>%
  select(c("Glottocode","Feature","Code_Name")) %>% na.omit()

wals_data <- as.data.frame(pivot_wider(wals,
                                  names_from = Feature,
                                  values_from = Code_Name,
                                  values_fill = "?",
                                  values_fn = function(x){if (length(unique(x))==1){unique(x)} else ("?")}))

# sanity checks
expect_false(any(duplicated(wals_data$Glottocode)))

######### collect and compile typological data from autotyp ######### 
autotyp_languages <- read.csv("../../raw_data/autotyp_v.1.1.1/languages.csv",header=TRUE, sep=",")
autotyp_parameters <- read.csv("../../raw_data/autotyp_v.1.1.1/parameters.csv",header=TRUE, sep=",")
autotyp_parameters <- autotyp_parameters %>% filter(dataset %in% c("Alienability",
                                                                   "AlignmentForDefaultPredicatesPerLanguage",
                                                                   "ClauseWordOrder",
                                                                   "Clusivity",
                                                                   "Gender",
                                                                   "GrammaticalMarkersPerLanguage",
                                                                   "LocusOfMarkingPerLanguage",
                                                                   "MaximallyInflectedVerbAgreementAggregatedByMarkerPositionBinned4",
                                                                   "MaximallyInflectedVerbSynthesis", 
                                                                   "MaximallyInflectedVerbInflectionAndAgreementCountsByPosition",
                                                                   "MaximallyInflectedVerbInflectionCategoriesAggregatedPresence",
                                                                   "MorphologyPerLanguage",
                                                                   "NPStructurePerLanguage",
                                                                   "NumeralClassifiers",
                                                                   "PredicateClassesSemanticsPerLanguage"))

autotyp_parameters$Feature <- apply(autotyp_parameters, 1, function(x) paste("autotyp",x[11],x[2],sep="_"))
autotyp_parameters <- autotyp_parameters %>% select(c("ID","Feature"))

names(autotyp_languages)[1]<-"Language_ID"
names(autotyp_parameters)[1]<-"Parameter_ID"

autotyp_values <- read.csv("../../raw_data/autotyp_v.1.1.1/values.csv",header=TRUE, sep=",")
autotyp_values <- autotyp_values %>% filter(Parameter_ID%in%unique(autotyp_parameters$Parameter_ID))

autotyp_merged <- autotyp_languages %>% merge(autotyp_values, by = "Language_ID", all = TRUE) %>% 
  merge(autotyp_parameters, by = "Parameter_ID", all = TRUE) %>% 
  filter(Glottocode!="") %>% 
  select(c("Glottocode","Feature","Value")) %>%
  na.omit()

autotyp_data<-pivot_wider(autotyp_merged, 
                       names_from = Feature,
                       values_from = Value, 
                       values_fill = "?", 
                       values_fn = function(x){if (length(unique(x))==1){unique(x)} else ("?")})

######### collect and compile typological data from lexibank ######### 
lexibank_lgs <- read.csv("../../raw_data/lexibank_v.1.0/languages.csv",header=TRUE, sep=",")
lexibank_lexicon <- read.csv("../../raw_data/lexibank_v.1.0/lexicon-values.csv",header=TRUE, sep=",")
lexibank_phonology <- read.csv("../../raw_data/lexibank_v.1.0/phonology-values.csv",header=TRUE, sep=",")

lexibank_lexicon <- lexibank_lexicon %>% filter(Parameter_ID%in%c("concepts","forms","senses")==F) %>% select(Language_ID,Value,Parameter_ID)
lexibank_phonology <- lexibank_phonology %>% filter(Parameter_ID%in%c("concepts","forms","senses","forms_with_sounds")==F) %>% select(Language_ID,Value,Parameter_ID)
lexibank_both <- rbind(lexibank_lexicon,lexibank_phonology)

names(lexibank_lgs)[1]<-"Language_ID"

lexibank_merged <- lexibank_both %>% merge(lexibank_lgs,by="Language_ID") %>% select(Glottocode,Value,Parameter_ID)
lexibank_merged$Parameter_ID<-paste("lexibank_",lexibank_merged$Parameter_ID,sep="")

lexibank_data<-pivot_wider(lexibank_merged, 
                           names_from = Parameter_ID,
                           values_from = Value, 
                           values_fill = "?", 
                           values_fn = function(x){if (length(unique(x))==1){unique(x)} else ("?")})

######### collect and compile typological data from PHOIBLE ######### 

phoible_data <- local({
  # load the PHOIBLE phoneme list, only take the privileged sources
  phoible <- read_csv("../../raw_data/phoible_v.20230416.csv", # url("https://raw.githubusercontent.com/phoible/dev/master/data/phoible.csv") after the last commit as of January 2024 (https://github.com/phoible/dev/commit/7030ae02863f0e1ddaf67f0f950c0ea1477cd4ee)
                      col_types=c(InventoryID='i', 
                                  Marginal='l', 
                                  .default='c')) %>%
    filter(!is.na(Glottocode)) %>%
    # select features
    select(
      InventoryID,
      Glottocode,
      # list of all distinctive features
      tone,
      stress,
      syllabic,
      short,
      long,
      consonantal,
      sonorant,
      continuant,
      delayedRelease,
      approximant,
      tap,
      trill,
      nasal,
      lateral,
      labial,
      round,
      labiodental,
      coronal,
      anterior,
      distributed,
      strident,
      dorsal,
      high,
      low,
      front,
      back,
      tense,
      retractedTongueRoot, 
      advancedTongueRoot,
      periodicGlottalSource,
      epilaryngealSource,
      spreadGlottis,
      constrictedGlottis,
      fortis,
      lenis,
      raisedLarynxEjective,
      loweredLarynxImplosive,
      click
    )
  
  # if there is more than one inventory per glottocode, select the largest
  counts <- phoible %>% 
    group_by(InventoryID, Glottocode) %>% 
    summarize(phonemes = n())
  
  biggest_inventories_by_glottocode <- counts %>% 
    group_by(Glottocode) %>% 
    slice_max(order_by = phonemes, n = 1)
  
  phoible <- phoible %>% 
    filter(InventoryID %in% biggest_inventories_by_glottocode$InventoryID) %>%
    select(-InventoryID)
  
  # binarized aggregated features based on distinctive features
  phoible_aggregated_features <- transmute(phoible,
                                           PHOIBLE_BAF_affricates = (continuant %in% '-' & delayedRelease %in% '+' & click %in% '-' & nasal %in% '-') | (continuant %in% '-,+' & delayedRelease %in% '-,+' & click %in% '-' & coronal %in% '+' & nasal %in% '-') | (continuant %in% '-,+' & delayedRelease %in% c('+','-,+') & click %in% '-' & dorsal %in% '+' & nasal %in% '-'), 
                                           PHOIBLE_BAF_alveolarApproximant = consonantal %in% '-' & trill %in% '-' & lateral %in% '-' & coronal %in% '+' & dorsal %in% '-' & strident %in% '-', 
                                           PHOIBLE_BAF_aspirated = spreadGlottis %in% c('+','-,+') & consonantal %in% '+' & periodicGlottalSource %in% '-',
                                           PHOIBLE_BAF_aspiratedFricatives = sonorant %in% '-' & continuant %in% '+' & spreadGlottis %in% '+' & consonantal %in% '+' & periodicGlottalSource %in% '-',
                                           PHOIBLE_BAF_aspiratedStops = consonantal %in% '+' & continuant %in% '-' & sonorant %in% '-' & spreadGlottis %in% c('+','-,+') & periodicGlottalSource %in% '-',
                                           PHOIBLE_BAF_bilabials = consonantal %in% '+' & labial %in% '+' & labiodental %in% '-',
                                           PHOIBLE_BAF_breathyVowels = spreadGlottis %in% '+' & periodicGlottalSource %in% '+' & consonantal %in% c('-','0') & syllabic %in% '+',
                                           PHOIBLE_BAF_breathyConsonants = spreadGlottis %in% '+' & periodicGlottalSource %in% '+' & consonantal %in% '+' & syllabic %in% '-',
                                           PHOIBLE_BAF_clicks = click %in% c('+','-,+','+,-'),
                                           PHOIBLE_BAF_coronalTrills = trill %in% c('+', '-,+', '-,-,+', '+,-') & coronal %in% '+',
                                           PHOIBLE_BAF_creakyVowels = constrictedGlottis %in% '+' & periodicGlottalSource %in% '+' & consonantal %in% c('-','0') & syllabic %in% '+',
                                           PHOIBLE_BAF_creakyConsonants = constrictedGlottis %in% '+' & periodicGlottalSource %in% '+' & consonantal %in% '+' & syllabic %in% '-',
                                           PHOIBLE_BAF_dentalAlveolarVoicedLateralApproximant = lateral %in% c('+','-,-,+','-,+','-,+,-','+,-','0,-,+') & approximant %in% c('+','+,-','-,+','-,-,+') & dorsal %in% c('-','-,+','+,-') & anterior %in% '+',
                                           PHOIBLE_BAF_fricatives = sonorant %in% '-' & continuant %in% '+' & consonantal %in% '+',
                                           PHOIBLE_BAF_frontRoundVowels = syllabic %in% '+' & consonantal %in% '-' & front %in% '+' & round %in% '+',
                                           PHOIBLE_BAF_glides = syllabic %in% '-' & consonantal %in% '-' & approximant %in% '+' & continuant %in% '+' & strident %in% '0',
                                           PHOIBLE_BAF_glottalizedResonants = consonantal %in% '+' & sonorant %in% '+' & constrictedGlottis %in% c('+','+,-'), 
                                           PHOIBLE_BAF_interdentals = (consonantal %in% '+' & continuant %in% '+' & tap %in% '-' & trill %in% '-' & lateral %in% '-' & coronal %in% '+' & anterior %in% '+' & distributed %in% '+' & strident %in% '-' & dorsal %in% '-') | (consonantal %in% '+' & continuant %in% '-' & sonorant %in% '-' & delayedRelease %in% '+' & coronal %in% '+' & anterior %in% '+' & distributed %in% '+' & strident %in% '-' & dorsal %in% '-'),
                                           PHOIBLE_BAF_labializedConsonants = consonantal %in% '+' & labial %in% c('+', '-,+') & round %in% c('+','-,+','-,-,+'),
                                           PHOIBLE_BAF_labialVelars = consonantal %in% '+' & labial %in% '+' & dorsal %in% '+' & high %in% '+' & click %in% '-',
                                           PHOIBLE_BAF_labiodentals = labiodental %in% c('+','-,+','+,-','+,+,-'),
                                           PHOIBLE_BAF_lateralObstruents = sonorant %in% '-' & lateral %in% c('+','-,+') & continuant %in% c('+','-,+') & click %in% '-',
                                           PHOIBLE_BAF_laterals = lateral %in% c('+','-,-,+','-,+','-,+,-','+,-','0,-,+'),
                                           PHOIBLE_BAF_liquids = consonantal %in% '+' & approximant %in% '+',
                                           PHOIBLE_BAF_longConsonants = consonantal %in% '+' & long %in% c('+','+,-','-,+'), 
                                           PHOIBLE_BAF_longVowels = syllabic %in% '+' & consonantal %in% '-' & long %in% c('+','+,-','-,+'), 
                                           PHOIBLE_BAF_loweredLarynxImplosives = loweredLarynxImplosive %in% c('+','-,+','+,-'),
                                           PHOIBLE_BAF_nasalizedVowels = syllabic %in% '+' & consonantal %in% '-' & nasal %in% c('+', '-,+', '+,-'),
                                           PHOIBLE_BAF_nasals = sonorant %in% '+' & approximant %in% '-',
                                           PHOIBLE_BAF_palatalAndPalatalizedObstruents = consonantal %in% '+' & high %in% c('+','-,+') & front %in% c('+','-,+') & sonorant %in% '-',
                                           PHOIBLE_BAF_palatalAndPalatalizedSonorants = consonantal %in% '+' & high %in% '+' & front %in% c('+','-,+') & sonorant %in% '+',
                                           PHOIBLE_BAF_palatalLaterals = consonantal %in% '+' & high %in% '+' & front %in% '+' & sonorant %in% '+' & continuant %in% '+' & tap %in% '-' & lateral %in% '+',
                                           PHOIBLE_BAF_pharyngealsEpiglottals = (low %in% c('+','+,-') & syllabic %in% '-') | (epilaryngealSource %in% '+' & syllabic %in% '-') | (retractedTongueRoot %in% c('+','-,+','-,-,+') & syllabic %in% '-'),
                                           PHOIBLE_BAF_raisedLarynxEjectives = raisedLarynxEjective %in% c('+','-,+','-,-,+'),
                                           PHOIBLE_BAF_retroflexConsonants = coronal %in% '+' & anterior %in% '-' & distributed %in% '-' & click %in% '-',
                                           PHOIBLE_BAF_stridents = strident %in% c('+','-,-,+','-,+','-,+,-','+,-','0,-,+','0,0,-,+'),
                                           PHOIBLE_BAF_taps = tap %in% c('+','-,+','-,-,+'),
                                           PHOIBLE_BAF_tone = tone %in% '+',
                                           PHOIBLE_BAF_uvularTrills = trill %in% c('+', '-,+', '-,-,+', '+,') & coronal %in% '-' & labial %in% '-',
                                           PHOIBLE_BAF_uvulars = consonantal %in% '+' & coronal %in% '-' & dorsal %in% '+' & low %in% '-' & back %in% c('+','+,-','-,+'),
                                           PHOIBLE_BAF_velarNasals = dorsal %in% '+' & nasal %in% '+' & consonantal %in% '+' & coronal %in% '-' & sonorant %in% '+' & back %in% '-' & click %in% '-'
  )
  
  
  phoible <- phoible %>% select(Glottocode) %>%
    cbind(phoible_aggregated_features)
  
  # compute presence per glottocode
  phoible_data <- phoible %>% 
    group_by(Glottocode) %>%
    summarize_all(list(
      presence = ~any(.)
    ))
})

######### combine all data sets into one data frame ######### 

typology_merged <- wals_data %>% 
  full_join(autotyp_data, by= "Glottocode") %>% 
  full_join(lexibank_data, by= "Glottocode") %>% 
  full_join(phoible_data, by= "Glottocode") %>%
  filter(!is.na(Glottocode))

typology_merged <- apply(typology_merged,2,function(x) as.character(x))

typology_merged[is.na(typology_merged)]<-"?"

colnames(typology_merged)[1] <- "glottocode"

write.csv(typology_merged,"../../curated_data/TypLinkInd/compiled_external_input_features.csv",row.names = F)
