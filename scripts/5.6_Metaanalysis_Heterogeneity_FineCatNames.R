# this script attempts to explain the heterogeneity in paths
# focusing on phenological responses to temperature
# by specific trait types and taxon (as categorical variable)

library(sTraitChange)
library(tidyverse)
library(ape)
library(magrittr)

# 1. data read-in  and prepare --------------------------------------------
Coefs_Aut <- readRDS(file = './output_fSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')


## subset to only retain the studies with windur max 52
Coefs_Aut <- Coefs_Aut[Coefs_Aut$WinDur <= 51, ]
length(unique(Coefs_Aut$ID))  ## 202 (from initial 210)


## read in the trait data and merge
traits <- read.csv('./data/speciesTraits.csv')

traits_sub <- subset(traits, select = c(Species, GenLength_y_IUCN))
Coefs_Aut_sp <- merge(Coefs_Aut, traits_sub, by = 'Species', all.x = TRUE)


# read in the phylo tree
vert_tree <- read.tree('./data/phylogenies_100/vert1.tre')
plot(vert_tree)
is.ultrametric(vert_tree)
# get the matrix of phylogenetic correlations
Mat_phylo <- ape::vcv.phylo(vert_tree, corr = TRUE)



# to be able to compare the models prepare the dataset as for the analyses with phylogeny directly
## first update the species names to correspond to the ones on phylogeny
Coefs_Aut_sp <- Coefs_Aut_sp %>%
  mutate(Sp_phylo = case_when(
    Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
    Species == 'Thalasseus sandvicensis' ~ 'Sterna sandvicensis',
    Species == 'Setophaga caerulescens' ~ 'Dendroica caerulescens',
    Species == 'Thalassarche melanophris' ~ 'Thalassarche melanophrys',
    Species == 'Ichthyaetus audouinii' ~ 'Larus audouinii',
    Species == 'Stercorarius maccormicki' ~ 'Catharacta maccormicki',
    TRUE ~ Species))


Coefs_Aut_sp$Species <- unlist(lapply(1:nrow(Coefs_Aut_sp), FUN = function(x){
  binary <- strsplit(as.character(Coefs_Aut_sp$Sp_phylo[x]), " ")
  Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")
}))

# need a copy with the sp names, to be used to account for phylogeny (apart from that
# we also account for among-sp variation by including species as random intercept
# in the model)
Coefs_Aut_sp$Sp_phylo <- Coefs_Aut_sp$Species

# check that all species that are on the tree are also in the dataset
for(i in 1:(nrow(Coefs_Aut_sp))){
  check <- match(Coefs_Aut_sp$Species[i], vert_tree$tip.label)
  if (is.na(check)){
    print(paste(Coefs_Aut_sp$Species[i], 'is not in the phylogeny', sep = ' '))}
}

# other way around
for(i in 1:(nrow(Mat_phylo))){
  check <- match(rownames(Mat_phylo)[i], Coefs_Aut_sp$Species)
  if (is.na(check)){
    print(paste(rownames(Mat_phylo)[i], 'is not in the dataset', sep = ' '))}
}

# correcting for the Oppel's study the entry for weather qauality: it is "2"
Coefs_Aut_sp$WeathQ[Coefs_Aut_sp$WeathQ == ''] <-  "2"
Coefs_Aut_sp %<>% mutate(across(where(is_character), as_factor))



# 2. Heterogeneity due to variety of phenological traits ---------------------

Coefs_phen <- droplevels(Coefs_Aut_sp %>%
  filter(Trait_Categ == 'Phenological'))

unique(Coefs_phen$Trait)
table(Coefs_phen$Trait)

## check the different traits to see whether some of them can be pulled together
unique(Coefs_phen$ID[Coefs_phen$Trait == 'MaleArrivalDate'])   # 186
unique(Coefs_phen$ID[Coefs_phen$Trait == 'FemaleArrivalDate'])  # 185
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'MaleArrivalDate'])   # Moller_1
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'FemaleArrivalDate'])  # Moller_1
unique(Coefs_phen$ID[Coefs_phen$Trait == 'ArrivalDateFemales'])   #  421
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'ArrivalDateFemales'])   # Schaub_et_al
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'ArrivalDateMales'])   # Schaub_et_al
unique(Coefs_phen$ID[Coefs_phen$Trait == 'ArrivalDateMales'])   #  420
# all arrival dates merge under ArrivalDate

unique(Coefs_phen$ID[Coefs_phen$Trait == 'InitiationDate'])   #  430
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'InitiationDate'])   # Oppel_et_al
Coefs_phen[Coefs_phen$Study_Authors == 'Oppel_et_al', ]
Coefs_phen$Trait[Coefs_phen$Study_Authors == 'Oppel_et_al'] <- 'FirstLayDate'

# reclassify into a broader trait type that encompassess similar
# traits (e.g. arrival of females and arrivla of males)
Coefs_phen_type <- Coefs_phen %>%
  mutate(TraitType = as.factor(case_when(Trait %in% c('ArrivalDateFemales',
                                            'ArrivalDateMales',
                                            'FemaleArrivalDate',
                                            'MaleArrivalDate') ~ 'ArrivalDate',
                               Trait %in% c('BreedingDate',
                                            'LayingDateAllBroods',
                                            'MeanBreedingDate',
                                            'LayingDate',
                                            'NestInitiationDate',
                                            'NestDate') ~ 'OnsetBreeding',
                               Trait %in% c('StartOfLaying') ~ 'FirstLayDate',
                               Trait %in% c('AntlerCastDate', 'RutEndDate',
                                            'OestrusDate') ~
                                 'RutDate',
                              .default = as.character(Trait))))

table(Coefs_phen_type$TraitType)
levels(Coefs_phen_type$TraitType)

#Coefs_phen$TraitType <- as.character(Coefs_phen$Trait)

# check
#unique(Coefs_phen$ID[Coefs_phen$Trait == 'NestInitiationDate'])   # 472 471
#unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'NestInitiationDate'])   # Ross_et_al
#Coefs_phen[Coefs_phen$Trait == 'NestInitiationDate', ]  ## bird


# check what is the study with breedingDate
#unique(Coefs_phen$ID[Coefs_phen$Trait == 'BreedingDate'])   #  136
#unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'BreedingDate']) ## Tarwater & Beissinger
# I checked and this is the study by Corey Tarwater on Forpus passerinus, so
# it can be changed to "OnsetBreeding"
#Coefs_phen$Trait[Coefs_phen$Trait == 'BreedingDate'] <- 'OnsetBreeding'


# check LayingDateAllBroods
#unique(Coefs_phen$ID[Coefs_phen$Trait == 'LayingDateAllBroods'])   #  14
#unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'LayingDateAllBroods'])  # Schroeder_et_al
## all broods of the sparrow were included, can potentially be included under the "Onset breeding"


# check MeanBreedingDate
#unique(Coefs_phen$ID[Coefs_phen$Trait == 'MeanBreedingDate'])   #  126
#unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'MeanBreedingDate'])  # Tarwater&Beissinger
## in fact these two studies by Tarwater & Beissinger are simply from different locations;
## so renaming also as "OnsetBreeding"
#Coefs_phen$Trait[Coefs_phen$Trait == 'MeanBreedingDate'] <- 'OnsetBreeding'


# # check FirstLayDate
# unique(Coefs_phen$ID[Coefs_phen$Trait == 'FirstLayDate'])   #   556 559 558 557 560 430
# unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'FirstLayDate'])  # Jenouvrier_et_al Oppel_et_al

# check StartOfLaying
#unique(Coefs_phen$ID[Coefs_phen$Trait == 'StartOfLaying'])   #   486
#unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'StartOfLaying'])  # Kruuk_et_al
# study on the superb fairy wren. I checked and So it is basically firstLaydate (mean of
# first laying dates) CHANGE
#Coefs_phen$Trait[Coefs_phen$Trait == 'StartOfLaying'] <- 'FirstLayDate'

# check Oestrus date - this is only related to rutting, so can be reclassified under
# that type of traits
#unique(Coefs_phen$ID[Coefs_phen$Trait == 'OestrusDate'])   #   312
#unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'OestrusDate'])  # Moyes et al

# check ParturitionDate
#unique(Coefs_phen$ID[Coefs_phen$Trait == 'ParturitionDate'])   #   561 403 402
#unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'ParturitionDate'])  # Renaud_et_al Massot_et_al
# for Ovis canadensis and zootoca vivipara


# Also replace "LayingDate", "NestDate" and "NestInitiationDate" with "OnsetBreeding"
#Coefs_phen$Trait[Coefs_phen$Trait == 'LayingDate'] <- 'OnsetBreeding'
#Coefs_phen$Trait[Coefs_phen$Trait == 'NestInitiationDate'] <- 'OnsetBreeding'
#Coefs_phen$Trait[Coefs_phen$Trait == 'NestDate'] <- 'OnsetBreeding'

# also include "LayingDateAllBroods" under "OnsetBreeding"
#Coefs_phen$Trait[Coefs_phen$Trait == 'LayingDateAllBroods'] <- 'OnsetBreeding'


# 3. Heterogeneity in CZ --------------------------------------------
tre_sub <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% Coefs_phen_type$Species))
Mat_phylo <- ape::vcv.phylo(tre_sub, corr = TRUE)
diag(Mat_phylo)
# check that all species that are on the tree are also in the dataset
for(i in 1:(nrow(Coefs_phen_type))){
  check <- match(Coefs_phen_type$Species[i], tre_sub$tip.label)
  if (is.na(check)){
    print(paste(Coefs_phen_type$Species[i], 'is not in the phylogeny', sep = ' '))}
}

# other way around
for(i in 1:(nrow(Mat_phylo))){
  check <- match(rownames(Mat_phylo)[i], Coefs_phen_type$Species)
  if (is.na(check)){
    print(paste(rownames(Mat_phylo)[i], 'is not in the dataset', sep = ' '))}
}

meta_Phen_CZ_bySpeTrait <- fit_meta_phylo(data_MA = Coefs_phen_type,
                                  Type_EfS = 'Trait_mean<-det_Clim', COV = 'Pvalue',
                                  Cov_fact = 'TraitType',  DD = 'n_effectGR', simpleSEM = TRUE,
                                  A = Mat_phylo,
                                  Trait = FALSE, des.matrix = 'identity')
meta_Phen_CZ_bySpeTrait$data
# a significant effect of the type of measure we use. This is in a way good
# may be good to plot this in a supplement figure, or at least give as a table

CZ_bySpecTrait <- meta_Phen_CZ_bySpeTrait$data[[1]] %>%
  dplyr::select(., Levels_Covar, Estimate, SError,
                EfS_Low, EfS_Upper, Chi2, pval_Covar) %>%
  dplyr::rename(., Variable = Levels_Covar,
                CI.lb = EfS_Low, CI.ub = EfS_Upper,
                pval = pval_Covar) %>%
  dplyr::mutate(dplyr::across(.cols= c(Estimate, SError, dplyr::starts_with('CI')),
                              .fns = ~format(round(.x, 3), nsmall = 3, scientific = FALSE)),
                Chi2 = format(round(Chi2, 2), nsmall = 2, scientific = FALSE))

CZ_bySpecTrait$`p value` <- numeric(length = nrow(CZ_bySpecTrait))
for(i in 1:nrow(CZ_bySpecTrait)){
  if (CZ_bySpecTrait$pval[i] < 0.0001){
    CZ_bySpecTrait$`p value`[i] <- '<0.0001'
  } else {
    CZ_bySpecTrait$`p value`[i] <- format(round(CZ_bySpecTrait$pval[i], 4), nsmall = 4, scientific = FALSE)

  }
}

CZ_bySpecTrait$pval <- NULL

save_xlsx(table = CZ_bySpecTrait, table_name = './tables/Explain_CZHeterogen_SpecificTraits')

# heterogeneity in each path - explain by taxon
meta_phen_HeterTax <- fit_all_meta(data_MA = Coefs_phen,
                              Demog_rate = NULL,
                              Trait_categ = 'Phenological',
                              Clim = 'Temperature',
                              Cov_fact ='Taxon',
                              COV =  'Pvalue',
                              sel = 'Temp_Phen_TaxonExplainHet',
                              A = Mat_phylo,
                              folder_name = './output_all/',
                              colr = c('red', 'black', 'blue'),
                              DD = 'n_effectGR',
                              simpleSEM = TRUE,
                              des.matrix = 'identity',
                              all_Relations = c('Trait_mean<-det_Clim',
                                                'GR<-Trait_mean',
                                                'Ind_GR<-det_Clim',
                                                'GR<-det_Clim'))
meta_phen_HeterTax$meta_res

CZ_byTaxon<- meta_phen_HeterTax$meta_res[[1]] %>%
  dplyr::filter(Relation == 'Trait_mean<-det_Clim') %>%
  dplyr::select(., Levels_Covar, Estimate, SError,
                EfS_Low, EfS_Upper, Chi2, pval_Covar) %>%
  dplyr::rename(., Variable = Levels_Covar,
                CI.lb = EfS_Low, CI.ub = EfS_Upper,
                pval = pval_Covar) %>%
  dplyr::mutate(dplyr::across(.cols= c(Estimate, SError, dplyr::starts_with('CI')),
                              .fns = ~format(round(.x, 3), nsmall = 3, scientific = FALSE)),
                Chi2 = format(round(Chi2, 2), nsmall = 2, scientific = FALSE))

CZ_byTaxon$`p value` <- numeric(length = nrow(CZ_byTaxon))
for(i in 1:nrow(CZ_byTaxon)){
  if (CZ_byTaxon$pval[i] < 0.0001){
    CZ_byTaxon$`p value`[i] <- '<0.0001'
  } else {
    CZ_byTaxon$`p value`[i] <- format(round(CZ_byTaxon$pval[i], 4), nsmall = 4, scientific = FALSE)

  }
}

CZ_byTaxon$pval <- NULL

save_xlsx(table = CZ_byTaxon, table_name = './tables/Explain_CZHeterogen_byTaxon')

# save also for ZG and CG results
meta_Phen_ZG_byTaxon <- fit_meta_phylo(data_MA = Coefs_phen,
                                       Type_EfS = 'GR<-Trait_mean', COV = 'Pvalue',
                                       Cov_fact = 'Taxon',  DD = 'n_effectGR', simpleSEM = TRUE,
                                       A = Mat_phylo,
                                       Trait = FALSE, des.matrix = 'identity')
meta_Phen_ZG_byTaxon$data


ZG_byTaxon <- meta_phen_HeterTax$meta_res[[1]] %>%
  dplyr::filter(Relation == 'GR<-Trait_mean') %>%
  dplyr::select(., Levels_Covar, Estimate, SError,
                EfS_Low, EfS_Upper, Chi2, pval_Covar) %>%
  dplyr::rename(., Variable = Levels_Covar,
                CI.lb = EfS_Low, CI.ub = EfS_Upper,
                pval = pval_Covar) %>%
  dplyr::mutate(dplyr::across(.cols= c(Estimate, SError, dplyr::starts_with('CI')),
                              .fns = ~format(round(.x, 3), nsmall = 3, scientific = FALSE)),
                Chi2 = format(round(Chi2, 2), nsmall = 2, scientific = FALSE))

ZG_byTaxon$`p value` <- numeric(length = nrow(ZG_byTaxon))
for(i in 1:nrow(ZG_byTaxon)){
  if (ZG_byTaxon$pval[i] < 0.0001){
    ZG_byTaxon$`p value`[i] <- '<0.0001'
  } else {
    ZG_byTaxon$`p value`[i] <- format(round(ZG_byTaxon$pval[i], 4), nsmall = 4, scientific = FALSE)

  }
}

ZG_byTaxon$pval <- NULL

save_xlsx(table = ZG_byTaxon, table_name = './tables/Explain_ZGHeterogen_byTaxon')


CG_byTaxon <- meta_phen_HeterTax$meta_res[[1]] %>%
  dplyr::filter(Relation == 'GR<-det_Clim') %>%
  dplyr::select(., Levels_Covar, Estimate, SError,
                EfS_Low, EfS_Upper, Chi2, pval_Covar) %>%
  dplyr::rename(., Variable = Levels_Covar,
                CI.lb = EfS_Low, CI.ub = EfS_Upper,
                pval = pval_Covar) %>%
  dplyr::mutate(dplyr::across(.cols= c(Estimate, SError, dplyr::starts_with('CI')),
                              .fns = ~format(round(.x, 3), nsmall = 3, scientific = FALSE)),
                Chi2 = format(round(Chi2, 2), nsmall = 2, scientific = FALSE))

CG_byTaxon$`p value` <- numeric(length = nrow(CG_byTaxon))
for(i in 1:nrow(CG_byTaxon)){
  if (CG_byTaxon$pval[i] < 0.0001){
    CG_byTaxon$`p value`[i] <- '<0.0001'
  } else {
    CG_byTaxon$`p value`[i] <- format(round(CG_byTaxon$pval[i], 4), nsmall = 4, scientific = FALSE)

  }
}

CG_byTaxon$pval <- NULL
save_xlsx(table = CG_byTaxon, table_name = './tables/Explain_CGHeterogen_byTaxon')


