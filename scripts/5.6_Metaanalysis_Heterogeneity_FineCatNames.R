# this script attempts to explain the heterogeneity in paths
# focusing on phenological responses to temperature
# by specific trait types

library(ggpubr)
library(metafor)
library(patchwork)
library(multcomp)
library(tidyverse)
library(magrittr)

# 1. data read-in  and prepare --------------------------------------------
Coefs_Aut <- readRDS(file = './output_fSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')


## subset to only retain the studies with windur max 52
Coefs_Aut <- Coefs_Aut[Coefs_Aut$WinDur <= 51, ]
length(unique(Coefs_Aut$ID))  ## 202 (from initial 210)


## read in the trait data and merge
traits <- read.csv('./data-raw/speciesTraits.csv')

## some preparation of the trait data prior to use
## adding two alternative diet classif + fixing the migratory mode
traits$Diet.5Cat_Elton <- trimws(traits$Diet.5Cat_Elton, 'both')
traits_proc <- traits %>%
  dplyr::mutate(., DietPS = dplyr::case_when(
    Diet.5Cat_Elton %in% c('Herbivore', 'PlantSeed') ~ 'primary consumer',
    Diet.5Cat_Elton %in% c('Invertebrate', 'Insectivore') ~ 'secondary consumer, invert',
    Diet.5Cat_Elton %in% c('Carnivore', 'VertFish', 'VertFishScav') ~ 'secondary consumer, vert',
    Diet.5Cat_Elton %in% c('Omnivore', 'VertInvertEggs', 'InvertFish') ~ 'secondary consumer, omnivore'),
    Diet_HCO = dplyr::case_when(
      Diet.5Cat_Elton %in% c('Herbivore', 'PlantSeed') ~ 'herbivore',
      Diet.5Cat_Elton %in% c('Omnivore') ~ 'omnivore',
      Diet.5Cat_Elton %in% c('Invertebrate', 'Insectivore', 'Carnivore',
                             'VertFish', 'VertFishScav', 'VertInvertEggs', 'InvertFish') ~ 'carnivore'),
    Migrat = dplyr::case_when(
      Migratory.mode_Sibly %in% c('nonmigrant', 'resident', 'nonmigrant???') ~ 'resident',
      Migratory.mode_Sibly %in% c('migrant') ~ 'migrant',
      TRUE ~ 'unknown'))

traits_sub <- subset(traits_proc, select = c(Species, GenLength_y_IUCN, concern_IUCN, reprod_rate))
Coefs_Aut_sp <- merge(Coefs_Aut, traits_sub, by = 'Species', all.x = TRUE)


# read in the phylo tree
vert_tree <- read.tree('./data/phylogenies/vert1.tre')
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
# after checking the original papers:
# For 185 & 186 (studies by Moeller on Hirundo rustica): mean arival date of females and Males;
# for studies by Schaub et al: it is the First arival date for males and females, so I will not merge them.
# but I will rename the traits from Shcaub et al. as "FirstArrivalDateMale..."
Coefs_phen$Trait <- as.character(Coefs_phen$Trait)
Coefs_phen$Trait[Coefs_phen$Trait == 'ArrivalDateFemales'] <- 'FirstArrivalDateFemales'
Coefs_phen$Trait[Coefs_phen$Trait == 'ArrivalDateMales'] <- 'FirstArrivalDateMales'


unique(Coefs_phen$ID[Coefs_phen$Trait == 'NestInitiationDate'])   # 472 471
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'NestInitiationDate'])   # Ross_et_al
Coefs_phen[Coefs_phen$Trait == 'NestInitiationDate', ]  ## bird


unique(Coefs_phen$ID[Coefs_phen$Trait == 'InitiationDate'])   #  430
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'InitiationDate'])   # Oppel_et_al
Coefs_phen[Coefs_phen$Study_Authors == 'Oppel_et_al', ]
Coefs_phen$Trait[Coefs_phen$Study_Authors == 'Oppel_et_al'] <- 'FirstLayDate'




# check what is the study with breedingDate
unique(Coefs_phen$ID[Coefs_phen$Trait == 'BreedingDate'])   #  136
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'BreedingDate']) ## Tarwater & Beissinger
# I checked and this is the study by Corey Tarwater on Forpus passerinus, so
# it can be changed to "OnsetBreeding"
Coefs_phen$Trait[Coefs_phen$Trait == 'BreedingDate'] <- 'OnsetBreeding'


# check LayingDateAllBroods
unique(Coefs_phen$ID[Coefs_phen$Trait == 'LayingDateAllBroods'])   #  14
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'LayingDateAllBroods'])  # Schroeder_et_al
## keep as is, because here all broods of the sparrow were included...


# check MeanBreedingDate
unique(Coefs_phen$ID[Coefs_phen$Trait == 'MeanBreedingDate'])   #  126
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'MeanBreedingDate'])  # Tarwater&Beissinger
## in fact these two studies by Tarwater & Beissinger are simply from different locations;
## so renaming also as "OnsetBreeding"
Coefs_phen$Trait[Coefs_phen$Trait == 'MeanBreedingDate'] <- 'OnsetBreeding'

# check OnsetBreeding
unique(Coefs_phen$ID[Coefs_phen$Trait == 'OnsetBreeding'])   #  443
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'OnsetBreeding'])  # Rodel et al
# it is a rabbit, so makes sense to keep OnsetBreeding

# check FirstLayDate
unique(Coefs_phen$ID[Coefs_phen$Trait == 'FirstLayDate'])   #   556 559 558 557 560 430
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'FirstLayDate'])  # Jenouvrier_et_al Oppel_et_al

# check StartOfLaying
unique(Coefs_phen$ID[Coefs_phen$Trait == 'StartOfLaying'])   #   486
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'StartOfLaying'])  # Kruuk_et_al
# study on the superb fairy wren. I checked and So it is basically firstLaydate (mean of
# first laying dates) CHANGE
Coefs_phen$Trait[Coefs_phen$Trait == 'StartOfLaying'] <- 'FirstLayDate'

# check Oestrus date
unique(Coefs_phen$ID[Coefs_phen$Trait == 'OestrusDate'])   #   312
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'OestrusDate'])  # Moyes et al

# check ParturitionDate
unique(Coefs_phen$ID[Coefs_phen$Trait == 'ParturitionDate'])   #   561 403 402
unique(Coefs_phen$Study_Authors[Coefs_phen$Trait == 'ParturitionDate'])  # Renaud_et_al Massot_et_al
# for Ovis canadensis and zootoca vivipara


# Also replace "LayingDate", "NestDate" and "NestInitiationDate" with "OnsetBreeding"
Coefs_phen$Trait[Coefs_phen$Trait == 'LayingDate'] <- 'OnsetBreeding'
Coefs_phen$Trait[Coefs_phen$Trait == 'NestInitiationDate'] <- 'OnsetBreeding'
Coefs_phen$Trait[Coefs_phen$Trait == 'NestDate'] <- 'OnsetBreeding'

# also include "LayingDateAllBroods" under "OnsetBreeding"
Coefs_phen$Trait[Coefs_phen$Trait == 'LayingDateAllBroods'] <- 'OnsetBreeding'


# 3. Heterogeneity in CZ --------------------------------------------
tre_sub <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% Coefs_phen$Species))
Mat_phylo <- ape::vcv.phylo(tre_sub, corr = TRUE)
diag(Mat_phylo)
# check that all species that are on the tree are also in the dataset
for(i in 1:(nrow(Coefs_phen))){
  check <- match(Coefs_phen$Species[i], tre_sub$tip.label)
  if (is.na(check)){
    print(paste(Coefs_phen$Species[i], 'is not in the phylogeny', sep = ' '))}
}

# other way around
for(i in 1:(nrow(Mat_phylo))){
  check <- match(rownames(Mat_phylo)[i], Coefs_phen$Species)
  if (is.na(check)){
    print(paste(rownames(Mat_phylo)[i], 'is not in the dataset', sep = ' '))}
}

Coefs_phen$Trait <- as.factor(Coefs_phen$Trait)
meta_Phen_CZ_bySpeTrait <- fit_meta_phylo(data_MA = Coefs_phen,
                                  Type_EfS = 'Trait_mean<-det_Clim', COV = 'Pvalue',
                                  Cov_fact = 'Trait',  DD = 'n_effectGR', simpleSEM = TRUE,
                                  A = Mat_phylo,
                                  Trait = FALSE, des.matrix = 'identity')
meta_Phen_CZ_bySpeTrait$data
# Okay, significant effect of the type of measure we use. This is in a way good
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



meta_Phen_CZ_byTaxon <- fit_meta_phylo(data_MA = Coefs_phen,
                                       Type_EfS = 'Trait_mean<-det_Clim', COV = 'Pvalue',
                                       Cov_fact = 'Taxon',  DD = 'n_effectGR', simpleSEM = TRUE,
                                       A = Mat_phylo,
                                       Trait = FALSE, des.matrix = 'identity')
meta_Phen_CZ_byTaxon$data
## also difference by taxon

CZ_byTaxon<- meta_Phen_CZ_byTaxon$data[[1]] %>%
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



# 4. Heterogeneity in ZG and CG --------------------------------------------
meta_Phen_ZG_byTaxon <- fit_meta_phylo(data_MA = Coefs_phen,
                                       Type_EfS = 'GR<-Trait_mean', COV = 'Pvalue',
                                       Cov_fact = 'Taxon',  DD = 'n_effectGR', simpleSEM = TRUE,
                                       A = Mat_phylo,
                                       Trait = FALSE, des.matrix = 'identity')
meta_Phen_ZG_byTaxon$data


ZG_byTaxon<- meta_Phen_ZG_byTaxon$data[[1]] %>%
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


meta_Phen_CG_byTaxon <- fit_meta_phylo(data_MA = Coefs_phen,
                                        Type_EfS = 'GR<-det_Clim', COV = 'Pvalue',
                                        Cov_fact = 'Taxon',  DD = 'n_effectGR', simpleSEM = TRUE,
                                        A = Mat_phylo,
                                        Trait = FALSE, des.matrix = 'identity')
meta_Phen_CG_byTaxon$data



CG_byTaxon<- meta_Phen_CG_byTaxon$data[[1]] %>%
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


