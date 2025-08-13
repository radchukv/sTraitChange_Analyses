## Script with meta-analysis of the SEM results on temperature
## accounting for phylogeny
## using a subset of SEM studies with winDur limited to max 52 weeks
## including covariates only in the needed models
## (i.e. where response is related to climate, and
## hence to the climatic window)
## per Taxon separately

library(metafor)
library(ggplot2)
library(tidyverse)
library(ape)
library(corrplot)
library(sTraitChange)
library(magrittr)

# 0. Data read  ---------------------------------------------------

Coefs_Aut <- readRDS(file = './output_fSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')

## subset to only retain the studies with windur max 52
Coefs_Aut <- Coefs_Aut[Coefs_Aut$WinDur <= 51, ]
length(unique(Coefs_Aut$ID))  ## 199 (from initial 210)


## read in the trait data and merge
traits <- read.csv('./data/speciesTraits.csv')
traits_sub <- subset(traits, select = c(Species, GenLength_y_IUCN))
Coefs_Aut_sp <- merge(Coefs_Aut, traits_sub, by = 'Species', all.x = TRUE)

Coefs_Aut_sp %<>%
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

table(Coefs_Aut_sp$TraitType)

# 1. Prepare data for birds -----------------------------------------------

Coefs_Aut_Bird <- Coefs_Aut_sp %>%
  filter(Taxon == "Bird")



# to be able to compare the models prepare the dataset as for the analyses with phylogeny directly
## first update the species names to correspond to the ones on phylogeny
Coefs_Aut_Bird <- Coefs_Aut_Bird %>%
  mutate(Sp_phylo = case_when(
    Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
    Species == 'Thalasseus sandvicensis' ~ 'Sterna sandvicensis',
    Species == 'Setophaga caerulescens' ~ 'Dendroica caerulescens',
    Species == 'Thalassarche melanophris' ~ 'Thalassarche melanophrys',
    Species == 'Ichthyaetus audouinii' ~ 'Larus audouinii',
    Species == 'Stercorarius maccormicki' ~ 'Catharacta maccormicki',
    TRUE ~ Species))


Coefs_Aut_Bird$Species <- unlist(lapply(1:nrow(Coefs_Aut_Bird), FUN = function(x){
  binary <- strsplit(as.character(Coefs_Aut_Bird$Sp_phylo[x]), " ")
  Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")
}))

# read in the phylo tree
vert_tree <- read.tree('./data/phylogenies_100/vert1.tre')
plot(vert_tree)
is.ultrametric(vert_tree)

pruned_bird <- keep.tip(vert_tree, Coefs_Aut_Bird$Species)
pruned_bird
# get the matrix of phylogenetic correlations
Mat_phylo <- ape::vcv.phylo(pruned_bird, corr = TRUE)


# need a copy with the sp names, to be used to account for phylogeny (apart from that
# we also account for among-sp variation by including species as random intercept
# in the model)
Coefs_Aut_Bird$Sp_phylo <- Coefs_Aut_Bird$Species

# check that all species that are on the tree are also in the dataset
for(i in 1:(nrow(Coefs_Aut_Bird))){
  check <- match(Coefs_Aut_Bird$Species[i], pruned_bird$tip.label)
  if (is.na(check)){
    print(paste(Coefs_Aut_Bird$Species[i], 'is not in the phylogeny', sep = ' '))}
}

# other way around
for(i in 1:(nrow(Mat_phylo))){
  check <- match(rownames(Mat_phylo)[i], Coefs_Aut_Bird$Species)
  if (is.na(check)){
    print(paste(rownames(Mat_phylo)[i], 'is not in the dataset', sep = ' '))}
}

# correcting for the Oppel's study the entry for weather qauality: it is "2"
Coefs_Aut_Bird$WeathQ[Coefs_Aut_Bird$WeathQ == ''] <-  "2"
Coefs_Aut_Bird %<>% mutate(across(where(is_character), as_factor))

# I. BIRDS: Phenology     -------------------


## test function
Coefs_phenClim <- subset(Coefs_Aut_Bird, Relation == 'Trait_mean<-det_Clim' &
                           Trait_Categ == 'Phenological')
tre_sub <- drop.tip(pruned_bird, which(!pruned_bird$tip.label %in% Coefs_phenClim$Species))
Mat_phen_climtr <- ape::vcv.phylo(tre_sub, corr = TRUE)
corrplot(Mat_phen_climtr)

check <- fit_meta_phylo(data_MA = Coefs_phenClim,
                  Type_EfS = 'Trait_mean<-det_Clim',
                  Cov_fact = NULL, COV = NULL,
                  DD = 'n_effectGR', simpleSEM = TRUE, A = Mat_phen_climtr)



## fitting the meta-analyses
# Always have to prepare the subtree for a specific combi of trait and clim!! (here done in the test)
meta_Phen_Cov <- fit_all_meta(data_MA = Coefs_Aut_Bird,
                               Demog_rate = NULL,
                               Trait_categ = 'Phenological',
                               Clim = 'Temperature',
                               Cov_fact = 'WeathQ',
                               COV = 'Pvalue',
                               sel = 'Temp_Phen_Cov_Birds',
                               folder_name = './output_all/',
                               colr = c('black', 'red'),
                               DD = 'n_effectGR',
                               simpleSEM = TRUE,
                               A = Mat_phen_climtr,
                               all_Relations = c('Trait_mean<-det_Clim',
                                                 'GR<-det_Clim',
                                                 'Ind_GR<-det_Clim',
                                                 'Tot_GR<-det_Clim'))


meta_Phen <- fit_all_meta(data_MA = Coefs_Aut_Bird,
                          Demog_rate = NULL,
                          Trait_categ = 'Phenological',
                          Clim = 'Temperature',
                          Cov_fact = NULL,
                          COV = NULL,
                          sel = 'Temp_Phen_Birds',
                          A = Mat_phen_climtr,
                          folder_name = './output_all/',
                          colr = c('black'),
                          DD = 'n_effectGR',
                          simpleSEM = TRUE,
                          all_Relations = c('GR<-Pop_mean',
                                            'GR<-Trait_mean'))


# saving the results also with the random variance estimates, for the SI
meta_Phen$meta_res
colnames(meta_Phen_Cov$meta_res[[1]])[1] <- 'Variable'
phen_Birds_res <-
  bind_rows(meta_Phen$meta_res, meta_Phen_Cov$meta_res) %>%
  filter(Variable == 'intrcpt') %>%
  select(Relation, Species.SD, ID.SD, Location.SD)


phen_Birds_res$Species.SD <- format(phen_Birds_res$Species.SD, scientific = TRUE)
phen_Birds_res$ID.SD <- format(phen_Birds_res$ID.SD,  scientific = TRUE)
phen_Birds_res$Location.SD <- format(phen_Birds_res$Location.SD,  scientific = TRUE)

## save this for the SI
save_xlsx(table = phen_Birds_res, table_name = './tables/random_sd_PhenTempBirds')


hist(Coefs_Aut_Bird$GenLength_y_IUCN)


## checking the effect of Generation length
# have to update the tree because generation lentgh is not available for all species, so some
# rows are dropped inside the function
GenLength_phen <- Coefs_phenClim[! is.na(Coefs_phenClim$GenLength_y_IUCN), ]
tre_GL_phenClim <- drop.tip(pruned_bird, which(!pruned_bird$tip.label %in% GenLength_phen$Species))
# plot(tre_GL_phenClim)
# is.ultrametric(tre_GL_phenClim)

Mat_GL_climtr <- ape::vcv.phylo(tre_GL_phenClim, corr = TRUE)
corrplot(Mat_GL_climtr)

meta_Phen_Glen <- fit_all_meta(data_MA = Coefs_Aut_Bird,
                               Clim = 'Temperature',
                               Demog_rate = NULL,
                               Trait_categ = 'Phenological',
                               Cov_fact = NULL,
                               COV = 'GenLength_y_IUCN',
                               sel = 'Temp_Phen_GenL_Birds',
                               A = Mat_GL_climtr,
                               folder_name = './output_all/',
                               colr = c('black'),
                               DD = 'n_effectGR',
                               simpleSEM = TRUE,
                               all_Relations = c('Ind_GR<-det_Clim', 'Tot_GR<-det_Clim',
                                                 'GR<-Trait_mean'))
meta_Phen_Glen$meta_res


# 2. Prepare data for mammals -----------------------------------------------

Coefs_Aut_Mam <- Coefs_Aut_sp %>%
  filter(Taxon == "Mammal")

Coefs_Aut_Mam$Sp_phylo <- Coefs_Aut_Mam$Species

Coefs_Aut_Mam$Species <- unlist(lapply(1:nrow(Coefs_Aut_Mam), FUN = function(x){
  binary <- strsplit(as.character(Coefs_Aut_Mam$Sp_phylo[x]), " ")
  Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")
}))
length(unique(Coefs_Aut_Mam$Species))  #9

# prun the vertebrate tree to mammals only

pruned_mam <- keep.tip(vert_tree, Coefs_Aut_Mam$Species)
pruned_mam
plot(pruned_mam)
# get the matrix of phylogenetic correlations
Mat_phylo <- ape::vcv.phylo(pruned_mam, corr = TRUE)


# need a copy with the sp names, to be used to account for phylogeny (apart from that
# we also account for among-sp variation by including species as random intercept
# in the model)
Coefs_Aut_Mam$Sp_phylo <- Coefs_Aut_Mam$Species

# check that all species that are on the tree are also in the dataset
for(i in 1:(nrow(Coefs_Aut_Mam))){
  check <- match(Coefs_Aut_Mam$Species[i], pruned_mam$tip.label)
  if (is.na(check)){
    print(paste(Coefs_Aut_Mam$Species[i], 'is not in the phylogeny', sep = ' '))}
}

# other way around
for(i in 1:(nrow(Mat_phylo))){
  check <- match(rownames(Mat_phylo)[i], Coefs_Aut_Mam$Species)
  if (is.na(check)){
    print(paste(rownames(Mat_phylo)[i], 'is not in the dataset', sep = ' '))}
}

# correcting for the Oppel's study the entry for weather qauality: it is "2"
Coefs_Aut_Mam$WeathQ[Coefs_Aut_Mam$WeathQ == ''] <-  "2"
Coefs_Aut_Mam %<>% mutate(across(where(is_character), as_factor))

# II. MAMMALS: Phenology     -------------------

## test function
Coefs_phenClim_mam <- subset(Coefs_Aut_Mam, Relation == 'Trait_mean<-det_Clim' &
                           Trait_Categ == 'Phenological')
tre_sub <- drop.tip(pruned_mam, which(!pruned_mam$tip.label %in% Coefs_phenClim_mam$Species))
Mat_phen_climtr_mam <- ape::vcv.phylo(tre_sub, corr = TRUE)
corrplot(Mat_phen_climtr_mam)

check <- fit_meta_phylo(data_MA = Coefs_phenClim_mam,
                        Type_EfS = 'Trait_mean<-det_Clim',
                        Cov_fact = NULL, COV = NULL,
                        DD = 'n_effectGR', simpleSEM = TRUE, A = Mat_phen_climtr_mam)



## fitting the meta-analyses
# Always have to prepare the subtree for a specific combi of trait and clim!! (here done in the test)
meta_Phen_Cov_Mam <- fit_all_meta(data_MA = Coefs_Aut_Mam,
                              Demog_rate = NULL,
                              Trait_categ = 'Phenological',
                              Clim = 'Temperature',
                              Cov_fact = 'WeathQ',
                              COV = 'Pvalue',
                              sel = 'Temp_Phen_Cov_Mammals',
                              folder_name = './output_all/',
                              colr = c('black', 'red'),
                              DD = 'n_effectGR',
                              simpleSEM = TRUE,
                              A = Mat_phen_climtr_mam,
                              all_Relations = c('Trait_mean<-det_Clim',
                                                'GR<-det_Clim',
                                                'Ind_GR<-det_Clim',
                                                'Tot_GR<-det_Clim'))


meta_Phen_mam <- fit_all_meta(data_MA = Coefs_Aut_Mam,
                          Demog_rate = NULL,
                          Trait_categ = 'Phenological',
                          Clim = 'Temperature',
                          Cov_fact = NULL,
                          COV = NULL,
                          sel = 'Temp_Phen_Mammals',
                          A = Mat_phen_climtr_mam,
                          folder_name = './output_all/',
                          colr = c('black'),
                          DD = 'n_effectGR',
                          simpleSEM = TRUE,
                          all_Relations = c('GR<-Pop_mean',
                                            'GR<-Trait_mean'))


# saving the results also with the random variance estimates, for the SI
meta_Phen_mam$meta_res
colnames(meta_Phen_Cov_Mam$meta_res[[1]])[1] <- 'Variable'
phen_Mams_res <-
  bind_rows(meta_Phen_mam$meta_res, meta_Phen_Cov_Mam$meta_res) %>%
  filter(Variable == 'intrcpt') %>%
  select(Relation, Species.SD, ID.SD, Location.SD)


phen_Mams_res$Species.SD <- format(phen_Mams_res$Species.SD, scientific = TRUE)
phen_Mams_res$ID.SD <- format(phen_Mams_res$ID.SD,  scientific = TRUE)
phen_Mams_res$Location.SD <- format(phen_Mams_res$Location.SD,  scientific = TRUE)

## save this for the SI
save_xlsx(table = phen_Mams_res, table_name = './tables/random_sd_PhenTempMammals')


# III. Summary across groups -------------------------------------------

## binding all across-studies effect sizes together
colnames(meta_Phen_Cov$meta_res[[1]])[1] <- 'Variable'
ef_Birds <- rbind(meta_Phen$meta_res[[1]], meta_Phen_Cov$meta_res[[1]]) %>%
  dplyr::mutate(., Taxon = 'Bird')

colnames(meta_Phen_Cov_Mam$meta_res[[1]])[1] <- 'Variable'
ef_Mam <- rbind(meta_Phen_mam$meta_res[[1]], meta_Phen_Cov_Mam$meta_res[[1]]) %>%
  dplyr::mutate(., Taxon = 'Mammal')

ef_all <- rbind(ef_Birds, ef_Mam)
ef_all %<>%
  dplyr::mutate(., col = dplyr::case_when(EfS_Low < 0 & EfS_Upper < 0 ~ 'orange',
                                          EfS_Low > 0 & EfS_Upper > 0 ~ 'lightblue',
                                          TRUE ~ 'darkgrey'),
                col_var = dplyr::case_when(Variable == 'intrcpt' & pval_Covar < 0.05  ~ 'Intercept significant',
                                           Variable == 'intrcpt' & pval_Covar >= 0.05 ~ 'Intercept non-significant',
                                           Variable == 'Pvalue' & pval_Covar < 0.05 ~ 'Pvalue significant',
                                           Variable == 'Pvalue' & pval_Covar >= 0.05 ~ 'Pvalue non-significant',
                                           Variable == 'WeathQ2' & pval_Covar < 0.05 ~ 'Weather quality significant',
                                           Variable == 'WeathQ2' & pval_Covar >= 0.05 ~ 'Weather quality non-significant'),
                REL = dplyr::case_when(Relation == 'Trait_mean<-det_Clim' ~ 'CZ',
                                       Relation == 'GR<-Trait_mean' ~ 'ZG',
                                       Relation == 'GR<-det_Clim' ~ 'CG',
                                       Relation == 'GR<-Pop_mean' ~ 'PG',
                                       Relation == 'Ind_GR<-det_Clim' ~ 'CZG',
                                       Relation == 'Tot_GR<-det_Clim' ~ 'TotalCG'))


## sort the factor so as to have the coefficients from bottom to the top of the scheme
ef_all$REL <- factor(ef_all$REL, levels = c('CZ', 'ZG','CG', 'PG', 'CZG', 'TotalCG'))

## sort the data frme to have the levels in the desired order
ord <- order(ef_all$Taxon, ef_all$REL)
ef_all <- ef_all[ord, ]

## save ef_all
saveRDS(object = ef_all, file = './output_all/all_efSizes_temperature_Phen_ByTaxon.RDS')
# ef_all <- readRDS(file = './output_all/all_efSizes_temperature_Phen_ByTaxon.RDS')

# IV. Extract heterogeneity metrics -------------------------------------------

heter_data_phen1 <- data.table::rbindlist(lapply(meta_Phen_Cov$data_meta[[1]]$heter_mod, function(x){x[1, -4]}))
heter_data_phen1$name <- meta_Phen_Cov$data_meta[[1]]$names

heter_data_phen2 <- data.table::rbindlist(lapply(meta_Phen$data_meta[[1]]$heter_mod, function(x){x[1, -4]}))
heter_data_phen2$name <- meta_Phen$data_meta[[1]]$names

heter_data_phen <- rbind(heter_data_phen1, heter_data_phen2)
heter_data_phen$Taxon <- 'Bird'

heter_data_Mam <- data.table::rbindlist(lapply(meta_Phen_Cov_Mam$data_meta[[1]]$heter_mod, function(x){x[1, -4]}))
heter_data_Mam$name <- meta_Phen_Cov_Mam$data_meta[[1]]$names

heter_data_Mam2 <- data.table::rbindlist(lapply(meta_Phen_mam$data_meta[[1]]$heter_mod, function(x){x[1, -4]}))
heter_data_Mam2$name <- meta_Phen_mam$data_meta[[1]]$names

heter_data_phen_Mam <- rbind(heter_data_Mam, heter_data_Mam2)
heter_data_phen_Mam$Taxon <- 'Mammal'

heter_data_bothT <- rbind(heter_data_phen, heter_data_phen_Mam)
heter_data_bothT$Climate <- 'Temperature'
## recode the Relations name
heter_data_bothT %<>%
  dplyr::mutate(Relation =
                  dplyr::recode_factor(name, 'Trait_mean<-det_Clim' = 'CZ',
                                       'GR<-Trait_mean' = 'ZG',
                                       'GR<-det_Clim' = 'CG',
                                       'GR<-Pop_mean' = 'PG',
                                       'Ind_GR<-det_Clim' = 'CZG',
                                       'Tot_GR<-det_Clim' = 'TotalCG'))
heter_data_bothT <- heter_data_bothT[order(heter_data_bothT$Taxon,
                                           heter_data_bothT$Relation), ]



heter_data_bothT$Q <- format(round(heter_data_bothT$Q, 1), scientific = FALSE)
heter_data_bothT$I2 <- format(round(heter_data_bothT$I2, 3),  scientific = FALSE)
heter_data_bothT$pval <- numeric(length = nrow(heter_data_bothT))
for(i in 1:nrow(heter_data_bothT)){
  if(heter_data_bothT$Qp[i] < 0.001){
    heter_data_bothT$pval[i] <- '<0.001'
  } else {
    heter_data_bothT$pval[i] <- format(round(heter_data_bothT$Qp[i], 3), nsmall = 3, scientific = FALSE)
  }
}
heter_data_bothT$Qp <- NULL
heter_data_bothT$name <- NULL

## save this file for the SI
save_xlsx(table = heter_data_bothT, table_name = './tables/Heterogen_MetaAnalyses_TemperaturePhenology_perTaxon')


#  V. Extract all effect sizes, also indirect ones, estimated with meta-analyses --------
## I have to extract the Coefs that were estimated within the meta-analysis, ie.
## 'Ind_DemRate<-det_Clim', 'Ind_GR<-det_Clim', 'Tot_DemRate<-det_Clim', 'Tot_GR<-det_Clim'
## and then merge that data with the Coefs_Aut and save on the disk

# test the function
two <- dplyr::bind_rows(lapply(list('Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'),
                               function(x){extr_coefs(obj = meta_Phen_Cov, Type_EfS = x)}))

## cycling through all the subsets and all estimated path coefficients to get them all into one dataset
efSizes_estim <- dplyr::bind_rows(lapply(list(meta_Phen_Cov, meta_Phen_Cov_Mam), function(x){
  lapply(list('Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'),
         function(y){extr_coefs(obj = x, Type_EfS = y)})
}))


Coefs_Aut_subst <- Coefs_Aut_sp %>%
  dplyr::rename(.,  SError = Std.Error)

## col names that we want to exclude
names(Coefs_Aut_subst)[! names(Coefs_Aut_subst) %in% names(efSizes_estim)]

## add a column P.Value to efSizes_estim, because this column is of interest for meta-analyses
## (even if not available for the indirect path)
efSizes_estim$P.Value = rep(NA, nrow(efSizes_estim))

Coefs_Aut_subst <- Coefs_Aut_subst %>%
  dplyr::select(., -c(names(Coefs_Aut_subst)[! names(Coefs_Aut_subst) %in% names(efSizes_estim)]))


names(efSizes_estim)[! names(efSizes_estim) %in% names(Coefs_Aut_subst)]
efSizes_estim <- efSizes_estim %>%
  dplyr::select(., -c(names(efSizes_estim)[! names(efSizes_estim) %in% names(Coefs_Aut_subst)]))
Coef_all_T <- rbind(Coefs_Aut_subst, efSizes_estim)

# save all the effect sizes
saveRDS(object = Coef_all_T, file = './output_all/PathCoefs_Temp_Phen_PerTaxon_AlsoEstimatedRelations.RDS')
#check <- readRDS('./output_all/PathCoefs_Temp_AlsoEstimatedRelations.RDS')

