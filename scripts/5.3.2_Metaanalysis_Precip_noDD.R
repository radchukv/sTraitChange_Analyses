## Script with meta-analysis of the SEM results on precipitation
## SEM models fitted without the effect of population size on pop
# growth rate
## accounting for phylogeny
## using a subset of SEM studies with winDur limited to max 52 weeks
## including covariates only in the needed models (i.e.
## where response is related to climate, and hence to the climatic window)

library(tidyverse)
library(metafor)
library(ggplot2)
library(magrittr)
library(ape)
library(corrplot)
library(sTraitChange)


## Read in the data
Coefs_Aut <- readRDS(file = './output_fSEM_noDD_precip/PathCoefs_allMods_Precip_Weights_noDD_Autocor.RDS')
table(Coefs_Aut$Taxon)

## subset to only retain the studies with windur max 52
Coefs_Aut <- Coefs_Aut[Coefs_Aut$WinDur <= 51, ]
length(unique(Coefs_Aut$ID))  ## 209


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


# read in the phylo tree
vert_tree <- read.tree('./data/phylogenies_100/vert1.tre')
plot(vert_tree)
is.ultrametric(vert_tree)

# drop the atlantic yellow-nosed as it is not in the precip dataset - no precipitatio
# for that island
vert_tree <- drop.tip(vert_tree, "Thalassarche_chlororhynchos")
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


Coefs_Aut_sp %<>% mutate(across(where(is_character), as_factor))


# I. Phenology     -------------------

## test function
Coefs_phenClim <- subset(Coefs_Aut_sp, Relation == 'Trait_mean<-det_Clim' &
                           Trait_Categ == 'Phenological')
tre_sub <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% Coefs_phenClim$Species))
# plot(tre_sub)

Mat_phen_climtr <- ape::vcv.phylo(tre_sub, corr = TRUE)
corrplot(Mat_phen_climtr)

check <- fit_meta_phylo(data_MA = Coefs_phenClim,
                        Type_EfS = 'Trait_mean<-det_Clim',
                        Cov_fact = "WeathQ", COV = "Pvalue",
                        DD = 'none', simpleSEM = TRUE,
                        A = Mat_phen_climtr)


## some quick exploration for the pValue
Coefs_Aut_sp %>%
  dplyr::filter(.,  Relation == 'Trait_mean<-det_Clim') %>%
  ggplot(., aes(x = Pvalue, y = Estimate)) + geom_point() +
  facet_grid(rows = vars(Trait_Categ), cols= vars(Demog_rate_Categ)) +
  geom_smooth(method = 'lm', col ='red') + theme_bw()


## split by trait categ only
Coefs_Aut_sp %>%
  dplyr::filter(.,  Relation == 'Trait_mean<-det_Clim') %>%
  ggplot(., aes(x = Pvalue, y = Estimate)) + geom_point() +
  facet_grid(rows = vars(Trait_Categ)) +
  geom_smooth(method = 'lm', col ='red') + theme_bw()  ## seems ot be relevant for phenol traits only

## fitting the meta-analyses
meta_Phen_Cov <- fit_all_meta(data_MA = Coefs_Aut_sp,
                              Demog_rate = NULL,
                              Trait_categ = 'Phenological',
                              Clim = 'Precipitation',
                              Cov_fact = 'WeathQ',
                              COV = 'Pvalue',
                              sel = 'Precip_Phen_Cov_noDD',
                              A = Mat_phen_climtr,
                              folder_name = './output_all_noDD/',
                              colr = c('black', 'red'),
                              DD = 'none',
                              simpleSEM = TRUE,
                              all_Relations = c('Trait_mean<-det_Clim',
                                                'GR<-det_Clim',
                                                'Ind_GR<-det_Clim',
                                                'Tot_GR<-det_Clim'))


meta_Phen <- fit_all_meta(data_MA = Coefs_Aut_sp,
                          Demog_rate = NULL,
                          Trait_categ = 'Phenological',
                          Clim = 'Precipitation',
                          Cov_fact = NULL,
                          COV = NULL,
                          sel = 'Precip_Phen_noDD',
                          A = Mat_phen_climtr,
                          folder_name = './output_all_noDD/',
                          colr = c('black'),
                          DD = 'none',
                          simpleSEM = TRUE,
                          all_Relations = c('GR<-Trait_mean'))



# II. Morphology -------------------------------------------

## fitting the meta-analyses
# prepare the tree
Coefs_morphClim <- subset(Coefs_Aut_sp, Relation == 'Trait_mean<-det_Clim' &
                            Trait_Categ == 'Morphological')
tre_sub <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% Coefs_morphClim$Species))
# plot(tre_sub)
# is.ultrametric(tre_sub)

Mat_morph_climtr <- ape::vcv.phylo(tre_sub, corr = TRUE)
corrplot(Mat_morph_climtr)

meta_Morph_Cov <- fit_all_meta(data_MA = Coefs_Aut_sp,
                               Demog_rate = NULL,
                               Trait_categ = 'Morphological',
                               Clim = 'Precipitation',
                               Cov_fact = 'WeathQ',
                               COV = 'Pvalue',
                               sel = 'Precip_Morph_Cov_noDD',
                               A = Mat_morph_climtr,
                               folder_name = './output_all_noDD/',
                               colr = c('black', 'red'),
                               DD = 'none',
                               simpleSEM = TRUE,
                               all_Relations = c('Trait_mean<-det_Clim',
                                                 'GR<-det_Clim',
                                                 'Ind_GR<-det_Clim',
                                                 'Tot_GR<-det_Clim'))



meta_Morph <- fit_all_meta(data_MA = Coefs_Aut_sp,
                           Demog_rate = NULL,
                           Trait_categ = 'Morphological',
                           Clim = 'Precipitation',
                           Cov_fact = NULL,
                           COV = NULL,
                           sel = 'Precip_Morph_noDD',
                           A = Mat_morph_climtr,
                           folder_name = './output_all_noDD/',
                           colr = c('black'),
                           DD = 'none',
                           simpleSEM = TRUE,
                           all_Relations = c('GR<-Trait_mean'))



# III. Summary across groups  ---------------------------------------------

## binding all across-studies effect sizes together
## binding all across-studies effect sizes together
colnames(meta_Phen_Cov$meta_res[[1]])[1] <- 'Variable'
ef_Phen <- rbind(meta_Phen$meta_res[[1]], meta_Phen_Cov$meta_res[[1]]) %>%
  dplyr::mutate(., Trait_Categ = 'Phenology')

colnames(meta_Morph_Cov$meta_res[[1]])[1] <- 'Variable'
ef_Morph <- rbind(meta_Morph$meta_res[[1]], meta_Morph_Cov$meta_res[[1]]) %>%
  dplyr::mutate(., Trait_Categ = 'Morphology')

ef_all <- rbind(ef_Phen, ef_Morph)
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
ord <- order(ef_all$Trait_Categ, ef_all$REL)
ef_all <- ef_all[ord, ]

## save ef_all
saveRDS(object = ef_all, file = './output_all_noDD/all_efSizes_noDD_precip.RDS')



# IV. Extract all effect sizes, also indirect ones, estimated with meta-analyses --------

meta_Phen_Cov$data_meta[[1]]$data_EfS[meta_Phen_Cov$data_meta[[1]]$names == 'Ind_GR<-det_Clim']

## so, i have to extract the Coefs that were estimated within the meta-analysis, ie.
## 'Ind_DemRate<-det_Clim', 'Ind_GR<-det_Clim', 'Tot_DemRate<-det_Clim', 'Tot_GR<-det_Clim'
## and then merge that data with the Coefs_Aut and save on the disk

## cycling through all the subsets and all estimated path coefficients to get them all into one dataset
efSizes_estim <- dplyr::bind_rows(lapply(list(meta_Phen_Cov, meta_Morph_Cov), function(x){
  lapply(list('Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'),
         function(y){extr_coefs(obj = x, Type_EfS = y)})
}))

Coefs_Aut_subst <- Coefs_Aut_sp %>%
  dplyr::rename(.,  SError = Std.Error)

## col names that we want to exclude
names(Coefs_Aut_subst)[! names(Coefs_Aut_subst) %in% names(efSizes_estim)]

## add a column P.Value to efSizes_estim, because this column is of interest for meta-analyses (even if not available for the indirect path)
efSizes_estim$P.Value = rep(NA, nrow(efSizes_estim))

Coefs_Aut_subst <- Coefs_Aut_subst %>%
  dplyr::select(., -c(names(Coefs_Aut_subst)[! names(Coefs_Aut_subst) %in% names(efSizes_estim)]))


names(efSizes_estim)[! names(efSizes_estim) %in% names(Coefs_Aut_subst)]
efSizes_estim <- efSizes_estim %>%
  dplyr::select(., -c(names(efSizes_estim)[! names(efSizes_estim) %in% names(Coefs_Aut_subst)]))
Coef_all_T <- rbind(Coefs_Aut_subst, efSizes_estim)

# save all the effect sizes
saveRDS(object = Coef_all_T, file = './output_all_noDD/PathCoefs_noDD_Precip_AlsoEstimatedRelations.RDS')
