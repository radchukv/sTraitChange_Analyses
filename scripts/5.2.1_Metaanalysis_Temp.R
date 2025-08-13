## Script with meta-analysis of the SEM results on temperature
## accounting for phylogeny
## using a subset of SEM studies with winDur limited to max 52 weeks
## including covariates only in the needed models
## (i.e. where response is related to climate, and
## hence to the climatic window)

library(metafor)
library(ggplot2)
library(tidyverse)
library(ape)
library(corrplot)
library(sTraitChange)
library(magrittr)

# 0. Data read and prep ---------------------------------------------------

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

# I. Phenology     -------------------


## test function
Coefs_phenClim <- subset(Coefs_Aut_sp, Relation == 'Trait_mean<-det_Clim' &
                           Trait_Categ == 'Phenological')
tre_sub <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% Coefs_phenClim$Species))
Mat_phen_climtr <- ape::vcv.phylo(tre_sub, corr = TRUE)
corrplot(Mat_phen_climtr)

check <- fit_meta_phylo(data_MA = Coefs_phenClim,
                  Type_EfS = 'Trait_mean<-det_Clim',
                  Cov_fact = NULL, COV = NULL,
                  DD = 'n_effectGR', simpleSEM = TRUE, A = Mat_phen_climtr)

check_cov <- fit_meta_phylo(data_MA = Coefs_phenClim, Type_EfS = 'Trait_mean<-det_Clim',
                  Cov_fact = 'Taxon', COV = NULL,
                  DD = 'n_effectGR', simpleSEM = TRUE,
                  A = Mat_phen_climtr, des.matrix = 'identity')


## some quick exploration for the pValue

## split by trait categ only
dat_lab <- data.frame(x = c(0,0), y = c(1.2,1.2),
                      label = c('a', 'b'),
                      Trait_Categ = c('Morphological', 'Phenological'))
pdf('./output_all/Relation_PValue&TempEstimate_ByTraitCateg_ClimTraitRelation.pdf',
    width  = 10)
Coefs_Aut_sp %>%
  dplyr::filter(.,  Relation == 'Trait_mean<-det_Clim') %>%
  ggplot(., aes(x = Pvalue, y = Estimate)) + geom_point() +
  facet_wrap(. ~ Trait_Categ) +
  geom_smooth(method = 'lm', col ='red') + theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  xlab(bquote(P[Delta*AICc])) +
  geom_text(data = dat_lab, aes(x = x, y = y, label = label),
            fontface = 'bold', size = 6)
dev.off()


## fitting the meta-analyses
# Always have to prepare the subtree for a specific combi of trait and clim!! (here done in the test)
meta_Phen_Cov <- fit_all_meta(data_MA = Coefs_Aut_sp,
                               Demog_rate = NULL,
                               Trait_categ = 'Phenological',
                               Clim = 'Temperature',
                               Cov_fact = 'WeathQ',
                               COV = 'Pvalue',
                               sel = 'Temp_Phen_Cov',
                               folder_name = './output_all/',
                               colr = c('black', 'red'),
                               DD = 'n_effectGR',
                               simpleSEM = TRUE,
                               A = Mat_phen_climtr,
                               all_Relations = c('Trait_mean<-det_Clim',
                                                 'GR<-det_Clim',
                                                 'Ind_GR<-det_Clim',
                                                 'Tot_GR<-det_Clim'))

## for sensitivity analysis - no p value
meta_Phen_noP <- fit_all_meta(data_MA = Coefs_Aut_sp,
                              Demog_rate = NULL,
                              Trait_categ = 'Phenological',
                              Clim = 'Temperature',
                              Cov_fact = 'WeathQ',
                              COV = NULL,
                              sel = 'Temp_Phen_noP',
                              A = Mat_phen_climtr,
                              folder_name = './output_all/',
                              colr = c('black', 'red'),
                              DD = 'n_effectGR',
                              simpleSEM = TRUE,
                              all_Relations = c('Trait_mean<-det_Clim',
                                                'GR<-det_Clim',
                                                'Ind_GR<-det_Clim',
                                                'Tot_GR<-det_Clim'))


meta_Phen <- fit_all_meta(data_MA = Coefs_Aut_sp,
                          Demog_rate = NULL,
                          Trait_categ = 'Phenological',
                          Clim = 'Temperature',
                          Cov_fact = NULL,
                          COV = NULL,
                          sel = 'Temp_Phen',
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
phen_res <-
  bind_rows(meta_Phen$meta_res, meta_Phen_Cov$meta_res) %>%
  filter(Variable == 'intrcpt') %>%
  select(Relation, Species.SD, ID.SD, Location.SD)


phen_res$Species.SD <- format(phen_res$Species.SD, scientific = TRUE)
phen_res$ID.SD <- format(phen_res$ID.SD,  scientific = TRUE)
phen_res$Location.SD <- format(phen_res$Location.SD,  scientific = TRUE)

## save this for the SI
save_xlsx(table = phen_res, table_name = './tables/random_sd_PhenTemp')


hist(Coefs_Aut_sp$GenLength_y_IUCN)


## checking the effect of Generation length
# have to update the tree because generation lentgh is not available for all species, so some
# rows are dropped inside the function
GenLength_phen <- Coefs_phenClim[! is.na(Coefs_phenClim$GenLength_y_IUCN), ]
tre_GL_phenClim <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% GenLength_phen$Species))
# plot(tre_GL_phenClim)
# is.ultrametric(tre_GL_phenClim)

Mat_GL_climtr <- ape::vcv.phylo(tre_GL_phenClim, corr = TRUE)
corrplot(Mat_GL_climtr)

meta_Phen_Glen <- fit_all_meta(data_MA = Coefs_Aut_sp,
                               Clim = 'Temperature',
                               Demog_rate = NULL,
                               Trait_categ = 'Phenological',
                               Cov_fact = NULL,
                               COV = 'GenLength_y_IUCN',
                               sel = 'Temp_Phen_GenL',
                               A = Mat_GL_climtr,
                               folder_name = './output_all/',
                               colr = c('black'),
                               DD = 'n_effectGR',
                               simpleSEM = TRUE,
                               all_Relations = c('Ind_GR<-det_Clim', 'Tot_GR<-det_Clim',
                                                 'GR<-Trait_mean'))



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
                              Clim = 'Temperature',
                              Cov_fact = 'WeathQ',
                              COV = 'Pvalue',
                              sel = 'Temp_Morph_Cov',
                              A = Mat_morph_climtr,
                              folder_name = './output_all/',
                              colr = c('black', 'red'),
                              DD = 'n_effectGR',
                              simpleSEM = TRUE,
                              all_Relations = c('Trait_mean<-det_Clim',
                                                'GR<-det_Clim',
                                                'Ind_GR<-det_Clim',
                                                'Tot_GR<-det_Clim'))

meta_Morph_noP <- fit_all_meta(data_MA = Coefs_Aut_sp,
                               Demog_rate = NULL,
                               Trait_categ = 'Morphological',
                               Clim = 'Temperature',
                               Cov_fact = 'WeathQ',
                               COV = NULL,
                               sel = 'Temp_Morph_noP',
                               A = Mat_morph_climtr,
                               folder_name = './output_all/',
                               colr = c('black', 'red'),
                               DD = 'n_effectGR',
                               simpleSEM = TRUE,
                               all_Relations = c('Trait_mean<-det_Clim',
                                                 'GR<-det_Clim',
                                                 'Ind_GR<-det_Clim',
                                                 'Tot_GR<-det_Clim'))

meta_Morph <- fit_all_meta(data_MA = Coefs_Aut_sp,
                          Demog_rate = NULL,
                          Trait_categ = 'Morphological',
                          Clim = 'Temperature',
                          Cov_fact = NULL,
                          COV = NULL,
                          sel = 'Temp_Morph',
                          A = Mat_morph_climtr,
                          folder_name = './output_all/',
                          colr = c('black'),
                          DD = 'n_effectGR',
                          simpleSEM = TRUE,
                          all_Relations = c('GR<-Pop_mean',
                                            'GR<-Trait_mean'))

hist(Coefs_Aut_sp$GenLength_y_IUCN)

## checking the effect of Generation length
# here I should prepare a new tree, with a subset of species (as gen length i snot available for all, and
# function then removes those spp for which no gen length is available)
GenLength_morph <- Coefs_morphClim[! is.na(Coefs_morphClim$GenLength_y_IUCN), ]
tre_GL_morphClim <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% GenLength_morph$Species))
# plot(tre_GL_phenClim)
# is.ultrametric(tre_GL_phenClim)

Mat_GL_climMorph <- ape::vcv.phylo(tre_GL_morphClim, corr = TRUE)
corrplot(Mat_GL_climMorph)
meta_Morph_Glen <- fit_all_meta(data_MA = Coefs_Aut_sp,
                               Clim = 'Temperature',
                               Demog_rate = NULL,
                               Trait_categ = 'Morphological',
                               Cov_fact = NULL,
                               COV = 'GenLength_y_IUCN',
                               sel = 'Temp_Morph_GenL',
                               A = Mat_GL_climMorph,
                               folder_name = './output_all/',
                               colr = c('black'),
                               DD = 'n_effectGR',
                               simpleSEM = TRUE,
                               all_Relations = c('Ind_GR<-det_Clim', 'Tot_GR<-det_Clim',
                                                 'GR<-Trait_mean'))


# III. Summary across groups -------------------------------------------

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
saveRDS(object = ef_all, file = './output_all/all_efSizes_temperature_cor.RDS')
# ef_all <- readRDS(file = './output_all/all_efSizes_temperature.RDS')

# IV. Extract heterogeneity metrics -------------------------------------------

heter_data_phen1 <- data.table::rbindlist(lapply(meta_Phen_Cov$data_meta[[1]]$heter_mod, function(x){x[1, -4]}))
heter_data_phen1$name <- meta_Phen_Cov$data_meta[[1]]$names

heter_data_phen2 <- data.table::rbindlist(lapply(meta_Phen$data_meta[[1]]$heter_mod, function(x){x[1, -4]}))
heter_data_phen2$name <- meta_Phen$data_meta[[1]]$names

heter_data_phen <- rbind(heter_data_phen1, heter_data_phen2)
heter_data_phen$Trait_Category <- 'Phenological'

heter_data_morph1 <- data.table::rbindlist(lapply(meta_Morph_Cov$data_meta[[1]]$heter_mod, function(x){x[1, -4]}))
heter_data_morph1$name <- meta_Morph_Cov$data_meta[[1]]$names

heter_data_morph2 <- data.table::rbindlist(lapply(meta_Morph$data_meta[[1]]$heter_mod, function(x){x[1, -4]}))
heter_data_morph2$name <- meta_Morph$data_meta[[1]]$names

heter_data_morph <- rbind(heter_data_morph1, heter_data_morph2)
heter_data_morph$Trait_Category <- 'Morphological'

heter_data_temp <- rbind(heter_data_morph, heter_data_phen)
heter_data_temp$Climate <- 'Temperature'
## recode the Relations name
heter_data_temp %<>%
  dplyr::mutate(Relation =
                  dplyr::recode_factor(name, 'Trait_mean<-det_Clim' = 'CZ',
                                       'GR<-Trait_mean' = 'ZG',
                                       'GR<-det_Clim' = 'CG',
                                       'GR<-Pop_mean' = 'PG',
                                       'Ind_GR<-det_Clim' = 'CZG',
                                       'Tot_GR<-det_Clim' = 'TotalCG'))
heter_data_temp <- heter_data_temp[order(heter_data_temp$Trait_Category,
                                         heter_data_temp$Relation), ]



heter_data_temp$Q <- format(round(heter_data_temp$Q, 1), scientific = FALSE)
heter_data_temp$I2 <- format(round(heter_data_temp$I2, 3),  scientific = FALSE)
heter_data_temp$pval <- numeric(length = nrow(heter_data_temp))
for(i in 1:nrow(heter_data_temp)){
  if(heter_data_temp$Qp[i] < 0.001){
    heter_data_temp$pval[i] <- '<0.001'
  } else {
    heter_data_temp$pval[i] <- format(round(heter_data_temp$Qp[i], 3), nsmall = 3, scientific = FALSE)
  }
}
heter_data_temp$Qp <- NULL
heter_data_temp$name <- NULL

## save this file for the SI
save_xlsx(table = heter_data_temp, table_name = './tables/Heterogen_MetaAnalyses_Temperature')

# V. Summary across groups: noP -------------------------------------------

## binding all across-studies effect sizes together
ef_Phen_noP <- meta_Phen_noP$meta_res[[1]] %>%
  dplyr::mutate(., Trait_Categ = 'Phenology')

ef_Morph_noP <- meta_Morph_noP$meta_res[[1]] %>%
  dplyr::mutate(., Trait_Categ = 'Morphology')

ef_all_noP <- rbind(ef_Phen_noP, ef_Morph_noP)
ef_all_noP %<>%
  dplyr::mutate(., col = dplyr::case_when(EfS_Low < 0 & EfS_Upper < 0 ~ 'orange',
                                          EfS_Low > 0 & EfS_Upper > 0 ~ 'lightblue',
                                          TRUE ~ 'darkgrey'),
                col_var = dplyr::case_when(Levels_Covar == 'intrcpt' & pval_Covar < 0.05  ~ 'Intercept significant',
                                           Levels_Covar == 'intrcpt' & pval_Covar >= 0.05 ~ 'Intercept non-significant',
                                           Levels_Covar == 'Pvalue' & pval_Covar < 0.05  ~ 'Pvalue significant',
                                           Levels_Covar == 'Pvalue' & pval_Covar >= 0.05 ~ 'Pvalue non-significant',
                                           Levels_Covar == 'WeathQ2' & pval_Covar < 0.05 ~ 'Weather quality significant',
                                           Levels_Covar == 'WeathQ2' & pval_Covar >= 0.05 ~ 'Weather quality non-significant'),
                REL = dplyr::case_when(Relation == 'Trait_mean<-det_Clim' ~ 'CZ',
                                       Relation == 'Ind_GR<-det_Clim' ~ 'CZG',
                                       Relation == 'Tot_GR<-det_Clim' ~ 'TotalCG',
                                       Relation == 'GR<-det_Clim' ~ 'CG'))


## sort the factor so as to have the coefficients from bottom to the top of the scheme
ef_all_noP$REL <- factor(ef_all_noP$REL, levels = c('CZ', 'CG', 'CZG', 'TotalCG'))

## sort the data frme to have the levels in the desired order
ord_noP <- order(ef_all_noP$Trait_Categ, ef_all_noP$REL)
ef_all_noP <- ef_all_noP[ord_noP, ]

## save ef_all
saveRDS(object = ef_all_noP, file = './output_all/all_efSizes_noPvalue_temperature.RDS')
# ef_all_noP <- readRDS(file = './output_all/all_efSizes_noPvalue_temperature.RDS')

#    VII. Extract all effect sizes, also indirect ones, estimated with meta-analyses --------
## I have to extract the Coefs that were estimated within the meta-analysis, ie.
## 'Ind_DemRate<-det_Clim', 'Ind_GR<-det_Clim', 'Tot_DemRate<-det_Clim', 'Tot_GR<-det_Clim'
## and then merge that data with the Coefs_Aut and save on the disk

# test the function
two <- dplyr::bind_rows(lapply(list('Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'),
                               function(x){extr_coefs(obj = meta_Phen_Cov, Type_EfS = x)}))

## cycling through all the subsets and all estimated path coefficients to get them all into one dataset
efSizes_estim <- dplyr::bind_rows(lapply(list(meta_Phen_Cov, meta_Morph_Cov), function(x){
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
saveRDS(object = Coef_all_T, file = './output_all/PathCoefs_Temp_AlsoEstimatedRelations.RDS')
#check <- readRDS('./output_all/PathCoefs_Temp_AlsoEstimatedRelations.RDS')

