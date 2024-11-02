## Script with meta-analysis of the SEM results on temperature
## SEM models fitted without the effect of population size on pop
# growth rate
## accounting for phylogeny
## using a subset of SEM studies with winDur limited to max 52 weeks
## including covariates only in the needed models (i.e.
## where response is related to climate, and hence to the climatic window)

library(metafor)
library(vegan)
library(ggplot2)
library(magrittr)
library(ape)
library(corrplot)



## read in the data
Coefs_Aut <- readRDS(file = './output_fSEM_noDD_temp/PathCoefs_allMods_Temp_Weights_noDD_Autocor.RDS')
table(Coefs_Aut$Taxon)

## subset to only retain the studies with windur max 52
Coefs_Aut <- Coefs_Aut[Coefs_Aut$WinDur <= 51, ]
length(unique(Coefs_Aut$ID))  ## 205

## read in the trait data and merge
traits <- read.csv('./data/speciesTraits.csv')

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
# best is probably to not turn itno ultrametric but then calculate proper correlation
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



## fitting the meta-analyses
# Always have to prepare the subtree for a specific combi of trait and clim!! (here done in the test)
meta_Phen_Cov <- fit_all_meta(data_MA = Coefs_Aut_sp,
                              Demog_rate = NULL,
                              Trait_categ = 'Phenological',
                              Clim = 'Temperature',
                              Cov_fact = 'WeathQ',
                              COV = 'Pvalue',
                              sel = 'Temp_Phen_Cov_noDD',
                              folder_name = './output_all_noDD/',
                              colr = c('black', 'red'),
                              DD = 'none',
                              simpleSEM = TRUE,
                              A = Mat_phen_climtr,
                              all_Relations = c('Trait_mean<-det_Clim',
                                                'Ind_GR<-det_Clim',
                                                'Tot_GR<-det_Clim'))
meta_Phen <- fit_all_meta(data_MA = Coefs_Aut_sp,
                          Demog_rate = NULL,
                          Trait_categ = 'Phenological',
                          Clim = 'Temperature',
                          Cov_fact = NULL,
                          COV = NULL,
                          sel = 'Temp_Phen_noDD',
                          A = Mat_phen_climtr,
                          folder_name = './output_all_noDD/',
                          colr = c('black'),
                          DD = 'none',
                          simpleSEM = TRUE,
                          all_Relations = c('GR<-det_Clim',
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
                               sel = 'Temp_Morph_Cov_noDD',
                               A = Mat_morph_climtr,
                               folder_name = './output_all_noDD/',
                               colr = c('black', 'red'),
                               DD = 'none',
                               simpleSEM = TRUE,
                               all_Relations = c('Trait_mean<-det_Clim',
                                                 'Ind_GR<-det_Clim',
                                                 'Tot_GR<-det_Clim'))


meta_Morph <- fit_all_meta(data_MA = Coefs_Aut_sp,
                           Demog_rate = NULL,
                           Trait_categ = 'Morphological',
                           Clim = 'Temperature',
                           Cov_fact = NULL,
                           COV = NULL,
                           sel = 'Temp_Morph_noDD',
                           A = Mat_morph_climtr,
                           folder_name = './output_all_noDD/',
                           colr = c('black'),
                           DD = 'none',
                           simpleSEM = TRUE,
                           all_Relations = c('GR<-det_Clim',
                                             'GR<-Trait_mean'))


# III. Summary across groups      -----------------------------------------

## binding all rows together
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



## binding all across-studies effect sizes together
ef_Phen <- rbind(meta_Phen$meta_res[[1]], meta_Phen_Cov$meta_res[[1]]) %>%
  dplyr::mutate(., Trait_Categ = 'Phenology')

ef_Morph <- rbind(meta_Morph$meta_res[[1]], meta_Morph_Cov$meta_res[[1]]) %>%
  dplyr::mutate(., Trait_Categ = 'Morphology')

ef_all <- rbind(ef_Phen, ef_Morph)
ef_all %<>%
  dplyr::mutate(., col = dplyr::case_when(EfS_Low < 0 & EfS_Upper < 0 ~ 'orange',
                                          EfS_Low > 0 & EfS_Upper > 0 ~ 'lightblue',
                                          TRUE ~ 'darkgrey'),
                col_var = dplyr::case_when(Variable == 'intrcpt' & EfS_Low < 0 & EfS_Upper < 0 ~ 'Intercept significant',
                                           Variable == 'intrcpt' & EfS_Low > 0 & EfS_Upper > 0 ~ 'Intercept significant',
                                           (Variable == 'intrcpt' & EfS_Low > 0 & EfS_Upper < 0) | (Variable == 'intrcpt' & EfS_Low <0 & EfS_Upper > 0) ~ 'Intercept non-significant',
                                           Variable == 'Pvalue' & EfS_Low < 0 & EfS_Upper < 0 ~ 'Pvalue significant',
                                           Variable == 'Pvalue' & EfS_Low > 0 & EfS_Upper > 0 ~ 'Pvalue significant',
                                           (Variable == 'Pvalue' & EfS_Low > 0 & EfS_Upper < 0) | (Variable == 'Pvalue' &  EfS_Low < 0 & EfS_Upper > 0) ~ 'Pvalue non-significant',
                                           Variable == 'WeathQ2' & EfS_Low < 0 & EfS_Upper < 0 ~ 'Weather quality significant',
                                           Variable == 'WeathQ2' & EfS_Low > 0 & EfS_Upper > 0 ~ 'Weather quality significant',
                                           (Variable == 'WeathQ2' & EfS_Low > 0 & EfS_Upper < 0) | ( Variable == 'WeathQ2' & EfS_Low < 0 & EfS_Upper > 0) ~ 'Weather quality non-significant'),
                REL = dplyr::case_when(Relation == 'Trait_mean<-det_Clim' ~ 'CZ',
                                       Relation == 'GR<-Trait_mean' ~ 'ZG',
                                       Relation == 'GR<-det_Clim' ~ 'CG',
                                       Relation == 'Ind_GR<-det_Clim' ~ 'CZG',
                                       Relation == 'Tot_GR<-det_Clim' ~ 'TotalCG'))


## sort the factor so as to have the coefficients from bottom to the top of the scheme
ef_all$REL <- factor(ef_all$REL, levels = c('CZ', 'ZG','CG', 'CZG', 'TotalCG'))

## sort the data frme to have the levels in the desired order
ord <- order(ef_all$Trait_Categ, ef_all$REL)
ef_all <- ef_all[ord, ]

# ef_all$yAx <- rep(c(1:11), 2) ## did not run it here, so may need to run it before plotting in the Fig. script

## save ef_all
saveRDS(object = ef_all, file = './output_all_noDD/all_efSizes_noDD_temp.RDS')



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
saveRDS(object = Coef_all_T, file = './output_all_noDD/PathCoefs_noDD_Temp_AlsoEstimatedRelations.RDS')


