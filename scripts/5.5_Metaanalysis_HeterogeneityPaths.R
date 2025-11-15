# this script analyses heterogeneity in paths derived with SEMs that
# focused on the effects of temperature on phenological traits
# attempts to explain the heterogeneity in these paths by
# species characteristics and latitude,
# according to the hypotheses outlined in the manuscript

library(metafor)
library(tidyverse)
library(magrittr)
library(sTraitChange)
library(corrplot)
library(ape)
library(patchwork)

# 0. data read-in  and prepare --------------------------------------------
## here read in the files with all thepath coefficeints, also those for indirect and total paths
Coefs_Aut_prec <- readRDS(file = './output_all/PathCoefs_Precip_AlsoEstimatedRelations.RDS')
Coefs_Aut_prec$Climate <- 'Precipitation'

Coefs_Aut <- readRDS(file = './output_all/PathCoefs_Temp_AlsoEstimatedRelations.RDS')
Coefs_Aut$Climate <- 'Temperature'

Coefs <- rbind(Coefs_Aut, Coefs_Aut_prec)


## data on traits
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



traits_proc <- traits_proc %>%
  mutate(Sp_phylo = case_when(
    Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
    Species == 'Thalasseus sandvicensis' ~ 'Sterna sandvicensis',
    Species == 'Setophaga caerulescens' ~ 'Dendroica caerulescens',
    Species == 'Thalassarche melanophris' ~ 'Thalassarche melanophrys',
    Species == 'Ichthyaetus audouinii' ~ 'Larus audouinii',
    Species == 'Stercorarius maccormicki' ~ 'Catharacta maccormicki',
    TRUE ~ Species))


traits_proc$Species <- unlist(lapply(1:nrow(traits_proc), FUN = function(x){
  binary <- strsplit(as.character(traits_proc$Sp_phylo[x]), " ")
  Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")
}))


Coefs_Aut_sp <- merge(Coefs, subset(traits_proc, select = -Sp_phylo), by = 'Species')
Coefs_Aut_sp %<>%
  dplyr::mutate(., Blood = dplyr::case_when(
    Taxon %in% c('Bird', 'Mammal') ~ 'endotherm',
    Taxon %in% c('Reptilia', 'Fish') ~ 'ectotherm'
  ),
  abs_lat = abs(Latitude)) %>%
  mutate(across(where(is_character), as_factor))

# need a copy with the sp names, to be used to account for phylogeny (apart from that
# we also account for among-sp variation by including species as random intercept
# in the model)
Coefs_Aut_sp$Sp_phylo <- Coefs_Aut_sp$Species

# quick check of how many species have generation time > 2 years
# for the rebuttal letter
gen_time_check <- Coefs_Aut_sp %>%
  distinct(., Species, .keep_all = TRUE)
#%>% filter(GenLength_y_IUCN.x <= 3) %>%
# summarise(n())# 5 species

pdf('./plots_ms/Fig_GenTimeDistribution_rebutt.pdf', width = 9)
ggplot(gen_time_check, aes(x = GenLength_y_IUCN.x)) +
  geom_histogram() +
  geom_vline(
    data = . %>% summarise(med = median(GenLength_y_IUCN.x, na.rm = TRUE)),
    mapping = aes(xintercept = med),
    col = 'red', lwd = 2) + theme_bw() +
  xlab("Generation time, years") + theme(axis.title = element_text(size = 20))
dev.off()


# read in the phylo tree
vert_tree <- read.tree('./data/phylogenies_100/vert1.tre')
plot(vert_tree)
is.ultrametric(vert_tree)
# get the matrix of phylogenetic correlations
Mat_phylo <- ape::vcv.phylo(vert_tree, corr = TRUE)


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

# chen rossii is not in the dataframe
vert_tree$tip.label[! vert_tree$tip.label %in% Coefs_Aut_sp$Species]
# drop it from the tree
vert_tree <- drop.tip(vert_tree, 'Chen_rossii')

plot(vert_tree)
is.ultrametric(vert_tree)
# get the matrix of phylogenetic correlations
Mat_phylo <- ape::vcv.phylo(vert_tree, corr = TRUE)

# check again
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


# 1.  CZ  --------------------------------------------------------
CZ_Phen_T <- droplevels(Coefs_Aut_sp %>%
   dplyr::filter(., Relation == 'Trait_mean<-det_Clim' & Trait_Categ == 'Phenological' &
                   Climate == 'Temperature'))

hist(CZ_Phen_T$Estimate)

## am afraid we o not have that nice a spread across the latitude to test thi shypothesis
ggplot(CZ_Phen_T, aes(x = Latitude, y = Estimate)) + geom_point() +
  facet_grid(rows = vars(Climate))

ggplot(CZ_Phen_T, aes(x = Estimate)) + geom_histogram() +
  facet_grid(rows = vars(Climate), cols = vars(Taxon)) ## not much to see here
table(CZ_Phen_T$Taxon)
ggplot(CZ_Phen_T, aes(x = GenLength_y_IUCN.y, y = Estimate)) + geom_point() +
  facet_grid(rows = vars(Climate))

ggplot(CZ_Phen_T, aes(y = Estimate, x = Migrat)) + geom_boxplot() +
  facet_grid(rows = vars(Climate))
ggplot(CZ_Phen_T, aes(y = Estimate, x = Diet_HCO)) + geom_boxplot() +
  facet_grid(rows = vars(Climate))


# 1.1. CZ - overall model ---------------------------------------------

# 1. use the dataset where the rows with no data for migratory mode are dropped
# 2. center the quantitative predictors
CZ_PhenT_all <- droplevels(CZ_Phen_T %>%
  filter(! is.na(GenLength_y_IUCN.y) &  Migrat != 'unknown') %>%
  mutate(Hemisphere = if_else(Latitude > 0, 'N', 'S')))
table(CZ_PhenT_all$Blood)  # only endotherms! only 3 ectotherms are left, So we can't look at thermoreg as predictor

CZ_PhenT_scale <- CZ_PhenT_all %>%
  mutate(abs_lat = scale(abs_lat, scale = FALSE),
         GenLength_y_IUCN.y = scale(GenLength_y_IUCN.y, scale = FALSE),
         Pvalue = scale(Pvalue, scale = FALSE))




# including phylogeny
tre_sub <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% CZ_PhenT_scale$Species))
Mat_phen_climtr <- ape::vcv.phylo(tre_sub, corr = TRUE)
corrplot(Mat_phen_climtr)

# run the model
mod_CZ_PhenT_all <- rma.mv(Estimate ~ Diet_HCO + Migrat +
                             GenLength_y_IUCN.y + abs_lat + Pvalue,
                           V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                       ~1|Location, ~1|Sp_phylo),
                           R = list(Sp_phylo = Mat_phen_climtr),
                           data = CZ_PhenT_scale, method = 'ML')
mod_CZ_PhenT_all
anova(mod_CZ_PhenT_all, btt = c(2:7))  ##
# Test of Moderators (coefficients 2:7):
#  QM(df = 6) = 17.2202, p-val = 0.0085

tab_spSpecific_uni(mod_mv = mod_CZ_PhenT_all,
                   table_name = './tables/CZ_PhenT_allScaled',
                   explanatory = c('Diet_HCO', 'Migrat',
                                   'GenLength_y_IUCN.y', 'abs_lat', 'Pvalue'),
                   interact_fac = NULL)


# 1.2. CZ - univariate tests ------------------------------------------
## fitting per each single variable
## 1. Absolute Latitude

mod_CZ_PhenT_AbsLat <- rma.mv(Estimate ~ abs_lat + Pvalue,
                     V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                 ~1|Location, ~1|Sp_phylo),
                     R = list(Sp_phylo = Mat_phen_climtr),
                     data = CZ_PhenT_all, method = 'ML')
summary(mod_CZ_PhenT_AbsLat) ## marginally signif

tab_spSpecific_uni(mod_mv = mod_CZ_PhenT_AbsLat,
                   table_name = './tables/CZ_PhenT_AbsLat',
                   explanatory = c('abs_lat', 'Pvalue'),
                   interact_fac = NULL)
plot_CZ_PhenT_AbsLat <- plot_uni_spSpec(data_allEstim = CZ_PhenT_all,
                               mod_mv = mod_CZ_PhenT_AbsLat,
                               lOut = 10, xLab = 'Absolute latitude',
                               yLab = 'Effect of temperature on \n phenology (CZ)',
                               pdf_basename = './output_all/PlotCZ_PhenT_byAbsLat_SD',
                               byHemisphere = TRUE,
                               miny = min(CZ_PhenT_all$Estimate) - 0.1,
                               maxy = max(CZ_PhenT_all$Estimate) + 0.1)

plot_CZ_PhenT_AbsLat
## 2. Generation length
mod_CZ_PhenT_Gen <- rma.mv(Estimate ~ GenLength_y_IUCN.y + Pvalue,
                     V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                 ~1|Location,  ~1|Sp_phylo),
                     R = list(Sp_phylo = Mat_phen_climtr),
                     data = CZ_PhenT_all, method = 'ML')
summary(mod_CZ_PhenT_Gen)

tab_spSpecific_uni(mod_mv = mod_CZ_PhenT_Gen,
                   table_name = './tables/CZ_PhenT_GenTime',
                   explanatory = c('GenLength_y_IUCN.y'),
                   interact_fac = NULL)

## plotting predictions
plot_CZ_PhenT_Gen <- plot_uni_spSpec(data_allEstim = CZ_PhenT_all,
                                    mod_mv = mod_CZ_PhenT_Gen,
                                    lOut = 10, xLab = 'Generation time (years)',
                                    yLab = 'CZ estimate',
                                    pdf_basename = './output_all/PlotCZ_PhenT_byGen_SD',
                                    byHemisphere = FALSE)

plot_CZ_PhenT_Gen


## 3. Migratory mode
CZ_PhenT_migr <- CZ_Phen_T %>%
  dplyr::filter(Migrat != 'unknown')

# including phylogeny
tre_sub_migr <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% CZ_PhenT_migr$Species))
Mat_phen_climtr_migr <- ape::vcv.phylo(tre_sub_migr, corr = TRUE)
corrplot(Mat_phen_climtr_migr)

mod_CZ_PhenT_Migr <- rma.mv(Estimate ~ Migrat + Pvalue,
                     V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                 ~1|Location, ~1|Sp_phylo),
                     R = list(Sp_phylo = Mat_phen_climtr_migr),
                     data = CZ_PhenT_migr, method = 'ML')
summary(mod_CZ_PhenT_Migr)  ## nonsignif

stats::anova(mod_CZ_PhenT_Migr, btt = c(2))
# non-signif

tab_spSpecific_uni(mod_mv = mod_CZ_PhenT_Migr,
                   table_name = './tables/CZ_PhenT_Migrat',
                   explanatory = c('Migrat'),
                   interact_fac = NULL)


# 4. Diet
table(CZ_Phen_T$Diet_HCO)
tre_sub_diet <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% CZ_Phen_T$Species))
Mat_phen_climtr_d <- ape::vcv.phylo(tre_sub_diet, corr = TRUE)
corrplot(Mat_phen_climtr_d)
mod_CZ_PhenT_Diet <- rma.mv(Estimate ~ Diet_HCO + Pvalue,
                             V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                         ~1|Location, ~1|Sp_phylo),
                            R = list(Sp_phylo = Mat_phen_climtr_d),
                             data = CZ_Phen_T, method = 'ML')
summary(mod_CZ_PhenT_Diet)
stats::anova(mod_CZ_PhenT_Diet, btt = c(2,3))
## non-signif

tab_spSpecific_uni(mod_mv = mod_CZ_PhenT_Diet,
                   table_name = './tables/CZ_PhenT_Diet',
                   explanatory = c('Diet_HCO', 'Pvalue'),
                   interact_fac = NULL)


# plot
head(predict(mod_CZ_PhenT_Diet,  addx=TRUE))

PhenCZ_TDiet <- matrix(c(c(0, 1, 0), c(0, 0, 1), rep(mean(CZ_Phen_T$Pvalue, na.rm = TRUE), 3)),
                       ncol = 3, byrow = FALSE)
Pred_PhenCZ_TDiet <- as.data.frame(predict(mod_CZ_PhenT_Diet, newmods = PhenCZ_TDiet,  addx=TRUE))

Pred_PhenCZ_TDiet %<>%
  dplyr::rename(., Pvalue = X.Pvalue,
                Estimate = pred) %>%
  dplyr::mutate(., Diet_HCO = c("carnivore", "herbivore", "omnivore"),
                SD = (ci.ub - ci.lb) / 4,
                Est_PlSD = Estimate + SD,
                Est_MinSD = Estimate - SD)

## adding the counts in each category
counts_Phen_Diet <- as.data.frame(table(CZ_Phen_T$Diet_HCO))
counts_Phen_Diet$x <- c(0.5, 0.7, 1)
counts_Phen_Diet$Estimate <- rep(1, 2, 3)
counts_Phen_Diet %<>%
  dplyr::rename(Diet_HCO = Var1)

## plotting the predictions
plot_CZ_PhenT_Diet <- ggplot(CZ_Phen_T, aes(x = Diet_HCO, y = Estimate)) +
  geom_point(position = position_dodge(0.3), size = 3, alpha = 0.2,
             show.legend = FALSE) +
  geom_errorbar(data = Pred_PhenCZ_TDiet,
                aes(ymin = Est_MinSD, ymax = Est_PlSD), width = 0.2,
                position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_point(data = Pred_PhenCZ_TDiet,
             aes(x = Diet_HCO, y = Estimate),
             position = position_dodge(width = 0.3), alpha = 0.8,
             size = 5, pch = 21) +  # , show.legend = FALSE
  theme_bw() + xlab('Diet') + ylab('CZ Estimate') +
  geom_text(counts_Phen_Diet, mapping = aes(label = Freq),
            position = position_dodge(width = 0.3), fontface = 'bold',
            cex = rel(5), show.legend = FALSE)

pdf('./output_all/PlotCZ_PhenT_byDiet_SD.pdf')
plot_CZ_PhenT_Diet
dev.off()


# 2.  ZG  --------------------------------------------------------
ZG_Phen_T <- droplevels(Coefs_Aut_sp %>%
                          dplyr::filter(., Relation == 'GR<-Trait_mean' &
                                          Trait_Categ == 'Phenological' &
                                          Climate == 'Temperature'))

hist(ZG_Phen_T$Estimate)

## am afraid we o not have that nice a spread across the latitude to test thi shypothesis
ggplot(ZG_Phen_T, aes(x = abs_lat, y = Estimate)) + geom_point()
table(ZG_Phen_T$Taxon)

ggplot(ZG_Phen_T, aes(x = GenLength_y_IUCN.y, y = Estimate)) + geom_point()
ggplot(ZG_Phen_T, aes(y = Estimate, x = Migrat)) + geom_boxplot()

# 2.1. ZG - overall model ---------------------------------------------

# 1. use the dataset where the rows with no data for migratory mode are dropped
# 2. center the quantitative predictors
ZG_PhenT_all <- droplevels(ZG_Phen_T %>%
                             filter(! is.na(GenLength_y_IUCN.y) &  Migrat != 'unknown') %>%
                             mutate(Hemisphere = if_else(Latitude > 0, 'N', 'S')))
table(ZG_PhenT_all$Blood)  # only endotherms! So we can't look at thermoreg as predictor

ZG_PhenT_scale <- ZG_PhenT_all %>%
  mutate(abs_lat = scale(abs_lat, scale = FALSE),
         GenLength_y_IUCN.y = scale(GenLength_y_IUCN.y, scale = FALSE),
         Pvalue = scale(Pvalue, scale = FALSE))

# including phylogeny
tre_sub <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% ZG_PhenT_scale$Species))
Mat_phen_climtr <- ape::vcv.phylo(tre_sub, corr = TRUE)
corrplot(Mat_phen_climtr)


# here Pvalue is not needed as a covariate
mod_ZG_PhenT_all <- rma.mv(Estimate ~ Diet_HCO + Migrat +
                             GenLength_y_IUCN.y + abs_lat,
                           V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                       ~1|Location, ~1|Sp_phylo),
                           R = list(Sp_phylo = Mat_phen_climtr),
                           data = ZG_PhenT_scale, method = 'ML',
                           control = list(optimizer = 'uobyqa'))
mod_ZG_PhenT_all
anova(mod_ZG_PhenT_all, btt = c(2:6))  ##
# Test of Moderators (coefficients 2:6):
# QM(df = 5) = 9.1296, p-val = 0.1040


tab_spSpecific_uni(mod_mv = mod_ZG_PhenT_all,
                   table_name = './tables/ZG_PhenT_allScaled',
                   explanatory = c('Diet_HCO', 'Migrat',
                                   'GenLength_y_IUCN.y', 'abs_lat'),
                   interact_fac = NULL)


# 2.2. ZG - univariate tests --------------------------------------------
## fitting per each single variable
## 1.1. Absolute Latitude
mod_ZG_PhenT_AbsLat <- rma.mv(Estimate ~ abs_lat + Pvalue,
                              V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                          ~1|Location, ~1|Sp_phylo),
                              R = list(Sp_phylo = Mat_phen_climtr),
                              data = ZG_PhenT_all, method = 'ML')
summary(mod_ZG_PhenT_AbsLat) ## non-signif

tab_spSpecific_uni(mod_mv = mod_ZG_PhenT_AbsLat,
                   table_name = './tables/ZG_PhenT_AbsLat',
                   explanatory = c('abs_lat', 'Pvalue'),
                   interact_fac = NULL)
plot_ZG_PhenT_AbsLat <- plot_uni_spSpec(data_allEstim = ZG_PhenT_all,
                                                mod_mv = mod_ZG_PhenT_AbsLat,
                                                lOut = 10, xLab = 'Absolute latitude',
                                                yLab = 'Effect of phenology on\n population growth rate (ZG)',
                                                pdf_basename = './output_all/PlotZG_PhenT_byAbsLat_SD',
                                          byHemisphere = TRUE,
                                        miny = min(CZ_PhenT_all$Estimate) - 0.1,
                                        maxy = max(CZ_PhenT_all$Estimate) + 0.1)

plot_ZG_PhenT_AbsLat
## 2. Generation length
mod_ZG_PhenT_Gen <- rma.mv(Estimate ~ GenLength_y_IUCN.y + Pvalue,
                           V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                       ~1|Location, ~1|Sp_phylo),
                           R = list(Sp_phylo = Mat_phen_climtr),
                           data = ZG_PhenT_all, method = 'ML')
summary(mod_ZG_PhenT_Gen) ## marginal

tab_spSpecific_uni(mod_mv = mod_ZG_PhenT_Gen,
                   table_name = './tables/ZG_PhenT_GenTime',
                   explanatory = c('GenLength_y_IUCN.y', 'Pvalue'),
                   interact_fac = NULL)

## 3. Migratory mode
ZG_PhenT_migr <- droplevels(ZG_PhenT_all %>%
  dplyr::filter(Migrat != 'unknown'))

# including phylogeny
tre_sub_migr <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% ZG_PhenT_migr$Species))
Mat_phen_climtr_migr <- ape::vcv.phylo(tre_sub_migr, corr = TRUE)
corrplot(Mat_phen_climtr_migr)

mod_ZG_PhenT_Migr <- rma.mv(Estimate ~ Migrat + Pvalue,
                            V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                        ~1|Location, ~1|Sp_phylo),
                            R = list(Sp_phylo = Mat_phen_climtr_migr),
                            data = ZG_PhenT_migr, method = 'ML')
summary(mod_ZG_PhenT_Migr)  ## nonsignif


tab_spSpecific_uni(mod_mv = mod_ZG_PhenT_Migr,
                   table_name = './tables/ZG_PhenT_Migrat',
                   explanatory = c('Migrat', 'Pvalue'),
                   interact_fac = NULL)


## 4. Thermoreg
table(ZG_PhenT_all$Blood) # only 3 ectotherms vs 86 enodtherms... does not make sense

## 5. Diet
table(ZG_PhenT_all$Diet_HCO)
mod_ZG_PhenT_Diet <- rma.mv(Estimate ~ Diet_HCO + Pvalue,
                            V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                        ~1|Location, ~1|Sp_phylo),
                            R = list(Sp_phylo = Mat_phen_climtr),
                            data = ZG_PhenT_all, method = 'ML')
summary(mod_ZG_PhenT_Diet)
stats::anova(mod_ZG_PhenT_Diet, btt = c(2,3))
## non-signif

tab_spSpecific_uni(mod_mv = mod_ZG_PhenT_Diet,
                   table_name = './tables/ZG_PhenT_Diet',
                   explanatory = c('Diet_HCO', 'Pvalue'),
                   interact_fac = NULL)


# 3.  CZG --------------------------------------------------------
CZG_Phen_T <- droplevels(Coefs_Aut_sp %>%
                        dplyr::filter(., Relation == 'Ind_GR<-det_Clim' &
                                        Trait_Categ == 'Phenological' &
                                        Climate == 'Temperature'))

hist(CZG_Phen_T$Estimate)
ggplot(CZG_Phen_T, aes(x = Latitude, y = Estimate)) + geom_point()

ggplot(CZG_Phen_T, aes(x = GenLength_y_IUCN.y, y = Estimate)) + geom_point()
ggplot(CZG_Phen_T, aes(y = Estimate, x = Migrat)) + geom_boxplot()

# 3.1. CZG - overall model ---------------------------------------------

# 1. use the dataset where the rows with no data for migratory mode are dropped
# 2. center the quantitative predictors
CZG_PhenT_all <- droplevels(CZG_Phen_T %>%
                             filter(! is.na(GenLength_y_IUCN.y) &  Migrat != 'unknown') %>%
                              mutate(Hemisphere = if_else(Latitude > 0, 'N', 'S')))
table(CZG_PhenT_all$Blood)  # only endotherms! So we can't look at thermoreg as predictor

CZG_PhenT_scale <- CZG_PhenT_all %>%
  mutate(abs_lat = scale(abs_lat, scale = FALSE),
         GenLength_y_IUCN.y = scale(GenLength_y_IUCN.y, scale = FALSE),
         Pvalue = scale(Pvalue, scale = FALSE))


# including phylogeny
tre_sub <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% CZG_PhenT_scale$Species))
Mat_phen_climtr <- ape::vcv.phylo(tre_sub, corr = TRUE)
corrplot(Mat_phen_climtr)

mod_CZG_PhenT_all <- rma.mv(Estimate ~ Diet_HCO + Migrat +
                             GenLength_y_IUCN.y + abs_lat + Pvalue,
                           V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                       ~1|Location, ~1|Sp_phylo),
                           R = list(Sp_phylo = Mat_phen_climtr),
                           data = CZG_PhenT_scale, method = 'ML',
                           control = list(optimizer = 'uobyqa'))
mod_CZG_PhenT_all
anova(mod_CZG_PhenT_all, btt = c(2:7))  ##
# Test of Moderators (coefficients 2:7):
# QM(df = 6) = 6.1427, p-val = 0.4074

tab_spSpecific_uni(mod_mv = mod_CZG_PhenT_all,
                   table_name = './tables/CZG_PhenT_allScaled',
                   explanatory = c('Diet_HCO', 'Migrat',
                                   'GenLength_y_IUCN.y', 'abs_lat', 'Pvalue'),
                   interact_fac = NULL)



# 3.2. CZG - univariate tests -------------------------------------------

## 1. Absolute Latitude
mod_CZG_PhenT_AbsLat <- rma.mv(Estimate ~ abs_lat + Pvalue,
                        V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                    ~1|Location, ~1|Sp_phylo),
                        R = list(Sp_phylo = Mat_phen_climtr),
                        data = CZG_PhenT_all, method = 'ML')
summary(mod_CZG_PhenT_AbsLat)  # non-signif

tab_spSpecific_uni(mod_mv = mod_CZG_PhenT_AbsLat, table_name = './tables/CZG_PhenT_AbsLat',
                   explanatory = c('abs_lat'),
                   interact_fac = NULL)

##  plot
plot_CZG_PhenT_AbsLat <- plot_uni_spSpec(data_allEstim = CZG_PhenT_all,
                                        mod_mv = mod_CZG_PhenT_AbsLat,
                                        lOut = 10, xLab = 'Absolute latitude',
                                        yLab = 'Phenology-mediated effect of temperature \n on population growth rate (CZG)',
                                        pdf_basename = './output_all/PlotCZG_PhenT_byAbsLat',
                                        byHemisphere = TRUE,
                                        miny = min(CZ_PhenT_all$Estimate) - 0.1,
                                        maxy = max(CZ_PhenT_all$Estimate) + 0.1)
plot_CZG_PhenT_AbsLat

## 2. Generation length
mod_CZG_PhenT_Gen <- rma.mv(Estimate ~ GenLength_y_IUCN.y + Pvalue,
                     V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                 ~1|Location, ~1|Sp_phylo),
                     R = list(Sp_phylo = Mat_phen_climtr),
                     data = CZG_PhenT_all, method = 'ML')
summary(mod_CZG_PhenT_Gen)

tab_spSpecific_uni(mod_mv = mod_CZG_PhenT_Gen,
                   table_name = './tables/CZG_PhenT_GenTime',
                   explanatory = c('GenLength_y_IUCN.y'),
                   interact_fac = NULL)

## 3. Migratory mode
CZG_PhenT_migr <- CZG_PhenT_all %>%
  dplyr::filter(Migrat != 'unknown')
# including phylogeny
tre_sub_migr <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% CZG_PhenT_migr$Species))
Mat_phen_climtr_migr <- ape::vcv.phylo(tre_sub_migr, corr = TRUE)
corrplot(Mat_phen_climtr_migr)

mod_CZG_PhenT_Migr <- rma.mv(Estimate ~ Migrat + Pvalue,
                      V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                  ~1|Location, ~1|Sp_phylo),
                      R = list(Sp_phylo = Mat_phen_climtr_migr),
                      data = CZG_PhenT_migr, method = 'ML')
summary(mod_CZG_PhenT_Migr)

tab_spSpecific_uni(mod_mv = mod_CZG_PhenT_Migr,
                   table_name = './tables/CZG_PhenT_Migrat',
                   explanatory = c('Migrat'),
                   interact_fac = NULL)


## 4. Thermoreg
table(CZG_PhenT_migr$Blood)  ## only 2 ectotherms, is not very meaningful

## 5. Diet
table(CZG_PhenT_all$Diet_HCO)
mod_CZG_PhenT_Diet <- rma.mv(Estimate ~ Diet_HCO + Pvalue,
                            V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                        ~1|Location, ~1|Sp_phylo),
                            R = list(Sp_phylo = Mat_phen_climtr),
                            data = CZG_PhenT_all, method = 'ML')
summary(mod_CZG_PhenT_Diet)  ## p-val = 0.2679
stats::anova(mod_CZG_PhenT_Diet, btt = c(2,3))
## non-signif

tab_spSpecific_uni(mod_mv = mod_CZG_PhenT_Diet,
                   table_name = './tables/CZG_PhenT_Diet',
                   explanatory = c('Diet_HCO', 'Pvalue'),
                   interact_fac = NULL)


# 4.  CG  --------------------------------------------------------
CG_Phen_T <- droplevels(Coefs_Aut_sp %>%
                         dplyr::filter(., Relation == 'GR<-det_Clim'
                                       & Trait_Categ == 'Phenological' &
                                         Climate == 'Temperature'))

hist(CG_Phen_T$Estimate)
ggplot(CG_Phen_T, aes(x = Latitude, y = Estimate)) + geom_point()
ggplot(CG_Phen_T, aes(x = GenLength_y_IUCN.y, y = Estimate)) + geom_point()
ggplot(CG_Phen_T, aes(y = Estimate, x = Migrat)) + geom_boxplot()

# 4.1. CG - overall model ---------------------------------------------

# 1. use the dataset where the rows with no data for migratory mode are dropped
# 2. center the quantitative predictors
CG_PhenT_all <- droplevels(CG_Phen_T %>%
                              filter(! is.na(GenLength_y_IUCN.y) &  Migrat != 'unknown') %>%
                             mutate(Hemisphere = if_else(Latitude > 0, 'N', 'S')))
table(CG_PhenT_all$Blood)  # only endotherms! So we can't look at thermoreg as predictor

CG_PhenT_scale <- CG_PhenT_all %>%
  mutate(abs_lat = scale(abs_lat, scale = FALSE),
         GenLength_y_IUCN.y = scale(GenLength_y_IUCN.y, scale = FALSE),
         Pvalue = scale(Pvalue, scale = FALSE))

# including phylogeny
tre_sub <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% CG_PhenT_scale$Species))
Mat_phen_climtr <- ape::vcv.phylo(tre_sub, corr = TRUE)
corrplot(Mat_phen_climtr)

# run the model
mod_CG_PhenT_all <- rma.mv(Estimate ~ Diet_HCO + Migrat +
                              GenLength_y_IUCN.y + abs_lat + Pvalue,
                            V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                        ~1|Location, ~1|Sp_phylo),
                           R = list(Sp_phylo = Mat_phen_climtr),
                            data = CG_PhenT_scale, method = 'ML',
                            control = list(optimizer = 'uobyqa'))
mod_CG_PhenT_all
anova(mod_CG_PhenT_all, btt = c(2:7))  ##
# Test of Moderators (coefficients 2:7):
# QM(df = 6) = 21.2497, p-val = 0.0017


tab_spSpecific_uni(mod_mv = mod_CG_PhenT_all,
                   table_name = './tables/CG_PhenT_allScaled',
                   explanatory = c('Diet_HCO', 'Migrat',
                                   'GenLength_y_IUCN.y', 'abs_lat', 'Pvalue'),
                   interact_fac = NULL)
# here again absolute latitude is significant

# 4.2. CG - univariate tests -------------------------------------------
## fitting a model with one explanatory at a time

## 1.1. Absolute Latitude
mod_CG_PhenT_AbsLat <- rma.mv(Estimate ~ abs_lat + Pvalue,
                         V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                     ~1|Location, ~1|Sp_phylo),
                         R = list(Sp_phylo = Mat_phen_climtr),
                         data = CG_PhenT_all, method = 'ML',
                         control = list(optimizer = 'uobyqa'))
summary(mod_CG_PhenT_AbsLat)

tab_spSpecific_uni(mod_mv = mod_CG_PhenT_AbsLat,
                   table_name = './tables/CG_PhenT_AbsLat',
                   explanatory = c('abs_lat', 'Pvalue'),
                   interact_fac = NULL)
## signif. latitude again

## checking the plot
plot_CG_PhenT_AbsLat <- plot_uni_spSpec(data_allEstim = CG_PhenT_all,
                                                mod_mv = mod_CG_PhenT_AbsLat,
                                                lOut = 10,
                                                xLab = 'Absolute latitude',
                                                yLab = 'Direct effect of temperature on \n population growth rate (CG)',
                                                pdf_basename = './output_all/PlotCG_PhenT_byAbsLat',
                                        byHemisphere = TRUE,
                                        miny = min(CZ_PhenT_all$Estimate) - 0.1,
                                        maxy = max(CZ_PhenT_all$Estimate) + 0.1)

plot_CG_PhenT_AbsLat
## 2. Generation time
mod_CG_PhenT_Gen <- rma.mv(Estimate ~ GenLength_y_IUCN.y + Pvalue,
                      V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                  ~1|Location, ~1|Sp_phylo),
                      R = list(Sp_phylo = Mat_phen_climtr),
                      data = CG_PhenT_all, method = 'ML',
                      control = list(optimizer = 'BFGS'))
summary(mod_CG_PhenT_Gen)
## non-signif

tab_spSpecific_uni(mod_mv = mod_CG_PhenT_Gen,
                   table_name = './tables/CG_PhenT_GenTime',
                   explanatory = c('GenLength_y_IUCN.y', 'Pvalue'),
                   interact_fac = NULL)


## 3. Migratory mode
CG_PhenT_migr <- droplevels(CG_PhenT_all %>%
  dplyr::filter(Migrat != 'unknown'))
table(CG_PhenT_migr$Migrat)

# including phylogeny
tre_sub_migr <- drop.tip(vert_tree, which(!vert_tree$tip.label %in% CG_PhenT_migr$Species))
Mat_phen_climtr_migr <- ape::vcv.phylo(tre_sub_migr, corr = TRUE)
corrplot(Mat_phen_climtr_migr)

mod_CG_PhenT_Migr <- rma.mv(Estimate ~ Migrat + Pvalue,
                       V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                   ~1|Location, ~1|Sp_phylo),
                       R = list(Sp_phylo = Mat_phen_climtr_migr),
                       data = CG_PhenT_migr, method = 'ML',
                       control = list(optimizer = 'BFGS'))
summary(mod_CG_PhenT_Migr)

tab_spSpecific_uni(mod_mv = mod_CG_PhenT_Migr,
                   table_name = './tables/CG_PhenT_Migrat',
                   explanatory = c('Migrat', 'Pvalue'),
                   interact_fac = NULL)

## 4. Thermoreg
table(CG_PhenT_all$Blood) # only 3 ectothemrs, won't work


## 5. Diet
table(CG_PhenT_all$Diet_HCO)
mod_CG_PhenT_Diet <- rma.mv(Estimate ~ Diet_HCO + Pvalue,
                             V = SError^2, random = list(~ 1|Species, ~1|ID,
                                                         ~1|Location, ~1|Sp_phylo),
                            R = list(Sp_phylo = Mat_phen_climtr),
                             data = CG_PhenT_all, method = 'ML')
summary(mod_CG_PhenT_Diet)  # non-signif
stats::anova(mod_CG_PhenT_Diet, btt = c(2,3))
## non-signif

tab_spSpecific_uni(mod_mv = mod_CG_PhenT_Diet,
                   table_name = './tables/CG_PhenT_Diet',
                   explanatory = c('Diet_HCO', 'Pvalue'),
                   interact_fac = NULL)


# 5.  Overall Plots -------------------------------------------------------

# Figure with two panels onlY: effect of latitude and generation time
pdf('./plots_ms/FigS6_3Panels_GenT_Diet&AbsLat_SD.pdf',
    width = 12, height = 5)
plot_CZ_PhenT_AbsLat[[2]]  + plot_CZ_PhenT_Gen[[2]] +
  plot_CZ_PhenT_Diet +
  plot_layout(guides = 'collect', ncol = 3, nrow = 1) +
  plot_annotation(tag_levels = "a") &
  theme(plot.margin = unit(c(0.3, 0.3, 0, 0), 'line'),
        plot.tag = element_text(face = 'bold'),
        legend.position = 'bottom')
dev.off()

# plot of the associations between abs_latitude and each of the path coefficients in the SEM
pdf('./plots_ms/Fig4_Latitude_explains_CZ_CZG_CG_based_SD.pdf',
    width = 8, height = 8)
plot_CZ_PhenT_AbsLat +
  plot_ZG_PhenT_AbsLat +
  plot_CZG_PhenT_AbsLat +
  plot_CG_PhenT_AbsLat +
  plot_layout(guides = 'collect', axis_titles = "collect", ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = "a") &
  theme(plot.margin = unit(c(0.3, 0.3, 0, 0), 'line'),
        plot.tag = element_text(face = 'bold'),
        legend.position = 'bottom')
dev.off()
