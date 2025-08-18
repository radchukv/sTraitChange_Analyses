## Script to run multiple meta-analyses on possible posterior
## phylogenies, to account for the uncertainty in phylogenetic trees

library(ggplot2)
library(tidyverse)
library(sTraitChange)
library(magrittr)
library(grafify)
# read in the results of SEMs on temperature
Coefs_Aut <- readRDS(file = './output_fSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')

## subset to only retain the studies with windur max 52
Coefs_Aut <- Coefs_Aut[Coefs_Aut$WinDur <= 51, ]
length(unique(Coefs_Aut$ID))  ## 202 (from initial 210)


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

Coefs_Aut_sp$WeathQ[Coefs_Aut_sp$WeathQ == ''] <-  "2"
Coefs_Aut_sp %<>% mutate(across(where(is_character), as_factor))


# 1. Prepare data for birds -----------------------------------------------

Coefs_Aut_Bird <- Coefs_Aut_sp %>%
  filter(Taxon == "Bird")


# read in the results of SEMs on precip
Coefs_Aut_precip <- readRDS(file = './output_fSEM_precip/PathCoefs_allMods_Precip_Weights_DD_Autocor.RDS')
## subset to only retain the studies with windur max 52
Coefs_Aut_precip <- Coefs_Aut_precip[Coefs_Aut_precip$WinDur <= 51, ]
length(unique(Coefs_Aut_precip$ID))  ## 207 (from initial 210)
unique(Coefs_Aut_precip$WeathQ)

Coefs_Aut_precip <- merge(Coefs_Aut_precip, traits_sub, by = 'Species', all.x = TRUE)
Coefs_Aut_precip %<>%
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

table(Coefs_Aut_precip$TraitType)
Coefs_Aut_precip %<>% mutate(across(where(is_character), as_factor))

Coefs_Aut_precip_Bird <- Coefs_Aut_precip %>%
  filter(Taxon == "Bird")


# 2. Bird analyses --------------------------------------------------------


# test the function
test_phenT <- fit_all_acountphylo(data_MA = Coefs_Aut_Bird, vtree_folder = "./data/phylogenies_100",
                              ind = 2, Trait_categ = "Phenological",
                            Clim = "Temperature")
test_morphT <- fit_all_acountphylo(data_MA = Coefs_Aut_Bird, vtree = "./data/phylogenies_100",
                            ind = 2, Trait_categ = "Morphological",
                            Clim = "Temperature")
# applying hte function to precip
test_phenP <- fit_all_acountphylo(data_MA = Coefs_Aut_precip_Bird, vtree = "./data/phylogenies_100",
                             ind = 2, Trait_categ = "Phenological",
                             Clim = "Precipitation")
test_morphP <- fit_all_acountphylo(data_MA = Coefs_Aut_precip_Bird, vtree = "./data/phylogenies_100",
                              ind = 2, Trait_categ = "Morphological",
                              Clim = "Precipitation")


ef_all <- rbind(test_phenT, test_morphT, test_phenP, test_morphP)


# apply the function to 100 selected trees
## !!!!!!!! ATTENTION: takes long time to run !!!!!!!!!!!!
## if desired, run overnight (on the machine with specs as in Readme)
## or on cluster
## Alternatively: read in the output file that is saved in
## output_all (see below)
start.time <- Sys.time()
print(start.time)
for(i in 2:100){

  print(i)
  phenT <- fit_all_acountphylo(data_MA = Coefs_Aut_Bird, vtree = "./data/phylogenies_100",
                                    ind = i, Trait_categ = "Phenological",
                                    Clim = "Temperature")
  morphT <- fit_all_acountphylo(data_MA = Coefs_Aut_Bird, vtree = "./data/phylogenies_100",
                                     ind = i, Trait_categ = "Morphological",
                                     Clim = "Temperature")
  phenP <- fit_all_acountphylo(data_MA = Coefs_Aut_precip_Bird, vtree = "./data/phylogenies_100",
                                    ind = i, Trait_categ = "Phenological",
                                    Clim = "Precipitation")
  morphP <- fit_all_acountphylo(data_MA = Coefs_Aut_precip_Bird, vtree = "./data/phylogenies_100",
                                     ind = i, Trait_categ = "Morphological",
                                     Clim = "Precipitation")
  ef_all <- rbind(ef_all, phenT, morphT, phenP, morphP)
end.time <- Sys.time()
print(end.time)
}


# save the file for the future use
saveRDS(ef_all, file = "./output_all/ef_all_difPhylo_Birds.rds")

# read the file in
ef_all <- readRDS("./output_all/ef_all_difPhylo_Birds.rds")


# from here on some cleaning will be needed

length(unique(ef_all$phylo)) # 100 replicates
# check in what proportion of those random phylogeny-based analyses the phylo-corrected model
# fits the data better than the model without accounting for phylogeny

ef_all_phylo <- ef_all %>%
  mutate(Phylobetter = ifelse(AIC_phylo < AIC_nophyl, "Yes", "No")) %>%
  filter(Variable == 'intrcpt') %>%
  group_by(Trait_Categ, Clim, Relation, phylo) %>%
  slice(1L) %>%
  ungroup() %>%
  mutate(lambda = Phylo.SD /(Species.SD + Phylo.SD))
100*2*2*6  # 2400 ok
nrow(ef_all_phylo)

nrow(ef_all_phylo[ef_all_phylo$Phylobetter == 'Yes', ]) /nrow(ef_all_phylo)  # 0.04208754
# in n0.04208754  of the 100 replicates the model with phylocorrection fits
# better than the model without accounting for phylogeny

# checking the lambda
ef_all_phylo %>%
  dplyr::filter(.,  Relation == 'Trait_mean<-det_Clim') %>%
ggplot(., aes(x = lambda)) + geom_histogram() +
  facet_grid(Trait_Categ ~ Clim, scales = 'free_y') + theme_bw()

ef_all_phylo %>%
  dplyr::filter(.,  Relation == 'GR<-Trait_mean') %>%
  ggplot(., aes(x = lambda)) + geom_histogram() +
  facet_grid(Trait_Categ ~ Clim, scales = 'free_y') + theme_bw()

ef_all_phylo %>%
  dplyr::filter(.,  Relation == 'Ind_GR<-det_Clim') %>%
  ggplot(., aes(x = lambda)) + geom_histogram() +
  facet_grid(Trait_Categ ~ Clim, scales = 'free_y') + theme_bw()

ef_all_phylo %>%
  dplyr::filter(.,  Relation == 'GR<-det_Clim') %>%
  ggplot(., aes(x = lambda)) + geom_histogram() +
  facet_grid(Trait_Categ ~ Clim, scales = 'free_y') + theme_bw()

# out of curiousity
ef_all_phylo %>%
  dplyr::filter(.,  Relation == 'GR<-Pop_mean') %>%
  ggplot(., aes(x = lambda)) + geom_histogram() +
  facet_grid(Trait_Categ ~ Clim, scales = 'free_y') + theme_bw()

quantile(ef_all_phylo$lambda, probs = c(0.05, 0.1, 0.2, 0.5, 0.7, 0.9, 0.95), na.rm = TRUE)


# 3. Prepare data for mammals -----------------------------------------------

Coefs_Aut_m <- Coefs_Aut_sp %>%
  filter(Taxon == "Mammal")


# subset the data on precip for reptiles only
Coefs_Aut_precip_m <- Coefs_Aut_precip %>%
  filter(Taxon == "Mammal")


# 4. Mammal analyses --------------------------------------------------------


# run for the first phylo treee
test_phenT_m <- fit_all_acountphylo(data_MA = Coefs_Aut_m, vtree_folder = "./data/phylogenies_100",
                                  ind = 2, Trait_categ = "Phenological",
                                  Clim = "Temperature")
test_morphT_m <- fit_all_acountphylo(data_MA = Coefs_Aut_m, vtree = "./data/phylogenies_100",
                                   ind = 2, Trait_categ = "Morphological",
                                   Clim = "Temperature")
# applying hte function to precip
test_phenP_m <- fit_all_acountphylo(data_MA = Coefs_Aut_precip_m, vtree = "./data/phylogenies_100",
                                  ind = 2, Trait_categ = "Phenological",
                                  Clim = "Precipitation")
test_morphP_m <- fit_all_acountphylo(data_MA = Coefs_Aut_precip_m, vtree = "./data/phylogenies_100",
                                   ind = 2, Trait_categ = "Morphological",
                                   Clim = "Precipitation")


ef_all_m <- rbind(test_phenT_m, test_morphT_m, test_phenP_m, test_morphP_m)


# apply the function to 100 selected trees
## !!!!!!!! ATTENTION: takes long time to run !!!!!!!!!!!!
## if desired, run overnight (on the machine with specs as in Readme)
## or on cluster
## Alternatively: read in the output file that is saved in
## output_all (see below)
start.time <- Sys.time()
print(start.time)
for(i in 2:100){

  print(i)
  phenT <- fit_all_acountphylo(data_MA = Coefs_Aut_m, vtree = "./data/phylogenies_100",
                               ind = i, Trait_categ = "Phenological",
                               Clim = "Temperature")
  morphT <- fit_all_acountphylo(data_MA = Coefs_Aut_m, vtree = "./data/phylogenies_100",
                                ind = i, Trait_categ = "Morphological",
                                Clim = "Temperature")
  phenP <- fit_all_acountphylo(data_MA = Coefs_Aut_precip_m, vtree = "./data/phylogenies_100",
                               ind = i, Trait_categ = "Phenological",
                               Clim = "Precipitation")
  morphP <- fit_all_acountphylo(data_MA = Coefs_Aut_precip_m, vtree = "./data/phylogenies_100",
                                ind = i, Trait_categ = "Morphological",
                                Clim = "Precipitation")
  ef_all_m <- rbind(ef_all_m, phenT, morphT, phenP, morphP)
  end.time <- Sys.time()
  print(end.time)
}


# save the file for the future use
saveRDS(ef_all_m, file = "./output_all/ef_all_difPhylo_Mams.rds")

# read the file in
ef_all_m <- readRDS("./output_all/ef_all_difPhylo_Mams.rds")


# from here on some cleaning will be needed

length(unique(ef_all_m$phylo)) # 100 replicates
# check in what proportion of those random phylogeny-based analyses the phylo-corrected model
# fits the data better than the model without accounting for phylogeny


# 5. Plots across all phylo analyses ------------------------------------------
# read the dtaa on all vertebrates
ef_all <- readRDS("./output_all/ef_all_difPhylo1.rds")



# keep 100 replicates only
ef_all_phylo <- ef_all %>%
  mutate(Phylobetter = ifelse(AIC_phylo < AIC_nophyl, "Yes", "No")) %>%
  filter(Variable == 'intrcpt') %>%
  group_by(Trait_Categ, Clim, Relation, phylo) %>%
  slice(1L) %>%
  ungroup() %>%
  mutate(lambda = Phylo.SD /(Species.SD + Phylo.SD),
         Analysis = 'full',
         REL = dplyr::case_when(Relation == 'Trait_mean<-det_Clim' ~ 'CZ',
                                Relation == 'GR<-Trait_mean' ~ 'ZG',
                                Relation == 'GR<-det_Clim' ~ 'CG',
                                Relation == 'GR<-Pop_mean' ~ 'PG',
                                Relation == 'Ind_GR<-det_Clim' ~ 'CZG',
                                Relation == 'Tot_GR<-det_Clim' ~ 'TotalCG'))
ef_all_phylo$REL <- factor(ef_all_phylo$REL, levels = c('CZ', 'ZG','CG', 'PG', 'CZG', 'TotalCG'))
100*2*2*6  # 2400 ok
nrow(ef_all_phylo)
nrow(ef_all_phylo[ef_all_phylo$Phylobetter == 'Yes', ]) /nrow(ef_all_phylo)  # 0.02609428


# read the file on birds
ef_all_b <- readRDS("./output_all/ef_all_difPhylo_Birds.rds")


length(unique(ef_all_b$phylo)) # 100 replicates
# check in what proportion of those random phylogeny-based analyses the phylo-corrected model
# fits the data better than the model without accounting for phylogeny

ef_all_phylo_b <- ef_all_b %>%
  mutate(Phylobetter = ifelse(AIC_phylo < AIC_nophyl, "Yes", "No")) %>%
  filter(Variable == 'intrcpt') %>%
  group_by(Trait_Categ, Clim, Relation, phylo) %>%
  slice(1L) %>%
  ungroup() %>%
  mutate(lambda = Phylo.SD /(Species.SD + Phylo.SD),
         Analysis = 'bird',
         REL = dplyr::case_when(Relation == 'Trait_mean<-det_Clim' ~ 'CZ',
                                Relation == 'GR<-Trait_mean' ~ 'ZG',
                                Relation == 'GR<-det_Clim' ~ 'CG',
                                Relation == 'GR<-Pop_mean' ~ 'PG',
                                Relation == 'Ind_GR<-det_Clim' ~ 'CZG',
                                Relation == 'Tot_GR<-det_Clim' ~ 'TotalCG'))
100*2*2*6  # 2400 ok
nrow(ef_all_phylo_b)

nrow(ef_all_phylo_b[ef_all_phylo_b$Phylobetter == 'Yes', ]) /nrow(ef_all_phylo_b)  # 0.04208754
ef_all_phylo_b$REL <- factor(ef_all_phylo_b$REL, levels = c('CZ', 'ZG','CG', 'PG', 'CZG', 'TotalCG'))

# read the data for mammals
# read the file in
ef_all_m <- readRDS("./output_all/ef_all_difPhylo_Mams.rds")
length(unique(ef_all_m$phylo)) # 100 replicates

ef_all_phylo_m <- ef_all_m %>%
  mutate(Phylobetter = ifelse(AIC_phylo < AIC_nophyl, "Yes", "No")) %>%
  filter(Variable == 'intrcpt') %>%
  group_by(Trait_Categ, Clim, Relation, phylo) %>%
  slice(1L) %>%
  ungroup() %>%
  mutate(lambda = Phylo.SD /(Species.SD + Phylo.SD),
         Analysis = 'mammal',
         REL = dplyr::case_when(Relation == 'Trait_mean<-det_Clim' ~ 'CZ',
                                Relation == 'GR<-Trait_mean' ~ 'ZG',
                                Relation == 'GR<-det_Clim' ~ 'CG',
                                Relation == 'GR<-Pop_mean' ~ 'PG',
                                Relation == 'Ind_GR<-det_Clim' ~ 'CZG',
                                Relation == 'Tot_GR<-det_Clim' ~ 'TotalCG'))
100*2*2*6  # 2400 ok
nrow(ef_all_phylo_m)

nrow(ef_all_phylo_m[ef_all_phylo_m$Phylobetter == 'Yes', ]) /nrow(ef_all_phylo_m)  #0.0008417508
ef_all_phylo_m$REL <- factor(ef_all_phylo_m$REL, levels = c('CZ', 'ZG','CG', 'PG', 'CZG', 'TotalCG'))



# plots
# 1. for the full dataset only
plot_lambda_tempPhen <- ef_all_phylo %>%
  dplyr::filter(.,  Trait_Categ == 'Phenological' & Clim == 'Temperature') %>%
  ggplot(., aes(x = lambda)) + geom_histogram() +
  facet_wrap(~ REL, scales = 'free_y') +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = rel(1.4)))

pdf('./plots_ms/FigS3_Hist_lambda_TempPhen_AllRelations.pdf', width = 10)
plot_lambda_tempPhen
dev.off()
# summary
sum_stats <- ef_all_phylo %>%
  dplyr::filter(.,  Trait_Categ == 'Phenological' & Clim == 'Temperature') %>%
  group_by(REL) %>%
  summarise(mean = mean(lambda, na.rm = TRUE), min = min(lambda, na.rm = TRUE),
            max = max(lambda, na.rm = TRUE), med = median(lambda, na.rm = TRUE))
sum_stats


# 2. for all datasets, on CZ only
ef <- rbind(ef_all_phylo, ef_all_phylo_b, ef_all_phylo_m)
# here still set up levels for Analysis
ef$Analysis <- factor(ef$Analysis, levels = c('bird', 'mammal', 'full'))

plot_CZ_AllTraitClim_analyses <- ef %>%
  dplyr::filter(.,  REL == 'CZ') %>%
  ggplot(., aes(x = lambda)) + geom_histogram(aes(fill = Analysis), alpha = 0.6) +
  geom_vline(
    data = . %>%
      group_by(Trait_Categ, Clim, Analysis) %>%
      summarise(mean = mean(lambda, na.rm = TRUE)),
    mapping = aes(xintercept = mean, col = Analysis),
     lwd = 1) +
  geom_text(data = . %>%
              group_by(Trait_Categ, Clim, Analysis) %>%
              summarise(mean = mean(lambda, na.rm = TRUE)),
            mapping = aes(x = mean +0.1, y = c(20, 5, 15, 5, 20, 15, 100, 70, 200, 100, 70, 20),
                          label = paste(round(mean, 2)), col = Analysis),
            size = 9) +
  facet_grid(Trait_Categ ~ Clim, scales = 'free_y') + theme_bw() +
  theme_bw() + ylab("Count") + xlab("Lambda") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = rel(1.1)),
        legend.position = 'bottom',
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1))) +
  scale_colour_grafify(palette = 'muted') +
  scale_fill_grafify(palette = 'muted')

pdf('./plots_ms/FigS22_Lambda_CZ_acrossAllsets.pdf')
print(plot_CZ_AllTraitClim_analyses)
dev.off()

# 3. for all datasets, on ZG only
plot_ZG_AllTraitClim_analyses <- ef %>%
  dplyr::filter(.,  REL == 'ZG') %>%
  ggplot(., aes(x = lambda)) + geom_histogram(aes(fill = Analysis), alpha = 0.6) +
  geom_vline(
    data = . %>%
      group_by(Trait_Categ, Clim, Analysis) %>%
      summarise(mean = mean(lambda, na.rm = TRUE)),
    mapping = aes(xintercept = mean, col = Analysis),
    lwd = 1) +
  geom_text(data = . %>%
              group_by(Trait_Categ, Clim, Analysis) %>%
              summarise(mean = mean(lambda, na.rm = TRUE)),
            mapping = aes(x = mean +0.1, y = c(25, 50, 75, 50, 25, 75, 50, 75, 100, 50, 75, 100),
                          label = paste(round(mean, 2)), col = Analysis),
            size = 9) +
  facet_grid(Trait_Categ ~ Clim, scales = 'free_y') + theme_bw() +
  theme_bw() + ylab("Count") + xlab("Lambda") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = rel(1.1)),
        legend.position = 'bottom',
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1))) +
  scale_colour_grafify(palette = 'muted') +
  scale_fill_grafify(palette = 'muted')


pdf('./plots_ms/FigS23_Lambda_ZG_acrossAllsets.pdf')
print(plot_ZG_AllTraitClim_analyses)
dev.off()


# 4. for all datasets, on CZG only
plot_CZG_AllTraitClim_analyses <- ef %>%
  dplyr::filter(.,  REL == 'CZG') %>%
  ggplot(., aes(x = lambda)) + geom_histogram(aes(fill = Analysis), alpha = 0.6) +
  geom_vline(
    data = . %>%
      group_by(Trait_Categ, Clim, Analysis) %>%
      summarise(mean = mean(lambda, na.rm = TRUE)),
    mapping = aes(xintercept = mean, col = Analysis),
    lwd = 1) +
  geom_text(data = . %>%
              group_by(Trait_Categ, Clim, Analysis) %>%
              summarise(mean = mean(lambda, na.rm = TRUE)),
            mapping = aes(x = mean +0.1, y = c(25, 50, 75, 50, 25, 75, 50, 75, 100, 50, 75, 100),
                          label = paste(round(mean, 2)), col = Analysis),
            size = 9) +
  facet_grid(Trait_Categ ~ Clim, scales = 'free_y') + theme_bw() +
  theme_bw() + ylab("Count") + xlab("Lambda") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = rel(1.1)),
        legend.position = 'bottom',
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1))) +
  scale_colour_grafify(palette = 'muted') +
  scale_fill_grafify(palette = 'muted')


pdf('./plots_ms/FigS24_Lambda_CZG_acrossAllsets.pdf')
print(plot_CZG_AllTraitClim_analyses)
dev.off()

# 5. for all datasets, on CG only
plot_CG_AllTraitClim_analyses <- ef %>%
  dplyr::filter(.,  REL == 'CG') %>%
  ggplot(., aes(x = lambda)) + geom_histogram(aes(fill = Analysis), alpha = 0.6) +
  geom_vline(
    data = . %>%
      group_by(Trait_Categ, Clim, Analysis) %>%
      summarise(mean = mean(lambda, na.rm = TRUE)),
    mapping = aes(xintercept = mean, col = Analysis),
    lwd = 1) +
  geom_text(data = . %>%
              group_by(Trait_Categ, Clim, Analysis) %>%
              summarise(mean = mean(lambda, na.rm = TRUE)),
            mapping = aes(x = mean +0.1, y = c(25, 50, 75, 50, 25, 75, 50, 120, 200, 50, 75, 100),
                          label = paste(round(mean, 2)), col = Analysis),
            size = 9) +
  facet_grid(Trait_Categ ~ Clim, scales = 'free_y') + theme_bw() +
  theme_bw() + ylab("Count") + xlab("Lambda") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = rel(1.1)),
        legend.position = 'bottom',
        axis.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1))) +
  scale_colour_grafify(palette = 'muted') +
  scale_fill_grafify(palette = 'muted')


pdf('./plots_ms/FigS25_Lambda_CG_acrossAllsets.pdf')
print(plot_CG_AllTraitClim_analyses)
dev.off()

# summarize pagles lambda for phen + temperature in ef (per dataset)
sum_stats_ef <- ef %>%
  dplyr::filter(.,  Trait_Categ == 'Phenological' & Clim == 'Temperature') %>%
  group_by(REL, Analysis) %>%
  summarise(mean = mean(lambda, na.rm = TRUE), min = min(lambda, na.rm = TRUE),
            max = max(lambda, na.rm = TRUE), med = median(lambda, na.rm = TRUE))
sum_stats_ef
