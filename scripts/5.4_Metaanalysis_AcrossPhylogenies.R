## Script to run multiple meta-analyses on possible posterior
## phylogenies, to account for the uncertainty in phylogenetic trees

library(ggplot2)
library(tidyverse)
library(sTraitChange)
library(magrittr)

# read in the results of SEMs on temperature
Coefs_Aut <- readRDS(file = './output_fSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')

## subset to only retain the studies with windur max 52
Coefs_Aut <- Coefs_Aut[Coefs_Aut$WinDur <= 51, ]
length(unique(Coefs_Aut$ID))  ## 202 (from initial 210)


## read in the trait data and merge
traits <- read.csv('./data/speciesTraits.csv')
traits_sub <- subset(traits, select = c(Species, GenLength_y_IUCN))
Coefs_Aut_sp <- merge(Coefs_Aut, traits_sub, by = 'Species', all.x = TRUE)
Coefs_Aut_sp$WeathQ[Coefs_Aut_sp$WeathQ == ''] <-  "2"
Coefs_Aut_sp %<>% mutate(across(where(is_character), as_factor))



# read in the results of SEMs on precip
Coefs_Aut_precip <- readRDS(file = './output_fSEM_precip/PathCoefs_allMods_Precip_Weights_DD_Autocor.RDS')
## subset to only retain the studies with windur max 52
Coefs_Aut_precip <- Coefs_Aut_precip[Coefs_Aut_precip$WinDur <= 51, ]
length(unique(Coefs_Aut_precip$ID))  ## 207 (from initial 210)
unique(Coefs_Aut_precip$WeathQ)

Coefs_Aut_precip <- merge(Coefs_Aut_precip, traits_sub, by = 'Species', all.x = TRUE)
Coefs_Aut_precip %<>% mutate(across(where(is_character), as_factor))


# test the function
test_phenT <- fit_all_acountphylo(data_MA = Coefs_Aut_sp, vtree_folder = "./data/phylogenies_100",
                              ind = 2, Trait_categ = "Phenological",
                            Clim = "Temperature")
test_morphT <- fit_all_acountphylo(data_MA = Coefs_Aut_sp, vtree = "./data/phylogenies_100",
                            ind = 2, Trait_categ = "Morphological",
                            Clim = "Temperature")
# applying hte function to precip
test_phenP <- fit_all_acountphylo(data_MA = Coefs_Aut_precip, vtree = "./data/phylogenies_100",
                             ind = 2, Trait_categ = "Phenological",
                             Clim = "Precipitation")
test_morphP <- fit_all_acountphylo(data_MA = Coefs_Aut_precip, vtree = "./data/phylogenies_100",
                              ind = 2, Trait_categ = "Morphological",
                              Clim = "Precipitation")


ef_all <- rbind(test_phenT, test_morphT, test_phenP, test_morphP)


# apply the function to 107 selected trees
## !!!!!!!! ATTENTION: takes long time to run !!!!!!!!!!!!
## if desired, run overnight (on the machine with specs as in Readme)
## or on cluster
## Alternatively: read in the output file that is saved in
## output_all (see below)
start.time <- Sys.time()
print(start.time)
for(i in 2:108){

  print(i)
  phenT <- fit_all_acountphylo(data_MA = Coefs_Aut_sp, vtree = "./data/phylogenies_100",
                                    ind = i, Trait_categ = "Phenological",
                                    Clim = "Temperature")
  morphT <- fit_all_acountphylo(data_MA = Coefs_Aut_sp, vtree = "./data/phylogenies_100",
                                     ind = i, Trait_categ = "Morphological",
                                     Clim = "Temperature")
  phenP <- fit_all_acountphylo(data_MA = Coefs_Aut_precip, vtree = "./data/phylogenies_100",
                                    ind = i, Trait_categ = "Phenological",
                                    Clim = "Precipitation")
  morphP <- fit_all_acountphylo(data_MA = Coefs_Aut_precip, vtree = "./data/phylogenies_100",
                                     ind = i, Trait_categ = "Morphological",
                                     Clim = "Precipitation")
  ef_all <- rbind(ef_all, phenT, morphT, phenP, morphP)
end.time <- Sys.time()
print(end.time)
}


# save the file for the future use
saveRDS(ef_all, file = "./output_all/ef_all_difPhylo1.rds")

# read the file in
ef_1 <- readRDS("./output_all/ef_all_difPhylo.rds")
ef_2 <- readRDS("./output_all/ef_all_difPhylo1.rds")
ef_2 %<>%
  filter(phylo != 2)

ef_all <- rbind(ef_1, ef_2)
length(unique(ef_all$phylo)) # 108 replicates
# check in what proportion of those random phylogeny-based analyses the phylo-corrected model
# fits the data better than the model without accounting for phylogeny

ef_all_phylo <- ef_all %>%
  mutate(Phylobetter = ifelse(AIC_phylo < AIC_nophyl, "Yes", "No")) %>%
  group_by(Trait_Categ, Clim, Relation, phylo) %>%
  slice(1L) %>%
  ungroup()
108*2*2*6  # 2592 ok
nrow(ef_all_phylo)

nrow(ef_all_phylo[ef_all_phylo$Phylobetter == 'Yes', ]) /nrow(ef_all_phylo)  # 0
# in none of the 108 replicates the model with phylocorrection fits
# better than the model without accounting for phylogeny

