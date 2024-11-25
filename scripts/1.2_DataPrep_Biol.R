#### This script prepares biological data for sliding window analyses using climwin

#devtools::install_github('radchukv/sTraitChange')
library(sTraitChange)
library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyverse)

# read in the dataset
biol_dat <- read.csv('./data/Data_14_08_2024.csv', header = T)

biol_dat %<>% mutate(across(where(is_character), as_factor))

levels(biol_dat$Country)

length(unique(biol_dat$ID))  ## 217
## alaskan study with many NAs, exclude
biol_dat <- droplevels(subset(biol_dat, ! ID %in% c(228, 231, 232, 234)))


biol_dat$Unit_trait[biol_dat$Unit_trait == 'KG'] <- 'kg'
biol_dat <- droplevels(biol_dat)

length(unique(biol_dat$ID))     ## 213

## some Demog_rate_Categ have white spaces... Correcting
biol_dat$Demog_rate_Categ <- trimws(biol_dat$Demog_rate_Categ)

# some exploration for stats
biol_start <- biol_dat %>%
  group_by(ID) %>%
  slice_head()
min(biol_start$Year)
max(biol_start$Year)
hist(biol_start$Year)
median(biol_start$Year)  ## 1994

biol_end <- biol_dat %>%
  group_by(ID) %>%
  slice_tail()
min(biol_end$Year)
max(biol_end$Year)
hist(biol_end$Year)
median(biol_end$Year)  #2014

## keep only European studies and US - for now (two separate datasets, because climwin are diff.)
biol_eu <- droplevels(subset(biol_dat, ! Country %in% c('Antarctica', 'Australia',
                                                        'Canada', 'Falkland Islands',
                                                        'Greenland', 'Mexico',
                                                        'New Zealand', 'South Africa',
                                                        'South Atlantic Ocean',
                                                        'South Georgia', 'Svalbard',
                                                        'Taiwan', 'USA', 'Venezuela')))

length(unique(biol_eu$ID))        # 110
length(levels(biol_eu$Species))   ## 30
length(levels(biol_eu$Location))  ## 40
levels(biol_eu$Taxon)  # 3


# standardize units
levels(biol_eu$Unit_trait)

## use the funciton to converting the dates that are not Julian to Julian dates
biol_eu_stand <- convert_JulianDay(biol_data = biol_eu)

## US
biol_us <- droplevels(subset(biol_dat, Country %in% c('USA')))

length(unique(biol_us$ID))        ## 46
length(levels(biol_us$Species))   ## 10
length(levels(biol_us$Location))  ## 15
levels(biol_us$Taxon)             ## 4
table(biol_us$Taxon)   ## reptiles lead
table(biol_us$Trait_Categ)  ## mostly morphol

# standardize units
levels(biol_us$Unit_trait)  ## so convertion of dates to Julian not needed

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                               prepare data for analyses                                      ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 1. EUrope no sea birds
eu_noSea <- prep_subset(data = biol_eu_stand, Seabird = FALSE)
nrow(eu_noSea$Sel[[1]])    # 87
length(unique(eu_noSea$subdata[[1]]$ID))   #87

max(eu_noSea$subdata[[1]]$Trait_mean[eu_noSea$subdata[[1]]$Trait_Categ == 'Phenological'], na.rm = T)
hist(eu_noSea$subdata[[1]]$Trait_mean[eu_noSea$subdata[[1]]$Trait_Categ == 'Phenological'])


# 2. US no sea birds
## exclude sea birds
us_noSea <- prep_subset(data = biol_us, Seabird = FALSE)
nrow(us_noSea$Sel[[1]])    ## 44
length(unique(us_noSea$subdata[[1]]$ID))   ## 44

max(us_noSea$subdata[[1]]$Trait_mean[us_noSea$subdata[[1]]$Trait_Categ == 'Phenological'], na.rm = T)
hist(us_noSea$subdata[[1]]$Trait_mean[us_noSea$subdata[[1]]$Trait_Categ == 'Phenological'])


## 3. Sea birds
## keep only sea birds
all_Sea <- prep_subset(data = biol_dat, Seabird = TRUE)
nrow(all_Sea$Sel[[1]])  # 46
length(unique(all_Sea$subdata[[1]]$ID))      # 46

max(all_Sea$subdata[[1]]$Trait_mean[all_Sea$subdata[[1]]$Trait_Categ == 'Phenological'], na.rm = T)
hist(all_Sea$subdata[[1]]$Trait_mean[all_Sea$subdata[[1]]$Trait_Categ == 'Phenological'])

## 4. Other studies (no sea birds) for which we use world data - pay attention, this list has to be updated?
biol_w <- droplevels(subset(biol_dat, ID %in%
                              c(unique(subset(biol_dat, Country == 'Canada')$ID),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Wheelwright_et_al']),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Cheng_et_al']),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Nater_et_al']),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Veiberg_et_al']),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Jenouvrier_et_al']),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Tarwater&Beissinger']),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Veran&Beissinger']))))
length(unique(biol_w$ID))        ## 35
length(levels(biol_w$Species))   ## 22
length(levels(biol_w$Location))  ## 15
levels(biol_w$Taxon)             ## 3

# no need to standardize units
levels(biol_w$Unit_trait)

rest_w <- prep_subset(data = biol_w, Seabird = FALSE)
nrow(rest_w$Sel[[1]])                     ## 27
length(unique(rest_w$subdata[[1]]$ID))    ## 27

max(rest_w$subdata[[1]]$Trait_mean[rest_w$subdata[[1]]$Trait_Categ == 'Phenological'], na.rm = T)
hist(rest_w$subdata[[1]]$Trait_mean[rest_w$subdata[[1]]$Trait_Categ == 'Phenological'])


## 5. Australia
biol_au <- droplevels(subset(biol_dat, Country %in% c('Australia')))
biol_au_noSea <- prep_subset(biol_au, Seabird = FALSE)
nrow(biol_au_noSea$Sel[[1]])                      ## 9
length(unique(biol_au_noSea$subdata[[1]]$ID))     ## 9

max(biol_au_noSea$subdata[[1]]$Trait_mean[biol_au_noSea$subdata[[1]]$Trait_Categ == 'Phenological'], na.rm = T)
hist(biol_au_noSea$subdata[[1]]$Trait_mean[biol_au_noSea$subdata[[1]]$Trait_Categ == 'Phenological'])

