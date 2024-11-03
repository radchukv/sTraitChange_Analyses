#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## this script merges climwin results with demographic data,
## to have full dataset for SEM analyses


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                 1. For temperature data                      ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source('./scripts/1.2_DataPrep_Biol.R')

## getting the IDs for analyses
files <- list.files('./output_climwin_temp/', pattern = 'Rand.RDS', full.names = TRUE)

testIDs <- unlist(lapply(files, FUN = function(x){
  as.numeric(strsplit(strsplit(x, 'temp//')[[1]][2], split = '_')[[1]][1])
}))
length(testIDs)    ## 213
length(unique(eu_noSea$Sel[[1]]$ID))    ##87


## prepare the data (merging unique climwin studies with the respective
## biological data) for EU
temp_eu_SEM <- prep_SEM_input(prep_subset_climwin = eu_noSea,
                              out_for_SEM = 'output_forSEM_temp',
                              oneGrid = FALSE,
                              explanYear = TRUE,
                              endWindow = 104,
                              RefMon = NA,
                              selIDs = unique(eu_noSea$Sel[[1]]$ID)[which (unique(eu_noSea$Sel[[1]]$ID) %in%
                                                                             testIDs)])
length(unique(temp_eu_SEM$ID))  # 87
length(unique(eu_noSea$subdata[[1]]$ID))  #87
## check which ones are not there
# unique(eu_noSea$subdata[[1]]$ID)[which (! unique(eu_noSea$subdata[[1]]$ID) %in% unique(temp_eu_SEM$ID))]

## add the continent
temp_eu_SEM$Continent <- 'Europe'


## prepare the data for US
temp_us_SEM <- prep_SEM_input(prep_subset_climwin = us_noSea,
                              out_for_SEM = 'output_forSEM_temp',
                              oneGrid = FALSE,
                              explanYear = TRUE,
                              endWindow = 104,
                              RefMon = NA,
                              selIDs = unique(us_noSea$Sel[[1]]$ID)[which (unique(us_noSea$Sel[[1]]$ID) %in%
                                                                             testIDs)])

length(unique(temp_us_SEM$ID))   # 44
length(unique(us_noSea$subdata[[1]]$ID))    # 44

temp_us_SEM$Continent <- 'North America'


## prepare the data for seabirds
temp_seaB_SEM <- prep_SEM_input(prep_subset_climwin = all_Sea,
                                out_for_SEM = 'output_forSEM_temp',
                                oneGrid = FALSE,
                                explanYear = TRUE,
                                endWindow = 104,
                                RefMon = NA,
                                selIDs = unique(all_Sea$Sel[[1]]$ID)[which (unique(all_Sea$Sel[[1]]$ID) %in%
                                                                              testIDs)])
length(unique(temp_seaB_SEM$ID))    # 46
length(unique(all_Sea$subdata[[1]]$ID))   # 46
# check if all are there
# unique(all_Sea$subdata[[1]]$ID)[which (! unique(all_Sea$subdata[[1]]$ID) %in% unique(temp_seaB_SEM$ID))]


## assigning the continent based on the proximity to the mainland
unique(temp_seaB_SEM$Country)

IDs_perIslandCountry <- temp_seaB_SEM %>%
  dplyr::group_by(., Country, ID) %>%
  dplyr::group_by(., Country) %>%
  dplyr::summarise(., numID = dplyr::n_distinct(ID))
temp_seaB_SEM$Continent[temp_seaB_SEM$Country %in%  c('New Zealand', 'Australia')] <- 'Australia'
temp_seaB_SEM$Continent[temp_seaB_SEM$Country %in%  c('Sweden', 'Norway', 'Germany',
                                                      'Spain', 'Denmark', 'Scotland')] <- 'Europe'
temp_seaB_SEM$Continent[temp_seaB_SEM$Country %in% c('Falkland Islands', 'South Georgia',
                                                     'South Atlantic Ocean')] <- 'South America'
temp_seaB_SEM$Continent[temp_seaB_SEM$Country %in% c('Canada', 'Greenland', 'Mexico', 'USA')] <- 'North America'
temp_seaB_SEM$Continent[temp_seaB_SEM$Country == 'Antarctica'] <- 'Antarctica'


## prepare the data for other studies for which we used world data..
temp_restW_SEM <- prep_SEM_input(prep_subset_climwin = rest_w,
                                 out_for_SEM = 'output_forSEM_temp',
                                 oneGrid = FALSE,
                                 explanYear = TRUE,
                                 endWindow = 104,
                                 RefMon = NA,
                                 selIDs = unique(rest_w$Sel[[1]]$ID)[which (unique(rest_w$Sel[[1]]$ID) %in%
                                                                              testIDs)])
length(unique(temp_restW_SEM$ID))  # 27
length(unique(rest_w$subdata[[1]]$ID))      # 27
## unique(rest_w$subdata[[1]]$ID)[! unique(rest_w$subdata[[1]]$ID) %in%  unique(temp_restW_SEM$ID)]


## add continents, depending on the countries
unique(temp_restW_SEM$Country)
temp_restW_SEM$Continent[temp_restW_SEM$Country == 'Canada'] <- 'North America'
temp_restW_SEM$Continent[temp_restW_SEM$Country == 'Venezuela'] <- 'South America'
temp_restW_SEM$Continent[temp_restW_SEM$Country == 'Taiwan'] <- 'Asia'
temp_restW_SEM$Continent[temp_restW_SEM$Country == 'Svalbard'] <- 'Europe'
temp_restW_SEM$Continent[temp_restW_SEM$Country == 'South Africa'] <- 'Africa'


## prepare the data for Australia
temp_au_SEM <- prep_SEM_input(prep_subset_climwin = biol_au_noSea,
                              out_for_SEM = 'output_forSEM_temp',
                              oneGrid = FALSE,
                              explanYear = TRUE,
                              endWindow = 104,
                              RefMon = NA,
                              selIDs = unique(biol_au_noSea$Sel[[1]]$ID)[which (unique(biol_au_noSea$Sel[[1]]$ID) %in%
                                                                                  testIDs)])

length(unique(temp_au_SEM$ID))  #9
length(unique(biol_au_noSea$subdata[[1]]$ID))  # 9

temp_au_SEM$Continent <- 'Australia'


## bind all the SEM input datasets
temp_SEM <- rbind(temp_eu_SEM, temp_us_SEM, temp_seaB_SEM,
                  temp_restW_SEM, temp_au_SEM)
length(unique(temp_SEM$ID))    # 213
temp_SEM$Demog_rate <- trimws(temp_SEM$Demog_rate, 'both')


## check for the diff
length(unique(biol_dat$ID)[which(! unique(biol_dat$ID) %in% unique(temp_SEM$ID))])
unique(biol_dat$ID)[which(! unique(biol_dat$ID) %in% unique(temp_SEM$ID))]

## for now also fixing the Demog_rateCateg here (later will be fixed directly in the
## biol_dat) so that there are no white spaces
temp_SEM$Demog_rate_Categ <- trimws(temp_SEM$Demog_rate_Categ)
temp_SEM$Demog_rate <- trimws(temp_SEM$Demog_rate)
temp_SEM$Trait <- trimws(temp_SEM$Trait)

## add the variable differentiating the age class for the morphological trait (diff. effects for morphology)
# test_grepl <- grepl(pattern = 'Fledging|Subadult|Juv|Young', x= temp_SEM$Trait)
# sum(test_grepl)
# test_grep <- grep(pattern = 'Fledging|Subadult|Juv|Young', x= temp_SEM$Trait)
# length(test_grep)

temp_SEM$Trait_ageClass[temp_SEM$Trait_Categ == 'Morphological'] <- 'Adult'
temp_SEM$Trait_ageClass[temp_SEM$Trait_Categ == 'Morphological' &
                          grepl(pattern = 'Birth|Chick|Neonatal', temp_SEM$Trait)] <- 'Neonatal'
temp_SEM$Trait_ageClass[temp_SEM$Trait_Categ == 'Morphological' &
                          grepl(pattern = 'Fledging|Subadult|Juv|Young', temp_SEM$Trait)] <- 'Young'

## save it, for faster access later on
saveRDS(object = temp_SEM, file = './output_forSEM_temp/all_SEM.RDS')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`#
####                 2. For precipitation data                      ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## getting the IDs for analyses
files_precip <- list.files('./output_climwin_precip/', pattern = 'Rand.RDS', full.names = TRUE)

testIDs_precip <- unlist(lapply(files_precip, FUN = function(x){
  as.numeric(strsplit(strsplit(x, 'precip//')[[1]][2], split = '_')[[1]][1])
}))
length(testIDs_precip)     ## 307
length(unique(eu_noSea$Sel[[1]]$ID))     ##  87


## prepare the data (merging unique climwin studies with the respective
## biological data) for EU
precip_eu_SEM <- prep_SEM_input(prep_subset_climwin = eu_noSea,
                                out_for_SEM = 'output_forSEM_precip',
                                oneGrid = FALSE,
                                explanYear = TRUE,
                                endWindow = 104,
                                RefMon = NA,
                                selIDs = unique(eu_noSea$Sel[[1]]$ID)[which (unique(eu_noSea$Sel[[1]]$ID) %in%
                                                                               testIDs_precip)])

length(unique(precip_eu_SEM$ID))   ## 87
length(unique(eu_noSea$subdata[[1]]$ID))  ## 87

# unique(eu_noSea$subdata[[1]]$ID)[which (! unique(eu_noSea$subdata[[1]]$ID) %in% unique(precip_eu_SEM$ID))]

## add the continent
precip_eu_SEM$Continent <- 'Europe'

## prepare the data for US
precip_us_SEM <- prep_SEM_input(prep_subset_climwin = us_noSea,
                                out_for_SEM = 'output_forSEM_precip',
                                oneGrid = FALSE,
                                explanYear = TRUE,
                                endWindow = 104,
                                RefMon = NA,
                                selIDs = unique(us_noSea$Sel[[1]]$ID)[which (unique(us_noSea$Sel[[1]]$ID) %in%
                                                                               testIDs_precip)])

length(unique(precip_us_SEM$ID))  ##44
length(unique(us_noSea$subdata[[1]]$ID))    ##44

precip_us_SEM$Continent <- 'North America'
## check what is missing
# unique(us_noSea$subdata[[1]]$ID)[which (! unique(us_noSea$subdata[[1]]$ID) %in% unique(precip_us_SEM$ID))]


## prepare the data for seabirds
precip_seaB_SEM <- prep_SEM_input(prep_subset_climwin = all_Sea,
                                  out_for_SEM = 'output_forSEM_precip',
                                  oneGrid = FALSE,
                                  explanYear = TRUE,
                                  endWindow = 104,
                                  RefMon = NA,
                                  selIDs = unique(all_Sea$Sel[[1]]$ID)[which (unique(all_Sea$Sel[[1]]$ID) %in%
                                                                                testIDs_precip)])
length(unique(precip_seaB_SEM$ID))   ##45
length(unique(all_Sea$subdata[[1]]$ID))    ##46
## difference of 3 studies - figure out why
unique(all_Sea$subdata[[1]]$ID)[which (! unique(all_Sea$subdata[[1]]$ID) %in%
                                         unique(precip_seaB_SEM$ID))]
## 430
unique(all_Sea$subdata[[1]]$Study_Authors[all_Sea$subdata[[1]]$ID %in% c(430)])
## Oppel study - have to be excluded, no precip data from the sea


## assigning the continent based on the proximity to the mainland
unique(precip_seaB_SEM$Country)
## count number of ids per country

IDs_perIslandCountry <- precip_seaB_SEM %>%
  dplyr::group_by(., Country, ID) %>%
  dplyr::group_by(., Country) %>%
  dplyr::summarise(., numID = dplyr::n_distinct(ID))
precip_seaB_SEM$Continent[precip_seaB_SEM$Country %in%  c('New Zealand', 'Australia')] <- 'Australia'
precip_seaB_SEM$Continent[precip_seaB_SEM$Country %in%  c('Sweden', 'Norway', 'Germany',
                                                          'Spain', 'Denmark', 'Scotland')] <- 'Europe'
precip_seaB_SEM$Continent[precip_seaB_SEM$Country %in% c('Falkland Islands', 'South Georgia')] <-
  'South America'
precip_seaB_SEM$Continent[precip_seaB_SEM$Country %in% c('Canada', 'Greenland', 'Mexico', 'USA')] <- 'North America'
precip_seaB_SEM$Continent[precip_seaB_SEM$Country == 'Antarctica'] <- 'Antarctica'


## prepare the data for other studies for which we used world data..
precip_restW_SEM <- prep_SEM_input(prep_subset_climwin = rest_w,
                                   out_for_SEM = 'output_forSEM_precip',
                                   oneGrid = FALSE,
                                   explanYear = TRUE,
                                   endWindow = 104,
                                   RefMon = NA,
                                   selIDs = unique(rest_w$Sel[[1]]$ID)[which (unique(rest_w$Sel[[1]]$ID) %in%
                                                                                testIDs_precip)])
length(unique(precip_restW_SEM$ID))  ##27
length(unique(rest_w$subdata[[1]]$ID))    ##27
# check the difference
unique(rest_w$subdata[[1]]$ID)[which (! unique(rest_w$subdata[[1]]$ID) %in%
                                         unique(precip_restW_SEM$ID))]
## in fact 300 is former 321 (just overwrite here)
precip_restW_SEM$ID[precip_restW_SEM$ID == 321] <- 300


## add continents, depending on the countries
unique(precip_restW_SEM$Country)
precip_restW_SEM$Continent[precip_restW_SEM$Country == 'Canada'] <- 'North America'
precip_restW_SEM$Continent[precip_restW_SEM$Country == 'Venezuela'] <- 'South America'
precip_restW_SEM$Continent[precip_restW_SEM$Country == 'Taiwan'] <- 'Asia'
precip_restW_SEM$Continent[precip_restW_SEM$Country == 'Svalbard'] <- 'Europe'
precip_restW_SEM$Continent[precip_restW_SEM$Country == 'South Africa'] <- 'Africa'


## prepare the data for Australia
precip_au_SEM <- prep_SEM_input(prep_subset_climwin = biol_au_noSea,
                                out_for_SEM = 'output_forSEM_precip',
                                oneGrid = FALSE,
                                explanYear = TRUE,
                                endWindow = 104,
                                RefMon = NA,
                                selIDs = unique(biol_au_noSea$Sel[[1]]$ID)[which (unique(biol_au_noSea$Sel[[1]]$ID) %in%
                                                                                    testIDs_precip)])

length(unique(precip_au_SEM$ID))      ##9
length(unique(biol_au_noSea$subdata[[1]]$ID))  ##9

precip_au_SEM$Continent <- 'Australia'


## bind all the SEM input datasets
precip_SEM <- rbind(precip_eu_SEM, precip_us_SEM, precip_seaB_SEM,
                    precip_restW_SEM, precip_au_SEM)
length(unique(precip_SEM$ID))   ##212
length(unique(biol_dat$ID))  # 213 - okay, since Oppel study had to be ommited


## for now also fixing the Demog_rateCateg here (later will be fixed directly in the
## biol_dat) so that there are no white spaces
precip_SEM$Demog_rate_Categ <- trimws(precip_SEM$Demog_rate_Categ)
precip_SEM$Demog_rate <- trimws(precip_SEM$Demog_rate)
precip_SEM$Trait <- trimws(precip_SEM$Trait)

## add the variable differentiating the age class for the morphological trait (diff. effects for morphology)
precip_SEM$Trait_ageClass[precip_SEM$Trait_Categ == 'Morphological'] <- 'Adult'
precip_SEM$Trait_ageClass[precip_SEM$Trait_Categ == 'Morphological' &
                            grepl(pattern = 'Birth|Chick|Neonatal', precip_SEM$Trait)] <- 'Neonatal'
precip_SEM$Trait_ageClass[precip_SEM$Trait_Categ == 'Morphological' &
                            grepl(pattern = 'Fledging|Subadult|Juv|Young', precip_SEM$Trait)] <- 'Young'


## save it, for faster access later on
saveRDS(object = precip_SEM, file = './output_forSEM_precip/all_SEM.RDS')
