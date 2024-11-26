## This script runs sliding window analyses (using climwin package)
# across all studies in the dataset
# !!!!!!! Attention: it will not run out of the box
# as climate data are not shared in the repo (becasue it is too heavy
# for GitHub)  and this script depends on that data
library(sTraitChange)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                            I. T analyses                                       ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source('./scripts/1.1_DataPrep_Clim.R') ## !!! THIS WILL Not run unless the user downloads all the required
                                        ## climate data on their own (too big size to share via GitHub!!)
source('./scripts/1.2_DataPrep_Biol.R')
#~~~~~~~~~~~~~~~~~~~~~~~~
## test for SST analyses  - DO NOT RUN, just a test
unique(all_Sea$Sel[[1]]$ID)
test_seaT <- climwin_proc(biol_data = all_Sea$subdata[[1]],
             clim_data = rot_SST,
             ID = 430, randwin = TRUE,
             seednum = 1302, repeats = 200,
             cinterval = 'week',
             plot_check = TRUE, stat = 'mean',
             out_clim = 'output_climwin_temp',
             oneGrid = FALSE, explanYear = TRUE,
             startWindow = 0, endWindow = 104)

## test Europe no sea bird
test_climwin <- climwin_proc(biol_data = eu_noSea$subdata[[1]],
               clim_data = meanT,
               ID = 21, randwin = TRUE,
               seednum = 1302, repeats = 200,
               cinterval = 'week',
               plot_check = TRUE,
               out_clim = 'output_climwin_temp',
               stat = 'mean',
               oneGrid = FALSE, explanYear = TRUE,
               startWindow = 0, endWindow = 104, RefMon = NA)

## 1. Europe
## check for the difference between the run and not run studies
# files <- list.files('./output_climwin_temp/', pattern = 'Rand.RDS', full.names = TRUE)
#
# testIDs <- unlist(lapply(files, FUN = function(x){
#   as.numeric(strsplit(strsplit(x, 'temp/')[[1]][2], split = '_')[[1]][1])
# }))  ## on Mac: temp//, on PC: temp/
# length(testIDs)
# unique(eu_noSea$Sel[[1]]$ID)[which (! unique(eu_noSea$Sel[[1]]$ID) %in%
#                                       testIDs)]

eu_noSea_temp <- lapply(unique(eu_noSea$Sel[[1]]$ID), FUN = function(x){
  climwin_proc(biol_data = eu_noSea$subdata[[1]],
               clim_data = meanT,
               ID = x, randwin = TRUE,
               seednum = 1302, repeats = 200,
               cinterval = 'week',
               plot_check = TRUE,
               out_clim = 'output_climwin_temp',
               stat = 'mean',
               oneGrid = FALSE, explanYear = TRUE,
               startWindow = 0, endWindow = 104, RefMon = NA)})


## 2. USA
unique(us_noSea$Sel[[1]]$ID)[which (! unique(us_noSea$Sel[[1]]$ID) %in%
                                      testIDs)]

US_noSea_temp <- lapply(unique(us_noSea$Sel[[1]]$ID), FUN = function(x){
  climwin_proc(biol_data = us_noSea$subdata[[1]],
               clim_data = meanT_world,
               ID = x, randwin = TRUE,
               seednum = 1302, repeats = 200,
               cinterval = 'week',
               plot_check = TRUE,
               out_clim = 'output_climwin_temp',
               stat = 'mean',
               oneGrid = FALSE, explanYear = TRUE,
               startWindow = 0, endWindow = 104, RefMon = NA)})

## 3. Seabirds
unique(all_Sea$Sel[[1]]$ID)[which (! unique(all_Sea$Sel[[1]]$ID) %in%
                                     testIDs)]
Sea_temp <- lapply(unique(all_Sea$Sel[[1]]$ID), FUN = function(x){
  climwin_proc(biol_data = all_Sea$subdata[[1]],
               clim_data = rot_SST,
               ID = x, randwin = TRUE,
               seednum = 1302, repeats = 200,
               cinterval = 'week',
               plot_check = TRUE,
               out_clim = 'output_climwin_temp',
               stat = 'mean',
               oneGrid = FALSE, explanYear = TRUE,
               startWindow = 0, endWindow = 104, RefMon = NA)})

## 4. other studies for which we use world data
unique(rest_w$Sel[[1]]$ID)[which (! unique(rest_w$Sel[[1]]$ID) %in%
                                    testIDs)]

others_w <- lapply(unique(rest_w$Sel[[1]]$ID), FUN = function(x){
  climwin_proc(biol_data = rest_w$subdata[[1]],
               clim_data = meanT_world,
               ID = x, randwin = TRUE,
               seednum = 1302, repeats = 200,
               cinterval = 'week',
               plot_check = TRUE,
               out_clim = 'output_climwin_temp',
               stat = 'mean',
               oneGrid = FALSE, explanYear = TRUE,
               startWindow = 0, endWindow = 104, RefMon = NA)})


## 5. Australia
unique(biol_au_noSea$Sel[[1]]$ID)[which (! unique(biol_au_noSea$Sel[[1]]$ID) %in%
                                           testIDs)]

au_temp <- lapply(unique(biol_au_noSea$Sel[[1]]$ID), FUN = function(x){
  climwin_proc(biol_data = biol_au_noSea$subdata[[1]],
               clim_data = meanT_AU,
               ID = x, randwin = TRUE,
               seednum = 1302, repeats = 200,
               cinterval = 'week',
               plot_check = TRUE,
               out_clim = 'output_climwin_temp',
               stat = 'mean',
               oneGrid = FALSE, explanYear = TRUE,
               startWindow = 0, endWindow = 104, RefMon = NA)})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                            II. Precip analyses                                  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## 1. Europe
eu_noSea_precip <- lapply(unique(eu_noSea$Sel[[1]]$ID), FUN = function(x){
  climwin_proc(biol_data = eu_noSea$subdata[[1]],
               clim_data = meanP,
               ID = x, randwin = TRUE,
               seednum = 1302, repeats = 200,
               cinterval = 'week',
               plot_check = TRUE,
               out_clim = 'output_climwin_precip',
               stat = 'sum',
               oneGrid = FALSE, explanYear = TRUE,
               startWindow = 0, endWindow = 104, RefMon = NA)})


## 2. US
spring_precip_us <- lapply(unique(us_noSea$Sel[[1]]$ID), FUN = function(x){
  climwin_proc(biol_data = us_noSea$subdata[[1]],
               clim_data = totP_US,
               ID = x, randwin = TRUE,
               seednum = 1302, repeats = 200,
               cinterval = 'week',
               plot_check = TRUE,
               out_clim = 'output_climwin_precip',
               stat = 'sum',
               oneGrid = FALSE, explanYear = TRUE,
               startWindow = 0, endWindow = 104, RefMon = NA)})

## 3. Seabirds
Sea_precip <- lapply(unique(all_Sea$Sel[[1]]$ID), FUN = function(x){
  climwin_proc(biol_data = all_Sea$subdata[[1]],
               clim_data = rot_prec,
               ID = x, randwin = TRUE,
               seednum = 1302, repeats = 200,
               cinterval = 'week',
               plot_check = TRUE,
               out_clim = 'output_climwin_precip',
               stat = 'sum',
               oneGrid = FALSE, explanYear = TRUE,
               startWindow = 0, endWindow = 104,
               RefMon = NA, weatherVar = NA)})


## 4. other studies for which we use world data
others_w <- lapply(unique(rest_w$Sel[[1]]$ID), FUN = function(x){
  climwin_proc(biol_data = rest_w$subdata[[1]],
               clim_data = rot_prec,
               ID = x, randwin = TRUE,
               seednum = 1302, repeats = 200,
               cinterval = 'week',
               plot_check = TRUE,
               out_clim = 'output_climwin_precip',
               stat = 'sum',
               oneGrid = FALSE, explanYear = TRUE,
               startWindow = 0, endWindow = 104, RefMon = NA)})

## 5. Australia
au_precip <- lapply(unique(biol_au_noSea$Sel[[1]]$ID), FUN = function(x){
  climwin_proc(biol_data = biol_au_noSea$subdata[[1]],
               clim_data = totP_AU,
               ID = x, randwin = TRUE,
               seednum = 1302, repeats = 200,
               cinterval = 'week',
               plot_check = TRUE,
               out_clim = 'output_climwin_precip',
               stat = 'sum',
               oneGrid = FALSE, explanYear = TRUE,
               startWindow = 0, endWindow = 104, RefMon = NA)})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                III. Analyse climwin results                       ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## test analyse_climwin - single dataset - DO NOT RUN
t_anal <- analyse_climwin(ID = 300, biol_data = rest_w$subdata[[1]],
                          out_clim = 'output_climwin_temp',
                          randwin = TRUE, metric = 'AIC',
                          MinDur = 1, MaxDur = 40,
                          deltaThresh = -7, test_winDur = FALSE,
                          out_for_SEM = 'output_forSEM_temp',
                          oneGrid = FALSE, explanYear = TRUE,
                          endWindow = 104, RefMon = NA)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                III.1. Climwin results Temperature                       ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## if not all studies were analysed yet, get the IDs directly from the
## files that are in the folder

### Analyse temperature datasets
## 1. EU no sea

## check first if all IDs in the dataset are actually present in the folder
files <- list.files('./output_climwin_temp/', pattern = 'Rand.RDS', full.names = TRUE)

testIDs <- unlist(lapply(files, FUN = function(x){
  as.numeric(strsplit(strsplit(x, 'temp//')[[1]][2], split = '_')[[1]][1])
}))
length(testIDs)  ## 213
length(unique(eu_noSea$Sel[[1]]$ID))    ## 87
length(unique(eu_noSea$Sel[[1]]$ID)[which (unique(eu_noSea$Sel[[1]]$ID) %in%
                               testIDs)])

eu_Anal <- lapply(unique(eu_noSea$Sel[[1]]$ID)[which (unique(eu_noSea$Sel[[1]]$ID) %in%
                                                        testIDs)],
                    FUN = function(x){
                      analyse_climwin(ID = x, biol_data = eu_noSea$subdata[[1]],
                                      out_clim = 'output_climwin_temp',
                                      randwin = TRUE, metric = 'AIC',
                                      MinDur = 1, MaxDur = 40,
                                      deltaThresh = -7, test_winDur = FALSE,
                                      out_for_SEM = 'output_forSEM_temp',
                                      oneGrid = FALSE, explanYear = TRUE,
                                      endWindow = 104, RefMon = NA)
                    })

length(eu_Anal)  ## 87

## 2. US no sea
length(unique(us_noSea$Sel[[1]]$ID))    ## 44
length(unique(us_noSea$Sel[[1]]$ID)[which (unique(us_noSea$Sel[[1]]$ID) %in%
                                      testIDs)])  ## 44
us_Anal <- lapply(unique(us_noSea$Sel[[1]]$ID)[which (unique(us_noSea$Sel[[1]]$ID) %in%
                                                        testIDs)],
                  FUN = function(x){
                    analyse_climwin(ID = x, biol_data = us_noSea$subdata[[1]],
                                    out_clim = 'output_climwin_temp',
                                    randwin = TRUE, metric = 'AIC',
                                    MinDur = 1, MaxDur = 40,
                                    deltaThresh = -7, test_winDur = FALSE,
                                    out_for_SEM = 'output_forSEM_temp',
                                    oneGrid = FALSE, explanYear = TRUE,
                                    endWindow = 104, RefMon = NA)
                  })
length(us_Anal)  ##  44

## 3. seabirds
length(unique(all_Sea$Sel[[1]]$ID)[which (unique(all_Sea$Sel[[1]]$ID) %in%
                                            testIDs)])    ## 46
length(unique(all_Sea$Sel[[1]]$ID))     ## 46
length(unique(all_Sea$subdata[[1]]$ID))  ## 46
sea_Anal <- lapply(unique(all_Sea$Sel[[1]]$ID)[which (unique(all_Sea$Sel[[1]]$ID) %in%
                                                         testIDs)],
                  FUN = function(x){
                    analyse_climwin(ID = x, biol_data = all_Sea$subdata[[1]],
                                    out_clim = 'output_climwin_temp',
                                    randwin = TRUE, metric = 'AIC',
                                    MinDur = 1, MaxDur = 40,
                                    deltaThresh = -7, test_winDur = FALSE,
                                    out_for_SEM = 'output_forSEM_temp',
                                    oneGrid = FALSE, explanYear = TRUE,
                                    endWindow = 104, RefMon = NA)
                  })

length(sea_Anal)   ## 46


## 4. other studies with world climatic data
length(unique(rest_w$Sel[[1]]$ID)[which (unique(rest_w$Sel[[1]]$ID) %in%
                                            testIDs)])  ## 27
length(unique(rest_w$Sel[[1]]$ID))  ## 27
length(unique(rest_w$subdata[[1]]$ID))
world_Anal <- lapply(unique(rest_w$Sel[[1]]$ID)[which (unique(rest_w$Sel[[1]]$ID) %in%
                                                        testIDs)],
                   FUN = function(x){
                     analyse_climwin(ID = x, biol_data = rest_w$subdata[[1]],
                                     out_clim = 'output_climwin_temp',
                                     randwin = TRUE, metric = 'AIC',
                                     MinDur = 1, MaxDur = 40,
                                     deltaThresh = -7, test_winDur = FALSE,
                                     out_for_SEM = 'output_forSEM_temp',
                                     oneGrid = FALSE, explanYear = TRUE,
                                     endWindow = 104, RefMon = NA)
                   })

length(world_Anal)

## 5. Australia no sea
length(unique(biol_au_noSea$Sel[[1]]$ID)[which (unique(biol_au_noSea$Sel[[1]]$ID) %in%
                                            testIDs)])  ## 9
length(unique(biol_au_noSea$Sel[[1]]$ID))  ## 9
length(unique(biol_au_noSea$subdata[[1]]$ID))
Au_Anal <- lapply(unique(biol_au_noSea$Sel[[1]]$ID)[which (unique(biol_au_noSea$Sel[[1]]$ID) %in%
                                                        testIDs)],
                  FUN = function(x){
                    analyse_climwin(ID = x, biol_data = biol_au_noSea$subdata[[1]],
                                    out_clim = 'output_climwin_temp',
                                    randwin = TRUE, metric = 'AIC',
                                    MinDur = 1, MaxDur = 40,
                                    deltaThresh = -7, test_winDur = FALSE,
                                    out_for_SEM = 'output_forSEM_temp',
                                    oneGrid = FALSE, explanYear = TRUE,
                                    endWindow = 104, RefMon = NA)
                  })
length(Au_Anal)  ## 9

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                III.2. Climwin results Precipitation                       ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Analyse precipitation datasets
## 1. EU no sea

## check first if all IDs in the dataset are actually present in the folder
files_prec <- list.files('./output_climwin_precip/', pattern = 'Rand.RDS', full.names = TRUE)

testIDs_prec <- unlist(lapply(files_prec, FUN = function(x){
  as.numeric(strsplit(strsplit(x, 'precip//')[[1]][2], split = '_')[[1]][1])
}))
length(testIDs_prec)
length(unique(eu_noSea$Sel[[1]]$ID))    ## 87
length(unique(eu_noSea$Sel[[1]]$ID)[which (unique(eu_noSea$Sel[[1]]$ID) %in%
                                             testIDs_prec)])  ## 87
eu_Anal_p <- lapply(unique(eu_noSea$Sel[[1]]$ID)[which (unique(eu_noSea$Sel[[1]]$ID) %in%
                                                          testIDs_prec)],
                  FUN = function(x){
                    analyse_climwin(ID = x, biol_data = eu_noSea$subdata[[1]],
                                    out_clim = 'output_climwin_precip',
                                    randwin = TRUE, metric = 'AIC',
                                    MinDur = 1, MaxDur = 40,
                                    deltaThresh = -7, test_winDur = FALSE,
                                    out_for_SEM = 'output_forSEM_precip',
                                    oneGrid = FALSE, explanYear = TRUE,
                                    endWindow = 104, RefMon = NA)
                  })
length(eu_Anal_p)

## 2. US no sea
length(unique(us_noSea$Sel[[1]]$ID))    ## 44
length(unique(us_noSea$Sel[[1]]$ID)[which (unique(us_noSea$Sel[[1]]$ID) %in%
                                             testIDs_prec)]) ##44
us_Anal_p <- lapply(unique(us_noSea$Sel[[1]]$ID)[which (unique(us_noSea$Sel[[1]]$ID) %in%
                                                          testIDs_prec)],
                  FUN = function(x){
                    analyse_climwin(ID = x, biol_data = us_noSea$subdata[[1]],
                                    out_clim = 'output_climwin_precip',
                                    randwin = TRUE, metric = 'AIC',
                                    MinDur = 1, MaxDur = 40,
                                    deltaThresh = -7, test_winDur = FALSE,
                                    out_for_SEM = 'output_forSEM_precip',
                                    oneGrid = FALSE, explanYear = TRUE,
                                    endWindow = 104, RefMon = NA)
                  })

length(us_Anal_p)  ##44


## 3. seabirds
length(unique(all_Sea$Sel[[1]]$ID))    ## 46
length(unique(all_Sea$Sel[[1]]$ID)[which (unique(all_Sea$Sel[[1]]$ID) %in%
                                            testIDs_prec)])  ## 45 - the ID 430 by Oppel et al is missing - that is fine
## because we cannot run it since no precipitation data available  over ocean..
sea_Anal_p <- lapply(unique(all_Sea$Sel[[1]]$ID)[which (unique(all_Sea$Sel[[1]]$ID) %in%
                                                          testIDs_prec)],
                   FUN = function(x){
                     analyse_climwin(ID = x, biol_data = all_Sea$subdata[[1]],
                                     out_clim = 'output_climwin_precip',
                                     randwin = TRUE, metric = 'AIC',
                                     MinDur = 1, MaxDur = 40,
                                     deltaThresh = -7, test_winDur = FALSE,
                                     out_for_SEM = 'output_forSEM_precip',
                                     oneGrid = FALSE, explanYear = TRUE,
                                     endWindow = 104, RefMon = NA)
                   })
length(sea_Anal_p)  ## 45

## 4. other studies with world climatic data
length(unique(rest_w$Sel[[1]]$ID))    ## 27
length(unique(rest_w$Sel[[1]]$ID)[which (unique(rest_w$Sel[[1]]$ID) %in%
                                           testIDs_prec)]) ## 27
world_Anal_p <- lapply(unique(rest_w$Sel[[1]]$ID)[which (unique(rest_w$Sel[[1]]$ID) %in%
                                                           testIDs_prec)],
                     FUN = function(x){
                       analyse_climwin(ID = x, biol_data = rest_w$subdata[[1]],
                                       out_clim = 'output_climwin_precip',
                                       randwin = TRUE, metric = 'AIC',
                                       MinDur = 1, MaxDur = 40,
                                       deltaThresh = -7, test_winDur = FALSE,
                                       out_for_SEM = 'output_forSEM_precip',
                                       oneGrid = FALSE, explanYear = TRUE,
                                       endWindow = 104, RefMon = NA)
                     })

length(world_Anal_p) ## 27

## 5. Australia no sea
length(unique(biol_au_noSea$Sel[[1]]$ID)[which (unique(biol_au_noSea$Sel[[1]]$ID) %in%
                                                  testIDs_prec)]) ## 9
length(unique(biol_au_noSea$Sel[[1]]$ID)) ## 9
length(unique(biol_au_noSea$subdata[[1]]$ID))
Au_Anal_p <- lapply(unique(biol_au_noSea$Sel[[1]]$ID)[which (unique(biol_au_noSea$Sel[[1]]$ID) %in%
                                                               testIDs_prec)],
                  FUN = function(x){
                    analyse_climwin(ID = x, biol_data = biol_au_noSea$subdata[[1]],
                                    out_clim = 'output_climwin_precip',
                                    randwin = TRUE, metric = 'AIC',
                                    MinDur = 1, MaxDur = 40,
                                    deltaThresh = -7, test_winDur = FALSE,
                                    out_for_SEM = 'output_forSEM_precip',
                                    oneGrid = FALSE, explanYear = TRUE,
                                    endWindow = 104, RefMon = NA)
                  })

length(Au_Anal_p) ## 9
