## This script runs sliding window analyses (using climwin package)
# across all studies in the dataset
# !!!!!!! Attention: it will not run out of the box
# as climate data are not shared in the repo - too heavy


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                            I. T analyses                                       ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source('./scripts/1.1_DataPrep_Clim.R')
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
## check for the difference betweeen the run and not run studies
files <- list.files('./output_climwin_temp/', pattern = 'Rand.RDS', full.names = TRUE)

testIDs <- unlist(lapply(files, FUN = function(x){
  as.numeric(strsplit(strsplit(x, 'temp/')[[1]][2], split = '_')[[1]][1])
}))  ## on Mac: temp//, on PC: temp/
length(testIDs)
unique(eu_noSea$Sel[[1]]$ID)[which (! unique(eu_noSea$Sel[[1]]$ID) %in%
                                      testIDs)]  ##422, 423 - Schaub
# ## + check what to re=run for new recruitment cat
# recr <- c(80, 83, 92, 151, 2, 319, 554, 330, 331, 332, 408, 409, 410, 411, 412, 413,
#           422, 423, 424, 425, 432, 445, 461, 464, 465, 312, 313, 314, 315, 539, 541,
#           543, 547, 549, 552, 188, 565, 581, 582)
# unique(eu_noSea$Sel[[1]]$ID[eu_noSea$Sel[[1]]$Demog_rate_Categ ==
#                               'Recruitment'])[
#                                 which (unique(eu_noSea$Sel[[1]]$ID[
#                                   eu_noSea$Sel[[1]]$Demog_rate_Categ ==
#                                     'Recruitment']) %in% recr)] ## for these to re-run climwin


## run additional from Schroeder and Avilles
IDs_to_run <- c(499, 514)

IDs_subs <- c(14:17)

eu_noSea_temp_add <- lapply(IDs_to_run, FUN = function(x){
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

# eu_noSea_temp_sub <- lapply(unique(eu_noSea$Sel[[1]]$ID[eu_noSea$Sel[[1]]$Demog_rate_Categ ==
#                                                           'Recruitment'])[
#                                                             which (unique(eu_noSea$Sel[[1]]$ID[
#                                                               eu_noSea$Sel[[1]]$Demog_rate_Categ ==
#                                                                 'Recruitment']) %in% recr)], FUN = function(x){
#                                                                   climwin_proc(biol_data = eu_noSea$subdata[[1]],
#                                                                                clim_data = meanT,
#                                                                                ID = x, randwin = TRUE,
#                                                                                seednum = 1302, repeats = 200,
#                                                                                cinterval = 'week',
#                                                                                plot_check = TRUE,
#                                                                                out_clim = 'output_climwin_temp',
#                                                                                stat = 'mean',
#                                                                                oneGrid = FALSE, explanYear = TRUE,
#                                                                                startWindow = 0, endWindow = 104, RefMon = NA)})

## 2. USA
unique(us_noSea$Sel[[1]]$ID)[which (! unique(us_noSea$Sel[[1]]$ID) %in%
                                      testIDs)] #567, 555 - Lany et al - run!
# unique(us_noSea$Sel[[1]]$ID[us_noSea$Sel[[1]]$Demog_rate_Categ ==
#                               'Recruitment'])[
#                                 which (unique(us_noSea$Sel[[1]]$ID[
#                                   us_noSea$Sel[[1]]$Demog_rate_Categ ==
#                                     'Recruitment']) %in% recr)]  ## these + 567, 555


## run additional from Blumstein & Martin
IDs_to_run <- c(18:20, 23:26)


US_noSea_temp_add <- lapply(IDs_to_run, FUN = function(x){
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

# US_noSea_temp_sub <- lapply(c(567, 555, unique(us_noSea$Sel[[1]]$ID[us_noSea$Sel[[1]]$Demog_rate_Categ ==
#                                                                       'Recruitment'])[
#                                                                         which (unique(us_noSea$Sel[[1]]$ID[
#                                                                           us_noSea$Sel[[1]]$Demog_rate_Categ ==
#                                                                             'Recruitment']) %in% recr)]), FUN = function(x){
#                                                                               climwin_proc(biol_data = us_noSea$subdata[[1]],
#                                                                                            clim_data = meanT_world,
#                                                                                            ID = x, randwin = TRUE,
#                                                                                            seednum = 1302, repeats = 200,
#                                                                                            cinterval = 'week',
#                                                                                            plot_check = TRUE,
#                                                                                            out_clim = 'output_climwin_temp',
#                                                                                            stat = 'mean',
#                                                                                            oneGrid = FALSE, explanYear = TRUE,
#                                                                                            startWindow = 0, endWindow = 104, RefMon = NA)})

## 3. Seabirds
unique(all_Sea$Sel[[1]]$ID)[which (! unique(all_Sea$Sel[[1]]$ID) %in%
                                     testIDs)]  ## 2 - new Reed, 124, 125 (new Teplitsky) - run
# unique(all_Sea$Sel[[1]]$ID[all_Sea$Sel[[1]]$Demog_rate_Categ ==
#                              'Recruitment'])[
#                                which (unique(all_Sea$Sel[[1]]$ID[
#                                  all_Sea$Sel[[1]]$Demog_rate_Categ ==
#                                    'Recruitment']) %in% recr)]  ## 2 - only 2 to run
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

# Sea_temp_sub <- lapply(c(2), FUN = function(x){
#   climwin_proc(biol_data = all_Sea$subdata[[1]],
#                clim_data = rot_SST,
#                ID = x, randwin = TRUE,
#                seednum = 1302, repeats = 200,
#                cinterval = 'week',
#                plot_check = TRUE,
#                out_clim = 'output_climwin_temp',
#                stat = 'mean',
#                oneGrid = FALSE, explanYear = TRUE,
#                startWindow = 0, endWindow = 104, RefMon = NA)})

## check for 430 climwin_proc (whether NYears carry over to the next step in the output dataset)
check_430 <- climwin_proc(biol_data = all_Sea$subdata[[1]],
                          clim_data = rot_SST,
                          ID = 430, randwin = TRUE,
                          seednum = 1302, repeats = 200,
                          cinterval = 'week',
                          plot_check = TRUE,
                          out_clim = 'output_climwin_temp',
                          stat = 'mean',
                          oneGrid = FALSE, explanYear = TRUE,
                          startWindow = 0, endWindow = 104, RefMon = NA)

## 4. other studies for which we use world data
unique(rest_w$Sel[[1]]$ID)[which (! unique(rest_w$Sel[[1]]$ID) %in%
                                    testIDs)]  ##581, 582 - Aubry, run
# unique(rest_w$Sel[[1]]$ID[rest_w$Sel[[1]]$Demog_rate_Categ ==
#                             'Recruitment'])[
#                               which (unique(rest_w$Sel[[1]]$ID[
#                                 rest_w$Sel[[1]]$Demog_rate_Categ ==
#                                   'Recruitment']) %in% recr)]  ## run all of these


## run additional for 565
IDs_to_run_w <- c(565)
others_w_add <- lapply(IDs_to_run_w, FUN = function(x){
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

# others_w_sub <- lapply(unique(rest_w$Sel[[1]]$ID[rest_w$Sel[[1]]$Demog_rate_Categ ==
#                                                    'Recruitment'])[
#                                                      which (unique(rest_w$Sel[[1]]$ID[
#                                                        rest_w$Sel[[1]]$Demog_rate_Categ ==
#                                                          'Recruitment']) %in% recr)], FUN = function(x){
#                                                            climwin_proc(biol_data = rest_w$subdata[[1]],
#                                                                         clim_data = meanT_world,
#                                                                         ID = x, randwin = TRUE,
#                                                                         seednum = 1302, repeats = 200,
#                                                                         cinterval = 'week',
#                                                                         plot_check = TRUE,
#                                                                         out_clim = 'output_climwin_temp',
#                                                                         stat = 'mean',
#                                                                         oneGrid = FALSE, explanYear = TRUE,
#                                                                         startWindow = 0, endWindow = 104, RefMon = NA)})

## 5. Australia
unique(biol_au_noSea$Sel[[1]]$ID)[which (! unique(biol_au_noSea$Sel[[1]]$ID) %in%
                                           testIDs)] ##OK
# unique(biol_au_noSea$Sel[[1]]$ID[biol_au_noSea$Sel[[1]]$Demog_rate_Categ ==
#                                    'Recruitment'])[
#                                      which (unique(biol_au_noSea$Sel[[1]]$ID[
#                                        biol_au_noSea$Sel[[1]]$Demog_rate_Categ ==
#                                          'Recruitment']) %in% recr)]  ## 552 only


## run additional for 549
IDs_to_run_au <- c(549)
au_temp_add <- lapply(IDs_to_run_au, FUN = function(x){
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

# au_temp_sub <- lapply(c(552), FUN = function(x){
#   climwin_proc(biol_data = biol_au_noSea$subdata[[1]],
#                clim_data = meanT_AU,
#                ID = x, randwin = TRUE,
#                seednum = 1302, repeats = 200,
#                cinterval = 'week',
#                plot_check = TRUE,
#                out_clim = 'output_climwin_temp',
#                stat = 'mean',
#                oneGrid = FALSE, explanYear = TRUE,
#                startWindow = 0, endWindow = 104, RefMon = NA)})


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


# eu_noSea_precip_sub <- lapply(unique(eu_noSea$Sel[[1]]$ID[eu_noSea$Sel[[1]]$Demog_rate_Categ ==
#                                                             'Recruitment'])[
#                                                               which (unique(eu_noSea$Sel[[1]]$ID[
#                                                                 eu_noSea$Sel[[1]]$Demog_rate_Categ ==
#                                                                   'Recruitment']) %in% recr)], FUN = function(x){
#                                                                     climwin_proc(biol_data = eu_noSea$subdata[[1]],
#                                                                                  clim_data = meanP,
#                                                                                  ID = x, randwin = TRUE,
#                                                                                  seednum = 1302, repeats = 200,
#                                                                                  cinterval = 'week',
#                                                                                  plot_check = TRUE,
#                                                                                  out_clim = 'output_climwin_precip',
#                                                                                  stat = 'sum',
#                                                                                  oneGrid = FALSE, explanYear = TRUE,
#                                                                                  startWindow = 0, endWindow = 104, RefMon = NA)})






## 2. US

## run additional from Blumstein& Martin

us_noSea_prec_add <- lapply(IDs_to_run, FUN = function(x){
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

# US_noSea_precip_sub <- lapply(c(567, 555, unique(us_noSea$Sel[[1]]$ID[us_noSea$Sel[[1]]$Demog_rate_Categ ==
#                                                                         'Recruitment'])[
#                                                                           which (unique(us_noSea$Sel[[1]]$ID[
#                                                                             us_noSea$Sel[[1]]$Demog_rate_Categ ==
#                                                                               'Recruitment']) %in% recr)]), FUN = function(x){
#                                                                                 climwin_proc(biol_data = us_noSea$subdata[[1]],
#                                                                                              clim_data = totP_US,
#                                                                                              ID = x, randwin = TRUE,
#                                                                                              seednum = 1302, repeats = 200,
#                                                                                              cinterval = 'week',
#                                                                                              plot_check = TRUE,
#                                                                                              out_clim = 'output_climwin_precip',
#                                                                                              stat = 'sum',
#                                                                                              oneGrid = FALSE, explanYear = TRUE,
#                                                                                              startWindow = 0, endWindow = 104, RefMon = NA)})

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

# Sea_precip_sub <- lapply(c(2), FUN = function(x){
#   climwin_proc(biol_data = all_Sea$subdata[[1]],
#                clim_data = rot_prec,
#                ID = x, randwin = TRUE,
#                seednum = 1302, repeats = 200,
#                cinterval = 'week',
#                plot_check = TRUE,
#                out_clim = 'output_climwin_precip',
#                stat = 'sum',
#                oneGrid = FALSE, explanYear = TRUE,
#                startWindow = 0, endWindow = 104, RefMon = NA)})

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

# others_w_precip_sub <- lapply(unique(rest_w$Sel[[1]]$ID[rest_w$Sel[[1]]$Demog_rate_Categ ==
#                                                           'Recruitment'])[
#                                                             which (unique(rest_w$Sel[[1]]$ID[
#                                                               rest_w$Sel[[1]]$Demog_rate_Categ ==
#                                                                 'Recruitment']) %in% recr)], FUN = function(x){
#                                                                   climwin_proc(biol_data = rest_w$subdata[[1]],
#                                                                                clim_data = rot_prec,
#                                                                                ID = x, randwin = TRUE,
#                                                                                seednum = 1302, repeats = 200,
#                                                                                cinterval = 'week',
#                                                                                plot_check = TRUE,
#                                                                                out_clim = 'output_climwin_precip',
#                                                                                stat = 'sum',
#                                                                                oneGrid = FALSE, explanYear = TRUE,
#                                                                                startWindow = 0, endWindow = 104, RefMon = NA)})


## run additional for 565
others_w_prec_add <- lapply(IDs_to_run_w, FUN = function(x){
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

## run additional for 549
au_prec_add <- lapply(IDs_to_run_au, FUN = function(x){
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

# au_precip_sub <- lapply(c(552), FUN = function(x){
#   climwin_proc(biol_data = biol_au_noSea$subdata[[1]],
#                clim_data = meanT_AU,
#                ID = x, randwin = TRUE,
#                seednum = 1302, repeats = 200,
#                cinterval = 'week',
#                plot_check = TRUE,
#                out_clim = 'output_climwin_precip',
#                stat = 'sum',
#                oneGrid = FALSE, explanYear = TRUE,
#                startWindow = 0, endWindow = 104, RefMon = NA)})

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

seaT_anal <- analyse_climwin(ID = 430, biol_data = all_Sea$subdata[[1]],
                          out_clim = 'output_climwin_temp',
                          randwin = TRUE, metric = 'AIC',
                          MinDur = 1, MaxDur = 40,
                          deltaThresh = -7, test_winDur = FALSE,
                          out_for_SEM = 'output_forSEM_temp',
                          oneGrid = FALSE, explanYear = TRUE,
                          endWindow = 104)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                III.1. Climwin results Temperature                       ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## if not all studies were analysed yet, get the IDs directly from the
## files that are in the folder
## getting the IDs for analyses
# files <- list.files('./output_climwin_temp/', full.names = TRUE)
#
# testIDs <- unlist(lapply(files, FUN = function(x){
#   as.numeric(strsplit(strsplit(x, 'temp/')[[1]][2], split = '_')[[1]][1])
# })
#        )
#
# Test_anal <- lapply(testIDs,
#                         FUN = function(x){
#                           analyse_climwin(ID = x, biol_data = eu_noSea$subdata[[1]],
#                                           out_clim = 'output_climwin_temp',
#                                           randwin = TRUE, metric = 'AIC',
#                                           MinDur = 1, MaxDur = 40,
#                                           deltaThresh = -7, test_winDur = FALSE,
#                                           out_for_SEM = 'output_forSEM_temp',
#                                           oneGrid = FALSE, explanYear = TRUE,
#                                           endWindow = 104, RefMon = NA)
#                         })

### Analyse temperature datasets
## 1. EU no sea

## check first if all IDs in the dataset are actually present in the folder
files <- list.files('./output_climwin_temp/', pattern = 'Rand.RDS', full.names = TRUE)

testIDs <- unlist(lapply(files, FUN = function(x){
  as.numeric(strsplit(strsplit(x, 'temp//')[[1]][2], split = '_')[[1]][1])
}))
length(testIDs)  ## 283  ##300  ## 308
length(unique(eu_noSea$Sel[[1]]$ID))    ## 132  ##144  ## 133
length(unique(eu_noSea$Sel[[1]]$ID)[which (unique(eu_noSea$Sel[[1]]$ID) %in%
                               testIDs)]) ## 144   ## 133

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

length(eu_Anal)  ## 133

## 2. US no sea
length(unique(us_noSea$Sel[[1]]$ID))    ## 44  ## 47
length(unique(us_noSea$Sel[[1]]$ID)[which (unique(us_noSea$Sel[[1]]$ID) %in%
                                      testIDs)])  ##44  ##47
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
length(us_Anal)  ##  47

## 3. seabirds
length(unique(all_Sea$Sel[[1]]$ID)[which (unique(all_Sea$Sel[[1]]$ID) %in%
                                            testIDs)])  ## 58  ## 54
length(unique(all_Sea$Sel[[1]]$ID)) ## 58  ## 54
length(unique(all_Sea$subdata[[1]]$ID))  ## 68
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

length(sea_Anal)  ##   ## 54

## check - why data for some studies (apparently rest of the wrold?) contains NYears if it should not
check_187 <- analyse_climwin(ID = 187, biol_data = rest_w$subdata[[1]],
                out_clim = 'output_climwin_temp',
                randwin = TRUE, metric = 'AIC',
                MinDur = 1, MaxDur = 40,
                deltaThresh = -7, test_winDur = FALSE,
                out_for_SEM = 'output_forSEM_temp',
                oneGrid = FALSE, explanYear = TRUE,
                endWindow = 104, RefMon = NA)

## 4. other studies with world climatic data
length(unique(rest_w$Sel[[1]]$ID)[which (unique(rest_w$Sel[[1]]$ID) %in%
                                            testIDs)])  ## 33
length(unique(rest_w$Sel[[1]]$ID))  ## 33
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
                                            testIDs)])  ##11
length(unique(biol_au_noSea$Sel[[1]]$ID))  ## 11
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
length(Au_Anal)

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
length(testIDs_prec)  ## 283  ##304
length(unique(eu_noSea$Sel[[1]]$ID))    ## 132  ## 144  ## 133
length(unique(eu_noSea$Sel[[1]]$ID)[which (unique(eu_noSea$Sel[[1]]$ID) %in%
                                             testIDs_prec)])  ##144   ## 133
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
length(unique(us_noSea$Sel[[1]]$ID))    ## 44  ##47
length(unique(us_noSea$Sel[[1]]$ID)[which (unique(us_noSea$Sel[[1]]$ID) %in%
                                             testIDs_prec)]) ##44   ##47
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

length(us_Anal_p)  ##47


## 3. seabirds
length(unique(all_Sea$Sel[[1]]$ID))    ## 58   ## 54
length(unique(all_Sea$Sel[[1]]$ID)[which (unique(all_Sea$Sel[[1]]$ID) %in%
                                            testIDs_prec)])  ## 57 - the ID 430 by Oppel et al is missing - that is fine   ## 53
## because we cannot run it since the precipitation is not over ocean..
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
length(sea_Anal_p)  ## 57   ## 53

## 4. other studies with world climatic data
length(unique(rest_w$Sel[[1]]$ID))    ## 33
length(unique(rest_w$Sel[[1]]$ID)[which (unique(rest_w$Sel[[1]]$ID) %in%
                                           testIDs_prec)]) ## 33
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

length(world_Anal_p) ## 33
## 5. Australia no sea
length(unique(biol_au_noSea$Sel[[1]]$ID)[which (unique(biol_au_noSea$Sel[[1]]$ID) %in%
                                                  testIDs_prec)]) ##11
length(unique(biol_au_noSea$Sel[[1]]$ID)) ##11
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

length(Au_Anal_p) ##11
