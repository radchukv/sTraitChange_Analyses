#### This script prepares climate data for climwin analyses
# ATTENTION: it will not run out of the box as no climate
# data are shared in the repo. In order to run this script
# the user has to download the climate data from the respective
# sources mentioned in Supplementary materials (Table S6) and
# place them in the respective directory on the path (./data)

#devtools::install_github('radchukv/sTraitChange')
#devtools::install_github("LiamDBailey/climwin", ref = "devel")
library(climwin)
library(sTraitChange)
library(dplyr)
library(ggplot2)
library(magrittr)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        Weather EUROPE                        ####
meanT <- raster::stack('./data/tg_ens_mean_0.1deg_reg_v18.0e.nc', varname = 'tg')
meanP <- raster::stack('./data/rr_ens_mean_0.1deg_reg_v18.0e.nc', varname = 'rr')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        Weather USA                           ####
# precip
precip_CPC_US <- list.files('./data/weather/precip_USA/', full.names = T)
totP_US <- raster::stack(lapply(1:length(precip_CPC_US), FUN = function(x){
  raster::stack(precip_CPC_US[x])}))
## need to bring all the coordinate systems to the common denominator...
## climatic model use longitude coordinates that run from 0 to 360, rather than -180 to 180
## see here: http://r-sig-geo.2731867.n2.nabble.com/Raster-projection-alignment-issues-td7587564.html
## needs rotation
totP_US <- raster::rotate(totP_US)

temp_CPC <- list.files('./data/weather/temp_world/', full.names = T)


## only read min_CPC, to assign the names of the layers for the mean stack
min_CPC <- grep('tmin', temp_CPC, value = TRUE)
minT <- raster::stack(lapply(1:length(min_CPC), FUN = function(x){
  raster::stack(min_CPC[x])}))


## read the prepared mean T stack
meanT_world <- raster::stack('./data/weather/temp_world/meanT_world_CPC.nc', varname = 'variable') ## does not preserve the names of the layers
names(meanT_world) <- names(minT)
# pdf('./output_climwin_temp/T_world.pdf')
# raster::plot(meanT_US[[1]])
# dev.off()


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        World: SST                            ####
SST <- list.files('./data/weather/SST', full.names = T)
daily_SST_single <- grep('day.mean', SST, value = TRUE)
all_SST <- raster::stack(lapply(1:length(daily_SST_single), FUN = function(x){
  raster::stack(daily_SST_single[x])}))


## read the prepared rotated SST stack
rot_SST <- raster::stack('./data/weather/SST/SST_rotated.nc') ## does not preserve the names of the layers
names(rot_SST) <- names(all_SST)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        World: Precipit                            ####
precip_W <- list.files('./data/weather/prec_world', full.names = T)
Precip_W_single <- grep('precip.', precip_W, value = T)
world_P <- raster::stack(lapply(1:length(Precip_W_single), FUN = function(x){
  raster::stack(Precip_W_single[x])}))

## read the prepared rotated world precip stack
rot_prec <- raster::stack('./data/weather/prec_world/Prec_rotated.nc') ## does not preserve the names of the layers
names(rot_prec) <- names(world_P)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        Weather Australia                           ####
# precip
precip_AU <- list.files('./data/weather/Australia/Precipit_AU/', full.names = T)
totP_AU <- raster::stack(lapply(1:length(precip_AU), FUN = function(x){
  raster::stack(precip_AU[x])}))
totP_AU

## temp
temp_AU <- list.files('./data/weather/Australia/Temp_AU/', full.names = T)

## minT - for names of the layers
min_T_AU <- grep('tmin', temp_AU, value = TRUE)
minT_AU <- raster::stack(lapply(1:length(min_T_AU), FUN = function(x){
  raster::stack(min_T_AU[x])}))

## read the prepared mean T
meanT_AU <- raster::stack('./data/weather/Australia/Temp_AU/meanT_AU.nc') ## does not preserve the names of the layers
names(meanT_AU) <- names(minT_AU)


