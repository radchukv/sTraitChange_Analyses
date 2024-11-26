#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## this script fits SEMs to the data
library(sTraitChange)
library(ggplot2)
library(tidyverse)
library(magrittr)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                 1. For temperature data                      ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

temp_SEM <- readRDS(file = './output_forSEM_temp/all_SEM.RDS')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####              1.2. SEMs DD, weights, standard, autocor                         ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## test the function ID = 1
ID1_st_Aut <- fit_SEM(biol_data = temp_SEM, ID = 7,
                                 out_SEM = NULL,
                                 DD = 'n_effectGR',
                                 weight = TRUE, correlation = TRUE,
                                 standardize = TRUE, Trait = FALSE,
                                 simpleSEM = TRUE)


## fit SEMs for all studies

#  3:
# Error in nlme::gls(GR ~ det_Clim + Pop_mean +
# Trait_mean, correlation = nlme::corAR1(form = ~Year |  :
#    false convergence (8)                                                                                    false convergence (8)

## 221:: Error in glsEstimate(object, control = control) :
## computed "gls" fit is singular, rank 4
## Called from: glsEstimate(object, control = control)
## this was originally in the set - not sure this dataset at all is still in the dataset, and for starters would anyways try to fit everything that's possible

## 550:
# `coef<-.corARMA`(`*tmp*`, value = value[parMap[, i]]) :
#   Coefficient matrix not invertible
fitted_SEM_Autoc <- lapply(unique(temp_SEM$ID)[- which(unique(temp_SEM$ID) %in%
                                                         c(3, 221, 550))],
                             FUN = function(x){
                               fit_SEM(biol_data = temp_SEM, ID = x,
                                       out_SEM = NULL, DD = 'n_effectGR', # out_SEM  = 'output_fSEM_temp'
                                       weight = TRUE, correlation = TRUE,
                                       standardize = TRUE, Trait = FALSE,
                                       simpleSEM = TRUE)
                             })
length(fitted_SEM_Autoc)  ## 210


# getting the results of SEM
Cstat_Aut <- extract_res_SEM(list_fitSEM = fitted_SEM_Autoc, stat_extr = 'Cstat')
R2_Aut <- extract_res_SEM(list_fitSEM = fitted_SEM_Autoc, stat_extr = 'R2')
Coefs_Aut <- extract_res_SEM(list_fitSEM = fitted_SEM_Autoc, stat_extr = 'coefs')



## read-in the previously saved Cstat_Aut
#Cstat_Aut <- readRDS(file = './output_fSEM_temp/Cstat_allMods_Temp_Weights_DD_Autocorr.RDS')


## checking how Cstat p values look like
pdf('./output_SEM_all/Cstat_allMods_Temp_Weights_DD_Autocor.pdf')
hist(Cstat_Aut$P.Value, col = 'grey', breaks = 20, xlab = 'p value for Cstat', main = '')
dev.off()


ggplot(Cstat_Aut, aes(P.Value)) + geom_histogram() +
  facet_grid(rows = vars(Continent), scales = 'free') +
  theme_bw() + theme(axis.title = element_text(size = rel(2)),
                     axis.text = element_text(size = rel(1.5)),
                     strip.text = element_text(size = rel(1.7)))

sum(Cstat_Aut$P.Value > 0.05) ## 162 studies pass the Cstat GOF test
sum(Cstat_Aut$P.Value > 0.05) / nrow(Cstat_Aut)  ## 77%
## per trait category
sum(Cstat_Aut$P.Value[Cstat_Aut$Trait_Categ == 'Phenological'] > 0.05) / nrow(Cstat_Aut[Cstat_Aut$Trait_Categ == 'Phenological', ])  ##71%
sum(Cstat_Aut$P.Value[Cstat_Aut$Trait_Categ == 'Morphological'] > 0.05) / nrow(Cstat_Aut[Cstat_Aut$Trait_Categ == 'Morphological', ])  ##82%

## a figure for P value of Cstat split by trait category, for Supplement
summary_Cstat <- Cstat_Aut %>%
  dplyr::group_by(., Trait_Categ) %>%
  dplyr::summarise(., Median = median(P.Value)) %>%
  dplyr::mutate(., lab = c('a', 'b'))


##  description of the all data we got (not excluding the studies with bad GOF)
table(Cstat_Aut$Trait_Categ, Cstat_Aut$Continent)
table(Cstat_Aut$Trait_Categ, Cstat_Aut$Taxon)

table(Cstat_Aut$Taxon)  ## unbalanced
table(Cstat_Aut$Taxon, Cstat_Aut$Continent)

table(Cstat_Aut$Country)
sp <- as.data.frame(table(Cstat_Aut$Species))


colnames(sp) <- c('Species', 'Freq')
pdf('./output_SEM_all/data_descript_species_Temp_Autocor.pdf', width = 8)
ggplot(sp, aes(Species, Freq)) +
  geom_bar(stat = 'identity') +
  theme_bw() + theme(legend.position = 'none',
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text = element_text(size = rel(1.1), colour = 'black'),
                     axis.text.x = element_text(angle = 90, vjust = 0.5),
                     axis.title = element_text(size  = rel(1.5))) + xlab('Species')
dev.off()

mean(sp$Freq)  ## 2.8 - mean frequency per sp
median(sp$Freq)  ## 1
sp[sp$Freq >= 5, ]

hist(sp$Freq, breaks = 20)


## write out the Cstat file
saveRDS(object = Cstat_Aut, file = './output_fSEM_temp/Cstat_allMods_Temp_Weights_DD_Autocorr.RDS')
# Cstat_Aut <- readRDS(file = './output_fSEM_temp/Cstat_allMods_Temp_Weights_DD_Autocorr.RDS')

## looking at R2

## reading in the file saved previously
#R2_Aut <- readRDS(file = './output_fSEM_temp/R2_allMods_Temp_Weights_DD_Autocorr.RDS')

## renaming the levels of the response, for a nicer figure labels
R2_Aut_proc <- R2_Aut %>%
  dplyr::mutate(., Response = dplyr::recode(Response,
                                 Trait_mean = 'Trait'))
R2_Aut_proc$Response <- factor(R2_Aut_proc$Response,
                               levels = c('Trait', 'GR'))
summ_coefs_Aut <- R2_Aut_proc %>%
  dplyr::group_by(., Trait_Categ, Response) %>%
  dplyr::summarise(Mean = mean(R.squared, na.rm = T),
            Median = median(R.squared, na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(label = letters[1:4],
                y = c(21, 21, 15, 15))

summ_coefs_Aut


plot_R2_Tr <- ggplot(R2_Aut_proc, aes(x = R.squared)) +
  facet_grid(Trait_Categ ~ Response, scale = 'free_y') +
  geom_histogram() + theme_bw() +
  geom_vline(data = summ_coefs_Aut, aes(xintercept = Median),
             col = 'red', lwd = 1, lty = 2) +
  theme(strip.text = element_text(size = rel(1.3)),
        strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_text(data = summ_coefs_Aut, aes(label  = round(Median, 2),
            x = round(Median, 2) + 0.08, y = 13),
            col = 'red', size = 5) +
  geom_text(data = summ_coefs_Aut, aes(x = 0, y = y, label = label),
            fontface = 'bold', size = 5) + xlab(expression('R'^2)) +
  ylab('Frequency')

pdf('./plots_ms/FigS7_R2_allMods_Temp_Weights_DD_Autocor_DiffTraitCateg.pdf',
    width = 10, height = 10)
print(plot_R2_Tr)
dev.off()


summ_coefs_byResp_Aut <- R2_Aut_proc %>%
  dplyr::group_by(., Response) %>%
  dplyr::summarise(Mean = mean(R.squared, na.rm = T),
            Median = median(R.squared, na.rm = T))

ggplot(R2_Aut_proc, aes(x = R.squared)) +
  facet_wrap(vars(Response)) +
  geom_histogram() + theme_bw() +
  geom_vline(data = summ_coefs_byResp_Aut, aes(xintercept = Median),
             col = 'red', lwd = 1.4) + theme(strip.text.y = element_text(size = rel(0.8)))

## write out the R2 file
saveRDS(object = R2_Aut, file = './output_fSEM_temp/R2_allMods_Temp_Weights_DD_Autocorr.RDS')


## looking at the path coefficients
# merging predictor and response in one
Coefs_Aut$Relation <- paste(Coefs_Aut$Response, Coefs_Aut$Predictor, sep = '<-')

## split by traits and dem. rates
pdf('./output_SEM_all/PathCoefs_Boxplot_byTraitDemRate_Temp_Weights_DD_Autocor.pdf',
    width = 9)
ggplot(Coefs_Aut, aes(y = Estimate, x = Demog_rate_Categ, col = Trait_Categ)) +
  facet_wrap(vars(Relation)) + geom_boxplot() + theme_bw() +
  geom_hline(yintercept = 0, col = 'black', lty = 2, lwd= 0.7) +
  theme(axis.title = element_text(size = rel(1.4)),
        axis.text = element_text(size = rel(1.1)),
        strip.text = element_text(size = rel(1.4)),
        legend.position = 'bottom', legend.text = element_text(size = rel(1.1))) +
  xlab('Demographic rate')
dev.off()


## write out the coefs file
saveRDS(object = Coefs_Aut, file = './output_fSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####              1.3. SEMs no DD, weights, standard, autocor                         ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


## test the function ID = 418
SEM_ID1_st_Aut_noDD <- fit_SEM(biol_data = temp_SEM, ID = 7,
                               out_SEM = NULL,  # output_fSEM_noDD_temp
                               DD = 'none', weight = TRUE,
                               correlation = TRUE,
                               standardize = TRUE, Trait = FALSE,
                               simpleSEM = TRUE)



# [- which(unique(temp_SEM$ID) %in%
## c(221, 412, 413, 573, 581, 582))]
fitted_SEM_noDD <- lapply(unique(temp_SEM$ID),
                           FUN = function(x){
                             fit_SEM(biol_data = temp_SEM, ID = x,
                                     out_SEM = NULL, # 'output_fSEM_noDD_temp'
                                     DD = 'none',
                                     weight = TRUE, correlation = TRUE,
                                     standardize = TRUE, Trait = FALSE,
                                     simpleSEM = TRUE)
                           })
length(fitted_SEM_noDD)  ##213


# getting the results of SEM
Cstat_Aut_noDD <- extract_res_SEM(list_fitSEM = fitted_SEM_noDD, stat_extr = 'Cstat')
R2_Aut_noDD <- extract_res_SEM(list_fitSEM = fitted_SEM_noDD, stat_extr = 'R2')
Coefs_Aut_noDD <- extract_res_SEM(list_fitSEM = fitted_SEM_noDD, stat_extr = 'coefs')


## read-in the previously saved Cstat_Aut
#Cstat_Aut_noDD <- readRDS(file = './output_fSEM_noDD_temp/Cstat_allMods_Temp_Weights_noDD_Autocorr.RDS')

## write out the Cstat file
saveRDS(object = Cstat_Aut_noDD, file = './output_fSEM_noDD_temp/Cstat_allMods_Temp_Weights_noDD_Autocorr.RDS')


## looking at R2

## reading in the file saved previously
#R2_Aut_noDD <- readRDS(file = './output_fSEM_noDD_temp/R2_allMods_Temp_Weights_noDD_Autocorr.RDS')

## renaming the levels of the response, for a nicer figure labels
R2_Aut_noDD_proc <- R2_Aut_noDD %>%
  dplyr::mutate(., Response = dplyr::recode(Response,
                                            Demog_rate_mean = 'Demographic rate',
                                            Trait_mean = 'Trait'))
R2_Aut_noDD_proc$Response <- factor(R2_Aut_noDD_proc$Response, levels = c('Trait',
                                                      'Demographic rate', 'GR'))


## checking for trait category only
summ_coefs_Trait_noDD <- R2_Aut_noDD_proc %>%
  dplyr::group_by(., Response, Trait_Categ) %>%
  dplyr::summarise(Mean = mean(R.squared, na.rm = T),
                   Median = median(R.squared, na.rm = T))
summ_coefs_Trait_noDD$label <- letters[1:4]

summ_coefs_Trait_noDD
plot_R2_Trait_noDD <- ggplot(R2_Aut_noDD_proc, aes(x = R.squared)) +
  facet_grid(rows = vars(Response), cols = vars(Trait_Categ)) +
  geom_histogram() + theme_bw() +
  geom_vline(data = summ_coefs_Trait_noDD, aes(xintercept = Median),
             col = 'red', lwd = 1, lty = 2) +
  theme(strip.text = element_text(size = rel(1.3)),
        strip.background = element_rect(fill = 'white')) +
  geom_text(data = summ_coefs_Trait_noDD, mapping = aes(label  = round(Median, 2),
                                                   x = round(Median, 2) + 0.08, y = 5),
            col = 'red', size = 4) +
  geom_text(data = summ_coefs_Trait_noDD, mapping = aes(x = 0.03, y = 29, label = label),
            fontface = 'bold') + xlab(expression('R'^2))

plot_R2_Trait_noDD

## overall statistics across the trait categories
summ_coefs_byResp_Aut_noDD <- R2_Aut_noDD_proc %>%
  dplyr::group_by(., Response) %>%
  dplyr::summarise(Mean = mean(R.squared, na.rm = T),
                   Median = median(R.squared, na.rm = T))

ggplot(R2_Aut_noDD_proc, aes(x = R.squared)) +
  facet_wrap(vars(Response)) +
  geom_histogram() + theme_bw() +
  geom_vline(data = summ_coefs_byResp_Aut_noDD, aes(xintercept = Median),
             col = 'red', lwd = 1.4) + theme(strip.text.y = element_text(size = rel(0.8)))

## write out the R2 file
saveRDS(object = R2_Aut_noDD, file = './output_fSEM_noDD_temp/R2_allMods_Temp_Weights_noDD_Autocorr.RDS')


## looking at the path coefficients
# merging predictor and response in one
Coefs_Aut_noDD$Relation <- paste(Coefs_Aut_noDD$Response, Coefs_Aut_noDD$Predictor, sep = '<-')
table(Coefs_Aut_noDD$Demog_rate_Categ)

## plot, split by traits and dem. rates
pdf('./output_SEM_all/PathCoefs_Boxplot_byTraitDemRate_Temp_Weights_noDD_Autocor.pdf',
    width = 9)
ggplot(Coefs_Aut_noDD, aes(y = Estimate, x = Demog_rate_Categ, col = Trait_Categ)) +
  facet_wrap(vars(Relation)) + geom_boxplot() + theme_bw() +
  geom_hline(yintercept = 0, col = 'black', lty = 2, lwd= 0.7) +
  theme(axis.title = element_text(size = rel(1.4)),
        axis.text = element_text(size = rel(1.1)),
        strip.text = element_text(size = rel(1.4)),
        legend.position = 'bottom', legend.text = element_text(size = rel(1.1))) +
  xlab('Demographic rate')
dev.off()


## write out the coefs file
saveRDS(object = Coefs_Aut_noDD, file = './output_fSEM_noDD_temp/PathCoefs_allMods_Temp_Weights_noDD_Autocor.RDS')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                 2. For precipitation data                      ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

precip_SEM <- readRDS(file = './output_forSEM_precip/all_SEM.RDS')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####              2.1. SEMs DD, weights, standard, autocor                         ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## fit SEMs for all studies
## [- which(unique(precip_SEM$ID) %in%  c(221))]

fitted_SEM_Autoc_precip <- lapply(unique(precip_SEM$ID)[- which(unique(precip_SEM$ID) %in%  c(221))],
                           FUN = function(x){
                             fit_SEM(biol_data = precip_SEM, ID = x,
                                     out_SEM = NULL, DD = 'n_effectGR',  # 'output_fSEM_precip'
                                     weight = TRUE, correlation = TRUE,
                                     standardize = TRUE, Trait = FALSE,
                                     simpleSEM = TRUE)
                           })
length(fitted_SEM_Autoc_precip)  ## 211


# getting the results of SEM
Cstat_Aut_precip <- extract_res_SEM(list_fitSEM = fitted_SEM_Autoc_precip, stat_extr = 'Cstat')
R2_Aut_precip <- extract_res_SEM(list_fitSEM = fitted_SEM_Autoc_precip, stat_extr = 'R2')
Coefs_Aut_precip <- extract_res_SEM(list_fitSEM = fitted_SEM_Autoc_precip, stat_extr = 'coefs')

## comparing the IDs for the studies with SEM and those for which climwin was run
#unique(biol_dat$ID)[which (! unique(biol_dat$ID) %in% unique(Cstat_Aut_precip$ID))]
## [1] 221 430 - the study by Oppel excluded, no precip data for the sea (430)
# and 221 did not fit - convergence issues

## looking at the Cstat
pdf('./output_SEM_all/Cstat_allMods_Precip_Weights_DD_Autocor.pdf')
hist(Cstat_Aut_precip$P.Value, col = 'grey', breaks = 20, xlab = 'p value for Cstat', main = '')
dev.off()


ggplot(Cstat_Aut_precip, aes(P.Value)) + geom_histogram() +
  facet_grid(rows = vars(Continent), scales = 'free') +
  theme_bw() + theme(axis.title = element_text(size = rel(2)),
                     axis.text = element_text(size = rel(1.5)),
                     strip.text = element_text(size = rel(1.7)))
sum(Cstat_Aut_precip$P.Value > 0.05) ## 165 studies pass the Cstat GOF test
sum(Cstat_Aut_precip$P.Value > 0.05) / nrow(Cstat_Aut_precip)  ## 78%


## write out the Cstat file
saveRDS(object = Cstat_Aut_precip, file = './output_fSEM_precip/Cstat_allMods_Precip_Weights_DD_Autocorr.RDS')


## looking at R2
## loadign the previously saved file
#R2_Aut_precip <- readRDS(file = './output_fSEM_precip/R2_allMods_Precip_Weights_DD_Autocorr.RDS')

R2_Aut_precip %<>%
  dplyr::mutate(., Response = dplyr::recode(Response,
                                            Trait_mean = 'Trait'))
R2_Aut_precip$Response <- factor(R2_Aut_precip$Response,
                               levels = c('Trait', 'GR'))
summ_coefs_Aut_precip <- R2_Aut_precip %>%
  dplyr::group_by(., Trait_Categ, Response) %>%
  dplyr::summarise(Mean = mean(R.squared, na.rm = T),
                   Median = median(R.squared, na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(label = letters[1:4],
                y = rep(17, 4))

summ_coefs_Aut_precip


plot_R2_Tr_precip <- ggplot(R2_Aut_precip, aes(x = R.squared)) +
  facet_grid(Trait_Categ ~ Response, scale = 'free_y') +
  geom_histogram() + theme_bw() +
  geom_vline(data = summ_coefs_Aut_precip, aes(xintercept = Median),
             col = 'red', lwd = 1, lty = 2) +
  theme(strip.text = element_text(size = rel(1.3)),
        strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_text(data = summ_coefs_Aut_precip, aes(label  = round(Median, 2),
                                       x = round(Median, 2) + 0.08, y = 13),
            col = 'red', size = 5) +
  geom_text(data = summ_coefs_Aut_precip, aes(x = 0, y = y, label = label),
            fontface = 'bold', size = 5) + xlab(expression('R'^2)) +
  ylab('Frequency')

pdf('./plots_ms/FigS10_R2_allMods_Precip_Weights_DD_Autocor_DiffTraitCateg.pdf',
    width = 10, height = 10)
print(plot_R2_Tr_precip)
dev.off()


summ_coefs_byResp_Aut_precip <- R2_Aut_precip %>%
  dplyr::group_by(., Response) %>%
  dplyr::summarise(Mean = mean(R.squared, na.rm = T),
                   Median = median(R.squared, na.rm = T))

ggplot(R2_Aut_precip, aes(x = R.squared)) +
  facet_wrap(vars(Response)) +
  geom_histogram() + theme_bw() +
  geom_vline(data = summ_coefs_byResp_Aut_precip, aes(xintercept = Median),
             col = 'red', lwd = 1.4) + theme(strip.text.y = element_text(size = rel(0.8)))

## write out the R2 file
saveRDS(object = R2_Aut_precip, file = './output_fSEM_precip/R2_allMods_Precip_Weights_DD_Autocorr.RDS')


## looking at the path coefficients
# merging predictor and response in one
Coefs_Aut_precip$Relation <- paste(Coefs_Aut_precip$Response,
                                   Coefs_Aut_precip$Predictor, sep = '<-')
table(Coefs_Aut_precip$Demog_rate_Categ)

## split by traits and dem. rates
pdf('./output_SEM_all/PathCoefs_Boxplot_byTraitDemRate_Precip_Weights_DD_Autocor.pdf',
    width = 9)
ggplot(Coefs_Aut_precip, aes(y = Estimate, x = Demog_rate_Categ, col = Trait_Categ)) +
  facet_wrap(vars(Relation)) + geom_boxplot() + theme_bw() +
  geom_hline(yintercept = 0, col = 'black', lty = 2, lwd= 0.7) +
  theme(axis.title = element_text(size = rel(1.4)),
        axis.text = element_text(size = rel(1.1)),
        strip.text = element_text(size = rel(1.4)),
        legend.position = 'bottom', legend.text = element_text(size = rel(1.1))) +
  xlab('Demographic rate')
dev.off()


## write out the coefs file
saveRDS(object = Coefs_Aut_precip, file = './output_fSEM_precip/PathCoefs_allMods_Precip_Weights_DD_Autocor.RDS')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####              2.2. SEMs no DD, weights, standard, autocor                         ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


## fit SEMs for all studies
# Currently fitting SEM for study 581 for species Chen caerulescens caerulescens in La Perouse Bay for LayingDate
# Error in `coef<-.corARMA`(`*tmp*`, value = value[parMap[, i]]) :
#   Coefficient matrix not invertible

# Currently fitting SEM for study 582 for species Chen caerulescens caerulescens in La Perouse Bay for YoungBodyCondition
# Error in nlme::gls(GR ~ det_Clim + Trait_mean, correlation = nlme::corAR1(form = ~Year |  :
#     false convergence (8)
#     Called from: nlme::gls(GR ~ det_Clim + Trait_mean, correlation = nlme::corAR1(form = ~Year |
#                                                                                                                                                           ID), method = "REML", data = dat)

fitted_SEM_Autoc_precip_noDD <- lapply(unique(precip_SEM$ID)[- which(unique(precip_SEM$ID) %in%  c(581, 582))],
                                  FUN = function(x){
                                    fit_SEM(biol_data = precip_SEM, ID = x,
                                            out_SEM = NULL,  #'output_fSEM_noDD_precip'
                                            DD = 'none',
                                            weight = TRUE, correlation = TRUE,
                                            standardize = TRUE, Trait = FALSE,
                                            simpleSEM = TRUE)
                                  })
length(fitted_SEM_Autoc_precip_noDD)  ## 210


# getting the results of SEM
Cstat_Aut_precip_noDD <- extract_res_SEM(list_fitSEM = fitted_SEM_Autoc_precip_noDD, stat_extr = 'Cstat')
R2_Aut_precip_noDD <- extract_res_SEM(list_fitSEM = fitted_SEM_Autoc_precip_noDD, stat_extr = 'R2')
Coefs_Aut_precip_noDD <- extract_res_SEM(list_fitSEM = fitted_SEM_Autoc_precip_noDD, stat_extr = 'coefs')


## write out the Cstat file
saveRDS(object = Cstat_Aut_precip_noDD, file = './output_fSEM_noDD_precip/Cstat_allMods_Precip_Weights_noDD_Autocorr.RDS')


## looking at R2
summ_coefs_Aut_precip_noDD <- R2_Aut_precip_noDD %>%
  dplyr::group_by(., Response, Trait_Categ) %>%
  dplyr::summarise(Mean = mean(R.squared, na.rm = T),
                   Median = median(R.squared, na.rm = T))
summ_coefs_Aut_precip_noDD

p_Aut_noDD <- ggplot(R2_Aut_precip_noDD, aes(x = R.squared)) +
  facet_grid(rows = vars(Response), cols = vars(Trait_Categ)) +
  geom_histogram() + theme_bw() +
  geom_vline(data = summ_coefs_Aut_precip_noDD, aes(xintercept = Median),
             col = 'red', lwd = 1.4) +
  theme(strip.text = element_text(size = rel(1.3))) +
  geom_text(data = summ_coefs_Aut_precip_noDD,
            aes(label  = round(summ_coefs_Aut_precip_noDD$Median, 2)),
            x = round(summ_coefs_Aut_precip_noDD$Median, 2) + 0.12, y = rep(10, 4),
            col = rep('red', 4), size = rep(4, 4))


print(p_Aut_noDD)

## again, R2 explained for DR and GR is much higher if DD is included (so here is lower)

 summ_coefs_byResp_Aut_precip_noDD <- R2_Aut_precip_noDD %>%
  dplyr::group_by(., Response) %>%
  dplyr::summarise(Mean = mean(R.squared, na.rm = T),
                   Median = median(R.squared, na.rm = T))

ggplot(R2_Aut_precip_noDD, aes(x = R.squared)) +
  facet_wrap(vars(Response)) +
  geom_histogram() + theme_bw() +
  geom_vline(data = summ_coefs_byResp_Aut_precip_noDD, aes(xintercept = Median),
             col = 'red', lwd = 1.4) + theme(strip.text.y = element_text(size = rel(0.8)))

## write out the R2 file
saveRDS(object = R2_Aut_precip_noDD, file = './output_fSEM_noDD_precip/R2_allMods_Precip_Weights_noDD_Autocorr.RDS')


## looking at the path coefficients
Coefs_Aut_precip_noDD$Relation <- paste(Coefs_Aut_precip_noDD$Response,
                                        Coefs_Aut_precip_noDD$Predictor, sep = '<-')
table(Coefs_Aut_precip_noDD$Relation)

## split by traits and dem. rates
pdf('./output_SEM_all/PathCoefs_Boxplot_byTraitDemRate_Precip_Weights_noDD_Autocor.pdf',
    width = 9)
ggplot(Coefs_Aut_precip_noDD, aes(y = Estimate, x = Demog_rate_Categ, col = Trait_Categ)) +
  facet_wrap(vars(Relation)) + geom_boxplot() + theme_bw() +
  geom_hline(yintercept = 0, col = 'black', lty = 2, lwd= 0.7) +
  theme(axis.title = element_text(size = rel(1.4)),
        axis.text = element_text(size = rel(1.1)),
        strip.text = element_text(size = rel(1.4)),
        legend.position = 'bottom', legend.text = element_text(size = rel(1.1))) +
  xlab('Demographic rate')
dev.off()


## write out the coefs file
saveRDS(object = Coefs_Aut_precip_noDD, file = './output_fSEM_noDD_precip/PathCoefs_allMods_Precip_Weights_noDD_Autocor.RDS')



# 3. Some summary for text ------------------------------------------------

R2_Aut_proc$Climate <- 'Temperature'
R2_Aut_precip$Climate <- 'Precipitation'
all_R2 <- rbind(R2_Aut_proc, R2_Aut_precip)


sum_R2 <- all_R2 %>%
  dplyr::group_by(Response, Trait_Categ) %>%
  dplyr::summarise(Median = median(R.squared))
