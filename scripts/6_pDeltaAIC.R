# The script to assess relations between PDeltaAIC, study duration and the CZ coefficients

# I. Number of years and PDeltaAIC -------------------------------------------

library(ggplot2)
library(tidyverse)
library(magrittr)

temp_SEM <- readRDS(file = './output_forSEM_temp/all_SEM.RDS')
precip_SEM <- readRDS(file = './output_forSEM_precip/all_SEM.RDS')

temp_SEM$Climate <- 'Temperature'
precip_SEM$Climate <- 'Precipitation'

all <- rbind(temp_SEM, precip_SEM)
ind <- all %>%
  dplyr::distinct(., ID, Climate, .keep_all = T)

## and separate ind for temp and precip, to be able to fit the models separately
temp_ind <- temp_SEM %>%
  dplyr::distinct(., ID, .keep_all = T)
precip_ind <- precip_SEM %>%
  dplyr::distinct(., ID, .keep_all = T)


ggplot(ind, aes(NYears, Pvalue)) +
  geom_point(alpha = 0.6) +
  theme_bw() + ylab('PdeltaAICc') +
  facet_wrap(. ~ Climate)



## 1. models for temperature dataset

## a simple linear model with trait category as a covariate
mod_full_temp <- lm(Pvalue ~ NYears*Trait_Categ, data = temp_ind)
summary(mod_full_temp)

## test significance
mod_noInt_temp <- lm(Pvalue ~ NYears + Trait_Categ, data = temp_ind)
anova(mod_noInt_temp, mod_full_temp, test = 'Chisq')  ## non-signif. interaction


summary(mod_noInt_temp)
## an increase in the study duration by one year corresponds to the decline of pDeltaAIC by approx. 0.01 (0.008)

## test significance of Trait Categ
mod_NYrs_temp <- lm(Pvalue ~ NYears, data = temp_ind)
anova(mod_NYrs_temp, mod_noInt_temp)  ## signif. effect of trait categ
## Model 1: Pvalue ~ NYears
# Model 2: Pvalue ~ NYears + Trait_Categ
# Res.Df    RSS Df Sum of Sq      F   Pr(>F)
# 1    211 19.251
# 2    210 18.534  1   0.71689 8.1227 0.004807 **


## test significance of NYears
mod_TrCateg_temp <- lm(Pvalue ~ Trait_Categ, data = temp_ind)
anova(mod_TrCateg_temp, mod_noInt_temp)  ## signif. effect of  Nyears

summary(mod_TrCateg_temp)
hist(temp_ind$Pvalue)


## 2. models for precipitation dataset

## a simple linear model with trait category as a covariate
mod_full_precip <- lm(Pvalue ~ NYears*Trait_Categ, data = precip_ind)
summary(mod_full_precip)

## test significance
mod_noInt_precip <- lm(Pvalue ~ NYears + Trait_Categ, data = precip_ind)
anova(mod_noInt_precip, mod_full_precip)  ## non-signif. interaction


summary(mod_noInt_precip)
## an increase in the study duration by one year corresponds to the decline of pDeltaAIC by approx. 0.003

## test significance of Trait Categ
mod_NYrs_precip <- lm(Pvalue ~ NYears, data = precip_ind)
anova(mod_NYrs_precip, mod_noInt_precip)  ## non-signif. effect of trait categ
summary(mod_NYrs_precip)
## Model 1: Pvalue ~ NYears
## Model 2: Pvalue ~ NYears + Trait_Categ
##   Res.Df    RSS Df Sum of Sq      F  Pr(>F)
# 1    210 17.920
# 2    209 17.608  1   0.31138 3.6958 0.05591 .

## test significance of NYears
mod_TrCateg_precip <- lm(Pvalue ~ Trait_Categ, data = precip_ind)
anova(mod_TrCateg_precip, mod_noInt_precip)  ## non-signif. effect of  Nyears

# Model 1: Pvalue ~ Trait_Categ
# Model 2: Pvalue ~ NYears + Trait_Categ
# Res.Df    RSS Df Sum of Sq      F  Pr(>F)
# 1    210 17.933
# 2    209 17.608  1   0.32495 3.8569 0.05087 .

summary(mod_TrCateg_precip)
hist(precip_ind$Pvalue)

## change the trait category levels order, so that it matches other analyses
ind$Trait_Categ <- factor(ind$Trait_Categ, levels = c('Morphological', 'Phenological'))
ind$Climate <- factor(ind$Climate, levels = c('Precipitation', 'Temperature'))

dat_lab <- data.frame(x = rep(9,4), y = rep(1,4),
                      label = letters[1:4],
                      Trait_Categ = rep(levels(ind$Trait_Categ), 2),
                      Climate = rep(levels(ind$Climate), each = 2),
                      slope = c(coef(mod_NYrs_precip)[[2]], coef(mod_NYrs_precip)[[2]],
                                coef(mod_noInt_temp)[[2]], coef(mod_noInt_temp)[[2]]),
                      intercept = c(coef(mod_NYrs_precip)[[1]], coef(mod_NYrs_precip)[[1]],
                      coef(mod_noInt_temp)[[1]] + coef(mod_noInt_temp)[[3]],
                      coef(mod_noInt_temp)[[1]]))
dat_lab %<>%
  dplyr::mutate(round_slope = round(slope, 3))
dat_lab$expres <- sapply(as.numeric(dat_lab$round_slope), function(x){
  deparse(bquote(paste(beta, ' = ', .(x))))
})

pdf('./plots_ms/FigS12_PdeltaAICc_vs_NumYears_byTraitCateg_Climate.pdf',
    width = 10, height = 10)
ggplot(ind, aes(NYears, Pvalue)) + geom_point(alpha = 0.6) +
  theme_bw() + ylab(bquote(P[Delta*AICc])) +
  facet_grid(Climate ~Trait_Categ) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid.minor = element_blank()) +
  xlab('Study duration (years)') +
  geom_abline(data = dat_lab, aes(slope = slope, intercept = intercept),
              linetype = c('dashed', 'dashed', 'solid', 'solid'), color = 'red', lwd = 1.2) +
  geom_text(data = dat_lab, aes(x = x, y = y, label = label),
            size = 6, fontface = 'bold') +
  geom_text(data = dat_lab, aes(x = 40, y = 0.7, label = expres), parse = TRUE)
dev.off()



# II. Relation between PDeltaAIC and estimates --------------------------------

## loading temperature data
Coefs_Aut <- readRDS(file = './output_fSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')
table(Coefs_Aut$Taxon)
Coefs_Aut$Climate <- 'Temperature'
## loading precipitation data
Coefs_Aut_precip <- readRDS(file = './output_fSEM_precip/PathCoefs_allMods_Precip_Weights_DD_Autocor.RDS')
table(Coefs_Aut_precip$Taxon)
Coefs_Aut_precip$Climate <- 'Precipitation'

## subset to only retain the studies with windur max 52
Coefs_Aut <- Coefs_Aut[Coefs_Aut$WinDur <= 51, ]
length(unique(Coefs_Aut$ID))  ## 202
Coefs_Aut_precip <- Coefs_Aut_precip[Coefs_Aut_precip$WinDur <= 51, ]
length(unique(Coefs_Aut_precip$ID))  ## 210

## replace IceBird with Seabird
Coefs_Aut$BirdType[Coefs_Aut$BirdType == 'IceBird'] <- 'Seabird'
Coefs_Aut_precip$BirdType[Coefs_Aut_precip$BirdType == 'IceBird'] <- 'Seabird'
Coefs_all <- bind_rows(Coefs_Aut, Coefs_Aut_precip)


## plot

## set the trait category levels order, so that it matches other analyses
Coefs_all$Trait_Categ <- factor(Coefs_all$Trait_Categ, levels = c('Morphological', 'Phenological'))
Coefs_all$Climate <- factor(Coefs_all$Climate, levels = c('Precipitation', 'Temperature'))

## split of the relaiton between PDeltaAIC and estimates by trait categ only
dat_lab <- data.frame(x = rep(0, 4), y = rep(1.2, 4),
                      label = letters[1:4],
                      Trait_Categ = rep(levels(Coefs_all$Trait_Categ), 2),
                      Climate = rep(levels(Coefs_all$Climate), each = 2))

CZ_only <- Coefs_all %>%
  filter(Relation == 'Trait_mean<-det_Clim')

pdf('./plots_ms/FigS13_Relation_PValue&Estimate_ByTraitCateg_Climate.pdf',
    width  = 10, height = 10)
CZ_only %>%
  dplyr::filter(.,  Relation == 'Trait_mean<-det_Clim') %>%
  ggplot(., aes(x = Pvalue, y = Estimate)) + geom_point() +
  facet_grid(Climate ~ Trait_Categ) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  theme_bw() +
  xlab(bquote(P[Delta*AICc])) +
  geom_text(data = dat_lab, aes(x = x, y = y, label = label),
            fontface = 'bold', size = 6) +
  theme(panel.grid.major = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12))
dev.off()
