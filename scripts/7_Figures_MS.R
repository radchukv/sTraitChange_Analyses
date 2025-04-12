## The script to produce the figures included in the MS

library(ggplot2)
library(patchwork)
library(tidyverse)
library(magrittr)
library(sTraitChange)
library(ggtext) ## needed for element_markdown() for the labs



# 1. Figure for SEMs, with temp and precip  --------

## read in the effect sizes on temperature
ef_T <- readRDS(file = './output_all/all_efSizes_temperature.RDS')
ef_T$Climatic_var <- 'Temperature'

## effect sizes precipitation
ef_P <- readRDS(file = './output_all/all_efSizes_Precip.RDS')
ef_P$Climatic_var <- 'Precipitation'

ef_all <- rbind(ef_T, ef_P)


## here save ef_all with stats for reporting in the MS - this has to be re-run YET !!!!
ef_all_perModelstats <- ef_all %>%
  dplyr::filter(., Variable == 'intrcpt') %>%
  dplyr::select(., Climatic_var, Trait_Categ, REL,
                AIC_EfS_Covar, Species.SD, ID.SD,
                Location.SD, Phylo.SD, Chi2, pval_Covar)

ef_all_perModelstats <- ef_all_perModelstats %>%
  dplyr::rename(., Climate = Climatic_var) %>%
  dplyr::mutate(AIC_EfS_Covar = format(round(AIC_EfS_Covar, 3), nsmall = 3, scientific = FALSE),
                dplyr::across(.cols = tidyr::ends_with('SD'),
                              .fns = ~ format(round(.x, 4), nsmall = 4, scientific = FALSE)))

save_xlsx(table = ef_all_perModelstats, table_name = './tables/StatsPerModel_AllGroups')

## and saving the stats per parameter
ef_all_perParstats <- ef_all %>%
  dplyr::select(., Climatic_var, Trait_Categ, REL,
                Variable, Estimate, SError,
                EfS_Low, EfS_Upper, Chi2, pval_Covar) %>%
  dplyr::rename(., Climate = Climatic_var, `Trait category` = Trait_Categ,
                CI.lb = EfS_Low, CI.ub = EfS_Upper,
                pval = pval_Covar) %>%
  dplyr::mutate(dplyr::across(.cols= c(Estimate, SError, dplyr::starts_with('CI')),
                              .fns = ~format(round(.x, 3), nsmall = 3, scientific = FALSE)),
                Chi2 = format(round(Chi2, 2), nsmall = 2, scientific = FALSE))

ef_all_perParstats$`p value` <- numeric(length = nrow(ef_all_perParstats))
for(i in 1:nrow(ef_all_perParstats)){
  if (ef_all_perParstats$pval[i] < 0.0001){
    ef_all_perParstats$`p value`[i] <- '<0.0001'
  } else {
    ef_all_perParstats$`p value`[i] <- format(round(ef_all_perParstats$pval[i], 4), nsmall = 4, scientific = FALSE)

  }
}

ef_all_perParstats$pval <- NULL
save_xlsx(table = ef_all_perParstats, table_name = './tables/StatsPerParam_AllGroups')



## continue with data prep for the plot
ef_all <- droplevels(ef_all %>%
  dplyr::filter(., Variable == 'intrcpt'))

ef_all$Climatic_var <- factor(ef_all$Climatic_var, levels =c('Temperature', 'Precipitation'))


## sort the data frme to have the levels in the desired order
ord <- order(ef_all$Trait_Categ, ef_all$REL)
ef_all <- ef_all[ord, ]

ef_all$yAx <- rep(rep(1:6, each = 2), 2)

## setting the Lower estimate to 0 for Z->G, for precipitation and Morp  as p value is 0.059 (will have to separate it clearly into
## marginally-signif. later on)
ef_all$EfS_Low[ef_all$Trait_Categ == 'Phenology' & ef_all$REL == 'ZG' &
                 ef_all$Climatic_var == 'Precipitation'] <- 0


ef_all <- ef_all %>%
  dplyr::mutate(., COL = paste(col_var, Climatic_var, sep = ', '),
                Col = dplyr::case_when(
                  COL == 'Intercept significant, Temperature' ~ 'Signif, temperature',
                  COL == "Intercept significant, Precipitation" ~ 'Signif, precipitation',
                  COL == "Intercept non-significant, Temperature" ~ 'Non-signif, temperature',
                  COL == "Intercept non-significant, Precipitation" ~ 'Non-signif, precipitation'),
                REL_Clim = paste(REL, Climatic_var, sep = '_'))
ef_all$REL_Clim <- factor(ef_all$REL_Clim, levels = c("CZ_Temperature", "CZ_Precipitation",
                                                      "ZG_Temperature", "ZG_Precipitation",
                                                      "CG_Temperature", "CG_Precipitation",
                                                      "PG_Temperature", "PG_Precipitation",
                                                      "CZG_Temperature", "CZG_Precipitation",
                                                      "TotalCG_Temperature", "TotalCG_Precipitation"))



## read in the data with all the effect sizes for the  original studies
allEs_T <- readRDS(file = './output_all/PathCoefs_Temp_AlsoEstimatedRelations.RDS')
allEs_T$Climatic_var <- 'Temperature'

## abd also precip
allEs_P <- readRDS(file = './output_all/PathCoefs_Precip_AlsoEstimatedRelations.RDS')
allEs_P$Climatic_var <- 'Precipitation'

## recode Relation to same levels as REL
all_Es <- rbind(allEs_T, allEs_P)

all_Es %<>%
  dplyr::mutate(REL = dplyr::case_when(Relation == 'Trait_mean<-det_Clim' ~ 'CZ',
                                       Relation == 'GR<-Trait_mean' ~ 'ZG',
                                       Relation == 'GR<-det_Clim' ~ 'CG',
                                       Relation == 'GR<-Pop_mean' ~ 'PG',
                                       Relation == 'Ind_GR<-det_Clim' ~ 'CZG',
                                       Relation == 'Tot_GR<-det_Clim' ~ 'TotalCG'),
                Trait_Categ = dplyr::recode(Trait_Categ, "Morphological" =  "Morphology",
                                         "Phenological" = "Phenology"),
                REL_Clim = paste(REL, Climatic_var, sep = '_'))
## setting the levels of the RelClim the way I need
all_Es$REL_Clim <- factor(all_Es$REL_Clim, levels = c("CZ_Temperature", "CZ_Precipitation",
                                                      "ZG_Temperature", "ZG_Precipitation",
                                                      "CG_Temperature", "CG_Precipitation",
                                                      "PG_Temperature", "PG_Precipitation",
                                                      "CZG_Temperature", "CZG_Precipitation",
                                                      "TotalCG_Temperature", "TotalCG_Precipitation"))
all_Es$REL <- factor(all_Es$REL, levels = c("CZ", "ZG", "CG", "PG", "CZG", "TotalCG"))

## labels for the y axis
axisLab <- data.frame(x = rep(-2.3, length(unique(all_Es$REL_Clim))),
                      y = c(0, NA, 1, NA, 2, NA, 4, NA,  4.8, NA, 6, NA),
                      label = rep(levels(all_Es$REL),  each = 2),
                      Climatic_var = rep(unique(all_Es$Climatic_var),
                                         length(levels(all_Es$REL))),
                      REL_Clim = factor(levels(all_Es$REL_Clim), levels= levels(all_Es$REL_Clim)))
axisLab <- na.omit(axisLab)

## add the number of studies for each combi on the panels
unique_studiesT <- allEs_T %>%
  dplyr::distinct(., ID, .keep_all = TRUE)
tab_Tr_Temp <- table(unique_studiesT$Trait_Categ)
tab_Tr_Temp

unique_studiesP <- allEs_P %>%
  dplyr::distinct(., ID, .keep_all = TRUE)
tab_Tr_Prec <- table(unique_studiesP$Trait_Categ)
tab_Tr_Prec



# temperature
tab_T <- as.data.frame(tab_Tr_Temp)
colnames(tab_T) <- c('Trait_Categ', 'Count')
tab_T$Trait_Categ <- as.character(tab_T$Trait_Categ)
tab_T$Trait_Categ[tab_T$Trait_Categ == 'Morphological'] <- 'Morphology'
tab_T$Trait_Categ[tab_T$Trait_Categ == 'Phenological'] <- 'Phenology'
tab_T$Climatic_var <- 'Temperature'
tab_T$x <- rep(1.8, 2)
tab_T$y <- rep(25, 2)

## precipitation
tab_P <- as.data.frame(tab_Tr_Prec)
colnames(tab_P) <- c('Trait_Categ', 'Count')
tab_P$Trait_Categ <- as.character(tab_P$Trait_Categ)
tab_P$Trait_Categ[tab_P$Trait_Categ == 'Morphological'] <- 'Morphology'
tab_P$Trait_Categ[tab_P$Trait_Categ == 'Phenological'] <- 'Phenology'
tab_P$x <- rep(1.8, 2)
tab_P$y <- rep(5, 2)
tab_P$Climatic_var <- 'Precipitation'

tab_clim <- rbind(tab_T, tab_P)
tab_clim$REL <- 'TotalCG'
tab_clim$REL_Clim <- 'TotalCG_Precipitation'
tab_clim$REL_Clim <- factor(tab_clim$REL_Clim, levels  =
                              c("CZ_Temperature", "CZ_Precipitation",
                                "ZG_Temperature", "ZG_Precipitation",
                                "CG_Temperature", "CG_Precipitation",
                                "PG_Temperature", "PG_Precipitation",
                                "CZG_Temperature", "CZG_Precipitation",
                                "TotalCG_Temperature", "TotalCG_Precipitation"))

label_facet <- data.frame(lab = c('a', 'b'),
                          Trait_Categ = c('Temperature', 'Precipitation'),
                          x = rep(-1.8, 2), y = rep(12, 2))


## applying a function to get the effect sizes per a combination of dem. and trait category
p_Phen <- plot_hist_points(data_allEstim = all_Es,
                             Ef_sizesEstim = ef_all,
                             dataAxis = axisLab,
                             Traitdem = NULL,
                             tabC = tab_clim, annot = NULL,
                             Trait_categ = 'Phenology')
p_Morph <- plot_hist_points(data_allEstim = all_Es,
                             Ef_sizesEstim = ef_all,
                             dataAxis = axisLab,
                             Traitdem = NULL,
                             tabC = tab_clim,
                            Trait_categ = 'Morphology')



pdf('./plots_ms/FigS9_DistributionEffectSizes_PerTrait&Climate.pdf', width = 9, height = 7)
p_Morph + p_Phen +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = 'bottom',
        plot.margin = unit(c(0.3, 0.3, 0, 0), 'line'),
        plot.tag = element_text(face = 'bold'))
dev.off()


# 1.1. Additional analyses: correlation between CZ and ZG -----------------
cz_ef_temp <- subset(all_Es, Relation == 'Trait_mean<-det_Clim' & Climatic_var == 'Temperature')
# this dataset contains phenol and morph studies in response to temperature
cz_ef_prec <- subset(all_Es, Relation == 'Trait_mean<-det_Clim' & Climatic_var == 'Precipitation')
# this dataset contains phenol and morph studies in response to precip

zg_ef_temp <- subset(all_Es, Relation == 'GR<-Trait_mean' & Climatic_var == 'Temperature')
zg_ef_prec <- subset(all_Es, Relation == 'GR<-Trait_mean' & Climatic_var == 'Precipitation')

cz_zg_temp <- merge(cz_ef_temp, zg_ef_temp, by = 'ID')
phen_cz_zg_temp <- subset(cz_zg_temp, Trait_Categ.x == 'Phenology')
cor.test(x = phen_cz_zg_temp$Estimate.x, phen_cz_zg_temp$Estimate.y)
# positive significant correlation between slopes: meaning if the cz slopes are
# positive then zg wil also be positive. OR: if ZG are negative, zg will also be negative
# this of course is reflected in Fig. 3 in the panel of CZ vs ZG
# r  = 0.3269118; p-value = 0.001382, df= 91

# morphology and temp:
morph_cz_zg_temp <- subset(cz_zg_temp, Trait_Categ.x == 'Morphology')
cor.test(x = morph_cz_zg_temp$Estimate.x, morph_cz_zg_temp$Estimate.y)
# non-signif: r = 0.04719467, p = 0.626, df = 107

# look at phenology relation with precipitation
cz_zg_prec <- merge(cz_ef_prec, zg_ef_prec, by = 'ID')
phen_cz_zg_prec <- subset(cz_zg_prec, Trait_Categ.x == 'Phenology')
cor.test(x = phen_cz_zg_prec$Estimate.x, phen_cz_zg_prec$Estimate.y)
# non-significant, but an indication of the negative weak correlation
# r =-0.1437449 , df = 93, p = 0.1646

# and morphology with precip:
morph_cz_zg_prec <- subset(cz_zg_prec, Trait_Categ.x == 'Morphology')
cor.test(x = morph_cz_zg_prec$Estimate.x, morph_cz_zg_prec$Estimate.y)
# non-signif, positive: 0.1817777, p = 0.05186, df = 113

# 2. SA: Figure sensitivity to DD  --------

## read in the effect sizes on temperature without DD
ef_T_noDD <- readRDS(file = './output_all_noDD/all_efSizes_noDD_temp.RDS')
ef_T_noDD$Climatic_var <- 'Temperature'

## effect sizes precipitation
ef_P_noDD <- readRDS(file = './output_all_noDD/all_efSizes_noDD_precip.RDS')
ef_P_noDD$Climatic_var <- 'Precipitation'

ef_all_noDD <- rbind(ef_T_noDD, ef_P_noDD)


## continue with data prep for the plot
ef_all_noDD <- droplevels(ef_all_noDD %>%
                       dplyr::filter(., Variable == 'intrcpt')) %>%
                      mutate(across(where(is_character), as_factor))



## sort the data frme to have the levels in the desired order
ord <- order(ef_all_noDD$Trait_Categ, ef_all_noDD$REL)
ef_all_noDD <- ef_all_noDD[ord, ]


ef_all_noDD <- ef_all_noDD %>%
  dplyr::mutate(., COL = paste(col_var, Climatic_var, sep = ', '),
                Col = dplyr::case_when(
                  COL == 'Intercept significant, Temperature' ~ 'Signif, temperature',
                  COL == "Intercept significant, Precipitation" ~ 'Signif, precipitation',
                  COL == "Intercept non-significant, Temperature" ~ 'Non-signif, temperature',
                  COL == "Intercept non-significant, Precipitation" ~ 'Non-signif, precipitation'),
                REL_Clim = paste(REL, Climatic_var, sep = '_'))
ef_all_noDD$REL_Clim <- factor(ef_all_noDD$REL_Clim, levels = c("CZ_Temperature", "CZ_Precipitation",
                                                      "ZG_Temperature", "ZG_Precipitation",
                                                      "CG_Temperature", "CG_Precipitation",
                                                      "PG_Temperature", "PG_Precipitation",
                                                      "CZG_Temperature", "CZG_Precipitation",
                                                      "TotalCG_Temperature", "TotalCG_Precipitation"))

## labels
labs_df <- data.frame(x = rep(-2.3, 6), y = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5),
                      Trait_Categ = rep(c('Morphology'),  6),
                      label =  as.character(levels(ef_all$REL)))
tab_T$x <- rep(1.2, 2)
tab_T$y <- rep(12, 2)
tab_P$x <- rep(1.2, 2)
tab_P$y <- rep(11.5, 2)
tab_T$Trait_Categ <- as.factor(tab_T$Trait_Categ)

axisLab <- data.frame(x = c(-3.3), y = 6, label = 'Relation',
                      Trait_Categ = c('Morphology'))
label_facet <- data.frame(lab = c('a', 'b'),
                          Trait_Categ = levels(tab_T$Trait_Categ),
                          x = rep(-1.8, 2), y = rep(12, 2))

## merge ef_all with and without DD
ef_all_noDD$DD <- 'No'
ef_all$DD <- 'Yes'

ef_all$yAx <- rep(1:12, 2)
ef_all_noDD$yAx <- rep(c(1:6, 9:12), 2)

ef_all_DD <- rbind(ef_all_noDD, ef_all)


## for these plots consider adding either shading fro each par-r or the lines separating each par-r, so that the points can be linked to the labels easily
Fig_DD <- ggplot(ef_all_DD, aes(x = Estimate, y = yAx)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey', lwd = 1.1) +
  geom_errorbar(width=.1,
                aes(xmin = EfS_Low, xmax = EfS_Upper, colour = Col, linetype = DD), lwd = 0.4) +
  geom_point(aes(shape = as.factor(DD), color = Col), alpha = 0.8, cex= 4) +  ## did not succeed with jittering both the lines and the points
  labs(x = 'Effect size', y = 'Relation') +
  facet_grid(. ~ Trait_Categ) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size = 12),
        plot.margin = margin(0.2, 0.2, 1.5, 4.5, unit = 'line'),
        legend.position = 'bottom',
        strip.background = element_blank(),
        legend.box="vertical", legend.margin=margin(),
        legend.spacing.x  = unit(0, 'line'),
        legend.spacing.y  = unit(0, 'line'),
        legend.box.margin = margin(0,0,0,0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  +
  geom_text(data = labs_df, aes(x = x, y =y, label = label), hjust=1) +
  coord_cartesian(clip = 'off', xlim = c(-2, 1.5)) +
  scale_colour_manual(values = c("Signif, temperature" = "brown2",
                                 "Non-signif, temperature" = "lightsalmon1",
                                 "Signif, precipitation" = "royalblue1",
                                 "Non-signif, precipitation" = "lightblue1"),
                      name = 'Significance') +
  scale_shape_manual(values = c(21, 4), name = 'DD') +
  scale_linetype_manual(values = c('dotted', 'solid'), name = 'DD') +
  guides(color=guide_legend(nrow=2, byrow = TRUE)) +
  geom_text(data = tab_T, aes(x = x, y= y, label = Count), col = "brown2") +
  geom_text(data = tab_P, aes(x = x, y= y, label = Count), col = "royalblue1") +
  geom_text(data = label_facet, aes(x = x, y = y, label = lab), fontface = 'bold') +
  geom_text(data = axisLab, aes(x = x, y = y, label= label), angle = 90)

## output the figure
pdf('plots_ms/FigS13_SensitivityDD.pdf', width = 8, height = 6)
print(Fig_DD)
dev.off()



# 3. SA: Figure sensitivity to P value  --------
## read in the effect sizes on temperature without DD
ef_T_noP <- readRDS(file = './output_all/all_efSizes_noPvalue_temperature.RDS')
ef_T_noP$Climatic_var <- 'Temperature'

## effect sizes precipitation
ef_P_noP <- readRDS(file = './output_all/all_efSizes_noPvalue_Precip.RDS')
ef_P_noP$Climatic_var <- 'Precipitation'

ef_all_noP <- rbind(ef_T_noP, ef_P_noP)
colnames(ef_all_noP)[1] <- 'Variable'


## continue with data prep for the plot
ef_all_noP <- droplevels(ef_all_noP %>%
                            dplyr::filter(., Variable == 'intrcpt')) %>%
  mutate(across(where(is_character), as_factor))


## sort the data frme to have the levels in the desired order
ord <- order(ef_all_noP$Trait_Categ, ef_all_noP$REL)
ef_all_noP <- ef_all_noP[ord, ]


ef_all_noP <- ef_all_noP %>%
  dplyr::mutate(., COL = paste(col_var, Climatic_var, sep = ', '),
                Col = dplyr::case_when(
                  COL == 'Intercept significant, Temperature' ~ 'Signif, temperature',
                  COL == "Intercept significant, Precipitation" ~ 'Signif, precipitation',
                  COL == "Intercept non-significant, Temperature" ~ 'Non-signif, temperature',
                  COL == "Intercept non-significant, Precipitation" ~ 'Non-signif, precipitation'),
                REL_Clim = paste(REL, Climatic_var, sep = '_'))
ef_all_noP$REL_Clim <- factor(ef_all_noP$REL_Clim, levels = c("CZ_Temperature", "CZ_Precipitation",
                                                                "CZG_Temperature", "CZG_Precipitation",
                                                                "TotalCG_Temperature", "TotalCG_Precipitation"))

## subset ef_all to only have the same relations that are relevant to the sensitivity analysis without Pval
ef_all_subP <- droplevels(subset(ef_all, REL %in% unique(ef_all_noP$REL)))

## labels
labs_df <- data.frame(x = rep(-2.3, 3), y = c(1.5, 3.5, 5.5),
                      Trait_Categ = rep(c('Morphology'),  3),
                      label =  as.character(levels(ef_all_subP$REL)))
tab_T$x <- rep(1.2, 2)
tab_T$y <- rep(6, 2)
tab_P$x <- rep(1.2, 2)
tab_P$y <- rep(5.5, 2)

axisLab <- data.frame(x = c(-3.3), y = 3, label = 'Relation',
                      Trait_Categ = c('Morphology'))
label_facet <- data.frame(lab = c('a', 'b'),
                          Trait_Categ = levels(tab_T$Trait_Categ),
                          x = rep(-1.8, 2), y = rep(6, 2))

## merge ef_all with and without DD

ef_all_noP$PdeltaAICc <- 'No'
ef_all_subP$PdeltaAICc <- 'Yes'
ef_all_subP$DD <- NULL

ef_all_subP$yAx <- rep(1:6, 2)
ef_all_noP$yAx <- rep(1:6, 2)

ef_all_P <- rbind(ef_all_noP, ef_all_subP)


Fig_Pval <- ggplot(ef_all_P, aes(x = Estimate, y = yAx)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey', lwd = 1.1) +
  geom_errorbar(width=.1,
                aes(xmin = EfS_Low, xmax = EfS_Upper, colour = Col, linetype = PdeltaAICc), lwd = 0.4) +
  geom_point(aes(shape = as.factor(PdeltaAICc), color = Col), alpha = 0.8, cex = 4) +
  labs(x = 'Effect size', y = 'Relation') +
  facet_grid(. ~ Trait_Categ) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size = 12),
        plot.margin = margin(0.2, 0.2, 1.5, 4.5, unit = 'line'),
        legend.position = 'bottom',
        strip.background = element_blank(),
        legend.box="vertical", legend.margin=margin(),
        legend.spacing.x  = unit(0, 'line'),
        legend.spacing.y  = unit(0, 'line'),
        legend.box.margin = margin(0,0,0,0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  +
  geom_text(data = labs_df, aes(x = x, y =y, label = label), hjust=1) +
  coord_cartesian(clip = 'off', xlim = c(-2, 1.5)) +
  scale_colour_manual(values = c("Signif, temperature" = "brown2",
                                 "Non-signif, temperature" = "lightsalmon1",
                                 "Signif, precipitation" = "royalblue1",
                                 "Non-signif, precipitation" = "lightblue1"),
                      name = 'Significance') +
  scale_shape_manual(values = c(21, 4), name = 'PdeltaAICc') +
  scale_linetype_manual(values = c('dotted', 'solid'), name = 'PdeltaAICc') +
  guides(color=guide_legend(nrow=2, byrow = TRUE)) +
  geom_text(data = tab_T, aes(x = x, y= y, label = Count), col = "brown2") +
  geom_text(data = tab_P, aes(x = x, y= y, label = Count), col = "royalblue1") +
  geom_text(data = label_facet, aes(x = x, y = y, label = lab), fontface = 'bold') +
  geom_text(data = axisLab, aes(x = x, y = y, label= label), angle = 90)



pdf('./plots_ms/FigS12_Sensitivity_PdeltaAICc.pdf', width = 8, height = 5)
print(Fig_Pval)
dev.off()



# 4.  Results Plots ------------------
## temperature

# 1. loading the raw data from SEMs
## 1.1. read in the raw data for temperature
temp_SEM <- readRDS(file = './output_forSEM_temp/all_SEM.RDS')
head(temp_SEM)

## add the years where they are simply omitted as rows
allYrs_T <- do.call('rbind', lapply(X = unique(temp_SEM$ID), FUN = function(x){
  subs <- droplevels(temp_SEM[temp_SEM$ID == x, ])
  full_NA <- data.frame(Year = seq(min(subs$Year), max(subs$Year), by = 1),
                        ID = x)
}))

consec_yrs_T <- merge(allYrs_T, temp_SEM, by = c('ID','Year'), all= T)

## calculate GR
temp_GR <- consec_yrs_T %>%
  dplyr::mutate(., Pop_mean_lag = c(Pop_mean[-1], NA)) %>%
  dplyr::mutate(., GR = log(Pop_mean_lag / Pop_mean)) %>%
  dplyr::filter(., !is.na(GR) & !is.na(Trait_mean) &
                  !is.na(Demog_rate_mean) & !is.na(Pop_mean))


temp_GRRes <- split(temp_GR, temp_GR$ID) %>%
  purrr::map(., ~lm(Clim ~ Year, data = .)) %>%
  purrr::map2(.x = ., .y = split(temp_GR, f= temp_GR$ID),
              .f = ~broom::augment_columns(x = .x, data = .y)) %>%
  dplyr::bind_rows() %>%
  dplyr::select(-.rownames)


temp_std <- temp_GRRes %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(Trait_SE = Trait_SE / sd(Trait_mean, na.rm = T),
                Demog_rate_SE = Demog_rate_SE /sd(Demog_rate_mean, na.rm = T),
                det_Clim = as.numeric(scale(`.resid`)),
                Trait_mean = scale(Trait_mean),
                Demog_rate_mean = scale(Demog_rate_mean),
                Pop_mean = scale(Pop_mean),
                GR = scale(GR)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Climatic_var = 'Temperature')

temp_std$GR <- as.numeric(temp_std$GR[,1])


# 1.2 for precipitation
precip_SEM <- readRDS(file = './output_forSEM_precip/all_SEM.RDS')
head(precip_SEM)

## add the years where they are simply omitted as rows
allYrs_P <- do.call('rbind', lapply(X = unique(precip_SEM$ID), FUN = function(x){
  subs <- droplevels(precip_SEM[precip_SEM$ID == x, ])
  full_NA <- data.frame(Year = seq(min(subs$Year), max(subs$Year), by = 1),
                        ID = x)
}))

consec_yrs_P <- merge(allYrs_P, precip_SEM, by = c('ID','Year'), all= T)

## calculate GR
precip_GR <- consec_yrs_P %>%
  dplyr::mutate(., Pop_mean_lag = c(Pop_mean[-1], NA)) %>%
  dplyr::mutate(., GR = log(Pop_mean_lag / Pop_mean)) %>%
  dplyr::filter(., !is.na(GR) & !is.na(Trait_mean) &
                  !is.na(Demog_rate_mean) & !is.na(Pop_mean))


precip_GRRes <- split(precip_GR, precip_GR$ID) %>%
  purrr::map(., ~lm(Clim ~ Year, data = .)) %>%
  purrr::map2(.x = ., .y = split(precip_GR, f= precip_GR$ID),
              .f = ~broom::augment_columns(x = .x, data = .y)) %>%
  dplyr::bind_rows() %>%
  dplyr::select(-.rownames)


precip_std <- precip_GRRes %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(Trait_SE = Trait_SE / sd(Trait_mean, na.rm = T),
                Demog_rate_SE = Demog_rate_SE /sd(Demog_rate_mean, na.rm = T),
                det_Clim = as.numeric(scale(`.resid`)),
                Trait_mean = scale(Trait_mean),
                Demog_rate_mean = scale(Demog_rate_mean),
                Pop_mean = scale(Pop_mean),
                GR = scale(GR)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Climatic_var = 'Precipitation')
precip_std$GR <- as.numeric(precip_std$GR[,1])




## 2. now prepare effect sizes from each study
# do I really need to keep P.Value in here? CHECK!!!
# 2.1. temperature
wide_temp_all <- allEs_T %>%
  dplyr::select(ID, Estimate, SError, P.Value, Relation) %>%
  tidyr::pivot_wider(id_cols = ID, names_from = Relation,
                     values_from = c(Estimate, SError, P.Value), names_sep = '/')

metaD_temp_all <- allEs_T %>%
  dplyr::distinct(ID, .keep_all = TRUE) %>%
  dplyr::select(-c(Estimate, SError, Relation,
                   P.Value, Pvalue, Count, Nyears,
                   WinDur, Ref.day, Ref.month,
                   WindowClose, deltaAIC, Trait_ageClass,
                   WeathQ, GenLength_y_IUCN, Sp_phylo))

wide_tempES_all <- (merge(wide_temp_all, metaD_temp_all, by = 'ID') %>%
                  dplyr::mutate(Climatic_var = 'Temperature'))


## 2.2 precipitation
wide_precip_all <- allEs_P %>%
  dplyr::select(ID, Estimate, SError, P.Value, Relation) %>%
  tidyr::pivot_wider(id_cols = ID, names_from = Relation,
                     values_from = c(Estimate, SError, P.Value), names_sep = '/')

metaD_prec_all <- allEs_P %>%
  dplyr::distinct(ID, .keep_all = TRUE) %>%
  dplyr::select(-c(Estimate, SError, Relation,
                   P.Value, Pvalue, Count, Nyears,
                   WinDur, Ref.day, Ref.month,
                   WindowClose, deltaAIC, Trait_ageClass,
                   WeathQ, Sp_phylo, GenLength_y_IUCN))

wide_precES_all <- (merge(wide_precip_all, metaD_prec_all, by = 'ID') %>%
                      dplyr::mutate(Climatic_var = 'Precipitation'))

## 3. and load the maeta-analytic across-study ES estimates (meta-analyses not split by the sign of the clim effect)
# 3.1 Temperature
globES_T <-  ef_T %>%
  dplyr::filter(Variable == 'intrcpt') %>%
  dplyr::mutate(Trait_Categ = dplyr::recode(Trait_Categ,
                                     'Phenology' = 'Phenological',
                                     'Morphology' = 'Morphological'))

# 3.2. Precipitation
globES_P <- ef_P %>%
  dplyr::filter(Variable == 'intrcpt') %>%
  dplyr::mutate(Trait_Categ = dplyr::recode(Trait_Categ,
                                     'Phenology' = 'Phenological',
                                     'Morphology' = 'Morphological'))


# 4.1. Temp + Phen --------------------------------------------------------

PhenT_CZ <- plot_concept(Trait_categ = 'Phenological',
                         raw_dat = temp_std,
                         GlobES_dat = globES_T,
                         ES_dat = wide_tempES_all,
                         path = 'CZ',
                         xvar_raw = 'det_Clim',
                         yvar_raw = 'Trait_mean',
                         slope_ES = 'Estimate/Trait_mean<-det_Clim',
                         ylab = 'Phenology, Z',
                         xlab = 'Temperature, C')
PhenT_ZG <- plot_concept(Trait_categ = 'Phenological',
                         raw_dat = temp_std,
                         GlobES_dat = globES_T,
                         ES_dat = wide_tempES_all,
                         path = 'ZG',
                         xvar_raw = 'Trait_mean',
                         yvar_raw = 'GR',
                         slope_ES = 'Estimate/GR<-Trait_mean',
                         ylab = 'Population growth rate, G',
                         xlab = 'Phenology, Z | Temperature, C')
PhenT_CG <- plot_concept(Trait_categ = 'Phenological',
                         raw_dat = temp_std,
                         GlobES_dat = globES_T,
                         ES_dat = wide_tempES_all,
                         path = 'CG',
                         xvar_raw = 'det_Clim',
                         yvar_raw = 'GR',
                         slope_ES = 'Estimate/GR<-det_Clim',
                         ylab = 'Population growth rate, G',
                         xlab = 'Temperature, C | Phenology, Z')
PhenT_CZGvsCG <- pl_conc_DirInd(Trait_categ = 'Phenological',
                                GlobES_dat = globES_T,
                                ES_dat = wide_tempES_all,
                                xlab = 'Phenology-mediated effect <br> of temperature on G (CZ)',
                                ylab = 'Direct effect of temperature on G (CG)')$plot

## plot of CZ vs ZG
##1. subset for the needed trait category
ES_datPhen <- subset(wide_tempES_all, Trait_Categ == 'Phenological')


## original, all points in black
## rectangles for expectation
rect <- data.frame(xmin = c(0, -Inf), xmax = c(Inf, 0),
                   ymin = c(0, -Inf), ymax = c(Inf, 0),
                   id = c(1, 1))
rect$`Estimate/Trait_mean<-det_Clim` <- c(-1, 1)
rect$`Estimate/GR<-Trait_mean` <- c(-1, 1.5)
PhenT_CZvsZG_bw <- ggplot(ES_datPhen,
                           aes(x = `Estimate/Trait_mean<-det_Clim`,
                               y = `Estimate/GR<-Trait_mean`)) +
  geom_rect(data = rect, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax,
                fill = factor(id)), alpha = 0.45) +
  geom_point(alpha = 0.45, cex= 6) +
  geom_hline(yintercept = 0, col = 'black', lty = 3) +
  geom_vline(xintercept = 0, col = 'black', lty = 3) +
  xlab('Effect of temperature on phenology (CZ)') +
  ylab('Effect of phenology on G (ZG)') +
  ylim(min(ES_datPhen$`Estimate/GR<-Trait_mean`) - 0.1,
       max(ES_datPhen$`Estimate/GR<-Trait_mean`) + 0.1) +
  scale_fill_manual(values= c("#1B9E77","#D95F02" )) +
  theme_bw() + theme(legend.position = 'none',
                     strip.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.text = element_text(size  =12),
                     panel.grid.major = element_blank(),
                     axis.title = element_text(size = 24),
                     axis.text = element_text(size = 20))

  PhenT_CZvsZG_bw

pdf('./output_all/CZvsZG_colourExpect.pdf', height = 7.7, width = 7.7)
  PhenT_CZvsZG_bw
dev.off()

PhenT_CZvsZG_bw_nobg <- ggplot(ES_datPhen,
                            aes(x = `Estimate/Trait_mean<-det_Clim`,
                                y = `Estimate/GR<-Trait_mean`)) +
    geom_point(alpha = 0.45, cex= 6) +
    geom_hline(yintercept = 0, col = 'black', lty = 3) +
    geom_vline(xintercept = 0, col = 'black', lty = 3) +
    xlab('Effect of temperature on phenology (CZ)') +
    ylab('Effect of phenology on G (ZG)') +
    ylim(min(ES_datPhen$`Estimate/GR<-Trait_mean`) - 0.1,
         max(ES_datPhen$`Estimate/GR<-Trait_mean`) + 0.1) +
    theme_bw() + theme(legend.position = 'none',
                       strip.background = element_blank(),
                       panel.grid.minor = element_blank(),
                       strip.text = element_text(size  =12),
                       panel.grid.major = element_blank(),
                       axis.title = element_text(size = 24),
                       axis.text = element_text(size = 20))

pdf('./output_all/CZvsZG_bw.pdf', height = 7, width = 7)
PhenT_CZvsZG_bw_nobg
dev.off()

PhenT_CZvsZG_col <- ggplot(ES_datPhen,
       aes(x = `Estimate/Trait_mean<-det_Clim`,
           y = `Estimate/GR<-Trait_mean`, col = `Estimate/Ind_GR<-det_Clim`)) +
  geom_point(alpha = 0.9, cex= 4) +
  geom_hline(yintercept = 0, col = 'black', lty = 3) +
  geom_vline(xintercept = 0, col = 'black', lty = 3) +
  xlab('Effect of climate on trait (CZ)') +
  ylab('Effect of trait on GR (ZG)') +
  ylim(min(ES_datPhen$`Estimate/GR<-Trait_mean`) - 0.1,
       max(ES_datPhen$`Estimate/GR<-Trait_mean`) + 0.1) +
  theme_bw() + theme(legend.position = 'bottom',
                     strip.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.text = element_text(size  =12),
                     panel.grid.major = element_blank(),
                     axis.title = element_text(size = 24),
                     axis.text = element_text(size = 20),
                     legend.title = element_text(size = 16),
                     legend.text = element_text(size = 16)
                     ) +
  scale_colour_gradient2(low = 'red', mid = 'lightgrey',
                         high = 'blue', midpoint = 0,
                         name = 'CZG estimate')


pdf('./output_all/CZvsZG_colourMagnCZG.pdf')
PhenT_CZvsZG_col
dev.off()


## plotting with two colours: one for CZG >=0 and the other for <0
ES_datPhenT_CZG <- ES_datPhen %>%
  dplyr::mutate(CZG_magn = dplyr::case_when(
    `Estimate/Ind_GR<-det_Clim` < 0 ~ 'Negative',
    `Estimate/Ind_GR<-det_Clim` >= 0 ~ 'Non-negative',
    TRUE ~ 'other'
  )) %>%
  dplyr::mutate(CZG_magn= factor(CZG_magn)) %>%
  dplyr::mutate(CZG_magn = forcats::fct_relevel(CZG_magn, c("Non-negative", "Negative")))

PhenT_CZvsZG_CZGNonneg <- ggplot(ES_datPhenT_CZG,
                               aes(x = `Estimate/Trait_mean<-det_Clim`,
                                   y = `Estimate/GR<-Trait_mean`, col = CZG_magn)) +
  geom_point(alpha = 0.9, cex = 5) +
  geom_hline(yintercept = 0, col = 'black', lty = 3) +
  geom_vline(xintercept = 0, col = 'black', lty = 3) +
  xlab('Effect of temperature <br> on phenology (CZ)') +
  ylab('Effect of phenology on G (ZG)') +
  ylim(min(ES_datPhen$`Estimate/GR<-Trait_mean`) - 0.1,
       max(ES_datPhen$`Estimate/GR<-Trait_mean`) + 0.1) +
  theme_bw() + theme(legend.position = 'bottom',
                     strip.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.text = element_text(size =12),
                     panel.grid.major = element_blank(),
                     axis.title = element_text(size = 20),
                     axis.text = element_text(size = 15),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 20),
                     axis.title.x = element_markdown(),
                     axis.title.y = element_markdown()
                     )  +
  scale_color_brewer(palette = 'Dark2', name = 'Sign of the CZG')
PhenT_CZvsZG_CZGNonneg

pdf('./output_all/CZvsZG_colourSignCZG.pdf')
print(PhenT_CZvsZG_CZGNonneg)
dev.off()


# combining the range of the CZG magnitude (point fill) with the positive/negative
# for outline of the point
PhenT_CZvsZG_mult <- ggplot(ES_datPhenT_CZG,
                                 aes(x = `Estimate/Trait_mean<-det_Clim`,
                                     y = `Estimate/GR<-Trait_mean`, col = CZG_magn,
                                     fill = `Estimate/Ind_GR<-det_Clim`)) +
  geom_point(alpha = 0.9, cex = 4, pch = 21, stroke = 2) +
  geom_hline(yintercept = 0, col = 'black', lty = 3, lwd = 1.5) +
  geom_vline(xintercept = 0, col = 'black', lty = 3, lwd = 1.5) +
  xlab('Effect of temperature <br> on phenology (CZ)') +
  ylab('Effect of phenology on G (ZG)') +
  ylim(min(ES_datPhen$`Estimate/GR<-Trait_mean`) - 0.1,
       max(ES_datPhen$`Estimate/GR<-Trait_mean`) + 0.1) +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size =12),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.key.width=unit(1.5,"cm"),
        axis.title.x = element_markdown())  +
  scale_color_brewer(palette = 'Dark2', name = 'Sign of the CZG') +
  scale_fill_gradient2(low = 'red', mid = 'lightgrey',
                         high = 'blue', midpoint = 0,
                         name = 'CZG estimate')
PhenT_CZvsZG_mult

pdf('./output_all/Fig_CZvsZG_multGuides.pdf', height = 8, width = 7.7)
print(PhenT_CZvsZG_mult)
dev.off()

# Saving the two other panels (on CZ and ZG, now with specific x labels)
pdf('./plots_ms/Fig3_CZPan&ZGPan.pdf', height = 12, width = 7)
PhenT_ZG + PhenT_CZ +
  plot_layout(nrow = 2)
dev.off()

# Saving the plot CG as a regression (and not as a density plot)
pdf('./plots_ms/Fig3_CG_Pan.pdf', height = 6, width = 6)
PhenT_CG
dev.off()


## plotting the density plots for CZG and CG
CZG_dens <- ggplot(subset(wide_tempES_all, Trait_Categ ==
                            'Phenological'), aes(x = `Estimate/Ind_GR<-det_Clim`)) +
  geom_density(fill = 'grey', color = 'darkgrey') +
  geom_vline(data = subset(globES_T, REL == 'CZG' & Trait_Categ == 'Phenological'),
             aes(xintercept = Estimate), lty = 2, lwd = 1.3) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        axis.title.x = element_markdown()) +
  xlab('Phenology-mediated effect of <br> temperature on G (CZG)') +
  ylab('Density')

CG_dens <- ggplot(subset(wide_tempES_all, Trait_Categ ==
                            'Phenological'), aes(x = `Estimate/GR<-det_Clim`)) +
  geom_density(fill = 'grey', color = 'darkgrey') +
  geom_vline(data = subset(globES_T, REL == 'CG' & Trait_Categ == 'Phenological'),
             aes(xintercept = Estimate), lty = 2, lwd = 1.3) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20)) +
  xlab('Direct effect of temperature on G (CG)') +
  ylab('Density')

# Saving the two panels on densities
pdf('./output_all/Fig_CZG&CG_Pan_Densities.pdf', height = 12, width = 7)
CZG_dens + CG_dens +
  plot_layout(nrow = 2)
dev.off()


layout_man <- "
AB
#D
CD
E#
"

# this part of the script will have to be deleted if this is not
# the format of the figure that is included in the MS
pdf('./output_all/Fig_4panels_Concept_PhenT_ColouredByCZGSign.pdf',
    height = 13, width = 13)
PhenT_CZ +
  PhenT_ZG +
  PhenT_CZvsZG_CZGNonneg +
  PhenT_CZGvsCG+
  guide_area() +
  plot_annotation(tag_levels = 'a') +
  plot_layout(design = layout_man,
              widths = c(4,4), heights = c(4, 0.1, 4, 0.2),
              guides = 'collect'
              ) &
  theme(plot.tag = element_text(face = 'bold', size = 20))
dev.off()


# binomial test to assess whehter the proportion of non-negative is
# higher than by chance only
table(ES_datPhenT_CZG$CZG_magn)
binom.test(x = table(ES_datPhenT_CZG$CZG_magn)['Non-negative'],
           n = table(ES_datPhenT_CZG$CZG_magn)['Non-negative'] + table(ES_datPhenT_CZG$CZG_magn)['Negative'],
           p = 0.5,
           alternative = 'greater')
#  number of successes = 56, number of trials = 93, p-value = 0.03069
# 95 percent confidence interval: 0.511718 1.00
# probability of success: 0.6021505

# 4.2. Temp + Morph -------------------------------------------------------

MorphT_CZ <- plot_concept(Trait_categ = 'Morphological',
                                raw_dat = temp_std,
                                GlobES_dat = globES_T,
                                ES_dat = wide_tempES_all,
                                path = 'CZ',
                                xvar_raw = 'det_Clim',
                                yvar_raw = 'Trait_mean',
                                slope_ES = 'Estimate/Trait_mean<-det_Clim',
                                expression(paste('Morphology, ', Z[mo])),
                                xlab = expression(paste('Temperature, ', C[te])))
MorphT_ZG <- plot_concept(Trait_categ = 'Morphological',
                                raw_dat = temp_std,
                                GlobES_dat = globES_T,
                                ES_dat = wide_tempES_all,
                                path = 'ZG',
                                xvar_raw = 'Trait_mean',
                                yvar_raw = 'GR',
                                slope_ES = 'Estimate/GR<-Trait_mean',
                                ylab = 'Population growth rate, G',
                                xlab = expression(paste('Morphology, ', Z[mo], ' | Temperature, ', C[te])))
MorphT_CG <- plot_concept(Trait_categ = 'Morphological',
                         raw_dat = temp_std,
                         GlobES_dat = globES_T,
                         ES_dat = wide_tempES_all,
                         path = 'CG',
                         xvar_raw = 'det_Clim',
                         yvar_raw = 'GR',
                         slope_ES = 'Estimate/GR<-det_Clim',
                         ylab = 'Population growth rate, G',
                         xlab = 'Climate, C | Trait, Z')
MorphT_CZGvsCG <- pl_conc_DirInd(Trait_categ = 'Morphological',
                                GlobES_dat = globES_T,
                                ES_dat = wide_tempES_all,
                                xlab = 'Morphology-mediated  effect <br> of temperature on G (C<sub>te</sub>Z<sub>mo</sub>G)',
                                ylab = 'Direct effect of temperature on G (C<sub>te</sub>G)')$plot

## plot of CZ vs ZG
ES_datMorph <- subset(wide_tempES_all, Trait_Categ == 'Morphological')

MorphT_CZvsZG <- ggplot(ES_datMorph,
                       aes(x = `Estimate/Trait_mean<-det_Clim`,
                           y = `Estimate/GR<-Trait_mean`)) +
  geom_point(colour = 'black', alpha = 0.45) +
  geom_hline(yintercept = 0, col = 'black', lty = 3) +
  geom_vline(xintercept = 0, col = 'black', lty = 3) +
  xlab('Effect of climate on trait (CZ)') +
  ylab('Effect of trait on GR (ZG)') +
  theme_bw() + theme(legend.position = 'bottom',
                     strip.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.text = element_text(size  =12),
                     panel.grid.major = element_blank())


## plot with color for positive vs negative CZG per study
## plotting with two colours: one for CZG >=0 and the other for <0
ES_datMorphT_CZG <- ES_datMorph %>%
  dplyr::mutate(CZG_magn = dplyr::case_when(
    `Estimate/Ind_GR<-det_Clim` <= 0 ~ 'Negative',
    `Estimate/Ind_GR<-det_Clim` > 0 ~ 'Non-negative',
    TRUE ~ 'other'
  )) %>%
  dplyr::mutate(CZG_magn= factor(CZG_magn)) %>%
  dplyr::mutate(CZG_magn = forcats::fct_relevel(CZG_magn, c("Non-negative", "Negative")))

MorphT_CZvsZG_CZGNonneg <- ggplot(ES_datMorphT_CZG,
                                 aes(x = `Estimate/Trait_mean<-det_Clim`,
                                     y = `Estimate/GR<-Trait_mean`, col = CZG_magn)) +
  geom_point(alpha = 0.9, cex = 4) +
  geom_hline(yintercept = 0, col = 'black', lty = 3) +
  geom_vline(xintercept = 0, col = 'black', lty = 3) +
  xlab('Effect of temperature <br> on morphology (C<sub>te</sub>Z<sub>mo</sub>)') +
  ylab('Effect of morphology on G (Z<sub>mo</sub>G)') +
  ylim(min(ES_datPhen$`Estimate/GR<-Trait_mean`) - 0.1,
       max(ES_datPhen$`Estimate/GR<-Trait_mean`) + 0.1) +
  theme_bw() + theme(legend.position = 'bottom',
                     strip.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.text = element_text(size =12),
                     panel.grid.major = element_blank(),
                     axis.title = element_text(size = 20),
                     axis.text = element_text(size = 15),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 20),
                     axis.title.x = element_markdown(),
                     axis.title.y = element_markdown())  +
  scale_color_brewer(palette = 'Dark2', name = 'Sign of the CZG')
MorphT_CZvsZG_CZGNonneg


## final plot for the SI
layout_man <- "
AB
#D
CD
E#
"

pdf('./plots_ms/FigS17_MorphT_CZvsZG_ColouredByCZGSign.pdf',
    height = 13, width = 13)
MorphT_CZ +
  MorphT_ZG +
  MorphT_CZvsZG_CZGNonneg +
  MorphT_CZGvsCG +
  guide_area() +
  plot_annotation(tag_levels = 'a') +
  plot_layout(design = layout_man,
              widths = c(4,4), heights = c(4, 0.1, 4, 0.2),  ## not ideal but the best I could do to not have legend overlap the plot
              guides = 'collect'
  ) &
  theme(plot.tag = element_text(face = 'bold', size = 20))
dev.off()



## stats for studies with non-negaitve CZG being higher than by random
table(ES_datMorphT_CZG$CZG_magn)
binom.test(x = table(ES_datMorphT_CZG$CZG_magn)['Non-negative'],
           n = table(ES_datMorphT_CZG$CZG_magn)['Non-negative'] + table(ES_datMorphT_CZG$CZG_magn)['Negative'],
           p = 0.5,
           alternative = 'greater')
#  number of successes = 53, number of trials = 109, p-value = 0.6491
# 95 percent confidence interval: 0.4039488 1.00
# probability of success: 0.4862385


# 4.3. Precip + Phen ------------------------------------------------------
PhenP_CZ <- plot_concept(Trait_categ = 'Phenological',
                                raw_dat = precip_std,
                                GlobES_dat = globES_P,
                                ES_dat = wide_precES_all,
                                path = 'CZ',
                                xvar_raw = 'det_Clim',
                                yvar_raw = 'Trait_mean',
                                slope_ES = 'Estimate/Trait_mean<-det_Clim',
                                ylab = expression(paste('Phenology, ', Z[ph])),
                                xlab = expression(paste('Precipitation, ', C[pr])))
PhenP_ZG <- plot_concept(Trait_categ = 'Phenological',
                                raw_dat = precip_std,
                                GlobES_dat = globES_P,
                                ES_dat = wide_precES_all,
                                path = 'ZG',
                                xvar_raw = 'Trait_mean',
                                yvar_raw = 'GR',
                                slope_ES = 'Estimate/GR<-Trait_mean',
                                ylab = 'Population growth rate, G',
                                xlab = expression(paste('Phenology, ', Z[ph], ' | Precipitation, ', C[pr])))
PhenP_CG <- plot_concept(Trait_categ = 'Phenological',
                                raw_dat = precip_std,
                                GlobES_dat = globES_P,
                                ES_dat = wide_precES_all,
                                path = 'CG',
                                xvar_raw = 'det_Clim',
                                yvar_raw = 'GR',
                                slope_ES = 'Estimate/GR<-det_Clim',
                                ylab = 'Population growth rate, G',
                                xlab = 'Climate, C')
PhenP_CZGvsCG <- pl_conc_DirInd(Trait_categ = 'Phenological',
                                       GlobES_dat = globES_P,
                                       ES_dat = wide_precES_all,
                                       xlab = 'Phenology-mediated effect of <br> precipitation on G (C<sub>pr</sub>Z<sub>ph</sub>G)',
                                       ylab = 'Direct effect of precipitation on G (C<sub>pr</sub>G)')$plot


## plot of CZ vs ZG
ES_datP_phen <- subset(wide_precES_all, Trait_Categ == 'Phenological')


PhenP_CZvsZG <- ggplot(ES_datP_phen,
                        aes(x = `Estimate/Trait_mean<-det_Clim`,
                            y = `Estimate/GR<-Trait_mean`)) +
  geom_point(colour = 'black', alpha = 0.45) +
  geom_hline(yintercept = 0, col = 'black', lty = 3) +
  geom_vline(xintercept = 0, col = 'black', lty = 3) +
  xlab('Effect of climate on trait (CZ)') +
  ylab('Effect of trait on GR (ZG)') +
  theme_bw() + theme(legend.position = 'bottom',
                     strip.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.text = element_text(size  =12),
                     panel.grid.major = element_blank())

## plot with color for positive vs negative CZG per study
## plotting with two colours: one for CZG >-0 and the other for <0
ES_datPhenP_CZG <- ES_datP_phen %>%
  dplyr::mutate(CZG_magn = dplyr::case_when(
    `Estimate/Ind_GR<-det_Clim` <= 0 ~ 'Negative',
    `Estimate/Ind_GR<-det_Clim` > 0 ~ 'Non-negative',
    TRUE ~ 'other'
  )) %>%
  dplyr::mutate(CZG_magn= factor(CZG_magn)) %>%
  dplyr::mutate(CZG_magn = forcats::fct_relevel(CZG_magn, c("Non-negative", "Negative")))

PhenP_CZvsZG_CZGNonneg <- ggplot(ES_datPhenP_CZG,
                                  aes(x = `Estimate/Trait_mean<-det_Clim`,
                                      y = `Estimate/GR<-Trait_mean`, col = CZG_magn)) +
  geom_point(alpha = 0.9, cex = 4) +
  geom_hline(yintercept = 0, col = 'black', lty = 3) +
  geom_vline(xintercept = 0, col = 'black', lty = 3) +
  xlab('Effect of precipitation <br> on phenology (C<sub>pr</sub>Z<sub>ph</sub>)') +
  ylab('Effect of phenology on G (Z<sub>ph</sub>G)') +
  ylim(min(ES_datPhen$`Estimate/GR<-Trait_mean`) - 0.1,
       max(ES_datPhen$`Estimate/GR<-Trait_mean`) + 0.1) +
  theme_bw() + theme(legend.position = 'bottom',
                     strip.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.text = element_text(size =12),
                     panel.grid.major = element_blank(),
                     axis.title = element_text(size = 20),
                     axis.text = element_text(size = 15),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 20),
                     axis.title.x = element_markdown(),
                     axis.title.y = element_markdown())  +
  scale_color_brewer(palette = 'Dark2', name = 'Sign of the CZG')
PhenP_CZvsZG_CZGNonneg


## final plot for the SI
layout_man <- "
AB
#D
CD
E#
"

pdf('./plots_ms/FigS18_PhenP_CZvsZG_ColouredByCZGSign.pdf',
    height = 13, width = 13)
PhenP_CZ +
  PhenP_ZG +
  PhenP_CZvsZG_CZGNonneg +
  PhenP_CZGvsCG +
  guide_area() +
  plot_annotation(tag_levels = 'a') +
  plot_layout(design = layout_man,
              widths = c(4,4), heights = c(4, 0.1, 4, 0.2),  ## not ideal but the best I could do to not have legend overlap the plot
              guides = 'collect'
  ) &
  theme(plot.tag = element_text(face = 'bold', size = 20),
        #legend.position = 'bottom'
  )
dev.off()


## stats for studies with non-negaitve CZG being higher than by random
table(ES_datPhenP_CZG$CZG_magn)
binom.test(x = table(ES_datPhenP_CZG$CZG_magn)['Non-negative'],
           n = table(ES_datPhenP_CZG$CZG_magn)['Non-negative'] + table(ES_datPhenP_CZG$CZG_magn)['Negative'],
           p = 0.5,
           alternative = 'greater')
#  number of successes = 41, number of trials = 95, p-value = 0.9247
# 95 percent confidence interval: 0.3452896 1.00
# probability of success: 0.4315789


# 4.4. Precip  + Morph ----------------------------------------------------
MorphP_CZ <- plot_concept(Trait_categ = 'Morphological',
                                 raw_dat = precip_std,
                                 GlobES_dat = globES_P,
                                 ES_dat = wide_precES_all,
                                 path = 'CZ',
                                 xvar_raw = 'det_Clim',
                                 yvar_raw = 'Trait_mean',
                                 slope_ES = 'Estimate/Trait_mean<-det_Clim',
                                 ylab = expression(paste('Morphology, ', Z[mo])),
                                 xlab = expression(paste('Precipitation, ', C[pr])))
MorphP_ZG <- plot_concept(Trait_categ = 'Morphological',
                                 raw_dat = precip_std,
                                 GlobES_dat = globES_P,
                                 ES_dat = wide_precES_all,
                                 path = 'ZG',
                                 xvar_raw = 'Trait_mean',
                                 yvar_raw = 'GR',
                                 slope_ES = 'Estimate/GR<-Trait_mean',
                                 ylab = 'Population growth rate, G',
                                 xlab = expression(paste('Morphology, ', Z[mo], ' | Precipitation, ', C[pr])))
MorphP_CG <- plot_concept(Trait_categ = 'Morphological',
                                 raw_dat = precip_std,
                                 GlobES_dat = globES_P,
                                 ES_dat = wide_precES_all,
                                 path = 'CG',
                                 xvar_raw = 'det_Clim',
                                 yvar_raw = 'GR',
                                 slope_ES = 'Estimate/GR<-det_Clim',
                                 ylab = 'Population growth rate, G',
                                 xlab = 'Climate, C')
MorphP_CZGvsCG <- pl_conc_DirInd(Trait_categ = 'Morphological',
                                        GlobES_dat = globES_P,
                                        ES_dat = wide_precES_all,
                                        xlab = 'Morphology-mediated effect of <br> precipitation on G (C<sub>pr</sub>Z<sub>mo</sub>G)',
                                        ylab = 'Direct effect of precipitation on G (C<sub>pr</sub>G)')$plot
## plot of CZ vs ZG
ES_datP_morph <- subset(wide_precES_all, Trait_Categ == 'Morphological')

## maybe pimp up this plot with either density distributions on y axis or Mean + SD over the all data points (negative and positive in x-axis)
## depends on how the final conceptual plot will look like - if that one contains the density distributions, here they should also be used
MorphP_CZvsZG <- ggplot(ES_datP_morph,
                       aes(x = `Estimate/Trait_mean<-det_Clim`,
                           y = `Estimate/GR<-Trait_mean`)) +
  geom_point(colour = 'black', alpha = 0.45) +
  geom_hline(yintercept = 0, col = 'black', lty = 3) +
  geom_vline(xintercept = 0, col = 'black', lty = 3) +
  xlab('Effect of climate on trait (CZ)') +
  ylab('Effect of trait on GR (ZG)') +
  theme_bw() + theme(legend.position = 'bottom',
                     strip.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.text = element_text(size  =12),
                     panel.grid.major = element_blank())


## plot with color for positive vs negative CZG per study
## plotting with two colours: one for CZG >-0 and the other for <0
ES_datMorphP_CZG <- ES_datP_morph %>%
  dplyr::mutate(CZG_magn = dplyr::case_when(
    `Estimate/Ind_GR<-det_Clim` <= 0 ~ 'Negative',
    `Estimate/Ind_GR<-det_Clim` > 0 ~ 'Non-negative',
    TRUE ~ 'other'
  )) %>%
  dplyr::mutate(CZG_magn= factor(CZG_magn)) %>%
  dplyr::mutate(CZG_magn = forcats::fct_relevel(CZG_magn, c("Non-negative", "Negative")))

MorphP_CZvsZG_CZGNonneg <- ggplot(ES_datMorphP_CZG,
                                 aes(x = `Estimate/Trait_mean<-det_Clim`,
                                     y = `Estimate/GR<-Trait_mean`, col = CZG_magn)) +
  geom_point(alpha = 0.9, cex = 4) +
  geom_hline(yintercept = 0, col = 'black', lty = 3) +
  geom_vline(xintercept = 0, col = 'black', lty = 3) +
  xlab('Effect of precipitation <br> on morphology (C<sub>pr</sub>Z<sub>mo</sub>)') +
  ylab('Effect of morphology on G (Z<sub>mo</sub>G)') +
  ylim(min(ES_datPhen$`Estimate/GR<-Trait_mean`) - 0.1,
       max(ES_datPhen$`Estimate/GR<-Trait_mean`) + 0.1) +
  theme_bw() + theme(legend.position = 'bottom',
                     strip.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.text = element_text(size =12),
                     panel.grid.major = element_blank(),
                     axis.title = element_text(size = 20),
                     axis.text = element_text(size = 15),
                     legend.title = element_text(size = 20),
                     legend.text = element_text(size = 20),
                     axis.title.x = element_markdown(),
                     axis.title.y = element_markdown())  +
  scale_color_brewer(palette = 'Dark2', name = 'Sign of the CZG')
MorphP_CZvsZG_CZGNonneg


## final plot for the SI
layout_man <- "
AB
#D
CD
E#
"

pdf('./plots_ms/FigS20_MorphP_CZvsZG_ColouredByCZGSign.pdf',
    height = 13, width = 13)
MorphP_CZ +
  MorphP_ZG +
  MorphP_CZvsZG_CZGNonneg +
  MorphP_CZGvsCG +
  guide_area() +
  plot_annotation(tag_levels = 'a') +
  plot_layout(design = layout_man,
              widths = c(4,4), heights = c(4, 0.1, 4, 0.2),  ## not ideal but the best I could do to not have legend overlap the plot
              guides = 'collect'
  ) &
  theme(plot.tag = element_text(face = 'bold', size = 20))
dev.off()


## stats for studies with non-negaitve CZG being higher than by random
table(ES_datMorphP_CZG$CZG_magn)
binom.test(x = table(ES_datMorphP_CZG$CZG_magn)['Non-negative'],
           n = table(ES_datMorphP_CZG$CZG_magn)['Non-negative'] + table(ES_datMorphP_CZG$CZG_magn)['Negative'],
           p = 0.5,
           alternative = 'greater')
#  number of successes = 59, number of trials = 115, p-value = 0.4261
# 95 percent confidence interval: 0.4324683 1.00
# probability of success: 0.5130435



# 5. Forest plots for SI --------------------------------------

# 5.1. For temperature
ef_T$yAx <- rep(c(1:12), 2)
labs_df <- data.frame(x = rep(-2.3, 6), y = c(2, 4, 5, 6, 8, 11),
                      Trait_Categ = rep(c('Morphology'),  2),
                      label =  as.character(levels(ef_all$REL)))


# modifying tab_T somewhat to produce the plot with intercept only
tab_T$x <- rep(1, 2)
tab_T$y <- rep(12, 2)
tab_T <- tab_T[order(tab_T$Trait_Categ), ]
tab_T$Lab <- letters[1:2]

pdf('./plots_ms/FigS10_Distribution_EfSizes_Covariates_perTraitCateg_Temperature.pdf', width = 9)
ggplot(ef_T, aes(x = Estimate, y = yAx)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'lightgrey', lwd = 1.1) +
  geom_errorbar(width=.1,
                aes(xmin = EfS_Low, xmax = EfS_Upper, colour = col_var)) +
  geom_point(aes(colour = col_var), cex= 4) +
  labs(x = 'Effect size', y = '') +
  facet_grid(cols = vars(Trait_Categ)) +
  theme_bw() + theme(axis.title.y = element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     plot.margin = margin(0.2, 0.2, 1.5, 4, unit = 'line'),
                     legend.position = 'bottom',
                     strip.background = element_blank(),
                     strip.text = element_text(size = 14),
                     axis.text = element_text(color = 'black', size = 11),
                     legend.text = element_text(size = 12),
                     axis.title = element_text(size = 14))  +
  geom_text(data = labs_df, aes(x = x, y =y, label = label), hjust=1,
            size = 4.5) +
  coord_cartesian(clip = 'off', xlim = c(-2, 1.5)) +
  scale_colour_manual(values = c("Intercept non-significant" = "gray",
                                 "Intercept significant" = "black",
                                 "Pvalue non-significant" = "orange",
                                 "Pvalue significant" = "red",
                                 "Weather quality non-significant" = "lightblue",
                                 "Weather quality significant" = "blue")) +
  guides(colour=guide_legend('Legend')) +
  geom_text(data = tab_T, aes(x = x, y= y, label = paste(Count, 'studies'))) +
  geom_text(data = tab_T, aes(y = y, x = -1.9, label= Lab),
            fontface = 'bold', size = 6) +
  guides(color= guide_legend(ncol=2, title = 'Effect', size = 14))
dev.off()


## 5.2. for Precipitation
ef_P$yAx <- rep(c(1:12), 2)

labs_df <- data.frame(x = rep(-2.3, 6), y = c(2, 4, 5, 6, 8, 11),
                      Trait_Categ = rep(c('Morphology'),  2),
                      label =  as.character(levels(ef_all$REL)))

# modifying the x and y coordinates in the summary table on precip such as to
# show the number of studies in this plot
tab_P$x <- rep(1, 2)
tab_P$y <- rep(12, 2)
tab_P$Trait_Categ <- as.factor(tab_P$Trait_Categ)
tab_P <- tab_P[order(tab_P$Trait_Categ), ]
tab_P$Lab <- letters[1:2]

pdf('./plots_ms/FigS11_Distribution_EfSizes_Covariates_perTraitCateg_Precip.pdf', width = 9)
ggplot(ef_P, aes(x = Estimate, y = yAx)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'lightgrey', lwd = 1.1) +
  geom_errorbar(width=.1,
                aes(xmin = EfS_Low, xmax = EfS_Upper, colour = col_var)) +
  geom_point(aes(colour = col_var), cex = 4) +
  labs(x = 'Effect size', y = '') +
  facet_grid(cols = vars(Trait_Categ)) +
  theme_bw() + theme(axis.title.y = element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     plot.margin = margin(0.2, 0.2, 1.5, 4, unit = 'line'),
                     legend.position = 'bottom',
                     strip.background = element_blank(),
                     strip.text = element_text(size = 14),
                     axis.text = element_text(color = 'black', size = 11),
                     legend.text = element_text(size = 12),
                     axis.title = element_text(size = 14))  +
  geom_text(data = labs_df, aes(x = x, y =y, label = label), hjust=1) +
  coord_cartesian(clip = 'off', xlim = c(-2, 1.5)) +
  scale_colour_manual(values = c("Intercept non-significant" = "gray",
                                 "Intercept significant" = "black",
                                 "Pvalue non-significant" = "orange",
                                 "Pvalue significant" = "red",
                                 "Weather quality non-significant" = "lightblue",
                                 "Weather quality significant" = "blue")) +
  guides(colour=guide_legend('Legend')) +
  geom_text(data = tab_P, aes(x = x, y= y, label = paste(Count, 'studies'))) +
  geom_text(data = tab_P, aes(y = y, x = -1.9, label= Lab),
            fontface = 'bold', size = 6) +
  guides(color= guide_legend(ncol=2, title = 'Effect', size = 14))
dev.off()

