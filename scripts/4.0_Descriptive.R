## This script provides some descriptive stats of the whole dataset
## some of these plots were produced for a ppt and are not relevant for the MS (
## stylized / colourful versions)

library(metafor)
library(ggplot2)
library(tidyverse)


source('./scripts/1.2_DataPrep_Biol.R')


## reading in the data per study returned from merging all the climwin outputs prepared for SEM
temp_SEM <- readRDS(file = './output_forSEM_temp/all_SEM.RDS')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                      Describing the trend in pop size                                     ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## simple linear models of pop over time
fit_lm <- function(df){lm(Pop_mean ~ Year, data = df)}
fit_lm_scale <- function(df){lm(pop_scale ~ Year, data = df)}

trendsPop <- temp_SEM %>%
  dplyr::group_by(., ID) %>%
  dplyr::filter(., ! is.na(Pop_mean)) %>%
  dplyr::mutate(., pop_scale = scale(Pop_mean)) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    LM = purrr::map(purrr::map(data, fit_lm),
                    broom::tidy),
    LM_std = purrr::map(purrr::map(data, fit_lm_scale), broom::tidy)) %>%
  tidyr::unnest(cols= c(LM, LM_std), names_sep = '_') %>%
  dplyr::filter(., LM_term == 'Year') %>%
  dplyr::select(., -c(LM_term, LM_std_term,
                      LM_std_p.value, LM_std_statistic)) %>%
  dplyr::mutate(., Trend = dplyr::case_when(LM_p.value >= 0.05 ~ 'stable',
                                            LM_p.value < 0.05 & LM_estimate > 0 ~ 'positive',
                                            LM_p.value < 0.05 & LM_estimate < 0 ~ 'negative'))

hist(trendsPop$LM_std_estimate, breaks = 20, col = 'lightgrey',
     xlab = 'Population trend over time', cex.lab = 1.4,
     main ='')
abline(v = median(trendsPop$LM_std_estimate), lwd = 2, col = 'blue')


## merging Coefs_Aut and LinMod, to have the pop. trend in the same dataset
dat <- merge(biol_dat, subset(trendsPop, select = c(ID, LM_std_estimate, LM_std_std.error, Trend)), by = c('ID'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                1. Descriptive stat of the complete dataset                                  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# directly here rename that typo, thhalassarsche melanophrys to be correct
dat$Species[dat$Species == 'Thalassarche melanophrys'] <- 'Thalassarche melanophris'
dat <- droplevels(dat)
data_descript <- droplevels(dat %>%
  group_by(, ID) %>%
  mutate(Nyears = n()) %>%
  ungroup() %>%
  dplyr::distinct(., ID, .keep_all = T))

# add the continent
data_descript %<>%
  mutate(Continent = case_when(ID %in% unique(eu_noSea$biol_NY[[1]]$ID) ~ 'Europe',
                   ID %in% unique(us_noSea$biol_NY[[1]]$ID) ~  'North America',
                   ID %in% unique(biol_au_noSea$biol_NY[[1]]$ID) ~ 'Australia',
                   Country %in%  c('New Zealand', 'Australia') ~ 'Australia',
                   Country %in%  c('Sweden', 'Norway', 'Germany',
                                 'Spain', 'Denmark', 'Scotland', 'Svalbard') ~ 'Europe',
                   Country %in% c('Falkland Islands', 'South Georgia',
                                  'Venezuela', 'South Atlantic Ocean') ~ 'South America',
                   Country %in% c('Canada', 'Greenland', 'Mexico', 'USA') ~ 'North America',
                   Country == 'Antarctica' ~ 'Antarctica',
                   Country == 'Taiwan' ~ 'Asia',
                   Country == 'South Africa' ~ 'Africa',
                   .default  = 'other'))

data_descr <- data_descript %>%
  mutate(Taxon = recode(Taxon, Reptilia = 'Reptile'))

# summaries
nrow(data_descr)  ## 213 studies
length(unique(data_descr$Species))  # 73 species
table(data_descr$Taxon)
table(data_descr$Trait_Categ)  ##  97 Phenology and  116 Morphology
table(data_descr$Demog_rate_Categ)
table(data_descr$Trait_Categ, data_descr$Taxon)
length(unique(data_descr$Species[data_descr$Taxon== 'Bird'])) # 53
length(unique(data_descr$Species[data_descr$Taxon== 'Mammal'])) # 10
length(unique(data_descr$Species[data_descr$Taxon== 'Reptile'])) # 7
length(unique(data_descr$Species[data_descr$Taxon== 'Fish'])) # 3


# across countries
tab_country <- table(data_descr$Country)
pdf('./output_all/Descriptive_Countries.pdf')
par(mar = c(9, 4, 2, 2))
barplot(tab_country, col = 'chocolate2', cex.axis = 1.5, cex.names = 0.8, las = 2)
mtext(side = 2, 'Number of studies', cex = 1.2, line = 2.4)
mtext(side = 1, 'Country', cex = 1.2, line = 7.5)
dev.off()


# per population trend category
tab_Trend <- table(data_descr$Trend)
pdf('./output_all/Descriptive_PopTrend.pdf')
par(mar = c(4, 4, 2, 2))
barplot(tab_Trend, col = 'chocolate2', cex.axis = 1.5, cex.names = 1.2)
mtext(side = 2, 'Number of studies', cex = 1.2, line = 2.6)
mtext(side = 1, 'Population trend over years', cex = 1.2, line = 2.5)
dev.off()
table(data_descr$Trait_Categ,  data_descr$Trend)

# 1b. Heterogeneity in the exact morphological traits --------------------
# keeping unique studies on morphology
morph <- droplevels(data_descr %>%
                            filter(Trait_Categ == 'Morphological'))

unique(morph$Trait)
morph_tab <- as.data.frame(table(morph$Trait))

morph_tab[order(morph_tab$Freq, decreasing = TRUE),]

BMass <- morph %>%
  dplyr::filter(Trait %in% grep("Mass", morph_tab$Var1, value = TRUE))
nrow(BMass)  # 65



# 2. Map of the studies: Fig 2 in the MS ---------------------------------------

table(data_descr$Trait_Categ, data_descr$Taxon)
tab_Trait_Tax <- table(data_descr$Trait_Categ, data_descr$Taxon)  ## Dem. rate by trait category: in each cell min 30 studies, so perhaps it is OK to run analyses on these cells

tab_Trait_Tax <- as.data.frame(tab_Trait_Tax)
names(tab_Trait_Tax) <- c('Trait', 'Taxon', 'Count')

# change the order of the levels
tab_Trait_Tax$Taxon <- factor(tab_Trait_Tax$Taxon, levels = c('Bird', 'Fish', 'Mammal','Reptile'))

plot_TaxByTraitCat <- ggplot(tab_Trait_Tax, aes(x = Trait, y = Count, fill = Taxon)) +
  geom_bar(stat= 'identity', position = position_dodge(),
           alpha = 0.7) +
  scale_fill_brewer(palette = 'Dark2') +
  geom_text(aes(label=Count, y= c(90, 55, 16, 16, 16, 48, 16, 16)), vjust=1.6, color="black",
            position = position_dodge(0.9), size= rel(4)) +
  theme_bw() + theme(axis.text = element_text(size = rel(1),
                                              colour = 'black'),
                     axis.title = element_text(size = rel(1)),
                     legend.position = 'none',
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     legend.text = element_text(size= rel(0.8)),
                     legend.title = element_text(size = rel(0.8),
                                                 hjust = 0.5),
                     plot.margin = unit(c(0,0,0,0), "cm"),
                     legend.key.size = unit(0.8,"line"),
                     legend.spacing.x  = unit(0, 'line'),
                     legend.spacing.y  = unit(0, 'line'),
                     legend.box.margin = margin(-12,-10,-1,-10),
                     legend.background = element_rect(fill= "transparent")) +
  guides(fill=guide_legend(ncol=2))

pdf('./output_all/Descriptive_Barplot_TraitByTaxonPPT.pdf', width = 10)
plot_TaxByTraitCat
dev.off()


# map of all studies
mp <- NULL
mapWorld <- borders("world", colour="white", fill="gray70") # create a layer of borders
mp <- ggplot() + mapWorld + theme_bw()


mapSubs <- mapWorld$data[mapWorld$data$long >-140 & mapWorld$data$long <= 180, ]
# Create a base plot with gpplot2
p <- ggplot() + coord_fixed() +
  xlab("") + ylab("")


base_world <- p + geom_polygon(data = mapSubs, aes(x=long, y=lat, group=group),
                               colour="white", fill="gray70") + theme_bw()


## correcting the study that is from Newfoundland but has
## a positive longitude
data_descr$Longitude[data_descr$Species == 'Rangifer tarandus caribou'] <- -57

# change the order of the levels in the taxon, so that the fish is of orange colour
data_descr$Taxon <- factor(data_descr$Taxon, levels = c('Bird', 'Fish', 'Mammal','Reptile'))

map_tax <-
  base_world +
  geom_point(data = data_descr,
             aes(x = Longitude, y = Latitude, colour = Taxon, fill = Taxon),  size = 2,
             pch=21,  alpha=I(0.7), stroke = 1.2) +
  scale_color_brewer(palette = 'Dark2') +
  scale_fill_brewer(palette = 'Dark2') +
  theme(legend.position = c(0.09,0.16), legend.direction = 'vertical',
        legend.background = element_rect(color = 'white', fill="white", linewidth = 1, linetype="solid"),
        legend.text = element_text(size = rel(1.4)), legend.title =
          element_text(size = rel(1.4)), legend.key =
          element_blank()) +
  guides(color = guide_legend('Taxon', override.aes = list(size=5)),
         fill = guide_legend('Taxon')) +
  xlab('Longitude') + ylab('Latitude')

## map without the inset on number of studies yet
pdf('./output_all/Fig_Studies_onMap.pdf', width= 10, height = 6)
print(map_tax)
dev.off()


vpTax_perTraitCat <- grid::viewport(width = 0.32, height = 0.25,
                                    x = 0.675, y = 0.725, just = c('left', 'bottom'))

# final map with the inset on the number of studies per taxon within each
# trait category
pdf('./plots_ms/Fig2_Studies_onMap_PerTraitCat.pdf', width= 10, height = 6)
print(map_tax)
print(plot_TaxByTraitCat, vp = vpTax_perTraitCat)
dev.off()



# 3. Study duration histogram: SI figs ---------------------------------------


## histogram of the studies duration
sum_dur <- data_descr %>%
  dplyr::group_by(., Trait_Categ) %>%
  dplyr::summarise(., Median = median(Nyears)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Lab = c('a','b'), x = 5, y = 37)

pdf('./plots_ms/FigS1_StudyDuration.pdf')
ggplot(data_descr, aes(x = Nyears)) + geom_histogram() +
  facet_wrap(vars(Trait_Categ), nrow = 2) +
  geom_vline(data = sum_dur, mapping = aes(xintercept = Median), linetype = 2, color = 'red') +
  geom_text(data = sum_dur, aes(x = Median + 2, y = 20, label = Median), color = 'red') +
  theme_bw() + theme(strip.background = element_blank(),
                     strip.text = element_text(size = 12),
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 10, color = 'black')) +
  xlab('Study duration (years)') + ylab('Count') +
  geom_text(data = sum_dur, aes(x = x, y = y, label = Lab),
            fontface = 'bold', size = 5)
dev.off()

# For PPT, colours changed and enlarged font
pdf('./output_all/Descriptive_StudyDurationPPT.pdf', width = 12)
ggplot(data_descr, aes(x = Nyears)) + geom_histogram() +
  facet_wrap(vars(Trait_Categ), nrow = 1) +
  geom_vline(data = sum_dur, mapping = aes(xintercept = Median), linetype = 2, lwd = 3, color = 'darkorange2') +
  geom_text(data = sum_dur, aes(x = Median + 5, y = 20, label = Median), color = 'darkorange2', size = 10) +
  theme_bw() + theme(strip.background = element_blank(),
                     strip.text = element_text(size = 23),
                     axis.title = element_text(size = 23),
                     axis.text = element_text(size = 23, color = 'black')) +
  xlab('Study duration (years)') + ylab('Count') +
  geom_text(data = sum_dur, aes(x = x, y = y, label = Lab),
            fontface = 'bold', size = 10)
dev.off()


# number of studies per species fro the ppt
table(data_descr$Species)
tab_Sp <- as.data.frame(table(data_descr$Species))
names(tab_Sp) <- c('Species', 'Count')
tab_Sp <- tab_Sp[order(tab_Sp$Count), ]

#  - at least three studies per species
tab_Sp_sel <- tab_Sp[tab_Sp$Count > 2, ]
tab_Sp_sel$Species <- factor(tab_Sp_sel$Species, levels = tab_Sp_sel$Species)
pdf('./output_all/Descriptive_Barplot_bySpeciesPPT.pdf', height = 7)
ggplot(tab_Sp_sel, aes(x = Species, y = Count)) +
  geom_bar(stat= 'identity', position = position_dodge(),
           fill = 'darkorange2') +
  theme_bw() + theme(axis.text = element_text(size = 18,
                                              colour = 'black'),
                     axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 1),
                     axis.title = element_text(size = 23),
                     legend.position = 'bottom',
                     legend.text = element_text(size= 21),
                     legend.title = element_text(size = 23),
                     plot.margin = unit(c(2, 1, 10, 5), 'pt')) +
  guides(fill=guide_legend(ncol=2))
dev.off()


# and keeping only species for which we have at least five studies
tab_Sp_sel5 <- tab_Sp_sel[tab_Sp_sel$Count > 4, ]
pdf('./output_all/Descriptive_Barplot_bySpecies_min5Studies_PPT.pdf', height = 7)
ggplot(tab_Sp_sel5, aes(x = Species, y = Count)) +
  geom_bar(stat= 'identity', position = position_dodge(),
           fill = 'darkorange2') +
  theme_bw() + theme(axis.text = element_text(size = 18,
                                              colour = 'black'),
                     axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 1),
                     axis.title = element_text(size = 23),
                     legend.position = 'bottom',
                     legend.text = element_text(size= 21),
                     legend.title = element_text(size = 23),
                     plot.margin = unit(c(2, 1, 10, 5), 'pt')) +
  guides(fill=guide_legend(ncol=2))
dev.off()


# 4. Histograms of popsize per study --------------------------------------

## plot histograms of the original pop size data in one file, with a facet per study
## abbreviations for study and species per facet
dat_abbrev <- temp_SEM %>%
  dplyr::group_by(., ID) %>%
  dplyr::distinct(., ID, Study_Authors, Species, .keep_all = TRUE)


abbrev_d <- dat_abbrev %>%
  dplyr::group_by(., ID) %>%
  tidyr::nest(.data =., cols = ID) %>%
  dplyr::mutate(
    fLet = purrr::map_chr(dat_abbrev$Species, ~substr(., start = 1, stop = 1)),
    sName = purrr::map_chr(dat_abbrev$Species, ~substr(strsplit(as.character(.), split = ' ')[[1]][2], 1, 4)),
    label = paste(Study_Authors, paste(fLet, sName, sep = '.'), sep = ': ')
  ) %>%
  tidyr::unnest(cols = c(cols))


hist_pop <- ggplot(temp_SEM, aes(Pop_mean)) + geom_histogram(bins = 20) +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_text(size = 3, data = abbrev_d,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = colours()[505])) +
  theme(legend.position = 'none')


print(hist_pop)


ggforce::n_pages(hist_pop)  ## calculates number of pages in a plot to know via how many to iterate
pdf('./output_all/Pop_hist.pdf')
for(i in 1:11){
  p <- ggplot(temp_SEM, aes(Pop_mean)) +
    geom_histogram(bins = 20) +
    theme_bw() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    theme(strip.background = element_rect(colour = 'black', fill = 'white')) +
    geom_text(size = 2, data = abbrev_d,
              aes(y = Inf, x = -Inf, label = label, vjust = 1, hjust = 0,
                  colour = colours()[505]), fontface ='bold') + theme(legend.position = 'none')
  print(p)
}
dev.off()


# 5. Prep supplementary data file -----------------------------------------

# adding the traits also

## data on traits
traits <- read.csv('./data/speciesTraits.csv')

## some preparation of the trait data prior to use
## adding two alternative diet classif + fixing the migratory mode
traits$Diet.5Cat_Elton <- trimws(traits$Diet.5Cat_Elton, 'both')
traits_proc <- traits %>%
  dplyr::mutate(., DietPS = dplyr::case_when(
    Diet.5Cat_Elton %in% c('Herbivore', 'PlantSeed') ~ 'primary consumer',
    Diet.5Cat_Elton %in% c('Invertebrate', 'Insectivore') ~ 'secondary consumer, invert',
    Diet.5Cat_Elton %in% c('Carnivore', 'VertFish', 'VertFishScav') ~ 'secondary consumer, vert',
    Diet.5Cat_Elton %in% c('Omnivore', 'VertInvertEggs', 'InvertFish') ~ 'secondary consumer, omnivore'),
    Diet_HCO = dplyr::case_when(
      Diet.5Cat_Elton %in% c('Herbivore', 'PlantSeed') ~ 'herbivore',
      Diet.5Cat_Elton %in% c('Omnivore') ~ 'omnivore',
      Diet.5Cat_Elton %in% c('Invertebrate', 'Insectivore', 'Carnivore',
                             'VertFish', 'VertFishScav', 'VertInvertEggs', 'InvertFish') ~ 'carnivore'),
    Migrat = dplyr::case_when(
      Migratory.mode_Sibly %in% c('nonmigrant', 'resident', 'nonmigrant???') ~ 'resident',
      Migratory.mode_Sibly %in% c('migrant') ~ 'migrant',
      TRUE ~ 'unknown'))


dat_traits <- merge(data_descr, subset(traits_proc, select = c(Species, GenLength_y_IUCN, Migrat, Diet_HCO)),
                    by = 'Species', all.x = TRUE)
dat_suppl <- dat_traits %>%
  dplyr::select(Species, ID, Study_Authors, Journal, Year_pub, Title, Taxon, Location, Longitude, Latitude,
         Country, Trait, Trait_Categ, Unit_trait, Nyears, GenLength_y_IUCN, Migrat, Diet_HCO) %>%
  rename(Duration = Nyears, Diet = Diet_HCO, GenTime_y_IUCN = GenLength_y_IUCN)

dat_suppl <- dat_suppl[order(dat_suppl$ID), ]


# reclassify into a broader trait type that encompassess similar
# traits (e.g. arrival of females and arrivla of males)
dat_s_phen_type <- dat_suppl %>%
  mutate(TraitType = as.factor(case_when(Trait %in% c('ArrivalDateFemales',
                                                      'ArrivalDateMales',
                                                      'FemaleArrivalDate',
                                                      'MaleArrivalDate',
                                                      'SettlementDateYearOld') ~ 'ArrivalDate',
                                         Trait %in% c('BreedingDate',
                                                      'LayingDateAllBroods',
                                                      'MeanBreedingDate',
                                                      'LayingDate',
                                                      'NestInitiationDate',
                                                      'NestDate') ~ 'OnsetBreeding',
                                         Trait %in% c('StartOfLaying') ~ 'FirstLayDate',
                                         Trait %in% c('AntlerCastDate', 'RutEndDate',
                                                      'OestrusDate', 'RutStartDate') ~
                                           'RutDate',
                                         Trait %in% c('BirthDate') ~ 'ParturitionDate',
                                         .default = as.character(Trait))))
test <- droplevels(subset(dat_s_phen_type, Trait_Categ == 'Phenological'))
table(test$TraitType)
levels(dat_s_phen_type$TraitType)

table(test$TraitType)/sum(table(test$TraitType))

# ArrivalDate   EmergenceDate    FirstLayDate   Fledging_Date
# 0.06185567      0.02061856      0.07216495      0.01030928
# HatchingDate   OnsetBreeding ParturitionDate         RutDate
# 0.01030928      0.74226804      0.04123711      0.04123711


save_xlsx(table = dat_s_phen_type, table_name = './tables/SupplementaryData_AllStudies')

