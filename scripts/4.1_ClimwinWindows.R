## This script describes results from sliding window analyses,
## e.g. window durations and the p deltaAIC values acrorss studies

library(ggplot2)
library(patchwork)
library(tidyr)

# I. Data prep ------------------------------------------------------------
## read in the data for temperature
Cstat_Aut <- readRDS(file = './output_fSEM_temp/Cstat_allMods_Temp_Weights_DD_Autocorr.RDS')


#make open and close columns
Cstat_Aut$WinDur <- Cstat_Aut$WinDur + 1 ## because dfor cases where the window is open and closed in the same week, the duration of the window is 1 week (not 0)
Cstat_Aut$refdate <- as.Date(with(Cstat_Aut, paste(Ref.day, Ref.month, '2020', sep="-")), "%d-%m-%Y")
Cstat_Aut$WinCloseDate <- Cstat_Aut$refdat - as.difftime(Cstat_Aut$WindowClose, unit = 'weeks')
Cstat_Aut$WinOpenDate <- Cstat_Aut$WinCloseDate - as.difftime(Cstat_Aut$WinDur, unit = 'weeks')
Cstat_Aut$Climate <- 'Temperature'


## read in the data on precip
Cstat_Aut_precip <- readRDS(file = './output_fSEM_precip/Cstat_allMods_Precip_Weights_DD_Autocorr.RDS')
Cstat_Aut_precip$WinDur <- Cstat_Aut_precip$WinDur + 1
Cstat_Aut_precip$refdate <- as.Date(with(Cstat_Aut_precip, paste(Ref.day, Ref.month, '2020', sep="-")), "%d-%m-%Y")
Cstat_Aut_precip$WinCloseDate <- Cstat_Aut_precip$refdat - as.difftime(Cstat_Aut_precip$WindowClose, unit = 'weeks')
Cstat_Aut_precip$WinOpenDate <- Cstat_Aut_precip$WinCloseDate - as.difftime(Cstat_Aut_precip$WinDur, unit = 'weeks')
Cstat_Aut_precip$Climate <- 'Precipitation'


## checking how it looks for the complete dataset
nrow(Cstat_Aut_precip)   ## 211

## binding the data from both datasets
bothClimate <- rbind(Cstat_Aut, Cstat_Aut_precip)


# II. DeltaAIC,  pDeltaAIC & WinDur exploratory plots ---------------------


## some summary stats for the supplementary material
bothClim <- bothClimate %>%
  dplyr::filter(., WinDur <= 51) %>%
  dplyr::mutate(., dataset = 'Reduced')

ecdf(bothClim$deltaAIC[bothClim$Climate == 'Temperature'])(-5)  ##0.82
ecdf(bothClim$deltaAIC[bothClim$Climate == 'Precipitation'])(-7)  ## 0.66

ecdf(bothClim$Pvalue[bothClim$Climate == 'Temperature'])(0.05)  ## 0.21
ecdf(bothClim$Pvalue[bothClim$Climate == 'Precipitation'])(0.05)  ## 0.07

ecdf(bothClim$Pvalue[bothClim$Climate == 'Temperature'])(0.1)  ##  0.29
ecdf(bothClim$Pvalue[bothClim$Climate == 'Precipitation'])(0.1)  ##  0.15


## plot of p value vs climate variable and trait category
Clim_summ <- bothClim %>%
  dplyr::group_by(Climate, Trait_Categ) %>%
  dplyr::summarize(Median = median(Pvalue)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(label = letters[1:4])

pl_Pval_byTrait_Cl <- ggplot(bothClim, aes(x = Pvalue)) + geom_histogram() +
  facet_grid(Climate ~Trait_Categ) +
  geom_vline(data = Clim_summ, aes(xintercept = Median),
             col = 'red', lwd = 1.3) +
  geom_text(data = Clim_summ,
            aes(x = Median + 0.06, y = 20,
                label = round(Median, 2)),
            col = 'red', size = 5) +
  xlab(bquote(P[Delta*AICc])) + ylab('Frequency') +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(color = 'black',
                                 size = 10)) +
  geom_text(data = Clim_summ, aes(x = 0, y = 35, label = label),
            fontface = 'bold', size = 6)

pdf('./plots_ms/FigS2_PDeltaAIC_byTrait&Clim.pdf',
    height = 10, width = 10)
print(pl_Pval_byTrait_Cl)
dev.off()

Clim_summ_NumSt <- bothClim %>%
  dplyr::group_by(Climate, Trait_Categ) %>%
  dplyr::mutate(PvLess =
                  dplyr::case_when(Pvalue < 0.05 ~ 1,
                                   TRUE ~ 0)) %>%
  dplyr::mutate(AllSt = sum(PvLess))  %>%
  dplyr::add_count() %>%
  dplyr::mutate(Prop = AllSt / n) %>%
  dplyr::distinct(., Climate, Trait_Categ, .keep_all = TRUE)


## and now a histogram plot of deltaAICc
sum_deltaAIC <- bothClim %>%
  dplyr::group_by(Climate, Trait_Categ) %>%
  dplyr::summarize(Median = median(deltaAIC)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(label = letters[1:4],
                x = c(-10.5, -14.5, -10.5, -18),
                xlable = c(-30, -105, -30, -105))

pl_deltaAIC_byTrait_Cl <- ggplot(bothClim, aes(x = deltaAIC)) + geom_histogram() +
  facet_grid(Climate ~Trait_Categ, scales = 'free_x') +
  geom_vline(data = sum_deltaAIC, aes(xintercept = Median),
             col = 'red', lwd = 1.3) +
  geom_text(data = sum_deltaAIC,
            aes(x = x, y = 20, label = round(Median, 1)),
            col = 'red', size = 5) +
  xlab(bquote(Delta*AICc)) + ylab('Frequency') +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(color = 'black',
                                 size = 10)) +
  geom_text(data = sum_deltaAIC, aes(x = xlable, y = 41, label = label),
            fontface = 'bold', size = 6)

pdf('./output_all/DeltaAICc_byTrait&Clim.pdf',
    height = 10, width = 10)
print(pl_deltaAIC_byTrait_Cl)
dev.off()



## Supplementary Fig:  histograms of windur for precip and temp

nrow(bothClimate)  ##421
nrow(bothClim[bothClim$Climate == 'Temperature', ]) ## 202

## some stasts for the MS
ecdf(bothClim$WinDur[bothClim$Climate == 'Temperature'])(1)  ## 36%
ecdf(bothClim$WinDur[bothClim$Climate == 'Precipitation'])(1)  ## 31%


bothClimate$dataset <-  'Full'
bothCl <- dplyr::bind_rows(bothClimate, bothClim)

## summary stats
summar_winDur <- bothCl %>%
  dplyr::group_by(., Climate, dataset) %>%
  dplyr::summarise(., Median = median(WinDur), mean = mean(WinDur)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(xval = c(6, 4.5, 6, 5), x = rep(-2.5, 4),
                y = rep(100, 4), lab_txt = letters[1:4])

plot_winDur <- ggplot(bothCl, aes(x = WinDur)) + geom_histogram(binwidth = 1) +
  facet_grid(rows = vars(Climate), cols = vars(dataset), scales = 'free_x') +
  theme_bw() + theme(strip.background = element_blank(),
                     strip.text = element_text(size = 12)) +
  geom_vline(data = summar_winDur, aes(xintercept = Median),
             color = 'red', lwd = 1, linetype = 2) +
  geom_text(data = summar_winDur,
            aes(x = xval, y = 75, label = round(Median, 1)),
            color = 'red', size = 4.5) +
  labs(x = 'Window duration (weeks)', y = 'Frequency') +
  geom_text(data = summar_winDur,
            aes(x = x, y = y, label = lab_txt),
            fontface = 'bold', size = 5)

pdf('./plots_ms/FigS16_Histogram_WinDur_Full&Reduced.pdf')
print(plot_winDur)
dev.off()


# III. P value for the Cstat  ---------------------------------------------
## supplementary for the p value on Fisher's Cstat
summ_Cstat <- bothClim %>%
  dplyr::group_by(Climate, Trait_Categ) %>%
  dplyr::summarise(Median = median(P.Value)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(lab = letters[1:4], x = 0,
                y = c(26, 26, 19, 19))

pl_Cstat_pVal <- ggplot(bothClim, aes(x = P.Value)) +
  geom_histogram() +
  facet_grid(Climate ~ Trait_Categ, scales = 'free_y') +
  geom_vline(data = summ_Cstat, aes(xintercept = Median),
             col = 'red', lwd = 1.2) +
  geom_text(data = summ_Cstat,
            aes(x = Median + 0.07, y = 15,
                label = round(Median, 2)),
            col ='red', size = 5) + theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(color = 'black',
                                 size = 10)) +
  xlab('P value associated with Fishers C statistics') +
  ylab('Frequency') +
  geom_text(data = summ_Cstat, aes(x = x, y = y, label = lab),
            fontface = 'bold', size = 6)


pdf('./plots_ms/FigS6_Cstat_Pval_Clim&TraitCateg.pdf',
    height = 10, width = 10)
print(pl_Cstat_pVal)
dev.off()

# some stats for the SI
prop_More05 <- bothClim %>%
  dplyr::group_by(Trait_Categ) %>%
  dplyr::mutate(More05 = dplyr::case_when(
    P.Value >= 0.05 ~ 1,
    P.Value < 0.05 ~ 0)) %>%
  dplyr::summarize(Sum_more05 = sum(More05),
            Tot = dplyr::n()) %>%
  dplyr::mutate(prop = Sum_more05 / Tot)

# IV. Further description of window duration and Cstat ---------------------------------------
## some description of the window analyses
hist(Cstat_Aut$WinDur)
median(Cstat_Aut$WinDur)

## window duration
pdf('./output_all/WinDur_allMods_Temp_Weights_DD_Autocor.pdf')
hist(Cstat_Aut$WinDur, breaks = 30, col = 'grey',
     xlab = 'Climatic window duration', main = '')
dev.off()

median_winDur_Aut <- Cstat_Aut %>%
  dplyr::group_by(., Continent) %>%
  dplyr::summarize(Median = median(WinDur))


ggplot(Cstat_Aut, aes(WinDur)) + geom_histogram() +
  facet_grid(rows = vars(Continent), scales = 'free') +
  theme_bw() + theme(axis.title = element_text(size = rel(2)),
                     axis.text = element_text(size = rel(1.5)),
                     strip.text = element_text(size = rel(1.1))) +
  geom_vline(data = median_winDur_Aut, aes(xintercept = Median),
             col = 'red', lwd = 1.4) +
  geom_text(data = median_winDur_Aut, aes(label=round(Median, 2)),
            x = round(median_winDur_Aut$Median, 2) + 5,
            y = rep(0.8, nrow(median_winDur_Aut)),
            col = rep('red',  nrow(median_winDur_Aut)),
            size = rep(6, nrow(median_winDur_Aut)))


pdf('./output_all/DeltaAIC_allMods_Temp_Weights_DD_Autocor.pdf')
hist(Cstat_Aut$deltaAIC, breaks = 40, xlim = c(min(Cstat_Aut$deltaAIC, na.rm = TRUE) - 10, 0),
     col = 'grey', xlab = 'deltaAIC to the model with intercept only', main = '')
abline(v = median(Cstat_Aut$deltaAIC, na.rm = T), col = 'blue', lwd = 2)
dev.off()

pdf('./output_all/DeltaAIC_allMods_byTaxon_Temp_Weights_DD_Autocor.pdf')
ggplot(Cstat_Aut, aes(deltaAIC)) + geom_histogram() +
  facet_grid(rows = vars(Taxon), scales = 'free') +
  theme_bw() + theme(axis.title = element_text(size = rel(2)),
                     axis.text = element_text(size = rel(1.5)),
                     strip.text = element_text(size = rel(1.7)))
dev.off()
