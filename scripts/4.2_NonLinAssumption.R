## checking nonlinearity assumption underlying different relations

library(spaMM)
library(tidyverse)

# I. Nonlinearity in Temp-trait relation -------------------------------
## read in the raw data for temperature (for each year climate values, traits, dem. rates etc.)
temp_SEM <- readRDS(file = './output_forSEM_temp/all_SEM.RDS')
head(temp_SEM)


## before fitting the model have to perform all the preparations and standardizations as done for SEMs
## add the years where they are  omitted as rows
allYrs_T <- do.call('rbind', lapply(X = unique(temp_SEM$ID), FUN = function(x){
  subs <- droplevels(temp_SEM[temp_SEM$ID == x, ])
  full_NA <- data.frame(Year = seq(min(subs$Year), max(subs$Year), by = 1),
                        ID = x)
}))

consec_yrs_T <- merge(allYrs_T, temp_SEM, by = c('ID','Year'), all= T)


# 1. Phenological studies only --------------------------------------------

phen_temp <- subset(consec_yrs_T, Trait_Categ == 'Phenological')

## calculate GR
temp_GR_phen <- phen_temp %>%
  dplyr::mutate(., Pop_mean_lag = c(Pop_mean[-1], NA)) %>%
  dplyr::mutate(., GR = log(Pop_mean_lag / Pop_mean)) %>%
  dplyr::filter(., !is.na(GR) & !is.na(Trait_mean) &
                  !is.na(Demog_rate_mean) & !is.na(Pop_mean))


temp_GRRes_phen <- split(temp_GR_phen, temp_GR_phen$ID) %>%
  purrr::map(., ~lm(Clim ~ Year, data = .)) %>%
  purrr::map2(.x = ., .y = split(temp_GR_phen, f= temp_GR_phen$ID),
              .f = ~broom::augment_columns(x = .x, data = .y)) %>%
  dplyr::bind_rows()

quantile(temp_GR_phen$GR, probs = c(0.01, 0.05, 0.2, 0.5, 0.75, 0.95, 0.99))

temp_std_phen <- temp_GRRes_phen %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(Trait_SE = Trait_SE / sd(Trait_mean, na.rm = T),
                Demog_rate_SE = Demog_rate_SE /sd(Demog_rate_mean, na.rm = T),
                det_Clim = as.numeric(scale(`.resid`)),
                Trait_mean = scale(Trait_mean),
                Demog_rate_mean = scale(Demog_rate_mean),
                Pop_mean = scale(Pop_mean),
                GR = scale(GR)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Trait_mean2 = Trait_mean^2,
                det_Clim2 = det_Clim^2)

temp_std_phen$GR <- as.numeric(temp_std_phen$GR[,1])

glmm_Traitphen_quad <- fitme(Trait_mean ~ det_Clim + det_Clim2 +
                          (det_Clim|ID) +
                          AR1(1|Year),
                        data = temp_std_phen, method = 'REML')
summary(glmm_Traitphen_quad)

## fitting the model without the det_Clim2 as fixed effect, to assess the signif. of this quadratic effect ACROSS the studies
glmm_Traitphen_lin_REML <- fitme(Trait_mean ~ det_Clim  +
                                 (det_Clim|ID) +
                                 AR1(1|Year),
                               data = temp_std_phen, method = 'REML')
summary(glmm_Traitphen_lin_REML)

# and ML
glmm_Traitphen_quad_ML <- fitme(Trait_mean ~ det_Clim + det_Clim2 +
                               (det_Clim|ID) +
                               AR1(1|Year),
                             data = temp_std_phen, method = 'ML')
summary(glmm_Traitphen_quad_ML)
glmm_Traitphen_lin_ML <- fitme(Trait_mean ~ det_Clim  +
                               (det_Clim|ID) +
                               AR1(1|Year),
                             data = temp_std_phen, method = 'ML')
summary(glmm_Traitphen_lin_ML)
anova(glmm_Traitphen_lin_ML, glmm_Traitphen_quad_ML, boot.repl = 0)
# Chi2:= -0.1764523  df = 1 and p = 1
AIC(glmm_Traitphen_quad_ML)
AIC(glmm_Traitphen_lin_ML)


# and if we looked at the models with random intercepts only (no random-slope model)
glmm_Traitphen_quad_ML_RI <- fitme(Trait_mean ~ det_Clim + det_Clim2 +
                                  (1|ID) +
                                  AR1(1|Year),
                                data = temp_std_phen, method = 'ML')
summary(glmm_Traitphen_quad_ML_RI)
glmm_Traitphen_lin_ML_RI <- fitme(Trait_mean ~ det_Clim  +
                                 (1|ID) +
                                 AR1(1|Year),
                               data = temp_std_phen, method = 'ML')
summary(glmm_Traitphen_lin_ML_RI)
anova(glmm_Traitphen_lin_ML_RI, glmm_Traitphen_quad_ML_RI, boot.repl = 0)
## Chi2 = 4.761957  df = 1 p = 0.02909545
AIC(glmm_Traitphen_lin_ML_RI)
AIC(glmm_Traitphen_quad_ML_RI)

# 2. Morphological studies only --------------------------------------------

morph_temp <- subset(consec_yrs_T, Trait_Categ == 'Morphological')
nrow(morph_temp)
length(unique(morph_temp$ID))

## calculate GR
temp_GR_morph <- morph_temp %>%
  dplyr::mutate(., Pop_mean_lag = c(Pop_mean[-1], NA)) %>%
  dplyr::mutate(., GR = log(Pop_mean_lag / Pop_mean)) %>%
  dplyr::filter(., !is.na(GR) & !is.na(Trait_mean) &
                  !is.na(Demog_rate_mean) & !is.na(Pop_mean))


temp_GRRes_morph <- split(temp_GR_morph, temp_GR_morph$ID) %>%
  purrr::map(., ~lm(Clim ~ Year, data = .)) %>%
  purrr::map2(.x = ., .y = split(temp_GR_morph, f= temp_GR_morph$ID),
              .f = ~broom::augment_columns(x = .x, data = .y)) %>%
  dplyr::bind_rows()

quantile(temp_GRRes_morph$GR, probs = c(0.01, 0.05, 0.2, 0.5, 0.75, 0.95, 0.99))

temp_std_morph <- temp_GRRes_morph %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(Trait_SE = Trait_SE / sd(Trait_mean, na.rm = T),
                Demog_rate_SE = Demog_rate_SE /sd(Demog_rate_mean, na.rm = T),
                det_Clim = as.numeric(scale(`.resid`)),
                Trait_mean = scale(Trait_mean),
                Demog_rate_mean = scale(Demog_rate_mean),
                Pop_mean = scale(Pop_mean),
                GR = scale(GR)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Trait_mean2 = Trait_mean^2,
                det_Clim2 = det_Clim^2)

temp_std_morph$GR <- as.numeric(temp_std_morph$GR[,1])

glmm_Traitmorph_quad <- fitme(Trait_mean ~ det_Clim + det_Clim2 +
                               (det_Clim|ID) +
                               AR1(1|Year),
                             data = temp_std_morph, method = 'REML')
summary(glmm_Traitmorph_quad)

## fitting the model without the det_Clim2 as fixed effect, to assess how signif. is this quadratic effect ACROSS the studies
glmm_Traitmorph_lin_REML <- fitme(Trait_mean ~ det_Clim  +
                                  (det_Clim|ID) +
                                  AR1(1|Year),
                                data = temp_std_morph, method = 'REML')
summary(glmm_Traitmorph_lin_REML)


# and with ML, for LRT
glmm_Traitmorph_quad_ML <- fitme(Trait_mean ~ det_Clim + det_Clim2 +
                                  (det_Clim|ID) +
                                  AR1(1|Year),
                                data = temp_std_morph, method = 'ML')
summary(glmm_Traitmorph_quad_ML)

glmm_Traitmorph_lin_ML <- fitme(Trait_mean ~ det_Clim  +
                                 (det_Clim|ID) +
                                 AR1(1|Year),
                               data = temp_std_morph, method = 'ML')
summary(glmm_Traitmorph_lin_ML)
anova(glmm_Traitmorph_lin_ML, glmm_Traitmorph_quad_ML, boot.repl = 0)
## Chi2 = 0.9215953, df = 1, p value = 0.3370564
AIC(glmm_Traitmorph_quad_ML)
AIC(glmm_Traitmorph_lin_ML)


# and if we looked at the models with mixed intercepts only (no random-slope model)
glmm_Traitmorph_quad_ML_RI <- fitme(Trait_mean ~ det_Clim + det_Clim2 +
                                     (1|ID) +
                                     AR1(1|Year),
                                   data = temp_std_morph, method = 'ML')
summary(glmm_Traitmorph_quad_ML_RI)
glmm_Traitmorph_lin_ML_RI <- fitme(Trait_mean ~ det_Clim  +
                                    (1|ID) +
                                    AR1(1|Year),
                                  data = temp_std_morph, method = 'ML')
summary(glmm_Traitmorph_lin_ML_RI)
anova(glmm_Traitmorph_lin_ML_RI, glmm_Traitmorph_quad_ML_RI, boot.repl = 0)
## signif: Chi2 =   4.861175,  df = 1 ,p = 0.02746763

AIC(glmm_Traitmorph_lin_ML_RI)
AIC(glmm_Traitmorph_quad_ML_RI)



# II. Nonlinearity in Precip-trait relation -------------------------------

## read in the raw data for precipitaion (for each year climate values, traits, dem. rates etc.)
precip_SEM <- readRDS(file = './output_forSEM_precip/all_SEM.RDS')
head(precip_SEM)

## before fitting the model have to perform all the preparations and standardizations as done for SEMs
## add the years where they are simply omitted as rows
allYrs_P <- do.call('rbind', lapply(X = unique(precip_SEM$ID), FUN = function(x){
  subs <- droplevels(precip_SEM[precip_SEM$ID == x, ])
  full_NA <- data.frame(Year = seq(min(subs$Year), max(subs$Year), by = 1),
                        ID = x)
}))

consec_yrs_P <- merge(allYrs_P, precip_SEM, by = c('ID','Year'), all= T)


# 1. Phenological studies only --------------------------------------------

phen_prec <- subset(consec_yrs_P, Trait_Categ == 'Phenological')

## calculate GR
prec_GR_phen <- phen_prec %>%
  dplyr::mutate(., Pop_mean_lag = c(Pop_mean[-1], NA)) %>%
  dplyr::mutate(., GR = log(Pop_mean_lag / Pop_mean)) %>%
  dplyr::filter(., !is.na(GR) & !is.na(Trait_mean) &
                  !is.na(Demog_rate_mean) & !is.na(Pop_mean))


prec_GRRes_phen <- split(prec_GR_phen, prec_GR_phen$ID) %>%
  purrr::map(., ~lm(Clim ~ Year, data = .)) %>%
  purrr::map2(.x = ., .y = split(prec_GR_phen, f= prec_GR_phen$ID),
              .f = ~broom::augment_columns(x = .x, data = .y)) %>%
  dplyr::bind_rows()

quantile(prec_GRRes_phen$GR, probs = c(0.01, 0.05, 0.2, 0.5, 0.75, 0.95, 0.99))

prec_std_phen <- prec_GRRes_phen %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(Trait_SE = Trait_SE / sd(Trait_mean, na.rm = T),
                Demog_rate_SE = Demog_rate_SE /sd(Demog_rate_mean, na.rm = T),
                det_Clim = as.numeric(scale(`.resid`)),
                Trait_mean = scale(Trait_mean),
                Demog_rate_mean = scale(Demog_rate_mean),
                Pop_mean = scale(Pop_mean),
                GR = scale(GR)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Trait_mean2 = Trait_mean^2,
                det_Clim2 = det_Clim^2)

prec_std_phen$GR <- as.numeric(prec_std_phen$GR[,1])

glmm_P_Traitphen_quad <- fitme(Trait_mean ~ det_Clim + det_Clim2 +
                               (det_Clim|ID) +
                               AR1(1|Year),
                             data = prec_std_phen, method = 'REML')
summary(glmm_P_Traitphen_quad) ## hmm precip really does not have an effect on phenology

## fitting the model without the det_Clim2 as fixed effect, to assess how signif. is this quadratic effect ACROSS the studies
glmm_P_Traitphen_lin_REML <- fitme(Trait_mean ~ det_Clim  +
                                   (det_Clim|ID) +
                                   AR1(1|Year),
                                 data = prec_std_phen, method = 'REML')
summary(glmm_P_Traitphen_lin_REML)


# and with ML, for LRT
glmm_P_Traitphen_quad_ML <- fitme(Trait_mean ~ det_Clim + det_Clim2 +
                                  (det_Clim|ID) +
                                  AR1(1|Year),
                                data = prec_std_phen, method = 'ML')
summary(glmm_P_Traitphen_quad_ML)
glmm_P_Traitphen_lin_ML <- fitme(Trait_mean ~ det_Clim  +
                                 (det_Clim|ID) +
                                 AR1(1|Year),
                               data = prec_std_phen, method = 'ML')
summary(glmm_P_Traitphen_lin_ML)
anova(glmm_P_Traitphen_lin_ML, glmm_P_Traitphen_quad_ML, boot.repl = 0)
## Chi2 = 0.5280198, df = 1, p value = 0.4674406
AIC(glmm_P_Traitphen_quad_ML)
AIC(glmm_P_Traitphen_lin_ML)



# and if we looked at the models with mixed intercepts only (no random- slope model)
glmm_P_Traitphen_quad_ML_RI <- fitme(Trait_mean ~ det_Clim + det_Clim2 +
                                     (1|ID) +
                                     AR1(1|Year),
                                   data = prec_std_phen, method = 'ML')
summary(glmm_P_Traitphen_quad_ML_RI)
glmm_P_Traitphen_lin_ML_RI <- fitme(Trait_mean ~ det_Clim  +
                                    (1|ID) +
                                    AR1(1|Year),
                                  data = prec_std_phen, method = 'ML')
summary(glmm_P_Traitphen_lin_ML_RI)
anova(glmm_P_Traitphen_lin_ML_RI, glmm_P_Traitphen_quad_ML_RI, boot.repl = 0)
## non-signif. Chi2 = 0.853069  df = 1, p =  0.3556856


# 2. Morphological studies only --------------------------------------------

morph_prec <- subset(consec_yrs_P, Trait_Categ == 'Morphological')
nrow(morph_prec)
length(unique(morph_prec$ID))

## calculate GR
prec_GR_morph <- morph_prec %>%
  dplyr::mutate(., Pop_mean_lag = c(Pop_mean[-1], NA)) %>%
  dplyr::mutate(., GR = log(Pop_mean_lag / Pop_mean)) %>%
  dplyr::filter(., !is.na(GR) & !is.na(Trait_mean) &
                  !is.na(Demog_rate_mean) & !is.na(Pop_mean))


prec_GRRes_morph <- split(prec_GR_morph, prec_GR_morph$ID) %>%
  purrr::map(., ~lm(Clim ~ Year, data = .)) %>%
  purrr::map2(.x = ., .y = split(prec_GR_morph, f= prec_GR_morph$ID),
              .f = ~broom::augment_columns(x = .x, data = .y)) %>%
  dplyr::bind_rows()

quantile(prec_GRRes_morph$GR, probs = c(0.01, 0.05, 0.2, 0.5, 0.75, 0.95, 0.99))

prec_std_morph <- prec_GRRes_morph %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(Trait_SE = Trait_SE / sd(Trait_mean, na.rm = T),
                Demog_rate_SE = Demog_rate_SE /sd(Demog_rate_mean, na.rm = T),
                det_Clim = as.numeric(scale(`.resid`)),
                Trait_mean = scale(Trait_mean),
                Demog_rate_mean = scale(Demog_rate_mean),
                Pop_mean = scale(Pop_mean),
                GR = scale(GR)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Trait_mean2 = Trait_mean^2,
                det_Clim2 = det_Clim^2)

prec_std_morph$GR <- as.numeric(prec_std_morph$GR[,1])

glmm_P_Traitmorph_quad <- fitme(Trait_mean ~ det_Clim + det_Clim2 +
                                (det_Clim|ID) +
                                AR1(1|Year),
                              data = prec_std_morph, method = 'REML')
summary(glmm_P_Traitmorph_quad) ## gmmmm no really an effect of temp on morphology... Much heterogeneity..

## fitting the model without the det_Clim2 as fixed effect, to assess how signif. is this quadratic effect ACROSS the studies
glmm_P_Traitmorph_lin_REML <- fitme(Trait_mean ~ det_Clim  +
                                    (det_Clim|ID) +
                                    AR1(1|Year),
                                  data = prec_std_morph, method = 'REML')
summary(glmm_P_Traitmorph_lin_REML)

## and with ML, for LRT
glmm_P_Traitmorph_quad_ML <- fitme(Trait_mean ~ det_Clim + det_Clim2 +
                                   (det_Clim|ID) +
                                   AR1(1|Year),
                                 data = prec_std_morph, method = 'ML')
summary(glmm_P_Traitmorph_quad_ML)

glmm_P_Traitmorph_lin_ML <- fitme(Trait_mean ~ det_Clim  +
                                  (det_Clim|ID) +
                                  AR1(1|Year),
                                data = prec_std_morph, method = 'ML')
summary(glmm_P_Traitmorph_lin_ML)
anova(glmm_P_Traitmorph_lin_ML, glmm_P_Traitmorph_quad_ML, boot.repl = 0)
## Chi2 = 0.02349689, df =1,  p value = 0.878172
AIC(glmm_P_Traitmorph_quad_ML)
AIC(glmm_P_Traitmorph_lin_ML)



# and if we looked at the models with mixed intercepts only (no random slope model)
glmm_P_Traitmorph_quad_ML_RI <- fitme(Trait_mean ~ det_Clim + det_Clim2 +
                                      (1|ID) +
                                      AR1(1|Year),
                                    data = prec_std_morph, method = 'ML')
summary(glmm_P_Traitmorph_quad_ML_RI)
glmm_P_Traitmorph_lin_ML_RI <- fitme(Trait_mean ~ det_Clim  +
                                     (1|ID) +
                                     AR1(1|Year),
                                   data = prec_std_morph, method = 'ML')
summary(glmm_P_Traitmorph_lin_ML_RI)
anova(glmm_P_Traitmorph_lin_ML_RI, glmm_P_Traitmorph_quad_ML_RI, boot.repl = 0)
## non-signif: Chi2 =  1.906926, df =  1, p = 0.1673051
AIC(glmm_P_Traitmorph_lin_ML_RI)
AIC(glmm_P_Traitmorph_quad_ML_RI)


# III. Nonlinearity in temperature-driven GR-trait  relations ---------------------------

# 1. Phenology + Temp------------------------------------------------------------

## plotting all the data together
temp_std_phen$ID <- as.factor(temp_std_phen$ID)

pdf('./output_all/shape_phenTraitGR_raw.pdf')
for(i in seq(0, 100, by = 10)){

  print(ggplot(subset(temp_std_phen, ID %in% unique(temp_std_phen$ID)[i:(i+9)]),
               aes(Trait_mean, GR, colour=ID)) +
          ylim(-2, 2) +
          # geom_point() +
          geom_smooth(size=1.5, linetype="dotted", se=FALSE,
                      method="gam",
                      formula = y ~ s(x, bs = "cs", k = 4)) +
          theme_classic() + theme(legend.position = 'bottom'))
}
dev.off()


# a test with trait as predictor of GR  + random trait slope
test_phenTGR_TraitOnly <- test_nonlin(data = temp_std_phen,
                                  formula_full  = 'GR ~ Trait_mean + Trait_mean2 +
                                  (Trait_mean|ID)',
                                  formula_null  = 'GR ~ Trait_mean +
                                  (Trait_mean|ID)')
test_phenTGR_TraitOnly
# Chi2 =  0.2555225, df =  1, p =  0.6132131

# and if using rand intercept only
test_phenTGR_TraitOnly_int <- test_nonlin(data = temp_std_phen,
                                     formula_full  = 'GR ~ Trait_mean + Trait_mean2 +
                                  (1|ID)',
                                     formula_null  = 'GR ~ Trait_mean +
                                  (1|ID)')

test_phenTGR_TraitOnly_int  ## also non-signif.
## Chi2 = 0.5251621, df =  1, p = 0.468648



# 2. Morphology + Temp------------------------------------------------------------

## plotting all the data together
temp_std_morph$ID <- as.factor(temp_std_morph$ID)
pdf('./output_all/shape_morphTraitGR_raw.pdf')
for(i in seq(0, 120, by = 10)){

  print(ggplot(subset(temp_std_morph, ID %in% unique(temp_std_morph$ID)[i:(i+9)]),
               aes(Trait_mean, GR, colour=ID)) +
          ylim(-2, 2) +
          # geom_point() +
          geom_smooth(size=1.5, linetype="dotted", se=FALSE,
                      method="gam",
                      formula = y ~ s(x, bs = "cs", k = 4)) +
          theme_classic() + theme(legend.position = 'bottom'))
}
dev.off()


# test for the non-linear effect of traits on GR (not including climate)
## using random trait slope
test_morphTGR_TraitOnly <- test_nonlin(data = temp_std_morph,
                                     formula_full  = 'GR ~ Trait_mean + Trait_mean2 +
                                  (Trait_mean|ID)',
                                     formula_null  = 'GR ~ Trait_mean +
                                  (Trait_mean|ID)')
test_morphTGR_TraitOnly
# Chi2 = 1.937405 , df = 1, p value = 0.1639505

# a simpler test, with trait only as predictor  + random intercept only
test_morphTGR_TraitOnly_int <- test_nonlin(data = temp_std_morph,
                                      formula_full  = 'GR ~ Trait_mean + Trait_mean2 +
                                  (1|ID)',
                                      formula_null  = 'GR ~ Trait_mean +
                                  (1|ID)')
test_morphTGR_TraitOnly_int
# Chi2 = 1.869467, df=  1, p = 0.1715362



# 3. Phenology + Precip------------------------------------------------------------
# a simpler test, with trait only as predictor  + random trait slope
test_phenPGR_TraitOnly <- test_nonlin(data = prec_std_phen,
                                      formula_full  = 'GR ~ Trait_mean + Trait_mean2 +
                                  (Trait_mean|ID)',
                                      formula_null  = 'GR ~ Trait_mean +
                                  (Trait_mean|ID)')
test_phenPGR_TraitOnly
# chi2_LR = p_v 0.2842072, df = 1,  p_value = 0.5939568

# a simpler test, with trait only as predictor  + random intercept only
test_phenPGR_TraitOnly_int <- test_nonlin(data = prec_std_phen,
                                          formula_full  = 'GR ~ Trait_mean + Trait_mean2 +
                                  (1|ID)',
                                          formula_null  = 'GR ~ Trait_mean +
                                  (1|ID)')
test_phenPGR_TraitOnly_int
# Chi2 = 0.5620872 , df = 1, p= 0.4534205


# 4. Morphology + Precip ------------------------------------------------------------

# a simpler test, with trait only as predictor  + random trait slope
test_morphPGR_TraitOnly <- test_nonlin(data = prec_std_morph,
                                      formula_full  = 'GR ~ Trait_mean + Trait_mean2 +
                                  (Trait_mean|ID)',
                                      formula_null  = 'GR ~ Trait_mean +
                                  (Trait_mean|ID)')
test_morphPGR_TraitOnly
# chi2_LR = 1.937405, df =   1, p = 0.1639505

# a simpler test, with trait only as predictor  + random intercept only
test_morphPGR_TraitOnly_int <- test_nonlin(data = prec_std_morph,
                                          formula_full  = 'GR ~ Trait_mean + Trait_mean2 +
                                  (1|ID)',
                                          formula_null  = 'GR ~ Trait_mean +
                                  (1|ID)')
test_morphPGR_TraitOnly_int
# Chi2 = 1.869467, df =   1, p =  0.1715362
