
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sTraitChange_Analyses

<!-- badges: start -->
<!-- badges: end -->

The goal of sTraitChange_Analyses project is to provide the full
workflow as well as the data needed to reproduce the analyses for the
manuscript ‘Changes in phenology mediate vertebrate population responses
to climate globally’ by Radchuk et al. (submitted). The analyses in this
repository rely heavily on the functions available from
[sTraitChange](https://github.com/radchukv/sTraitChange) package.

To be able to run these analyses you will have to install
[sTraitChange](https://github.com/radchukv/sTraitChange) package.

## Structure

This project contains the following folders:  
- data: datasets used in the analyses;  
- output_SEM_all: general outputs from the analyses based on structural
equation models (SEMs);  
- output_forSEM_precip and output_forSEM_temp: the prepared datasets
ready for running SEMs using, respectively precipitation and temperature
as climate variables;  
- output_fSEM_precip and output_fSEM_temp: the output from the SEM
analyses, that are used later on in meta-analyses, obtained, using
respectively precipitation and temperature as climate variables;  
- output_fSEM_noDD_precip and output_fSEM_noDD_temp: the output from the
SEM analyses conducted without population size included in SEM
(sensitivity analyses: for details, see the manuscript). Some of these
generated results are intermediate and are later used in meta-analyses.
These results are obtained, using respectively precipitation and
temperature as climate variable;  
- output_all_noDD: general output from the analyses conducted when
population size is not included in the SEM;  
- output_all: general output obtained from the complete workflow;  
- plots_ms: plots used in the manuscript;  
- tables: tables produced during the analyses;  
- scripts: all scripts used in the workflow (more details below).

## Scripts

The scripts are numbered in the logical order in which they are used in
the analyses. The scripts pertaining to sliding window analyses will not
run out of the box because relevant climate datasets could not be shared
in the repo (due to its prohibitively big size). In order to run those
scripts the user will have to download the climate data from the
respective sources mentioned in Supplementary materials (Table S6) and
place them in the respective directory on the path (./data).

Brief description of what each script does:  
- **1.1_DataPrep_Clim.R**: prepares climate data for sliding window
analyses using `climwin` package.  
- **1.2_DataPrep_Biol.R**: prepares biological data for sliding window
analyses using `climwin` package.  
- **2_climwin.R**: performs sliding window analyses (will not run
directly, see above!);  
- **3.0_PrepData_forSEM.R**: prepares the data for SEM analyses;  
- **3.1_SEM.R**: performs SEM analyses;  
- **4.0_Descriptive.R**: provides description of the total dataset;  
- **4.1_ClimwinWindows.R**: describes results from sliding window
analyses;  
- **4.2_NonLinAssumption.R**: assess non-linearity assumptions;  
- **5.1_Phylo_prep.R**: prepares a mega-tree for performing
meta-analyses that account for phylogenetic relatedness;  
- **5.2.1_Metaanalysis_Temp.R**: meta-analyses focusing on temperature
as climate variable;  
- **5.2.2_Metaanalysis_Precip.R**: meta-analyses focusing on
precipitation as climate variable;  
- **5.3.1_Metaanalysis_Temp_noDD.R** &
**5.3.2_Metaanalysis_Precip_noDD.R**: meta-analyses focusing,
respectively on temperature and precipitation as climate variable and
performed without the population size in SEM (sensitivity analyses);  
- **5.4_Metaanalysis_AcrossPhylogenies.R**: a script running
meta-analyses across 100 posterior phylogenies;  
- **5.5_Metaanalysis_HeterogeneityPaths.R**: a script explaining the
heterogeneity in the SEM paths with hypothesized predictors;  
- **5.6_Metaanalysis_Heterogeneity_FineCatNames.R**: assessing the
impact of the exact measured trait on heterogeneity in the phenotypic
responses to climate;  
- **6_pDeltaAIC.R**: exploration of the potential impact that pDeltaAIC
has on the results;  
- **7_Figures_MS.R**: Figures used in the manuscript.

## Support

For any further information, please contact <radchuk@izw-berlin.de>.
