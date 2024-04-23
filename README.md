# Endo_Stoch_Demo
This repository houses analysis and manuscript for the Stochastic  Demography of Fungal Endophyte - Grass Populations
### Repository Authors: 
Josh Fowler and Tom Miller
### Project Co-Authors: 
Jenn Rudgers, Ken Whitney, and Shaun Ziegler

## README file last updated: 
Apr. 23, 2024

### Project Overview:
Fungal endophytes are widespread symbionts of grasses that have been shown to provide a variety of context dependent benefits under environmental stress such as drought or salinity tolerance. Context-dependence may make interactions seem unpredictable as environmental conditions vary between years, but it provides a distinct mechanism by which symbionts may act as mutualists beyond influencing mean population growth rates, buffering their hosts from the detrimental effects of variance. We are quantifying the relative importance of mean and variance effects of microbial symbionts using long-term demographic data from experimental plots in Indiana at Lilly-Dickey Woods. 

The experiment, started in 2007, comprises 10-18 plots for each of 7 species of grass hosts. Half the plots were planted with 20 endophyte-infected plants, and half with 20 endophyte-free plants. The plots are censused annually for plant survival, plant size (measured as number of tillers) and reproduction (measured as number of flowering tillers and counts of seed and spikelets).

We estimate individual level vital rates from this experimental data. From these vital rate estimates, we build stochastic matrix projection models to make population projections which allow us to assess the contributions of endophyte partnership on both the mean and variance components of host fitness. 


Raw and cleaned data can be found at the following EDI data repository: [https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sev.343.2]{https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sev.343.2}

Certain figures in the publication have been compiled outside of R. 

## Repository Folder Description:
This repository is set up with three folders:
### Analyses 
holds scripts to cobble data and run analyses
File Name  | Description
------------- | -------------
2022-07-13_README_for_year_random_effects_in_MPM.txt | notes about meeting between Josh and Tom clarifying how we are lining up the year random effects of the vital rate models within the population model. summary is that surv/growth occur in year t->t1, and reproduction occurs in year t1, but this reproduction informs new recruits in the following year.
LTREB_endodemog_2021_Plant_locations_and_maps | script to compile maps of surviving plants in the plots during census based on XY coordinates
LTREB_endodemog_2022_Plant_locations_and_maps | script to compile maps of surviving plants in the plots during census based on XY coordinates
MPM_analysis.R  | Script to analyse matrix population model with and without endo effects on mean and variance. Includes plots of mean and variance effects. This script also performs stochastic simulations for simulation experiment determining the contrbution of mean and variance effects to stochastic growth rates.  Also includes a few plots at the end visualizing growth rates.
MPM_analysis_TOM.R  | Script to analyse matrix population model with and without endo effects on mean and variance, and to run simulations of increasing variance, including by sampling different subsets of observed years, or by sampling from the posterior distribution of year effects.
MPM_analysis_decomposition_quadXorigin.R | Script to run stochastic growth rate simulations for matrix population models, as above in MPM_analysis_TOM.R, but using the quadraticXorigin interaction vital rate models following reviewer suggestions.
MPM_analysis_quadXorigin.R | Script to analyse matrix population model as above in MPM_analysis.R, but using the quadraticXorigin interaction vital rate models following reviewer suggestions.
MPM_analysis_vital_rate_decomposition.R | Script to perform vital rate level decomposition analysis of endo effects on mean and variance. 
MPM_climate_analysis.R | Script to run climate-explicit version of matrix population models, including either 3 month or 12 month drought indices. This is used to examine endophyte-specific relationships between population growth and these climate drivers, and assess the relative sensitivity to climate. 
MPM_functions.R | Includes functions used in MPM_analysis scripts above.
endo_spp_grow_fert.stan | Stan model for growth and no. of inflorescences vital rate model for multiple species with zero-truncated negative binomial distribution with year and plot random effects (underestimates large sizes, and over-estimates variance and skew compared to Poisson Inverse Gaussian distribution).
endo_spp_grow_fert_PIG.stan | Stan model for growth and no. of inflorescences vital rate model for multiple species with Poisson Inverse Gaussian distribution, with year and plot random effects.
endo_spp_grow_fert_PIG-quadraticXorigin.stan | Stan model for growth and no. of inflorescences vital rate model for multiple species with Poisson Inverse Gaussian distribution, with year and plot random effects, incorporating quadratic size effects specific to recruit and original plants.
endo_spp_grow_fert_noplot.stan | Stan model for growth and no. of inflorescences vital rate model for multiple species with negative binomial distribution with only year random effects
endo_spp_s_to_s.stan | Stan model for seed to seedling transition  (germination) vital rate model for multiple species with Binomial distribution with year and plot random effects
endo_spp_spike_poisson.stan | Stan model for spikelets per inflorescence vital rate model for multiple species with Poisson distribution with year and plot random effects
endo_spp_spike_NB.stan | Stan model for spikelets per inflorescence vital rate model for multiple species with Negative Binomial distribution with year and plot random effects
endo_spp_spike_NB-quadraticXorigin.stan | Stan model for spikelets per inflorescence vital rate model for multiple species with Negative Binomial distribution with year and plot random effects for updated quadXorigina models.
endo_spp_surv_flw.stan | Stan model for survival and flowering status vital rate models for multiple species with Bernoulli distribution with year and plot random effects.
endo_spp_surv_flw-quadraticXorigin.stan | Stan model for survival and flowering status vital rate models for multiple species with Bernoulli distribution with year and plot random effects.
climate_endo_spp_grow_fert_PIG.stan | climate versions of vr models, incorporating 3 month or 12 month SPEI drought index as predictors.
climate_endo_spp_grow_fert_PIG_quadXorigin.stan | climate-explicit vital rate models including quadraditic size effects…
climate_endo_spp_s_to_s.stan | climate versions of vr models
climate_endo_spp_spike_nb.stan | climate versions of vr models
climate_endo_spp_spike_nb_quadXorigin.stan | climate-explicit vital rate models including quadraditic size effects…
climate_endo_spp_surv_flw.stan | climate versions of vr models
climate_endo_spp_surv_flw_quadXorigin.stan | climate-explicit vital rate models including quadraditic size effects…
climate_seedling_grow_PIG.stan | climate versions of vr models
climate_seedling_surv.stan | climate versions of vr models
endo_stoch_demo_matrix_compile.R | Script to compile average matrixes for project about life history and reproductive output with Robin Snyder.
endodemog_data_processing.R | Script that cleans legacy experimental data (2007-2018) and merges this legacy with ongoing field data (2019-2020). Data is stored in the Dropbox folder "EndodemogData". Legacy data manipulation involves pulling reproductive data out of the spreadsheets and merging this with size and survival data. For recent data, these are cleaned and also stored in the Dropbox sub-folder "Field Data", which include data as collected in the field, as well as cleaned data that correct some missing fields and assign tag not found status to recruits or previously found plants if possible. This script also downloads and merges weather station data with the demographic data.
seed_means.R | Script to run seed means model and visualize model diagnostics
seed_to_seedling.R | Script to run germination model and visualize model diagnostics
seedling_grow_PIG.stan | Stan model for first year seedling growth vital rate model for multiple species with Poisson Inverse Gaussian distribution with year and plot random effects. Our matrix model assumes a reproductive delay where first year plants of 1 tiller size do not reproduce.
seedling_surv.stan | Stan model for first year seedling survival vital rate model for multiple species with Bernoulli distribution with year and plot random effects. Our matrix model assumes a reproductive delay where first year plants of 1 tiller size do not reproduce.
stochastic_lambda_analysis.R | Test script to run stochastic population growth simulations and life table response experiment, replaced by MPM_analysis.R 
vital_rate_analysis.R | Script to run survival, growth, and reproductive vital rate models, and visualize model diagnostics
climate_explicit_vital_rate_analysis.R | Script to run survival, growth, and fertility vital rate models with climate drivers, and visualize model diagnostics
vital_rate_fit_figures.R | Script to plot model diagnostic plots for all species and all vital rates. Uses saved model outputs and saved y_rep outputs for each vital rate from vital_rate_analysis.R.
vital_rate_fit_figures_quadXorigin.R | Script to plot model diagnostic plots for all species and all vital rates for updated vital rate models following reviewer suggestions. Uses saved model outputs and saved y_rep outputs for each vital rate from vital_rate_analysis.R.
climate_vital_rate_fit_figures.R | Script to plot model diagnostic plots for all species and all vital rates for climate-explicit version of analyses. Uses saved model outputs and saved y_rep outputs for each vital rate from ckunate_explicit_vital_rate_analysis.R.
life_history_analysis.R  |  Script to calculate life history metrics from matrices, and perform phylogenetic mixed-effects model regressions between slow-fast life history traits and the magnitude of symbiont-mediated variance buffering effects. Includes multiple version of this analysis, including early attempts using phylogenetic corrections.
seed_mean.stan | Stan model for mean seed per spikelet estimates for multiple species with normal distribution and no random effects
seed_means.R | Script to analyse mean seed per spikelet
seed_to_seedling.R | Script to analysis seed-to-seedling transition vital rate model.
seedling_grow.stan | vital rate model for seedling growth with negative binomial distribution
seedling_grow.stan | vital rate model for seedling growth with the PIG distribution
seedling_surv.stan | vital rate model for seedling survival
stochastic_lambda_analysis.R | old version of stochastic growth rate decomposition
stroma_observations_summary.R | Script to summarize the frequency of stroma observations in the plots, and of the fidelity of symbiont status from opportunistic endophyte scoring over the course of the experiment.



### Manuscript 
holds drafts and figures of manuscript. Final manuscript is in folder submission_to_EcolLetters

File Name  | Description
------------- | -------------
EndoStochDemo.bib | bib file containing references
Endo_Stoch_Demo_EcolLetters_revision2.tex | tex file for compiling final submission to journal.
Other misc. files | Various style and log files for latex, based around formatting (sn-jnl.csl is necessary to compile), getting the appendix TOC to compile requires running first with an overall TOC then removing this line and re-compiling.


### Meeting Notes (last updated 2021)
holds weekly-ish Josh and Tom's meeting notes, including weekly goals and discussions beyond this project

File Name  | Description
------------- | -------------
MeetingUpdates.tex | Latex document for recording weekly notes and goals
MeetingUpdates.pdf | PDF document output of MeetingUpdates.tex
Other misc. files | Various style and log files for latex



