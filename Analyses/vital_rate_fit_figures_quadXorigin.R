## Title: Grass endophyte vital rate fit visualizations and plots for all vital rates and all species
## Authors: Joshua and Tom
#############################################################


library(tidyverse)
library(RColorBrewer)
library(viridis)
library(rstan)
library(StanHeaders)
library(bayesplot)
library(devtools)
library(moments)
library(patchwork)
library(actuar) # for PIG distribution 


invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }
Lkurtosis=function(x) log(kurtosis(x)); 

quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}
#############################################################################################
####### Reading in the demographic data ------------------
#############################################################################################

#  data are prepared in the endodemog_data_processing.R file, 
# source("Analyses/endodemog_data_processing.R")
species_factor_key <- c("AGPE" = 1, "ELRI" = 2, "ELVI" = 3, 
                        "FESU" = 4, "LOAR" = 5, "POAL" = 6, 
                        "POSY" = 7)
LTREB_full <- read_csv("~/Dropbox/EndodemogData/Fulldataplusmetadata/LTREB_full.csv")
LTREB_repro1 <- read_csv("~/Dropbox/EndodemogData/Fulldataplusmetadata/LTREB_repro1.csv")

## Clean up the main data frame for NA's, other small data entry errors
LTREB_data_forsurv <- LTREB_full %>% 
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(endo_01)) %>%   # There are a few LOAR that don't have a plot level endo assigned
  filter(origin_01 == 1 & year_t != birth | origin_01 == 0)  # filtering out first year germinants (including those that are bigger than 1 tiller)
dim(LTREB_data_forsurv)

LTREB_surv_seedling <- LTREB_full %>% 
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(endo_01)) %>%  # There are a few LOAR that don't have a plot level endo assigned
  filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
  filter(logsize_t == 0) # this is filtering out the plants that are "recruits" but are larger than 1 tiller
dim(LTREB_surv_seedling)

# I want to look at these "seedlings" that are bigger than 1 tiller, We are dropping these from our seedling model, as they are most likely missed plants from the previous year
LTREB_surv_big_seedling <- LTREB_full %>% 
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(endo_01)) %>%  # There are a few LOAR that don't have a plot level endo assigned
  filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
  filter(logsize_t != 0) 
dim(LTREB_surv_big_seedling)
table(LTREB_surv_big_seedling$year_t,LTREB_surv_big_seedling$size_t, LTREB_surv_big_seedling$species)

LTREB_data_forflw <- LTREB_full %>% 
  filter(!is.na(FLW_STAT_T1)) %>% 
  filter(!is.na(logsize_t1)) %>% 
  filter(!is.na(endo_01))
dim(LTREB_data_forflw)


LTREB_data_forgrow <- LTREB_full %>%
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(size_t1)) %>% 
  filter(!is.na(endo_01)) %>% 
  filter(origin_01 == 1 & year_t != birth | origin_01 == 0)  # filtering out first year germinants (including those that are bigger than 1 tiller)
dim(LTREB_data_forgrow)

LTREB_grow_seedling <- LTREB_full %>%
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(size_t1)) %>% 
  filter(!is.na(endo_01)) %>% 
  filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
  filter(logsize_t == 0) # this is filtering out the plants that are "recruits" but are larger than 1 tiller

dim(LTREB_grow_seedling)


LTREB_data_forfert <- LTREB_full %>% 
  filter(!is.na(FLW_COUNT_T1)) %>% 
  filter(FLW_COUNT_T1 > 0) %>% 
  filter(!is.na(logsize_t1))
dim(LTREB_data_forfert)

LTREB_data_forspike <- LTREB_full %>%
  dplyr::select(-FLW_COUNT_T, -FLW_STAT_T, -SPIKE_A_T, -SPIKE_B_T, -SPIKE_C_T, -SPIKE_D_T, -SPIKE_AGPE_MEAN_T, -census_month, -year, -spei1, -spei3, -spei12, -spei24, -annual_temp, -annual_precip, -endo_status_from_check, -plot_endo_for_check, -endo_mismatch, -dist_a, -dist_b) %>% 
  filter(!is.na(FLW_STAT_T1)) %>% 
  filter(FLW_STAT_T1>0) %>% 
  reshape2::melt(id.var = c("plot_fixed" ,   "plot_index",         "pos"         ,           "id",
                  "species"       ,         "species_index"  ,        "endo_01",
                  "endo_index"  ,           "origin_01"       ,       "birth" ,
                  "year_t1"         ,       "year_t1_index"       ,   "surv_t1" ,
                  "size_t1"         ,       "logsize_t1"       ,
                  "year_t",
                  "year_t_index"     ,      "size_t"           ,      "logsize_t"  ,
                  "FLW_COUNT_T1"      ,      "FLW_STAT_T1"),
       value.name = "spike_count_t1") %>% 
  rename(spikelet_id = variable) %>% 
  filter(!is.na(spike_count_t1)) %>% 
  mutate(spike_count_t1 = as.integer(spike_count_t1))


# Create data lists to be used for the Stan model

surv_data_list <- list(y = LTREB_data_forsurv$surv_t1,
                       logsize = LTREB_data_forsurv$logsize_t,
                       origin_01 = LTREB_data_forsurv$origin_01,
                       endo_01 = as.integer(LTREB_data_forsurv$endo_01),
                       endo_index = as.integer(LTREB_data_forsurv$endo_index),
                       spp = as.integer(LTREB_data_forsurv$species_index),
                       year_t = as.integer(LTREB_data_forsurv$year_t_index),
                       plot = as.integer(LTREB_data_forsurv$plot_index),
                       N = nrow(LTREB_data_forsurv),
                       nSpp = length(unique(LTREB_data_forsurv$species_index)),
                       nYear = max(unique(LTREB_data_forsurv$year_t_index)),
                       nPlot = max(unique(LTREB_data_forsurv$plot_index)),
                       nEndo =   length(unique(LTREB_data_forsurv$endo_01)))
str(surv_data_list)

seed_surv_data_list <- list(y = LTREB_surv_seedling$surv_t1,
                            logsize = LTREB_surv_seedling$logsize_t,
                            origin_01 = LTREB_surv_seedling$origin_01,
                            endo_01 = as.integer(LTREB_surv_seedling$endo_01),
                            endo_index = as.integer(LTREB_surv_seedling$endo_index),
                            spp = as.integer(LTREB_surv_seedling$species_index),
                            year_t = as.integer(LTREB_surv_seedling$year_t_index),
                            plot = as.integer(LTREB_surv_seedling$plot_index),
                            N = nrow(LTREB_surv_seedling),
                            nSpp = length(unique(LTREB_surv_seedling$species_index)),
                            nYear = max(unique(LTREB_surv_seedling$year_t_index)),
                            nPlot = max(unique(LTREB_surv_seedling$plot_index)),
                            nEndo =   length(unique(LTREB_surv_seedling$endo_01)))
str(seed_surv_data_list)


flw_data_list <- list(y = LTREB_data_forflw$FLW_STAT_T1,
                      logsize = LTREB_data_forflw$logsize_t1,
                      origin_01 = LTREB_data_forflw$origin_01,
                      endo_01 = as.integer(LTREB_data_forflw$endo_01),
                      endo_index = as.integer(LTREB_data_forflw$endo_index),
                      spp = as.integer(LTREB_data_forflw$species_index),
                      year_t = as.integer(LTREB_data_forflw$year_t_index),
                      plot = as.integer(LTREB_data_forflw$plot_index),
                      N = nrow(LTREB_data_forflw),
                      nSpp = length(unique(LTREB_data_forflw$species_index)),
                      nYear = max(unique(LTREB_data_forflw$year_t_index)),
                      nPlot = length(unique(LTREB_data_forflw$plot_index)),
                      nEndo =   length(unique(LTREB_data_forflw$endo_01)))
str(flw_data_list)

grow_data_list <- list(y = as.integer(LTREB_data_forgrow$size_t1),
                       logsize = LTREB_data_forgrow$logsize_t,
                       origin_01 = as.integer(LTREB_data_forgrow$origin_01),
                       endo_01 = as.integer(LTREB_data_forgrow$endo_01),
                       endo_index = as.integer(LTREB_data_forgrow$endo_index),
                       spp = as.integer(LTREB_data_forgrow$species_index),
                       year_t = as.integer(LTREB_data_forgrow$year_t_index),
                       plot = as.integer(LTREB_data_forgrow$plot_index),
                       N = nrow(LTREB_data_forgrow),
                       nSpp = length(unique(LTREB_data_forgrow$species_index)),
                       nYear = max(unique(LTREB_data_forgrow$year_t_index)),
                       nPlot = max(unique(LTREB_data_forgrow$plot_index)),
                       nEndo =   length(unique(LTREB_data_forgrow$endo_01)))
str(grow_data_list)
seed_grow_data_list <- list(y = as.integer(LTREB_grow_seedling$size_t1),
                            logsize = LTREB_grow_seedling$logsize_t,
                            origin_01 = as.integer(LTREB_grow_seedling$origin_01),
                            endo_01 = as.integer(LTREB_grow_seedling$endo_01),
                            endo_index = as.integer(LTREB_grow_seedling$endo_index),
                            spp = as.integer(LTREB_grow_seedling$species_index),
                            year_t = as.integer(LTREB_grow_seedling$year_t_index),
                            plot = as.integer(LTREB_grow_seedling$plot_index),
                            N = nrow(LTREB_grow_seedling),
                            nSpp = length(unique(LTREB_grow_seedling$species_index)),
                            nYear = max(unique(LTREB_grow_seedling$year_t_index)),
                            nPlot = max(unique(LTREB_grow_seedling$plot_index)),
                            nEndo =   length(unique(LTREB_grow_seedling$endo_01)))
str(seed_grow_data_list)



fert_data_list <- list(y = as.integer(LTREB_data_forfert$FLW_COUNT_T1),
                       logsize = LTREB_data_forfert$logsize_t1,
                       origin_01 = LTREB_data_forfert$origin_01,
                       endo_01 = as.integer(LTREB_data_forfert$endo_01),
                       endo_index = as.integer(LTREB_data_forfert$endo_index),
                       spp = as.integer(LTREB_data_forfert$species_index),
                       year_t = as.integer(LTREB_data_forfert$year_t_index),
                       plot = as.integer(LTREB_data_forfert$plot_index),
                       N = nrow(LTREB_data_forfert),
                       nSpp = length(unique(LTREB_data_forfert$species_index)),
                       nYear = max(unique(LTREB_data_forfert$year_t_index)),
                       nPlot = max(unique(LTREB_data_forfert$plot_index)),
                       nEndo =   length(unique(LTREB_data_forfert$endo_01)))
str(fert_data_list)



spike_data_list <- list(nYear = max(unique(LTREB_data_forspike$year_t_index)),
                        nPlot = max(unique(LTREB_data_forspike$plot_index)),
                        nSpp = length(unique(LTREB_data_forspike$species)),
                        nEndo=length(unique(LTREB_data_forspike$endo_01)),
                        N = nrow(LTREB_data_forspike),
                        year_t = as.integer(LTREB_data_forspike$year_t_index),
                        plot = as.integer(LTREB_data_forspike$plot_index),
                        spp = as.integer(LTREB_data_forspike$species_index),
                        y = LTREB_data_forspike$spike_count_t1,
                        logsize = LTREB_data_forspike$logsize_t1,
                        endo_01 = LTREB_data_forspike$endo_01,
                        origin_01 = LTREB_data_forspike$origin_01)
str(spike_data_list)

LTREB_data_for_seedmeans <- LTREB_repro1 %>%  
  mutate(seed_per_spike = seed/spikelets) %>% 
  mutate(SEED_PER_SPIKE= case_when(species != "AGPE" ~ seed_per_spike,
                                   species == "AGPE" & tillerid_fixed == "multitillermean" ~ seed, # AGPE has some of its seeds data already recorded as seed/spike
                                   species == "AGPE" & tillerid_fixed != "multitillermean" ~ seed_per_spike)) %>% 
  mutate(species_index = as.integer(recode(species, !!!species_factor_key))) %>%
  mutate(endo_index = as.integer(as.factor(endo_01+1)))  %>% 
  filter(!is.na(SEED_PER_SPIKE), SEED_PER_SPIKE >0)

dim(LTREB_data_for_seedmeans)
# View(LTREB_data_for_seedmeans)

# Creating data list to be passed to the model
seed_mean_data_list <- list(seed = LTREB_data_for_seedmeans$SEED_PER_SPIKE,
                            endo_01 = LTREB_data_for_seedmeans$endo_01,
                            endo_index = LTREB_data_for_seedmeans$endo_index,
                            year = LTREB_data_for_seedmeans$year,
                            plot = LTREB_data_for_seedmeans$plot_fixed,
                            species = LTREB_data_for_seedmeans$species_index,
                            N = length(na.omit(LTREB_data_for_seedmeans$SEED_PER_SPIKE)),
                            nYear = length(unique(LTREB_data_for_seedmeans$year)),
                            nPlot = length(unique(LTREB_data_for_seedmeans$plot_fixed)),
                            nSpecies = length(unique(LTREB_data_for_seedmeans$species_index)),
                            nEndo =   length(unique(LTREB_data_for_seedmeans$endo_01)))
str(seed_mean_data_list)

s_to_s_data_list <- read_rds("Analyses/s_to_s_data_list.rds")

#############################################################################################
####### Read in matrix population functions ------------------
#############################################################################################
# we can use functions from this to set up our parameters for the plots with mean and variance effects
source("Analyses/MPM_functions.R")


#############################################################################################
####### Read in Stan vital rate model outputs ------------------
#############################################################################################

# The stan objects for each vital rate
# surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_surv.rds")
# surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_surv_woseedling.rds")
# grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_grow_PIG_10000iterations.rds")
# grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow_PIG.rds")
# flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_flw.rds")
# fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_PIG.rds")
# spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike_year_plot_nb.rds")
# seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_seed_mean.rds")
# stos_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_s_to_s.rds") 

# surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_surv.rds")
# surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_surv_woseedling_quad.rds")
# grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_grow_PIG_10000iterations.rds")
# grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow_PIG_quad.rds")
# flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_flw_quad.rds")
# fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_PIG_quad.rds")
# spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike_year_plot_nb_quad.rds")
# seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_seed_mean.rds")
# stos_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_s_to_s.rds") 


surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_surv.rds")
surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_surv_woseedling_quadXorigin.rds")
grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_grow_PIG_10000iterations.rds")
grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow_PIG_quadXorigin.rds")
flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_flw_quadXorigin.rds")
fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_PIG_quadXorigin.rds")
spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike_year_plot_nb_quadXorigin.rds")
seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_seed_mean.rds")
stos_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_s_to_s.rds") 

# Pulling out the actual parameters
surv_par <- rstan::extract(surv_fit, pars =quote_bare(beta0,betasize,betasize_2,betaendo,
                                                      tau_year, tau_plot, sigma_year, sigma0, sigmaendo))
surv_sdlg_par <- rstan::extract(surv_fit_seedling, pars =quote_bare(beta0,betaendo,
                                                                    tau_year, tau_plot, sigma_year, sigma0, sigmaendo))
grow_par <- rstan::extract(grow_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo,
                                                       tau_year, tau_plot,
                                                       sigma, sigma_year, sigma0, sigmaendo))
grow_sdlg_par <- rstan::extract(grow_fit_seedling, pars = quote_bare(beta0,betaendo,
                                                                     tau_year, tau_plot,
                                                                     sigma, sigma_year, sigma0, sigmaendo))
flow_par <- rstan::extract(flw_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo,
                                                      tau_year, tau_plot, sigma_year, sigma0, sigmaendo))
fert_par <- rstan::extract(fert_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo,
                                                       tau_year, tau_plot, sigma_year, sigma0, sigmaendo))
spike_par <- rstan::extract(spike_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo,
                                                         tau_year, tau_plot,
                                                         phi, sigma_year, sigma0, sigmaendo))
seed_par <- rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
recruit_par <- rstan::extract(stos_fit, pars = quote_bare(beta0,betaendo,
                                                          tau_year, tau_plot, sigma_year, sigma0, sigmaendo))

# Saved y_rep values from the script "vital_rate_analysis.R", "seed_means.R" and "seed_to_seedling.R"
y_s_sim <- readRDS(file = "yrep_survivalmodel_quadXorigin.rds")
y_seed_s_sim <- readRDS(file = "yrep_seedlingsurvivalmodel.rds")
y_f_sim <- readRDS(file = "yrep_floweringmodel_quadXorigin.rds")
y_g_sim <- readRDS(file = "yrep_growthPIGmodel_quadXorigin.rds")
y_seed_g_sim <- readRDS(file = "yrep_seedlinggrowthPIGmodel.rds")
y_fert_sim <- readRDS(file = "yrep_fertilityPIGmodel_quadXorigin.rds")
y_spike_sim <- readRDS(file = "yrep_spikeletNBmodel_quadXorigin.rds")
y_seedmean_sim <- readRDS(file = "yrep_seedmeanmodel.rds")
y_recruit_sim <- readRDS(file = "yrep_stosmodel.rds")

#########################################################################################################
# Plots for all species all vital rates and model fits ------------------------------
#########################################################################################################
# Set color scheme based on analine blue
endophyte_color_scheme <- c("#fdedd3","#f3c8a8", "#5a727b", "#4986c7", "#181914",  "#163381")
color_scheme_set(endophyte_color_scheme)
color_scheme_view()

twotone_endophyte_color_scheme <- c("#FDC800", "#0082EA")

# And creating a color palette for each year
yearcount = length(unique(LTREB_full$year_t))
yearcolors<- colorRampPalette(brewer.pal(8,"Dark2"))(yearcount)
# scales::show_col(yearcolors)
species_list <- c("A. perennans", "E. villosus", "E. virginicus", "F. subverticillata", "L. arundinaceae", "P. alsodes", "P. sylvestris")
species_code_list <- c("AGPE", "ELRI", "ELVI", "FESU", "LOAR", "POAL", "POSY")

species_colors <- c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")


## Plot for all models vital rate fits #####
# # survival
# surv_densplot <- ppc_dens_overlay(surv_data_list$y, y_s_sim) + theme_classic() + labs(title = "Adult Survival", x = "Survival status", y = "Density")
# # surv_densplot
# # ggsave(surv_densplot, filename = "surv_densplot.png", width = 4, height = 4)
# 
# mean_s_plot <-   ppc_stat(surv_data_list$y, y_s_sim, stat = "mean") + legend_none()+labs(title = "Mean")
# sd_s_plot <- ppc_stat(surv_data_list$y, y_s_sim, stat = "sd")+ legend_none()+labs(title = "SD")
# skew_s_plot <- ppc_stat(surv_data_list$y, y_s_sim, stat = "skewness")+ legend_none()+labs(title = "Skew")
# kurt_s_plot <- ppc_stat(surv_data_list$y, y_s_sim, stat = "Lkurtosis")+ legend_none()+labs(title = "Kurtosis")
# surv_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Survival")
# # surv_moments
# # ggsave(surv_moments, filename = "surv_momentplot.png", width = 4, height = 4)
# 
# #seedliing survival
# seedsurv_densplot <- ppc_dens_overlay(seed_surv_data_list$y, y_seed_s_sim) + theme_classic() + labs(title = "Seedling Survival", x = "Survival status", y = "Density")
# # seedsurv_densplot
# # ggsave(seedsurv_densplot, filename = "seedsurv_densplot.png", width = 4, height = 4)
# 
# mean_s_plot <-   ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "mean") + legend_none()+labs(title = "Mean")
# sd_s_plot <- ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "sd")+ legend_none()+labs(title = "SD")
# skew_s_plot <- ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "skewness")+ legend_none()+labs(title = "Skew")
# kurt_s_plot <- ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "Lkurtosis")+ legend_none()+labs(title = "Kurtosis")
# seedsurv_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Seedling Survival")
# # seedsurv_moments
# # ggsave(seedsurv_moments, filename = "seedsurv_momentsplot.png", width = 4, height = 4)
# 
# # Flowering
# flw_densplot <- ppc_dens_overlay(flw_data_list$y, y_f_sim) + theme_classic() + labs(title = "Flowering", x = "Flowering status", y = "Density")
# # flw_densplot
# # ggsave(flw_densplot, filename = "flw_densplot.png", width = 4, height = 4)
# 
# mean_f_plot <-   ppc_stat(flw_data_list$y, y_f_sim, stat = "mean") + legend_none()+labs(title = "Mean")
# sd_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "sd")+ legend_none()+labs(title = "SD")
# skew_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "skewness")+ legend_none()+labs(title = "Skew")
# kurt_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "Lkurtosis")+ legend_none()+labs(title = "Kurtosis")
# flw_moments <- mean_f_plot+sd_f_plot+skew_f_plot+kurt_f_plot +plot_annotation(title = "Flowering")
# # flw_moments
# # ggsave(flw_moments, filename = "flw_momentsplot.png", width = 4, height = 4)
# 
# # Growth (with the poisson inverse gaussian distribution)
# grow_densplot <- ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(0,60) + theme_classic() + labs(title = "Adult Growth ", x = "No. of Tillers", y = "Density")
# # grow_densplot
# # ggsave(grow_densplot, filename = "grow_densplot.png", width = 4, height = 4)
# 
# mean_g_plot <-   ppc_stat(grow_data_list$y, y_g_sim, stat = "mean") + legend_none()+labs(title = "Mean")
# sd_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "sd")+ legend_none()+labs(title = "SD")
# skew_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "skewness")+ legend_none()+labs(title = "Skew")
# kurt_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "Lkurtosis")+ legend_none()+labs(title = "Kurtosis")
# grow_moments <- mean_g_plot+sd_g_plot+skew_g_plot+kurt_g_plot+ plot_annotation(title = "Growth")
# # grow_moments
# # ggsave(grow_moments, filename = "grow_momentsplot.png", width = 4, height = 4)
# 
# # seedling growth (with poisson inverse gaussian distribution)
# seedgrow_densplot <- ppc_dens_overlay(seed_grow_data_list$y, y_seed_g_sim) + xlim(0,60) + theme_classic() + labs(title = "Seedling Growth", x = "No. of Tillers", y = "Density")
# # seedgrow_densplot
# # ggsave(seedgrow_densplot, filename = "seed_grow_densplot.png", width = 4, height = 4)
# 
# mean_seed_g_plot <-   ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "mean")+ legend_none()+labs(title = "Mean")
# sd_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "sd")+ legend_none()+labs(title = "SD")
# skew_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "skewness")+ legend_none()+labs(title = "Skew")
# kurt_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "Lkurtosis")+ legend_none()+labs(title = "Kurtosis")
# seedgrow_moments <- mean_seed_g_plot+sd_seed_g_plot+skew_seed_g_plot+kurt_seed_g_plot+ plot_annotation(title = "Seedling Growth")
# # seedgrow_moments
# # ggsave(seedgrow_moments, filename = "seedgrow_momentsplot.png", width = 4, height = 4)
# 
# # Fertility
# fert_densplot <- ppc_dens_overlay(fert_data_list$y, y_fert_sim) + xlim(0,40) + theme_classic() + labs(title = "Fertility", x = "No. of Reproductive Tillers", y = "Density")
# # fert_densplot
# # ggsave(fert_densplot, filename = "fert_densplot_withPIG.png", width = 4, height = 4)
# 
# # This doesn't fit all the moments quite as well, but it's likely good enough
# mean_fert_plot <-   ppc_stat(fert_data_list$y, y_fert_sim, stat = "mean")+ legend_none()+labs(title = "Mean")
# sd_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "sd")+ legend_none()+labs(title = "SD")
# skew_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "skewness")+ legend_none()+labs(title = "Skew")
# kurt_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "Lkurtosis")+ legend_none()+labs(title = "Kurtosis")
# fert_moments <- mean_fert_plot+sd_fert_plot+skew_fert_plot+kurt_fert_plot + plot_annotation(title = "Fertility")
# # fert_moments
# # ggsave(fert_moments, filename = "fert_momentsplot_withPIG.png", width = 4, height = 4)
# 
# # Spikelets per infl
# spike_densplot <- ppc_dens_overlay(spike_data_list$y, y_spike_sim) + xlim(0,250) + theme_classic() + labs(title = "Spikelet Count", x = "Spikelets/Inflorescence", y = "Density")
# # spike_densplot
# # ggsave(spike_densplot, filename = "spike_densplot.png", width = 4, height = 4)
# 
# mean_spike_plot <-   ppc_stat(spike_data_list$y, y_spike_sim, stat = "mean")+ legend_none()+labs(title = "Mean")
# sd_spike_plot <- ppc_stat(spike_data_list$y, y_spike_sim, stat = "sd")+ legend_none()+labs(title = "SD")
# skew_spike_plot <- ppc_stat(spike_data_list$y, y_spike_sim, stat = "skewness")+ legend_none()+labs(title = "Skew")
# kurt_spike_plot <- ppc_stat(spike_data_list$y, y_spike_sim, stat = "Lkurtosis")+ legend_none()+labs(title = "Kurtosis")
# spike_moments <- mean_spike_plot+sd_spike_plot+skew_spike_plot+kurt_spike_plot+ plot_annotation(title = "Spikelets per Infl.")
# # spike_moments
# # ggsave(spike_moments, filename = "spike_momentsplot.png", width = 4, height = 4)
# 
# # Seed mean per spikelet
# seedmean_densplot <-ppc_dens_overlay(seed_mean_data_list$seed, y_seedmean_sim)+ theme_classic() + labs(title = "Seed Means", x = "Seeds per Spikelet.", y = "Density") 
# # ggsave(seedmean_densplot, filename = "seedmean_densplot.png", width = 4, height = 4)
# 
# mean_sm_plot <-   ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "mean")+ legend_none()+labs(title = "Mean")
# sd_sm_plot <- ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "sd")+ legend_none()+labs(title = "SD")
# skew_sm_plot <- ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "skewness")+ legend_none()+labs(title = "Skew")
# kurt_sm_plot <- ppc_stat(seed_mean_data_list$seed, y_seedmean_sim, stat = "Lkurtosis")+ legend_none()+labs(title = "Kurtosis")
# seedmean_moments <- mean_sm_plot + sd_sm_plot + skew_sm_plot + kurt_sm_plot
# # seedmean_moments
# # ggsave(seedmean_moments, filename = "seedmean_moments.png", width = 4, height = 4)
# 
# # Seed to Seedling
# stos_densplot <- ppc_dens_overlay(s_to_s_data_list$tot_recruit_t1, y_recruit_sim) +xlim(0,40) + theme_classic() + labs(title = "Recruitment", x = "Successful Germination", y = "Density") 
# # ggsave(stos_densplot, filename = "stos_densplot.png", width = 4, height = 4)
# 
# mean_stos_plot <-   ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "mean")+ legend_none()+labs(title = "Mean")
# sd_stos_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "sd")+ legend_none()+labs(title = "SD")
# skew_stos_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "skewness")+ legend_none()+labs(title = "Skew")
# kurt_stos_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "Lkurtosis")+ legend_none()+labs(title = "Kurtosis")
# stos_moments <- mean_stos_plot+sd_stos_plot+skew_stos_plot+kurt_stos_plot
# # stos_moments
# # ggsave(stos_moments, filename = "stos_moments.png", width = 4, height = 4)
# 
# 
# # all together
# ## Plot for all models fits and moments
# fm_layout <- "
# AB
# CD
# EG
# HI
# J#
# "
# fitsandmoments_plot <- wrap_plots(A = (surv_densplot | surv_moments), B = (seedsurv_densplot | seedsurv_moments), 
#   C = (grow_densplot | grow_moments), D = (seedgrow_densplot | seedgrow_moments),
#   E =(flw_densplot | flw_moments), G = (fert_densplot | fert_moments),
#   H = (spike_densplot | spike_moments), I = (seedmean_densplot | seedmean_moments),
#   J = (stos_densplot | stos_moments), design = fm_layout, guides = "collect") +
#   plot_annotation(title = "Vital rate fits and moments with 500 posterior draws")
# ggsave(fitsandmoments_plot, filename = "fitsandmoments_plot.png", width = 12, height = 15)
# 
# 


### Plot for vital rate fits vs simulated data for each species separately. ####

# First making a function to filter our vital rate data and the simulations we already have by species
filterspecies_densplot <- function(data, data_list, sim, species = NA, x, y, xlim = NA){
  require(tidyverse); require(bayesplot); require(patchwork)
  filtered_data_list <- lapply(data_list, function(x) x[data$species == species])
  
  sim_list <- as.list(data.frame(t(sim)))
  filtered_sim_list <- lapply(sim_list, function(x) x[data$species == species])
  
  filtered_sim <- matrix(unlist(filtered_sim_list), ncol = length(filtered_data_list$y))
  min_x <- min(filtered_sim)
  
  densplot <- ppc_dens_overlay(filtered_data_list$y, filtered_sim) + theme_classic() + xlim(min_x,xlim) + labs(subtitle = species, x = x, y = y)

  mean_plot <-   ppc_stat(filtered_data_list$y, filtered_sim, stat = "mean") + legend_none()+labs(x = "Mean")
  sd_plot <- ppc_stat(filtered_data_list$y, filtered_sim, stat = "sd")+ legend_none()+labs(x = "SD")
  skew_plot <- ppc_stat(filtered_data_list$y, filtered_sim, stat = "skewness")+ legend_none()+labs(x = "Skew")
  kurt_plot <- ppc_stat(filtered_data_list$y, filtered_sim, stat = "Lkurtosis")+ legend_none()+labs(x = "Kurtosis")
  moments <- mean_plot+sd_plot+skew_plot+kurt_plot + plot_layout(nrow = 1)
  # 
  densplot_combo <- densplot + moments
  return(densplot_combo)
}


# Adult Survival
AGPE_surv_densplot <- filterspecies_densplot(data = LTREB_data_forsurv, data_list = surv_data_list, sim = y_s_sim, 
                                             species = "AGPE", x = "Survival Status", y = "Density")
ELRI_surv_densplot <- filterspecies_densplot(data = LTREB_data_forsurv, data_list = surv_data_list, sim = y_s_sim, 
                                             species = "ELRI", x = "Survival Status", y = "Density")
ELVI_surv_densplot <- filterspecies_densplot(data = LTREB_data_forsurv, data_list = surv_data_list, sim = y_s_sim, 
                                             species = "ELVI", x = "Survival Status", y = "Density")
FESU_surv_densplot <- filterspecies_densplot(data = LTREB_data_forsurv, data_list = surv_data_list, sim = y_s_sim, 
                                             species = "FESU", x = "Survival Status", y = "Density")
LOAR_surv_densplot <- filterspecies_densplot(data = LTREB_data_forsurv, data_list = surv_data_list, sim = y_s_sim, 
                                             species = "LOAR", x = "Survival Status", y = "Density")
POAL_surv_densplot <- filterspecies_densplot(data = LTREB_data_forsurv, data_list = surv_data_list, sim = y_s_sim, 
                                             species = "POAL", x = "Survival Status", y = "Density")
POSY_surv_densplot <- filterspecies_densplot(data = LTREB_data_forsurv, data_list = surv_data_list, sim = y_s_sim, 
                                             species = "POSY", x = "Survival Status", y = "Density")

fm_layout <- "
A
B
C
D
E
G
H"
survbyspecies_densplot <- wrap_plots(A = (AGPE_surv_densplot), B = (ELRI_surv_densplot), C = (ELVI_surv_densplot), D = (FESU_surv_densplot), E = (LOAR_surv_densplot), G = (POAL_surv_densplot), H = (POSY_surv_densplot),
design = fm_layout, guides = "collect") + plot_annotation(title = "Adult Survival", subtitle = "Vital rate fits and moments with 500 posterior draws")

ggsave(survbyspecies_densplot, filename = "survbyspecies_densplot_quadXorigin.png", width = 8, height = 10)


# Seedling Survival
AGPE_seedsurv_densplot <- filterspecies_densplot(data = LTREB_surv_seedling, data_list = seed_surv_data_list, sim = y_seed_s_sim, 
                                             species = "AGPE", x = "Survival Status", y = "Density")
ELRI_seedsurv_densplot <- filterspecies_densplot(data = LTREB_surv_seedling, data_list = seed_surv_data_list, sim = y_seed_s_sim, 
                                             species = "ELRI", x = "Survival Status", y = "Density")
ELVI_seedsurv_densplot <- filterspecies_densplot(data = LTREB_surv_seedling, data_list = seed_surv_data_list, sim = y_seed_s_sim, 
                                             species = "ELVI", x = "Survival Status", y = "Density")
FESU_seedsurv_densplot <- filterspecies_densplot(data = LTREB_surv_seedling, data_list = seed_surv_data_list, sim = y_seed_s_sim, 
                                             species = "FESU", x = "Survival Status", y = "Density")
LOAR_seedsurv_densplot <- filterspecies_densplot(data = LTREB_surv_seedling, data_list = seed_surv_data_list, sim = y_seed_s_sim, 
                                             species = "LOAR", x = "Survival Status", y = "Density")
POAL_seedsurv_densplot <- filterspecies_densplot(data = LTREB_surv_seedling, data_list = seed_surv_data_list, sim = y_seed_s_sim, 
                                             species = "POAL", x = "Survival Status", y = "Density")
POSY_seedsurv_densplot <- filterspecies_densplot(data = LTREB_surv_seedling, data_list = seed_surv_data_list, sim = y_seed_s_sim, 
                                             species = "POSY", x = "Survival Status", y = "Density")


seedsurvbyspecies_densplot <- wrap_plots(A = (AGPE_seedsurv_densplot), B = (ELRI_seedsurv_densplot), C = (ELVI_seedsurv_densplot), D = (FESU_seedsurv_densplot), E = (LOAR_seedsurv_densplot), G = (POAL_seedsurv_densplot), H = (POSY_seedsurv_densplot),
                                     design = fm_layout, guides = "collect") + plot_annotation(title = "Seedling Survival", subtitle = "Vital rate fits and moments with 500 posterior draws")

ggsave(seedsurvbyspecies_densplot, filename = "seedsurvbyspeciesXdensplot.png", width = 8, height = 10)

# Flowering

AGPE_flw_densplot <- filterspecies_densplot(data = LTREB_data_forflw, data_list = flw_data_list, sim = y_f_sim, 
                                             species = "AGPE", x = "Flowering Status", y = "Density")
ELRI_flw_densplot <- filterspecies_densplot(data = LTREB_data_forflw, data_list = flw_data_list, sim = y_f_sim, 
                                             species = "ELRI", x = "Flowering Status", y = "Density")
ELVI_flw_densplot <- filterspecies_densplot(data = LTREB_data_forflw, data_list = flw_data_list, sim = y_f_sim, 
                                             species = "ELVI", x = "Flowering Status", y = "Density")
FESU_flw_densplot <- filterspecies_densplot(data = LTREB_data_forflw, data_list = flw_data_list, sim = y_f_sim, 
                                             species = "FESU", x = "Flowering Status", y = "Density")
LOAR_flw_densplot <- filterspecies_densplot(data = LTREB_data_forflw, data_list = flw_data_list, sim = y_f_sim, 
                                             species = "LOAR", x = "Flowering Status", y = "Density")
POAL_flw_densplot <- filterspecies_densplot(data = LTREB_data_forflw, data_list = flw_data_list, sim = y_f_sim, 
                                             species = "POAL", x = "Flowering Status", y = "Density")
POSY_flw_densplot <- filterspecies_densplot(data = LTREB_data_forflw, data_list = flw_data_list, sim = y_f_sim, 
                                             species = "POSY", x = "Flowering Status", y = "Density")

flwbyspecies_densplot <- wrap_plots(A = (AGPE_flw_densplot), B = (ELRI_flw_densplot), C = (ELVI_flw_densplot), D = (FESU_flw_densplot), E = (LOAR_flw_densplot), G = (POAL_flw_densplot), H = (POSY_flw_densplot),
                                     design = fm_layout, guides = "collect") + plot_annotation(title = "Flowering Probability", subtitle = "Vital rate fits and moments with 500 posterior draws")

ggsave(flwbyspecies_densplot, filename = "flwbyspecies_densplot_quadXorigin.png", width = 8, height = 10)



# Growth

AGPE_grow_densplot <- filterspecies_densplot(data = LTREB_data_forgrow, data_list = grow_data_list, sim = y_g_sim, 
                                            species = "AGPE", x = "Adult Growth", y = "Density", xlim = 50)
ELRI_grow_densplot <- filterspecies_densplot(data = LTREB_data_forgrow, data_list = grow_data_list, sim = y_g_sim, 
                                            species = "ELRI", x = "Adult Growth", y = "Density", xlim = 30)
ELVI_grow_densplot <- filterspecies_densplot(data = LTREB_data_forgrow, data_list = grow_data_list, sim = y_g_sim, 
                                            species = "ELVI", x = "Adult Growth", y = "Density", xlim = 30)
FESU_grow_densplot <- filterspecies_densplot(data = LTREB_data_forgrow, data_list = grow_data_list, sim = y_g_sim, 
                                            species = "FESU", x = "Adult Growth", y = "Density", xlim = 40)
LOAR_grow_densplot <- filterspecies_densplot(data = LTREB_data_forgrow, data_list = grow_data_list, sim = y_g_sim, 
                                            species = "LOAR", x = "Adult Growth", y = "Density", xlim = 50)
POAL_grow_densplot <- filterspecies_densplot(data = LTREB_data_forgrow, data_list = grow_data_list, sim = y_g_sim, 
                                            species = "POAL", x = "Adult Growth", y = "Density", xlim = 60)
POSY_grow_densplot <- filterspecies_densplot(data = LTREB_data_forgrow, data_list = grow_data_list, sim = y_g_sim, 
                                            species = "POSY", x = "Adult Growth", y = "Density", xlim = 50)


growbyspecies_densplot <- wrap_plots(A = (AGPE_grow_densplot), B = (ELRI_grow_densplot), C = (ELVI_grow_densplot), D = (FESU_grow_densplot), E = (LOAR_grow_densplot), G = (POAL_grow_densplot), H = (POSY_grow_densplot),
                                    design = fm_layout, guides = "collect") + plot_annotation(title = "Adult Growth", subtitle = "Vital rate fits and moments with 500 posterior draws")

ggsave(growbyspecies_densplot, filename = "growbyspecies_densplott_quadXorigin.png", width = 8, height = 10)



# Seedling Growth

AGPE_seedgrow_densplot <- filterspecies_densplot(data = LTREB_grow_seedling, data_list = seed_grow_data_list, sim = y_seed_g_sim, 
                                             species = "AGPE", x = "Seedling Growth", y = "Density", xlim = 20)
ELRI_seedgrow_densplot <- filterspecies_densplot(data = LTREB_grow_seedling, data_list = seed_grow_data_list, sim = y_seed_g_sim, 
                                             species = "ELRI", x = "Seedling Growth", y = "Density", xlim = 10)
ELVI_seedgrow_densplot <- filterspecies_densplot(data = LTREB_grow_seedling, data_list = seed_grow_data_list, sim = y_seed_g_sim, 
                                             species = "ELVI", x = "Seedling Growth", y = "Density", xlim = 10)
FESU_seedgrow_densplot <- filterspecies_densplot(data = LTREB_grow_seedling, data_list = seed_grow_data_list, sim = y_seed_g_sim, 
                                             species = "FESU", x = "Seedling Growth", y = "Density", xlim = 20)
LOAR_seedgrow_densplot <- filterspecies_densplot(data = LTREB_grow_seedling, data_list = seed_grow_data_list, sim = y_seed_g_sim, 
                                             species = "LOAR", x = "Seedling Growth", y = "Density",  xlim = 30)
POAL_seedgrow_densplot <- filterspecies_densplot(data = LTREB_grow_seedling, data_list = seed_grow_data_list, sim = y_seed_g_sim, 
                                             species = "POAL", x = "Seedling Growth", y = "Density",  xlim = 20)
POSY_seedgrow_densplot <- filterspecies_densplot(data = LTREB_grow_seedling, data_list = seed_grow_data_list, sim = y_seed_g_sim, 
                                             species = "POSY", x = "Seedling Growth", y = "Density",  xlim = 20)


seedgrowbyspecies_densplot <- wrap_plots(A = (AGPE_seedgrow_densplot), B = (ELRI_seedgrow_densplot), C = (ELVI_seedgrow_densplot), D = (FESU_seedgrow_densplot), E = (LOAR_seedgrow_densplot), G = (POAL_seedgrow_densplot), H = (POSY_seedgrow_densplot),
                                     design = fm_layout, guides = "collect") + plot_annotation(title = "Seedling Growth", subtitle = "Vital rate fits and moments with 500 posterior draws")

ggsave(seedgrowbyspecies_densplot, filename = "seedgrowbyspecies_densplot.png", width = 8, height = 10)

# Fertility
AGPE_fert_densplot <- filterspecies_densplot(data = LTREB_data_forfert, data_list = fert_data_list, sim = y_fert_sim, 
                                            species = "AGPE", x = "No. of Reproductive Tillers", y = "Density")
ELRI_fert_densplot <- filterspecies_densplot(data = LTREB_data_forfert, data_list = fert_data_list, sim = y_fert_sim, 
                                            species = "ELRI", x = "No. of Reproductive Tillers", y = "Density")
ELVI_fert_densplot <- filterspecies_densplot(data = LTREB_data_forfert, data_list = fert_data_list, sim = y_fert_sim, 
                                            species = "ELVI", x = "No. of Reproductive Tillers", y = "Density")
FESU_fert_densplot <- filterspecies_densplot(data = LTREB_data_forfert, data_list = fert_data_list, sim = y_fert_sim, 
                                            species = "FESU", x = "No. of Reproductive Tillers", y = "Density")
LOAR_fert_densplot <- filterspecies_densplot(data = LTREB_data_forfert, data_list = fert_data_list, sim = y_fert_sim, 
                                            species = "LOAR", x = "No. of Reproductive Tillers", y = "Density")
POAL_fert_densplot <- filterspecies_densplot(data = LTREB_data_forfert, data_list = fert_data_list, sim = y_fert_sim, 
                                            species = "POAL", x = "No. of Reproductive Tillers", y = "Density")
POSY_fert_densplot <- filterspecies_densplot(data = LTREB_data_forfert, data_list = fert_data_list, sim = y_fert_sim, 
                                            species = "POSY", x = "No. of Reproductive Tillers", y = "Density")

fertbyspecies_densplot <- wrap_plots(A = (AGPE_fert_densplot), B = (ELRI_fert_densplot), C = (ELVI_fert_densplot), D = (FESU_fert_densplot), E = (LOAR_fert_densplot), G = (POAL_fert_densplot), H = (POSY_fert_densplot),
                                    design = fm_layout, guides = "collect") + plot_annotation(title = "Fertility", subtitle = "Vital rate fits and moments with 500 posterior draws")

ggsave(fertbyspecies_densplot, filename = "fertbyspecies_densplot_quadXorigin.png", width = 8, height = 10)


# Spikelet Count
AGPE_spike_densplot <- filterspecies_densplot(data = LTREB_data_forspike, data_list = spike_data_list, sim = y_spike_sim, 
                                             species = "AGPE", x = "Spikelets/Inflorescence", y = "Density")
ELRI_spike_densplot <- filterspecies_densplot(data = LTREB_data_forspike, data_list = spike_data_list, sim = y_spike_sim, 
                                             species = "ELRI", x = "Spikelets/Inflorescence", y = "Density")
ELVI_spike_densplot <- filterspecies_densplot(data = LTREB_data_forspike, data_list = spike_data_list, sim = y_spike_sim, 
                                             species = "ELVI", x = "Spikelets/Inflorescence", y = "Density")
FESU_spike_densplot <- filterspecies_densplot(data = LTREB_data_forspike, data_list = spike_data_list, sim = y_spike_sim, 
                                             species = "FESU", x = "Spikelets/Inflorescence", y = "Density")
LOAR_spike_densplot <- filterspecies_densplot(data = LTREB_data_forspike, data_list = spike_data_list, sim = y_spike_sim, 
                                             species = "LOAR", x = "Spikelets/Inflorescence", y = "Density")
POAL_spike_densplot <- filterspecies_densplot(data = LTREB_data_forspike, data_list = spike_data_list, sim = y_spike_sim, 
                                             species = "POAL", x = "Spikelets/Inflorescence", y = "Density")
POSY_spike_densplot <- filterspecies_densplot(data = LTREB_data_forspike, data_list = spike_data_list, sim = y_spike_sim, 
                                             species = "POSY", x = "Spikelets/Inflorescence", y = "Density")

spikebyspecies_densplot <- wrap_plots(A = (AGPE_spike_densplot), B = (ELRI_spike_densplot), C = (ELVI_spike_densplot), D = (FESU_spike_densplot), E = (LOAR_spike_densplot), G = (POAL_spike_densplot), H = (POSY_spike_densplot),
                                     design = fm_layout, guides = "collect") + plot_annotation(title = "Spikelet Production", subtitle = "Vital rate fits and moments with 500 posterior draws")

ggsave(spikebyspecies_densplot, filename = "spikebyspecies_densplot_quadXorigin.png", width = 8, height = 10)

# Spikelet Count
# Need to rename the outcome variable because I had used seed in the past, and this function, I'm using 'y' for every other vital rate.
seed_mean_data_list$y <- seed_mean_data_list$seed

AGPE_seedmean_densplot <- filterspecies_densplot(data = LTREB_data_for_seedmeans, data_list = seed_mean_data_list, sim = y_seedmean_sim, 
                                              species = "AGPE", x = "Seeds/Spikelets", y = "Density")
ELRI_seedmean_densplot <- filterspecies_densplot(data = LTREB_data_for_seedmeans, data_list = seed_mean_data_list, sim = y_seedmean_sim, 
                                              species = "ELRI", x = "Seeds/Spikelets", y = "Density")
ELVI_seedmean_densplot <- filterspecies_densplot(data = LTREB_data_for_seedmeans, data_list = seed_mean_data_list, sim = y_seedmean_sim, 
                                              species = "ELVI", x = "Seeds/Spikelets", y = "Density")
FESU_seedmean_densplot <- filterspecies_densplot(data = LTREB_data_for_seedmeans, data_list = seed_mean_data_list, sim = y_seedmean_sim, 
                                              species = "FESU", x = "Seeds/Spikelets", y = "Density")
LOAR_seedmean_densplot <- filterspecies_densplot(data = LTREB_data_for_seedmeans, data_list = seed_mean_data_list, sim = y_seedmean_sim, 
                                              species = "LOAR", x = "Seeds/Spikelets", y = "Density")
POAL_seedmean_densplot <- filterspecies_densplot(data = LTREB_data_for_seedmeans, data_list = seed_mean_data_list, sim = y_seedmean_sim, 
                                              species = "POAL", x = "Seeds/Spikelets", y = "Density")
POSY_seedmean_densplot <- filterspecies_densplot(data = LTREB_data_for_seedmeans, data_list = seed_mean_data_list, sim = y_seedmean_sim, 
                                              species = "POSY", x = "Seeds/Spikelets", y = "Density")

seedmeanbyspecies_densplot <- wrap_plots(A = (AGPE_seedmean_densplot), B = (ELRI_seedmean_densplot), C = (ELVI_seedmean_densplot), D = (FESU_seedmean_densplot), E = (LOAR_seedmean_densplot), G = (POAL_seedmean_densplot), H = (POSY_seedmean_densplot),
                                      design = fm_layout, guides = "collect") + plot_annotation(title = "Seed Production", subtitle = "Vital rate fits and moments with 500 posterior draws")

ggsave(seedmeanbyspecies_densplot, filename = "seedmeanbyspecies_densplot.png", width = 8, height = 10)


# Recruitment
# Need to create a dataframe version of the s_to_s list
rev_species_factor_key <- setNames(names(species_factor_key), species_factor_key)

s_to_s_df <- data.frame(s_to_s_data_list) %>% 
  mutate(species = as.character(recode(spp, !!!rev_species_factor_key)))

s_to_s_data_list$y <- s_to_s_data_list$tot_recruit_t1

AGPE_stos_densplot <- filterspecies_densplot(data = s_to_s_df, data_list = s_to_s_data_list, sim = y_recruit_sim, 
                                              species = "AGPE", x = "Successful Germination", y = "Density", xlim = 30)
ELRI_stos_densplot <- filterspecies_densplot(data = s_to_s_df, data_list = s_to_s_data_list, sim = y_recruit_sim, 
                                              species = "ELRI", x = "Successful Germination", y = "Density", xlim = 30)
ELVI_stos_densplot <- filterspecies_densplot(data = s_to_s_df, data_list = s_to_s_data_list, sim = y_recruit_sim, 
                                              species = "ELVI", x = "Successful Germination", y = "Density", xlim = 20)
FESU_stos_densplot <- filterspecies_densplot(data = s_to_s_df, data_list = s_to_s_data_list, sim = y_recruit_sim, 
                                              species = "FESU", x = "Successful Germination", y = "Density", xlim = 20)
LOAR_stos_densplot <- filterspecies_densplot(data = s_to_s_df, data_list = s_to_s_data_list, sim = y_recruit_sim, 
                                              species = "LOAR", x = "Successful Germination", y = "Density", xlim = 20)
POAL_stos_densplot <- filterspecies_densplot(data = s_to_s_df, data_list = s_to_s_data_list, sim = y_recruit_sim, 
                                              species = "POAL", x = "Successful Germination", y = "Density", xlim = 10)
POSY_stos_densplot <- filterspecies_densplot(data = s_to_s_df, data_list = s_to_s_data_list, sim = y_recruit_sim, 
                                              species = "POSY", x = "Successful Germination", y = "Density", xlim = 30)

stosbyspecies_densplot <- wrap_plots(A = (AGPE_stos_densplot), B = (ELRI_stos_densplot), C = (ELVI_stos_densplot), D = (FESU_stos_densplot), E = (LOAR_stos_densplot), G = (POAL_stos_densplot), H = (POSY_stos_densplot),
                                      design = fm_layout, guides = "collect") + plot_annotation(title = "Recruitment", subtitle = "Vital rate fits and moments with 500 posterior draws")

ggsave(stosbyspecies_densplot, filename = "stosbyspecies_densplot.png", width = 8, height = 10)



## Plot for size-specific moments for growth model minus seedlings, which have only one size ####
## could plot these for other vital rates if desired

# Function for looking at binned size_t fits, particularly important for the growth kernel as this determines the transitions through the matrix model
# plots the mean, sd, skew and kertosis of the posteriors (grey) as well as the mean of the posteriors for each moment (black) and the data (red) for size bins
size_moments_ppc <- function(data,y_name,sim, n_bins, title = NA){
  require(tidyverse)
  require(patchwork)
  data$y_name <- data[[y_name]]
  bins <- data %>%
    ungroup() %>% 
    arrange(logsize_t) %>% 
    mutate(size_bin = cut_number(logsize_t, n_bins)) %>% 
    group_by(size_bin)  %>% 
    dplyr::summarize(mean_t1 = mean(y_name),
                     sd_t1 = sd(y_name),
                     skew_t1 = skewness(y_name),
                     kurt_t1 = Lkurtosis(y_name),
                     bin_mean = mean(logsize_t),
                     bin_n = n())
  sim_moments <- bind_cols(enframe(data$logsize_t), as_tibble(t(sim))) %>%
    rename(logsize_t = value) %>%
    arrange(logsize_t) %>%
    mutate(size_bin = cut_number(logsize_t, n_bins)) %>%
    pivot_longer(., cols = starts_with("V"), names_to = "post_draw", values_to = "sim") %>%
    group_by(size_bin, post_draw) %>%
    summarize( Mean = mean((sim)),
               SD = sd((sim)),
               Skew = skewness((sim)),
               Kurtosis = Lkurtosis((sim)),
               bin_mean = mean(logsize_t),
               bin_n = n())
  sim_medians <- sim_moments %>%
    group_by(size_bin, bin_mean) %>%
    summarize(median_mean_sim = median(Mean),
              median_sd_sim = median(SD),
              median_skew_sim = median(Skew),
              median_kurt_sim = median(Kurtosis))
  meanplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = Mean), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_mean_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = mean_t1), shape = 1, color = "firebrick2") + xlab("log(size_t)") + theme_classic()
  sdplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = SD), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_sd_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = sd_t1), shape = 1, color = "firebrick2") + xlab("log(size_t)") + theme_classic()
  skewplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = Skew), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_skew_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = skew_t1), shape = 1, color = "firebrick2") + xlab("log(size_t)") + theme_classic()
  kurtplot <- ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = Kurtosis), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_kurt_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = kurt_t1), shape = 1, color = "firebrick2") + xlab("log(size_t)") + theme_classic()
  size_ppc_plot <- meanplot+ sdplot+skewplot+ kurtplot+plot_annotation(title = title)
  return(size_ppc_plot)
}



size_moments_ppc_by_species <- function(data,y_name,sim, n_bins, Species = NA, title = NA){
  require(tidyverse)
  require(patchwork)
  data$y_name <- data[[y_name]]
  bins <- data %>%
    ungroup() %>%
    arrange(logsize_t) %>%
    mutate(size_bin = cut_number(logsize_t, n_bins)) %>%
    group_by(species, size_bin)  %>%
    dplyr::summarize(mean_t1 = mean(y_name),
                     sd_t1 = sd(y_name),
                     skew_t1 = skewness(y_name),
                     kurt_t1 = Lkurtosis(y_name),
                     bin_mean = mean(logsize_t),
                     bin_n = n()) 
  sim_moments <- bind_cols(species = data$species, logsize_t = data$logsize_t, as_tibble(t(sim))) %>%
    arrange(logsize_t) %>%
    mutate(size_bin = cut_number(logsize_t, n_bins)) %>%
    pivot_longer(., cols = starts_with("V"), names_to = "post_draw", values_to = "sim") %>%
    group_by(species, size_bin, post_draw) %>%
    summarize( Mean = mean((sim)),
               SD = sd((sim)),
               Skew = skewness((sim)),
               Kurtosis = Lkurtosis((sim)),
               bin_mean = mean(logsize_t),
               bin_n = n()) 
  sim_medians <- sim_moments %>%
    group_by(species, size_bin, bin_mean) %>%
    summarize(median_mean_sim = median(Mean),
              median_sd_sim = median(SD),
              median_skew_sim = median(Skew),
              median_kurt_sim = median(Kurtosis)) 
  bins <- bins %>% filter(species == Species)
  sim_moments <- sim_moments %>% filter(species == Species)
  sim_medians <- sim_medians %>% filter(species == Species)
  
  meanplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = Mean), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_mean_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = mean_t1), shape = 1, color = "firebrick2") + xlab("log(size_t)") + theme_classic()
  sdplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = SD), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_sd_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = sd_t1), shape = 1, color = "firebrick2") + xlab("log(size_t)") + theme_classic()
  skewplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = Skew), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_skew_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = skew_t1), shape = 1, color = "firebrick2") + xlab("log(size_t)") + theme_classic()
  kurtplot <- ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = Kurtosis), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_kurt_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = kurt_t1), shape = 1, color = "firebrick2") + xlab("log(size_t)") + theme_classic()
  size_ppc_plot <- meanplot+ sdplot+skewplot+ kurtplot+plot_annotation(title = title, subtitle = Species)
  return(size_ppc_plot)
  # return(list(meanplot = meanplot, sdplot = sdplot, skewplot = skewplot, kurtplot = kurtplot, size_ppc_plot))
}



##### Growth ####
# all together

PIG_growth_size_ppc <- size_moments_ppc(data = LTREB_data_forgrow,
                                        y_name = "size_t1",
                                        sim = y_g_sim, 
                                        n_bins = 5, 
                                        title = "Adult Growth")
# PIG_growth_size_ppc
# ggsave(PIG_growth_size_ppc, filename = "PIG_growth_size_ppc_quadXorigin.png", width = 4, height = 4)

# by species
AGPE_growth_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forgrow,
                                                    y_name = "size_t1",
                                                    sim = y_g_sim,
                                                    n_bins = 5,
                                                    title = "Adult Growth",
                                                    Species = "AGPE")
# AGPE_growth_size_ppc
ggsave(AGPE_growth_size_ppc, filename = "AGPE_growth_size_ppc_quadXorigin.png", width = 4, height = 4)

ELRI_growth_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forgrow,
                                                    y_name = "size_t1",
                                                    sim = y_g_sim,
                                                    n_bins = 5,
                                                    title = "Adult Growth",
                                                    Species = "ELRI")
# ELRI_growth_size_ppc
ggsave(ELRI_growth_size_ppc, filename = "ELRI_growth_size_ppc_quadXorigin.png", width = 4, height = 4)

ELVI_growth_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forgrow,
                                                    y_name = "size_t1",
                                                    sim = y_g_sim,
                                                    n_bins = 5,
                                                    title = "Adult Growth",
                                                    Species = "ELVI")
# ELVI_growth_size_ppc
ggsave(ELVI_growth_size_ppc, filename = "ELVI_growth_size_ppc_quadXorigin.png", width = 4, height = 4)

FESU_growth_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forgrow,
                                                    y_name = "size_t1",
                                                    sim = y_g_sim,
                                                    n_bins = 5,
                                                    title = "Adult Growth",
                                                    Species = "FESU")
# FESU_growth_size_ppc
ggsave(FESU_growth_size_ppc, filename = "FESU_growth_size_ppc_quadXorigin.png", width = 4, height = 4)

LOAR_growth_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forgrow,
                                                    y_name = "size_t1",
                                                    sim = y_g_sim,
                                                    n_bins = 5,
                                                    title = "Adult Growth",
                                                    Species = "LOAR")
# LOAR_growth_size_ppc
ggsave(LOAR_growth_size_ppc, filename = "LOAR_growth_size_ppc_quadXorigin.png", width = 4, height = 4)

POAL_growth_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forgrow,
                                                    y_name = "size_t1",
                                                    sim = y_g_sim,
                                                    n_bins = 5,
                                                    title = "Adult Growth",
                                                    Species = "POAL")
# POAL_growth_size_ppc
ggsave(POAL_growth_size_ppc, filename = "POAL_growth_size_ppc_quadXorigin.png", width = 4, height = 4)

POSY_growth_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forgrow,
                                                    y_name = "size_t1",
                                                    sim = y_g_sim,
                                                    n_bins = 5,
                                                    title = "Adult Growth",
                                                    Species = "POSY")
# POSY_growth_size_ppc
ggsave(POSY_growth_size_ppc, filename = "POSY_growth_size_ppc_quadXorigin.png", width = 4, height = 4)





##### adult survival ####
# all together

surv_size_ppc <- size_moments_ppc(data = LTREB_data_forsurv,
                                  y_name = "surv_t1",
                                  sim = y_s_sim, 
                                  n_bins = 4, 
                                  title = "Adult Survival")
# surv_size_ppc
# ggsave(surv_size_ppc, filename = "surv_size_ppc_quadXorigin.png", width = 4, height = 4)

# by species
AGPE_surv_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forsurv,
                                                  y_name = "surv_t1",
                                                  sim = y_s_sim,
                                                  n_bins = 4,
                                                  title = "Adult Survival",
                                                  Species = "AGPE")
# AGPE_surv_size_ppc
ggsave(AGPE_surv_size_ppc, filename = "AGPE_surv_size_ppc_quadXorigin.png", width = 4, height = 4)

ELRI_surv_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forsurv,
                                                  y_name = "surv_t1",
                                                  sim = y_s_sim,
                                                  n_bins = 4,
                                                  title = "Adult Survival",
                                                  Species = "ELRI")
# ELRI_surv_size_ppc
ggsave(ELRI_surv_size_ppc, filename = "ELRI_surv_size_ppc_quadXorigin.png", width = 4, height = 4)

ELVI_surv_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forsurv,
                                                  y_name = "surv_t1",
                                                  sim = y_s_sim,
                                                  n_bins = 4,
                                                  title = "Adult Survival",
                                                  Species = "ELVI")
# ELVI_surv_size_ppc
ggsave(ELVI_surv_size_ppc, filename = "ELVI_surv_size_ppc_quadXorigin.png", width = 4, height = 4)

FESU_surv_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forsurv,
                                                  y_name = "surv_t1",
                                                  sim = y_s_sim,
                                                  n_bins = 4,
                                                  title = "Adult Survival",
                                                  Species = "FESU")
# FESU_surv_size_ppc
ggsave(FESU_surv_size_ppc, filename = "FESU_surv_size_ppc_quadXorigin.png", width = 4, height = 4)

LOAR_surv_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forsurv,
                                                  y_name = "surv_t1",
                                                  sim = y_s_sim,
                                                  n_bins = 4,
                                                  title = "Adult Survival",
                                                  Species = "LOAR")
# LOAR_surv_size_ppc
ggsave(LOAR_surv_size_ppc, filename = "LOAR_surv_size_ppc_quadXorigin.png", width = 4, height = 4)

POAL_surv_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forsurv,
                                                  y_name = "surv_t1",
                                                  sim = y_s_sim,
                                                  n_bins = 4,
                                                  title = "Adult Survival",
                                                  Species = "POAL")
# POAL_surv_size_ppc
ggsave(POAL_surv_size_ppc, filename = "POAL_surv_size_ppc_quadXorigin.png", width = 4, height = 4)

POSY_surv_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forsurv,
                                                  y_name = "surv_t1",
                                                  sim = y_s_sim,
                                                  n_bins = 4,
                                                  title = "Adult Survival",
                                                  Species = "POSY")
# POSY_surv_size_ppc
ggsave(POSY_surv_size_ppc, filename = "POSY_surv_size_ppc_quadXorigin.png", width = 4, height = 4)




##### Flowering####
# all together

flw_size_ppc <- size_moments_ppc(data = LTREB_data_forflw,
                                 y_name = "FLW_STAT_T1",
                                 sim = y_f_sim, 
                                 n_bins = 5, 
                                 title = "Flowering")
# surv_size_ppc
# ggsave(surv_size_ppc, filename = "surv_size_ppc_quadXorigin.png", width = 4, height = 4)

# by species
AGPE_flw_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forflw,
                                                 y_name = "FLW_STAT_T1",
                                                 sim = y_f_sim,
                                                 n_bins = 2,
                                                 title = "Flowering",
                                                 Species = "AGPE")
# AGPE_flw_size_ppc
ggsave(AGPE_flw_size_ppc, filename = "AGPE_flw_size_ppc_quadXorigin.png", width = 4, height = 4)

ELRI_flw_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forflw,
                                                 y_name = "FLW_STAT_T1",
                                                 sim = y_f_sim,
                                                 n_bins = 2,
                                                 title = "Flowering",
                                                 Species = "ELRI")
# ELRI_flw_size_ppc
ggsave(ELRI_flw_size_ppc, filename = "ELRI_flw_size_ppc_quadXorigin.png", width = 4, height = 4)

ELVI_flw_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forflw,
                                                 y_name = "FLW_STAT_T1",
                                                 sim = y_f_sim,
                                                 n_bins = 2,
                                                 title = "Flowering",
                                                 Species = "ELVI")
# ELVI_flw_size_ppc
ggsave(ELVI_flw_size_ppc, filename = "ELVI_flw_size_ppc_quadXorigin.png", width = 4, height = 4)

FESU_flw_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forflw,
                                                 y_name = "FLW_STAT_T1",
                                                 sim = y_f_sim,
                                                 n_bins = 2,
                                                 title = "Flowering",
                                                 Species = "FESU")
# FESU_flw_size_ppc
ggsave(FESU_flw_size_ppc, filename = "FESU_flw_size_ppc_quadXorigin.png", width = 4, height = 4)

LOAR_flw_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forflw,
                                                 y_name = "FLW_STAT_T1",
                                                 sim = y_f_sim,
                                                 n_bins = 2,
                                                 title = "Flowering",
                                                 Species = "LOAR")
# LOAR_flw_size_ppc
ggsave(LOAR_flw_size_ppc, filename = "LOAR_flw_size_ppc_quadXorigin.png", width = 4, height = 4)

POAL_flw_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forflw,
                                                 y_name = "FLW_STAT_T1",
                                                 sim = y_f_sim,
                                                 n_bins = 2,
                                                 title = "Flowering",
                                                 Species = "POAL")
# POAL_flw_size_ppc
ggsave(POAL_flw_size_ppc, filename = "POAL_flw_size_ppc_quadXorigin.png", width = 4, height = 4)

POSY_flw_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forflw,
                                                 y_name = "FLW_STAT_T1",
                                                 sim = y_f_sim,
                                                 n_bins = 2,
                                                 title = "Flowering",
                                                 Species = "POSY")
# POSY_flw_size_ppc
ggsave(POSY_flw_size_ppc, filename = "POSY_flw_size_ppc_quadXorigin.png", width = 4, height = 4)






##### Inflorescence production ####
# all together

fert_size_ppc <- size_moments_ppc(data = LTREB_data_forfert,
                                  y_name = "FLW_COUNT_T1",
                                  sim = y_fert_sim, 
                                  n_bins = 5, 
                                  title = "Infl. Production")
# fert_size_ppc
# ggsave(fert_size_ppc, filename = "fert_size_ppc_quadXorigin.png", width = 4, height = 4)

# by species
AGPE_fert_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forfert,
                                                  y_name = "FLW_COUNT_T1",
                                                  sim = y_fert_sim,
                                                  n_bins = 5,
                                                  title = "Infl. Production",
                                                  Species = "AGPE")
# AGPE_fert_size_ppc
ggsave(AGPE_fert_size_ppc, filename = "AGPE_fert_size_ppc_quadXorigin.png", width = 4, height = 4)

ELRI_fert_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forfert,
                                                  y_name = "FLW_COUNT_T1",
                                                  sim = y_fert_sim,
                                                  n_bins = 5,
                                                  title = "Infl. Production",
                                                  Species = "ELRI")
# ELRI_fert_size_ppc
ggsave(ELRI_fert_size_ppc, filename = "ELRI_fert_size_ppc_quadXorigin.png", width = 4, height = 4)

ELVI_fert_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forfert,
                                                  y_name = "FLW_COUNT_T1",
                                                  sim = y_fert_sim,
                                                  n_bins = 5,
                                                  title = "Infl. Production",
                                                  Species = "ELVI")
# ELVI_fert_size_ppc
ggsave(ELVI_fert_size_ppc, filename = "ELVI_fert_size_ppc_quadXorigin.png", width = 4, height = 4)

FESU_fert_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forfert,
                                                  y_name = "FLW_COUNT_T1",
                                                  sim = y_fert_sim,
                                                  n_bins = 5,
                                                  title = "Infl. Production",
                                                  Species = "FESU")
# FESU_fert_size_ppc
ggsave(FESU_fert_size_ppc, filename = "FESU_fert_size_ppc_quadXorigin.png", width = 4, height = 4)

LOAR_fert_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forfert,
                                                  y_name = "FLW_COUNT_T1",
                                                  sim = y_fert_sim,
                                                  n_bins = 5,
                                                  title = "Infl. Production",
                                                  Species = "LOAR")
# LOAR_fert_size_ppc
ggsave(LOAR_fert_size_ppc, filename = "LOAR_fert_size_ppc_quadXorigin.png", width = 4, height = 4)

POAL_fert_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forfert,
                                                  y_name = "FLW_COUNT_T1",
                                                  sim = y_fert_sim,
                                                  n_bins = 5,
                                                  title = "Infl. Production",
                                                  Species = "POAL")
# POAL_fert_size_ppc
ggsave(POAL_fert_size_ppc, filename = "POAL_fert_size_ppc_quadXorigin.png", width = 4, height = 4)

POSY_fert_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forfert,
                                                  y_name = "FLW_COUNT_T1",
                                                  sim = y_fert_sim,
                                                  n_bins = 5,
                                                  title = "Infl. Production",
                                                  Species = "POSY")
# POSY_fert_size_ppc
ggsave(POSY_fert_size_ppc, filename = "POSY_fert_size_ppc_quadXorigin.png", width = 4, height = 4)





##### Inflorescence production ####
# all together

spike_size_ppc <- size_moments_ppc(data = LTREB_data_forspike,
                                   y_name = "spike_count_t1",
                                   sim = y_spike_sim, 
                                   n_bins = 5, 
                                   title = "Spikelets/Infl.")
# spike_size_ppc
# ggsave(fert_size_ppc, filename = "fert_size_ppc_quadXorigin.png", width = 4, height = 4)

# by species
AGPE_spike_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forspike,
                                                   y_name = "spike_count_t1",
                                                   sim = y_spike_sim,
                                                   n_bins = 5,
                                                   title = "Spikelets/Infl.",
                                                   Species = "AGPE")
# AGPE_spike_size_ppc
ggsave(AGPE_spike_size_ppc, filename = "AGPE_spike_size_ppc_quadXorigin.png", width = 4, height = 4)

ELRI_spike_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forspike,
                                                   y_name = "spike_count_t1",
                                                   sim = y_spike_sim,
                                                   n_bins = 5,
                                                   title = "Spikelets/Infl.",
                                                   Species = "ELRI")
# ELRI_spike_size_ppc
ggsave(ELRI_spike_size_ppc, filename = "ELRI_spike_size_ppc_quadXorigin.png", width = 4, height = 4)

ELVI_spike_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forspike,
                                                   y_name = "spike_count_t1",
                                                   sim = y_spike_sim,
                                                   n_bins = 5,
                                                   title = "Spikelets/Infl.",
                                                   Species = "ELVI")
# ELVI_spike_size_ppc
ggsave(ELVI_spike_size_ppc, filename = "ELVI_spike_size_ppc_quadXorigin.png", width = 4, height = 4)

FESU_spike_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forspike,
                                                   y_name = "spike_count_t1",
                                                   sim = y_spike_sim,
                                                   n_bins = 5,
                                                   title = "Spikelets/Infl.",
                                                   Species = "FESU")
# FESU_spike_size_ppc
ggsave(FESU_spike_size_ppc, filename = "FESU_spike_size_ppc_quadXorigin.png", width = 4, height = 4)

LOAR_spike_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forspike,
                                                   y_name = "spike_count_t1",
                                                   sim = y_spike_sim,
                                                   n_bins = 5,
                                                   title = "Spikelets/Infl.",
                                                   Species = "LOAR")
# LOAR_spike_size_ppc
ggsave(LOAR_spike_size_ppc, filename = "LOAR_spike_size_ppc_quadXorigin.png", width = 4, height = 4)

POAL_spike_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forspike,
                                                   y_name = "spike_count_t1",
                                                   sim = y_spike_sim,
                                                   n_bins = 5,
                                                   title = "Spikelets/Infl.",
                                                   Species = "POAL")
# POAL_spike_size_ppc
ggsave(POAL_spike_size_ppc, filename = "POAL_spike_size_ppc_quadXorigin.png", width = 4, height = 4)

POSY_spike_size_ppc <- size_moments_ppc_by_species(data = LTREB_data_forspike,
                                                   y_name = "spike_count_t1",
                                                   sim = y_spike_sim,
                                                   n_bins = 5,
                                                   title = "Spikelets/Infl.",
                                                   Species = "POSY")
# POSY_spike_size_ppc
ggsave(POSY_spike_size_ppc, filename = "POSY_spike_size_ppc_quadXorigin.png", width = 4, height = 4)


size_ppc_layout <- "
AB
CD
EH
H#
"


size_ppc_plot <- wrap_plots(A = wrap_elements(AGPE_spike_size_ppc), B =  wrap_elements(ELRI_spike_size_ppc),
                            C = wrap_elements(ELVI_spike_size_ppc), D = wrap_elements(FESU_spike_size_ppc),
                            E = wrap_elements(LOAR_spike_size_ppc), G = wrap_elements(POAL_spike_size_ppc), 
                            H = wrap_elements(POSY_spike_size_ppc), I = plot_spacer(), design = size_ppc_layout) + plot_annotation(title = "Size specific vital rate moments")






# a plot putting all the cross species plots together

size_ppc_layout <- "
AB
CD
E#
"


size_ppc_plot <- wrap_plots(A = wrap_elements(PIG_growth_size_ppc), B =  wrap_elements(surv_size_ppc),
                            C = wrap_elements(flw_size_ppc), D = wrap_elements(PIG_fert_size_ppc),
                            E = wrap_elements(PIG_spike_size_ppc), design = size_ppc_layout) + plot_annotation(title = "Size specific vital rate moments")


# size_ppc_plot
ggsave(size_ppc_plot, filename = "size_ppc_plott_quad_origin.png", width = 10, height = 13)




## Traceplots for select parameters from all vital rates for AGPE
surv_trace <- mcmc_trace(surv_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Survival")
seedsurv_trace <- mcmc_trace(surv_fit_seedling, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Seedling Survival")
flw_trace <- mcmc_trace(flw_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Flowering")
grow_trace <- mcmc_trace(grow_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Growth")
seedgrow_trace <- mcmc_trace(grow_fit_seedling, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Seedling Growth")
fert_trace <- mcmc_trace(fert_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Fertility")
spike_trace <- mcmc_trace(spike_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Spikelets per inflorescence")
seedmean_trace <- mcmc_trace(seedmean_fit, pars = c("beta0[1]", "betaendo[1]"))+ggtitle("Mean Seeds per spikelet")
stos_trace <- mcmc_trace(stos_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Recruitment")

vr_traceplots <- (surv_trace)/
                 (seedsurv_trace)/
                 (grow_trace)/
                 (seedgrow_trace)/
                 (flw_trace)/
                 (fert_trace)/
                 (spike_trace)/
                 (seedmean_trace)/
                 (stos_trace) + plot_annotation(title = "Traceplots for select parameters from all vital rates (AGPE)")
ggsave(vr_traceplots, filename = "vr_traceplots.png", width = 25, height = 20)


## Plots for mean and year variance endophyte effects on all vital rates, all species ####
sigmayear_species_key <- c("sigma_year[1,1]" = "AGPE E-", "sigma_year[1,2]" = "AGPE E+", "sigma_year[2,1]" = "ELRI E-",  "sigma_year[2,2]" = "ELRI E+", "sigma_year[3,1]" = "ELVI E-", "sigma_year[3,2]" = "ELVI E+", "sigma_year[4,1]" = "FESU E-","sigma_year[4,2]" = "FESU E+", "sigma_year[5,1]" = "LOAR E-", "sigma_year[5,2]" = "LOAR E+", "sigma_year[6,1]" = "POAL E-", "sigma_year[6,2]" = "POAL E+", "sigma_year[7,1]" = "POSY E-","sigma_year[7,2]" = "POSY E+")
mean_species_key <- c("betaendo[1]" = "AGPE", "betaendo[2]" = "ELRI", "betaendo[3]" = "ELVI", "betaendo[4]" = "FESU", "betaendo[5]" = "LOAR", "betaendo[6]" = "POAL", "betaendo[7]" = "POSY")
sd_species_key <- c("sigmaendo[1]" = "AGPE", "sigmaendo[2]" = "ELRI", "sigmaendo[3]" = "ELVI", "sigmaendo[4]" = "FESU", "sigmaendo[5]" = "LOAR", "sigmaendo[6]" = "POAL", "sigmaendo[7]" = "POSY")

#posteriors of variance for all vital rates
surv_sigmayear_posteriors <-  mcmc_areas(surv_fit, prob = 0.8, regex_pars = c("sigma_year"))+labs(title = "Adult Survival", subtitle = "Vital rate SD with 80% credible intervals") + scale_y_discrete(labels = sigmayear_species_key)
seedsurv_sigmayear_posteriors <-  mcmc_areas(surv_fit_seedling, prob = 0.8, regex_pars = c("sigma_year"))+labs(title = "Seedling Survival", subtitle = "Vital rate SD with 80% credible intervals") + scale_y_discrete(labels = sigmayear_species_key)
grow_sigmayear_posteriors <- mcmc_areas(grow_fit, prob = 0.8, regex_pars = c("sigma_year"))+labs(title = "Adult Growth", subtitle = "Vital rate SD with 80% credible intervals") + scale_y_discrete(labels = sigmayear_species_key)
seedgrow_sigmayear_posteriors <- mcmc_areas(grow_fit_seedling, prob = 0.8, regex_pars = c("sigma_year"))+labs(title = "Seedling Growth", subtitle = "Vital rate SD with 80% credible intervals") + scale_y_discrete(labels = sigmayear_species_key)
flw_sigmayear_posteriors <- mcmc_areas(flw_fit, prob = 0.8, regex_pars = c("sigma_year"))+labs(title = "Flowering", subtitle = "Vital rate SD with 80% credible intervals") + scale_y_discrete(labels = sigmayear_species_key)
fert_sigmayear_posteriors <- mcmc_areas(fert_fit, prob = 0.8, regex_pars = c("sigma_year"))+labs(title = "Fertility", subtitle = "Vital rate SD with 80% credible intervals") + scale_y_discrete(labels = sigmayear_species_key)
spike_sigmayear_posteriors <- mcmc_areas(spike_fit, prob = 0.8, regex_pars = c("sigma_year"))+labs(title = "Spikelets per infl.", subtitle = "Vital rate SD with 80% credible intervals") + scale_y_discrete(labels = sigmayear_species_key)
stos_sigmayear_posteriors <- mcmc_areas(stos_fit, prob = 0.8, regex_pars = c("sigma_year"))+labs(title = "Germination", subtitle = "Vital rate SD with 80% credible intervals") + scale_y_discrete(labels = sigmayear_species_key)

vital_rate_sd_posteriors <- surv_sigmayear_posteriors +seedsurv_sigmayear_posteriors +grow_sigmayear_posteriors + seedgrow_sigmayear_posteriors +flw_sigmayear_posteriors + fert_sigmayear_posteriors + spike_sigmayear_posteriors +stos_sigmayear_posteriors+
  plot_layout(ncol = 1)
ggsave(vital_rate_sd_posteriors, filename = "vital_rate_sd_posteriors.png", width = 12, height = 30)


# effect of endophyte on mean (betaendo)
surv_endomean_posteriors <-  mcmc_areas(surv_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Adult Survival", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
seedsurv_endomean_posteriors <-  mcmc_areas(surv_fit_seedling, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Seedling Survival", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
grow_endomean_posteriors <- mcmc_areas(grow_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Adult Growth", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
seedgrow_endomean_posteriors <- mcmc_areas(grow_fit_seedling, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Seedling Growth", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
flw_endomean_posteriors <- mcmc_areas(flw_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Flowering", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
fert_endomean_posteriors <- mcmc_areas(fert_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Fertility", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
spike_endomean_posteriors <- mcmc_areas(spike_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Spikelets per infl.", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
seedmean_endomean_posteriors <- mcmc_areas(seedmean_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Mean seeds per spikelet", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
stos_endomean_posteriors <- mcmc_areas(stos_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Germination", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)

# effect of endophyte on year standard deviation (seed mean does not have an endophyte effect on variance)
surv_endosd_posteriors <-  mcmc_areas(surv_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Adult Survival", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
seedsurv_endosd_posteriors <-  mcmc_areas(surv_fit_seedling, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Seedling Survival", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
grow_endosd_posteriors <- mcmc_areas(grow_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Adult Growth", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
seedgrow_endosd_posteriors <- mcmc_areas(grow_fit_seedling, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Seedling Growth", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
flw_endosd_posteriors <- mcmc_areas(flw_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Flowering", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
fert_endosd_posteriors <- mcmc_areas(fert_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Fertility", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
spike_endosd_posteriors <- mcmc_areas(spike_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Spikelets per infl.", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)
stos_endosd_posteriors <- mcmc_areas(stos_fit, prob = 0.8, regex_pars = c("sigmaendo"))+labs(title = "Germination", subtitle = "Endophyte effect on SD with 80% credible intervals") + scale_y_discrete(labels = sd_species_key)



endomeanandvar_posteriors <- (surv_endomean_posteriors + surv_endosd_posteriors)/
                       (seedsurv_endomean_posteriors + seedsurv_endosd_posteriors)/
                       (grow_endomean_posteriors + grow_endosd_posteriors)/
                       (seedgrow_endomean_posteriors + seedgrow_endosd_posteriors)/
                       (flw_endomean_posteriors + flw_endosd_posteriors)/
                       (fert_endomean_posteriors + fert_endosd_posteriors)/
                       (spike_endomean_posteriors + spike_endosd_posteriors)/
                       (seedmean_endomean_posteriors + plot_spacer())/
                       (stos_endomean_posteriors + stos_endosd_posteriors) + plot_annotation(title = "Endophyte effect on mean and interannual SD across vital rates by species")
  
ggsave(endomeanandvar_posteriors, filename = "endomeanandvar_posteriors.png", width = 12, height = 30)

# making plots with the histogram of E+ and E- variance for each vital rate
#surv
dimnames(surv_par$sigma_year) <- list(Draw = paste0("i",1:dim(surv_par$sigma_year)[1]), Species = species_list, Endo = c("S-", "S+"))
surv_sigmayear_cube <- cubelyr::as.tbl_cube(surv_par$sigma_year)
surv_sigmayear_df <- as_tibble(surv_sigmayear_cube)  %>% 
  rename(sigma_year = `surv_par$sigma_year`)

surv_sigmayear_hist <- ggplot(data = surv_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( twotone_endophyte_color_scheme[1], twotone_endophyte_color_scheme[2]))+
  facet_wrap(~Species, scales = "free", ncol = 4)+
  labs(title = "Survival", x = "Interannual SD")+
  theme_minimal()+
  theme(strip.text = element_text(face = "italic"))
# surv_sigmayear_hist
ggsave(surv_sigmayear_hist, filename = "surv_sigmayear_hist_quadXorigin.png", width = 7, height = 3)


#seedling surv
dimnames(surv_sdlg_par$sigma_year) <- list(Draw = paste0("i",1:dim(surv_sdlg_par$sigma_year)[1]), Species = species_list, Endo = c("S-", "S+"))
seedsurv_sigmayear_cube <- cubelyr::as.tbl_cube(surv_sdlg_par$sigma_year)
seedsurv_sigmayear_df <- as_tibble(seedsurv_sigmayear_cube)  %>% 
  rename(sigma_year = `surv_sdlg_par$sigma_year`)

seedsurv_sigmayear_hist <- ggplot(data = seedsurv_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( twotone_endophyte_color_scheme[1], twotone_endophyte_color_scheme[2]))+
  facet_wrap(~Species, scales = "free", ncol = 4)+
  labs(title = "Seedling Survival", x = "Interannual SD")+
  theme_minimal()+
  theme(strip.text = element_text(face = "italic"))
# seedsurv_sigmayear_hist
ggsave(seedsurv_sigmayear_hist, filename = "seedsurv_sigmayear_hist.png", width = 7, height = 3)


#grow
dimnames(grow_par$sigma_year) <- list(Draw = paste0("i",1:dim(grow_par$sigma_year)[1]), Species = species_list, Endo = c("S-", "S+"))
grow_sigmayear_cube <- cubelyr::as.tbl_cube(grow_par$sigma_year)
grow_sigmayear_df <- as_tibble(grow_sigmayear_cube)  %>% 
  rename(sigma_year = `grow_par$sigma_year`)

grow_sigmayear_hist <- ggplot(data = grow_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( twotone_endophyte_color_scheme[1], twotone_endophyte_color_scheme[2]))+
  facet_wrap(~Species, scales = "free", ncol = 4)+
  labs(title = "Growth", x = "Interannual SD")+
  theme_minimal()+
  theme(strip.text = element_text(face = "italic"))
# grow_sigmayear_hist
ggsave(grow_sigmayear_hist, filename = "grow_sigmayear_hist_quadXorigin.png", width = 7, height = 3)


#seedling grow
dimnames(grow_sdlg_par$sigma_year) <- list(Draw = paste0("i",1:dim(grow_sdlg_par$sigma_year)[1]), Species = species_list, Endo = c("S-", "S+"))
seedgrow_sigmayear_cube <- cubelyr::as.tbl_cube(grow_sdlg_par$sigma_year)
seedgrow_sigmayear_df <- as_tibble(seedgrow_sigmayear_cube)  %>% 
  rename(sigma_year = `grow_sdlg_par$sigma_year`)

seedgrow_sigmayear_hist <- ggplot(data = seedgrow_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( twotone_endophyte_color_scheme[1], twotone_endophyte_color_scheme[2]))+
  facet_wrap(~Species, scales = "free", ncol = 4)+
  labs(title = "Seeding Growth", x = "Interannual SD")+
  theme_minimal()+
  theme(strip.text = element_text(face = "italic"))
# seedgrow_sigmayear_hist
ggsave(seedgrow_sigmayear_hist, filename = "seedgrow_sigmayear_hist.png", width = 7, height = 3)


#flw
dimnames(flow_par$sigma_year) <- list(Draw = paste0("i",1:dim(flow_par$sigma_year)[1]), Species = species_list, Endo = c("S-", "S+"))
flow_sigmayear_cube <- cubelyr::as.tbl_cube(flow_par$sigma_year)
flow_sigmayear_df <- as_tibble(flow_sigmayear_cube)  %>% 
  rename(sigma_year = `flow_par$sigma_year`)

flow_sigmayear_hist <- ggplot(data = flow_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( twotone_endophyte_color_scheme[1], twotone_endophyte_color_scheme[2]))+
  facet_wrap(~Species, scales = "free", ncol = 4)+
  labs(title = "Flowering", x = "Interannual SD")+
  theme_minimal()+
  theme(strip.text = element_text(face = "italic"))
# flow_sigmayear_hist
ggsave(flow_sigmayear_hist, filename = "flow_sigmayear_hist_quadXorigin.png", width = 7, height = 3)


#fert
dimnames(fert_par$sigma_year) <- list(Draw = paste0("i",1:dim(fert_par$sigma_year)[1]), Species = species_list, Endo = c("S-", "S+"))
fert_sigmayear_cube <- cubelyr::as.tbl_cube(fert_par$sigma_year)
fert_sigmayear_df <- as_tibble(fert_sigmayear_cube)  %>% 
  rename(sigma_year = `fert_par$sigma_year`)

fert_sigmayear_hist <- ggplot(data = fert_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( twotone_endophyte_color_scheme[1], twotone_endophyte_color_scheme[2]))+
  facet_wrap(~Species, scales = "free", ncol = 4)+
  labs(title = "# of Flw Tillers", x = "Interannual SD")+
  theme_minimal()+
  theme(strip.text = element_text(face = "italic"))
# fert_sigmayear_hist
ggsave(fert_sigmayear_hist, filename = "fert_sigmayear_hist_quadXorigin.png", width = 7, height = 3)


#spike
dimnames(spike_par$sigma_year) <- list(Draw = paste0("i",1:dim(spike_par$sigma_year)[1]), Species = species_list, Endo = c("S-", "S+"))
spike_sigmayear_cube <- cubelyr::as.tbl_cube(spike_par$sigma_year)
spike_sigmayear_df <- as_tibble(spike_sigmayear_cube)  %>% 
  rename(sigma_year = `spike_par$sigma_year`)

spike_sigmayear_hist <- ggplot(data = spike_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( twotone_endophyte_color_scheme[1], twotone_endophyte_color_scheme[2]))+
  facet_wrap(~Species, scales = "free", ncol = 4)+
  labs(title = "Spikelets", x = "Interannual SD")+
  theme_minimal()+
  theme(strip.text = element_text(face = "italic"))
# spike_sigmayear_hist
ggsave(spike_sigmayear_hist, filename = "spike_sigmayear_hist_quadXorigin.png", width = 7, height = 3)


#germination
dimnames(recruit_par$sigma_year) <- list(Draw = paste0("i",1:dim(recruit_par$sigma_year)[1]), Species = species_list, Endo = c("S-", "S+"))
recruit_sigmayear_cube <- cubelyr::as.tbl_cube(recruit_par$sigma_year)
recruit_sigmayear_df <- as_tibble(recruit_sigmayear_cube)  %>% 
  rename(sigma_year = `recruit_par$sigma_year`)

recruit_sigmayear_hist <- ggplot(data = recruit_sigmayear_df)+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( twotone_endophyte_color_scheme[1], twotone_endophyte_color_scheme[2]))+
  facet_wrap(~Species, scales = "free", ncol = 4)+
  labs(title = "Recruitment", x = "Interannual SD")+
  theme_minimal()+
  theme(strip.text = element_text(face = "italic"))
# recruit_sigmayear_hist
ggsave(recruit_sigmayear_hist, filename = "recruit_sigmayear_hist.png", width = 7, height = 3)

# posterior histograms of E+ and E- variance

endo_sigmayear_histograms <- surv_sigmayear_hist + seedsurv_sigmayear_hist + grow_sigmayear_hist + seedgrow_sigmayear_hist + flow_sigmayear_hist + fert_sigmayear_hist + spike_sigmayear_hist + recruit_sigmayear_hist+
  plot_layout( nrow  = 2, guides = "collect")+
  plot_annotation(title = "Endophyte effect on interannual SD across vital rates by species")

ggsave(endo_sigmayear_histograms, filename = "endo_sigmayear_histograms_quadXorigin.png", width = 20, height = 15)


####### Making pop out histogram panels for figure 2 (C and D)

FESUsurv_sigmayear_hist <- ggplot(data = filter(surv_sigmayear_df, Species == "F. subverticillata"))+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( twotone_endophyte_color_scheme[1], twotone_endophyte_color_scheme[2]))+
  facet_wrap(~Species, scales = "free", ncol = 4)+
  labs(x = "Interannual SD", y = "Count", fill = "Symbiont Status")+
  theme_minimal()+
  theme(strip.text = element_text(face = "italic"),
        axis.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)))
# FESUsurv_sigmayear_hist
ggsave(FESUsurv_sigmayear_hist, filename = "FESUsurv_sigmayear_hist_quadXorigin.png", width = 6, height = 4)


POALfert_sigmayear_hist <- ggplot(data = filter(fert_sigmayear_df, Species == "P. alsodes"))+
  geom_histogram(aes(x = sigma_year, fill = Endo), alpha = .9, bins = 250, position = "identity") +
  scale_fill_manual(values = c( twotone_endophyte_color_scheme[1], twotone_endophyte_color_scheme[2]))+
  facet_wrap(~Species, scales = "free", ncol = 4)+
  labs(x = "Interannual SD", y = "Count", fill = "Symbiont Status")+
  theme_minimal()+
  theme(strip.text = element_text(face = "italic"),
        axis.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)))
# POALfert_sigmayear_hist
ggsave(POALfert_sigmayear_hist, filename = "POALfert_sigmayear_hist_quadXorigin.png", width = 6, height = 4)




######## Plots of all posteriors from vital rate models #########
parameter_key <- c("beta0[1,1]" = expression(paste(beta[0][h][o],"- AGPE; Orig.")), "beta0[2,1]" = expression(paste(beta[0][h][o],"- ELRI; Rec.")), "beta0[3,1]" = expression(paste(beta[0][h][o],"- ELRI; Orig.")), "beta0[4,1]" = expression(paste(beta[0][h][o],"- FESU; Orig.")), "beta0[5,1]" = expression(paste(beta[0][h][o],"- LOAR; Orig.")), "beta0[6,1]" = expression(paste(beta[0][h][o],"- POAL: Orig.")), "beta0[7,1]" = expression(paste(beta[0][h],"- POSY; Orig.")),
                   "beta0[1,2]" = expression(paste(beta[0][h][o],"- AGPE; Rec.")), "beta0[2,2]" = expression(paste(beta[0][h][o],"- ELRI; Rec.")), "beta0[3,2]" = expression(paste(beta[0][h][o],"- ELRI; Rec.")), "beta0[4,2]" = expression(paste(beta[0][h][o],"- FESU; Rec.")), "beta0[5,2]" = expression(paste(beta[0][h][o],"- LOAR; Rec.")), "beta0[6,2]" = expression(paste(beta[0][h][o],"- POAL: Rec.")), "beta0[7,2]" = expression(paste(beta[0][h],"- POSY; Orig.")),
                                      
                   "betaendo[1]" = expression(paste(beta[1][h],"- AGPE Endo")), "betaendo[2]" = expression(paste(beta[1][h], "- ELRI Endo")), "betaendo[3]" = expression(paste(beta[1][h], "- ELVI Endo")), "betaendo[4]" = expression(paste( beta[1][h], "- FESU Endo")), "betaendo[5]" = expression(paste( beta[1][h], "- LOAR Endo")), "betaendo[6]" = expression(paste( beta[1][h], "- POAL Endo")), "betaendo[7]" = expression(paste( beta[1][h], "- POSY Endo")),
                   
                   "betasize[1,1]" = expression(paste(beta[2][h][o],"- AGPE Size; Orig.")), "betasize[2,1]" = expression(paste(beta[2][h][o], "- ELRI Size; Orig.")), "betasize[3,1]" = expression(paste(beta[2][h][o], "- ELVI Size; Orig.")), "betasize[4,1]" = expression(paste( beta[2][h][o], "- FESU Size; Orig.")), "betasize[5,1]" = expression(paste( beta[2][h][o], "- LOAR Size; Orig.")), "betasize[6,1]" = expression(paste( beta[2][h][o], "- POAL Size; Orig.")), "betasize[7,1]" = expression(paste( beta[2][h][o], "- POSY Size; Orig.")),
                   "betasize[1,2]" = expression(paste(beta[2][h][o],"- AGPE Size; Rec.")), "betasize[2,2]" = expression(paste(beta[2][h][o], "- ELRI Size; Rec.")), "betasize[3,2]" = expression(paste(beta[2][h][o], "- ELVI Size; Rec.")), "betasize[4,2]" = expression(paste( beta[2][h][o], "- FESU Size; Rec.")), "betasize[5,2]" = expression(paste( beta[2][h][o], "- LOAR Size; Rec.")), "betasize[6,2]" = expression(paste( beta[2][h][o], "- POAL Size; Rec.")), "betasize[7,2]" = expression(paste( beta[2][h][o], "- POSY Size; Rec.")),
                   
                   "betasize_2[1,1]" = expression(paste(beta[3][h][o],"- AGPE Quadratic; Orig.")), "betasize_2[2,1]" = expression(paste(beta[3][h][o], "- ELRI Quadratic; Orig.")), "betasize_2[3,1]" = expression(paste(beta[3][h][o], "- ELVI Quadratic; Orig.")), "betasize_2[4,1]" = expression(paste( beta[3][h][o], "- FESU Quadratic; Orig.")), "betasize_2[5,1]" = expression(paste( beta[3][h][o], "- LOAR Quadratic; Orig.")), "betasize_2[6,1]" = expression(paste( beta[3][h][o], "- POAL Quadratic; Orig.")), "betasize_2[7,1]" = expression(paste( beta[3][h][o], "- POSY Quadratic; Orig.")),
                   "betasize_2[1,2]" = expression(paste(beta[3][h][o],"- AGPE Quadratic; Rec.")), "betasize_2[2,2]" = expression(paste(beta[3][h][o], "- ELRI Quadratic; Rec.")), "betasize_2[3,2]" = expression(paste(beta[3][h][o], "- ELVI Quadratic; Rec.")), "betasize_2[4,2]" = expression(paste( beta[3][h][o], "- FESU Quadratic; Rec.")), "betasize_2[5,2]" = expression(paste( beta[3][h][o], "- LOAR Quadratic; Rec.")), "betasize_2[6,2]" = expression(paste( beta[3][h][o], "- POAL Quadratic; Rec.")), "betasize_2[7,2]" = expression(paste( beta[3][h][o], "- POSY Quadratic; Rec.")),
                   
                   "sigma_plot" = expression(paste(sigma[rho], "- Plot SD")),
                   "sigma_year[1,1]" = expression(paste(sigma[tau], "- AGPE S- Year SD")), "sigma_year[1,2]" = expression(paste(sigma[tau], "- AGPE S+ Year SD")),
                   "sigma_year[2,1]" = expression(paste(sigma[tau], "- ELRI S- Year SD")), "sigma_year[2,2]" = expression(paste(sigma[tau], "- ELRI S+ Year SD")),
                   "sigma_year[3,1]" = expression(paste(sigma[tau], "- ELVI S- Year SD")), "sigma_year[3,2]" = expression(paste(sigma[tau], "- ELVI S+ Year SD")),
                   "sigma_year[4,1]" = expression(paste(sigma[tau], "- FESU S- Year SD")), "sigma_year[4,2]" = expression(paste(sigma[tau], "- FESU S+ Year SD")),
                   "sigma_year[5,1]" = expression(paste(sigma[tau], "- LOAR S- Year SD")), "sigma_year[5,2]" = expression(paste(sigma[tau], "- LOAR S+ Year SD")),
                   "sigma_year[6,1]" = expression(paste(sigma[tau], "- POAL S- Year SD")), "sigma_year[6,2]" = expression(paste(sigma[tau], "- POAL S+ Year SD")),
                   "sigma_year[7,1]" = expression(paste(sigma[tau], "- POSY S- Year SD")), "sigma_year[7,2]" = expression(paste(sigma[tau], "- POSY S+ Year SD")))
sigma_year_scale <- c("sigma_plot", 
                      "sigma_year[1,1]", "sigma_year[1,2]",
                      "sigma_year[2,1]", "sigma_year[2,2]",
                      "sigma_year[3,1]", "sigma_year[3,2]",
                      "sigma_year[4,1]", "sigma_year[4,2]",
                      "sigma_year[5,1]", "sigma_year[5,2]",
                      "sigma_year[6,1]", "sigma_year[6,2]",
                      "sigma_year[7,1]", "sigma_year[7,2]")

# Surv
surv_beta_posteriors <-  mcmc_areas(surv_fit, prob = 0.8, regex_pars = c("beta0", "betasize", "betaendo", "betaorigin"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=rev) + theme(axis.text.y = element_text(size = rel(.8)))
surv_sigma_posteriors <-  mcmc_areas(surv_fit, prob = 0.8, regex_pars = c("sigma_plot", "sigma_year"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=sigma_year_scale)+ theme(axis.text.y = element_text(size = rel(.8)))

surv_posteriors_plot <- surv_beta_posteriors + surv_sigma_posteriors + plot_annotation(title = "Adult Survival", subtitle = "Posterior mean with 80% credible intervals")

ggsave(surv_posteriors_plot, filename = "surv_posteriors_plot.png", width = 6, height = 8)

#Sdlg Surv
seedsurv_beta_posteriors <-  mcmc_areas(surv_fit_seedling, prob = 0.8, regex_pars = c("beta0", "betasize", "betaendo", "betaorigin"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=rev) + theme(axis.text.y = element_text(size = rel(.8)))
seedsurv_sigma_posteriors <-  mcmc_areas(surv_fit_seedling, prob = 0.8, regex_pars = c("sigma_plot", "sigma_year"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=sigma_year_scale)+ theme(axis.text.y = element_text(size = rel(.8)))

seedsurv_posteriors_plot <- seedsurv_beta_posteriors + seedsurv_sigma_posteriors + plot_annotation(title = "Seedling Survival", subtitle = "Posterior mean with 80% credible intervals")

ggsave(seedsurv_posteriors_plot, filename = "seedsurv_posteriors_plot.png", width = 6, height = 8)

# Grow
grow_beta_posteriors <-  mcmc_areas(grow_fit, prob = 0.8, regex_pars = c("beta0", "betasize", "betaendo", "betaorigin"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=rev) + theme(axis.text.y = element_text(size = rel(.8)))
grow_sigma_posteriors <-  mcmc_areas(grow_fit, prob = 0.8, regex_pars = c("sigma_plot", "sigma_year"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=sigma_year_scale)+ theme(axis.text.y = element_text(size = rel(.8)))

grow_posteriors_plot <- grow_beta_posteriors + grow_sigma_posteriors + plot_annotation(title = "Adult Growth", subtitle = "Posterior mean with 80% credible intervals")

ggsave(grow_posteriors_plot, filename = "grow_posteriors_plot.png", width = 6, height = 8)

#Sdlg Growth
seedgrow_beta_posteriors <-  mcmc_areas(grow_fit_seedling, prob = 0.8, regex_pars = c("beta0", "betasize", "betaendo", "betaorigin"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=rev) + theme(axis.text.y = element_text(size = rel(.8)))
seedgrow_sigma_posteriors <-  mcmc_areas(grow_fit_seedling, prob = 0.8, regex_pars = c("sigma_plot", "sigma_year"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=sigma_year_scale)+ theme(axis.text.y = element_text(size = rel(.8)))

seedgrow_posteriors_plot <- seedgrow_beta_posteriors + seedgrow_sigma_posteriors + plot_annotation(title = "Seedling Growth", subtitle = "Posterior mean with 80% credible intervals")

ggsave(seedgrow_posteriors_plot, filename = "seedgrow_posteriors_plot.png", width = 6, height = 8)



# Flowering
flw_beta_posteriors <-  mcmc_areas(flw_fit, prob = 0.8, regex_pars = c("beta0", "betasize", "betaendo", "betaorigin"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=rev) + theme(axis.text.y = element_text(size = rel(.8)))
flw_sigma_posteriors <-  mcmc_areas(flw_fit, prob = 0.8, regex_pars = c("sigma_plot", "sigma_year"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=sigma_year_scale)+ theme(axis.text.y = element_text(size = rel(.8)))

flw_posteriors_plot <- flw_beta_posteriors + flw_sigma_posteriors + plot_annotation(title = "Flowering Probability", subtitle = "Posterior mean with 80% credible intervals")

ggsave(flw_posteriors_plot, filename = "flw_posteriors_plot.png", width = 6, height = 8)

# Fertility
fert_beta_posteriors <-  mcmc_areas(fert_fit, prob = 0.8, regex_pars = c("beta0", "betasize", "betaendo", "betaorigin"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=rev) + theme(axis.text.y = element_text(size = rel(.8)))
fert_sigma_posteriors <-  mcmc_areas(fert_fit, prob = 0.8, regex_pars = c("sigma_plot", "sigma_year"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=sigma_year_scale)+ theme(axis.text.y = element_text(size = rel(.8)))

fert_posteriors_plot <- fert_beta_posteriors + fert_sigma_posteriors + plot_annotation(title = "Infl. Production", subtitle = "Posterior mean with 80% credible intervals")

ggsave(fert_posteriors_plot, filename = "fert_posteriors_plot.png", width = 6, height = 8)


# Spikelet
spike_beta_posteriors <-  mcmc_areas(spike_fit, prob = 0.8, regex_pars = c("beta0", "betasize", "betaendo", "betaorigin"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=rev) + theme(axis.text.y = element_text(size = rel(.8)))
spike_sigma_posteriors <-  mcmc_areas(spike_fit, prob = 0.8, regex_pars = c("sigma_plot", "sigma_year"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=sigma_year_scale)+ theme(axis.text.y = element_text(size = rel(.8)))

spike_posteriors_plot <- spike_beta_posteriors + spike_sigma_posteriors + plot_annotation(title = "Spikelets/Infl.", subtitle = "Posterior mean with 80% credible intervals")

ggsave(spike_posteriors_plot, filename = "spike_posteriors_plot.png", width = 6, height = 8)

# Seed means
seedmean_beta_posteriors <-  mcmc_areas(seedmean_fit, prob = 0.8, regex_pars = c("beta0", "betasize", "betaendo", "betaorigin"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=rev) + theme(axis.text.y = element_text(size = rel(.8)))
#seedmean_sigma_posteriors <-  mcmc_areas(seedmean_fit, prob = 0.8, regex_pars = c("sigma_plot", "sigma_year", "sigma0"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=sigma_year_scale)+ theme(axis.text.y = element_text(size = rel(.8)))

seedmean_posteriors_plot <- seedmean_beta_posteriors + plot_annotation(title = "Seeds/Spikelet", subtitle = "Posterior mean with 80% credible intervals")

ggsave(seedmean_posteriors_plot, filename = "seedmean_posteriors_plot.png", width = 6, height = 8)

# Seed to seedling
stos_beta_posteriors <-  mcmc_areas(stos_fit, prob = 0.8, regex_pars = c("beta0", "betasize", "betaendo", "betaorigin"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=rev) + theme(axis.text.y = element_text(size = rel(.8)))
stos_sigma_posteriors <-  mcmc_areas(stos_fit, prob = 0.8, regex_pars = c("sigma_plot", "sigma_year"))+labs(title = "", subtitle = "") + theme_minimal() + scale_y_discrete(labels = parameter_key, limits=sigma_year_scale)+ theme(axis.text.y = element_text(size = rel(.8)))

stos_posteriors_plot <- stos_beta_posteriors + stos_sigma_posteriors + plot_annotation(title = "Recruitment", subtitle = "Posterior mean with 80% credible intervals")

ggsave(stos_posteriors_plot, filename = "stos_posteriors_plot.png", width = 6, height = 8)

seedsurv_endomean_posteriors <-  mcmc_areas(surv_fit_seedling, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Seedling Survival", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
grow_endomean_posteriors <- mcmc_areas(grow_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Adult Growth", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
seedgrow_endomean_posteriors <- mcmc_areas(grow_fit_seedling, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Seedling Growth", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
flw_endomean_posteriors <- mcmc_areas(flw_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Flowering", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
fert_endomean_posteriors <- mcmc_areas(fert_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Fertility", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
spike_endomean_posteriors <- mcmc_areas(spike_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Spikelets per infl.", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
seedmean_endomean_posteriors <- mcmc_areas(seedmean_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Mean seeds per spikelet", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)
stos_endomean_posteriors <- mcmc_areas(stos_fit, prob = 0.8, regex_pars = c("betaendo"))+labs(title = "Germination", subtitle = "Endophyte effect on mean with 80% credible intervals") + scale_y_discrete(labels = mean_species_key)




## Plots for all vital rates, all species by size with data #####
max_size <- LTREB_full %>% 
  dplyr::filter(origin_01 == 1) %>% 
  dplyr::select(species,species_index, size_t) %>% 
  filter(!is.na(size_t)) %>% 
  group_by(species, species_index) %>% 
  summarise(actual_max_size = max(size_t),
            max_size = quantile(size_t,probs=0.975))

n_post_draws <- 500
post_draws <- sample.int(7500, n_post_draws)

x_seq_length <- 100
x_seq <- array(dim= c(x_seq_length,7))
for(s in 1:7){
 x_seq[,s] <-  seq(from = 1, to = filter(max_size, species_index == s)$actual_max_size, length.out = 100)
}

surv_iter <- grow_iter <- flw_iter <- fert_iter <- spike_iter <-array(dim = c(length(x_seq[,1]),2,7, n_post_draws))
surv_mean <- grow_mean <- flw_mean <- fert_mean <- spike_mean<-array(dim = c(length(x_seq[,1]),2,7,3))

survyear_iter <- growyear_iter <- flwyear_iter <- spikeyear_iter <- fertyear_iter <-array(dim = c(length(x_seq[,1]),2,7, (length(unique(LTREB_full$year_t_index))),n_post_draws))
survyear_mean <- growyear_mean <- flwyear_mean <- spikeyear_mean <- fertyear_mean <-array(dim = c(length(x_seq[,1]),2,7,(length(unique(LTREB_full$year_t_index))),3))


sx<-function(x,params){
  invlogit(params$surv_int + params$surv_slope*log(x) + params$surv_slope_2*(log(x)^2))
}
gx <- function(x,params){
  exp(params$grow_int + params$grow_slope*log(x)+ params$grow_slope_2*(log(x)^2))
}
flwx <- function(x,params){
  invlogit(params$flow_int + params$flow_slope*log(x)+ params$flow_slope_2*(log(x)^2))
}
fertx <- function(x,params){
  exp(params$fert_int + params$fert_slope*log(x)+ params$fert_slope_2*(log(x)^2))
}

spikex <- function(x,params){
  exp(params$spike_int + params$spike_slope*log(x)+ params$spike_slope_2*(log(x)^2))
}

for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){

surv_iter[,e,s,i] <- sx(make_params_quadXorigin(species=s,
                      endo_mean=(e-1),
                      endo_var=(e-1),
                      original = 0, # should be =1 to represent recruit
                      draw=post_draws[i],
                      max_size=max_size,
                      rfx=F,
                      surv_par=surv_par,
                      surv_sdlg_par = surv_sdlg_par,
                      grow_par=grow_par,
                      grow_sdlg_par = grow_sdlg_par,
                      flow_par=flow_par,
                      fert_par=fert_par,
                      spike_par=spike_par,
                      seed_par=seed_par,
                      recruit_par=recruit_par), x = x_seq[,s])
grow_iter[,e,s,i] <- gx(make_params_quadXorigin(species=s,
                                    endo_mean=(e-1),
                                    endo_var=(e-1),
                                    original = 0, # should be =1 to represent recruit
                                    draw=post_draws[i],
                                    max_size=max_size,
                                    rfx=F,
                                    surv_par=surv_par,
                                    surv_sdlg_par = surv_sdlg_par,
                                    grow_par=grow_par,
                                    grow_sdlg_par = grow_sdlg_par,
                                    flow_par=flow_par,
                                    fert_par=fert_par,
                                    spike_par=spike_par,
                                    seed_par=seed_par,
                                    recruit_par=recruit_par), x = x_seq[,s])
flw_iter[,e,s,i] <- flwx(make_params_quadXorigin(species=s,
                                    endo_mean=(e-1),
                                    endo_var=(e-1),
                                    original = 0, # should be =1 to represent recruit
                                    draw=post_draws[i],
                                    max_size=max_size,
                                    rfx=F,
                                    surv_par=surv_par,
                                    surv_sdlg_par = surv_sdlg_par,
                                    grow_par=grow_par,
                                    grow_sdlg_par = grow_sdlg_par,
                                    flow_par=flow_par,
                                    fert_par=fert_par,
                                    spike_par=spike_par,
                                    seed_par=seed_par,
                                    recruit_par=recruit_par), x = x_seq[,s])
fert_iter[,e,s,i] <- fertx(make_params_quadXorigin(species=s,
                                    endo_mean=(e-1),
                                    endo_var=(e-1),
                                    original = 0, # should be =1 to represent recruit
                                    draw=post_draws[i],
                                    max_size=max_size,
                                    rfx=F,
                                    surv_par=surv_par,
                                    surv_sdlg_par = surv_sdlg_par,
                                    grow_par=grow_par,
                                    grow_sdlg_par = grow_sdlg_par,
                                    flow_par=flow_par,
                                    fert_par=fert_par,
                                    spike_par=spike_par,
                                    seed_par=seed_par,
                                    recruit_par=recruit_par), x = x_seq[,s])
spike_iter[,e,s,i] <- spikex(make_params_quadXorigin(species=s,
                                       endo_mean=(e-1),
                                       endo_var=(e-1),
                                       original = 0, # should be =1 to represent recruit
                                       draw=post_draws[i],
                                       max_size=max_size,
                                       rfx=F,
                                       surv_par=surv_par,
                                       surv_sdlg_par = surv_sdlg_par,
                                       grow_par=grow_par,
                                       grow_sdlg_par = grow_sdlg_par,
                                       flow_par=flow_par,
                                       fert_par=fert_par,
                                       spike_par=spike_par,
                                       seed_par=seed_par,
                                       recruit_par=recruit_par), x = x_seq[,s])

    }
  }
}

for(x in 1:length(x_seq[,1])){
  for(e in 1:2){
    for(s in 1:7){
surv_mean[x,e,s,1] <- mean(surv_iter[x,e,s,], na.rm = T)
surv_mean[x,e,s,2:3] <- quantile(surv_iter[x,e,s,], probs = c(.1,.9), na.rm = T)

grow_mean[x,e,s,1] <- mean(grow_iter[x,e,s,], na.rm = T)
grow_mean[x,e,s,2:3] <- quantile(grow_iter[x,e,s,], probs = c(.1,.9), na.rm = T)

flw_mean[x,e,s,1] <- mean(flw_iter[x,e,s,], na.rm = T)
flw_mean[x,e,s,2:3] <- quantile(flw_iter[x,e,s,], probs = c(.1,.9), na.rm = T)

fert_mean[x,e,s,1] <- mean(fert_iter[x,e,s,], na.rm = T)
fert_mean[x,e,s,2:3] <- quantile(fert_iter[x,e,s,], probs = c(.1,.9), na.rm = T)

spike_mean[x,e,s,1] <- mean(spike_iter[x,e,s,], na.rm = T)
spike_mean[x,e,s,2:3] <- quantile(spike_iter[x,e,s,], probs = c(.1,.9), na.rm = T)
    }
  }
}

dimnames(surv_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list, c("mean","twenty","eighty"))
dimnames(grow_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list, c("mean","twenty","eighty"))
dimnames(flw_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list, c("mean","twenty","eighty"))
dimnames(fert_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list, c("mean","twenty","eighty"))
dimnames(spike_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list, c("mean","twenty","eighty"))

# Now the same thing for each year specific vital rate
for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      for(y in 1:14){
      
      survyear_iter[,e,s,y,i] <- sx(make_params_quadXorigin(species=s,
                                                endo_mean=(e-1),
                                                endo_var=(e-1),
                                                original = 0, # should be =1 to represent recruit
                                                draw=post_draws[i],
                                                max_size=max_size,
                                                rfx=T,
                                                year = y,
                                                repro_offset = 0,
                                                surv_par=surv_par,
                                                surv_sdlg_par = surv_sdlg_par,
                                                grow_par=grow_par,
                                                grow_sdlg_par = grow_sdlg_par,
                                                flow_par=flow_par,
                                                fert_par=fert_par,
                                                spike_par=spike_par,
                                                seed_par=seed_par,
                                                recruit_par=recruit_par), x = x_seq[,s])
      growyear_iter[,e,s,y,i] <- gx(make_params_quadXorigin(species=s,
                                          endo_mean=(e-1),
                                          endo_var=(e-1),
                                          original = 0, # should be =1 to represent recruit
                                          draw=post_draws[i],
                                          max_size=max_size,
                                          rfx=T,
                                          year = y,
                                          repro_offset = 0,
                                          surv_par=surv_par,
                                          surv_sdlg_par = surv_sdlg_par,
                                          grow_par=grow_par,
                                          grow_sdlg_par = grow_sdlg_par,
                                          flow_par=flow_par,
                                          fert_par=fert_par,
                                          spike_par=spike_par,
                                          seed_par=seed_par,
                                          recruit_par=recruit_par), x = x_seq[,s])
      flwyear_iter[,e,s,y,i] <- flwx(make_params_quadXorigin(species=s,
                                           endo_mean=(e-1),
                                           endo_var=(e-1),
                                           original = 0, # should be =1 to represent recruit
                                           draw=post_draws[i],
                                           max_size=max_size,
                                           rfx=T,
                                           year = y,
                                           repro_offset = 0,
                                           surv_par=surv_par,
                                           surv_sdlg_par = surv_sdlg_par,
                                           grow_par=grow_par,
                                           grow_sdlg_par = grow_sdlg_par,
                                           flow_par=flow_par,
                                           fert_par=fert_par,
                                           spike_par=spike_par,
                                           seed_par=seed_par,
                                           recruit_par=recruit_par), x = x_seq[,s])
      fertyear_iter[,e,s,y,i] <- fertx(make_params_quadXorigin(species=s,
                                             endo_mean=(e-1),
                                             endo_var=(e-1),
                                             original = 0, # should be =1 to represent recruit
                                             draw=post_draws[i],
                                             max_size=max_size,
                                             rfx=T,
                                             year = y,
                                             repro_offset = 0,
                                             surv_par=surv_par,
                                             surv_sdlg_par = surv_sdlg_par,
                                             grow_par=grow_par,
                                             grow_sdlg_par = grow_sdlg_par,
                                             flow_par=flow_par,
                                             fert_par=fert_par,
                                             spike_par=spike_par,
                                             seed_par=seed_par,
                                             recruit_par=recruit_par), x = x_seq[,s])
      spikeyear_iter[,e,s,y,i] <- spikex(make_params_quadXorigin(species=s,
                                                   endo_mean=(e-1),
                                                   endo_var=(e-1),
                                                   original = 0, # should be =1 to represent recruit
                                                   draw=post_draws[i],
                                                   max_size=max_size,
                                                   rfx=T,
                                                   year = y,
                                                   repro_offset = 0,
                                                   surv_par=surv_par,
                                                   surv_sdlg_par = surv_sdlg_par,
                                                   grow_par=grow_par,
                                                   grow_sdlg_par = grow_sdlg_par,
                                                   flow_par=flow_par,
                                                   fert_par=fert_par,
                                                   spike_par=spike_par,
                                                   seed_par=seed_par,
                                                   recruit_par=recruit_par), x = x_seq[,s])
      
      }
    }
  }
}
        
for(x in 1:length(x_seq[,1])){
  for(e in 1:2){
    for(s in 1:7){
      for(y in 1:14){
      survyear_mean[x,e,s,y,1] <- mean(survyear_iter[x,e,s,y,], na.rm = T)
      survyear_mean[x,e,s,y,2:3] <- quantile(survyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)

      growyear_mean[x,e,s,y,1] <- mean(growyear_iter[x,e,s,y,], na.rm = T)
      growyear_mean[x,e,s,y,2:3] <- quantile(growyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)

      flwyear_mean[x,e,s,y,1] <- mean(flwyear_iter[x,e,s,y,], na.rm = T)
      flwyear_mean[x,e,s,y,2:3] <- quantile(flwyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)

      fertyear_mean[x,e,s,y,1] <- mean(fertyear_iter[x,e,s,y,], na.rm = T)
      fertyear_mean[x,e,s,y,2:3] <- quantile(fertyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)
      
      spikeyear_mean[x,e,s,y,1] <- mean(spikeyear_iter[x,e,s,y,], na.rm = T)
      spikeyear_mean[x,e,s,y,2:3] <- quantile(spikeyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)
      }
    }
  }
}

dimnames(survyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list,c(2008:2021), c("mean","twenty","eighty"))
dimnames(growyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list,c(2008:2021), c("mean","twenty","eighty"))
dimnames(flwyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list,c(2008:2021), c("mean","twenty","eighty"))
dimnames(fertyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list,c(2008:2021), c("mean","twenty","eighty"))
dimnames(spikeyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list,c(2008:2021), c("mean","twenty","eighty"))

#Now I'm gonna make these into  tidy dataframes for plotting

dimnames(x_seq) <- list(paste0("x_size", 1:length(x_seq[,1])), paste0(species_code_list))

x_seq_df <- as_tibble(x_seq) %>%
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = -no_row,
               names_to = "Species",
               values_to = "x_seq") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  mutate(log_x_seq = log(x_seq))


surv_mean_df <- as_tibble(surv_mean)  %>%   
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("surv") ) %>% 
  separate(name, c("Endo", "Species","quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "surv") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)


survyear_mean_df <- as_tibble(survyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = -no_row,
               values_to = c("surv")) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "surv") %>%
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)

grow_mean_df <- as_tibble(grow_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("size_t1") )  %>% 
  separate(name, c("Endo", "Species","quantile"), "\\.") %>%   
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "size_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)

growyear_mean_df <- as_tibble(growyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("size_t1") ) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "size_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)

flw_mean_df <- as_tibble(flw_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("flw_t1") ) %>% 
  separate(name, c("Endo", "Species","quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "flw_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)

flwyear_mean_df <- as_tibble(flwyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("flw_t1") ) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "flw_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)

fert_mean_df <- as_tibble(fert_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("fert_t1") ) %>% 
  separate(name, c("Endo", "Species","quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "fert_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)

fertyear_mean_df <- as_tibble(fertyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("fert_t1") ) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "fert_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)

spike_mean_df <- as_tibble(spike_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("spike_t1") ) %>% 
  separate(name, c("Endo", "Species","quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "spike_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df) 


spikeyear_mean_df <- as_tibble(spikeyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("spike_t1") ) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "spike_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)


# Bin Data by size and then by year for plotting

bin_by_size_t <- function(df_raw, vr, nbins){
  require(tidyverse)
  df_raw <- ungroup(df_raw)
  size_bin_df <- df_raw %>% 
    rename(Endo = endo_01, Year = year_t1, Species = species) %>%
    mutate(size_bin = cut(logsize_t, breaks = nbins)) %>%
    group_by(size_bin, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_vr = mean({{vr}},na.rm=T),
              samplesize = n()) %>% 
    mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
    mutate(Endo = case_when(Endo == 0 ~ "S-", 
                           Endo == 1 ~ "S+"))
  
  return(size_bin_df)
}
bin_by_size_t1 <- function(df_raw, vr, nbins){
  require(tidyverse)
  df_raw <- ungroup(df_raw)
  size_bin_df <- df_raw %>% 
    rename(Endo = endo_01, Year = year_t1, Species = species) %>%
    mutate(size_bin = cut(logsize_t1, breaks = nbins)) %>%
    group_by(size_bin, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t1),na.rm=T),
              mean_vr = mean({{vr}},na.rm=T),
              samplesize = n()) %>% 
    mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
    mutate(Endo = case_when(Endo == 0 ~ "S-", 
                           Endo == 1 ~ "S+"))
  
  return(size_bin_df)
}


bin_by_year_size_t <- function(df_raw, vr, nbins){
  require(tidyverse)
  df_raw <- ungroup(df_raw)
  size_bin_df <- df_raw %>% 
    rename(Endo = endo_01, Year = year_t1, Species = species) %>%
    mutate(size_bin = cut(logsize_t, breaks = nbins)) %>%
    group_by(size_bin, Year, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_vr = mean({{vr}},na.rm=T),
              samplesize = n()) %>% 
    mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
    mutate(Endo = case_when(Endo == 0 ~ "S-", 
                            Endo == 1 ~ "S+"))
  
  return(size_bin_df)
}
bin_by_year_size_t1 <- function(df_raw, vr, nbins){
  require(tidyverse)
  df_raw <- ungroup(df_raw)
  size_bin_df <- df_raw %>% 
    rename(Endo = endo_01, Year = year_t1, Species = species) %>%
    mutate(size_bin = cut(logsize_t1, breaks = nbins)) %>%
    group_by(size_bin, Year, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t1),na.rm=T),
              mean_vr = mean({{vr}},na.rm=T),
              samplesize = n()) %>% 
    mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
    mutate(Endo = case_when(Endo == 0 ~ "S-", 
                            Endo == 1 ~ "S+")) %>% 

  return(size_bin_df)
}

neat_names <- function(df_raw, vr){
  require(tidyverse)
  df_raw <- ungroup(df_raw)
  renamed_df <- df_raw %>% 
    rename(Endo = endo_01, Year = year_t1, Species = species) %>%
    mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
    mutate(Endo = case_when(Endo == 0 ~ "S-", 
                            Endo == 1 ~ "S+"))
}


surv_sizebin <- bin_by_size_t(filter(LTREB_data_forsurv,origin_01==0),vr = surv_t1, nbins = 20)
# seedsurv_sizebin <- bin_by_size_t(filter(LTREB_surv_seedling,origin_01==0),vr = surv_t1, nbins = 40)
grow_sizebin <- bin_by_size_t(filter(LTREB_data_forgrow,origin_01==0),vr = size_t1, nbins = 20)
# seedgrow_sizebin <- bin_by_size_t(filter(LTREB_grow_seedling,origin_01==0), vr = size_t1, nbins = 40)
flw_sizebin <- bin_by_size_t1(filter(LTREB_data_forflw,origin_01==0), vr = FLW_STAT_T1,nbins = 20)
fert_sizebin <- bin_by_size_t1(filter(LTREB_data_forfert,origin_01==0),vr = FLW_COUNT_T1, nbins = 20)
spike_sizebin <- bin_by_size_t1(filter(LTREB_data_forspike,origin_01==0), vr = spike_count_t1, nbins = 20)

surv_yearsizebin <- bin_by_year_size_t(filter(LTREB_data_forsurv,origin_01==0),vr = surv_t1, nbins = 20)
# seedsurv_yearsizebin <- bin_by_year_size_t(filter(LTREB_surv_seedling,origin_01==0),vr = surv_t1, nbins = 40)
grow_yearsizebin <- bin_by_year_size_t(filter(LTREB_data_forgrow,origin_01==0),vr = size_t1, nbins = 20)
# seedgrow_yearsizebin <- bin_by_year_size_t(filter(LTREB_grow_seedling,origin_01==0), vr = size_t1, nbins = 40)
flw_yearsizebin <- bin_by_year_size_t1(filter(LTREB_data_forflw,origin_01==0), vr = FLW_STAT_T1,nbins = 20)
fert_yearsizebin <- bin_by_year_size_t1(filter(LTREB_data_forfert,origin_01==0),vr = FLW_COUNT_T1, nbins = 20)
spike_yearsizebin <- bin_by_year_size_t1(filter(LTREB_data_forspike,origin_01==0), vr = spike_count_t1, nbins = 20)


surv_neatdata <- neat_names(LTREB_data_forsurv)
seedsurv_neatdata <- neat_names(LTREB_surv_seedling)
grow_neatdata <- neat_names(LTREB_data_forgrow)
seedgrow_neatdata <- neat_names(LTREB_grow_seedling)
flw_neatdata <- neat_names(LTREB_data_forflw)
fert_neatdata <- neat_names(LTREB_data_forfert)
spike_neatdata <- neat_names(LTREB_data_forspike)


#The plots
surv_meanplot <- ggplot()+
  geom_point(data = surv_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo, col = Species), alpha = .7) +
  geom_ribbon(data = surv_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, alpha = Endo, fill = Species))+
  geom_line(data = surv_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = Species)) +
  scale_color_manual(values = species_colors)+  scale_fill_manual(values = species_colors)+
  scale_shape_manual( values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + scale_alpha_manual( values = c(.1,.3))+
  facet_wrap(~Species, scales = "free", ncol = 2) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none", col = "none", alpha = guide_legend(order = 1), linetype = guide_legend(order = 1), shape = guide_legend(order = 1), size = guide_legend(order = 2))+
  labs(title = "Adult Survival",subtitle = "Data from Original plants only", y = "Survival Probability", x = expression("log(# of tillers)"[" year t"]), size = "Sample Size")
# surv_meanplot
ggsave(surv_meanplot, filename = "surv_meanplot_quadXorigin.png", width = 8, height = 12)

surv_yearplot <- ggplot()+
  geom_point(data = surv_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .4) +
  # geom_ribbon(data = survyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = survyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = as.factor(Year), group = Year), alpha = .8) +
  scale_shape_manual( values = c(1,19)) + scale_linetype_manual(values = c(2,1)) +
  scale_color_manual( values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free",ncol = 4) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none",shape = guide_legend(order = 2, ncol = 1),linetype = guide_legend(order = 2, ncol = 1), color = guide_legend(order = 1, ncol = 2), size = guide_legend(order = 3))+
  labs(title = "Adult Survival",subtitle = "Data from Original plants only", y = "Survival Probability", x = expression("log(# of tillers)"[" year t"]), size = "Sample Size", col = "Year")
# surv_yearplot
ggsave(surv_yearplot, filename = "surv_yearplot_quadXorigin.png", width = 16, height = 17)

grow_meanplot <- ggplot()+
  # geom_point(data = grow_neatdata, aes(x = logsize_t, y = size_t1, shape = Endo, col = Species), alpha = .1) +
  geom_point(data = grow_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo, col = Species), alpha = .7) +
  geom_ribbon(data = grow_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, alpha = Endo, fill = Species))+
  geom_line(data = grow_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = Species)) +
  scale_color_manual(values = species_colors)+  scale_fill_manual(values = species_colors)+
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + scale_alpha_manual(values = c(.1,.3))+
  facet_wrap(~Species, scales = "free", ncol = 2) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none", col = "none", alpha = guide_legend(order = 1), linetype = guide_legend(order = 1), shape = guide_legend(order = 1), size = guide_legend(order = 2))+
  labs(title = "Adult Growth",subtitle = "Data from Original plants only", y = "# of Tillers", x = expression("log(# of tillers)"[" year t"]), size = "Sample Size")
# grow_meanplot
ggsave(grow_meanplot, filename = "grow_meanplot_quadXorigin.png", width = 8, height = 12)

grow_yearplot <- ggplot()+
  # geom_point(data = grow_neatdata, aes(x = logsize_t, y = size_t1, col = as.factor(Year), shape = Endo), alpha = .4) +
  geom_point(data = grow_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .4) +
  # geom_ribbon(data = growyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = growyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = as.factor(Year),group = Year), alpha = .8) +
  scale_shape_manual( values = c(1,19)) + scale_linetype_manual(values = c(2,1)) +
  scale_color_manual( values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free", ncol = 4) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none",shape = guide_legend(order = 2, ncol = 1),linetype = guide_legend(order = 2, ncol = 1), color = guide_legend(order = 1, ncol = 2), size = guide_legend(order = 3))+
  labs(title = "Adult Growth",subtitle = "Data from Original plants only",  y = "# of Tillers", x = expression("log(# of tillers)"[" year t"]), col = "Year",size = "Sample Size")
# grow_yearplot
ggsave(grow_yearplot, filename = "grow_yearplot_quadXorigin.png", width = 16, height = 17)

flw_meanplot <- ggplot()+
  geom_point(data = flw_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo, col = Species), alpha = .7) +
  geom_ribbon(data = flw_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, alpha = Endo, fill = Species))+
  geom_line(data = flw_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = Species)) +
  scale_color_manual(values = species_colors)+  scale_fill_manual(values = species_colors)+
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + scale_alpha_manual(values = c(.1,.3))+  
  facet_wrap(~Species, scales = "free", ncol = 2) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none", col = "none", alpha = guide_legend(order = 1), linetype = guide_legend(order = 1), shape = guide_legend(order = 1), size = guide_legend(order = 2))+
  labs(title = "Flowering",subtitle = "Data from Original plants only", y = "Flowering Probability", x = expression("log(# of tillers)"[" year t+1"]), size = "Sample Size")
# flw_meanplot
ggsave(flw_meanplot, filename = "flw_meanplot_quadXorigin.png", width = 8, height = 12)


flw_yearplot <- ggplot()+
  geom_point(data = flw_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .4) +
  # geom_ribbon(data = flwyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = flwyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col= as.factor(Year), group = Year)) +
  scale_shape_manual( values = c(1,19)) + scale_linetype_manual(values = c(2,1)) +
  scale_color_manual( values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free", ncol = 4) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none",shape = guide_legend(order = 2, ncol = 1),linetype = guide_legend(order = 2, ncol = 1), color = guide_legend(order = 1, ncol = 2), size = guide_legend(order = 3))+
  labs(title = "Flowering",subtitle = "Data from Original plants only", y = "Flowering Probability", x = expression("log(# of tillers)"[" year t+1"]), col = "Year", size = "Sample Size")
# flw_yearplot
ggsave(flw_yearplot, filename = "flw_yearplot_quadXorigin.png", width = 16, height = 17)

fert_meanplot <- ggplot()+
  geom_point(data = fert_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo, col = Species), alpha = .7) +
  geom_ribbon(data = fert_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, alpha = Endo, fill = Species))+
  geom_line(data = fert_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, color = Species)) +
  scale_color_manual(values = species_colors)+  scale_fill_manual(values = species_colors)+
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + scale_alpha_manual(values = c(.1,.3))+  
  facet_wrap(~Species, scales = "free", ncol = 2) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none", col = "none", alpha = guide_legend(order = 1), linetype = guide_legend(order = 1), shape = guide_legend(order = 1), size = guide_legend(order = 2))+
  labs(title = "Fertility",subtitle = "Data from Original plants only", y = "# of Repro. Tillers", x = expression("log(# of tillers)"[" year t+1"]), size = "Sample Size")
# fert_meanplot
ggsave(fert_meanplot, filename = "fert_meanplot_quadXorigin.png", width = 8, height = 12)

fert_yearplot <- ggplot()+
  geom_point(data = fert_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .4) +
  # geom_ribbon(data = fertyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = fertyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = as.factor(Year), group = Year)) +
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + 
  scale_color_manual(values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free", ncol = 4) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none",shape = guide_legend(order = 2, ncol = 1),linetype = guide_legend(order = 2, ncol = 1), color = guide_legend(order = 1, ncol = 2), size = guide_legend(order = 3))+
  labs(title = "Fertility",subtitle = "Data from Original plants only",  y = "# of Repro. Tillers", x = expression("log(# of tillers)"[" year t+1"]), col = "Year", size = "Sample Size")
# fert_yearplot
ggsave(fert_yearplot, filename = "fert_yearplot_quadXorigin.png",width = 16, height = 17)


spike_meanplot <- ggplot()+
  # geom_point(data = spike_neatdata, aes(x = logsize_t, y = spike_count_t1, shape = Endo, col = Species), alpha = .4) +
  geom_point(data = spike_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo, col = Species), alpha = .7) +
  geom_ribbon(data = spike_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, alpha = Endo, fill = Species))+
  geom_line(data = spike_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, color = Species)) +
  scale_color_manual(values = species_colors)+  scale_fill_manual(values = species_colors)+
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + scale_alpha_manual(values = c(.1,.3))+  
  facet_wrap(~Species, scales = "free", ncol = 2) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none", col = "none", alpha = guide_legend(order = 1), linetype = guide_legend(order = 1), shape = guide_legend(order = 1), size = guide_legend(order = 2))+
  labs(title = "Spikes/Infl.",subtitle = "Data from Original plants only", y = "Spikelet/Infl.", x = expression("log(# of tillers)"[" year t+1"]), size = "Sample Size")
# spike_meanplot
ggsave(spike_meanplot, filename = "spike_meanplot_quadXorigin.png", width = 8, height = 12)

spike_yearplot <- ggplot()+
  # geom_point(data = spike_neatdata, aes(x = logsize_t, y = spike_count_t1, col = as.factor(Year), shape = Endo), alpha = .4) +
  geom_point(data = spike_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .4) +
  # geom_ribbon(data = spikeyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = spikeyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = as.factor(Year), group = Year)) +
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + 
  scale_color_manual(values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free", ncol = 4) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none",shape = guide_legend(order = 2, ncol = 1),linetype = guide_legend(order = 2, ncol = 1), color = guide_legend(order = 1, ncol = 2), size = guide_legend(order = 3))+
  labs(title = "Spikes/Infl.",subtitle = "Data from Original plants only", y = "Spikelet/Infl.", x = expression("log(# of tillers)"[" year t+1"]), col = "Year", size = "Sample Size")
# spike_yearplot
ggsave(spike_yearplot, filename = "spike_yearplot_quadXorigin.png",width = 16, height = 17)


# meaneffect_fitplot <- surv_meanplot+grow_meanplot+flw_meanplot+fert_meanplot+spike_meanplot + 
#   plot_layout( nrow  = 1) + plot_annotation(title = "Endophyte effect on mean of size-structured vital rates", subtitle = "with 80% credible intervals")& 
#   theme(text = element_text(size = 16))
# ggsave(meaneffect_fitplot, filename = "meaneffect_fitplot.png", width = 20, height = 18 )
# 
# 
# vareffect_fitplot <- surv_yearplot+grow_yearplot+flw_yearplot+fert_yearplot+spike_yearplot +
#   plot_layout( nrow  = 1) + plot_annotation(title = "Endophyte effect on interannual variance of size-structured vital rates", subtitle = "with 80% credible intervals")& 
#   theme(text = element_text(size = 16))
# ggsave(vareffect_fitplot, filename = "vareffect_fitplot.png", width = 30, height = 25 )
# 
# meanvareffect_fitplot <-  surv_meanplot + surv_yearplot+grow_meanplot +grow_yearplot +flw_meanplot+flw_yearplot+fert_meanplot +fert_yearplot +spike_meanplot +spike_yearplot + 
#   plot_layout( nrow  = 1,widths = c(1,2,1,2,1,2,1,2,1,2),guides = "collect")
# ggsave(meanvareffect_fitplot, filename = "meanvareffect_fitplot.png", width = 40, height = 25 )  


# plot of just AGPE growth and FESU survival
AGPEgrow_meanplot <- ggplot()+
  geom_point(data = filter(grow_sizebin, Species == "A. perennans"), aes(x = mean_size, y = mean_vr, size = samplesize, fill = Endo), shape  = 21, alpha = .8) +
  geom_ribbon(data = filter(grow_mean_df, Species == "A. perennans"), aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = filter(grow_mean_df, Species == "A. perennans"), aes(x = log_x_seq, y = mean, linetype = Endo, col = Endo)) +
  scale_color_manual(values = twotone_endophyte_color_scheme)+  scale_fill_manual(values = twotone_endophyte_color_scheme)+
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + scale_alpha_manual(values = c(.3,.5))+
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 9),
        strip.background = element_blank(), 
        strip.text.x = element_blank()) + 
  guides(lwd = guide_legend(order = 1))+
  labs(y = "# of tillers in year t+1", x = "log(# of tillers in year t)", size = "Sample Size",  fill = "Symbiont Status", color = "Symbiont Status", linetype = "Symbiont Status")
AGPEgrow_meanplot
ggsave(AGPEgrow_meanplot, filename = "AGPEgrow_meanplot_quadXorigin.png", width = 4.5, height = 4)
  
AGPEgrow_yearplot <- ggplot()+
  geom_point(data = filter(grow_yearsizebin, Species == "A. perennans"), aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .5) +
  geom_line(data = filter(growyear_mean_df, Species == "A. perennans"), aes(x = log_x_seq, y = mean, linetype = Endo, col = Year)) +
  scale_shape_manual(values = c(1,19))+ 
  scale_color_manual(values = yearcolors)+ # this is using the set of colors above, but you could also supply hex codes
  scale_linetype_manual(values = c(2,1))+
  facet_wrap(~Species + Endo, scales = "free", ncol = 2) + 
  theme_classic() + theme(strip.background = element_blank(), strip.text.x = element_blank()) + labs(title = "Adult Growth", y = "# of tillers in year t+1", x = "log(# of tillers in year t)")
AGPEgrow_yearplot

FESUsurv_meanplot <- ggplot()+
  geom_point(data = filter(surv_sizebin, Species == "F. subverticillata"), aes(x = mean_size, y = mean_vr, size = samplesize, fill = Endo), shape  = 21, alpha = .8) +
  geom_ribbon(data = filter(surv_mean_df, Species == "F. subverticillata"), aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = filter(surv_mean_df, Species == "F. subverticillata"), aes(x = log_x_seq, y = mean, linetype = Endo, col = Endo)) +
  scale_color_manual(values = twotone_endophyte_color_scheme)+  scale_fill_manual(values = twotone_endophyte_color_scheme)+
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + scale_alpha_manual(values = c(.3,.5))+
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 9),
        strip.background = element_blank(), 
        strip.text.x = element_blank()) + 
  guides(lwd = guide_legend(order = 1))+
  labs( y = "Survival Probability", x = "log(# of tillers in year t)", size = "Sample Size", fill = "Symbiont Status", color = "Symbiont Status", linetype = "Symbiont Status")
FESUsurv_meanplot
ggsave(FESUsurv_meanplot, filename = "FESUsurv_meanplot_quadXorigin.png", width = 4.5, height = 4)

FESUsurv_yearplot <- ggplot()+
  geom_point(data = filter(surv_yearsizebin, Species == "F. subverticillata"), aes(x = mean_size, y = mean_vr, size = samplesize, fill = Endo),shape = 21, alpha = .8) +
  geom_line(data = filter(survyear_mean_df, Species == "F. subverticillata"), aes(x = log_x_seq, y = mean, col = Endo, group = Year)) +
  scale_shape_manual(values = c(1,19))+ 
  scale_fill_manual(values = twotone_endophyte_color_scheme)+   scale_color_manual(values = twotone_endophyte_color_scheme)+
  scale_linetype_manual(values = c(2,1))+
  facet_wrap(~Species + Endo, scales = "free", ncol = 2) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        strip.background = element_blank(), 
        strip.text.x = element_blank()) + 
  guides(fill = "none",lwd = guide_legend(order = 1))+
  labs( y = "Survival Probability", x = "log(# of tillers in year t)", linetype = "Symbiont Status", fill = "Year", size = "Sample Size")
FESUsurv_yearplot
ggsave(FESUsurv_yearplot, filename = "FESUsurv_yearplot_quadXorigin.png", width = 10, height = 7)

POALfert_meanplot <- ggplot()+
  geom_point(data = filter(fert_sizebin, Species == "P. alsodes"), aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo), alpha = .5) +
  geom_ribbon(data = filter(fert_mean_df, Species == "P. alsodes"), aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = filter(fert_mean_df, Species == "P. alsodes"), aes(x = log_x_seq, y = mean, linetype = Endo)) +
  scale_shape_manual(values = c(1,19))+   scale_linetype_manual(values = c(2,1))+ scale_fill_manual(values = c( endophyte_color_scheme[3], endophyte_color_scheme[5]))+ 
  facet_wrap(~Species, scales = "free", ncol = 1) + 
  theme_classic() + theme(strip.background = element_blank())+ labs(title = "Adult Fertility", y = "# of repro. tillers", x = "log(# of tillers in year t+1)", col = "Year", fill = "Year", lwd = "Sample Size")
POALfert_meanplot
ggsave(POALfert_meanplot, filename = "POALfert_meanplot_quadXorigin.png", width = 3.5, height = 4)


POALfert_yearplot <- ggplot()+
  geom_point(data = filter(fert_yearsizebin, Species == "P. alsodes"), aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo, col = as.numeric(Year)), alpha = .6) +
  # geom_ribbon(data = fertyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = filter(fertyear_mean_df, Species == "P. alsodes"), aes(x = log_x_seq, y = mean, linetype = Endo, col = as.numeric(Year), group = Year), alpha = .8) +
  scale_shape_manual(values = c(1,19))+ 
  scale_color_gradient2(low = "grey10", mid = "grey90", high = species_colors[6], n.breaks = 14, midpoint = 2009, guide = "legend")+
  scale_linetype_manual(values = c(2,1))+
  facet_wrap(~Species + Endo, ncol = 2) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.background = element_blank(), 
        strip.text.x = element_blank()) + 
  guides(fill = "none", lwd = guide_legend(order = 1))+
  # lims(y = c(0,50))+
  labs( y = "# of repro. tillers", x = "log(# of tillers in year t+1)", col = "Year", fill = "Year", size = "Sample Size")
POALfert_yearplot
ggsave(POALfert_yearplot, filename = "POALfert_yearplot_quadXorigin.png",  width = 8, height = 7)



AGPE_growplot <- AGPEgrow_meanplot+AGPEgrow_yearplot+plot_layout( nrow  = 1,
                                                                 widths = c(1,2),
                                                                 guides = "collect")
ggsave(AGPE_growplot, filename = "AGPE_growplot.png", height = 4, width = 8)

FESU_survplot <- FESUsurv_meanplot+FESUsurv_yearplot+plot_layout( nrow  = 1,
                                                                  widths = c(1,2),
                                                                  guides = "collect")
ggsave(FESU_survplot, filename = "FESU_survplot.png", height = 4, width = 8)

POAL_fertplot <- POALfert_meanplot+POALfert_yearplot+plot_layout( nrow  = 1,
                                                                  widths = c(1,2),
                                                                  guides = "collect")
ggsave(POAL_fertplot, filename = "POAL_fertplot.png", height = 4, width = 8)



######## Remaking the size plots specific to recruit plants to demonstrate that this is contributing to some of the apparent bad fit in certain vital rates.#####
max_size <- LTREB_full %>% 
  dplyr::filter(origin_01 == 1) %>%
  dplyr::select(species,species_index, size_t) %>% 
  filter(!is.na(size_t)) %>% 
  group_by(species, species_index) %>% 
  summarise(actual_max_size = max(size_t),
            max_size = quantile(size_t,probs=0.975))
x_seq_length <- 100
x_seq <- array(dim= c(x_seq_length,7))
for(s in 1:7){
  x_seq[,s] <-  seq(from = 1, to = filter(max_size, species_index == s)$actual_max_size, length.out = 100)
}

recruit_surv_iter <- recruit_grow_iter <- recruit_flw_iter <- recruit_fert_iter <- recruit_spike_iter <-array(dim = c(length(x_seq[,1]),2,7, n_post_draws))
recruit_surv_mean <- recruit_grow_mean <- recruit_flw_mean <- recruit_fert_mean <- recruit_spike_mean<-array(dim = c(length(x_seq[,1]),2,7,3))

recruit_survyear_iter <- recruit_growyear_iter <- recruit_flwyear_iter <- recruit_spikeyear_iter <- recruit_fertyear_iter <-array(dim = c(length(x_seq[,1]),2,7, (length(unique(LTREB_full$year_t_index))),n_post_draws))
recruit_survyear_mean <- recruit_growyear_mean <- recruit_flwyear_mean <- recruit_spikeyear_mean <- recruit_fertyear_mean <-array(dim = c(length(x_seq[,1]),2,7,(length(unique(LTREB_full$year_t_index))),3))

for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      
      recruit_surv_iter[,e,s,i] <- sx(make_params_quadXorigin(species=s,
                                          endo_mean=(e-1),
                                          endo_var=(e-1),
                                          original = 1, # should be =1 to represent recruit
                                          draw=post_draws[i],
                                          max_size=max_size,
                                          rfx=F,
                                          surv_par=surv_par,
                                          surv_sdlg_par = surv_sdlg_par,
                                          grow_par=grow_par,
                                          grow_sdlg_par = grow_sdlg_par,
                                          flow_par=flow_par,
                                          fert_par=fert_par,
                                          spike_par=spike_par,
                                          seed_par=seed_par,
                                          recruit_par=recruit_par), x = x_seq[,s])
      recruit_grow_iter[,e,s,i] <- gx(make_params_quadXorigin(species=s,
                                          endo_mean=(e-1),
                                          endo_var=(e-1),
                                          original = 1, # should be =1 to represent recruit
                                          draw=post_draws[i],
                                          max_size=max_size,
                                          rfx=F,
                                          surv_par=surv_par,
                                          surv_sdlg_par = surv_sdlg_par,
                                          grow_par=grow_par,
                                          grow_sdlg_par = grow_sdlg_par,
                                          flow_par=flow_par,
                                          fert_par=fert_par,
                                          spike_par=spike_par,
                                          seed_par=seed_par,
                                          recruit_par=recruit_par), x = x_seq[,s])
      recruit_flw_iter[,e,s,i] <- flwx(make_params_quadXorigin(species=s,
                                           endo_mean=(e-1),
                                           endo_var=(e-1),
                                           original = 1, # should be =1 to represent recruit
                                           draw=post_draws[i],
                                           max_size=max_size,
                                           rfx=F,
                                           surv_par=surv_par,
                                           surv_sdlg_par = surv_sdlg_par,
                                           grow_par=grow_par,
                                           grow_sdlg_par = grow_sdlg_par,
                                           flow_par=flow_par,
                                           fert_par=fert_par,
                                           spike_par=spike_par,
                                           seed_par=seed_par,
                                           recruit_par=recruit_par), x = x_seq[,s])
      recruit_fert_iter[,e,s,i] <- fertx(make_params_quadXorigin(species=s,
                                             endo_mean=(e-1),
                                             endo_var=(e-1),
                                             original = 1, # should be =1 to represent recruit
                                             draw=post_draws[i],
                                             max_size=max_size,
                                             rfx=F,
                                             surv_par=surv_par,
                                             surv_sdlg_par = surv_sdlg_par,
                                             grow_par=grow_par,
                                             grow_sdlg_par = grow_sdlg_par,
                                             flow_par=flow_par,
                                             fert_par=fert_par,
                                             spike_par=spike_par,
                                             seed_par=seed_par,
                                             recruit_par=recruit_par), x = x_seq[,s])
      recruit_spike_iter[,e,s,i] <- spikex(make_params_quadXorigin(species=s,
                                               endo_mean=(e-1),
                                               endo_var=(e-1),
                                               original = 1, # should be =1 to represent recruit
                                               draw=post_draws[i],
                                               max_size=max_size,
                                               rfx=F,
                                               surv_par=surv_par,
                                               surv_sdlg_par = surv_sdlg_par,
                                               grow_par=grow_par,
                                               grow_sdlg_par = grow_sdlg_par,
                                               flow_par=flow_par,
                                               fert_par=fert_par,
                                               spike_par=spike_par,
                                               seed_par=seed_par,
                                               recruit_par=recruit_par), x = x_seq[,s])
      
    }
  }
}

for(x in 1:length(x_seq[,1])){
  for(e in 1:2){
    for(s in 1:7){
      recruit_surv_mean[x,e,s,1] <- mean(recruit_surv_iter[x,e,s,], na.rm = T)
      recruit_surv_mean[x,e,s,2:3] <- quantile(recruit_surv_iter[x,e,s,], probs = c(.1,.9), na.rm = T)
      
      recruit_grow_mean[x,e,s,1] <- mean(recruit_grow_iter[x,e,s,], na.rm = T)
      recruit_grow_mean[x,e,s,2:3] <- quantile(recruit_grow_iter[x,e,s,], probs = c(.1,.9), na.rm = T)
      
      recruit_flw_mean[x,e,s,1] <- mean(recruit_flw_iter[x,e,s,], na.rm = T)
      recruit_flw_mean[x,e,s,2:3] <- quantile(recruit_flw_iter[x,e,s,], probs = c(.1,.9), na.rm = T)
      
      recruit_fert_mean[x,e,s,1] <- mean(recruit_fert_iter[x,e,s,], na.rm = T)
      recruit_fert_mean[x,e,s,2:3] <- quantile(recruit_fert_iter[x,e,s,], probs = c(.1,.9), na.rm = T)
      
      recruit_spike_mean[x,e,s,1] <- mean(recruit_spike_iter[x,e,s,], na.rm = T)
      recruit_spike_mean[x,e,s,2:3] <- quantile(recruit_spike_iter[x,e,s,], probs = c(.1,.9), na.rm = T)
    }
  }
}

dimnames(recruit_surv_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list, c("mean","twenty","eighty"))
dimnames(recruit_grow_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list, c("mean","twenty","eighty"))
dimnames(recruit_flw_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list, c("mean","twenty","eighty"))
dimnames(recruit_fert_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list, c("mean","twenty","eighty"))
dimnames(recruit_spike_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list, c("mean","twenty","eighty"))


        
# Now the same thing for each year specific vital rate
for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      for(y in 1:14){
      
      recruit_survyear_iter[,e,s,y,i] <- sx(make_params_quadXorigin(species=s,
                                                endo_mean=(e-1),
                                                endo_var=(e-1),
                                                original = 1, # should be =1 to represent recruit
                                                draw=post_draws[i],
                                                max_size=max_size,
                                                rfx=T,
                                                year = y,
                                                repro_offset = 0,
                                                surv_par=surv_par,
                                                surv_sdlg_par = surv_sdlg_par,
                                                grow_par=grow_par,
                                                grow_sdlg_par = grow_sdlg_par,
                                                flow_par=flow_par,
                                                fert_par=fert_par,
                                                spike_par=spike_par,
                                                seed_par=seed_par,
                                                recruit_par=recruit_par), x = x_seq[,s])
      recruit_growyear_iter[,e,s,y,i] <- gx(make_params_quadXorigin(species=s,
                                          endo_mean=(e-1),
                                          endo_var=(e-1),
                                          original = 1, # should be =1 to represent recruit
                                          draw=post_draws[i],
                                          max_size=max_size,
                                          rfx=T,
                                          year = y,
                                          repro_offset = 0,
                                          surv_par=surv_par,
                                          surv_sdlg_par = surv_sdlg_par,
                                          grow_par=grow_par,
                                          grow_sdlg_par = grow_sdlg_par,
                                          flow_par=flow_par,
                                          fert_par=fert_par,
                                          spike_par=spike_par,
                                          seed_par=seed_par,
                                          recruit_par=recruit_par), x = x_seq[,s])
      recruit_flwyear_iter[,e,s,y,i] <- flwx(make_params_quadXorigin(species=s,
                                           endo_mean=(e-1),
                                           endo_var=(e-1),
                                           original = 1, # should be =1 to represent recruit
                                           draw=post_draws[i],
                                           max_size=max_size,
                                           rfx=T,
                                           year = y,
                                           repro_offset = 0,
                                           surv_par=surv_par,
                                           surv_sdlg_par = surv_sdlg_par,
                                           grow_par=grow_par,
                                           grow_sdlg_par = grow_sdlg_par,
                                           flow_par=flow_par,
                                           fert_par=fert_par,
                                           spike_par=spike_par,
                                           seed_par=seed_par,
                                           recruit_par=recruit_par), x = x_seq[,s])
      recruit_fertyear_iter[,e,s,y,i] <- fertx(make_params_quadXorigin(species=s,
                                             endo_mean=(e-1),
                                             endo_var=(e-1),
                                             original = 1, # should be =1 to represent recruit
                                             draw=post_draws[i],
                                             max_size=max_size,
                                             rfx=T,
                                             year = y,
                                             repro_offset = 0,
                                             surv_par=surv_par,
                                             surv_sdlg_par = surv_sdlg_par,
                                             grow_par=grow_par,
                                             grow_sdlg_par = grow_sdlg_par,
                                             flow_par=flow_par,
                                             fert_par=fert_par,
                                             spike_par=spike_par,
                                             seed_par=seed_par,
                                             recruit_par=recruit_par), x = x_seq[,s])
      recruit_spikeyear_iter[,e,s,y,i] <- spikex(make_params_quadXorigin(species=s,
                                                   endo_mean=(e-1),
                                                   endo_var=(e-1),
                                                   original = 1, # should be =1 to represent recruit
                                                   draw=post_draws[i],
                                                   max_size=max_size,
                                                   rfx=T,
                                                   year = y,
                                                   repro_offset = 0,
                                                   surv_par=surv_par,
                                                   surv_sdlg_par = surv_sdlg_par,
                                                   grow_par=grow_par,
                                                   grow_sdlg_par = grow_sdlg_par,
                                                   flow_par=flow_par,
                                                   fert_par=fert_par,
                                                   spike_par=spike_par,
                                                   seed_par=seed_par,
                                                   recruit_par=recruit_par), x = x_seq[,s])
      
      }
    }
  }
}
        
for(x in 1:length(x_seq[,1])){
  for(e in 1:2){
    for(s in 1:7){
      for(y in 1:14){
        recruit_survyear_mean[x,e,s,y,1] <- mean(recruit_survyear_iter[x,e,s,y,], na.rm = T)
      recruit_survyear_mean[x,e,s,y,2:3] <- quantile(recruit_survyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)

      recruit_growyear_mean[x,e,s,y,1] <- mean(recruit_growyear_iter[x,e,s,y,], na.rm = T)
      recruit_growyear_mean[x,e,s,y,2:3] <- quantile(recruit_growyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)

      recruit_flwyear_mean[x,e,s,y,1] <- mean(recruit_flwyear_iter[x,e,s,y,], na.rm = T)
      recruit_flwyear_mean[x,e,s,y,2:3] <- quantile(recruit_flwyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)

      recruit_fertyear_mean[x,e,s,y,1] <- mean(recruit_fertyear_iter[x,e,s,y,], na.rm = T)
      recruit_fertyear_mean[x,e,s,y,2:3] <- quantile(recruit_fertyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)
      
      recruit_spikeyear_mean[x,e,s,y,1] <- mean(recruit_spikeyear_iter[x,e,s,y,], na.rm = T)
      recruit_spikeyear_mean[x,e,s,y,2:3] <- quantile(recruit_spikeyear_iter[x,e,s,y,], probs = c(.1,.9), na.rm = T)
      }
    }
  }
}

dimnames(recruit_survyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list,c(2008:2021), c("mean","twenty","eighty"))
dimnames(recruit_growyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list,c(2008:2021), c("mean","twenty","eighty"))
dimnames(recruit_flwyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list,c(2008:2021), c("mean","twenty","eighty"))
dimnames(recruit_fertyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list,c(2008:2021), c("mean","twenty","eighty"))
dimnames(recruit_spikeyear_mean) <- list(paste0("size", 1:length(x_seq[,1])), c("S-","S+"), species_code_list,c(2008:2021), c("mean","twenty","eighty"))



#Now I'm gonna make these into  tidy dataframes for plotting

dimnames(x_seq) <- list(paste0("x_size", 1:length(x_seq[,1])), paste0(species_code_list))

x_seq_df <- as_tibble(x_seq) %>%
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = -no_row,
               names_to = "Species",
               values_to = "x_seq") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  mutate(log_x_seq = log(x_seq))


recruit_surv_mean_df <- as_tibble(recruit_surv_mean)  %>%   
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("surv") ) %>% 
  separate(name, c("Endo", "Species","quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "surv") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)


recruit_survyear_mean_df <- as_tibble(recruit_survyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = -no_row,
               values_to = c("surv")) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "surv") %>%
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)



recruit_grow_mean_df <- as_tibble(recruit_grow_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("size_t1") )  %>% 
  separate(name, c("Endo", "Species","quantile"), "\\.") %>%   
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "size_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)

recruit_growyear_mean_df <- as_tibble(recruit_growyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("size_t1") ) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "size_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)

recruit_flw_mean_df <- as_tibble(recruit_flw_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("flw_t1") ) %>% 
  separate(name, c("Endo", "Species","quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "flw_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)

recruit_flwyear_mean_df <- as_tibble(recruit_flwyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("flw_t1") ) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "flw_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)

recruit_fert_mean_df <- as_tibble(recruit_fert_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("fert_t1") ) %>% 
  separate(name, c("Endo", "Species","quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "fert_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)

recruit_fertyear_mean_df <- as_tibble(recruit_fertyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("fert_t1") ) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "fert_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)

recruit_spike_mean_df <- as_tibble(recruit_spike_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("spike_t1") ) %>% 
  separate(name, c("Endo", "Species","quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species"), names_from = "quantile", values_from = "spike_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df) 

recruit_spikeyear_mean_df <- as_tibble(recruit_spikeyear_mean) %>%    
  mutate(no_row = row_number()) %>% 
  pivot_longer(cols = starts_with("S"),
               values_to = c("spike_t1") ) %>% 
  separate(name, c("Endo", "Species", "Year", "quantile"), "\\.") %>% 
  pivot_wider(id_cols = c("no_row", "Endo", "Species", "Year"), names_from = "quantile", values_from = "spike_t1") %>% 
  mutate(Species = recode(Species, AGPE = species_list[1], ELRI = species_list[2], ELVI = species_list[3], FESU = species_list[4], LOAR = species_list[5], POAL = species_list[6], POSY = species_list[7])) %>% 
  left_join(x_seq_df)
# Bin Data by size and then by year for plotting for only recruit data


recruit_surv_sizebin <- bin_by_size_t(filter(LTREB_data_forsurv, origin_01 == 1),vr = surv_t1, nbins = 20)
recruit_grow_sizebin <- bin_by_size_t(filter(LTREB_data_forgrow, origin_01 == 1),vr = size_t1, nbins = 20)
recruit_flw_sizebin <- bin_by_size_t1(filter(LTREB_data_forflw, origin_01 == 1), vr = FLW_STAT_T1,nbins = 20)
recruit_fert_sizebin <- bin_by_size_t1(filter(LTREB_data_forfert, origin_01 == 1),vr = FLW_COUNT_T1, nbins = 20)
recruit_spike_sizebin <- bin_by_size_t1(filter(LTREB_data_forspike, origin_01 ==1), vr = spike_count_t1, nbins = 20)

recruit_surv_yearsizebin <- bin_by_year_size_t(filter(LTREB_data_forsurv,origin_01==1),vr = surv_t1, nbins = 20)
recruit_grow_yearsizebin <- bin_by_year_size_t(filter(LTREB_data_forgrow,origin_01==1),vr = size_t1, nbins = 20)
recruit_flw_yearsizebin <- bin_by_year_size_t1(filter(LTREB_data_forflw,origin_01==1), vr = FLW_STAT_T1,nbins = 20)
recruit_fert_yearsizebin <- bin_by_year_size_t1(filter(LTREB_data_forfert,origin_01==1),vr = FLW_COUNT_T1, nbins = 20)
recruit_spike_yearsizebin <- bin_by_year_size_t1(filter(LTREB_data_forspike,origin_01==1), vr = spike_count_t1, nbins = 20)


recruit_surv_neatdata <- neat_names(filter(LTREB_data_forsurv, origin_01 == 1))
recruit_grow_neatdata <- neat_names(filter(LTREB_data_forgrow, origin_01 == 1))
recruit_flw_neatdata <- neat_names(filter(LTREB_data_forflw, origin_01 == 1))
recruit_fert_neatdata <- neat_names(filter(LTREB_data_forfert, origin_01 == 1))
recruit_spike_neatdata <- neat_names(filter(LTREB_data_forspike, origin_01 == 1))


#The plots
recruit_surv_meanplot <- ggplot()+
  geom_point(data = recruit_surv_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo, col = Species), alpha = .7) +
  geom_ribbon(data = recruit_surv_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, alpha = Endo, fill = Species))+
  geom_line(data = recruit_surv_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = Species)) +
  scale_color_manual(values = species_colors)+  scale_fill_manual(values = species_colors)+
  scale_shape_manual( values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + scale_alpha_manual( values = c(.1,.3))+
  facet_wrap(~Species, scales = "free", ncol = 2) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none", col = "none", alpha = guide_legend(order = 1), linetype = guide_legend(order = 1), shape = guide_legend(order = 1), size = guide_legend(order = 2))+
  labs(title = "Adult Survival",subtitle = "Data from Recruit plants only", y = "Survival Probability", x = expression("log(# of tillers)"[" year t"]), size = "Sample Size")
# recruit_surv_meanplot
ggsave(recruit_surv_meanplot, filename = "recruit_surv_meanplot_quadXorigin.png", width = 8, height = 12)



recruit_surv_yearplot <- ggplot()+
  geom_point(data = recruit_surv_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .4) +
  # geom_ribbon(data = recruit_survyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = recruit_survyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = as.factor(Year), group = Year), alpha = .8) +
  scale_shape_manual( values = c(1,19)) + scale_linetype_manual(values = c(2,1)) +
  scale_color_manual( values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free",ncol = 4) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none",shape = guide_legend(order = 2, ncol = 1),linetype = guide_legend(order = 2, ncol = 1), color = guide_legend(order = 1, ncol = 2), size = guide_legend(order = 3))+
  labs(title = "Adult Survival", subtitle = "Data from Recruit plants only", y = "Survival Probability", x = expression("log(# of tillers)"[" year t"]), size = "Sample Size", col = "Year")
# surv_yearplot
ggsave(recruit_surv_yearplot, filename = "recruit_surv_yearplot_quadXorigin.png", width = 16, height = 17)


recruit_grow_meanplot <- ggplot()+
  # geom_point(data = recruit_grow_neatdata, aes(x = logsize_t, y = size_t1, shape = Endo, col = Species), alpha = .1) +
  geom_point(data = recruit_grow_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo, col = Species), alpha = .7) +
  geom_ribbon(data = recruit_grow_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, alpha = Endo, fill = Species))+
  geom_line(data = recruit_grow_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = Species)) +
  scale_color_manual(values = species_colors)+  scale_fill_manual(values = species_colors)+
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + scale_alpha_manual(values = c(.1,.3))+
  facet_wrap(~Species, scales = "free", ncol = 2) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none", col = "none", alpha = guide_legend(order = 1), linetype = guide_legend(order = 1), shape = guide_legend(order = 1), size = guide_legend(order = 2))+
  labs(title = "Adult Growth", subtitle = "Data from Recruit plants only", y = "# of Tillers", x = expression("log(# of tillers)"[" year t"]), size = "Sample Size")
# recruit_grow_meanplot
ggsave(recruit_grow_meanplot, filename = "recruit_grow_meanplot_quadXorigin.png", width = 8, height = 12)



recruit_grow_yearplot <- ggplot()+
  # geom_point(data = recruit_grow_neatdata, aes(x = logsize_t, y = size_t1, col = as.factor(Year), shape = Endo), alpha = .4) +
  geom_point(data = recruit_grow_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .4) +
  # geom_ribbon(data = recruit_growyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = recruit_growyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = as.factor(Year),group = Year), alpha = .8) +
  scale_shape_manual( values = c(1,19)) + scale_linetype_manual(values = c(2,1)) +
  scale_color_manual( values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free", ncol = 4) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none",shape = guide_legend(order = 2, ncol = 1),linetype = guide_legend(order = 2, ncol = 1), color = guide_legend(order = 1, ncol = 2), size = guide_legend(order = 3))+
  labs(title = "Adult Growth", subtitle = "Data from Recruit plants only",  y = "# of Tillers", x = expression("log(# of tillers)"[" year t"]), col = "Year",size = "Sample Size")
# recruit_grow_yearplot
ggsave(recruit_grow_yearplot, filename = "recruit_grow_yearplot_quadXorigin.png", width = 16, height = 17)




recruit_flw_meanplot <- ggplot()+
  geom_point(data = recruit_flw_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo, col = Species), alpha = .7) +
  geom_ribbon(data = recruit_flw_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, alpha = Endo, fill = Species))+
  geom_line(data = recruit_flw_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = Species)) +
  scale_color_manual(values = species_colors)+  scale_fill_manual(values = species_colors)+
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + scale_alpha_manual(values = c(.1,.3))+  
  facet_wrap(~Species, scales = "free", ncol = 2) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none", col = "none", alpha = guide_legend(order = 1), linetype = guide_legend(order = 1), shape = guide_legend(order = 1), size = guide_legend(order = 2))+
  labs(title = "Flowering", subtitle = "Data from Recruit plants only", y = "Flowering Probability", x = expression("log(# of tillers)"[" year t+1"]), size = "Sample Size")
# recruit_flw_meanplot
ggsave(recruit_flw_meanplot, filename = "recruit_flw_meanplot_quadXorigin.png", width = 8, height = 12)


recruit_flw_yearplot <- ggplot()+
  geom_point(data = recruit_flw_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .4) +
  # geom_ribbon(data = flwyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = recruit_flwyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col= as.factor(Year), group = Year)) +
  scale_shape_manual( values = c(1,19)) + scale_linetype_manual(values = c(2,1)) +
  scale_color_manual( values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free", ncol = 4) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none",shape = guide_legend(order = 2, ncol = 1),linetype = guide_legend(order = 2, ncol = 1), color = guide_legend(order = 1, ncol = 2), size = guide_legend(order = 3))+
  labs(title = "Flowering",subtitle = "Data from Recruit plants only", y = "Flowering Probability", x = expression("log(# of tillers)"[" year t+1"]), col = "Year", size = "Sample Size")
# flw_yearplot
ggsave(recruit_flw_yearplot, filename = "recruit_flw_yearplot_quadXorigin.png", width = 16, height = 17)

recruit_fert_meanplot <- ggplot()+
  geom_point(data = recruit_fert_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo, col = Species), alpha = .7) +
  geom_ribbon(data = recruit_fert_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, alpha = Endo, fill = Species))+
  geom_line(data = recruit_fert_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, color = Species)) +
  scale_color_manual(values = species_colors)+  scale_fill_manual(values = species_colors)+
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + scale_alpha_manual(values = c(.1,.3))+  
  facet_wrap(~Species, scales = "free", ncol = 2) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none", col = "none", alpha = guide_legend(order = 1), linetype = guide_legend(order = 1), shape = guide_legend(order = 1), size = guide_legend(order = 2))+
  labs(title = "Fertility",subtitle = "Data from Recruit plants only", y = "# of Repro. Tillers", x = expression("log(# of tillers)"[" year t+1"]), size = "Sample Size")
# fert_meanplot
ggsave(recruit_fert_meanplot, filename = "recruit_fert_meanplot_quadXorigin.png", width = 8, height = 12)

recruit_fert_yearplot <- ggplot()+
  geom_point(data = recruit_fert_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .4) +
  # geom_ribbon(data = fertyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = recruit_fertyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = as.factor(Year), group = Year)) +
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + 
  scale_color_manual(values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free", ncol = 4) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none",shape = guide_legend(order = 2, ncol = 1),linetype = guide_legend(order = 2, ncol = 1), color = guide_legend(order = 1, ncol = 2), size = guide_legend(order = 3))+
  labs(title = "Fertility", subtitle = "Data from Recruit plants only", y = "# of Repro. Tillers", x = expression("log(# of tillers)"[" year t+1"]), col = "Year", size = "Sample Size")
# fert_yearplot
ggsave(recruit_fert_yearplot, filename = "recruit_fert_yearplot_quadXorigin.png",width = 16, height = 17)


recruit_spike_meanplot <- ggplot()+
  # geom_point(data = spike_neatdata, aes(x = logsize_t, y = spike_count_t1, shape = Endo, col = Species), alpha = .4) +
  geom_point(data = recruit_spike_sizebin, aes(x = mean_size, y = mean_vr, size = samplesize, shape = Endo, col = Species), alpha = .7) +
  geom_ribbon(data = recruit_spike_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, alpha = Endo, fill = Species))+
  geom_line(data = recruit_spike_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, color = Species)) +
  scale_color_manual(values = species_colors)+  scale_fill_manual(values = species_colors)+
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + scale_alpha_manual(values = c(.1,.3))+  
  facet_wrap(~Species, scales = "free", ncol = 2) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none", col = "none", alpha = guide_legend(order = 1), linetype = guide_legend(order = 1), shape = guide_legend(order = 1), size = guide_legend(order = 2))+
  labs(title = "Spikes/Infl.",subtitle = "Data from Recruit plants only", y = "Spikelet/Infl.", x = expression("log(# of tillers)"[" year t+1"]), size = "Sample Size")
# spike_meanplot
ggsave(recruit_spike_meanplot, filename = "recruit_spike_meanplot_quadXorigin.png", width = 8, height = 12)

recruit_spike_yearplot <- ggplot()+
  # geom_point(data = spike_neatdata, aes(x = logsize_t, y = spike_count_t1, col = as.factor(Year), shape = Endo), alpha = .4) +
  geom_point(data = recruit_spike_yearsizebin, aes(x = mean_size, y = mean_vr, size = samplesize, col = as.factor(Year), shape = Endo), alpha = .4) +
  # geom_ribbon(data = spikeyear_mean_df, aes(x = log_x_seq, ymin = twenty, ymax = eighty, fill = Endo), alpha = .3)+
  geom_line(data = recruit_spikeyear_mean_df, aes(x = log_x_seq, y = mean, linetype = Endo, col = as.factor(Year), group = Year)) +
  scale_shape_manual(values = c(1,19)) + scale_linetype_manual(values = c(2,1)) + 
  scale_color_manual(values = yearcolors)+
  facet_wrap(~Species + Endo, scales = "free", ncol = 4) + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = rel(1.5)),
        legend.position = "right",
        legend.justification = "center",
        legend.margin = margin(10,0,10,0),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(.8)),
        axis.text = element_text(size = rel(1)),
        title = element_text(size = rel(2))) + 
  guides(fill = "none",shape = guide_legend(order = 2, ncol = 1),linetype = guide_legend(order = 2, ncol = 1), color = guide_legend(order = 1, ncol = 2), size = guide_legend(order = 3))+
  labs(title = "Spikes/Infl.",subtitle = "Data from Recruit plants only", y = "Spikelet/Infl.", x = expression("log(# of tillers)"[" year t+1"]), col = "Year", size = "Sample Size")
# spike_yearplot
ggsave(recruit_spike_yearplot, filename = "recruit_spike_yearplot_quadXorigin.png",width = 16, height = 17)












 ######## making a heatmap of each vital rate for degree of buffering #####



# Pulling out the actual parameters
surv_par <- rstan::extract(surv_fit, pars =quote_bare(beta0,betasize,betaendo,
                                                      tau_year, tau_plot, sigma_year, sigma0, sigmaendo))
surv_sdlg_par <- rstan::extract(surv_fit_seedling, pars =quote_bare(beta0,betaendo,
                                                                    tau_year, tau_plot, sigma_year, sigma0, sigmaendo))
grow_par <- rstan::extract(grow_fit, pars = quote_bare(beta0,betasize,betaendo,
                                                       tau_year, tau_plot,
                                                       sigma, sigma_year, sigma0, sigmaendo))
grow_sdlg_par <- rstan::extract(grow_fit_seedling, pars = quote_bare(beta0,betaendo,
                                                                     tau_year, tau_plot,
                                                                     sigma, sigma_year, sigma0, sigmaendo))
flow_par <- rstan::extract(flw_fit, pars = quote_bare(beta0,betasize,betaendo,
                                                      tau_year, tau_plot, sigma_year, sigma0, sigmaendo))
fert_par <- rstan::extract(fert_fit, pars = quote_bare(beta0,betasize,betaendo,
                                                       tau_year, tau_plot, sigma_year, sigma0, sigmaendo))
spike_par <- rstan::extract(spike_fit, pars = quote_bare(beta0,betasize,betaendo,
                                                         tau_year, tau_plot,
                                                         phi, sigma_year, sigma0, sigmaendo))
seed_par <- rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
recruit_par <- rstan::extract(stos_fit, pars = quote_bare(beta0,betaendo,
                                                          tau_year, tau_plot, sigma_year, sigma0, sigmaendo))


# First I'm making dataframes for the endophytes effects on mean and variance
#surv
dimnames(surv_par$sigmaendo) <- list(Draw = paste0("i",1:dim(surv_par$sigma_year)[1]), Species = species_list)
surv_sigmaendo_cube <- cubelyr::as.tbl_cube(surv_par$sigmaendo)
surv_sigmaendo_df <- as_tibble(surv_sigmaendo_cube)  %>% 
  rename(estimate = `surv_par$sigmaendo`) %>% 
  mutate(vital_rate = "Survival", effect = "Standard Deviation")
dimnames(surv_par$betaendo) <- list(Draw = paste0("i",1:dim(surv_par$sigma_year)[1]), Species = species_list)
surv_betaendo_cube <- cubelyr::as.tbl_cube(surv_par$betaendo)
surv_betaendo_df <- as_tibble(surv_betaendo_cube)  %>% 
  rename(estimate = `surv_par$betaendo`) %>% 
  mutate(vital_rate = "Survival", effect = "Mean")
#seedling survival
dimnames(surv_sdlg_par$sigmaendo) <- list(Draw = paste0("i",1:dim(surv_sdlg_par$sigma_year)[1]), Species = species_list)
surv_sdlg_sigmaendo_cube <- cubelyr::as.tbl_cube(surv_sdlg_par$sigmaendo)
surv_sdlg_sigmaendo_df <- as_tibble(surv_sdlg_sigmaendo_cube)  %>% 
  rename(estimate = `surv_sdlg_par$sigmaendo`) %>% 
  mutate(vital_rate = "Seedling Survival", effect = "Standard Deviation")
dimnames(surv_sdlg_par$betaendo) <- list(Draw = paste0("i",1:dim(surv_sdlg_par$sigma_year)[1]), Species = species_list)
surv_sdlg_betaendo_cube <- cubelyr::as.tbl_cube(surv_sdlg_par$betaendo)
surv_sdlg_betaendo_df <- as_tibble(surv_sdlg_betaendo_cube)  %>% 
  rename(estimate = `surv_sdlg_par$betaendo`) %>% 
  mutate(vital_rate = "Seedling Survival", effect = "Mean")
#grow
dimnames(grow_par$sigmaendo) <- list(Draw = paste0("i",1:dim(grow_par$sigmaendo)[1]), Species = species_list)
grow_sigmaendo_cube <- cubelyr::as.tbl_cube(grow_par$sigmaendo)
grow_sigmaendo_df <- as_tibble(grow_sigmaendo_cube)  %>% 
  rename(estimate = `grow_par$sigmaendo`)%>% 
  mutate(vital_rate = "Growth", effect = "Standard Deviation")
dimnames(grow_par$betaendo) <- list(Draw = paste0("i",1:dim(grow_par$betaendo)[1]), Species = species_list)
grow_betaendo_cube <- cubelyr::as.tbl_cube(grow_par$betaendo)
grow_betaendo_df <- as_tibble(grow_betaendo_cube)  %>% 
  rename(estimate = `grow_par$betaendo`)%>% 
  mutate(vital_rate = "Growth", effect = "Mean")
#seedling grow
dimnames(grow_sdlg_par$sigmaendo) <- list(Draw = paste0("i",1:dim(grow_sdlg_par$sigmaendo)[1]), Species = species_list)
seedgrow_sigmaendo_cube <- cubelyr::as.tbl_cube(grow_sdlg_par$sigmaendo)
seedgrow_sigmaendo_df <- as_tibble(seedgrow_sigmaendo_cube)  %>% 
  rename(estimate = `grow_sdlg_par$sigmaendo`)%>% 
  mutate(vital_rate = "Seedling Growth", effect = "Standard Deviation")
dimnames(grow_sdlg_par$betaendo) <- list(Draw = paste0("i",1:dim(grow_sdlg_par$betaendo)[1]), Species = species_list)
seedgrow_betaendo_cube <- cubelyr::as.tbl_cube(grow_sdlg_par$betaendo)
seedgrow_betaendo_df <- as_tibble(seedgrow_betaendo_cube)  %>% 
  rename(estimate = `grow_sdlg_par$betaendo`)%>% 
  mutate(vital_rate = "Seedling Growth", effect = "Mean")
#flw
dimnames(flow_par$sigmaendo) <- list(Draw = paste0("i",1:dim(flow_par$sigmaendo)[1]), Species = species_list)
flow_sigmaendo_cube <- cubelyr::as.tbl_cube(flow_par$sigmaendo)
flow_sigmaendo_df <- as_tibble(flow_sigmaendo_cube)  %>% 
  rename(estimate = `flow_par$sigmaendo`)%>% 
  mutate(vital_rate = "Flowering", effect = "Standard Deviation")
dimnames(flow_par$betaendo) <- list(Draw = paste0("i",1:dim(flow_par$betaendo)[1]), Species = species_list)
flow_betaendo_cube <- cubelyr::as.tbl_cube(flow_par$betaendo)
flow_betaendo_df <- as_tibble(flow_betaendo_cube)  %>% 
  rename(estimate = `flow_par$betaendo`)%>% 
  mutate(vital_rate = "Flowering", effect = "Mean")
#fert
dimnames(fert_par$sigmaendo) <- list(Draw = paste0("i",1:dim(fert_par$sigmaendo)[1]), Species = species_list)
fert_sigmaendo_cube <- cubelyr::as.tbl_cube(fert_par$sigmaendo)
fert_sigmaendo_df <- as_tibble(fert_sigmaendo_cube)  %>% 
  rename(estimate = `fert_par$sigmaendo`)%>% 
  mutate(vital_rate = "Inflorescence Production", effect = "Standard Deviation")
dimnames(fert_par$betaendo) <- list(Draw = paste0("i",1:dim(fert_par$betaendo)[1]), Species = species_list)
fert_betaendo_cube <- cubelyr::as.tbl_cube(fert_par$betaendo)
fert_betaendo_df <- as_tibble(fert_betaendo_cube)  %>% 
  rename(estimate = `fert_par$betaendo`)%>% 
  mutate(vital_rate = "Inflorescence Production", effect = "Mean")
#spike
dimnames(spike_par$sigmaendo) <- list(Draw = paste0("i",1:dim(spike_par$sigmaendo)[1]), Species = species_list)
spike_sigmaendo_cube <- cubelyr::as.tbl_cube(spike_par$sigmaendo)
spike_sigmaendo_df <- as_tibble(spike_sigmaendo_cube)  %>% 
  rename(estimate = `spike_par$sigmaendo`)%>% 
  mutate(vital_rate = "Spikelets/Infl.", effect = "Standard Deviation")
dimnames(spike_par$betaendo) <- list(Draw = paste0("i",1:dim(spike_par$betaendo)[1]), Species = species_list)
spike_betaendo_cube <- cubelyr::as.tbl_cube(spike_par$betaendo)
spike_betaendo_df <- as_tibble(spike_betaendo_cube)  %>% 
  rename(estimate = `spike_par$betaendo`)%>% 
  mutate(vital_rate = "Spikelets/Infl.", effect = "Mean")
#germination
dimnames(recruit_par$sigmaendo) <- list(Draw = paste0("i",1:dim(recruit_par$sigmaendo)[1]), Species = species_list)
recruit_sigmaendo_cube <- cubelyr::as.tbl_cube(recruit_par$sigmaendo)
recruit_sigmaendo_df <- as_tibble(recruit_sigmaendo_cube)  %>% 
  rename(estimate = `recruit_par$sigmaendo`)%>% 
  mutate(vital_rate = "Recruitment", effect = "Standard Deviation")
dimnames(recruit_par$betaendo) <- list(Draw = paste0("i",1:dim(recruit_par$betaendo)[1]), Species = species_list)
recruit_betaendo_cube <- cubelyr::as.tbl_cube(recruit_par$betaendo)
recruit_betaendo_df <- as_tibble(recruit_betaendo_cube)  %>% 
  rename(estimate = `recruit_par$betaendo`)%>% 
  mutate(vital_rate = "Recruitment", effect = "Mean")

#Combining all of those into one dataframe
endo_vr_effects_df <- surv_sigmaendo_df %>% 
  rbind(surv_betaendo_df,
        surv_sdlg_sigmaendo_df,surv_sdlg_betaendo_df,
        grow_sigmaendo_df,grow_betaendo_df,
        seedgrow_sigmaendo_df,seedgrow_betaendo_df,
        flow_sigmaendo_df, flow_betaendo_df,
        fert_sigmaendo_df,fert_betaendo_df,
        spike_sigmaendo_df, spike_betaendo_df,
        recruit_sigmaendo_df,recruit_betaendo_df)
# calculating the average effects for each species and vital rate
endo_vr_effects_summary <- endo_vr_effects_df %>% 
  group_by(Species, vital_rate, effect) %>% 
  summarize(average_effect = mean(estimate))

endo_vr_effects_overallsd <- endo_vr_effects_df %>% 
  group_by(effect) %>% 
  summarize(sd = sd(estimate)) %>% 
  pivot_wider(names_from = c(effect), values_from = c(sd))

endo_vr_effects_standardized <- endo_vr_effects_summary %>% 
  mutate(stand_effect = case_when(effect == "Mean" ~ average_effect/endo_vr_effects_overallsd$Mean,
                                  effect == "Standard Deviation" ~ average_effect/endo_vr_effects_overallsd$`Standard Deviation`))
# some overall summary numbers for manuscript
summary_endo_vr_effects_standardized <- endo_vr_effects_standardized %>% 
  group_by(effect) %>% 
  summarize(mean_average = mean(average_effect),
            mean_stand = mean(stand_effect))

#now we can make a heat map based on those means
vr_order <- c("Survival","Seedling Survival", "Growth", "Seedling Growth", "Flowering", "Inflorescence Production", "Spikelets/Infl.", "Recruitment")
meanvar_effect_heatmap <- ggplot()+
  geom_tile(data = endo_vr_effects_summary, aes(x = Species, y = factor(vital_rate, levels=vr_order), fill = average_effect), color = "lightgrey")+
  scale_fill_gradient2(low = "#ef8a62", high = "#67a9cf")+
  facet_wrap(~effect)+
  labs(x = "Species", y = "Vital Rate", fill = "Avg. Endophyte Effect")+  
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

meanvar_effect_heatmap

#Version filling by the standardized effect size
stand_effect_heatmap <- ggplot()+
  geom_tile(data = endo_vr_effects_standardized, aes(x = Species, y = factor(vital_rate, levels=vr_order), fill = stand_effect), color = "lightgrey")+
  geom_text(data = endo_vr_effects_standardized, aes(x = Species, y = factor(vital_rate, levels = vr_order), label = round(stand_effect, digits = 2)), size = 3)+
  scale_fill_gradient2(low = "#ef8a62", high = "#67a9cf")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 10))+
  facet_wrap(~effect)+
  labs(x = "Species", y = "Vital Rate", fill = "Standardized Effect")+  
  theme_minimal()+
  theme(text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1))+
  guides(fill = "none")
stand_effect_heatmap

ggsave(stand_effect_heatmap, file = "stand_effect_heatmap_text_quadXorigin.png", width = 7.5, height = 5)


par(mfrow=c(4,1),mar=c(2,4,0,0))
hist(filter(endo_vr_effects_standardized, effect == "Mean")$average_effect)
hist(filter(endo_vr_effects_standardized, effect == "Mean")$stand_effect)
hist(filter(endo_vr_effects_standardized, effect == "Variance")$average_effect)
hist(filter(endo_vr_effects_standardized, effect == "Variance")$stand_effect)
#pulling out 100 random draws from each vital rate/species/effect
thinned_endo_vr_effects_df <- endo_vr_effects_df %>% 
  group_by(vital_rate, effect, Species) %>% 
  slice_sample(n=100)

meanvar_effect_pointmap <- ggplot()+
  geom_tile(data = filter(endo_vr_effects_summary), aes(x = Species, y = factor(vital_rate, levels=vr_order), fill = average_effect), color = "lightgrey")+
  geom_jitter(data = filter(thinned_endo_vr_effects_df), aes(x = Species, y = factor(vital_rate, levels=vr_order), color = estimate), alpha = .8)+
  scale_fill_gradient2(low = "#ef8a62", high = "#67a9cf")+
  scale_color_gradient2(low = "#ef8a62",high="#67a9cf")+
  facet_wrap(~effect)+
  labs(x = "Species", y = "Vital Rate", fill = "Avg. Endophyte Effect")+  
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
meanvar_effect_pointmap
