## Title: Climate explicit Grass endophyte vital rate models with a Bayesian framework
## Purpose: Assembles data lists and runs vital rate kernels written in STAN with mixed effects, 
## and does visualisation of posterior predictive checks
## Authors: Joshua and Tom
#############################################################

library(tidyverse)
library(rstan)
library(StanHeaders)
# library(shinystan)
library(bayesplot)
library(devtools)
library(moments)
library(gridExtra)

library(RColorBrewer)


invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }
Lkurtosis=function(x) log(kurtosis(x)); 


#############################################################################################
####### Data manipulation to prepare data for Stan models------------------
#############################################################################################

#  data are prepared in the endodemog_data_processing.R file, 
# source("Analyses/endodemog_data_processing.R")
LTREB_full <- read_csv("~/Dropbox/EndodemogData/Fulldataplusmetadata/LTREB_full.csv")


#############################################################################################
####### Preparing data lists for vital rate kernels ------------------
#############################################################################################
# I'm gonna start with just one species
## Clean up the main data frame for NA's, other small data entry errors
LTREB_data_forsurv <- LTREB_full %>% 
  # filter(species == "FESU", year_t1 < 2019) %>%
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(endo_01)) %>%   # There are a few LOAR that don't have a plot level endo assigned
  filter(origin_01 == 1 & year_t != birth | origin_01 == 0) %>%    # filtering out first year germinants (including those that are bigger than 1 tiller)
  mutate(origin_index = origin_01+1) %>% 
  filter(!is.na(spei12)) # I think I can try to get more recent climate data, but there are NA's for the 2020 and 2019 AGPE data right now
dim(LTREB_data_forsurv)

LTREB_surv_means <- LTREB_data_forsurv %>% 
  group_by(species, endo_01, year_t1) %>% 
  summarize(mean_surv = mean(surv_t1),
            spei12 = mean(as.numeric(spei12)),
            annual_precip = mean(as.numeric(annual_precip)),
            annual_temp = mean(as.numeric(annual_temp)),
            count = n())
# 
# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = year_t1, y = mean_surv, color = as.factor(endo_01), size = count))+ facet_wrap(~species)+
#   theme_classic()
# # 
# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = year_t1, y = spei12, color = endo_01, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = year_t1, y = annual_precip, color = endo_01, size = count))+ facet_wrap(~species)

# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = year_t1, y = annual_temp, color = endo_01, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = spei12, y = mean_surv, color = as.factor(endo_01), size = count), alpha = .6)+ facet_wrap(~species) +
#   theme_classic()
# 
# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = annual_precip, y = mean_surv, color = as.factor(endo_01), size = count), alpha = .6)+ facet_wrap(~species) +
#   theme_classic()
# 
# ggplot(data = LTREB_surv_means)+
#   geom_point(aes(x = annual_temp, y = mean_surv, color = as.factor(endo_01), size = count), alpha = .6)+ facet_wrap(~species) +
#   theme_classic()
# 
# 
# # 
LTREB_surv_seedling <- LTREB_full %>%
  filter(!is.na(surv_t1)) %>%
  filter(!is.na(logsize_t)) %>%
  filter(!is.na(endo_01)) %>%  # There are a few LOAR that don't have a plot level endo assigned
  filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
  filter(logsize_t == 0) # this is filtering out the plants that are "recruits" but are larger than 1 tiller
dim(LTREB_surv_seedling)
# 
LTREB_seedlingsurv_means <- LTREB_surv_seedling %>%
  group_by(species, endo_01, year_t1) %>%
  summarize(mean_surv = mean(surv_t1),
            spei12 = mean(as.numeric(spei12)),
            annual_precip = mean(as.numeric(annual_precip)),
            annual_temp = mean(as.numeric(annual_temp)),
            endo = as.character(endo_01),
            count = n())

# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = (year_t1), y = mean_surv, color = species, size = count))+ facet_wrap(~species) +
#   theme_classic()

# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = year_t1, y = spei12, color = species, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = year_t1, y = spei12, color = endo_01, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = year_t1, y = annual_precip, color = endo_01, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = year_t1, y = annual_temp, color = endo_01, size = count))+ facet_wrap(~species)
# 
# 
# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = spei12, y = mean_surv, color = as.factor(endo_01), size = count))+ facet_wrap(~species) +
#   theme_classic()


# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = annual_precip, y = mean_surv, color = as.factor(endo_01), size = count), alpha = .6)+ facet_wrap(~species) +
#   theme_classic()
# 
# ggplot(data = LTREB_seedlingsurv_means)+
#   geom_point(aes(x = annual_temp, y = mean_surv, color = as.factor(endo_01), size = count), alpha = .6)+ facet_wrap(~species) +
#   theme_classic()


# 
# # I want to look at these "seedlings" that are bigger than 1 tiller, We are dropping these from our seedling model, as they are most likely missed plants from the previous year
# LTREB_surv_big_seedling <- LTREB_full %>% 
#   filter(!is.na(surv_t1)) %>% 
#   filter(!is.na(logsize_t)) %>% 
#   filter(!is.na(endo_01)) %>%  # There are a few LOAR that don't have a plot level endo assigned
#   filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
#   filter(logsize_t != 0) 
# dim(LTREB_surv_big_seedling)
# table(LTREB_surv_big_seedling$year_t,LTREB_surv_big_seedling$size_t, LTREB_surv_big_seedling$species)
# 
LTREB_data_forflw <- LTREB_full %>%
  filter(!is.na(FLW_STAT_T1)) %>%
  filter(!is.na(logsize_t1)) %>%
  filter(!is.na(endo_01)) %>% 
  mutate(origin_index = origin_01+1)

# dim(LTREB_data_forflw)


LTREB_flw_means <- LTREB_data_forflw %>%
  group_by(species, endo_01, year_t1) %>%
  summarize(mean_flw = mean(FLW_STAT_T),
            spei12 = mean(as.numeric(spei12)),
            annual_precip = mean(as.numeric(annual_precip)),
            annual_temp = mean(as.numeric(annual_temp)),
            endo  = as.character(endo_01),
            count = n())

# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = year_t1, y = mean_flw, color = as.factor(endo_01), size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = year_t1, y = spei12, color = species, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = year_t1, y = spei12, color = endo_01, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = year_t1, y = annual_precip, color = endo_01, size = count))+ facet_wrap(~species)
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = year_t1, y = annual_temp, color = endo_01, size = count))+ facet_wrap(~species)
# 
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = spei12, y = mean_flw, color = as.factor(endo_01), size = count))+ facet_wrap(~species) +
#   theme_classic()
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = annual_precip, y = mean_flw, color = as.factor(endo_01), size = count))+ facet_wrap(~species) +
#   theme_classic()
# 
# ggplot(data = LTREB_flw_means)+
#   geom_point(aes(x = annual_temp, y = mean_flw, color = as.factor(endo_01), size = count))+ facet_wrap(~species) +
#   theme_classic()

# 
LTREB_data_forgrow <- LTREB_full %>%
  filter(!is.na(logsize_t)) %>%
  filter(!is.na(size_t1)) %>%
  filter(!is.na(endo_01)) %>%
  filter(origin_01 == 1 & year_t != birth | origin_01 == 0) %>%  # filtering out first year germinants (including those that are bigger than 1 tiller)
  mutate(origin_index = origin_01+1)

dim(LTREB_data_forgrow)

# 
LTREB_grow_seedling <- LTREB_full %>%
  filter(!is.na(logsize_t)) %>%
  filter(!is.na(size_t1)) %>%
  filter(!is.na(endo_01)) %>%
  filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
  filter(logsize_t == 0) # this is filtering out the plants that are "recruits" but are larger than 1 tiller

# dim(LTREB_grow_seedling)
# 
# 
LTREB_data_forfert <- LTREB_full %>%
  filter(!is.na(FLW_COUNT_T1)) %>%
  filter(FLW_COUNT_T1 > 0) %>%
  filter(!is.na(logsize_t1)) %>% 
  mutate(origin_index = origin_01+1)

# dim(LTREB_data_forfert)


LTREB_data_forspike <- LTREB_full %>%
  dplyr::select(-FLW_COUNT_T, -FLW_STAT_T, -SPIKE_A_T, -SPIKE_B_T, -SPIKE_C_T, -SPIKE_D_T, -SPIKE_AGPE_MEAN_T, -census_month, -year, -spei1,-spei24, -annual_temp, -annual_precip, -endo_status_from_check, -plot_endo_for_check, -endo_mismatch, -dist_a, -dist_b) %>% 
  filter(!is.na(FLW_STAT_T1)) %>% 
  filter(FLW_STAT_T1>0) %>% 
  reshape2::melt(id.var = c("plot_fixed" ,   "plot_index",         "pos"         ,           "id",
                  "species"       ,         "species_index"  ,        "endo_01",
                  "endo_index"  ,           "origin_01"       ,       "birth" , "spei12", "spei3",
                  "year_t1"         ,       "year_t1_index"       ,   "surv_t1" ,
                  "size_t1"         ,       "logsize_t1"       ,
                  "year_t",
                  "year_t_index"     ,      "size_t"           ,      "logsize_t"  ,
                  "FLW_COUNT_T1"      ,      "FLW_STAT_T1"),
       value.name = "spike_count_t1") %>% 
  rename(spikelet_id = variable) %>% 
  filter(!is.na(spike_count_t1), spike_count_t1 > 0) %>% 
  mutate(spike_count_t1 = as.integer(spike_count_t1)) %>% 
  mutate(origin_index = origin_01+1)

# 
# # ggplot(LTREB_data_forspike)+
#   geom_histogram(aes(x=spike_count_t))+
#   facet_grid(year_t~species)
## I don't think there are enough data to fit year variances
## so I am just going to fit fixed effects of size and endo

# rm(LTREB_full)


# Create data lists to be used for the Stan model with 12month and 3month SPEI values
surv_data_list <- list(y = LTREB_data_forsurv$surv_t1,
                       logsize = LTREB_data_forsurv$logsize_t,
                       origin_01 = LTREB_data_forsurv$origin_01,
                       origin_index= LTREB_data_forsurv$origin_index,
                       endo_01 = as.integer(LTREB_data_forsurv$endo_01),
                       endo_index = as.integer(LTREB_data_forsurv$endo_index),
                       spp = as.integer(LTREB_data_forsurv$species_index),
                       year_t = as.integer(LTREB_data_forsurv$year_t_index),
                       plot = as.integer(as.factor(LTREB_data_forsurv$plot_index)),
                       N = nrow(LTREB_data_forsurv),
                       nSpp = length(unique(LTREB_data_forsurv$species_index)),
                       nYear = max(unique(LTREB_data_forsurv$year_t_index)),
                       nPlot = length(unique(LTREB_data_forsurv$plot_index)),
                       nEndo =   length(unique(LTREB_data_forsurv$endo_01)),
                       nOrigin =   length(unique(LTREB_data_forsurv$origin_01)))
surv_spei12_data_list <- append(surv_data_list, list(spei = as.numeric(LTREB_data_forsurv$spei12),
                                                     spei_nl = as.numeric(LTREB_data_forsurv$spei12^2)))
surv_spei3_data_list <- append(surv_data_list, list(spei = as.numeric(LTREB_data_forsurv$spei3),
                                                     spei_nl = as.numeric(LTREB_data_forsurv$spei3^2)))

str(surv_data_list);str(surv_spei12_data_list);str(surv_spei3_data_list)


seed_surv_data_list <- list(y = LTREB_surv_seedling$surv_t1,
                            logsize_t = LTREB_surv_seedling$logsize_t,
                            origin_01 = LTREB_surv_seedling$origin_01,
                            endo_01 = as.integer(LTREB_surv_seedling$endo_01),
                            endo_index = as.integer(LTREB_surv_seedling$endo_index),
                            spp = as.integer(LTREB_surv_seedling$species_index),
                            year_t = as.integer(LTREB_surv_seedling$year_t_index),
                            plot = as.integer(LTREB_surv_seedling$plot_index),
                            N = nrow(LTREB_surv_seedling),
                            nSpp = length(unique(LTREB_surv_seedling$species_index)),
                            nYear = max(unique(LTREB_surv_seedling$year_t_index)),
                            nPlot = length(unique(LTREB_data_forsurv$plot_index)),
                            nEndo =   length(unique(LTREB_surv_seedling$endo_01)))
seed_surv_spei12_data_list <- append(seed_surv_data_list, list(spei = as.numeric(LTREB_surv_seedling$spei12),
                                                               spei_nl = as.numeric(LTREB_surv_seedling$spei12^2)))
seed_surv_spei3_data_list <- append(seed_surv_data_list, list(spei = as.numeric(LTREB_surv_seedling$spei3),
                                                               spei_nl = as.numeric(LTREB_surv_seedling$spei3^2)))

str(seed_surv_data_list);str(seed_surv_spei12_data_list);str(seed_surv_spei3_data_list)


flw_data_list <- list(y = LTREB_data_forflw$FLW_STAT_T1,
                      logsize = LTREB_data_forflw$logsize_t1,
                      origin_01 = LTREB_data_forflw$origin_01,
                      origin_index= LTREB_data_forflw$origin_index,
                      endo_01 = as.integer(LTREB_data_forflw$endo_01),
                      endo_index = as.integer(LTREB_data_forflw$endo_index),
                      spp = as.integer(LTREB_data_forflw$species_index),
                      year_t = as.integer(LTREB_data_forflw$year_t_index),
                      plot = as.integer(LTREB_data_forflw$plot_index),
                      N = nrow(LTREB_data_forflw),
                      nSpp = length(unique(LTREB_data_forflw$species_index)),
                      nYear = max(unique(LTREB_data_forflw$year_t_index)),
                      nPlot = length(unique(LTREB_data_forflw$plot_index)),
                      nEndo =   length(unique(LTREB_data_forflw$endo_01)),
                      nOrigin =   length(unique(LTREB_data_forflw$origin_01)))
flw_spei12_data_list <- append(flw_data_list, list(spei = as.numeric(LTREB_data_forflw$spei12),
                                                   spei_nl = as.numeric(LTREB_data_forflw$spei12^2)))
flw_spei3_data_list <- append(flw_data_list, list(spei = as.numeric(LTREB_data_forflw$spei3),
                                                   spei_nl = as.numeric(LTREB_data_forflw$spei3^2)))

str(flw_data_list);str(flw_spei12_data_list);str(flw_spei3_data_list)

# 
grow_data_list <- list(y = as.integer(LTREB_data_forgrow$size_t1),
                       logsize = LTREB_data_forgrow$logsize_t,
                       origin_01 = as.integer(LTREB_data_forgrow$origin_01),
                       origin_index= LTREB_data_forgrow$origin_index,
                       endo_01 = as.integer(LTREB_data_forgrow$endo_01),
                       endo_index = as.integer(LTREB_data_forgrow$endo_index),
                       spp = as.integer(LTREB_data_forgrow$species_index),
                       year_t = as.integer(LTREB_data_forgrow$year_t_index),
                       plot = as.integer(LTREB_data_forgrow$plot_index),
                       N = nrow(LTREB_data_forgrow),
                       nSpp = length(unique(LTREB_data_forgrow$species_index)),
                       nYear = max(unique(LTREB_data_forgrow$year_t_index)),
                       nPlot = length(unique(LTREB_data_forgrow$plot_index)),
                       nEndo =   length(unique(LTREB_data_forgrow$endo_01)),
                       nOrigin =   length(unique(LTREB_data_forgrow$origin_01)))
grow_spei12_data_list <- append(grow_data_list, list(spei = as.numeric(LTREB_data_forgrow$spei12),
                                                     spei_nl = as.numeric(LTREB_data_forgrow$spei12^2)))
grow_spei3_data_list <- append(grow_data_list, list(spei = as.numeric(LTREB_data_forgrow$spei3),
                                                     spei_nl = as.numeric(LTREB_data_forgrow$spei3^2)))
str(grow_data_list);str(grow_spei12_data_list);str(grow_spei3_data_list)

#
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
seed_grow_spei12_data_list <- append(seed_grow_data_list, list(spei = as.numeric(LTREB_grow_seedling$spei12),
                                                               spei_nl = as.numeric(LTREB_grow_seedling$spei12^2)))
seed_grow_spei3_data_list <- append(seed_grow_data_list, list(spei = as.numeric(LTREB_grow_seedling$spei3),
                                                               spei_nl = as.numeric(LTREB_grow_seedling$spei3^2)))
str(seed_grow_data_list);str(seed_grow_spei12_data_list);str(seed_grow_spei3_data_list)

fert_data_list <- list(y = as.integer(LTREB_data_forfert$FLW_COUNT_T1),
                       logsize = LTREB_data_forfert$logsize_t1,
                       origin_01 = LTREB_data_forfert$origin_01,
                       origin_index= LTREB_data_forfert$origin_index,
                       endo_01 = as.integer(LTREB_data_forfert$endo_01),
                       endo_index = as.integer(LTREB_data_forfert$endo_index),
                       spp = as.integer(LTREB_data_forfert$species_index),
                       year_t = as.integer(LTREB_data_forfert$year_t_index),
                       plot = as.integer(LTREB_data_forfert$plot_index),
                       N = nrow(LTREB_data_forfert),
                       nSpp = length(unique(LTREB_data_forfert$species_index)),
                       nYear = max(unique(LTREB_data_forfert$year_t_index)),
                       nPlot = max(unique(LTREB_data_forfert$plot_index)),
                       nEndo =   length(unique(LTREB_data_forfert$endo_01)),
                       nOrigin =   length(unique(LTREB_data_forfert$origin_01)))
fert_spei12_data_list <- append(fert_data_list, list(spei = as.numeric(LTREB_data_forfert$spei12),
                                                     spei_nl = as.numeric(LTREB_data_forfert$spei12^2)))
fert_spei3_data_list <- append(fert_data_list, list(spei = as.numeric(LTREB_data_forfert$spei3),
                                                     spei_nl = as.numeric(LTREB_data_forfert$spei3^2)))
str(fert_data_list);str(fert_spei12_data_list);str(fert_spei3_data_list);

spike_data_list <- list(nYear = max(unique(LTREB_data_forspike$year_t_index)),
                        nPlot = max(unique(LTREB_data_forspike$plot_index)),
                        nSpp = length(unique(LTREB_data_forspike$species)),
                        nEndo=length(unique(LTREB_data_forspike$endo_01)),
                        nOrigin =   length(unique(LTREB_data_forspike$origin_01)),
                        N = nrow(LTREB_data_forspike),
                        year_t = as.integer(LTREB_data_forspike$year_t_index),
                        plot = as.integer(LTREB_data_forspike$plot_index),
                        spp = as.integer(LTREB_data_forspike$species_index),
                        y = LTREB_data_forspike$spike_count_t1,
                        logsize = LTREB_data_forspike$logsize_t1,
                        endo_01 = LTREB_data_forspike$endo_01,
                        origin_01 = LTREB_data_forspike$origin_01,
                        origin_index= LTREB_data_forspike$origin_index)

spike_spei12_data_list <- append(spike_data_list, list(spei = as.numeric(LTREB_data_forspike$spei12),
                                                       spei_nl = as.numeric(LTREB_data_forspike$spei12^2)))
spike_spei3_data_list <- append(spike_data_list, list(spei = as.numeric(LTREB_data_forspike$spei3),
                                                       spei_nl = as.numeric(LTREB_data_forspike$spei3^2)))
str(spike_data_list);str(spike_spei12_data_list);str(spike_spei3_data_list)
#########################################################################################################
# Stan model runs------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
mcmc_pars <- list(
  warmup = 2500, 
  iter = 5000, 
  thin = 1, 
  chains = 3
)
# Fitting the model with 12 month spei, no quadratic terms
# The model with squared climate terrms gives max_treedepth warnings but otherwise converges okay but without non-linear terms does not have warning.
sm_surv_spei12 <- stan(file = "Analyses/climate_endo_spp_surv_flw.stan", data = surv_spei12_data_list,
                iter = mcmc_pars$iter,
                warmup = mcmc_pars$warmup,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin)
saveRDS(sm_surv_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/climate_spei12_endo_spp_surv_woseedling_linear.rds")
# Fitting the model with 3 month spei, no quadratic terms as well to explore how growing season climate affects vital rates
sm_surv_spei3 <- stan(file = "Analyses/climate_endo_spp_surv_flw.stan", data = surv_spei3_data_list,
                iter = mcmc_pars$iter,
                warmup = mcmc_pars$warmup,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin)
saveRDS(sm_surv_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/climate_spei3_endo_spp_surv_woseedling_linear.rds")



# Quad X origin models
sm_surv_spei12 <- stan(file = "Analyses/climate_endo_spp_surv_flw_quadXorigin.stan", data = surv_spei12_data_list,
                       iter = mcmc_pars$iter,
                       warmup = mcmc_pars$warmup,
                       chains = mcmc_pars$chains, 
                       thin = mcmc_pars$thin)
saveRDS(sm_surv_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/climate_spei12_endo_spp_surv_woseedling_linear_quadXorigin.rds")
# Fitting the model with 3 month spei, no quadratic terms as well to explore how growing season climate affects vital rates
sm_surv_spei3 <- stan(file = "Analyses/climate_endo_spp_surv_flw_quadXorigin.stan", data = surv_spei3_data_list,
                      iter = mcmc_pars$iter,
                      warmup = mcmc_pars$warmup,
                      chains = mcmc_pars$chains, 
                      thin = mcmc_pars$thin)
saveRDS(sm_surv_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/climate_spei3_endo_spp_surv_woseedling_linear_quadXorigin.rds")



#Seedlling survival
# fit with low ESS, so Ran for two times iter 
sm_seed_surv_spei12 <- stan(file = "Analyses/climate_seedling_surv.stan", data = seed_surv_spei12_data_list,
                     iter = mcmc_pars$iter*2,
                     warmup = mcmc_pars$warmup*2,
                     chains = mcmc_pars$chains, 
                     thin = mcmc_pars$thin)
saveRDS(sm_seed_surv_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_spei12.rds")

# spei3 gave maximum reedepth warnings but otherwise is okay
sm_seed_surv_spei3 <- stan(file = "Analyses/climate_seedling_surv.stan", data = seed_surv_spei3_data_list,
                     iter = mcmc_pars$iter,
                     warmup = mcmc_pars$warmup,
                     chains = mcmc_pars$chains, 
                     thin = mcmc_pars$thin)
saveRDS(sm_seed_surv_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_spei3.rds")


# This is looking at the difference between linear and polynomial spei12
sm_seed_surv_linear <- readRDS(file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_withoutquadraticterm.rds")
sm_seed_surv_nonlinear <- readRDS(file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_withtightpriors.rds")

quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}
surv_sdlg_parlinear <- rstan::extract(sm_seed_surv_linear, pars =quote_bare(beta0,betaendo,betaspei_endo,
                                                                tau_year, tau_plot))
surv_sdlg_parnonlinear <- rstan::extract(sm_seed_surv_nonlinear, pars =quote_bare(beta0,betaendo,betaspei_endo,betaspei_nl_endo,
                                                                     tau_year, tau_plot))
linear_s_seed_eminus <- matrix(NA, 7,100)
linear_s_seed_eplus <- matrix(NA, 7,100)

nonlinear_s_seed_eminus <- matrix(NA, 7,100)
nonlinear_s_seed_eplus <- matrix(NA, 7,100)

for(s in 1:7){
linear_s_seed_eminus[s,] <- invlogit(mean(surv_sdlg_parlinear$beta0[,s])+mean(surv_sdlg_parlinear$betaendo[,s])*0+mean(surv_sdlg_parlinear$betaspei_endo[,s,1])*seq(-1,2, length.out = 100))#+mean(surv_sdlg_parlinear$betaspei_nl_endo[,s,1]*(seq(-1,2, length.out = 100)^2)))
linear_s_seed_eplus[s,] <- invlogit(mean(surv_sdlg_parlinear$beta0[,s])+mean(surv_sdlg_parlinear$betaendo[,s])*1+mean(surv_sdlg_parlinear$betaspei_endo[,s,2])*seq(-1,2, length.out = 100))#+mean(surv_sdlg_parlinear$betaspei_nl_endo[,s,2]*(seq(-1,2, length.out = 100)^2)))
nonlinear_s_seed_eminus[s,] <- invlogit(mean(surv_sdlg_parnonlinear$beta0[,s])+mean(surv_sdlg_parnonlinear$betaendo[,s])*0+mean(surv_sdlg_parnonlinear$betaspei_endo[,s,1])*seq(-1,2, length.out = 100)+mean(surv_sdlg_parnonlinear$betaspei_nl_endo[,s,1])*(seq(-1,2, length.out = 100)^2))
nonlinear_s_seed_eplus[s,] <- invlogit(mean(surv_sdlg_parnonlinear$beta0[,s])+mean(surv_sdlg_parnonlinear$betaendo[,s])*1+mean(surv_sdlg_parnonlinear$betaspei_endo[,s,2])*seq(-1,2, length.out = 100)+mean(surv_sdlg_parnonlinear$betaspei_nl_endo[,s,2])*(seq(-1,2, length.out = 100)^2))
}


# plotting the predicted linear and nonlinear fits vs the seedling survival data
pred <- as_tibble(t(linear_s_seed_eminus)) %>% 
  rename("lin_Eminus_"= contains("V")) %>% 
  cbind(spei = seq(-1,2,length.out=100)) %>% 
  pivot_longer(cols = contains("lin")) %>% 
  separate(name, c("effect", "endo", "species"))

pred2 <- as_tibble(t(linear_s_seed_eplus)) %>% 
  rename("lin_Eplus_"= contains("V")) %>% 
  cbind(spei = seq(-1,2,length.out=100)) %>% 
  pivot_longer(cols = contains("lin")) %>% 
  separate(name, c("effect", "endo", "species"))

pred3 <- as_tibble(t(nonlinear_s_seed_eminus)) %>% 
  rename("nonlin_Eminus_"= contains("V")) %>% 
  cbind(spei = seq(-1,2,length.out=100)) %>% 
  pivot_longer(cols = contains("nonlin")) %>% 
  separate(name, c("effect", "endo", "species"))

pred4 <- as_tibble(t(nonlinear_s_seed_eplus)) %>% 
  rename("nonlin_Eplus_"= contains("V")) %>% 
  cbind(spei = seq(-1,2,length.out=100)) %>% 
  pivot_longer(cols = contains("nonlin")) %>% 
  separate(name, c("effect", "endo", "species"))

pred0 <- pred %>% 
  full_join(pred2) %>% 
  full_join(pred3) %>% 
  full_join(pred4) %>% 
  mutate(species = case_when(species == 1 ~ "AGPE",
                             species == 2 ~ "ELRI",
                             species == 3 ~ "ELVI",
                             species == 4 ~ "FESU",
                             species == 5 ~ "LOAR",
                             species == 6 ~ "POAL",
                             species == 7 ~ "POSY")) %>% 
  mutate(endo = case_when(endo == "Eminus" ~ "0",
                          endo == "Eplus" ~ "1"))

seedlingsurv_climate_plot <- ggplot(data = LTREB_seedlingsurv_means)+
  geom_point(aes(x = spei12, y = mean_surv, color = species, shape = endo, size = count))+ 
  geom_line(data = pred0, aes(x = spei, y = value, color = species, lty = endo))+
  # geom_line(data = filter(pred0,effect == "lin"), aes(x = spei, y = value, color = species, lty = endo))+
  # facet_wrap(~species)+
  facet_wrap(~species+effect) +
  scale_color_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84")) +
  scale_shape_manual(values = c(16,1))+
  ylab("Mean Seedling Survival") + xlab("SPEI (12 month)") +
  theme_classic()
seedlingsurv_climate_plot
ggsave(seedlingsurv_climate_plot, filename = "seedlingsurv_climate_plot.png", width = 8, height = 5)





# running the flowering model with just linear SPEI term fits without errors
# fitting the fowering model with the polynomial SPEI term gives max_treedepth errors
sm_flw_spei12 <- stan(file = "Analyses/climate_endo_spp_surv_flw.stan", data = flw_spei12_data_list,
                iter = mcmc_pars$iter,
                warmup = mcmc_pars$warmup,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin)
saveRDS(sm_flw_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw_spei12.rds")
# sm_flw_spei3 had max treedepth warning
sm_flw_spei3 <- stan(file = "Analyses/climate_endo_spp_surv_flw.stan", data = flw_spei3_data_list,
               iter = mcmc_pars$iter,
               warmup = mcmc_pars$warmup,
               chains = mcmc_pars$chains, 
               thin = mcmc_pars$thin)
saveRDS(sm_flw_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw_spei3.rds")


# This is looking at the difference between linear and polynomial spei12
sm_flw_linear <- readRDS(file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw_spei12.rds")

quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}
flw_parlinear <- rstan::extract(sm_flw_linear, pars =quote_bare(beta0,betaendo,betaspei_endo,
                                                                            tau_year, tau_plot))

linear_f_eminus <- matrix(NA, 7,100)
linear_f_eplus <- matrix(NA, 7,100)


for(s in 1:7){
  linear_f_eminus[s,] <- invlogit(mean(flw_parlinear$beta0[,s])+mean(flw_parlinear$betaendo[,s])*0+mean(flw_parlinear$betaspei_endo[,s,1])*seq(-1,2, length.out = 100))#+mean(flw_parlinear$betaspei_nl_endo[,s,1]*(seq(-1,2, length.out = 100)^2)))
  linear_f_eplus[s,] <- invlogit(mean(flw_parlinear$beta0[,s])+mean(flw_parlinear$betaendo[,s])*1+mean(flw_parlinear$betaspei_endo[,s,2])*seq(-1,2, length.out = 100))#+mean(flw_parlinear$betaspei_nl_endo[,s,2]*(seq(-1,2, length.out = 100)^2)))
}


# plotting the predicted linear and nonlinear fits vs the seedling survival data
pred <- as_tibble(t(linear_f_eminus)) %>% 
  rename("lin_Eminus_"= contains("V")) %>% 
  cbind(spei = seq(-1,2,length.out=100)) %>% 
  pivot_longer(cols = contains("lin")) %>% 
  separate(name, c("effect", "endo", "species"))

pred2 <- as_tibble(t(linear_f_eplus)) %>% 
  rename("lin_Eplus_"= contains("V")) %>% 
  cbind(spei = seq(-1,2,length.out=100)) %>% 
  pivot_longer(cols = contains("lin")) %>% 
  separate(name, c("effect", "endo", "species"))

pred0 <- pred %>% 
  full_join(pred2) %>% 
  mutate(species = case_when(species == 1 ~ "AGPE",
                             species == 2 ~ "ELRI",
                             species == 3 ~ "ELVI",
                             species == 4 ~ "FESU",
                             species == 5 ~ "LOAR",
                             species == 6 ~ "POAL",
                             species == 7 ~ "POSY")) %>% 
  mutate(endo = case_when(endo == "Eminus" ~ "0",
                          endo == "Eplus" ~ "1"))

flw_climate_plot <- ggplot(data = LTREB_flw_means)+
  # geom_point(aes(x = spei12, y = mean_flw, color = species, shape = endo, size = count))+
  geom_line(data = pred0, aes(x = spei, y = value, color = species, lty = endo))+
  # geom_line(data = filter(pred0,effect == "lin"), aes(x = spei, y = value, color = species, lty = endo))+
  # facet_wrap(~species)+
  facet_wrap(~species+effect) +
  scale_color_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84")) +
  scale_shape_manual(values = c(16,1))+
  ylab("Mean Flowering Prob.") + xlab("SPEI (12 month)") +
  theme_classic()
flw_climate_plot
ggsave(seedlingsurv_climate_plot, filename = "seedlingsurv_climate_plot.png", width = 8, height = 5)



# fitting the fowering model with the polynomial SPEI term gives max_treedepth errors
sm_flw_spei12 <- stan(file = "Analyses/climate_endo_spp_surv_flw_quadXorigin.stan", data = flw_spei12_data_list,
                      iter = mcmc_pars$iter,
                      warmup = mcmc_pars$warmup,
                      chains = mcmc_pars$chains, 
                      thin = mcmc_pars$thin)
saveRDS(sm_flw_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw_spei12_quadXorigin.rds")
# sm_flw_spei3 had max treedepth warning
sm_flw_spei3 <- stan(file = "Analyses/climate_endo_spp_surv_flw_quadXorigin.stan", data = flw_spei3_data_list,
                     iter = mcmc_pars$iter,
                     warmup = mcmc_pars$warmup,
                     chains = mcmc_pars$chains, 
                     thin = mcmc_pars$thin)
saveRDS(sm_flw_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw_spei3_quadXorigin.rds")



# Running the growth model with the PIG with climate effects
# fits with no errors including polynomial term
# running without the polynomial term , no errors
sm_grow_spei12 <- stan(file = "Analyses/climate_endo_spp_grow_fert_PIG.stan", data = grow_spei12_data_list,
               iter = mcmc_pars$iter,
               warmup = mcmc_pars$warmup,
               chains = mcmc_pars$chains, 
               thin = mcmc_pars$thin)
saveRDS(sm_grow_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_spei12.rds")
# saveRDS(sm_grow, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_nonlinear.rds")

# 3 month spei fits with no warrnings
sm_grow_spei3 <- stan(file = "Analyses/climate_endo_spp_grow_fert_PIG.stan", data = grow_spei3_data_list,
                       iter = mcmc_pars$iter,
                       warmup = mcmc_pars$warmup,
                       chains = mcmc_pars$chains, 
                       thin = mcmc_pars$thin)
saveRDS(sm_grow_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_spei3.rds")


sm_grow_spei12 <- stan(file = "Analyses/climate_endo_spp_grow_fert_PIG_quadXorigin.stan", data = grow_spei12_data_list,
                       iter = mcmc_pars$iter,
                       warmup = mcmc_pars$warmup,
                       chains = mcmc_pars$chains, 
                       thin = mcmc_pars$thin)
saveRDS(sm_grow_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_spei12_quadXorigin.rds")
# saveRDS(sm_grow, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_nonlinear.rds")

# 3 month spei fits with no warrnings
sm_grow_spei3 <- stan(file = "Analyses/climate_endo_spp_grow_fert_PIG_quadXorigin.stan", data = grow_spei3_data_list,
                      iter = mcmc_pars$iter,
                      warmup = mcmc_pars$warmup,
                      chains = mcmc_pars$chains, 
                      thin = mcmc_pars$thin)
saveRDS(sm_grow_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_spei3_quadXorigin.rds")


# sm_grow_spei3 <- readRDS(file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_spei3_quadXorigin.rds")


# fitting linear model with 5000 iterations gives low bulk and tail ESS, rerunning with 10000 (this is what we had to do for non-climate model as well)
# ftting linear model for 10000 iteration fits with no errors or warnings
# fitting non-linear model for 10000 iterations fits with no errors or warnings
sm_seedgrow_spei12 <- stan(file = "Analyses/climate_seedling_grow_PIG.stan", data = seed_grow_spei12_data_list,
                iter = mcmc_pars$iter*2,
                warmup = mcmc_pars$warmup*2,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin)
saveRDS(sm_seedgrow_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/climate_seedling_grow_spei12_10000iterations.rds")
# 3 month spei fits with no errors or warnings
sm_seedgrow_spei3 <- stan(file = "Analyses/climate_seedling_grow_PIG.stan", data = seed_grow_spei3_data_list,
                    iter = mcmc_pars$iter*2,
                    warmup = mcmc_pars$warmup*2,
                    chains = mcmc_pars$chains, 
                    thin = mcmc_pars$thin)
saveRDS(sm_seedgrow_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/climate_seedling_grow_spei3_10000iterations.rds")
 # sm_seedgrow_linear <- readRDS(file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_grow_linear.rds")

# linear fert model fits with no warnings or errors
# nonlinear fert model fits with no warnings or errors
sm_fert_spei12 <- stan(file = "Analyses/climate_endo_spp_grow_fert_PIG.stan", data =fert_spei12_data_list,
                iter = mcmc_pars$iter,
                warmup = mcmc_pars$warmup,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin)
saveRDS(sm_fert_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_spei12.rds")
# saveRDS(sm_fert, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_linear.rds")
# saveRDS(sm_fert, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_nonlinear.rds")

sm_fert_spei3 <- stan(file = "Analyses/climate_endo_spp_grow_fert_PIG.stan", data =fert_spei3_data_list,
                iter = mcmc_pars$iter,
                warmup = mcmc_pars$warmup,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin)
saveRDS(sm_fert_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_spei3.rds")



sm_fert_spei12 <- stan(file = "Analyses/climate_endo_spp_grow_fert_PIG_quadXorigin.stan", data =fert_spei12_data_list,
                       iter = mcmc_pars$iter,
                       warmup = mcmc_pars$warmup,
                       chains = mcmc_pars$chains, 
                       thin = mcmc_pars$thin)
saveRDS(sm_fert_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_spei12_quadXorigin.rds")
# saveRDS(sm_fert, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_linear.rds")
# saveRDS(sm_fert, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_nonlinear.rds")

sm_fert_spei3 <- stan(file = "Analyses/climate_endo_spp_grow_fert_PIG_quadXorigin.stan", data =fert_spei3_data_list,
                      iter = mcmc_pars$iter,
                      warmup = mcmc_pars$warmup,
                      chains = mcmc_pars$chains, 
                      thin = mcmc_pars$thin)
saveRDS(sm_fert_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_spei3_quadXorigin.rds")


#fitting the spike data as negative binomial with just linear spei term
# fitting the spide data with nonlinear term gives tree-depth warnings
sm_spike_nb_spei12 <- stan(file = "Analyses/climate_endo_spp_spike_nb.stan", data = spike_spei12_data_list,
                    iter = mcmc_pars$iter,
                    warmup = mcmc_pars$warmup,
                    chains = mcmc_pars$chains,
                    thin = mcmc_pars$thin)
saveRDS(sm_spike_nb_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_spei12.rds")

# saveRDS(sm_spike_nb, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_linear.rds")
# saveRDS(sm_spike_nb, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_nonlinear.rds")
# 3 month spei fitss max treedepth warnings
sm_spike_nb_spei3 <- stan(file = "Analyses/climate_endo_spp_spike_nb.stan", data = spike_spei3_data_list,
                    iter = mcmc_pars$iter,
                    warmup = mcmc_pars$warmup,
                    chains = mcmc_pars$chains,
                    thin = mcmc_pars$thin)
saveRDS(sm_spike_nb_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_spei3.rds")




sm_spike_nb_spei12 <- stan(file = "Analyses/climate_endo_spp_spike_nb_quadXorigin.stan", data = spike_spei12_data_list,
                           iter = mcmc_pars$iter,
                           warmup = mcmc_pars$warmup,
                           chains = mcmc_pars$chains,
                           thin = mcmc_pars$thin)
saveRDS(sm_spike_nb_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_spei12_quadXorigin.rds")

# saveRDS(sm_spike_nb, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_linear.rds")
# saveRDS(sm_spike_nb, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_nonlinear.rds")
# 3 month spei fitss max treedepth warnings
sm_spike_nb_spei3 <- stan(file = "Analyses/climate_endo_spp_spike_nb_quadXorigin.stan", data = spike_spei3_data_list,
                          iter = mcmc_pars$iter,
                          warmup = mcmc_pars$warmup,
                          chains = mcmc_pars$chains,
                          thin = mcmc_pars$thin)
saveRDS(sm_spike_nb_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_spei3_quadXorigin.rds")



#########################################################################################################
# Linear Model Diagnostics ------------------------------
#########################################################################################################
# Function for looking at binned size_t fits, particularly important for the growth kernel as this determines the transitions through the matrix model
size_moments_ppc <- function(data,y_name,sim, n_bins, title = NA){
  require(tidyverse)
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
    summarize( mean_sim = mean((sim)),
               sd_sim = sd((sim)),
               skew_sim = skewness((sim)),
               kurt_sim = Lkurtosis((sim)),
               bin_mean = mean(logsize_t),
               bin_n = n())
  sim_medians <- sim_moments %>%
    group_by(size_bin, bin_mean) %>%
    summarize(median_mean_sim = median(mean_sim),
              median_sd_sim = median(sd_sim),
              median_skew_sim = median(skew_sim),
              median_kurt_sim = median(kurt_sim))
  meanplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = mean_sim), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_mean_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = mean_t1), shape = 1, color = "firebrick2") +
    theme_classic()
  sdplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = sd_sim), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_sd_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = sd_t1), shape = 1, color = "firebrick2") + theme_classic()
  skewplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = skew_sim), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_skew_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = skew_t1), shape = 1, color = "firebrick2") + theme_classic()
  kurtplot <- ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = kurt_sim), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_kurt_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = kurt_t1), shape = 1, color = "firebrick2") + theme_classic()
  size_ppc_plot <- grid.arrange(meanplot, sdplot,skewplot, kurtplot, top = title)
  return(size_ppc_plot)
}


#### survival ppc ####
# surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_surv_woseedling_linear.rds")
# surv_fit_spei12 <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_spei12_endo_spp_surv_woseedling_linear.rds")
# surv_fit_spei3 <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_spei3_endo_spp_surv_woseedling_linear.rds")


surv_fit_spei12 <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_spei12_endo_spp_surv_woseedling_linear_quadXorigin.rds")
surv_fit_spei3 <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_spei3_endo_spp_surv_woseedling_linear_quadXorigin.rds")



predS_spei12<- rstan::extract(surv_fit_spei12, pars = c("p"))$p # extract the linear predictor
predS_spei3<- rstan::extract(surv_fit_spei3, pars = c("p"))$p # extract the linear predictor
n_post_draws <- 500
post_draws <- sample.int(dim(predS_spei12)[1], n_post_draws) # draw samples from the posterior of the linear predictor
y_s_sim_spei12 <- y_s_sim_spei3 <-  matrix(NA,n_post_draws,length(surv_data_list$y))
for(i in 1:n_post_draws){
  y_s_sim_spei12[i,] <- rbinom(n=length(surv_data_list$y), size=1, prob = invlogit(predS_spei12[post_draws[i],]))
  y_s_sim_spei3[i,] <- rbinom(n=length(surv_data_list$y), size=1, prob = invlogit(predS_spei3[post_draws[i],]))
  
}
saveRDS(y_s_sim_spei12, file = "yrep_climatesurvivalmodel_linear_spei12_quadXorigin.rds")
saveRDS(y_s_sim_spei3, file = "yrep_climatesurvivalmodel_linear_spei3_quadXorigin.rds")

y_s_sim_spei12 <- readRDS(file = "yrep_climatesurvivalmodel_linear_spei12_quadXorigin.rds")
y_s_sim_spei3 <- readRDS(file = "yrep_climatesurvivalmodel_linear_spei3_quadXorigin.rds")

# ppc_dens_overlay(surv_data_list$y, y_s_sim)
surv_densplot_spei12 <- ppc_dens_overlay(surv_spei12_data_list$y, y_s_sim_spei12) + theme_classic() + labs(title = "Adult Survival (12 month SPEI)", x = "Survival status", y = "Density")
surv_densplot_spei12
surv_densplot_spei3 <- ppc_dens_overlay(surv_spei3_data_list$y, y_s_sim_spei3) + theme_classic() + labs(title = "Adult Survival (3 month SPEI)", x = "Survival status", y = "Density")
surv_densplot_spei3

ggsave(surv_densplot_spei12, filename = "climate_surv_densplot_spei12.png", width = 4, height = 4)
ggsave(surv_densplot_spei12, filename = "climate_surv_densplot_spei12.png", width = 4, height = 4)

traceplot(surv_fit_spei12, pars = "sigma0")
traceplot(surv_fit_spei3, pars = "sigma0")


mean_s_spei12_plot <-   ppc_stat(surv_spei12_data_list$y, y_s_sim_spei12, stat = "mean")
sd_s_spei12_plot <- ppc_stat(surv_spei12_data_list$y, y_s_sim_spei12, stat = "sd")
skew_s_spei12_plot <- ppc_stat(surv_spei12_data_list$y, y_s_sim_spei12, stat = "skewness")
kurt_s_spei12_plot <- ppc_stat(surv_spei12_data_list$y, y_s_sim_spei12, stat = "Lkurtosis")
surv_moments_spei12 <- grid.arrange(mean_s_spei12_plot,sd_s_spei12_plot,skew_s_spei12_plot,kurt_s_spei12_plot,  top = "Survival (12 month SPEI)")
# ggsave(surv_moments_spei12, filename = "climate_surv_spei12_momentplot.png", width = 4, height = 4)

mean_s_spei3_plot <-   ppc_stat(surv_spei3_data_list$y, y_s_sim_spei3, stat = "mean")
sd_s_spei3_plot <- ppc_stat(surv_spei3_data_list$y, y_s_sim_spei3, stat = "sd")
skew_s_spei3_plot <- ppc_stat(surv_spei3_data_list$y, y_s_sim_spei3, stat = "skewness")
kurt_s_spei3_plot <- ppc_stat(surv_spei3_data_list$y, y_s_sim_spei3, stat = "Lkurtosis")
surv_moments_spei3 <- grid.arrange(mean_s_spei3_plot,sd_s_spei3_plot,skew_s_spei3_plot,kurt_s_spei3_plot,  top = "Survival (12 month SPEI)")
# ggsave(surv_moments_spei3, filename = "climate_surv_spei3_momentplot.png", width = 4, height = 4)


# now we want to look at how the the model is fitting across sizes

surv_size_spei12_ppc <- size_moments_ppc(data = LTREB_data_forsurv,
                                  y_name = "surv_t1",
                                  sim = y_s_sim_spei12, 
                                  n_bins = 4, 
                                  title = "Survival")
# ggsave(surv_size_spei12_ppc, filename = "climate_surv_size_spei12_ppc.png", width = 4, height = 4)
surv_size_spei3_ppc <- size_moments_ppc(data = LTREB_data_forsurv,
                                         y_name = "surv_t1",
                                         sim = y_s_sim_spei3, 
                                         n_bins = 4, 
                                         title = "Survival")
# ggsave(surv_size_spei3_ppc, filename = "climate_surv_size_spei3_ppc.png", width = 4, height = 4)


#### seedling survival ppc ####
surv_seed_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_withoutquadraticterm.rds")
predseedS <- rstan::extract(surv_seed_fit, pars = c("p"))$p
n_post_draws <- 500
post_draws <- sample.int(dim(predseedS)[1], n_post_draws)
y_seed_s_sim <- matrix(NA,n_post_draws,length(seed_surv_data_list$y))
for(i in 1:n_post_draws){
  y_seed_s_sim[i,] <- rbinom(n=length(seed_surv_data_list$y), size=1, prob = invlogit(predseedS[post_draws[i],]))
}
saveRDS(y_seed_s_sim, file = "yrep_climate_seedlingsurvivalmodel_linear.rds")
y_seed_s_sim <- readRDS(file = "yrep_climate_seedlingsurvivalmodel_linear.rds")
# ppc_dens_overlay(seed_surv_data_list$y, y_s_sim)
seedsurv_densplot <- ppc_dens_overlay(seed_surv_data_list$y, y_seed_s_sim) + theme_classic() + labs(title = "Seedling Survival", x = "Survival status", y = "Density")
# seedsurv_densplot
# ggsave(seedsurv_densplot, filename = "climate_seedsurv_densplot.png", width = 4, height = 4)

mean_s_plot <-   ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "mean")
sd_s_plot <- ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "sd")
skew_s_plot <- ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "skewness")
kurt_s_plot <- ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "Lkurtosis")
climate_seedsurv_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Seedling Survival")
# climate_seedsurv_moments
# ggsave(climate_seedsurv_moments, filename = "climate_seedsurv_momentsplot.png", width = 4, height = 4)


#### flowering ppc ####
flow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw.rds")
predF <- rstan::extract(flow_fit, pars = c("p"))$p
n_post_draws <- 500
post_draws <- sample.int(dim(predF)[1], n_post_draws)
y_f_sim <- matrix(NA,n_post_draws,length(flw_data_list$y))
for(i in 1:n_post_draws){
  y_f_sim[i,] <- rbinom(n=length(flw_data_list$y), size=1, prob = invlogit(predF[post_draws[i],]))
}
saveRDS(y_f_sim, file = "yrep_climate_floweringmodel_linear.rds")
# y_f_sim <- readRDS(file = "yrep_climate_floweringmodel_linear.rds")

# ppc_dens_overlay(flw_data_list$y, y_f_sim)
flw_densplot <- ppc_dens_overlay(flw_data_list$y, y_f_sim) + theme_classic() + labs(title = "Flowering", x = "Flowering status", y = "Density")
flw_densplot
ggsave(flw_densplot, filename = "climate_flw_densplot.png", width = 4, height = 4)


mean_f_plot <-   ppc_stat(flw_data_list$y, y_f_sim, stat = "mean")
sd_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "sd")
skew_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "skewness")
kurt_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "Lkurtosis")
flw_moments <- mean_f_plot+sd_f_plot+skew_f_plot+kurt_f_plot +plot_annotation(title = "Flowering")
flw_moments
# ggsave(flw_moments, filename = "climate_flw_momentsplot.png", width = 4, height = 4)

# now we want to look at how the the model is fitting across sizes
flw_size_ppc <- size_moments_ppc(data = LTREB_data_forflw,
                                 y_name = "FLW_STAT_T1",
                                 sim = y_f_sim, 
                                 n_bins = 2, 
                                 title = "Flowering")
# ggsave(flw_size_ppc, filename = "climate_flw_size_ppc.png", width = 4, height = 4)

#### Growth ppc ####
#growth PIG distribution
grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_linear.rds")
# pairs(grow_fit, pars = c("beta0"))
grow_par <- rstan::extract(grow_fit, pars = c("predG","beta0","betasize","betaendo","betaorigin","tau_year","tau_plot","theta", "sigma"))
predG <- grow_par$predG
theta <- grow_par$theta
n_post_draws <- 500
post_draws <- sample.int(dim(predG)[1], n_post_draws)
y_g_sim   <-  matrix(NA,n_post_draws,length(grow_data_list$y))

# simulate data
for(i in 1:n_post_draws){
  ## sample growth data (zero-truncated PIG)
  for(j in 1:length(grow_data_list$y)){
    # probability without truncation
    prob_v <- dpois( 1:1000, 
                     lambda = (predG[i,j] * theta[i,j]) )
    # probability for trunctation (denominator)
    prob_t <- (1 - dpois(0, lambda = (predG[i,j] * theta[i,j]) ) )
    
    y_g_sim[i,j] <- sample( x=1:1000, 
                            size = 1, replace = T,
                            prob = prob_v / prob_t )
  }
}
saveRDS(y_g_sim, file = "yrep_climate_growthPIGmodel_linear.rds")
y_g_sim <- readRDS(file = "yrep_climate_growthPIGmodel_linear.rds")

# Posterior predictive check
ppc_dens_overlay(grow_data_list$y, y_g_sim)
ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(0,60) + ggtitle("growth w/o seedling PIG")
ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(50,120)

grow_densplot <- ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(0,60) + theme_classic() + labs(title = "Growth w/o seedling PIG", x = "No. of Tillers", y = "Density")
grow_densplot
ggsave(grow_densplot, filename = "climate_grow_densplot.png", width = 4, height = 4)



mean_g_plot <-   ppc_stat(grow_data_list$y, y_g_sim, stat = "mean")
sd_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "sd")
skew_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "skewness")
kurt_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "Lkurtosis")
grow_moments <- mean_g_plot+sd_g_plot+skew_g_plot+kurt_g_plot+ plot_annotation(title = "Growth PIG")
grow_moments
ggsave(grow_moments, filename = "climate_grow_momentsplot.png", width = 4, height = 4)


# now we want to look at how the the growth model is fitting especially across sizes, as this determines the transitions through the matrix model
PIG_growth_size_ppc <- size_moments_ppc(data = LTREB_data_forgrow,
                                        y_name = "size_t1",
                                        sim = y_g_sim, 
                                        n_bins = 6, 
                                        title = "Growth PIG")
ggsave(PIG_growth_size_ppc, filename = "PIG_climate_growth_size_pcc.png", width = 4, height = 4)


#### Seedling Growth ppc ####
#seedling growth PIG distribution
s_grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_seedling_grow_linear_10000iterations.rds")
# pairs(grow_fit, pars = c("beta0"))
s_grow_par <- rstan::extract(s_grow_fit, pars = c("predG","beta0","betaendo","tau_year","tau_plot","theta", "sigma"))
predseedG <- s_grow_par$predG
seed_theta <- s_grow_par$theta
n_post_draws <- 500
post_draws <- sample.int(dim(predseedG)[1], n_post_draws)
y_seed_g_sim   <-  matrix(NA,n_post_draws,length(seed_grow_data_list$y))

# simulate data
for(i in 1:n_post_draws){
  ## sample growth data (zero-truncated PIG)
  for(j in 1:length(seed_grow_data_list$y)){
    # probability without truncation
    prob_v <- dpois( 1:1000, 
                     lambda = (predseedG[i,j] * seed_theta[i,j]) )
    # probability for trunctation (denominator)
    prob_t <- (1 - dpois(0, lambda = (predseedG[i,j] * seed_theta[i,j]) ) )
    
    y_seed_g_sim[i,j] <- sample( x=1:1000, 
                                 size = 1, replace = T,
                                 prob = prob_v / prob_t )
  }
}

saveRDS(y_seed_g_sim, file = "yrep_climate_seedlinggrowthPIGmodel_linear.rds")
y_seed_g_sim <- readRDS(file = "yrep_climate_seedlinggrowthPIGmodel_linear.rds")

# Posterior predictive check
ppc_dens_overlay(seed_grow_data_list$y, y_seed_g_sim)
ppc_dens_overlay(seed_grow_data_list$y, y_seed_g_sim) + xlim(0,20) + ggtitle("seedling growth w/ PIG")
ppc_dens_overlay(seed_grow_data_list$y, y_seed_g_sim) + xlim(50,120)

seedgrow_densplot <- ppc_dens_overlay(seed_grow_data_list$y, y_seed_g_sim) + xlim(0,60) + theme_classic() + labs(title = "Seedling Growth with PIG", x = "No. of Tillers", y = "Density")
seedgrow_densplot
ggsave(seedgrow_densplot, filename = "seed_grow_densplot.png", width = 4, height = 4)



mean_seed_g_plot <-   ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "mean")
sd_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "sd")
skew_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "skewness")
kurt_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "Lkurtosis")
seedgrow_moments <- mean_seed_g_plot+sd_seed_g_plot+skew_seed_g_plot+kurt_seed_g_plot+ plot_annotation(title = "Seedling Growth PIG")
seedgrow_moments
ggsave(seedgrow_moments, filename = "seedgrow_momentsplot.png", width = 4, height = 4)


#### Fertility ppc ####
# Fert fit with Poisson inverse Gaussian

fert_fit_pig <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_linear.rds")
fert_par <- rstan::extract(fert_fit_pig, pars = c("predG","beta0","betasize","betaendo","betaorigin","tau_year","tau_plot","theta", "sigma"))
predG <- fert_par$predG
theta <- fert_par$theta
n_post_draws <- 500
post_draws <- sample.int(dim(predG)[1], n_post_draws)
y_fert_sim   <-  matrix(NA,n_post_draws,length(fert_data_list$y))

# simulate data
for(i in 1:n_post_draws){
  ## sample fertilty data (zero-truncated PIG)
  for(j in 1:length(fert_data_list$y)){
    # probability without truncation
    prob_v <- dpois( 1:1000, 
                     lambda = (predG[i,j] * theta[i,j]) )
    # probability for trunctation (denominator)
    prob_t <- (1 - dpois(0, lambda = (predG[i,j] * theta[i,j]) ) )
    
    y_fert_sim[i,j] <- sample( x=1:1000, 
                               size = 1, replace = T,
                               prob = prob_v / prob_t )
  }
}
saveRDS(y_fert_sim, file = "yrep_climate_fertilityPIGmodel_linear.rds")
y_fert_sim <- readRDS(file = "yrep_climate_fertilityPIGmodel_linear.rds")


# Posterior predictive check
ppc_dens_overlay(fert_data_list$y, y_fert_sim)
ppc_dens_overlay(fert_data_list$y, y_fert_sim) + xlim(50,120)
fert_densplot <- ppc_dens_overlay(fert_data_list$y, y_fert_sim) + xlim(0,40) + ggtitle("Fertility with PIG")
fert_densplot
ggsave(fert_densplot, filename = "fert_densplot_withPIG.png", width = 4, height = 4)

# This doesn't fit all the moments quite as well, but it's likely good enough
mean_fert_plot <-   ppc_stat(fert_data_list$y, y_fert_sim, stat = "mean")
sd_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "sd")
skew_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "skewness")
kurt_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "Lkurtosis")
fert_moments <- mean_fert_plot+sd_fert_plot+skew_fert_plot+kurt_fert_plot + plot_annotation(title = "Fertility with PIG")
fert_moments
ggsave(fert_moments, filename = "fert_momentsplot_withPIG.png", width = 4, height = 4)

# Now we can look at the binned size fit for fertility
fert_size_ppc <- size_moments_ppc(data = LTREB_data_forfert,
                                  y_name = "FLW_COUNT_T1",
                                  sim = y_fert_sim, 
                                  n_bins = 3, 
                                  title = "Inflorescence Count")
ggsave(fert_size_ppc, filename = "fert_size_ppc.png", width = 4, height = 4)


#### spikelet ppc ####
# looking at the negative binomial fit, fits really nicely
spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_linear.rds")
spike_par <- rstan::extract(spike_fit, pars = c("lambda","beta0","betasize","betaendo","betaorigin","tau_year","tau_plot", "phi", "od"))
predSpike <- spike_par$lambda
phiSpike <- spike_par$phi
odSpike <- spike_par$od

n_post_draws <- 500
post_draws <- sample.int(dim(predSpike)[1], n_post_draws)
y_spike_sim <- matrix(NA,n_post_draws,length(spike_data_list$y))
for(i in 1:n_post_draws){
  for(j in 1:length(spike_data_list$y)){
    y_spike_sim[i,j] <- sample(x = 1:max(spike_data_list$y), size = 1, replace = T, prob = dnbinom(1:max(spike_data_list$y), mu = exp(predSpike[post_draws[i],j]), size = odSpike[post_draws[i],j]))
  }
}
saveRDS(y_spike_sim, file = "yrep_climate_spikeletNBmodel_linear.rds")
y_spike_sim <- readRDS(file = "yrep_climate_spikeletNBmodel_linear.rds")

ppc_dens_overlay(spike_data_list$y, y_spike_sim)
ppc_dens_overlay(spike_data_list$y, y_spike_sim) + xlim(0,250) + ggtitle("Spikelet Count")
spike_densplot <- ppc_dens_overlay(spike_data_list$y, y_spike_sim) + xlim(0,250) + ggtitle("Spikelet Count")
spike_densplot
ggsave(spike_densplot, filename = "spike_densplot.png", width = 4, height = 4)

mean_spike_plot <-   ppc_stat(spike_data_list$y, y_spike_sim, stat = "mean")
sd_spike_plot <- ppc_stat(spike_data_list$y, y_spike_sim, stat = "sd")
skew_spike_plot <- ppc_stat(spike_data_list$y, y_spike_sim, stat = "skewness")
kurt_spike_plot <- ppc_stat(spike_data_list$y, y_spike_sim, stat = "Lkurtosis")
spike_moments <- mean_spike_plot+sd_spike_plot+skew_spike_plot+kurt_spike_plot+ plot_annotation(title = "Spikelets per Infl. ZTNB")
spike_moments
ggsave(spike_moments, filename = "spike_momentsplot.png", width = 4, height = 4)

