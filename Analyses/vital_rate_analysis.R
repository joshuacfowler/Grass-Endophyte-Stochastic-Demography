## Title: Grass endophyte population model with a bayesian framework
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
library(patchwork)


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
# Curious about lifespan

max_ages <- LTREB_full %>% 
  mutate(age = year_t1 - birth,
         age_at_death = case_when(surv_t1 == 0 ~ year_t1 - birth)) 

max_survival <- LTREB_full %>% 
  mutate(age = year_t1 - birth,
         age_at_death = case_when(surv_t1 == 0 ~ year_t1 - birth)) %>% 
  group_by(species) %>% 
  summarize(max_age = max(age, na.rm = T),
            mean_age = mean(age, na.rm = T),
            max_survival = max(age_at_death, na.rm = T),
            mean_survival = mean(age_at_death, na.rm = T),
            median_survival = median(age_at_death, na.rm = T))
# max_survival
mean_ageplot <- ggplot(data = max_ages)+
  geom_histogram(aes(x = age)) +
  facet_wrap(~species) + theme_classic()
# mean_ageplot
# ggsave(mean_ageplot, filename = "~/Documents/mean_ageplot.png", width = 4, height = 4)

# Getting some summary numbers
# indiv_plants <- LTREB_full %>% 
#   group_by(id) %>% 
#   summarize(obs = n())
# adul_plant <- LTREB_full %>% 
#   filter(!is.na(size_t))
# nrow(LTREB_full)
# nrow(indiv_plants)

## Clean up the main data frame for NA's, other small data entry errors
LTREB_data_forsurv <- LTREB_full %>% 
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(endo_01)) %>%  # There are a few LOAR that don't have a plot level endo assigned
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

# Looking at the data coverage across years, 
# and if there are differences between Shaun data and the data since 2018
# LTREB_data_forspike %>% 
#   filter(species == "AGPE") %>% 
# ggplot()+
#   geom_histogram(aes(x=spike_count_t1))+
#   facet_grid(year_t1~species)
# 
# LTREB_data_forspike %>% 
#   filter(species == "ELVI") %>% 
#   ggplot()+
#   geom_histogram(aes(x=spike_count_t1))+
#   facet_grid(year_t1~species)
# 
# LTREB_data_forspike %>% 
#   filter(species == "ELRI") %>% 
#   ggplot()+
#   geom_histogram(aes(x=spike_count_t1))+
#   facet_grid(year_t1~species)
# 
# LTREB_data_forspike %>% 
#   filter(species == "FESU") %>% 
#   ggplot()+
#   geom_histogram(aes(x=spike_count_t1))+
#   facet_grid(year_t1~species)
# 
# LTREB_data_forspike %>% 
#   filter(species == "POAL") %>% 
#   ggplot()+
#   geom_histogram(aes(x=spike_count_t1))+
#   facet_grid(year_t1~species)
# 
# LTREB_data_forspike %>% 
#   filter(species == "POSY") %>% 
#   ggplot()+
#   geom_histogram(aes(x=spike_count_t1))+
#   facet_grid(year_t1~species)
# 
# LTREB_data_forspike %>% 
#   filter(species == "LOAR") %>% 
#   ggplot()+
#   geom_histogram(aes(x=spike_count_t1))+
#   facet_grid(year_t1~species)
# 


# 
# pre_2018_data <- filter(LTREB_data_forspike, year_t1<2018 & species == "AGPE")$spike_count_t1
# post_2018_data <- filter(LTREB_data_forspike, year_t1>=2018 & species == "AGPE")$spike_count_t1
# 
# t.test(pre_2018_data, post_2018_data, two.sided = T)
# wilcox.test(pre_2018_data, post_2018_data)
# 
# pre_2018_data <- filter(LTREB_data_forspike, year_t1<2018 & species == "ELVI")$spike_count_t1
# post_2018_data <- filter(LTREB_data_forspike, year_t1>=2018 & species == "ELVI")$spike_count_t1
# 
# t.test(pre_2018_data, post_2018_data, two.sided = T)
# wilcox.test(pre_2018_data, post_2018_data)
# 
# pre_2018_data <- filter(LTREB_data_forspike, year_t1<2018 & species == "ELRI")$spike_count_t1
# post_2018_data <- filter(LTREB_data_forspike, year_t1>=2018 & species == "ELRI")$spike_count_t1
# 
# t.test(pre_2018_data, post_2018_data, two.sided = T)
# wilcox.test(pre_2018_data, post_2018_data)
# 
# pre_2018_data <- filter(LTREB_data_forspike, year_t1<2018 & species == "POAL")$spike_count_t1
# post_2018_data <- filter(LTREB_data_forspike, year_t1>=2018 & species == "POAL")$spike_count_t1
# 
# t.test(pre_2018_data, post_2018_data, two.sided = T)
# wilcox.test(pre_2018_data, post_2018_data)
# 
# pre_2018_data <- filter(LTREB_data_forspike, year_t1<2018 & species == "FESU")$spike_count_t1
# post_2018_data <- filter(LTREB_data_forspike, year_t1>=2018 & species == "FESU")$spike_count_t1
# 
# t.test(pre_2018_data, post_2018_data, two.sided = T)
# wilcox.test(pre_2018_data, post_2018_data)
# 
# #Note to self, look at number of rows that have flower data but not spikelet data
# dim(LTREB_data_forfert) #There are 5072 rows that have flower status = 1 and have inflorescence counts
# 
# flow_but_not_spike <- LTREB_full %>%
#   filter(FLW_STAT_T1 == 1 & is.na(SPIKE_A_T1) & is.na(SPIKE_AGPE_MEAN_T1)) # There are 330 rows that have that have flow status without spikelet counts.
# 
# #Some of these are one off plants where the infloresce wasn't expanded (especially some of the AGPE in 2021), but others were potentially just forgotten. There are a few years and species that have signifiicant chunks missing, and so I'm going to look at the raw data for those years.
# table(data = flow_but_not_spike$species, flow_but_not_spike$year_t1) #Look at AGPE 2014, FESU 2012, LOAR 2011, POAL 2010, POSY 2010-2011

#And revisit the data merging for each species



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

#########################################################################################################
# Stan model runs------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
mcmc_pars <- list(
  iter = 5000, 
  warmup = 2500, 
  thin = 1, 
  chains = 3
)


sm_surv <- stan(file = "Analyses/endo_spp_surv_flw.stan", data = surv_data_list,
                iter = mcmc_pars$iter,
                warmup = mcmc_pars$warmup,
                chains = mcmc_pars$chains, 
                thin = mcmc_pars$thin)
saveRDS(sm_surv, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_surv_woseedling.rds")

sm_seed_surv <- stan(file = "Analyses/seedling_surv.stan", data = seed_surv_data_list,
                     iter = mcmc_pars$iter,
                     warmup = mcmc_pars$warmup,
                     chains = mcmc_pars$chains, 
                     thin = mcmc_pars$thin)
saveRDS(sm_seed_surv, file = "~/Dropbox/EndodemogData/Model_Runs/endo_seedling_surv.rds")

sm_flw <- stan(file = "Analyses/endo_spp_surv_flw.stan", data = flw_data_list,
               iter = mcmc_pars$iter,
               warmup = mcmc_pars$warmup,
               chains = mcmc_pars$chains, 
               thin = mcmc_pars$thin)
saveRDS(sm_flw, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_flw.rds")

# Negative binomial growth model, we are using the pig for better fit
# sm_grow <- stan(file = "Analyses/endo_spp_grow_fert.stan", data = grow_data_list,
#                 iter = mcmc_pars$iter,
#                 warmup = mcmc_pars$warmup,
#                 chains = mcmc_pars$chains, 
#                 thin = mcmc_pars$thin)
# saveRDS(sm_grow, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow.rds")

sm_grow_pig <- stan(file = "Analyses/endo_spp_grow_fert_PIG.stan", data = grow_data_list,
                    iter = mcmc_pars$iter,
                    warmup = mcmc_pars$warmup,
                    chains = mcmc_pars$chains, 
                    thin = mcmc_pars$thin)
saveRDS(sm_grow_pig, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow_PIG.rds")

# Negative binomial seedling growth model. We are going with the PIG cause it fits better and this model has some sampling issues
# sm_seed_grow <- stan(file = "Analyses/seedling_grow.stan", data = seed_grow_data_list,
#                 iter = mcmc_pars$iter,
#                 warmup = mcmc_pars$warmup,
#                 chains = mcmc_pars$chains, 
#                 thin = mcmc_pars$thin,
#                 control = list(max_treedepth = 15))
# saveRDS(sm_seed_grow, file = "~/Dropbox/EndodemogData/Model_Runs/endo_seedling_grow.rds")

#This fits well, but had to run for a longer time (10000 iterations) to get good effective sample size from Stan.
sm_seed_grow_pig <- stan(file = "Analyses/seedling_grow_PIG.stan", data = seed_grow_data_list,
                    iter = mcmc_pars$iter*2,
                    warmup = mcmc_pars$warmup*2,
                    chains = mcmc_pars$chains, 
                    thin = mcmc_pars$thin)
saveRDS(sm_seed_grow_pig, file = "~/Dropbox/EndodemogData/Model_Runs/endo_seedling_grow_PIG_10000iterations.rds")


# 
# sm_fert <- stan(file = "Analyses/endo_spp_grow_fert.stan", data = fert_data_list,
#                iter = mcmc_pars$iter,
#                warmup = mcmc_pars$warmup,
#                chains = mcmc_pars$chains, 
#                thin = mcmc_pars$thin)
# # saveRDS(sm_fert, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert.rds")
# 
# sm_fert_nc <- stan(file = "Analyses/endo_spp_grow_fert_nc.stan", data = fert_data_list,
#                 iter = mcmc_pars$iter,
#                 warmup = mcmc_pars$warmup,
#                 chains = mcmc_pars$chains, 
#                 thin = mcmc_pars$thin)
# # saveRDS(sm_fert_nc, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_nc.rds")
# 
# 
# sm_fert_noplot <- stan(file = "Analyses/endo_spp_grow_fert_noplot.stan", data = fert_data_list,
#                 iter = mcmc_pars$iter,
#                 warmup = mcmc_pars$warmup,
#                 chains = mcmc_pars$chains, 
#                 thin = mcmc_pars$thin,
#                 control = list(adapt_delta = .99)) # run without this, there is a small number of divergent transitions
# # saveRDS(sm_fert_noplot, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_noplot.rds")

sm_fert_pig <- stan(file = "Analyses/endo_spp_grow_fert_PIG.stan", data = fert_data_list,
                    iter = mcmc_pars$iter,
                    warmup = mcmc_pars$warmup,
                    chains = mcmc_pars$chains, 
                    thin = mcmc_pars$thin)
saveRDS(sm_fert_pig, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_pig.rds")


# Fitting spikelet data as a poisson, this converges without errors mostly sampling with the increased treedepth
# sm_spike_pois <- stan(file = "Analyses/endo_spp_spike_poisson.stan", data = spike_data_list,
#                  iter = mcmc_pars$iter,
#                  warmup = mcmc_pars$warmup,
#                  chains = mcmc_pars$chains, 
#                  thin = mcmc_pars$thin,
#                  control = list(max_treedepth = 15))
# # saveRDS(sm_spike_pois, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike_year_plot_poisson.rds")

#fitting the spike data as negative binomial because the poisson fits well for mean, but not sd
sm_spike_nb <- stan(file = "Analyses/endo_spp_spike_nb.stan", data = spike_data_list,
                      iter = mcmc_pars$iter,
                      warmup = mcmc_pars$warmup,
                      chains = mcmc_pars$chains, 
                      thin = mcmc_pars$thin)
saveRDS(sm_spike_nb, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike_year_plot_nb.rds")

#########################################################################################################
# Model Diagnostics ------------------------------
#########################################################################################################
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
  size_ppc_plot <- meanplot+ sdplot+skewplot+ kurtplot+plot_annotation(title = title)
  return(size_ppc_plot)
}
endophyte_color_scheme <- c("#fdedd3","#f3c8a8", "#5a727b", "#4986c7", "#181914",  "#163381")
color_scheme_set(endophyte_color_scheme)
# color_scheme_view()

#### survival ppc ####
surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_surv_woseedling.rds")
predS <- rstan::extract(surv_fit, pars = c("p"))$p
n_post_draws <- 500
post_draws <- sample.int(dim(predS)[1], n_post_draws)
y_s_sim <- matrix(NA,n_post_draws,length(surv_data_list$y))
for(i in 1:n_post_draws){
  y_s_sim[i,] <- rbinom(n=length(surv_data_list$y), size=1, prob = invlogit(predS[post_draws[i],]))
}
saveRDS(y_s_sim, file = "yrep_survivalmodel.rds")
y_s_sim <- readRDS(file = "yrep_survivalmodel.rds")
# ppc_dens_overlay(surv_data_list$y, y_s_sim)
surv_densplot <- ppc_dens_overlay(surv_data_list$y, y_s_sim) + theme_classic() + labs(title = "Adult Survival", x = "Survival status", y = "Density")
surv_densplot
# ggsave(surv_densplot, filename = "surv_densplot.png", width = 4, height = 4)

mean_s_plot <-   ppc_stat(surv_data_list$y, y_s_sim, stat = "mean")
sd_s_plot <- ppc_stat(surv_data_list$y, y_s_sim, stat = "sd")
skew_s_plot <- ppc_stat(surv_data_list$y, y_s_sim, stat = "skewness")
kurt_s_plot <- ppc_stat(surv_data_list$y, y_s_sim, stat = "Lkurtosis")
surv_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Survival")
surv_moments
# ggsave(surv_moments, filename = "surv_momentplot.png", width = 4, height = 4)
# now we want to look at how the the model is fitting across sizes

surv_size_ppc <- size_moments_ppc(data = LTREB_data_forsurv,
                                         y_name = "surv_t1",
                                         sim = y_s_sim, 
                                         n_bins = 4, 
                                         title = "Survival")
surv_size_ppc
# ggsave(surv_size_ppc, filename = "surv_size_ppc.png", width = 4, height = 4)

#### seedling survival ppc ####
surv_seed_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_surv.rds")
predseedS <- rstan::extract(surv_seed_fit, pars = c("p"))$p
n_post_draws <- 500
post_draws <- sample.int(dim(predseedS)[1], n_post_draws)
y_seed_s_sim <- matrix(NA,n_post_draws,length(seed_surv_data_list$y))
for(i in 1:n_post_draws){
  y_seed_s_sim[i,] <- rbinom(n=length(seed_surv_data_list$y), size=1, prob = invlogit(predseedS[post_draws[i],]))
}
saveRDS(y_seed_s_sim, file = "yrep_seedlingsurvivalmodel.rds")
y_seed_s_sim <- readRDS(file = "yrep_seedlingsurvivalmodel.rds")
# ppc_dens_overlay(seed_surv_data_list$y, y_s_sim)
seedsurv_densplot <- ppc_dens_overlay(seed_surv_data_list$y, y_seed_s_sim) + theme_classic() + labs(title = "Seedling Survival", x = "Survival status", y = "Density")
seedsurv_densplot
# ggsave(seedsurv_densplot, filename = "seedsurv_densplot.png", width = 4, height = 4)

mean_s_plot <-   ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "mean")
sd_s_plot <- ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "sd")
skew_s_plot <- ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "skewness")
kurt_s_plot <- ppc_stat(seed_surv_data_list$y, y_seed_s_sim, stat = "Lkurtosis")
seedsurv_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Seedling Survival")
seedsurv_moments
# ggsave(seedsurv_moments, filename = "seedsurv_momentsplot.png", width = 4, height = 4)

#### flowering ppc ####
flow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_flw.rds")
predF <- rstan::extract(flow_fit, pars = c("p"))$p
n_post_draws <- 500
post_draws <- sample.int(dim(predF)[1], n_post_draws)
y_f_sim <- matrix(NA,n_post_draws,length(flw_data_list$y))
for(i in 1:n_post_draws){
  y_f_sim[i,] <- rbinom(n=length(flw_data_list$y), size=1, prob = invlogit(predF[post_draws[i],]))
}
saveRDS(y_f_sim, file = "yrep_floweringmodel.rds")
y_f_sim <- readRDS(file = "yrep_floweringmodel.rds")

# ppc_dens_overlay(flw_data_list$y, y_f_sim)
flw_densplot <- ppc_dens_overlay(flw_data_list$y, y_f_sim) + theme_classic() + labs(title = "Flowering", x = "Flowering status", y = "Density")
flw_densplot
ggsave(flw_densplot, filename = "flw_densplot.png", width = 4, height = 4)


mean_f_plot <-   ppc_stat(flw_data_list$y, y_f_sim, stat = "mean")
sd_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "sd")
skew_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "skewness")
kurt_f_plot <- ppc_stat(flw_data_list$y, y_f_sim, stat = "Lkurtosis")
flw_moments <- mean_f_plot+sd_f_plot+skew_f_plot+kurt_f_plot +plot_annotation(title = "Flowering")
flw_moments
ggsave(flw_moments, filename = "flw_momentsplot.png", width = 4, height = 4)

# now we want to look at how the the model is fitting across sizes
flw_size_ppc <- size_moments_ppc(data = LTREB_data_forflw,
                                  y_name = "FLW_STAT_T1",
                                  sim = y_f_sim, 
                                  n_bins = 2, 
                                  title = "Flowering")
ggsave(flw_size_ppc, filename = "flw_size_ppc.png", width = 4, height = 4)


#### growth ppc ####
#This is the ZTNB model, it doesn't fit the variance very well
grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow.rds")
# pairs(grow_fit, pars = c("beta0"))
grow_par <- rstan::extract(grow_fit, pars = c("lambda","beta0","betasize","betaendo","betaorigin","tau_year","tau_plot", "phi", "od"))
predG <- grow_par$lambda
phiG <- grow_par$phi
odG <- grow_par$od

n_post_draws <- 100
post_draws <- sample.int(dim(predG)[1], n_post_draws)
y_g_sim <- matrix(NA,n_post_draws,length(grow_data_list$y))
for(i in 1:n_post_draws){
  for(j in 1:length(grow_data_list$y)){
    y_g_sim[i,j] <- sample(x = 1:max(grow_data_list$y), size = 1, replace = T, prob = dnbinom(1:max(grow_data_list$y), mu = exp(predG[post_draws[i],j]), size = odG[post_draws[i],j])/(1-dnbinom(0, mu = exp(predG[post_draws[i],j]), size = odG[post_draws[i],j])))
  }
}
ppc_dens_overlay(grow_data_list$y, y_g_sim)
ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(0,60) + ggtitle("growth w/o seedling")
ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(50,120)

mean_g_plot <-   ppc_stat(grow_data_list$y, y_g_sim, stat = "mean")
sd_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "sd")
skew_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "skewness")
kurt_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "Lkurtosis")
grid.arrange(mean_g_plot,sd_g_plot,skew_g_plot,kurt_g_plot,  top = "Growth ZTNB")


# now we want to look at how the the growth model is fitting especially across sizes, as this determines the transitions through the matrix model
# data bins

ZTNB_growth_size_ppc <- size_moments_ppc(data = LTREB_data_forgrow,
                                    y_name = "size_t1",
                                    sim = y_g_sim, 
                                    n_bins = 6, 
                                    title = "Growth ZTNB")
# ggsave(ZTNB_growth_size_ppc, filename = "ZTNB_growth_size_pcc.png", width = 4, height = 4)

#growth PIG distribution
grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow_PIG.rds")
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
saveRDS(y_g_sim, file = "yrep_growthPIGmodel.rds")
y_g_sim <- readRDS(file = "yrep_growthPIGmodel.rds")

# Posterior predictive check
ppc_dens_overlay(grow_data_list$y, y_g_sim)
ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(0,60) + ggtitle("growth w/o seedling PIG")
ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(50,120)

grow_densplot <- ppc_dens_overlay(grow_data_list$y, y_g_sim) + xlim(0,60) + theme_classic() + labs(title = "Growth w/o seedling PIG", x = "No. of Tillers", y = "Density")
grow_densplot
ggsave(grow_densplot, filename = "grow_densplot.png", width = 4, height = 4)



mean_g_plot <-   ppc_stat(grow_data_list$y, y_g_sim, stat = "mean")
sd_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "sd")
skew_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "skewness")
kurt_g_plot <- ppc_stat(grow_data_list$y, y_g_sim, stat = "Lkurtosis")
grow_moments <- mean_g_plot+sd_g_plot+skew_g_plot+kurt_g_plot+ plot_annotation(title = "Growth PIG")
grow_moments
ggsave(grow_moments, filename = "grow_momentsplot.png", width = 4, height = 4)


# now we want to look at how the the growth model is fitting especially across sizes, as this determines the transitions through the matrix model
PIG_growth_size_ppc <- size_moments_ppc(data = LTREB_data_forgrow,
                                         y_name = "size_t1",
                                         sim = y_g_sim, 
                                         n_bins = 6, 
                                         title = "Growth PIG")
ggsave(PIG_growth_size_ppc, filename = "PIG_growth_size_pcc.png", width = 4, height = 4)


#### seedling growth ppc ####
# negative binomial
s_grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_grow.rds")
# pairs(grow_fit, pars = c("beta0"))
s_grow_par <- rstan::extract(s_grow_fit, pars = c("lambda","beta0","betaendo","tau_year","tau_plot", "phi", "od"))
predseedG <- s_grow_par$lambda
phiseedG <- s_grow_par$phi
odseedG <- s_grow_par$od

n_post_draws <- 100
post_draws <- sample.int(dim(predseedG)[1], n_post_draws)
y_seed_g_sim <- matrix(NA,n_post_draws,length(seed_grow_data_list$y))
for(i in 1:n_post_draws){
  for(j in 1:length(seed_grow_data_list$y)){
    y_seed_g_sim[i,j] <- sample(x = 1:max(seed_grow_data_list$y), size = 1, replace = T, prob = dnbinom(1:max(seed_grow_data_list$y), mu = exp(predseedG[post_draws[i],j]), size = odseedG[post_draws[i],j])/(1-dnbinom(0, mu = exp(predseedG[post_draws[i],j]), size = odseedG[post_draws[i],j])))
  }
}
ppc_dens_overlay(seed_grow_data_list$y, y_seed_g_sim)
ppc_dens_overlay(seed_grow_data_list$y, y_seed_g_sim) + xlim(0,20) + ggtitle("seedling growth")

mean_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "mean") + ggtitle("mean")
sd_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "sd")+ ggtitle("sd")
skew_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "skewness")+ ggtitle("skew")
kurt_seed_g_plot <- ppc_stat(seed_grow_data_list$y, y_seed_g_sim, stat = "Lkurtosis")+ ggtitle("kurt")
grid.arrange(mean_seed_g_plot,sd_seed_g_plot,skew_seed_g_plot,kurt_seed_g_plot, top = "Seedling Growth")


#seedling growth PIG distribution
s_grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_grow_PIG_10000iterations.rds")
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

saveRDS(y_seed_g_sim, file = "yrep_seedlinggrowthPIGmodel.rds")
y_seed_g_sim <- readRDS(file = "yrep_seedlinggrowthPIGmodel.rds")

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

#### fertility ppc ####
# This fit is maybe not super great, but it's good enough for now probably.
fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_noplot.rds")
fert_par <- rstan::extract(fert_fit, pars = c("lambda","beta0","betasize","betaendo","betaorigin","tau_year", "phi", "od"))
predFert <- fert_par$lambda
phiFert <- fert_par$phi
odFert <- fert_par$od

n_post_draws <- 500
post_draws <- sample.int(dim(predFert)[1], n_post_draws)
y_fert_sim <- matrix(NA,n_post_draws,length(fert_data_list$y))
for(i in 1:n_post_draws){
  for(j in 1:length(fert_data_list$y)){
    y_fert_sim[i,j] <- sample(x = 1:max(fert_data_list$y), size = 1, replace = T, prob = dnbinom(1:max(fert_data_list$y), mu = exp(predFert[post_draws[i],j]), size = odFert[post_draws[i],j])/(1-dnbinom(0, mu = exp(predFert[post_draws[i],j]), size = odFert[post_draws[i],j])))
  }
}
ppc_dens_overlay(fert_data_list$y, y_fert_sim)
ppc_dens_overlay(fert_data_list$y, y_fert_sim) + xlim(0,20)
ppc_dens_overlay(fert_data_list$y, y_fert_sim) + xlim(20,80)

fert_densplot <- ppc_dens_overlay(fert_data_list$y, y_fert_sim) + xlim(0,30) + theme_classic() + labs(title = "Panicles", x = "No. of Panicles", y = "Density")
fert_densplot
ggsave(fert_densplot, filename = "fert_densplot.png", width = 4, height = 4)




mean_fert_plot <-   ppc_stat(fert_data_list$y, y_fert_sim, stat = "mean")
sd_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "sd")
skew_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "skewness")
kurt_fert_plot <- ppc_stat(fert_data_list$y, y_fert_sim, stat = "Lkurtosis")
grid.arrange(mean_fert_plot,sd_fert_plot,skew_fert_plot,kurt_fert_plot)

# Now we can look at the binned size fit for fertility
fert_size_ppc <- size_moments_ppc(data = LTREB_data_forfert,
                                        y_name = "FLW_COUNT_T1",
                                        sim = y_fert_sim, 
                                        n_bins = 6, 
                                        title = "Inflorescence Count")

mcmc_trace(fert_fit,pars=c("betaendo[1]","betaendo[2]","betaendo[3]"
                           ,"betaendo[4]","betaendo[5]","betaendo[6]"
                           ,"betaendo[7]"))


# Fert fit with Poisson inverse Gaussian

fert_fit_pig <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_pig.rds")
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
saveRDS(y_fert_sim, file = "yrep_fertilityPIGmodel.rds")
y_fert_sim <- readRDS(file = "yrep_fertilityPIGmodel.rds")


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
# Looking at the poisson fit
spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike_year_plot_poisson.rds")
spike_pars <- rstan::extract(spike_fit, pars = c("lambda"))
predSpike <- spike_pars$lambda
n_post_draws <- 500
post_draws <- sample.int(dim(predSpike)[1], n_post_draws)
y_spike_sim <- matrix(NA,n_post_draws,length(spike_data_list$y))
for(i in 1:n_post_draws){
  y_spike_sim[i,] <- rpois(n=length(spike_data_list$y), lambda = exp(predSpike[post_draws[i],]))
}
# ppc_dens_overlay(spike_data_list$y, y_spike_sim)
# ppc_dens_overlay(spike_data_list$y, y_spike_sim) + xlim(0,250)

spike_densplot <- ppc_dens_overlay(spike_data_list$y, y_spike_sim) + xlim(0,30) + theme_classic() + labs(title = "Panicles", x = "No. of Panicles", y = "Density")
spike_densplot
ggsave(spike_densplot, filename = "spike_densplot.png", width = 4, height = 4)

# Fit isn't super great for the poisson, but probably close enough for the mean. (Also this is just fit as a poisson, but could look at neg binom and zero truncate)

mean_spike_plot <-   ppc_stat(spike_data_list$y, y_spike_sim, stat = "mean")
sd_spike_plot <- ppc_stat(spike_data_list$y, y_spike_sim, stat = "sd")
skew_spike_plot <- ppc_stat(spike_data_list$y, y_spike_sim, stat = "skewness")
kurt_spike_plot <- ppc_stat(spike_data_list$y, y_spike_sim, stat = "Lkurtosis")
grid.arrange(mean_spike_plot,sd_spike_plot,skew_spike_plot,kurt_spike_plot)


mcmc_trace(spike_fit,pars=c("betaendo[1]","betaendo[2]","betaendo[3]"
                           ,"betaendo[4]","betaendo[5]","betaendo[6]"
                           ,"betaendo[7]"))

# looking at the negative binomial fit, fits really nicely
spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike_year_plot_nb.rds")
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
saveRDS(y_spike_sim, file = "yrep_spikeletNBmodel.rds")
y_spike_sim <- readRDS(file = "yrep_spikeletNBmodel.rds")

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

#########################################################################################################
# Plots for all species all vital rates and model fits ------------------------------
#########################################################################################################

## Plot for all models vital rate fits


fits_plot <- surv_densplot + seedsurv_densplot+
             grow_densplot + seedgrow_densplot+
             flw_densplot+ fert_densplot+
             spike_densplot+plot_spacer() + 
             plot_layout(ncol = 2) + plot_annotation(title = "Vital rate fits with 500 posterior draws")
  
fits_plot
ggsave(fits_plot, filename = "fits_plot.png", width = 18, height = 20)


## Plot for all models moments

moments_plot <- surv_moments + seedsurv_moments+
                grow_moments + seedgrow_moments+
                flw_moments+ fert_moments+
                spike_moments+ plot_spacer()+
                plot_layout(ncol = 2) + plot_annotation(title = "Vital rate moments with distribution of posterior draws")
ggsave(moments_plot, filename = "moments_plot.png", width = 18, height = 20)


## Plot for all models fits and moments
fitsandmoments_plot <- (surv_densplot + surv_moments)/
                      (seedsurv_densplot + seedsurv_moments)/
                      (grow_densplot + grow_moments)/ 
                      (seedgrow_densplot + seedgrow_moments)/
                       (flw_densplot + flw_moments)/
                        (fert_densplot + fert_moments)/
                       (spike_densplot + spike_moments)+ 
                      plot_annotation(title = "Vital rate fits and moments with 500 posterior draws")
ggsave(fitsandmoments_plot, filename = "fitsandmoments_plot.png", width = 18, height = 20)

## Plot for size-specific moments for growth model minus seedlings, which have only one size
## could plot these for other vital rates if desired
size_ppc_plot <- (PIG_growth_size_ppc+ plot_annotation(title = 'Growth'))+
                  plot_annotation(title = "Size specific vital rate moments")
size_ppc_plot
ggsave(size_ppc_plot, filename = "size_ppc_plot.png", width = 12, height = 20)


## Plot of traceplots select parameters for all models

surv_trace <- mcmc_trace(surv_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Survival")
seedsurv_trace <- mcmc_trace(surv_seed_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Seedling Survival")
flw_trace <- mcmc_trace(flow_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Flowering")
grow_trace <- mcmc_trace(grow_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Growth")
seedgrow_trace <- mcmc_trace(s_grow_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Seedling Growth")
fert_trace <- mcmc_trace(fert_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Fertility")
spike_trace <- mcmc_trace(spike_fit, pars = c("beta0[1]", "betaendo[1]", "sigmaendo[1]"))+ggtitle("Spikelets per inflorescence")

vr_traceplots <- (surv_trace)/
                  (seedsurv_trace)/
                  (grow_trace)/
                  (seedgrow_trace)/
                  (flw_trace)/
                  (fert_trace)/
                  (spike_trace) + plot_annotation(title = "Traceplots for select parameters from all vital rates")
ggsave(vr_traceplots, filename = "vr_traceplots.png", width = 25, height = 20)

######## Plots of endophyte mean and variance effects for all vital rates ####
# max size for each species to create a sequence of size values
max_size <- LTREB_full %>% 
  dplyr::select(species,size_t) %>% 
  filter(!is.na(size_t)) %>% 
  group_by(species) %>% 
  summarise(actual_max_size = max(size_t),
            max_size = quantile(size_t,probs=0.975))

# Our parameter estimates are in the _fit stan objects for each vital rate
surv_fit
surv_par <- rstan::extract(surv_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                      tau_year, tau_plot))
n_post_draws <- 100
post_draws <- sample.int(dim(surv_par[[1]])[1], n_post_draws)
source("Analyses/MPM_functions.R")
surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_surv.rds")
surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_surv_woseedling.rds")
grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_seedling_grow_PIG_10000iterations.rds")
grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_grow_PIG.rds")
flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_flw.rds")
fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_fert_PIG.rds")
spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike_year_plot_nb.rds")
seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/seed_mean.rds")
stos_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_s_to_s.rds") 

surv_par <- rstan::extract(surv_fit, pars =quote_bare(beta0,betasize,betaendo,betaorigin,
                                                      tau_year, tau_plot))
surv_sdlg_par <- rstan::extract(surv_fit_seedling, pars =quote_bare(beta0,betaendo,
                                                                    tau_year, tau_plot))
grow_par <- rstan::extract(grow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                       tau_year, tau_plot,
                                                       sigma))
grow_sdlg_par <- rstan::extract(grow_fit_seedling, pars = quote_bare(beta0,betaendo,
                                                                     tau_year, tau_plot,
                                                                     sigma))
flow_par <- rstan::extract(flw_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                      tau_year, tau_plot))
fert_par <- rstan::extract(fert_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                       tau_year, tau_plot))
spike_par <- rstan::extract(spike_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,
                                                         tau_year, tau_plot,
                                                         phi))
seed_par <- rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
recruit_par <- rstan::extract(stos_fit, pars = quote_bare(beta0,betaendo,
                                                          tau_year, tau_plot))

params <- make_params(species=1,
                     endo_mean=(2-1),
                     endo_var=(2-1),
                     original = 1, # should be =1 to represent recruit
                     draw=post_draws[1],
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
                     recruit_par=recruit_par)
for(s in 1:7)
surv_int <- surv_par$beta0[,1] + beta_endo






