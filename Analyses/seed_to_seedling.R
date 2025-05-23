## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates seed to seedling recruitment kernel written in STAN, 
## and does visualisation of posterior predictive checks
## Authors: Joshua and Tom
#############################################################

library(tidyverse)
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)
library(moments)
library(patchwork)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }
Lkurtosis=function(x) log(kurtosis(x)); 


#############################################################################################
####### Data manipulation to prepare data as lists for Stan models------------------
#############################################################################################

source("Analyses/endodemog_data_processing.R")

# Loading in our vital rate model samples to estimate seed production
sm_seed_means <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_seed_mean.rds")
sm_spikelet <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_spike_year_plot_nb.rds")


# I'm going to use LTREB_full to generate estimates of seed production and then calculate our numbers of recruits
seed_mean_pars <- rstan::extract(sm_seed_means, pars = c("mu_seed","beta0", "betaendo", "sigma0"))
spikelet_pars <- rstan::extract(sm_spikelet, pars = c("lambda", "phi",  "od", "beta0", "betasize", "betaendo", "betaorigin", 
                                                      "tau_plot", "sigma_plot", 
                                                      "tau_year", "sigma_year"))


seedmean_prob <- c(NA)
spikeperinf_prob <- c(NA)
spikeperinf_pred <- c(NA)
spike_od <- c(NA)
for(i in 1:nrow(LTREB_full)){
      seedmean_prob[i] <- sample(seed_mean_pars$beta0[,LTREB_full$species_index[i]], size = 1) +
                            sample(seed_mean_pars$betaendo[,LTREB_full$species_index[i]], size = 1)*LTREB_full$endo_01[i]
      
      spikeperinf_prob[i] <- exp(sample(spikelet_pars$beta0[,LTREB_full$species_index[i]], size = 1) +
                               sample(spikelet_pars$betasize[,LTREB_full$species_index[i]], size = 1)*LTREB_full$logsize_t1[i] +
                                 sample(spikelet_pars$betaendo[,LTREB_full$species_index[i]], size = 1)*LTREB_full$endo_01[i] +
                                   sample(spikelet_pars$betaorigin, size = 1)*LTREB_full$origin_01[i] +
                                       sample(spikelet_pars$tau_plot[,LTREB_full$plot_index[i]], size = 1) +
                                         sample(spikelet_pars$tau_year[,LTREB_full$species_index[i],LTREB_full$endo_index[i],LTREB_full$year_t_index[i]], size = 1))
      spike_od[i] <-   exp(sample(spikelet_pars$phi[,LTREB_full$species_index[i]], size = 1))

}
LTREB_full$spikeperinf_prob_t1 <- spikeperinf_prob
LTREB_full$spikeperinf_pred_t1 <- rnbinom(n = nrow(LTREB_full), mu = spikeperinf_prob, size = spike_od)
LTREB_full$seedmean_pred_t1 <- seedmean_prob


#My spikelet predictions are higher than seen in the data:
# max(filter(LTREB_full, species == "AGPE")$SPIKE_A_T1, na.rm = T)
# max(filter(LTREB_full, species == "AGPE")$spikeperinf_pred_t1, na.rm = T)
# 
# max(filter(LTREB_full, species == "ELVI")$SPIKE_A_T1, na.rm = T)
# max(filter(LTREB_full, species == "ELVI")$spikeperinf_pred_t1, na.rm = T)

# plot(LTREB_full$logsize_t1, LTREB_full$spikeperinf_pred, col = LTREB_full$species_index)
# plot( LTREB_full$spikeperinf_pred_t1, LTREB_full$SPIKE_A_T1)
# abline(0,1)
# # 
# hist(LTREB_full$spikeperinf_pred_t1)
# hist(LTREB_full$SPIKE_A_T1)
# # hist(LTREB_full$SPIKE_AGPE_MEAN_T1)


# Here we want to generate a estimate of the seed production for each plant. We are using the our reproductive data and our seed mean model. Where we don't have spikelet data, we use the spikelet model to fill in our estimates
# There are 12 rows where we have size data but not flower data, most of which are 1 or 2 tiller plants. We probably forgot to record the flowering data, and many are likely to be non-flowering. Given the small number, I'm dropping these.
# Then I want to get a value for the total number of seeds for each plot for each year.
LTREB_annual_seed_data <- LTREB_full %>% 
  mutate(endo_01 = case_when(plot_index == 30 | plot_index == 23 ~ 1,
                             TRUE ~ endo_01),
         endo_index = case_when(plot_index == 30 | plot_index == 23 ~ 2,
                             TRUE ~ endo_index)) %>% # fixing some typos for plot endophyte status for two LOAR plots
  rowwise() %>%
  mutate(spike_t1_mean = mean(c(SPIKE_A_T1, SPIKE_B_T1, SPIKE_C_T1, SPIKE_AGPE_MEAN_T1), na.rm = T)) %>% # getting a spikelet mean for each plant where we actually have data
  mutate(seed_est_t1 = case_when( is.na(SPIKE_A_T1) & !is.na(FLW_STAT_T1) ~ round(FLW_STAT_T1*FLW_COUNT_T1*spikeperinf_pred_t1*seedmean_pred_t1),
                                  !is.na(SPIKE_A_T1) & !is.na(FLW_STAT_T1) ~ round(FLW_STAT_T1*FLW_COUNT_T1*spike_t1_mean*seedmean_pred_t1))) %>%  # This is where we fill in the spikelet and seed mean predictions to get our seed estimates.
  mutate(year_t1_factor = factor(year_t1_index, levels = min(year_t1_index):max(year_t1_index)),
         plot_factor = factor(plot_index, levels = min(plot_index):max(plot_index))) %>% 
  group_by(species, species_index, plot_index, year_t, year_t_index, year_t1, year_t1_index, endo_01, endo_index) %>% 
  summarize(seeds_n_rows = n(),
            seeds_surviving_plants_t1 = sum(surv_t1, na.rm = T),
            tot_seed_t1 = sum(seed_est_t1, na.rm = T)) %>% 
  ungroup() %>% 
  complete(nesting(species, species_index, plot_index, endo_01, endo_index), nesting(year_t1, year_t1_index, year_t, year_t_index), 
           fill =  list(tot_seed_t1 = 0)) %>% 
  mutate(tot_seed_t = dplyr::lag(tot_seed_t1, order_by = plot_index)) %>% 
  mutate(tot_seed_t = case_when(year_t == 2007 ~ 0,
                                year_t != 2007 ~ tot_seed_t))

LTREB_annual_recruit_data <- LTREB_full %>% 
  filter(origin_01 == "1" & year_t == birth) %>% 
  mutate(endo_01 = case_when(plot_index == 30 | plot_index == 23 ~ 1,
                             TRUE ~ endo_01),
         endo_index = case_when(plot_index == 30 | plot_index == 23 ~ 2,
                                TRUE ~ endo_index)) %>% # fixing some typos for plot endophyte status for two LOAR plots
  mutate(surv_t = case_when(!is.na(size_t) ~ 1,
                           is.na(size_t) ~ 0)) %>% 
  filter(surv_t == 1) %>%       # There are a few plants which were TNF or other issues that weren't actually recruits, so I'm dropping these (~22 plants)
  filter(size_t <= 4) %>%   # dropping all the plants bigger than 4 tillers. There are ~13000 recruits that are size 1 tiller, and 682 "recruits" that are bigger. Some of these we have notes saying they arer probablly older. For currrent data collection, we are leaviing the birth year as NA, instead of adding a definite year. Filtering out at 4 tillers drops 142 plants, which allows for the occasional recruit that is actually bigger than 1 tiller.
  group_by(species, species_index, plot_index, year_t, year_t_index, year_t1, year_t1_index, endo_01, endo_index) %>% 
  summarize(recruits_n_rows = n(),
            recruits_surviving_plants_t1 = sum(surv_t1, na.rm = T),
            tot_recruits_t = sum(surv_t, na.rm = T)) %>% 
  ungroup() %>% 
  complete(nesting(species, species_index, plot_index, endo_01, endo_index), nesting( year_t1, year_t1_index, year_t, year_t_index), 
           fill = list(tot_recruits_t = 0)) %>% 
  mutate(tot_recruits_t1 = dplyr::lead(tot_recruits_t, order_by = plot_index))

# getting the spei data for each species and year 
# LOAR only goes up to 2019
LTREB_spei <- LTREB_full %>% 
  ungroup() %>% 
  dplyr::select(species, species_index, year_t1, year_t1_index, year_t, year_t_index, spei12, spei3) %>% 
  unique()

# These are the problem recruit id's
# View(filter(LTREB_annual_recruit_data, surv_t == 0))
# Should also consider dropping the plants that are big and aren't actually recruits. In new data, i'm trying to leave the birth year NA, but some of the previous have assigned the birth year, with a note saying it's probably an older plant. There aren't a ton, but maybe worth dropping the really big plants
# unique(LTREB_annual_recruit_data$size_t)
# table(LTREB_annual_recruit_data$size_t)
     

  
LTREB_s_to_s_data <- LTREB_annual_seed_data %>% 
  left_join(LTREB_annual_recruit_data, by = c("species", "species_index", "plot_index", "year_t", "year_t_index", "year_t1", "year_t1_index", "endo_01", "endo_index")) %>% 
  left_join(LTREB_spei,  by = c("species", "species_index", "year_t", "year_t_index", "year_t1", "year_t1_index")) %>% 
  mutate(tot_recruits_t1 = case_when(is.na(tot_recruits_t1) ~ 0, 
                                      TRUE ~ tot_recruits_t1))  %>%      # Giving a zero value of recruits for 2007 and 2008. This pretty true, I think, but there may be some instance where the birth years of those earlly recruits are just weird because the Original plants all have birth year 2007
  filter(tot_seed_t >= tot_recruits_t1) %>%   # There are ~119 rows where there are recruits, but no seeds from previous year, or more recruits than seeds. Think about a seed bank
  filter(!is.na(spei12)) # this drops the spei values for LOAR where we don't have any actual data. This is actually good because they were being labeled as 0 values for seed and recruits

dim(LTREB_s_to_s_data)

# Create data lists to be used for the Stan model
# We are using the same predicted seed values for the climate-implicit and -explicit versions of the models
s_to_s_data_list <- list(tot_recruit_t1 = LTREB_s_to_s_data$tot_recruits_t1,
                         tot_seed_t = LTREB_s_to_s_data$tot_seed_t,
                         endo_01 = as.integer(LTREB_s_to_s_data$endo_01),
                         endo_index = as.integer(LTREB_s_to_s_data$endo_index),
                         year_t = as.integer(LTREB_s_to_s_data$year_t_index),
                         plot = as.integer(LTREB_s_to_s_data$plot_index),
                         spp = LTREB_s_to_s_data$species_index,
                         N = nrow(LTREB_s_to_s_data),
                         nYear = as.integer(max(unique(LTREB_s_to_s_data$year_t_index))),
                         nPlot = length(unique(LTREB_s_to_s_data$plot_index)),
                         nSpp = length(unique(LTREB_s_to_s_data$species_index)),
                         nEndo = length(unique(LTREB_s_to_s_data$endo_01)))
s_to_s_spei12_data_list <- append(s_to_s_data_list, list(spei = as.numeric(LTREB_s_to_s_data$spei12),
                                                               spei_nl = as.numeric(LTREB_s_to_s_data$spei12^2)))
s_to_s_spei3_data_list <- append(s_to_s_data_list, list(spei = as.numeric(LTREB_s_to_s_data$spei3),
                                                              spei_nl = as.numeric(LTREB_s_to_s_data$spei3^2)))


str(s_to_s_data_list);str(s_to_s_spei12_data_list);str(s_to_s_spei3_data_list)

# saveRDS(s_to_s_data_list, file = "Analyses/s_to_s_data_list.rds")
# saveRDS(s_to_s_spei12_data_list, file = "Analyses/s_to_s_spei12_data_list.rds")
# saveRDS(s_to_s_spei3_data_list, file = "Analyses/s_to_s_spei3_data_list.rds")

s_to_s_data_list <- read_rds(file = "Analyses/s_to_s_data_list.rds")
s_to_s_spei12_data_list <- read_rds(file = "Analyses/s_to_s_spei12_data_list.rds")
s_to_s_spei3_data_list <- read_rds(file = "Analyses/s_to_s_spei3_data_list.rds")


#########################################################################################################
# Stan model for seed to seedling recruitment rate ------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
mcmc_pars <- list(
  warmup = 2500, 
  iter =5000, 
  thin = 1, 
  chains = 3
)

# Stan model -------------

## Run the model by calling stan()
## Save the outputs as rds files
# sm_s_to_s <- stanc(file = "Analyses/endo_spp_s_to_s.stan")
# Fits with max_treedepth warning, so increased max_treedepth
sm_s_to_s <- stan(file = "Analyses/endo_spp_s_to_s.stan", data = s_to_s_data_list,
                     iter = mcmc_pars$iter,
                     warmup = mcmc_pars$warmup,
                     chains = mcmc_pars$chains, 
                     thin = mcmc_pars$thin,
                  control = list(max_treedepth = 15))
saveRDS(sm_s_to_s, file = "~/Dropbox/EndodemogData/Model_Runs/endo_spp_s_to_s.rds")

print(sm_s_to_s, pars = c("sigmaendo"))
traceplot(sm_s_to_s, pars = c("beta0"))


#running the climate-explicit version with just the linear climate effect
# running this gives maximum treedepth warnings
sm_s_to_s_spei12 <- stan(file = "Analyses/climate_endo_spp_s_to_s.stan", data = s_to_s_spei12_data_list,
                  iter = mcmc_pars$iter,
                  warmup = mcmc_pars$warmup,
                  chains = mcmc_pars$chains, 
                  thin = mcmc_pars$thin,
                  control = list(max_treedepth = 15))
saveRDS(sm_s_to_s_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_s_to_s_spei12.rds")

sm_s_to_s_spei3 <- stan(file = "Analyses/climate_endo_spp_s_to_s.stan", data = s_to_s_spei3_data_list,
                         iter = mcmc_pars$iter,
                         warmup = mcmc_pars$warmup,
                         chains = mcmc_pars$chains, 
                         thin = mcmc_pars$thin,
                         control = list(max_treedepth = 15))
saveRDS(sm_s_to_s_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_s_to_s_spei3.rds")

#########################################################################################################
# Model Diagnostics ------------------------------
#########################################################################################################

s_to_s_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_s_to_s.rds") 
s_to_s_par <- rstan::extract(s_to_s_fit, pars = c("p"))$p
predRecruit<- s_to_s_par

n_post_draws <- 500
post_draws <- sample.int(dim(predRecruit)[1], n_post_draws)
y_recruit_sim <- matrix(NA,n_post_draws,length(s_to_s_data_list$tot_recruit_t1))

for(i in 1:n_post_draws){
  y_recruit_sim[i,] <- rbinom(n = length(s_to_s_data_list$tot_recruit_t1), size = s_to_s_data_list$tot_seed_t, prob = invlogit(predRecruit[post_draws[i],]))
}
saveRDS(y_recruit_sim, file = "yrep_stosmodel.rds")
y_recruit_sim <- read_rds(file = "yrep_stosmodel.rds")

ppc_dens_overlay(s_to_s_data_list$tot_recruit_t1, y_recruit_sim)
ppc_dens_overlay(s_to_s_data_list$tot_recruit_t1, y_recruit_sim) +xlim(0,30) # seems to be fitting okay
stos_densplot <- ppc_dens_overlay(s_to_s_data_list$tot_recruit_t1, y_recruit_sim) +xlim(0,40) + labs(title = "Recruitment", x = "Successful Germination", y = "Density") 
ggsave(stos_densplot, filename = "stos_densplot.png", width = 4, height = 4)

# overall mean looks okay, but not great

mean_stos_plot <-   ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "mean")
sd_stos_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "sd")
skew_stos_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "skewness")
kurt_stos_plot <- ppc_stat(s_to_s_data_list$tot_recruit_t1, y_recruit_sim, stat = "Lkurtosis")
stos_moments <- mean_stos_plot+sd_stos_plot+skew_stos_plot+kurt_stos_plot
stos_moments
ggsave(stos_moments, filename = "stos_moments.png", width = 4, height = 4)


# PPC for the climate-explicit recruitment model
s_to_s_fit_spei12 <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_s_to_s_spei12.rds") 
s_to_s_fit_spei3 <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_s_to_s_spei3.rds") 

s_to_s_par_spei12 <- rstan::extract(s_to_s_fit_spei12, pars = c("p"))$p
s_to_s_par_spei3 <- rstan::extract(s_to_s_fit_spei3, pars = c("p"))$p

predRecruit_spei12<- s_to_s_par_spei12
predRecruit_spei3<- s_to_s_par_spei3

n_post_draws <- 500
post_draws <- sample.int(dim(predRecruit_spei12)[1], n_post_draws)
y_recruit_sim_spei12 <- y_recruit_sim_spei3 <- matrix(NA,n_post_draws,length(s_to_s_spei12_data_list$tot_recruit_t1))

for(i in 1:n_post_draws){
  y_recruit_sim_spei12[i,] <- rbinom(n = length(s_to_s_spei12_data_list$tot_recruit_t1), size = s_to_s_spei12_data_list$tot_seed_t, prob = invlogit(predRecruit_spei12[post_draws[i],]))
  y_recruit_sim_spei3[i,] <- rbinom(n = length(s_to_s_spei3_data_list$tot_recruit_t1), size = s_to_s_spei3_data_list$tot_seed_t, prob = invlogit(predRecruit_spei3[post_draws[i],]))
  
}
saveRDS(y_recruit_sim_spei12, file = "yrep_climate_stosmodel_spei12.rds")
saveRDS(y_recruit_sim_spei3, file = "yrep_climate_stosmodel_spei3.rds")

y_recruit_sim_spei12 <- read_rds(file = "yrep_climate_stosmodel_spei12.rds")
y_recruit_sim_spei3 <- read_rds(file = "yrep_climate_stosmodel_spei3.rds")

ppc_dens_overlay(s_to_s_spei12_data_list$tot_recruit_t1, y_recruit_sim_spei12)+xlim(0,30)# seems to be fitting okay
stos_spei12_densplot <- ppc_dens_overlay(s_to_s_spei12_data_list$tot_recruit_t1, y_recruit_sim_spei12) +xlim(0,30) + labs(title = "Recruitment (12 month SPEI)", x = "Successful Germination", y = "Density") 
ggsave(stos_spei12_densplot, filename = "climate_stos_spei12_densplot.png", width = 4, height = 4)

stos_spei3_densplot <- ppc_dens_overlay(s_to_s_spei12_data_list$tot_recruit_t1, y_recruit_sim_spei12) +xlim(0,30) + labs(title = "Recruitment (12 month SPEI)", x = "Successful Germination", y = "Density") 
ggsave(stos_spei3_densplot, filename = "climate_stos_spei3_densplot.png", width = 4, height = 4)

# overall mean looks okay.

mean_stos_spei12_plot <-   ppc_stat(s_to_s_spei12_data_list$tot_recruit_t1, y_recruit_sim_spei12, stat = "mean")
sd_stos_spei12_plot <- ppc_stat(s_to_s_spei12_data_list$tot_recruit_t1, y_recruit_sim_spei12, stat = "sd")
skew_stos_spei12_plot <- ppc_stat(s_to_s_spei12_data_list$tot_recruit_t1, y_recruit_sim_spei12, stat = "skewness")
kurt_stos_spei12_plot <- ppc_stat(s_to_s_spei12_data_list$tot_recruit_t1, y_recruit_sim_spei12, stat = "Lkurtosis")
stos_spei12_moments <- mean_stos_spei12_plot+sd_stos_spei12_plot+skew_stos_spei12_plot+kurt_stos_spei12_plot
# stos_spei12_moments
ggsave(stos_spei12_moments, filename = "climate_stos_spei12_moments.png", width = 4, height = 4)


mean_stos_spei3_plot <-   ppc_stat(s_to_s_spei3_data_list$tot_recruit_t1, y_recruit_sim_spei3, stat = "mean")
sd_stos_spei3_plot <- ppc_stat(s_to_s_spei3_data_list$tot_recruit_t1, y_recruit_sim_spei3, stat = "sd")
skew_stos_spei3_plot <- ppc_stat(s_to_s_spei3_data_list$tot_recruit_t1, y_recruit_sim_spei3, stat = "skewness")
kurt_stos_spei3_plot <- ppc_stat(s_to_s_spei3_data_list$tot_recruit_t1, y_recruit_sim_spei3, stat = "Lkurtosis")
stos_spei3_moments <- mean_stos_spei3_plot+sd_stos_spei3_plot+skew_stos_spei3_plot+kurt_stos_spei3_plot
# stos_spei3_moments
ggsave(stos_spei3_moments, filename = "climate_stos_spei3_moments.png", width = 4, height = 4)


# making a plot for total estimate seed production in each plot
s_to_s_data_df <- as_tibble(s_to_s_data_list) %>% 
  mutate(spp  = case_when(spp == 1 ~ "AGPE",
                          spp == 2 ~ "ELRI",
                          spp == 3 ~ "ELVI",
                          spp == 4 ~ "FESU",
                          spp == 5 ~ "LOAR",
                          spp == 6 ~ "POAL",
                          spp == 7 ~ "POSY",)) %>% 
  filter(tot_seed_t>0)

s_to_s_summary <- s_to_s_data_df %>% 
  group_by(spp, endo_01) %>% 
  summarize(mean_seed = mean(tot_seed_t),
            sd_seed = sd(tot_seed_t),
            mean_recruit = mean(tot_recruit_t1),
            sd_recruit = sd(tot_recruit_t1))

seed_hist <- ggplot(s_to_s_data_df)+
  geom_histogram(aes(x = tot_seed_t, fill = as.factor(endo_01)), alpha = .6)+
  geom_vline(data = s_to_s_summary, aes(xintercept = mean_seed, color = as.factor(endo_01)))+
  facet_wrap(~spp, scales = "free")+
  theme_classic()
ggsave(seed_hist, filename = "seed_hist.png")

seed_series <- ggplot(s_to_s_data_df)+
  geom_point(aes(x = year_t, y = tot_seed_t, color = as.factor(endo_01)), alpha = .6)+
  geom_hline(data = s_to_s_summary, aes(yintercept = mean_seed, color = as.factor(endo_01)))+
  facet_wrap(~spp, scales = "free")+
  theme_classic()
ggsave(seed_series, filename = "seed_series.png")


recruit_hist <- ggplot(s_to_s_data_df)+
  geom_histogram(aes(x = tot_recruit_t1, fill = as.factor(endo_01)), alpha = .6)+
  geom_vline(data = s_to_s_summary, aes(xintercept = mean_recruit, color = as.factor(endo_01)))+
  facet_wrap(~spp, scales = "free")+
  theme_classic()
ggsave(recruit_hist, filename = "recruit_hist.png")


recruit_series <- ggplot(s_to_s_data_df)+
  geom_point(aes(x = year_t, y = tot_recruit_t1, color = as.factor(endo_01)), alpha = .6)+
  geom_hline(data = s_to_s_summary, aes(yintercept = mean_recruit, color = as.factor(endo_01)))+
  facet_wrap(~spp, scales = "free")+
  theme_classic()
ggsave(recruit_series, filename = "recruit_series.png")

  
