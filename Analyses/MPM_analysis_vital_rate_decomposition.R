## Title: Grass endophyte population model with a bayesian framework
## Purpose: Builds matrix population model from vital rate estimates
## and runs stochastic LTRE of vital rate contributions
## Authors: Joshua and Tom
#############################################################

library(tidyverse)
library(scales)
library(bayesplot)
library(popbio)
# library(countreg)
library(RColorBrewer)
library(actuar)
library(rstan)
library(patchwork) # for putting plots together

# library(gridExtra)
# library(grid)
# library(cowplot) # for pulling legend from ggplots
# library(cubelyr) # for working between lists of matrixes and dataframes

euclidean <- function(a, b) sqrt(sum((a - b)^2))

quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

#############################################################################################
####### Read in Data and creating size bins------------------
#############################################################################################

# source("Analyses/endodemog_data_processing.R")
tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"
joshpath <- "~/Dropbox/EndodemogData/"
path<-joshpath

LTREB_full <-  read_csv(paste0(joshpath,"Fulldataplusmetadata/LTREB_full.csv"))
max_size <- LTREB_full %>% 
  dplyr::select(species,species_index, size_t) %>% 
  filter(!is.na(size_t)) %>% 
  group_by(species, species_index) %>% 
  summarise(actual_max_size = max(size_t),
            max_size = quantile(size_t,probs=0.975),
            max_size_99 = quantile(size_t,probs=0.99)) # The mean and sd effects plots look basically identical with either max size

# ggplot(data = LTREB_full)+
#   geom_histogram(aes(size_t)) +
#   geom_vline(data = max_size, aes(xintercept = max_size))+
#   geom_vline(data = max_size, aes(xintercept = max_size_99), col = "red")+
#   facet_wrap(~species, scales = "free")

#############################################################################################
####### Read in matrix population functions ------------------
#############################################################################################

source("Analyses/MPM_functions.R")


#############################################################################################
####### Read in Stan vital rate model outputs ------------------
#############################################################################################

surv_fit_seedling <- read_rds(paste0(path,"/Model_Runs/endo_seedling_surv.rds"))
surv_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_surv_woseedling.rds"))
grow_fit_seedling <- read_rds(paste0(path,"/Model_Runs/endo_seedling_grow_PIG_10000iterations.rds"))
grow_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_grow_PIG.rds"))
flw_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_flw.rds"))
fert_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_fert_PIG.rds"))
spike_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_spike_year_plot_nb.rds"))
seedmean_fit <- read_rds(paste0(path,"/Model_Runs/seed_mean.rds"))
stos_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_s_to_s.rds")) 

# surv_fit_seedling <- readRDS(url("https://www.dropbox.com/s/vf1mju5u4c4fs3t/endo_seedling_surv.rds?dl=1"))
# surv_fit <- readRDS(url("https://www.dropbox.com/s/00bor35inv5dypd/endo_spp_surv_woseedling.rds?dl=1"))
# grow_fit_seedling <- readRDS(url("https://www.dropbox.com/s/m0mw5z29slpm4p7/endo_seedling_grow_PIG_10000iterations.rds?dl=1"))
# grow_fit <- readRDS(url("https://www.dropbox.com/s/0ze8aooi9axj3oq/endo_spp_grow_PIG.rds?dl=1"))
# flw_fit <- readRDS(url("https://www.dropbox.com/s/ej65pn5k0km0z9c/endo_spp_flw.rds?dl=1"))
# fert_fit <- readRDS(url("https://www.dropbox.com/s/pk4x1j97kazu6pb/endo_spp_fert_pig.rds?dl=1"))
# spike_fit <- readRDS(url("https://www.dropbox.com/s/pjgui0n9tng6427/endo_spp_spike_year_plot_nb.rds?dl=1"))
# seedmean_fit <- readRDS(url("https://www.dropbox.com/s/3ma5yc8iusu8bh0/endo_spp_seed_mean.rds?dl=1"))
# stos_fit <- readRDS(url("https://www.dropbox.com/s/nf50hd76iw3hucw/endo_spp_s_to_s.rds?dl=1"))

surv_par <- rstan::extract(surv_fit, pars =quote_bare(beta0,betasize,betaendo,betaorigin,sigma_year,
                                                      tau_year, tau_plot))
surv_sdlg_par <- rstan::extract(surv_fit_seedling, pars =quote_bare(beta0,betaendo,sigma_year,
                                                                    tau_year, tau_plot))
grow_par <- rstan::extract(grow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,sigma_year,
                                                       tau_year, tau_plot,
                                                       sigma))
grow_sdlg_par <- rstan::extract(grow_fit_seedling, pars = quote_bare(beta0,betaendo,sigma_year,
                                                                     tau_year, tau_plot,
                                                                     sigma))
flow_par <- rstan::extract(flw_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,sigma_year,
                                                      tau_year, tau_plot))
fert_par <- rstan::extract(fert_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,sigma_year,
                                                       tau_year, tau_plot))
spike_par <- rstan::extract(spike_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,sigma_year,
                                                         tau_year, tau_plot,
                                                         phi))
seed_par <- rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
recruit_par <- rstan::extract(stos_fit, pars = quote_bare(beta0,betaendo,sigma_year,tau_year, tau_plot))

# endophyte effects on vital rate mean and variance ---------------------------


## SET-UP
# make the list of parameters and calculate mean lambdas
n_draws <- 500 # the means are the same whether we do 500 or 1000 draws
post_draws <- sample.int(7500,size=n_draws) # The models except for seedling growth have 7500 iterations. That one has more (15000 iterations) to help it converge.
n_spp <- length(unique(LTREB_full$species))
n_endo <- 2
# observation years
years_obs <- max(unique(LTREB_full$year_t_index))-1 # we are sampling years 2-14

# VITAL RATE decomposition analysis --------------------------------------------------
## stochastic simulation of lambda_S, with mean/var effects on/off in combination for vital rates (grouped by Surv+Growth vs Reproduction)
## We will sample matrices for observed years
## We added arguments LTRE_U and LTRE_F to turn on and off endophyte effects on mean and variance in the make_params function

## store lambdaS output: 8 species (7+mean),16 scenarios
VRdecomp_lambdaS_obs<-array(NA,dim=c(16,n_spp+1,n_draws))
VRdecomp_lambdaS_obs_extreme2_em<-VRdecomp_lambdaS_obs_extreme6_em<-array(NA,dim=c(16,n_spp+1,n_draws))

# save_VRdecomp_lambda_obs <- array(NA,dim=c(years_obs,16,7,n_draws))
# save_VRdecomp_lambda_obs_extreme2_em <- array(NA,dim=c(2,16,7,n_draws))
# save_VRdecomp_lambda_obs_extreme6_em <- array(NA,dim=c(6,16,7,n_draws))


for(d in 1:n_draws){
  ## list of transition years that we observed
  A_t_obs <-list()
  for(s in 1:n_spp){
    U.VminusMminusF.VminusMminus_list <- U.VplusMplusRVplusVplus_list <- U.VplusMminusF.VminusMminus_list <- U.VplusMplusF.VminusMminus_list <- U.VplusMplusF.VplusMminus_list <- U.VplusMminusF.VplusMminus_list <- U.VplusMminusF.VminusMplus_list <- U.VplusMplusF.VminusMplus_list  <- U.VplusMminusF.VplusMplus_list <- U.VminusMplusF.VminusMminus_list <- U.VminusMplusF.VplusMminus_list <- U.VminusMplusF.VplusMplus_list <- U.VminusMminusF.VplusMplus_list <- U.VminusMminusF.VplusMminus_list <- U.VminusMminusF.VminusMplus_list <- U.VminusMplusF.VminusMplus_list <- list()
    # 1: sample observed years
    for(y in 1:years_obs){ ## 13 transitions matrices (2008-09 through 2020-21)
      U.VminusMminusF.VminusMminus_list[[y]] <- bigmatrix(make_params(species=s,
                                                LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                endo_mean_U=0,
                                                endo_var_U=0,
                                                endo_mean_F=0,
                                                endo_var_F=0,
                                                original = 1, # should be =1 to represent recruit
                                                draw=post_draws[d],
                                                max_size=max_size,
                                                rfx=T,
                                                year=y+1,
                                                surv_par=surv_par,
                                                surv_sdlg_par = surv_sdlg_par,
                                                grow_par=grow_par,
                                                grow_sdlg_par = grow_sdlg_par,
                                                flow_par=flow_par,
                                                fert_par=fert_par,
                                                spike_par=spike_par,
                                                seed_par=seed_par,
                                                recruit_par=recruit_par),
                                    extension = 100)$MPMmat
      U.VplusMplusRVplusVplus_list[[y]] <- bigmatrix(make_params(species=s,
                                               LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                               endo_mean_U=1,
                                               endo_var_U=1,
                                               endo_mean_F=1,
                                               endo_var_F=1,
                                               original = 1, # should be =1 to represent recruit
                                               draw=post_draws[d],
                                               max_size=max_size,
                                               rfx=T,
                                               year=y+1,
                                               surv_par=surv_par,
                                               surv_sdlg_par = surv_sdlg_par,
                                               grow_par=grow_par,
                                               grow_sdlg_par = grow_sdlg_par,
                                               flow_par=flow_par,
                                               fert_par=fert_par,
                                               spike_par=spike_par,
                                               seed_par=seed_par,
                                               recruit_par=recruit_par),
                                   extension = 100)$MPMmat
      U.VplusMminusF.VminusMminus_list[[y]] <- bigmatrix(make_params(species=s,
                                                                     LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                                     endo_mean_U=0,
                                                                     endo_var_U=1,
                                                                     endo_mean_F=0,
                                                                     endo_var_F=0,
                                                                 original = 1, # should be =1 to represent recruit
                                                                 draw=post_draws[d],
                                                                 max_size=max_size,
                                                                 rfx=T,
                                                                 year=y+1,
                                                                 surv_par=surv_par,
                                                                 surv_sdlg_par = surv_sdlg_par,
                                                                 grow_par=grow_par,
                                                                 grow_sdlg_par = grow_sdlg_par,
                                                                 flow_par=flow_par,
                                                                 fert_par=fert_par,
                                                                 spike_par=spike_par,
                                                                 seed_par=seed_par,
                                                                 recruit_par=recruit_par),
                                                     extension = 100)$MPMmat
      U.VplusMplusF.VminusMminus_list[[y]] <- bigmatrix(make_params(species=s,
                                                                    LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                                    endo_mean_U=1,
                                                                    endo_var_U=1,
                                                                    endo_mean_F=0,
                                                                    endo_var_F=0,
                                                                     original = 1, # should be =1 to represent recruit
                                                                     draw=post_draws[d],
                                                                     max_size=max_size,
                                                                     rfx=T,
                                                                     year=y+1,
                                                                     surv_par=surv_par,
                                                                     surv_sdlg_par = surv_sdlg_par,
                                                                     grow_par=grow_par,
                                                                     grow_sdlg_par = grow_sdlg_par,
                                                                     flow_par=flow_par,
                                                                     fert_par=fert_par,
                                                                     spike_par=spike_par,
                                                                     seed_par=seed_par,
                                                                     recruit_par=recruit_par),
                                                         extension = 100)$MPMmat
      U.VplusMplusF.VplusMminus_list[[y]]  <- bigmatrix(make_params(species=s,
                                                              LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                              endo_mean_U=1,
                                                              endo_var_U=1,
                                                              endo_mean_F=1,
                                                              endo_var_F=0,
                                                              original = 1, # should be =1 to represent recruit
                                                              draw=post_draws[d],
                                                              max_size=max_size,
                                                              rfx=T,
                                                              year=y+1,
                                                              surv_par=surv_par,
                                                              surv_sdlg_par = surv_sdlg_par,
                                                              grow_par=grow_par,
                                                              grow_sdlg_par = grow_sdlg_par,
                                                              flow_par=flow_par,
                                                              fert_par=fert_par,
                                                              spike_par=spike_par,
                                                              seed_par=seed_par,
                                                              recruit_par=recruit_par),
                                                  extension = 100)$MPMmat
      U.VplusMminusF.VplusMminus_list[[y]] <- bigmatrix(make_params(species=s,
                                                               LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                               endo_mean_U=0,
                                                               endo_var_U=1,
                                                               endo_mean_F=0,
                                                               endo_var_F=1,
                                                               original = 1, # should be =1 to represent recruit
                                                               draw=post_draws[d],
                                                               max_size=max_size,
                                                               rfx=T,
                                                               year=y+1,
                                                               surv_par=surv_par,
                                                               surv_sdlg_par = surv_sdlg_par,
                                                               grow_par=grow_par,
                                                               grow_sdlg_par = grow_sdlg_par,
                                                               flow_par=flow_par,
                                                               fert_par=fert_par,
                                                               spike_par=spike_par,
                                                               seed_par=seed_par,
                                                               recruit_par=recruit_par),
                                                   extension = 100)$MPMmat
      U.VplusMminusF.VminusMplus_list[[y]] <- bigmatrix(make_params(species=s,
                                                              LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                              endo_mean_U=0,
                                                              endo_var_U=1,
                                                              endo_mean_F=1,
                                                              endo_var_F=0,
                                                              original = 1, # should be =1 to represent recruit
                                                              draw=post_draws[d],
                                                              max_size=max_size,
                                                              rfx=T,
                                                              year=y+1,
                                                              surv_par=surv_par,
                                                              surv_sdlg_par = surv_sdlg_par,
                                                              grow_par=grow_par,
                                                              grow_sdlg_par = grow_sdlg_par,
                                                              flow_par=flow_par,
                                                              fert_par=fert_par,
                                                              spike_par=spike_par,
                                                              seed_par=seed_par,
                                                              recruit_par=recruit_par),
                                                  extension = 100)$MPMmat
      U.VplusMplusF.VminusMplus_list[[y]] <- bigmatrix(make_params(species=s,
                                                             LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                             endo_mean_U=1,
                                                             endo_var_U=1,
                                                             endo_mean_F=1,
                                                             endo_var_F=0,
                                                             original = 1, # should be =1 to represent recruit
                                                             draw=post_draws[d],
                                                             max_size=max_size,
                                                             rfx=T,
                                                             year=y+1,
                                                             surv_par=surv_par,
                                                             surv_sdlg_par = surv_sdlg_par,
                                                             grow_par=grow_par,
                                                             grow_sdlg_par = grow_sdlg_par,
                                                             flow_par=flow_par,
                                                             fert_par=fert_par,
                                                             spike_par=spike_par,
                                                             seed_par=seed_par,
                                                             recruit_par=recruit_par),
                                                 extension = 100)$MPMmat
      U.VplusMminusF.VplusMplus_list[[y]] <- bigmatrix(make_params(species=s,
                                                             LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                             endo_mean_U=0,
                                                             endo_var_U=1,
                                                             endo_mean_F=1,
                                                             endo_var_F=1,
                                                             original = 1, # should be =1 to represent recruit
                                                             draw=post_draws[d],
                                                             max_size=max_size,
                                                             rfx=T,
                                                             year=y+1,
                                                             surv_par=surv_par,
                                                             surv_sdlg_par = surv_sdlg_par,
                                                             grow_par=grow_par,
                                                             grow_sdlg_par = grow_sdlg_par,
                                                             flow_par=flow_par,
                                                             fert_par=fert_par,
                                                             spike_par=spike_par,
                                                             seed_par=seed_par,
                                                             recruit_par=recruit_par),
                                                 extension = 100)$MPMmat
      U.VminusMplusF.VminusMminus_list[[y]] <- bigmatrix(make_params(species=s,
                                                             LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                             endo_mean_U=1,
                                                             endo_var_U=0,
                                                             endo_mean_F=0,
                                                             endo_var_F=0,
                                                             original = 1, # should be =1 to represent recruit
                                                             draw=post_draws[d],
                                                             max_size=max_size,
                                                             rfx=T,
                                                             year=y+1,
                                                             surv_par=surv_par,
                                                             surv_sdlg_par = surv_sdlg_par,
                                                             grow_par=grow_par,
                                                             grow_sdlg_par = grow_sdlg_par,
                                                             flow_par=flow_par,
                                                             fert_par=fert_par,
                                                             spike_par=spike_par,
                                                             seed_par=seed_par,
                                                             recruit_par=recruit_par),
                                                 extension = 100)$MPMmat
      U.VminusMplusF.VplusMminus_list[[y]] <- bigmatrix(make_params(species=s,
                                                               LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                               endo_mean_U=1,
                                                               endo_var_U=0,
                                                               endo_mean_F=0,
                                                               endo_var_F=1,
                                                               original = 1, # should be =1 to represent recruit
                                                               draw=post_draws[d],
                                                               max_size=max_size,
                                                               rfx=T,
                                                               year=y+1,
                                                               surv_par=surv_par,
                                                               surv_sdlg_par = surv_sdlg_par,
                                                               grow_par=grow_par,
                                                               grow_sdlg_par = grow_sdlg_par,
                                                               flow_par=flow_par,
                                                               fert_par=fert_par,
                                                               spike_par=spike_par,
                                                               seed_par=seed_par,
                                                               recruit_par=recruit_par),
                                                   extension = 100)$MPMmat
      U.VminusMplusF.VplusMplus_list[[y]] <-  bigmatrix(make_params(species=s,
                                                           LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                           endo_mean_U=1,
                                                           endo_var_U=0,
                                                           endo_mean_F=1,
                                                           endo_var_F=1,
                                                           original = 1, # should be =1 to represent recruit
                                                           draw=post_draws[d],
                                                           max_size=max_size,
                                                           rfx=T,
                                                           year=y+1,
                                                           surv_par=surv_par,
                                                           surv_sdlg_par = surv_sdlg_par,
                                                           grow_par=grow_par,
                                                           grow_sdlg_par = grow_sdlg_par,
                                                           flow_par=flow_par,
                                                           fert_par=fert_par,
                                                           spike_par=spike_par,
                                                           seed_par=seed_par,
                                                           recruit_par=recruit_par),
                                               extension = 100)$MPMmat
      U.VminusMminusF.VplusMplus_list[[y]] <-  bigmatrix(make_params(species=s,
                                                                LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                                endo_mean_U=0,
                                                                endo_var_U=0,
                                                                endo_mean_F=1,
                                                                endo_var_F=1,
                                                                original = 1, # should be =1 to represent recruit
                                                                draw=post_draws[d],
                                                                max_size=max_size,
                                                                rfx=T,
                                                                year=y+1,
                                                                surv_par=surv_par,
                                                                surv_sdlg_par = surv_sdlg_par,
                                                                grow_par=grow_par,
                                                                grow_sdlg_par = grow_sdlg_par,
                                                                flow_par=flow_par,
                                                                fert_par=fert_par,
                                                                spike_par=spike_par,
                                                                seed_par=seed_par,
                                                                recruit_par=recruit_par),
                                                    extension = 100)$MPMmat
      U.VminusMminusF.VplusMminus_list[[y]] <-  bigmatrix(make_params(species=s,
                                                                 LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                                 endo_mean_U=0,
                                                                 endo_var_U=0,
                                                                 endo_mean_F=0,
                                                                 endo_var_F=1,
                                                                 original = 1, # should be =1 to represent recruit
                                                                 draw=post_draws[d],
                                                                 max_size=max_size,
                                                                 rfx=T,
                                                                 year=y+1,
                                                                 surv_par=surv_par,
                                                                 surv_sdlg_par = surv_sdlg_par,
                                                                 grow_par=grow_par,
                                                                 grow_sdlg_par = grow_sdlg_par,
                                                                 flow_par=flow_par,
                                                                 fert_par=fert_par,
                                                                 spike_par=spike_par,
                                                                 seed_par=seed_par,
                                                                 recruit_par=recruit_par),
                                                     extension = 100)$MPMmat
      U.VminusMminusF.VminusMplus_list[[y]] <-  bigmatrix(make_params(species=s,
                                                                LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                                endo_mean_U=0,
                                                                endo_var_U=0,
                                                                endo_mean_F=1,
                                                                endo_var_F=0,
                                                                original = 1, # should be =1 to represent recruit
                                                                draw=post_draws[d],
                                                                max_size=max_size,
                                                                rfx=T,
                                                                year=y+1,
                                                                surv_par=surv_par,
                                                                surv_sdlg_par = surv_sdlg_par,
                                                                grow_par=grow_par,
                                                                grow_sdlg_par = grow_sdlg_par,
                                                                flow_par=flow_par,
                                                                fert_par=fert_par,
                                                                spike_par=spike_par,
                                                                seed_par=seed_par,
                                                                recruit_par=recruit_par),
                                                    extension = 100)$MPMmat
      U.VminusMplusF.VminusMplus_list[[y]] <-  bigmatrix(make_params(species=s,
                                                               LTRE = T, # This turns on the mean and variance terms specific to survival vs reproduction
                                                               endo_mean_U=1,
                                                               endo_var_U=0,
                                                               endo_mean_F=1,
                                                               endo_var_F=0,
                                                               original = 1, # should be =1 to represent recruit
                                                               draw=post_draws[d],
                                                               max_size=max_size,
                                                               rfx=T,
                                                               year=y+1,
                                                               surv_par=surv_par,
                                                               surv_sdlg_par = surv_sdlg_par,
                                                               grow_par=grow_par,
                                                               grow_sdlg_par = grow_sdlg_par,
                                                               flow_par=flow_par,
                                                               fert_par=fert_par,
                                                               spike_par=spike_par,
                                                               seed_par=seed_par,
                                                               recruit_par=recruit_par),
                                                   extension = 100)$MPMmat
      
    }#y loop
    ## store matrix lists in a list
    A_t_obs[[s]]<-list(U.VminusMminusF.VminusMminus = U.VminusMminusF.VminusMminus_list,
                       U.VplusMplusRVplusVplus = U.VplusMplusRVplusVplus_list,
                       U.VplusMminusF.VminusMminus = U.VplusMminusF.VminusMminus_list,
                       U.VplusMplusF.VminusMminus = U.VplusMplusF.VminusMminus_list,
                       U.VplusMplusF.VplusMminus = U.VplusMplusF.VplusMminus_list,
                       U.VplusMminusF.VplusMminus = U.VplusMminusF.VplusMminus_list, 
                       U.VplusMminusF.VminusMplus = U.VplusMminusF.VminusMplus_list,
                       U.VplusMplusF.VminusMplus = U.VplusMplusF.VminusMplus_list,
                       U.VplusMminusF.VplusMplus = U.VplusMminusF.VplusMplus_list,
                       U.VminusMplusF.VminusMminus = U.VminusMplusF.VminusMminus_list,
                       U.VminusMplusF.VplusMminus = U.VminusMplusF.VplusMminus_list,
                       U.VminusMplusF.VplusMplus = U.VminusMplusF.VplusMplus_list,
                       U.VminusMminusF.VplusMplus = U.VminusMminusF.VplusMplus_list,
                       U.VminusMminusF.VplusMminus = U.VminusMminusF.VplusMminus_list,
                       U.VminusMminusF.VminusMplus = U.VminusMminusF.VminusMplus_list,
                       U.VminusMplusF.VminusMplus = U.VminusMplusF.VminusMplus_list)
    
    ## get lambda by year for E+ and E-
    lambda_t<-matrix(NA,2,years_obs)

    for(i in 1:years_obs){
      lambda_t[1,i]<-lambda(A_t_obs[[s]][[1]][[i]])
      lambda_t[2,i]<-lambda(A_t_obs[[s]][[2]][[i]])
    }
 
    
    ## selects extreme years based on E- only
    toptwo_em <- lambda_t[1,]%in%sort(lambda_t[1,])[c(1,years_obs)]
    topsix_em <- lambda_t[1,]%in%sort(lambda_t[1,])[c(1:3,(years_obs-2):years_obs)]

    for(e in 1:16){
      VRdecomp_lambdaS_obs[e,s,d]<-lambdaSim(mat_list = A_t_obs[[s]][[e]],max_yrs = 500)$lambdaS

      VRdecomp_lambdaS_obs_extreme2_em[e,s,d]<-lambdaSim(mat_list = A_t_obs[[s]][[e]][toptwo_em],max_yrs = 500)$lambdaS
      VRdecomp_lambdaS_obs_extreme6_em[e,s,d]<-lambdaSim(mat_list = A_t_obs[[s]][[e]][topsix_em],max_yrs = 500)$lambdaS
      # We could save the annual lambda values if we chose
      # save_VRdecomp_lambda_obs[,e,s,d] <- sapply(A_t_obs[[s]][[e]], FUN = lambda)
      # 
      # save_VRdecomp_lambda_obs_extreme2_em[,e,s,d] <- sapply(A_t_obs[[s]][[e]][toptwo_em], FUN = lambda)
      # save_VRdecomp_lambda_obs_extreme6_em[,e,s,d] <- sapply(A_t_obs[[s]][[e]][topsix_em], FUN = lambda)

    }
    
  }#s loop
  # Calculating cross species means
  VRdecomp_lambdaS_obs[1,8,d] <- mean(VRdecomp_lambdaS_obs[1,1:7,d]) # species mean U.VminusMminusF.VminusMminus
  VRdecomp_lambdaS_obs[2,8,d] <- mean(VRdecomp_lambdaS_obs[2,1:7,d]) # species mean U.VplusMplusRVplusVplus
  VRdecomp_lambdaS_obs[3,8,d] <- mean(VRdecomp_lambdaS_obs[3,1:7,d]) # species mean U.VplusMminusF.VminusMminus
  VRdecomp_lambdaS_obs[4,8,d] <- mean(VRdecomp_lambdaS_obs[4,1:7,d]) # species mean U.VplusMplusF.VminusMminus
  VRdecomp_lambdaS_obs[5,8,d] <- mean(VRdecomp_lambdaS_obs[5,1:7,d]) # species mean U.VplusMplusF.VplusMminus 
  VRdecomp_lambdaS_obs[6,8,d] <- mean(VRdecomp_lambdaS_obs[6,1:7,d]) # species mean  U.VplusMminusF.VplusMminus
  VRdecomp_lambdaS_obs[7,8,d] <- mean(VRdecomp_lambdaS_obs[7,1:7,d]) # species mean U.VplusMminusF.VminusMplus
  VRdecomp_lambdaS_obs[8,8,d] <- mean(VRdecomp_lambdaS_obs[8,1:7,d]) # species mean U.VplusMplusF.VminusMplus 
  VRdecomp_lambdaS_obs[9,8,d] <- mean(VRdecomp_lambdaS_obs[9,1:7,d]) # species mean U.VplusMminusF.VplusMplus 
  VRdecomp_lambdaS_obs[10,8,d] <- mean(VRdecomp_lambdaS_obs[10,1:7,d]) # species mean U.VminusMplusF.VminusMminus
  VRdecomp_lambdaS_obs[11,8,d] <- mean(VRdecomp_lambdaS_obs[11,1:7,d]) # species mean  U.VminusMplusF.VplusMminus
  VRdecomp_lambdaS_obs[12,8,d] <- mean(VRdecomp_lambdaS_obs[12,1:7,d]) # species mean U.VminusMplusF.VplusMplus 
  VRdecomp_lambdaS_obs[13,8,d] <- mean(VRdecomp_lambdaS_obs[13,1:7,d]) # species mean U.VminusMminusF.VplusMplus
  VRdecomp_lambdaS_obs[14,8,d] <- mean(VRdecomp_lambdaS_obs[14,1:7,d]) # species mean U.VminusMminusF.VplusMminus
  VRdecomp_lambdaS_obs[15,8,d] <- mean(VRdecomp_lambdaS_obs[15,1:7,d]) # species mean  U.VminusMminusF.VminusMplus 
  VRdecomp_lambdaS_obs[16,8,d] <- mean(VRdecomp_lambdaS_obs[16,1:7,d]) # species mean U.VminusMplusF.VminusMplus
  
  
  VRdecomp_lambdaS_obs_extreme2_em[1,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[1,1:7,d]) # species mean U.VminusMminusF.VminusMminus
  VRdecomp_lambdaS_obs_extreme2_em[2,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[2,1:7,d]) # species mean U.VplusMplusRVplusVplus
  VRdecomp_lambdaS_obs_extreme2_em[3,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[3,1:7,d]) # species mean U.VplusMminusF.VminusMminus
  VRdecomp_lambdaS_obs_extreme2_em[4,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[4,1:7,d]) # species mean U.VplusMplusF.VminusMminus
  VRdecomp_lambdaS_obs_extreme2_em[5,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[5,1:7,d]) # species mean U.VplusMplusF.VplusMminus 
  VRdecomp_lambdaS_obs_extreme2_em[6,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[6,1:7,d]) # species mean  U.VplusMminusF.VplusMminus
  VRdecomp_lambdaS_obs_extreme2_em[7,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[7,1:7,d]) # species mean U.VplusMminusF.VminusMplus
  VRdecomp_lambdaS_obs_extreme2_em[8,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[8,1:7,d]) # species mean U.VplusMplusF.VminusMplus 
  VRdecomp_lambdaS_obs_extreme2_em[9,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[9,1:7,d]) # species mean U.VplusMminusF.VplusMplus 
  VRdecomp_lambdaS_obs_extreme2_em[10,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[10,1:7,d]) # species mean U.VminusMplusF.VminusMminus
  VRdecomp_lambdaS_obs_extreme2_em[11,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[11,1:7,d]) # species mean  U.VminusMplusF.VplusMminus
  VRdecomp_lambdaS_obs_extreme2_em[12,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[12,1:7,d]) # species mean U.VminusMplusF.VplusMplus 
  VRdecomp_lambdaS_obs_extreme2_em[13,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[13,1:7,d]) # species mean U.VminusMminusF.VplusMplus
  VRdecomp_lambdaS_obs_extreme2_em[14,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[14,1:7,d]) # species mean U.VminusMminusF.VplusMminus
  VRdecomp_lambdaS_obs_extreme2_em[15,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[15,1:7,d]) # species mean  U.VminusMminusF.VminusMplus 
  VRdecomp_lambdaS_obs_extreme2_em[16,8,d] <- mean(VRdecomp_lambdaS_obs_extreme2_em[16,1:7,d]) # species mean U.VminusMplusF.VminusMplus
  
  VRdecomp_lambdaS_obs_extreme6_em[1,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[1,1:7,d]) # species mean U.VminusMminusF.VminusMminus
  VRdecomp_lambdaS_obs_extreme6_em[2,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[2,1:7,d]) # species mean U.VplusMplusRVplusVplus
  VRdecomp_lambdaS_obs_extreme6_em[3,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[3,1:7,d]) # species mean U.VplusMminusF.VminusMminus
  VRdecomp_lambdaS_obs_extreme6_em[4,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[4,1:7,d]) # species mean U.VplusMplusF.VminusMminus
  VRdecomp_lambdaS_obs_extreme6_em[5,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[5,1:7,d]) # species mean U.VplusMplusF.VplusMminus 
  VRdecomp_lambdaS_obs_extreme6_em[6,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[6,1:7,d]) # species mean  U.VplusMminusF.VplusMminus
  VRdecomp_lambdaS_obs_extreme6_em[7,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[7,1:7,d]) # species mean U.VplusMminusF.VminusMplus
  VRdecomp_lambdaS_obs_extreme6_em[8,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[8,1:7,d]) # species mean U.VplusMplusF.VminusMplus 
  VRdecomp_lambdaS_obs_extreme6_em[9,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[9,1:7,d]) # species mean U.VplusMminusF.VplusMplus 
  VRdecomp_lambdaS_obs_extreme6_em[10,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[10,1:7,d]) # species mean U.VminusMplusF.VminusMminus
  VRdecomp_lambdaS_obs_extreme6_em[11,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[11,1:7,d]) # species mean  U.VminusMplusF.VplusMminus
  VRdecomp_lambdaS_obs_extreme6_em[12,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[12,1:7,d]) # species mean U.VminusMplusF.VplusMplus 
  VRdecomp_lambdaS_obs_extreme6_em[13,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[13,1:7,d]) # species mean U.VminusMminusF.VplusMplus
  VRdecomp_lambdaS_obs_extreme6_em[14,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[14,1:7,d]) # species mean U.VminusMminusF.VplusMminus
  VRdecomp_lambdaS_obs_extreme6_em[15,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[15,1:7,d]) # species mean  U.VminusMminusF.VminusMplus 
  VRdecomp_lambdaS_obs_extreme6_em[16,8,d] <- mean(VRdecomp_lambdaS_obs_extreme6_em[16,1:7,d]) # species mean U.VminusMplusF.VminusMplus
  
}

# # Saving all of the simulations
saveRDS(VRdecomp_lambdaS_obs, file = paste0(path,"/Model_Runs/MPM_output/VRdecomp_lambdaS_obs.rds"))
saveRDS(VRdecomp_lambdaS_obs_extreme2_em, file = paste0(path,"/Model_Runs/MPM_output/VRdecomp_lambdaS_obs_extreme2_em.rds"))
saveRDS(VRdecomp_lambdaS_obs_extreme6_em, file = paste0(path,"/Model_Runs/MPM_output/VRdecomp_lambdaS_obs_extreme6_em.rds"))

VRdecomp_lambdaS_obs <- readRDS( file = paste0(path,"/Model_Runs/MPM_output/VRdecomp_lambdaS_obs.rds"))
VRdecomp_lambdaS_obs_extreme2_em <- readRDS( file = paste0(path,"/Model_Runs/MPM_output/VRdecomp_lambdaS_obs_extreme2_em.rds"))
VRdecomp_lambdaS_obs_extreme6_em <- readRDS( file = paste0(path,"/Model_Runs/MPM_output/VRdecomp_lambdaS_obs_extreme6_em.rds"))




# calculate cross species mean and the posterior means of the contributions 
VRdecomp_lambdaS_obs_diff <- array(NA,dim=c(3,8,8,7))

# VRdecomp_lambdaS_obs_diff <- array(NA,dim=c(3,16,8,7))
for(s in 1:8){
  # eplus-eminus
  VRdecomp_lambdaS_obs_diff[1,1,s,1] = mean(VRdecomp_lambdaS_obs[2,s,]-VRdecomp_lambdaS_obs[1,s,], na.rm = T) # eplus-eminus
  VRdecomp_lambdaS_obs_diff[1,1,s,2:7] = quantile(VRdecomp_lambdaS_obs[2,s,]-VRdecomp_lambdaS_obs[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[2,1,s,1] = mean(VRdecomp_lambdaS_obs_extreme2_em[2,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], na.rm = T) # eplus-eminus for extreme2_em
  VRdecomp_lambdaS_obs_diff[2,1,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme2_em[2,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[3,1,s,1] = mean(VRdecomp_lambdaS_obs_extreme6_em[2,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], na.rm = T) # eplus-eminus for extreme6_em
  VRdecomp_lambdaS_obs_diff[3,1,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme6_em[2,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  

  # eplus mean only - eminus
  VRdecomp_lambdaS_obs_diff[1,2,s,1] = mean(VRdecomp_lambdaS_obs[16,s,]-VRdecomp_lambdaS_obs[1,s,], na.rm = T) # eplus mean only - eminus
  VRdecomp_lambdaS_obs_diff[1,2,s,2:7] = quantile(VRdecomp_lambdaS_obs[16,s,]-VRdecomp_lambdaS_obs[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[2,2,s,1] = mean(VRdecomp_lambdaS_obs_extreme2_em[16,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], na.rm = T) # eplus mean only - eminus for extreme2_em
  VRdecomp_lambdaS_obs_diff[2,2,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme2_em[16,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[3,2,s,1] = mean(VRdecomp_lambdaS_obs_extreme6_em[16,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], na.rm = T)# eplus mean only - eminus for extreme6_em
  VRdecomp_lambdaS_obs_diff[3,2,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme6_em[16,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  
  # eplus var only - eminus
  VRdecomp_lambdaS_obs_diff[1,3,s,1] = mean(VRdecomp_lambdaS_obs[6,s,]-VRdecomp_lambdaS_obs[1,s,], na.rm = T) #  eplus var only - eminus
  VRdecomp_lambdaS_obs_diff[1,3,s,2:7] = quantile(VRdecomp_lambdaS_obs[6,s,]-VRdecomp_lambdaS_obs[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[2,3,s,1] = mean(VRdecomp_lambdaS_obs_extreme2_em[6,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], na.rm = T) #  eplus var only - eminusfor extreme2_em
  VRdecomp_lambdaS_obs_diff[2,3,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme2_em[6,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[3,3,s,1] = mean(VRdecomp_lambdaS_obs_extreme6_em[6,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], na.rm = T)#  eplus var only - eminus for extreme6_em
  VRdecomp_lambdaS_obs_diff[3,3,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme6_em[6,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  
  
  # mean var interaction
  VRdecomp_lambdaS_obs_diff[1,4,s,] = VRdecomp_lambdaS_obs_diff[1,1,s,]-VRdecomp_lambdaS_obs_diff[1,2,s,]-VRdecomp_lambdaS_obs_diff[1,3,s,]# mean variance interaction
  VRdecomp_lambdaS_obs_diff[2,4,s,] = VRdecomp_lambdaS_obs_diff[2,1,s,]-VRdecomp_lambdaS_obs_diff[2,2,s,]-VRdecomp_lambdaS_obs_diff[2,3,s,]# mean variance interaction for extreme2_em
  VRdecomp_lambdaS_obs_diff[3,4,s,] = VRdecomp_lambdaS_obs_diff[3,1,s,]-VRdecomp_lambdaS_obs_diff[3,2,s,]-VRdecomp_lambdaS_obs_diff[3,3,s,]# mean variance interaction for extreme6_em
  
  # eplus SURVmean only - eminus
  VRdecomp_lambdaS_obs_diff[1,5,s,1] = mean(VRdecomp_lambdaS_obs[10,s,]-VRdecomp_lambdaS_obs[1,s,], na.rm = T) 
  VRdecomp_lambdaS_obs_diff[1,5,s,2:7] = quantile(VRdecomp_lambdaS_obs[10,s,]-VRdecomp_lambdaS_obs[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[2,5,s,1] = mean(VRdecomp_lambdaS_obs_extreme2_em[10,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], na.rm = T) # for extreme2_em
  VRdecomp_lambdaS_obs_diff[2,5,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme2_em[10,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[3,5,s,1] = mean(VRdecomp_lambdaS_obs_extreme6_em[10,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], na.rm = T) # for extreme6_em
  VRdecomp_lambdaS_obs_diff[3,5,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme6_em[10,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  
  # eplus SURV var only - eminus
  VRdecomp_lambdaS_obs_diff[1,6,s,1] = mean(VRdecomp_lambdaS_obs[3,s,]-VRdecomp_lambdaS_obs[1,s,], na.rm = T) 
  VRdecomp_lambdaS_obs_diff[1,6,s,2:7] = quantile(VRdecomp_lambdaS_obs[3,s,]-VRdecomp_lambdaS_obs[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[2,6,s,1] = mean(VRdecomp_lambdaS_obs_extreme2_em[3,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], na.rm = T) # for extreme2_em
  VRdecomp_lambdaS_obs_diff[2,6,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme2_em[3,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[3,6,s,1] = mean(VRdecomp_lambdaS_obs_extreme6_em[3,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], na.rm = T) # for extreme6_em
  VRdecomp_lambdaS_obs_diff[3,6,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme6_em[3,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)

  # eplus Repro mean only - eminus
  VRdecomp_lambdaS_obs_diff[1,7,s,1] = mean(VRdecomp_lambdaS_obs[15,s,]-VRdecomp_lambdaS_obs[1,s,], na.rm = T) 
  VRdecomp_lambdaS_obs_diff[1,7,s,2:7] = quantile(VRdecomp_lambdaS_obs[15,s,]-VRdecomp_lambdaS_obs[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[2,7,s,1] = mean(VRdecomp_lambdaS_obs_extreme2_em[15,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], na.rm = T) # for extreme2_em
  VRdecomp_lambdaS_obs_diff[2,7,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme2_em[15,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[3,7,s,1] = mean(VRdecomp_lambdaS_obs_extreme6_em[15,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], na.rm = T) # for extreme6_em
  VRdecomp_lambdaS_obs_diff[3,7,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme6_em[15,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  
  # eplus Repro Var only - eminus
  VRdecomp_lambdaS_obs_diff[1,8,s,1] = mean(VRdecomp_lambdaS_obs[14,s,]-VRdecomp_lambdaS_obs[1,s,], na.rm = T) 
  VRdecomp_lambdaS_obs_diff[1,8,s,2:7] = quantile(VRdecomp_lambdaS_obs[14,s,]-VRdecomp_lambdaS_obs[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[2,8,s,1] = mean(VRdecomp_lambdaS_obs_extreme2_em[14,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], na.rm = T) # for extreme2_em
  VRdecomp_lambdaS_obs_diff[2,8,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme2_em[14,s,]-VRdecomp_lambdaS_obs_extreme2_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  VRdecomp_lambdaS_obs_diff[3,8,s,1] = mean(VRdecomp_lambdaS_obs_extreme6_em[14,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], na.rm = T) # for extreme6_em
  VRdecomp_lambdaS_obs_diff[3,8,s,2:7] = quantile(VRdecomp_lambdaS_obs_extreme6_em[14,s,]-VRdecomp_lambdaS_obs_extreme6_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  
}
################################################################
##### Plot of stochastic lambda contributions
################################################################

# Set color scheme based on analine blue
endophyte_color_scheme <- c("#fdedd3","#f3c8a8", "#5a727b", "#4986c7", "#181914",  "#163381")
color_scheme_set(endophyte_color_scheme)
# color_scheme_view()

simulation_color_scheme <- c( "#163381", "#4986c7", "#bdc9e1", "#f1eef6")
# scales::show_col(simulation_color_scheme)

# And creating a color palette for each year
yearcount = length(unique(LTREB_full$year_t))
yearcolors<- colorRampPalette(brewer.pal(8,"Dark2"))(yearcount)
# scales::show_col(yearcolors)
species_list <- c("A. perennans", "E. villosus", "E. virginicus", "F. subverticillata", "L. arundinaceae", "P. alsodes", "P. sylvestris")
species_code_list <- c("AGPE", "ELRI", "ELVI", "FESU", "LOAR", "POAL", "POSY")

dimnames(VRdecomp_lambdaS_obs_diff) <- list(Scenario = c("Ambient variability","Most variability (2 best and worst years)","More variability (6 best and worst years)"), Contribution = c("Full Effect","Mean only","Variance only","Interaction", "Surv&Growth - Mean only", "Surv&Growth - Var only", "Reproduction - Mean only", "Reproduction - Var only"), Species = c(species_list, "Species Mean"), Quantile = c("mean","fifth","twelvepointfive","twentyfifth","seventyfifth","eightysevenpointfive","ninetyfifth"))
VRdecomp_lambdaS_obs_diff_cube <- cubelyr::as.tbl_cube(VRdecomp_lambdaS_obs_diff)
VRdecomp_lambdaS_obs_diff_df <- as_tibble(VRdecomp_lambdaS_obs_diff_cube) %>% 
  pivot_wider(names_from = "Quantile", values_from = VRdecomp_lambdaS_obs_diff) %>% 
  mutate(vr_or_full = case_when(Contribution %in% c("Full Effect","Mean only","Variance only","Interaction") ~ "Full",
                                Contribution %in% c("Surv&Growth - Mean only", "Surv&Growth - Var only", "Reproduction - Mean only", "Reproduction - Var only") ~ "Vital Rate")) %>% 
  filter(Contribution %in% c("Full Effect","Surv&Growth - Mean only", "Surv&Growth - Var only", "Reproduction - Mean only", "Reproduction - Var only"))


x_levels <- c("Full Effect", "Surv&Growth - Mean only", "Surv&Growth - Var only", "Reproduction - Mean only", "Reproduction - Var only")
VRdecomp_contributions_obs_plot <- ggplot(data = VRdecomp_lambdaS_obs_diff_df) +
  geom_hline(yintercept = 0, col = "black") +
  geom_linerange(aes(x = Contribution, ymin = fifth, ymax = ninetyfifth, group = Scenario), color = "grey38",position = position_dodge(width = .7), lwd = .8) +
  geom_linerange(aes(x = Contribution, ymin = twelvepointfive, ymax = eightysevenpointfive, group = Scenario),color = "grey38", position = position_dodge(width = .7), lwd = 1.5) +
  geom_linerange(aes(x = Contribution, ymin = twentyfifth, ymax = seventyfifth, group = Scenario),color = "grey38", position = position_dodge(width = .7), lwd = 2.4) +
  
  geom_linerange(aes(x = Contribution, ymin = fifth, ymax = ninetyfifth, color = Scenario),position = position_dodge(width = .7), lwd = .4) +
  geom_linerange(aes(x = Contribution, ymin = twelvepointfive, ymax = eightysevenpointfive, color = Scenario),position = position_dodge(width = .7), lwd = 1) +
  geom_linerange(aes(x = Contribution, ymin = twentyfifth, ymax = seventyfifth, color = Scenario),position = position_dodge(width = .7), lwd = 2) +
  
  geom_point(aes(y = mean, x = Contribution, group = Scenario), color = "grey38", position = position_dodge(width = .7), size = 3.9) +
  geom_point(aes(y = mean, x = Contribution, color = Scenario), position = position_dodge(width = .7), size = 3.5) +
  facet_wrap(~Species, nrow = 2, scales = "free") + coord_flip() +
  # scale_shape_manual(values = c(16,18,15,17))+
  scale_fill_manual(values = c(simulation_color_scheme), labels = ~ stringr::str_wrap(.x, width = 25))+
  scale_color_manual(values = c(simulation_color_scheme), labels = ~ stringr::str_wrap(.x, width = 25))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  scale_x_discrete(limits = rev(x_levels))+
  labs( x = "", y = expression(paste("Symbiosis effect on", " ", lambda["s"])))+
  guides(pch = "none")+
  theme(panel.background = element_blank(),
        plot.background = element_rect(color = "white"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(linewidth = .5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = rel(1.2), face = "bold"),
        axis.text.x = element_text(size = rel(1), face = "bold"),
        axis.title = element_text(size = rel(1.5)),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.3), face = "italic"),
        legend.key = element_rect(fill = NA),
        legend.title=element_text(size=rel(1.2)), 
        legend.text=element_text(size=rel(1)),
        legend.key.size = unit(1.1,"cm"))
# VRdecomp_contributions_obs_plot
ggsave(VRdecomp_contributions_obs_plot, filename = "VRdecomp_contributions_obs_plot.png", width = 16, height = 10)

