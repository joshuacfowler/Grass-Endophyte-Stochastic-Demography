## Title: Grass endophyte population model with a bayesian framework
## Purpose: Builds matrix population model from vital rate estimates
## and runs stochastic LTRE
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


# endophyte effects on lambda mean and variance ---------------------------

## There are two possible ways of running the analysis:
## 1. Sample observation years (10-12 of these depending how you count)
## 2. Sample many hypothetical years from posterior of year effects

## SET-UP
# make the list of parameters and calculate mean lambdas
n_draws <- 500 # the means are the same whether we do 500 or 1000 draws
post_draws <- sample.int(7500,size=n_draws) # The models except for seedling growth have 7500 iterations. That one has more (15000 iterations) to help it converge.
n_spp <- length(unique(LTREB_full$species))
n_endo <- 2
# observation years
years_obs <- max(unique(LTREB_full$year_t_index))-1 # we are sampling years 2-14
lambda_year_obs <- array(dim = c(years_obs,(n_spp+1),n_endo,n_draws))
# sampled years
n_years_samp <- 50
lambda_year_samp <- array(dim = c(n_years_samp,(n_spp+1),n_endo,n_draws))
lambda_year_samp_extr <- array(dim = c(n_years_samp,(n_spp+1),n_endo,n_draws))

# store mean and variance outputs for observed and sampled years
lambda_mean_obs <- lambda_mean_samp <- array(dim = c((n_spp+1),n_endo,n_draws))
lambda_sd_obs <- lambda_sd_samp <- array(dim = c((n_spp+1),n_endo,n_draws))
lambda_cv_obs <- lambda_cv_samp <- array(dim = c((n_spp+1),n_endo,n_draws))

for(i in 1:n_draws){
  for(e in 1:n_endo){
    for(s in 1:n_spp){
      for(y in 1:years_obs){ # loops 1-13, which samples years 2 (2008-2009) through 14 (2020-2021)
        # 1. Sample observation years, calculate mean and variance
        lambda_year_obs[y,s,e,i] <- lambda(bigmatrix(make_params(species=s,
                                                             endo_mean=(e-1),
                                                             endo_var=(e-1),
                                                             original = 1, # should be =1 to represent recruit
                                                             draw=post_draws[i],
                                                             max_size=max_size,
                                                             rfx=T,
                                                             year=y+1, ## this pulls the t to t1 growth-surv transition and t1 reproduction -- see make_params()
                                                             surv_par=surv_par,
                                                             surv_sdlg_par = surv_sdlg_par,
                                                             grow_par=grow_par,
                                                             grow_sdlg_par = grow_sdlg_par,
                                                             flow_par=flow_par,
                                                             fert_par=fert_par,
                                                             spike_par=spike_par,
                                                             seed_par=seed_par,
                                                             recruit_par=recruit_par),
                                                 extension = 100)$MPMmat) # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
      }
      lambda_mean_obs[s,e,i] <- mean(lambda_year_obs[,s,e,i])
      lambda_sd_obs[s,e,i] <- sd(lambda_year_obs[,s,e,i])
      lambda_cv_obs[s,e,i] <- lambda_sd_obs[s,e,i]/lambda_mean_obs[s,e,i]
    
      # 2.Sample many hypothetical years from fitted year variances  and for hypothetical years with increased variance
      for(yi in 1:n_years_samp){
        lambda_year_samp[yi,s,e,i] <- lambda(bigmatrix(make_params(species=s,
                                                                  endo_mean=(e-1),
                                                                  endo_var=(e-1),
                                                                  original = 1, # should be =1 to represent recruit
                                                                  draw=post_draws[i],
                                                                  max_size=max_size,
                                                                  rfx=T,
                                                                  samp=T,
                                                                  #year=y+1, ## this pulls the t to t1 growth-surv transition and t1 reproduction -- see make_params()
                                                                  surv_par=surv_par,
                                                                  surv_sdlg_par = surv_sdlg_par,
                                                                  grow_par=grow_par,
                                                                  grow_sdlg_par = grow_sdlg_par,
                                                                  flow_par=flow_par,
                                                                  fert_par=fert_par,
                                                                  spike_par=spike_par,
                                                                  seed_par=seed_par,
                                                                  recruit_par=recruit_par),
                                                      extension = 100)$MPMmat)
        
      }
      lambda_mean_samp[s,e,i] <- mean(lambda_year_samp[,s,e,i])
      lambda_sd_samp[s,e,i] <- sd(lambda_year_samp[,s,e,i])
      lambda_cv_samp[s,e,i] <- lambda_sd_samp[s,e,i]/lambda_mean_samp[s,e,i]

    }
    #calculating overall species means
    lambda_mean_obs[8,e,i] <- mean(lambda_mean_obs[1:7,e,i])
    lambda_sd_obs[8,e,i] <- sd(lambda_mean_obs[1:7,e,i])
    lambda_cv_obs[8,e,i] <- lambda_sd_obs[8,e,i]/lambda_mean_obs[8,e,i]
    
    lambda_mean_samp[8,e,i] <- mean(lambda_mean_samp[1:7,e,i])
    lambda_sd_samp[8,e,i] <- sd(lambda_mean_samp[1:7,e,i])
    lambda_cv_samp[8,e,i] <- lambda_sd_samp[8,e,i]/lambda_mean_samp[8,e,i]
   
  }
}

# Saving all of the simulations
# saveRDS(lambda_year_obs, file = paste0(path,"/Model_Runs/MPM_output/lambda_year_obs.rds"))
# 
# saveRDS(lambda_mean_obs, file = paste0(path,"/Model_Runs/MPM_output/lambda_mean_obs.rds"))
# saveRDS(lambda_sd_obs, file = paste0(path,"/Model_Runs/MPM_output/lambda_sd_obs.rds"))
# saveRDS(lambda_cv_obs, file = paste0(path,"/Model_Runs/MPM_output/lambda_cv_obs.rds"))
# 
# saveRDS(lambda_year_samp, file = paste0(path,"/Model_Runs/MPM_output/lambda_year_samp.rds"))
# 
# saveRDS(lambda_mean_samp, file = paste0(path,"/Model_Runs/MPM_output/lambda_mean_samp.rds"))
# saveRDS(lambda_sd_samp, file = paste0(path,"/Model_Runs/MPM_output/lambda_sd_samp.rds"))
# saveRDS(lambda_cv_samp, file = paste0(path,"/Model_Runs/MPM_output/lambda_cv_samp.rds"))
# 



# reading back in the simulations
 lambda_year_obs <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambda_year_obs.rds"))
 
 lambda_mean_obs <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambda_mean_obs.rds"))
 lambda_sd_obs <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambda_sd_obs.rds"))
 lambda_cv_obs <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambda_cv_obs.rds"))
 
 lambda_year_samp <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambda_year_samp.rds"))
 
 lambda_mean_samp <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambda_mean_samp.rds"))
 lambda_sd_samp <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambda_sd_samp.rds"))
 lambda_cv_samp <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambda_cv_samp.rds"))

 
# Turning the sampled lambdas into a dataframe
dimnames(lambda_mean_samp) <- list(species = paste0("s",1:8), Endo = paste0("e",1:2), Iteration= paste0("i",1:n_draws))
dimnames(lambda_sd_samp) <- list(species = paste0("s",1:8), Endo = paste0("e",1:2), Iteration= paste0("i",1:n_draws))
dimnames(lambda_cv_samp) <- list(species = paste0("s",1:8), Endo = paste0("e",1:2), Iteration= paste0("i",1:n_draws))

lambda_mean_samp_cube <- cubelyr::as.tbl_cube(lambda_mean_samp)
lambda_sd_samp_cube <- cubelyr::as.tbl_cube(lambda_sd_samp)
lambda_cv_samp_cube <- cubelyr::as.tbl_cube(lambda_cv_samp)

lambda_mean_samp_df <- as_tibble(lambda_mean_samp_cube) %>% 
  pivot_wider(names_from = Endo, values_from = lambda_mean_samp) %>% 
  mutate(lambda_diff = e2-e1) %>% 
  mutate(Species = case_when(species == "s1" ~ "Agrostis perennans",
                             species == "s2" ~ "Elymus villosus",
                             species == "s3" ~ "Elymus virginicus",
                             species == "s4" ~ "Festuca subverticillata",
                             species == "s5" ~ "Lolium arundinaceum",
                             species == "s6" ~ "Poa alsodes",
                             species == "s7" ~ "Poa sylvestris",
                             species == "s8" ~ "Species Mean"))  %>% 
  mutate(sampling = "samp")

lambda_sd_samp_df <- as_tibble(lambda_sd_samp_cube) %>% 
  pivot_wider(names_from = Endo, values_from = lambda_sd_samp) %>% 
  mutate(lambda_diff = e2-e1) %>% 
  mutate(Species = case_when(species == "s1" ~ "Agrostis perennans",
                             species == "s2" ~ "Elymus villosus",
                             species == "s3" ~ "Elymus virginicus",
                             species == "s4" ~ "Festuca subverticillata",
                             species == "s5" ~ "Lolium arundinaceum",
                             species == "s6" ~ "Poa alsodes",
                             species == "s7" ~ "Poa sylvestris",
                             species == "s8" ~ "Species Mean"))  %>% 
  mutate(sampling = "samp")

lambda_cv_samp_df <- as_tibble(lambda_cv_samp_cube) %>% 
  pivot_wider(names_from = Endo, values_from = lambda_cv_samp) %>% 
  mutate(lambda_diff = e2-e1) %>% 
  mutate(Species = case_when(species == "s1" ~ "Agrostis perennans",
                             species == "s2" ~ "Elymus villosus",
                             species == "s3" ~ "Elymus virginicus",
                             species == "s4" ~ "Festuca subverticillata",
                             species == "s5" ~ "Lolium arundinaceum",
                             species == "s6" ~ "Poa alsodes",
                             species == "s7" ~ "Poa sylvestris",
                             species == "s8" ~ "Species Mean"))  %>% 
  mutate(sampling = "samp")


# Turning the obs lambdas into a dataframe
dimnames(lambda_mean_obs) <- list(species = paste0("s",1:8), Endo = paste0("e",1:2), Iteration= paste0("i",1:n_draws))
dimnames(lambda_sd_obs) <- list(species = paste0("s",1:8), Endo = paste0("e",1:2), Iteration= paste0("i",1:n_draws))
dimnames(lambda_cv_obs) <- list(species = paste0("s",1:8), Endo = paste0("e",1:2), Iteration= paste0("i",1:n_draws))

lambda_mean_obs_cube <- cubelyr::as.tbl_cube(lambda_mean_obs)
lambda_sd_obs_cube <- cubelyr::as.tbl_cube(lambda_sd_obs)
lambda_cv_obs_cube <- cubelyr::as.tbl_cube(lambda_cv_obs)

lambda_mean_obs_df <- as_tibble(lambda_mean_obs_cube) %>% 
  pivot_wider(names_from = Endo, values_from = lambda_mean_obs) %>% 
  mutate(lambda_diff = e2-e1) %>% 
  mutate(Species = case_when(species == "s1" ~ "Agrostis perennans",
                             species == "s2" ~ "Elymus villosus",
                             species == "s3" ~ "Elymus virginicus",
                             species == "s4" ~ "Festuca subverticillata",
                             species == "s5" ~ "Lolium arundinaceum",
                             species == "s6" ~ "Poa alsodes",
                             species == "s7" ~ "Poa sylvestris",
                             species == "s8" ~ "Species Mean")) %>% 
  mutate(sampling = "obs")
lambda_sd_obs_df <- as_tibble(lambda_sd_obs_cube) %>% 
  pivot_wider(names_from = Endo, values_from = lambda_sd_obs) %>% 
  mutate(lambda_diff = e2-e1) %>% 
  mutate(Species = case_when(species == "s1" ~ "Agrostis perennans",
                             species == "s2" ~ "Elymus villosus",
                             species == "s3" ~ "Elymus virginicus",
                             species == "s4" ~ "Festuca subverticillata",
                             species == "s5" ~ "Lolium arundinaceum",
                             species == "s6" ~ "Poa alsodes",
                             species == "s7" ~ "Poa sylvestris",
                             species == "s8" ~ "Species Mean")) %>% 
  mutate(sampling = "obs")
lambda_cv_obs_df <- as_tibble(lambda_cv_obs_cube) %>% 
  pivot_wider(names_from = Endo, values_from = lambda_cv_obs) %>% 
  mutate(lambda_diff = e2-e1) %>% 
  mutate(Species = case_when(species == "s1" ~ "Agrostis perennans",
                             species == "s2" ~ "Elymus villosus",
                             species == "s3" ~ "Elymus virginicus",
                             species == "s4" ~ "Festuca subverticillata",
                             species == "s5" ~ "Lolium arundinaceum",
                             species == "s6" ~ "Poa alsodes",
                             species == "s7" ~ "Poa sylvestris",
                             species == "s8" ~ "Species Mean")) %>% 
  mutate(sampling = "obs")

# Combining the dataframes
lambda_mean_df <- rbind(lambda_mean_obs_df, lambda_mean_samp_df)
lambda_sd_df <- rbind(lambda_sd_obs_df, lambda_sd_samp_df)
lambda_cv_df <- rbind(lambda_cv_obs_df, lambda_cv_samp_df)


# Plots of endo effects on mean, sd and cv

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

meanlambda_plot <- ggplot(data = filter(lambda_mean_df, sampling == "obs")) +
  geom_hline(yintercept = 0, col = "black") + 
  # geom_linerange(data = lambda_mean_diff_df, aes(x = species, y = mean, ymin = 0, ymax = mean, color = species)) + 
  geom_point(aes(y = lambda_diff, x = species, color = sampling, group = sampling), position = position_jitterdodge(dodge.width = 0.75,jitter.width=0.2), alpha = .2) +
  stat_summary(aes(y = lambda_diff, x = species, fill = sampling,  group = sampling),positio = position_dodge(width= 0.75),pch = 21, col = "black", fun = mean,geom = "point", size = 2) +
  # geom_point(data = lambda_mean_diff_df, aes(y = mean, x = species, color = species), lwd = 2) +
  # scale_color_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
  # scale_fill_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
  coord_flip(ylim = c(-.5,.5)) +
  labs(y = expression(paste("Effect on ", bar(lambda))),
       x = "")+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = NA))
# meanlambda_plot

sdlambda_plot <- ggplot(data = lambda_sd_df) +
  geom_hline(yintercept = 0, col = "black") + 
  # geom_linerange(data = lambda_mean_diff_df, aes(x = species, y = mean, ymin = 0, ymax = mean, color = species)) + 
  geom_point(aes(y = lambda_diff, x = species, color = sampling, group = sampling), position = position_jitterdodge(dodge.width = 0.75,jitter.width=0.2), alpha = .2) +
  stat_summary(aes(y = lambda_diff, x = species, fill = sampling,  group = sampling),positio = position_dodge(width= 0.75),pch = 21, col = "black", fun = mean,geom = "point", size = 2) +
  # geom_point(data = lambda_mean_diff_df, aes(y = mean, x = species, color = species), lwd = 2) +
  # scale_color_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
  # scale_fill_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
  coord_flip(ylim = c(-.3,.2)) +
  labs(y = expression(paste("Effect on ", "SD(",lambda,")")),
       x = "")+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = NA))
# sdlambda_plot

cvlambda_plot <- ggplot(data = lambda_cv_df) +
  geom_hline(yintercept = 0, col = "black") + 
  # geom_linerange(data = lambda_mean_diff_df, aes(x = species, y = mean, ymin = 0, ymax = mean, color = species)) + 
  geom_point(aes(y = lambda_diff, x = species, color = sampling, group = sampling), position = position_jitterdodge(dodge.width = 0.75,jitter.width=0.2), alpha = .2) +
  stat_summary(aes(y = lambda_diff, x = species, fill = sampling,  group = sampling),positio = position_dodge(width= 0.75),pch = 21, col = "black", fun = mean,geom = "point", size = 2) +
  # geom_point(data = lambda_mean_diff_df, aes(y = mean, x = species, color = species), lwd = 2) +
  # scale_color_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
  # scale_fill_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
  coord_flip(ylim = c(-1,.8)) +
  labs(y = expression(paste("Effect on ", "CV(",lambda,")")),
       x = "")+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = NA))
# cvlambda_plot

altogether <- meanlambda_plot + sdlambda_plot + cvlambda_plot + 
  plot_layout(nrow = 1) +
  plot_annotation(title = "Observed vs Modeled Year effects")
# ggsave(altogether, filename = "endo_effects_obsVSsampled.png", width = 12, height = 6)

# I want to see the vital rate correlations
LTREB_full %>% 
  # filter(FLW_COUNT_T1<=size_t1) %>% 
  mutate(growth = size_t1-size_t,
         prop_flow = FLW_COUNT_T1) %>% 
  group_by(species, year_t1) %>% 
  # summarize(cor = cor(growth, prop_flow, use= "complete.obs")) %>% 
  summarize(mean_growth = mean(growth, na.rm = T),
            mean_flow = mean(prop_flow,na.rm = T),
            mean_surv = mean(surv_t1, na.rm = T)) %>% 
  ggplot()+
  # geom_histogram(aes(x = cor))+
  # facet_wrap(~species)
  geom_smooth(aes(y = mean_flow, x= mean_growth),col = "darkgrey", method = "glm")+#, method.args=list(family = "binomial") )+
  geom_point(aes(y = mean_flow, x = mean_growth, col = as.factor(year_t1)), alpha = .4)+
  facet_wrap(~species, scale = "free") +
  ylab("Mean Survival")+
  xlab("Mean Growth")

#some summary stats for manuscript
lambda_sd_summary <- lambda_sd_df %>% 
  filter(sampling == "obs") %>% 
  group_by(species, Species) %>% 
  summarize(mean_e1 = mean(e1),
         mean_e2 = mean(e2),
         lambda_diff = mean_e2-mean_e1,
         perc_diff = (lambda_diff/mean_e1)*100)
  

# decomposition analysis --------------------------------------------------
## stochastic simulation of lambda_S, with mean/var effects on/off in combination
## same options as above: we could randomly sample just the years we observed, 
## or simulate hypothetical years. In the latter case we need worry about correlations. 
euclidean <- function(a, b) sqrt(sum((a - b)^2))

## store lambdaS output: 8 species (7+mean),4 scenarios
n_draws<-500
lambdaS_obs<-array(NA,dim=c(4,n_spp+1,n_draws))
lambdaS_obs_extreme2<-lambdaS_obs_extreme4<-lambdaS_obs_extreme6<-array(NA,dim=c(4,n_spp+1,n_draws))
lambdaS_obs_extreme2_em<-lambdaS_obs_extreme4_em<-lambdaS_obs_extreme6_em<-array(NA,dim=c(4,n_spp+1,n_draws))

lambdaS_samp<-array(NA,dim=c(4,n_spp+1,n_draws))
lambdaS_samp_extreme10<-lambdaS_samp_extreme20<-lambdaS_samp_extreme30<-array(NA,dim=c(4,n_spp+1,n_draws))
lambdaS_samp_extreme10_em<-lambdaS_samp_extreme20_em<-lambdaS_samp_extreme30_em<-array(NA,dim=c(4,n_spp+1,n_draws))

save_lambda_obs <- array(NA,dim=c(years_obs,4,7,n_draws))
save_lambda_samp <- array(NA,dim=c(n_years_samp,4,7,n_draws))
save_lambda_obs_extreme2 <- save_lambda_obs_extreme2_em <- array(NA,dim=c(2,4,7,n_draws))
save_lambda_obs_extreme4 <- save_lambda_obs_extreme4_em <- array(NA,dim=c(4,4,7,n_draws))
save_lambda_obs_extreme6 <- save_lambda_obs_extreme6_em <- array(NA,dim=c(6,4,7,n_draws))
save_lambda_samp_extreme10 <- save_lambda_samp_extreme10_em <- array(NA,dim=c(10,4,7,n_draws))
save_lambda_samp_extreme20 <- save_lambda_samp_extreme20_em <- array(NA,dim=c(20,4,7,n_draws))
save_lambda_samp_extreme30 <- save_lambda_samp_extreme30_em <- array(NA,dim=c(30,4,7,n_draws))

for(d in 1:n_draws){
## list of transition years that we observed
A_t_obs <- A_t_samp <-list()
  for(s in 1:n_spp){
    eminus_list <- eplus_list <- eplus_mean_only_list <- eplus_var_only_list <- list()
    eminus_list_samp <- eplus_list_samp <- eplus_mean_only_list_samp <- eplus_var_only_list_samp <- list()
    # 1: sample observed years
    for(y in 1:years_obs){ ## 13 transitions matrices (2008-09 through 2020-21)
      eminus_list[[y]] <- bigmatrix(make_params(species=s,
                                                endo_mean=0,
                                                endo_var=0,
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
      eplus_list[[y]] <- bigmatrix(make_params(species=s,
                                               endo_mean=1,
                                               endo_var=1,
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
      eplus_mean_only_list[[y]] <- bigmatrix(make_params(species=s,
                                                          endo_mean=1,
                                                          endo_var=0,
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
      eplus_var_only_list[[y]] <- bigmatrix(make_params(species=s,
                                                         endo_mean=0,
                                                         endo_var=1,
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
    # 2: Sample many hypothetical years from fitted year variances
    for(yi in 1:n_years_samp){
      eminus_list_samp[[yi]] <- bigmatrix(make_params(species=s,
                                                endo_mean=0,
                                                endo_var=0,
                                                original = 1, # should be =1 to represent recruit
                                                draw=post_draws[d],
                                                max_size=max_size,
                                                rfx=T,
                                                samp=T,
                                                # year=y+1,
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
      eplus_list_samp[[yi]] <- bigmatrix(make_params(species=s,
                                               endo_mean=1,
                                               endo_var=1,
                                               original = 1, # should be =1 to represent recruit
                                               draw=post_draws[d],
                                               max_size=max_size,
                                               rfx=T,
                                               samp=T,
                                               # year=y+1,
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
      eplus_mean_only_list_samp[[yi]] <- bigmatrix(make_params(species=s,
                                                         endo_mean=1,
                                                         endo_var=0,
                                                         original = 1, # should be =1 to represent recruit
                                                         draw=post_draws[d],
                                                         max_size=max_size,
                                                         rfx=T,
                                                         samp=T,
                                                         # year=y+1,
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
      eplus_var_only_list_samp[[yi]] <- bigmatrix(make_params(species=s,
                                                        endo_mean=0,
                                                        endo_var=1,
                                                        original = 1, # should be =1 to represent recruit
                                                        draw=post_draws[d],
                                                        max_size=max_size,
                                                        rfx=T,
                                                        samp=T,
                                                        # year=y+1,
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
    }#yi loop
    ## store matrix lists in a list
    A_t_obs[[s]]<-list(eminus=eminus_list,
                   eplus_mean_only=eplus_mean_only_list,
                   eplus_var_only=eplus_var_only_list,
                   eplus=eplus_list)
    A_t_samp[[s]]<-list(eminus=eminus_list_samp,
                       eplus_mean_only=eplus_mean_only_list_samp,
                       eplus_var_only=eplus_var_only_list_samp,
                       eplus=eplus_list_samp)

    ## get lambda by year for E+ and E-
    lambda_t<-matrix(NA,2,years_obs)
    lambda_t_samp <- matrix(NA,2,n_years_samp)
    
    for(i in 1:years_obs){
      lambda_t[1,i]<-lambda(A_t_obs[[s]][[1]][[i]])
      lambda_t[2,i]<-lambda(A_t_obs[[s]][[4]][[i]])
    }
    for(i in 1:n_years_samp){
      lambda_t_samp[1,i]<-lambda(A_t_samp[[s]][[1]][[i]])
      lambda_t_samp[2,i]<-lambda(A_t_samp[[s]][[4]][[i]])
    }
    
    dist <- dist_samp <- c()
    for(i in 1:years_obs){
      dist[i] <- euclidean(c(lambda_t[1,i],lambda_t[2,i]),
                           c(mean(lambda_t[1,]),mean(lambda_t[2,])))
    }
    for(i in 1:n_years_samp){
      dist_samp[i] <- euclidean(c(lambda_t_samp[1,i],lambda_t_samp[2,i]),
                           c(mean(lambda_t_samp[1,]),mean(lambda_t_samp[2,])))
    }
    toptwo <- dist%in%rev(sort(dist))[1:2]
    topfour <- dist%in%rev(sort(dist))[1:4]
    topsix <- dist%in%rev(sort(dist))[1:6]
    top10_samp <- dist_samp%in%rev(sort(dist_samp))[1:10]
    top20_samp <- dist_samp%in%rev(sort(dist_samp))[1:20]
    top30_samp <- dist_samp%in%rev(sort(dist_samp))[1:30]
    
    ## alternative version that selects extreme years based on E- only
    toptwo_em <- lambda_t[1,]%in%sort(lambda_t[1,])[c(1,years_obs)]
    topfour_em <- lambda_t[1,]%in%sort(lambda_t[1,])[c(1:2,(years_obs-1):years_obs)]
    topsix_em <- lambda_t[1,]%in%sort(lambda_t[1,])[c(1:3,(years_obs-2):years_obs)]
    top10_samp_em <- lambda_t_samp[1,]%in%sort(lambda_t_samp[1,])[c(1:5,(n_years_samp-4):n_years_samp)]
    top20_samp_em <- lambda_t_samp[1,]%in%sort(lambda_t_samp[1,])[c(1:10,(n_years_samp-9):n_years_samp)]
    top30_samp_em <- lambda_t_samp[1,]%in%sort(lambda_t_samp[1,])[c(1:15,(n_years_samp-14):n_years_samp)]
    
    for(e in 1:4){
      lambdaS_obs[e,s,d]<-lambdaSim(mat_list = A_t_obs[[s]][[e]],max_yrs = 500)$lambdaS
      lambdaS_obs_extreme2[e,s,d]<-lambdaSim(mat_list = A_t_obs[[s]][[e]][toptwo],max_yrs = 500)$lambdaS
      lambdaS_obs_extreme4[e,s,d]<-lambdaSim(mat_list = A_t_obs[[s]][[e]][topfour],max_yrs = 500)$lambdaS
      lambdaS_obs_extreme6[e,s,d]<-lambdaSim(mat_list = A_t_obs[[s]][[e]][topsix],max_yrs = 500)$lambdaS
      lambdaS_obs_extreme2_em[e,s,d]<-lambdaSim(mat_list = A_t_obs[[s]][[e]][toptwo_em],max_yrs = 500)$lambdaS
      lambdaS_obs_extreme4_em[e,s,d]<-lambdaSim(mat_list = A_t_obs[[s]][[e]][topfour_em],max_yrs = 500)$lambdaS
      lambdaS_obs_extreme6_em[e,s,d]<-lambdaSim(mat_list = A_t_obs[[s]][[e]][topsix_em],max_yrs = 500)$lambdaS
      
      save_lambda_obs[,e,s,d] <- sapply(A_t_obs[[s]][[e]], FUN = lambda)
      save_lambda_obs_extreme2[,e,s,d] <- sapply(A_t_obs[[s]][[e]][toptwo], FUN = lambda)
      save_lambda_obs_extreme4[,e,s,d] <- sapply(A_t_obs[[s]][[e]][topfour], FUN = lambda)
      save_lambda_obs_extreme6[,e,s,d] <- sapply(A_t_obs[[s]][[e]][topsix], FUN = lambda)
      save_lambda_obs_extreme2_em[,e,s,d] <- sapply(A_t_obs[[s]][[e]][toptwo_em], FUN = lambda)
      save_lambda_obs_extreme4_em[,e,s,d] <- sapply(A_t_obs[[s]][[e]][topfour_em], FUN = lambda)
      save_lambda_obs_extreme6_em[,e,s,d] <- sapply(A_t_obs[[s]][[e]][topsix_em], FUN = lambda)
      
      lambdaS_samp[e,s,d]<-lambdaSim(mat_list = A_t_samp[[s]][[e]],max_yrs = 500)$lambdaS
      lambdaS_samp_extreme10[e,s,d]<-lambdaSim(mat_list = A_t_samp[[s]][[e]][top10_samp],max_yrs = 500)$lambdaS
      lambdaS_samp_extreme20[e,s,d]<-lambdaSim(mat_list = A_t_samp[[s]][[e]][top20_samp],max_yrs = 500)$lambdaS
      lambdaS_samp_extreme30[e,s,d]<-lambdaSim(mat_list = A_t_samp[[s]][[e]][top30_samp],max_yrs = 500)$lambdaS
      lambdaS_samp_extreme10_em[e,s,d]<-lambdaSim(mat_list = A_t_samp[[s]][[e]][top10_samp_em],max_yrs = 500)$lambdaS
      lambdaS_samp_extreme20_em[e,s,d]<-lambdaSim(mat_list = A_t_samp[[s]][[e]][top20_samp_em],max_yrs = 500)$lambdaS
      lambdaS_samp_extreme30_em[e,s,d]<-lambdaSim(mat_list = A_t_samp[[s]][[e]][top30_samp_em],max_yrs = 500)$lambdaS
      
      save_lambda_samp[,e,s,d] <- sapply(A_t_samp[[s]][[e]], FUN = lambda)
      save_lambda_samp_extreme10[,e,s,d] <- sapply(A_t_samp[[s]][[e]][top10_samp], FUN = lambda)
      save_lambda_samp_extreme20[,e,s,d] <- sapply(A_t_samp[[s]][[e]][top20_samp], FUN = lambda)
      save_lambda_samp_extreme30[,e,s,d] <- sapply(A_t_samp[[s]][[e]][top30_samp], FUN = lambda)
      save_lambda_samp_extreme10_em[,e,s,d] <- sapply(A_t_samp[[s]][[e]][top10_samp_em], FUN = lambda)
      save_lambda_samp_extreme20_em[,e,s,d] <- sapply(A_t_samp[[s]][[e]][top20_samp_em], FUN = lambda)
      save_lambda_samp_extreme30_em[,e,s,d] <- sapply(A_t_samp[[s]][[e]][top30_samp_em], FUN = lambda)
    }
    
  }#s loop
# Calculating cross species means
lambdaS_obs[1,8,d] <- mean(lambdaS_obs[1,1:7,d]) # species mean eminus
lambdaS_obs[2,8,d] <- mean(lambdaS_obs[2,1:7,d]) # species mean eplus mean only
lambdaS_obs[3,8,d] <- mean(lambdaS_obs[3,1:7,d]) # species mean eplus var only
lambdaS_obs[4,8,d] <- mean(lambdaS_obs[4,1:7,d]) # species mean eplus

lambdaS_obs_extreme2[1,8,d] <- mean(lambdaS_obs_extreme2[1,1:7,d]) # species mean eminus
lambdaS_obs_extreme2[2,8,d] <- mean(lambdaS_obs_extreme2[2,1:7,d]) # species mean eplus mean only
lambdaS_obs_extreme2[3,8,d] <- mean(lambdaS_obs_extreme2[3,1:7,d]) # species mean eplus var only
lambdaS_obs_extreme2[4,8,d] <- mean(lambdaS_obs_extreme2[4,1:7,d]) # species mean eplus
lambdaS_obs_extreme4[1,8,d] <- mean(lambdaS_obs_extreme4[1,1:7,d]) # species mean eminus
lambdaS_obs_extreme4[2,8,d] <- mean(lambdaS_obs_extreme4[2,1:7,d]) # species mean eplus mean only
lambdaS_obs_extreme4[3,8,d] <- mean(lambdaS_obs_extreme4[3,1:7,d]) # species mean eplus var only
lambdaS_obs_extreme4[4,8,d] <- mean(lambdaS_obs_extreme4[4,1:7,d]) # species mean eplus
lambdaS_obs_extreme6[1,8,d] <- mean(lambdaS_obs_extreme6[1,1:7,d]) # species mean eminus
lambdaS_obs_extreme6[2,8,d] <- mean(lambdaS_obs_extreme6[2,1:7,d]) # species mean eplus mean only
lambdaS_obs_extreme6[3,8,d] <- mean(lambdaS_obs_extreme6[3,1:7,d]) # species mean eplus var only
lambdaS_obs_extreme6[4,8,d] <- mean(lambdaS_obs_extreme6[4,1:7,d]) # species mean eplus

lambdaS_obs_extreme2_em[1,8,d] <- mean(lambdaS_obs_extreme2_em[1,1:7,d]) # species mean eminus
lambdaS_obs_extreme2_em[2,8,d] <- mean(lambdaS_obs_extreme2_em[2,1:7,d]) # species mean eplus mean only
lambdaS_obs_extreme2_em[3,8,d] <- mean(lambdaS_obs_extreme2_em[3,1:7,d]) # species mean eplus var only
lambdaS_obs_extreme2_em[4,8,d] <- mean(lambdaS_obs_extreme2_em[4,1:7,d]) # species mean eplus
lambdaS_obs_extreme4_em[1,8,d] <- mean(lambdaS_obs_extreme4_em[1,1:7,d]) # species mean eminus
lambdaS_obs_extreme4_em[2,8,d] <- mean(lambdaS_obs_extreme4_em[2,1:7,d]) # species mean eplus mean only
lambdaS_obs_extreme4_em[3,8,d] <- mean(lambdaS_obs_extreme4_em[3,1:7,d]) # species mean eplus var only
lambdaS_obs_extreme4_em[4,8,d] <- mean(lambdaS_obs_extreme4_em[4,1:7,d]) # species mean eplus
lambdaS_obs_extreme6_em[1,8,d] <- mean(lambdaS_obs_extreme6_em[1,1:7,d]) # species mean eminus
lambdaS_obs_extreme6_em[2,8,d] <- mean(lambdaS_obs_extreme6_em[2,1:7,d]) # species mean eplus mean only
lambdaS_obs_extreme6_em[3,8,d] <- mean(lambdaS_obs_extreme6_em[3,1:7,d]) # species mean eplus var only
lambdaS_obs_extreme6_em[4,8,d] <- mean(lambdaS_obs_extreme6_em[4,1:7,d]) # species mean eplus

lambdaS_samp[1,8,d] <- mean(lambdaS_samp[1,1:7,d]) # species mean eminus
lambdaS_samp[2,8,d] <- mean(lambdaS_samp[2,1:7,d]) # species mean eplus mean only
lambdaS_samp[3,8,d] <- mean(lambdaS_samp[3,1:7,d]) # species mean eplus var only
lambdaS_samp[4,8,d] <- mean(lambdaS_samp[4,1:7,d]) # species mean eplus

lambdaS_samp_extreme10[1,8,d] <- mean(lambdaS_samp_extreme10[1,1:7,d]) # species mean eminus
lambdaS_samp_extreme10[2,8,d] <- mean(lambdaS_samp_extreme10[2,1:7,d]) # species mean eplus mean only
lambdaS_samp_extreme10[3,8,d] <- mean(lambdaS_samp_extreme10[3,1:7,d]) # species mean eplus var only
lambdaS_samp_extreme10[4,8,d] <- mean(lambdaS_samp_extreme10[4,1:7,d]) # species mean eplus
lambdaS_samp_extreme20[1,8,d] <- mean(lambdaS_samp_extreme20[1,1:7,d]) # species mean eminus
lambdaS_samp_extreme20[2,8,d] <- mean(lambdaS_samp_extreme20[2,1:7,d]) # species mean eplus mean only
lambdaS_samp_extreme20[3,8,d] <- mean(lambdaS_samp_extreme20[3,1:7,d]) # species mean eplus var only
lambdaS_samp_extreme20[4,8,d] <- mean(lambdaS_samp_extreme20[4,1:7,d]) # species mean eplus
lambdaS_samp_extreme30[1,8,d] <- mean(lambdaS_samp_extreme30[1,1:7,d]) # species mean eminus
lambdaS_samp_extreme30[2,8,d] <- mean(lambdaS_samp_extreme30[2,1:7,d]) # species mean eplus mean only
lambdaS_samp_extreme30[3,8,d] <- mean(lambdaS_samp_extreme30[3,1:7,d]) # species mean eplus var only
lambdaS_samp_extreme30[4,8,d] <- mean(lambdaS_samp_extreme30[4,1:7,d]) # species mean eplus

lambdaS_samp_extreme10_em[1,8,d] <- mean(lambdaS_samp_extreme10_em[1,1:7,d]) # species mean eminus
lambdaS_samp_extreme10_em[2,8,d] <- mean(lambdaS_samp_extreme10_em[2,1:7,d]) # species mean eplus mean only
lambdaS_samp_extreme10_em[3,8,d] <- mean(lambdaS_samp_extreme10_em[3,1:7,d]) # species mean eplus var only
lambdaS_samp_extreme10_em[4,8,d] <- mean(lambdaS_samp_extreme10_em[4,1:7,d]) # species mean eplus
lambdaS_samp_extreme20_em[1,8,d] <- mean(lambdaS_samp_extreme20_em[1,1:7,d]) # species mean eminus
lambdaS_samp_extreme20_em[2,8,d] <- mean(lambdaS_samp_extreme20_em[2,1:7,d]) # species mean eplus mean only
lambdaS_samp_extreme20_em[3,8,d] <- mean(lambdaS_samp_extreme20_em[3,1:7,d]) # species mean eplus var only
lambdaS_samp_extreme20_em[4,8,d] <- mean(lambdaS_samp_extreme20_em[4,1:7,d]) # species mean eplus
lambdaS_samp_extreme30_em[1,8,d] <- mean(lambdaS_samp_extreme30_em[1,1:7,d]) # species mean eminus
lambdaS_samp_extreme30_em[2,8,d] <- mean(lambdaS_samp_extreme30_em[2,1:7,d]) # species mean eplus mean only
lambdaS_samp_extreme30_em[3,8,d] <- mean(lambdaS_samp_extreme30_em[3,1:7,d]) # species mean eplus var only
lambdaS_samp_extreme30_em[4,8,d] <- mean(lambdaS_samp_extreme30_em[4,1:7,d]) # species mean eplus
}#end d loop

# # Saving all of the simulations
saveRDS(lambdaS_obs, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_obs.rds"))
saveRDS(lambdaS_obs_extreme2, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_obs_extreme2.rds"))
saveRDS(lambdaS_obs_extreme4, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_obs_extreme4.rds"))
saveRDS(lambdaS_obs_extreme6, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_obs_extreme6.rds"))
saveRDS(lambdaS_obs_extreme2_em, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_obs_extreme2_em.rds"))
saveRDS(lambdaS_obs_extreme4_em, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_obs_extreme4_em.rds"))
saveRDS(lambdaS_obs_extreme6_em, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_obs_extreme6_em.rds"))

saveRDS(save_lambda_obs, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_obs.rds"))
saveRDS(save_lambda_obs_extreme2, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_obs_extreme2.rds"))
saveRDS(save_lambda_obs_extreme4, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_obs_extreme4.rds"))
saveRDS(save_lambda_obs_extreme6, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_obs_extreme6.rds"))
saveRDS(save_lambda_obs_extreme2_em, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_obs_extreme2_em.rds"))
saveRDS(save_lambda_obs_extreme4_em, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_obs_extreme4_em.rds"))
saveRDS(save_lambda_obs_extreme6_em, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_obs_extreme6_em.rds"))

saveRDS(lambdaS_samp, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_samp.rds"))
saveRDS(lambdaS_samp_extreme10, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_samp_extreme10.rds"))
saveRDS(lambdaS_samp_extreme20, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_samp_extreme20.rds"))
saveRDS(lambdaS_samp_extreme30, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_samp_extreme30.rds"))
saveRDS(lambdaS_samp_extreme10_em, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_samp_extreme10_em.rds"))
saveRDS(lambdaS_samp_extreme20_em, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_samp_extreme20_em.rds"))
saveRDS(lambdaS_samp_extreme30_em, file = paste0(path,"/Model_Runs/MPM_output/lambdaS_samp_extreme30_em.rds"))

saveRDS(save_lambda_samp, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_samp.rds"))
saveRDS(save_lambda_samp_extreme10, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_samp_extreme10.rds"))
saveRDS(save_lambda_samp_extreme20, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_samp_extreme20.rds"))
saveRDS(save_lambda_samp_extreme30, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_samp_extreme30.rds"))
saveRDS(save_lambda_samp_extreme10_em, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_samp_extreme10_em.rds"))
saveRDS(save_lambda_samp_extreme20_em, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_samp_extreme20_em.rds"))
saveRDS(save_lambda_samp_extreme30_em, file = paste0(path,"/Model_Runs/MPM_output/save_lambda_samp_extreme30_em.rds"))
# read in the saved outputs
lambdaS_obs <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_obs.rds"))
lambdaS_obs_extreme2 <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_obs_extreme2.rds"))
lambdaS_obs_extreme4 <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_obs_extreme4.rds"))
lambdaS_obs_extreme6 <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_obs_extreme6.rds"))
lambdaS_obs_extreme2_em <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_obs_extreme2_em.rds"))
lambdaS_obs_extreme4_em <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_obs_extreme4_em.rds"))
lambdaS_obs_extreme6_em <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_obs_extreme6_em.rds"))

save_lambda_obs <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_obs.rds"))
save_lambda_obs_extreme2 <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_obs_extreme2.rds"))
save_lambda_obs_extreme4 <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_obs_extreme4.rds"))
save_lambda_obs_extreme6 <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_obs_extreme6.rds"))
save_lambda_obs_extreme2_em <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_obs_extreme2_em.rds"))
save_lambda_obs_extreme4_em <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_obs_extreme4_em.rds"))
save_lambda_obs_extreme6_em <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_obs_extreme6_em.rds"))

lambdaS_samp <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_samp.rds"))
lambdaS_samp_extreme10 <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_samp_extreme10.rds"))
lambdaS_samp_extreme20 <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_samp_extreme20.rds"))
lambdaS_samp_extreme30 <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_samp_extreme30.rds"))
lambdaS_samp_extreme10_em <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_samp_extreme10_em.rds"))
lambdaS_samp_extreme20_em <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_samp_extreme20_em.rds"))
lambdaS_samp_extreme30_em <- read_rds(paste0(path,"/Model_Runs/MPM_output/lambdaS_samp_extreme30_em.rds"))

save_lambda_samp <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_samp.rds"))
save_lambda_samp_extreme10 <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_samp_extreme10.rds"))
save_lambda_samp_extreme20 <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_samp_extreme20.rds"))
save_lambda_samp_extreme30 <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_samp_extreme30.rds"))
save_lambda_samp_extreme10_em <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_samp_extreme10_em.rds"))
save_lambda_samp_extreme20_em <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_samp_extreme20_em.rds"))
save_lambda_samp_extreme30_em <- read_rds(paste0(path,"/Model_Runs/MPM_output/save_lambda_samp_extreme30_em.rds"))

## look at mean and variance of lambda in observed and extreme samples
cv<-function(x){sd(x)/mean(x)}

trt_order <- c("All Years", "6 Extr. Years", "4 Extr. Years", "2 Extr. Years")
mean_lambdaT_obs <- as.data.frame.table(apply(save_lambda_obs[,c(1,4),,],c(2,3,4),mean))
mean_lambdaT_obs_extreme2 <- as.data.frame.table(apply(save_lambda_obs_extreme2[,c(1,4),,],c(2,3,4),mean))
mean_lambdaT_obs_extreme4 <- as.data.frame.table(apply(save_lambda_obs_extreme4[,c(1,4),,],c(2,3,4),mean))
mean_lambdaT_obs_extreme6 <- as.data.frame.table(apply(save_lambda_obs_extreme6[,c(1,4),,],c(2,3,4),mean))

names(mean_lambdaT_obs)<-names(mean_lambdaT_obs_extreme2)<-names(mean_lambdaT_obs_extreme4)<-names(mean_lambdaT_obs_extreme6)<-c("Endo","Spp","Draw")
mean_lambdaT_obs$Trt<-"All Years"
mean_lambdaT_obs_extreme2$Trt<-"2 Extr. Years"
mean_lambdaT_obs_extreme4$Trt<-"4 Extr. Years"
mean_lambdaT_obs_extreme6$Trt<-"6 Extr. Years"
mean_lambdaT_combo<-rbind(mean_lambdaT_obs,mean_lambdaT_obs_extreme2,
                          mean_lambdaT_obs_extreme4,mean_lambdaT_obs_extreme6)
names(mean_lambdaT_combo)[4]<-"lambda"
mean_lambdaT_combo <- as_tibble(mean_lambdaT_combo) %>% 
  mutate(Species = case_when(Spp == "A" ~ "AGPE",Spp == "B" ~ "ELRI",Spp == "C" ~ "ELVI",Spp == "D" ~ "FESU",Spp == "E" ~ "LOAR",Spp == "F" ~ "POAL",Spp == "G" ~ "POSY"),
         Endo = case_when(Endo == "A" ~ "E-", Endo == "B" ~ "E+"))

ggplot(data = filter(mean_lambdaT_combo, Trt != "4 Extr. Years"))+
  geom_boxplot(aes(y=lambda,x = factor(Trt, levels = trt_order),fill = Trt))+
  facet_wrap(~Species, scales = "free_y") +
  scale_fill_manual(values = simulation_color_scheme)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text())

sd_lambdaT_obs <- as.data.frame.table(apply(save_lambda_obs[,c(1,4),,],c(2,3,4),sd))
sd_lambdaT_obs_extreme2 <- as.data.frame.table(apply(save_lambda_obs_extreme2[,c(1,4),,],c(2,3,4),sd))
sd_lambdaT_obs_extreme4 <- as.data.frame.table(apply(save_lambda_obs_extreme4[,c(1,4),,],c(2,3,4),sd))
sd_lambdaT_obs_extreme6 <- as.data.frame.table(apply(save_lambda_obs_extreme6[,c(1,4),,],c(2,3,4),sd))

names(sd_lambdaT_obs)<-names(sd_lambdaT_obs_extreme2)<-names(sd_lambdaT_obs_extreme4)<-names(sd_lambdaT_obs_extreme6)<-c("Endo","Spp","Draw")
sd_lambdaT_obs$Trt<-"All Years"
sd_lambdaT_obs_extreme2$Trt<-"2 Extr. Years"
sd_lambdaT_obs_extreme4$Trt<-"4 Extr. Years"
sd_lambdaT_obs_extreme6$Trt<-"6 Extr. Years"
sd_lambdaT_combo<-rbind(sd_lambdaT_obs,sd_lambdaT_obs_extreme2,
                          sd_lambdaT_obs_extreme4,sd_lambdaT_obs_extreme6)
names(sd_lambdaT_combo)[4]<-"sd_lambda"
sd_lambdaT_combo <- as_tibble(sd_lambdaT_combo) %>% 
  mutate(Species = case_when(Spp == "A" ~ "AGPE",Spp == "B" ~ "ELRI",Spp == "C" ~ "ELVI",Spp == "D" ~ "FESU",Spp == "E" ~ "LOAR",Spp == "F" ~ "POAL",Spp == "G" ~ "POSY"),
         Endo = case_when(Endo == "A" ~ "E-", Endo == "B" ~ "E+"))

ggplot(data = filter(sd_lambdaT_combo, Trt != "4 Extr. Years"))+
  geom_boxplot(aes(y=sd_lambda,x = factor(Trt, levels = trt_order),fill = Trt))+
  facet_wrap(~Species, scales = "free_y") +
  scale_fill_manual(values = simulation_color_scheme)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text())

## same now but for Em years
mean_lambdaT_obs_extreme2_em <- as.data.frame.table(apply(save_lambda_obs_extreme2_em[,c(1,4),,],c(2,3,4),mean))
mean_lambdaT_obs_extreme4_em <- as.data.frame.table(apply(save_lambda_obs_extreme4_em[,c(1,4),,],c(2,3,4),mean))
mean_lambdaT_obs_extreme6_em <- as.data.frame.table(apply(save_lambda_obs_extreme6_em[,c(1,4),,],c(2,3,4),mean))

names(mean_lambdaT_obs_extreme2_em)<-names(mean_lambdaT_obs_extreme4_em)<-names(mean_lambdaT_obs_extreme6_em)<-c("Endo","Spp","Draw")

mean_lambdaT_obs_extreme2_em$Trt<-"2 Extr. Years"
mean_lambdaT_obs_extreme4_em$Trt<-"4 Extr. Years"
mean_lambdaT_obs_extreme6_em$Trt<-"6 Extr. Years"
mean_lambdaT_combo_em<-rbind(mean_lambdaT_obs,mean_lambdaT_obs_extreme2_em,
                          mean_lambdaT_obs_extreme4_em,mean_lambdaT_obs_extreme6_em)
names(mean_lambdaT_combo_em)[4]<-"lambda"
mean_lambdaT_combo_em <- as_tibble(mean_lambdaT_combo_em) %>% 
  mutate(Species = case_when(Spp == "A" ~ "AGPE",Spp == "B" ~ "ELRI",Spp == "C" ~ "ELVI",Spp == "D" ~ "FESU",Spp == "E" ~ "LOAR",Spp == "F" ~ "POAL",Spp == "G" ~ "POSY"),
         Endo = case_when(Endo == "A" ~ "E-", Endo == "B" ~ "E+"))

sim_mean_boxplot <- ggplot(data = filter(mean_lambdaT_combo_em, Trt != "4 Extr. Years"))+
  geom_boxplot(aes(y=lambda,x = factor(Trt, levels = trt_order),fill = Trt))+
  facet_wrap(~Species, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = simulation_color_scheme)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y = expression(paste("mean(",lambda["t"], ")")), fill = "Scenario")
  

sd_lambdaT_obs_extreme2_em <- as.data.frame.table(apply(save_lambda_obs_extreme2_em[,c(1,4),,],c(2,3,4),sd))
sd_lambdaT_obs_extreme4_em <- as.data.frame.table(apply(save_lambda_obs_extreme4_em[,c(1,4),,],c(2,3,4),sd))
sd_lambdaT_obs_extreme6_em <- as.data.frame.table(apply(save_lambda_obs_extreme6_em[,c(1,4),,],c(2,3,4),sd))

names(sd_lambdaT_obs_extreme2_em)<-names(sd_lambdaT_obs_extreme4_em)<-names(sd_lambdaT_obs_extreme6_em)<-c("Endo","Spp","Draw")
sd_lambdaT_obs_extreme2_em$Trt<-"2 Extr. Years"
sd_lambdaT_obs_extreme4_em$Trt<-"4 Extr. Years"
sd_lambdaT_obs_extreme6_em$Trt<-"6 Extr. Years"
sd_lambdaT_combo_em<-rbind(sd_lambdaT_obs,sd_lambdaT_obs_extreme2_em,
                        sd_lambdaT_obs_extreme4_em,sd_lambdaT_obs_extreme6_em)
names(sd_lambdaT_combo_em)[4]<-"sd_lambda"
sd_lambdaT_combo_em <- as_tibble(sd_lambdaT_combo_em) %>% 
  mutate(Species = case_when(Spp == "A" ~ "AGPE",Spp == "B" ~ "ELRI",Spp == "C" ~ "ELVI",Spp == "D" ~ "FESU",Spp == "E" ~ "LOAR",Spp == "F" ~ "POAL",Spp == "G" ~ "POSY"),
         Endo = case_when(Endo == "A" ~ "E-", Endo == "B" ~ "E+"))

sim_sd_boxplot <- ggplot(data = filter(sd_lambdaT_combo_em, Trt != "4 Extr. Years"))+
  geom_boxplot(aes(y=sd_lambda,x = factor(Trt, levels = trt_order),fill = Trt))+
  facet_wrap(~Species, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = simulation_color_scheme)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = .5))+
  labs(y = expression(paste("sd(",lambda["t"], ")")), fill = "Scenario")

sim_boxplots <- sim_mean_boxplot+sim_sd_boxplot+
  plot_layout(nrow = 2, guides = "collect")+ plot_annotation(tag_levels = "A")
ggsave(sim_boxplots, filename = "sim_boxplots.png", width = 10, height = 6)
# calc some summary stats for manuscript
mean_lambdaT_em_summary <- mean_lambdaT_combo_em %>% 
  group_by(Trt) %>% 
  summarize(mean_lambda = mean(lambda))
100-mean_lambdaT_em_summary[4,2]/mean_lambdaT_em_summary[1,2]*100

sd_lambdaT_em_summary <- sd_lambdaT_combo_em %>% 
  group_by(Trt) %>% 
  summarize(mean_sd_lambda = mean(sd_lambda))
100-sd_lambdaT_em_summary[4,2]/sd_lambdaT_em_summary[1,2]*100

# calculate cross species mean and the posterior means
lambdaS_obs_diff <- lambdaS_samp_diff <- array(NA,dim=c(4,4,8,7))
for(s in 1:8){
  # eplus-eminus
  lambdaS_obs_diff[1,1,s,1] = mean(lambdaS_obs[4,s,]-lambdaS_obs[1,s,], na.rm = T) # eplus-eminus
  lambdaS_obs_diff[1,1,s,2:7] = quantile(lambdaS_obs[4,s,]-lambdaS_obs[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_obs_diff[2,1,s,1] = mean(lambdaS_obs_extreme2_em[4,s,]-lambdaS_obs_extreme2_em[1,s,], na.rm = T) # eplus-eminus for extreme2_em
  lambdaS_obs_diff[2,1,s,2:7] = quantile(lambdaS_obs_extreme2_em[4,s,]-lambdaS_obs_extreme2_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_obs_diff[3,1,s,1] = mean(lambdaS_obs_extreme4_em[4,s,]-lambdaS_obs_extreme4_em[1,s,], na.rm = T) # eplus-eminus for extreme4_em
  lambdaS_obs_diff[3,1,s,2:7] = quantile(lambdaS_obs_extreme4_em[4,s,]-lambdaS_obs_extreme4_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_obs_diff[4,1,s,1] = mean(lambdaS_obs_extreme6_em[4,s,]-lambdaS_obs_extreme6_em[1,s,], na.rm = T) # eplus-eminus for extreme6_em
  lambdaS_obs_diff[4,1,s,2:7] = quantile(lambdaS_obs_extreme6_em[4,s,]-lambdaS_obs_extreme6_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  
  lambdaS_samp_diff[1,1,s,1] = mean(lambdaS_samp[4,s,]-lambdaS_samp[1,s,], na.rm = T) # eplus-eminus
  lambdaS_samp_diff[1,1,s,2:7] = quantile(lambdaS_samp[4,s,]-lambdaS_samp[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_samp_diff[2,1,s,1] = mean(lambdaS_samp_extreme10_em[4,s,]-lambdaS_samp_extreme10_em[1,s,], na.rm = T) # eplus-eminus for extreme2_em
  lambdaS_samp_diff[2,1,s,2:7] = quantile(lambdaS_samp_extreme10_em[4,s,]-lambdaS_samp_extreme10_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_samp_diff[3,1,s,1] = mean(lambdaS_samp_extreme20_em[4,s,]-lambdaS_samp_extreme20_em[1,s,], na.rm = T) # eplus-eminus for extreme4_em
  lambdaS_samp_diff[3,1,s,2:7] = quantile(lambdaS_samp_extreme20_em[4,s,]-lambdaS_samp_extreme20_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_samp_diff[4,1,s,1] = mean(lambdaS_samp_extreme30_em[4,s,]-lambdaS_samp_extreme30_em[1,s,], na.rm = T) # eplus-eminus for extreme6_em
  lambdaS_samp_diff[4,1,s,2:7] = quantile(lambdaS_samp_extreme30_em[4,s,]-lambdaS_samp_extreme30_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  # eplus mean only - eminus
  lambdaS_obs_diff[1,2,s,1] = mean(lambdaS_obs[2,s,]-lambdaS_obs[1,s,], na.rm = T) # eplus mean only - eminus
  lambdaS_obs_diff[1,2,s,2:7] = quantile(lambdaS_obs[2,s,]-lambdaS_obs[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_obs_diff[2,2,s,1] = mean(lambdaS_obs_extreme2_em[2,s,]-lambdaS_obs_extreme2_em[1,s,], na.rm = T) # eplus mean only - eminus for extreme2_em
  lambdaS_obs_diff[2,2,s,2:7] = quantile(lambdaS_obs_extreme2_em[2,s,]-lambdaS_obs_extreme2_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_obs_diff[3,2,s,1] = mean(lambdaS_obs_extreme4_em[2,s,]-lambdaS_obs_extreme4_em[1,s,], na.rm = T) # eplus mean only - eminus for extreme4_em
  lambdaS_obs_diff[3,2,s,2:7] = quantile(lambdaS_obs_extreme4_em[2,s,]-lambdaS_obs_extreme4_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_obs_diff[4,2,s,1] = mean(lambdaS_obs_extreme6_em[2,s,]-lambdaS_obs_extreme6_em[1,s,], na.rm = T)# eplus mean only - eminus for extreme6_em
  lambdaS_obs_diff[4,2,s,2:7] = quantile(lambdaS_obs_extreme6_em[2,s,]-lambdaS_obs_extreme6_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  
  lambdaS_samp_diff[1,2,s,1] = mean(lambdaS_samp[2,s,]-lambdaS_samp[1,s,], na.rm = T) # eplus mean only - eminus
  lambdaS_samp_diff[1,2,s,2:7] = quantile(lambdaS_samp[2,s,]-lambdaS_samp[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_samp_diff[2,2,s,1] = mean(lambdaS_samp_extreme10_em[2,s,]-lambdaS_samp_extreme10_em[1,s,], na.rm = T) # eplus mean only - eminus for extreme2_em
  lambdaS_samp_diff[2,2,s,2:7] = quantile(lambdaS_samp_extreme10_em[2,s,]-lambdaS_samp_extreme10_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_samp_diff[3,2,s,1] = mean(lambdaS_samp_extreme20_em[2,s,]-lambdaS_samp_extreme20_em[1,s,], na.rm = T) # eplus mean only - eminus for extreme4_em
  lambdaS_samp_diff[3,2,s,2:7] = quantile(lambdaS_samp_extreme20_em[2,s,]-lambdaS_samp_extreme20_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_samp_diff[4,2,s,1] = mean(lambdaS_samp_extreme30_em[2,s,]-lambdaS_samp_extreme30_em[1,s,], na.rm = T)# eplus mean only - eminus for extreme6_em
  lambdaS_samp_diff[4,2,s,2:7] = quantile(lambdaS_samp_extreme30_em[2,s,]-lambdaS_samp_extreme30_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  # eplus var only - eminus
  lambdaS_obs_diff[1,3,s,1] = mean(lambdaS_obs[3,s,]-lambdaS_obs[1,s,], na.rm = T) #  eplus var only - eminus
  lambdaS_obs_diff[1,3,s,2:7] = quantile(lambdaS_obs[3,s,]-lambdaS_obs[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_obs_diff[2,3,s,1] = mean(lambdaS_obs_extreme2_em[3,s,]-lambdaS_obs_extreme2_em[1,s,], na.rm = T) #  eplus var only - eminusfor extreme2_em
  lambdaS_obs_diff[2,3,s,2:7] = quantile(lambdaS_obs_extreme2_em[3,s,]-lambdaS_obs_extreme2_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_obs_diff[3,3,s,1] = mean(lambdaS_obs_extreme4_em[3,s,]-lambdaS_obs_extreme4_em[1,s,], na.rm = T) #  eplus var only - eminus for extreme4_em
  lambdaS_obs_diff[3,3,s,2:7] = quantile(lambdaS_obs_extreme4_em[3,s,]-lambdaS_obs_extreme4_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_obs_diff[4,3,s,1] = mean(lambdaS_obs_extreme6_em[3,s,]-lambdaS_obs_extreme6_em[1,s,], na.rm = T)#  eplus var only - eminus for extreme6_em
  lambdaS_obs_diff[4,3,s,2:7] = quantile(lambdaS_obs_extreme6_em[3,s,]-lambdaS_obs_extreme6_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  
  lambdaS_samp_diff[1,3,s,1] = mean(lambdaS_samp[3,s,]-lambdaS_samp[1,s,], na.rm = T) #  eplus var only - eminus
  lambdaS_samp_diff[1,3,s,2:7] = quantile(lambdaS_samp[3,s,]-lambdaS_samp[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_samp_diff[2,3,s,1] = mean(lambdaS_samp_extreme10_em[3,s,]-lambdaS_samp_extreme10_em[1,s,], na.rm = T) #  eplus var only - eminusfor extreme2_em
  lambdaS_samp_diff[2,3,s,2:7] = quantile(lambdaS_samp_extreme10_em[3,s,]-lambdaS_samp_extreme10_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_samp_diff[3,3,s,1] = mean(lambdaS_samp_extreme20_em[3,s,]-lambdaS_samp_extreme20_em[1,s,], na.rm = T) #  eplus var only - eminus for extreme4_em
  lambdaS_samp_diff[3,3,s,2:7] = quantile(lambdaS_samp_extreme20_em[3,s,]-lambdaS_samp_extreme20_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  lambdaS_samp_diff[4,3,s,1] = mean(lambdaS_samp_extreme30_em[3,s,]-lambdaS_samp_extreme30_em[1,s,], na.rm = T)#  eplus var only - eminus for extreme6_em
  lambdaS_samp_diff[4,3,s,2:7] = quantile(lambdaS_samp_extreme30_em[3,s,]-lambdaS_samp_extreme30_em[1,s,], probs = c(0.025, 0.125,0.25,0.75,0.875,0.975), na.rm = T)
  # mean var interaction
  lambdaS_obs_diff[1,4,s,] = lambdaS_obs_diff[1,1,s,]-lambdaS_obs_diff[1,2,s,]-lambdaS_obs_diff[1,3,s,]# mean variance interaction
  lambdaS_obs_diff[2,4,s,] = lambdaS_obs_diff[2,1,s,]-lambdaS_obs_diff[2,2,s,]-lambdaS_obs_diff[2,3,s,]# mean variance interaction for extreme2_em
  lambdaS_obs_diff[3,4,s,] = lambdaS_obs_diff[3,1,s,]-lambdaS_obs_diff[3,2,s,]-lambdaS_obs_diff[3,3,s,]# mean variance interaction for extreme4_em
  lambdaS_obs_diff[4,4,s,] = lambdaS_obs_diff[4,1,s,]-lambdaS_obs_diff[4,2,s,]-lambdaS_obs_diff[4,3,s,]# mean variance interaction for extreme6_em
  
  lambdaS_samp_diff[1,4,s,] = lambdaS_samp_diff[1,1,s,]-lambdaS_samp_diff[1,2,s,]-lambdaS_samp_diff[1,3,s,]# mean variance interaction
  lambdaS_samp_diff[2,4,s,] = lambdaS_samp_diff[2,1,s,]-lambdaS_samp_diff[2,2,s,]-lambdaS_samp_diff[2,3,s,]# mean variance interaction for extreme2_em
  lambdaS_samp_diff[3,4,s,] = lambdaS_samp_diff[3,1,s,]-lambdaS_samp_diff[3,2,s,]-lambdaS_samp_diff[3,3,s,]# mean variance interaction for extreme4_em
  lambdaS_samp_diff[4,4,s,] = lambdaS_samp_diff[4,1,s,]-lambdaS_samp_diff[4,2,s,]-lambdaS_samp_diff[4,3,s,]# mean variance interaction for extreme6_em
}

# calculating the certainty of a total effect above zero
dim(lambdaS_obs)
sum(lambdaS_obs[1,8,]>0)/500*100

post_contrib <- array(NA, dim = c(4,8,500))
prob_pos <- array(NA, dim = c(8))

for(s in 1:8){
  for(i in 1:n_draws){
  post_contrib[1,s,i] <-  lambdaS_obs[4,s,i] - lambdaS_obs[1,s,i] # eplus-eminus
  post_contrib[2,s,i] <-  lambdaS_obs[2,s,i] - lambdaS_obs[1,s,i] # eplus mean only -eminus
  post_contrib[3,s,i] <-  lambdaS_obs[3,s,i] - lambdaS_obs[1,s,i] # eplus var only-eminus
  post_contrib[4,s,i] <-  post_contrib[1,s,i] - post_contrib[2,s,i] - post_contrib[3,s,i] # mean-variance interaction
  
  prob_pos[s] <- sum(post_contrib[1,s,]>0)/500*100
  
}
}
################################################################
##### Plot of stochastic lambda contributions
################################################################
dimnames(lambdaS_obs_diff) <- list(Scenario = c("Ambient variability","Most variability (2 best and worst years)","4 Most Extreme Years","More variability (6 best and worst years)"), Contribution = c("Full Effect","Mean only","Variance only","Interaction"), Species = c(species_list, "Species Mean"), Quantile = c("mean","fifth","twelvepointfive","twentyfifth","seventyfifth","eightysevenpointfive","ninetyfifth"))
lambdaS_obs_diff_cube <- cubelyr::as.tbl_cube(lambdaS_obs_diff)
lambdaS_obs_diff_df <- as_tibble(lambdaS_obs_diff_cube) %>% 
  pivot_wider(names_from = "Quantile", values_from = lambdaS_obs_diff) %>% 
  mutate(Sampling = "observed")

dimnames(lambdaS_samp_diff) <- list(Scenario = c("All Years","2 Most Extreme Years","4 Most Extreme Years","6 Most Extreme Years"), Contribution = c("Full Effect","Mean only","Variance only","Interaction"), Species = c(species_list, "Species Mean"), Quantile = c("mean","fifth","twelvepointfive","twentyfifth","seventyfifth","eightysevenpointfive","ninetyfifth"))
lambdaS_samp_diff_cube <- cubelyr::as.tbl_cube(lambdaS_samp_diff)
lambdaS_samp_diff_df <- as_tibble(lambdaS_samp_diff_cube) %>% 
  pivot_wider(names_from = "Quantile", values_from = lambdaS_samp_diff) %>% 
  mutate(Sampling = "sampled")

lambdaS_diff_df <- lambdaS_obs_diff_df %>% 
  rbind(lambdaS_samp_diff_df)

x_levels <- c( "Interaction", "Variance only", "Mean only", "Full Effect")
contributions_obs_plot <- ggplot(data = filter(lambdaS_obs_diff_df, Scenario != "4 Most Extreme Years")) +
  geom_hline(yintercept = 0, col = "black") +
  geom_linerange(aes(x = Contribution, ymin = fifth, ymax = ninetyfifth, group = Scenario), color = "grey38",position = position_dodge(width = .7), lwd = .8) +
  geom_linerange(aes(x = Contribution, ymin = twelvepointfive, ymax = eightysevenpointfive, group = Scenario),color = "grey38", position = position_dodge(width = .7), lwd = 1.5) +
  geom_linerange(aes(x = Contribution, ymin = twentyfifth, ymax = seventyfifth, group = Scenario),color = "grey38", position = position_dodge(width = .7), lwd = 2.4) +
  
  geom_linerange(aes(x = Contribution, ymin = fifth, ymax = ninetyfifth, color = Scenario),position = position_dodge(width = .7), lwd = .4) +
  geom_linerange(aes(x = Contribution, ymin = twelvepointfive, ymax = eightysevenpointfive, color = Scenario),position = position_dodge(width = .7), lwd = 1) +
  geom_linerange(aes(x = Contribution, ymin = twentyfifth, ymax = seventyfifth, color = Scenario),position = position_dodge(width = .7), lwd = 2) +
  
  geom_point(aes(y = mean, x = Contribution, group = Scenario, pch = Contribution), color = "grey38", position = position_dodge(width = .7), size = 3.9) +
  geom_point(aes(y = mean, x = Contribution, color = Scenario, pch = Contribution), position = position_dodge(width = .7), size = 3.5) +
  facet_wrap(~Species, nrow = 2, scales = "free") + coord_flip() +
  scale_shape_manual(values = c(16,18,15,17))+
  scale_fill_manual(values = c(simulation_color_scheme), labels = ~ stringr::str_wrap(.x, width = 25))+
  scale_color_manual(values = c(simulation_color_scheme), labels = ~ stringr::str_wrap(.x, width = 25))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  scale_x_discrete(limits = x_levels)+
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
contributions_obs_plot
ggsave(contributions_obs_plot, filename = "contributions_obs_plot.png", width = 13, height = 8)


contributions_samp_plot <- ggplot(data = lambdaS_samp_diff_df) +
  geom_hline(yintercept = 0, col = "black") +
  geom_linerange(aes(x = Contribution, ymin = twentyfifth, ymax = seventyfifth, color = Scenario),position = position_dodge(width = .6), lwd = 2) +
  geom_linerange(aes(x = Contribution, ymin = twelvepointfive, ymax = eightysevenpointfive, color = Scenario),position = position_dodge(width = .6), lwd = 1) +
  geom_linerange(aes(x = Contribution, ymin = fifth, ymax = ninetyfifth, color = Scenario),position = position_dodge(width = .6)) +
  geom_point(aes(y = mean, x = Contribution, fill = Scenario, pch = Contribution), position = position_dodge(width = .6), lwd = 3) +
  facet_wrap(~Species, nrow = 2, scales = "free") + coord_flip() +
  scale_shape_manual(values = c(21, 22, 23,24))+
  scale_fill_manual(values = c(simulation_color_scheme))+
  scale_color_manual(values = c(simulation_color_scheme))+
  scale_x_discrete(limits = x_levels)+
  labs(title = "Sampled Years", x = "", y = expression(paste("Endophyte effect on", " ", lambda["s"])))+
  theme(panel.background = element_blank(),
        plot.background = element_rect(color = "white"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(size = .5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = rel(1.5), face = "bold"),
        axis.text.x = element_text(size = rel(1.5), face = "bold"),
        axis.title = element_text(size = rel(1.5)),
        strip.background = element_blank())
# contributions_samp_plot
# ggsave(contributions_samp_plot, filename = "contributions_samp_plot.tiff", width = 15, height = 6)

contributions_plot <- contributions_obs_plot+contributions_samp_plot + plot_layout(ncol = 1)
ggsave(contributions_plot, filename = "contributions_plot.tiff", width = 15, height = 20)


# just the species mean obs contribution plot

mean_contributions_obs_plot <- ggplot(data = filter(lambdaS_obs_diff_df, Species == "Species Mean" & Scenario != "4 Most Extreme Years")) +
  geom_hline(yintercept = 0, col = "black") +
  geom_linerange(aes(x = Contribution, ymin = fifth, ymax = ninetyfifth, group = Scenario), color = "grey38",position = position_dodge(width = .7), lwd = .8) +
  geom_linerange(aes(x = Contribution, ymin = twelvepointfive, ymax = eightysevenpointfive, group = Scenario),color = "grey38", position = position_dodge(width = .7), lwd = 1.5) +
  geom_linerange(aes(x = Contribution, ymin = twentyfifth, ymax = seventyfifth, group = Scenario),color = "grey38", position = position_dodge(width = .7), lwd = 2.4) +
  
  geom_linerange(aes(x = Contribution, ymin = fifth, ymax = ninetyfifth, color = Scenario),position = position_dodge(width = .7), lwd = .4) +
  geom_linerange(aes(x = Contribution, ymin = twelvepointfive, ymax = eightysevenpointfive, color = Scenario),position = position_dodge(width = .7), lwd = 1) +
  geom_linerange(aes(x = Contribution, ymin = twentyfifth, ymax = seventyfifth, color = Scenario),position = position_dodge(width = .7), lwd = 2) +
  
  geom_point(aes(y = mean, x = Contribution, group = Scenario, pch = Contribution), color = "grey38", position = position_dodge(width = .7), size = 3.9) +
  geom_point(aes(y = mean, x = Contribution, color = Scenario, pch = Contribution), position = position_dodge(width = .7), size = 3.5) +
  facet_wrap(~Species, nrow = 2, scales = "free") + coord_flip() +
  scale_shape_manual(values = c(16,18,15,17))+
  scale_fill_manual(values = c(simulation_color_scheme), labels = ~ stringr::str_wrap(.x, width = 25))+
  scale_color_manual(values = c(simulation_color_scheme), labels = ~ stringr::str_wrap(.x, width = 25))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  scale_x_discrete(limits = x_levels)+
  labs( x = "", y = expression(paste("Symbiosis effect on", " ", lambda["s"])))+
  guides(pch = "none")+
  theme(panel.background = element_blank(),
        plot.background = element_rect(color = "white"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(size = .5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = rel(1.2), face = "bold"),
        axis.text.x = element_text(size = rel(1), face = "bold"),
        axis.title = element_text(size = rel(1.5)),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.title=element_text(size=rel(.7)),
        legend.text=element_text(size=rel(.6)),
        legend.key.size = unit(.7,"cm"))
mean_contributions_obs_plot
ggsave(mean_contributions_obs_plot, filename = "mean_contributions_obs_plot.png", width = 5, height =4)



# some calculations for manuscript
meanspecies_calcs <- lambdaS_obs_diff_df %>% 
  filter(Scenario == "Ambient variability" | Scenario == "Most variability (2 best and worst years)") %>%
  select(Species,Scenario,Contribution,Sampling,mean) %>% 
  pivot_wider(names_from = Contribution, values_from = mean) %>% 
  mutate(fullminusintx = `Full Effect` - `Interaction`,
         percent_var_offullminusintx = (`Variance only`/fullminusintx)*100,
         percent_mean_offullminusintx = (`Mean only`/fullminusintx)*100,
         percent_var = (`Variance only`/`Full Effect`)*100,
         percent_mean = (`Mean only`/`Full Effect`)*100,
         percent_intx = (`Interaction`/`Full Effect`)*100,
         percent_varofmean= (`Variance only`/`Mean only`)*100,
         meandivbyvar = (`Mean only`/`Variance only`))

meanspecies_change <- meanspecies_calcs %>% 
  dplyr::select(Species,Scenario, `Full Effect`, `Mean only`, `Variance only`, Interaction) %>% 
  group_by(Species) %>% 
  summarize(perc_change_full = (max(`Full Effect`)- min(`Full Effect`))/min(`Full Effect`)*100,
            perc_change_var = (max(`Variance only`)- min(`Variance only`))/min(`Variance only`)*100)




## one random draw
par(mfrow=c(1,2))
barplot(apply(lambdaS_mat,c(1,2),mean,na.rm=T),beside=T,main="normal")
barplot(apply(lambdaS_mat_extreme,c(1,2),mean,na.rm=T),beside=T,main="extreme")


normal<-apply(lambdaS_mat,c(1,2),mean,na.rm=T)
extreme<-apply(lambdaS_mat_extreme,c(1,2),mean,na.rm=T)

barplot(rbind(normal[c(1,4),]),beside=T)
barplot(rbind(extreme[c(1,4),]),beside=T)

tot.endo<-normal[4,]-normal[1,]
mean.endo<-normal[2,]-normal[1,]
var.endo<-normal[3,]-normal[1,]
intx.endo<-tot.endo-(mean.endo+var.endo)


(tot.endo)
colSums(rbind(mean.endo,var.endo,intx.endo))

tot.endo.ex<-extreme[4,]-extreme[1,]
mean.endo.ex<-extreme[2,]-extreme[1,]
var.endo.ex<-extreme[3,]-extreme[1,]
intx.endo.ex<-tot.endo.ex-(mean.endo.ex+var.endo.ex)

barplot(rbind(tot.endo,mean.endo,var.endo,intx.endo),beside=T,main="normal")
barplot(rbind(tot.endo.ex,mean.endo.ex,var.endo.ex,intx.endo.ex),beside=T,main="extreme")






## if we wanted to skew the sampling to worst and best years
lambdaS_mat_extreme<-array(NA,dim=c(4,7,n_draws))
## get E+ and E- lambdas by year 

## E- lambda_t
s<-2
lambda_t<-matrix(NA,2,13)
for(i in 1:13){
  lambda_t[1,i]<-lambda(A_t_obs[[s]][[1]][[i]])
  lambda_t[2,i]<-lambda(A_t_obs[[s]][[4]][[i]])
}

plot(lambda_t[1,],lambda_t[2,])

dist<-c()
for(i in 1:13){
  dist[i] <- euclidean(c(lambda_t[1,i],lambda_t[2,i]),
                       c(mean(lambda_t[1,]),mean(lambda_t[2,])))
}
topsix <- dist%in%rev(sort(dist))[1:6]
points(lambda_t[1,topsix],lambda_t[2,topsix],col="blue",pch=16)

topsix_em <- lambda_t[1,]%in%sort(lambda_t[1,])[c(1:3,(years_obs-2):years_obs)]
points(lambda_t[1,topsix_em],lambda_t[2,topsix_em],pch="X")

x<-lambda_t[1,]
sort(x)
rank(x)
order(x)
which.max(lambda_t[1,])
sort(lambda_t[1,])
rank(lambda_t[1,])

order(runif(5))

euc.cent <- c(mean(lambda_t[1,]),mean(lambda_t[2,]))
points(x=euc.cent[1],y=euc.cent[2],pch=16,col="red")
#euclidian center
euclidean <- function(a, b) sqrt(sum((a - b)^2))
dist<-c()
for(i in 1:13){
  dist[i] <- euclidean(c(lambda_t[1,i],lambda_t[2,i]),
          c(mean(lambda_t[1,]),mean(lambda_t[2,])))
}
topsix <- dist%in%rev(sort(dist))[1:6]
# getting ranks of E+ and E- lambdas to try to get a top three and bottom three separately
recentered <-lambda_t-c(mean(lambda_t[1,]),mean(lambda_t[2,])) #recentering the points around (0,0)
positive <- (atan2(x=recentered[1,], y = recentered[2,])) # calculating the angle in radians between the points and the x axis where negattive values are less than the mean. 


points(lambda_t[1,topsix],lambda_t[2,topsix],col="blue",pch=16)
points(lambda_t[1,topsix & positive>0], lambda_t[2,topsix & positive>0], col = "green", pch = 16)
#euclidian center of the extreme points
euc.cent_extreme <- c(c(mean(lambda_t[1,topsix]),mean(lambda_t[2,topsix])))
points(x=euc.cent_extreme[1],y=euc.cent_extreme[2],pch=16,col="orange")


for(e in 1:4){
  lambdaS_mat_extreme[e,s,d]<-lambdaSim(mat_list = A_t_obs[[s]][[e]][topfour],max_yrs = 500)$lambdaS
}

# getting ranks of E+ and E- lambdas to try to get a top three and bottom three separately
recentered <-lambda_t-c(mean(lambda_t[1,]),mean(lambda_t[2,])) #recentering the points around (0,0)

positive <- (atan2(x=recentered[1,], y = recentered[2,])) # calculating the angle in radians between the points and the x axis where negattive values are less than the mean. 


topsix
euclidean(c(-1,-1),c(0,0))

origin_dist <- c()
for(i in 1:13){
  origin_dist[i] <- euclidean(c(recentered[1,i],recentered[2,i]),
                              c(mean(recentered[1,]),mean(recentered[2,])))
}
  


# trying to figure out how to sample from the full posterior just 5 and 95% quantiles
qnorm(p = .95, mean = 0, sd = 1)

plot(sample(x = dnorm(x = -1:1, mean = 0, sd = 2), size = 100, replace = T))
plot(sample(x = dnorm(x = qnorm(p = .95, mean = 0, sd = 1):10, mean = 0, sd = 2),size = 100, replace = T), col = "red")
# Could do this in the vital rates, but maybe better to just increase the sd by 10 percent
rfx_surv <- sample(qnorm(p=c(seq(.9,1,.001),seq(0,.1,.001)),mean = 0, sd=surv_par$sigma_year[draw,species,(endo_var+1)]),size = 1) # sample the tenth and ninetieth percentiles

