## Purpose: compile and save matrix population model from vital rate estimates for stochastic grass endophyte population model
### Compiling these to share matrices with Robin Snyder and to make a life cycle graph
## Authors: Joshua Fowler and Tom Miller
#############################################################

library(tidyverse)
library(scales)
library(bayesplot)
library(popbio)
# library(countreg)
library(actuar)
library(rstan)
library(Rage)
library(scales)


quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

#############################################################################################
####### Read in Data and creating size bins------------------
#############################################################################################
tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"
joshpath <- "~/Dropbox/EndodemogData/"
path<-joshpath
#source("Analyses/endodemog_data_processing.R")
#instead of running processing script, read in LTREB_full, which that script creates
LTREB_full <- read.csv(paste0(path,"Fulldataplusmetadata/LTREB_full.csv"))

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


# surv_fit_seedling <- read_rds(paste0(path,"/Model_Runs/endo_seedling_surv.rds"))
# surv_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_surv_woseedling.rds"))
# grow_fit_seedling <- read_rds(paste0(path,"/Model_Runs/endo_seedling_grow_PIG_10000iterations.rds"))
# grow_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_grow_PIG.rds"))
# flw_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_flw.rds"))
# fert_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_fert_PIG.rds"))
# spike_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_spike_year_plot_nb.rds"))
# seedmean_fit <- read_rds(paste0(path,"/Model_Runs/seed_mean.rds"))
# stos_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_s_to_s.rds")) 

surv_fit_seedling <- read_rds(paste0(path,"/Model_Runs/endo_seedling_surv.rds"))
surv_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_surv_woseedling_quadXorigin.rds"))
grow_fit_seedling <- read_rds(paste0(path,"/Model_Runs/endo_seedling_grow_PIG_10000iterations.rds"))
grow_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_grow_PIG_quadXorigin.rds"))
flw_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_flw_quadXorigin.rds"))
fert_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_fert_PIG_quadXorigin.rds"))
spike_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_spike_year_plot_nb_quadXorigin.rds"))
seedmean_fit <- read_rds(paste0(path,"/Model_Runs/seed_mean.rds"))
stos_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_s_to_s.rds")) 

surv_par <- rstan::extract(surv_fit, pars =quote_bare(beta0,betasize,betasize_2,betaendo,
                                                      tau_year, tau_plot))
surv_sdlg_par <- rstan::extract(surv_fit_seedling, pars =quote_bare(beta0,betaendo,
                                                                    tau_year, tau_plot))
grow_par <- rstan::extract(grow_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo,
                                                       tau_year, tau_plot,
                                                       sigma))
grow_sdlg_par <- rstan::extract(grow_fit_seedling, pars = quote_bare(beta0,betaendo,
                                                                     tau_year, tau_plot,
                                                                     sigma))
flow_par <- rstan::extract(flw_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo,
                                                      tau_year, tau_plot))
fert_par <- rstan::extract(fert_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo,
                                                       tau_year, tau_plot, 
                                                       theta, # individual-level deviates
                                                       sigma)) # Inverse Gaussian Shape parameter))
spike_par <- rstan::extract(spike_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo,
                                                         tau_year, tau_plot,
                                                         phi))
seed_par <- rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
recruit_par <- rstan::extract(stos_fit, pars = quote_bare(beta0,betaendo,
                                                          tau_year, tau_plot))
# dim(surv_par$tau_year)


#############################################################################################
####### Run the MPM ------------------
#############################################################################################

# make the list of parameters and calculate mean lambdas
n_draws <- 500# the means are the same whether we do 500 or 1000 draws
post_draws <- sample.int(7500,size=n_draws) # The models except for seedling growth have 7500 iterations. That one has more (15000 iterations) to help it converge.

# Running and saving the annual matrices for each species

spp_vec = c("Agrostis_perennans",
             "Elymus_villosus",
             "Elymus_virginicus",
             "Festuca_subverticillata",
             "Lolium_arundinaceum",
             "Poa_alsodes",
             "Poa_sylvestris")
year_vec <- c(2009:2021) # 13 years because of reproduction measured in year t1; needs to be before growth, so no year 1
endo_vec <- c("Eminus", "Eplus")

iter.list <- list()
year.list <- list()
spp.list <- list()
endo.list <- list()

fert_iter.list <- list()
fert_year.list <- list()
fert_spp.list <- list()
fert_endo.list <- list()
# 
# for(e in 1:2){
#   for(s in 1:7){
#     for(y in 1:13){
#         for(i in 1:length(post_draws)){
#       name <- paste0("iter", i)
#           iter.list[[name]] <- bigmatrix(make_params_quadXorigin(species=s,
#                                                                endo_mean=(e-1),
#                                                                endo_var=(e-1),
#                                                                original = 1, # should be =1 to represent recruit
#                                                                draw=post_draws[i],
#                                                                max_size=max_size,
#                                                                rfx=T,
#                                                                year=y+1,
#                                                                surv_par=surv_par,
#                                                                surv_sdlg_par = surv_sdlg_par,
#                                                                grow_par=grow_par,
#                                                                grow_sdlg_par = grow_sdlg_par,
#                                                                flow_par=flow_par,
#                                                                fert_par=fert_par,
#                                                                spike_par=spike_par,
#                                                                seed_par=seed_par,
#                                                                recruit_par=recruit_par),
#                                                    quadratic = 1,
#                                                    extension = 100)$MPMmat # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
#       name <- paste0("y",year_vec[y])
#       year.list[[name]] <- iter.list
#     }
#     name <- paste0(spp_vec[s])
#     spp.list[[name]] <- year.list
#   }
#   name <- paste0(endo_vec[e])
#   endo.list[[name]] <- spp.list
#   }
# }
# saveRDS(endo.list, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/GrassEndo_list_of_matrices.rds")
# 
# saveRDS(endo.list, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/AGPE_oneiter_list_of_matrices.rds")
# GrassEndo_list_of_matrices <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/GrassEndo_list_of_matrices.rds")
# 
# 
# 
# 
# for(e in 1:2){
#   for(s in 1:7){
#     for(y in 1:13){
#       for(i in 1:length(post_draws)){
#         name <- paste0("iter", i)
#         fert_iter.list[[name]] <- bigmatrix(make_params(species=s,
#                                                         endo_mean=(e-1),
#                                                         endo_var=(e-1),
#                                                         original = 1, # should be =1 to represent recruit
#                                                         draw=post_draws[i],
#                                                         max_size=max_size,
#                                                         rfx=T,
#                                                         year=y+1,
#                                                         surv_par=surv_par,
#                                                         surv_sdlg_par = surv_sdlg_par,
#                                                         grow_par=grow_par,
#                                                         grow_sdlg_par = grow_sdlg_par,
#                                                         flow_par=flow_par,
#                                                         fert_par=fert_par,
#                                                         spike_par=spike_par,
#                                                         seed_par=seed_par,
#                                                         recruit_par=recruit_par),
#                                             extension = 100)$Fmat # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
#       }
#       name <- paste0("y",year_vec[y])
#       fert_year.list[[name]] <- fert_iter.list
#       
#     }
#     name <- paste0(spp_vec[s])
#     fert_spp.list[[name]] <- fert_year.list
#     
#   }
#   name <- paste0(endo_vec[e])
#   fert_endo.list[[name]] <- fert_spp.list
#   
# }
# 
# 
# saveRDS(fert_endo.list, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/GrassEndo_list_of_matrices_fert.rds")
# 
# GrassEndo_list_of_matrices_fert <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/GrassEndo_list_of_matrices_fert.rds")
# 
# # GrassEndo_list_of_matrices_fert <- fert_endo.list
# 


## get posterior mean parameter esimates -- rather than take the mean of the matrices
## if one had to design the worst way to do this, it would be as follows!

#1-surv_par
surv_par_means<-list()
surv_par_means$beta0 <- apply(surv_par$beta0,c(2,3),mean)
surv_par_means$betasize <- apply(surv_par$betasize,c(2,3),mean)
surv_par_means$betasize_2 <- apply(surv_par$betasize_2,c(2,3),mean)
surv_par_means$betaendo <- apply(surv_par$betaendo,c(2),mean)
surv_par_means$tau_year <- apply(surv_par$tau_year,c(2,3,4),mean)
surv_par_means$tau_plot <- apply(surv_par$tau_plot,c(2),mean)

#2-surv_sdlg_par
surv_sdlg_par_means<-list()
surv_sdlg_par_means$beta0<-apply(surv_sdlg_par$beta0,2,mean)
surv_sdlg_par_means$betaendo<-apply(surv_sdlg_par$betaendo,2,mean)
surv_sdlg_par_means$tau_year <- apply(surv_sdlg_par$tau_year,c(2,3,4),mean)
surv_sdlg_par_means$tau_plot <- apply(surv_sdlg_par$tau_plot,c(2),mean)

#3-grow_par
grow_par_means<-list()
grow_par_means$beta0 <- apply(grow_par$beta0,c(2,3),mean)
grow_par_means$betasize <- apply(grow_par$betasize,c(2,3),mean)
grow_par_means$betasize_2 <- apply(grow_par$betasize_2,c(2,3),mean)
grow_par_means$betaendo <- apply(grow_par$betaendo,c(2),mean)
grow_par_means$tau_year <- apply(grow_par$tau_year,c(2,3,4),mean)
grow_par_means$tau_plot <- apply(grow_par$tau_plot,c(2),mean)
grow_par_means$sigma <- mean(grow_par$sigma)

#4-grow_sdlg_par
grow_sdlg_par_means<-list()
grow_sdlg_par_means$beta0<-apply(grow_sdlg_par$beta0,2,mean)
grow_sdlg_par_means$betaendo<-apply(grow_sdlg_par$betaendo,2,mean)
grow_sdlg_par_means$tau_year <- apply(grow_sdlg_par$tau_year,c(2,3,4),mean)
grow_sdlg_par_means$tau_plot <- apply(grow_sdlg_par$tau_plot,c(2),mean)
grow_sdlg_par_means$sigma <- mean(grow_sdlg_par$sigma)


#5-flow_par
flow_par_means<-list()
flow_par_means$beta0 <- apply(flow_par$beta0,c(2,3),mean)
flow_par_means$betasize <- apply(flow_par$betasize,c(2,3),mean)
flow_par_means$betasize_2 <- apply(flow_par$betasize_2,c(2,3),mean)
flow_par_means$betaendo <- apply(flow_par$betaendo,c(2),mean)
flow_par_means$tau_year <- apply(flow_par$tau_year,c(2,3,4),mean)
flow_par_means$tau_plot <- apply(flow_par$tau_plot,c(2),mean)

#6-fert_par
fert_par_means<-list()
fert_par_means$beta0 <- apply(fert_par$beta0,c(2,3),mean)
fert_par_means$betasize <- apply(fert_par$betasize,c(2,3),mean)
fert_par_means$betasize_2 <- apply(fert_par$betasize_2,c(2,3),mean)
fert_par_means$betaendo <- apply(fert_par$betaendo,c(2),mean)
fert_par_means$tau_year <- apply(fert_par$tau_year,c(2,3,4),mean)
fert_par_means$tau_plot <- apply(fert_par$tau_plot,c(2),mean)
fert_par_means$sigma <- mean(fert_par$sigma)

#7-spike_par
spike_par_means<-list()
spike_par_means$beta0 <- apply(spike_par$beta0,c(2,3),mean)
spike_par_means$betasize <- apply(spike_par$betasize,c(2,3),mean)
spike_par_means$betasize_2 <- apply(spike_par$betasize_2,c(2,3),mean)
spike_par_means$betaendo <- apply(spike_par$betaendo,c(2),mean)
spike_par_means$tau_year <- apply(spike_par$tau_year,c(2,3,4),mean)
spike_par_means$tau_plot <- apply(spike_par$tau_plot,c(2),mean)

#8-seed_par
seed_par_means<-list()
seed_par_means$beta0 <- apply(seed_par$beta0,c(2),mean)
seed_par_means$betaendo <- apply(seed_par$betaendo,c(2),mean)

#9-recruit_par
recruit_par_means<-list()
recruit_par_means$beta0 <- apply(recruit_par$beta0,c(2),mean)
recruit_par_means$betaendo <- apply(recruit_par$betaendo,c(2),mean)
recruit_par_means$tau_year <- apply(recruit_par$tau_year,c(2,3,4),mean)
recruit_par_means$tau_plot <- apply(recruit_par$tau_plot,c(2),mean)



##### rewriting the fertility function to return the number of inflorescencees, not the number of recruits#####
fx<-function(x, params,spei=0, quadratic){
  xb<-pmin(x,params$max_size) # any predicted plants larger than max size will set to be max size
  flw <- invlogit(params$flow_int + params$flow_spei*spei + params$flow_slope*log(xb) + params$flow_slope_2*(log(xb)^2)*quadratic)
  fert_mean <- exp(params$fert_int + params$fert_spei*spei + params$fert_slope*log(xb) + params$fert_slope_2*(log(xb)^2)*quadratic)
  
  fert <- fert_mean/(1-exp(   (-1*params$fert_sigma)* (sqrt(1 + ((2*fert_mean)/params$fert_sigma)   )   -   1)  ))
  # spike <- exp(params$spike_int + params$spike_spei*spei + params$spike_slope*log(xb))
  # seeds_per_spike <- params$seeds_per_spike
  # recruits_per_seed <- invlogit(params$recruits_per_seed + params$recruits_spei*spei)
  seedlings <- flw * fert# * spike * seeds_per_spike * recruits_per_seed
  return(seedlings)
}


#simulation to check expected value of ZT-PIG
# xb<-pmin(x,params$max_size) 
# fert_mean = exp(params$fert_int + params$fert_spei*spei + params$fert_slope*log(xb) + params$fert_slope_2*(log(xb)^2)*quadratic)
# fert_pois <- fert_pig <- fert_ztpig <- c()
# for(i in 1:length(fert_mean)){
# fert_pois[i] <- mean(rpois(fert_mean[i], n = 1000000))
# fert_pig[i]<- mean(rpoisinvgauss(mean=fert_mean[i],shape=(fert_mean[i]*params$fert_sigma),n = 1000000))
# fert_ztpig[i] <- mean(rpoisinvgauss(mean = fert[i], shape = (fert[i]*params$fert_sigma), n = 1000000))
# }
# 
lambda <- 2; sigma <-  8

# confirming the math with simulation for poisson and zero-truncated poisson
mean(rpois(n = 100000, lambda = lambda))
mean(rztpois(n = 100000, lambda = lambda))
2/(1-exp(-2))


# now trying this for the zero truncated PIG
theta <- mean(rinvgauss(n = 1000, mean = 1, shape = sigma))


  
  prob <- dpois(1:1000, lambda=lambda * theta)
  prob_trunc <- (1- dpois(0, lambda=lambda * theta))

fert_sim <- mean(sample(x = 1:1000, size = 1000000, replace = T, prob = prob/prob_trunc))
fert_sim

lambda_trunc <- lambda/(1-exp(   (-1*sigma)* (sqrt(1 + ((2*lambda)/sigma)   )   -   1)  ))

lambda_trunc


# hist(rpoisinvgauss(n = 100000, mean = lambda_trunc, shape = lambda_trunc*sigma))
# theta <- rinvgauss(n = 1000, mean = 1, sigma)
# 
# lambda_sim <- dpois(1:1000, lambda=lambda * theta)
# lambda_0 <- (1- dpois(0, lambda=lambda * theta))
# for(i in 1:1000){
# fert_sim[i] <- sample(x = 1:1000, size = 1, replace = T, prob = lambda_sim/lambda_0)
# }


# another way using the dpoisinvgauss function
# fert_pig <- dpoisinvgauss(x=x,mean=(fert_mean),shape=((fert_mean)*params$fert_sigma))
# fert_Zero <- dpoisinvgauss(x=0,mean=(fert_mean),shape=((fert_mean)*params$fert_sigma))
# 
# fert_ztpig <- fert_pig/(1-fert_Zero)


## re-write make_params to use posterior mean params

make_params_quadXorigin <- function(species,endo_mean,endo_var,LTRE=F, endo_mean_U, endo_mean_F, endo_var_U, endo_var_F, original=0,rfx=F,spei=F,year=NULL,repro_offset = 1, max_size,samp=F,samp_extreme=0,
                                    surv_par_means,surv_sdlg_par_means,grow_par_means,grow_sdlg_par_means,flow_par_means,fert_par_means,spike_par_means,seed_par_means,recruit_par_means){
  if(LTRE==F){
    endo_mean_U <- endo_mean_F <- endo_mean
    endo_var_U <- endo_var_F <- endo_var
  }
  if(LTRE==T){
    endo_mean_U <- endo_mean_U
    endo_var_U <- endo_var_U
    endo_mean_F <- endo_mean_F
    endo_var_F <- endo_var_F
  }
  
  if(rfx==F){rfx_surv <- rfx_surv_sdlg <- rfx_grow <- rfx_grow_sdlg <- rfx_flow <- rfx_fert <- rfx_spike <- rfx_rct <-  0}
  
  if(rfx==T & samp==F){
    ## timing and survival and growth (size_t / y_t1) is meant to line up with reproduction (size_t1 / y_t1)
    rfx_surv <- surv_par_means$tau_year[species,(endo_var_U+1),(year)]; 
    rfx_surv_sdlg <-surv_sdlg_par_means$tau_year[species,(endo_var_U+1),(year)];
    rfx_grow <- grow_par_means$tau_year[species,(endo_var_U+1),(year)];
    rfx_grow_sdlg <- grow_sdlg_par_means$tau_year[species,(endo_var_U+1),(year)];
    rfx_flow <- flow_par_means$tau_year[species,(endo_var_F+1),(year-repro_offset)]; # fitting 
    rfx_fert <- fert_par_means$tau_year[species,(endo_var_F+1),(year-repro_offset)]; 
    rfx_spike <- spike_par_means$tau_year[species,(endo_var_F+1),(year-repro_offset)];
    rfx_rct <- recruit_par_means$tau_year[species,(endo_var_F+1),year];
  }
  
  if(spei == F){spei_surv <- spei_surv_sdlg <- spei_grow <- spei_grow_sdlg <- spei_flow <- spei_fert <- spei_spike <- spei_rct <-  0}

  
  
  params <- c()
  #survival
  params$surv_int <- surv_par_means$beta0[species,original+1] + 
    endo_mean_U * surv_par_means$betaendo[species] + 
    rfx_surv
  params$surv_spei <- spei_surv
  params$surv_slope <- surv_par_means$betasize[species,original+1]
  params$surv_slope_2 <- surv_par_means$betasize_2[species,original+1]
  
  # seedling survival
  params$surv_sdlg_int <- surv_sdlg_par_means$beta0[species] + 
    endo_mean_U * surv_sdlg_par_means$betaendo[species] + 
    rfx_surv_sdlg
  params$surv_sdlg_spei <- spei_surv_sdlg
  
  #growth
  params$grow_int <- grow_par_means$beta0[species,original+1] + 
    endo_mean_U * grow_par_means$betaendo[species] + 
    rfx_grow
  params$grow_spei <- spei_grow
  params$grow_slope <- grow_par_means$betasize[species,original+1] 
  params$grow_slope_2 <- grow_par_means$betasize_2[species,original+1]  
  
  params$grow_sigma <- grow_par_means$sigma
  # seedling growth
  params$grow_sdlg_int <- grow_sdlg_par$beta0[species] + 
    endo_mean_U * grow_sdlg_par$betaendo[species] + 
    rfx_grow_sdlg
  params$grow_sdlg_spei <- spei_grow_sdlg
  params$grow_sdlg_sigma <- grow_sdlg_par$sigma
  
  #flowering
  params$flow_int <- flow_par_means$beta0[species,original+1] + 
    endo_mean_F * flow_par_means$betaendo[species] + 
    + rfx_flow
  params$flow_spei <- spei_flow
  params$flow_slope <- flow_par_means$betasize[species,original+1]  
  params$flow_slope_2 <- flow_par_means$betasize_2[species,original+1]  
  
  #fertility
  params$fert_int <- fert_par_means$beta0[species,original+1] +
    endo_mean_F * fert_par_means$betaendo[species] +
    rfx_fert
  params$fert_spei <- spei_fert
  params$fert_slope <- fert_par_means$betasize[species,original+1]
  params$fert_slope_2 <- fert_par_means$betasize_2[species,original+1]
  params$fert_sigma <- fert_par_means$sigma
  
  
  #spikelets
  params$spike_int <- spike_par_means$beta0[species,original+1]  +
    endo_mean_F * spike_par_means$betaendo[species] +
    rfx_spike
  params$spike_spei <- spei_spike
  params$spike_slope <- spike_par_means$betasize[species,original+1]  
  params$spike_slope_2 <- spike_par_means$betasize_2[species,original+1]  
  
  
  #seeds per spikelet
  params$seeds_per_spike <- seed_par_means$beta0[species] + 
    endo_mean_F * seed_par_means$betaendo[species]
  #recruits per seed
  params$recruits_per_seed <- recruit_par_means$beta0[species] + 
    endo_mean_F * recruit_par_means$betaendo[species] +
    rfx_rct
  params$recruits_spei <- spei_rct
  
  #tack on max size
  params$max_size <- max_size$max_size[species]
  
  return(params)
}


## this code is a little convoluted, but Josh and I confirmed that y iteration 1
## generates the 2008-2009 transition year, and the reproductive rates can be compared 
## to the 2008 data. y iteration 13 generates the 2020-2021 transition matrix and can
## be compared wth 2020 flowering data. 

for(e in 1:2){
  for(s in 1:7){
    for(y in 1:13){
      for(i in 1){
        name <- paste0("iter", i)
        fert_iter.list[[name]] <- bigmatrix(make_params_quadXorigin(species=s,
                                                        endo_mean=(e-1),
                                                        endo_var=(e-1),
                                                        original = 1, # should be =1 to represent recruit
                                                        max_size=max_size,
                                                        rfx=T,
                                                        year=y+1,
                                                        surv_par_means=surv_par_means,
                                                        surv_sdlg_par_means = surv_sdlg_par_means,
                                                        grow_par_means=grow_par_means,
                                                        grow_sdlg_par_means = grow_sdlg_par_means,
                                                        flow_par_means=flow_par_means,
                                                        fert_par_means=fert_par_means,
                                                        spike_par_means=spike_par_means,
                                                        seed_par_means=seed_par_means,
                                                        recruit_par_means=recruit_par_means),
                                            quadratic = 1,
                                            extension = 100)$Fmat # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
      }
      name <- paste0("y",year_vec[y])
      fert_year.list[[name]] <- fert_iter.list
      
    }
    name <- paste0(spp_vec[s])
    fert_spp.list[[name]] <- fert_year.list
    
  }
  name <- paste0(endo_vec[e])
  fert_endo.list[[name]] <- fert_spp.list
  
}



## compare predicted fertility (# infs) and observed -- by species, year, and endo
## prep the data - pull out recruits only, drop 2007--2008, and drop the birth year (since these are in the recruit bin)
LTREB_full %>% 
  filter(origin_01 == 1) %>% 
  filter(year_t!=2007) %>% 
  filter(year_t!=birth) -> LTREB_recruits

## for every row, pull the expected value from the corresponding matrix
LTREB_recruits$expected_infs<-NA
for(i in 1:nrow(LTREB_recruits)){
  row_dat<-LTREB_recruits[i,]
  ## skip if anything I need is NA
  if(is.na(row_dat$endo_01) | is.na(row_dat$species_index) | is.na(row_dat$year_t_index) | is.na(row_dat$size_t)){next}
  LTREB_recruits$expected_infs[i] <- fert_endo.list[[(row_dat$endo_01+1)]][[row_dat$species_index]][[(row_dat$year_t_index-1)]]$iter1[1,(row_dat$size_t+1)]
}

## compare predicted and observed
pdf(paste0(path,"/Model_Runs/MPM_output/predicted_observed_inflorescences.pdf"))
plot(jitter(LTREB_recruits$expected_infs),
     jitter(LTREB_recruits$FLW_COUNT_T),col=alpha("black",0.25),
     xlab="Predicted inflorescneces",ylab="Observed inflorescences")
abline(0,1,col="red",lwd=3)
dev.off()

## compare across years
LTREB_recruits %>% 
  group_by(year_t,species_index) %>% 
  summarise(obs_mean_inf = mean(FLW_COUNT_T,na.rm=T),
            pred_mean_inf = mean(expected_infs,na.rm=T)) -> infs_by_year

pdf(paste0(path,"/Model_Runs/MPM_output/inflorescences_by_year.pdf"),width=4,height=7)
par(mfrow=c(4,2),mar=c(5,4,1,1))
for(i in 1:7){
plot(infs_by_year$year_t[infs_by_year$species_index==i],infs_by_year$obs_mean_inf[infs_by_year$species_index==i],
     type="n",xlab="Year",ylab="Inflorescences")
lines(infs_by_year$year_t[infs_by_year$species_index==i],infs_by_year$pred_mean_inf[infs_by_year$species_index==i])
points(infs_by_year$year_t[infs_by_year$species_index==i],infs_by_year$obs_mean_inf[infs_by_year$species_index==i])
title(spp_vec[i])
}
dev.off()

## write out matrices
saveRDS(fert_endo.list, file = paste0(path,"/Model_Runs/MPM_output/GrassEndo_list_of_mean_matrices.rds"))







# Now we can also get the growth/surv matrix
T_iter.list <- list()
T_year.list <- list()
T_spp.list <- list()
T_endo.list <- list()

for(e in 1:2){
  for(s in 1:7){
    for(y in 1:13){
      for(i in 1){
        name <- paste0("iter", i)
        T_iter.list[[name]] <- bigmatrix(make_params_quadXorigin(species=s,
                                                                 endo_mean=(e-1),
                                                                 endo_var=(e-1),
                                                                 original = 1, # should be =1 to represent recruit
                                                                 max_size=max_size,
                                                                 rfx=T,
                                                                 year=y+1,
                                                                 surv_par_means=surv_par_means,
                                                                 surv_sdlg_par_means = surv_sdlg_par_means,
                                                                 grow_par_means=grow_par_means,
                                                                 grow_sdlg_par_means = grow_sdlg_par_means,
                                                                 flow_par_means=flow_par_means,
                                                                 fert_par_means=fert_par_means,
                                                                 spike_par_means=spike_par_means,
                                                                 seed_par_means=seed_par_means,
                                                                 recruit_par_means=recruit_par_means),
                                         quadratic = 1,
                                         extension = 100)$Tmat # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
      }
      name <- paste0("y",year_vec[y])
      T_year.list[[name]] <- T_iter.list
      
    }
    name <- paste0(spp_vec[s])
    T_spp.list[[name]] <- T_year.list
    
  }
  name <- paste0(endo_vec[e])
  T_endo.list[[name]] <- T_spp.list
  
}



## write out matrices
saveRDS(T_endo.list, file = paste0(path,"/Model_Runs/MPM_output/GrassEndo_mean_list_of_Pmatrices.rds"))

