## Purpose: compile and save matrix population model from vital rate estimates for stochastic grass endophyte population model
### Compiling these to share matrices with Robin Snyder
## Authors: Joshua Fowler and Tom Miller
#############################################################

library(tidyverse)
library(scales)
library(bayesplot)
library(popbio)
# library(countreg)
library(actuar)
library(rstan)



quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

#############################################################################################
####### Read in Data and creating size bins------------------
#############################################################################################

source("Analyses/endodemog_data_processing.R")

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
tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"
joshpath <- "~/Dropbox/EndodemogData/"
path<-joshpath

surv_fit_seedling <- read_rds(paste0(path,"/Model_Runs/endo_seedling_surv.rds"))
surv_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_surv_woseedling.rds"))
grow_fit_seedling <- read_rds(paste0(path,"/Model_Runs/endo_seedling_grow_PIG_10000iterations.rds"))
grow_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_grow_PIG.rds"))
flw_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_flw.rds"))
fert_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_fert_PIG.rds"))
spike_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_spike_year_plot_nb.rds"))
seedmean_fit <- read_rds(paste0(path,"/Model_Runs/seed_mean.rds"))
stos_fit <- read_rds(paste0(path,"/Model_Runs/endo_spp_s_to_s.rds")) 


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

for(e in 1:2){
  for(s in 1:7){
    for(y in 1:13){
        for(i in 1:length(post_draws)){
      name <- paste0("iter", i)
          iter.list[[name]] <- bigmatrix(make_params(species=s,
                                                               endo_mean=(e-1),
                                                               endo_var=(e-1),
                                                               original = 1, # should be =1 to represent recruit
                                                               draw=post_draws[i],
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
                                                   extension = 100)$MPMmat # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
      name <- paste0("y",year_vec[y])
      year.list[[name]] <- iter.list
    }
    name <- paste0(spp_vec[s])
    spp.list[[name]] <- year.list
  }
  name <- paste0(endo_vec[e])
  endo.list[[name]] <- spp.list
  }
}
saveRDS(endo.list, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/GrassEndo_list_of_matrices.rds")

saveRDS(endo.list, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/AGPE_oneiter_list_of_matrices.rds")
GrassEndo_list_of_matrices <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/GrassEndo_list_of_matrices.rds")



for(e in 1:2){
  for(s in 1:7){
    for(y in 1:13){
      for(i in 1:length(post_draws)){
        name <- paste0("iter", i)
        fert_iter.list[[name]] <- bigmatrix(make_params(species=s,
                                                        endo_mean=(e-1),
                                                        endo_var=(e-1),
                                                        original = 1, # should be =1 to represent recruit
                                                        draw=post_draws[i],
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


saveRDS(fert_endo.list, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/GrassEndo_list_of_matrices_fert.rds")

GrassEndo_list_of_matrices_fert <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/GrassEndo_list_of_matrices_fert.rds")

# GrassEndo_list_of_matrices_fert <- fert_endo.list



AGPE_oneiter_list_of_matrices <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/AGPE_oneiter_list_of_matrices.rds")

#pulling out just agrostis perennans
AGPE_list_of_matrices <- list()

AGPE_list_of_matrices[["Eminus"]] <- GrassEndo_list_of_matrices$Eminus$Agrostis_perennans
AGPE_list_of_matrices[["Eplus"]] <- GrassEndo_list_of_matrices$Eplus$Agrostis_perennans


saveRDS(AGPE_list_of_matrices, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/AGPE_list_of_matrices.rds")
 


# Making mean matrices for each year/endo/species

GrassEndo_mean_list_of_matrices <- list()

GrassEndo_mean_list_of_matrices_fert <- list()

for(e in 1:2){
  for(s in 1:7){
    for(y in 1:13){
      mean_matrix <- popbio::mean.list(GrassEndo_list_of_matrices[[e]][[s]][[y]])
      name <- paste0("y",year_vec[y])
      year.list[[name]]  <-mean_matrix 
    }
    name <- paste0(spp_vec[s])
    spp.list[[name]]  <- year.list
  }
  name <- paste0(endo_vec[e])
  endo.list[[name]] <- spp.list
}

saveRDS(endo.list, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/GrassEndo_mean_list_of_matrices.rds")
GrassEndo_mean_list_of_matrices <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/GrassEndo_mean_list_of_matrices.rds")




for(e in 1:2){
  for(s in 1:7){
    for(y in 1:13){
      mean_matrix <- popbio::mean.list(GrassEndo_list_of_matrices_fert[[e]][[s]][[y]])
      name <- paste0("y",year_vec[y])
      year.list[[name]]  <-mean_matrix 
    }
    name <- paste0(spp_vec[s])
    spp.list[[name]]  <- year.list
  }
  name <- paste0(endo_vec[e])
  endo.list[[name]] <- spp.list
}

saveRDS(endo.list, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/GrassEndo_mean_list_of_matrices_fert.rds")
GrassEndo_mean_list_of_matrices_fert <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/GrassEndo_mean_list_of_matrices_fert.rds")




View(endo.list)


image(GrassEndo_list_of_matrices$Eminus$Agrostis_perennans$y2009$iter3)
image(GrassEndo_mean_list_of_matrices$Eminus$Agrostis_perennans$y2009)


image(GrassEndo_mean_list_of_matrices$Eminus$Elymus_villosus$y2010)

