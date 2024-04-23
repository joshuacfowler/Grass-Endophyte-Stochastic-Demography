## Title: Grass endophyte population model with a bayesian framework
## Purpose: Builds matrix population model from vital rate estimates with SPEI effect
## Authors: Joshua and Tom
#############################################################


library(tidyverse)
library(scales)
library(bayesplot)
library(popbio)
# library(countreg)
library(actuar)
library(rstan)
library(patchwork) # for putting plots together


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

# source("Analyses/endodemog_data_processing.R")
LTREB_full <- read_csv(paste0(joshpath,"Fulldataplusmetadata/LTREB_full.csv"))

max_size <- LTREB_full %>% 
  dplyr::select(species,species_index, size_t) %>% 
  filter(!is.na(size_t)) %>% 
  group_by(species, species_index) %>% 
  summarise(actual_max_size = max(size_t),
            max_size = quantile(size_t,probs=0.975),
            max_size_99 = quantile(size_t,probs=0.99)) # The mean and sd effects plots look basically identical with either max size
# We use these min max values of SPEI for the relationship between growth rates and climate
spei12_minmax <- LTREB_full %>% 
  dplyr::select(species,species_index, spei12) %>% 
  filter(!is.na(spei12)) %>% 
  group_by(species, species_index) %>% 
  summarise(min_spei = min(spei12),
            max_spei = max(spei12))
spei3_minmax <- LTREB_full %>% 
  dplyr::select(species,species_index, spei3) %>% 
  filter(!is.na(spei3)) %>% 
  group_by(species, species_index) %>% 
  summarise(min_spei = min(spei3),
            max_spei = max(spei3))
  
#############################################################################################
####### Read in matrix population functions ------------------
#############################################################################################

source("Analyses/MPM_functions.R")

#############################################################################################
####### Read in Stan vital rate model outputs ------------------
#############################################################################################

# spei12_surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_spei12.rds")
# spei3_surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_spei3.rds")
# 
# spei12_surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_spei12_endo_spp_surv_woseedling_linear.rds")
# spei3_surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_spei3_endo_spp_surv_woseedling_linear.rds")
# 
# spei12_grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_seedling_grow_spei12_10000iterations.rds")
# spei3_grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_seedling_grow_spei3_10000iterations.rds")
# 
# spei12_grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_spei12.rds")
# spei3_grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_spei3.rds")
# 
# spei12_flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw_spei12.rds")
# spei3_flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw_spei3.rds")
# 
# spei12_fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_spei12.rds")
# spei3_fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_spei3.rds")
# 
# spei12_spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_spei12.rds")
# spei3_spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_spei3.rds")
# 
# seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_seed_mean.rds") # doesn't include spei predictor
# 
# spei12_stos_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_s_to_s_spei12.rds") 
# spei3_stos_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_s_to_s_spei3.rds") 
# 
# 

# versions of models include linear and quadratic size effects and more flexible origin/recruit intercepts
spei12_surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_spei12.rds")
spei3_surv_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_seedling_surv_spei3.rds")

spei12_surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_spei12_endo_spp_surv_woseedling_linear_quadXorigin.rds")
spei3_surv_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_spei3_endo_spp_surv_woseedling_linear_quadXorigin.rds")

spei12_grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_seedling_grow_spei12_10000iterations.rds")
spei3_grow_fit_seedling <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_seedling_grow_spei3_10000iterations.rds")

spei12_grow_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_spei12_quadXorigin.rds")
spei3_grow_fit <- readRDS("~/Dropbox/EndodemogData/Model_Runs/climate_endo_grow_PIG_spei3_quadXorigin.rds")

spei12_flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw_spei12_quadXorigin.rds")
spei3_flw_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_flw_spei3_quadXorigin.rds")

spei12_fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_spei12_quadXorigin.rds")
spei3_fert_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_fert_PIG_spei3_quadXorigin.rds")

spei12_spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_spei12_quadXorigin.rds")
spei3_spike_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_spike_spei3_quadXorigin.rds")

seedmean_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/endo_spp_seed_mean.rds") # doesn't include spei predictor

spei12_stos_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_s_to_s_spei12.rds") 
spei3_stos_fit <- read_rds("~/Dropbox/EndodemogData/Model_Runs/climate_endo_spp_s_to_s_spei3.rds") 


# Pulling out the actual parameters
#survival
spei12_surv_par <- rstan::extract(spei12_surv_fit, pars =quote_bare(beta0,betasize,betasize_2,betaendo,betaspei_endo,tau_year, tau_plot, sigma_year))
spei3_surv_par <- rstan::extract(spei3_surv_fit, pars =quote_bare(beta0,betasize,betasize_2,betaendo,betaspei_endo,tau_year, tau_plot, sigma_year))
# seedling survival
spei12_surv_sdlg_par <- rstan::extract(spei12_surv_fit_seedling, pars =quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))
spei3_surv_sdlg_par <- rstan::extract(spei3_surv_fit_seedling, pars =quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))

# growth
spei12_grow_par <- rstan::extract(spei12_grow_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma, sigma_year))
spei3_grow_par <- rstan::extract(spei3_grow_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma, sigma_year))

# seedling growth
spei12_grow_sdlg_par <- rstan::extract(spei12_grow_fit_seedling, pars = quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma, sigma_year))
spei3_grow_sdlg_par <- rstan::extract(spei3_grow_fit_seedling, pars = quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma, sigma_year))

#flowering
spei12_flow_par <- rstan::extract(spei12_flw_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))
spei3_flow_par <- rstan::extract(spei3_flw_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))

# number of flower tillers
spei12_fert_par <- rstan::extract(spei12_fert_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))
spei3_fert_par <- rstan::extract(spei3_fert_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))

# spikelets/infl
spei12_spike_par <- rstan::extract(spei12_spike_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))
spei3_spike_par <- rstan::extract(spei3_spike_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))

# seed mean/spikelet
seed_par <- rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
# seedling recruitment
spei12_recruit_par <- rstan::extract(spei12_stos_fit, pars = quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))
spei3_recruit_par <- rstan::extract(spei3_stos_fit, pars = quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))


# endophyte effects on lambda mean across spei values ---------------------------


## SET-UP
# make the list of parameters and calculate mean lambdas
n_draws <- 500# the means are the same whether we do 500 or 1000 draws
post_draws <- sample.int(7500,size=n_draws) # The models except for seedling growth have 7500 iterations. That one has more (15000 iterations) to help it converge.
n_spp <- length(unique(LTREB_full$species))
n_endo <- 2
# stepsize of spei values
spei_steps <- 10
spei12_range <- spei3_range <- array(dim = c(n_spp,spei_steps))
for(s in 1:n_spp){
spei12_range[s,] <- seq(from = spei12_minmax$min_spei[s], to = spei12_minmax$max_spei[s], length.out = 10)
spei3_range[s,] <- seq(from = spei3_minmax$min_spei[s], to = spei3_minmax$max_spei[s], length.out = 10)
}
lambda_spei12 <- lambda_spei3 <- array(dim = c(spei_steps,(n_spp+1),n_endo,n_draws))


for(i in 1:n_draws){
  for(e in 1:n_endo){
    for(s in 1:n_spp){
      for(c in 1:spei_steps){
        # 1. fit MPM with 12 Month SPEI effect, without random effects
        lambda_spei12[c,s,e,i] <- lambda(bigmatrix(make_params_quadXorigin(species=s,
                                                                 endo_mean=(e-1),
                                                                 endo_var=(e-1),
                                                                 original = 1, # should be =1 to represent recruit
                                                                 draw=post_draws[i],
                                                                 max_size=max_size,
                                                                 rfx=F,
                                                                 spei = T,
                                                                 surv_par=spei12_surv_par,
                                                                 surv_sdlg_par = spei12_surv_sdlg_par,
                                                                 grow_par=spei12_grow_par,
                                                                 grow_sdlg_par = spei12_grow_sdlg_par,
                                                                 flow_par=spei12_flow_par,
                                                                 fert_par=spei12_fert_par,
                                                                 spike_par=spei12_spike_par,
                                                                 seed_par=seed_par,
                                                                 recruit_par=spei12_recruit_par),
                                                     spei = spei12_range[s,c],
                                                     quadratic = 1,
                                                     extension = 100)$MPMmat) # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
        # 1. fit MPM with 12 Month SPEI effect, without random effects
        lambda_spei3[c,s,e,i] <- lambda(bigmatrix(make_params_quadXorigin(species=s,
                                                               endo_mean=(e-1),
                                                               endo_var=(e-1),
                                                               original = 1, # should be =1 to represent recruit
                                                               draw=post_draws[i],
                                                               max_size=max_size,
                                                               rfx=F,
                                                               spei = T,
                                                               surv_par=spei3_surv_par,
                                                               surv_sdlg_par = spei3_surv_sdlg_par,
                                                               grow_par=spei3_grow_par,
                                                               grow_sdlg_par = spei3_grow_sdlg_par,
                                                               flow_par=spei3_flow_par,
                                                               fert_par=spei3_fert_par,
                                                               spike_par=spei3_spike_par,
                                                               seed_par=seed_par,
                                                               recruit_par=spei3_recruit_par),
                                                   spei = spei3_range[s,c],
                                                   quadratic = 1,
                                                   extension = 100)$MPMmat) # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
        
      } # endo of spei loop
    } # endo of species loop
  }# endo of endo loop
} # end of iteration loop

# calculating overall species mean
for(i in 1:n_draws){
  for(e in 1:n_endo){
    lambda_spei12[,8,e,i] <- rowMeans(lambda_spei12[,1:7,e,i])
    lambda_spei3[,8,e,i] <- rowMeans(lambda_spei3[,1:7,e,i])
    
  }
}

saveRDS(lambda_spei12, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_spei12.rds")
saveRDS(lambda_spei3, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_spei3.rds")

lambda_spei12 <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_spei12.rds")
lambda_spei3 <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_spei3.rds")

dimnames(lambda_spei12) <-dimnames(lambda_spei3) <- list(spei = paste0("spei",1:spei_steps),species = paste0("s",1:(n_spp+1)), Endo = paste0("e",1:2), Iteration= paste0("i",1:n_draws))
lambda_spei12_cube <- cubelyr::as.tbl_cube(lambda_spei12)
lambda_spei3_cube <- cubelyr::as.tbl_cube(lambda_spei3)

lambda_spei12_df <- as_tibble(lambda_spei12_cube) %>% 
  mutate(Species = case_when(species == "s1" ~ "Agrostis perennans",
                             species == "s2" ~ "Elymus villosus",
                             species == "s3" ~ "Elymus virginicus",
                             species == "s4" ~ "Festuca subverticillata",
                             species == "s5" ~ "Lolium arundinaceum",
                             species == "s6" ~ "Poa alsodes",
                             species == "s7" ~ "Poa sylvestris",
                             species == "s8" ~ "Species Mean"),
         Endo = case_when(Endo == "e1" ~ "S-",
                          Endo == "e2" ~ "S+")) 
lambda_spei3_df <- as_tibble(lambda_spei3_cube) %>% 
  mutate(Species = case_when(species == "s1" ~ "Agrostis perennans",
                             species == "s2" ~ "Elymus villosus",
                             species == "s3" ~ "Elymus virginicus",
                             species == "s4" ~ "Festuca subverticillata",
                             species == "s5" ~ "Lolium arundinaceum",
                             species == "s6" ~ "Poa alsodes",
                             species == "s7" ~ "Poa sylvestris",
                             species == "s8" ~ "Species Mean"),
         Endo = case_when(Endo == "e1" ~ "S-",
                          Endo == "e2" ~ "S+")) 

dimnames(spei12_range) <- dimnames(spei3_range) <- list(species = paste0("s",1:(n_spp)), spei = paste0("spei",1:spei_steps))
spei12_range_cube <- cubelyr::as.tbl_cube(spei12_range)
spei3_range_cube <- cubelyr::as.tbl_cube(spei3_range)

spei12_range_df <- as_tibble(spei12_range_cube) %>% 
  rename(spei_value = spei12_range)
spei3_range_df <- as_tibble(spei3_range_cube) %>% 
  rename(spei_value = spei3_range)

lambda_spei12_df <- lambda_spei12_df %>% 
  left_join(spei12_range_df) %>% 
  filter(species != "s8")

lambda_spei3_df <- lambda_spei3_df %>% 
  left_join(spei3_range_df) %>% 
  filter(species != "s8")

lambda_spei12_mean <- lambda_spei12_df %>% 
  left_join(spei12_range_df) %>% 
  group_by(spei_value,species,Endo,Species) %>% 
  summarize(lambda_spei_mean = mean(lambda_spei12)) %>% 
  left_join(spei12_range_df) 

lambda_spei3_mean <- lambda_spei3_df %>% 
  left_join(spei3_range_df) %>% 
  group_by(spei_value,species,Endo,Species) %>% 
  summarize(lambda_spei_mean = mean(lambda_spei3)) %>% 
  left_join(spei3_range_df) 


# reading in lambda_mean and lambda_var with 500 post draws from dropbox, derived from MPM_analysis script
lambda_mean <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_mean_quadXorigin.rds")

lambda_hold <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_hold_quadXorigin.rds")
lambda_var <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_var_quadXorigin.rds")

# Mean endophyte difference and quantiles
lambda_means <- matrix(NA,8,2)
lambda_mean_diff <- matrix(NA,8,7)
for(s in 1:8){
  lambda_means[s,1] <- mean(lambda_mean[s,1,])
  lambda_means[s,2] <- mean(lambda_mean[s,2,])
  lambda_mean_diff[s,1] = mean(lambda_mean[s,2,] - lambda_mean[s,1,])
  lambda_mean_diff[s,2:7] = quantile(lambda_mean[s,2,] - lambda_mean[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
}

# Calculating endophyte effect on sd and variance
lambda_sds <- matrix(NA,8,2)
lambda_vars <- matrix(NA,8,2)

lambda_cv <- (lambda_var)/(lambda_mean) # this is the "variance penalty" term not the true CV right now

lambda_cvs <- matrix(NA,8,2)

lambda_sd_diff <- matrix(NA,8,7)
lambda_var_diff <- matrix(NA,8,7)
lambda_cv_diff <-  matrix(NA,8,7)
for(s in 1:8){
  lambda_sds[s,1] <- mean(lambda_var[s,1,])
  lambda_sds[s,2] <- mean(lambda_var[s,2,])
  
  lambda_vars[s,1] <- mean(lambda_var[s,1,])^2
  lambda_vars[s,2] <- mean(lambda_var[s,2,])^2
  # 
  lambda_cvs[s,1] <- mean(lambda_cv[s,1,])
  lambda_cvs[s,2] <- mean(lambda_cv[s,2,])
  
  lambda_sd_diff[s,1] = mean(lambda_var[s,2,]) - mean(lambda_var[s,1,])
  lambda_sd_diff[s,2:7] = quantile(lambda_var[s,2,] - lambda_var[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
  
  lambda_var_diff[s,1] = mean(lambda_var[s,2,]^2 - lambda_var[s,1,]^2)
  lambda_var_diff[s,2:7] = quantile(lambda_var[s,2,]^2 - lambda_var[s,1,]^2,probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
  
  lambda_cv_diff[s,1] = mean(lambda_cv[s,2,] - lambda_cv[s,1,])
  lambda_cv_diff[s,2:7] = quantile(lambda_cv[s,2,] - lambda_cv[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
  
}
##### Making a plot
# Set color scheme based on analine blue
endophyte_color_scheme <- c("#fdedd3","#f3c8a8", "#5a727b", "#4986c7", "#181914",  "#163381")
color_scheme_set(endophyte_color_scheme)
# color_scheme_view()
# And creating a color palette for each year
# yearcount = length(unique(LTREB_full$year_t))
# yearcolors<- colorRampPalette(brewer.pal(8,"Dark2"))(yearcount)
# scales::show_col(yearcolors)
species_list <- c("A. perennans", "E. villosus", "E. virginicus", "F. subverticillata", "L. arundinaceae", "P. alsodes", "P. sylvestris")
species_code_list <- c("AGPE", "ELRI", "ELVI", "FESU", "LOAR", "POAL", "POSY")


spei12_lambda_plot <- ggplot(data = lambda_spei12_df)+
  geom_path(aes(x = spei_value, y = lambda_spei12, group = interaction(Endo,Iteration), color = Endo), alpha = .05, lwd = .4)+
  geom_line(data = lambda_spei12_mean, aes(x = spei_value, y = lambda_spei_mean, group = Endo, color = Endo),lwd = 1)+
  facet_wrap(~Species, scales = "free", nrow = 2)+
  scale_color_manual(values = endophyte_color_scheme[c(2,6)])+
  theme_classic()+ theme(strip.background = element_blank(),
                         strip.text = element_text(face = "italic")) + 
  labs(x = "12-month SPEI", y = expression("Pop. Growth Rate " (lambda)))
spei12_lambda_plot
ggsave(spei12_lambda_plot, filename = "spei12_lamda_plot.png", width = 6, height = 6 )

spei3_lambda_plot <- ggplot(data = lambda_spei3_df)+
  geom_path(aes(x = spei_value, y = lambda_spei3, group = interaction(Endo,Iteration), color = Endo), alpha = .05, lwd = .4)+
  geom_line(data = lambda_spei3_mean, aes(x = spei_value, y = lambda_spei_mean, group = Endo, color = Endo),lwd = 1)+
  facet_wrap(~Species, scales = "free", nrow = 2)+
  scale_color_manual(values = endophyte_color_scheme[c(2,6)])+
  theme_classic()+ theme(strip.background = element_blank(),
                         strip.text = element_text(face = "italic")) + 
  labs(x = "3-month SPEI", y = expression("Pop. Growth Rate " (lambda)))
spei3_lambda_plot
ggsave(spei3_lambda_plot, filename = "spei3_lamda_plot.png", width = 6, height = 6 )

spei_combo_lambda_plot <- spei3_lambda_plot + spei12_lambda_plot +
  plot_layout(nrow = 2, guides = "collect") + plot_annotation(tag_levels = "A")
ggsave(spei_combo_lambda_plot, filename = "spei_combo_lambda_plot.png", width = 8, height = 10)

# calculating the slope for E+ and E- for each species
lambda_spei3_slopes <- lambda_spei3_mean %>% 
  filter(spei == "spei1" | spei == "spei10") %>% 
  pivot_wider(names_from = c(spei), values_from = c(spei_value, lambda_spei_mean)) %>% 
  mutate(change_lambda = lambda_spei_mean_spei10-lambda_spei_mean_spei1,
         change_spei = spei_value_spei10 - spei_value_spei1,
         slope = change_lambda/change_spei) %>% 
  arrange(Species)
lambda_spei3_slopes_diff <- lambda_spei3_slopes %>% 
  ungroup() %>% 
  dplyr::select(Species, Endo, slope) %>% 
  pivot_wider(names_prefix = expr("slope"), names_from = Endo, values_from = slope) %>% 
  mutate(slope_ratio = abs(`slopeS-`/`slopeS+`))  # interpretation:  ratio>1 means that S- is more sensitive than S+, in line with expectation of drought contributing to buffering
lambda_spei3_slopes_diff$sd_effect <- lambda_sd_diff[1:7,1]
lambda_spei3_slopes_diff$cv_effect <- lambda_cv_diff[1:7,1]
lambda_spei3_slopes_diff$mean_effect <- lambda_mean_diff[1:7,1]

ggplot(data = lambda_spei3_slopes_diff)+
  geom_point(aes(x = sd_effect, y = slope_ratio))
ggplot(data = lambda_spei3_slopes_diff)+
  geom_point(aes(x = cv_effect, y = slope_ratio))
ggplot(data = lambda_spei3_slopes_diff)+
  geom_point(aes(x = mean_effect, y = slope_ratio))


lambda_spei12_slopes <- lambda_spei12_mean %>% 
  filter(spei == "spei1" | spei == "spei10") %>% 
  pivot_wider(names_from = c(spei), values_from = c(spei_value, lambda_spei_mean)) %>% 
  mutate(change_lambda = lambda_spei_mean_spei10-lambda_spei_mean_spei1,
         change_spei = spei_value_spei10 - spei_value_spei1,
         slope = change_lambda/change_spei) %>% 
  arrange(Species)
lambda_spei12_slopes_diff <- lambda_spei12_slopes %>% 
  ungroup() %>% 
  dplyr::select(Species, Endo, slope) %>% 
  pivot_wider(names_prefix = expr("slope"), names_from = Endo, values_from = slope) %>% 
  mutate(slope_ratio = abs(`slopeS-`/`slopeS+`))
lambda_spei12_slopes_diff$sd_effect <- lambda_sd_diff[1:7,1]
lambda_spei12_slopes_diff$cv_effect <- lambda_cv_diff[1:7,1]
lambda_spei12_slopes_diff$mean_effect <- lambda_mean_diff[1:7,1]


ggplot(data = lambda_spei12_slopes_diff)+
  geom_point(aes(x = sd_effect, y = slope_ratio))
ggplot(data = lambda_spei12_slopes_diff)+
  geom_point(aes(x = cv_effect, y = slope_ratio))
ggplot(data = lambda_spei12_slopes_diff)+
  geom_point(aes(x = mean_effect, y = slope_ratio))


######################################################################
########## Calculating some posterior certainties for manuscript #####
######################################################################


spei12_surv_par <- rstan::extract(spei12_surv_fit, pars =quote_bare(beta0,betasize,betasize_2,betaendo,betaspei_endo,tau_year, tau_plot, sigma_year))
spei3_surv_par <- rstan::extract(spei3_surv_fit, pars =quote_bare(beta0,betasize,betasize_2,betaendo,betaspei_endo,tau_year, tau_plot, sigma_year))
# seedling survival
spei12_surv_sdlg_par <- rstan::extract(spei12_surv_fit_seedling, pars =quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))
spei3_surv_sdlg_par <- rstan::extract(spei3_surv_fit_seedling, pars =quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))

# growth
spei12_grow_par <- rstan::extract(spei12_grow_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma, sigma_year))
spei3_grow_par <- rstan::extract(spei3_grow_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma, sigma_year))

# seedling growth
spei12_grow_sdlg_par <- rstan::extract(spei12_grow_fit_seedling, pars = quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma, sigma_year))
spei3_grow_sdlg_par <- rstan::extract(spei3_grow_fit_seedling, pars = quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma, sigma_year))

#flowering
spei12_flow_par <- rstan::extract(spei12_flw_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))
spei3_flow_par <- rstan::extract(spei3_flw_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))

# number of flower tillers
spei12_fert_par <- rstan::extract(spei12_fert_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))
spei3_fert_par <- rstan::extract(spei3_fert_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))

# spikelets/infl
spei12_spike_par <- rstan::extract(spei12_spike_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))
spei3_spike_par <- rstan::extract(spei3_spike_fit, pars = quote_bare(beta0,betasize,betasize_2,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))

# seed mean/spikelet
seed_par <- rstan::extract(seedmean_fit, pars = quote_bare(beta0,betaendo)) #no plot or year effect
# seedling recruitment
spei12_recruit_par <- rstan::extract(spei12_stos_fit, pars = quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))
spei3_recruit_par <- rstan::extract(spei3_stos_fit, pars = quote_bare(beta0,betaendo, betaspei_endo, tau_year, tau_plot, sigma_year))






# Histograms of the slope estimates
surv_post3 <- as.array(spei3_surv_par$betaspei_endo)
surv_post12 <- as.array(spei12_surv_par$betaspei_endo)

surv_sdlg_post3 <- as.array(spei3_surv_sdlg_par$betaspei_endo)
surv_sdlg_post12 <- as.array(spei12_surv_sdlg_par$betaspei_endo)

flow_post3 <- as.array(spei3_flow_par$betaspei_endo)
flow_post12 <- as.array(spei12_flow_par$betaspei_endo)

grow_post3 <- as.array(spei3_grow_par$betaspei_endo)
grow_post12 <- as.array(spei12_grow_par$betaspei_endo)

grow_sdlg_post3 <- as.array(spei3_grow_sdlg_par$betaspei_endo)
grow_sdlg_post12 <- as.array(spei12_grow_sdlg_par$betaspei_endo)

fert_post3 <- as.array(spei3_fert_par$betaspei_endo)
fert_post12 <- as.array(spei12_fert_par$betaspei_endo)

spike_post3 <- as.array(spei3_spike_par$betaspei_endo)
spike_post12 <- as.array(spei12_spike_par$betaspei_endo)

recruit_post3 <- as.array(spei3_recruit_par$betaspei_endo)
recruit_post12 <- as.array(spei12_recruit_par$betaspei_endo)


dimnames(surv_post3) <- dimnames(surv_post12) <- dimnames(surv_sdlg_post3) <- dimnames(flow_post3) <- dimnames(flow_post12) <- dimnames(fert_post3) <- dimnames(fert_post12) <- dimnames(grow_post3) <- dimnames(grow_post12) <- dimnames(spike_post3) <- dimnames(spike_post12) <- dimnames(recruit_post3) <- dimnames(recruit_post12) <- list(Iteration= paste0("i",1:7500), species = paste0("s",1:(n_spp)), Endo = paste0("e",1:2))
dimnames(surv_sdlg_post12) <- dimnames(grow_sdlg_post3) <- dimnames(grow_sdlg_post12) <- list(Iteration= paste0("i",1:15000), species = paste0("s",1:(n_spp)), Endo = paste0("e",1:2))
  
surv_post3_df <- as_tibble(cubelyr::as.tbl_cube(surv_post3)) %>% 
  rename(value = surv_post3) %>% mutate(spei = "spei3",vital_rate = "Survival")
surv_post12_df <- as_tibble(cubelyr::as.tbl_cube(surv_post12)) %>% 
  rename(value = surv_post12) %>% mutate(spei = "spei12",vital_rate = "Survival")

surv_sdlg_post3_df <- as_tibble(cubelyr::as.tbl_cube(surv_sdlg_post3)) %>% 
  rename(value = surv_sdlg_post3) %>% mutate(spei = "spei3",vital_rate = "Seedling Survival")
surv_sdlg_post12_df <- as_tibble(cubelyr::as.tbl_cube(surv_sdlg_post12)) %>% 
  rename(value = surv_sdlg_post12) %>% mutate(spei = "spei12",vital_rate = "Seedling Survival")

flow_post3_df <- as_tibble(cubelyr::as.tbl_cube(flow_post3)) %>% 
  rename(value = flow_post3) %>% mutate(spei = "spei3",vital_rate = "Flowering")
flow_post12_df <- as_tibble(cubelyr::as.tbl_cube(flow_post12)) %>% 
  rename(value = flow_post12) %>% mutate(spei = "spei12",vital_rate = "Flowering")

fert_post3_df <- as_tibble(cubelyr::as.tbl_cube(fert_post3)) %>% 
  rename(value = fert_post3) %>% mutate(spei = "spei3",vital_rate = "Fertility")
fert_post12_df <- as_tibble(cubelyr::as.tbl_cube(fert_post12)) %>% 
  rename(value = fert_post12) %>% mutate(spei = "spei12",vital_rate = "Fertility")

grow_post3_df <- as_tibble(cubelyr::as.tbl_cube(grow_post3)) %>% 
  rename(value = grow_post3) %>% mutate(spei = "spei3",vital_rate = "Growth")
grow_post12_df <- as_tibble(cubelyr::as.tbl_cube(grow_post12)) %>% 
  rename(value = grow_post12) %>% mutate(spei = "spei12",vital_rate = "Growth")

grow_sdlg_post3_df <- as_tibble(cubelyr::as.tbl_cube(grow_sdlg_post3)) %>% 
  rename(value = grow_sdlg_post3) %>% mutate(spei = "spei3",vital_rate = "Seedling Growth")
grow_sdlg_post12_df <- as_tibble(cubelyr::as.tbl_cube(grow_sdlg_post12)) %>% 
  rename(value = grow_sdlg_post12) %>% mutate(spei = "spei12",vital_rate = "Seedling Growth")

spike_post3_df <- as_tibble(cubelyr::as.tbl_cube(spike_post3)) %>% 
  rename(value = spike_post3) %>% mutate(spei = "spei3",vital_rate = "Spike")
spike_post12_df <- as_tibble(cubelyr::as.tbl_cube(spike_post12)) %>% 
  rename(value = spike_post12) %>% mutate(spei = "spei12",vital_rate = "Spike")

recruit_post3_df <- as_tibble(cubelyr::as.tbl_cube(recruit_post3)) %>% 
  rename(value = recruit_post3) %>% mutate(spei = "spei3",vital_rate = "Recruit")
recruit_post12_df <- as_tibble(cubelyr::as.tbl_cube(recruit_post12)) %>% 
  rename(value = recruit_post12) %>% mutate(spei = "spei12",vital_rate = "Recruit")


spei_posterior <- rbind(surv_post3_df, surv_post12_df,
                        surv_sdlg_post3_df, surv_sdlg_post12_df,
                        flow_post3_df, flow_post12_df,
                        fert_post3_df, fert_post12_df,
                        grow_post3_df, grow_post12_df,
                        grow_sdlg_post3_df, grow_sdlg_post12_df,
                        spike_post3_df, spike_post12_df,
                        recruit_post3_df, recruit_post12_df) %>% 
  mutate(if_positive = case_when(value>0 ~ TRUE, value<0 ~ FALSE))


# calculating the probability of positive slopes for the MS
spei_post_summary <- spei_posterior %>% 
  group_by(vital_rate, species, Endo, spei) %>% 
  summarize(iterations = n(),
            mean = mean(value),
            if_positive = sum(value>0)/iterations)



