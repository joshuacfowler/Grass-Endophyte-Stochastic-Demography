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
library(actuar)
library(rstan)
library(patchwork) # for putting plots together
library(ggrepel) # for labels that don't overlap


# library(gridExtra)
# library(grid)
# library(cowplot) # for pulling legend from ggplots
# library(cubelyr) # for working between lists of matrixes and dataframes


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

# surv_fit_seedling <- readRDS(url("https://www.dropbox.com/s/vf1mju5u4c4fs3t/endo_seedling_surv.rds?dl=1"))
# surv_fit <- readRDS(url("https://www.dropbox.com/s/00bor35inv5dypd/endo_spp_surv_woseedling.rds?dl=1"))
# grow_fit_seedling <- readRDS(url("https://www.dropbox.com/s/m0mw5z29slpm4p7/endo_seedling_grow_PIG_10000iterations.rds?dl=1"))
# grow_fit <- readRDS(url("https://www.dropbox.com/s/0ze8aooi9axj3oq/endo_spp_grow_PIG.rds?dl=1"))
# flw_fit <- readRDS(url("https://www.dropbox.com/s/ej65pn5k0km0z9c/endo_spp_flw.rds?dl=1"))
# fert_fit <- readRDS(url("https://www.dropbox.com/s/pk4x1j97kazu6pb/endo_spp_fert_pig.rds?dl=1"))
# spike_fit <- readRDS(url("https://www.dropbox.com/s/pjgui0n9tng6427/endo_spp_spike_year_plot_nb.rds?dl=1"))
# seedmean_fit <- readRDS(url("https://www.dropbox.com/s/3ma5yc8iusu8bh0/endo_spp_seed_mean.rds?dl=1"))
# stos_fit <- readRDS(url("https://www.dropbox.com/s/nf50hd76iw3hucw/endo_spp_s_to_s.rds?dl=1"))

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
# plot(surv_par$tau_year[,1,1,], grow_par$tau_year[,1,1,], col = levels(as_factor(grow_par$tau_year[1,1,1,])))
#############################################################################################
####### Run the MPM ------------------
#############################################################################################

# make the list of parameters and calculate mean lambdas
n_draws <- 500 # the means are the same whether we do 500 or 1000 draws
post_draws <- sample.int(7500,size=n_draws) # The models except for seedling growth have 7500 iterations. That one has more (15000 iterations) to help it converge.

lambda_mean <- array(dim = c(8,2,n_draws))
for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      lambda_mean[s,e,i] <- lambda(bigmatrix(make_params(species=s,
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
                                                         recruit_par=recruit_par), 
                                             extension = 100)$MPMmat) # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
    lambda_mean[8,e,i] <- mean(lambda_mean[1:7,e,i])
  }
}
# saving lambda_mean with 500 post draws to dropbox
saveRDS(lambda_mean, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_mean.rds")
lambda_mean <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_mean.rds")



# Mean endophyte difference and quantiles
lambda_means <- matrix(NA,8,2)
lambda_mean_diff <- matrix(NA,8,7)
for(s in 1:8){
  lambda_means[s,1] <- mean(lambda_mean[s,1,])
  lambda_means[s,2] <- mean(lambda_mean[s,2,])
  lambda_mean_diff[s,1] = mean(lambda_mean[s,2,] - lambda_mean[s,1,])
  lambda_mean_diff[s,2:7] = quantile(lambda_mean[s,2,] - lambda_mean[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
}

# look at the lambda values for a general gut check
lambda_means

#calculating the percent increase for each species
lambda_means_percentdiff <- 100-(lambda_means[,1]/lambda_means[,2]*100)

## now do variance in lambda 

lambda_hold <- array(dim = c(13,7,2,n_draws)) # 13 years because of reproduction measured in year t1; needs to be before growth, so no year 1
lambda_var <- array(dim = c(8,2,n_draws))
for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      for(y in 1:13){
        
        lambda_hold[y,s,e,i] <- lambda(bigmatrix(make_params(species=s,
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
                                                 extension = 100)$MPMmat) # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
      }
      lambda_var[s,e,i] <- sd(lambda_hold[,s,e,i]) # we calulate the standard deviation here
    }
    lambda_var[8,e,i] <- mean(lambda_var[1:7,e,i])
  }
}
#saving the yearly lambdas and the sd of lambdas to dropbox
saveRDS(lambda_hold, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_hold.rds")
saveRDS(lambda_var, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_var.rds")
lambda_hold <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_hold.rds")
lambda_var <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_var.rds")



#Calculationg endophyte effect on sd and variance
lambda_sds <- matrix(NA,8,2)
lambda_vars <- matrix(NA,8,2)

lambda_cv <- (lambda_var^2)/(2*lambda_mean^2) # this is the "variance penalty" term not the true CV right now

lambda_cvs <- matrix(NA,8,2)

lambda_sd_diff <- matrix(NA,8,7)
lambda_var_diff <- matrix(NA,8,7)
lambda_cv_diff <-  matrix(NA,8,7)
for(s in 1:8){
  lambda_sds[s,1] <- mean(lambda_var[s,1,])
  lambda_sds[s,2] <- mean(lambda_var[s,2,])
  
  lambda_vars[s,1] <- mean(lambda_var[s,1,])^2
  lambda_vars[s,2] <- mean(lambda_var[s,2,])^2
  
  lambda_cvs[s,1] <- mean(lambda_cv[s,1,])
  lambda_cvs[s,2] <- mean(lambda_cv[s,2,])
  
  lambda_sd_diff[s,1] = mean(lambda_var[s,2,] - lambda_var[s,1,])
  lambda_sd_diff[s,2:7] = quantile(lambda_var[s,2,] - lambda_var[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
  
  lambda_var_diff[s,1] = mean(lambda_var[s,2,]^2 - lambda_var[s,1,]^2)
  lambda_var_diff[s,2:7] = quantile(lambda_var[s,2,]^2 - lambda_var[s,1,]^2,probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
  
  lambda_cv_diff[s,1] = mean(lambda_cv[s,2,] - lambda_cv[s,1,])
  lambda_cv_diff[s,2:7] = quantile(lambda_cv[s,2,]^2 - lambda_cv[s,1,]^2,probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
  
}

# calcualating the percent different in standard deviation
lambda_sds_percentdiff <- 100 - (lambda_sds[,1]/lambda_sds[,2])*100
# calcualating the percent different in variance
lambda_vars_percentdiff <- 100 - (lambda_vars[,1]/lambda_vars[,2])*100
lambda_cvs_percentdiff <- 100 - (lambda_cvs[,1]/lambda_cvs[,2])*100

# Making plots of yearly lambdas
yearly_lambda <- array(dim = c(13,7,2))

for(y in 1:13){
  for(s in 1:7){
    for( e in 1:2){
      yearly_lambda[y,s,e] <- mean(lambda_hold[y,s,e,])
    }
  }
}
dimnames(yearly_lambda) <- list(Year = c(1:13), Species = paste0("s",1:7), Endo = paste0("e",1:2))
yearly_lambda_cube <- cubelyr::as.tbl_cube(yearly_lambda)
yearly_lambda_df <- as_tibble(yearly_lambda_cube)

dimnames(lambda_hold) <- list(Year = c(1:13), Species = paste0("s",1:7), Endo = paste0("e",1:2), Iterations = paste0("i",1:n_draws))
yearly_lambda_iter_cube <- cubelyr::as.tbl_cube(lambda_hold)
yearly_lambda_iter_df <- as_tibble(yearly_lambda_iter_cube)

ggplot(data = yearly_lambda_df)+
  geom_hline(yintercept = 1)+
  geom_path(aes(x = Year, y = yearly_lambda, color = Endo))+
  facet_wrap(~Species)+
  theme_minimal()

ggplot(data = yearly_lambda_df)+
  geom_vline(xintercept = 1)+
  geom_histogram(aes(x = yearly_lambda, fill = Endo),alpha = .5, position = "identity", bins = 10)+
  facet_wrap(~Species)+
  theme_minimal()

y_df <- yearly_lambda_df %>% 
  pivot_wider(id_cols = c(Year,Species), names_from = Endo, values_from = yearly_lambda)
ggplot(data = y_df)+
  geom_point(aes(x = e1, y = e2),alpha = .5)+
  facet_wrap(~Species)+
  theme_minimal()
################################################################
##### Plot of mean and variance endo effect on lambda ##########
################################################################

lambda_mean_diff_df <- as_tibble(lambda_mean_diff)  %>% 
  rename( "mean" = V1, fifth = V2, twelfthpointfive = V3, twentyfifth = V4, seventyfifth = V5, eightyseventhpointfive = V6, ninetyfifth = V7) %>% 
  mutate(rownames = row.names(.)) %>% 
  mutate(species = case_when(rownames == 1 ~ "Agrostis perennans",
                             rownames == 2 ~ "Elymus villosus",
                             rownames == 3 ~ "Elymus virginicus",
                             rownames == 4 ~ "Festuca subverticillata",
                             rownames == 5 ~ "Lolium arundinaceum",
                             rownames == 6 ~ "Poa alsodes",
                             rownames == 7 ~ "Poa sylvestris",
                             rownames == 8 ~ "Species Mean"))

  
ggplot(data = lambda_mean_diff_df) +
  geom_hline(yintercept = 0, col = "black") + 
  geom_linerange(aes(y = mean, x = species, ymin = 0, ymax = mean), color = "white", lwd =4)+
  geom_point(aes(y = mean, x = species, color = species), lwd = 4) +
  geom_linerange(aes(y = mean, x = species, ymin = twentyfifth, ymax = seventyfifth, color = species), lwd = 2) +
  geom_linerange(aes(y = mean, x = species, ymin = twelfthpointfive, ymax = eightyseventhpointfive, color = species), lwd = 1) +
  geom_linerange(aes(y = mean, x = species, ymin = fifth, ymax = ninetyfifth, color = species)) +
  scale_color_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
  coord_flip()+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = NA))
  
# Version with raw posterior draws
dimnames(lambda_mean) <- list(Species = paste0("s",1:8), Endo = paste0("e",1:2), Iteration= paste0("i",1:n_draws))
lambda_mean_cube <- cubelyr::as.tbl_cube(lambda_mean)
lambda_mean_df <- as_tibble(lambda_mean_cube) %>% 
  pivot_wider(names_from = Endo, values_from = lambda_mean) %>% 
  mutate(lambda_diff = e2-e1) %>% 
  mutate(species = case_when(Species == "s1" ~ "Agrostis perennans",
                           Species == "s2" ~ "Elymus villosus",
                           Species == "s3" ~ "Elymus virginicus",
                           Species == "s4" ~ "Festuca subverticillata",
                           Species == "s5" ~ "Lolium arundinaceum",
                           Species == "s6" ~ "Poa alsodes",
                           Species == "s7" ~ "Poa sylvestris",
                           Species == "s8" ~ "Species Mean"))

meanlambda_plot <- ggplot(data = lambda_mean_df) +
  geom_hline(yintercept = 0, col = "black") + 
  # geom_linerange(data = lambda_mean_diff_df, aes(x = species, y = mean, ymin = 0, ymax = mean, color = species)) + 
  geom_jitter( aes(y = lambda_diff, x = species, color = species), width = .2, alpha = .2) +
  stat_summary(aes(y = lambda_diff, x = species), fun = mean,geom = "point", size = 3) +
  geom_point(data = lambda_mean_diff_df, aes(y = mean, x = species, color = species), lwd = 2) +
  scale_color_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
  coord_flip(ylim = c(-.5,.5)) +
  labs(y = expression(paste("Symbiosis Effect on ", bar(lambda))),
       x = "",
       color = "Species")+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = NA),
        axis.line = element_line(color = "grey"),
        axis.text.y = element_text(face = "italic"),
        legend.position = "none")
meanlambda_plot
ggsave(meanlambda_plot, filename = "meanlambda_plot.png", width = 4, height = 5)

# now for the effect on variance plot
lambda_cv_diff_df <- as_tibble(lambda_cv_diff)  %>% 
  rename( "mean" = V1, fifth = V2, twelfthpointfive = V3, twentyfifth = V4, seventyfifth = V5, eightyseventhpointfive = V6, ninetyfifth = V7) %>% 
  mutate(rownames = row.names(.)) %>% 
  mutate(species = case_when(rownames == 1 ~ "Agrostis perennans",
                             rownames == 2 ~ "Elymus villosus",
                             rownames == 3 ~ "Elymus virginicus",
                             rownames == 4 ~ "Festuca subverticillata",
                             rownames == 5 ~ "Lolium arundinaceum",
                             rownames == 6 ~ "Poa alsodes",
                             rownames == 7 ~ "Poa sylvestris",
                             rownames == 8 ~ "Species Mean"))


# ggplot(data = lambda_var_diff_df) +
#   geom_hline(yintercept = 0, col = "black") + 
#   geom_linerange(aes(y = mean, x = species, ymin = 0, ymax = mean), color = "white", lwd =4) +
#   geom_point(aes(y = mean, x = species, color = species), lwd = 4) +
#   geom_linerange(aes(y = mean, x = species, ymin = twentyfifth, ymax = seventyfifth, color = species), lwd = 2) +
#   geom_linerange(aes(y = mean, x = species, ymin = twelfthpointfive, ymax = eightyseventhpointfive, color = species), lwd = 1) +
#   geom_linerange(aes(y = mean, x = species, ymin = fifth, ymax = ninetyfifth, color = species)) +
#   scale_color_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
#   coord_flip()+
#   theme(panel.background = element_rect(fill = "white"),
#         panel.grid = element_line(color = NA))


# Version with raw posterior draws
dimnames(lambda_cv) <- list(Species = paste0("s",1:8), Endo = paste0("e",1:2), Iteration= paste0("i",1:n_draws))
lambda_cv_cube <- cubelyr::as.tbl_cube(lambda_cv)
lambda_cv_df <- as_tibble(lambda_cv_cube) %>% 
  pivot_wider(names_from = Endo, values_from = lambda_cv) %>% 
  mutate(lambda_diff = e2-e1) %>% 
  mutate(species = case_when(Species == "s1" ~ "Agrostis perennans",
                             Species == "s2" ~ "Elymus villosus",
                             Species == "s3" ~ "Elymus virginicus",
                             Species == "s4" ~ "Festuca subverticillata",
                             Species == "s5" ~ "Lolium arundinaceum",
                             Species == "s6" ~ "Poa alsodes",
                             Species == "s7" ~ "Poa sylvestris",
                             Species == "s8" ~ "Species Mean"))

lambdavar_plot <- ggplot(data = lambda_cv_df) +
  geom_hline(yintercept = 0, col = "black") + 
  # geom_linerange(data = lambda_var_diff_df, aes(x = species, y = mean, ymin = 0, ymax = mean, color = species))+
  geom_jitter( aes(y = lambda_diff, x = species, color = species), width = .2, alpha = .2) +
  stat_summary(aes(y = lambda_diff, x = species), fun = mean,geom = "point", size = 3) +
  geom_point(data = lambda_cv_diff_df, aes(y = mean, x = species, color = species), lwd = 2) +
  scale_color_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
  coord_flip(ylim = c(-.75,.25)) +  #There's a few iterations out at -1.5
  labs(y = expression(paste("Symbiosis Effect on ", "CV(",lambda,")")),
       x = "",
       color = "Species")+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = NA),
        axis.line = element_line(color = "grey"),
        axis.text.y = element_text(face = "italic"),
        legend.position = "none")
lambdavar_plot
ggsave(lambdavar_plot, filename = "lambdavar_plot_coefficientofvariation.png", width = 4, height = 5)

endo_lambdaeffects_plot <-  meanlambda_plot +lambdavar_plot + plot_layout(nrow = 1, guides = "collect")
endo_lambdaeffects_plot
ggsave(endo_lambdaeffects_plot, filename = "endo_lambdaeffects_plot.png", width = 12, height = 6)



# a plot of variance by mean effects
lambda_var_join <- lambda_cv_df %>% 
  rename( varlambda_diff = lambda_diff, var_e1 = e1, var_e2 = e2)

lambda_mean_join <- lambda_mean_df %>% 
  rename( meanlambda_diff = lambda_diff, mean_e1 = e1, mean_e2 = e2)

lambda_join_df <- lambda_mean_join %>% 
  full_join(lambda_var_join, by = c("Species", "Iteration", "species"))

summarylambda_join_df <- lambda_join_df %>% 
  group_by(Species, species) %>% 
  summarize(avg_meandiff = mean(meanlambda_diff),
            avg_vardiff = mean(varlambda_diff))

agpe_meanvar <- ggplot(data = filter(lambda_join_df, species == "Agrostis perennans"))+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  geom_point(aes(x = meanlambda_diff, y = varlambda_diff, color = species), alpha = .3)+
  scale_color_manual(values = c("#dbdb42"))+
  labs(title = "AGPE")+
  theme(panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none")
# agpe_meanvar

elri_meanvar <- ggplot(data = filter(lambda_join_df, species == "Elymus villosus"))+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  geom_point(aes(x = meanlambda_diff, y = varlambda_diff, color = species), alpha = .3)+
  scale_color_manual(values = c("#b8e3a0"))+
  labs(title = "ELRI")+
  theme(panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none")
# elri_meanvar

elvi_meanvar <- ggplot(data = filter(lambda_join_df, species == "Elymus virginicus"))+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  geom_point(aes(x = meanlambda_diff, y = varlambda_diff, color = species), alpha = .3)+
  scale_color_manual(values = c("#7fcdbb"))+
  labs(title = "ELVI")+
  theme(panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none")
# elvi_meanvar
  
fesu_meanvar <- ggplot(data = filter(lambda_join_df, species == "Festuca subverticillata"))+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  geom_point(aes(x = meanlambda_diff, y = varlambda_diff, color = species), alpha = .3)+
  scale_color_manual(values = c("#41b6c4"))+
  labs(title = "FESU")+
  theme(panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none")
# fesu_meanvar

loar_meanvar <- ggplot(data = filter(lambda_join_df, species == "Lolium arundinaceum"))+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  geom_point(aes(x = meanlambda_diff, y = varlambda_diff, color = species), alpha = .3)+
  scale_color_manual(values = c("#1d91c0"))+
  labs(title = "LOAR")+
  theme(panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none")
# loar_meanvar

poal_meanvar <- ggplot(data = filter(lambda_join_df, species == "Poa alsodes"))+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  geom_point(aes(x = meanlambda_diff, y = varlambda_diff, color = species), alpha = .3)+
  scale_color_manual(values = c("#225ea8"))+
  labs(title = "POAL")+
  theme(panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none")
# poal_meanvar

posy_meanvar <- ggplot(data = filter(lambda_join_df, species == "Poa sylvestris"))+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  geom_point(aes(x = meanlambda_diff, y = varlambda_diff, color = species), alpha = .3)+
  scale_color_manual(values = c("#0c2c84"))+
  labs(title = "POSY")+
  theme(panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none")
# posy_meanvar

sppmean_meanvar <- ggplot(data = filter(lambda_join_df, species == "Species Mean"))+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  geom_point(aes(x = meanlambda_diff, y = varlambda_diff, color = species), alpha = .3)+
  scale_color_manual(values = c("#A9A9A9"))+
  labs(title = "Species Mean")+
  theme(panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none")
# sppmean_meanvar
meanvar_biplot <- ggplot(data = lambda_join_df) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  # geom_point(aes(x = meanlambda_diff, y = varlambda_diff, color = species), alpha = .5) +
  geom_point(data = summarylambda_join_df, aes(x = avg_meandiff, y = avg_vardiff, color = species), lwd  = 3) +
  geom_label_repel(data = summarylambda_join_df, aes(x = avg_meandiff, y = (avg_vardiff), label = species), fontface = "italic", lwd  = 3) +
  xlim(-.05,.17) + ylim(-.26,.05)+
  # xlim(-.05,.2) + ylim(-.08,.05)+
  scale_color_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
  xlab(expression(paste("Effect on ", bar(lambda)))) +
  ylab(expression(paste("Effect on ", "CV(",lambda,")"))) +
  theme(panel.background = element_blank(), 
        axis.line = element_blank(),
        legend.position = "none",)
meanvar_biplot
ggsave(meanvar_biplot, filename = "meanvar_biplot.png", width = 5, height = 5)

meanvar_biplot_post <- meanvar_biplot + (agpe_meanvar+ elri_meanvar + elvi_meanvar +fesu_meanvar + loar_meanvar + poal_meanvar + posy_meanvar + sppmean_meanvar )+
  plot_layout(ncol = 2, nrow = 1, widths = unit(c(9, 1), c("cm", "null")), heights = unit(c(9, 1), c('cm', 'null')))
meanvar_biplot_post
ggsave(meanvar_biplot_post, filename = "meanvar_biplot_post.png", width = 10, height = 5)


meanvar_biplot_with_
#############################################################################################
####### Stochastic lambda simulations ------------------
#############################################################################################

# Make a list of year specific transition matrices for with and without variance and mean effects of endophytes
lambdaS_out <- array(dim = c(8,4,n_draws))
for(i in 1:length(post_draws)){
  for(s in 1:7){
    eminus_list <- eplus_list <- eplus__mean_only_list <- eplus__var_only_list <- list()
    for(y in 1:11){ ## why is y looped 11 times here?
      eminus_list[[y]] <- bigmatrix(make_params(species=s,
                                                endo_mean=0,
                                                endo_var=0,
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
                                    extension = 100)$MPMmat
      eplus_list[[y]] <- bigmatrix(make_params(species=s,
                                               endo_mean=1,
                                               endo_var=1,
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
                                   extension = 100)$MPMmat
      eplus__mean_only_list[[y]] <- bigmatrix(make_params(species=s,
                                                          endo_mean=1,
                                                          endo_var=0,
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
                                              extension = 100)$MPMmat
      eplus__var_only_list[[y]] <- bigmatrix(make_params(species=s,
                                                         endo_mean=0,
                                                         endo_var=1,
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
                                             extension = 100)$MPMmat
    }
    lambdaS_out[s,1,i] <- lambdaSim(eminus_list)$lambdaS
    lambdaS_out[s,2,i] <- lambdaSim(eplus__mean_only_list)$lambdaS
    lambdaS_out[s,3,i] <- lambdaSim(eplus__var_only_list)$lambdaS
    lambdaS_out[s,4,i] <- lambdaSim(eplus_list)$lambdaS
  }
  lambdaS_out[8,1,i] <- mean(lambdaS_out[1:7,1,i]) # species mean eminus
  lambdaS_out[8,2,i] <- mean(lambdaS_out[1:7,2,i]) # species mean eplus mean only
  lambdaS_out[8,3,i] <- mean(lambdaS_out[1:7,3,i]) # species mean eplus var only
  lambdaS_out[8,4,i] <- mean(lambdaS_out[1:7,4,i]) # species mean eplus
}

saveRDS(lambdaS_out, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambdaS_out.rds")
lambdaS_out <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambdaS_out.rds")


lambdaS_diff <- lambdaS_diff_mean_only <- lambdaS_diff_var_only <- lambdaS_diff_mean_var_inter <- matrix(NA,8,7)
for(s in 1:8){
  lambdaS_diff[s,1] = mean(lambdaS_out[s,4,] - lambdaS_out[s,1,], na.rm = T) # eplus - eminus
  lambdaS_diff[s,2:7] = quantile(lambdaS_out[s,4,] - lambdaS_out[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95), na.rm = T)
  
  lambdaS_diff_mean_var_inter[s,1] = mean(lambdaS_out[s,4,] - lambdaS_out[s,1,] -(lambdaS_out[s,2,]-lambdaS_out[s,1,]) - (lambdaS_out[s,3,]-lambdaS_out[s,1,]), na.rm = T) # Trying to calcuallte the contribution of the mean-variance interaction
  lambdaS_diff_mean_var_inter[s,2:7] = quantile(lambdaS_out[s,4,] - lambdaS_out[s,1,]-(lambdaS_out[s,2,]-lambdaS_out[s,1,]) - (lambdaS_out[s,3,]-lambdaS_out[s,1,]),probs=c(0.05,0.125,0.25,0.75,0.875,0.95), na.rm = T)
  
  lambdaS_diff_mean_only[s,1] = mean(lambdaS_out[s,2,] - lambdaS_out[s,1,], na.rm = T) # eplus mean only - eminus
  lambdaS_diff_mean_only[s,2:7] = quantile(lambdaS_out[s,2,] - lambdaS_out[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95), na.rm = T)
  
  lambdaS_diff_var_only[s,1] = mean(lambdaS_out[s,3,] - lambdaS_out[s,1,], na.rm = T) # eplus var only - eminus
  lambdaS_diff_var_only[s,2:7] = quantile(lambdaS_out[s,3,] - lambdaS_out[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95), na.rm = T)
}
# sanity check base r plot
lambdaS_effect <- rowMeans(lambdaS_out[8,2:4,]-lambdaS_out[8,1,], na.rm = T)
barplot(lambdaS_effect,col=rainbow(3),
        names.arg = c("Mean effect","Variance effect","Total effect"))
################################################################
##### Plot of stochastic lambda contributions
################################################################
# I need to check the indexing of the variance endophyte effect because it seems like it is backwards... (total effect is less than mean + var)
# I am relabeling them here for ESA, but I need to fix this in the vital rate models.
lambdaS_diff_df <- as_tibble(lambdaS_diff) %>%  
  rename( "mean" = V1, fifth = V2, twelfthpointfive = V3, twentyfifth = V4, seventyfifth = V5, eightyseventhpointfive = V6, ninetyfifth = V7) %>% 
  mutate(rownames = row.names(.)) %>% 
  mutate(species = case_when(rownames == 1 ~ "Agrostis perennans",
                             rownames == 2 ~ "Elymus villosus",
                             rownames == 3 ~ "Elymus virginicus",
                             rownames == 4 ~ "Festuca subverticillata",
                             rownames == 5 ~ "Lolium arundinaceum",
                             rownames == 6 ~ "Poa alsodes",
                             rownames == 7 ~ "Poa sylvestris",
                             rownames == 8 ~ "Species Mean")) %>% 
  mutate(contribution = "Total")
lambdaS_diff_mean_var_inter_df <- as_tibble(lambdaS_diff_mean_var_inter) %>%  
  rename( "mean" = V1, fifth = V2, twelfthpointfive = V3, twentyfifth = V4, seventyfifth = V5, eightyseventhpointfive = V6, ninetyfifth = V7) %>% 
  mutate(rownames = row.names(.)) %>% 
  mutate(species = case_when(rownames == 1 ~ "Agrostis perennans",
                             rownames == 2 ~ "Elymus villosus",
                             rownames == 3 ~ "Elymus virginicus",
                             rownames == 4 ~ "Festuca subverticillata",
                             rownames == 5 ~ "Lolium arundinaceum",
                             rownames == 6 ~ "Poa alsodes",
                             rownames == 7 ~ "Poa sylvestris",
                             rownames == 8 ~ "Species Mean")) %>% 
  mutate(contribution = "Interaction")
lambdaS_diff_mean_only_df <- as_tibble(lambdaS_diff_mean_only) %>%  
  rename( "mean" = V1, fifth = V2, twelfthpointfive = V3, twentyfifth = V4, seventyfifth = V5, eightyseventhpointfive = V6, ninetyfifth = V7) %>% 
  mutate(rownames = row.names(.)) %>% 
  mutate(species = case_when(rownames == 1 ~ "Agrostis perennans",
                             rownames == 2 ~ "Elymus villosus",
                             rownames == 3 ~ "Elymus virginicus",
                             rownames == 4 ~ "Festuca subverticillata",
                             rownames == 5 ~ "Lolium arundinaceum",
                             rownames == 6 ~ "Poa alsodes",
                             rownames == 7 ~ "Poa sylvestris",
                             rownames == 8 ~ "Species Mean")) %>% 
  mutate(contribution = "Mean")
lambdaS_diff_var_only_df <- as_tibble(lambdaS_diff_var_only) %>%  
  rename( "mean" = V1, fifth = V2, twelfthpointfive = V3, twentyfifth = V4, seventyfifth = V5, eightyseventhpointfive = V6, ninetyfifth = V7) %>% 
  mutate(rownames = row.names(.)) %>% 
  mutate(species = case_when(rownames == 1 ~ "Agrostis perennans",
                             rownames == 2 ~ "Elymus villosus",
                             rownames == 3 ~ "Elymus virginicus",
                             rownames == 4 ~ "Festuca subverticillata",
                             rownames == 5 ~ "Lolium arundinaceum",
                             rownames == 6 ~ "Poa alsodes",
                             rownames == 7 ~ "Poa sylvestris",
                             rownames == 8 ~ "Species Mean")) %>% 
  mutate(contribution = "Var")

lambdaS_diff_all_df <- bind_rows(lambdaS_diff_df, lambdaS_diff_mean_var_inter_df, lambdaS_diff_mean_only_df, lambdaS_diff_var_only_df)
lambdaS_diff_all_df$contribution <- factor(lambdaS_diff_all_df$contribution, levels = c("Var", "Mean", "Interaction", "Total"))
x_labels <- c(expression(paste("Var(", bar(lambda), ")")), expression(bar(lambda)),"Interaction", "Total")
byspp_stochplot <- ggplot(data = lambdaS_diff_all_df) +
  geom_hline(yintercept = 0, col = "black") + 
  geom_linerange(aes(y = mean, x = contribution, ymin = twentyfifth, ymax = seventyfifth, color = species), lwd = 2) +
  geom_linerange(aes(y = mean, x = contribution, ymin = twelfthpointfive, ymax = eightyseventhpointfive, color = species), lwd = 1) +
  geom_linerange(aes(y = mean, x = contribution, ymin = fifth, ymax = ninetyfifth, color = species)) +
  geom_point(aes(y = mean, x = contribution, fill = species, pch = contribution), lwd = 4) + 
  facet_wrap(~species, nrow = 2, scales = "free") + coord_flip() + 
  scale_shape_manual(values = c(21, 22, 23,1)) +
  scale_fill_manual(values = c("#59A1BC", "#4E816D", "#EFAD3A", "#397BB7", "#9E78A1", "#E04D55", "#9D5251", "#5C5C5D")) +
  scale_color_manual(values = c("#59A1BC", "#4E816D", "#EFAD3A", "#397BB7", "#9E78A1", "#E04D55", "#9D5251", "#5C5C5D")) +
  xlab("") + scale_x_discrete(labels = x_labels) +
  ylab(expression(paste("Endophyte effect on", " ", lambda["s"]))) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(size = .5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = rel(1.5), face = "bold"),
        axis.text.x = element_text(size = rel(1.5), face = "bold"),
        axis.title = element_text(size = rel(1.5)),
        strip.background = element_blank(),
        legend.position="none")
byspp_stochplot
ggsave(byspp_stochplot, filename = "~/Documents/byspp_stochplot.tiff", width = 10, height = 5)

# Just the mean effect
sppmean_stochplot <- ggplot(data = subset(lambdaS_diff_all_df, species == "Species Mean")) +
  geom_hline(yintercept = 0, col = "black") + 
  geom_linerange(aes(y = mean, x = contribution, ymin = twentyfifth, ymax = seventyfifth, color = species), lwd = 2) +
  geom_linerange(aes(y = mean, x = contribution, ymin = twelfthpointfive, ymax = eightyseventhpointfive, color = species), lwd = 1) +
  geom_linerange(aes(y = mean, x = contribution, ymin = fifth, ymax = ninetyfifth, color = species)) +
  geom_point(aes(y = mean, x = contribution, fill = species, pch = contribution), lwd = 4) + 
  facet_wrap(~species, nrow = 2, scales = "free") + coord_flip() + 
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c( "#5C5C5D")) +
  scale_color_manual(values = c("#5C5C5D")) +
  xlab("Species Mean Contributions") + scale_x_discrete(labels = x_labels) +
  ylab(expression(paste("Endophyte effect on", " ", lambda["s"]))) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(size = .5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = rel(2), face = "bold"),
        axis.text.x = element_text(size = rel(2), face = "bold"),
        axis.title = element_text(size = rel(2)),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position="none")
sppmean_stochplot
ggsave(sppmean_stochplot, filename = "sppmean_stochplot.tiff",  width = 6, height = 4.5, bg = "transparent")


lambdaS_percent <- lambdaS_diff_all_df %>% 
  group_by(species) %>% 
  summarize(percent_var = mean[contribution == "Var"]/mean[contribution == "Total"],
            percent_mean = mean[contribution == "Mean"]/mean[contribution == "Total"],
            percent_inter = mean[contribution == "Interaction"]/mean[contribution == "Total"])







################################################################
########### getting yearly lambda values #####
eminus_yearly <- eplus_yearly <- array(dim = c(1,11,n_draws)) # just doing this for one one species for now


  for(y in 1:11){
    for(i in 1:length(post_draws)){
    eminus_yearly[1,y,i] <- lambda(bigmatrix(make_params(species=4,
                                       endo_mean=0,
                                       endo_var=0,
                                       original = 1, # should be =1 to represent recruit
                                       draw=post_draws[i],
                                       max_size=max_size,
                                       rfx=T,
                                       year=y,surv_par=surv_par,
                                       surv_sdlg_par = surv_sdlg_par,
                                       grow_par=grow_par,
                                       grow_sdlg_par = grow_sdlg_par,
                                       flow_par=flow_par,
                                       fert_par=fert_par,
                                       spike_par=spike_par,
                                       seed_par=seed_par,
                                       recruit_par=recruit_par))$MPMmat)
    eplus_yearly[1,y,i] <- lambda(bigmatrix(make_params(species=4,
            endo_mean=1,
            endo_var=1,
            original = 1, # should be =1 to represent recruit
            draw=post_draws[i],
            max_size=max_size,
            rfx=T,
            year=y,surv_par=surv_par,
            surv_sdlg_par = surv_sdlg_par,
            grow_par=grow_par,
            grow_sdlg_par = grow_sdlg_par,
            flow_par=flow_par,
            fert_par=fert_par,
            spike_par=spike_par,
            seed_par=seed_par,
            recruit_par=recruit_par))$MPMmat)

}}

saveRDS(eplus_yearly, file = "~/Documents/eplus_yearly.rds")
saveRDS(eminus_yearly, file = "~/Documents/eminus_yearly.rds")


mean_eplus_yearly <- matrix(NA, 11,7)
mean_eminus_yearly <- matrix(NA, 11,7)

for(y in 1:11){
mean_eplus_yearly[y,1] <- mean(eplus_yearly[1,y,])
mean_eplus_yearly[y,2:7] <- quantile(eplus_yearly[1,y,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95), na.rm = T)

mean_eminus_yearly[y,1] <- mean(eminus_yearly[1,y,])
mean_eplus_yearly[y,2:7] <- quantile(eminus_yearly[1,y,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95), na.rm = T)

}
spp_names <- as_tibble(names(species_factor_key))

focal_spp <- "FESU"
focal_index <- which(spp_names==focal_spp)
focal_month <- c(9,7,7,5,7,5,5,NA)[focal_index]


# read in PPT and Temp
climate_data <- read_csv(file = "~/Dropbox/EndodemogData/PRISMClimateData_BrownCo.csv") %>% 
  mutate(year = year(Date), month = month(Date), day = day(Date)) %>% 
  rename(ppt = `ppt (mm)`, tmean = `tmean (degrees C)`) %>% 
  mutate(site_lat = 39.235900000000, site_long = -86.218100000000) %>% 
  mutate(census_month = focal_month, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>% 
  filter(climate_year != 2006) %>% 
  group_by(climate_year) %>% 
  summarize(`Cumulative PPT (mm)` = sum(ppt),
            `Mean Temp. (C˚)` = mean(tmean)) %>% 
  filter(climate_year <= max(LTREB_full$year_t))
climate_lambdas <- bind_cols(Eminus = mean_eminus_yearly[,1], Eplus = mean_eplus_yearly[,1]) %>% 
  bind_cols(climate_data[2:12,])

climate_lambdaplot <- ggplot(data = climate_lambdas) +
  geom_point(aes(x = `Cumulative PPT (mm)`, y = Eplus), size = rel(4), color = "#397BB7") +
  geom_point(aes(x = `Cumulative PPT (mm)`, y = Eminus), size = rel(4), color = "#397BB7", pch = 21) +
  geom_smooth(aes(x = `Cumulative PPT (mm)`, y = Eplus), linetype = "solid", color = "#397BB7", method = "glm", se = FALSE) +
  geom_smooth(aes(x = `Cumulative PPT (mm)`, y = Eminus), linetype = "dashed", color = "#397BB7", method = "glm", se = FALSE) +
  ylab(expression(paste("Annual", " ", lambda))) + scale_linetype_manual("",breaks = c("E+","E-"), values = c("solid", "dashed")) + 
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(size = .5, colour = "black"),
        axis.line.y = element_line(size = .5, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = rel(2), face = "bold"),
        axis.text.x = element_text(size = rel(2), face = "bold"),
        axis.title = element_text(size = rel(2)))


climate_lambdaplot
ggsave(climate_lambdaplot, filename = "~/Documents/climate_lambdaplot.tiff",width = 4, height = 4, bg = "transparent")


# Get surv intercepts
eminus_surv_int <- eplus_surv_int <- matrix(NA, 11,length(post_draws))
for(y in 1:11){
  for(i in 1:length(post_draws)){
eminus_surv_int[y,i] <- invlogit(make_params(species=4,
                        endo_mean=0,
                        endo_var=0,
                        original = 1, # should be =1 to represent recruit
                        draw=post_draws[i],
                        max_size=max_size,
                        rfx=T,
                        year=y,surv_par=surv_par,
                        surv_sdlg_par = surv_sdlg_par,
                        grow_par=grow_par,
                        grow_sdlg_par = grow_sdlg_par,
                        flow_par=flow_par,
                        fert_par=fert_par,
                        spike_par=spike_par,
                        seed_par=seed_par,
                        recruit_par=recruit_par)$surv_int)
eplus_surv_int[y,i] <- invlogit(make_params(species=4,
                               endo_mean=1,
                               endo_var=1,
                               original = 1, # should be =1 to represent recruit
                               draw=post_draws[i],
                               max_size=max_size,
                               rfx=T,
                               year=y,surv_par=surv_par,
                               surv_sdlg_par = surv_sdlg_par,
                               grow_par=grow_par,
                               grow_sdlg_par = grow_sdlg_par,
                               flow_par=flow_par,
                               fert_par=fert_par,
                               spike_par=spike_par,
                               seed_par=seed_par,
                               recruit_par=recruit_par)$surv_int)
}}

mean_surv_eplus_yearly <- matrix(NA, 11,7)
mean_surv_eminus_yearly <- matrix(NA, 11,7)

for(y in 1:11){
  mean_surv_eplus_yearly[y,1] <- mean(eplus_surv_int[y,])
  mean_surv_eplus_yearly[y,2:7] <- quantile(eplus_surv_int[y,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95), na.rm = T)
  
  mean_surv_eminus_yearly[y,1] <- mean(eminus_surv_int[y,])
  mean_surv_eplus_yearly[y,2:7] <- quantile(eminus_surv_int[y,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95), na.rm = T)
  
}

climate_surv <- bind_cols(Eminus = mean_surv_eminus_yearly[,1], Eplus = mean_surv_eplus_yearly[,1]) %>%
  bind_cols(climate_data[2:12,])

  
  
  
climate_survplot <- ggplot(data = climate_surv) +
  geom_point(aes(x = `Cumulative PPT (mm)`, y = Eplus), color = "red") +
  geom_point(aes(x = `Cumulative PPT (mm)`, y = Eminus), color = "darkred") +
  geom_smooth(aes(x = `Cumulative PPT (mm)`, y = Eplus), color = "red", method = "glm", se = FALSE) +
  geom_smooth(aes(x = `Cumulative PPT (mm)`, y = Eminus), color = "darkred", method = "glm", se = FALSE)
  # geom_point(aes(x = `Mean Temp. (C˚)`, y = Eplus), color = "red") +
  # geom_point(aes(x = `Mean Temp. (C˚)`, y = Eminus), color = "darkred") +
  # geom_smooth(aes(x = `Mean Temp. (C˚)`, y = Eplus), color = "red", method = "glm", se = FALSE) +
  # geom_smooth(aes(x = `Mean Temp. (C˚)`, y = Eminus), color = "darkred", method = "glm", se = FALSE)

climate_survplot
ggsave(climate_survplot, filename = "~/Documents/climate_survplot.tiff",width = 6, height = 4.5, bg = "transparent")


  
  
  
######## Making vital rate  yearly figures #####
params_yearly <- array(dim = c(17,7,2,11,length(post_draws)))
for(i in 1:length(post_draws)){
  for(y in 1:11){
    for(e in 1:2){
      for(s in 1:7){
        params_yearly[,s,e,y,i] <- unlist(make_params(species=s,
                                              endo_mean=(e-1),
                                              endo_var=(e-1),
                                              original = 1, # should be =1 to represent recruit
                                              draw=post_draws[i],
                                              max_size=max_size,
                                              rfx=T,
                                              year=y,
                                              surv_par=surv_par,
                                              surv_sdlg_par = surv_sdlg_par,
                                              grow_par=grow_par,
                                              grow_sdlg_par = grow_sdlg_par,
                                              flow_par=flow_par,
                                              fert_par=fert_par,
                                              spike_par=spike_par,
                                              seed_par=seed_par,
                                              recruit_par=recruit_par))
      }
    }
  }
}
  
mean_params_yearly <- array(dim = c(17,7,2,11))
for(p in 1:17){
  for(s in 1:7){
    for(e in 1:2){
      for(y in 1:11){
      mean_params_yearly[p,s,e,y] <- mean(params_yearly[p,s,e,y,])
      }
    }
  }
}
x_seq <- array(dim = c(100,7), dimnames = list(size_t = paste0("size", 1:100), Species = c("AGPE", "ELRI", "ELVI", "FESU", "LOAR", "POAL", "POSY")))
for(s in 1:7){
x_seq[,s] <- seq(length = 100, from = 1, to = max_size$actual_max_size[s])
}
x_seq_df <- as_tibble(x_seq, rownames = "row_seq") %>% 
  pivot_longer(-row_seq, names_to = c("Species"), values_to = c("x_seq"))

surv_mod <- grow_mod <- flw_mod <- fert_mod <- array(dim = c(100,11,2,7), dimnames = list(size_t = paste0("size", 1:100), Year = paste0("y", 1:11), Endo = paste0("e",1:2), Species = c("AGPE", "ELRI", "ELVI", "FESU", "LOAR", "POAL", "POSY")))
for(s in 1:7){
  for(e in 1:2){
    for(y in 1:11){
surv_mod[,y,e,s] <- invlogit(mean_params_yearly[1,s,e,y] + mean_params_yearly[2,s,e,y]*log(x_seq[,s]))
grow_mod[,y,e,s] <- exp(mean_params_yearly[4,s,e,y] + mean_params_yearly[5,s,e,y]*log(x_seq[,s]))
flw_mod[,y,e,s] <- invlogit(mean_params_yearly[9,s,e,y] + mean_params_yearly[10,s,e,y]*log(x_seq[,s]))
fert_mod[,y,e,s] <- exp(mean_params_yearly[11,s,e,y] + mean_params_yearly[12,s,e,y]*log(x_seq[,s]))
    }
  }
}
surv_fit_df <- as_tibble(surv_mod, rownames = "row_seq") %>% 
  pivot_longer( -row_seq, names_to = c("Year", "Endo", "Species"), names_pattern = ("([^.]+).([^.]+).([^.]+)"), values_to = "prob_surv_t1") %>% 
  left_join(x_seq_df, by = c("Species", "row_seq")) %>% 
  mutate(Endo = case_when(Endo == "e1" ~ 0,
                          Endo == "e2" ~ 1))
grow_fit_df <- as_tibble(grow_mod, rownames = "row_seq") %>% 
  pivot_longer( -row_seq, names_to = c("Year", "Endo", "Species"), names_pattern = ("([^.]+).([^.]+).([^.]+)"), values_to = "size_t1") %>% 
  left_join(x_seq_df, by = c("Species", "row_seq")) %>% 
  mutate(Endo = case_when(Endo == "e1" ~ 0,
                          Endo == "e2" ~ 1))

flw_fit_df <- as_tibble(flw_mod, rownames = "row_seq") %>% 
  pivot_longer( -row_seq, names_to = c("Year", "Endo", "Species"), names_pattern = ("([^.]+).([^.]+).([^.]+)"), values_to = "prob_flw_t") %>% 
  left_join(x_seq_df, by = c("Species", "row_seq")) %>% 
  mutate(Endo = case_when(Endo == "e1" ~ 0,
                          Endo == "e2" ~ 1))
fert_fit_df <- as_tibble(fert_mod, rownames = "row_seq") %>% 
  pivot_longer( -row_seq, names_to = c("Year", "Endo", "Species"), names_pattern = ("([^.]+).([^.]+).([^.]+)"), values_to = "no_flw_t") %>% 
  left_join(x_seq_df, by = c("Species", "row_seq")) %>% 
  mutate(Endo = case_when(Endo == "e1" ~ 0,
                          Endo == "e2" ~ 1))
# Bin Data by year

bin_by_size <- function(df_raw, nbins){
  require(tidyverse)
  df_raw <- ungroup(df_raw)
  surv_bin_byyear <- df_raw %>% 
    filter(!is.na(surv_t1)) %>% 
    filter(!is.na(logsize_t)) %>% 
    filter(!is.na(endo_01)) %>%   # There are a few LOAR that don't have a plot level endo assigned
    filter(origin_01 == 1 & year_t != birth | origin_01 == 0) %>% 
    dplyr::select(surv_t1, logsize_t, year_t, endo_01, species) %>% 
    rename(Endo = endo_01, Year = year_t, Species = species) %>%
    mutate(size_bin = cut(logsize_t, breaks = nbins)) %>%
    group_by(size_bin, Year, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_surv = mean(surv_t1,na.rm=T),
              samplesize = n())
    
  surv_sdlg_bin_byyear <-  df_raw %>% 
    filter(!is.na(surv_t1)) %>%
    filter(!is.na(logsize_t)) %>% 
    filter(!is.na(endo_01)) %>%  # There are a few LOAR that don't have a plot level endo assigned
    filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
    filter(logsize_t == 0) %>% 
    dplyr::select(surv_t1, logsize_t, year_t, endo_01, species) %>%
    rename(Endo = endo_01, Year = year_t, Species = species) %>%
    group_by(logsize_t, Year, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_sdlg_surv = mean(surv_t1,na.rm=T),
              samplesize = n())
    
  grow_bin_byyear <-   df_raw %>% 
    filter(!is.na(logsize_t)) %>% 
    filter(!is.na(size_t1)) %>% 
    filter(!is.na(endo_01)) %>% 
    filter(origin_01 == 1 & year_t != birth | origin_01 == 0) %>% 
    dplyr::select(size_t1, logsize_t, year_t, endo_01, species) %>%
    rename(Endo = endo_01, Year = year_t, Species = species) %>% 
    mutate(size_bin = cut(logsize_t, breaks = nbins)) %>%
    group_by(size_bin, Year, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_size_t1 = mean(size_t1,na.rm=T),
              samplesize = n())
  
  grow_sdlg_bin_byyear <-df_raw %>% 
    filter(!is.na(logsize_t)) %>% 
    filter(!is.na(size_t1)) %>% 
    filter(!is.na(endo_01)) %>% 
    filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
    filter(logsize_t == 0) %>% 
    dplyr::select(size_t1, logsize_t, year_t, endo_01, species) %>%
    rename(Endo = endo_01,Year = year_t, Species = species) %>%
    group_by(logsize_t, Year, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_size_t1 = mean(size_t1,na.rm=T),
              samplesize = n())
  
  
  flw_bin_byyear <- df_raw %>% 
    filter(!is.na(FLW_STAT_T)) %>% 
    filter(!is.na(logsize_t)) %>% 
    filter(!is.na(endo_01)) %>% 
    dplyr::select(FLW_STAT_T, logsize_t, year_t, endo_01, species) %>%
    rename(Endo = endo_01, Year = year_t, Species = species) %>%
    mutate(size_bin = cut(logsize_t, breaks = nbins)) %>%
    group_by(size_bin, Year, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_flw = mean(FLW_STAT_T,na.rm=T),
              samplesize = n())
  
  fert_bin_byyear <- df_raw %>% 
    filter(!is.na(FLW_COUNT_T)) %>% 
    filter(FLW_COUNT_T > 0) %>% 
    filter(!is.na(logsize_t)) %>% 
    dplyr::select(FLW_COUNT_T, logsize_t, year_t, endo_01, species) %>%
    rename(Endo = endo_01,Year = year_t, Species = species) %>%
    mutate(size_bin = cut(logsize_t, breaks = nbins)) %>%
    group_by(size_bin, Year, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_fert = mean(FLW_COUNT_T,na.rm=T),
              samplesize = n())

    size_bins <- list(surv_bin_byyear = surv_bin_byyear, 
                      surv_sdlg_bin_byyear = surv_sdlg_bin_byyear,
                      grow_bin_byyear = grow_bin_byyear, 
                      grow_sldg_bin_byyear = grow_sdlg_bin_byyear, 
                      flw_bin_byyear = flw_bin_byyear, 
                      fert_bin_byyear = fert_bin_byyear)
  return(size_bins)
}
size_bin_data <- bin_by_size(LTREB_full, nbins = 10)

# surv plot
surv_dataplot <- ggplot(data = surv_fit_df) +
  geom_line( aes(x = log(x_seq), y = prob_surv_t1, group = Year, color = Species, linetype = as.factor(Endo))) +
  geom_point(data = size_bin_data$surv_bin_byyear, aes(x = mean_size, y = mean_surv, color = as.factor(Species), size = samplesize, shape = as.factor(Endo))) + 
  facet_wrap(~Species + Endo, ncol = 2, scales = "free")  +
  scale_fill_manual(values = c("#59A1BC", "#4E816D", "#EFAD3A", "#397BB7", "#9E78A1", "#E04D55", "#9D5251")) +
  scale_color_manual(values = c("#59A1BC", "#4E816D", "#EFAD3A", "#397BB7", "#9E78A1", "#E04D55", "#9D5251")) +
  scale_shape_manual(values = c(19, 1),limits = factor(c(1,0)), name = "Endo Status", labels = c("E-", "E+"))+ 
  scale_linetype_discrete(name = "Endo Status",limits = factor(c(1,0)), labels = c("E+", "E-"))+ 
  ggtitle("Yearly Prob. of Survival")+ xlab("log(size_t)") + ylab("Prob. Surv.") +
  guides(color = FALSE, size = FALSE, shape = FALSE, linetype = FALSE)+
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        # panel.grid.major.y = element_line(size = .5, colour = "gray"),
        # panel.grid.major.x = element_line(size = .5, colour = "gray"),
        # panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(size = .5, colour = "black"),
        axis.line.y = element_line(size = .5, colour = "black"),
        # axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size = rel(2), face = "bold"),
        # axis.text.x = element_text(size = rel(2), face = "bold"),
        axis.title = element_text(size = rel(2)),
        strip.background = element_rect(fill = "white", color = "lightgray"))

surv_dataplot
ggsave(surv_dataplot, filename = "~/Documents/surv_dataplot.tiff",width = 6, height = 20, bg = "white")

FESUsurv_dataplot <- ggplot(data = subset(surv_fit_df, Species == "FESU")) +
  geom_line( aes(x = log(x_seq), y = prob_surv_t1, group = Year, color = Species, linetype = as.factor(Endo))) +
  geom_point(data = subset(size_bin_data$surv_bin_byyear, Species == "FESU"), aes(x = mean_size, y = mean_surv, color = as.factor(Species), size = samplesize, shape = as.factor(Endo))) + 
  facet_wrap(~Species + Endo, ncol = 2, scales = "free")  +
  scale_fill_manual(values = c( "#397BB7")) +
  scale_color_manual(values = c( "#397BB7")) +
  scale_shape_manual(values = c(19, 1),limits = factor(c(1,0)), name = "Endo Status", labels = c("E+", "E-"))+ 
  scale_linetype_discrete(name = "Endo Status",limits = factor(c(1,0)), labels = c("E+", "E-"))+ 
  ggtitle("Yearly Prob. of Survival")+ xlab("log(size_t)") + ylab("Prob. Surv.") +
  guides(color = FALSE, size = FALSE, shape = FALSE, linetype = FALSE)+
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        # panel.grid.major.y = element_line(size = .5, colour = "gray"),
        # panel.grid.major.x = element_line(size = .5, colour = "gray"),
        # panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(size = .5, colour = "black"),
        axis.line.y = element_line(size = .5, colour = "black"),
        # axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size = rel(2), face = "bold"),
        # axis.text.x = element_text(size = rel(2), face = "bold"),
        axis.title = element_text(size = rel(2)),
        strip.background = element_rect(fill = "white", color = "lightgray"))

FESUsurv_dataplot
ggsave(FESUsurv_dataplot, filename = "~/Documents/FESUsurv_dataplot.png",width = 5, height = 4, bg = "white")

# growth plot
grow_dataplot <- ggplot(data = grow_fit_df) +
  geom_line( aes(x = log(x_seq), y = size_t1, group = Year, color = Species, linetype = as.factor(Endo))) +
  geom_point(data = size_bin_data$grow_bin_byyear, aes(x = mean_size, y = mean_size_t1, color = as.factor(Species), size = samplesize, shape = as.factor(Endo))) +
  facet_wrap(~Species + Endo, ncol = 2, scales = "free")  +
  scale_fill_manual(values = c("#59A1BC", "#4E816D", "#EFAD3A", "#397BB7", "#9E78A1", "#E04D55", "#9D5251")) +
  scale_color_manual(values = c("#59A1BC", "#4E816D", "#EFAD3A", "#397BB7", "#9E78A1", "#E04D55", "#9D5251")) +
  scale_shape_manual(values = c(19, 1), name = "Endo Status", labels = c("E+", "E-"))+ 
  scale_linetype_discrete(name = "Endo Status",limits = factor(c(1,0)), labels = c("E+", "E-"))+ 
  ggtitle("Yearly Growth")+ xlab("log(size_t)") + ylab("Size_t1") +
  guides(color = FALSE, size = FALSE, shape = FALSE, linetype = FALSE)+
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        # panel.grid.major.y = element_line(size = .5, colour = "gray"),
        # panel.grid.major.x = element_line(size = .5, colour = "gray"),
        # panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(size = .5, colour = "black"),
        axis.line.y = element_line(size = .5, colour = "black"),
        # axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size = rel(2), face = "bold"),
        # axis.text.x = element_text(size = rel(2), face = "bold"),
        axis.title = element_text(size = rel(2)),
        strip.background = element_rect(fill = "white", color = "lightgray"))

grow_dataplot
ggsave(grow_dataplot, filename = "~/Documents/grow_dataplot.tiff",width = 6, height = 20, bg = "white")


# flw plot
flw_dataplot <- ggplot(data = flw_fit_df) +
  geom_line( aes(x = log(x_seq), y = prob_flw_t, group = Year, color = Species, linetype = as.factor(Endo))) +
  geom_point(data = size_bin_data$flw_bin_byyear, aes(x = mean_size, y = mean_flw, color = as.factor(Species), size = samplesize, shape = as.factor(Endo))) + 
  facet_wrap(~Species + Endo, ncol = 2, scales = "free")  +
  scale_fill_manual(values = c("#59A1BC", "#4E816D", "#EFAD3A", "#397BB7", "#9E78A1", "#E04D55", "#9D5251")) +
  scale_color_manual(values = c("#59A1BC", "#4E816D", "#EFAD3A", "#397BB7", "#9E78A1", "#E04D55", "#9D5251")) +
  scale_shape_manual(values = c(19, 1),limits = factor(c(1,0)), name = "Endo Status", labels = c("E+", "E-"))+ 
  scale_linetype_discrete(name = "Endo Status",limits = factor(c(1,0)), labels = c("E+", "E-"))+ 
  ggtitle("Yearly Prob. of Flowering")+ xlab("log(size_t)") + ylab("Prob. Flw.") +
  guides(color = FALSE, size = FALSE, shape = FALSE, linetype = FALSE)+
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        # panel.grid.major.y = element_line(size = .5, colour = "gray"),
        # panel.grid.major.x = element_line(size = .5, colour = "gray"),
        # panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(size = .5, colour = "black"),
        axis.line.y = element_line(size = .5, colour = "black"),
        # axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size = rel(2), face = "bold"),
        # axis.text.x = element_text(size = rel(2), face = "bold"),
        axis.title = element_text(size = rel(2)),
        strip.background = element_rect(fill = "white", color = "lightgray"))

flw_dataplot
ggsave(flw_dataplot, filename = "~/Documents/flw_dataplot.tiff",width = 6, height = 20, bg = "white")

# Fert plot
fert_dataplot <- ggplot(data = fert_fit_df) +
  geom_line( aes(x = log(x_seq), y = no_flw_t, group = Year, color = Species, linetype = as.factor(Endo))) +
  geom_point(data = size_bin_data$fert_bin_byyear, aes(x = mean_size, y = mean_fert, color = as.factor(Species), size = samplesize, shape = as.factor(Endo))) + 
  facet_wrap(~Species + Endo, ncol = 2, scales = "free")  +
  scale_fill_manual(values = c("#59A1BC", "#4E816D", "#EFAD3A", "#397BB7", "#9E78A1", "#E04D55", "#9D5251")) +
  scale_color_manual(values = c("#59A1BC", "#4E816D", "#EFAD3A", "#397BB7", "#9E78A1", "#E04D55", "#9D5251")) +
  scale_shape_manual(values = c(19, 1),limits = factor(c(1,0)), name = "Endo Status", labels = c("E+", "E-"))+ 
  scale_linetype_discrete(name = "Endo Status",limits = factor(c(1,0)), labels = c("E+", "E-"))+ 
  ggtitle("Yearly Flowering No.")+ xlab("log(size_t)") + ylab("No. Flw Tillers") +
  guides()+
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        # panel.grid.major.y = element_line(size = .5, colour = "gray"),
        # panel.grid.major.x = element_line(size = .5, colour = "gray"),
        # panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(size = .5, colour = "black"),
        axis.line.y = element_line(size = .5, colour = "black"),
        # axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size = rel(2), face = "bold"),
        # axis.text.x = element_text(size = rel(2), face = "bold"),
        axis.title = element_text(size = rel(2)),
        strip.background = element_rect(fill = "white", color = "lightgray"))

legend <- cowplot::get_legend(fert_dataplot + theme(legend.title = element_text(size = rel(2)),
                                                    legend.text = element_text(size = rel(1.5))))
grid.newpage()
grid.draw(legend)  

fert_dataplot <- fert_dataplot + guides(size = FALSE, color = FALSE, linetype = FALSE, shape = FALSE)
ggsave(fert_dataplot, filename = "~/Documents/fert_dataplot.tiff",width = 6, height = 20, bg = "white")



vr_plot <- grid.arrange(surv_dataplot, grow_dataplot, flw_dataplot, fert_dataplot, legend, nrow = 1)       
ggsave(vr_plot, filename = "~/Documents/vr_plot.tiff",width = 20, height = 20, bg = "white")



##### Making a mean endophyte effect vital rate plot ####

params_mean <- array(dim = c(17,7,2,length(post_draws)))
for(i in 1:length(post_draws)){
    for(e in 1:2){
      for(s in 1:7){
        params_mean[,s,e,i] <- unlist(make_params(species=s,
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
                                                      recruit_par=recruit_par))
      }
    }
  }


mean_params_mean <- array(dim = c(17,7,2, 7))
for(p in 1:17){
  for(s in 1:7){
    for(e in 1:2){
        mean_params_mean[p,s,e,1] <- mean(params_mean[p,s,e,])
        mean_params_mean[p,s,e,2:7] <- quantile(params_mean[p,s,e,], c(0.05,0.125,0.25,0.75,0.875,0.95))
      }
    }
  }

x_seq <- array(dim = c(100,7), dimnames = list(size_t = paste0("size", 1:100), Species = c("AGPE", "ELRI", "ELVI", "FESU", "LOAR", "POAL", "POSY")))
for(s in 1:7){
  x_seq[,s] <- seq(length = 100, from = 1, to = max_size$actual_max_size[s])
}
x_seq_df <- as_tibble(x_seq, rownames = "row_seq") %>% 
  pivot_longer(-row_seq, names_to = c("Species"), values_to = c("x_seq"))

surv_mod <- grow_mod <- flw_mod <- fert_mod <- array(dim = c(100,2,7,7), dimnames = list(size_t = paste0("size", 1:100), Endo = paste0("e",1:2), Species = c("AGPE", "ELRI", "ELVI", "FESU", "LOAR", "POAL", "POSY"), iteration = c("mean", "five","twelvepointfive", "twentyfive", "seventyfive", "eightysevenpointfive", "ninetyfive")))
for(s in 1:7){
  for(e in 1:2){
    for(i in 1:7)
      surv_mod[,e,s,i] <- invlogit(mean_params_mean[1,s,e,i] + mean_params_mean[2,s,e,i]*log(x_seq[,s]))
      grow_mod[,e,s,i] <- exp(mean_params_mean[4,s,e,i] + mean_params_mean[5,s,e,i]*log(x_seq[,s]))
      flw_mod[,e,s,i] <- invlogit(mean_params_mean[9,s,e,i] + mean_params_mean[10,s,e,i]*log(x_seq[,s]))
      fert_mod[,e,s,i] <- exp(mean_params_mean[11,s,e,i] + mean_params_mean[12,s,e,i]*log(x_seq[,s]))
    }
  }

surv_fit_df <- as_tibble(surv_mod, rownames = "row_seq") %>% 
  pivot_longer( -row_seq, names_to = c( "Endo", "Species", "Iterations"), names_pattern = ("([^.]+).([^.]+).([^.]+)"), values_to = "prob_surv_t1") %>% 
  left_join(x_seq_df, by = c("Species", "row_seq")) %>% 
  mutate(Endo = case_when(Endo == "e1" ~ 0,
                          Endo == "e2" ~ 1)) %>% 
  pivot_wider(names_from = "Iterations", values_from = prob_surv_t1)

grow_fit_df <- as_tibble(grow_mod, rownames = "row_seq") %>% 
  pivot_longer( -row_seq, names_to = c( "Endo", "Species", "Iterations"), names_pattern = ("([^.]+).([^.]+).([^.]+)"), values_to = "prob_surv_t1") %>% 
  left_join(x_seq_df, by = c("Species", "row_seq")) %>% 
  mutate(Endo = case_when(Endo == "e1" ~ 0,
                          Endo == "e2" ~ 1)) %>% 
  pivot_wider(names_from = "Iterations", values_from = prob_surv_t1)


flw_fit_df <- as_tibble(flw_mod, rownames = "row_seq") %>% 
  pivot_longer( -row_seq, names_to = c( "Endo", "Species", "Iterations"), names_pattern = ("([^.]+).([^.]+).([^.]+)"), values_to = "prob_surv_t1") %>% 
  left_join(x_seq_df, by = c("Species", "row_seq")) %>% 
  mutate(Endo = case_when(Endo == "e1" ~ 0,
                          Endo == "e2" ~ 1)) %>% 
  pivot_wider(names_from = "Iterations", values_from = prob_surv_t1)

fert_fit_df <- as_tibble(fert_mod, rownames = "row_seq") %>% 
  pivot_longer( -row_seq, names_to = c( "Endo", "Species", "Iterations"), names_pattern = ("([^.]+).([^.]+).([^.]+)"), values_to = "prob_surv_t1") %>% 
  left_join(x_seq_df, by = c("Species", "row_seq")) %>% 
  mutate(Endo = case_when(Endo == "e1" ~ 0,
                          Endo == "e2" ~ 1)) %>% 
  pivot_wider(names_from = "Iterations", values_from = prob_surv_t1)

  
  
# Bin Data by Endo

bin_by_size_mean <- function(df_raw, nbins){
  require(tidyverse)
  df_raw <- ungroup(df_raw)
  surv_bin_mean <- df_raw %>% 
    filter(!is.na(surv_t1)) %>% 
    filter(!is.na(logsize_t)) %>% 
    filter(!is.na(endo_01)) %>%   # There are a few LOAR that don't have a plot level endo assigned
    filter(origin_01 == 1 & year_t != birth | origin_01 == 0) %>% 
    dplyr::select(surv_t1, logsize_t, year_t, endo_01, species) %>% 
    rename(Endo = endo_01, Year = year_t, Species = species) %>%
    mutate(size_bin = cut(logsize_t, breaks = nbins)) %>%
    group_by(size_bin, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_surv = mean(surv_t1,na.rm=T),
              samplesize = n())
  
  surv_sdlg_bin_mean <-  df_raw %>% 
    filter(!is.na(surv_t1)) %>%
    filter(!is.na(logsize_t)) %>% 
    filter(!is.na(endo_01)) %>%  # There are a few LOAR that don't have a plot level endo assigned
    filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
    filter(logsize_t == 0) %>% 
    dplyr::select(surv_t1, logsize_t, year_t, endo_01, species) %>%
    rename(Endo = endo_01, Year = year_t, Species = species) %>%
    group_by(logsize_t, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_sdlg_surv = mean(surv_t1,na.rm=T),
              samplesize = n())
  
  grow_bin_mean <-   df_raw %>% 
    filter(!is.na(logsize_t)) %>% 
    filter(!is.na(size_t1)) %>% 
    filter(!is.na(endo_01)) %>% 
    filter(origin_01 == 1 & year_t != birth | origin_01 == 0) %>% 
    dplyr::select(size_t1, logsize_t, year_t, endo_01, species) %>%
    rename(Endo = endo_01, Year = year_t, Species = species) %>% 
    mutate(size_bin = cut(logsize_t, breaks = nbins)) %>%
    group_by(size_bin, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_size_t1 = mean(size_t1,na.rm=T),
              samplesize = n())
  
  grow_sdlg_bin_mean <-df_raw %>% 
    filter(!is.na(logsize_t)) %>% 
    filter(!is.na(size_t1)) %>% 
    filter(!is.na(endo_01)) %>% 
    filter(origin_01 == 1 & year_t == birth) %>%  #filtering for recruits that are just germinated
    filter(logsize_t == 0) %>% 
    dplyr::select(size_t1, logsize_t, year_t, endo_01, species) %>%
    rename(Endo = endo_01,Year = year_t, Species = species) %>%
    group_by(logsize_t, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_size_t1 = mean(size_t1,na.rm=T),
              samplesize = n())
  
  
  flw_bin_mean <- df_raw %>% 
    filter(!is.na(FLW_STAT_T)) %>% 
    filter(!is.na(logsize_t)) %>% 
    filter(!is.na(endo_01)) %>% 
    dplyr::select(FLW_STAT_T, logsize_t, year_t, endo_01, species) %>%
    rename(Endo = endo_01, Year = year_t, Species = species) %>%
    mutate(size_bin = cut(logsize_t, breaks = nbins)) %>%
    group_by(size_bin, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_flw = mean(FLW_STAT_T,na.rm=T),
              samplesize = n())
  
  fert_bin_mean <- df_raw %>% 
    filter(!is.na(FLW_COUNT_T)) %>% 
    filter(FLW_COUNT_T > 0) %>% 
    filter(!is.na(logsize_t)) %>% 
    dplyr::select(FLW_COUNT_T, logsize_t, year_t, endo_01, species) %>%
    rename(Endo = endo_01,Year = year_t, Species = species) %>%
    mutate(size_bin = cut(logsize_t, breaks = nbins)) %>%
    group_by(size_bin, Endo, Species) %>%
    summarise(mean_size = mean((logsize_t),na.rm=T),
              mean_fert = mean(FLW_COUNT_T,na.rm=T),
              samplesize = n())
  
  size_bins <- list(surv_bin_mean = surv_bin_mean, 
                    surv_sdlg_bin_mean = surv_sdlg_bin_mean,
                    grow_bin_mean = grow_bin_mean, 
                    grow_sldg_bin_mean = grow_sdlg_bin_mean, 
                    flw_bin_mean = flw_bin_mean, 
                    fert_bin_mean = fert_bin_mean)
  return(size_bins)
}
size_bin_data_mean <- bin_by_size_mean(LTREB_full, nbins = 10)


FESUsurv_dataplot_mean <- ggplot(data = subset(surv_fit_df, Species == "FESU")) +
  geom_ribbon(aes(x = log(x_seq), ymin = five , ymax = ninetyfive, fill = as.factor(Endo), linetype = as.factor(Endo)), color = c("gray80"), alpha = .2) +
  geom_line( aes(x = log(x_seq), y = mean, color = Species, linetype = as.factor(Endo), lwd = 3)) +
  geom_point(data = subset(size_bin_data_mean$surv_bin_mean, Species == "FESU"), aes(x = mean_size, y = mean_surv, color = as.factor(Species), size = samplesize, shape = as.factor(Endo))) + 
  facet_wrap(~Species, ncol = 1, scales = "free")  +
  scale_fill_manual(values = c( "gray70", "gray50")) +
  scale_color_manual(values = c( "#397BB7")) +
  scale_shape_manual(values = c(19, 1),limits = factor(c(1,0)), name = "Endo Status", labels = c("E-", "E+"))+ 
  scale_linetype_discrete(name = "Endo Status",limits = factor(c(1,0)), labels = c("E-", "E+"))+ 
  ggtitle("Mean Survival")+ xlab("log(size_t)") + ylab("Prob. Surv.") +
  guides(color = FALSE, size = FALSE, shape = FALSE, linetype = FALSE, fill = FALSE)+
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        # panel.grid.major.y = element_line(size = .5, colour = "gray"),
        # panel.grid.major.x = element_line(size = .5, colour = "gray"),
        # panel.grid.minor.x = element_blank(),
        axis.line.x = element_line(size = .5, colour = "black"),
        axis.line.y = element_line(size = .5, colour = "black"),
        # axis.ticks.y = element_blank(),
        # axis.text.y = element_text(size = rel(2), face = "bold"),
        # axis.text.x = element_text(size = rel(2), face = "bold"),
        axis.title = element_text(size = rel(2)),
        strip.background = element_rect(fill = "white", color = "lightgray"))

FESUsurv_dataplot_mean
ggsave(FESUsurv_dataplot_mean, filename = "~/Documents/FESUsurv_dataplot_mean.png",width = 5, height = 4, bg = "white")

#### Some summary calculations for the manuscript ####
# getting the mean values as well as some posterior summaries from the mean effect
summary_lambda_mean <- lambda_mean_df %>% 
  group_by(species, Species) %>% 
  summarize(mean_lambda_diff = mean(lambda_diff),
            percent_endo_comparison = 100-(mean(e1)/mean(e2))*100,
            percent_diffgreaterthanzero = (sum(lambda_diff>0)/n())*100,
            draws = n())
  
summary_lambda_sd <- lambda_var_df %>% 
  group_by(species, Species) %>% 
  summarize(sd_lambda_diff = mean(lambda_diff),
            percent_endo_comparison = 100-(mean(e1)/mean(e2))*100,
            percent_difflessthanzero = (sum(lambda_diff<0)/n())*100,
            draws = n())

summary_lambda_cv <- lambda_cv_df  %>% 
  filter(sampling == "obs") %>%
  group_by(species, Species) %>% 
  summarize(cv_lambda_diff = mean(lambda_diff),
            percent_endo_comparison = 100-(mean(e1)/mean(e2))*100,
            percent_difflessthanzero = (sum(lambda_diff<0)/n())*100,
            draws = n())



