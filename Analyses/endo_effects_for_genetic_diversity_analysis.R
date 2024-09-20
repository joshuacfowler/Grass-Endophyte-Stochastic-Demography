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
                                                       tau_year, tau_plot))
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

lambda_mean <- array(dim = c(8,2,n_draws))
for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      lambda_mean[s,e,i] <- lambda(bigmatrix(make_params_quadXorigin(species=s,
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
                                                                     recruit_par=recruit_par), 
                                             quadratic = 1, # set this to one when using the vital rate models fit with a quadratic size-structure parameter
                                             extension = 100)$MPMmat) # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
      lambda_mean[8,e,i] <- mean(lambda_mean[1:7,e,i])
    }
  }
}
# saving lambda_mean with 500 post draws to dropbox
saveRDS(lambda_mean, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_mean_geneticanalysis.rds")
lambda_mean <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_mean_geneticanalysis.rds")


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


## now do variance in lambda 
# In this instance, we are calculating the effect on variance only for years 2007-2018, which is when samples where taken from plots for genetic diversity measurment

lambda_hold <- array(dim = c(8,7,2,n_draws)) # 8 years because of reproduction measured in year t1; needs to be before growth, so no year 1. This model should encompass transition years 2008-2009 to 2017-2018
lambda_var <- array(dim = c(8,2,n_draws))
for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      for(y in 1:8){
        
        lambda_hold[y,s,e,i] <- lambda(bigmatrix(make_params_quadXorigin(species=s,
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
                                                 quadratic = 1, # set this to one when using the vital rate models fit with a quadratic size-structure parameter
                                                 extension = 100)$MPMmat) # the extension parameter is used to fit the growth kernel to sizes larger than max size without losing probability density
      }
      lambda_var[s,e,i] <- sd(lambda_hold[,s,e,i]) # we calulate the standard deviation here
    }
    lambda_var[8,e,i] <- mean(lambda_var[1:7,e,i])
  }
}
#saving the yearly lambdas and the sd of lambdas to dropbox
# saveRDS(lambda_hold, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_hold.rds")
# saveRDS(lambda_var, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_var.rds")
# lambda_hold <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_hold.rds")
# lambda_var <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_var.rds")

# saveRDS(lambda_hold, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_hold_quadXorigin.rds")
saveRDS(lambda_var, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_var_geneticanalysis.rds")
lambda_var <- read_rds(file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/lambda_var_geneticanalysis.rds")



#Calculationg endophyte effect on sd and variance
lambda_sds <- matrix(NA,8,2)
lambda_vars <- matrix(NA,8,2)

lambda_cv <- (lambda_var)/(lambda_mean) 
# lambda_cv <- (lambda_var^2)/(2*lambda_mean^2) # This is the "variance penalty" not actually CV
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

lambda_mean_diff
lambda_cv_diff


####################################
####### plotting the effects #######
####################################


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


# ggplot(data = lambda_mean_diff_df) +
#   geom_hline(yintercept = 0, col = "black") + 
#   geom_linerange(aes(y = mean, x = species, ymin = 0, ymax = mean), color = "white", lwd =4)+
#   geom_point(aes(y = mean, x = species, color = species), lwd = 4) +
#   geom_linerange(aes(y = mean, x = species, ymin = twentyfifth, ymax = seventyfifth, color = species), lwd = 2) +
#   geom_linerange(aes(y = mean, x = species, ymin = twelfthpointfive, ymax = eightyseventhpointfive, color = species), lwd = 1) +
#   geom_linerange(aes(y = mean, x = species, ymin = fifth, ymax = ninetyfifth, color = species)) +
#   scale_color_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
#   coord_flip()+
#   theme(panel.background = element_rect(fill = "white"),
#         panel.grid = element_line(color = NA))

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
  stat_summary(aes(y = lambda_diff, x = species), fun = median,geom = "point", size = 3) +
  stat_summary(aes(y = lambda_diff, x = species, color = species), fun = median,geom = "point", size = 2) +
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
# ggsave(meanlambda_plot, filename = "meanlambda_plot_quadXorigin.png", width = 4, height = 5)

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
                             Species == "s8" ~ "Species Mean")) %>% 
  filter(lambda_diff>-60) # There is one really weird iteration

lambdavar_plot <- ggplot(data = lambda_cv_df) +
  geom_hline(yintercept = 0, col = "black") + 
  # geom_linerange(data = lambda_var_diff_df, aes(x = species, y = mean, ymin = 0, ymax = mean, color = species))+
  geom_jitter( aes(y = lambda_diff, x = species, color = species), width = .2, alpha = .2) +
  stat_summary(aes(y = lambda_diff, x = species), fun = median,geom = "point", size = 3) +
  stat_summary(aes(y = lambda_diff, x = species, color = species), fun = median,geom = "point", size = 2) +
  scale_color_manual(values = c("#dbdb42", "#b8e3a0", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84", "#A9A9A9")) +
  coord_flip(ylim = c(-1.5,.3)) +  #There's a few iterations out at -1.5
  labs(y = expression(paste("Symbiosis Effect on ", "CV(",lambda,")")),
       x = "",
       color = "Species")+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = NA),
        axis.line = element_line(color = "grey"),
        axis.text.y = element_text(face = "italic"),
        legend.position = "none")
lambdavar_plot
# ggsave(lambdavar_plot, filename = "lambdavar_plot_coefficientofvariation_quadXorigin.png", width = 4, height = 5)

endo_lambdaeffects_plot <-  meanlambda_plot +lambdavar_plot + plot_layout(nrow = 1, guides = "collect")
endo_lambdaeffects_plot
# ggsave(endo_lambdaeffects_plot, filename = "endo_lambdaeffects_plot_quadXorigin.png", width = 12, height = 6)



############################################################### 
#### Creating a dataframe to share with Maya Shamsid-Deen #####
###############################################################
# creating dataframes for sd and for variance

dimnames(lambda_var) <- list(Species = paste0("s",1:8), Endo = paste0("e",1:2), Iteration= paste0("i",1:n_draws))
lambda_sd_cube <- cubelyr::as.tbl_cube(lambda_var)
lambda_sd_df <- as_tibble(lambda_sd_cube) %>% 
  pivot_wider(names_from = Endo, values_from = lambda_var) %>% 
  mutate(lambda_diff = e2-e1) %>% 
  mutate(species = case_when(Species == "s1" ~ "Agrostis perennans",
                             Species == "s2" ~ "Elymus villosus",
                             Species == "s3" ~ "Elymus virginicus",
                             Species == "s4" ~ "Festuca subverticillata",
                             Species == "s5" ~ "Lolium arundinaceum",
                             Species == "s6" ~ "Poa alsodes",
                             Species == "s7" ~ "Poa sylvestris",
                             Species == "s8" ~ "Species Mean"))



lambda_variance <- lambda_var^2

dimnames(lambda_variance) <- list(Species = paste0("s",1:8), Endo = paste0("e",1:2), Iteration= paste0("i",1:n_draws))
lambda_variance_cube <- cubelyr::as.tbl_cube(lambda_variance)
lambda_variance_df <- as_tibble(lambda_variance_cube) %>% 
  pivot_wider(names_from = Endo, values_from = lambda_variance) %>% 
  mutate(lambda_diff = e2-e1) %>% 
  mutate(species = case_when(Species == "s1" ~ "Agrostis perennans",
                             Species == "s2" ~ "Elymus villosus",
                             Species == "s3" ~ "Elymus virginicus",
                             Species == "s4" ~ "Festuca subverticillata",
                             Species == "s5" ~ "Lolium arundinaceum",
                             Species == "s6" ~ "Poa alsodes",
                             Species == "s7" ~ "Poa sylvestris",
                             Species == "s8" ~ "Species Mean"))

# adding a column to keep track of each metric
lambda_mean_df$metric <- "mean"
lambda_cv_df$metric <- "CV"
lambda_sd_df$metric <-  "sd"
lambda_variance_df$metric <-  "variance"

lambda_effects_df <- bind_rows(lambda_mean_df, lambda_cv_df, lambda_sd_df, lambda_variance_df)


endo_effects_summary <- lambda_effects_df %>% 
  group_by(species, metric) %>% 
  dplyr::summarise(effect.average = mean(lambda_diff),
                   effect.sd = sd(lambda_diff),
                   effect.025 = quantile(lambda_diff, probs = .025),
                   effect.975 = quantile(lambda_diff, probs = .975)) %>% 
  filter(species != "Species Mean")

endo_effects_plot <- endo_effects_summary %>% 
  mutate(mean_var = case_when(metric == "mean" ~ "mean",
                              metric != "mean" ~ "var"))

ggplot(endo_effects_plot)+
  geom_bar(aes(y = effect.average, x = species, fill = metric),position = "dodge", stat = "identity") +
  facet_wrap(~mean_var)



endo_effects_summary_forMaya <- endo_effects_summary %>% 
  filter(metric %in% c("mean", "variance", "CV")) %>% 
  arrange( metric, species)


write_csv(endo_effects_summary_forMaya, file = "~/Dropbox/EndodemogData/Model_Runs/MPM_output/endo_effects_to_2018.csv")







################################################################################################
##### Looking at the proportion of surviving original plants vs recruitment in each plot #######
################################################################################################
endo_colors <- c("#dbdb42", "#0c2c84")

LTREB_plot_cleaned <- LTREB_full %>% 
  filter(!is.na(surv_t1), year_t1<=2018) %>% 
  mutate(endo_01 = case_when(plot_fixed %in% c(33,34,39,40) ~ 1,
                             plot_fixed %in% c(32,35) ~ 0, TRUE ~ endo_01)) %>%  # fixing a few plants that were labeled as E+ in E+ plots and vice versa for LOAR
  mutate(Endo = case_when(endo_01 == 0 ~ "E minus", endo_01 == 1 ~ "E plus",
                          is.na(endo_01) ~ "E plus"),
         Origin = case_when(origin_01 == 0 ~ "Original", origin_01 == 1 ~ "Recruit")) 

LTREB_plot_summary <- LTREB_plot_cleaned %>% 
  group_by(species, Endo, plot_fixed, year_t1) %>% 
  summarize(num_original = sum(origin_01==0),
            num_recruit = sum(origin_01==1),
            num_total = sum(!is.na(origin_01)),
            prop_original = num_original/num_total,
            prop_recruit = num_recruit/num_total)


originalcount_plot <- ggplot(LTREB_plot_summary)+
  geom_point(aes(x = year_t1, y = num_original, color = Endo))+
  scale_color_manual(values = c(endo_colors))+
  facet_wrap(species ~ plot_fixed) + theme_light()
# originalcount_plot
ggsave(originalcount_plot, filename = "GeneticDiversityPlots/originalcount_plot.png", width = 10, height = 12)



originaltrends_plot <- ggplot(LTREB_plot_summary)+
  geom_line(aes(x = year_t1, y = num_original, color = Endo, group = plot_fixed))+
  scale_color_manual(values = c(endo_colors))+
  facet_wrap(~species) + theme_light()
# originaltrends_plot
ggsave(originaltrends_plot, filename = "GeneticDiversityPlots/originaltrends_plot.png", width = 10, height = 12)




recruitcount_plot <- ggplot(LTREB_plot_summary)+
  geom_point(aes(x = year_t1, y = num_recruit, color = Endo))+
  scale_color_manual(values = c(endo_colors))+
  facet_wrap(species ~ plot_fixed) + theme_light()
# recruitcount_plot
ggsave(recruitcount_plot, filename = "GeneticDiversityPlots/recruitcount_plot.png", width = 10, height = 12)


recruittrends_plot <- ggplot(LTREB_plot_summary)+
  geom_line(aes(x = year_t1, y = num_recruit, color = Endo, group = plot_fixed))+
  scale_color_manual(values = c(endo_colors))+
  facet_wrap(~species) + theme_light()
# recruittrends_plot
ggsave(recruittrends_plot, filename = "GeneticDiversityPlots/recruittrends_plot.png", width = 10, height = 12)


turnoverprop_plot <- ggplot(LTREB_plot_cleaned)+
  geom_bar(aes(x = year_t1, fill = Origin), position = "fill")+
  facet_wrap(species ~ plot_fixed+Endo) + theme_light()
# turnoverprop_plot
ggsave(turnoverprop_plot, filename = "GeneticDiversityPlots/turnoverprop_plot.png", width = 10, height = 12)


LTREB_plot_summary_2018 <- LTREB_plot_summary %>% 
  filter(year_t1 == 2018)


write_csv(LTREB_plot_summary, "LTREB_plantorigin_summary.csv")
