## Title: Grass endophyte population model with a bayesian framework
## Purpose: functions for building matrix population model from vital rate estimates 
## Authors: Joshua and Tom
#############################################################

invlogit<-function(x){exp(x)/(1+exp(x))}

# Parameter assembly function ---------------------------------------------
make_params <- function(species,endo_mean,endo_var,LTRE=F, endo_mean_U, endo_mean_F, endo_var_U, endo_var_F, original=0,draw,rfx=F,spei=F,year=NULL,repro_offset = 1, max_size,samp=F,samp_extreme=0,
                        surv_par,surv_sdlg_par,grow_par,grow_sdlg_par,flow_par,fert_par,spike_par,seed_par,recruit_par){
  if(LTRE==F){
    endo_mean_U <- endo_mean_F <- endo_mean
    endo_var_U <- endo_var_F <- endo_var
  }
  if(LTRE==T){
    endo_mean_U <- endo_mean_U
    endo_var_U <- endo_var_U
    endo_mean_F <- endo_var_F
  }
  
  if(rfx==F){rfx_surv <- rfx_surv_sdlg <- rfx_grow <- rfx_grow_sdlg <- rfx_flow <- rfx_fert <- rfx_spike <- rfx_rct <-  0}
  
  if(rfx==T & samp==F){
    ## timing and survival and growth (size_t / y_t1) is meant to line up with reproduction (size_t1 / y_t1)
    rfx_surv <- surv_par$tau_year[draw,species,(endo_var_U+1),(year)]; 
    rfx_surv_sdlg <-surv_sdlg_par$tau_year[draw,species,(endo_var_U+1),(year)];
    rfx_grow <- grow_par$tau_year[draw,species,(endo_var_U+1),(year)];
    rfx_grow_sdlg <- grow_sdlg_par$tau_year[draw,species,(endo_var_U+1),(year)];
    rfx_flow <- flow_par$tau_year[draw,species,(endo_var_F+1),(year-repro_offset)]; # fitting 
    rfx_fert <- fert_par$tau_year[draw,species,(endo_var_F+1),(year-repro_offset)]; 
    rfx_spike <- spike_par$tau_year[draw,species,(endo_var_F+1),(year-repro_offset)];
    rfx_rct <- recruit_par$tau_year[draw,species,(endo_var_F+1),year];
  }
  
  if(rfx==T & samp==T){
    rfx_surv <- rnorm(1,mean=0,sd=(surv_par$sigma_year[draw,species,(endo_var_U+1)]+surv_par$sigma_year[draw,species,(endo_var_U+1)]*samp_extreme)) # including sample extreme as a way to increase the sd by a percentage 
    rfx_surv_sdlg <- rnorm(1,mean=0,sd=(surv_sdlg_par$sigma_year[draw,species,(endo_var_U+1)]+surv_sdlg_par$sigma_year[draw,species,(endo_var_U+1)]*samp_extreme)) #surv_sdlg_par$tau_year[draw,species,(endo_var+1),(year)];
    rfx_grow <- rnorm(1,mean=0,sd=(grow_par$sigma_year[draw,species,(endo_var_U+1)]+grow_par$sigma_year[draw,species,(endo_var_U+1)]*samp_extreme)) #grow_par$tau_year[draw,species,(endo_var+1),(year)];
    rfx_grow_sdlg <- rnorm(1,mean=0,sd=(grow_sdlg_par$sigma_year[draw,species,(endo_var_U+1)]+grow_sdlg_par$sigma_year[draw,species,(endo_var_U+1)]*samp_extreme)) #grow_sdlg_par$tau_year[draw,species,(endo_var+1),(year)];
    rfx_flow <- rnorm(1,mean=0,sd=(flow_par$sigma_year[draw,species,(endo_var_F+1)]+flow_par$sigma_year[draw,species,(endo_var_F+1)]*samp_extreme)) #flow_par$tau_year[draw,species,(endo_var+1),year-1]; # fitting 
    rfx_fert <- rnorm(1,mean=0,sd=(fert_par$sigma_year[draw,species,(endo_var_F+1)]+fert_par$sigma_year[draw,species,(endo_var_F+1)]*samp_extreme)) #fert_par$tau_year[draw,species,(endo_var+1),year-1]; 
    rfx_spike <- rnorm(1,mean=0,sd=(spike_par$sigma_year[draw,species,(endo_var_F+1)]+spike_par$sigma_year[draw,species,(endo_var_F+1)]*samp_extreme)) #spike_par$tau_year[draw,species,(endo_var+1),year-1];
    rfx_rct <- rnorm(1,mean=0,sd=(recruit_par$sigma_year[draw,species,(endo_var_F+1)]+recruit_par$sigma_year[draw,species,(endo_var_F+1)]*samp_extreme))#recruit_par$tau_year[draw,species,(endo_var+1),year];
  }  
  if(spei == F){spei_surv <- spei_surv_sdlg <- spei_grow <- spei_grow_sdlg <- spei_flow <- spei_fert <- spei_spike <- spei_rct <-  0}
  if(spei == T){
    spei_surv <- surv_par$betaspei_endo[draw,species,(endo_mean+1)]; 
    spei_surv_sdlg <-surv_sdlg_par$betaspei_endo[draw,species,(endo_mean+1)];
    spei_grow <- grow_par$betaspei_endo[draw,species,(endo_mean+1)];
    spei_grow_sdlg <- grow_sdlg_par$betaspei_endo[draw,species,(endo_mean+1)];
    spei_flow <- flow_par$betaspei_endo[draw,species,(endo_mean+1)]; # fitting 
    spei_fert <- fert_par$betaspei_endo[draw,species,(endo_mean+1)]; 
    spei_spike <- spike_par$betaspei_endo[draw,species,(endo_mean+1)];
    spei_rct <- recruit_par$betaspei_endo[draw,species,(endo_mean+1)];
  }

  
  params <- c()
  #survival
  params$surv_int <- surv_par$beta0[draw,species] + 
    endo_mean_U * surv_par$betaendo[draw,species] + 
    original * surv_par$betaorigin[draw] +
    rfx_surv
  params$surv_spei <- spei_surv
  params$surv_slope <- surv_par$betasize[draw,species]
  # seedling survival
  params$surv_sdlg_int <- surv_sdlg_par$beta0[draw,species] + 
    endo_mean_U * surv_sdlg_par$betaendo[draw,species] + 
    rfx_surv_sdlg
  params$surv_sdlg_spei <- spei_surv_sdlg
  #growth
  params$grow_int <- grow_par$beta0[draw,species] + 
    endo_mean_U * grow_par$betaendo[draw,species] + 
    original * grow_par$betaorigin[draw] +
    rfx_grow
  params$grow_spei <- spei_grow
  params$grow_slope <- grow_par$betasize[draw,species]  
  params$grow_sigma <- grow_par$sigma[draw] 
  # seedling growth
  params$grow_sdlg_int <- grow_sdlg_par$beta0[draw,species] + 
    endo_mean_U * grow_sdlg_par$betaendo[draw,species] + 
    rfx_grow_sdlg
  params$grow_sdlg_spei <- spei_grow_sdlg
  params$grow_sdlg_sigma <- grow_sdlg_par$sigma[draw] 
  
  #flowering
  params$flow_int <- flow_par$beta0[draw,species] + 
    endo_mean_F * flow_par$betaendo[draw,species] + 
    original * flow_par$betaorigin[draw] +
    spei_flow + rfx_flow
  params$flow_spei <- spei_flow
  params$flow_slope <- flow_par$betasize[draw,species]  
  #fertility
  params$fert_int <- fert_par$beta0[draw,species] +
   endo_mean_F * fert_par$betaendo[draw,species] +
   original * fert_par$betaorigin[draw] +
   rfx_fert
  params$fert_spei <- spei_fert
  params$fert_slope <- fert_par$betasize[draw,species]

  #spikelets
  params$spike_int <- spike_par$beta0[draw,species]  +
    endo_mean_F * spike_par$betaendo[draw,species] +
    original * spike_par$betaorigin[draw] +
    rfx_spike
  params$spike_spei <- spei_spike
  params$spike_slope <- spike_par$betasize[draw,species]  
  #seeds per spikelet
  params$seeds_per_spike <- seed_par$beta0[draw,species] + 
    endo_mean_F * seed_par$betaendo[draw,species]
  #recruits per seed
  params$recruits_per_seed <- recruit_par$beta0[draw,species] + 
    endo_mean_F * recruit_par$betaendo[draw,species] +
    rfx_rct
  params$recruits_spei <- spei_rct
  
  #tack on max size
  params$max_size <- max_size$max_size[species]
  
  return(params)
}
# Vital rate functions ----------------------------------------------------
sx<-function(x,params,spei=0){
  xb<-pmin(x,params$max_size) # any predicted plants larger than max size will set to be max size
  invlogit(params$surv_int + params$surv_spei*spei + params$surv_slope*log(xb))
}
sx_sdlg <- function(params, spei=0){
  invlogit(params$surv_sdlg_int+params$surv_sdlg_spei*spei)
}

#gxy <- function(x,y,params){
#  grow_mean <- params$grow_int + params$grow_slope*log(x)
#  
#  grow<-dpoisinvgauss(x=y,mean=exp(grow_mean),shape=(exp(grow_mean)*params$grow_sigma))
#  grow<-ifelse(is.nan(grow) | is.infinite(grow),0,grow)
#  
#  truncLower<-dpoisinvgauss(x=0,mean=exp(grow_mean), shape=(exp(grow_mean)*params$grow_sigma))
#  # truncLower<-sum(ifelse(is.nan(truncLower) | is.infinite(truncLower),0,truncLower))
#  
#  truncUpper<-sum(dpoisinvgauss(x=params$max_size:10000,mean=exp(grow_mean),shape=(exp(grow_mean)*params$grow_sigma)))
#  # truncUpper<-sum(ifelse(is.nan(truncUpper) | is.infinite(truncUpper),0,truncUpper))
#  return(grow/(1-(truncLower+truncUpper)))
#}

gxy <- function(x,y,params, spei=0){
  xb<-pmin(x,params$max_size)
  grow_mean <- params$grow_int + params$grow_spei*spei + params$grow_slope*log(xb)
  grow<-dpoisinvgauss(x=y,mean=exp(grow_mean),shape=(exp(grow_mean)*params$grow_sigma))
  truncZero<-dpoisinvgauss(x=0,mean=exp(grow_mean),shape=(exp(grow_mean)*params$grow_sigma))
  return(grow/(1-truncZero))
}

gxy_sdlg <- function(x,y,params, spei=0){
  grow_mean <- params$grow_sdlg_int + params$grow_sdlg_spei*spei
  grow<-dpoisinvgauss(x=y,mean=exp(grow_mean),shape=(exp(grow_mean)*params$grow_sigma))
  truncZero<-dpoisinvgauss(x=0,mean=exp(grow_mean),shape=(exp(grow_mean)*params$grow_sigma))
  return(grow/(1-truncZero))
}

pxy<-function(x,y,params,spei=0){
  sx(x,params,spei) * gxy(x,y,params,spei)
}

fx<-function(x, params,spei=0){
  xb<-pmin(x,params$max_size) # any predicted plants larger than max size will set to be max size
  flw <- invlogit(params$flow_int + params$flow_spei*spei + params$flow_slope*log(xb))
  fert <- exp(params$fert_int + params$fert_spei*spei + params$fert_slope*log(xb))
  spike <- exp(params$spike_int + params$spike_spei*spei + params$spike_slope*log(xb))
  seeds_per_spike <- params$seeds_per_spike
  recruits_per_seed <- invlogit(params$recruits_per_seed + params$recruits_spei*spei)
  seedlings <- flw * fert * spike * seeds_per_spike * recruits_per_seed
  return(seedlings)
}

# Bigmatrix function ------------------------------------------------------
# This includes a reproductive delay till the first tiller
bigmatrix<-function(params,spei=0,extension=0){   
  matdim<-params$max_size+extension #extension here adds sizes larger than max size which accounts for probability density lost at predicted sizes larger than our upper bound
  y <- 1:matdim 
  #fertility transition
  Fmat <- matrix(0,matdim+1,matdim+1)
  Fmat[1,2:(matdim+1)]<-fx(x = y, params, spei)
  
  #growth/survival transition
  Tmat <-matrix(0,matdim+1,matdim+1)
  Tmat[2:(matdim+1),2:(matdim+1)] <- t(outer(y,y,pxy,params,spei))
  # surviving seedlings emerge into continuous population
  Tmat[2:(matdim+1),1] <- gxy_sdlg(x=1,y=y, params = params,spei)*sx_sdlg(params = params,spei)
  MPMmat<-Tmat + Fmat
  return(list(MPMmat = MPMmat, Tmat = Tmat, Fmat = Fmat))
}

# lambdaS function##########################################################
lambdaSim<-function(mat_list, ## a list of transition matrices, each corresponding to a study year
                    max_yrs=1000 ## how many years the simulation runs (arbitrarily large)
){
  ## grab the dimension of the projection matrix
  matdim<-dim(mat_list[[1]])[1]
  ## grab the number of study years / matrices we have available
  n_years <- length(mat_list)
  ## vector that will hold year-by-year growth rates
  rtracker <- rep(0,max_yrs)
  ## initial vector of population structure -- note that this sums to one, which will be convenient
  n0 <- rep(1/matdim,matdim)
  for(t in 1:max_yrs){ #Start loop
    ## for each year, randomly sample one of the matrices
    A_t <- mat_list[[sample.int(n=n_years,size=1)]]
    ## project the population one step forward
    n0 <- A_t %*% n0
    ## total population size after one year of growth
    N  <- sum(n0)
    ## calculate r as log(N_t+1 / N_t), note that here N_t=1
    rtracker[t]<-log(N)
    ## rescale population vector to sum to one, so the same trick works again next time step
    n0 <-n0/N
  }
  #discard initial values (to get rid of transient)
  burnin    <- round(max_yrs*0.1)
  #Finish and return
  log_lambdaS <- mean(rtracker[-c(1:burnin)])
  lambdaS<-exp(log_lambdaS)
  return(list(log_lambdaS=log_lambdaS,lambdaS=lambdaS))
}
