
data {
    // indices
    int<lower=0> nYear;                       // number of years
    int<lower=0> nPlot;                         // number of plots
    int<lower=0> nSpp;                        // number of host species
    int<lower=0> nEndo;                        // number of endophyte levels
    int<lower=0> nOrigin;                    // number of plant origin levels
    // vital rate data
    int<lower=0> N;                       // number of observations for surv model
    int<lower=0, upper=nYear> year_t[N];         // year of observation for surv model
    int<lower=0> plot[N];                   // plot of observation for surv model
    int<lower=0, upper=nSpp> spp[N];         // year of observation for surv model
    int<lower=0, upper=1> y[N];      // plant survival at time t+1 or flowering at time t+1
    vector<lower=0>[N] logsize;             // plant size at time t for surv model or time t+1 for flowering
    int<lower=0,upper=1> endo_01[N];            // plant endophyte status for surv model
    int<lower=1,upper=2> origin_index[N];          // plant origin status for surv model
}

transformed data {
  vector<lower=0>[N] logsize_2;
  for(n in 1:N){
  logsize_2[n] = logsize[n]^2;
  }
}
parameters {
    // vr params
    real beta0[nSpp,nOrigin];                  // predictor parameters as grand means 
    real betasize[nSpp,nOrigin];                  //   spp specific size slope 
    real betasize_2[nSpp,nOrigin];                  //   spp specific size slope - quadratic term
    vector[nSpp] betaendo;                  // spp specific endophyt effect 
    //vector[nSpp] betaorigin;               // spp specific origin effect --> I think we don't need this origin effect in this model since I just made the size effects completely crossed with origin and species
    
    real tau_year[nSpp,nEndo,nYear];      // random year effect, unique to species and endo

    vector[nSpp] sigma0;                 // year variance
    vector[nSpp] sigmaendo;              // endo effect on variance

    vector[nPlot] tau_plot;        // random plot effect
    real<lower=0> sigma_plot;          // plot variance effect
}

transformed parameters {
    real p[N];                           
    real sigma_year[nSpp,nEndo];

    // surv Linear Predictor
    for(n in 1:N){
    p[n] = beta0[spp[n],origin_index[n]] + betasize[spp[n],origin_index[n]]*logsize[n] + betasize_2[spp[n],origin_index[n]]*logsize_2[n] + betaendo[spp[n]]*endo_01[n]
    + tau_year[spp[n],(endo_01[n]+1),year_t[n]] 
    + tau_plot[plot[n]]
    ;
    }
    
    // endo effect on variance
    for(s in 1:nSpp){
      for(d in 1:nEndo){
        sigma_year[s,d] = exp(sigma0[s] + sigmaendo[s]*(d-1));
      }
    }
}

model {
    // priors
    //this is plot variance
      tau_plot ~ normal(0,sigma_plot);
      sigma_plot ~ normal(0, 1);

      
    //fixed effect priors
      to_vector(beta0[,1]) ~ normal(0,5);
      to_vector(beta0[,2]) ~ normal(0,5);

      betaendo ~ normal(0,5);
      to_vector(betasize[,1]) ~ normal(0,5); // this should set up the prior for the original plants size effect across species
      to_vector(betasize[,2]) ~ normal(0,5); // this should set up the prior for the recruit plants size effect across species

      to_vector(betasize_2[,1]) ~ normal(0,5);
      to_vector(betasize_2[,2]) ~ normal(0,5);

      //betaorigin ~ normal(0,5);      
      sigma0 ~ normal(0,1);
      sigmaendo ~ normal(0,1);
      
    //species endo year priors
    for(s in 1:nSpp){
          to_vector(tau_year[s,1,]) ~ normal(0,sigma_year[s,1]); // sample year effects for each species for each endo status
          to_vector(tau_year[s,2,]) ~ normal(0,sigma_year[s,2]);
    }
    
    y ~ bernoulli_logit(p);
}
    
