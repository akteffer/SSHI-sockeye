data {
    int<lower=0> n;    //observations 
    int<lower=0> PRV[n];  //response (0: no PRV, 1: PRV)
    
    int<lower=0> n_zones;    
    matrix[n,n_zones] zones; //zone associated with each observation (design matrix)
   
    int<lower=0> n_farms;    
    matrix[n,n_farms] farms; //farm associated with each observation (design matrix)
   
    int<lower=0> n_samples; //number of fish-health-audit sampling events
    matrix[n,n_samples] samples;  //visit ID for each sample (design matrix)
    
    int<lower=0> n_covariates;
    matrix[n,n_covariates] X;
}
parameters {
// fixed effects    
    real beta_0;
    vector[n_covariates] beta; // fixed effects

// random effects
    vector[n_zones] z_prelim; // 
    vector[n_farms] f_int_prelim; //farm random intercept
    vector[n_farms] f_slope_prelim; //farm random slope (of time in SW)
    vector[n_samples] s_prelim; //the (random) sampling event mean

// zone standard deviation
    real<lower=0> sigma_Z;
    
// zone standard deviation
    real<lower=0> sigma_F_int;
    real<lower=0> sigma_F_slope;
    real<lower=-1,upper=1> cor_F;
    
// sampling event standard deviation
    real<lower=0> sigma_S;

}
transformed parameters {
    vector[n_samples] s;  //levels of the sampling-event-level random intercept
    matrix[n_farms,2] f; 
    matrix[2,2] vcov;
//    matrix[2,2] L;
    vector[n_zones] z;  //levels of the zone-level random slope
   
    s = sigma_S*s_prelim;    //scale the random effects
    
    //use the magic of the Cholesky decomposition to implement the "centred" random effect model for correlated random effects
    vcov[1,1] = sigma_F_int*sigma_F_int;
    vcov[1,2] = cor_F*sigma_F_int*sigma_F_slope;
    vcov[2,1] = cor_F*sigma_F_int*sigma_F_slope;
    vcov[2,2] = sigma_F_slope*sigma_F_slope;
    f = append_col(f_int_prelim,f_slope_prelim)*cholesky_decompose(vcov)';
   
    z = sigma_Z*z_prelim;    //scale the random effects
}
model {
    vector[n] mu;
    
    //assemble the mean of the linear predictor, for this iteration
    mu = beta_0 + X*beta; 
    mu = mu + samples*s + farms*f[,1];     //sampling-event-level random effect
    for(i in 1:n){
        mu[i] = mu[i] + X[i,1]*(zones[i]*z) + X[i,1]*(farms[i]*f[,2]);    //add in random time-since-entry component
    }
    for(i in 1:n){
        target += binomial_lpmf(PRV[i] | 1, 1/(1 + exp(-mu[i])));//inv_logit(mu[i]));//   //final likelihood
    }
    
    //spatial random effect
    target += normal_lpdf(z_prelim | 0, 1);  
    
    //farm random effect
    target += normal_lpdf(f_int_prelim | 0, 1);
    target += normal_lpdf(f_slope_prelim | 0, 1);
     
    //sampling-event random effect
    target += normal_lpdf(s_prelim | 0, 1); 
    
    //standard deviation(s) & correlation(s) for random effects
    target += uniform_lpdf(sigma_S | 0,10);
    
    target += uniform_lpdf(sigma_F_int | 0,10);
    target += uniform_lpdf(sigma_F_slope | 0,10);
    target += uniform_lpdf(cor_F | -1,1);
    
    target += uniform_lpdf(sigma_Z | 0,10);
    
    //fixed effects
    target += uniform_lpdf(beta_0 | -25, 25);//beta_0 ~ uniform(-10,10);//
    target += uniform_lpdf(beta[1] | -25, 25);//beta_0 ~ uniform(-10,10);//
   
    //print("post-sampling target: ", target());
}
//generated quantities {




