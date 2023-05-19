library(postpack)
library(tidyr)
library(dplyr)
library(Hmisc)
library(pals)
library(ggplot2)
library(gplots)
library(ggpubr)
library(scales)

#### freshwater analysis ####
d_fw<-readRDS("sockFWdatLW.rds")

pathogens<-c("arena2","c_b_cys","ce_sha","Circo.virus","de_sal","fl_psy",
             "ic_hof","ic_mul","IcD","ihnv","lo_sal","my_arc","nu_sal",
             "pa_kab","pa_min","pa_ther","prv","pspv","Qin","Rhabdo3_virus",
             "rlo","sp_des","te_bry","te_mar")

# remove these models due to insufficient positives
remove<-c() 

spsu_load<-NULL
warning_list<-list()
converge_list<-list()
count<-0
post_probs<-data.frame(matrix(NA,nrow=150,ncol=3))
ppath_list<-list()


for(j in 1:length(pathogens)){
	path_list<-list()
  path<-pathogens[j]

d_fw$pathogen<-d_fw[,match(path,colnames(d_fw))]
# remove any NAs (already accounted for length and mass)
dat2<-d_fw %>%
	filter(!is.na(jaz),
	        !is.na(pathogen),
          !is.na(Length.x),
          !is.na(Mass.x),
	        jaz!="UFR+GStr") 

### step 1: data ###
jags_data = list(
	n_yrs = length(unique(dat2$ocean_yr)),  
	n_obs = nrow(dat2),
	n_jaz = length(unique(dat2$jaz)),  
	yr = as.numeric(as.factor(dat2$ocean_yr)),
	jaz = as.numeric(as.factor(dat2$jaz)),
	Mass = arm::rescale(log(dat2$Mass.x)),
	length = arm::rescale(log(dat2$Length.x)),
	pathogen = dat2$pathogen
	)
### step 2: model ###
jags_model = function() {
  # error variability priors
  sig_yr ~ dunif(0,2)
  sig_resid ~ dunif(0,2)
  sig_jaz ~ dunif(0,2)
  sig_B2 ~ dunif(0,2)

  # power transform
  pow ~ dunif(0,1)
  
  # the random effect for year: influences the intercepts only
  for (y in 1:n_yrs) {
  	yr_eff[y] ~ dnorm(0, 1/sig_yr^2)
  }
  
  # random effect for JAZ: slope and intercept
  for (s in 1:n_jaz){
    jaz_eff[s] ~ dnorm(0, 1/sig_jaz^2)
    #b2[s] ~ dnorm(B2,1/sig_B2^2)
  }
  
  #uninformed priors for betas
  b0 ~ dnorm(0, 1E-6)
  b1 ~ dnorm(0, 1E-6)
  b2 ~ dnorm(0, 1E-6)

  for (q in 1:n_obs) {
    Mass[q] ~ dnorm(Mass_hat[q], 1/sig_resid^2)
    Mass_hat[q] <- b0 + b1*length[q] + b2*(pathogen[q])^pow + jaz_eff[jaz[q]] + yr_eff[yr[q]] #
    residual[q] <- Mass[q] - Mass_hat[q]
  }
}

jags_file = "model.txt"
postpack::write_model(jags_model, jags_file)

### step 3: initial values ###
fit = with(jags_data, lme4::lmer(Mass ~ length + pathogen + (1|jaz) + (1|yr) )) 
coefs = lme4::fixef(fit)
coef_names = names(coefs)

# extract which coefficients are what
int = stringr::str_detect(coef_names, "Intercept")
slope = stringr::str_detect(coef_names, "^length$")
pathogen = stringr::str_detect(coef_names, "^pathogen$")

# obtain the coefficient estimates in the parameterization estimated by JAGS model
b0_fit = unname(coefs[int]) 
b1_fit = unname(coefs[slope]) 
b2_fit = unname(coefs[pathogen])
sig_yr_fit = 0.1
sig_jaz_fit = 0.1
sig_mon_fit = 0.1
sig_resid_fit = summary(fit)$sigma

jags_inits = function(nc) {
  inits = list()
  for (c in 1:nc) {
    inits[[c]] = list(
 
      b0 = rnorm(1, b0_fit, abs(b0_fit) * 0.1),
      b1 = rnorm(1, b1_fit, abs(b1_fit) * 0.2),
      b2 = rnorm(1, b2_fit, abs(b2_fit) * 0.2),
      #B2 = rnorm(1, b2_fit, 1e-3),
      pow=0.05,
      sig_B2 = runif(1,0,1),
      sig_yr = runif(1, sig_yr_fit * 0.5, sig_yr_fit * 1.5),
      sig_jaz = runif(1, sig_jaz_fit * 0.5, sig_jaz_fit * 1.5),
      sig_resid = runif(1, sig_resid_fit * 0.5, sig_resid_fit * 1.5)
    )
  }
  return(inits)
}

### step 4: parameters to monitor ###
jags_params = c("b0", "b1","b2","pow", "sig_resid", "sig_jaz",  "Mass_hat", "residual","jaz_eff","yr_eff")

### step 5: dimensions ###
jags_dims = c(
  ni = 10000,  # number of post-burn-in samples per chain
  nb = 10000,  # number of burn-in samples
  nt = 1,     # thinning rate
  nc = 3      # number of chains
)

### step 6: mcmc ###
post = jagsUI::jags.basic(
  data = jags_data,
  model.file = jags_file,
  inits = jags_inits(jags_dims["nc"]),
  parameters.to.save = jags_params,
  n.adapt = 1000,
  n.iter = sum(jags_dims[c("ni", "nb")]),
  n.thin = jags_dims["nt"],
  n.burnin = jags_dims["nb"],
  n.chains = jags_dims["nc"],
  parallel = T
)

if( length(warnings())>0){
  warning_list[[length(warning_list)+1]]<-warnings()
  names(warning_list)[[length(warning_list)]]<-pathogens[j]
  assign("last.warning", NULL, envir = baseenv())  #this step clears the warnings
}

### step 7: convergence diagnosis ###
diag<-t(postpack::post_summ(post, c("B|b","sig"), neff = T, Rhat = T)[c("Rhat", "neff"),])
converge_list[[j]]<-diag

### step 8: inference ###
b_est = postpack::post_summ(post, "B|b",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.975))
b2_est = postpack::post_summ(post, "b2",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.975))

pwr<-post_summ(post,"pow")

spsu_load<-cbind(spsu_load,b2_est[,1])
count<-count+1

post_path<-post_subset(post,"b2",matrix=T)

post_probs[count,1:3]<-c(path,"spsu",length(which(post_path<0))/nrow(post_path))

ppath_list[[count]]<-post_path
}

names(converge_list)<-pathogens

names(ppath_list)<-pathogens

# now adding the coinfection variables
mod_ind<-c("richness","rib")

# use this to loop
for(j in 1:length(mod_ind)){
  path_list<-list()
  
  path<-mod_ind[j]

  d_fw$pathogen<-d_fw[,match(mod_ind[j],colnames(d_fw))]

  # remove any NAs (already accounted for length and mass)
  dat2<-d_fw %>%
    filter(!is.na(jaz),
           !is.na(pathogen),
           !is.na(Length.x),
           !is.na(Mass.x),
           jaz!="UFR+GStr") 

  ### step 1: data ###
  jags_data = list(
    n_yrs = length(unique(dat2$ocean_yr)),  
    n_obs = nrow(dat2),
    n_jaz = length(unique(dat2$jaz)),  
    yr = as.numeric(as.factor(dat2$ocean_yr)),
    jaz = as.numeric(as.factor(dat2$jaz)),
    Mass = arm::rescale(log(dat2$Mass.x)),
    length = arm::rescale(log(dat2$Length.x)),
    pathogen = dat2$pathogen
  )
  ### step 2: model ###
  jags_model = function() {
    # error variability priors
    sig_yr ~ dunif(0,2)
    sig_resid ~ dunif(0,2)
    sig_prior ~ dunif(0,2)
    sig_jaz ~ dunif(0,2)
    sig_B2 ~ dunif(0,2)

    # the random effect for year: influences the intercepts only
    for (y in 1:n_yrs) {
      yr_eff[y] ~ dnorm(0, 1/sig_yr^2)
    }
    
    # random effect for JAZ: influences intercept only
    for (s in 1:n_jaz){
      jaz_eff[s] ~ dnorm(0, 1/sig_jaz^2)
      #b2[s] ~ dnorm(B2,1/sig_B2^2)
    }
    
    #uninformed priors for betas
    b0 ~ dnorm(0, 1E-6)
    b1 ~ dnorm(0, 1E-6)
    b2 ~ dnorm(0, 1E-6)
    
    for (q in 1:n_obs) {
      Mass[q] ~ dnorm(Mass_hat[q], 1/sig_resid^2)
      Mass_hat[q] <- b0 + b1*length[q] + b2*(pathogen[q]) + jaz_eff[jaz[q]] + yr_eff[yr[q]] 
      residual[q] <- Mass[q] - Mass_hat[q]
    }
  }
  
  jags_file = "model.txt"
  postpack::write_model(jags_model, jags_file)
  
  ### step 3: initial values ###
  fit = with(jags_data, lme4::lmer(Mass ~ length + pathogen + (1|jaz) + (1|yr) )) 
  coefs = lme4::fixef(fit)
  coef_names = names(coefs)
  
  # extract which coefficients are what
  int = stringr::str_detect(coef_names, "Intercept")
  slope = stringr::str_detect(coef_names, "^length$")
  pathogen = stringr::str_detect(coef_names, "^pathogen$")
  
  # obtain the coefficient estimates in the parameterization estimated by JAGS model
  b0_fit = unname(coefs[int]) 
  b1_fit = unname(coefs[slope]) 
  b2_fit = unname(coefs[pathogen])
  sig_yr_fit = 0.1
  sig_jaz_fit = 0.1
  sig_mon_fit = 0.1
  sig_resid_fit = summary(fit)$sigma
  
  jags_inits = function(nc) {
    inits = list()
    for (c in 1:nc) {
      inits[[c]] = list(
        b0 = rnorm(1, b0_fit, abs(b0_fit) * 0.1),
        b1 = rnorm(1, b1_fit, abs(b1_fit) * 0.2),
        b2 = rnorm(1, b2_fit, abs(b2_fit) * 0.2),
        #B2 = rnorm(1, b2_fit, 1e-3),
        sig_B2 = runif(1,0,1),
        sig_yr = runif(1, sig_yr_fit * 0.5, sig_yr_fit * 1.5),
        sig_jaz = runif(1, sig_jaz_fit * 0.5, sig_jaz_fit * 1.5),
        sig_resid = runif(1, sig_resid_fit * 0.5, sig_resid_fit * 1.5)
      )
    }
    return(inits)
  }
  
  ### step 4: parameters to monitor ###
  jags_params = c("b0", "b1","b2", "sig_resid", "sig_jaz",  "Mass_hat", "residual","jaz_eff","yr_eff")
  
  ### step 5: dimensions ###
  jags_dims = c(
    ni = 10000,  # number of post-burn-in samples per chain
    nb = 10000,  # number of burn-in samples
    nt = 1,     # thinning rate
    nc = 3      # number of chains
  )
  
  ### step 6: mcmc ###
  post = jagsUI::jags.basic(
    data = jags_data,
    model.file = jags_file,
    inits = jags_inits(jags_dims["nc"]),
    parameters.to.save = jags_params,
    n.adapt = 1000,
    n.iter = sum(jags_dims[c("ni", "nb")]),
    n.thin = jags_dims["nt"],
    n.burnin = jags_dims["nb"],
    n.chains = jags_dims["nc"],
    parallel = T
  )
  
  if( length(warnings())>0){
    warning_list[[length(warning_list)+1]]<-warnings()
    names(warning_list)[[length(warning_list)]]<-mod_ind[j]
    assign("last.warning", NULL, envir = baseenv())  #this step clears the warnings
  }
    
  count<-count+1
  ### step 7: convergence diagnosis ###
  diag<-t(postpack::post_summ(post, c("B|b","sig"), neff = T, Rhat = T)[c("Rhat", "neff"),])
  converge_list[[count]]<-diag

  ### step 8: inference ###
  b_est = postpack::post_summ(post, "B|b",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.975))
  b2_est = postpack::post_summ(post, "b2",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.975))
  
  spsu_load<-cbind(spsu_load,b2_est[,1])
  
  post_path<-post_subset(post,"b2",matrix=T)
  
  post_probs[count,1:3]<-c(path,"spsu",length(which(post_path<0))/nrow(post_path))
  
  ppath_list[[count]]<-post_path
  
}

names(converge_list)[25:26]<-mod_ind

names(ppath_list)[25:26]<-mod_ind

pathogens<-c("arena2","c_b_cys","ce_sha","Circo.virus","de_sal","fl_psy","ic_hof","ic_mul",
             "IcD","ihnv","lo_sal","my_arc","nu_sal","pa_kab","pa_min",
             "pa_ther","prv","pspv","Qin","Rhabdo3_virus","rlo",
             "sp_des","te_bry","te_mar","richness","rib")

colnames(spsu_load)<-pathogens

# collect Rhat and ess
conv<-matrix(NA,nrow=26,ncol=6)
for(i in 1:length(names(converge_list))){
  conv[i,1]<-converge_list[[i]][2,1]
  conv[i,2]<-converge_list[[i]][2,2]
  conv[i,3]<-converge_list[[i]][3,1]
  conv[i,4]<-converge_list[[i]][3,2]
  conv[i,5]<-converge_list[[i]][1,1]
  conv[i,6]<-converge_list[[i]][1,2]

}

colnames(conv)<-c("Rhat b0","ess b0","Rhat b1", "ess b1", "Rhat B2", "ess B2")
rownames(conv)<-names(converge_list)

post_probs<-post_probs[1:26,]

sockeye_lw<-list(spsu_load,post_probs,converge_list,conv,ppath_list)
names(sockeye_lw)<-c("estimates","posterior_probabilities","convergence_data","convergence_summary","samples")
#saveRDS(sockeye_lw,"sockeye_FW_norandslope_lw_results_jan26_2023.rds")

#### marine analysis ####
d_sw<-readRDS("sockSWdatLW.rds")

pathogens<-c("arena2","c_b_cys","ce_sha","Circo.virus","de_sal","fa_mar",
              "fl_psy","ic_hof","ic_mul","IcD",
              "ku_thy","lo_sal","mo_vis","my_arc","nu_sal","pa_kab","pa_min",
              "pa_pse","pa_ther","prv","pspv","Qin","Rhabdo3_virus","rlo",
              "sch","smallUK","sp_des","te_bry","te_mar","ven")

spsu_load<-NULL
warning_list<-list()
converge_list<-list()
count<-0
post_probs<-data.frame(matrix(NA,nrow=150,ncol=3))
ppath_list<-list()


for(j in 1:length(pathogens)){
  path_list<-list()
  path<-pathogens[j]
  
  d_sw$pathogen<-d_sw[,match(path,colnames(d_sw))]
  # remove any NAs (already accounted for length and mass)
  dat2<-d_sw %>%
    filter(!is.na(jaz),
           !is.na(pathogen),
           !is.na(Length),
           !is.na(Mass),
           jaz!="UFR+GStr",
           !is.na(temp_devs)) 
  
  ### step 1: data ###
  jags_data = list(
    n_yrs = length(unique(dat2$ocean_yr)),  
    n_obs = nrow(dat2),
    n_jaz = length(unique(dat2$jaz)),  
    yr = as.numeric(as.factor(dat2$ocean_yr)),
    jaz = as.numeric(as.factor(dat2$jaz)),
    Mass = arm::rescale(log(dat2$Mass)),
    length = arm::rescale(log(dat2$Length)),
    pathogen = dat2$pathogen,
    sst = arm::rescale(dat2$temp_devs)
      )
  ### step 2: model ###
  jags_model = function() {
    # error variability priors
    sig_yr ~ dunif(0,2)
    sig_resid ~ dunif(0,2)
    sig_jaz ~ dunif(0,2)
    sig_B2 ~ dunif(0,2)
    
    # power transform
    pow ~ dunif(0,1)
    
    # the random effect for year: influences the intercepts only
    for (y in 1:n_yrs) {
      yr_eff[y] ~ dnorm(0, 1/sig_yr^2)
    }
    
    # random effect for JAZ: slope and intercept
    for (s in 1:n_jaz){
      jaz_eff[s] ~ dnorm(0, 1/sig_jaz^2)
      b2[s] ~ dnorm(B2,1/sig_B2^2)
    }
    
    #uninformed priors for betas
    b0 ~ dnorm(0, 1E-6)
    b1 ~ dnorm(0, 1E-6)
    B2 ~ dnorm(0, 1E-6)
    b3 ~ dnorm(0, 1E-6)
    
    for (q in 1:n_obs) {
      Mass[q] ~ dnorm(Mass_hat[q], 1/sig_resid^2)
      Mass_hat[q] <- b0 + b1*length[q] + b2[jaz[q]]*(pathogen[q])^pow + b3*sst[q] + jaz_eff[jaz[q]] + yr_eff[yr[q]] #
      residual[q] <- Mass[q] - Mass_hat[q]
    }
  }
  
  jags_file = "model.txt"
  postpack::write_model(jags_model, jags_file)
  
  ### step 3: initial values ###
  fit = with(jags_data, lme4::lmer(Mass ~ length + pathogen + sst + (1|jaz) + (1|yr) )) 
  coefs = lme4::fixef(fit)
  coef_names = names(coefs)
  
  # extract which coefficients are what
  int = stringr::str_detect(coef_names, "Intercept")
  slope = stringr::str_detect(coef_names, "^length$")
  pathogen = stringr::str_detect(coef_names, "^pathogen$")
  sstdev = stringr::str_detect(coef_names, "^sst$")
  
  # obtain the coefficient estimates in the parameterization estimated by JAGS model
  b0_fit = unname(coefs[int]) 
  b1_fit = unname(coefs[slope]) 
  b2_fit = unname(coefs[pathogen])
  b3_fit = unname(coefs[sstdev])
  sig_yr_fit = 0.1
  sig_jaz_fit = 0.1
  sig_resid_fit = summary(fit)$sigma
  
  jags_inits = function(nc) {
    inits = list()
    for (c in 1:nc) {
      inits[[c]] = list(
        
        b0 = rnorm(1, b0_fit, abs(b0_fit) * 0.1),
        b1 = rnorm(1, b1_fit, abs(b1_fit) * 0.2),
        b2 = rnorm(jags_data$n_jaz, b2_fit, abs(b2_fit) * 0.2),
        B2 = rnorm(1, b2_fit, 1e-3),
        b3 = rnorm(1, b3_fit, abs(b3_fit) * 0.2),
        pow=0.05,
        sig_B2 = runif(1,0,1),
        sig_yr = runif(1, sig_yr_fit * 0.5, sig_yr_fit * 1.5),
        sig_jaz = runif(1, sig_jaz_fit * 0.5, sig_jaz_fit * 1.5),
        sig_resid = runif(1, sig_resid_fit * 0.5, sig_resid_fit * 1.5)
      )
    }
    return(inits)
  }
  
  ### step 4: parameters to monitor ###
  jags_params = c("b0", "b1","b2","b3","pow","B2", "sig_resid", "sig_jaz",  "Mass_hat", "residual","jaz_eff","yr_eff")
  
  ### step 5: dimensions ###
  jags_dims = c(
    ni = 10000,  # number of post-burn-in samples per chain
    nb = 10000,  # number of burn-in samples
    nt = 1,     # thinning rate
    nc = 3      # number of chains
  )
  
  ### step 6: mcmc ###
  post = jagsUI::jags.basic(
    data = jags_data,
    model.file = jags_file,
    inits = jags_inits(jags_dims["nc"]),
    parameters.to.save = jags_params,
    n.adapt = 1000,
    n.iter = sum(jags_dims[c("ni", "nb")]),
    n.thin = jags_dims["nt"],
    n.burnin = jags_dims["nb"],
    n.chains = jags_dims["nc"],
    parallel = T
  )
  
  if( length(warnings())>0){
    warning_list[[length(warning_list)+1]]<-warnings()
    names(warning_list)[[length(warning_list)]]<-pathogens[j]
    assign("last.warning", NULL, envir = baseenv())  #this step clears the warnings
  }
  
  ### step 7: convergence diagnosis ###
  diag<-t(postpack::post_summ(post, c("B|b","sig"), neff = T, Rhat = T)[c("Rhat", "neff"),])
  converge_list[[j]]<-diag
  
  ### step 8: inference ###
  b_est = postpack::post_summ(post, "B|b",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.975))
  b2_est = postpack::post_summ(post, "B2|b2",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.975))
  b3_est = postpack::post_summ(post, "b3",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.95))
  
  
  pwr<-post_summ(post,"pow")
  
  spsu_load<-cbind(spsu_load,b2_est[,1])
  count<-count+1
  
  post_path<-post_subset(post,"B2",matrix=T)
  post_sst<-post_subset(post,"b3",matrix=T)
  
  
  post_probs[count,1:4]<-c(path,"spsu",length(which(post_path<0))/nrow(post_path),
                           length(which(post_sst<0))/nrow(post_sst))
  
  ppath_list[[count]]<-post_path
}

  names(converge_list)<-pathogens
  
  names(ppath_list)<-pathogens

  
  # now adding the coinfection variables
  mod_ind<-c("richness","rib")
  
  # use this to loop
  for(j in 1:length(mod_ind)){
    path_list<-list()
    
    path<-mod_ind[j]
    
    d_sw$pathogen<-d_sw[,match(mod_ind[j],colnames(d_sw))]
    
    # remove any NAs (already accounted for length and mass)
    dat2<-d_sw %>%
      filter(!is.na(jaz),
             !is.na(pathogen),
             !is.na(Length),
             !is.na(Mass),
             jaz!="UFR+GStr",
             !is.na(temp_devs)) 
    
    ### step 1: data ###
    jags_data = list(
      n_yrs = length(unique(dat2$ocean_yr)),  
      n_obs = nrow(dat2),
      n_jaz = length(unique(dat2$jaz)),  
      yr = as.numeric(as.factor(dat2$ocean_yr)),
      jaz = as.numeric(as.factor(dat2$jaz)),
      Mass = arm::rescale(log(dat2$Mass)),
      length = arm::rescale(log(dat2$Length)),
      pathogen = dat2$pathogen,
      sst = arm::rescale(dat2$temp_devs)
    )
    ### step 2: model ###
    jags_model = function() {
      # error variability priors
      sig_yr ~ dunif(0,2)
      sig_resid ~ dunif(0,2)
      sig_jaz ~ dunif(0,2)
      sig_B2 ~ dunif(0,2)
      
      # the random effect for year: influences the intercepts only
      for (y in 1:n_yrs) {
        yr_eff[y] ~ dnorm(0, 1/sig_yr^2)
      }
      
      # random effect for JAZ: slope and intercept
      for (s in 1:n_jaz){
        jaz_eff[s] ~ dnorm(0, 1/sig_jaz^2)
        b2[s] ~ dnorm(B2,1/sig_B2^2)
      }
      
      #uninformed priors for betas
      b0 ~ dnorm(0, 1E-6)
      b1 ~ dnorm(0, 1E-6)
      B2 ~ dnorm(0, 1E-6)
      b3 ~ dnorm(0, 1E-6)
      
      for (q in 1:n_obs) {
        Mass[q] ~ dnorm(Mass_hat[q], 1/sig_resid^2)
        Mass_hat[q] <- b0 + b1*length[q] + b2[jaz[q]]*(pathogen[q]) + b3*sst[q] + jaz_eff[jaz[q]] + yr_eff[yr[q]] #
        residual[q] <- Mass[q] - Mass_hat[q]
      }
    }
    
    jags_file = "model.txt"
    postpack::write_model(jags_model, jags_file)
    
    ### step 3: initial values ###
    fit = with(jags_data, lme4::lmer(Mass ~ length + pathogen + sst + (1|jaz) + (1|yr) )) 
    coefs = lme4::fixef(fit)
    coef_names = names(coefs)
    
    # extract which coefficients are what
    int = stringr::str_detect(coef_names, "Intercept")
    slope = stringr::str_detect(coef_names, "^length$")
    pathogen = stringr::str_detect(coef_names, "^pathogen$")
    sstdev = stringr::str_detect(coef_names, "^sst$")
    
    # obtain the coefficient estimates in the parameterization estimated by JAGS model
    b0_fit = unname(coefs[int]) 
    b1_fit = unname(coefs[slope]) 
    b2_fit = unname(coefs[pathogen])
    b3_fit = unname(coefs[sstdev])
    sig_yr_fit = 0.1
    sig_jaz_fit = 0.1
    sig_resid_fit = summary(fit)$sigma
    
    jags_inits = function(nc) {
      inits = list()
      for (c in 1:nc) {
        inits[[c]] = list(
          
          b0 = rnorm(1, b0_fit, abs(b0_fit) * 0.1),
          b1 = rnorm(1, b1_fit, abs(b1_fit) * 0.2),
          b2 = rnorm(jags_data$n_jaz, b2_fit, abs(b2_fit) * 0.2),
          B2 = rnorm(1, b2_fit, 1e-3),
          b3 = rnorm(1, b3_fit, abs(b3_fit) * 0.2),
          sig_B2 = runif(1,0,1),
          sig_yr = runif(1, sig_yr_fit * 0.5, sig_yr_fit * 1.5),
          sig_jaz = runif(1, sig_jaz_fit * 0.5, sig_jaz_fit * 1.5),
          sig_resid = runif(1, sig_resid_fit * 0.5, sig_resid_fit * 1.5)
        )
      }
      return(inits)
    }
    
    ### step 4: parameters to monitor ###
    jags_params = c("b0", "b1","b2","b3","B2", "sig_resid", "sig_jaz",  "Mass_hat", "residual","jaz_eff","yr_eff")
    
    ### step 5: dimensions ###
    jags_dims = c(
      ni = 10000,  # number of post-burn-in samples per chain
      nb = 10000,  # number of burn-in samples
      nt = 1,     # thinning rate
      nc = 3      # number of chains
    )
    
    ### step 6: mcmc ###
    post = jagsUI::jags.basic(
      data = jags_data,
      model.file = jags_file,
      inits = jags_inits(jags_dims["nc"]),
      parameters.to.save = jags_params,
      n.adapt = 1000,
      n.iter = sum(jags_dims[c("ni", "nb")]),
      n.thin = jags_dims["nt"],
      n.burnin = jags_dims["nb"],
      n.chains = jags_dims["nc"],
      parallel = T
    )
    
    if( length(warnings())>0){
      warning_list[[length(warning_list)+1]]<-warnings()
      names(warning_list)[[length(warning_list)]]<-mod_ind[j]
      assign("last.warning", NULL, envir = baseenv())  #this step clears the warnings
    }
    count<-count+1
    ### step 7: convergence diagnosis ###
    diag<-t(postpack::post_summ(post, c("B|b","sig"), neff = T, Rhat = T)[c("Rhat", "neff"),])
    converge_list[[count]]<-diag
    
    ### step 8: inference ###
    b_est = postpack::post_summ(post, "B|b",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.975))
    b2_est = postpack::post_summ(post, "B2|b2",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.975))
    b3_est = postpack::post_summ(post, "b3",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.95))
    
    spsu_load<-cbind(spsu_load,b2_est[,1])
    
    
    post_path<-post_subset(post,"B2",matrix=T)
    post_sst<-post_subset(post,"b3",matrix=T)
    
    
    post_probs[count,1:4]<-c(path,"spsu",length(which(post_path<0))/nrow(post_path),
                             length(which(post_sst<0))/nrow(post_sst))
    
    ppath_list[[count]]<-post_path
    
  }
  
  names(converge_list)[31:32]<-mod_ind
  
  names(ppath_list)[31:32]<-mod_ind
  
  pathogens<-c("arena2","c_b_cys","ce_sha","Circo.virus","de_sal","fa_mar",
               "fl_psy","ic_hof","ic_mul","IcD",
               "ku_thy","lo_sal","mo_vis","my_arc","nu_sal","pa_kab","pa_min",
               "pa_pse","pa_ther","prv","pspv","Qin","Rhabdo3_virus","rlo",
               "sch","smallUK","sp_des","te_bry","te_mar","ven","richness","rib")
  
  colnames(spsu_load)<-pathogens
  
  # collect Rhat and ess
  conv<-matrix(NA,nrow=32,ncol=8)
  for(i in 1:length(names(converge_list))){
    conv[i,1]<-converge_list[[i]][2,1]
    conv[i,2]<-converge_list[[i]][2,2]
    conv[i,3]<-converge_list[[i]][3,1]
    conv[i,4]<-converge_list[[i]][3,2]
    conv[i,5]<-converge_list[[i]][1,1]
    conv[i,6]<-converge_list[[i]][1,2]
    conv[i,7]<-converge_list[[i]][11,1]
    conv[i,8]<-converge_list[[i]][11,2]
    
  }
  
  colnames(conv)<-c("Rhat b0","ess b0","Rhat b1", "ess b1", "Rhat B2", "ess B2","Rhat b3", "ess b3")
  rownames(conv)<-names(converge_list)
  
  post_probs<-post_probs[1:32,]
  
  colnames(post_probs)<-c("taxa","season","pathogen_pp","sst_pp")
  
  sockeye_lw<-list(spsu_load,post_probs,converge_list,conv,ppath_list)
  names(sockeye_lw)<-c("estimates","posterior_probabilities","convergence_data","convergence_summary","samples")
  saveRDS(sockeye_lw,"sockeye_SW_lw_results_jan24_2023.rds")
  
