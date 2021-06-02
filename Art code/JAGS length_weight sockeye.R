library(postpack)
library(StatonMisc) # I don't think you actually need this package for this code (and might be true for some others)
library(tidyr)
library(dplyr)
library(Hmisc)
library(pals)
library(ggplot2)
library(gplots)
library(ggpubr)


d_sw<-readRDS("C:/Users/artie/OneDrive/PBS_postdoc/sockeye/sockeye_dataset_Nov27.rds")
sockCU<-read.csv("C:/Users/artie/OneDrive/PBS_postdoc/sockeye/sockeye_MGLvsCU_working_Oct20.csv")
d_sw$jaz<-sockCU$JAZ_ACRO[match(d_sw$Stock,sockCU$df_ID)]

# remove any fish without length or mass
d_sw<-d_sw %>%
	filter(!is.na(Length),
				 !is.na(Mass))

# checking for outliers
plot(log(d_sw$Mass)~log(d_sw$Length)) # looks like there are some
d_sw$Length[which(d_sw$Length<50)]<-d_sw$Length[which(d_sw$Length<50)]*10 # move decimal place for some

fit<-lm(log(d_sw$Mass)~log(d_sw$Length))			
d_sw$lwres<-fit$residuals
d_sw$month<-substr(d_sw$mo_yr,6,7)


# check pathogen prevalence to decide on what to look at
prevalence<-function(x){100*(length(which(x>0))/length(which(!is.na(x))))}
prev<-NULL
for(i in 39:80){
	a<-append(colnames(d_sw)[i],prevalence(d_sw[,i]))
	prev<-rbind(prev,a)
}

# here are all the pathogens we want to model
# Amy will have to decide upon this ultimately
pathogens<-c("arena2","c_b_cys","ce_sha","de_sal","fa_mar","fl_psy","ic_hof","ic_mul",
						 "ku_thy","lo_sal","mo_vis","my_arc","pa_kab","pa_min",
						 "pa_pse","pa_ther","prv","pspv","rlo",
						 "sch","smallUK","sp_des","te_bry","te_mar","ven")


# use this to loop
model_index<-paste0(rep(c("prev","prev","load","load"),25),sep="_",
                    rep(pathogens,4)[order(rep(pathogens,4))],sep="_",
                    rep(c("SpSu","FaWi"),25))

# remove these models due to insufficient positives
remove<-c("load_de_sal_FaWi","load_ku_thy_FaWi","load_mo_vis_FaWi","load_prv_FaWi","load_mo_vis_SpSu",
          "prev_mo_vis_FaWi","load_te_bry_FaWi")

mod_ind<-model_index[!model_index %in% remove]

spsu_prev<-NULL
fawi_prev<-NULL
spsu_load<-NULL
fawi_load<-NULL

warning_list<-list()
converge_list<-list()
count<-0
sst_results<-data.frame(matrix(NA,nrow=150,ncol=6))
post_probs<-data.frame(matrix(NA,nrow=150,ncol=5))

for(j in 1:length(mod_ind)){
	path_list<-list()
    
	  path<-substr(mod_ind[j],6,nchar(mod_ind[j])-5)
	  seas<-substr(mod_ind[j],nchar(mod_ind[j])-3,nchar(mod_ind[j]))
	  param<-substr(mod_ind[j],1,4)
    	
d_sw$pathogen<-d_sw[,match(path,colnames(d_sw))]
d_sw$present<-ifelse(d_sw$pathogen>0,1,0)
d_sw$month<-substr(d_sw$mo_yr,6,7)

# remove any NAs (already accounted for length and mass)
dat<-d_sw %>%
	filter(!is.na(jaz)) %>% 
  filter(season2==seas) %>% 
	filter(!is.na(pathogen)) %>%
  filter(!is.na(Length)) %>%
  filter(!is.na(Mass)) %>%
  filter(!is.na(temp_devs))

# need this step if doing load models
#if(param=="prev"){dat2<-dat}else{dat2<-dat %>% filter(pathogen>0)}
dat2<-dat

# I am dropping empty stocks so that the code works but this does mean that the same numbers won't
# correspond across models - this only matters if we decide to go in and look at stock effects
dat2$jaz<- droplevels(dat2$jaz)
# same would be true for years - those won't always match up across models

### step 1: data ###
jags_data = list(
	n_yrs = length(unique(dat2$Year)),  # I assume this is ocean year, not brood year
	n_obs = nrow(dat2),
	n_jaz = length(unique(dat2$jaz)), # gonna leave coding in as "jaz" because am lazy
	n_mon = length(unique(dat2$month)),
	yr = as.numeric(as.factor(dat2$Year)),
	jaz = as.numeric(as.factor(dat2$jaz)),
	mon = as.numeric(as.factor(dat2$month)),
	Mass = log(dat2$Mass),
	length = log(dat2$Length),
	pathogen = arm::rescale(if(param=="prev"){dat2$present}else{log(dat2$pathogen + min(dat2$pathogen[dat2$pathogen>0])/2)}),
	sst = arm::rescale(dat2$temp_devs)
	)
### step 2: model ###
jags_model = function() {
	# error variability priors
	sig_yr ~ dunif(0,2)
	sig_resid ~ dunif(0,2)
	sig_prior ~ dunif(0,2)
	sig_jaz ~ dunif(0,2)
	sig_mon ~ dunif(0,2)
	sig_B1 ~ dunif(0,2)
	
	# the random effect for year: influences the intercepts only
	for (y in 1:n_yrs) {
		yr_eff[y] ~ dnorm(0, 1/sig_yr^2)
	}
	# random effect for JAZ: influences intercept only
	for (s in 1:n_jaz){
		jaz_eff[s] ~ dnorm(0, 1/sig_jaz^2)
	}
	
	# random effect for month: intercept and slope for length 
	for(m in 1:n_mon){
	  mon_eff[m] ~ dnorm(0, 1/sig_mon^2)
	  b1[m] ~ dnorm(B1, 1/sig_B1^2)
	}
	
	#uninformed priors for betas
	b0 ~ dnorm(0, 1E-6)
	B1 ~ dnorm(0, 1E-6)
	b2 ~ dnorm(0, 1E-6)
	b3 ~ dnorm(0, 1E-6) 
	
	for (q in 1:n_obs) {
		Mass[q] ~ dnorm(Mass_hat[q], 1/sig_resid^2)
		Mass_hat[q] <- b0 + b1[mon[q]]*length[q] + b2*pathogen[q] + b3*sst[q] + yr_eff[yr[q]] + jaz_eff[jaz[q]] + mon_eff[mon[q]]
		residual[q] <- Mass[q] - Mass_hat[q]
	}
}

jags_file = "model.txt"
postpack::write_model(jags_model, jags_file)



### step 3: initial values ###
fit = with(jags_data, lme4::lmer(Mass ~ length + pathogen  + sst + (1|jaz) + (1|yr) + (1|mon)))
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
sig_mon_fit = 0.1
sig_resid_fit = summary(fit)$sigma

jags_inits = function(nc) {
	inits = list()
	for (c in 1:nc) {
		inits[[c]] = list(
			b0 = rnorm(1, b1_fit, abs(b1_fit) * 0.2),
			b1 = rnorm(jags_data$n_mon, b1_fit, abs(b1_fit) * 0.2),
			b2 = rnorm(1, b2_fit, abs(b2_fit) * 0.2),
			b3 = rnorm(1, b3_fit, abs(b3_fit) * 0.2),
			B1 = rnorm(1, b1_fit, 1e-3),
			sig_B1 = runif(1,0,1),
			sig_yr = runif(1, sig_yr_fit * 0.5, sig_yr_fit * 1.5),
			sig_jaz = runif(1, sig_jaz_fit * 0.5, sig_jaz_fit * 1.5),
			sig_mon = runif(1, sig_mon_fit * 0.5, sig_mon_fit * 1.5),
			sig_resid = runif(1, sig_resid_fit * 0.5, sig_resid_fit * 1.5)
		)
	}
	return(inits)
}

### step 4: parameters to monitor ###
jags_params = c("b0", "B1","b1","b2","b3", "sig_resid", "sig_yr", "sig_jaz", "yr_eff", "Mass_hat", "residual")

### step 5: dimensions ###
jags_dims = c(
	ni = 10000,  # number of post-burn-in samples per chain
	nb = 10000,  # number of burn-in samples
	nt = 2,     # thinning rate
	nc = 2      # number of chains
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
	parallel = F
)

 if( length(warnings())>0){
 	warning_list[[length(warning_list)+1]]<-warnings()
 	names(warning_list)[[length(warning_list)]]<-mod_ind[j]
 	assign("last.warning", NULL, envir = baseenv())  #this step clears the warnings
 }

### step 7: convergence diagnosis ###
diag<-t(postpack::post_summ(post, c("B|b", "sig"), ess = T, Rhat = T)[c("Rhat", "ess"),])
converge_list[[j]]<-diag

#postpack::diag_plots(post,"b")

### step 8: inference ###
b_est = postpack::post_summ(post, "b",p_summ=c(0.5,0.025,0.25,0.4,0.6,0.75,0.95))
b2_est = postpack::post_summ(post, "b2",p_summ=c(0.5,0.025,0.25,0.4,0.6,0.75,0.95))
b3_est = postpack::post_summ(post, "b3",p_summ=c(0.5,0.025,0.25,0.4,0.6,0.75,0.95))

mod_type = paste0(seas,sep="_",param)

if(mod_type=="SpSu_prev"){spsu_prev<-cbind(spsu_prev,b2_est)}else{
  if(mod_type=="SpSu_load"){spsu_load<-cbind(spsu_load,b2_est)}else{
    if(mod_type=="FaWi_prev"){fawi_prev<-cbind(fawi_prev,b2_est)}else{
      fawi_load<-cbind(fawi_load,b2_est)}}}

count<-count+1
sst_results[count,1:6]<-c(path,seas,param,b3_est[1,1],b3_est[4,1],b3_est[7,1])

post_path<-post_subset(post,"b2",matrix=T)
post_sst<-post_subset(post,"b3",matrix=T)

post_probs[count,1:5]<-c(path,seas,param,length(which(post_path<0))/nrow(post_path),
                         length(which(post_sst<0))/nrow(post_sst))



} # end of j loop

names(converge_list)<-mod_ind


colnames(spsu_prev)<-pathogens
colnames(fawi_prev)<-pathogens[-11]
colnames(spsu_load)<-pathogens[-11]
colnames(fawi_load)<-pathogens[-c(4,9,11,17,23)]




for_dwplot<-function(x){
	dw1<-as.data.frame(t(x))
	dw1$grp<-as.factor(rownames(dw1))
	dw1$grp<-factor(dw1$grp,levels=dw1$grp[order(dw1[,3])])
	colnames(dw1)[3:9]<-c("median","lci","mlci","clci","cuci","muci","uci")
	return(dw1)
}

ssp<-for_dwplot(spsu_prev)
ssl<-for_dwplot(spsu_load)
fwp<-for_dwplot(fawi_prev) # remove mo_vis
fwl<-for_dwplot(fawi_load)

names(post_probs)<-c("pathogen","season","parameter","path_pp","sst_pp")
post_probs2<-post_probs[1:count,]

names(sst_results)<-c("pathogen","season","parameter","sst_mean","lci","uci")
sst_results2<-sst_results[1:count,]

sockeye_lw<-list(ssp,ssl,fwp,fwl,converge_list,sst_results2,post_probs2)
names(sockeye_lw)<-c("ss_prev","ss_load","fw_prev","fw_load","convergence","sst_results","post_probs")
#saveRDS(sockeye_lw,"sockeye_lw_results_nov27_monthRE_sst_zeros.rds")

dw_plot_all<-function(x){
  ggplot(x,aes(x=mean,y=grp))+
    geom_vline(xintercept=0, color=2) + 
    geom_point(size=4) +
    labs(x="effect size",y=NULL) +
    geom_segment(aes(y=grp,yend=grp,x=lci, xend=uci),size=1) +
    geom_segment(aes(y=grp,yend=grp,x=mlci, xend=muci),size=2) +
    theme(legend.position="none", plot.title=element_text(size=10,face="bold")) +
    theme_bw()
}

ss_prev<-dw_plot_all(ssp) +ggtitle("spring-summer prevalence") 
ss_load<-dw_plot_all(ssl) +ggtitle("spring-summer load") 
fw_prev<-dw_plot_all(fwp) +ggtitle("fall-winter prevalence") 
fw_load<-dw_plot_all(fwl) +ggtitle("fall-winter load") 

ggarrange(ss_prev,ss_load,fw_prev,fw_load,nrow=1)



