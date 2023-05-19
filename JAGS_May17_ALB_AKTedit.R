library(postpack)
library(StatonMisc) # I don't think you actually need this package for this code (and might be true for some others)
library(tidyr)
library(dplyr)
library(Hmisc)
library(pals)
library(ggplot2)
library(gplots)
library(ggpubr)
library(scales)

pathogen_names<-read.csv("~/Documents.nosync/DFO PDF/Data/SSHI-sockeye/data/pathogen_names.csv")

# new data source
d_sw<-readRDS("~/Documents.nosync/DFO PDF/Data/SSHI-sockeye/sockeye_dataset_Aug27AB.rds")
sockCU<-read.csv("~/Documents.nosync/DFO PDF/Data/SSHI-sockeye/data/sockeye_MGLvsCU_working_Oct20.csv")
d_sw$jaz<-sockCU$JAZ_ACRO[match(d_sw$Stock,sockCU$df_ID)]
fraser<-c("FRCany+GStr","LFR+GStr","LILL+GStr","LTh+GStr","MFR+GStr","NTh+GStr",
          "STh+GStr","UFR+GStr")

# remove any fish without length or mass, not from fraser, not from spring-summer
d_sw<-d_sw %>%
	filter(!is.na(Length),
				 !is.na(Mass),
				 jaz %in% fraser,
				 Latitude<=51.5,
				 season2=="SpSu") 

# replace -999 with NA
d_sw[d_sw==-999]<-NA

# merge in IcD and other new agents
d<-read.csv("data/ONNE metadata no LOD_6.17.2021.csv")
d2<-d %>%
  filter(SWFW=="SW")
d_sw <- merge(d_sw, d2[,c(1,37,48,55,64,65,71,76,78,79,81,82,84,87,94)], by="Unique", all.x = TRUE)
d_sw <- d_sw %>% relocate(c(83:96), .before = season2) #move in line with other agents

# checking for outliers
plot(log(d_sw$Mass)~log(d_sw$Length)) # looks like there are some
d_sw$Length[which(d_sw$Length<50)]<-d_sw$Length[which(d_sw$Length<50)]*10 # move decimal place for some

# create residuals
fit<-lm(log(d_sw$Mass)~log(d_sw$Length))			
d_sw$lwres<-fit$residuals

# there are a bunch of ultra-skinny fish with few pathogen detections (may be biasing a number of relationships)
# they seem abnormally skinny (perhaps measurement error?) - here am removing anything outside of 3 SDs from the mean to deal with them
sdx3<- which(d_sw$lwres > (sd(d_sw$lwres) * 3) | d_sw$lwres < -(sd(d_sw$lwres) * 3))
d_sw<-d_sw[-sdx3,]

## Calculate richness - sum of agents per fish w/ positive detections
# SW
rich_sw =
  data.frame(
    d_sw %>% 
      select(c(28:84)) %>%
      mutate(richness = rowSums(.>0,na.rm=TRUE))) 

# Calculate RIB - Artie Bass et al. composite richness+load metric
#SW
ribdata_sw <- as.matrix(d_sw[,c(28:84)])
rib1_sw <- data.frame(apply(ribdata_sw,2,function(x) {x/max(x, na.rm=TRUE)}))
rib2_sw <- data.frame(rowSums(rib1_sw, na.rm=TRUE))
richrib_sw <- cbind(rich_sw,rib2_sw)
colnames(richrib_sw)[59] <- "rib"
richrib_sw$rib[which(richrib_sw$rib==-Inf)] <- 0
d_sw <- cbind(d_sw, richrib_sw[,58:59])
head(d_sw)
d_sw <- d_sw %>% relocate(c("richness","rib"), .before = season2) #move these in line with agents

# check pathogen prevalence to decide on what to look at
prevalence<-function(x){100*(length(na.omit(x[x>0]))/length(na.omit(x)))}
prev<-NULL
for(i in 28:86){
  a<-append(colnames(d_sw)[i],prevalence(d_sw[,i]))
  prev<-rbind(prev,a)
}
prev

# Positive detections - just the raw integers
posdet<-function(x){length(na.omit(x))}
posdet.df <- NULL
for(i in 28:86){
  a<-append(colnames(d_sw)[i],posdet(d_sw[,i]))
  posdet.df<-rbind(posdet.df,a)
}
posdet.df


# Here are the pathogens we want to model based on prevalence
pathogens<-c("richness",  "rib",  "c_b_cys",  "IcD",  "lo_sal", "pa_ther",  "my_arc",  "pa_min",
             "pspv",  "ce_sha",  "ic_mul",  "fa_mar",  "arena2",  "te_mar", "pa_kab",  "sch",  "ven",  
             "pa_pse",  "ic_hof",  "rlo",  "te_bry",  "sp_des",  "smallUK",
             "fl_psy",  "prv",  "de_sal")

# Get the model components ready
model_index<-paste0("load",sep="_",pathogens,sep="_","SpSu")
mod_ind<-model_index[!model_index %in% remove]
spsu_load<-NULL
fawi_load<-NULL
warning_list<-list()
converge_list<-list()
count<-0
sst_results<-data.frame(matrix(NA,nrow=150,ncol=11))
post_probs<-data.frame(matrix(NA,nrow=150,ncol=4))
ppath_list<-list()
psst_list<-list()

pdf(file="data/sock_lwr_Dec22_sw.pdf",width=11,height=8.5,paper="USr") #this did not work

#########################################################################
#Run it
for(j in 1:length(mod_ind)){
	path_list<-list()
    
	  path<-substr(mod_ind[j],6,nchar(mod_ind[j])-5)
	  seas<-substr(mod_ind[j],nchar(mod_ind[j])-3,nchar(mod_ind[j]))
	  param<-substr(mod_ind[j],1,4)
    	
d_sw$pathogen<-d_sw[,match(path,colnames(d_sw))]
d_sw$present<-ifelse(d_sw$pathogen>0,1,0)
d_sw$month<-substr(d_sw$mo_yr,6,7)

# remove any NAs (already accounted for length and mass)
dat2<-d_sw %>%
	filter(!is.na(jaz)) %>% 
  filter(season2==seas) %>% 
	filter(!is.na(pathogen)) %>%
  filter(!is.na(Length)) %>%
  filter(!is.na(Mass)) %>%
  filter(!is.na(temp_devs))

## NOTE - I did not drop empty stocks as it tripped the code
# I am dropping empty stocks so that the code works but this does mean that the same numbers won't
# correspond across models - this only matters if we decide to go in and look at stock effects
#dat2$jaz<- droplevels(dat2$jaz)
# same would be true for years - those won't always match up across models

### step 1: data ###
jags_data = list(
	n_yrs = length(unique(dat2$ocean_yr)),  
	n_obs = nrow(dat2),
	n_jaz = length(unique(dat2$jaz)),   # gonna leave coding in as "jaz" because am lazy
	yr = as.numeric(as.factor(dat2$ocean_yr)),
	jaz = as.numeric(as.factor(dat2$jaz)),
	Mass = log(dat2$Mass),
	length = log(dat2$fl),
	pathogen = dat2$pathogen,
	sst = arm::rescale(dat2$temp_devs)
	)

### step 2: model ###
jags_model = function() {
	# error variability priors
	sig_yr ~ dunif(0,2)
	sig_resid ~ dunif(0,2)
	sig_prior ~ dunif(0,2)
	sig_jaz ~ dunif(0,2)
	sig_B2 ~ dunif(0,2)
	sig_B0 ~ dunif(0,2)

	# power transform
	pow ~ dunif(0,2)
	
	# the random effect for year: influences the intercepts only
	for (y in 1:n_yrs) {
		yr_eff[y] ~ dnorm(0, 1/sig_yr^2)
	}
	# random effect for JAZ: influences intercept only
	for (s in 1:n_jaz){
	  #jaz_eff[s] ~ dnorm(0, 1/sig_jaz^2)
	  b0[s] ~ dnorm(B0,1/sig_B0^2)
	  b2[s] ~ dnorm(B2,1/sig_B2^2)
	  	}
	
	#uninformed priors for betas
	B0 ~ dnorm(0, 1E-6)
	b1 ~ dnorm(0, 1E-6)
	B2 ~ dnorm(0, 1E-6)
	b3 ~ dnorm(0, 1E-6) 
	
	for (q in 1:n_obs) {
		Mass[q] ~ dnorm(Mass_hat[q], 1/sig_resid^2)
		Mass_hat[q] <- b0[jaz[q]] + b1*length[q] + b2[jaz[q]]*(pathogen[q])^pow + b3*sst[q] + yr_eff[yr[q]] 
		residual[q] <- Mass[q] - Mass_hat[q]
	}
}

jags_file = "model.txt"
postpack::write_model(jags_model, jags_file)

### step 3: initial values ###
# put that log in there to make it something that will work for initial value
fit = with(jags_data, lme4::lmer(Mass ~ length + pathogen  + sst + (1|jaz) + (1|yr) ))
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
		  b0 = rnorm(jags_data$n_jaz, b0_fit, abs(b0_fit) * 0.2),
		  b1 = rnorm(1, b1_fit, abs(b1_fit) * 0.2),
		  b2 = rnorm(jags_data$n_jaz, b2_fit, abs(b2_fit) * 0.2),
		  b3 = rnorm(1, b3_fit, abs(b3_fit) * 0.2),
		  B2 = rnorm(1, b2_fit, 1e-3),
		  B0 = rnorm(1, b0_fit, 1e-3),
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
jags_params = c("b0", "b1","B0","B2","b2","b3","pow", "sig_resid", "sig_yr", "sig_jaz", "yr_eff", "Mass_hat", "residual")

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
	parallel = F
)

 if( length(warnings())>0){
 	warning_list[[length(warning_list)+1]]<-warnings()
 	names(warning_list)[[length(warning_list)]]<-mod_ind[j]
 	assign("last.warning", NULL, envir = baseenv())  #this step clears the warnings
 }

### step 7: convergence diagnosis ###
diag<-t(postpack::post_summ(post, c("B|b", "sig"), neff = T, Rhat = T)[c("Rhat", "neff"),])
converge_list[[j]]<-diag

### step 8: inference ###
b_est = postpack::post_summ(post, "B|b", probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.95))
b2_est = postpack::post_summ(post, "B2|b2",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.95))
b3_est = postpack::post_summ(post, "B3|b3",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.95))

pwr<-post_summ(post,"pow")

if(seas=="SpSu"){spsu_load<-cbind(spsu_load,b2_est[,1])}else{
  fawi_load<-cbind(fawi_load,b2_est[,1])}

count<-count+1
sst_results[count,1:11]<-c(path,seas,b3_est[,1])

post_path<-post_subset(post,"B2",matrix=T)
post_sst<-post_subset(post,"b3",matrix=T)

post_probs[count,1:4]<-c(path,seas,length(which(post_path<0))/nrow(post_path),
                         length(which(post_sst<0))/nrow(post_sst))

ppath_list[[count]]<-post_path
psst_list[[count]]<-post_sst


bf<-post_subset(post,"b1",matrix=T)
bf<-as.data.frame(bf)

b0<-post_subset(post,"b0",matrix=T)
b0<-as.data.frame(b0)
colnames(b0)<-c(as.character(unique(dat2$jaz)))

b1<-post_subset(post,"B2|b2",matrix=T)
b1<-as.data.frame(b1)
colnames(b1)<-c("mean",as.character(unique(dat2$jaz)))

grps<-length(unique(dat2$jaz))

rand<-sample(0:30000,1000,replace = F)

samps<-data.frame(pivot_longer(b0[rand,],cols=1:grps,names_to = "cwt1",values_to = "int"),
                  pivot_longer(b1[rand,2:(grps+1)],cols=1:grps,names_to = "cwt",values_to = "slp"))
samps$flb1<-rep(bf[rand,1],grps)
samps$new_int<-samps$int-(mean(jags_data$Mass)-mean(jags_data$length)*samps$flb1)

samps2<-samps[,2:6]; samps2$pathogen<-NA; samps2$lw_res<-NA
datx<-data.frame(int=NA,cwt=dat2$jaz,
                 slp=NA,flb1=NA,new_int=NA,
                 pathogen=dat2$pathogen,lw_res=dat2$lwres)

samps3<-rbind(samps2,datx)

p<-ggplot(samps3,aes(x=pathogen^pwr[3],y=lw_res,fill=cwt)) + 
  geom_abline(aes(intercept = new_int, slope = slp),
              color="gray30",alpha=0.1,size=0.2) +
  geom_point(pch=21,size=2) +
  facet_wrap(~cwt,ncol=3) + theme_classic() + labs(tag="C)") +
  xlab(paste0("pathogen load (copies RNA)",sep=" ^ ",round(pwr[3],2))) + 
  ylab("mass deviation") +
  theme(legend.position = "none",strip.text.x = element_text(size = 8))

hexcodes<-hue_pal()(grps)

g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
stripr2<-stripr[-which(unlist(lapply(g$grobs[stripr],grepl,pattern="zero")))]

to_ord<-order(substr(g$layout$name[stripr2],11,11),substr(g$layout$name[stripr2],9,9))
fills <- hexcodes
k <- 1
for (i in stripr2[to_ord]) {
  z <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[z]]$gp$fill <- fills[k]
  k <- k+1
}
p3<-g

b2dat<-pivot_longer(b1,names_to = "JAZ",values_to = "est", cols=2:(grps+1))

b2<-pivot_longer(b1,col=1:ncol(b1),names_to ="stock",values_to = "est")
b2$est<-as.numeric(b2$est)

p4<-ggplot(data=b2[which(b2$stock!="mean"),],aes(x=est, group=stock, color=stock)) +
  geom_density(size=1.2,adjust=5) + theme_classic() +
  geom_density(data=b2[b2$stock=="mean",],col=1,size=2,adjust=5) +
  geom_vline(xintercept = 0,color=2) + coord_cartesian(xlim=c(min(b2_est[4,]),max(b2_est[9,]))) +
  theme(legend.position = "none", legend.key.size = unit(0.8,"lines")) +
  xlab("slope coefficient for pathogen load")

g1<-annotate_figure(ggarrange(p3,p4),
                    top=paste0(pathogen_names[match(path,pathogen_names$abbrev),2],sep="  ",seas))

print(g1)

} # end of j loop
dev.off()
###############################################################################################

# Explore the data
names(converge_list)<-mod_ind
colnames(spsu_load)<-pathogens

for_dwplot<-function(x){
	dw1<-as.data.frame(t(x))
	dw1$grp<-as.factor(rownames(dw1))
	dw1$grp<-factor(dw1$grp,levels=dw1$grp[order(dw1[,3])])
	colnames(dw1)[3:9]<-c("median","lci","mlci","clci","cuci","muci","uci")
	return(dw1)
}
ssl<-for_dwplot(spsu_load)
names(post_probs)<-c("pathogen","season","path_pp","sst_pp")
post_probs2<-post_probs[1:count,]
names(sst_results)<-c("pathogen","season","parameter","sst_mean","lci","uci")
sst_results2<-sst_results[1:count,]

# collect Rhat and neff
conv<-matrix(NA,nrow=26,ncol=8)
for(i in 1:length(names(converge_list))){
  conv[i,1]<-converge_list[[i]]["B0",1]
  conv[i,2]<-converge_list[[i]]["B0",2]
  conv[i,3]<-converge_list[[i]]["b1",1]
  conv[i,4]<-converge_list[[i]]["b1",2]
  conv[i,5]<-converge_list[[i]]["B2",1]
  conv[i,6]<-converge_list[[i]]["B2",2]
  conv[i,7]<-converge_list[[i]]["b3",1]
  conv[i,8]<-converge_list[[i]]["b3",2]
}

colnames(conv)<-c("Rhat B0","neff B0","Rhat b1", "neff b1", 
                  "Rhat B2", "neff B2", "Rhat b3", "neff b3")
rownames(conv)<-names(converge_list)
sockeye_lw<-list(ssl,converge_list,sst_results2,post_probs2,conv)
names(sockeye_lw)<-c("ss_load","convergence","sst_results","post_probs","convergence_summary")

saveRDS(sockeye_lw,"sockeye_lw_results_sw_Dec23_51.5.rds")

# Plot
output <- sockeye_lw$ss_load
ggplot(output, aes(x=median, y=reorder(grp,-median))) +
  geom_point() +
  geom_errorbarh(aes(xmax = uci, xmin = lci), size=.2, height = 0) +
  geom_errorbarh(aes(xmax = muci, xmin = mlci), size=.5, height = 0) +
  geom_errorbarh(aes(xmax = cuci, xmin = clci), size=1.2, height = 0) +
  #xlim(-0.2,0.2) +
  geom_vline(xintercept = 0, lty = 2, size=.25) +
  labs(x="Effect size", y = "Marine-detected agents", title = "Marine") +
  theme(legend.position = "none", axis.text.y = element_text(face="italic"))


