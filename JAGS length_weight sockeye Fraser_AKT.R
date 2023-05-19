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


d_sw<-readRDS("data/sockeye_dataset_Nov27.rds")
sockCU<-read.csv("data/sockeye_MGLvsCU_working_Oct20.csv")
d_sw$jaz<-sockCU$JAZ_ACRO[match(d_sw$Stock,sockCU$df_ID)]

fraser<-c("FRCany+GStr","LFR+GStr","LILL+GStr","LTh+GStr","MFR+GStr","NTh+GStr",
          "STh+GStr","UFR+GStr")

# remove any fish without length or mass, not from fraser, not from spring-summer
d_sw<-d_sw %>%
	filter(!is.na(Length),
				 !is.na(Mass),
				 jaz %in% fraser,
				 season2=="SpSu")


# checking for outliers
plot(log(d_sw$Mass)~log(d_sw$Length)) # looks like there are some
d_sw$Length[which(d_sw$Length<50)]<-d_sw$Length[which(d_sw$Length<50)]*10 # move decimal place for some

fit<-lm(log(d_sw$Mass)~log(d_sw$Length))			
d_sw$lwres<-fit$residuals

## Calculate richness
# SW
rich_sw =
  data.frame(
    d_sw %>% 
      select(39:80) %>%
      mutate(richness = rowSums(.>0,na.rm=TRUE))) 
#rich_sw$rich_std <- scale(rich_sw$richness)

# Calculate RIB
#SW
ribdata_sw <- as.matrix(d_sw[,39:80])
rib1_sw <- data.frame(apply(ribdata_sw,2,function(x) {x/max(x, na.rm=TRUE)}))
rib2_sw <- data.frame(rowSums(rib1_sw, na.rm=TRUE))
richrib_sw <- cbind(rich_sw,rib2_sw)
colnames(richrib_sw)[44] <- "rib"
richrib_sw$rib[which(richrib_sw$rib==-Inf)] <- 0
#richrib_sw$rib_std <- scale(na.omit(richrib_sw$rib))
d_sw <- cbind(d_sw, richrib_sw[,43:44])
head(d_sw)

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
						 "ku_thy","lo_sal","my_arc","pa_kab","pa_min",
						 "pa_pse","pa_ther","prv","pspv","rlo",
						 "sch","smallUK","sp_des","te_bry","te_mar","ven","richness","rib")

model_index<-paste0("load",sep="_",pathogens,sep="_","SpSu")

# remove these models due to insufficient positives
#remove<-c("load_mo_vis_SpSu") # "prev_mo_vis_FaWi",

mod_ind<-model_index#[!model_index %in% remove]

spsu_load<-NULL
fawi_load<-NULL
warning_list<-list()
converge_list<-list()
count<-0
sst_results<-data.frame(matrix(NA,nrow=150,ncol=11))
post_probs<-data.frame(matrix(NA,nrow=150,ncol=4))
ppath_list<-list()
psst_list<-list()

pdf(file="sock_lwr_July14_21_AKT.pdf",width=11,height=8.5,paper="USr")

### Run it
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
    filter(!is.na(jaz),
           season2==seas,
           !is.na(pathogen),
           !is.na(Length),
           !is.na(Mass),
           jaz!="UFR+GStr") 
  # I am dropping empty stocks so that the code works but this does mean that the same numbers won't
  # correspond across models - this only matters if we decide to go in and look at stock effects
  dat2$jaz<- droplevels(dat2$jaz)
  # same would be true for years - those won't always match up across models
  
  ### step 1: data ###
  jags_data = list(
    n_yrs = length(unique(dat2$Year)),  # I assume this is ocean year, not brood year
    n_obs = nrow(dat2),
    n_jaz = length(unique(dat2$jaz)), # gonna leave coding in as "jaz" because am lazy
    yr = as.numeric(as.factor(dat2$Year)),
    jaz = as.numeric(as.factor(dat2$jaz)),
    Mass = log(dat2$Mass),
    length = log(dat2$Length),
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
    #sig_B3 ~ dunif(0,2)
    
    # power transform
    pow ~ dunif(0,2)
    
    # the random effect for year: influences the intercepts only
    for (y in 1:n_yrs) {
      yr_eff[y] ~ dnorm(0, 1/sig_yr^2)
    }
    # random effect for JAZ: influences intercept only
    for (s in 1:n_jaz){
      jaz_eff[s] ~ dnorm(0, 1/sig_jaz^2)
      b2[s] ~ dnorm(B2,1/sig_B2^2)
      #b3[s] ~ dnorm(B3,1/sig_B3^2)
    }
    
    #uninformed priors for betas
    b0 ~ dnorm(0, 1E-6)
    b1 ~ dnorm(0, 1E-6)
    B2 ~ dnorm(0, 1E-6)
    b3 ~ dnorm(0, 1E-6) 
    
    for (q in 1:n_obs) {
      Mass[q] ~ dnorm(Mass_hat[q], 1/sig_resid^2)
      Mass_hat[q] <- b0 + b1*length[q] + b2[jaz[q]]*(pathogen[q])^pow + b3*sst[q] + yr_eff[yr[q]] + jaz_eff[jaz[q]]
      residual[q] <- Mass[q] - Mass_hat[q]
    }
  }
  
  jags_file = "model.txt"
  postpack::write_model(jags_model, jags_file)
  
  
  
  ### step 3: initial values ###
  # put that log in there to make it something that will work for initial value
  fit = with(jags_data, lme4::lmer(Mass ~ length + log(pathogen+0.1)  + sst + (1|jaz) + (1|yr) ))
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
        b0 = rnorm(1, b0_fit, abs(b0_fit) * 0.1),
        b1 = rnorm(1, b1_fit, abs(b1_fit) * 0.2),
        b2 = rnorm(jags_data$n_jaz, b2_fit, abs(b2_fit) * 0.2),
        b3 = rnorm(1, b3_fit, abs(b3_fit) * 0.2),
        B2 = rnorm(1, b2_fit, 1e-3),
        #B3 = rnorm(1, b3_fit, 1e-3),
        pow=0.05,
        sig_B2 = runif(1,0,1),
        #sig_B3 = runif(1,0,1),
        sig_yr = runif(1, sig_yr_fit * 0.5, sig_yr_fit * 1.5),
        sig_jaz = runif(1, sig_jaz_fit * 0.5, sig_jaz_fit * 1.5),
        sig_resid = runif(1, sig_resid_fit * 0.5, sig_resid_fit * 1.5)
      )
    }
    return(inits)
  }
  
  ### step 4: parameters to monitor ###
  jags_params = c("b0", "b1","B2","b2","b3","pow", "sig_resid", "sig_yr", "sig_jaz", "yr_eff", "Mass_hat", "residual")
  
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
  diag<-t(postpack::post_summ(post, c("B|b", "sig"), neff = T, Rhat = T)[c("Rhat", "neff"),])
  converge_list[[j]]<-diag
  
  #postpack::diag_plots(post,"b")
  
  ### step 8: inference ###
  b_est = postpack::post_summ(post, "B|b",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.95))
  b2_est = postpack::post_summ(post, "B2|b2",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.95))
  b3_est = postpack::post_summ(post, "B3|b3",probs=c(0.5,0.025,0.25,0.4,0.6,0.75,0.95))
  
  power<-post_summ(post,"pow")
  
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
  
  
  p3<-ggplot(data=dat2,
             aes(x=pathogen^power[3],y=lwres,color=jaz,grp=jaz)) +
    geom_point() + geom_smooth(method="lm",se=T, color=1) + facet_wrap(~jaz,ncol=4,scales = "free") +
    #coord_cartesian(ylim=c(-0.25,0.25)) +
    theme_classic() + theme(legend.position = "none",strip.text.x = element_text(size = 8)) +
    ylab("raw residuals from length x mass") + xlab(paste("load, (copies RNA) ^",round(power[3],2)))
  
  b1<-post_subset(post,"B2|b2",matrix=T)
  
  b1<-as.data.frame(b1)
  
  colnames(b1)<-c("mean",as.character(unique(dat2$jaz)))
  
  b2<-pivot_longer(b1,col=1:ncol(b1),names_to ="stock",values_to = "est")
  b2$est<-as.numeric(b2$est)
  
  p4<-ggplot(data=b2[which(b2$stock!="mean"),],aes(x=est, group=stock, color=stock)) +
    geom_density(size=1.2,adjust=5) + theme_classic() +
    geom_density(data=b2[b2$stock=="mean",],col=1,size=2,adjust=5) +
    geom_vline(xintercept = 0,color=2) + coord_cartesian(xlim=c(min(b2_est[4,]),max(b2_est[9,]))) +
    theme(legend.position = "top", legend.key.size = unit(0.8,"lines")) +
    xlab("Effect size")
  
  g1<-annotate_figure(ggarrange(p3,p4),top=paste0(path,sep="  ",seas))
  print(g1)
  
} # end of j loop
dev.off()


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

sockeye_lw<-list(ssl,converge_list,sst_results2,post_probs2)
names(sockeye_lw)<-c("ss_load","convergence","sst_results","post_probs")
saveRDS(sockeye_lw,"sockeye_lw_results_July14.rds")

output <- sockeye_lw$ss_load
jpeg(filename='figs/Fig_MD_SW.jpg', 
     width=500, height=500, quality=300)
ggplot(output, aes(x=median, y=reorder(grp,-median))) +
  geom_point() +
  geom_errorbarh(aes(xmax = uci, xmin = lci), size=.2, height = 0) +
  geom_errorbarh(aes(xmax = muci, xmin = mlci), size=.5, height = 0) +
  geom_errorbarh(aes(xmax = cuci, xmin = clci), size=1.2, height = 0) +
  #xlim(-0.2,0.2) +
  geom_vline(xintercept = 0, lty = 2, size=.25) +
  labs(x="Effect size", y = "Freshwater-detected agents", title = "Freshwater") +
  theme(legend.position = "none", axis.text.y = element_text(face="italic"))
dev.off()

# dw_plot_all<-function(x){
#   ggplot(x,aes(x=mean,y=grp))+
#     geom_vline(xintercept=0, color=2) + 
#     geom_point(size=4) +
#     labs(x="effect size",y=NULL) +
#     geom_segment(aes(y=grp,yend=grp,x=lci, xend=uci),size=1) +
#     geom_segment(aes(y=grp,yend=grp,x=mlci, xend=muci),size=2) +
#     theme(legend.position="none", plot.title=element_text(size=10,face="bold")) +
#     theme_bw()
# }
# 
# ss_load<-dw_plot_all(ssl) +ggtitle("spring-summer load") 
# fw_load<-dw_plot_all(fwl) +ggtitle("fall-winter load") 
# 
# ggarrange(ss_load,fw_load,nrow=1)

# collect Rhat and neff
conv<-matrix(NA,nrow=44,ncol=8)
for(i in 1:length(names(converge_list))){
  conv[i,1]<-converge_list[[i]][3,1]
  conv[i,2]<-converge_list[[i]][3,2]
  conv[i,3]<-converge_list[[i]][4,1]
  conv[i,4]<-converge_list[[i]][4,2]
  conv[i,5]<-converge_list[[i]][1,1]
  conv[i,6]<-converge_list[[i]][1,2]
  conv[i,7]<-converge_list[[i]][2,1]
  conv[i,8]<-converge_list[[i]][2,2]
}

colnames(conv)<-c("Rhat b0","neff b0","Rhat b1", "neff b1", "Rhat B2", "neff B2", "Rhat B3", "neff B3")
rownames(conv)<-names(converge_list)


