# CHILKO ANALYSIS
# SURVIVAL RELATIVE TO INFECTION BURDENS

## This analysis is to focus on the Chilko stock in freshwater as most samples came from this 
# group and eliminating stock as a cofactor can increase our ability to detct pathogen related
# impacts. 


#### Load packages and set directory
#setwd("~/Documents.nosync/DFO PDF/Data/SSHI-sockeye")
library(lme4)
library(rstanarm) # https://mc-stan.org/users/documentation/case-studies/tutorial_rstanarm.html
library(ggplot2)
library(plotrix)
library(tidyverse)
library(gridExtra)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(shinystan)
library(data.table)
library(base)
library(ggpubr)
library(dplyr)

  
##### FW SR_resid data 
inf_agt_resid_data_fw <- read.csv("data/REDUCED_ONNE_agents only_FW_220126.csv")
head(inf_agt_resid_data_fw)
fw_rr <- read.csv("data/RR_ONNE_cumulative metrics_FW_220126.csv")
head(fw_rr)
inf_agt_resid_data_fw <- rbind(inf_agt_resid_data_fw, fw_rr)

# Data cleaning
## Standardize and incorporate into FW df
inf_std_fw <- plyr::ddply(inf_agt_resid_data_fw, c("agent"),function(x) {
  scaled_prev_fw <- scale(x$prev)
  scaled_load_fw <- scale(x$mean_load)
  xx <- data.frame(scaled_prev_fw, scaled_load_fw)
})
inf_agt_resid_data_fw$prev_std <- inf_std_fw[,2]
inf_agt_resid_data_fw$load_std <- inf_std_fw[,3]
dim(inf_agt_resid_data_fw)

# Replace RIB and richness with re-scaled values
inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="rib"),14] <- scale(inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="rib"),10])
inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="richness"),14] <- scale(inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="richness"),10])
inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="prv"),14] <- scale(inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="prv"),10])
inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="pspv"),14] <- scale(inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="pspv"),10])
inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="Rhabdo3_virus"),14] <- scale(inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="Rhabdo3_virus"),10])
inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="rlo"),14] <- scale(inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="rlo"),10])
inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="sch"),14] <- scale(inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="sch"),10])
inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="smallUK"),14] <- scale(inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="smallUK"),10])
inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="sp_des"),14] <- scale(inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="sp_des"),10])
inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="te_bry"),14] <- scale(inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="te_bry"),10])
inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="te_mar"),14] <- scale(inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="te_mar"),10])
inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="ven"),14] <- scale(inf_agt_resid_data_fw[which(inf_agt_resid_data_fw$agent=="ven"),10])

# Check on scaling
temp <- inf_agt_resid_data_fw[inf_agt_resid_data_fw$agent=="rib",]
ggplot(temp) +
  geom_point(aes(x=prev, y=prev_std), pch=21)
ggplot(inf_agt_resid_data_fw) +
  geom_point(aes(x=prev, y=prev_std, col=agent), pch=21) +
  stat_smooth(aes(x=prev, y=prev_std, col=agent), method="lm", se=F, linewidth=.5)

## Isolate chilko fish in fw
dim(inf_agt_resid_data_fw)
chilkofw <- inf_agt_resid_data_fw[inf_agt_resid_data_fw$Stock_Analysis=="Chilko",]
dim(chilkofw)  
head(chilkofw)  

## Bring in Chilko SR data
chSRdata <- read.csv("data/survival_indices_ONNE_chilkoFW.csv")
head(chSRdata)

## Merge FW Chilko pathogen data with FW SR data (long)
head(chilkofw)
chilkofw <- merge(chilkofw, chSRdata, by="brood_year", all.x=T)
names(chilkofw)
chilko <- chilkofw[,c(1,3,4,8,9,10,11,13,14,15,18,19)]
names(chilko)
chilko <- chilko[!duplicated(chilko), ]
head(chilko)

# Bring in temp deviation data (Chilko Lake dock, Chilcotin R)
temp.dev.chilko <- read.csv("data/FR_temp_dev.csv")

# Merge SST data with FR FW data
temp <- read.csv("data/master_brood_table_covar_210528.csv") #use brood table
temp2 <- temp[temp$stock_name=="Chilko",] #only chilko
temp2$Year <- temp2$BY+2
temp3 <- temp2[,c(35,3,33)]
temp4 <- merge(temp.dev.chilko, temp3, by = "Year", all.x=T)

#merge temp devs - use Chilcotin temps as more compete data set reflective of chilko lake conditions
temp5 <- temp4[,c(1,4,6)] #use only chilcotin R and SST
chilko <- merge(chilko, temp5, by = "Year")
head(chilko)
unique(chilko$Year)

# Remove mig year 2008 - insufficient data
chilko <- chilko[!chilko$Year==2008,]

#Create residual dataframes and remove missing values
## SR_resid_chfw
chilko_resid <- chilko[chilko$survival_index=="SR_resid_chfw",]
head(chilko_resid)
unique(chilko_resid$Year)

## RSM
chilko_rsm <- chilko[chilko$survival_index=="RSM",]
chilko_rsm
unique(chilko_rsm$Year)

# Create objects for Chilko analysis
years <- unique(chilko$Year)
brdyears <- unique(chilko$brood_year)
agents <- unique(chilko$agent)

## Prevalence Chilko FW
## Plot total detections of each agent by year 
samp.chilkofw<-chilko_resid %>% 
  group_by(agent, Year) %>%
  summarise(posdet = sum(posdet))
ggplot(data=samp.chilkofw, aes(x=factor(Year), y=posdet, fill=agent))+
  geom_bar(stat="identity")+
  labs(fill="Agents") +
  coord_flip()+
  xlab("Sampling year")+
  ylab("Positive detections")


# Prevalence Chilko by SR_resid_cov
## use SR_cov_resid
ggplot(chilko,aes(prev, value, color=survival_index, shape=factor(Year)))+
  geom_smooth(aes(prev, value, group=survival_index), method = "lm", se=F, linewidth=.2)+
  geom_point() +
  scale_shape_manual(values=1:nlevels(factor(chilko$Year))) +
  facet_wrap(~ agent,nrow=5, scales = "free")+
  xlab("FW prevalence")+
  ylab("SR residual")+
  theme_bw()

# Plot prev by year
ggplot(data=chilko_resid, aes(x=brood_year, y=prev_std, col=agent)) +
  geom_line()+
  geom_point()

## Are the 2 metrics correlated?
temp <- chilko[,c(1,11,12)]
temp2 <- temp %>% distinct()
temp3 <- spread(temp2, survival_index, value) #go long to wide
ggplot(data=temp3, aes(x=RSM, y=SR_resid_chfw, col=factor(Year))) +
  geom_point()


#######################################################################################
# STAN Approach for Multi-level Modeling

# FW PREVALENCE - INDEPENDENT MODELS by AGENT
## Stock-specific metric
### Create files for each agent
for(i in unique(chilko_resid$agent)) {
  nam <- paste("df.chilko", i, sep = ".")
  assign(nam, chilko_resid[chilko_resid$agent==i,])
}

### Loop for STAN independent models by agent
for(i in agents){
  data <- subset(chilko_resid, agent==i)
  nam <- paste("mod.chilko", i, sep = ".")
  assign(nam, stan_glm(value ~ 0 + prev_std + chilcotinjundev,
                        data = data,
                        adapt_delta=0.99))
}

## Derive coefficient estimates and save in .csv file
coefs_stan_chilko <- matrix(NA,
                     nrow = length(agents),
                     ncol = 5,
                     dimnames = list(agents,c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.chilko",i, sep="."))
  ind_coef <- summary(model, 
                      pars = c("prev_std"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan_chilko[i,] <- ind_coef[1,c(4:8)]
}
write.csv(coefs_stan_chilko, file="data/prev_coefs_stan_ChilkoFW_220906_tempcov_SRresidch.csv")


# Load estimates from file (if not running full model) and assign rownames
coefs_stan_chilko <- read.csv("data/prev_coefs_stan_ChilkoFW_220906_tempcov_SRresidch.csv")
rownames(coefs_stan_chilko) <- coefs_stan_chilko[,1]
coefs_stan_chilko <- coefs_stan_chilko[,-1]  

# Plot effect size per agent
coefs_order <- coefs_stan_chilko[order(-coefs_stan_chilko[,3]),]

jpeg(filename='figs/Fig_prev_coefs_stan_chilko_221206_tempcov_SRresidch.jpg', 
     width=480, height=500, quality=75)
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))
plotCI(x = coefs_order[,3],
       y = seq(1,length(agents)),
       li = (coefs_order[,1]),
       ui = (coefs_order[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-6,4),
       pch = 16,
       scol = "grey")
plotCI(x = coefs_order[,3],
       y = seq(1,length(agents)),
       li = (coefs_order[,2]),
       ui = (coefs_order[,4]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "grey")
text(rep(-6,length(agents)), 
     seq(1,length(agents)), 
     labels = rownames(coefs_order), 
     pos = 4,
     font = 2,
     cex=0.95)
axis(1, at = c(-8, -6, -4, -2, 0, 2, 4, 6, 8))
abline(v = 0, lty = 2)
box(col="grey")	
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Smolt migrants per spawner",3,line=0.25)
dev.off()

####################################################
# Examine fit
plot(mod.chilko.c_b_cys, prob = 0.5)
plot(mod.chilko.arena2, prob = 0.5, pars = "beta")
plot(mod.chilko.arena2, "hist", pars = "sigma")
pp_check(mod.chilko.arena2, plotfun = "error_binned")


####################################################
## Derive posterior estimates 

### Slopes
coefs_stan_stk_slp.ch <- matrix(NA,
                                nrow = length(agents),
                                ncol = 5,
                                dimnames = list(agents,c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.chilko.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("prev_std"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan_stk_slp.ch[i,] <- ind_coef[1,c(4:8)]
}

### Temperature slopes
coefs_stan_stk_slp_temp.ch <- matrix(NA,
                                nrow = length(agents),
                                ncol = 5,
                                dimnames = list(agents,c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.chilko.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("chilcotinjundev"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan_stk_slp_temp.ch[i,] <- ind_coef[1,c(4:8)]
}

### Sigmas
coefs_stan_stk_sig.ch <- matrix(NA,
                                nrow = length(agents),
                                ncol = 5,
                                dimnames = list(agents, c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.chilko.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("sigma"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_sig.ch", i, sep = ".")
  coefs_stan_stk_sig.ch[i,] <- ind_coef[1,c(4:8)]
}

### Rhat and Neff
sims.ch <-as.matrix(mod.chilko.c_b_cys) #extract parameter names from any agent model 
dim(sims.ch)
para_name.ch <- c(colnames(sims.ch), "mean_PPD", "log-posterior")
para_name.ch

#### Rhat
coefs_stan_stk_rhat.ch <- matrix(NA,
                                 nrow = length(para_name.ch),
                                 ncol = 21,
                                 dimnames = list(para_name.ch, agents))
for(i in agents){
  model<-get(paste("mod.chilko",i, sep="."))
  ind_coef <- as.matrix(summary(model, 
                                digits = 3) [,"Rhat"])
  nam <- paste("coefs_stan_stk_rhat.ch", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_rhat.ch <- cbind(ind_coef[c(1:dim(ind_coef)[1]),1], paste(i))))
}

#### Neff
coefs_stan_stk_neff.ch <- matrix(NA,
                                 nrow = length(para_name.ch),
                                 ncol = 2,
                                 dimnames = list(para_name.ch,c("Neff","agent")))
for(i in agents){
  model<-get(paste("mod.chilko",i, sep="."))
  ind_coef <- as.matrix(summary(model) [,"n_eff"])
  nam <- paste("coefs_stan_stk_neff.ch", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_neff.ch <- cbind(ind_coef[c(1:dim(ind_coef)[1]),1], paste(i), paste("Neff"))))
}




## Derive proportions of posterior draws <0 per model - FW prevalence
param.ch<-colnames(data.frame(mod.chilko.c_b_cys)) #create object of parameters in model
param.prop0.ch <- matrix(NA,
                         ncol = 1,
                         nrow = length(agents),
                         dimnames = list(agents,"prop0"))

for (i in agents){
  model<-as.matrix(get(paste("mod.chilko.",i, sep="")))
  model2<-model[2001:4000,1]
  param.prop0.ch[i,] <- (sum(model2 < 0))/2000
}
write.csv(param.prop0.ch, file="data/Percent post draws >0_prev_chilko_220915_SRresidch.csv")
propzero.ch <- read.csv("data/Percent post draws >0_prev_chilko_220915_SRresidch.csv")

jpeg(filename='figs/Fig_percent post draws >0_prev_chilko_220915_SRresidch.jpg', 
     width=480, height=500, quality=75)
ggplot(propzero.ch) +
  geom_bar(stat="identity", aes(reorder(X, prop0), prop0), fill="gray", col="black") +
  ylim(0,1)+
  labs(x="Agent", y="Proportion of total", title="Chilko Posterior estimates <0") +
  coord_flip()
dev.off()




#################### REPEAT WITH RSM ###################################################
#######################################################################################
# STAN Approach for Multi-level Modeling

# FW PREVALENCE - INDEPENDENT MODELS by AGENT
## Stock-specific metric
### Create files for each agent
for(i in unique(chilko_rsm$agent)) {
  nam <- paste("df.chilko", i, sep = ".")
  assign(nam, chilko_rsm[chilko_rsm$agent==i,])
}

### Loop for STAN independent models by agent
for(i in agents){
  data <- subset(chilko_rsm, agent==i)
  nam <- paste("mod.chilko", i, sep = ".")
  assign(nam, stan_glm(value ~ 0 + prev_std + early_sst_stnd,
                       data = data,
                       adapt_delta=0.99))
}

## Derive coefficient estimates and save in .csv file
coefs_stan_chilko.rsm <- matrix(NA,
                            nrow = length(agents),
                            ncol = 5,
                            dimnames = list(agents,c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.chilko",i, sep="."))
  ind_coef <- summary(model, 
                      pars = c("prev_std"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan_chilko.rsm[i,] <- ind_coef[1,c(4:8)]
}
write.csv(coefs_stan_chilko.rsm, file="data/prev_coefs_stan_ChilkoFW_220906_tempcov_RSM.csv")


# Load estimates from file (if not running full model) and assign rownames
coefs_stan_chilko.rsm <- read.csv("data/prev_coefs_stan_ChilkoFW_220906_tempcov_RSM.csv")
rownames(coefs_stan_chilko.rsm) <- coefs_stan_chilko.rsm[,1]
coefs_stan_chilko.rsm <- coefs_stan_chilko.rsm[,-1]  


# Plot effect size per agent
coefs_order <- coefs_stan_chilko.rsm[order(-coefs_stan_chilko.rsm[,3]),]

jpeg(filename='figs/Fig_prev_coefs_stan_chilko_221206_tempcov_RSM.jpg', 
     width=480, height=500, quality=75)
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))
plotCI(x = coefs_order[,3],
       y = seq(1,length(agents)),
       li = (coefs_order[,1]),
       ui = (coefs_order[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-3,4),
       pch = 16,
       scol = "grey")
plotCI(x = coefs_order[,3],
       y = seq(1,length(agents)),
       li = (coefs_order[,2]),
       ui = (coefs_order[,4]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "grey")
text(rep(-3,length(agents)), 
     seq(1,length(agents)), 
     labels = rownames(coefs_order), 
     pos = 4,
     font = 2,
     cex=0.95)
axis(1, at = c(-2, 0, 2, 4))
abline(v = 0, lty = 2)
box(col="grey")	
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Returning adults per smolt migrant",3,line=0.25)
dev.off()

####################################################
# Examine fit
plot(mod.chilko.c_b_cys, prob = 0.5)
plot(mod.chilko.arena2, prob = 0.5, pars = "beta")
plot(mod.chilko.arena2, "hist", pars = "sigma")
pp_check(mod.chilko.arena2, plotfun = "error_binned")


####################################################
## Derive posterior estimates 

### Slopes
coefs_stan_stk_slp.chrsm <- matrix(NA,
                                nrow = length(agents),
                                ncol = 5,
                                dimnames = list(agents,c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.chilko.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("prev_std"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan_stk_slp.chrsm[i,] <- ind_coef[1,c(4:8)]
}

### Temperature slopes
coefs_stan_stk_slp_temp.chrsm <- matrix(NA,
                                     nrow = length(agents),
                                     ncol = 5,
                                     dimnames = list(agents,c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.chilko.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("chilcotinjundev"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan_stk_slp_temp.chrsm[i,] <- ind_coef[1,c(4:8)]
}

### Sigmas
coefs_stan_stk_sig.chrsm <- matrix(NA,
                                nrow = length(agents),
                                ncol = 5,
                                dimnames = list(agents, c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.chilko.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("sigma"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_sig.ch", i, sep = ".")
  coefs_stan_stk_sig.chrsm[i,] <- ind_coef[1,c(4:8)]
}

### Rhat and Neff
sims.ch <-as.matrix(mod.chilko.c_b_cys) #extract parameter names from any agent model 
dim(sims.ch)
para_name.ch <- c(colnames(sims.ch), "mean_PPD", "log-posterior")
para_name.ch

#### Rhat
coefs_stan_stk_rhat.chrsm <- matrix(NA,
                                 nrow = length(para_name.ch),
                                 ncol = 21,
                                 dimnames = list(para_name.ch, agents))
for(i in agents){
  model<-get(paste("mod.chilko",i, sep="."))
  ind_coef <- as.matrix(summary(model, 
                                digits = 3) [,"Rhat"])
  nam <- paste("coefs_stan_stk_rhat.ch", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_rhat.chrsm <- cbind(ind_coef[c(1:dim(ind_coef)[1]),1], paste(i))))
}

#### Neff
coefs_stan_stk_neff.chrsm <- matrix(NA,
                                 nrow = length(para_name.ch),
                                 ncol = 2,
                                 dimnames = list(para_name.ch,c("Neff","agent")))
for(i in agents){
  model<-get(paste("mod.chilko",i, sep="."))
  ind_coef <- as.matrix(summary(model) [,"n_eff"])
  nam <- paste("coefs_stan_stk_neff.chrsm", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_neff.chrsm <- cbind(ind_coef[c(1:dim(ind_coef)[1]),1], paste(i), paste("Neff"))))
}




## Derive proportions of posterior draws <0 per model - FW prevalence
param.ch<-colnames(data.frame(mod.chilko.c_b_cys)) #create object of parameters in model
param.prop0.ch <- matrix(NA,
                         ncol = 1,
                         nrow = length(agents),
                         dimnames = list(agents,"prop0"))

for (i in agents){
  model<-as.matrix(get(paste("mod.chilko.",i, sep="")))
  model2<-model[2001:4000,1]
  param.prop0.ch[i,] <- (sum(model2 < 0))/2000
}
write.csv(param.prop0.ch, file="data/Percent post draws >0_prev_chilko_220915_RSM.csv")
propzero.ch <- read.csv("data/Percent post draws >0_prev_chilko_220915_RSM.csv")

jpeg(filename='figs/Fig_percent post draws >0_prev_chilko_220915_RSM.jpg', 
     width=480, height=500, quality=75)
ggplot(propzero.ch) +
  geom_bar(stat="identity", aes(reorder(X, prop0), prop0), fill="gray", col="black") +
  ylim(0,1)+
  labs(x="Agent", y="Proportion of total", title="Chilko Posterior estimates <0") +
  coord_flip()
dev.off()

###########################################################################






##Post hoc analysis - prevalence in FW
## Prevalence Chilko FW
## Plot total SW detections of each agent by year - note variable prevalence across agents and years
chilko_resid <- chilko[chilko$fw_metric=="SR_resid_chfw",]
samp.chilkofw<-chilko_resid %>% 
  group_by(agent, Year) %>%
  summarise(posdet = sum(posdet))

#All agents
ggplot(data=samp.chilkofw, aes(x=reorder(agent, posdet), y=posdet, fill=factor(Year)))+
  geom_bar(stat="identity")+
  labs(fill="Sampling year") +
  coord_flip()+
  xlab("Infectious agents")+
  ylab("Positive detections")

#ARENA2
samp.chilkofw.arena <- samp.chilkofw[samp.chilkofw$agent=="arena2",]
ggplot(data=samp.chilkofw.arena, aes(x=reorder(agent, posdet), y=posdet, fill=factor(Year)))+
  geom_bar(stat="identity")+
  labs(fill="Sampling year") +
  coord_flip()+
  xlab("Infectious agents")+
  ylab("Positive detections")



vec3 <- unique(prev_beta_ch$X)
blank.agents.ch <- setdiff(vec1, vec3) 
tempch <- data.frame(
  "X" = blank.agents.ch,
  "lower" = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
  "X25" = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
  "mid" = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
  "X75" = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
  "upper" = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
samps.ch <- read.csv("figs/Assay_table_Chilko_ONNE2021_220916.csv")
samps.ch$SWFW <- "CH"
prev.ch <- samps.ch[,c(1,5)]
prev_beta_ch <- merge(prev_beta_ch, prev.ch, by="X", all.x=TRUE)

#CH
temp <- read.csv("data/Percent post draws >0_prev_chilko_220915.csv")
prev_neg_ch <- data.frame(temp[,1:2])
tempch <- data.frame(
  "X" = blank.agents.ch,
  "prop0" = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
prev_neg_ch <- rbind(prev_neg_ch,tempch)

#CH
prev_beta_ch <-prev_beta_ch[order(-prev_beta_ch$mid),]
colnames(prev_beta_ch) <- c("agent","lwr.ch","midlwr.ch","mid.ch","midupp.ch","upp.ch","CH_prev")
prev_neg_ch <- prev_neg_ch[order(match(prev_neg_ch$X, prev_beta_ch$agent)),]
colnames(prev_neg_ch) <- c("agent","prop0.ch")

all.ch.red <- all.ch[!all.ch$CH_prev<1,]
all.ch.red <- all.ch.red[!all.ch.red$prop0.ch=="NA",]
all.ch.red <- all.ch.red[!all.ch.red$agent=="NA",]
#remove from fw: Circo.virus,ihnv,pa_kab,pch_sal,prv,te_mar,ven,smallUK

plot.chilko <- all.ch.red %>%
  ggplot(aes(x=mid.ch, y=reorder(plot.agent,-mid.ch), alpha = CH_prev)) +
  geom_point() +
  geom_errorbarh(aes(xmax = upp.ch, xmin = lwr.ch), size=.5, height = 0) +
  geom_errorbarh(aes(xmax = midupp.ch, xmin = midlwr.ch), size=1, height = 0) +
  #xlim(-1,1) +
  geom_vline(xintercept = 0, lty = 2, size=.25) +
  labs(x="", y = "", title = "Chilko Freshwater") +
  theme(legend.position = "none", axis.text.y = element_text(face="italic")) +
  geom_text(aes(x=1.5, y=reorder(plot.agent,-mid.ch), label = round((100*prop0.ch),0)),inherit.aes = FALSE) +
  geom_text(aes(x=-1.5, y=reorder(plot.agent,-mid.ch), label = CH_prev, alpha=CH_prev))

jpeg(filename='figs/Fig_slope_prev_Chilko_220916.jpg', 
     width=500, height=500, quality=300)
chilko.plot
dev.off()




#######################################################################################
# STAN Approach for Multi-level Modeling

# FW PREVALENCE - INDEPENDENT MODELS by AGENT
## Stock-specific metric
### Create files for each agent
for(i in unique(chilko_rsm$agent)) {
  nam <- paste("df.chilko", i, sep = ".")
  assign(nam, chilko_rsm[chilko_rsm$agent==i,])
}

### Loop for STAN independent models by agent
for(i in agents){
  data <- subset(chilko_rsm, agent==i)
  nam <- paste("mod.chilko", i, sep = ".")
  assign(nam, stan_glm(value ~ 0 + prev_std + chilcotinjundev,
                       data = data,
                       adapt_delta=0.99))
}

## Derive coefficient estimates and save in .csv file
coefs_stan_chilko <- matrix(NA,
                            nrow = length(agents),
                            ncol = 5,
                            dimnames = list(agents,c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.chilko",i, sep="."))
  ind_coef <- summary(model, 
                      pars = c("prev_std"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan_chilko[i,] <- ind_coef[1,c(4:8)]
}
write.csv(coefs_stan_chilko, file="data/prev_coefs_stan_ChilkoFW_220906_tempcov_RSM.csv")


# Load estimates from file (if not running full model) and assign rownames
coefs_stan_chilko <- read.csv("data/prev_coefs_stan_ChilkoFW_220906_tempcov_RSM.csv")
rownames(coefs_stan_chilko) <- coefs_stan_chilko[,1]
coefs_stan_chilko <- coefs_stan_chilko[,-1]  

# Plot effect size per agent
coefs_order <- coefs_stan_chilko[order(-coefs_stan_chilko[,3]),]

jpeg(filename='figs/Fig_prev_coefs_stan_chilko_221206_tempcov_RSM.jpg', 
     width=480, height=500, quality=75)
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))
plotCI(x = coefs_order[,3],
       y = seq(1,length(agents)),
       li = (coefs_order[,1]),
       ui = (coefs_order[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-10,10),
       pch = 16,
       scol = "grey")
plotCI(x = coefs_order[,3],
       y = seq(1,length(agents)),
       li = (coefs_order[,2]),
       ui = (coefs_order[,4]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "grey")
text(rep(-10,length(agents)), 
     seq(1,length(agents)), 
     labels = rownames(coefs_order), 
     pos = 4,
     font = 2,
     cex=0.95)
axis(1, at = c(-8, -6, -4, -2, 0, 2, 4, 6, 8))
abline(v = 0, lty = 2)
box(col="grey")	
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Prevalence",3,line=0.25)
dev.off()

