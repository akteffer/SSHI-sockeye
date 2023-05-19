# Chilko SW prevalence ~ Marine survival


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
inf_agt_resid_data <- read.csv("data/REDUCED_ONNE_agents only_SW_220126.csv")
head(inf_agt_resid_data)
sw_rr <- read.csv("data/RR_ONNE_cumulative metrics_SW_220126.csv")
head(sw_rr)
inf_agt_resid_data <- rbind(inf_agt_resid_data, sw_rr)

# Data cleaning
# Standardize and incorporate into SW df
inf_std <- plyr::ddply(inf_agt_resid_data, c("agent"),function(x) {
  scaled_prev <- scale(x$prev)
  scaled_load <- scale(x$mean_load)
  xx <- data.frame(scaled_prev, scaled_load)
})
inf_agt_resid_data$prev_std <- inf_std[,2]
inf_agt_resid_data$load_std <- inf_std[,3]

# Replace RIB and richness with re-scaled values
inf_agt_resid_data[which(inf_agt_resid_data$agent=="rib"),14] <- scale(inf_agt_resid_data[which(inf_agt_resid_data$agent=="rib"),10])
inf_agt_resid_data[which(inf_agt_resid_data$agent=="richness"),14] <- scale(inf_agt_resid_data[which(inf_agt_resid_data$agent=="richness"),10])
inf_agt_resid_data[which(inf_agt_resid_data$agent=="prv"),14] <- scale(inf_agt_resid_data[which(inf_agt_resid_data$agent=="prv"),10])
inf_agt_resid_data[which(inf_agt_resid_data$agent=="pspv"),14] <- scale(inf_agt_resid_data[which(inf_agt_resid_data$agent=="pspv"),10])
inf_agt_resid_data[which(inf_agt_resid_data$agent=="Rhabdo3_virus"),14] <- scale(inf_agt_resid_data[which(inf_agt_resid_data$agent=="Rhabdo3_virus"),10])
inf_agt_resid_data[which(inf_agt_resid_data$agent=="re_sal"),14] <- scale(inf_agt_resid_data[which(inf_agt_resid_data$agent=="re_sal"),10])
inf_agt_resid_data[which(inf_agt_resid_data$agent=="rlo"),14] <- scale(inf_agt_resid_data[which(inf_agt_resid_data$agent=="rlo"),10])
inf_agt_resid_data[which(inf_agt_resid_data$agent=="sch"),14] <- scale(inf_agt_resid_data[which(inf_agt_resid_data$agent=="sch"),10])
inf_agt_resid_data[which(inf_agt_resid_data$agent=="smallUK"),14] <- scale(inf_agt_resid_data[which(inf_agt_resid_data$agent=="smallUK"),10])
inf_agt_resid_data[which(inf_agt_resid_data$agent=="sp_des"),14] <- scale(inf_agt_resid_data[which(inf_agt_resid_data$agent=="sp_des"),10])
inf_agt_resid_data[which(inf_agt_resid_data$agent=="te_bry"),14] <- scale(inf_agt_resid_data[which(inf_agt_resid_data$agent=="te_bry"),10])
inf_agt_resid_data[which(inf_agt_resid_data$agent=="te_mar"),14] <- scale(inf_agt_resid_data[which(inf_agt_resid_data$agent=="te_mar"),10])
inf_agt_resid_data[which(inf_agt_resid_data$agent=="ven"),14] <- scale(inf_agt_resid_data[which(inf_agt_resid_data$agent=="ven"),10])


# Check on scaling
#temp <- inf_agt_resid_data[inf_agt_resid_data$agent=="ven",]
#ggplot(temp) +
#  geom_point(aes(x=prev, y=prev_std), pch=21)
#ggplot(inf_agt_resid_data) +
#  geom_point(aes(x=prev, y=prev_std, col=agent), pch=21)+
#  stat_smooth(aes(x=prev, y=prev_std, col=agent), method="lm", se=F, size=.5)

## Isolate chilko fish in sw
dim(inf_agt_resid_data)
chilkosw <- inf_agt_resid_data[inf_agt_resid_data$Stock_Analysis=="Chilko",]
chilkosw <- chilkosw[chilkosw$metric=="SR_cov_resid",]
dim(chilkosw)  
head(chilkosw)  

## Bring in Chilko SR data
chSRdata <- read.csv("data/survival_indices_ONNE_chilkoFW.csv")
chSRdata <- chSRdata[chSRdata$survival_index=="RSM",]
head(chSRdata)

## Merge SW Chilko pathogen data with FW SR data (long)
head(chilkosw)
chilkosw <- merge(chilkosw, chSRdata, by="brood_year", all.x=T)
names(chilkosw)
chilko <- chilkosw[,c(1,3,4,8,9,10,11,13,14,15,19,7)]
names(chilko)
colnames(chilko) <- c("brood_year","Stock_Analysis" ,"Year" ,    
                        "N" , "posdet" ,"prev" ,"mean_load" ,   
                       "agent", "prev_std","load_std", "RSM","SR_cov_resid")
names(chilko)
chilko <- chilko[!duplicated(chilko), ]
head(chilko)
dim(chilko)

# Bring in temp deviation data (Chilko Lake dock, Chilcotin R)
#temp.dev.chilko <- read.csv("data/FR_temp_dev.csv")

# Merge SST data 
temp <- read.csv("data/master_brood_table_covar_210528.csv") #use brood table
temp2 <- temp[temp$stock_name=="Chilko",] #only chilko
temp2$Year <- temp2$BY+2
temp3 <- temp2[,c(35,33)]
chilko <- merge(chilko, temp3, by = "Year", all.x=T)
head(chilko)
unique(chilko$Year)

# Remove mig year 2008 - insufficient data
chilkoW <- chilko[!chilko$Year==2008,]
head(chilkoW)

# Make it long
chilko <- gather(chilkoW, metric, value, RSM:SR_cov_resid, factor_key=T)

# Create objects for Chilko analysis
years <- unique(chilko$Year)
brdyears <- unique(chilko$brood_year)
agents <- unique(chilko$agent)

# Prevalence Chilko by SR_resid_cov
## 
jpeg(filename='figs/Fig_prev_coefs_stan_chilko_SW_RSM_SRcovresid_230220.jpg', 
     width=1000, height=1000, quality=300)
ggplot(chilko,aes(prev, value, color=metric, shape=factor(Year)))+
  geom_smooth(aes(prev, value, group=metric), method = "lm", se=F, linewidth=.2)+
  geom_point() +
  scale_shape_manual(values=1:nlevels(factor(chilko$Year))) +
  facet_wrap(~ agent,nrow=5, scales = "free")+
  xlab("SW prevalence")+
  ylab("SR residual")+
  theme_bw()
dev.off()

# Plot prev by year
ggplot(data=chilko_resid, aes(x=brood_year, y=prev_std, col=agent)) +
  geom_line()+
  geom_point()

## Prevalence Chilko FW
## Plot total detections of each agent by year 
chilko_rsm <- chilko[chilko$metric=="RSM",]
chilko_resid <- chilko[chilko$metric=="SR_cov_resid",]

samp.chilkosw<-chilko_rsm %>% 
  group_by(agent, Year) %>%
  summarise(posdet = sum(posdet))
ggplot(data=samp.chilkosw, aes(x=factor(Year), y=posdet, fill=agent))+
  geom_bar(stat="identity")+
  labs(fill="Agents") +
  coord_flip()+
  xlab("Sampling year")+
  ylab("Positive detections")



#######################################################################################
# STAN Approach for Multi-level Modeling - RSM
# SW PREVALENCE - INDEPENDENT MODELS by AGENT
## Stock-specific metric

# Create objects for Chilko analysis
years <- unique(chilko$Year)
brdyears <- unique(chilko$brood_year)
agents <- unique(chilko$agent)

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
write.csv(coefs_stan_chilko, file="data/prev_coefs_stan_ChilkoSW_RSM.csv")


# Load estimates from file (if not running full model) and assign rownames
coefs_stan_chilko <- read.csv("data/prev_coefs_stan_ChilkoSW_RSM.csv")
rownames(coefs_stan_chilko) <- coefs_stan_chilko[,1]
coefs_stan_chilko <- coefs_stan_chilko[,-1]  

# Plot effect size per agent
coefs_order <- coefs_stan_chilko[order(-coefs_stan_chilko[,3]),]

jpeg(filename='figs/Fig_prev_coefs_stan_chilkoSW_RSM.jpg', 
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
       xlim = c(-3,3),
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
axis(1, at = c(-8, -6, -4, -2, 0, 2, 4, 6, 8))
abline(v = 0, lty = 2)
box(col="grey")	
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Marine survival vs SW prev",3,line=0.25)
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
                      regex_pars = c("early_sst_stnd"),
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
                                 nrow = length(agents),
                                 ncol = 5,
                                 dimnames = list(agents, para_name.ch))
for(i in agents){
  model<-get(paste("mod.chilko",i, sep="."))
  ind_coef <- as.matrix(summary(model, 
                                digits = 3) [,"Rhat"])
  coefs_stan_stk_rhat.ch[i,] <- ind_coef[c(1:5),1]
}

#### Neff
coefs_stan_stk_neff.ch <- matrix(NA,
                                 nrow = length(agents),
                                 ncol = 5,
                                 dimnames = list(agents, para_name.ch))
for(i in agents){
  model<-get(paste("mod.chilko",i, sep="."))
  ind_coef <- as.matrix(summary(model) [,"n_eff"])
  coefs_stan_stk_neff.ch[i,] <- ind_coef[c(1:5),1]
}

## Matrices to DFs add label for metric
coefs_stan_stk_slp.ch <- data.frame(coefs_stan_stk_slp.ch)
coefs_stan_stk_slp_temp.ch <- data.frame(coefs_stan_stk_slp_temp.ch)
coefs_stan_stk_sig.ch <- data.frame(coefs_stan_stk_sig.ch)
coefs_stan_stk_rhat.ch <- data.frame(coefs_stan_stk_rhat.ch)
coefs_stan_stk_neff.ch <- data.frame(coefs_stan_stk_neff.ch)
coefs_stan_stk_slp.ch$metric <- "prev_std"
coefs_stan_stk_slp_temp.ch$metric <- "SST"
coefs_stan_stk_sig.ch$metric <- "sigma"
coefs_stan_stk_rhat.ch$metric <- "Rhat"
coefs_stan_stk_neff.ch$metric <- "Neff"

## Merge files for Chilko marine survival and save as csv
chilko_betas_sigma <- rbind(coefs_stan_stk_slp.ch,coefs_stan_stk_slp_temp.ch,coefs_stan_stk_sig.ch)
chilko_rhat_neff <- rbind(coefs_stan_stk_rhat.ch,coefs_stan_stk_neff.ch)
write.csv(chilko_betas_sigma, file="data/Chilko_MarineSurvival_betassigma_230222.csv")
write.csv(chilko_rhat_neff, file="data/Chilko_MarineSurvival_rhatneff_230222.csv")

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
write.csv(param.prop0.ch, file="data/Percent post draws >0_prev_chilko_RSM.csv")
propzero.ch <- read.csv("data/Percent post draws >0_prev_chilko_RSM.csv")

jpeg(filename='figs/Fig_percent post draws >0_prev_chilko_RSM.jpg', 
     width=480, height=500, quality=75)
ggplot(propzero.ch) +
  geom_bar(stat="identity", aes(reorder(X, prop0), prop0), fill="gray", col="black") +
  ylim(0,1)+
  labs(x="Agent", y="Proportion of total", title="Chilko Posterior estimates <0") +
  coord_flip()
dev.off()





#######################################################################################
# STAN Approach for Multi-level Modeling - SR_cov_resid
# SW PREVALENCE - INDEPENDENT MODELS by AGENT
## Stock-specific metric

# Create objects for Chilko analysis
years <- unique(chilko$Year)
brdyears <- unique(chilko$brood_year)
agents <- unique(chilko$agent)

### Create files for each agent
for(i in unique(chilko_resid$agent)) {
  nam <- paste("df.chilko", i, sep = ".")
  assign(nam, chilko_resid[chilko_resid$agent==i,])
}

### Loop for STAN independent models by agent
for(i in agents){
  data <- subset(chilko_resid, agent==i)
  nam <- paste("mod.chilko", i, sep = ".")
  assign(nam, stan_glm(value ~ 0 + prev_std + early_sst_stnd,
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
write.csv(coefs_stan_chilko, file="data/prev_coefs_stan_ChilkoSW_SRcovresid.csv")


# Load estimates from file (if not running full model) and assign rownames
coefs_stan_chilko <- read.csv("data/prev_coefs_stan_ChilkoSW_SRcovresid.csv")
rownames(coefs_stan_chilko) <- coefs_stan_chilko[,1]
coefs_stan_chilko <- coefs_stan_chilko[,-1]  

# Plot effect size per agent
coefs_order <- coefs_stan_chilko[order(-coefs_stan_chilko[,3]),]

jpeg(filename='figs/Fig_prev_coefs_stan_chilkoSW_SRcovresid.jpg', 
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
       xlim = c(-3,3),
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
axis(1, at = c(-8, -6, -4, -2, 0, 2, 4, 6, 8))
abline(v = 0, lty = 2)
box(col="grey")	
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("SRcovresid vs SW prev",3,line=0.25)
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
write.csv(param.prop0.ch, file="data/Percent post draws >0_prev_chilko_SWprev_SRcovresid.csv")
propzero.ch <- read.csv("data/Percent post draws >0_prev_chilko_SWprev_SRcovresid.csv")

jpeg(filename='figs/Fig_percent post draws >0_prev_chilko_SWprev_SRcovresid.jpg', 
     width=480, height=500, quality=75)
ggplot(propzero.ch) +
  geom_bar(stat="identity", aes(reorder(X, prop0), prop0), fill="gray", col="black") +
  ylim(0,1)+
  labs(x="Agent", y="Proportion of total", title="Chilko Posterior estimates <0") +
  coord_flip()
dev.off()


