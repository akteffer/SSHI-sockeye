# SSHI ANALYSIS SOCKEYE STAGE II
# A.K. Teffer
#### Sockeye salmon productivity versus infection profiles

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

#### Read in data, clean, standardize - stock-specific metrics
##### SW
inf_agt_resid_data <- read.csv("data/ONNE productivity infection analysis_stockspecific_210519.csv")
head(inf_agt_resid_data)
##### FW
inf_agt_resid_data_fw <- read.csv("data/ONNE productivity infection analysis_stockspecific_FW_210519.csv")
head(inf_agt_resid_data_fw)

# Data cleaning
# Standardize and incorporate into SW df
inf_std <- plyr::ddply(inf_agt_resid_data, c("agent"),function(x) {
  scaled_prev <- scale(x$prev)
  scaled_load <- scale(x$mean_load)
  xx <- data.frame(scaled_prev, scaled_load)
})
inf_agt_resid_data$prev_std <- inf_std[,2]
inf_agt_resid_data$load_std <- inf_std[,3]
# Add Stock column, remove smallUK, add plot.agent column (full name)
inf_agt_resid_data$Stock <- inf_agt_resid_data$Stock_Analysis
inf_agt_resid_data$plot.agent <- inf_agt_resid_data$agent
inf_agt_resid_data$plot.agent <- recode(inf_agt_resid_data$plot.agent, "ic_mul" = "I. multifiliis", 
                                 "te_mar" = "T. maritinum",
                                 "pa_ther" = "P. theridion",
                                 "fl_psy" = "F. psychrophilum",
                                 "sch" = "Ca. S. salmonis",
                                 "te_bry" = "T. bryosalmonae",
                                 "pa_kab" = "P. kabatai",
                                 "c_b_cys" = "Ca. B. cysticola",
                                 "pa_min" = "P. minibicornis",
                                 "arena2" = "SPAV-2",
                                 "fa_mar" = "F. margolisi",
                                 "my_arc" = "M. arcticus",
                                 "ven" = "VENV",
                                 "ic_hof" = "I. hoferi",
                                 "lo_sal" = "L. salmonae",
                                 "rlo" = "RLO",
                                 "sp_des" = "S. destruens",
                                 "ku_thy" = "K. thyrsites",
                                 "prv" = "PRV",
                                 "pspv" = "PSPV",
                                 "ce_sha" = "C. shasta",
                                 "pa_pse" = "P. pseudobranchicola",
                                 "de_sal" = "D. salmonis")
dim(inf_agt_resid_data)

# Standardize and incorporate into FW df
inf_std_fw <- plyr::ddply(inf_agt_resid_data_fw, c("agent"),function(x) {
  scaled_prev_fw <- scale(x$prev)
  scaled_load_fw <- scale(x$mean_load)
  xx <- data.frame(scaled_prev_fw, scaled_load_fw)
})
inf_agt_resid_data_fw$prev_std <- inf_std_fw[,2]
inf_agt_resid_data_fw$load_std <- inf_std_fw[,3]
dim(inf_agt_resid_data_fw)

# Create objects for analysis
agents <- unique(inf_agt_resid_data$agent)
stocks <- unique(inf_agt_resid_data$Stock)
years <- unique(inf_agt_resid_data$Year)
brdyears <- unique(inf_agt_resid_data$brood_year)

## Plot total detections by stock
## SW - When were stocks sampled?
sample.totals.sw<-aggregate(inf_agt_resid_data$N, by=list(Category=inf_agt_resid_data$Stock, Year=inf_agt_resid_data$Year), FUN=sum)
jpeg(filename='figs/Fig_Total samples by stock year_SW.jpg', 
     width=480, height=500, quality=75)
ggplot(data=sample.totals.sw, aes(x=reorder(Category, x), y=x, fill=factor(Year))) +
  geom_col() +
  coord_flip() +
  xlab("Stocks")+
  ylab("Samples")
dev.off()

## FW - When were stocks sampled?
sample.totals.fw<-aggregate(inf_agt_resid_data_fw$N, by=list(Category=inf_agt_resid_data_fw$Stock, Year=inf_agt_resid_data_fw$Year), FUN=sum)
jpeg(filename='figs/Fig_Total samples by stock year_FW.jpg', 
     width=480, height=500, quality=75)
ggplot(data=sample.totals.fw, aes(x=reorder(Category, x), y=x, fill=factor(Year))) +
  geom_col() +
  coord_flip() +
  xlab("Stocks")+
  ylab("Samples")
dev.off()

## Plot total SW detections of each agent by year - note variable prevalence across agents and years
samplesperagent.sw<-inf_agt_resid_data %>% 
  group_by(agent, Year) %>%
  summarise(posdet = sum(posdet))
#plot
jpeg(filename='figs/Fig_Total agent detections by year_SW.jpg', 
     width=480, height=500, quality=75)
ggplot(data=samplesperagent.sw, aes(x=reorder(agent, posdet), y=posdet, fill=factor(Year)))+
  geom_bar(stat="identity")+
  labs(fill="Sampling year") +
  coord_flip()+
  xlab("Infectious agents")+
  ylab("Positive detections")+
  scale_x_discrete(labels=c("ic_mul" = expression(italic("I. multifiliis")), 
                            "te_mar" = expression(italic("T. maritinum")),
                            "pa_ther" = expression(italic("P. theridion")),
                            "fl_psy" = expression(italic("F. psychrophilum")),
                            "sch" = expression(italic("Ca. S. salmonis")),
                            "te_bry" = expression(italic("T. bryosalmonae")),
                            "pa_kab" = expression(italic("P. kabatai")),
                            "c_b_cys" = expression(italic("Ca. B. cysticola")),
                            "pa_min" = expression(italic("P. minibicornis")),
                            "arena2" = "SPAV-2",
                            "fa_mar" = expression(italic("F. margolisi")),
                            "my_arc" = expression(italic("M. arcticus")),
                            "ven" = "VENV",
                            "ic_hof" = expression(italic("I. hoferi")),
                            "lo_sal" = expression(italic("L. salmonae")),
                            "rlo" = "RLO",
                            "sp_des" = expression(italic("S. destruens")),
                            "ku_thy" = expression(italic("K. thyrsites")),
                            "prv" = "PRV",
                            "pspv" = "PSPV",
                            "ce_sha" = expression(italic("C. shasta")),
                            "pa_pse" = expression(italic("P. pseudobranchicola")),
                            "de_sal" = expression(italic("D. salmonis"))))
dev.off()

## Plot total FW detections of each agent by year - note variable prevalence across agents and years
samplesperagent.fw<-inf_agt_resid_data_fw %>% 
  group_by(agent, Year) %>%
  summarise(posdet = sum(posdet)) 
#plot
jpeg(filename='figs/Fig_Total agent detections by year_FW.jpg', 
     width=480, height=500, quality=75)
ggplot(data=samplesperagent.fw, aes(x=reorder(agent, posdet), y=posdet, fill=factor(Year)))+
  geom_bar(stat="identity") +
  labs(fill="Sampling year") +
  coord_flip()+
  xlab("Infectious agents")+
  ylab("Positive detections")+
  scale_x_discrete(labels=c("ic_mul" = expression(italic("I. multifiliis")), 
                            "te_mar" = expression(italic("T. maritinum")),
                            "pa_ther" = expression(italic("P. theridion")),
                            "fl_psy" = expression(italic("F. psychrophilum")),
                            "sch" = expression(italic("Ca. S. salmonis")),
                            "te_bry" = expression(italic("T. bryosalmonae")),
                            "pa_kab" = expression(italic("P. kabatai")),
                            "c_b_cys" = expression(italic("Ca. B. cysticola")),
                            "pa_min" = expression(italic("P. minibicornis")),
                            "arena2" = "SPAV-2",
                            "fa_mar" = expression(italic("F. margolisi")),
                            "my_arc" = expression(italic("M. arcticus")),
                            "ven" = "VENV",
                            "ic_hof" = expression(italic("I. hoferi")),
                            "lo_sal" = expression(italic("L. salmonae")),
                            "rlo" = "RLO",
                            "sp_des" = expression(italic("S. destruens")),
                            "ku_thy" = expression(italic("K. thyrsites")),
                            "prv" = "PRV",
                            "pspv" = "PSPV",
                            "ce_sha" = expression(italic("C. shasta")),
                            "pa_pse" = expression(italic("P. pseudobranchicola")),
                            "de_sal" = expression(italic("D. salmonis"))))
dev.off()


# DATA EXPLORATION AND VISUALIZATION
## Agent diversity by stock
### SW
samplesperagent.swst<-inf_agt_resid_data %>% 
  group_by(agent, Stock_Analysis) %>%
  summarise(det = sum(posdet)) %>%
  mutate(freq = det / sum(det)) #proportion of pos detections of total pos detecttions in this Stock

jpeg(filename='figs/Fig_detections by stock_proportion.jpg', 
     width=500, height=300, quality=75)
ggplot(data=samplesperagent.swst, aes(x=reorder(Stock_Analysis, det), y=det, fill=agent))+
  geom_col( position="fill") +
  labs(fill="Sampling year") +
  coord_flip() +
  xlab("Stocks")+
  ylab("Proportion of detections")+
  scale_y_continuous(labels = scales::percent)
dev.off()

### FW
samplesperagent.fwst<-inf_agt_resid_data_fw %>% 
  group_by(agent, Stock_Analysis) %>%
  summarise(det = sum(posdet)) %>%
  mutate(freq = det / sum(det)) #proportion of pos detections of total pos detecttions in this Stock

jpeg(filename='figs/Fig_detections by stock_proportion_FW.jpg', 
     width=430, height=300, quality=75)
ggplot(data=samplesperagent.fwst, aes(x=reorder(Stock_Analysis, det), y=det, fill=agent))+
  geom_col(position="fill") +
  labs(fill="Sampling year") +
  coord_flip() +
  xlab("Stocks")+
  ylab("Proportion of detections")+
  scale_y_continuous(labels = scales::percent)
dev.off()

### By sampling year in FW (only 5: 2011-2015)
samplesperagent.fwsy<-inf_agt_resid_data_fw %>% 
  group_by(agent, Year) %>%
  summarise(det = sum(posdet)) %>%
  mutate(freq = det / sum(det)) #proportion of pos detections of total pos detecttions in this Stock
ggplot(data=samplesperagent.fwsy, aes(x=factor(Year), y=det, fill=agent))+
  geom_col() +
  labs(fill="Sampling year") +
  coord_flip() +
  xlab("Year")+
  ylab("Proportion of detections")+
  scale_y_continuous(labels = scales::percent)


# Plot raw data by: 
## Prevalence SW
jpeg(filename='figs/Fig_Raw data by year_prev_stspec.jpg', 
     width=480, height=500, quality=75)
ggplot(inf_agt_resid_data,aes(prev, resid_value, color=Stock, shape=factor(Year)))+
  geom_smooth(aes(prev, resid_value, group=Stock), method = "lm", se=F, size=.2)+
  geom_point()+
  scale_shape_manual(values=1:nlevels(factor(inf_agt_resid_data$Year))) +
  facet_wrap(~ agent,nrow=5, scales = "free")+
  xlab("SW prevalence")+
  ylab("SR residual")+
  theme_bw()
dev.off()

## Load SW
jpeg(filename='figs/Fig_Raw data by year_load_stspec.jpg', 
     width=480, height=500, quality=75)
ggplot(inf_agt_resid_data,aes(log10(mean_load), resid_value, color=Stock_Analysis, shape=factor(Year)))+
  geom_smooth(aes(log10(mean_load), resid_value, group=Stock), method = "lm", se=F, size=.2)+
  geom_point()+
  scale_shape_manual(values=1:nlevels(factor(inf_agt_resid_data$Year))) +
  facet_wrap(~ agent,nrow=5, scales="free")+
  xlab("SW log10 load")+
  ylab("SR residual")+
  theme_bw()
dev.off()


## Prevalence FW
jpeg(filename='figs/Fig_Raw data by year_prev_stspec_FW.jpg', 
     width=480, height=500, quality=75)
ggplot(inf_agt_resid_data_fw,aes(prev, resid_value, color=Stock_Analysis, shape=factor(Year)))+
  geom_smooth(aes(prev, resid_value, group=Stock_Analysis), method = "lm", se=F, size=.2)+
  geom_point()+
  scale_shape_manual(values=1:nlevels(factor(inf_agt_resid_data_fw$Year))) +
  facet_wrap(~ agent,nrow=5, scales = "free")+
  xlab("FW prevalence")+
  ylab("SR residual")+
  theme_bw()
dev.off()

## Load FW - not enough for model
jpeg(filename='figs/Fig_Raw data by year_load_stspec_FW.jpg', 
     width=480, height=500, quality=75)
ggplot(inf_agt_resid_data_fw,aes(log10(mean_load), resid_value, color=Stock_Analysis, shape=factor(Year)))+
  geom_smooth(aes(log10(mean_load), resid_value, group=Stock_Analysis), method = "lm", se=F, size=.2)+
  geom_point()+
  scale_shape_manual(values=1:nlevels(factor(inf_agt_resid_data_fw$Year))) +
  facet_wrap(~ agent,nrow=5, scales="free")+
  xlab("FW log10 load")+
  ylab("SR residual")+
  theme_bw()
dev.off()

## Load SW among years and stocks - box plots
ggplot(inf_agt_resid_data) +
#  geom_violin(aes(y=log10(mean_load), x=Stock_Analysis),draw_quantiles=c(0.25,0.5,0.75)) +
  geom_boxplot(aes(y=log10(mean_load), x=Stock_Analysis, col=Stock_Analysis)) +
  coord_flip() +
  facet_wrap(~ agent,nrow=5, scales="free")+
  xlab("Stocks")+
  ylab("Load")+
  theme_bw()





#######################################################################################
# STAN Approach for Multi-level Modeling

# SW PREVALENCE - INDEPENDENT MODELS by AGENT
## Stock-specific metric
### Create files for each agent
for(i in unique(inf_agt_resid_data$agent)) {
  nam <- paste("df", i, sep = ".")
  assign(nam, inf_agt_resid_data[inf_agt_resid_data$agent==i,])
}

### Loop for STAN independent models by agent
for(i in agents){
  data <- subset(inf_agt_resid_data, agent==i)
  nam <- paste("mod", i, sep = ".")
  assign(nam, stan_lmer(resid_value ~ 0 +  prev_std + (prev_std|Stock) +(1|Year), 
                        data = data,
                        adapt_delta=0.99,
                        REML = F))
}

# uninformed priors - for example:
prior_summary(mod.pa_ther)
plot(mod.pa_ther, "ess")
plot(mod.pa_ther, "trace")
summary(mod.c_b_cys)

## Shiny Stan launch
#my_sso.arena2 <- launch_shinystan(mod.arena2)
#my_sso.ic_mul <- launch_shinystan(mod.ic_mul)

## Derive coefficient estimates and save in .csv file
coefs_stan <- matrix(NA,
                     nrow = length(agents),
                     ncol = 5,
                     dimnames = list(agents,c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      pars = c("prev_std"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan[i,] <- ind_coef[1,c(4:8)]
}
write.csv(coefs_stan, file="data/prev_coefs_stan_stspec_210519.csv")

# Load estimates from file (if not running full model) and assign rownames
coefs_stan <- read.csv("data/prev_coefs_stan_stspec_210519.csv")
rownames(coefs_stan) <- coefs_stan[,1]
coefs_stan <- coefs_stan[,-1]  


# Plot effect size per agent
coefs_order <- coefs_stan[order(-coefs_stan[,3]),]

jpeg(filename='figs/Fig_prev_coefs_stan_stspec.jpg', 
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
       xlim = c(-2,2),
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
text(rep(-2,length(agents)), 
     seq(1,length(agents)), 
     labels = rownames(coefs_order), 
     pos = 4,
     font = 2,
     cex=0.95)
axis(1, at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5))
abline(v = 0, lty = 2)
box(col="grey")	
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Prevalence",3,line=0.25)
dev.off()

## Derive posterior estimates by stock
### Intercepts
coefs_stan_stk_int <- matrix(NA,
                             nrow = length(stocks),
                             ncol = 6,
                             dimnames = list(stocks,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\(\\Intercept) Stock\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_int", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_int <- cbind(ind_coef[c(1:17),c(4:8)], paste(i))))
}

### Slopes
coefs_stan_stk_slp <- matrix(NA,
                             nrow = length(stocks),
                             ncol = 6,
                             dimnames = list(stocks,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\prev_std Stock\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_slp", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_slp <- cbind(ind_coef[c(1:17),c(4:8)], paste(i))))
}

### Year intercepts
coefs_stan_stk_year <- matrix(NA,
                              nrow = length(years),
                              ncol = 6,
                              dimnames = list(years,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = "Year",
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_year", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_year <- cbind(ind_coef[c(1:8),c(4:8)], paste(i))))
}

### Sigmas
coefs_stan_stk_sig <- matrix(NA,
                             nrow = length(stocks),
                             ncol = 6,
                             dimnames = list(stocks,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("sigma","Sigma"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_sig", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_sig <- cbind(ind_coef[c(1:5),c(4:8)], paste(i))))
}

### Rhat and Neff
sims <-as.matrix(mod.arena2) #extract parameter names from any agent model 
dim(sims)
para_name <- c(colnames(sims), "mean_PPD", "log-posterior")
para_name

#### Rhat
coefs_stan_stk_rhat <- matrix(NA,
                              nrow = length(para_name),
                              ncol = 2,
                              dimnames = list(para_name,c("Rhat","agent")))
for(i in agents){
  model<-get(paste("mod",i, sep="."))
  ind_coef <- as.matrix(summary(model, 
                                digits = 3) [,"Rhat"])
  nam <- paste("coefs_stan_stk_rhat", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_rhat <- cbind(ind_coef[c(1:50),1], paste(i), paste("Rhat"))))
}

#### Neff
coefs_stan_stk_neff <- matrix(NA,
                              nrow = length(para_name),
                              ncol = 2,
                              dimnames = list(para_name,c("Neff","agent")))
for(i in agents){
  model<-get(paste("mod",i, sep="."))
  ind_coef <- as.matrix(summary(model) [,"n_eff"])
  nam <- paste("coefs_stan_stk_neff", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_neff <- cbind(ind_coef[c(1:50),1], paste(i), paste("Neff"))))
}

#### Rbind all posteriors and save as .csv file
post_int_slpintsig <- rbind(coefs_stan_stk_int.arena2,
                            coefs_stan_stk_int.c_b_cys,
                            coefs_stan_stk_int.ce_sha,
                            coefs_stan_stk_int.de_sal,
                            coefs_stan_stk_int.fa_mar,
                            coefs_stan_stk_int.fl_psy,
                            coefs_stan_stk_int.ic_hof,
                            coefs_stan_stk_int.ic_mul,
                            coefs_stan_stk_int.ku_thy,
                            coefs_stan_stk_int.lo_sal,
                            coefs_stan_stk_int.my_arc,
                            coefs_stan_stk_int.pa_kab,
                            coefs_stan_stk_int.pa_min,
                            coefs_stan_stk_int.pa_pse,
                            coefs_stan_stk_int.pa_ther,
                            coefs_stan_stk_int.prv,
                            coefs_stan_stk_int.pspv,
                            coefs_stan_stk_int.rlo,
                            coefs_stan_stk_int.sch,
                            coefs_stan_stk_int.smallUK,
                            coefs_stan_stk_int.sp_des,
                            coefs_stan_stk_int.te_bry,
                            coefs_stan_stk_int.te_mar,
                            coefs_stan_stk_int.ven,
                            coefs_stan_stk_slp.arena2,
                            coefs_stan_stk_slp.c_b_cys,
                            coefs_stan_stk_slp.ce_sha,
                            coefs_stan_stk_slp.de_sal,
                            coefs_stan_stk_slp.fa_mar,
                            coefs_stan_stk_slp.fl_psy,
                            coefs_stan_stk_slp.ic_hof,
                            coefs_stan_stk_slp.ic_mul,
                            coefs_stan_stk_slp.ku_thy,
                            coefs_stan_stk_slp.lo_sal,
                            coefs_stan_stk_slp.my_arc,
                            coefs_stan_stk_slp.pa_kab,
                            coefs_stan_stk_slp.pa_min,
                            coefs_stan_stk_slp.pa_pse,
                            coefs_stan_stk_slp.pa_ther,
                            coefs_stan_stk_slp.prv,
                            coefs_stan_stk_slp.pspv,
                            coefs_stan_stk_slp.rlo,
                            coefs_stan_stk_slp.sch,
                            coefs_stan_stk_slp.smallUK,
                            coefs_stan_stk_slp.sp_des,
                            coefs_stan_stk_slp.te_bry,
                            coefs_stan_stk_slp.te_mar,
                            coefs_stan_stk_slp.ven,
                            coefs_stan_stk_sig.arena2,
                            coefs_stan_stk_sig.c_b_cys,
                            coefs_stan_stk_sig.ce_sha,
                            coefs_stan_stk_sig.de_sal,
                            coefs_stan_stk_sig.fa_mar,
                            coefs_stan_stk_sig.fl_psy,
                            coefs_stan_stk_sig.ic_hof,
                            coefs_stan_stk_sig.ic_mul,
                            coefs_stan_stk_sig.ku_thy,
                            coefs_stan_stk_sig.lo_sal,
                            coefs_stan_stk_sig.my_arc,
                            coefs_stan_stk_sig.pa_kab,
                            coefs_stan_stk_sig.pa_min,
                            coefs_stan_stk_sig.pa_pse,
                            coefs_stan_stk_sig.pa_ther,
                            coefs_stan_stk_sig.prv,
                            coefs_stan_stk_sig.pspv,
                            coefs_stan_stk_sig.rlo,
                            coefs_stan_stk_sig.sch,
                            coefs_stan_stk_sig.smallUK,
                            coefs_stan_stk_sig.sp_des,
                            coefs_stan_stk_sig.te_bry,
                            coefs_stan_stk_sig.te_mar,
                            coefs_stan_stk_sig.ven,
                            coefs_stan_stk_year.arena2,
                            coefs_stan_stk_year.c_b_cys,
                            coefs_stan_stk_year.ce_sha,
                            coefs_stan_stk_year.de_sal,
                            coefs_stan_stk_year.fa_mar,
                            coefs_stan_stk_year.fl_psy,
                            coefs_stan_stk_year.ic_hof,
                            coefs_stan_stk_year.ic_mul,
                            coefs_stan_stk_year.ku_thy,
                            coefs_stan_stk_year.lo_sal,
                            coefs_stan_stk_year.my_arc,
                            coefs_stan_stk_year.pa_kab,
                            coefs_stan_stk_year.pa_min,
                            coefs_stan_stk_year.pa_pse,
                            coefs_stan_stk_year.pa_ther,
                            coefs_stan_stk_year.prv,
                            coefs_stan_stk_year.pspv,
                            coefs_stan_stk_year.rlo,
                            coefs_stan_stk_year.sch,
                            coefs_stan_stk_year.smallUK,
                            coefs_stan_stk_year.sp_des,
                            coefs_stan_stk_year.te_bry,
                            coefs_stan_stk_year.te_mar,
                            coefs_stan_stk_year.ven)

write.csv(post_int_slpintsig, file="data/Posterior distributions_Int Slp Sig_stspec_prev.csv")

#### Rbind all convergence parameters and save as .csv file
post_rhatneff_prev <- rbind(coefs_stan_stk_rhat.arena2,
                            coefs_stan_stk_rhat.c_b_cys,
                            coefs_stan_stk_rhat.ce_sha,
                            coefs_stan_stk_rhat.de_sal,
                            coefs_stan_stk_rhat.fa_mar,
                            coefs_stan_stk_rhat.fl_psy,
                            coefs_stan_stk_rhat.ic_hof,
                            coefs_stan_stk_rhat.ic_mul,
                            coefs_stan_stk_rhat.ku_thy,
                            coefs_stan_stk_rhat.lo_sal,
                            coefs_stan_stk_rhat.my_arc,
                            coefs_stan_stk_rhat.pa_kab,
                            coefs_stan_stk_rhat.pa_min,
                            coefs_stan_stk_rhat.pa_pse,
                            coefs_stan_stk_rhat.pa_ther,
                            coefs_stan_stk_rhat.prv,
                            coefs_stan_stk_rhat.pspv,
                            coefs_stan_stk_rhat.rlo,
                            coefs_stan_stk_rhat.sch,
                            coefs_stan_stk_rhat.smallUK,
                            coefs_stan_stk_rhat.sp_des,
                            coefs_stan_stk_rhat.te_bry,
                            coefs_stan_stk_rhat.te_mar,
                            coefs_stan_stk_rhat.ven,
                            coefs_stan_stk_neff.arena2,
                            coefs_stan_stk_neff.c_b_cys,
                            coefs_stan_stk_neff.ce_sha,
                            coefs_stan_stk_neff.de_sal,
                            coefs_stan_stk_neff.fa_mar,
                            coefs_stan_stk_neff.fl_psy,
                            coefs_stan_stk_neff.ic_hof,
                            coefs_stan_stk_neff.ic_mul,
                            coefs_stan_stk_neff.ku_thy,
                            coefs_stan_stk_neff.lo_sal,
                            coefs_stan_stk_neff.my_arc,
                            coefs_stan_stk_neff.pa_kab,
                            coefs_stan_stk_neff.pa_min,
                            coefs_stan_stk_neff.pa_pse,
                            coefs_stan_stk_neff.pa_ther,
                            coefs_stan_stk_neff.prv,
                            coefs_stan_stk_neff.pspv,
                            coefs_stan_stk_neff.rlo,
                            coefs_stan_stk_neff.sch,
                            coefs_stan_stk_neff.smallUK,
                            coefs_stan_stk_neff.sp_des,
                            coefs_stan_stk_neff.te_bry,
                            coefs_stan_stk_neff.te_mar,
                            coefs_stan_stk_neff.ven)

write.csv(post_rhatneff_prev, file="data/Posterior distributions_Rhat and Neff_stspec_prev.csv")

### Plot posteriors per agent model from files
post_all <- post_int_slpintsig
post_all <- read.csv("data/Posterior distributions_Int Slp Sig_stspec_prev.csv")
post_agents <- read.csv("data/prev_coefs_stan_stspec_210519.csv")
post_agents <- post_agents[order(match(post_agents[,1],propzero[,1])),]

## Plot Posteriors for all agents
jpeg(filename='figs/Fig_prev_coefs_stan_stspec_blue_210519.jpg', 
      width=480, height=500, quality=75)
ggplot(post_agents) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(X, -mid), ymax = X75, ymin = X25), size=1.5, col="blue") +
  geom_linerange(aes(x = X, ymax = upper, ymin = lower), col="blue") +
  geom_point(aes(x = X, y = mid), size = 3, col="blue")+
  labs(x ="Agents", y = "Effect size")+
  ggtitle("Sockeye Salmon slope estimates: \n95% & 50% credible intervals")+
  scale_x_discrete(labels=c("ic_mul" = "I. multifiliis", 
                            "te_mar" = "T. maritinum",
                            "pa_ther" = "P. theridion",
                            "fl_psy" = "F. psychrophilum",
                            "sch" = "Ca. S. salmonis",
                            "te_bry" = "T. bryosalmonae",
                            "pa_kab" = "P. kabatai",
                            "c_b_cys" = "Ca. B. cysticola",
                            "pa_min" = "P. minibicornis",
                            "arena2" = "Arenavirus",
                            "fa_mar" = "F. margolisi",
                            "my_arc" = "M. arcticus",
                            "ven" = "ENV",
                            "ic_hof" = "I. hoferi",
                            "lo_sal" = "L. salmonae",
                            "rlo" = "RLO",
                            "sp_des" = "S. destruens",
                            "ku_thy" = "K. thyrsites",
                            "prv" = "PRV",
                            "pspv" = "PSPV",
                            "ce_sha" = "C. shasta",
                            "pa_pse" = "P. pseudobranchicola",
                            "de_sal" = "D. salmonis"))+
  theme(axis.text.y = element_text(face = "italic"), plot.title = element_text(hjust = 0.5))+
  coord_flip()
dev.off()

## Plots per agent
#### Extract output from agent model
post_ic_mul <- post_all[post_all$X.1=="ic_mul",]

## Extract Posterior slopes by Stock - ic_mul
post_ic_mul_stockslp <- post_ic_mul[grep("prev_std Stock", post_ic_mul$X) ,]

#### Plot
jpeg(filename='figs/Fig_ic_mul_prev_stspec_bystock.jpg', 
     width=480, height=500, quality=75)
ggplot(post_ic_mul_stockslp) +
  geom_hline(yintercept = 0, linetype = "dashed", col="blue")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="black") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="black") +
  geom_point(aes(x = X, y = X50.), size = 3) +
  ylim(-0.8,0.8)+
  coord_flip()
dev.off()

## Extract Posterior intercepts for Stocks - ic_mul example
post_ic_mul_stockint <- post_ic_mul[c(1:18) ,]
## Plot
ggplot(post_ic_mul_stockint) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()

## Extract sigma values
post_ic_mul_stocksig <- post_ic_mul[grep("igma", post_ic_mul$X) ,]
## Plot
ggplot(post_ic_mul_stocksig) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = X, ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()

## Extract Posterior intercepts for years
post_ic_mul_year <- post_ic_mul[grep("b\\[\\(\\Intercept) Year", post_ic_mul$X) ,]
## Plot
ggplot(post_ic_mul_year) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()



## Calculate stock-specific Posterior estimates by adding 
## stock-specific draws to "global" (averaged across stocks) and then averaging

## Define generic parameter names
model<- mod.arena2 #dummy model
temp2 <- data.frame(model)
temp3 <- temp2[,grepl("prev",names(temp2))] #include only columns with "prev" in name
temp4 <- temp3[,-grep("igma",colnames(temp3))] #remove columns with "igma" in name
temp5 <- temp4[2001:4000,] #remove warm up iterations
temp6 <- data.frame(temp5[,1] + temp5[,2:18]) #add the stock-specific draws to global column
temp7 <- cbind(temp5[,1],temp6) #bind with global column
para_name2 <- colnames(temp5) #create an object with column names

## Create a matrix and loop
stk.spec.slope <- matrix(NA,
                         nrow = 18,
                         ncol = 6,
                         dimnames = list(para_name2,c("2.5","25","50","75","97.5","agent")))

for(i in agents){
  model<-get(paste("mod",i, sep="."))
  temp2 <- data.frame(model)
  temp3 <- temp2[,grepl("prev",names(temp2))] #include only columns with "prev" in name
  temp4 <- temp3[,-grep("igma",colnames(temp3))] #remove columns with "igma" in name
  temp5 <- temp4[2001:4000,] #remove warm up iterations
  temp6 <- data.frame(temp5[,1] + temp5[,2:18]) #add the stock-specific draws to global column
  temp7 <- cbind(temp5[,1],temp6) #bind with global column
  para_name2 <- colnames(temp5) #create an object with column names
  colnames(temp7) <- colnames(temp5) #assign names to columns
  nam <- paste("stk.spec.slope", i, sep = ".")
  temp8 <- as.matrix(apply(temp7, 2, quantile, probs = c(0.025,0.25,0.50,0.75,0.975)))
  temp9 <- t(temp8)
  temp10<- as.matrix(cbind(temp9, paste(i)))
  colnames(temp10) <- c("2.5","25","50","75","97.5","agent") #assign names to columns
  as.matrix(assign(nam, stk.spec.slope[i] <- temp10))
}

#### Rbind all posteriors and save as .csv file
stk.spec.slope.all <- rbind(stk.spec.slope.arena2,
                            stk.spec.slope.c_b_cys,
                            stk.spec.slope.ce_sha,
                            stk.spec.slope.de_sal,
                            stk.spec.slope.fa_mar,
                            stk.spec.slope.fl_psy,
                            stk.spec.slope.ic_hof,
                            stk.spec.slope.ic_mul,
                            stk.spec.slope.ku_thy,
                            stk.spec.slope.lo_sal,
                            stk.spec.slope.my_arc,
                            stk.spec.slope.pa_kab,
                            stk.spec.slope.pa_min,
                            stk.spec.slope.pa_pse,
                            stk.spec.slope.pa_ther,
                            stk.spec.slope.prv,
                            stk.spec.slope.pspv,
                            stk.spec.slope.rlo,
                            stk.spec.slope.sch,
                            stk.spec.slope.sp_des,
                            stk.spec.slope.te_bry,
                            stk.spec.slope.te_mar,
                            stk.spec.slope.ven)
write.csv(stk.spec.slope.all, file="data/Stock specific slopes_prev_stspec_210519.csv")

##### READ IN DATA FROM FILE
stspslp <- read.csv("data/Stock specific slopes_prev_stspec_210519.csv")
stspslp$stock <- substr(stspslp$X, 18, 28)
stspslp$stock <- substr(stspslp$stock, 1, nchar(stspslp$stock)-1)
stspslp$stock <- sub("^$", "Overall", stspslp$stock)

## Extract stock-specific slopes - ic_mul model
stk.spec.ic_mul <-stspslp[stspslp$agent=="ic_mul",]
## Plot
jpeg(filename='figs/Fig_stockspslope_ic_mul_stspec_210519.jpg', 
     width=480, height=500, quality=75)
ggplot(stk.spec.ic_mul) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.ic_mul[stk.spec.ic_mul$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.ic_mul[stk.spec.ic_mul$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.ic_mul[stk.spec.ic_mul$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size", title="I. multifiliis") +
  coord_flip()
dev.off()

## Plot stock-specific slopes - pa_ther
stk.spec.pa_ther <-stspslp[stspslp$agent=="pa_ther",]
## Plot
jpeg(filename='figs/Fig_SSHI ONNE_stock sp slope_pa_ther_210519.jpg', 
     width=480, height=500, quality=75)
ggplot(stk.spec.pa_ther) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.pa_ther[stk.spec.pa_ther$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.pa_ther[stk.spec.pa_ther$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.pa_ther[stk.spec.pa_ther$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size") +
  coord_flip()
dev.off()


## Plot stock-specific slopes - c_b_cys
stk.spec.c_b_cys <-stspslp[stspslp$agent=="c_b_cys",]
## Plot
jpeg(filename='figs/Fig_stockspslope_c_b_cys_stspec_210519.jpg', 
     width=480, height=500, quality=75)
ggplot(stk.spec.c_b_cys) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.c_b_cys[stk.spec.c_b_cys$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.c_b_cys[stk.spec.c_b_cys$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.c_b_cys[stk.spec.c_b_cys$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size", title="P. theridion") +
  coord_flip()
dev.off()


## Derive proportions of posterior draws <0 per model - prevalence
param<-colnames(data.frame(mod.arena2)) #create object of parameters in model
param.prop0 <- matrix(NA,
                      ncol = 1,
                      nrow = length(agents),
                      dimnames = list(agents,"prev_prop_neg"))

for (i in agents){
  model<-as.matrix(get(paste("mod.",i, sep="")))
  model2<-model[2001:4000,]
  para_name <- colnames(model2)
  temp <- as.matrix((colSums(model2 < 0))/2000)
  param.prop0[i,] <- temp[1,]
}
write.csv(param.prop0, file="data/Percent post draws >0_prev_stspec_210519.csv")




# plot
propzero <- read.csv("data/Percent post draws >0_prev_stspec_210519.csv")
propzero$plot.agent <- propzero$X
propzero$plot.agent <- recode(propzero$plot.agent, "ic_mul" = "I. multifiliis", 
                                        "te_mar" = "T. maritinum",
                                        "pa_ther" = "P. theridion",
                                        "fl_psy" = "F. psychrophilum",
                                        "sch" = "Ca. S. salmonis",
                                        "te_bry" = "T. bryosalmonae",
                                        "pa_kab" = "P. kabatai",
                                        "c_b_cys" = "Ca. B. cysticola",
                                        "pa_min" = "P. minibicornis",
                                        "arena2" = "SPAV-2",
                                        "fa_mar" = "F. margolisi",
                                        "my_arc" = "M. arcticus",
                                        "ven" = "ENV",
                                        "ic_hof" = "I. hoferi",
                                        "lo_sal" = "L. salmonae",
                                        "rlo" = "RLO",
                                        "sp_des" = "S. destruens",
                                        "ku_thy" = "K. thyrsites",
                                        "prv" = "PRV",
                                        "pspv" = "PSPV",
                                        "ce_sha" = "C. shasta",
                                        "pa_pse" = "P. pseudobranchicola",
                                        "de_sal" = "D. salmonis")

jpeg(filename='figs/Fig_prop<0_pa_ther_stspec.jpg', 
     width=480, height=500, quality=75)
ggplot(propzero) +
  geom_bar(stat="identity", aes(reorder(plot.agent, prev_prop_neg), prev_prop_neg), col="blue", fill="blue", width=0.5, alpha=0.5) +
  labs(x="Agent", y="Proportion of total") +
  ggtitle("Sockeye salmon: \nProportion of correlation coefficients < 0")+
  theme(axis.text.y = element_text(face = "italic"))+
  coord_flip()
dev.off()

### Violin plot?
ggplot(propzero, aes(plot.agent, prev_prop_neg)) + 
  geom_violin()+
  coord_flip()





##############################################
## Prevalence model with SST as fixed effect and stock-specific metric
##############################################
#### SST data
raw.clim <- read.csv("data/sst_yr_1_stock_anomalies.csv")
early.sst <- clim.wgt.avg(brood.table = brood_table,
                          env.data = raw.clim,
                          env.covar = "sst_anomaly",
                          type = "first_year",
                          out.covar = "early_sst") 
head(early.sst)
stock.ids2 <- brood_table[,1:2]
stock.ids <- stock.ids2[!duplicated(stock.ids2),]
early.sst <- merge(early.sst, stock.ids, by="Stock.ID")
names(early.sst) <- c("Stock.ID", "brood_year", "sst_anom", "Stock")

## Merge with pathogen data
dim(inf_agt_resid_data)
inf_agt_resid_data_sst <- merge(inf_agt_resid_data, early.sst, 
                               by = c("Stock", "brood_year"), all.x=TRUE)
head(inf_agt_resid_data_sst)

### Create files for each agent
for(i in unique(inf_agt_resid_data_sst$agent)) {
  nam <- paste("df.sst", i, sep = ".")
  assign(nam, inf_agt_resid_data_sst[inf_agt_resid_data_sst$agent==i,])
}

### Loop for STAN independent models by agent THIS OVER WRITES THE MODELS ABOVE!
for(i in agents){
  data <- subset(inf_agt_resid_data_sst, agent==i)
  nam <- paste("mod", i, sep = ".")
  assign(nam, stan_lmer(resid_value ~ 0 +  prev_std + sst_anom + (prev_std|Stock) +(1|Year), 
                        data = data,
                        adapt_delta=0.99,
                        REML = F))
}

## Derive coefficient estimates and save in .csv file
coefs_stan <- matrix(NA,
                     nrow = length(agents),
                     ncol = 5,
                     dimnames = list(agents,c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      pars = c("prev_std"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan[i,] <- ind_coef[1,c(4:8)]
}
write.csv(coefs_stan, file="data/prev_coefs_stan_stspec.sst.csv")

## Derive SST coefficient estimates and save in .csv file
coefs_stan_sst <- matrix(NA,
                     nrow = length(agents),
                     ncol = 5,
                     dimnames = list(agents,c("lower","25","mid","75","upper")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      pars = c("sst_anom"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan_sst[i,] <- ind_coef[1,c(4:8)]
}
write.csv(coefs_stan_sst, file="data/prev_coefs_stan_stspec.sst_SST coefs.csv")

# Load estimates from file (if not running full model) and assign rownames
coefs_stan <- read.csv("data/prev_coefs_stan_stspec.sst.csv")
rownames(coefs_stan) <- coefs_stan[,1]
coefs_stan <- coefs_stan[,-1]  

## Derive posterior estimates by stock
### Intercepts
coefs_stan_stk_int.sst <- matrix(NA,
                             nrow = length(stocks),
                             ncol = 6,
                             dimnames = list(stocks,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod",i, sep="."))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\(\\Intercept) Stock\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_int.sst", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_int.sst <- cbind(ind_coef[c(1:17),c(4:8)], paste(i))))
}

### Slopes
coefs_stan_stk_slp.sst <- matrix(NA,
                             nrow = length(stocks),
                             ncol = 6,
                             dimnames = list(stocks,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\prev_std Stock\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_slp.sst", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_slp.sst <- cbind(ind_coef[c(1:17),c(4:8)], paste(i))))
}

### Year intercepts
coefs_stan_stk_year.sst <- matrix(NA,
                              nrow = length(years),
                              ncol = 6,
                              dimnames = list(years,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = "Year",
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_year.sst", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_year.sst <- cbind(ind_coef[c(1:8),c(4:8)], paste(i))))
}

### Sigmas
coefs_stan_stk_sig.sst <- matrix(NA,
                             nrow = length(stocks),
                             ncol = 6,
                             dimnames = list(stocks,c("lower","25","mid","75","upper","agent")))
for(i in agents){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("sigma","Sigma"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_sig.sst", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_sig.sst <- cbind(ind_coef[c(1:5),c(4:8)], paste(i))))
}

### Rhat and Neff
sims <-as.matrix(mod.arena2) #extract parameter names from any agent model 
dim(sims)
para_name <- c(colnames(sims), "mean_PPD", "log-posterior")
para_name

#### Rhat
coefs_stan_stk_rhat.sst <- matrix(NA,
                              nrow = length(para_name),
                              ncol = 2,
                              dimnames = list(para_name,c("Rhat","agent")))
for(i in agents){
  model<-get(paste("mod",i, sep="."))
  ind_coef <- as.matrix(summary(model, 
                                digits = 3) [,"Rhat"])
  nam <- paste("coefs_stan_stk_rhat.sst", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_rhat.sst <- cbind(ind_coef[c(1:50),1], paste(i), paste("Rhat"))))
}

#### Neff
coefs_stan_stk_neff.sst <- matrix(NA,
                              nrow = length(para_name),
                              ncol = 2,
                              dimnames = list(para_name,c("Neff","agent")))
for(i in agents){
  model<-get(paste("mod",i, sep="."))
  ind_coef <- as.matrix(summary(model) [,"n_eff"])
  nam <- paste("coefs_stan_stk_neff.sst", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_neff.sst <- cbind(ind_coef[c(1:50),1], paste(i), paste("Neff"))))
}

#### Rbind all posteriors and save as .csv file
post_int_slpintsig.sst <- rbind(coefs_stan_stk_int.sst.arena2,
                            coefs_stan_stk_int.sst.c_b_cys,
                            coefs_stan_stk_int.sst.ce_sha,
                            coefs_stan_stk_int.sst.de_sal,
                            coefs_stan_stk_int.sst.fa_mar,
                            coefs_stan_stk_int.sst.fl_psy,
                            coefs_stan_stk_int.sst.ic_hof,
                            coefs_stan_stk_int.sst.ic_mul,
                            coefs_stan_stk_int.sst.ku_thy,
                            coefs_stan_stk_int.sst.lo_sal,
                            coefs_stan_stk_int.sst.my_arc,
                            coefs_stan_stk_int.sst.pa_kab,
                            coefs_stan_stk_int.sst.pa_min,
                            coefs_stan_stk_int.sst.pa_pse,
                            coefs_stan_stk_int.sst.pa_ther,
                            coefs_stan_stk_int.sst.pspv,
                            coefs_stan_stk_int.sst.rlo,
                            coefs_stan_stk_int.sst.sch,
                            coefs_stan_stk_int.sst.sp_des,
                            coefs_stan_stk_int.sst.te_bry,
                            coefs_stan_stk_int.sst.te_mar,
                            coefs_stan_stk_int.sst.ven,
                            coefs_stan_stk_slp.sst.arena2,
                            coefs_stan_stk_slp.sst.c_b_cys,
                            coefs_stan_stk_slp.sst.ce_sha,
                            coefs_stan_stk_slp.sst.de_sal,
                            coefs_stan_stk_slp.sst.fa_mar,
                            coefs_stan_stk_slp.sst.fl_psy,
                            coefs_stan_stk_slp.sst.ic_hof,
                            coefs_stan_stk_slp.sst.ic_mul,
                            coefs_stan_stk_slp.sst.ku_thy,
                            coefs_stan_stk_slp.sst.lo_sal,
                            coefs_stan_stk_slp.sst.my_arc,
                            coefs_stan_stk_slp.sst.pa_kab,
                            coefs_stan_stk_slp.sst.pa_min,
                            coefs_stan_stk_slp.sst.pa_pse,
                            coefs_stan_stk_slp.sst.pa_ther,
                            coefs_stan_stk_slp.sst.pspv,
                            coefs_stan_stk_slp.sst.rlo,
                            coefs_stan_stk_slp.sst.sch,
                            coefs_stan_stk_slp.sst.sp_des,
                            coefs_stan_stk_slp.sst.te_bry,
                            coefs_stan_stk_slp.sst.te_mar,
                            coefs_stan_stk_slp.sst.ven,
                            coefs_stan_stk_sig.sst.arena2,
                            coefs_stan_stk_sig.sst.c_b_cys,
                            coefs_stan_stk_sig.sst.ce_sha,
                            coefs_stan_stk_sig.sst.de_sal,
                            coefs_stan_stk_sig.sst.fa_mar,
                            coefs_stan_stk_sig.sst.fl_psy,
                            coefs_stan_stk_sig.sst.ic_hof,
                            coefs_stan_stk_sig.sst.ic_mul,
                            coefs_stan_stk_sig.sst.ku_thy,
                            coefs_stan_stk_sig.sst.lo_sal,
                            coefs_stan_stk_sig.sst.my_arc,
                            coefs_stan_stk_sig.sst.pa_kab,
                            coefs_stan_stk_sig.sst.pa_min,
                            coefs_stan_stk_sig.sst.pa_pse,
                            coefs_stan_stk_sig.sst.pa_ther,
                            coefs_stan_stk_sig.sst.pspv,
                            coefs_stan_stk_sig.sst.rlo,
                            coefs_stan_stk_sig.sst.sch,
                            coefs_stan_stk_sig.sst.sp_des,
                            coefs_stan_stk_sig.sst.te_bry,
                            coefs_stan_stk_sig.sst.te_mar,
                            coefs_stan_stk_sig.sst.ven,
                            coefs_stan_stk_year.sst.arena2,
                            coefs_stan_stk_year.sst.c_b_cys,
                            coefs_stan_stk_year.sst.ce_sha,
                            coefs_stan_stk_year.sst.de_sal,
                            coefs_stan_stk_year.sst.fa_mar,
                            coefs_stan_stk_year.sst.fl_psy,
                            coefs_stan_stk_year.sst.ic_hof,
                            coefs_stan_stk_year.sst.ic_mul,
                            coefs_stan_stk_year.sst.ku_thy,
                            coefs_stan_stk_year.sst.lo_sal,
                            coefs_stan_stk_year.sst.my_arc,
                            coefs_stan_stk_year.sst.pa_kab,
                            coefs_stan_stk_year.sst.pa_min,
                            coefs_stan_stk_year.sst.pa_pse,
                            coefs_stan_stk_year.sst.pa_ther,
                            coefs_stan_stk_year.sst.pspv,
                            coefs_stan_stk_year.sst.rlo,
                            coefs_stan_stk_year.sst.sch,
                            coefs_stan_stk_year.sst.sp_des,
                            coefs_stan_stk_year.sst.te_bry,
                            coefs_stan_stk_year.sst.te_mar,
                            coefs_stan_stk_year.sst.ven)

write.csv(post_int_slpintsig.sst, file="data/Posterior distributions_Int Slp Sig_stspec_prev_SST.csv")

#### Rbind all convergence parameters and save as .csv file
post_rhatneff_prev.sst <- rbind(coefs_stan_stk_rhat.sst.arena2,
                            coefs_stan_stk_rhat.sst.c_b_cys,
                            coefs_stan_stk_rhat.sst.ce_sha,
                            coefs_stan_stk_rhat.sst.de_sal,
                            coefs_stan_stk_rhat.sst.fa_mar,
                            coefs_stan_stk_rhat.sst.fl_psy,
                            coefs_stan_stk_rhat.sst.ic_hof,
                            coefs_stan_stk_rhat.sst.ic_mul,
                            coefs_stan_stk_rhat.sst.ku_thy,
                            coefs_stan_stk_rhat.sst.lo_sal,
                            coefs_stan_stk_rhat.sst.my_arc,
                            coefs_stan_stk_rhat.sst.pa_kab,
                            coefs_stan_stk_rhat.sst.pa_min,
                            coefs_stan_stk_rhat.sst.pa_pse,
                            coefs_stan_stk_rhat.sst.pa_ther,
                            coefs_stan_stk_rhat.sst.pspv,
                            coefs_stan_stk_rhat.sst.rlo,
                            coefs_stan_stk_rhat.sst.sch,
                            coefs_stan_stk_rhat.sst.sp_des,
                            coefs_stan_stk_rhat.sst.te_bry,
                            coefs_stan_stk_rhat.sst.te_mar,
                            coefs_stan_stk_rhat.sst.ven,
                            coefs_stan_stk_neff.sst.arena2,
                            coefs_stan_stk_neff.sst.c_b_cys,
                            coefs_stan_stk_neff.sst.ce_sha,
                            coefs_stan_stk_neff.sst.de_sal,
                            coefs_stan_stk_neff.sst.fa_mar,
                            coefs_stan_stk_neff.sst.fl_psy,
                            coefs_stan_stk_neff.sst.ic_hof,
                            coefs_stan_stk_neff.sst.ic_mul,
                            coefs_stan_stk_neff.sst.ku_thy,
                            coefs_stan_stk_neff.sst.lo_sal,
                            coefs_stan_stk_neff.sst.my_arc,
                            coefs_stan_stk_neff.sst.pa_kab,
                            coefs_stan_stk_neff.sst.pa_min,
                            coefs_stan_stk_neff.sst.pa_pse,
                            coefs_stan_stk_neff.sst.pa_ther,
                            coefs_stan_stk_neff.sst.pspv,
                            coefs_stan_stk_neff.sst.rlo,
                            coefs_stan_stk_neff.sst.sch,
                            coefs_stan_stk_neff.sst.sp_des,
                            coefs_stan_stk_neff.sst.te_bry,
                            coefs_stan_stk_neff.sst.te_mar,
                            coefs_stan_stk_neff.sst.ven)

write.csv(post_rhatneff_prev.sst, file="data/Posterior distributions_Rhat and Neff_stspec_prev_SST.csv")

## Derive proportions of posterior draws <0 per model - prevalence+SST
prop_less_than_0_prevsst <- matrix(NA,
                                ncol = 1,
                                nrow = length(agents),
                                dimnames = list(agents,"SSTagents_prop_neg"))
for (i in agents){
  model<-as.matrix(get(paste("mod.",i, sep="")))
  model2<-model[2001:4000,]
  para_name <- colnames(model2)
  temp <- as.matrix((colSums(model2 < 0))/2000)
  prop_less_than_0_prevsst[i,] <- temp[1,]
}
write.csv(prop_less_than_0_prevsst, file="data/Prop posterior neg_prevSST.csv")


### Plot posteriors per agent model from files
post_all.sst <- post_int_slpintsig.sst
post_all.sst <- read.csv("data/Posterior distributions_Int Slp Sig_stspec_prev.csv")
post_agents <- read.csv("data/prev_coefs_stan_stspec.sst.csv")

## REMOVE ARENA2 AND SMALLUK
post_all.sst <- post_all.sst[!(post_all.sst$X.1=="arena2" | post_all.sst$X.1=="smallUK"),]
post_agents <- post_agents[!(post_agents$X=="arena2" | post_agents$X=="smallUK"),]
  
## Plot Posteriors for all agents
jpeg(filename='figs/Fig_prev_coefs_stan_stspec_blue_sst.jpg', 
     width=480, height=500, quality=75)
ggplot(post_agents) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(X, -mid), ymax = X75, ymin = X25), size=1.5, col="darkblue") +
  geom_linerange(aes(x = X, ymax = upper, ymin = lower), col="darkblue") +
  geom_point(aes(x = X, y = mid), size = 3, col="darkblue")+
  labs(x ="Infectious agents", y = "Effect size", title="Prevalence beta w SST")+
  scale_x_discrete(labels=c("ic_mul" = "I. multifiliis", 
                            "te_mar" = "T. maritinum",
                            "pa_ther" = "P. theridion",
                            "fl_psy" = "F. psychrophilum",
                            "sch" = "Ca. S. salmonis",
                            "te_bry" = "T. bryosalmonae",
                            "pa_kab" = "P. kabatai",
                            "c_b_cys" = "Ca. B. cysticola",
                            "pa_min" = "P. minibicornis",
                            "arena2" = "Arenavirus",
                            "fa_mar" = "F. margolisi",
                            "my_arc" = "M. arcticus",
                            "ven" = "ENV",
                            "ic_hof" = "I. hoferi",
                            "lo_sal" = "L. salmonae",
                            "rlo" = "RLO",
                            "sp_des" = "S. destruens",
                            "ku_thy" = "K. thyrsites",
                            "prv" = "PRV",
                            "pspv" = "PSPV",
                            "ce_sha" = "C. shasta",
                            "pa_pse" = "P. pseudobranchicola",
                            "de_sal" = "D. salmonis"))+
  theme(axis.text.y = element_text(face = "italic"), plot.title = element_text(hjust = 0.5))+
  coord_flip()
dev.off()

## Read files PREV models w/ & w/o SST - remove ARENA2
post_agents <- read.csv("data/prev_coefs_stan_stspec.csv")
post_agents.noarena2 <- post_agents[post_agents$X!="arena2",]
post_agents.sst <- read.csv("data/prev_coefs_stan_stspec.sst.csv")
post_agents.sst.noarena2 <- post_agents.sst[post_agents.sst$X!="arena2",]
post_sst <- read.csv("data/prev_coefs_stan_stspec.sst_SST coefs.csv")
post_sst.noarena2 <- post_sst[post_sst$X!="arena2",]
post_agents.noarena2$sst_cofactor <- "No"
post_agents.sst.noarena2$sst_cofactor <- "Yes"
sst_impact <- rbind(post_agents.noarena2, post_agents.sst.noarena2)
post_agents.noarena2 <- post_agents.noarena2[order(-post_agents.noarena2$mid),] 
post_agents.sst.noarena2 <- post_agents.sst.noarena2[order(match(post_agents.sst.noarena2[,1],post_agents.noarena2[,1])),]
post_sst.noarena2 <- post_sst.noarena2[order(match(post_sst.noarena2[,1],post_agents.noarena2[,1])),]
post_agents <- post_agents.noarena2
post_agents.sst <- post_agents.sst.noarena2
post_sst <- post_sst.noarena2
post_sst$sst_cofactor <- "Yes"

#order prop0 files to match post_agents.noarena2
#merge with beta files
#table and plot

# Plot
agents.sst <- data.frame(unique(post_agents$X))
agents.sst$plot.agent <- agents.sst$unique.post_agents.X.
agents.sst$plot.agent <- recode(agents.sst$plot.agent, "ic_mul" = "I. multifiliis", 
                                        "te_mar" = "T. maritinum",
                                        "pa_ther" = "P. theridion",
                                        "fl_psy" = "F. psychrophilum",
                                        "sch" = "Ca. S. salmonis",
                                        "te_bry" = "T. bryosalmonae",
                                        "pa_kab" = "P. kabatai",
                                        "c_b_cys" = "Ca. B. cysticola",
                                        "pa_min" = "P. minibicornis",
                                        "arena2" = "SPAV-2",
                                        "fa_mar" = "F. margolisi",
                                        "my_arc" = "M. arcticus",
                                        "ven" = "ENV",
                                        "ic_hof" = "I. hoferi",
                                        "lo_sal" = "L. salmonae",
                                        "rlo" = "RLO",
                                        "sp_des" = "S. destruens",
                                        "ku_thy" = "K. thyrsites",
                                        "prv" = "PRV",
                                        "pspv" = "PSPV",
                                        "ce_sha" = "C. shasta",
                                        "pa_pse" = "P. pseudobranchicola",
                                        "de_sal" = "D. salmonis")
jpeg(filename='figs/Fig_Beta estimates with-without SST cofactor_stspec_fullname.jpg', 
     width=510, height=500, quality=300)
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))
#add coefficient estimates for the sst cofactor
plotCI(x = post_sst[,4], col="red",
       y = seq(1,nrow(agents.sst)),
       li = (post_sst[,2]),
       ui = (post_sst[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-2.5,1.5),
       pch = 16,
       scol = "red")
plotCI(x = post_sst[,4], col="red",
       y = seq(1,nrow(agents.sst)),
       li = (post_sst[,3]),
       ui = (post_sst[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "red")
#add coefficient estimates with sst cofactor
plotCI(x = post_agents.sst[,4], col="gray",
       y = seq(1,nrow(agents.sst)),
       li = (post_agents.sst[,2]),
       ui = (post_agents.sst[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-2.5,1.5),
       pch = 16,
       add = TRUE,
       scol = "gray")
plotCI(x = post_agents.sst[,4], col="gray",
       y = seq(1,nrow(agents.sst)),
       li = (post_agents.sst[,3]),
       ui = (post_agents.sst[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       lwd = 3,
       add = TRUE,
       scol = "gray")
#add coefficient estimates without SST cofactor
plotCI(x = post_agents[,4], col="black",
       y = seq(1,nrow(agents.sst))+.25,
       li = (post_agents[,2]),
       ui = (post_agents[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-2.5,1.5),
       pch = 16,
       add = TRUE,
       scol = "black")
plotCI(x = post_agents[,4], col="black",
       y = seq(1,nrow(agents.sst))+.25,
       li = (post_agents[,3]),
       ui = (post_agents[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "black")
text(rep(-2.5,nrow(agents.sst)), 
     seq(1,nrow(agents.sst)), 
     labels = (agents.sst[,2]), 
     pos = 4,
     font = 3,
     cex=0.95)
legend(x=.5, # x coordinate of the top left of the legend
  y=21,
  legend=c("Agent", "Agent + SST", "SST"), 
  pch=21,
  pt.bg=c("black","gray","red"))
axis(1, at = c(-2, -1.5,-1, -0.5, 0, 0.5, 1, 1.5))
abline(v = 0, lty = 2)
box(col="grey") 
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Prevalence",3,line=0.25)
dev.off()


## Merge into 1 df
post_sst$beta <- "SST"
post_agents.sst$beta <- "prev"
post_agents$beta <- "prev"
prev.fixed.beta <- rbind(post_agents,post_agents.sst,post_sst)








#################################################################################
#################################################################################
# SW LOAD - INDEPENDENT MODELS by AGENT

## Trim data set to agents with >3 detections across all years per stock
head(inf_agt_resid_data)
table.temp <- data.frame(inf_agt_resid_data %>%
  group_by(agent, Stock_Analysis) %>%
  summarise(yrcount = sum(!is.na(load_std))))
table.temp2 <- data.frame(table.temp %>%
                           group_by(agent) %>%
                           summarise(yrstcount = sum(yrcount>2)))
table.temp2[order(table.temp2$yrstcount),]
ggplot(table.temp, aes(agent, Stock_Analysis, fill=yrcount)) + 
  geom_tile()
#### Reduce data set for load STAN analysis down to AGENTS WITH AT LEAST 3 YEARS OF DETECTIONS IN AT LEAST 5 STOCKS
##### Keep: ic_hof,arena2,lo_sal,ce_sha,ven,pa_pse,sch,pa_kab,pspv,c_b_cys,my_arc,pa_min,pa_ther  
agents_keep <- c("ic_hof","lo_sal","ce_sha","ven","pa_pse","sch","pa_kab","pspv","c_b_cys","my_arc","pa_min","pa_ther")
#inf_agt_resid_data_reduced <- inf_agt_resid_data[inf_agt_resid_data$agent %in% agents_keep, ]


## Stock-specific SW metric 
### Loop for STAN independent models
for(i in agents_keep){
  data <- subset(inf_agt_resid_data, agents_keep==i)
  nam <- paste("mod.load", i, sep = ".")
  assign(nam, stan_lmer(resid_value ~ 0 + load_std + (load_std|Stock) +(1|Year), 
                        data = data,
                        adapt_delta=0.99,
                        REML = F))
}

## Loop to derive coefficient estimates
coefs_stan_l <- matrix(NA,
                     nrow = length(agents_keep),
                     ncol = 5,
                     dimnames = list(agents_keep,c("lower","25","mid","75","upper")))
for(i in agents_keep){
  model<-get(paste("mod.load",i, sep="."))
  ind_coef <- summary(model, 
                      pars = c("load_std"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan_l[i,] <- ind_coef[1,c(4:8)]
}
write.csv(coefs_stan_l, file="data/load_coefs_stan_stspec_reducedagents.csv")

# Load estimates from file and add rownames
coefs_stan <- read.csv("data/load_coefs_stan_stspec_reducedagents.csv")
rownames(coefs_stan) <- coefs_stan[,1]
coefs_stan <- coefs_stan[,-1] 

# Plot effect size of agents
coefs_order <- coefs_stan[order(-coefs_stan[,3]),]
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))
plotCI(x = coefs_order[,3],
       y = seq(1,length(agents_keep)),
       li = (coefs_order[,1]),
       ui = (coefs_order[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-1.75,2),
       pch = 16,
       scol = "grey")
plotCI(x = coefs_order[,3],
       y = seq(1,length(agents_keep)),
       li = (coefs_order[,2]),
       ui = (coefs_order[,4]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "grey")
text(rep(-1.75,length(agents_keep)), 
     seq(1,length(agents_keep)), 
     labels = rownames(coefs_order), 
     pos = 4,
     font = 2,
     cex=0.95)
axis(1, at = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))
abline(v = 0, lty = 2)
box(col="grey")	
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Intensity",3,line=0.25)


## Loop to derive posteriors by stock
### Intercepts
for(i in agents_keep){
  model<-get(paste("mod.load",i, sep="."))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\(\\Intercept) Stock\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  stocksl <- as.matrix(ind_coef[,1])
  stocksl1 <- rownames(stocksl)
  stocksl2 <- substr(stocksl1, 21, 32)
  stocksl3 <- substr(stocksl2, 1, nchar(stocksl2)-1)
  coefs_stan_stk_int_load <- matrix(NA,
                                    nrow = length(stocksl3),
                                    ncol = 6,
                                    dimnames = list(stocksl3,c("lower","25","mid","75","upper","agent")))
  stocks.mod <- length(stocksl3)
  nam <- paste("coefs_stan_stk_int_load", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_int_load <- cbind(ind_coef[1:stocks.mod,c(4:8)], paste(i))))
}


### Slopes
for(i in agents_keep){
  model<-get(paste("mod.load",i, sep="."))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\load_std Stock\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  stocksl <- as.matrix(ind_coef[,1])
  stocksl1 <- rownames(stocksl)
  stocksl2 <- substr(stocksl1, 21, 32)
  stocksl3 <- substr(stocksl2, 1, nchar(stocksl2)-1)
  coefs_stan_stk_slp_load <- matrix(NA,
                                    nrow = length(stocksl3),
                                    ncol = 6,
                                    dimnames = list(stocksl3,c("lower","25","mid","75","upper","agent")))
  stocks.mod <- length(stocksl3)
  nam <- paste("coefs_stan_stk_slp_load", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_slp_load <- cbind(ind_coef[1:stocks.mod,c(4:8)], paste(i))))
}

### Year intercepts
for(i in agents_keep){
  model<-get(paste("mod.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = "b\\[\\(\\Intercept) Year",
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  temp1 <- as.matrix(ind_coef[,1])
  temp2 <- rownames(temp1)
  temp3 <- substr(temp2, 20, 32)
  temp4 <- substr(temp3, 1, nchar(temp3)-1)
  coefs_stan_stk_year_load <- matrix(NA,
                                    nrow = length(temp4),
                                    ncol = 6,
                                    dimnames = list(temp4,c("lower","25","mid","75","upper","agent")))
  nam <- paste("coefs_stan_stk_year_load", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_year_load[i] <- cbind(ind_coef[c(1:length(temp4)),c(4:8)], paste(i))))
}

### Sigmas
for(i in agents_keep){
  model<-get(paste("mod.load",i, sep="."))
  ind_coef <- summary(model, 
                      regex_pars = c("sigma","Sigma"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  temp1 <- as.matrix(ind_coef[,1])
  temp2 <- rownames(temp1)
  coefs_stan_stk_sig_load <- matrix(NA,
                                    nrow = length(temp2),
                                    ncol = 6,
                                    dimnames = list(temp2,c("lower","25","mid","75","upper","agent")))
  nam <- paste("coefs_stan_stk_sig_load", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_sig_load[i] <- cbind(ind_coef[c(1:length(temp2)),c(4:8)], paste(i))))
}

### Rhat and Neff
#extract parameter names
sims<-as.matrix(mod.load.c_b_cys) #use any model to get parameter names, they should be the same
dim(sims)
para_name <- c(colnames(sims), "mean_PPD", "log-posterior")
para_name


for(i in agents_keep){
  model<-get(paste("mod.load",i, sep="."))
  ind_coef <- as.matrix(summary(model, 
                                digits = 3) [,"Rhat"])
  coefs_stan_stk_rhat_load <- matrix(NA,
                                     nrow = nrow(ind_coef),
                                     ncol = 2,
                                     dimnames = list(rownames(ind_coef),c("Rhat","agent")))
  nam <- paste("coefs_stan_stk_rhat_load", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_rhat_load[i] <- cbind(ind_coef[c(1:nrow(ind_coef)),1], paste(i), paste("Rhat"))))
}

for(i in agents_keep){
  model<-get(paste("mod.load",i, sep="."))
  ind_coef <- as.matrix(summary(model) [,"n_eff"])
  coefs_stan_stk_neff_load <- matrix(NA,
                                     nrow = nrow(ind_coef),
                                     ncol = 2,
                                     dimnames = list(rownames(ind_coef),c("Neff","agent")))
  nam <- paste("coefs_stan_stk_neff_load", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_neff_load[i] <- cbind(ind_coef[c(1:nrow(ind_coef)),1], paste(i), paste("Neff"))))
}

### rbind all files
#intercept
df <- matrix(ncol = 6, nrow = 0)
colnames(df) <- colnames(coefs_stan_stk_int_load.c_b_cys)
for(i in agents_keep){
  df1 <- get(paste("coefs_stan_stk_int_load",i, sep="."))
  df2 <- rbind(df1,df)
  df <- df2
}
#slopes
for(i in agents_keep){
  df1 <- get(paste("coefs_stan_stk_slp_load",i, sep="."))
  df2 <- rbind(df1,df)
  df <- df2
}
#sigmas
for(i in agents_keep){
  df1 <- get(paste("coefs_stan_stk_sig_load",i, sep="."))
  df2 <- rbind(df1,df)
  df <- df2
}
#year
for(i in agents_keep){
  df1 <- get(paste("coefs_stan_stk_year_load",i, sep="."))
  df2 <- rbind(df1,df)
  df <- df2
}
write.csv(df, file="data/PosteriorEst_IntSlpSigYr_load_stsp_reduced.csv")


#### Rbind all convergence parameters and write as csv
#rhat
df <- matrix(ncol = 3, nrow = 0)
for(i in agents_keep){
  df1 <- get(paste("coefs_stan_stk_rhat_load",i, sep="."))
  df2 <- rbind(df1,df)
  df <- df2
}
#neff
for(i in agents_keep){
  df1 <- get(paste("coefs_stan_stk_neff_load",i, sep="."))
  df2 <- rbind(df1,df)
  df <- df2
}
write.csv(df, file="data/PosteriorEst_RhatNeff_load_stsp_reduced.csv")


#Plot posteriors
## Plot from file
post_all_load <- read.csv("data/PosteriorEst_IntSlpSigYr_load_stsp_reduced.csv")
post_agents_load <- read.csv("data/load_coefs_stan_stspec_reducedagents.csv")

## Plot Posterior for all agents
jpeg(filename='figs/Fig_LOAD_beta estimates_reduced.jpg', 
     width=480, height=300, quality=75)
ggplot(post_agents_load) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(X, -mid), ymax = X75, ymin = X25), size=1.5, col="darkblue") +
  geom_linerange(aes(x = X, ymax = upper, ymin = lower), col="darkblue") +
  geom_point(aes(x = X, y = mid), size = 3, col="darkblue")+
  labs(x ="Infectious agents", y = "Global effect size", title="Intensity")+
  scale_x_discrete(labels=c("ic_mul" = "I. multifiliis", 
                            "te_mar" = "T. maritinum",
                            "pa_ther" = "P. theridion",
                            "fl_psy" = "F. psychrophilum",
                            "sch" = "Ca. S. salmonis",
                            "te_bry" = "T. bryosalmonae",
                            "pa_kab" = "P. kabatai",
                            "c_b_cys" = "Ca. B. cysticola",
                            "pa_min" = "P. minibicornis",
                            "arena2" = "Arenavirus",
                            "fa_mar" = "F. margolisi",
                            "my_arc" = "M. arcticus",
                            "ven" = "ENV",
                            "ic_hof" = "I. hoferi",
                            "lo_sal" = "L. salmonae",
                            "rlo" = "RLO",
                            "sp_des" = "S. destruens",
                            "ku_thy" = "K. thyrsites",
                            "prv" = "PRV",
                            "pspv" = "PSPV",
                            "ce_sha" = "C. shasta",
                            "pa_pse" = "P. pseudobranchicola",
                            "de_sal" = "D. salmonis"))+
  theme(axis.text.y = element_text(face = "italic"), plot.title = element_text(hjust = 0.5))+
  coord_flip()
dev.off()


########################
## Derive proportions of posterior draws <0 per model - load
prop_less_than_0_load <- matrix(NA,
                                ncol = 1,
                                nrow = length(agents_keep),
                                dimnames = list(agents_keep,"load_prop_neg"))
for (i in agents_keep){
  model<-as.matrix(get(paste("mod.load.",i, sep="")))
  model2<-model[2001:4000,]
  para_name <- colnames(model2)
  temp <- as.matrix((colSums(model2 < 0))/2000)
  prop_less_than_0_load[i,] <- temp[1,]
}
write.csv(prop_less_than_0_load, file="data/Prop posterior neg_load_reduced.csv")


#################################################################################
#################################################################################
## Calculate stock-specific agent estimates - LOAD - 
#create empty martix to bind to
stk.spec.slope.load <- matrix(NA,
                              nrow = 0,
                              ncol = 6,
                              dimnames = list(NULL,c("2.5","25","50","75","97.5","agent")))
## Create a matrix per agent model and loop
for(i in agents_keep){
  model<-get(paste("mod.load",i, sep="."))
  temp2 <- data.frame(model)
  temp3 <- temp2[,grepl("load",names(temp2))] #include only columns with "load" in name
  temp4 <- temp3[,-grep("igma",colnames(temp3))] #remove columns with "igma" in name
  temp5 <- temp4[2001:4000,] #remove warm up iterations
  temp6 <- data.frame(temp5[,1] + temp5[,2:ncol(temp5)]) #add the stock-specific draws to global column
  temp7 <- cbind(temp5[,1],temp6) #bind with global column
  colnames(temp7) <- colnames(temp5) #assign names to columns
  temp8 <- as.matrix(apply(temp7, 2, quantile, probs = c(0.025,0.25,0.50,0.75,0.975)))
  temp9 <- t(temp8) #transpose
  temp10<- as.matrix(cbind(temp9, paste(i)))
  colnames(temp10) <- c("2.5","25","50","75","97.5","agent") #assign names to columns
  stk.spec.slope.load <- rbind(stk.spec.slope.load, temp10)
}
write.csv(stk.spec.slope.load, file="data/Stock specific slopes_load_reduced.csv")


##### READ IN DATA FROM FILE
stspslp.load<-read.csv("data/Stock specific slopes_load_reduced.csv")
stspslp.load$stock <- substr(stspslp.load$X, 18, 28)
stspslp.load$stock <- substr(stspslp.load$stock, 1, nchar(stspslp.load$stock)-1)
stspslp.load$stock <- sub("^$", "Global", stspslp.load$stock)
head(stspslp.load)

## Plot stock-specific slopes - ENV
stk.spec.ven.load <-stspslp.load[stspslp.load$agent=="ven",]
#jpeg(filename='figs/Fig_SSHI ONNE stock sp slope_ven_load.jpg', 
#     width=480, height=500, quality=75)
ggplot(stk.spec.ven.load) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.ven.load[stk.spec.ven.load$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.ven.load[stk.spec.ven.load$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.ven.load[stk.spec.ven.load$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size", title="ENV") +
  coord_flip()
#dev.off()

## Plot stock-specific slopes - lo_sal
stk.spec.lo_sal.load <-stspslp.load[stspslp.load$agent=="lo_sal",]
#jpeg(filename='figs/Fig_SSHI ONNE stock sp slope_lo_sal_load.jpg', 
#     width=480, height=500, quality=75)
ggplot(stk.spec.lo_sal.load) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.lo_sal.load[stk.spec.lo_sal.load$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.lo_sal.load[stk.spec.lo_sal.load$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.lo_sal.load[stk.spec.lo_sal.load$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size", title="ENV") +
  coord_flip()

## Plot stock-specific slopes - pa_min
stk.spec.pa_min.load <-stspslp.load[stspslp.load$agent=="pa_min",]
#jpeg(filename='figs/Fig_SSHI ONNE stock sp slope_pa_min_load.jpg', 
#     width=480, height=500, quality=75)
ggplot(stk.spec.pa_min.load) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.pa_min.load[stk.spec.pa_min.load$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.pa_min.load[stk.spec.pa_min.load$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.pa_min.load[stk.spec.pa_min.load$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size", title="ENV") +
  coord_flip()






####################################################################################
####################################################################################
##SW LOAD WITH SST COVARIATE
##STAN MULTILEVEL MODELING

## Stock-specific SW metric 
### Loop for STAN independent models WITH SST
for(i in agents_keep){
  data <- subset(inf_agt_resid_data_sst, agents_keep==i)
  nam <- paste("mod.load.sst", i, sep = ".")
  assign(nam, stan_lmer(resid_value ~ 0 + load_std + sst_anom + (load_std|Stock) +(1|Year), 
                        data = data,
                        adapt_delta=0.99,
                        REML = F))
}

## Loop to derive coefficient estimates
#load agents
coefs_stan_ls <- matrix(NA,
                       nrow = length(agents_keep),
                       ncol = 5,
                       dimnames = list(agents_keep,c("lower","25","mid","75","upper")))
for(i in agents_keep){
  model<-get(paste("mod.load.sst",i, sep="."))
  ind_coef <- summary(model, 
                      pars = c("load_std"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan_ls[i,] <- ind_coef[1,c(4:8)]
}
write.csv(coefs_stan_ls, file="data/load_coefs_stan_stspec_reduced_agents_SST.csv")

#sst
coefs_stan_lsst <- matrix(NA,
                        nrow = length(agents_keep),
                        ncol = 5,
                        dimnames = list(agents_keep,c("lower","25","mid","75","upper")))
for(i in agents_keep){
  model<-get(paste("mod.load.sst",i, sep="."))
  ind_coef <- summary(model, 
                      pars = c("sst_anom"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan_lsst[i,] <- ind_coef[1,c(4:8)]
}
write.csv(coefs_stan_lsst, file="data/load_coefs_stan_stspec_reduced_sst.csv")

# Load estimates from file and add rownames
coefs_stan <- read.csv("data/load_coefs_stan_stspec_reduced_agents.csv")
rownames(coefs_stan) <- coefs_stan[,1]
coefs_stan <- coefs_stan[,-1] 

# Plot effect size of agents
coefs_order <- coefs_stan[order(-coefs_stan[,3]),]
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))
plotCI(x = coefs_order[,3],
       y = seq(1,length(agents_keep)),
       li = (coefs_order[,1]),
       ui = (coefs_order[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-1.75,2),
       pch = 16,
       scol = "grey")
plotCI(x = coefs_order[,3],
       y = seq(1,length(agents_keep)),
       li = (coefs_order[,2]),
       ui = (coefs_order[,4]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "grey")
text(rep(-1.75,length(agents_keep)), 
     seq(1,length(agents_keep)), 
     labels = rownames(coefs_order), 
     pos = 4,
     font = 2,
     cex=0.95)
axis(1, at = c(-1, -0.5, 0, 0.5, 1, 1.5, 2))
abline(v = 0, lty = 2)
box(col="grey")	
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Intensity",3,line=0.25)


## Loop to derive posteriors by stock
### Intercepts
for(i in agents_keep){
  model<-get(paste("mod.load.sst",i, sep="."))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\(\\Intercept) Stock\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  stocksl <- as.matrix(ind_coef[,1])
  stocksl1 <- rownames(stocksl)
  stocksl2 <- substr(stocksl1, 21, 32)
  stocksl3 <- substr(stocksl2, 1, nchar(stocksl2)-1)
  coefs_stan_stk_int_load.sst <- matrix(NA,
                                    nrow = length(stocksl3),
                                    ncol = 6,
                                    dimnames = list(stocksl3,c("lower","25","mid","75","upper","agent")))
  stocks.mod <- length(stocksl3)
  nam <- paste("coefs_stan_stk_int_load.sst", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_int_load.sst <- cbind(ind_coef[1:stocks.mod,c(4:8)], paste(i))))
}


### Slopes - agents betas
for(i in agents_keep){
  model<-get(paste("mod.load.sst",i, sep="."))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\load_std Stock\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  stocksl <- as.matrix(ind_coef[,1])
  stocksl1 <- rownames(stocksl)
  stocksl2 <- substr(stocksl1, 21, 32)
  stocksl3 <- substr(stocksl2, 1, nchar(stocksl2)-1)
  coefs_stan_stk_slp_load.sst.agent <- matrix(NA,
                                    nrow = length(stocksl3),
                                    ncol = 6,
                                    dimnames = list(stocksl3,c("lower","25","mid","75","upper","agent")))
  stocks.mod <- length(stocksl3)
  nam <- paste("coefs_stan_stk_slp_load.sst.agent", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_slp_load.sst.agent <- cbind(ind_coef[1:stocks.mod,c(4:8)], paste(i))))
}

### Year intercepts
for(i in agents_keep){
  model<-get(paste("mod.load.sst.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = "b\\[\\(\\Intercept) Year",
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  temp1 <- as.matrix(ind_coef[,1])
  temp2 <- rownames(temp1)
  temp3 <- substr(temp2, 20, 32)
  temp4 <- substr(temp3, 1, nchar(temp3)-1)
  coefs_stan_stk_year_load.sst <- matrix(NA,
                                     nrow = length(temp4),
                                     ncol = 6,
                                     dimnames = list(temp4,c("lower","25","mid","75","upper","agent")))
  nam <- paste("coefs_stan_stk_year_load.sst", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_year_load.sst[i] <- cbind(ind_coef[c(1:length(temp4)),c(4:8)], paste(i))))
}

### Sigmas
for(i in agents_keep){
  model<-get(paste("mod.load.sst",i, sep="."))
  ind_coef <- summary(model, 
                      regex_pars = c("sigma","Sigma"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  temp1 <- as.matrix(ind_coef[,1])
  temp2 <- rownames(temp1)
  coefs_stan_stk_sig_load.sst <- matrix(NA,
                                    nrow = length(temp2),
                                    ncol = 6,
                                    dimnames = list(temp2,c("lower","25","mid","75","upper","agent")))
  nam <- paste("coefs_stan_stk_sig_load.sst", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_sig_load.sst[i] <- cbind(ind_coef[c(1:length(temp2)),c(4:8)], paste(i))))
}

### Rhat and Neff
#extract parameter names
for(i in agents_keep){
  model<-get(paste("mod.load.sst",i, sep="."))
  ind_coef <- as.matrix(summary(model, 
                                digits = 3) [,"Rhat"])
  coefs_stan_stk_rhat_load.sst <- matrix(NA,
                                     nrow = nrow(ind_coef),
                                     ncol = 2,
                                     dimnames = list(rownames(ind_coef),c("Rhat","agent")))
  nam <- paste("coefs_stan_stk_rhat_load.sst", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_rhat_load.sst[i] <- cbind(ind_coef[c(1:nrow(ind_coef)),1], paste(i), paste("Rhat"))))
}

for(i in agents_keep){
  model<-get(paste("mod.load.sst",i, sep="."))
  ind_coef <- as.matrix(summary(model) [,"n_eff"])
  coefs_stan_stk_neff_load.sst <- matrix(NA,
                                     nrow = nrow(ind_coef),
                                     ncol = 2,
                                     dimnames = list(rownames(ind_coef),c("Neff","agent")))
  nam <- paste("coefs_stan_stk_neff_load.sst", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_neff_load.sst[i] <- cbind(ind_coef[c(1:nrow(ind_coef)),1], paste(i), paste("Neff"))))
}

### rbind all files
#intercept
df <- matrix(ncol = 6, nrow = 0)
colnames(df) <- colnames(coefs_stan_stk_int_load.sst.c_b_cys)
for(i in agents_keep){
  df1 <- get(paste("coefs_stan_stk_int_load.sst",i, sep="."))
  df2 <- rbind(df1,df)
  df <- df2
}
#slopes
for(i in agents_keep){
  df1 <- get(paste("coefs_stan_stk_slp_load.sst.agent",i, sep="."))
  df2 <- rbind(df1,df)
  df <- df2
}
#sigmas
for(i in agents_keep){
  df1 <- get(paste("coefs_stan_stk_sig_load.sst",i, sep="."))
  df2 <- rbind(df1,df)
  df <- df2
}
#year
for(i in agents_keep){
  df1 <- get(paste("coefs_stan_stk_year_load.sst",i, sep="."))
  df2 <- rbind(df1,df)
  df <- df2
}
write.csv(df, file="data/PosteriorEst_IntSlpSigYr_load_stsp_reduced_SST.csv")


#### Rbind all convergence parameters and write as csv
#rhat
df <- matrix(ncol = 3, nrow = 0)
for(i in agents_keep){
  df1 <- get(paste("coefs_stan_stk_rhat_load.sst",i, sep="."))
  df2 <- rbind(df1,df)
  df <- df2
}
#neff
for(i in agents_keep){
  df1 <- get(paste("coefs_stan_stk_neff_load.sst",i, sep="."))
  df2 <- rbind(df1,df)
  df <- df2
}
write.csv(df, file="data/PosteriorEst_RhatNeff_load_stsp_reduced_SST.csv")


#Plot posteriors
## Plot from file
post_all_load_sst <- read.csv("data/PosteriorEst_IntSlpSigYr_load_stsp_reduced_SST.csv")
post_agents_load_sst <- read.csv("data/load_coefs_stan_stspec_reduced_agents.csv")
post_anom_load_sst <- read.csv("data/load_coefs_stan_stspec_reduced_sst.csv")

## Plot Posterior for all agents
jpeg(filename='figs/Fig_LOAD_beta estimates_reduced_agents.jpg', 
     width=480, height=300, quality=75)
ggplot(post_agents_load_sst) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(X, -mid), ymax = X75, ymin = X25), size=1.5, col="darkblue") +
  geom_linerange(aes(x = X, ymax = upper, ymin = lower), col="darkblue") +
  geom_point(aes(x = X, y = mid), size = 3, col="darkblue")+
  labs(x ="Infectious agents", y = "Global effect size", title="Intensity")+
  scale_x_discrete(labels=c("ic_mul" = "I. multifiliis", 
                            "te_mar" = "T. maritinum",
                            "pa_ther" = "P. theridion",
                            "fl_psy" = "F. psychrophilum",
                            "sch" = "Ca. S. salmonis",
                            "te_bry" = "T. bryosalmonae",
                            "pa_kab" = "P. kabatai",
                            "c_b_cys" = "Ca. B. cysticola",
                            "pa_min" = "P. minibicornis",
                            "arena2" = "Arenavirus",
                            "fa_mar" = "F. margolisi",
                            "my_arc" = "M. arcticus",
                            "ven" = "ENV",
                            "ic_hof" = "I. hoferi",
                            "lo_sal" = "L. salmonae",
                            "rlo" = "RLO",
                            "sp_des" = "S. destruens",
                            "ku_thy" = "K. thyrsites",
                            "prv" = "PRV",
                            "pspv" = "PSPV",
                            "ce_sha" = "C. shasta",
                            "pa_pse" = "P. pseudobranchicola",
                            "de_sal" = "D. salmonis")) +
  theme(axis.text.y = element_text(face = "italic"), plot.title = element_text(hjust = 0.5))+
  coord_flip()
dev.off()


########################
## Derive proportions of posterior draws <0 per model - load
prop_less_than_0_load.sst <- matrix(NA,
                                ncol = 1,
                                nrow = length(agents_keep),
                                dimnames = list(agents_keep,"loadsst_prop_neg"))
for (i in agents_keep){
  model<-as.matrix(get(paste("mod.load.sst.",i, sep="")))
  model2<-model[2001:4000,]
  #para_name <- colnames(model2)
  temp <- as.matrix((colSums(model2 < 0))/2000)
  prop_less_than_0_load.sst[i,] <- temp[1,]
}
write.csv(prop_less_than_0_load.sst, file="data/Prop posterior neg_load_reduced_SST.csv")




#################################################################################
#################################################################################
## Calculate stock-specific agent estimates - LOAD - 
#create empty martix to bind to
stk.spec.slope.load.sst <- matrix(NA,
                              nrow = 0,
                              ncol = 6,
                              dimnames = list(NULL,c("2.5","25","50","75","97.5","agent")))
## Create a matrix per agent model and loop
for(i in agents_keep){
  model<-get(paste("mod.load.sst",i, sep="."))
  temp2 <- data.frame(model)
  temp3 <- temp2[,grepl("load",names(temp2))] #include only columns with "load" in name
  temp4 <- temp3[,-grep("igma",colnames(temp3))] #remove columns with "igma" in name
  temp5 <- temp4[2001:4000,] #remove warm up iterations
  temp6 <- data.frame(temp5[,1] + temp5[,2:ncol(temp5)]) #add the stock-specific draws to global column
  temp7 <- cbind(temp5[,1],temp6) #bind with global column
  colnames(temp7) <- colnames(temp5) #assign names to columns
  temp8 <- as.matrix(apply(temp7, 2, quantile, probs = c(0.025,0.25,0.50,0.75,0.975)))
  temp9 <- t(temp8) #transpose
  temp10<- as.matrix(cbind(temp9, paste(i)))
  colnames(temp10) <- c("2.5","25","50","75","97.5","agent") #assign names to columns
  stk.spec.slope.load.sst <- rbind(stk.spec.slope.load.sst, temp10)
}
write.csv(stk.spec.slope.load.sst, file="data/Stock specific slopes_load_reduced_SST.csv")


##### READ IN DATA FROM FILE
stspslp.load.sst <- read.csv("data/Stock specific slopes_load_reduced_SST.csv")
stspslp.load.sst$stock <- substr(stspslp.load.sst$X, 18, 28)
stspslp.load.sst$stock <- substr(stspslp.load.sst$stock, 1, nchar(stspslp.load.sst$stock)-1)
stspslp.load.sst$stock <- sub("^$", "Global", stspslp.load.sst$stock)
head(stspslp.load.sst)






#######################################################################################
#######################################################################################
# STAN Approach for Multi-level Modeling
# FW PREVALENCE - INDEPENDENT MODELS by AGENT
## Stock-specific metric
### Create files for each agent
agents.fw <- unique(inf_agt_resid_data_fw$agent)
stocks.fw <- unique(inf_agt_resid_data_fw$Stock_Analysis)
for(i in unique(inf_agt_resid_data_fw$agent)) {
  nam <- paste("df", i, sep = ".")
  assign(nam, inf_agt_resid_data_fw[inf_agt_resid_data_fw$agent==i,])
}

### Loop for STAN independent models by agent
for(i in agents.fw){
  data <- subset(inf_agt_resid_data_fw, agent==i)
  nam <- paste("mod.fw", i, sep = ".")
  assign(nam, stan_lmer(resid_value ~ 0 +  prev_std + (prev_std|Stock_Analysis) +(1|Year), 
                        data = data,
                        adapt_delta=0.99,
                        REML = F))
}

# uninformed priors - for example:
prior_summary(mod.fw.ic_mul)
plot(mod.fw.ic_mul, "ess")
plot(mod.fw.ic_mul, "trace")

## Shiny Stan launch
#my_sso.arena2 <- launch_shinystan(mod.fw.ic_mul)

## Derive coefficient estimates and save in .csv file
coefs_stan_fw <- matrix(NA,
                     nrow = length(agents.fw),
                     ncol = 5,
                     dimnames = list(agents.fw,c("lower","25","mid","75","upper")))
for(i in agents.fw){
  model<-get(paste("mod.fw.",i, sep=""))
  ind_coef_fw <- summary(model, 
                      pars = c("prev_std"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  coefs_stan_fw[i,] <- ind_coef_fw[1,c(4:8)]
}
write.csv(coefs_stan_fw, file="data/prev_coefs_stan_stspec_fw.csv")

# Load estimates from file (if not running full model) and assign rownames
coefs_stan_fw <- read.csv("data/prev_coefs_stan_stspec_fw.csv")
rownames(coefs_stan_fw) <- coefs_stan_fw[,1]
coefs_stan_fw <- coefs_stan_fw[,-1]  

# Plot effect size per agent
coefs_order <- coefs_stan_fw[order(-coefs_stan_fw[,3]),]

jpeg(filename='figs/Fig_prev_coefs_stan_stspec_fw.jpg', 
     width=480, height=250, quality=75)
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))
plotCI(x = coefs_order[,3],
       y = seq(1,length(agents.fw)),
       li = (coefs_order[,1]),
       ui = (coefs_order[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-1.5,1),
       pch = 16,
       scol = "grey")
plotCI(x = coefs_order[,3],
       y = seq(1,length(agents.fw)),
       li = (coefs_order[,2]),
       ui = (coefs_order[,4]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "grey")
text(rep(-1.5,length(agents.fw)), 
     seq(1,length(agents.fw)), 
     labels = rownames(coefs_order), 
     pos = 4,
     font = 2,
     cex=0.95)
axis(1, at = c(-1.5, -1, -0.5, 0, 0.5, 1))
abline(v = 0, lty = 2)
box(col="grey")	
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("FW Prevalence",3,line=0.25)
dev.off()

## Derive posterior estimates by stock
### Intercepts
coefs_stan_stk_int.fw <- matrix(NA,
                             nrow = length(stocks.fw),
                             ncol = 6,
                             dimnames = list(stocks.fw,c("lower","25","mid","75","upper","agent")))
for(i in agents.fw){
  model<-get(paste("mod.fw.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\(\\Intercept) Stock_Analysis\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_int.fw", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_int.fw <- cbind(ind_coef[c(1:14),c(4:8)], paste(i))))
}

### Slopes
coefs_stan_stk_slp.fw <- matrix(NA,
                             nrow = length(stocks.fw),
                             ncol = 6,
                             dimnames = list(stocks.fw,c("lower","25","mid","75","upper","agent")))
for(i in agents.fw){
  model<-get(paste("mod.fw.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("b\\[\\prev_std Stock_Analysis\\:"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_slp.fw", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_slp.fw <- cbind(ind_coef[c(1:14),c(4:8)], paste(i))))
}

### Year intercepts
years.fw <- unique(inf_agt_resid_data_fw$Year)
coefs_stan_stk_year.fw <- matrix(NA,
                              nrow = length(years.fw),
                              ncol = 6,
                              dimnames = list(years.fw,c("lower","25","mid","75","upper","agent")))
for(i in agents.fw){
  model<-get(paste("mod.fw.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = "Year",
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_year.fw", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_year.fw <- cbind(ind_coef[c(1:5),c(4:8)], paste(i))))
}

### Sigmas
coefs_stan_stk_sig.fw <- matrix(NA,
                             nrow = length(stocks.fw),
                             ncol = 6,
                             dimnames = list(stocks.fw,c("lower","25","mid","75","upper","agent")))
for(i in agents.fw){
  model<-get(paste("mod.fw.",i, sep=""))
  ind_coef <- summary(model, 
                      regex_pars = c("sigma","Sigma"),
                      probs = c(0.025,0.25,0.5,0.75, 0.975),
                      digits = 2)
  nam <- paste("coefs_stan_stk_sig.fw", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_sig.fw <- cbind(ind_coef[c(1:5),c(4:8)], paste(i))))
}

### Rhat and Neff
sims.fw <-as.matrix(mod.fw.c_b_cys) #extract parameter names from any agent model 
dim(sims.fw)
para_name.fw <- c(colnames(sims.fw), "mean_PPD", "log-posterior")
para_name.fw

#### Rhat
coefs_stan_stk_rhat.fw <- matrix(NA,
                              nrow = length(para_name.fw),
                              ncol = 2,
                              dimnames = list(para_name.fw,c("Rhat","agent")))
for(i in agents.fw){
  model<-get(paste("mod.fw",i, sep="."))
  ind_coef <- as.matrix(summary(model, 
                                digits = 3) [,"Rhat"])
  nam <- paste("coefs_stan_stk_rhat.fw", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_rhat.fw <- cbind(ind_coef[c(1:41),1], paste(i), paste("Rhat"))))
}

#### Neff
coefs_stan_stk_neff.fw <- matrix(NA,
                              nrow = length(para_name.fw),
                              ncol = 2,
                              dimnames = list(para_name.fw,c("Neff","agent")))
for(i in agents.fw){
  model<-get(paste("mod.fw",i, sep="."))
  ind_coef <- as.matrix(summary(model) [,"n_eff"])
  nam <- paste("coefs_stan_stk_neff.fw", i, sep = ".")
  as.matrix(assign(nam, coefs_stan_stk_neff.fw <- cbind(ind_coef[c(1:41),1], paste(i), paste("Neff"))))
}


## Derive proportions of posterior draws <0 per model - FW prevalence
param.fw<-colnames(data.frame(mod.fw.c_b_cys)) #create object of parameters in model
param.prop0.fw <- matrix(NA,
                      ncol = length(param.fw),
                      nrow = length(agents.fw),
                      dimnames = list(agents.fw,param.fw))

for (i in agents.fw){
  model<-as.matrix(get(paste("mod.fw.",i, sep="")))
  model2<-model[2001:4000,]
  param.prop0.fw[i,] <- (colSums(model2 < 0))/2000
}
write.csv(param.prop0.fw, file="data/Percent post draws >0_prev_stspec_FW.csv")
propzero.fw <- read.csv("data/Percent post draws >0_prev_stspec_FW.csv")

jpeg(filename='figs/Fig_prop<0_pa_ther_stspec_FW.jpg', 
     width=480, height=500, quality=75)
ggplot(propzero.fw) +
  geom_bar(stat="identity", aes(reorder(X, prev_std), prev_std), fill="gray", col="black") +
  ylim(0,1)+
  labs(x="Agent", y="Proportion of total", title="FW Posterior estimates <0") +
  coord_flip()
dev.off()



###############################################################################
## BRING ALL MODELING DATA TOGETHER
#beta estimates
prev_beta <- read.csv("data/prev_coefs_stan_stspec.csv")
prev_beta <- prev_beta[!prev_beta$X=="arena2",]
load_beta <- read.csv("data/load_coefs_stan_stspec_reducedagents.csv")
prev.sst_beta <- read.csv("data/prev_coefs_stan_stspec.sst.csv")
prev.sst_beta <- prev.sst_beta[!prev.sst_beta$X=="arena2",]
load.sst_beta <- read.csv("data/load_coefs_stan_stspec_reduced_agents_SST.csv")
anom.prev_beta <- read.csv("data/prev_coefs_stan_stspec.sst_SST coefs.csv")
anom.prev_beta <- anom.prev_beta[!anom.prev_beta$X=="arena2",]
anom.load_beta <- read.csv("data/load_coefs_stan_stspec_reduced_sst.csv")
prev_beta_fw <- read.csv("data/prev_coefs_stan_stspec_fw.csv")

#proportion less than 0
temp <- read.csv("data/Percent post draws >0_prev_stspec.csv")
prev_neg <- temp[,1:2]
prev_neg <- prev_neg[!prev_neg$X=="arena2",]
load_neg <- read.csv("data/Prop posterior neg_load_reduced.csv")
prev.sst_neg <- read.csv("data/Prop posterior neg_prevSST.csv")
prev.sst_neg <- prev.sst_neg[prev.sst_neg$X!="arena2" & prev.sst_neg$X!="smallUK",]
load.sst_neg <- read.csv("data/Prop posterior neg_load_reduced_SST.csv")
temp <- read.csv("data/Percent post draws >0_prev_stspec_FW.csv")
prev_neg_fw <- data.frame(temp[,1:2])
  
#sort all df to align with agent betas
prev_beta <-prev_beta[order(-prev_beta$mid),]
prev.sst_beta <- prev.sst_beta[order(match(prev.sst_beta$X, prev_beta$X)),]
anom.prev_beta <- anom.prev_beta[order(match(anom.prev_beta$X, prev_beta$X)),]
prev_neg <- prev_neg[order(match(prev_neg$X, prev_beta$X)),]

load_beta <-load_beta[order(-load_beta$mid),]
load.sst_beta <- load.sst_beta[order(match(load.sst_beta$X, load_beta$X)),]
anom.load_beta <- anom.load_beta[order(match(anom.load_beta$X, load_beta$X)),]
load_neg <- load_neg[order(match(load_neg$X, load_beta$X)),]

prev_beta_fw <-prev_beta_fw[order(-prev_beta_fw$mid),]
prev_neg_fw <- prev_neg_fw[order(match(prev_neg_fw$X, prev_beta_fw$X)),]

#Extract agent names into df - prevalence
agent.names <- data.frame(unique(inf_agt_resid_data[c("agent", "plot.agent")]))
agent.names <- agent.names[!agent.names$agent=="arena2",]
agent.names$plot.agent <- as.character(agent.names$plot.agent)
agent.names$plot.agent[agent.names$plot.agent == "VENV"] <- "ENV"
agent.names <- agent.names[order(match(agent.names$agent, prev_beta$X)),]

#Extract agent names into df - load
load.names1 <- data.frame(load_beta[,1])
colnames(load.names1) <- "agent"
load.names2 <- merge(agent.names, load.names1, by.x="agent", by.y="agent")
load.names <- load.names2[order(match(load.names2$agent, load_beta$X)),]

#Extract agent names into df - fw prev
fw_names1 <- data.frame(prev_beta_fw[,1])
colnames(fw_names1) <- "agent"
fw_names2 <- merge(agent.names, fw_names1, by.x="agent", by.y="agent")
fw.names <- fw_names2[order(match(fw_names2$agent, prev_beta_fw$X)),]

jpeg(filename='figs/Fig_Beta estimates_agent_agentSST_prev.jpg', 
     width=510, height=600, quality=300)
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))
#add coefficient estimates for the sst cofactor
plotCI(x = anom.prev_beta[,4], col="red",
       y = seq(1,nrow(prev_beta)),
       li = (anom.prev_beta[,2]),
       ui = (anom.prev_beta[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-2.5,1.5),
       pch = 16,
       scol = "red")
plotCI(x = anom.prev_beta[,4], col="red",
       y = seq(1,nrow(prev_beta)),
       li = (anom.prev_beta[,3]),
       ui = (anom.prev_beta[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "red")
#add coefficient estimates with sst cofactor
plotCI(x = prev.sst_beta[,4], col="gray",
       y = seq(1,nrow(prev_beta)),
       li = (prev.sst_beta[,2]),
       ui = (prev.sst_beta[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-2.5,1.5),
       pch = 16,
       add = TRUE,
       scol = "gray")
plotCI(x = prev.sst_beta[,4], col="gray",
       y = seq(1,nrow(prev_beta)),
       li = (prev.sst_beta[,3]),
       ui = (prev.sst_beta[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       lwd = 3,
       add = TRUE,
       scol = "gray")
#add coefficient estimates without SST cofactor
plotCI(x = prev_beta[,4], col="black",
       y = seq(1,nrow(prev_beta))+.25,
       li = (prev_beta[,2]),
       ui = (prev_beta[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-2.5,1.5),
       pch = 16,
       add = TRUE,
       scol = "black")
plotCI(x = prev_beta[,4], col="black",
       y = seq(1,nrow(prev_beta))+.25,
       li = (prev_beta[,3]),
       ui = (prev_beta[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "black")
text(rep(-2.5,nrow(prev_beta)), 
     seq(1,nrow(prev_beta)), 
     labels = (agent.names[,2]), 
     pos = 4,
     font = 3,
     cex=0.95)
legend(x=.5, # x coordinate of the top left of the legend
       y=21,
       legend=c("Agent", "Agent+SST", "SST"), 
       pch=21,
       pt.bg=c("black","gray","red"))
axis(1, at = c(-2, -1.5,-1, -0.5, 0, 0.5, 1, 1.5))
abline(v = 0, lty = 2)
box(col="grey") 
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Marine Prevalence",3,line=0.25)
dev.off()


## LOAD PLOT
jpeg(filename='figs/Fig_Beta estimates_agent_agentSST_load.jpg', 
     width=510, height=400, quality=300)
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))
#add coefficient estimates for the sst cofactor
plotCI(x = anom.load_beta[,4], col="red",
       y = seq(1,nrow(load_beta)),
       li = (anom.load_beta[,2]),
       ui = (anom.load_beta[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-2.5,1.5),
       pch = 16,
       scol = "red")
plotCI(x = anom.load_beta[,4], col="red",
       y = seq(1,nrow(load_beta)),
       li = (anom.load_beta[,3]),
       ui = (anom.load_beta[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "red")
#add coefficient estimates with sst cofactor
plotCI(x = load.sst_beta[,4], col="gray",
       y = seq(1,nrow(load_beta)),
       li = (load.sst_beta[,2]),
       ui = (load.sst_beta[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-2.5,1.5),
       pch = 16,
       add = TRUE,
       scol = "gray")
plotCI(x = load.sst_beta[,4], col="gray",
       y = seq(1,nrow(load_beta)),
       li = (load.sst_beta[,3]),
       ui = (load.sst_beta[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       lwd = 3,
       add = TRUE,
       scol = "gray")
#add coefficient estimates without SST cofactor
plotCI(x = load_beta[,4], col="black",
       y = seq(1,nrow(load_beta))+.25,
       li = (load_beta[,2]),
       ui = (load_beta[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-2.5,1.5),
       pch = 16,
       add = TRUE,
       scol = "black")
plotCI(x = load_beta[,4], col="black",
       y = seq(1,nrow(load_beta))+.25,
       li = (load_beta[,3]),
       ui = (load_beta[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       add = TRUE,
       lwd = 3,
       scol = "black")
text(rep(-2.5,nrow(load_beta)), 
     seq(1,nrow(load_beta)), 
     labels = (load.names[,2]), 
     pos = 4,
     font = 3,
     cex=0.95)
legend(x=.5, # x coordinate of the top left of the legend
       y=21,
       legend=c("Agent", "Agent+SST", "SST"), 
       pch=21,
       pt.bg=c("black","gray","red"))
axis(1, at = c(-2, -1.5,-1, -0.5, 0, 0.5, 1, 1.5))
abline(v = 0, lty = 2)
box(col="grey") 
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Marine Intensity",3,line=0.25)
dev.off()


### Freshwater prevalence fig ##
jpeg(filename='figs/Fig_Beta estimates_agent_agent_prev_FW.jpg', 
     width=510, height=300, quality=300)
par(mfrow=c(1,1), mar=c(3,1,1,1),oma=c(0.5,0.5,0.5,0.5))
#add coefficient estimates with sst cofactor
plotCI(x = prev_beta_fw[,4], col="black",
       y = seq(1,nrow(prev_beta_fw)),
       li = (prev_beta_fw[,2]),
       ui = (prev_beta_fw[,6]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "",
       xlim = c(-2.5,1.5),
       pch = 16,
       scol = "black")
plotCI(x = prev_beta_fw[,4], col="black",
       y = seq(1,nrow(prev_beta_fw)),
       li = (prev_beta_fw[,3]),
       ui = (prev_beta_fw[,5]),
       err = "x",
       sfrac = 0 ,
       gap = 0,
       pch = 16,
       lwd = 3,
       add = TRUE,
       scol = "black")
text(rep(-2.5,nrow(prev_beta_fw)), 
     seq(1,nrow(prev_beta_fw)), 
     labels = (fw.names[,2]), 
     pos = 4,
     font = 3,
     cex=0.95)
axis(1, at = c(-2, -1.5,-1, -0.5, 0, 0.5, 1, 1.5))
abline(v = 0, lty = 2)
box(col="grey") 
mtext("Effect size",1,line=2.2, cex=1.1)
mtext("Freshwater Prevalence",3,line=0.25)
dev.off()


#### Table for model parameters
prev_beta$term <- "prev.beta"
prev_beta$model <- "prev"
prev.sst_beta$term <- "prev.beta"
prev.sst_beta$model <- "prev+sst"
anom.prev_beta$term <- "sst.beta"
anom.prev_beta$model <- "prev+sst"
prev_beta_fw$term <- "prev.beta"
prev_beta_fw$model <- "prevfw"

prev.mod.param <- cbind(prev_beta[,c(1,2,4,6)], prop_neg = prev_neg[,2],"prev")
colnames(prev.mod.param) <- c("agent","prev2.5","prev50","prev97.5","prop.neg","model")
prev.sst.mod.param <- cbind(prev.sst_beta[,c(1,2,4,6)], anom.prev_beta[,c(2,4,6)],"prev+sst")
colnames(prev.sst.mod.param) <- c("agent","prev.sst2.5","prev.sst50","prev.sst97.5","sst2.5","sst50","sst97.5","model")
prev.all.mod <- cbind(prev.mod.param, prev.sst.mod.param[,-c(1)])
write.csv(prev.all.mod, file="data/Prev_model_results_table.csv")

load.mod.param <- cbind(load_beta[,c(1,2,4,6)], prop_neg = load_neg[,2],"load")
colnames(load.mod.param) <- c("agent","load2.5","load50","load97.5","prop.neg","model")
load.sst.mod.param <- cbind(load.sst_beta[,c(1,2,4,6)], anom.load_beta[,c(2,4,6)],"load+sst")
colnames(load.sst.mod.param) <- c("agent","load.sst2.5","load.sst50","load.sst97.5","sst2.5","sst50","sst97.5","model")
load.all.mod <- cbind(load.mod.param, load.sst.mod.param[,-c(1)])
write.csv(load.all.mod, file="data/Load_model_results_table.csv")

prev.mod.param_fw <- cbind(prev_beta_fw[,c(1,2,4,6)], prop_neg = prev_neg_fw[,2],"prev")
colnames(prev.mod.param_fw) <- c("agent","prev2.5","prev50","prev97.5","prop.neg","model")
write.csv(prev.mod.param_fw, file="data/Prev_FW_model_results_table.csv")

temp <- unique(inf_agt_resid_data[c("agent", "plot.agent")])
write.csv(temp, file="data/agent names.csv")


