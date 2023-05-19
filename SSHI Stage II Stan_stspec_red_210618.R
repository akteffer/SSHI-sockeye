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
inf_agt_resid_data <- read.csv("data/REDUCED_ONNE_agents only_SW_220126.csv")
head(inf_agt_resid_data)
sw_rr <- read.csv("data/RR_ONNE_cumulative metrics_SW_220126.csv")
head(sw_rr)
inf_agt_resid_data <- rbind(inf_agt_resid_data, sw_rr)

##### FW
inf_agt_resid_data_fw <- read.csv("data/REDUCED_ONNE_agents only_FW_220126.csv")
head(inf_agt_resid_data_fw)
fw_rr <- read.csv("data/RR_ONNE_cumulative metrics_FW_220126.csv")
head(fw_rr)
inf_agt_resid_data_fw <- rbind(inf_agt_resid_data_fw, fw_rr)

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

# Add Stock column, remove smallUK, add plot.agent column (full name)
inf_agt_resid_data$Stock <- inf_agt_resid_data$Stock_Analysis
inf_agt_resid_data <- merge(inf_agt_resid_data, agent.names, by = c("agent"))
head(inf_agt_resid_data)
dim(inf_agt_resid_data)


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
#temp <- inf_agt_resid_data_fw[inf_agt_resid_data_fw$agent=="rib",]
#ggplot(temp) +
#  geom_point(aes(x=prev, y=prev_std), pch=21)
#ggplot(inf_agt_resid_data_fw) +
#  geom_point(aes(x=prev, y=prev_std, col=agent), pch=21)+
#  stat_smooth(aes(x=prev, y=prev_std, col=agent), method="lm", se=F, size=.5)

inf_agt_resid_data_fw$Stock <- inf_agt_resid_data_fw$Stock_Analysis
inf_agt_resid_data_fw <- merge(inf_agt_resid_data_fw, agent.names, by = c("agent"))
head(inf_agt_resid_data_fw)
dim(inf_agt_resid_data_fw)

# Create objects for analysis
agents <- unique(inf_agt_resid_data$agent)
agents.fw <- unique(inf_agt_resid_data_fw$agent)
agents.rr <- c("richness","rib")
stocks <- unique(inf_agt_resid_data$Stock_Analysis)
years <- unique(inf_agt_resid_data$Year)
brdyears <- unique(inf_agt_resid_data$brood_year)


## Write csv's of final data frames for tables
write.csv(inf_agt_resid_data, file="data/inf_agent_resid_data_SW_230125.csv")
write.csv(inf_agt_resid_data_fw, file="data/inf_agent_resid_data_FW_230125.csv")



## Total assays run per sample
### SW
sw.N.assays <- data.frame(inf_agt_resid_data %>% 
                            group_by(agent) %>%
                            summarize(assays = sum(N)))
write.csv(sw.N.assays, file="figs/Final data frames and tables/SW_Nassays_230123.csv")
### FW
fw.N.assays <- data.frame(inf_agt_resid_data_fw %>% 
                            group_by(agent) %>%
                            summarize(assays = sum(N)))
write.csv(fw.N.assays, file="figs/Final data frames and tables/FW_Nassays_230123.csv")

### Data Exploration ###
## SW - When were stocks sampled?
sample.totals.sw<-aggregate(inf_agt_resid_data$N, by=list(Category=inf_agt_resid_data$Stock, Year=inf_agt_resid_data$Year), FUN=sum)
jpeg(filename='figs/Fig_Total samples by stock year_SW_230123.jpg', 
     width=480, height=500, quality=300)
ggplot(data=sample.totals.sw, aes(x=reorder(Category, x), y=x, fill=factor(Year))) +
  geom_col() +
  coord_flip() +
  xlab("Stocks")+
  ylab("Assays")
dev.off()

## FW - When were stocks sampled?
sample.totals.fw<-aggregate(inf_agt_resid_data_fw$N, by=list(Category=inf_agt_resid_data_fw$Stock, Year=inf_agt_resid_data_fw$Year), FUN=sum)
jpeg(filename='figs/Fig_Total samples by stock year_FW_230123.jpg', 
     width=480, height=500, quality=75)
ggplot(data=sample.totals.fw, aes(x=reorder(Category, x), y=x, fill=factor(Year))) +
  geom_col() +
  coord_flip() +
  xlab("Stocks")+
  ylab("Assays")
dev.off()

## Plot total SW detections of each agent by year - note variable prevalence across agents and years
samplesperagent.sw<-inf_agt_resid_data %>% 
  group_by(agent, Year) %>%
  summarise(posdet = sum(posdet))
jpeg(filename='figs/Fig_Total detections by year_SW_230123.jpg', 
     width=480, height=500, quality=75)
ggplot(data=samplesperagent.sw, aes(x=reorder(agent, posdet), y=posdet, fill=factor(Year)))+
  geom_bar(stat="identity")+
  labs(fill="Sampling year") +
  coord_flip()+
  xlab("Infectious agents")+
  ylab("Positive detections")
dev.off()

## Plot total FW detections of each agent by year - note variable prevalence across agents and years
samplesperagent.fw<-inf_agt_resid_data_fw %>% 
  group_by(agent, Year) %>%
  summarise(posdet = sum(posdet)) 
jpeg(filename='figs/Fig_Total detections by year_FW_230123.jpg', 
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
jpeg(filename='figs/Fig_detections by stock_proportion_SW_230123.jpg', 
     width=500, height=300, quality=300)
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
jpeg(filename='figs/Fig_detections by stock_proportion_FW_230123.jpg', 
     width=430, height=300, quality=75)
ggplot(data=samplesperagent.fwst, aes(x=reorder(Stock_Analysis, det), y=det, fill=agent))+
  geom_col(position="fill") +
  labs(fill="Sampling year") +
  coord_flip() +
  xlab("Stocks")+
  ylab("Proportion of detections")+
  scale_y_continuous(labels = scales::percent)
dev.off()

### By sampling year in FW (only 6: 2008, 2011-2015)
samplesperagent.fwsy<-inf_agt_resid_data_fw %>% 
  group_by(agent, Year) %>%
  summarise(det = sum(posdet)) %>%
  mutate(freq = det / sum(det)) #proportion of pos detections of total pos detecttions in this Stock
jpeg(filename='figs/Fig_detections by year_FW_230123.jpg', 
     width=430, height=300, quality=75)
ggplot(data=samplesperagent.fwsy, aes(x=factor(Year), y=det, fill=agent))+
  geom_col() +
  labs(fill="Sampling year") +
  coord_flip() +
  xlab("Year")+
  ylab("Proportion of detections")+
  scale_y_continuous(labels = scales::percent)
dev.off()

# Plot raw data by: 
## Prevalence SW
jpeg(filename='figs/Fig_SW_Raw data by year stock_prev_230123.jpg', 
     width=1200, height=700, quality=300)
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
jpeg(filename='figs/Fig_SW_Raw data by year stock_load_230123.jpg', 
     width=1100, height=700, quality=300)
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
jpeg(filename='figs/Fig_FW_Raw data by year stock_prev_230123.jpg', 
     width=700, height=700, quality=300)
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
jpeg(filename='figs/Fig_FW_Raw data by year stock_load_230123.jpg', 
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


#######################################################################################
#######################################################################################
# STAN Approach for Multi-level Modeling
#######################################################################################

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
prior_summary(mod.c_b_cys)
plot(mod.te_mar, "ess")
plot(mod.pa_ther, "trace")
summary(mod.c_b_cys)

## Shiny Stan launch
my_sso.Circo.virus <- launch_shinystan(mod.Circo.virus)
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
write.csv(coefs_stan, file="data/SW_prev_coefs_stan_stspec_230123.csv")
write.csv(coefs_stan, file="figs/Final data frames and tables/SW_prev_coefs_stan_stspec_230123.csv")

# Load estimates from file (if not running full model) and assign rownames
coefs_stan <- read.csv("data/SW_prev_coefs_stan_stspec_230123.csv")
rownames(coefs_stan) <- coefs_stan[,1]
coefs_stan <- coefs_stan[,-1]  

# Plot effect size per agent
coefs_order <- coefs_stan[order(-coefs_stan[,3]),]

jpeg(filename='figs/Fig_SW_prev_coefs_stan_stspec_220123.jpg', 
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
  as.matrix(assign(nam, coefs_stan_stk_int <- cbind(ind_coef[c(1:dim(ind_coef)[1]),c(4:8)], paste(i))))
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
  as.matrix(assign(nam, coefs_stan_stk_slp <- cbind(ind_coef[c(1:dim(ind_coef)[1]),c(4:8)], paste(i))))
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
  as.matrix(assign(nam, coefs_stan_stk_year <- cbind(ind_coef[c(1:dim(ind_coef)[1]),c(4:8)], paste(i))))
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
  as.matrix(assign(nam, coefs_stan_stk_sig <- cbind(ind_coef[c(1:dim(ind_coef)[1]),c(4:8)], paste(i))))
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
  as.matrix(assign(nam, coefs_stan_stk_rhat <- cbind(ind_coef[c(1:dim(ind_coef)[1]),1], paste(i), paste("Rhat"))))
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
  as.matrix(assign(nam, coefs_stan_stk_neff <- cbind(ind_coef[c(1:dim(ind_coef)[1]),1], paste(i), paste("Neff"))))
}

#### Rbind all posteriors and save as .csv file
post_int_slpintsig <- rbind(#coefs_stan_stk_int.arena1,
                            coefs_stan_stk_int.arena2,
                            #coefs_stan_stk_int.ascv,
                            coefs_stan_stk_int.c_b_cys,
                            coefs_stan_stk_int.ce_sha,
                            coefs_stan_stk_int.Circo.virus,
                            coefs_stan_stk_int.de_sal,
                            coefs_stan_stk_int.fa_mar,
                            coefs_stan_stk_int.fl_psy,
                            coefs_stan_stk_int.ic_hof,
                            coefs_stan_stk_int.ic_mul,
                            coefs_stan_stk_int.IcD,
                            coefs_stan_stk_int.ihnv,
                            coefs_stan_stk_int.ku_thy,
                            coefs_stan_stk_int.lo_sal,
                            #coefs_stan_stk_int.mo_vis,
                            coefs_stan_stk_int.my_arc,
                            coefs_stan_stk_int.nu_sal,
                            coefs_stan_stk_int.pa_kab,
                            coefs_stan_stk_int.pa_min,
                            coefs_stan_stk_int.pa_pse,
                            coefs_stan_stk_int.pa_ther,
                            #coefs_stan_stk_int.pch_sal,
                            #coefs_stan_stk_int.Picorna2_virus,
                            coefs_stan_stk_int.prv,
                            coefs_stan_stk_int.pspv,
                            coefs_stan_stk_int.Qin,
                            coefs_stan_stk_int.rlo,
                            #coefs_stan_stk_int.re_sal,
                            coefs_stan_stk_int.Rhabdo3_virus,
                            coefs_stan_stk_int.sch,
                            #coefs_stan_stk_int.smallUK,
                            coefs_stan_stk_int.sp_des,
                            coefs_stan_stk_int.te_bry,
                            coefs_stan_stk_int.te_mar,
                            coefs_stan_stk_int.ven,
                            coefs_stan_stk_int.richness,
                            coefs_stan_stk_int.rib,
                            #coefs_stan_stk_slp.arena1,
                            coefs_stan_stk_slp.arena2,
                            #coefs_stan_stk_slp.ascv,
                            coefs_stan_stk_slp.c_b_cys,
                            coefs_stan_stk_slp.ce_sha,
                            coefs_stan_stk_slp.Circo.virus,
                            coefs_stan_stk_slp.de_sal,
                            coefs_stan_stk_slp.fa_mar,
                            coefs_stan_stk_slp.fl_psy,
                            coefs_stan_stk_slp.ic_hof,
                            coefs_stan_stk_slp.ic_mul,
                            coefs_stan_stk_slp.IcD,
                            coefs_stan_stk_slp.ihnv,
                            coefs_stan_stk_slp.ku_thy,
                            coefs_stan_stk_slp.lo_sal,
                            #coefs_stan_stk_slp.mo_vis,
                            coefs_stan_stk_slp.my_arc,
                            coefs_stan_stk_slp.nu_sal,
                            coefs_stan_stk_slp.pa_kab,
                            coefs_stan_stk_slp.pa_min,
                            coefs_stan_stk_slp.pa_pse,
                            coefs_stan_stk_slp.pa_ther,
                            #coefs_stan_stk_slp.pch_sal,
                            #coefs_stan_stk_slp.Picorna2_virus,
                            coefs_stan_stk_slp.prv,
                            coefs_stan_stk_slp.pspv,
                            coefs_stan_stk_slp.Qin,
                            coefs_stan_stk_slp.rlo,
                            #coefs_stan_stk_slp.re_sal,
                            coefs_stan_stk_slp.Rhabdo3_virus,
                            coefs_stan_stk_slp.sch,
                            #coefs_stan_stk_slp.smallUK,
                            coefs_stan_stk_slp.sp_des,
                            coefs_stan_stk_slp.te_bry,
                            coefs_stan_stk_slp.te_mar,
                            coefs_stan_stk_slp.ven,
                            coefs_stan_stk_slp.richness,
                            coefs_stan_stk_slp.rib,
                            #coefs_stan_stk_sig.arena1,
                            coefs_stan_stk_sig.arena2,
                            #coefs_stan_stk_sig.ascv,
                            coefs_stan_stk_sig.c_b_cys,
                            coefs_stan_stk_sig.ce_sha,
                            coefs_stan_stk_sig.Circo.virus,
                            coefs_stan_stk_sig.de_sal,
                            coefs_stan_stk_sig.fa_mar,
                            coefs_stan_stk_sig.fl_psy,
                            coefs_stan_stk_sig.ic_hof,
                            coefs_stan_stk_sig.ic_mul,
                            coefs_stan_stk_sig.IcD,
                            coefs_stan_stk_sig.ihnv,
                            coefs_stan_stk_sig.ku_thy,
                            coefs_stan_stk_sig.lo_sal,
                            #coefs_stan_stk_sig.mo_vis,
                            coefs_stan_stk_sig.my_arc,
                            coefs_stan_stk_sig.nu_sal,
                            coefs_stan_stk_sig.pa_kab,
                            coefs_stan_stk_sig.pa_min,
                            coefs_stan_stk_sig.pa_pse,
                            coefs_stan_stk_sig.pa_ther,
                            #coefs_stan_stk_sig.pch_sal,
                            #coefs_stan_stk_sig.Picorna2_virus,
                            coefs_stan_stk_sig.prv,
                            coefs_stan_stk_sig.pspv,
                            coefs_stan_stk_sig.Qin,
                            coefs_stan_stk_sig.rlo,
                            #coefs_stan_stk_sig.re_sal,
                            coefs_stan_stk_sig.Rhabdo3_virus,
                            coefs_stan_stk_sig.sch,
                            #coefs_stan_stk_sig.smallUK,
                            coefs_stan_stk_sig.sp_des,
                            coefs_stan_stk_sig.te_bry,
                            coefs_stan_stk_sig.te_mar,
                            coefs_stan_stk_sig.ven,
                            coefs_stan_stk_sig.richness,
                            coefs_stan_stk_sig.rib,
                            #coefs_stan_stk_year.arena1,
                            coefs_stan_stk_year.arena2,
                            #coefs_stan_stk_year.ascv,
                            coefs_stan_stk_year.c_b_cys,
                            coefs_stan_stk_year.ce_sha,
                            coefs_stan_stk_year.Circo.virus,
                            coefs_stan_stk_year.de_sal,
                            coefs_stan_stk_year.fa_mar,
                            coefs_stan_stk_year.fl_psy,
                            coefs_stan_stk_year.ic_hof,
                            coefs_stan_stk_year.ic_mul,
                            coefs_stan_stk_year.IcD,
                            coefs_stan_stk_year.ihnv,
                            coefs_stan_stk_year.ku_thy,
                            coefs_stan_stk_year.lo_sal,
                            #coefs_stan_stk_year.mo_vis,
                            coefs_stan_stk_year.my_arc,
                            coefs_stan_stk_year.nu_sal,
                            coefs_stan_stk_year.pa_kab,
                            coefs_stan_stk_year.pa_min,
                            coefs_stan_stk_year.pa_pse,
                            coefs_stan_stk_year.pa_ther,
                            #coefs_stan_stk_year.pch_sal,
                            #coefs_stan_stk_year.Picorna2_virus,
                            coefs_stan_stk_year.prv,
                            coefs_stan_stk_year.pspv,
                            coefs_stan_stk_year.Qin,
                            coefs_stan_stk_year.rlo,
                            #coefs_stan_stk_year.re_sal,
                            coefs_stan_stk_year.Rhabdo3_virus,
                            coefs_stan_stk_year.sch,
                            #coefs_stan_stk_year.smallUK,
                            coefs_stan_stk_year.sp_des,
                            coefs_stan_stk_year.te_bry,
                            coefs_stan_stk_year.te_mar,
                            coefs_stan_stk_year.ven,
                            coefs_stan_stk_year.richness,
                            coefs_stan_stk_year.rib)

write.csv(post_int_slpintsig, file="data/SW_Posterior distributions_Int Slp Sig_stspec_prev_230123.csv")

#### Rbind all convergence parameters and save as .csv file
post_rhatneff_prev <- rbind(#coefs_stan_stk_rhat.arena1,
                            coefs_stan_stk_rhat.arena2,
                            #coefs_stan_stk_rhat.ascv,
                            coefs_stan_stk_rhat.c_b_cys,
                            coefs_stan_stk_rhat.ce_sha,
                            coefs_stan_stk_rhat.Circo.virus,
                            coefs_stan_stk_rhat.de_sal,
                            coefs_stan_stk_rhat.fa_mar,
                            coefs_stan_stk_rhat.fl_psy,
                            coefs_stan_stk_rhat.ic_hof,
                            coefs_stan_stk_rhat.ic_mul,
                            coefs_stan_stk_rhat.IcD,
                            coefs_stan_stk_rhat.ihnv,
                            coefs_stan_stk_rhat.ku_thy,
                            coefs_stan_stk_rhat.lo_sal,
                            #coefs_stan_stk_rhat.mo_vis,
                            coefs_stan_stk_rhat.my_arc,
                            coefs_stan_stk_rhat.nu_sal,
                            coefs_stan_stk_rhat.pa_kab,
                            coefs_stan_stk_rhat.pa_min,
                            coefs_stan_stk_rhat.pa_pse,
                            coefs_stan_stk_rhat.pa_ther,
                            #coefs_stan_stk_rhat.pch_sal,
                            #coefs_stan_stk_rhat.Picorna2_virus,
                            coefs_stan_stk_rhat.prv,
                            coefs_stan_stk_rhat.pspv,
                            coefs_stan_stk_rhat.Qin,
                            coefs_stan_stk_rhat.rlo,
                            #coefs_stan_stk_rhat.re_sal,
                            coefs_stan_stk_rhat.Rhabdo3_virus,
                            coefs_stan_stk_rhat.sch,
                            #coefs_stan_stk_rhat.smallUK,
                            coefs_stan_stk_rhat.sp_des,
                            coefs_stan_stk_rhat.te_bry,
                            coefs_stan_stk_rhat.te_mar,
                            coefs_stan_stk_rhat.ven,
                            coefs_stan_stk_rhat.richness,
                            coefs_stan_stk_rhat.rib,
                            #coefs_stan_stk_neff.arena1,
                            coefs_stan_stk_neff.arena2,
                            #coefs_stan_stk_neff.ascv,
                            coefs_stan_stk_neff.c_b_cys,
                            coefs_stan_stk_neff.ce_sha,
                            coefs_stan_stk_neff.Circo.virus,
                            coefs_stan_stk_neff.de_sal,
                            coefs_stan_stk_neff.fa_mar,
                            coefs_stan_stk_neff.fl_psy,
                            coefs_stan_stk_neff.ic_hof,
                            coefs_stan_stk_neff.ic_mul,
                            coefs_stan_stk_neff.IcD,
                            coefs_stan_stk_neff.ihnv,
                            coefs_stan_stk_neff.ku_thy,
                            coefs_stan_stk_neff.lo_sal,
                            #coefs_stan_stk_neff.mo_vis,
                            coefs_stan_stk_neff.my_arc,
                            coefs_stan_stk_neff.nu_sal,
                            coefs_stan_stk_neff.pa_kab,
                            coefs_stan_stk_neff.pa_min,
                            coefs_stan_stk_neff.pa_pse,
                            coefs_stan_stk_neff.pa_ther,
                            #coefs_stan_stk_neff.pch_sal,
                            #coefs_stan_stk_neff.Picorna2_virus,
                            coefs_stan_stk_neff.prv,
                            coefs_stan_stk_neff.pspv,
                            coefs_stan_stk_neff.Qin,
                            coefs_stan_stk_neff.rlo,
                            #coefs_stan_stk_neff.re_sal,
                            coefs_stan_stk_neff.Rhabdo3_virus,
                            coefs_stan_stk_neff.sch,
                            #coefs_stan_stk_neff.smallUK,
                            coefs_stan_stk_neff.sp_des,
                            coefs_stan_stk_neff.te_bry,
                            coefs_stan_stk_neff.te_mar,
                            coefs_stan_stk_neff.ven,
                            coefs_stan_stk_neff.richness,
                            coefs_stan_stk_neff.rib)

write.csv(post_rhatneff_prev, file="data/SW_Posterior distributions_Rhat and Neff_stspec_prev_230123.csv")

### Plot posteriors per agent model from files
post_all <- post_int_slpintsig
post_all <- read.csv("data/SW_Posterior distributions_Int Slp Sig_stspec_prev_230123.csv")
post_agents <- read.csv("data/SW_prev_coefs_stan_stspec_230123.csv")
propzero <- read.csv("data/SW_Percent post draws >0_prev_stspec_230123.csv")
post_agents <- post_agents[order(match(post_agents[,1],propzero[,1])),]

## Plot Posteriors for all agents
jpeg(filename='figs/Fig_SW_prev_coefs_stan_stspec_blue_230123.jpg', 
      width=480, height=500, quality=75)
ggplot(post_agents) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(X, -mid), ymax = X75, ymin = X25), size=1.5, col="blue") +
  geom_linerange(aes(x = X, ymax = upper, ymin = lower), col="blue") +
  geom_point(aes(x = X, y = mid), size = 3, col="blue") +
  labs(x ="Agents", y = "Effect size") +
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
                            "de_sal" = "D. salmonis",
                            "richness" = "Richness",
                            "rib" = "RIB",
                            "ihnv" = "IHNV"))+
  theme(axis.text.y = element_text(face = "italic"), plot.title = element_text(hjust = 0.5))+
  coord_flip()
dev.off()

## Plots per agent
#### Extract output from agent model
post_c_b_cys <- post_all[post_all$X.1=="c_b_cys",]

## Extract Posterior slopes by Stock - c_b_cys
post_c_b_cys_stockslp <- post_c_b_cys[grep("prev_std Stock", post_c_b_cys$X) ,]

#### Plot
ggplot(post_c_b_cys_stockslp) +
  geom_hline(yintercept = 0, linetype = "dashed", col="blue")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), linewidth=1.5, col="black") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="black") +
  geom_point(aes(x = X, y = X50.), size = 3) +
  ylim(-0.8,0.8)+
  coord_flip()


## Extract Posterior intercepts for Stocks - c_b_cys example
post_c_b_cys_stockint <- post_c_b_cys[grep("b\\[\\(\\Intercept) Stock", post_c_b_cys$X) ,]
## Plot
ggplot(post_c_b_cys_stockint) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()

## Extract sigma values
post_c_b_cys_stocksig <- post_c_b_cys[grep("igma", post_c_b_cys$X) ,]
## Plot
ggplot(post_c_b_cys_stocksig) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = X, ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()

## Extract Posterior intercepts for years
post_c_b_cys_year <- post_c_b_cys[grep("b\\[\\(\\Intercept) Year", post_c_b_cys$X) ,]
## Plot
ggplot(post_c_b_cys_year) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_linerange(aes(x = reorder(X, -X50.), ymax = X75., ymin = X25.), size=1.5, col="gray") +
  geom_linerange(aes(x = X, ymax = X97.5., ymin = X2.5.), col="gray") +
  geom_point(aes(x = X, y = X50.), size = 2) +
  coord_flip()



## Calculate stock-specific Posterior estimates by adding 
## stock-specific draws to "global" (averaged across stocks) and then averaging

## Create a matrix and loop
temp2 <- data.frame(mod.c_b_cys)
temp3 <- temp2[,grepl("prev",names(temp2))] #include only columns with "prev" in name
temp4 <- temp3[,-grep("igma",colnames(temp3))] #remove columns with "igma" in name
temp5 <- temp4[2001:4000,] 
para_name2 <- colnames(temp5) #create an object with column names

stk.spec.slope <- matrix(NA,
                         nrow = 10,
                         ncol = 6,
                         dimnames = list(para_name2,c("2.5","25","50","75","97.5","agent")))

for(i in agents){
  model<-get(paste("mod",i, sep="."))
  temp2 <- data.frame(model)
  temp3 <- temp2[,grepl("prev",names(temp2))] #include only columns with "prev" in name
  temp4 <- temp3[,-grep("igma",colnames(temp3))] #remove columns with "igma" in name
  temp5 <- temp4[2001:4000,] #remove warm up iterations
  temp6 <- data.frame(temp5[,1] + temp5[,2:dim(temp5)[2]]) #add the stock-specific draws to global column
  temp7 <- cbind(temp5[,1],temp6) #bind with global column
  colnames(temp7) <- colnames(temp5) #assign names to columns
  nam <- paste("stk.spec.slope", i, sep = ".")
  temp8 <- as.matrix(apply(temp7, 2, quantile, probs = c(0.025,0.25,0.50,0.75,0.975)))
  temp9 <- t(temp8)
  temp10<- as.matrix(cbind(temp9, paste(i)))
  colnames(temp10) <- c("2.5","25","50","75","97.5","agent") #assign names to columns
  as.matrix(assign(nam, stk.spec.slope[i] <- temp10))
}

#### Rbind all posteriors and save as .csv file
stk.spec.slope.all <- rbind(#stk.spec.slope.arena1,
                            stk.spec.slope.arena2,
                            #stk.spec.slope.ascv,
                            stk.spec.slope.c_b_cys,
                            stk.spec.slope.ce_sha,
                            stk.spec.slope.Circo.virus,
                            stk.spec.slope.de_sal,
                            stk.spec.slope.fa_mar,
                            stk.spec.slope.fl_psy,
                            stk.spec.slope.ic_hof,
                            stk.spec.slope.ic_mul,
                            stk.spec.slope.IcD,
                            stk.spec.slope.ihnv,
                            stk.spec.slope.ku_thy,
                            stk.spec.slope.lo_sal,
                            #stk.spec.slope.mo_vis,
                            stk.spec.slope.my_arc,
                            stk.spec.slope.nu_sal,
                            stk.spec.slope.pa_kab,
                            stk.spec.slope.pa_min,
                            stk.spec.slope.pa_pse,
                            stk.spec.slope.pa_ther,
                            #stk.spec.slope.pch_sal,
                            #stk.spec.slope.Picorna2_virus,
                            stk.spec.slope.prv,
                            stk.spec.slope.pspv,
                            stk.spec.slope.Qin,
                            stk.spec.slope.rlo,
                            #stk.spec.slope.re_sal,
                            stk.spec.slope.Rhabdo3_virus,
                            stk.spec.slope.sch,
                            #stk.spec.slope.smallUK,
                            stk.spec.slope.sp_des,
                            stk.spec.slope.te_bry,
                            stk.spec.slope.te_mar,
                            stk.spec.slope.ven,
                            stk.spec.slope.richness,
                            stk.spec.slope.rib,
                            stk.spec.slope.arena2,
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
                            stk.spec.slope.ven,
                            stk.spec.slope.richness,
                            stk.spec.slope.rib)
write.csv(stk.spec.slope.all, file="data/SW_Stock specific slopes_prev_stspec_230123.csv")


##### READ IN DATA FROM FILE
stspslp <- read.csv("data/SW_Stock specific slopes_prev_stspec_230123.csv")
stspslp$stock <- substr(stspslp$X, 18, 28)
stspslp$stock <- substr(stspslp$stock, 1, nchar(stspslp$stock)-1)
stspslp$stock <- sub("^$", "Overall", stspslp$stock)



## Plot stock-specific slopes - c_b_cys
stk.spec.c_b_cys <-stspslp[stspslp$agent=="c_b_cys",]
## Plot
jpeg(filename='figs/SW_stockspslope_c_b_cys_230123.jpg', 
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
  labs(x="Stock", y="Effect size", title="Ca. B. cysticola") +
  coord_flip()
dev.off()

## Extract stock-specific slopes - ic_mul model
stk.spec.ic_mul <-stspslp[stspslp$agent=="ic_mul",]
## Plot
jpeg(filename='figs/SW_stockspslope_ic_mul_230123.jpg', 
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

## Plot stock-specific slopes - richness
stk.spec.richness <-stspslp[stspslp$agent=="richness",]
## Plot
jpeg(filename='figs/SW_stockspslope_230123_richness.jpg', 
     width=480, height=500, quality=75)
ggplot(stk.spec.richness) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.richness[stk.spec.richness$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.richness[stk.spec.richness$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.richness[stk.spec.richness$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size", title="Richness") +
  coord_flip()
dev.off()

## Plot stock-specific slopes - rib
stk.spec.rib <-stspslp[stspslp$agent=="rib",]
## Plot
jpeg(filename='figs/SW_stockspslope_230123_RIB.jpg', 
     width=480, height=500, quality=75)
ggplot(stk.spec.rib) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.richness[stk.spec.richness$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.richness[stk.spec.richness$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.richness[stk.spec.richness$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size", title="RIB") +
  coord_flip()
dev.off()

## Plot stock-specific slopes - te_mar
stk.spec.te_mar <-stspslp[stspslp$agent=="te_mar",]
## Plot
jpeg(filename='figs/SW_stockspslope_230123_te_mar.jpg', 
     width=480, height=500, quality=75)
ggplot(stk.spec.te_mar) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.te_mar[stk.spec.te_mar$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.te_mar[stk.spec.te_mar$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.te_mar[stk.spec.te_mar$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size", title="T. maritinum") +
  coord_flip()
dev.off()

## Plot stock-specific slopes - fl_psy
stk.spec.fl_psy <-stspslp[stspslp$agent=="fl_psy",]
## Plot
jpeg(filename='figs/SW_stockspslope_230123_fl_psy.jpg', 
     width=480, height=500, quality=75)
ggplot(stk.spec.fl_psy) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.fl_psy[stk.spec.fl_psy$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.fl_psy[stk.spec.fl_psy$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.fl_psy[stk.spec.fl_psy$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size", title="F. psychrophilum") +
  coord_flip()
dev.off()

## Plot stock-specific slopes - Qin
stk.spec.Qin <-stspslp[stspslp$agent=="Qin",]
## Plot
jpeg(filename='figs/SW_stockspslope_230123_Qin.jpg', 
     width=480, height=500, quality=75)
ggplot(stk.spec.Qin) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.Qin[stk.spec.Qin$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.Qin[stk.spec.Qin$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.Qin[stk.spec.Qin$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size", title="Qin virus") +
  coord_flip()
dev.off()

## Plot stock-specific slopes - arena2
stk.spec.arena2 <-stspslp[stspslp$agent=="arena2",]
## Plot
jpeg(filename='figs/SW_stockspslope_230123_arena2.jpg', 
     width=480, height=500, quality=75)
ggplot(stk.spec.arena2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.arena2[stk.spec.arena2$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.arena2[stk.spec.arena2$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.arena2[stk.spec.arena2$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size", title="SPAV-2") +
  coord_flip()
dev.off()

## Plot stock-specific slopes - Circo.virus
stk.spec.Circo.virus <-stspslp[stspslp$agent=="Circo.virus",]
## Plot
jpeg(filename='figs/SW_stockspslope_230123_Circovirus.jpg', 
     width=480, height=500, quality=75)
ggplot(stk.spec.Circo.virus) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(x = reorder(stock, -X50), ymax = X75, ymin = X25), size=1.5, col="gray") +
  geom_linerange(aes(x = stock, ymax = X97.5, ymin = X2.5), col="gray") +
  geom_linerange(data=stk.spec.Circo.virus[stk.spec.Circo.virus$stock=="Global",], 
                 aes(x = stock, ymax = X75, ymin = X25), size=2, col="black") +
  geom_linerange(data=stk.spec.Circo.virus[stk.spec.Circo.virus$stock=="Global",], 
                 aes(x = stock, ymax = X2.5, ymin = X97.5), col="black") +
  geom_point(aes(x = stock, y = X50), size = 2) +
  geom_point(data=stk.spec.Circo.virus[stk.spec.Circo.virus$stock=="Global",], aes(x = stock, y = X50), size = 3) +
  labs(x="Stock", y="Effect size", title="Circo virus") +
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
write.csv(param.prop0, file="data/SW_Percent post draws >0_prev_stspec_230123.csv")

# plot
propzero <- read.csv("data/SW_Percent post draws >0_prev_stspec_230123.csv")
colnames(propzero) <- c("agent","prev_prop_neg")
propzero <- merge(propzero, agent.names, by = c("agent"))

jpeg(filename='figs/Fig_prop<0 posterior_stspec_230123.jpg', 
     width=480, height=500, quality=75)
ggplot(propzero) +
  geom_bar(stat="identity", aes(reorder(plot.agent, prev_prop_neg), prev_prop_neg), col="blue", fill="blue", width=0.5, alpha=0.5) +
  labs(x="Agent", y="Proportion of total") +
  ggtitle("Sockeye salmon: \nProportion of correlation coefficients < 0")+
  theme(axis.text.y = element_text(face = "italic"))+
  coord_flip()
dev.off()


################
# Create plot of all points in MCMC
################
param2.plot <- matrix(NA,
                      ncol = 2000,
                      nrow = length(agents),
                      dimnames = list(agents,c(1:2000)))

for (i in agents){
  model<-as.matrix(get(paste("mod.",i, sep="")))
  model2<-model[2001:4000,1]
  temp <- as.matrix(model2)
  param2.plot[i,] <- temp
}
write.csv(param2.plot, file="data/SW_Individual_est_prev_stspec_230123.csv")



#######################################################################################
#######################################################################################
# STAN Approach for Multi-level Modeling
# FW PREVALENCE - INDEPENDENT MODELS by AGENT

## Stock-specific metric
### Create files for each agent

agents.fw <- unique(inf_agt_resid_data_fw$agent)
rr <- c("richness","rib")
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
prior_summary(mod.fw.richness)
plot(mod.fw.ic_mul, "ess")
plot(mod.fw.ic_mul, "trace")
pairs(mod.fw.richness)
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
write.csv(coefs_stan_fw, file="data/FW_prev_coefs_stan_stspec_230123.csv")

# Load estimates from file (if not running full model) and assign rownames
coefs_stan_fw <- read.csv("data/FW_prev_coefs_stan_stspec_230123.csv")
rownames(coefs_stan_fw) <- coefs_stan_fw[,1]
coefs_stan_fw <- coefs_stan_fw[,-1]  

# Plot effect size per agent
coefs_order <- coefs_stan_fw[order(-coefs_stan_fw[,3]),]

jpeg(filename='figs/FW_fig_prev_coefs_stan_stspec_230123.jpg', 
     width=480, height=500, quality=75)
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
       xlim = c(-2,1.5),
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
text(rep(-2,length(agents.fw)), 
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
  as.matrix(assign(nam, coefs_stan_stk_int.fw <- cbind(ind_coef[c(1:dim(ind_coef)[1]),c(4:8)], paste(i))))
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
  as.matrix(assign(nam, coefs_stan_stk_slp.fw <- cbind(ind_coef[c(1:dim(ind_coef)[1]),c(4:8)], paste(i))))
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
  as.matrix(assign(nam, coefs_stan_stk_year.fw <- cbind(ind_coef[c(1:dim(ind_coef)[1]),c(4:8)], paste(i))))
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
  as.matrix(assign(nam, coefs_stan_stk_sig.fw <- cbind(ind_coef[c(1:dim(ind_coef)[1]),c(4:8)], paste(i))))
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
  as.matrix(assign(nam, coefs_stan_stk_rhat.fw <- cbind(ind_coef[c(1:dim(ind_coef)[1]),1], paste(i), paste("Rhat"))))
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
  as.matrix(assign(nam, coefs_stan_stk_neff.fw <- cbind(ind_coef[c(1:dim(ind_coef)[1]),1], paste(i), paste("Neff"))))
}


## Derive proportions of posterior draws <0 per model - FW prevalence
param.fw<-colnames(data.frame(mod.fw.c_b_cys)) #create object of parameters in model
param.prop0.fw <- matrix(NA,
                      ncol = 1,
                      nrow = length(agents.fw),
                      dimnames = list(agents.fw,"prop0"))

for (i in agents.fw){
  model<-as.matrix(get(paste("mod.fw.",i, sep="")))
  model2<-model[2001:4000,1]
  param.prop0.fw[i,] <- (sum(model2 < 0))/2000
}
write.csv(param.prop0.fw, file="data/FW_Percent post draws >0_prev_stspec_230123.csv")
propzero.fw <- read.csv("data/FW_Percent post draws >0_prev_stspec_230123.csv")

jpeg(filename='figs/FW_Fig_prop<0_pa_ther_stspec_230123.jpg', 
     width=480, height=500, quality=75)
ggplot(propzero.fw) +
  geom_bar(stat="identity", aes(reorder(X, prop0), prop0), fill="gray", col="black") +
  ylim(0,1)+
  labs(x="Agent", y="Proportion of total", title="FW Posterior estimates <0") +
  coord_flip()
dev.off()



###############################################################################
## BRING ALL MODELING DATA TOGETHER

#beta estimates
prev_beta_sw <- read.csv("data/SW_prev_coefs_stan_stspec_230123.csv")
prev_beta_fw <- read.csv("data/FW_prev_coefs_stan_stspec_230123.csv")
#de<-data.frame("ku_thy",NA,NA,NA,NA,NA)
#colnames(de) <- c("X","lower","X25","mid","X75","upper")
#prev_beta_fw <- rbind(prev_beta_fw, de)

#proportion less than 0
temp <- read.csv("data/SW_Percent post draws >0_prev_stspec_230123.csv")
prev_neg_sw <- temp[,1:2]
temp <- read.csv("data/FW_Percent post draws >0_prev_stspec_230123.csv")
prev_neg_fw <- data.frame(temp[,1:2])
  
#sort all df to align with agent betas
prev_beta_sw <-prev_beta_sw[order(-prev_beta_sw$mid),]
prev_neg_sw <- prev_neg_sw[order(match(prev_neg_sw$X, prev_beta_sw$X)),]

prev_beta_fw <-prev_beta_fw[order(match(prev_beta_fw$X, prev_beta$X)),]
prev_neg_fw <- prev_neg_fw[order(match(prev_neg_fw$X, prev_beta_fw$X)),]

#Extract agent names into df - prevalence
agent.names <- data.frame(unique(inf_agt_resid_data[c("agent", "plot.agent")]))
agent.names$plot.agent <- as.character(agent.names$plot.agent)
agent.names$plot.agent[agent.names$plot.agent == "VENV"] <- "ENV"
agent.names$plot.agent[agent.names$plot.agent == "SPAV-2"] <- "SPAV-2"
agent.names$plot.agent[agent.names$plot.agent == "cr_sal"] <- "C. salmonicida"
agent.names$plot.agent[agent.names$plot.agent == "smallUK"] <- "Putative RNA virus"
agent.names$plot.agent[agent.names$plot.agent == "ihnv"] <- "IHNV"
agent.names$plot.agent[agent.names$plot.agent == "mo_vis"] <- "M. viscosa"
agent.names$plot.agent[agent.names$plot.agent == "Qin"] <- "Qin virus"
agent.names$plot.agent[agent.names$plot.agent == "pch_sal"] <- "P. salmonis"
agent.names$plot.agent[agent.names$plot.agent == "nu_sal"] <- "N. salmonis"
agent.names$plot.agent[agent.names$plot.agent == "arena1"] <- "SPAV-1"
agent.names$plot.agent[agent.names$plot.agent == "IcD"] <- "Ichthyobodo sp."
agent.names$plot.agent[agent.names$plot.agent == "Picorna2_virus"] <- "Picornavirus2"
agent.names$plot.agent[agent.names$plot.agent == "Rhabdo3_virus"] <- "Rhabdovirus3"
agent.names$plot.agent[agent.names$plot.agent == "ascv"] <- "ASCV"
agent.names$plot.agent[agent.names$plot.agent == "re_sal"] <- "R. salmoninarum"
agent.names$plot.agent[agent.names$plot.agent == "Circo.virus"] <- "Circovirus"
agent.names <- agent.names[order(match(agent.names$agent, prev_beta$X)),]

#Extract agent names into df - fw prev
#fw_names1 <- data.frame(prev_beta_fw[,1])
#colnames(fw_names1) <- "agent"
#fw_names2 <- merge(agent.names, fw_names1, by.x="agent", by.y="agent")
#fw.names <- fw_names2[order(match(fw_names2$agent, prev_beta_fw$X)),]





