
### SSHI sockeye exploratory analyses (stage 1)
#### by B. Connors
##### revised by A. Teffer 

#Quick exploratory analyses to derive Fraser sockeye survival indices. All code and associated data can be found on Github [here](https://github.com/brendanmichaelconnors/SSHI-sockeye). 
#Fraser brood table provided to A. Teffer by PSC, recent missing estimates of spawners in PSC brood table were added from file provided to A. Teffer by DFO (M. Hawkshaw).
#Updated data from Steve Latham at PSC: (productiondata(by_stock)VERYprel-20mailout(Oct02_20).xlsx). Updated return data from Tracy Cone at DFO:(Sockeye all from 1950.xlsx) -> all integrated into "fraser_brood_table_210519.txt" in data folder by AKT
#To (re)generate the SST data you need to run the "sst_anomalies.R" file that is commented out below and make sure the ersst package is installed (see [here](https://github.com/michaelmalick/r-ersst)).


## load required packages and functions
library(ersst)     
library(plyr)
library(tidyverse)
library(viridis)
source("functions.R")
#source("sst_anomalies.R") 

## load Fraser sockeye brood table
raw_brood_table <- read.delim("data/fraser_brood_table_210528.txt")
raw_brood_table$total_effective_female_spawner <- as.factor(raw_brood_table$total_effective_female_spawner)
str(raw_brood_table)

## clean up brood table, add a few new columns
brood_table <- subset(raw_brood_table, BY < 2017) 
brood_table <- subset(brood_table, BY > 1949) 
brood_table <- subset(brood_table, stock_name != "Cultus")  
brood_table$efs <-as.numeric(levels(brood_table$total_effective_female_spawner)
)[brood_table$total_effective_female_spawner]

brood_table$efs_lag1 <-lag(brood_table$efs,1)
brood_table$efs_lag2 <-lag(brood_table$efs,2)
brood_table$efs_lag3 <-lag(brood_table$efs,3)
brood_table$age.51 <- 0
# brood_table$lnRS <-log(brood_table$recruits_no_jacks/brood_table$efs)
brood_table$lnRS <-log(brood_table$total_recruits/brood_table$efs)
brood_table <- brood_table[!is.na(brood_table$lnRS),] 
brood_table <- brood_table[brood_table$lnRS != "-Inf",] 
brood_table<-brood_table[!is.na(brood_table$Stock.ID),]

## add columns with ocean entry age proportions
brood_table$ocean_0 <- (brood_table$age.21+
                          brood_table$age.31+
                          brood_table$age.41+
                          brood_table$age.51)/brood_table$total_recruits

brood_table$ocean_1 <- (brood_table$age.32+
                          brood_table$age.42+
                          brood_table$age.52+
                          brood_table$age.62)/brood_table$total_recruits

brood_table$ocean_2 <- (brood_table$age.43+
                          brood_table$age.53+
                          brood_table$age.63)/brood_table$total_recruits
head(brood_table)

## derive ocean entry year SST anomaly indices
##  weighted average across year and stock specific age at ocean entry
raw.clim <- read.csv(file="data/sst_yr_1_stock_anomalies.csv",header=TRUE)
early.sst <- clim.wgt.avg(brood.table = brood_table,
                          env.data = raw.clim,
                          env.covar = "sst_anomaly",
                          type = "first_year",
                          out.covar = "early_sst") 
head(early.sst)

## derive pink salmon competitor indices
##  empirical  through 2015, extrapolated based on last four years of data 
##  for odd and even lines thereafter
raw.comp <- read.csv(file="data/pink_abundance_2021_05_28.csv",header=TRUE)

np.pink <- pink.wgt.avg(brood.table = brood_table,
                        pink.data = raw.comp,
                        pink.covar = "Total",
                        type = "second_year",
                        out.covar = "np_pinks")
head(np.pink)

## merge datasets, normalized covariates, export
master.1 <- merge(brood_table, early.sst, by=c("BY","Stock.ID"),all.x=T)
master.bt <- merge(master.1, np.pink, by=c("BY","Stock.ID"),all.x=T) 
master.bt <- master.bt[order(master.bt$Stock.ID),]
master.bt_w_cov <- plyr::ddply(master.bt, .(Stock.ID), transform,
                               early_sst_stnd = scale(early_sst)[ , 1]
                               ,np_pinks_stnd = scale(np_pinks)[ , 1]
)  
head(master.bt_w_cov)
write.csv(master.bt_w_cov, "data/master_brood_table_covar_210528.csv", row.names=FALSE)

## derive sockeye survival indices
##  lnRS = natural log of recruits per spawner
##  SR_resid = residuals from ricker (or larkin) fit 
##  SR_cov_resid = residuals from ricker (or larkin) fit with covariates  

larkin <- c("Chilko", "L.Shuswap", "Quesnel", "Stellako", "Gates", " Pitt", 
            " Scotch", "Seymour", " Birkenhead")

survial_indices <- plyr::ddply(master.bt_w_cov, c("stock_name"),function(x) {
  SR_fit <- lm(x$lnRS~x$efs,na.action=na.exclude)
  SR_resid <- scale(resid(SR_fit))
  SR_fit_cov <- lm(x$lnRS~x$efs+
                     x$np_pinks_stnd+
                     x$early_sst_stnd,na.action=na.exclude)
  SR_cov_resid <- scale(resid(SR_fit_cov))
  
  if(unique(x$stock_name) %in% larkin){
    SR_fit <- lm(x$lnRS~x$efs+
                   x$efs_lag1+
                   x$efs_lag2,na.action=na.exclude)
    SR_resid <- scale(resid(SR_fit))
    SR_fit_cov <- lm(x$lnRS~x$efs+
                       x$efs_lag1+
                       x$efs_lag2+
                       x$efs_lag3+
                       x$np_pinks_stnd+
                       x$early_sst_stnd,na.action=na.exclude)
    SR_cov_resid <- scale(resid(SR_fit_cov))                            
  }
  
  lnRS <- scale(x$lnRS)
  brood_year <- x$BY
  data.frame(brood_year,SR_resid,SR_cov_resid,lnRS)})

survial_indicesL<-gather(survial_indices,survival_index,value,SR_resid:lnRS)
head(survial_indices)
write.csv(survial_indices, file="data/survival_indices_ONNE.csv")

## what are magnitude and direction of covariate "effects"?
cov_effects <- plyr::ddply(master.bt_w_cov, c("stock_name"),function(x) {
  SR_fit_cov <- lm(x$lnRS~x$efs+
                     x$np_pinks_stnd+
                     x$early_sst_stnd,na.action=na.exclude)
  pink <- round(SR_fit_cov$coefficients[3],digits=3)
  sst <- round(SR_fit_cov$coefficients[4],digits=3)
  if(unique(x$stock_name) %in% larkin){
    SR_fit_cov <- lm(x$lnRS~x$efs+
                       x$efs_lag1+
                       x$efs_lag2+
                       x$efs_lag3+
                       x$np_pinks_stnd+
                       x$early_sst_stnd,na.action=na.exclude)
    pink <- round(SR_fit_cov$coefficients[6],digits=3)
    sst <- round(SR_fit_cov$coefficients[7],digits=3)}
  data.frame(pink,sst)})

cov_effectsL<-gather(cov_effects,covariate,value,pink:sst)
print(cov_effects)

ggplot(cov_effectsL,aes(x=value, color = covariate, fill = covariate))+
  geom_density(alpha=0.5)+
  geom_rug()+
  scale_x_continuous(limits=c(-1,0.5))+
  scale_fill_viridis(discrete = T)+
  ylab("")+
  xlab("Std. effect size")+
  theme_bw()

#These are just point estimates of effects from independent linear regressions, but they give a general sense of the correlations between covariates and sockeye survival. THe standardized effect size is the change in sockeye productivity (in log space) per standard deviation unit increase in the covariate.

## plot full survival index time series
ggplot(survial_indicesL, aes(brood_year, value,colour = survival_index)) +
  geom_line()+
  facet_wrap(~stock_name,nrow=4)+
  scale_colour_viridis_d()+
  xlab("Brood year")+
  ylab("Index value")+
  theme_bw()

#Note that the survival indices have been standardized (mean = zero,  std. dev. = one).

## plot survival indices for recent years with SSHI pathogen data
ggplot(survial_indicesL[survial_indicesL$brood_year>2007,], aes(brood_year, value, colour = survival_index)) +
  geom_line()+
  facet_wrap(~stock_name,nrow=4)+
  scale_x_continuous(breaks = c(2009,2011,2013,2015))+ 
  scale_colour_viridis_d()+
  xlab("Brood year")+
  ylab("Index value")+
  theme_bw()

#These index values are what can be related to prevalence/load profiles of out-migrating juveniles in stage 2 of the analyses. Looks like some reasonable contrast among years and stocks. My advice for stage 2 is to use the stock-recruit residuals without covariates as the base survival index, and do a sensitivity analysis with the other two. 

write.csv(survial_indicesL[survial_indicesL$brood_year>2005,], file="data/survival_indices_truncated_210603.csv")

## plot survival time series 
plot.SR_cov_resid <- survial_indicesL[survial_indicesL$survival_index=="SR_cov_resid",]

sr_plot <- ggplot(plot.SR_cov_resid, aes(brood_year, value,colour = stock_name)) +
  geom_line(alpha=0.5)+
  geom_smooth(aes(brood_year, value), method = "gam", col="black") +
  xlab("Brood year")+
  ylab("Residual stock recruitment")+
  labs(col="Stock")+
  theme_bw()

jpeg(filename='figs/Fig_SR_plot_230203.jpg', width=700, height=355, quality=300)
sr_plot
dev.off()

