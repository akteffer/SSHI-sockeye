
### SSHI sockeye exploratory analyses (stage 1)
#### by B. Connors
##### revised by A. Teffer for FW Chilko analysis

#Quick exploratory analyses to derive Fraser sockeye survival indices. All code and associated data can be found on Github [here](https://github.com/brendanmichaelconnors/SSHI-sockeye). 
#Fraser brood table provided to A. Teffer by PSC, recent missing estimates of spawners in PSC brood table were added from file provided to A. Teffer by DFO (M. Hawkshaw).
#Updated data from Steve Latham at PSC: (productiondata(by_stock)VERYprel-20mailout(Oct02_20).xlsx). Updated return data from Tracy Cone at DFO:(Sockeye all from 1950.xlsx) -> all integrated into "fraser_brood_table_210519.txt" in data folder by AKT
#To (re)generate the SST data you need to run the "sst_anomalies.R" file that is commented out below and make sure the ersst package is installed (see [here](https://github.com/michaelmalick/r-ersst)).

## Structure of FW analysis:
###So stage one is model like: 
  #migrants/spawner ~ spawners + spawners_lag1 + spawners_lag2+ spawners_lag3. 
  #And then stage two is model like: resid_value ~ 0 +  prev_std + (1|Year).
  #Then, as a sensitivity analysis, you could repeat the stage two model 
    #with addition of a temperature covariate for years with data.
  
## load required packages and functions
#library(ersst)     
library(plyr)
library(tidyverse)
library(viridis)
source("functions.R")
library(dplyr)
library(ggplot2)
#source("sst_anomalies.R") 

## load Fraser sockeye brood table
raw_brood_table <- read.delim("data/fraser_brood_table_210528.txt")
raw_brood_table$total_effective_female_spawner <- 
  as.factor(raw_brood_table$total_effective_female_spawner)
str(raw_brood_table_ch)

## Only Chilko included
raw_brood_table_ch <- raw_brood_table[raw_brood_table$stock_name=="Chilko",]
dim(raw_brood_table_ch)
unique(raw_brood_table_ch$BY)

## Append smolt migration data (leaving Chilko lake at fence)
returnschilkofw <- read.csv("data/Chilko Smolts all years_FW.csv")
returnschilkofw$BY <- returnschilkofw$Smolt.migration.year-2 #assign brood year as 2 years prior to migration year
head(returnschilkofw)
raw_brood_table_ch <- merge(raw_brood_table_ch, returnschilkofw, by = "BY")
head(raw_brood_table_ch)
dim(raw_brood_table_ch)

# Calculate RSM -> "marine" survival (RSM = returning adults per smolt out-migrant)
raw_brood_table_ch$RSM <-raw_brood_table_ch$total_recruits/raw_brood_table_ch$Total.smolts
raw_brood_table_ch <- raw_brood_table_ch[!is.na(raw_brood_table_ch$RSM),] 
raw_brood_table_ch <- raw_brood_table_ch[raw_brood_table_ch$RSM != "-Inf",] 
head(raw_brood_table_ch)

## clean up brood table, add a few new columns for EFS lag
brood_table_ch <- subset(raw_brood_table_ch, BY < 2017) 
brood_table_ch <- subset(raw_brood_table_ch, BY > 1954) 
brood_table_ch$efs <-as.numeric(
  levels(brood_table_ch$total_effective_female_spawner))[
    brood_table_ch$total_effective_female_spawner]

brood_table_ch$efs_lag1 <-lag(brood_table_ch$efs,1)
brood_table_ch$efs_lag2 <-lag(brood_table_ch$efs,2)
brood_table_ch$efs_lag3 <-lag(brood_table_ch$efs,3)
brood_table_ch$age.51 <- 0

## Calculate lnMS -> migrants per spawner including density dep (lnMS = log out-migrants per effective female spawner)
brood_table_ch$lnMS <-log(brood_table_ch$Total.smolts/brood_table_ch$efs)
#brood_table_ch <- brood_table_ch[!is.na(brood_table_ch$lnMS),] 
#brood_table_ch <- brood_table_ch[brood_table_ch$lnMS != "-Inf",] 
#brood_table_ch$lnMS <- as.numeric(brood_table_ch$lnMS)
str(brood_table_ch)


# Replace missing 2015 lnMS & RSM values in 2015 using average of relationships
brood_table_ch$lnMS[brood_table_ch$Smolt.migration.year==2015] <- NA
brood_table_ch$Total.smolts[brood_table_ch$Smolt.migration.year==2015] <- NA
migpersp <- mean(brood_table_ch$Total.smolts/brood_table_ch$efs, na.rm=T)
brood_table_ch$Total.smolts[brood_table_ch$Smolt.migration.year==2015] <- 
  migpersp*(brood_table_ch$efs[brood_table_ch$Smolt.migration.year==2015])
brood_table_ch$lnMS[brood_table_ch$Smolt.migration.year==2015] <- 
  log(migpersp)

brood_table_ch$RSM[brood_table_ch$Smolt.migration.year==2015] <- NA
brood_table_ch$RSM[brood_table_ch$Smolt.migration.year==2015] <- 
  raw_brood_table_ch$total_recruits[brood_table_ch$Smolt.migration.year==2015]/
  raw_brood_table_ch$Total.smolts[brood_table_ch$Smolt.migration.year==2015]

## add columns with ocean entry age proportions
brood_table_ch$ocean_0 <- (brood_table_ch$age.21+
                          brood_table_ch$age.31+
                          brood_table_ch$age.41+
                          brood_table_ch$age.51)/brood_table_ch$Total.smolts

brood_table_ch$ocean_1 <- (brood_table_ch$age.32+
                          brood_table_ch$age.42+
                          brood_table_ch$age.52+
                          brood_table_ch$age.62)/brood_table_ch$Total.smolts

brood_table_ch$ocean_2 <- (brood_table_ch$age.43+
                          brood_table_ch$age.53+
                          brood_table_ch$age.63)/brood_table_ch$Total.smolts
head(brood_table_ch)


## Export csv file of Chilko brood table
brood_table_ch <- brood_table_ch[brood_table_ch$BY>1960,]
write.csv(brood_table_ch, "data/master_brood_table_FWchilko_220901.csv", row.names=FALSE)



##  Derive sockeye survival indices
##  lnMS = natural log of FW migrants per spawner
##  RSM = returning adults per out-migrating smolt
##  SR_resid_chfw = residuals from larkin fit chilko in FW

larkin <- "Chilko"
survival_indices_chfw <- plyr::ddply(brood_table_ch, c("stock_name"), function(x) {
    SR_fit_chfw <- lm(x$lnMS ~ x$efs+ #accounts for density dependence and EFS
                       x$efs_lag1+
                       x$efs_lag2+
                       x$efs_lag3)
    SR_resid_chfw <- scale(resid(SR_fit_chfw))

    lnMS <- scale(x$lnMS)
    RSM <- scale(x$RSM)
    brood_year <- x$BY
    data.frame(brood_year,SR_resid_chfw,lnMS,RSM)})

survival_indices_chfwL<-gather(survival_indices_chfw,survival_index,value,SR_resid_chfw:RSM)
head(survival_indices_chfwL)
write.csv(survival_indices_chfwL, file="data/survival_indices_ONNE_chilkoFW.csv")


## plot full survival index time series
ggplot(survival_indices_chfwL, aes(brood_year, value, colour = survival_index)) +
  geom_line()+
  facet_wrap(~stock_name,nrow=4)+
  scale_colour_viridis_d()+
  xlab("Brood year")+
  ylab("Index value")+
  theme_bw()

#Note that the survival indices have been standardized (mean = zero,  std. dev. = one).
ggplot(survival_indices_chfwL[survival_indices_chfwL$brood_year>2008,], aes(brood_year, value, colour = survival_index)) +
  geom_line()+
  facet_wrap(~stock_name,nrow=4)+
  scale_x_continuous(breaks = c(2007,2009,2011,2013,2015))+ 
  scale_colour_viridis_d()+
  xlab("Brood year")+
  ylab("Index value")+
  theme_bw()


