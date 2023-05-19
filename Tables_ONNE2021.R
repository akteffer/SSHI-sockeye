## Tables for ONNE pathogens ~sockeye productivity/health manuscipt 2021
## A. K. Teffer
### Run code from Raw Data Filtering script for sw.data and fw.data

library(grid)
library(gridExtra)
library(ggplot2)
library(lattice)
library(dplyr)
library(tidyverse)

# Run code to extract data in figures script
unique(sw.data$Year)
unique(fw.data$Year)
all.data <- rbind(fw.data,sw.data)

#remove 2008 - not enough data and not used in modeling
all.data <- all.data[!all.data$Year==2008,]
fw.data <- fw.data[!fw.data$Year==2008,]
sw.data <- sw.data[!sw.data$Year==2008,]
unique(all.data$Year)
unique(fw.data$Year)
unique(sw.data$Year)

## Create tables
### Total samples per stock, year, SWFW
sw.samps.table <- table(sw.data$Stock_Analysis, sw.data$Year)
write.csv(sw.samps.table, file="figs/Table_SW_stockbyyear_230125.csv")
fw.samps.table <- table(fw.data$Stock_Analysis, fw.data$Year)
write.csv(fw.samps.table, file="figs/Table_FW_stockbyyear_230125.csv")
table(df$Year, df$SWFW)

## Modeling data - years with data (years stocks screened for each agent)
sw.mod.table.years <- table(inf_agt_resid_data$agent, inf_agt_resid_data$Stock_Analysis)
fw.mod.table.years <- table(inf_agt_resid_data_fw$agent, inf_agt_resid_data_fw$Stock_Analysis)

### Total assays per pathogen, year
assay.table = 
  data.frame(
  all.data %>%
      select(28:92) %>%  # replace to your needs
      summarise_each(funs(
        list(
          N = sum(!is.na(.)),
          perc.run = sum(!is.na(.))/length(df$Unique),
          posdet = sum(!is.na(.[.>0])),
          prev = sum(!is.na(.[.>0]))/sum(!is.na(.)),
          mean.load = mean(.[!is.na(.) & .!=0])
        )
      )
      )
  )
assay.table <- t(assay.table)
colnames(assay.table) <- c("N","perc.run","posdet","prev","mean.load")
write.csv(assay.table, file="figs/Assay_table_all_ONNE2021_221206.csv")

### SW
assay.table.sw = 
  data.frame(
    sw.mod.data %>%
      select(28:92) %>%  # replace to your needs
      summarise_each(funs(
        list(
          N = sum(!is.na(.)),
          perc.run = sum(!is.na(.))/length(sw.data$Unique),
          posdet = sum(!is.na(.[.>0])),
          prev = sum(!is.na(.[.>0]))/sum(!is.na(.)),
          mean.load = mean(.[!is.na(.) & .!=0])
        )
      )
      )
  )
assay.table.sw <- t(assay.table.sw)
colnames(assay.table.sw) <- c("N","perc.run","posdet","prev","mean.load")
write.csv(assay.table.sw, file="figs/Assay_table_SW_ONNE2021_230125.csv")

### FW
assay.table.fw = 
  data.frame(
    fw.mod.data %>%
      select(28:92) %>%  # replace to your needs
      summarise_each(funs(
        list(
          N = sum(!is.na(.)),
          perc.run = sum(!is.na(.))/length(fw.data$Unique),
          posdet = sum(!is.na(.[.>0])),
          prev = sum(!is.na(.[.>0]))/sum(!is.na(.)),
          mean.load = mean(.[!is.na(.) & .!=0])
        )
      )
      )
  )
assay.table.fw <- t(assay.table.fw)
colnames(assay.table.fw) <- c("N","perc.run","posdet","prev","mean.load")
write.csv(assay.table.fw, file="figs/Assay_table_FW_ONNE2021_230125.csv")


## Chilko SW
### SW
ch.mod.data <- sw.mod.data[sw.mod.data$Stock_Analysis=="Chilko",]
assay.table.ch = 
  data.frame(
    ch.mod.data %>%
      select(28:92) %>%  # replace to your needs
      summarise_each(funs(
        list(
          N = sum(!is.na(.)),
          perc.run = sum(!is.na(.))/length(sw.data$Unique),
          posdet = sum(!is.na(.[.>0])),
          prev = sum(!is.na(.[.>0]))/sum(!is.na(.)),
          mean.load = mean(.[!is.na(.) & .!=0])
        )
      )
      )
  )
assay.table.ch <- t(assay.table.ch)
colnames(assay.table.ch) <- c("N","perc.run","posdet","prev","mean.load")
write.csv(assay.table.ch, file="figs/Assay_table_Chilko_ONNE2021_230222.csv")


## total SW detections of each agent by year - note variable prevalence across agents and years
samplesperagent.sw <- inf_agt_resid_data %>% 
  group_by(agent, Year) %>%
  summarise(posdet = sum(posdet))

## total FW detections of each agent by year - note variable prevalence across agents and years
samplesperagent.fw<-inf_agt_resid_data_fw %>% 
  group_by(agent, Year) %>%
  summarise(posdet = sum(posdet)) 

## total SWCHILKO detections of each agent by year - note variable prevalence across agents and years
inf_agt_resid_data_ch <- inf_agt_resid_data[inf_agt_resid_data$Stock_Analysis=="Chilko",]
samplesperagent.ch <- inf_agt_resid_data_ch %>% 
  group_by(agent, Year) %>%
  summarise(posdet = sum(posdet))

## Richness and RIB ranges
names(df)
min(sw.data$richness, na.rm=T)
max(sw.data$richness, na.rm=T)
median(sw.data$richness, na.rm=T)
min(fw.data$richness, na.rm=T)
max(fw.data$richness, na.rm=T)
median(fw.data$richness, na.rm=T)

min(sw.data$rib, na.rm=T)
max(sw.data$rib, na.rm=T)
mean(sw.data$rib, na.rm=T)
min(fw.data$rib, na.rm=T)
max(fw.data$rib, na.rm=T)
mean(fw.data$rib, na.rm=T)

### Agent prevalence in data sets
prev.table <- as.data.frame(lapply(all.data[28:98], function(x) length(which(x>0)) / length(which(!is.na(x)))))
prev.table<-t(prev.table)
prev.table.sw <- as.data.frame(lapply(sw.data[28:98], function(x) length(which(x>0)) / length(which(!is.na(x)))))
prev.table.sw<-t(prev.table.sw)
prev.table.fw <- as.data.frame(lapply(fw.data[28:98], function(x) length(which(x>0)) / length(which(!is.na(x)))))
prev.table.fw<-t(prev.table.fw)
write.csv(prev.table, file="data/prev.table.all.csv")
write.csv(prev.table.sw, file="data/prev.table.sw.csv")
write.csv(prev.table.fw, file="data/prev.table.fw.csv")


## Agents we would not expect to see in FW
#  arena2, Pa_kab, sch, ven, pa_pse, ic_hof, ye_ruc
names(fw.data)
rare.agents.fw <- subset(fw.data, 
                        select = c(Unique,Date,Stock_Analysis,Latitude,Longitude,
                          arena2, pa_kab, sch, ven, pa_pse, ic_hof, ye_ruc))
head(rare.agents.fw)
write.csv(rare.agents.fw, file="data/weird.fw.agents.csv")

