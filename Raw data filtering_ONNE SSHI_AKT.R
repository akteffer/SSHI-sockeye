## ONNE metadata filtering code 
## 210519 - Amy K Teffer

## Set your working directory to where you will source the files from (change the file path below)
## I keep my csv data files in a "data" folder within the project folder, so the code can call from there to load it
library(kableExtra)
library(grid)
library(gridExtra)
library(ggplot2)
library(lattice)
library(dplyr)

#Bring in truncated resdiduals
trnc_resid<-read.csv("data/survival_indices_truncated_210603.csv", head=TRUE)
str(trnc_resid)
trnc_resid$Year <- trnc_resid$brood_year+2
names(trnc_resid) <- c("orderID", "Stock_Analysis", "brood_year", "metric", "resid_value", "Year") #rename columns

##create object with just SRR_resid metric (including covariate effects) to align with infection data
trnc_resid_srr <- trnc_resid[trnc_resid$metric=="SR_cov_resid",]
unique(trnc_resid$Stock_Analysis)

#Bring in pathogen data without LOD and remove samples with no stock ID (may be updated data available!)
all <- read.csv("data/ONNE metadata no LOD_6.17.2021.csv",header=TRUE)
dim(all)
min(all[,c(27:97)], na.rm=T)
# replace negative values with 0
#all[all==-999]<-0
all[,c(27:97)][all[,c(27:97)]<0] <-0

## Check against fitchip data 210720
fitchip <- read.csv("data/fraser_fitchip.csv",header=TRUE)
min(fitchip[,c(112:154)], na.rm=T)
fitchip[,c(112:154)][fitchip[,c(112:154)]<0] <-0
fitchip2 <- fitchip[,c(5,112:154)]
names(fitchip2)
colnames(fitchip2) <- c("Unique",                                            
                        "c_b_cys","sch","fl_psy","mo_vis","pch_sal","pisck_sal","re_sal","rlo","te_mar",
                        "hab_he_aka","hab_pseudo","ce_sha","de_sal","fa_mar","IcD","ic_hof","ic_mul",
                        "ku_thy","lo_sal","my_arc","my_ins", "na_sal", "nu_sal","pa_ther","pa_kab","pa_min",                      
                        "pa_pse", "sp_des","te_bry", "Circo.virus","ctv", "ven","ihnv","pspv","Picorna2_virus",                                     
                        "prv", "smallUK",  "Qin", "Rhabdo3_virus","arena1","arena2","ver","vhsv")

temp1 <- merge(all, fitchip2, by = "Unique") #merge 2 df's based on ID
temp2 <- temp1[,-grep("\\.x",colnames(temp1))] 
substrLeft <- function(x, n){
  substr(x, 0, nchar(x)-n)}
n_col <- ncol(temp2) #how many agents?
names(temp2)[57:n_col] <- substrLeft(names(temp2[57:n_col]),2) #trim names of columns
temp3 <- temp2[,-c(55,66:67)] #remove HABs and ye_ruc (no detections)
colnames(temp3) [55] <- "ye_ruc" #rename ye_ruc_glnA
temp3 = temp3[!duplicated(temp3$Unique),] #remove duplicates

nochange.ids <- unique(temp2$Unique)
nochange.df <- all[!all$Unique %in% nochange.ids,]
nochange.df2 <- nochange.df[-96]
colnames(nochange.df2) [96] <- "ye_ruc" #rename ye_ruc_glnA
temp.match <- intersect(colnames(nochange.df), colnames(temp3))
all1 <- rbind(nochange.df2, temp3)
dim(all1)
all1 = all1[!duplicated(all1$Unique),] #remove duplicates
write.csv(all1, file="data/match_fitchip_with_all.csv")

#take out fish with no stock assignment or from rare stock regions (northern, QCI, Transboundary)
all2 <- droplevels(all1[!(all1$Stock_Region==0) ,]) #remove any fish without stock assignment 
temp2<-droplevels(all2[-which(all2$Stock_Area=="Northern") ,]) #Northern removed
dim(temp2)
temp3<-droplevels(temp2[-which(temp2$Stock_Area=="QCI") ,]) #QCI removed
dim(temp3)
major<-droplevels(temp3[-which(temp3$Stock_Area=="TransBoundary") ,]) #Transboundary removed - major stocks left
dim(major)

## Create separate objects for SW and FW collected samples
sw.major<-droplevels(major[(major$SWFW=="SW") ,])#SW only
head(sw.major)
fw.major<-droplevels(major[(major$SWFW=="FW") ,])#FW only 
head(fw.major)

## Temporal limits: Reduce sampling period to spring-summer only
spsu1<-droplevels(sw.major[-which(sw.major$SEASON1=="Overwinter") ,]) #remove overwinter 
dim(spsu1)
spsu2<-droplevels(spsu1[-which(spsu1$SEASON1=="Winter") ,]) #remove winter 
dim(spsu2)
spsu3<-droplevels(spsu2[-which(spsu2$SEASON1=="Fall") ,]) # remove fall 
dim(spsu3) 

## Geographically limit: Remove samples from WCVI, 2018 and high latitudes (>51.5 lat) 
spsu4<-droplevels(spsu3[-which(spsu3$Zone=="WCVI") ,]) #remove WCVI 
dim(spsu4)
spsu5<-droplevels(spsu4[-which(spsu4$Latitude > 51.5) ,]) # remove high latitude samples
dim(spsu5)
spsu <- subset(spsu5, !is.na(DOY)) #remove rows with missing Date
dim(spsu)

## Change stock names to align with SR data in both SW and FW data
#SW
levels(spsu$Stock_Analysis)[levels(spsu$Stock_Analysis)=="Early Stuart"] <- "E.Stuart"
levels(spsu$Stock_Analysis)[levels(spsu$Stock_Analysis)=="Late Shuswap"] <- "L.Shuswap"
levels(spsu$Stock_Analysis)[levels(spsu$Stock_Analysis)=="Late Stuart"] <- "L.Stuart"
levels(spsu$Stock_Analysis)[levels(spsu$Stock_Analysis)=="Harrison-Widgeon"] <- "Harrison"
plot(spsu$DOY)
sw.data <- spsu #rename

#FW
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Early Stuart"] <- "E.Stuart"
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Late Shuswap"] <- "L.Shuswap"
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Late Stuart"] <- "L.Stuart"
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Harrison-Widgeon"] <- "Harrison"
#Remove samples with a DOY>200 (out of range of samples collected in FW)
fw.data <- subset(fw.major, DOY<200)
dim(fw.data)

## Calculate richness
# SW
rich_sw =
  data.frame(
    sw.data %>% 
      select(28:96) %>%
      mutate(richness = rowSums(.>0,na.rm=TRUE))) 

# FW
rich_fw =
  data.frame(
    fw.data %>% 
      select(28:96) %>%
      mutate(richness = rowSums(.>0,na.rm=TRUE))) 

# Calculate RIB
#SW
ribdata_sw <- as.matrix(sw.data[,28:96])
rib1_sw <- data.frame(apply(ribdata_sw,2,function(x) {x/max(x, na.rm=TRUE)}))
rib2_sw <- data.frame(rowSums(rib1_sw, na.rm=TRUE))
richrib_sw <- cbind(rich_sw,rib2_sw)
colnames(richrib_sw)[71] <- "rib"
richrib_sw$rib[which(richrib_sw$rib==-Inf)] <- 0
sw.data <- cbind(sw.data, richrib_sw[,70:71])
# FW
ribdata_fw <- as.matrix(fw.data[,28:96])
rib1_fw <- data.frame(apply(ribdata_fw,2,function(x){x/max(x, na.rm=TRUE)}))
rib2_fw <- data.frame(rowSums(rib1_fw, na.rm=TRUE))
richrib_fw <- cbind(rich_fw,rib2_fw)
colnames(richrib_fw)[71] <- "rib"
richrib_fw$rib[which(richrib_fw$rib==-Inf)] <- 0
fw.data <- cbind(fw.data, richrib_fw[,70:71])

## Examine output
head(sw.data)
dim(sw.data)
head(fw.data)
dim(fw.data)
data.frame(table(sw.data$Stock_Analysis))
data.frame(table(fw.data$Stock_Analysis))

## Merge into one df
all.data <- rbind(fw.data, sw.data)

#write csvs
write.csv(fw.data, file="data/fw.data.csv")
write.csv(sw.data, file="data/sw.data.csv")
write.csv(all.data, file="data/all.data.csv")

## Merge with SR data
resid.sw <- merge(trnc_resid_srr, sw.data, by = c("Stock_Analysis","Year"))
resid.fw <- merge(trnc_resid_srr, fw.data, by = c("Stock_Analysis","Year"))
resid.all <- rbind(resid.sw,resid.fw)

### Data Exploration ###
## Total agent detections all agents
agents <- names(resid.all[,32:100])
rr.agents <- names(resid.all[32:102])
agents.data <- resid.all[,32:102]
nRows <- dim(agents.data)[1]
calcStats <- function(x) {
  pos <- sum(agents.data[,x] > 0, na.rm=TRUE)
  c("positives" = pos, "proportion" = pos / nRows)
}
result <- t(as.data.frame(Map(calcStats, colnames(agents.data))))
result[order(result[,1]),]
write.csv(result[order(result[,1]),], file="figs/Final data frames and tables/ONNE_all_agent_detections_220126.csv")

## SW only
agents.sw <- names(resid.sw[,32:100])
rr.agents.sw <- names(resid.sw[32:102])
agents.data.sw <- resid.sw[,32:102]
nRows <- dim(agents.data.sw)[1]
calcStats <- function(x) {
  pos <- sum(agents.data.sw[,x] > 0, na.rm=TRUE)
  c("positives" = pos, "proportion" = pos / nRows)
}
result.sw <- t(as.data.frame(Map(calcStats, colnames(agents.data.sw))))
result.sw[order(result.sw[,1]),]
write.csv(result.sw[order(result.sw[,1]),], file="figs/Final data frames and tables/ONNE_all_agent_detections_220126_SW.csv")

## FW only
agents.fw <- names(resid.fw[,32:100])
rr.agents.fw <- names(resid.fw[32:102])
agents.data.fw <- resid.fw[,32:102]
nRows <- dim(agents.data.fw)[1]
calcStats <- function(x) {
  pos <- sum(agents.data.fw[,x] > 0, na.rm=TRUE)
  c("positives" = pos, "proportion" = pos / nRows)
}
result.fw <- t(as.data.frame(Map(calcStats, colnames(agents.data.fw))))
result.fw[order(result.fw[,1]),]
write.csv(result.fw[order(result.fw[,1]),], file="figs/Final data frames and tables/ONNE_all_agent_detections_220126_FW.csv")

N.resid.sw <- data.frame(resid.sw %>% 
  group_by(Year, Stock_Analysis, SWFW) %>%
  count(SWFW))
dim(N.resid.sw)
N.resid.fw <- data.frame(resid.fw %>% 
  group_by(Year, Stock_Analysis, SWFW) %>%
  count(SWFW))
dim(N.resid.fw)

resid.N.all <- merge(N.resid.fw, N.resid.sw, by = c("Year","Stock_Analysis"), all.y=TRUE, all.x=TRUE)
resid.N.all <- resid.N.all[,c(1,2,4,6)]
colnames(resid.N.all) <- c("Year","Stock","FW.N","SW.N")

resid.N.all.wide <- reshape(resid.N.all, idvar = "Stock", 
                            timevar = "Year", direction = "wide")
rownames (resid.N.all.wide) <- NULL
colnames(resid.N.all.wide) <- c("Stock","FW","SW","FW","SW",
                                "FW","SW","FW","SW",
                                "FW","SW","FW","SW",
                                "FW","SW","FW","SW",
                                "FW","SW","FW","SW")
write.csv(resid.N.all.wide, file="figs/Final data frames and tables/ONNE_agent_detectionsSWFW_Stock_SURV_220126.csv")

## ADD COLUMN TOTALS
resid.sw <- merge(trnc_resid_srr, sw.data, by = c("Stock_Analysis","Year"))
N.resid.sw <- data.frame(table(resid.sw$Stock_Analysis))
sw.N.st <- ggplot(data=N.resid.sw, aes(x=reorder(Var1, Freq), y=Freq)) +
  geom_col() +
  coord_flip() +
  xlab("Stocks")+
  ylab("Samples")+
  ggtitle("Marine-collected fish")+
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())
resid.fw <- merge(trnc_resid_srr, fw.data, by = c("Stock_Analysis","Year"))
N.resid.fw <- data.frame(table(resid.fw$Stock_Analysis))
fw.N.st <- ggplot(data=N.resid.fw, aes(x=reorder(Var1, Freq), y=Freq)) +
  geom_col() +
  coord_flip() +
  xlab("Stocks")+
  ylab("Samples") +
  ggtitle("Freshwater-collected samples")+
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())
N.resid.fw$source <- "Freshwater"
N.resid.sw$source <- "Marine"
N.resid.all <- rbind(N.resid.fw, N.resid.sw)


jpeg(filename='figs/Fig_Total samples by stock_SWandFW_220126.jpg', 
     width=480, height=400, quality=300)
grid.arrange(arrangeGrob(fw.N.st + theme(legend.position="none"), 
            sw.N.st + theme(legend.position="none"), 
            ncol = 2,
            left = textGrob("Stock", rot = 90, vjust = 1)),
            bottom = textGrob("Total samples", vjust = -0.5))
dev.off()

all.resid.data <- merge(trnc_resid_srr, all.data, by = c("Stock_Analysis","Year"))
samplesstyr.fw <- resid.fw %>% 
  group_by(Year, Stock_Analysis, SWFW) %>%
  count(SWFW)
samplesstyr.all<-all.resid.data %>% 
  group_by(Year, Stock_Analysis, SWFW) %>%
  count(Stock_Analysis)

## How many PRV detections?
dim(all1)
##arena1
prv.detects =
  data.frame(
    all1 %>% 
      group_by(Stock_Analysis, Year) %>%
        summarise(
          length(which(prv!="NA")), #samples
          length(which(prv>0)), #positive detections
          length(which(prv>0))/length(which(!is.na(prv))),  #calculates prevalence
          mean(prv[prv!=0], na.rm=TRUE),
          "prv"
      )
  )
names(prv.detects) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "agent") #rename columns
write.csv(prv.detects, file="data/prv detections sockeye dataset_AKT.csv")


###########################################################################################
### Create df for survival analysis ####
###########################################################################################
# NEED >1 DETECTION TO BE INCLUDED
result.sw[order(result.sw[,1]),]

## SW data by stock and sampling year
spsu.stock <-sw.data %>% group_by(Stock_Analysis, Year) #create object to be summarized by year and stock

##ASSESS INFECTION PARAMETERS BY AGENT:
##arena1
stock.arena1.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(arena1!="NA")), #samples
        length(which(arena1>0)), #positive detections
        length(which(arena1>0))/length(which(!is.na(arena1))),  #calculates prevalence
        mean(arena1[arena1!=0], na.rm=TRUE),
        (length(which(arena1>0)) / length(which(!is.na(arena1)))) * mean(arena1[arena1!=0], na.rm=TRUE),
        "arena1"
      )
  )
names(stock.arena1.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##arena2
stock.arena2.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(arena2!="NA")), #samples
        length(which(arena2>0)), #positive detections
        length(which(arena2>0))/length(which(!is.na(arena2))),  #calculates prevalence
        mean(arena2[arena2!=0], na.rm=TRUE),
        (length(which(arena2>0)) / length(which(!is.na(arena2)))) * mean(arena2[arena2!=0], na.rm=TRUE),
        "arena2"
      )
  )
names(stock.arena2.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##ascv
stock.ascv.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(ascv!="NA")), #samples
        length(which(ascv>0)), #positive detections
        length(which(ascv>0))/length(which(!is.na(ascv))),  #calculates prevalence
        mean(ascv[ascv!=0], na.rm=TRUE),
        (length(which(ascv>0)) / length(which(!is.na(ascv)))) * mean(ascv[ascv!=0], na.rm=TRUE),
        "ascv"
      )
  )
names(stock.ascv.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##c_b_cys
stock.c_b_cys.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(c_b_cys!="NA")), #samples
        length(which(c_b_cys>0)), #positive detections
        length(which(c_b_cys>0))/length(which(!is.na(c_b_cys))),  #calculates prevalence
        mean(c_b_cys[c_b_cys!=0], na.rm=TRUE),
        (length(which(c_b_cys>0)) / length(which(!is.na(c_b_cys)))) * mean(c_b_cys[c_b_cys!=0], na.rm=TRUE),
        "c_b_cys"
      )
  )
names(stock.c_b_cys.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: ce_sha
stock.ce_sha.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(ce_sha!="NA")), #samples
        length(which(ce_sha>0)), #positive detections
        length(which(ce_sha>0))/length(which(!is.na(ce_sha))),  #calculates prevalence
        mean(ce_sha[ce_sha!=0], na.rm=TRUE),
        (length(which(ce_sha>0)) / length(which(!is.na(ce_sha)))) * mean(ce_sha[ce_sha!=0], na.rm=TRUE),
        "ce_sha"
      )
  )
names(stock.ce_sha.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: Circo.virus
stock.Circo.virus.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(Circo.virus!="NA")), #samples
        length(which(Circo.virus>0)), #positive detections
        length(which(Circo.virus>0))/length(which(!is.na(Circo.virus))),  #calculates prevalence
        mean(Circo.virus[Circo.virus!=0], na.rm=TRUE),
        (length(which(Circo.virus>0)) / length(which(!is.na(Circo.virus)))) * mean(Circo.virus[Circo.virus!=0], na.rm=TRUE),
        "Circo.virus"
      )
  )
names(stock.Circo.virus.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: cr_sal
stock.cr_sal.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(cr_sal!="NA")), #samples
        length(which(cr_sal>0)), #positive detections
        length(which(cr_sal>0))/length(which(!is.na(cr_sal))),  #calculates prevalence
        mean(cr_sal[cr_sal!=0], na.rm=TRUE),
        (length(which(cr_sal>0)) / length(which(!is.na(cr_sal)))) * mean(cr_sal[cr_sal!=0], na.rm=TRUE),
        "cr_sal"
      )
  )
names(stock.cr_sal.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: de_sal
stock.de_sal.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(de_sal!="NA")), #samples
        length(which(de_sal>0)), #positive detections
        length(which(de_sal>0))/length(which(!is.na(de_sal))),  #calculates prevalence
        mean(de_sal[de_sal!=0], na.rm=TRUE),
        (length(which(de_sal>0)) / length(which(!is.na(de_sal)))) * mean(de_sal[de_sal!=0], na.rm=TRUE),
        "de_sal"
      )
  )
names(stock.de_sal.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: fa_mar
stock.fa_mar.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(fa_mar!="NA")), #samples
        length(which(fa_mar>0)), #positive detections
        length(which(fa_mar>0))/length(which(!is.na(fa_mar))),  #calculates prevalence
        mean(fa_mar[fa_mar!=0], na.rm=TRUE),
        (length(which(fa_mar>0)) / length(which(!is.na(fa_mar)))) * mean(fa_mar[fa_mar!=0], na.rm=TRUE),
        "fa_mar"
      )
  )
names(stock.fa_mar.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: fl_psy
stock.fl_psy.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(fl_psy!="NA")), #samples
        length(which(fl_psy>0)), #positive detections
        length(which(fl_psy>0))/length(which(!is.na(fl_psy))),  #calculates prevalence
        mean(fl_psy[fl_psy!=0], na.rm=TRUE),
        (length(which(fl_psy>0)) / length(which(!is.na(fl_psy)))) * mean(fl_psy[fl_psy!=0], na.rm=TRUE),
        "fl_psy"
      )
  )
names(stock.fl_psy.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: ic_hof
stock.ic_hof.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(ic_hof!="NA")), #samples
        length(which(ic_hof>0)), #positive detections
        length(which(ic_hof>0))/length(which(!is.na(ic_hof))),  #calculates prevalence
        mean(ic_hof[ic_hof!=0], na.rm=TRUE),
        (length(which(ic_hof>0)) / length(which(!is.na(ic_hof)))) * mean(ic_hof[ic_hof!=0], na.rm=TRUE),
        "ic_hof"
      )
  )
names(stock.ic_hof.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: ic_mul
stock.ic_mul.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(ic_mul!="NA")), #samples
        length(which(ic_mul>0)), #positive detections
        length(which(ic_mul>0))/length(which(!is.na(ic_mul))),  #calculates prevalence
        mean(ic_mul[ic_mul!=0], na.rm=TRUE),
        (length(which(ic_mul>0)) / length(which(!is.na(ic_mul)))) * mean(ic_mul[ic_mul!=0], na.rm=TRUE),
        "ic_mul"
      )
  )
names(stock.ic_mul.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: IcD
stock.IcD.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(IcD!="NA")), #samples
        length(which(IcD>0)), #positive detections
        length(which(IcD>0))/length(which(!is.na(IcD))),  #calculates prevalence
        mean(IcD[IcD!=0], na.rm=TRUE),
        (length(which(IcD>0)) / length(which(!is.na(IcD)))) * mean(IcD[IcD!=0], na.rm=TRUE),
        "IcD"
      )
  )
names(stock.IcD.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: ihnv
stock.ihnv.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(ihnv!="NA")), #samples
        length(which(ihnv>0)), #positive detections
        length(which(ihnv>0))/length(which(!is.na(ihnv))),  #calculates prevalence
        mean(ihnv[ihnv!=0], na.rm=TRUE),
        (length(which(ihnv>0)) / length(which(!is.na(ihnv)))) * mean(ihnv[ihnv!=0], na.rm=TRUE),
        "ihnv"
      )
  )
names(stock.ihnv.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: ku_thy
stock.ku_thy.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(ku_thy!="NA")), #samples
        length(which(ku_thy>0)), #positive detections
        length(which(ku_thy>0))/length(which(!is.na(ku_thy))),  #calculates prevalence
        mean(ku_thy[ku_thy!=0], na.rm=TRUE),
        (length(which(ku_thy>0)) / length(which(!is.na(ku_thy)))) * mean(ku_thy[ku_thy!=0], na.rm=TRUE),
        "ku_thy"
      )
  )
names(stock.ku_thy.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: lo_sal
stock.lo_sal.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(lo_sal!="NA")), #samples
        length(which(lo_sal>0)), #positive detections
        length(which(lo_sal>0))/length(which(!is.na(lo_sal))),  #calculates prevalence
        mean(lo_sal[lo_sal!=0], na.rm=TRUE),
        (length(which(lo_sal>0)) / length(which(!is.na(lo_sal)))) * mean(lo_sal[lo_sal!=0], na.rm=TRUE),
        "lo_sal"
      )
  )
names(stock.lo_sal.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: mo_vis
stock.mo_vis.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(mo_vis!="NA")), #samples
        length(which(mo_vis>0)), #positive detections
        length(which(mo_vis>0))/length(which(!is.na(mo_vis))),  #calculates prevalence
        mean(mo_vis[mo_vis!=0], na.rm=TRUE),
        (length(which(mo_vis>0)) / length(which(!is.na(mo_vis)))) * mean(mo_vis[mo_vis!=0], na.rm=TRUE),
        "mo_vis"
      )
  )
names(stock.mo_vis.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: my_arc
stock.my_arc.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(my_arc!="NA")), #samples
        length(which(my_arc>0)), #positive detections
        length(which(my_arc>0))/length(which(!is.na(my_arc))),  #calculates prevalence
        mean(my_arc[my_arc!=0], na.rm=TRUE),
        (length(which(my_arc>0)) / length(which(!is.na(my_arc)))) * mean(my_arc[my_arc!=0], na.rm=TRUE),
        "my_arc"
      )
  )
names(stock.my_arc.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: nu_sal
stock.nu_sal.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(nu_sal!="NA")), #samples
        length(which(nu_sal>0)), #positive detections
        length(which(nu_sal>0))/length(which(!is.na(nu_sal))),  #calculates prevalence
        mean(nu_sal[nu_sal!=0], na.rm=TRUE),
        (length(which(nu_sal>0)) / length(which(!is.na(nu_sal)))) * mean(nu_sal[nu_sal!=0], na.rm=TRUE),
        "nu_sal"
      )
  )
names(stock.nu_sal.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: pa_kab
stock.pa_kab.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(pa_kab!="NA")), #samples
        length(which(pa_kab>0)), #positive detections
        length(which(pa_kab>0))/length(which(!is.na(pa_kab))),  #calculates prevalence
        mean(pa_kab[pa_kab!=0], na.rm=TRUE),
        (length(which(pa_kab>0)) / length(which(!is.na(pa_kab)))) * mean(pa_kab[pa_kab!=0], na.rm=TRUE),
        "pa_kab"
      )
  )
names(stock.pa_kab.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: pa_min
stock.pa_min.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(pa_min!="NA")), #samples
        length(which(pa_min>0)), #positive detections
        length(which(pa_min>0))/length(which(!is.na(pa_min))),  #calculates prevalence
        mean(pa_min[pa_min!=0], na.rm=TRUE),
        (length(which(pa_min>0)) / length(which(!is.na(pa_min)))) * mean(pa_min[pa_min!=0], na.rm=TRUE),
        "pa_min"
      )
  )
names(stock.pa_min.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: pa_pse
stock.pa_pse.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(pa_pse!="NA")), #samples
        length(which(pa_pse>0)), #positive detections
        length(which(pa_pse>0))/length(which(!is.na(pa_pse))),  #calculates prevalence
        mean(pa_pse[pa_pse!=0], na.rm=TRUE),
        (length(which(pa_pse>0)) / length(which(!is.na(pa_pse)))) * mean(pa_pse[pa_pse!=0], na.rm=TRUE),
        "pa_pse"
      )
  )
names(stock.pa_pse.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: pa_ther
stock.pa_ther.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(pa_ther!="NA")), #samples
        length(which(pa_ther>0)), #positive detections
        length(which(pa_ther>0))/length(which(!is.na(pa_ther))),  #calculates prevalence
        mean(pa_ther[pa_ther!=0], na.rm=TRUE),
        (length(which(pa_ther>0)) / length(which(!is.na(pa_ther)))) * mean(pa_ther[pa_ther!=0], na.rm=TRUE),
        "pa_ther"
      )
  )
names(stock.pa_ther.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: pch_sal
stock.pch_sal.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(pch_sal!="NA")), #samples
        length(which(pch_sal>0)), #positive detections
        length(which(pch_sal>0))/length(which(!is.na(pch_sal))),  #calculates prevalence
        mean(pch_sal[pch_sal!=0], na.rm=TRUE),
        (length(which(pch_sal>0)) / length(which(!is.na(pch_sal)))) * mean(pch_sal[pch_sal!=0], na.rm=TRUE),
        "pch_sal"
      )
  )
names(stock.pch_sal.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: Picorna2_virus
stock.Picorna2_virus.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(Picorna2_virus!="NA")), #samples
        length(which(Picorna2_virus>0)), #positive detections
        length(which(Picorna2_virus>0))/length(which(!is.na(Picorna2_virus))),  #calculates prevalence
        mean(Picorna2_virus[Picorna2_virus!=0], na.rm=TRUE),
        (length(which(Picorna2_virus>0)) / length(which(!is.na(Picorna2_virus)))) * mean(Picorna2_virus[Picorna2_virus!=0], na.rm=TRUE),
        "Picorna2_virus"
      )
  )
names(stock.Picorna2_virus.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: prv
stock.prv.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(prv!="NA")), #samples
        length(which(prv>0)), #positive detections
        length(which(prv>0))/length(which(!is.na(prv))),  #calculates prevalence
        mean(prv[prv!=0], na.rm=TRUE),
        (length(which(prv>0)) / length(which(!is.na(prv)))) * mean(prv[prv!=0], na.rm=TRUE),
        "prv"
      )
  )
names(stock.prv.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: pspv
stock.pspv.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(pspv!="NA")), #samples
        length(which(pspv>0)), #positive detections
        length(which(pspv>0))/length(which(!is.na(pspv))),  #calculates prevalence
        mean(pspv[pspv!=0], na.rm=TRUE),
        (length(which(pspv>0)) / length(which(!is.na(pspv)))) * mean(pspv[pspv!=0], na.rm=TRUE),
        "pspv"
      )
  )
names(stock.pspv.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: Qin
stock.Qin.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(Qin!="NA")), #samples
        length(which(Qin>0)), #positive detections
        length(which(Qin>0))/length(which(!is.na(Qin))),  #calculates prevalence
        mean(Qin[Qin!=0], na.rm=TRUE),
        (length(which(Qin>0)) / length(which(!is.na(Qin)))) * mean(Qin[Qin!=0], na.rm=TRUE),
        "Qin"
      )
  )
names(stock.Qin.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: re_sal
stock.re_sal.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(re_sal!="NA")), #samples
        length(which(re_sal>0)), #positive detections
        length(which(re_sal>0))/length(which(!is.na(re_sal))),  #calculates prevalence
        mean(re_sal[re_sal!=0], na.rm=TRUE),
        (length(which(re_sal>0)) / length(which(!is.na(re_sal)))) * mean(re_sal[re_sal!=0], na.rm=TRUE),
        "re_sal"
      )
  )
names(stock.re_sal.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: Rhabdo3_virus
stock.Rhabdo3_virus.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(Rhabdo3_virus!="NA")), #samples
        length(which(Rhabdo3_virus>0)), #positive detections
        length(which(Rhabdo3_virus>0))/length(which(!is.na(Rhabdo3_virus))),  #calculates prevalence
        mean(Rhabdo3_virus[Rhabdo3_virus!=0], na.rm=TRUE),
        (length(which(Rhabdo3_virus>0)) / length(which(!is.na(Rhabdo3_virus)))) * mean(Rhabdo3_virus[Rhabdo3_virus!=0], na.rm=TRUE),
        "Rhabdo3_virus"
      )
  )
names(stock.Rhabdo3_virus.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: rlo
stock.rlo.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(rlo!="NA")), #samples
        length(which(rlo>0)), #positive detections
        length(which(rlo>0))/length(which(!is.na(rlo))),  #calculates prevalence
        mean(rlo[rlo!=0], na.rm=TRUE),
        (length(which(rlo>0)) / length(which(!is.na(rlo)))) * mean(rlo[rlo!=0], na.rm=TRUE),
        "rlo"
      )
  )
names(stock.rlo.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: sch
stock.sch.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(sch!="NA")), #samples
        length(which(sch>0)), #positive detections
        length(which(sch>0))/length(which(!is.na(sch))),  #calculates prevalence
        mean(sch[sch!=0], na.rm=TRUE),
        (length(which(sch>0)) / length(which(!is.na(sch)))) * mean(sch[sch!=0], na.rm=TRUE),
        "sch"
      )
  )
names(stock.sch.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: smallUK
stock.smallUK.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(smallUK!="NA")), #samples
        length(which(smallUK>0)), #positive detections
        length(which(smallUK>0))/length(which(!is.na(smallUK))),  #calculates prevalence
        mean(smallUK[smallUK!=0], na.rm=TRUE),
        (length(which(smallUK>0)) / length(which(!is.na(smallUK)))) * mean(smallUK[smallUK!=0], na.rm=TRUE),
        "smallUK"
      )
  )
names(stock.smallUK.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: sp_des
stock.sp_des.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(sp_des!="NA")), #samples
        length(which(sp_des>0)), #positive detections
        length(which(sp_des>0))/length(which(!is.na(sp_des))),  #calculates prevalence
        mean(sp_des[sp_des!=0], na.rm=TRUE),
        (length(which(sp_des>0)) / length(which(!is.na(sp_des)))) * mean(sp_des[sp_des!=0], na.rm=TRUE),
        "sp_des"
      )
  )
names(stock.sp_des.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: te_bry
stock.te_bry.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(te_bry!="NA")), #samples
        length(which(te_bry>0)), #positive detections
        length(which(te_bry>0))/length(which(!is.na(te_bry))),  #calculates prevalence
        mean(te_bry[te_bry!=0], na.rm=TRUE),
        (length(which(te_bry>0)) / length(which(!is.na(te_bry)))) * mean(te_bry[te_bry!=0], na.rm=TRUE),
        "te_bry"
      )
  )
names(stock.te_bry.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: te_mar
stock.te_mar.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(te_mar!="NA")), #samples
        length(which(te_mar>0)), #positive detections
        length(which(te_mar>0))/length(which(!is.na(te_mar))),  #calculates prevalence
        mean(te_mar[te_mar!=0], na.rm=TRUE),
        (length(which(te_mar>0)) / length(which(!is.na(te_mar)))) * mean(te_mar[te_mar!=0], na.rm=TRUE),
        "te_mar"
      )
  )
names(stock.te_mar.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: ven
stock.ven.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(ven!="NA")), #samples
        length(which(ven>0)), #positive detections
        length(which(ven>0))/length(which(!is.na(ven))),  #calculates prevalence
        mean(ven[ven!=0], na.rm=TRUE),
        (length(which(ven>0)) / length(which(!is.na(ven)))) * mean(ven[ven!=0], na.rm=TRUE),
        "ven"
      )
  )
names(stock.ven.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: richness
stock.richness.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(richness!="NA")), #samples
        length(which(richness>0)), #positive detections
        mean(richness[richness!=0], na.rm=TRUE),#/6.285714, #we want the mean richness to be relative to data
        median(richness[richness!=0], na.rm=TRUE),
        (length(which(richness>0)) / length(which(!is.na(richness)))) * mean(richness[richness!=0], na.rm=TRUE),
        "richness"
      )
  )
names(stock.richness.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns

##For now, by agent: rib
stock.rib.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(rib!="NA")), #samples
        length(which(rib>0)), #positive detections
        mean(rib[rib!=0], na.rm=TRUE),#/0.3328581777, #we want the RIB to be relative
        median(rib[rib!=0], na.rm=TRUE),
        (length(which(rib>0)) / length(which(!is.na(rib)))) * mean(rib[rib!=0], na.rm=TRUE),
        "rib"
      )
  )
names(stock.rib.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload","agent") #rename columns


## Merge individual pathogen dfs into one file for evaluation independent of SR data
full.sw <- rbind(stock.arena1.sw, stock.arena2.sw, stock.ascv.sw,
  stock.c_b_cys.sw, stock.ce_sha.sw, stock.Circo.virus.sw, stock.cr_sal.sw, stock.de_sal.sw,
  stock.fa_mar.sw, stock.fl_psy.sw, stock.ic_hof.sw, stock.ic_mul.sw, stock.IcD.sw, stock.ihnv.sw,
  stock.ku_thy.sw, stock.lo_sal.sw, stock.mo_vis.sw, stock.my_arc.sw, stock.nu_sal.sw, stock.pa_kab.sw,
  stock.pa_min.sw, stock.pa_pse.sw, stock.pa_ther.sw, stock.pch_sal.sw, stock.Picorna2_virus.sw, stock.pspv.sw,
  stock.prv.sw, stock.Qin.sw, stock.re_sal.sw, stock.Rhabdo3_virus.sw, stock.rlo.sw, stock.sch.sw, stock.smallUK.sw, 
  stock.sp_des.sw, stock.te_bry.sw, stock.te_mar.sw, stock.ven.sw,
  stock.richness.sw, stock.rib.sw)

## Must be detected in >3 years and >5 fish to be included in survival analysis
#visualize agent detections across years
ggplot(full.sw) +  
  aes(agent, prev, color=as.factor(Year)) +
  geom_point() +
  ylim(0,1) +
  coord_flip()
sw.totals <- aggregate(full.sw$posdet, by=list(Category=full.sw$agent), FUN=sum)
sw.totals[order(sw.totals[,2]),]
temp <-
  full.sw %>% 
  group_by(Year, agent) %>%
  summarise(posdet = length(which(posdet>0)))
write.csv(temp, file="data/posdet per year sw.csv")
# Present in <3 years: Remove
  ## arena1, mo_vis, Picorna2_virus, re_sal, ascv, cr_sal, pch_sal

## Merge individual pathogen dfs with SR data
arena1.resid.sws <- merge(trnc_resid_srr, stock.arena1.sw, by = c("Stock_Analysis","Year"))
arena2.resid.sws <- merge(trnc_resid_srr, stock.arena2.sw, by = c("Stock_Analysis","Year"))
ascv.resid.sws <- merge(trnc_resid_srr, stock.ascv.sw, by = c("Stock_Analysis","Year"))
c_b_cys.resid.sws <- merge(trnc_resid_srr, stock.c_b_cys.sw, by = c("Stock_Analysis","Year"))
ce_sha.resid.sws <- merge(trnc_resid_srr, stock.ce_sha.sw, by = c("Stock_Analysis","Year"))
Circo.virus.resid.sws <- merge(trnc_resid_srr, stock.Circo.virus.sw, by = c("Stock_Analysis","Year"))
cr_sal.resid.sws <- merge(trnc_resid_srr, stock.cr_sal.sw, by = c("Stock_Analysis","Year"))
de_sal.resid.sws <- merge(trnc_resid_srr, stock.de_sal.sw, by = c("Stock_Analysis","Year"))
fa_mar.resid.sws <- merge(trnc_resid_srr, stock.fa_mar.sw, by = c("Stock_Analysis","Year"))
fl_psy.resid.sws <- merge(trnc_resid_srr, stock.fl_psy.sw, by = c("Stock_Analysis","Year"))
ic_hof.resid.sws <- merge(trnc_resid_srr, stock.ic_hof.sw, by = c("Stock_Analysis","Year"))
ic_mul.resid.sws <- merge(trnc_resid_srr, stock.ic_mul.sw, by = c("Stock_Analysis","Year"))
IcD.resid.sws <- merge(trnc_resid_srr, stock.IcD.sw, by = c("Stock_Analysis","Year"))
ihnv.resid.sws <- merge(trnc_resid_srr, stock.ihnv.sw, by = c("Stock_Analysis","Year"))
ku_thy.resid.sws <- merge(trnc_resid_srr, stock.ku_thy.sw, by = c("Stock_Analysis","Year"))
lo_sal.resid.sws <- merge(trnc_resid_srr, stock.lo_sal.sw, by = c("Stock_Analysis","Year"))
mo_vis.resid.sws <- merge(trnc_resid_srr, stock.mo_vis.sw, by = c("Stock_Analysis","Year"))
my_arc.resid.sws <- merge(trnc_resid_srr, stock.my_arc.sw, by = c("Stock_Analysis","Year"))
nu_sal.resid.sws <- merge(trnc_resid_srr, stock.nu_sal.sw, by = c("Stock_Analysis","Year"))
pa_kab.resid.sws <- merge(trnc_resid_srr, stock.pa_kab.sw, by = c("Stock_Analysis","Year"))
pa_min.resid.sws <- merge(trnc_resid_srr, stock.pa_min.sw, by = c("Stock_Analysis","Year"))
pa_pse.resid.sws <- merge(trnc_resid_srr, stock.pa_pse.sw, by = c("Stock_Analysis","Year"))
pa_ther.resid.sws <- merge(trnc_resid_srr, stock.pa_ther.sw, by = c("Stock_Analysis","Year"))
pch_sal.resid.sws <- merge(trnc_resid_srr, stock.pch_sal.sw, by = c("Stock_Analysis","Year"))
Picorna2_virus.resid.sws <- merge(trnc_resid_srr, stock.Picorna2_virus.sw, by = c("Stock_Analysis","Year"))
prv.resid.sws <- merge(trnc_resid_srr, stock.prv.sw, by = c("Stock_Analysis","Year"))
Qin.resid.sws <- merge(trnc_resid_srr, stock.Qin.sw, by = c("Stock_Analysis","Year"))
re_sal.resid.sws <- merge(trnc_resid_srr, stock.re_sal.sw, by = c("Stock_Analysis","Year"))
Rhabdo3_virus.resid.sws <- merge(trnc_resid_srr, stock.Rhabdo3_virus.sw, by = c("Stock_Analysis","Year"))
pspv.resid.sws <- merge(trnc_resid_srr, stock.pspv.sw, by = c("Stock_Analysis","Year"))
rlo.resid.sws <- merge(trnc_resid_srr, stock.rlo.sw, by = c("Stock_Analysis","Year"))
sch.resid.sws <- merge(trnc_resid_srr, stock.sch.sw, by = c("Stock_Analysis","Year"))
smallUK.resid.sws <- merge(trnc_resid_srr, stock.smallUK.sw, by = c("Stock_Analysis","Year"))
sp_des.resid.sws <- merge(trnc_resid_srr, stock.sp_des.sw, by = c("Stock_Analysis","Year"))
te_bry.resid.sws <- merge(trnc_resid_srr, stock.te_bry.sw, by = c("Stock_Analysis","Year"))
te_mar.resid.sws <- merge(trnc_resid_srr, stock.te_mar.sw, by = c("Stock_Analysis","Year"))
ven.resid.sws <- merge(trnc_resid_srr, stock.ven.sw, by = c("Stock_Analysis","Year"))
richness.resid.sws <- merge(trnc_resid_srr, stock.richness.sw, by = c("Stock_Analysis","Year"))
rib.resid.sws <- merge(trnc_resid_srr, stock.rib.sw, by = c("Stock_Analysis","Year"))

## Merge individual pathogen-SR data dfs into one df 
# Present in <3 years so Remove:
## arena1, mo_vis, Picorna2_virus, re_sal, ascv, cr_sal, pch_sal
## Note smallUK also removed after <4 samples per stock-agent-year combo applied, so remove
full.resid.sws <- rbind(
  arena2.resid.sws, c_b_cys.resid.sws, ce_sha.resid.sws, Circo.virus.resid.sws, 
  de_sal.resid.sws, fa_mar.resid.sws, fl_psy.resid.sws, ic_hof.resid.sws, ic_mul.resid.sws,
  IcD.resid.sws, ihnv.resid.sws, ku_thy.resid.sws, lo_sal.resid.sws, my_arc.resid.sws, 
  nu_sal.resid.sws, pa_kab.resid.sws, pa_min.resid.sws, pa_pse.resid.sws, pa_ther.resid.sws, 
  pspv.resid.sws, prv.resid.sws, Qin.resid.sws, rlo.resid.sws, Rhabdo3_virus.resid.sws, 
  sch.resid.sws, sp_des.resid.sws, te_bry.resid.sws, te_mar.resid.sws, ven.resid.sws)
rr.resid.sws  <- rbind(richness.resid.sws, rib.resid.sws)
write.csv(full.resid.sws, file="data/FULL_ONNE_agents only_SW_220126.csv")

# Cumulative metrics
write.csv(rr.resid.sws, file="data/RR_ONNE_cumulative metrics_SW_220126.csv")

## Remove prevalence values from estimates with <4 fish (N<10), but not for rib or richness
temp <- full.resid.sws
new.df.sw <- subset(temp, N > 9)

## Rename and save
inf_agt_resid_data <- new.df.sw
write.csv(inf_agt_resid_data, file="data/REDUCED_ONNE_agents only_SW_220126.csv")


## Data Exploration ###
# How many positive detections per year in SR dataset?
posi.det.total <- aggregate(posdet ~ agent, inf_agt_resid_data, sum)
posi.det.annual <- aggregate(posdet ~ Year + agent, inf_agt_resid_data, sum)
ggplot(data = posi.det.annual, aes(x=factor(Year), y=posdet, group=agent)) +
  geom_area(aes(fill=agent))

# How many positive detections overall in SR dataset? Ranked
posi.det.total[order(posi.det.total$posdet),]

# How many individual assays run?
sum(inf_agt_resid_data$N) #a shit ton




####################################################################################
### FW data cleaning
####################################################################################
result.fw
result.fw[order(result.fw[,1]),]

####All data in one table - Stock-specific metric - SW only
fw.stock <-fw.data %>% group_by(Stock_Analysis, Year) #create object to be summarized by year

##ASSESS INFECTION PARAMETERS BY AGENT:
##ae_sal
stock.arena2.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(arena2!="NA")), #samples
        length(which(arena2>0)), #positive detections
        length(which(arena2>0))/length(which(!is.na(arena2))),  #calculates prevalence
        mean(arena2[arena2!=0], na.rm=TRUE),
        (length(which(arena2>0)) / length(which(!is.na(arena2)))) * mean(arena2[arena2!=0], na.rm=TRUE)
      )
  )
names(stock.arena2.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##c_b_cys
stock.c_b_cys.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(c_b_cys!="NA")), #samples
        length(which(c_b_cys>0)), #positive detections
        length(which(c_b_cys>0))/length(which(!is.na(c_b_cys))),  #calculates prevalence
        mean(c_b_cys[c_b_cys!=0], na.rm=TRUE),
        (length(which(c_b_cys>0)) / length(which(!is.na(c_b_cys)))) * mean(ae_hyd[c_b_cys!=0], na.rm=TRUE)
      )
  )
names(stock.c_b_cys.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: ce_sha
stock.ce_sha.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(ce_sha!="NA")), #samples
        length(which(ce_sha>0)), #positive detections
        length(which(ce_sha>0))/length(which(!is.na(ce_sha))),  #calculates prevalence
        mean(ce_sha[ce_sha!=0], na.rm=TRUE),
        (length(which(ce_sha>0)) / length(which(!is.na(ce_sha)))) * mean(ce_sha[ce_sha!=0], na.rm=TRUE)
      )
  )
names(stock.ce_sha.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: cr_sal
stock.cr_sal.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(cr_sal!="NA")), #samples
        length(which(cr_sal>0)), #positive detections
        length(which(cr_sal>0))/length(which(!is.na(cr_sal))),  #calculates prevalence
        mean(cr_sal[cr_sal!=0], na.rm=TRUE),
        (length(which(cr_sal>0)) / length(which(!is.na(cr_sal)))) * mean(cr_sal[cr_sal!=0], na.rm=TRUE)
      )
  )
names(stock.cr_sal.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: Circo.virus
stock.Circo.virus.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(Circo.virus!="NA")), #samples
        length(which(Circo.virus>0)), #positive detections
        length(which(Circo.virus>0))/length(which(!is.na(Circo.virus))),  #calculates prevalence
        mean(Circo.virus[Circo.virus!=0], na.rm=TRUE),
        (length(which(Circo.virus>0)) / length(which(!is.na(Circo.virus)))) * mean(Circo.virus[Circo.virus!=0], na.rm=TRUE)
      )
  )
names(stock.Circo.virus.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: de_sal
stock.de_sal.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(de_sal!="NA")), #samples
        length(which(de_sal>0)), #positive detections
        length(which(de_sal>0))/length(which(!is.na(de_sal))),  #calculates prevalence
        mean(de_sal[de_sal!=0], na.rm=TRUE),
        (length(which(de_sal>0)) / length(which(!is.na(de_sal)))) * mean(de_sal[de_sal!=0], na.rm=TRUE)
      )
  )
names(stock.de_sal.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: fl_psy
stock.fl_psy.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(fl_psy!="NA")), #samples
        length(which(fl_psy>0)), #positive detections
        length(which(fl_psy>0))/length(which(!is.na(fl_psy))),  #calculates prevalence
        mean(fl_psy[fl_psy!=0], na.rm=TRUE),
        (length(which(fl_psy>0)) / length(which(!is.na(fl_psy)))) * mean(fl_psy[fl_psy!=0], na.rm=TRUE)
      )
  )
names(stock.fl_psy.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: ic_hof
stock.ic_hof.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(ic_hof!="NA")), #samples
        length(which(ic_hof>0)), #positive detections
        length(which(ic_hof>0))/length(which(!is.na(ic_hof))),  #calculates prevalence
        mean(ic_hof[ic_hof!=0], na.rm=TRUE),
        (length(which(ic_hof>0)) / length(which(!is.na(ic_hof)))) * mean(ic_hof[ic_hof!=0], na.rm=TRUE)
      )
  )
names(stock.ic_hof.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: ic_mul
stock.ic_mul.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(ic_mul!="NA")), #samples
        length(which(ic_mul>0)), #positive detections
        length(which(ic_mul>0))/length(which(!is.na(ic_mul))),  #calculates prevalence
        mean(ic_mul[ic_mul!=0], na.rm=TRUE),
        (length(which(ic_mul>0)) / length(which(!is.na(ic_mul)))) * mean(ic_mul[ic_mul!=0], na.rm=TRUE)
      )
  )
names(stock.ic_mul.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: IcD
stock.IcD.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(IcD!="NA")), #samples
        length(which(IcD>0)), #positive detections
        length(which(IcD>0))/length(which(!is.na(IcD))),  #calculates prevalence
        mean(IcD[IcD!=0], na.rm=TRUE),
        (length(which(IcD>0)) / length(which(!is.na(IcD)))) * mean(IcD[IcD!=0], na.rm=TRUE)
      )
  )
names(stock.IcD.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: ihnv
stock.ihnv.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(ihnv!="NA")), #samples
        length(which(ihnv>0)), #positive detections
        length(which(ihnv>0))/length(which(!is.na(ihnv))),  #calculates prevalence
        mean(ihnv[ihnv!=0], na.rm=TRUE),
        (length(which(ihnv>0)) / length(which(!is.na(ihnv)))) * mean(ihnv[ihnv!=0], na.rm=TRUE)
      )
  )
names(stock.ihnv.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: lo_sal
stock.lo_sal.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(lo_sal!="NA")), #samples
        length(which(lo_sal>0)), #positive detections
        length(which(lo_sal>0))/length(which(!is.na(lo_sal))),  #calculates prevalence
        mean(lo_sal[lo_sal!=0], na.rm=TRUE),
        (length(which(lo_sal>0)) / length(which(!is.na(lo_sal)))) * mean(lo_sal[lo_sal!=0], na.rm=TRUE)
      )
  )
names(stock.lo_sal.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: my_arc
stock.my_arc.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(my_arc!="NA")), #samples
        length(which(my_arc>0)), #positive detections
        length(which(my_arc>0))/length(which(!is.na(my_arc))),  #calculates prevalence
        mean(my_arc[my_arc!=0], na.rm=TRUE),
        (length(which(my_arc>0)) / length(which(!is.na(my_arc)))) * mean(my_arc[my_arc!=0], na.rm=TRUE)
      )
  )
names(stock.my_arc.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: pa_kab
stock.pa_kab.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(pa_kab!="NA")), #samples
        length(which(pa_kab>0)), #positive detections
        length(which(pa_kab>0))/length(which(!is.na(pa_kab))),  #calculates prevalence
        mean(pa_kab[pa_kab!=0], na.rm=TRUE),
        (length(which(pa_kab>0)) / length(which(!is.na(pa_kab)))) * mean(pa_kab[pa_kab!=0], na.rm=TRUE)
      )
  )
names(stock.pa_kab.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: pa_min
stock.pa_min.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(pa_min!="NA")), #samples
        length(which(pa_min>0)), #positive detections
        length(which(pa_min>0))/length(which(!is.na(pa_min))),  #calculates prevalence
        mean(pa_min[pa_min!=0], na.rm=TRUE),
        (length(which(pa_min>0)) / length(which(!is.na(pa_min)))) * mean(pa_min[pa_min!=0], na.rm=TRUE)
      )
  )
names(stock.pa_min.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: pa_pse
stock.pa_pse.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(pa_pse!="NA")), #samples
        length(which(pa_pse>0)), #positive detections
        length(which(pa_pse>0))/length(which(!is.na(pa_pse))),  #calculates prevalence
        mean(pa_pse[pa_pse!=0], na.rm=TRUE),
        (length(which(pa_pse>0)) / length(which(!is.na(pa_pse)))) * mean(pa_pse[pa_pse!=0], na.rm=TRUE)
      )
  )
names(stock.pa_pse.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: pa_ther
stock.pa_ther.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(pa_ther!="NA")), #samples
        length(which(pa_ther>0)), #positive detections
        length(which(pa_ther>0))/length(which(!is.na(pa_ther))),  #calculates prevalence
        mean(pa_ther[pa_ther!=0], na.rm=TRUE),
        (length(which(pa_ther>0)) / length(which(!is.na(pa_ther)))) * mean(pa_ther[pa_ther!=0], na.rm=TRUE)
      )
  )
names(stock.pa_ther.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: pch_sal
stock.pch_sal.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(pch_sal!="NA")), #samples
        length(which(pch_sal>0)), #positive detections
        length(which(pch_sal>0))/length(which(!is.na(pch_sal))),  #calculates prevalence
        mean(pch_sal[pch_sal!=0], na.rm=TRUE),
        (length(which(pch_sal>0)) / length(which(!is.na(pch_sal)))) * mean(pch_sal[pch_sal!=0], na.rm=TRUE)
      )
  )
names(stock.pch_sal.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: prv
stock.prv.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(prv!="NA")), #samples
        length(which(prv>0)), #positive detections
        length(which(prv>0))/length(which(!is.na(prv))),  #calculates prevalence
        mean(prv[prv!=0], na.rm=TRUE),
        (length(which(prv>0)) / length(which(!is.na(prv)))) * mean(prv[prv!=0], na.rm=TRUE)
      )
  )
names(stock.prv.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: pspv
stock.pspv.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(pspv!="NA")), #samples
        length(which(pspv>0)), #positive detections
        length(which(pspv>0))/length(which(!is.na(pspv))),  #calculates prevalence
        mean(pspv[pspv!=0], na.rm=TRUE),
        (length(which(pspv>0)) / length(which(!is.na(pspv)))) * mean(pspv[pspv!=0], na.rm=TRUE)
      )
  )
names(stock.pspv.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: Qin
stock.Qin.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(Qin!="NA")), #samples
        length(which(Qin>0)), #positive detections
        length(which(Qin>0))/length(which(!is.na(Qin))),  #calculates prevalence
        mean(Qin[Qin!=0], na.rm=TRUE),
        (length(which(Qin>0)) / length(which(!is.na(Qin)))) * mean(Qin[Qin!=0], na.rm=TRUE)
      )
  )
names(stock.Qin.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: rlo
stock.rlo.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(rlo!="NA")), #samples
        length(which(rlo>0)), #positive detections
        length(which(rlo>0))/length(which(!is.na(rlo))),  #calculates prevalence
        mean(rlo[rlo!=0], na.rm=TRUE),
        (length(which(rlo>0)) / length(which(!is.na(rlo)))) * mean(rlo[rlo!=0], na.rm=TRUE)
      )
  )
names(stock.rlo.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: sch
stock.sch.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(sch!="NA")), #samples
        length(which(sch>0)), #positive detections
        length(which(sch>0))/length(which(!is.na(sch))),  #calculates prevalence
        mean(sch[sch!=0], na.rm=TRUE),
        (length(which(sch>0)) / length(which(!is.na(sch)))) * mean(sch[sch!=0], na.rm=TRUE)
      )
  )
names(stock.sch.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: smallUK
stock.smallUK.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(smallUK!="NA")), #samples
        length(which(smallUK>0)), #positive detections
        length(which(smallUK>0))/length(which(!is.na(smallUK))),  #calculates prevalence
        mean(smallUK[smallUK!=0], na.rm=TRUE),
        (length(which(smallUK>0)) / length(which(!is.na(smallUK)))) * mean(smallUK[smallUK!=0], na.rm=TRUE)
      )
  )
names(stock.smallUK.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: sp_des
stock.sp_des.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(sp_des!="NA")), #samples
        length(which(sp_des>0)), #positive detections
        length(which(sp_des>0))/length(which(!is.na(sp_des))),  #calculates prevalence
        mean(sp_des[sp_des!=0], na.rm=TRUE),
        (length(which(sp_des>0)) / length(which(!is.na(sp_des)))) * mean(sp_des[sp_des!=0], na.rm=TRUE)
      )
  )
names(stock.sp_des.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: te_bry
stock.te_bry.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(te_bry!="NA")), #samples
        length(which(te_bry>0)), #positive detections
        length(which(te_bry>0))/length(which(!is.na(te_bry))),  #calculates prevalence
        mean(te_bry[te_bry!=0], na.rm=TRUE),
        (length(which(te_bry>0)) / length(which(!is.na(te_bry)))) * mean(te_bry[te_bry!=0], na.rm=TRUE)
      )
  )
names(stock.te_bry.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: te_mar
stock.te_mar.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(te_mar!="NA")), #samples
        length(which(te_mar>0)), #positive detections
        length(which(te_mar>0))/length(which(!is.na(te_mar))),  #calculates prevalence
        mean(te_mar[te_mar!=0], na.rm=TRUE),
        (length(which(te_mar>0)) / length(which(!is.na(te_mar)))) * mean(te_mar[te_mar!=0], na.rm=TRUE)
      )
  )
names(stock.te_mar.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: ven
stock.ven.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(ven!="NA")), #samples
        length(which(ven>0)), #positive detections
        length(which(ven>0))/length(which(!is.na(ven))),  #calculates prevalence
        mean(ven[ven!=0], na.rm=TRUE),
        (length(which(ven>0)) / length(which(!is.na(ven)))) * mean(ven[ven!=0], na.rm=TRUE)
      )
  )
names(stock.ven.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

stock.richness.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(richness!="NA")), #samples
        length(which(richness>0)), #positive detections
        mean(richness[richness!=0], na.rm=TRUE), #mean richness is in prevalence column!
        median(richness[richness!=0], na.rm=TRUE),
        (length(which(richness>0)) / length(which(!is.na(richness)))) * mean(richness[richness!=0], na.rm=TRUE)
      )
  )
names(stock.richness.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

stock.rib.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(rib!="NA")), #samples
        length(which(rib>0)), #positive detections
        mean(rib[rib!=0], na.rm=TRUE), #mean rib is in prevalence column!
        median(rib[rib!=0], na.rm=TRUE),
        (length(which(rib>0)) / length(which(!is.na(rib)))) * mean(rib[rib!=0], na.rm=TRUE)
      )
  )
names(stock.rib.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

## Merge individual pathogen dfs with SR data
arena2.resid.fws <- merge(trnc_resid_srr, stock.arena2.fw, by = c("Stock_Analysis","Year"))
c_b_cys.resid.fws <- merge(trnc_resid_srr, stock.c_b_cys.fw, by = c("Stock_Analysis","Year"))
ce_sha.resid.fws <- merge(trnc_resid_srr, stock.ce_sha.fw, by = c("Stock_Analysis","Year"))
Circo.virus.resid.fws <- merge(trnc_resid_srr, stock.Circo.virus.fw, by = c("Stock_Analysis","Year"))
cr_sal.resid.fws <- merge(trnc_resid_srr, stock.cr_sal.fw, by = c("Stock_Analysis","Year"))
de_sal.resid.fws <- merge(trnc_resid_srr, stock.de_sal.fw, by = c("Stock_Analysis","Year"))
fl_psy.resid.fws <- merge(trnc_resid_srr, stock.fl_psy.fw, by = c("Stock_Analysis","Year"))
ic_hof.resid.fws <- merge(trnc_resid_srr, stock.ic_hof.fw, by = c("Stock_Analysis","Year"))
ic_mul.resid.fws <- merge(trnc_resid_srr, stock.ic_mul.fw, by = c("Stock_Analysis","Year"))
IcD.resid.fws <- merge(trnc_resid_srr, stock.IcD.fw, by = c("Stock_Analysis","Year"))
ihnv.resid.fws <- merge(trnc_resid_srr, stock.ihnv.fw, by = c("Stock_Analysis","Year"))
lo_sal.resid.fws <- merge(trnc_resid_srr, stock.lo_sal.fw, by = c("Stock_Analysis","Year"))
my_arc.resid.fws <- merge(trnc_resid_srr, stock.my_arc.fw, by = c("Stock_Analysis","Year"))
pa_kab.resid.fws <- merge(trnc_resid_srr, stock.pa_kab.fw, by = c("Stock_Analysis","Year"))
pa_min.resid.fws <- merge(trnc_resid_srr, stock.pa_min.fw, by = c("Stock_Analysis","Year"))
pa_pse.resid.fws <- merge(trnc_resid_srr, stock.pa_pse.fw, by = c("Stock_Analysis","Year"))
pa_ther.resid.fws <- merge(trnc_resid_srr, stock.pa_ther.fw, by = c("Stock_Analysis","Year"))
pch_sal.resid.fws <- merge(trnc_resid_srr, stock.pch_sal.fw, by = c("Stock_Analysis","Year"))
prv.resid.fws <- merge(trnc_resid_srr, stock.prv.fw, by = c("Stock_Analysis","Year"))
Qin.resid.fws <- merge(trnc_resid_srr, stock.Qin.fw, by = c("Stock_Analysis","Year"))
pspv.resid.fws <- merge(trnc_resid_srr, stock.pspv.fw, by = c("Stock_Analysis","Year"))
rlo.resid.fws <- merge(trnc_resid_srr, stock.rlo.fw, by = c("Stock_Analysis","Year"))
sch.resid.fws <- merge(trnc_resid_srr, stock.sch.fw, by = c("Stock_Analysis","Year"))
smallUK.resid.fws <- merge(trnc_resid_srr, stock.smallUK.fw, by = c("Stock_Analysis","Year"))
sp_des.resid.fws <- merge(trnc_resid_srr, stock.sp_des.fw, by = c("Stock_Analysis","Year"))
te_bry.resid.fws <- merge(trnc_resid_srr, stock.te_bry.fw, by = c("Stock_Analysis","Year"))
te_mar.resid.fws <- merge(trnc_resid_srr, stock.te_mar.fw, by = c("Stock_Analysis","Year"))
ven.resid.fws <- merge(trnc_resid_srr, stock.ven.fw, by = c("Stock_Analysis","Year"))
richness.resid.fws <- merge(trnc_resid_srr, stock.richness.fw, by = c("Stock_Analysis","Year"))
rib.resid.fws <- merge(trnc_resid_srr, stock.rib.fw, by = c("Stock_Analysis","Year"))

##Add column with agent name
arena2.resid.fws$agent<-"arena2"
c_b_cys.resid.fws$agent<-"c_b_cys"
ce_sha.resid.fws$agent<-"ce_sha"
Circo.virus.resid.fws$agent<-"Circo.virus"
de_sal.resid.fws$agent<-"de_sal"
fl_psy.resid.fws$agent<-"fl_psy"
ic_hof.resid.fws$agent<-"ic_hof"
ic_mul.resid.fws$agent<-"ic_mul"
IcD.resid.fws$agent<-"IcD"
ihnv.resid.fws$agent<- "ihnv"
lo_sal.resid.fws$agent<-"lo_sal"
my_arc.resid.fws$agent<-"my_arc"
pa_kab.resid.fws$agent<-"pa_kab"
pa_min.resid.fws$agent<-"pa_min"
pa_pse.resid.fws$agent<-"pa_pse"
pa_ther.resid.fws$agent<-"pa_ther"
pch_sal.resid.fws$agent<-"pch_sal"
prv.resid.fws$agent<-"prv"
pspv.resid.fws$agent<-"pspv"
Qin.resid.fws$agent<-"Qin"
rlo.resid.fws$agent<-"rlo"
sch.resid.fws$agent<-"sch"
smallUK.resid.fws$agent<-"smallUK"
sp_des.resid.fws$agent<-"sp_des"
te_bry.resid.fws$agent<-"te_bry"
te_mar.resid.fws$agent<-"te_mar"
ven.resid.fws$agent<-"ven"
richness.resid.fws$agent<-"richness"
rib.resid.fws$agent<-"rib"

##Merge individual pathogen-SR data dfs into one 
# Present in <3 years: Remove
## Circo.virus,ihnv,pa_kab,pch_sal,prv,te_mar,ven,smallUK
full.resid.fws <- rbind(
  arena2.resid.fws, c_b_cys.resid.fws, ce_sha.resid.fws, de_sal.resid.fws,
  fl_psy.resid.fws, ic_hof.resid.fws, ic_mul.resid.fws, IcD.resid.fws, 
  lo_sal.resid.fws, my_arc.resid.fws, 
  pa_min.resid.fws, pa_pse.resid.fws, pa_ther.resid.fws,  pspv.resid.fws,
  Qin.resid.fws, rlo.resid.fws, sch.resid.fws,
  sp_des.resid.fws, te_bry.resid.fws)
rr.resid.fws <-  rbind(richness.resid.fws, rib.resid.fws)

## Remove Harrison from FW data - ocean type
# Pathogens only
full.resid.fws <- full.resid.fws[full.resid.fws$Stock_Analysis !=  "Harrison",]


## Must be detected in >3 years and >5 fish to be included in survival analysis
#visualize agent detections across years
ggplot(full.resid.fws) +  
  aes(agent, prev, color=as.factor(Year)) +
  geom_point() +
  ylim(0,1) +
  coord_flip()
fw.totals <- aggregate(full.resid.fws$posdet, by=list(Category=full.resid.fws$agent), FUN=sum)
fw.totals[order(fw.totals[,2]),]
temp <-
  full.resid.fws %>% 
  group_by(Year, agent) %>%
  summarise(posdet = length(which(posdet>0)))
write.csv(temp, file="data/posdet per year fw.csv")
# Present in <3 years: Remove
## Circo.virus,ihnv,pa_kab,pch_sal,prv,te_mar,ven,smallUK

write.csv(full.resid.fws, file="data/FULL_ONNE_agents only_FW_220126.csv")

# Cumulative metrics
rr.resid.fws <- rr.resid.fws[rr.resid.fws$Stock_Analysis !=  "Harrison",]
write.csv(rr.resid.fws, file="data/RR_ONNE_cumulative metrics_FW_220126.csv")

## Remove prevalence values from estimates with <10 fish (N<10), but not for rib or richness
temp <- full.resid.fws
new_df.fw <- subset(temp, N > 9) 

# Rename and save
inf_agt_resid_data_fw <- new_df.fw
write.csv(inf_agt_resid_data_fw, file="data/REDUCED_ONNE_agents only_FW_220126.csv")

## Print out dataframes in "final" folder for tables
dim(sw.data)
dim(fw.data)
write.csv(fw.data, file="figs/Final data frames and tables/FWdata_230123.csv")
write.csv(sw.data, file="figs/Final data frames and tables/SWdata_230123.csv")

# Agent Names file
# Agent names file
agent.names <- read.csv("data/agent.names_230124.csv")
agent.names <- agent.names[,-1]
agent.names




#########################################
####All data in one table - Stock-specific metric - chilko only
ch.data <- fw.data[fw.data$Stock_Analysis=="Chilko",]
ch.dataL<-gather(ch.data,agent,value,ae_hyd:rib)
head(ch.dataL)

##ASSESS INFECTION PARAMETERS BY AGENT:
ch.year <-ch.dataL %>% group_by(Year,agent) #create object to be summarized by year
ch.fw.prev =
  data.frame(
    ch.year %>% 
      summarise(
        length(which(value!="NA")), #samples
        length(which(value>0)), #positive detections
        length(which(value>0))/length(which(!is.na(value))),  #calculates prevalence
        mean(value[value!=0], na.rm=TRUE),
      )
  )
names(ch.fw.prev) <- c("Year", "Agent", "N", "posdet", "prev", "mean_load") #rename columns
head(ch.fw.prev)
unique(ch.fw.prev$Year)
write.csv(ch.fw.prev, file="data/chilko_fw_prev.csv")


#####################################
# Calculate range in RIB and Richness
head(sw.data)
head(fw.data)
RIB.range.sw <- c(min(sw.data$rib), max(sw.data$rib), mean(sw.data$rib))
RIB.range.fw <- c(min(fw.data$rib), max(fw.data$rib), mean(fw.data$rib))
Richness.range.sw <- c(min(sw.data$richness), max(sw.data$richness), median(sw.data$richness))
Richness.range.fw <- c(min(fw.data$richness), max(fw.data$richness), median(fw.data$richness))

