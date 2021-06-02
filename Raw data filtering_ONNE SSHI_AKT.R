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
trnc_resid<-read.csv("data/survival_indices_truncated_210519.csv", head=TRUE)
str(trnc_resid)
trnc_resid$Year <- trnc_resid$brood_year+2
names(trnc_resid) <- c("orderID", "Stock_Analysis", "brood_year", "metric", "resid_value", "Year") #rename columns
##create object with just SRR_resid metric to align with infection data
trnc_resid_srr <- trnc_resid[trnc_resid$metric=="SR_resid",]
levels(trnc_resid$Stock_Analysis)

#Bring in pathogen data without LOD and remove samples with no stock ID (may be updated data available!)
all <- read.csv("data/ONNE metadata no LOD_5.7.2021.csv",header=TRUE)
str(all)
all2 <- droplevels(all[!(all$Stock_Region=="") ,])#remove any fish without stock assignment for now
dim(all2)

#take out rare stock regions (northern, QCI, Transboundary)
temp2<-droplevels(all2[-which(all2$Stock_Area=="Northern") ,]) #Northern out
dim(temp2)
temp3<-droplevels(temp2[-which(temp2$Stock_Area=="QCI") ,]) #QCI out
dim(temp3)
major<-droplevels(temp3[-which(temp3$Stock_Area=="TransBoundary") ,]) #Transboundary out
dim(major)

## Divide up SW and FW collected samples
sw.major<-droplevels(major[(major$SWFW=="SW") ,])#SW only
head(sw.major)
fw.major<-droplevels(major[(major$SWFW=="FW") ,])#FW only 
head(fw.data)

## Temporally limit: Reduce sampling period to spring-summer only
spsu1<-droplevels(sw.major[-which(sw.major$SEASON1=="Overwinter") ,]) #remove winter 
dim(spsu1)
spsu2<-droplevels(spsu1[-which(spsu1$SEASON1=="Fall") ,]) # remove fall 
dim(spsu2) 
#spsu3<-droplevels(spsu2[-which(spsu2$Year=="2018") ,]) #remove 2018 - no SR data yet
#dim(spsu3)

## Geographically limit: Remove samples from WCVI, 2018 and high latitudes (>51.5 lat) 
spsu4<-droplevels(spsu2[-which(spsu2$Zone=="WCVI") ,]) #remove WCVI 
dim(spsu4)
spsu<-droplevels(spsu4[-which(spsu4$Latitude > 51.5) ,]) # remove high latitude samples
dim(spsu)

## Change stock names to align with SR data in both SW and FW data
levels(spsu$Stock_Analysis)[levels(spsu$Stock_Analysis)=="Early Stuart"] <- "E.Stuart"
levels(spsu$Stock_Analysis)[levels(spsu$Stock_Analysis)=="Late Shuswap"] <- "L.Shuswap"
levels(spsu$Stock_Analysis)[levels(spsu$Stock_Analysis)=="Late Stuart"] <- "L.Stuart"
levels(spsu$Stock_Analysis)[levels(spsu$Stock_Analysis)=="Harrison-Widgeon"] <- "Harrison"
sw.data <- spsu #rename
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Early Stuart"] <- "E.Stuart"
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Late Shuswap"] <- "L.Shuswap"
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Late Stuart"] <- "L.Stuart"
levels(fw.major$Stock_Analysis)[levels(fw.major$Stock_Analysis)=="Harrison-Widgeon"] <- "Harrison"
fw.data <- fw.major #rename

## Examine output
head(sw.data)
dim(sw.data)
head(fw.data)
dim(fw.data)
all.data <- rbind(fw.data, sw.data)
data.frame(table(sw.data$Stock_Analysis))
data.frame(table(fw.data$Stock_Analysis))

## Merge with SR data
resid.sw <- merge(trnc_resid_srr, sw.data, by = c("Stock_Analysis","Year"))
resid.fw <- merge(trnc_resid_srr, fw.data, by = c("Stock_Analysis","Year"))
resid.all <- rbind(resid.sw,resid.fw)

#Total agent detections all agents
agents <- names(resid.all[,32:101])
colSums(resid.all[,32:101])
all.agents <- resid.all[32:101]
nRows <- dim(all.agents)[1]
calcStats <- function(x) {
  pos <- sum(all.agents[,x] > 0, na.rm=TRUE)
  c("positives" = pos, "proportion" = pos / nRows)
}
result <- t(as.data.frame(Map(calcStats, colnames(all.agents))))
result
#write.csv(result[order(result[,1]),], file="data/ONNE_all_agent_detections_210519.csv")

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
colnames(resid.N.all.wide) <- c("Stock","FW","SW","FW","SW","FW","SW","FW","SW","FW","SW","FW","SW","FW","SW","FW","SW","FW","SW","FW","SW")

## create a table of stock and year totals by location
resid.N.all.wide[2,4,6,8,10,12,14,16,18,20] <- cell_spec(resid.N.all.wide[2,4,6,8,10,12,14,16,18,20], bold = T)
  
kbl(resid.N.all.wide, "html", align="c") %>%
  kable_paper(full_width = F) %>%
  collapse_rows(columns = 1, valign = "top") %>%
  column_spec (c(1,3,5,7,9,11,13,15,17,19,21), border_right = "0.5px solid gray") %>%
  kable_styling(font_size = 10, full_width=F, "striped") %>%
  add_header_above(c(" ", "2009" = 2, 
                     "2010" = 2, 
                     "2011" = 2, 
                     "2012" = 2, 
                     "2013" = 2, 
                     "2014" = 2, 
                     "2015" = 2, 
                     "2016" = 2,
                     "2017" = 2,
                     "2018" = 2)) %>%
  save_kable("test.html")


kbl(resid.N.all.wide, "html", align="c") %>%
  kable_paper(full_width = F) %>%
  collapse_rows(columns = 1, valign = "top") %>%
  column_spec (c(1,3,5,7,9,11,13,15,17), border_right = "0.5px solid gray") %>%
  kable_styling(font_size = 11, full_width=F, "striped") %>%
  add_header_above(c(" ", "2009" = 2, 
                     "2010" = 2, 
                     "2011" = 2, 
                     "2012" = 2, 
                     "2013" = 2, 
                     "2014" = 2, 
                     "2015" = 2, 
                     "2016" = 2)) %>%
  save_kable("test.html")



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
#sumrow <- as.data.frame(lapply(N.resid.fw[,Freq], func))
#sumrow
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

jpeg(filename='figs/Fig_Total samples by stock_SWFW_single plot.jpg', 
     width=400, height=400, quality=200)
ggplot(data=N.resid.all, aes(x=reorder(Var1, Freq), y=Freq, fill=source)) +
  geom_col(col="black") +
  coord_flip() +
  xlab("Stock")+
  ylab("Total samples") +
  labs(fill="Collection location")+
  scale_fill_manual(values=c("black","gray"))+
  theme(legend.position = c(0.7, 0.3))
dev.off()



jpeg(filename='figs/Fig_Total samples by stock_SWandFW.jpg', 
     width=480, height=400, quality=75)
grid.arrange(arrangeGrob(sw.N.st + theme(legend.position="none"), 
            fw.N.st + theme(legend.position="none"), 
            ncol = 2,
            left = textGrob("Stock", rot = 90, vjust = 1)),
            bottom = textGrob("Total samples", vjust = 1))
dev.off()

all.resid.data <- merge(trnc_resid_srr, all.data, by = c("Stock_Analysis","Year"))
samplesstyr.fw <- resid.fw %>% 
  group_by(Year, Stock_Analysis, SWFW) %>%
  count(SWFW)

names(all.resid.data)
samplesstyr.all<-all.resid.data %>% 
  group_by(Year, Stock_Analysis, SWFW) %>%
  count(Stock_Analysis)


## Total assays run per sample
### SW
sw.N.assays <- data.frame(inf_agt_resid_data %>% 
                            group_by(plot.agent) %>%
                            summarize(assays = sum(N)))
### FW
fw.N.assays <- data.frame(inf_agt_resid_data_fw %>% 
                            group_by(agent) %>%
                            summarize(assays = sum(N)))
dim(inf_agt_resid_data)




####All data in one table - Stock-specific metric - SW only
spsu.stock <-sw.data %>% group_by(Stock_Analysis, Year) #create object to be summarized by year

##ASSESS INFECTION PARAMETERS BY AGENT:
##ae_sal
stock.arena2.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(arena2!="NA")), #samples
        length(which(arena2>0)), #positive detections
        length(which(arena2>0))/length(which(!is.na(arena2))),  #calculates prevalence
        mean(arena2[arena2!=0], na.rm=TRUE),
        (length(which(arena2>0)) / length(which(!is.na(arena2)))) * mean(arena2[arena2!=0], na.rm=TRUE)
      )
  )
names(stock.arena2.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##c_b_cys
stock.c_b_cys.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(c_b_cys!="NA")), #samples
        length(which(c_b_cys>0)), #positive detections
        length(which(c_b_cys>0))/length(which(!is.na(c_b_cys))),  #calculates prevalence
        mean(c_b_cys[c_b_cys!=0], na.rm=TRUE),
        (length(which(c_b_cys>0)) / length(which(!is.na(c_b_cys)))) * mean(ae_hyd[c_b_cys!=0], na.rm=TRUE)
      )
  )
names(stock.c_b_cys.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: ce_sha
stock.ce_sha.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(ce_sha!="NA")), #samples
        length(which(ce_sha>0)), #positive detections
        length(which(ce_sha>0))/length(which(!is.na(ce_sha))),  #calculates prevalence
        mean(ce_sha[ce_sha!=0], na.rm=TRUE),
        (length(which(ce_sha>0)) / length(which(!is.na(ce_sha)))) * mean(ce_sha[ce_sha!=0], na.rm=TRUE)
      )
  )
names(stock.ce_sha.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: cr_sal
stock.cr_sal.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(cr_sal!="NA")), #samples
        length(which(cr_sal>0)), #positive detections
        length(which(cr_sal>0))/length(which(!is.na(cr_sal))),  #calculates prevalence
        mean(cr_sal[cr_sal!=0], na.rm=TRUE),
        (length(which(cr_sal>0)) / length(which(!is.na(cr_sal)))) * mean(cr_sal[cr_sal!=0], na.rm=TRUE)
      )
  )
names(stock.cr_sal.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: de_sal
stock.de_sal.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(de_sal!="NA")), #samples
        length(which(de_sal>0)), #positive detections
        length(which(de_sal>0))/length(which(!is.na(de_sal))),  #calculates prevalence
        mean(de_sal[de_sal!=0], na.rm=TRUE),
        (length(which(de_sal>0)) / length(which(!is.na(de_sal)))) * mean(de_sal[de_sal!=0], na.rm=TRUE)
      )
  )
names(stock.de_sal.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: fl_psy
stock.fl_psy.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(fl_psy!="NA")), #samples
        length(which(fl_psy>0)), #positive detections
        length(which(fl_psy>0))/length(which(!is.na(fl_psy))),  #calculates prevalence
        mean(fl_psy[fl_psy!=0], na.rm=TRUE),
        (length(which(fl_psy>0)) / length(which(!is.na(fl_psy)))) * mean(fl_psy[fl_psy!=0], na.rm=TRUE)
      )
  )
names(stock.fl_psy.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: fa_mar
stock.fa_mar.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(fa_mar!="NA")), #samples
        length(which(fa_mar>0)), #positive detections
        length(which(fa_mar>0))/length(which(!is.na(fa_mar))),  #calculates prevalence
        mean(fa_mar[fa_mar!=0], na.rm=TRUE),
        (length(which(fa_mar>0)) / length(which(!is.na(fa_mar)))) * mean(fa_mar[fa_mar!=0], na.rm=TRUE)
      )
  )
names(stock.fa_mar.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: ic_hof
stock.ic_hof.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(ic_hof!="NA")), #samples
        length(which(ic_hof>0)), #positive detections
        length(which(ic_hof>0))/length(which(!is.na(ic_hof))),  #calculates prevalence
        mean(ic_hof[ic_hof!=0], na.rm=TRUE),
        (length(which(ic_hof>0)) / length(which(!is.na(ic_hof)))) * mean(ic_hof[ic_hof!=0], na.rm=TRUE)
      )
  )
names(stock.ic_hof.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: ic_mul
stock.ic_mul.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(ic_mul!="NA")), #samples
        length(which(ic_mul>0)), #positive detections
        length(which(ic_mul>0))/length(which(!is.na(ic_mul))),  #calculates prevalence
        mean(ic_mul[ic_mul!=0], na.rm=TRUE),
        (length(which(ic_mul>0)) / length(which(!is.na(ic_mul)))) * mean(ic_mul[ic_mul!=0], na.rm=TRUE)
      )
  )
names(stock.ic_mul.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: ihnv
stock.ihnv.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(ihnv!="NA")), #samples
        length(which(ihnv>0)), #positive detections
        length(which(ihnv>0))/length(which(!is.na(ihnv))),  #calculates prevalence
        mean(ihnv[ihnv!=0], na.rm=TRUE),
        (length(which(ihnv>0)) / length(which(!is.na(ihnv)))) * mean(ihnv[ihnv!=0], na.rm=TRUE)
      )
  )
names(stock.ihnv.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: ku_thy
stock.ku_thy.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(ku_thy!="NA")), #samples
        length(which(ku_thy>0)), #positive detections
        length(which(ku_thy>0))/length(which(!is.na(ku_thy))),  #calculates prevalence
        mean(ku_thy[ku_thy!=0], na.rm=TRUE),
        (length(which(ku_thy>0)) / length(which(!is.na(ku_thy)))) * mean(ku_thy[ku_thy!=0], na.rm=TRUE)
      )
  )
names(stock.ku_thy.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: lo_sal
stock.lo_sal.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(lo_sal!="NA")), #samples
        length(which(lo_sal>0)), #positive detections
        length(which(lo_sal>0))/length(which(!is.na(lo_sal))),  #calculates prevalence
        mean(lo_sal[lo_sal!=0], na.rm=TRUE),
        (length(which(lo_sal>0)) / length(which(!is.na(lo_sal)))) * mean(lo_sal[lo_sal!=0], na.rm=TRUE)
      )
  )
names(stock.lo_sal.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: my_arc
stock.my_arc.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(my_arc!="NA")), #samples
        length(which(my_arc>0)), #positive detections
        length(which(my_arc>0))/length(which(!is.na(my_arc))),  #calculates prevalence
        mean(my_arc[my_arc!=0], na.rm=TRUE),
        (length(which(my_arc>0)) / length(which(!is.na(my_arc)))) * mean(my_arc[my_arc!=0], na.rm=TRUE)
      )
  )
names(stock.my_arc.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: pa_kab
stock.pa_kab.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(pa_kab!="NA")), #samples
        length(which(pa_kab>0)), #positive detections
        length(which(pa_kab>0))/length(which(!is.na(pa_kab))),  #calculates prevalence
        mean(pa_kab[pa_kab!=0], na.rm=TRUE),
        (length(which(pa_kab>0)) / length(which(!is.na(pa_kab)))) * mean(pa_kab[pa_kab!=0], na.rm=TRUE)
      )
  )
names(stock.pa_kab.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: pa_min
stock.pa_min.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(pa_min!="NA")), #samples
        length(which(pa_min>0)), #positive detections
        length(which(pa_min>0))/length(which(!is.na(pa_min))),  #calculates prevalence
        mean(pa_min[pa_min!=0], na.rm=TRUE),
        (length(which(pa_min>0)) / length(which(!is.na(pa_min)))) * mean(pa_min[pa_min!=0], na.rm=TRUE)
      )
  )
names(stock.pa_min.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: pa_pse
stock.pa_pse.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(pa_pse!="NA")), #samples
        length(which(pa_pse>0)), #positive detections
        length(which(pa_pse>0))/length(which(!is.na(pa_pse))),  #calculates prevalence
        mean(pa_pse[pa_pse!=0], na.rm=TRUE),
        (length(which(pa_pse>0)) / length(which(!is.na(pa_pse)))) * mean(pa_pse[pa_pse!=0], na.rm=TRUE)
      )
  )
names(stock.pa_pse.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: pa_ther
stock.pa_ther.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(pa_ther!="NA")), #samples
        length(which(pa_ther>0)), #positive detections
        length(which(pa_ther>0))/length(which(!is.na(pa_ther))),  #calculates prevalence
        mean(pa_ther[pa_ther!=0], na.rm=TRUE),
        (length(which(pa_ther>0)) / length(which(!is.na(pa_ther)))) * mean(pa_ther[pa_ther!=0], na.rm=TRUE)
      )
  )
names(stock.pa_ther.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: pch_sal
stock.pch_sal.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(pch_sal!="NA")), #samples
        length(which(pch_sal>0)), #positive detections
        length(which(pch_sal>0))/length(which(!is.na(pch_sal))),  #calculates prevalence
        mean(pch_sal[pch_sal!=0], na.rm=TRUE),
        (length(which(pch_sal>0)) / length(which(!is.na(pch_sal)))) * mean(pch_sal[pch_sal!=0], na.rm=TRUE)
      )
  )
names(stock.pch_sal.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: prv
stock.prv.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(prv!="NA")), #samples
        length(which(prv>0)), #positive detections
        length(which(prv>0))/length(which(!is.na(prv))),  #calculates prevalence
        mean(prv[prv!=0], na.rm=TRUE),
        (length(which(prv>0)) / length(which(!is.na(prv)))) * mean(prv[prv!=0], na.rm=TRUE)
      )
  )
names(stock.prv.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: pspv
stock.pspv.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(pspv!="NA")), #samples
        length(which(pspv>0)), #positive detections
        length(which(pspv>0))/length(which(!is.na(pspv))),  #calculates prevalence
        mean(pspv[pspv!=0], na.rm=TRUE),
        (length(which(pspv>0)) / length(which(!is.na(pspv)))) * mean(pspv[pspv!=0], na.rm=TRUE)
      )
  )
names(stock.pspv.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: rlo
stock.rlo.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(rlo!="NA")), #samples
        length(which(rlo>0)), #positive detections
        length(which(rlo>0))/length(which(!is.na(rlo))),  #calculates prevalence
        mean(rlo[rlo!=0], na.rm=TRUE),
        (length(which(rlo>0)) / length(which(!is.na(rlo)))) * mean(rlo[rlo!=0], na.rm=TRUE)
      )
  )
names(stock.rlo.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: sch
stock.sch.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(sch!="NA")), #samples
        length(which(sch>0)), #positive detections
        length(which(sch>0))/length(which(!is.na(sch))),  #calculates prevalence
        mean(sch[sch!=0], na.rm=TRUE),
        (length(which(sch>0)) / length(which(!is.na(sch)))) * mean(sch[sch!=0], na.rm=TRUE)
      )
  )
names(stock.sch.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: smallUK
stock.smallUK.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(smallUK!="NA")), #samples
        length(which(smallUK>0)), #positive detections
        length(which(smallUK>0))/length(which(!is.na(smallUK))),  #calculates prevalence
        mean(smallUK[smallUK!=0], na.rm=TRUE),
        (length(which(smallUK>0)) / length(which(!is.na(smallUK)))) * mean(smallUK[smallUK!=0], na.rm=TRUE)
      )
  )
names(stock.smallUK.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: sp_des
stock.sp_des.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(sp_des!="NA")), #samples
        length(which(sp_des>0)), #positive detections
        length(which(sp_des>0))/length(which(!is.na(sp_des))),  #calculates prevalence
        mean(sp_des[sp_des!=0], na.rm=TRUE),
        (length(which(sp_des>0)) / length(which(!is.na(sp_des)))) * mean(sp_des[sp_des!=0], na.rm=TRUE)
      )
  )
names(stock.sp_des.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: te_bry
stock.te_bry.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(te_bry!="NA")), #samples
        length(which(te_bry>0)), #positive detections
        length(which(te_bry>0))/length(which(!is.na(te_bry))),  #calculates prevalence
        mean(te_bry[te_bry!=0], na.rm=TRUE),
        (length(which(te_bry>0)) / length(which(!is.na(te_bry)))) * mean(te_bry[te_bry!=0], na.rm=TRUE)
      )
  )
names(stock.te_bry.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

##For now, by agent: te_mar
stock.te_mar.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(te_mar!="NA")), #samples
        length(which(te_mar>0)), #positive detections
        length(which(te_mar>0))/length(which(!is.na(te_mar))),  #calculates prevalence
        mean(te_mar[te_mar!=0], na.rm=TRUE),
        (length(which(te_mar>0)) / length(which(!is.na(te_mar)))) * mean(te_mar[te_mar!=0], na.rm=TRUE)
      )
  )
names(stock.te_mar.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

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
        
      )
  )
names(stock.ven.sw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

## Use loop to create individual agent DFs
library(tidyverse)
for(i in agents) {
  nam <- paste("111stock.", i, ".sw", sep = "")
  temp <- data.frame(
      spsu.stock %>% 
        group_by(Stock_Analysis, Year) %>%
        summarise(
          N = length(which(i!="NA")), #samples
          posdet = length(which(i>0)), #positive detections
          prev = length(which(i>0))/length(which(!is.na(i))),  #calculates prevalence
          meanload = mean(i[i!=0], na.rm=TRUE),
          agent = i,
          groups(Stock_Analysis, Year)
        )
    )
  assign(nam,temp)
}

iris.summary <- spsu.stock %>% 
  gather(variable, value, -Stock_Analysis, -Year) %>% 
  group_by(variable, Stock_Analysis, Year) %>% 
  summarize(
    N = length(variable), #samples
    posdet = length(which(variable != 0)), #positive detections
    prev = length(which(variable>0))/length(which(variable!="NA")),  #calculates prevalence
    meanload = mean(variable[variable!=0], na.rm=TRUE)
  )


categorical <- sw.data[,28:97] 
#stx <- levels(spsu.stock$Stock_Analysis)
#yrs <- unique(spsu.stock$Year)
library(dplyr)
for (i in colnames(categorical)) {
  nam <- paste("111stock.", i, ".sw", sep = "") 
  temp <- sw.data %>% 
    gather(Variable, Value, -Stock_Analysis, -Year) %>% 
    group_by(Stock_Analysis, Year) %>% 
      summarise(
        N=length(na.omit(categorical[[i]])), #samples
        postdet=sum(na.omit(categorical[[i]] > 0)), #positive detections
        meanload = mean(categorical[[i]]!=0, na.rm=TRUE),
        prev=sum(na.omit(categorical[[i]] > 0))/length(na.omit(categorical[[i]]))) 
  assign(nam,temp)
}

variable <- colnames(spsu.stock[,28:97]) 
purrr::map(variable, ~ spsu.stock %>%
             group_by_at(.x) %>%
             summarise(number = mean(mpg))) %>%
  set_names(variable) %>%
  bind_rows(., .id = 'variable')


temp <- spsu.stock %>% 
  group_by(Stock_Analysis, Year) %>% 
  summarise(
    N=length(na.omit(c_b_cys)), #samples
    postdet=sum(na.omit(c_b_cys > 0)), #positive detections
    meanload = mean(c_b_cys[c_b_cys!=0], na.rm=TRUE),
    prev=sum(na.omit(c_b_cys > 0))/length(na.omit(c_b_cys)))




## Merge raw pathogen data into one file

full.sw <- rbind(
  stock.arena2.sw, stock.c_b_cys.sw, stock.ce_sha.sw, stock.cr_sal.sw, stock.de_sal.sw,
  stock.fa_mar.sw, stock.fl_psy.sw, stock.ic_hof.sw, stock.ic_mul.sw, stock.ihnv.sw,
  stock.ku_thy.sw, stock.lo_sal.sw, stock.my_arc.sw, stock.pa_kab.sw,
  stock.pa_min.sw, stock.pa_pse.sw, stock.pa_ther.sw, stock.pspv.sw,
  stock.prv.sw, stock.rlo.sw, stock.sch.sw, stock.smallUK.sw, 
  stock.sp_des.sw, stock.te_bry.sw, stock.te_mar.sw, stock.ven.sw)

posi.full.sw <- aggregate(posdet ~ agent + Year, full.sw, sum)

#merge dataframes
arena2.resid.sws <- merge(trnc_resid_srr, stock.arena2.sw, by = c("Stock_Analysis","Year"))
c_b_cys.resid.sws <- merge(trnc_resid_srr, stock.c_b_cys.sw, by = c("Stock_Analysis","Year"))
ce_sha.resid.sws <- merge(trnc_resid_srr, stock.ce_sha.sw, by = c("Stock_Analysis","Year"))
cr_sal.resid.sws <- merge(trnc_resid_srr, stock.cr_sal.sw, by = c("Stock_Analysis","Year"))
de_sal.resid.sws <- merge(trnc_resid_srr, stock.de_sal.sw, by = c("Stock_Analysis","Year"))
fa_mar.resid.sws <- merge(trnc_resid_srr, stock.fa_mar.sw, by = c("Stock_Analysis","Year"))
fl_psy.resid.sws <- merge(trnc_resid_srr, stock.fl_psy.sw, by = c("Stock_Analysis","Year"))
ic_hof.resid.sws <- merge(trnc_resid_srr, stock.ic_hof.sw, by = c("Stock_Analysis","Year"))
ic_mul.resid.sws <- merge(trnc_resid_srr, stock.ic_mul.sw, by = c("Stock_Analysis","Year"))
ihnv.resid.sws <- merge(trnc_resid_srr, stock.ihnv.sw, by = c("Stock_Analysis","Year"))
ku_thy.resid.sws <- merge(trnc_resid_srr, stock.ku_thy.sw, by = c("Stock_Analysis","Year"))
lo_sal.resid.sws <- merge(trnc_resid_srr, stock.lo_sal.sw, by = c("Stock_Analysis","Year"))
my_arc.resid.sws <- merge(trnc_resid_srr, stock.my_arc.sw, by = c("Stock_Analysis","Year"))
pa_kab.resid.sws <- merge(trnc_resid_srr, stock.pa_kab.sw, by = c("Stock_Analysis","Year"))
pa_min.resid.sws <- merge(trnc_resid_srr, stock.pa_min.sw, by = c("Stock_Analysis","Year"))
pa_pse.resid.sws <- merge(trnc_resid_srr, stock.pa_pse.sw, by = c("Stock_Analysis","Year"))
pa_ther.resid.sws <- merge(trnc_resid_srr, stock.pa_ther.sw, by = c("Stock_Analysis","Year"))
prv.resid.sws <- merge(trnc_resid_srr, stock.prv.sw, by = c("Stock_Analysis","Year"))
pspv.resid.sws <- merge(trnc_resid_srr, stock.pspv.sw, by = c("Stock_Analysis","Year"))
rlo.resid.sws <- merge(trnc_resid_srr, stock.rlo.sw, by = c("Stock_Analysis","Year"))
sch.resid.sws <- merge(trnc_resid_srr, stock.sch.sw, by = c("Stock_Analysis","Year"))
smallUK.resid.sws <- merge(trnc_resid_srr, stock.smallUK.sw, by = c("Stock_Analysis","Year"))
sp_des.resid.sws <- merge(trnc_resid_srr, stock.sp_des.sw, by = c("Stock_Analysis","Year"))
te_bry.resid.sws <- merge(trnc_resid_srr, stock.te_bry.sw, by = c("Stock_Analysis","Year"))
te_mar.resid.sws <- merge(trnc_resid_srr, stock.te_mar.sw, by = c("Stock_Analysis","Year"))
ven.resid.sws <- merge(trnc_resid_srr, stock.ven.sw, by = c("Stock_Analysis","Year"))

##ADD COLUMN WITH AGENT NAME TO RBIND INTOLARGE FILE
arena2.resid.sws$agent<-"arena2"
c_b_cys.resid.sws$agent<-"c_b_cys"
ce_sha.resid.sws$agent<-"ce_sha"
cr_sal.resid.sws$agent<-"cr_sal"
de_sal.resid.sws$agent<-"de_sal"
fa_mar.resid.sws$agent<-"fa_mar"
fl_psy.resid.sws$agent<-"fl_psy"
ic_hof.resid.sws$agent<-"ic_hof"
ic_mul.resid.sws$agent<-"ic_mul"
ihnv.resid.sws$agent<- "ihnv"
ku_thy.resid.sws$agent<-"ku_thy"
lo_sal.resid.sws$agent<-"lo_sal"
my_arc.resid.sws$agent<-"my_arc"
pa_kab.resid.sws$agent<-"pa_kab"
pa_min.resid.sws$agent<-"pa_min"
pa_pse.resid.sws$agent<-"pa_pse"
pa_ther.resid.sws$agent<-"pa_ther"
prv.resid.sws$agent<-"prv"
pspv.resid.sws$agent<-"pspv"
rlo.resid.sws$agent<-"rlo"
sch.resid.sws$agent<-"sch"
smallUK.resid.sws$agent<-"smallUK"
sp_des.resid.sws$agent<-"sp_des"
te_bry.resid.sws$agent<-"te_bry"
te_mar.resid.sws$agent<-"te_mar"
ven.resid.sws$agent<-"ven"

##MERGE INTO LARGE FILE - only those w/ >30 detections stock-specific metrics
full.resid.sws <- rbind(
  arena2.resid.sws, c_b_cys.resid.sws, ce_sha.resid.sws, cr_sal.resid.sws, de_sal.resid.sws,
  fa_mar.resid.sws, fl_psy.resid.sws, ic_hof.resid.sws, ic_mul.resid.sws, ihnv.resid.sws,
  ku_thy.resid.sws, lo_sal.resid.sws, my_arc.resid.sws, pa_kab.resid.sws,
  pa_min.resid.sws, pa_pse.resid.sws, pa_ther.resid.sws, pspv.resid.sws,
  prv.resid.sws, rlo.resid.sws, sch.resid.sws, smallUK.resid.sws, 
  sp_des.resid.sws, te_bry.resid.sws, te_mar.resid.sws, ven.resid.sws)
write.csv(full.resid.sws, file="data/FULL_ONNE productivity infection analysis_stockspecific_210519.csv")

## How many positive detections?
posi.det.total <- aggregate(posdet ~ agent, full.resid.sws, sum)
posi.det.annual <- aggregate(posdet ~ Year + agent, full.resid.sws, sum)
ggplot(data = posi.det.annual, aes(x=Year, y=posdet, group=agent)) +
  geom_area(aes(fill=agent))


## Reduce data set to those agents with >30 detections
posi.det.total[order(posi.det.total$posdet),]
### smallUK 10, PRV 5, IHNV 8, CR_SAL 2 ---> remove from analysis due to insufficient data
## Going to leave them in for now - can remove later
inf_agt_resid_data <- full.resid.sws
#[full.resid.sws$agent !=  "smallUK" 
#                                     & full.resid.sws$agent != "ihnv"
#                                     & full.resid.sws$agent != "cr_sal",]
write.csv(inf_agt_resid_data, file="data/ONNE productivity infection analysis_stockspecific_210519.csv")

## Summarize agents
### SW
posi.inf.annual <- aggregate(posdet ~ Year + agent, inf_agt_resid_data, sum)
posi.inf <- aggregate(posdet ~ agent, inf_agt_resid_data, sum)
posi.inf[order(posi.inf$posdet),]







##################################################
### FW data cleaning

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

##For now, by agent: fa_mar
stock.fa_mar.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(fa_mar!="NA")), #samples
        length(which(fa_mar>0)), #positive detections
        length(which(fa_mar>0))/length(which(!is.na(fa_mar))),  #calculates prevalence
        mean(fa_mar[fa_mar!=0], na.rm=TRUE),
        (length(which(fa_mar>0)) / length(which(!is.na(fa_mar)))) * mean(fa_mar[fa_mar!=0], na.rm=TRUE)
      )
  )
names(stock.fa_mar.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

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

##For now, by agent: ku_thy
stock.ku_thy.fw =
  data.frame(
    fw.stock %>% 
      summarise(
        length(which(ku_thy!="NA")), #samples
        length(which(ku_thy>0)), #positive detections
        length(which(ku_thy>0))/length(which(!is.na(ku_thy))),  #calculates prevalence
        mean(ku_thy[ku_thy!=0], na.rm=TRUE),
        (length(which(ku_thy>0)) / length(which(!is.na(ku_thy)))) * mean(ku_thy[ku_thy!=0], na.rm=TRUE)
      )
  )
names(stock.ku_thy.fw) <- c("Stock_Analysis", "Year", "N", "posdet", "prev", "mean_load", "prevload") #rename columns

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


#merge dataframes
arena2.resid.fws <- merge(trnc_resid_srr, stock.arena2.fw, by = c("Stock_Analysis","Year"))
c_b_cys.resid.fws <- merge(trnc_resid_srr, stock.c_b_cys.fw, by = c("Stock_Analysis","Year"))
ce_sha.resid.fws <- merge(trnc_resid_srr, stock.ce_sha.fw, by = c("Stock_Analysis","Year"))
cr_sal.resid.fws <- merge(trnc_resid_srr, stock.cr_sal.fw, by = c("Stock_Analysis","Year"))
de_sal.resid.fws <- merge(trnc_resid_srr, stock.de_sal.fw, by = c("Stock_Analysis","Year"))
fa_mar.resid.fws <- merge(trnc_resid_srr, stock.fa_mar.fw, by = c("Stock_Analysis","Year"))
fl_psy.resid.fws <- merge(trnc_resid_srr, stock.fl_psy.fw, by = c("Stock_Analysis","Year"))
ic_hof.resid.fws <- merge(trnc_resid_srr, stock.ic_hof.fw, by = c("Stock_Analysis","Year"))
ic_mul.resid.fws <- merge(trnc_resid_srr, stock.ic_mul.fw, by = c("Stock_Analysis","Year"))
ihnv.resid.fws <- merge(trnc_resid_srr, stock.ihnv.fw, by = c("Stock_Analysis","Year"))
ku_thy.resid.fws <- merge(trnc_resid_srr, stock.ku_thy.fw, by = c("Stock_Analysis","Year"))
lo_sal.resid.fws <- merge(trnc_resid_srr, stock.lo_sal.fw, by = c("Stock_Analysis","Year"))
my_arc.resid.fws <- merge(trnc_resid_srr, stock.my_arc.fw, by = c("Stock_Analysis","Year"))
pa_kab.resid.fws <- merge(trnc_resid_srr, stock.pa_kab.fw, by = c("Stock_Analysis","Year"))
pa_min.resid.fws <- merge(trnc_resid_srr, stock.pa_min.fw, by = c("Stock_Analysis","Year"))
pa_pse.resid.fws <- merge(trnc_resid_srr, stock.pa_pse.fw, by = c("Stock_Analysis","Year"))
pa_ther.resid.fws <- merge(trnc_resid_srr, stock.pa_ther.fw, by = c("Stock_Analysis","Year"))
prv.resid.fws <- merge(trnc_resid_srr, stock.prv.fw, by = c("Stock_Analysis","Year"))
pspv.resid.fws <- merge(trnc_resid_srr, stock.pspv.fw, by = c("Stock_Analysis","Year"))
rlo.resid.fws <- merge(trnc_resid_srr, stock.rlo.fw, by = c("Stock_Analysis","Year"))
sch.resid.fws <- merge(trnc_resid_srr, stock.sch.fw, by = c("Stock_Analysis","Year"))
smallUK.resid.fws <- merge(trnc_resid_srr, stock.smallUK.fw, by = c("Stock_Analysis","Year"))
sp_des.resid.fws <- merge(trnc_resid_srr, stock.sp_des.fw, by = c("Stock_Analysis","Year"))
te_bry.resid.fws <- merge(trnc_resid_srr, stock.te_bry.fw, by = c("Stock_Analysis","Year"))
te_mar.resid.fws <- merge(trnc_resid_srr, stock.te_mar.fw, by = c("Stock_Analysis","Year"))
ven.resid.fws <- merge(trnc_resid_srr, stock.ven.fw, by = c("Stock_Analysis","Year"))

##ADD COLUMN WITH AGENT NAME TO RBIND INTOLARGE FILE
arena2.resid.fws$agent<-"arena2"
c_b_cys.resid.fws$agent<-"c_b_cys"
ce_sha.resid.fws$agent<-"ce_sha"
cr_sal.resid.fws$agent<-"cr_sal"
de_sal.resid.fws$agent<-"de_sal"
fa_mar.resid.fws$agent<-"fa_mar"
fl_psy.resid.fws$agent<-"fl_psy"
ic_hof.resid.fws$agent<-"ic_hof"
ic_mul.resid.fws$agent<-"ic_mul"
ihnv.resid.fws$agent<- "ihnv"
ku_thy.resid.fws$agent<-"ku_thy"
lo_sal.resid.fws$agent<-"lo_sal"
my_arc.resid.fws$agent<-"my_arc"
pa_kab.resid.fws$agent<-"pa_kab"
pa_min.resid.fws$agent<-"pa_min"
pa_pse.resid.fws$agent<-"pa_pse"
pa_ther.resid.fws$agent<-"pa_ther"
prv.resid.fws$agent<-"prv"
pspv.resid.fws$agent<-"pspv"
rlo.resid.fws$agent<-"rlo"
sch.resid.fws$agent<-"sch"
smallUK.resid.fws$agent<-"smallUK"
sp_des.resid.fws$agent<-"sp_des"
te_bry.resid.fws$agent<-"te_bry"
te_mar.resid.fws$agent<-"te_mar"
ven.resid.fws$agent<-"ven"

##MERGE INTO LARGE FILE - only those w/ >30 detections stock-specific metrics
full.resid.fws <- rbind(
  arena2.resid.fws, c_b_cys.resid.fws, ce_sha.resid.fws, cr_sal.resid.fws, de_sal.resid.fws,
  fa_mar.resid.fws, fl_psy.resid.fws, ic_hof.resid.fws, ic_mul.resid.fws, ihnv.resid.fws,
  ku_thy.resid.fws, lo_sal.resid.fws, my_arc.resid.fws, pa_kab.resid.fws,
  pa_min.resid.fws, pa_pse.resid.fws, pa_ther.resid.fws, pspv.resid.fws,
  prv.resid.fws, rlo.resid.fws, sch.resid.fws, smallUK.resid.fws, 
  sp_des.resid.fws, te_bry.resid.fws, te_mar.resid.fws, ven.resid.fws)
write.csv(full.resid.fws, file="data/FULL_ONNE productivity infection analysis_stockspecific_FW_210521.csv")

## How many positive detections?
posi.det.total.fw <- aggregate(posdet ~ agent, full.resid.fws, sum)
posi.det.annual.fw <- aggregate(posdet ~ Year + agent, full.resid.fws, sum)
ggplot(data = posi.det.annual.fw, aes(x=Year, y=posdet, group=agent)) +
  geom_area(aes(fill=agent))


## Reduce data set to those agents with >30 detections
posi.det.total.fw[order(posi.det.total.fw$posdet),]
###  agent posdet
#ku_thy      0
#cr_sal      1
#fa_mar      1
#prv      2
#ven      2
#sch      4
#te_mar      4
#ihnv      6
#sp_des      6
#arena2      8
#smallUK      8
#lo_sal     10
#pa_pse     10
#pa_kab     16
#pa_ther     17
#de_sal     19---> remove from analysis due to insufficient data <20 detections overall

## 
## Keep ALL for now for analysis - LW
inf_agt_resid_data_fw <- full.resid.fws
#[full.resid.fws$agent !=  "fa_mar" 
#& full.resid.fws$agent != "ku_thy"
 #                                       & full.resid.fws$agent != "cr_sal"
  #                                      & full.resid.fws$agent != "smallUK"
   #                                     & full.resid.fws$agent != "prv"
    #                                    & full.resid.fws$agent != "ven"
     #                                   & full.resid.fws$agent != "arena2"
      #                                  & full.resid.fws$agent != "sch"
       #                                 & full.resid.fws$agent != "te_mar"
        #                                & full.resid.fws$agent != "ihnv"
         #                               & full.resid.fws$agent != "sp_des"
          #                              & full.resid.fws$agent != "pa_pse"
           #                             & full.resid.fws$agent != "lo_sal"
            #                            & full.resid.fws$agent != "pa_kab"
             #                           & full.resid.fws$agent != "de_sal"
              #                          & full.resid.fws$agent != "pa_ther"
               #                         & full.resid.fws$agent != "ic_hof",]
inf_agt_resid_data_fw <- inf_agt_resid_data_fw[inf_agt_resid_data_fw$Stock_Analysis !=  "Harrison",]
write.csv(inf_agt_resid_data_fw, file="data/ONNE productivity infection analysis_stockspecific_FW_210519.csv")


## Summarize agents and total fish sampled
### FW
posi.inf.annual.fw <- aggregate(posdet ~ Year + agent, inf_agt_resid_data_fw, sum)
posi.inf.fw <- aggregate(posdet ~ agent, inf_agt_resid_data_fw, sum)
posi.inf.fw[order(posi.inf.fw$posdet),]
levels(inf_agt_resid_data$Stock_Analysis)


posi.inf.fw <- aggregate(posdet ~ agent, inf_agt_resid_data_fw, sum)

## REVIEW
head(inf_agt_resid_data)
head(inf_agt_resid_data_fw)
sum(inf_agt_resid_data$N)

head(sw.data)
head(fw.data)
length(sw.data$Unique)

## Tables
## How many fish per year/stock/source?
ts1 <- resid.all %>%
  group_by(Year, Stock_Analysis, SWFW) %>%
  count(SWFW) 
ts2 <- data.frame(ts1)
ts3 <- xtabs(n ~ Stock_Analysis + SWFW + Year, data=ts2)
#ts4 <- ts3[rowSums(ts3[,-1]) > 0, ]
ts5 <- addmargins(ts3)
ts6 <- ftable(ts5)
cont <- stats:::format.ftable(ts6, quote = FALSE)
#write.table(cont, sep = ",", file = "table.csv")

## Assays run per year 
library(tidyverse)
t0 <- resid.all %>%
  dplyr::select(Year, ae_hyd:ye_ruc_glnA) %>%
  group_by(Year) %>%
  summarise_each(funs(zero=sum(.==0), NAs=sum(is.na(.)), pos=sum(.>0)), ae_hyd:ye_ruc_glnA) %>%
  ungroup %>%
  as.data.frame()
t1 <- t(t0[,2:ncol(t0)])
colnames(t1) <- t0[,1] 
t2 <- data.frame(t1[order(rownames(t1)),])
t2$value <- seq(1:3)

##OK SO FROM HERE I CAN PLOT THESE AND ENSURE THAT THE VALUES TOTAL TO THE COUNT OF FISH SAMPLED






ts3 <- xtabs(n ~ Stock_Analysis + Year, data=ts2)
#ts4 <- ts3[rowSums(ts3[,-1]) > 0, ]
ts5 <- addmargins(ts3)
ts6 <- ftable(ts5)
cont <- stats:::format.ftable(ts6, quote = FALSE)
#write.table(cont, sep = ",", file = "table.csv")
