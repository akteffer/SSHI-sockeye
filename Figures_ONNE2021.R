### Figures ONNE 2021 pathogens ~ producivity/health
## A K. Teffer

### Run code from Raw Data Filtering script for sw.data and fw.data
### Also run code for SSHI STage I and II

library(lme4)
#library(rstanarm) # https://mc-stan.org/users/documentation/case-studies/tutorial_rstanarm.html
library(ggplot2)
library(plotrix)
library(tidyverse)
library(gridExtra)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
#library(shinystan)
library(data.table)
library(base)
library(ggpubr)
library(dplyr)
library(ggrepel)




###############################################################################
## BRING ALL MODELING DATA TOGETHER

# beta estimates
prev_beta_sw <- read.csv("data/SW_prev_coefs_stan_stspec_230123.csv")
prev_beta_fw <- read.csv("data/FW_prev_coefs_stan_stspec_230123.csv")
prev_beta_ch <- read.csv("data/prev_coefs_stan_ChilkoSW_RSM.csv")

# proportion less than 0
temp <- read.csv("data/SW_Percent post draws >0_prev_stspec_230123.csv")
prev_neg_sw <- temp[,1:2]
temp <- read.csv("data/FW_Percent post draws >0_prev_stspec_230123.csv")
prev_neg_fw <- data.frame(temp[,1:2])
temp <- read.csv("data/Percent post draws >0_prev_chilko_RSM.csv")
prev_neg_ch <- data.frame(temp[,1:2])

# Agent names file
agent.names <- read.csv("data/agent.names_230124.csv")
agent.names <- agent.names[,c(2:3)]

# Add extra agents to FW results with NA
vec1 <- unique(prev_beta_sw$X)
vec2 <- unique(prev_beta_fw$X)
blank.agents.fw <- setdiff(vec1, vec2) 
blank.agents.sw <- setdiff(vec2, vec1) 
tempfw <- data.frame(
  "X" = blank.agents.fw,
  "lower" = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
  "X25" = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
  "mid" = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
  "X75" = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
  "upper" = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
prev_beta_fw <- rbind(prev_beta_fw,tempfw)

## Bring in assays from Tables R script
#SW
samps.sw <- read.csv("figs/Assay_table_SW_ONNE2021_220126.csv")
samps.sw$SWFW <- "SW"
prev.sw <- samps.sw[,c(1,5)]
prev_beta_sw <- merge(prev_beta_sw, prev.sw, by="X", all.x=TRUE)
#FW
samps.fw <- read.csv("figs/Assay_table_FW_ONNE2021_220126.csv")
samps.fw$SWFW <- "FW"
prev.fw <- samps.fw[,c(1,5)]
prev_beta_fw <- merge(prev_beta_fw, prev.fw, by="X", all.x=TRUE)
#Chilko
samps.ch <- read.csv("figs/Assay_table_Chilko_ONNE2021_230222.csv")
samps.ch$SWFW <- "CH"
prev.ch <- samps.ch[,c(1,5)]
prev_beta_ch <- merge(prev_beta_ch, prev.ch, by="X", all.x=TRUE)


#proportion less than 0
#SW
prev_neg_sw <-read.csv("data/SW_Percent post draws >0_prev_stspec_230123.csv")
#FW
prev_neg_fw <- read.csv("data/FW_Percent post draws >0_prev_stspec_230123.csv")
tempfw <- data.frame(
  "X" = blank.agents.fw,
  "prop0" = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
prev_neg_fw <- rbind(prev_neg_fw,tempfw)
colnames(prev_neg_fw) <- c("X","prop0.fw")
#SW
prev_neg_ch <-read.csv("data/Percent post draws >0_prev_chilko_SWprev_SRcovresid.csv")

# Merge SW df's, FW df's, & agent names
names(prev_beta_sw)
names(prev_beta_fw)
names(prev_beta_ch)
names(prev_neg_sw)
names(prev_neg_fw)
names(prev_neg_ch)
names(agent.names)

#SW
prev_beta_sw <- merge(prev_beta_sw, prev_neg_sw, by="X")
colnames(prev_beta_sw) <- c("agent","lwr.sw","midlwr.sw","mid.sw","midupp.sw","upp.sw","SW_prev","prop0.sw")
#FW
prev_beta_fw <- merge(prev_beta_fw,prev_neg_fw, by="X")
colnames(prev_beta_fw) <- c("agent","lwr.fw","midlwr.fw","mid.fw","midupp.fw","upp.fw","FW_prev","prop0.fw")
#CHILKO
prev_beta_ch <- merge(prev_beta_ch, prev_neg_ch, by="X")
colnames(prev_beta_ch) <- c("agent","lwr.ch","midlwr.ch","mid.ch","midupp.ch","upp.ch","CH_prev","prop0.ch")
# SW + FW
prev_beta_all <- merge(prev_beta_sw,prev_beta_fw, by="agent")
#Add agent names
all <- merge(prev_beta_all,agent.names, by="agent")

## Prevalence as percent
all$SW_prev <- round(all$SW_prev*100,0)
all$FW_prev <- round(all$FW_prev*100,0)
all
prev_beta_ch$CH_prev <- round(prev_beta_ch$CH_prev*100,0)



### Mass Deviation results
# SW
MD.sw.full <- readRDS("sockeye_SW_lw_results_jan24_2023.rds")
temp <- MD.sw.full$estimates
temp2 <- as.data.frame(t(temp))
temp2$taxa <- rownames(temp2)
MD.sw.load <- temp2
MD.sw.post <- MD.sw.full$posterior_probabilities
MD.sw.post <- MD.sw.post[,c(1,3:4)]
colnames(MD.sw.post) <- c("taxa","pathogen_pp_sw","sst_pp")
MD.sw <- merge(MD.sw.load, MD.sw.post, by="taxa")
##SW: Agents over the line for B0 - all others look good (1.0-1.1) - Rhabdo3_virus, fa_mar
# Qin also has a very wide variance which likely indicates poor ESS
MD.sw<-MD.sw[!(MD.sw$taxa=="Rhabdo3_virus" | MD.sw$taxa=="fa_mar" | MD.sw$taxa=="Qin"),]

# FW - no random slope
MD.fw.full.nrs <- readRDS("sockeye_FW_norandslope_lw_results_jan26_2023.rds")
temp <- MD.fw.full.nrs$estimates
temp2 <- as.data.frame(t(temp))
temp2$taxa <- rownames(temp2)
MD.fw.load <- temp2
MD.fw.post <- MD.fw.full.nrs$posterior_probabilities
MD.fw.post <- MD.fw.post[,c(1,3)]
colnames(MD.fw.post) <- c("taxa","pathogen_pp_fw")
MD.fw <- merge(MD.fw.load, MD.fw.post, by="taxa")
##FW: ##Agents all look good (1.0-1.1) for rhats 
## Circo-virus has a wide CI, which indicates poor Neff and low sample size for Rhabdo
MD.fw<-MD.fw[!(MD.fw$taxa=="Circo.virus" | MD.fw$taxa=="Rhabdo3_virus"),]

# Save output files:
write.csv(MD.fw, file="data/MD.fw.csv")
write.csv(MD.sw, file="data/MD.sw.csv")


## Create plot that combines SW and MD effects - also FW?
#Merge md and surv datasets
all
md.all <- merge(MD.sw, MD.fw, by="taxa", all.x = T)
names(md.all)
#Define terms: "median.sw","lci.sw"=2.5%,"mlci.sw"=25%,"clci.sw"=40%,"cuci.sw"=60%,"muci.sw"=75%,"uci.sw"=97.5%
colnames(md.all) <- c("agent","mean.sw","sd.sw","median.sw","lci.sw","mlci.sw","clci.sw","cuci.sw",
                      "muci.sw","uci.sw","path_pp.sw","sst_pp_sw",
                      "mean.fw","sd.fw","median.fw","lci.fw","mlci.fw","clci.fw","cuci.fw",
                      "muci.fw","uci.fw","path_pp.fw")
md.all
md_surv <- merge(md.all, all, by="agent")
md_surv <- subset(md_surv, agent!="smallUK") #remove smallUK, not much data
write.csv(md_surv, file="data/md_surv_220126.csv")
write.csv(md_surv, file="figs/Final data frames and tables/md_surv_220126.csv")

##SW survival plot 
sw.surv.plot.df <- md_surv[order(md_surv$mid.sw),] #sort all by descending SW median beta (50th percentile)
sw.surv.plot.df2 <- subset(sw.surv.plot.df, SW_prev>1) #remove low prevalence agents (<1% prevalence)
temp <- md_surv[md_surv$agent %in% c('rib', 'richness'), ]
sw.surv.plot.df2 <- rbind(sw.surv.plot.df2, temp)
# Optional: Present in <3 years, Remove arena1, mo_vis, Picorna2_virus, re_sal, ascv, cr_sal, pch_sal, smallUK

sw.surv.plot <- sw.surv.plot.df2 %>%
  ggplot(aes(x=mid.sw, y=reorder(plot.agent,-mid.sw), alpha = SW_prev)) +
    geom_point() +
    geom_errorbarh(aes(xmax = upp.sw, xmin = lwr.sw), linewidth=.5, height = 0) +
    geom_errorbarh(aes(xmax = midupp.sw, xmin = midlwr.sw), linewidth=1, height = 0) +
    xlim(-1,1.5) +
    geom_vline(xintercept = 0, lty = 2, linewidth=.25) +
    labs(x="Effect size", y = "", title = "Marine Prevalence~Survival") +
    theme(legend.position = "none", axis.text.y = element_text(face="italic")) +
    geom_text(aes(x=1.5, y=reorder(plot.agent,-mid.sw), label = round((100*prop0.sw),0)), inherit.aes = FALSE)+
    geom_text(aes(x=-1, y=reorder(plot.agent,-mid.sw), label = SW_prev, alpha=100))

jpeg(filename='figs/Fig_SW_beta_prev_230124.jpg', width=500, height=600, quality=300)
sw.surv.plot
dev.off()

##FW survival plot 
fw.surv.plot.df <- subset(md_surv, !(is.na(lwr.fw))) #remove rows with no modeling data
fw.surv.plot.df2 <- subset(fw.surv.plot.df, FW_prev>1) #remove low prevalence agents (<1% prevalence)
fw.surv.plot.df3 <- fw.surv.plot.df2[order(fw.surv.plot.df2$mid.fw),] #sort all by descending SW median beta (50th percentile)
temp <- md_surv[md_surv$agent %in% c('rib', 'richness'), ]
fw.surv.plot.df3 <- rbind(fw.surv.plot.df3, temp)
# Optional: Present in <3 years?

fw.surv.plot <- fw.surv.plot.df3 %>%
  ggplot(aes(x=mid.fw, y=reorder(plot.agent,-mid.fw), alpha = FW_prev)) +
  geom_point() +
  geom_errorbarh(aes(xmax = upp.fw, xmin = lwr.fw), linewidth=.5, height = 0) +
  geom_errorbarh(aes(xmax = midupp.fw, xmin = midlwr.fw), linewidth=1, height = 0) +
  xlim(-1,1.5) +
  geom_vline(xintercept = 0, lty = 2, linewidth=.25) +
  labs(x="Effect size", y = "", title = "Freshwater Prevalence~Survival") +
  theme(legend.position = "none", axis.text.y = element_text(face="italic")) +
  geom_text(aes(x=1.5, y=reorder(plot.agent,-mid.fw), label = round((100*prop0.fw),0)), inherit.aes = FALSE)+
  geom_text(aes(x=-1, y=reorder(plot.agent,-mid.fw), label = FW_prev, alpha=100))

jpeg(filename='figs/Fig_FW_beta_prev_230124.jpg', width=500, height=600, quality=300)
fw.surv.plot
dev.off()

# Both survival plots
tiff(filename='figs/Fig_slope_prev_FWSW_230124.jpg', units="in", width=5, height=8, res=300)
figure <- ggarrange(fw.surv.plot,sw.surv.plot,nrow=2, heights=c(1,1.2))  
annotate_figure(figure, left = text_grob("Infectious agents", rot = 90))
dev.off()



##Marine Survival Chilko  plot - 
ch.surv.plot.df <- prev_beta_ch[order(prev_beta_ch$mid.ch),] #sort all by descending SW median beta (50th percentile)
ch.surv.plot.df2 <- subset(ch.surv.plot.df, CH_prev>1) #remove low prevalence agents (<1% prevalence)
temp <- md_surv[md_surv$agent %in% c('rib', 'richness'), ]
sw.surv.plot.df2 <- rbind(sw.surv.plot.df2, temp)
# Optional: Present in <3 years, Remove arena1, mo_vis, Picorna2_virus, re_sal, ascv, cr_sal, pch_sal, smallUK

sw.surv.plot <- sw.surv.plot.df2 %>%
  ggplot(aes(x=mid.sw, y=reorder(plot.agent,-mid.sw), alpha = SW_prev)) +
  geom_point() +
  geom_errorbarh(aes(xmax = upp.sw, xmin = lwr.sw), linewidth=.5, height = 0) +
  geom_errorbarh(aes(xmax = midupp.sw, xmin = midlwr.sw), linewidth=1, height = 0) +
  xlim(-1,1.5) +
  geom_vline(xintercept = 0, lty = 2, linewidth=.25) +
  labs(x="Effect size", y = "", title = "Marine Prevalence~Survival") +
  theme(legend.position = "none", axis.text.y = element_text(face="italic")) +
  geom_text(aes(x=1.5, y=reorder(plot.agent,-mid.sw), label = round((100*prop0.sw),0)), inherit.aes = FALSE)+
  geom_text(aes(x=-1, y=reorder(plot.agent,-mid.sw), label = SW_prev, alpha=100))

jpeg(filename='figs/Fig_SW_beta_prev_230124.jpg', width=500, height=600, quality=300)
sw.surv.plot
dev.off()



##SW MD plot 
sw.md.plot.df <- md_surv[order(md_surv$median.sw),] #sort all by descending SW median beta (50th percentile)
sw.md.plot.df2 <- subset(sw.md.plot.df, SW_prev>1) #remove low prevalence agents (<1% prevalence)
temp <- md_surv[md_surv$agent %in% c('rib', 'richness'), ]
sw.md.plot.df2 <- rbind(sw.md.plot.df2, temp)

sw.md.plot <- sw.md.plot.df2 %>%
  ggplot(aes(x=median.sw, y=reorder(plot.agent,-median.sw), alpha = SW_prev)) +
  geom_point() +
  geom_errorbarh(aes(xmax = uci.sw, xmin = lci.sw), linewidth=.5, height = 0) +
  geom_errorbarh(aes(xmax = muci.sw, xmin = mlci.sw), linewidth=1, height = 0) +
  #xlim(-1,1.5) +
  geom_vline(xintercept = 0, lty = 2, linewidth=.25) +
  labs(x="Effect size", y = "", title = "Marine Load~Mass Deviation") +
  theme(legend.position = "none", axis.text.y = element_text(face="italic")) +
  geom_text(aes(x=.2, y=reorder(plot.agent,-median.sw), label = round((100*as.numeric(path_pp.sw)),0)), inherit.aes = FALSE)+
  geom_text(aes(x=-.25, y=reorder(plot.agent,-median.sw), label = SW_prev, alpha=100))

jpeg(filename='figs/Fig_SW_beta_prev_MD_230124.jpg', width=500, height=600, quality=300)
sw.md.plot
dev.off()


##FW MD plot 
fw.md.plot.df <- subset(md_surv, !(is.na(median.fw))) #remove rows with no modeling data
fw.md.plot.df2 <- subset(fw.md.plot.df, FW_prev>1) #remove low prevalence agents (<1% prevalence)
fw.md.plot.df3 <- fw.md.plot.df2[order(fw.md.plot.df2$median.fw),] #sort all by descending SW median beta (50th percentile)
temp <- md_surv[md_surv$agent %in% c('rib', 'richness'), ]
fw.md.plot.df3 <- rbind(fw.md.plot.df3, temp)

fw.md.plot <- fw.md.plot.df3 %>%
  ggplot(aes(x=median.fw, y=reorder(plot.agent,-median.fw), alpha = FW_prev)) +
  geom_point() +
  geom_errorbarh(aes(xmax = uci.fw, xmin = lci.fw), linewidth=.5, height = 0) +
  geom_errorbarh(aes(xmax = muci.fw, xmin = mlci.fw), linewidth=1, height = 0) +
  #xlim(-1,1) +
  geom_vline(xintercept = 0, lty = 2, linewidth=.25) +
  labs(x="Effect size", y = "", title = "Freshwater Load~Mass Deviation") +
  theme(legend.position = "none", axis.text.y = element_text(face="italic")) +
  geom_text(aes(x=.2, y=reorder(plot.agent,-median.fw), label = round((100*as.numeric(path_pp.fw)),0)), inherit.aes = FALSE)+
  geom_text(aes(x=-.25, y=reorder(plot.agent,-median.fw), label = FW_prev, alpha=100))

jpeg(filename='figs/Fig_FW_beta_prev_MD_230124.jpg', width=500, height=600, quality=300)
fw.md.plot
dev.off()

# Both plots
tiff(filename='figs/Fig_slope_prev_MD_FWSW_230124.jpg', units="in", width=5, height=8, res=300)
figure <- ggarrange(fw.md.plot,sw.md.plot,nrow=2, heights=c(1,1.2))  
annotate_figure(figure, left = text_grob("Infectious agents", rot = 90))
dev.off()


## Seawater
#Remove Circo.virus and Qin - effect size distributions are too wide to have confidence
md_surv.plot<-md_surv#[!(md_surv$agent=="Circo.virus" | md_surv$agent=="Qin" | md_surv$agent=="nu_sal"),]
unique(md_surv.plot$agent)

# SW plot with agents >1% prevalence
temp1 <- subset(md_surv.plot, SW_prev>1)
temp2 <- md_surv.plot[md_surv.plot$agent %in% c('rib', 'richness'), ]
md_surv.plot.sw <- rbind(temp1, temp2)
sw_md_surv.plot2 <- 
  md_surv.plot.sw %>%
  ggplot() +
  geom_point(aes(x=mid.sw, y=median.sw, col=agent), size=.8, pch=21) +
  geom_errorbarh(aes(y=median.sw, xmin=midlwr.sw, xmax=midupp.sw, col=agent), linewidth=.25, height=0, alpha=0.7) +
  geom_errorbar(aes(x=mid.sw, y=median.sw, ymax = muci.sw, ymin = mlci.sw, col=agent), linewidth=.25, width=0, alpha=0.7) +
  ylim(-0.07,0.07) +
  xlim(-0.5,0.7) +
  annotate("rect", xmin=-Inf, xmax = 0, ymin=-Inf, ymax=Inf, alpha=.2) +
  annotate("rect", xmin=-Inf, xmax = Inf, ymin=-Inf, ymax=0, alpha=.2) +
  geom_text(aes(x=mid.sw, y=median.sw, label=agent, col=agent), size=4, show.legend = F)+
  xlab("Survival effect size")+
  ylab("Mass Deviation effect size") + 
  labs(title="Marine-sampled fish") +
  theme(legend.position = "none")

jpeg(filename='figs/Fig_MD_Surv_ONNE2021_230124_SW_CI.jpg', width=500, height=1000, quality=300)
sw_md_surv.plot2
dev.off()

## freshwater
temp1 <- subset(md_surv.plot, FW_prev>1)
temp2 <- subset(temp1, !lwr.fw=="NA")
temp3 <- md_surv.plot[md_surv.plot$agent %in% c('rib', 'richness'), ]
md_surv.plot.fw <- rbind(temp2, temp3)
min(md_surv.plot.fw$mlci.fw,na.rm=T)
max(md_surv.plot.fw$muci.fw,na.rm=T)
min(md_surv.plot.fw$midlwr.fw,na.rm=T)
max(md_surv.plot.fw$midupp.fw,na.rm=T)

fw_md_surv.plot2 <-
  md_surv.plot.fw %>%
  ggplot() +
  geom_point(aes(x=mid.fw, y=median.fw, col=agent), size=.5, pch=21) +
  geom_errorbarh(aes(y=median.fw, xmin=midlwr.fw, xmax=midupp.fw, col=agent), linewidth=.25, height=0, alpha=0.7) +
  geom_errorbar(aes(x=mid.fw, y=median.fw, ymax = muci.fw, ymin = mlci.fw, col=agent), linewidth=.25, width=0, alpha=0.7) +
  ylim(-0.07,0.07) +
  xlim(-0.5,0.7) +
  annotate("rect", xmin=-Inf, xmax = 0, ymin=-Inf, ymax=Inf, alpha=.2) +
  annotate("rect", xmin=-Inf, xmax = Inf, ymin=-Inf, ymax=0, alpha=.2) +
  xlab("Survival effect size")+
  ylab("Mass Deviation effect size")+ 
  labs(title="Freshwater-sampled fish") +
  geom_text(aes(x=mid.fw, y=median.fw, label=agent, col=agent), size=4, vjust=-.2) +
  theme(legend.position = "none")


jpeg(filename='figs/Fig_MD_Surv_ONNE2021_FW_230124.jpg', width=500, height=1000, quality=300)
fw_md_surv.plot2
dev.off()

tiff(filename='figs/Fig_MD_Surv_ONNE2021_230124_fwsw_CI.jpg', units="in",width=12, height=12, res=300)
ggarrange(fw_md_surv.plot2, sw_md_surv.plot2, nrow=1)
dev.off()




############################## Plot MD PP by Surv PP

## Seawater
#Remove Circo.virus and Qin - effect size distributions are too wide to have confidence
md_surv.plot<-md_surv#[!(md_surv$agent=="Circo.virus" | md_surv$agent=="Qin" | md_surv$agent=="nu_sal"),]
unique(md_surv.plot$agent)
#myDataFrame["rowName", "columnName"] <- value

# SW plot with agents >1% prevalence
temp1 <- subset(md_surv.plot, SW_prev>1)
temp2 <- md_surv.plot[md_surv.plot$agent %in% c('rib', 'richness'), ]
md_surv.plot.sw <- rbind(temp1, temp2)
md_surv.plot.sw$path_pp.sw <- as.numeric(md_surv.plot.sw$path_pp.sw)
md_surv.plot.sw$prop0.sw <- as.numeric(md_surv.plot.sw$prop0.sw)

sw_md_surv.plot2 <- 
  md_surv.plot.sw %>%
  ggplot() +
  geom_point(aes(x=-prop0.sw, y=-path_pp.sw), size=.8, pch=21) +
  #ylim(-0.07,0.07) +
  #xlim(-0.5,0.7) +
  annotate("rect", xmin=-.5, xmax = -1, ymin=0, ymax=-1, alpha=.2) +
  annotate("rect", xmin=0, xmax = -1, ymin=-.5, ymax=-1, alpha=.2) +
  geom_text(aes(x=-prop0.sw, y=-path_pp.sw, label=plot.agent), 
            col="black", vjust=-.4, size=3, show.legend = F)+
  xlab("Survival PP")+
  ylab("Mass Deviation PP") + 
  labs(title="Marine-sampled fish") +
  scale_y_continuous(breaks = c(-1, -.75,-.5,-.25,0),
                     labels = c("100","75","50","25","0")) +
  scale_x_continuous(breaks = c(-1, -.75,-.5,-.25,0),
                     labels = c("100","75","50","25","0")) +
  theme(legend.position = "none")

jpeg(filename='figs/Fig_MD_Surv_ONNE2021_230124_SW_CI_PP.jpg', width=500, height=500, quality=300)
sw_md_surv.plot2
dev.off()

## freshwater
temp1 <- subset(md_surv.plot, FW_prev>1)
temp2 <- subset(temp1, !lwr.fw=="NA")
temp3 <- md_surv.plot[md_surv.plot$agent %in% c('rib', 'richness'), ]
md_surv.plot.fw <- rbind(temp2, temp3)
min(md_surv.plot.fw$mlci.fw,na.rm=T)
max(md_surv.plot.fw$muci.fw,na.rm=T)
min(md_surv.plot.fw$midlwr.fw,na.rm=T)
max(md_surv.plot.fw$midupp.fw,na.rm=T)
md_surv.plot.fw$path_pp.fw <- as.numeric(md_surv.plot.fw$path_pp.fw)
md_surv.plot.fw$prop0.fw <- as.numeric(md_surv.plot.fw$prop0.fw)

fw_md_surv.plot2 <- 
  md_surv.plot.fw %>%
  ggplot() +
  geom_point(aes(x=-prop0.fw, y=-path_pp.fw), col="black", size=.8, pch=21) +
  annotate("rect", xmin=-.5, xmax = -1, ymin=0, ymax=-1, alpha=.2) +
  annotate("rect", xmin=0, xmax = -1, ymin=-.5, ymax=-1, alpha=.2) +
  geom_text(aes(x=-prop0.fw, y=-path_pp.fw, label=plot.agent), col="black", 
            vjust=-.5, size=3, show.legend = F)+ 
  xlab("Survival PP")+
  ylab("Mass Deviation PP") + 
  labs(title="Freshwater-sampled fish") +
  scale_y_continuous(breaks = c(-1, -.75,-.5,-.25,0),
                     labels = c("100","75","50","25","0")) +
  scale_x_continuous(breaks = c(-1, -.75,-.5,-.25,0),
                     labels = c("100","75","50","25","0")) +
  theme(legend.position = "none")


jpeg(filename='figs/Fig_MD_Surv_ONNE2021_FW_230124_PP.jpg', width=500, height=500, quality=300)
fw_md_surv.plot2
dev.off()

tiff(filename='figs/Fig_MD_Surv_ONNE2021_230124_fwsw_CI_PP.jpg', units="in",width=14, height=7, res=300)
ggarrange(fw_md_surv.plot2, sw_md_surv.plot2, nrow=1)
dev.off()







## Plot total SW detections of each agent by year - note variable prevalence across agents and years
samplesperagent.sw<-inf_agt_resid_data %>% 
  group_by(agent, Year) %>%
  summarise(posdet = sum(posdet)) 
#plot
jpeg(filename='figs/Fig_Agent detections per year_SW_relative_220126.jpg', 
     width=500, height=300, quality=300)
ggplot(data=samplesperagent.sw, aes(x=factor(Year), y=posdet, fill=agent))+
  geom_bar(stat="identity", position="fill")+
  labs(fill="Infectious agents") +
  coord_flip()+
  xlab("Sampling year")+
  ylab("Relative prevalence")+
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(limits = rev) +
  theme(legend.text = element_text(size=11))
dev.off()


## Plot prevalence of agents across years
prevyr <- data.frame(inf_agt_resid_data %>% 
  group_by(agent, Year) %>%
  summarise(
      posdet = sum(posdet), 
      N = sum(N), 
      prev=sum(posdet)/sum(N)))
data_ends.sw <- prevyr %>% filter(Year == 2018)
#plot
jpeg(filename='figs/Fig_Prevalence per year_ONNE2021_SW_220126.jpg', 
     width=350, height=300, quality=300)
ggplot(data=prevyr, aes(x=factor(Year), y=prev, col=agent, group=agent))+
  #geom_point() +
  geom_line() +
  theme(legend.position = "none") +
  geom_text_repel(aes(label=agent), data=data_ends.sw, size=2.5, nudge_x = 1, hjust="right", max.overlaps = 4) +
  xlab("Sampling year") +
  ylab("SW Prevalence") +
  scale_y_continuous(labels = scales::percent) 
dev.off()



## Plot prevalence of agents across years
prevyr_fw <- data.frame(inf_agt_resid_data_fw %>% 
                       group_by(agent, Year) %>%
                       summarise(
                         posdet = sum(posdet), 
                         N = sum(N), 
                         prev=sum(posdet)/sum(N)))
data_ends <- prevyr_fw %>% filter(Year == 2015)
#plot
jpeg(filename='figs/Fig_Prevalence per year_ONNE2021_fw_220126.jpg', 
     width=250, height=300, quality=300)
ggplot(data=prevyr_fw, aes(x=factor(Year), y=prev, col=agent, group=agent))+
  #geom_point() +
  geom_line() +
  theme(legend.position = "none") +
  xlab("Sampling year") +
  ylab("FW Prevalence") +
  geom_text_repel(aes(label=agent), data=data_ends, size=3, nudge_x = 2, hjust="right", max.overlaps = 2) +
  scale_y_continuous(labels = scales::percent) 
dev.off()


## plot full survival index time series
library(scales)
jpeg(filename='figs/Fig_SR 1950 to 2016_220105.jpg', 
     width=450, height=250, quality=300)
ggplot(survial_indices, aes(brood_year, SR_cov_resid, colour = stock_name)) +
  geom_line(alpha=.5)+
  annotate("rect", xmin=2006, xmax=2016, ymin=-5, ymax=Inf, alpha=.3) +
  xlab("Brood year")+
  ylab("SR residuals w covariates")+
  labs(col="Stock") +
  theme(legend.text = element_text(size=10)) +
  guides(color=guide_legend(ncol = 2)) +
  theme_bw()
dev.off()


## plot truncated survival index time series
jpeg(filename='figs/Fig_SR 2006 to 2016_220105.jpg', 
     width=450, height=450, quality=300)
ggplot(survial_indices[survial_indicesL$brood_year>2006,], aes(brood_year, SR_cov_resid, colour = stock_name)) +
  geom_line(alpha=.5)+
  xlab("Brood year")+
  ylab("SR residuals w covariates") +
  labs(col="Stock") +
  theme(legend.text = element_text(size=10), legend.position = c(.2,.25)) +
  guides(color=guide_legend(ncol = 2)) 
dev.off()


### Simple linear models of median - not actually legit statistically speaking
mod.sw <- lm(mid.sw ~ median.sw, data=md_surv)
summary(mod.sw)
mod.fw <- lm(mid.fw ~ median.fw, data=md_surv_fw)
summary(mod.fw)




head(brood_table)
## Plot prevalence of agents across years
tot.sp <- data.frame(brood_table %>% 
                       group_by(BY) %>%
                       summarise(sum(total_recruits)/sum(total_adult_spawner)))
#Plot spawners
jpeg(filename='figs/Fig_recperspawners_ONNE2021.jpg', 
     width=430, height=200, quality=300)
ggplot(tot.sp, aes(BY, sum.total_recruits..sum.total_adult_spawner.)) +
  geom_line() +
  xlab("Brood year") +
  ylab("Recruits/Spawners") +
  theme(legend.text = element_text(size=10), legend.position = c(.2,.25)) +
  guides(color=guide_legend(ncol = 2)) 
#theme_bw()
dev.off()


## Heat Map of Death
#create long version of md_surv
library(pheatmap)
prop0.swfw2 <- matrix(round(as.numeric(unlist(md_surv[,c(43,42,27,13)])),3),ncol=4)
row.names(prop0.swfw2) <- md_surv$plot.agent.y
colnames(prop0.swfw2) <- c("FW survival","SW survival","FW mass deviation","SW mass deviation")
prop0.swfw2 <- prop0.swfw2[order(rowSums(prop0.swfw2),decreasing=T),]
write.csv(prop0.swfw2, file="data/ONNE_SurvMDpostprobs_220126.csv")

jpeg(filename='figs/Fig_heatmap_all_ONNE2021_220126.jpg', 
     width=300, height=600, quality=300)
pheatmap(prop0.swfw2, display_numbers = T, colorRampPalette(c("#0072B2","white","#E69F00"))(100),
         cluster_rows = F, cluster_cols = F, angle_col=90, fontsize_number = 12,
         number_color=1, labels_row = as.expression(rownames(prop0.swfw2)))
dev.off()

# Map study area
# Survival
library(maps)
library(mapdata)
library(maptools)  #for shapefiles
library(scales)  #for transparency
samps_surv <- all.data[,c("Latitude", "Longitude")]   #my data for sampling sites, contains a column of "lat" and a column of "lon" with GPS points in decimal degrees
write.csv(samps_surv, file="data/LatLongONNE_Survival.csv")
samps_md <- d_sw[,c("Latitude", "Longitude")]   #my data for sampling sites, contains a column of "lat" and a column of "lon" with GPS points in decimal degrees
write.csv(samps_md, file="data/LatLongONNE_MD.csv")

jpeg(filename='figs/Fig_Survival_map_all_ONNE2021.jpg', 
     width=600, height=500, quality=300)
map("worldHires","Canada", xlim=c(-140,-115),ylim=c(48,52), col="gray90", fill=TRUE)  #plot the region of Canada I want
map("worldHires","usa", xlim=c(-140,-115),ylim=c(48,52), col="gray95", fill=TRUE, add=TRUE)  #add the adjacent parts of the US; can't forget my homeland
points(samps$Longitude, samps$Latitude, pch=21, col="black", bg="red", lwd=0.2, cex=0.2)  #plot my sample sites
dev.off()

