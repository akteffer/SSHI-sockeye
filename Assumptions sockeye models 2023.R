### Checking model assumptions and cutoffs
## Sockeye survival and mass deviation relative to pathogens
## A. Teffer

## Bring in data:
#Rhat and Neff

rhatneff_surv_mod <- read.csv("data/Posterior distributions_Rhat and Neff_stspec_prev_220126_REDUCED.csv")
colnames(rhatneff_surv_mod) <- c("param","value","agent","metric")
head(rhatneff_surv_mod)
rhat.surv <- rhatneff_surv_mod[rhatneff_surv_mod$metric=="Rhat",]
neff.surv <- rhatneff_surv_mod[rhatneff_surv_mod$metric=="Neff",]
rhat.surv <- rhat.surv[order(rhat.surv$value),]
neff.surv <- neff.surv[order(neff.surv$value),]

#range in rhat
min(rhat.surv$value)
max(rhat.surv$value)
plot(rhat.surv$value)

#range in neff
min(neff.surv$value)
max(neff.surv$value)
plot(neff.surv$value)

# Mass deviation models
#SW
MD.sw.full <- readRDS("sockeye_SW_lw_results_jan24_2023.rds")
MD.sw.convergence_data <- MD.sw.full$convergence_data
MD.sw.convergence_summary <- as.data.frame(MD.sw.full$convergence_summary)
MD.sw.convergence_summary$agent <- rownames(MD.sw.convergence_summary)
colnames(MD.sw.convergence_summary) <- c("Rhatb0","essb0","Rhatb1","essb1","Rhatb2","essb2","Rhatb3","essb3","agent")
#FW random slope
MD.fw.full <- readRDS("sockeye_FW_lw_results_jan24_2023.rds")
MD.fw.convergence_data <- MD.fw.full$convergence_data
MD.fw.convergence_summary <- as.data.frame(MD.fw.full$convergence_summary)
MD.fw.convergence_summary$agent <- rownames(MD.fw.convergence_summary)
colnames(MD.fw.convergence_summary) <- c("Rhatb0","essb0","Rhatb1","essb1","Rhatb2","essb2","agent")
#FW no random slope
MD.fw.full.nrs <- readRDS("sockeye_FW_norandslope_lw_results_jan26_2023.rds")
MD.fw.convergence_data <- MD.fw.full.nrs$convergence_data
MD.fw.convergence_summary <- as.data.frame(MD.fw.full.nrs$convergence_summary)
MD.fw.convergence_summary$agent <- rownames(MD.fw.convergence_summary)
colnames(MD.fw.convergence_summary) <- c("Rhatb0","essb0","Rhatb1","essb1","Rhatb2","essb2","agent")

#Art: There is a SST term in the marine model but not freshwater. b0 is intercept, b1 is length, B2 is 
#hyper parameter for JAZ random effect on pathogen slope (the estimate you really care about), b3 is SST 
#pathogen is power transformed in MD models



# Plot rhat summaries for FW random slope
ggplot(data= MD.fw.convergence_summary, aes(x=reorder(agent,Rhatb0), y=Rhatb0))+
  geom_bar(stat="identity") +
  ylim(0,1.5)+
  geom_hline(yintercept=1.1)+
  coord_flip()
ggplot(data= MD.fw.convergence_summary, aes(x=reorder(agent,Rhatb1), y=Rhatb1))+
  geom_bar(stat="identity") +
  ylim(0,1.5)+
  geom_hline(yintercept=1.1)+
  coord_flip()
ggplot(data= MD.fw.convergence_summary, aes(x=reorder(agent,Rhatb2), y=Rhatb2))+
  geom_bar(stat="identity") +
  ylim(0,1.5)+
  geom_hline(yintercept=1.1)+
  coord_flip()
##Agents over the line for B0 - all others look good (1.0-1.1)
# Circo.virus, nu_sal, Rhabdo3_virus, richness, arena2
## I think I will keep richness and arena2 with disclaimer

# Plot rhat summaries for FW NO random slope
ggplot(data= MD.fw.convergence_summary, aes(x=reorder(agent,Rhatb0), y=Rhatb0))+
  geom_bar(stat="identity") +
  ylim(0,1.5)+
  geom_hline(yintercept=1.1)+
  coord_flip()
ggplot(data= MD.fw.convergence_summary, aes(x=reorder(agent,Rhatb1), y=Rhatb1))+
  geom_bar(stat="identity") +
  ylim(0,1.5)+
  geom_hline(yintercept=1.1)+
  coord_flip()
ggplot(data= MD.fw.convergence_summary, aes(x=reorder(agent,Rhatb2), y=Rhatb2))+
  geom_bar(stat="identity") +
  ylim(0,1.5)+
  geom_hline(yintercept=1.1)+
  coord_flip()
##Agents over the line for B0 - all others look good (1.0-1.1)
# Circo.virus, nu_sal, Rhabdo3_virus, richness, arena2
## I think I will keep richness and arena2 with disclaimer

# Plot rhat summaries for SW
ggplot(data= MD.sw.convergence_summary, aes(x=reorder(agent,Rhatb0), y=Rhatb0))+
  geom_bar(stat="identity") +
  ylim(0,1.5)+
  geom_hline(yintercept=1.1)+
  coord_flip()
ggplot(data= MD.sw.convergence_summary, aes(x=reorder(agent,Rhatb1), y=Rhatb1))+
  geom_bar(stat="identity") +
  ylim(0,1.5)+
  geom_hline(yintercept=1.1)+
  coord_flip()
ggplot(data= MD.sw.convergence_summary, aes(x=reorder(agent,Rhatb2), y=Rhatb2))+
  geom_bar(stat="identity") +
  ylim(0,1.5)+
  geom_hline(yintercept=1.1)+
  coord_flip()
ggplot(data= MD.sw.convergence_summary, aes(x=reorder(agent,Rhatb3), y=Rhatb3))+
  geom_bar(stat="identity") +
  ylim(0,1.5)+
  geom_hline(yintercept=1.1)+
  coord_flip()
##Agents over the line for B0 - all others look good (1.0-1.1)
## Rhabdo3_virus, fa_mar


######################chilko###########################
rhatneff_surv_mod_ch <- read.csv("data/Chilko_MarineSurvival_rhatneff_230222.csv")
head(rhatneff_surv_mod_ch)
rhatneff_surv_mod_ch$prev_std <- as.numeric(rhatneff_surv_mod_ch$prev_std)
rhatneff_surv_mod_ch$early_sst_stnd <- as.numeric(rhatneff_surv_mod_ch$early_sst_stnd)
rhat.surv.ch <- rhatneff_surv_mod_ch[rhatneff_surv_mod_ch$metric=="Rhat",]
neff.surv.ch <- rhatneff_surv_mod_ch[rhatneff_surv_mod_ch$metric=="Neff",]
rhat.surv.ch <- rhat.surv.ch[order(rhat.surv.ch$prev_std),]
neff.surv.ch <- neff.surv.ch[order(neff.surv.ch$prev_std),]

#range in rhat
min(rhat.surv.ch$prev_std)
max(rhat.surv.ch$prev_std)
plot(rhat.surv.ch$prev_std)

#range in neff
min(neff.surv.ch$prev_std)
max(neff.surv.ch$prev_std)
plot(neff.surv.ch$prev_std)

# Rhat prev_std
ggplot(data= rhat.surv.ch, aes(x=reorder(X,prev_std), y=prev_std))+
  geom_bar(stat="identity") +
  ylim(0,1.5)+
  geom_hline(yintercept=1.1)+
  coord_flip()

# Rhat SST
ggplot(data= rhat.surv.ch, aes(x=reorder(X,early_sst_stnd), y=early_sst_stnd))+
  geom_bar(stat="identity") +
  ylim(0,1.5)+
  geom_hline(yintercept=1.1)+
  coord_flip()

# Neff
ggplot(data= neff.surv.ch, aes(x=reorder(X,prev_std), y=prev_std))+
  geom_bar(stat="identity") +
 # ylim(0,1.5)+
  #geom_hline(yintercept=1.1)+
  coord_flip()


