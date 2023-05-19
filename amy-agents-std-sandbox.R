library(tidyverse)

#### Read in data, clean, standardize - stock-specific metrics

##### SW

inf_agt_resid_data <- read.csv("data/REDUCED_ONNE productivity infection analysis_stockspecific_210630.csv")

head(inf_agt_resid_data)



##### FW

inf_agt_resid_data_fw <- read.csv("data/REDUCED_ONNE productivity infection analysis_stockspecific_FW_210630.csv")

str(inf_agt_resid_data_fw)



# Data cleaning

# Standardize and incorporate into SW df

inf_std <- plyr::ddply(inf_agt_resid_data, c("agent"),function(x) {
  
  scaled_prev <- scale(x$prev)
  
  scaled_load <- scale(x$mean_load)
  
  xx <- data.frame(scaled_prev, scaled_load)
  
})

inf_agt_resid_data$prev_std <- inf_std[,2]

inf_agt_resid_data$load_std <- inf_std[,3]

# Do not replace RIB and richness with scales values

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
                                        
                                        "de_sal" = "D. salmonis",
                                        
                                        "richness" = "Richness",
                                        
                                        "rib" = "RIB")

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



## Check on relationship between prevalence and standardized prevalence 

temp <- inf_agt_resid_data[inf_agt_resid_data$agent=="c_b_cys",]

ggplot(temp) +
  geom_point(aes(x=prev, y=prev_std), pch=21)+
  labs(x="prev", y="scaled prev")



#This is the code I use to calculate RIB if that helps:
  
ribdata_sw <- as.matrix(sw.data[,28:97]) #pull out raw agent values

rib1_sw <- data.frame(apply(ribdata_sw, 2, function(x) {x/max (x, na.rm=TRUE)})) #calculate the load divided by the max load for each agent

rib2_sw <- data.frame(rowSums(rib1_sw, na.rm=TRUE)) #sum across agents within each fish then bind this to the sw.data df



#This is the code I use to calculate N, total positive values, mean RIB within each stock-year group, median RIB, and prevalence including 0s: 
  
  stock.rib.sw =
  data.frame(
    spsu.stock %>% 
      summarise(
        length(which(rib!="NA")), #samples
        length(which(rib>0)), #positive detections
        mean(rib[rib!=0], na.rm=TRUE), 
        median(rib[rib!=0], na.rm=TRUE),
        (length(which(rib>0)) / length(which(!is.na(rib)))) * mean(rib[rib!=0], na.rm=TRUE),
        "rib"
      )
  )


