# SAFO MHC 2019
# Fin clips from brook trout in Western MA stream systems
library(ggplot2)
library(data.table)
library(reshape2)
library(ggpubr)

install.packages("devtools") 
devtools::install_github("BlakeRMills/MetBrewer")

setwd("~/Smith Conservation Fellow/Fasta Files")
mhc <- read.csv("MHC_AKT_211123.csv")
head(mhc)

# Create long version
names(mhc)
mhcl = melt(mhc, id.vars = c("Sample", "SampLoc", "brook"),
             measure.vars = c("V1", "V2", "V3","V4", "V5", "V6","V7", "V8"))
# Create df per brook
rb <- mhcl[mhcl$brook=="Roaring",]
fb <- mhcl[mhcl$brook=="Fourmile",]
wb <- mhcl[mhcl$brook=="West",]
ab <- mhcl[mhcl$brook=="Avery",]
ww <- mhcl[mhcl$brook=="West Whately",]

# Plots
ggplot(mhcl, aes(fill=variable, y=value, x=SampLoc)) + 
  geom_bar(position="fill", stat="identity")

# Pie charts by brook
pie.temp <-ggplot(rb, aes(x="", y=value, fill=variable)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  labs(title="Roaring Brook") +
  theme_void() +
  guides(fill=guide_legend(title="MHCII variant"))
pie.rb <- ggplot(rb, aes(x="", y=value, fill=variable)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  labs(title="Roaring Brook") +
  theme_void() +
  theme(legend.position = "none")
pie.fb <- ggplot(fb, aes(x="", y=value, fill=variable)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  labs(title="Fourmile Brook") +
  theme_void()+
  theme(legend.position = "none")
pie.wb <- ggplot(wb, aes(x="", y=value, fill=variable)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  labs(title="West Brook") +
  theme_void()+
  theme(legend.position = "none")
pie.ab <- ggplot(ab, aes(x="", y=value, fill=variable)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  labs(title="Avery Brook") +
  theme_void()+
  theme(legend.position = "none")
pie.ww <- ggplot(ww, aes(x="", y=value, fill=variable)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  labs(title="West Whately") +
  theme_void()+
  theme(legend.position = "none")

leg <- get_legend(pie.temp)
ggarrange(pie.rb, pie.fb, pie.wb, pie.ab, pie.ww, leg)

#print plot
jpeg("SAFO_MHCII_pie_brooks.jpg", units="in", width=6, height=6, res=300)
ggarrange(pie.rb, pie.fb, pie.wb, pie.ab, pie.ww, leg)
dev.off()

# by brook
jpeg("SAFO_MHCII_pie_rb.jpg", units="in", width=4, height=3, res=300)
ggarrange(pie.rb, leg)
dev.off()
jpeg("SAFO_MHCII_pie_fb.jpg", units="in", width=4, height=3, res=300)
ggarrange(pie.fb, leg)
dev.off()
jpeg("SAFO_MHCII_pie_wb.jpg", units="in", width=4, height=3, res=300)
ggarrange(pie.wb, leg)
dev.off()
jpeg("SAFO_MHCII_pie_ab.jpg", units="in", width=4, height=3, res=300)
ggarrange(pie.ab, leg)
dev.off()
jpeg("SAFO_MHCII_pie_ww.jpg", units="in", width=4, height=3, res=300)
ggarrange(pie.ww, leg)
dev.off()

jpeg("SAFO_MHCII_SampLoc_bar.jpg", units="in", width=5, height=1, res=300)
ggplot(mhcl, aes(fill=variable, y=value, x=SampLoc)) + 
  geom_bar(position="fill", stat="identity")
dev.off()

