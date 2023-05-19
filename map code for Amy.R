#load libraries
{
    library(ReacTran)
    library(magick)
    library(GA)
    library(mvtnorm)
    library(parallel)
    library(rgeos)
    library(geosphere)
    library(PBSmapping)
    library(sf)
    library(maptools)
    library(raster)
}
data("nepacLL")

agent = 'te_mar'
dat = read.csv('sockeye_agent_data_processed.csv')    #processed data from Amy

#  COMPARE DATA USED HERE WITH THOSE USED BY AMY TEFFER in SR ANALYSIS
#Amy.dat = read.csv('ONNE metadata no LOD_6.17.2021.csv') #data for cross referencing
#Amy.positives = Amy.dat$Unique[!is.na(Amy.dat$te_mar) & Amy.dat$te_mar>0]
#positives = dat$Unique[!is.na(dat$te_mar) & dat$te_mar>0]
#c(positives[!(positives%in%Amy.positives)], Amy.positives[!(Amy.positives%in%positives)])

#fish data processing
{
    dat = dat[dat$Year %in% years,]
    dat$DoY = as.numeric(format(as.Date(dat$Date, format='%d-%b-%y'), '%j'))
    dat = dat[dat$DoY>min.day & dat$DoY<max.day,] #use points up to the end of September (to reduce calculation of seaway distance for unnecessary points)
    dist.dat = read.csv('seaway distances - sockeye.csv')
    origin.points = data.frame(Longitude=-123.2060,Latitude=49.1491) #saltwater entry at the mouth of the Fraser River
    #subset points to omit those off the west coast of Vancouver Island 
    #(outside the scope of the linearised model)
    dat = dat[(dat$Latitude-7.3+1/3*dat$Longitude)>0,]# omit ECVI fish           & (dat$Latitude-90.2-1/3*dat$Longitude)>0,] #retain those southeast of the Fraser
    #...and omit stocks that don't directly migrate north
    dat = dat[!is.na(dat$Run_Timing),]
    slope=0.5; 
    dat$migration = trunc(sapply(1:nrow(dat),function(i){
        dist.dat$start_loc1[abs(dist.dat$Latitude-dat$Latitude[i])<0.0001 & abs(dist.dat$Longitude-dat$Longitude[i])<0.0001][1]/1000 #+ x.SW  = adjustment when origin was Mission
    }))
    #assign negative distance from origin for those fish that have moved southeast of the Fraser mouth
    dat$migration[dat$Latitude<((origin.points$Latitude-slope*origin.points$Longitude)+slope*dat$Longitude)] = -dat$migration[dat$Latitude<((origin.points$Latitude-slope*origin.points$Longitude)+slope*dat$Longitude)]
    #dat$set = as.factor(dat$set)
}
observations = as.numeric((dat[,agent]>0)[!is.na(dat[,agent])])
write.csv(dat, "cleaned_data.csv")


UTM.CRS = CRS('+proj=utm +zone=19 +datum=NAD83 +units=m +no_defs +ellps=WGS84')
LL.CRS = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')    #the '+ellps' argument doesn't "stick"
AM.zones = st_read(dsn='./AMD_zones', layer='DFO_FHZONES_POLY')
AM.zones.sp = spTransform(as(AM.zones[2:10,c('ZONE','geometry')],'Spatial'), LL.CRS)#proj.to.plot)

library(latex2exp)
library(sp)
library(sf)
library(bcmapsdata)
library(bcmaps)

use.cex = 0.8


######################################################################################
#plotting stuff and things
LL.CRS = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')    #the '+ellps' argument doesn't "stick"
AM.zones = st_read(dsn='./AMD_zones', layer='DFO_FHZONES_POLY')
AM.zones.sp = spTransform(as(AM.zones[2:10,],'Spatial'), LL.CRS)#proj.to.plot)
BC.Albers = CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 
                +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
BC = st_transform(get_layer("bc_bound_hres", class = "sf"), LL.CRS)#proj.to.plot)
BC = as(BC, 'Spatial')
farms.to.plot = farm.stocking$Facility.Name[farm.stocking$total.months>0]
DI.farms = farms[farms$Longitude>(-126) & farms$Longitude<(-124.4) & farms$Latitude>49.8,]


scaling = 0.75
#X11(display = "", width = 4.75, height=4.75)

# ZOOMED IN ##################################################################################################
#png("Figure 1 NEW.png", width = 6.75, height=6.25, units='in', res=400)
#pdf("Figure A1.pdf", width = 6.75, height=6.25)

samp.cex = 0.75
farm.cex = 1
par(mfrow=c(1,1), mar=c(4,4,0.5,0.5), cex=use.cex)
plot(BC, xlim=c(-127.0,-123.25), ylim=c(49,51), col='grey70', border=NA,
     xlab='longitude (°E)', ylab='latitude (°N)')
polygon(c(-127,-127,-126.1,-126.1),c(50.5,50.98,50.98,50.5), bg=NA, lty=2) #Broughton
text(-127-0.06,50.98+0.04, 'Broughton Archipelago', pos=4, font=1)
polygon(c(-126,-126,-124.6,-124.6),c(49.95,50.55,50.55,49.95), bg=NA, lty=2) #Discovery Islands
text(-126-0.06,50.55+0.04, 'Discovery Islands', pos=4, font=1)
axis(1); axis(2); box() #xlim=c(-133,-122), ylim=c(48,54)
#don't plot points that appear to be on land (dealt with in migration-distance code, but still present in raw data)
show.dat = dat[!(dat$Longitude<(-125.45) & dat$Latitude<50.2) & !(dat$Longitude<(-126.05) & dat$Latitude<50.4),]
points(show.dat$Longitude,show.dat$Latitude,pch=16, col='#1100aa', cex=samp.cex)#, col='#1100aa55', cex=0.75)
points(origin.points, pch=24, cex=1.5, col='white', bg='black', lwd=1.5)     
#abline(origin.points$Latitude-slope*origin.points$Longitude,slope) #cutoff line for neg/pos [southeast/northwest] direction
#plot(AM.zones.sp, add=TRUE, border='grey20', lty=1, lwd=1.5)
coordinates = t(sapply(1:length(AM.zones.sp@polygons), function(i){AM.zones.sp@polygons[[i]]@labpt}))
#INFECTED
points(dat[dat$te_mar>0,]$Longitude,dat[dat$te_mar>0,]$Latitude,
       pch=21, col='#FF0000', bg='white',cex=samp.cex, lwd=1)#, col='#1100aa55', cex=0.75)
#FARMS
points(farms$Longitude[farms$name%in%farms.to.plot],farms$Latitude[farms$name%in%farms.to.plot], pch=3, col='black', #'#ff8800FF', #
       cex=farm.cex, lwd=1)
#points(DI.farms$Longitude,DI.farms$Latitude, pch=16, col='red', cex=0.25)
# points(mean(hakDI$long),mean(hakDI$lat),pch=21, cex=1.5, col='white', bg='black', lwd=1.5)
# points(mean(hakJS$long),mean(hakJS$lat),pch=22, cex=1.5, col='white', bg='black', lwd=1.5)
#plot(BC, add=TRUE, col='grey50', border=NA)
text(-130,49.5,'Pacific Ocean', srt=-39, cex=2, col='grey80')
text(-125.5,49.55,'Vancouver Island', srt=-39, cex=1.5, col='white')
text(-123.75,47.9,'USA', cex=1, col='white')
text(-124,50.9,'British', col='white', cex=1.5)
text(-124,50.8,'Columbia', col='white', cex=1.5)

#dev.off()
############################################################################################################



# PROVINCIAL ##################################################################################################
#png("Figure 1.png", width = 4.75, height=4.75, units='in', res=400)
#pdf("Figure 1.pdf", width = 4.75, height=4.75)

par(mfrow=c(1,1), mar=c(4,4,0.5,0.5), cex=use.cex)
plot(sp.shore, xlim=c(-133.5,-122.5), ylim=c(48,54.75), col='grey80', border=NA,
     xlab='longitude (°E)', ylab='latitude (°N)')
plot(BC, add=TRUE, col='grey50', border=NA)
axis(1); axis(2); box() #xlim=c(-133,-122), ylim=c(48,54)
#don't plot points that appear to be on land (dealt with in migration-distance code, but still present in raw data)
show.dat = dat[!(dat$Longitude<(-125.45) & dat$Latitude<50.2) & !(dat$Longitude<(-126.05) & dat$Latitude<50.4),]
points(show.dat$Longitude,show.dat$Latitude,pch=16, col='#1100aa', cex=0.75)#, col='#1100aa55', cex=0.75)
points(origin.points, pch=24, cex=1.5, col='white', bg='black', lwd=1.5)     
#abline(origin.points$Latitude-slope*origin.points$Longitude,slope) #cutoff line for neg/pos [southeast/northwest] direction
plot(AM.zones.sp, add=TRUE, border='grey20', lty=1, lwd=1.5)
coordinates = t(sapply(1:length(AM.zones.sp@polygons), function(i){AM.zones.sp@polygons[[i]]@labpt}))
text(coordinates[AM.zones.sp$ZONE%in%paste('3-',1:5,sep=''),1]+0.15,
     coordinates[AM.zones.sp$ZONE%in%paste('3-',1:5,sep=''),2]+0.25,
     AM.zones.sp[AM.zones.sp$ZONE%in%paste('3-',1:5,sep=''),]$ZONE, col='grey20')
text(coordinates[!AM.zones.sp$ZONE%in%paste('3-',1:5,sep=''),1]-c(0,0.15,0.15,0.15),
     coordinates[!AM.zones.sp$ZONE%in%paste('3-',1:5,sep=''),2]+0.15,
     AM.zones.sp[!AM.zones.sp$ZONE%in%paste('3-',1:5,sep=''),]$ZONE, col='grey20')
points(farms$Longitude[farms$name%in%farms.to.plot],farms$Latitude[farms$name%in%farms.to.plot], pch=3, col='#ff8800FF', cex=1, lwd=1)
#points(DI.farms$Longitude,DI.farms$Latitude, pch=16, col='red', cex=0.25)
points(mean(hakDI$long),mean(hakDI$lat),pch=21, cex=1.5, col='white', bg='black', lwd=1.5)
points(mean(hakJS$long),mean(hakJS$lat),pch=22, cex=1.5, col='white', bg='black', lwd=1.5)
text(-130,49.5,'Pacific Ocean', srt=-39, cex=2, col='grey80')
text(-125.85,49.75,'Vancouver Island', srt=-39, cex=1, col='white')
text(-123.75,47.9,'USA', cex=1, col='white')
text(-125.5,54.0,'British', col='white', cex=1)
text(-125.5,53.7,'Columbia', col='white', cex=1)


#dev.off()
############################################################################################################
