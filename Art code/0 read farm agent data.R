library(viridis)

#read the data
#rm(list=ls())

#copy-number data
dat0 = read.csv('Full Export 19Oct18 AMD.csv')
dat0 = dat0[dat0$Prog == 'AMD',]
audit.details = read.csv('Audit sample data.csv')
farms = read.csv('salmon farms.csv')

#metadata - what are the microbes, which microbes were screened in the first and second batch of tests?
microbes = read.csv('microbes.csv')
microbes$FW[is.na(microbes$FW)] = -1
microbes$SW[is.na(microbes$SW)] = -1
microbes$test = as.character(microbes$test)
microbes = microbes[as.character(microbes$test)!='hkg',]

dat0 = dat0[order(as.character(dat0$Unique)),]
dat0 = cbind(1:nrow(dat0), dat0)
names(dat0)[1] = 'entry'

#add in the house-keeping gene data (hkg doesn't have copy number associated)
# datCT = read.csv('SSHI-CT-all.csv')[dat0[,1],]
# dat0$hkgCT = datCT$hkg
# dat0$hkgCT[abs(dat0$hkgCT)>200] = NA

dat = dat0#[,colnames(dat0)!='hkg'] #remove "house-keeping gene"
rm(dat0)
#rm(dat0.revised)

first.tests = c('ae_hyd','gy_sal','ipnv','isav7','isav8','mo_vis','nu_sal','omv','pmcv','sav','sp_sal')
second.tests = c('arena1','arena2','ascv','cov','ctv','mitov','ortho','reov','sgpx','smallUK','toti')

dat$Date = as.Date(as.character(dat$Date), format="%d-%b-%y")
dat[dat==-999] = NA     #remove '-999's

dat$numDate = as.numeric(dat$Date)
dat$DoY = as.numeric(format(dat$Date, '%j'))

dat$farm = unlist(sapply(as.character(dat$Fish), function(x){
    if(x %in% as.character(audit.details$geneticsassignednumber)){
        return(as.character(audit.details$Farm[as.character(audit.details$geneticsassignednumber)==x]))
    } else {
        return(NA)
    }
}))
dat$farm[dat$farm=='Baxter'] = 'Baxter Islets'
dat$farm[dat$farm=='Newcombe'] = 'Newcomb'
dat$farm[dat$farm=='Williamson Pass'] = 'Williamson'
dat$farm[dat$farm=='MacIntyre'] = 'McIntyre Lake'
dat$farm[dat$farm=='McCalls'] = 'McCall Islets'
dat$farm[dat$farm=='Doyle'] = 'Doyle Island'
dat$farm[dat$farm=='Cyrus Rocks'] = 'Cyrus Rock'
dat$farm[dat$farm=='Sir Edmond'] = 'Sir Edmund Bay'
dat$farm[dat$farm=='Dixon'] = 'Dixon Bay'
dat$farm[dat$farm=='Koskimo Bay'] = 'Koskimo'
dat$farm[dat$farm=='Doctor Islet'] = 'Doctor Islets'
dat$farm[dat$farm=='Hardwicke Island'] = 'Hardwicke'
dat$farm[dat$farm=='Sheep Pass'] = 'Sheep Passage'
dat$farm[dat$farm=='Saranac'] = 'Saranac Island'
dat$farm[dat$farm=='Swanson Island'] = 'Swanson'
dat$farm[dat$farm=='Midsummer Island'] = 'Midsummer'
dat$farm[dat$farm=='Chancellor'] = 'Chancellor Channel'
dat$farm[dat$farm=='Thurlow Point'] = 'Thurlow'
dat$farm[dat$farm=='Mt. Simmonds'] = 'Simmonds Point'
dat$farm[dat$farm=='Ahlstrom or Farm 13' & dat$Longitude<(-124)] = 'Ahlstrom'
dat$farm[dat$farm=='Ahlstrom or Farm 13' & dat$Longitude>(-124)] = 'Site 13'
dat$farm[dat$farm=='Lime Bay'] = 'Lime Point'
dat$farm[dat$farm=='Middle Point Bay'] = 'Middle Point'
dat$farm[dat$farm=='Conville Bay'] = 'Conville Point'   #lat/long indicate the Conville Point farm is more likely
dat$farm[dat$farm=='Bennett Point'] = 'Noo-la'  #name changed

#sum(is.na(dat$farm))
###############    ditch samples for which farm & stocking-date data are lacking ##############
dat = dat[!is.na(dat$farm),]

dat$sample.Date = unlist(sapply(as.character(dat$Fish), function(x){
    if(x %in% as.character(audit.details$geneticsassignednumber)){
        return(as.character(audit.details$coll_date[as.character(audit.details$geneticsassignednumber)==x]))
    } else {
        return(NA)
    }
}))
#convert the sampling date to R's date format.  Note that Excel incorrectly interpreted dates entered as d/m/Y as being in m/d/Y format ...leading to some sampling dates preceding saltwater entry dates
dat$sample.Date = sapply(dat$sample.Date, function(x){
        if(!is.na(x)){
            if(substr(x,3,3)=='/'){
                return(as.Date(x,format='%d/%m/%Y'))
            } else return(as.Date(x,format='%Y-%d-%m'))
        } else return(NA)
    })
dat$entry.Date = unlist(sapply(as.character(dat$Fish), function(x){
    if(x %in% as.character(audit.details$geneticsassignednumber)){
        return(as.character(audit.details$sw_entry[as.character(audit.details$geneticsassignednumber)==x]))
    } else {
        return(NA)
    }
}))
dat$entry.Date = sapply(dat$entry.Date, function(x){
    if(substr(x,3,3)=='/'){
        return(as.Date(x,format='%d/%m/%Y'))
    } else return(as.Date(x,format='%Y-%d-%m'))
})
dat$SWdays = as.numeric(dat$sample.Date-dat$entry.Date)


#omit data where days-in-saltwater information is absent or (apparently) in error
dat = dat[!is.na(dat$SWdays) & dat$SWdays>=0,]


#plot PRV vs time and a basic logistic model fit
mod = glm(cbind(prv>0,prv==0)~SWdays, family=binomial, data=dat)
step.size = 25
int.starts = seq(-100,1500,step.size)
props = sapply(int.starts, function(x){
    sum(dat$prv[dat$SWdays>=x & dat$SWdays<(x+step.size)]>0, na.rm=TRUE)/sum(dat$prv[dat$SWdays>=x & dat$SWdays<(x+step.size)]>(-1),na.rm=TRUE)
    })
counts = sapply(int.starts, function(x){sum(dat$prv[dat$SWdays>=x & dat$SWdays<(x+step.size)]>(-1),na.rm=TRUE)})

plot(dat$SWdays,(dat$prv>0)*1.05-0.025, 
     #pch=16 col='#00000033',
     pch='|', col='#000000', cex=0.5, lwd=2,
     xlab='days in salt water', ylab='PRV infection')
points(int.starts+step.size/2, props, pch=16, cex=sqrt(counts)/2, col='#00000055')
lines(0:max(dat$SWdays, na.rm=TRUE), 
      plogis(mod$coefficients[1] + (0:max(dat$SWdays, na.rm=TRUE))*mod$coefficients[2]),
      col='red')

# plot(dat$SWdays[!is.na(dat$SWdays)],log(dat$prv+1)[!is.na(dat$SWdays)], pch=16, col='#00000033',
#      #ylim = range((dat$prv)[!is.na(dat$SWdays)], na.rm=TRUE),
#      xlab='days in salt water', ylab='PRV load')
# 
# 
# 
# 
# 
# dat.farm = tapply(1:nrow(dat), dat$Stock, function(y){
#         x = dat[y,]
#         print(x)
#         return(c(mean(x$prv), sum(x$prv>0), length(x$prv), x$Latitude[1], x$Longitude[1]))
#     })
# 
# 
# 
# jit.sd = 0#0.05
# plot(dat$Longitude+rnorm(nrow(dat),0,jit.sd), dat$Latitude+rnorm(nrow(dat),0,jit.sd), pch=16,
#     col = magma(100, alpha=0.5)[trunc(100*sqrt(dat$prv)/max(sqrt(dat$prv), na.rm=TRUE))]
#     ,xlim=c(-127,-126)
#     ,ylim=c(50.25,51.25)
#     )
# points(farms$Longitude,farms$Latitude, pch=16, col='red')
# 
# plot(dat$Longitude+rnorm(nrow(dat),0,0.05), dat$Latitude+rnorm(nrow(dat),0,0.05), pch=16,
#      col = magma(100, alpha=0.5)[trunc(100*sqrt(dat$prv)/max(sqrt(dat$prv), na.rm=TRUE))]
# )
# 
# 
# 
# loc = tapply(1:nrow(dat), dat$STATION_NAME, function(x){dat[x,c('Longitude','Latitude')][1,]})
# locs = data.frame(location=c(0),long=c(0),lat=c(0))
# for(i in 1:length(loc)){
#     locs[i,] = c(names(loc)[i], loc[[i]][1], loc[[i]][2])  
# }
# 
# 
# farm.prv = tapply(dat$prv,dat$farm, function(x){mean(x, na.rm=TRUE)})
# hist((farm.prv)^(1/10), breaks=30)
