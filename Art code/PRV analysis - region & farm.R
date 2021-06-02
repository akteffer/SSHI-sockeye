##################################################################################################################
#              SETUP
##################################################################################################################

source('0 read farm agent data.R')

# dist.all = as.matrix(read.csv('distances in zone all.csv', row.names=1))
# for(zone in c('3.5', '3.4', '3.3', '3.1', '3.2', '2.4', '2.3')){
#     assign(paste('dist',zone,sep='.'), 
#            as.matrix(read.csv(paste('distances in zone ',zone,'.csv',sep=''), row.names=1)))
# }
# #distance correlations between simple and complex shoreline maps are >0.9999 
# #except for Aquaculture management zone 3.3 (0.9988)
# # dist.comp=dist.3.3; plot(dist.comp,dist.all[rownames(dist.comp),rownames(dist.comp)]); abline(0,1)
# # cor(c(dist.comp),c(dist.all[rownames(dist.comp),rownames(dist.comp)]))
# #...replace distances from coarse map for zone 3.3 from distances derived from detailed map
# dist.all[rownames(dist.3.3),rownames(dist.3.3)] = dist.3.3

#farms = read.csv('salmon farms.csv')
farms$ID = paste('F',farms$ref,sep='')
farms$X = farms$Longitude
farms$Y = farms$Latitude

dat$ID = sapply(dat$farm, function(x){farms$ID[as.character(farms$name)==x]})
                

dat$farmXdate = paste(dat$ID,dat$numDate,sep='.')
dat$zone = sapply(dat$ID, function(x){farms$zone[farms$ID==x]})
dat = dat[!is.na(dat$prv),]

#inter-sampling intervals
# hist(unique(dat$Date[dat$farm==unique(dat$farm)[7]]), breaks=100)
intervals = lapply(unique(dat$farm), function(x){diff(sort(unique(dat$numDate[dat$farm==x])))
    })
hist(unlist(intervals), breaks=50, xlim=c(0,900))


##################################################################################################################
#                    STAN
##################################################################################################################
#library(R2jags)
library(rstan)
#library(mgcv)
library(coda)
library(extraDistr)
library(rlist)
library(rstanarm)

Stan.data = list(
    n = as.integer(nrow(dat)),    #observations 
    PRV = as.integer(dat$prv>0),  #response (0: no PRV, 1: PRV)
    
    n_covariates = as.integer(1), #intercept and time-in-seawater effect
    X = as.matrix(data.frame(SWdays=dat$SWdays)),
    X.names = c('SWdays'),  
    #x.means = predictor.means,
    
    n_farms = as.integer(length(unique(dat$farm))),    
    farms = matrix(as.numeric(t(sapply(dat$farm, function(x){x==unique(dat$farm)}))), byrow=FALSE, nrow=nrow(dat)), #farm associated with each observation (design matrix)
    farm.names = unique(dat$farm),
    
    n_zones = as.integer(length(unique(dat$zone))),    
    zones = matrix(as.numeric(t(sapply(dat$zone, function(x){x==unique(dat$zone)}))), byrow=FALSE, nrow=nrow(dat)), #farm associated with each observation (design matrix)
    zone.names = unique(dat$zone),
    
    n_samples = as.integer(length(unique(dat$farmXdate))), #number of fish-health-audit sampling events
    samples = matrix(as.numeric(t(sapply(dat$farmXdate, function(x){x==unique(dat$farmXdate)}))), byrow=FALSE, nrow=nrow(dat)),  #visit ID for each sample (design matrix)
    sample_dates = unique(dat$farmXdate),
    sample_SWdays = sapply(unique(dat$farmXdate), function(x){dat$SWdays[x==dat$farmXdate][1]})
    # n_samples = as.integer(length(unique(dat$farm))), #number of fish-health-audit sampling events
    # samples = matrix(as.numeric(t(sapply(dat$farm, function(x){x==unique(dat$farm)}))), byrow=FALSE, nrow=nrow(dat)),  #visit ID for each sample (design matrix)
    # sample_dates = unique(dat$farm)
)

saveRDS(Stan.data,'Stan_data.rds')

Stan.inits = function(i=1){
    return(list(
        beta_0 = c(mod$coefficients)[1],
        beta = as.array(c(mod$coefficients)[2]),
        d_decompose = rep(0, length(unique(dat$ID))),
        
        s_prelim = rep(0, length(unique(dat$farmXdate))),
        f_int_prelim = rep(0, length(unique(dat$farm))),
        f_slope_prelim = rep(0, length(unique(dat$farm))),
        z_prelim = rep(0, length(unique(dat$zone))),
        
        sigma_Z = 0.1,
        sigma_F_int = 0.1,
        sigma_F_slope = 0.1,
        cor_F = 0,
        sigma_S = 5
    ))
}

fit = stan(
    file = "PRV_region&farm.stan",  # Stan program
    init = Stan.inits,
    data = Stan.data,    # named list of data
    chains = 3,             # number of Markov chains
    warmup = 2000,          # number of warmup iterations per chain
    iter = 4000,            # total number of iterations per chain
    cores = 3,              # number of cores 
    refresh = 100          # show progress every 'refresh' iterations
    ,thin = 1#0
    ,control = list(adapt_delta = 0.9, max_treedepth = 25)
)

saveRDS(fit,"StanFitPRV.rds") #save that fit!

##################################################################################################################
#                    diagnostic stuff to check out
##################################################################################################################
print(fit, pars = c('beta_0','beta[1]',
                    #'alpha_D','gamma_D','sigma_D',
                    'sigma_S', paste('s[',1:(as.integer(length(unique(dat$farmXdate)))),']',sep=''))
      ,digits_summary = 4)

pairs(fit, pars=c('d[1]','d[2]','alpha_D','sigma_D'))

pairs(fit, pars=c('beta_0','beta[1]',
                  'sigma_Z', 'sigma_S', 'sigma_F_int', 'sigma_F_slope', 'cor_F','s[1]'#, 'lp__'
))
pairs(fit, pars=c('z[1]','z[2]','z[3]','z[4]','z[5]','z[6]',
                  'sigma_Z'))

pairs(fit, pars=c('f[1,1]','f[2,1]','f[3,1]','f[4,1]','sigma_F_int'))
pairs(fit, pars=c('f[1,2]','f[2,2]','f[3,2]','f[4,2]','sigma_F_slope'))
i=2;pairs(fit, pars=c(paste('f[',1+i,',1]',sep=''),paste('f[',1+i,',2]',sep=''),paste('f[',2+i,',1]',sep=''),paste('f[',2+i,',2]',sep=''),
                      'sigma_F_int','sigma_F_slope','cor_F'))


pairs(fit, pars=c('sigma_S','s[1]','s[2]','s[3]'))
rstan::traceplot(fit, pars=c('s[1]',#'beta_0','beta[1]',
                  #'alpha_D',
                  #'gamma_D',
                  #'sigma_D',
                  'sigma_S'))
rstan::traceplot(fit, pars=c('sigma_S','s[1]','s[2]','s[3]','s[4]','s[5]','s[6]'))
rstan::traceplot(fit, pars=c('sigma_F_int', 'sigma_F_slope', 'cor_F')) 

hist(c(as.matrix(samples[,paste('s',c(1:252),'',sep='.')])))
pairs(fit, pars=c(names(fit)[1:10]))


##################################################################################################################
#                    model plotting
##################################################################################################################

#generate predictions, taking into account random effects
set.seed(1)
samples = rbind(as.data.frame(fit@sim$samples[[1]]), as.data.frame(fit@sim$samples[[2]]), as.data.frame(fit@sim$samples[[3]]))
seqdat = read.csv('seq_samples.csv')
rows = sample(1:dim(samples)[1],12000)
x.vals = 0:max(dat$SWdays, na.rm=TRUE)
preds = matrix(rep(NA,length(x.vals)*length(rows)), byrow = TRUE, ncol=length(x.vals))
preds.RE = preds
for(i in 1:length(rows)){
    #make the core prediction
    preds[i,] = plogis(samples[i,'beta_0'] + x.vals*samples[i,'beta.1.'])
    #randomly sample from among the farm visits, to incorporate the variation due to the random effects
    event = sample(Stan.data$sample_dates,1)
    farm = strsplit(event,'[.]')[[1]][1]
    preds.RE[i,] = plogis(samples[i,'beta_0'] #+ as.numeric(samples[i,grep('s[.]',names(samples))][which(Stan.data$sample_dates==event)])
                                      + as.numeric(samples[i,grep('f\\..+\\.1\\.',names(samples))][which(Stan.data$farm.names==farms$name[farms$ID==farm])]) + 
                              x.vals*(samples[i,'beta.1.'] 
                                      #+ as.numeric(samples[i,grep('z[.]',names(samples))][which(Stan.data$zone.names==farms$zone[farms$ID==farm])])
                                      + as.numeric(samples[i,grep('f\\..+\\.2',names(samples))][which(Stan.data$farm.names==farms$name[farms$ID==farm])])
                                      ))
}
CIs = apply(preds,2,function(x){quantile(x, probs=c(0.025,0.05,0.5,0.95,0.975))})    
CIs.RE = apply(preds.RE,2,function(x){quantile(x, probs=c(0.025,0.05,0.5,0.95,0.975))})    

# plot model, using Art's beautious colour scheme
plot(0:max(dat$SWdays, na.rm=TRUE), 
          plogis(mean(samples[,'beta_0']) + (0:max(dat$SWdays, na.rm=TRUE))*mean(samples[,'beta.1.'])),
    col='black', type='l', ylim=c(-0.025,1.025),
    xlab='days in salt water', ylab='PRV prevalence')

sampling.dates = tapply(dat$SWdays, dat$farmXdate, function(x){x[1]})
sampling.props = tapply(dat$prv, dat$farmXdate, function(x){x=x[!is.na(x)]; sum(x>0)/length(x)})
sampling.counts = tapply(dat$prv, dat$farmXdate, function(x){x=x[!is.na(x)]; length(x)})
points(sampling.dates, sampling.props, pch=16, cex=sqrt(sampling.counts)/2, 
       col='#00000055')
polygon(c(x.vals,rev(x.vals)), c(CIs[1,], rev(CIs[5,])), 
        col='#00000045', border=NA)
polygon(c(x.vals,rev(x.vals)), c(CIs.RE[1,], rev(CIs.RE[5,])), 
        col='#00000035', border=NA)





# plot model, pdf-style
pdf('PRV infection - audit samples.pdf', width=6, height=3)
par(mai=c(0.8,0.8,0.05,0.05))
plot(0:max(dat$SWdays, na.rm=TRUE), 
     plogis(mean(samples[,'beta_0']) + (0:max(dat$SWdays, na.rm=TRUE))*mean(samples[,'beta.1.'])),
     col='black', type='l', ylim=c(-0.025,1.025),
     xlab='days in salt water', ylab='PRV prevalence')
sampling.dates = tapply(dat$SWdays, dat$farmXdate, function(x){x[1]})
sampling.props = tapply(dat$prv, dat$farmXdate, function(x){x=x[!is.na(x)]; sum(x>0)/length(x)})
sampling.counts = tapply(dat$prv, dat$farmXdate, function(x){x=x[!is.na(x)]; length(x)})
points(sampling.dates, sampling.props, pch=16, cex=sqrt(sampling.counts)/2, 
       col='#00000055')
polygon(c(x.vals,rev(x.vals)), c(CIs[1,], rev(CIs[5,])), 
        col='#00000045', border=NA)
polygon(c(x.vals,rev(x.vals)), c(CIs.RE[1,], rev(CIs.RE[5,])), 
        col='#00000035', border=NA)
dev.off()



#old, perhaps useful plotting stuff

#sampling-event random effects
fitted.farms = sapply(Stan.data$sample_dates,function(x){strsplit(x,'[.]')[[1]][1]})
companies = sapply(fitted.farms,function(x){farms$company[farms$ID==x]})
hist(apply(samples[,grep('s.',names(samples))],2,mean), breaks=50)
hist(apply(samples[,grep('s.',names(samples))][companies=='Cermaq Canada'],2,mean), breaks=50)
hist(apply(samples[,grep('s.',names(samples))][companies=='MOWI Canada West'],2,mean), breaks=50)


#sampling-event random intercepts, against days in SW, by region
plot(Stan.data$sample_SWdays,
     apply(samples[,grep('s[.]',names(samples)) ],2,mean),
     col=sapply(unique(dat$farmXdate),function(x){
            farm=strsplit(x,'[.]')[[1]][1]
            which(farms$zone[farms$ID==farm]==unique(farms$zone))
         }), pch=16)


plot(dat$SWdays,log(dat$prv))
for(farm in unique(dat$farm)){
    points(dat$SWdays[dat$farm==farm], 
           log(dat$prv[dat$farm==farm]), pch=16, 
           col=rainbow(length(unique(dat$farm)))[which(farm==unique(dat$farm))])
}
