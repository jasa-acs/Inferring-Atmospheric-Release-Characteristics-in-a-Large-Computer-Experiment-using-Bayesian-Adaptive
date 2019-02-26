library(parallel)
library(mnormt)

out.name<-'ep'

source('code/emulator/predict_funcs.R')
# load('emu/flex_bassEmu.rda')
# mod: list of BASS models, one for each eof
# holdout.sim: index of 1800 (of 18000) holdouts from the entire process
# ho: index of 200 holdouts from BASS modeling
# X.sim: original input data
# x: X.sim without holdout.sim
# reduced.y.sim: from svd, the EOF weights (training data for BASS)
# eofs.sim: EOFs from svd

# for(i in 1:100){ # in case bass was run with small=F
#   mod[[i]]$curr.list<-mod[[i]]$des.basis<-mod[[i]]$cat.basis<-NULL
# }


################################################################################
## load model runs
# load('~/Google Drive/MARS/applications/flex/calib/dat/flex_dat.rda')
# Y.sim: output from 18000 model runs



################################################################################
## load observations
load('obs.rda')
# yobs: matrix of data for 137 locations and 12 time points (some missing)
# sites.obs: 150 sites used in obs
# sites.modelRuns: 137 sites used in model runs
# sites.obs.use: indes so that sites.obs[sites.obs.use] matches sites.modelRuns
# coordsLL: lat/lon coords for the 137 locations
# times.modelRuns: local times for 34 model run outputs for each location
# times.obs: same for obs, only 12
# Tr: transformation from model run time scale to obs time scale (average half-
#   hourly measurements to get hourly measurements) - note that 34/12 == 4658/1644
# tTr: transformation for times directly


labels <- sapply(1:4,function(i)
  as.expression(bquote(10^ .(i)))
)

sites.plot<-c(325,333,340,343,414,412,423,410,426,225,108,526,427,510,501)
ord<-order(sites.plot)
sites.plot<-sites.plot[ord]
which(sites.obs%in%sites.plot)
site.names<-c('1 mile from plant','4 miles from plant','Whalers Island','Avila Beach','mouth of See Canyon','midway up See Canyon','1 mile inland along Hwy 101','coastline extent of plume','inland from coastline extent','San Luis Obisbo','Morro Bay','Santa Maria','Pismo Beach','Arroyo Grande','Edna')[ord]



n.eofs<-100
eofs<-eofs.sim[,1:n.eofs]

## predict function
pred.em.bmars<-function(x,eofs,nc=1,mcmc.use=1){
  require(parallel)
  weights<-do.call(rbind,(mclapply(mod,function(i) predict(i,x,mcmc.use=mcmc.use),mc.cores = nc,mc.preschedule = T)))
  return(eofs%*%weights)
}

nc<-1
emuSample<-function(t1,t2){
  it<-sample(1:nmcmc.emu,1)
  betaMat<-bassBeta(mod,it,mb=mb)
  desBasis<-bassDesBasis(mod,t1,nc,it)
  catBasis<-raw2num(bassCat[,t2,,it],dd)
  weights<-weightsEOF(desBasis,catBasis,betaMat)
  eta.mean<-TrK%*%weights+Trmm
  eta<-eta.mean+TrK%*%rnorm(n.pc,0,s.bass[it,])+truncErrSamp(trunc.error) 
  # includes variance from BASS and from truncation
  return(list(it=it,betaMat=betaMat,desBasis=desBasis,catBasis=catBasis,
              weights=weights,eta.mean=eta.mean,eta=eta))
}

################################################################################
## get all 162 categorical combinations, build every categorical basis used

nmcmc.emu<-1000 # number of emulator MCMC samples
cats<-expand.grid(0:1,0:2,0:2,0:2,0:2)
cats<-data.frame(apply(cats,2,as.factor))
mb<-max(sapply(mod,function(m) max(m$nbasis)))+1 
# max number of basis functions plus intercept
dd<-c(mb,n.pc) # dimensions used in some predict functions
s.bass<-sqrt(do.call(cbind,lapply(mod,function(x) x$s2))) # matrix of bass std devs,

bassCat<-array(raw(0),dim=c(mb,nrow(cats),n.pc,nmcmc.emu))
for(it in 1:nmcmc.emu){
  bassCat[,,,it]<-array(as.raw(bassCatBasis(mod,cats,1,it)),dim = c(mb,n.pc,162))
  print(it)
}
# save(bassCat,file='bassCatRaw.rda')

# load('bassCatRaw.rda')


## get leftover error from using only 100 eofs (spatiotemporal)
TrK<-eofs.sim[,1:n.pc]
Trmm<-0
trunc.error<- Y.sim[,-holdout.sim][,-ho]-TrK%*%reduced.y.sim[1:n.pc,-ho] # deviations
#trunc.error<-trunc.error 
# same as distributing multiplication through previous line


truncErrSamp<-function(te.mat){ # a function to quickly sample one of the truncation 
  #   errors at each space and time, vectorized
  # te.mat is truncation error matrix, ncol is number of sims, nrow is length of EOFs
  #   (space time combinations)
  inds<-sample.int(ncol(te.mat),nrow(te.mat),replace=T) 
  samp<-te.mat[cbind(1:nrow(te.mat),inds)]
  return(samp)
}



for(i in 7:11)
  X.sim[,i]<-as.factor(X.sim[,i])

# pdf('flex_emu_fit4.pdf',width=10,height=6)
pdf(paste(out.name,'EmuPlots.pdf',sep=''),width=10,height=6)
par(mfrow=c(3,5),mar=c(0,0,0,0),oma=c(5.5,5.5,2,2))
locs<-which(sites.modelRuns%in%sites.plot)
times.disp<-paste(rep(c('05','06','07','08','09',10:21),each=2),':',rep(c('00','30'),17),sep='')

#YY.use<-YY.use+1#140 # choose one of the holdouts to calibrate to
mat<-matrix(1:4658,nrow=137)
ind<-ind<-seq(1,1800,100)+1
mcmc.use<-seq(1,1000,1)
pred<-array(dim=c(137,34,length(mcmc.use)))
k2<-0
aa<-aa2<-matrix(nrow=137*34,ncol=length(mcmc.use))
for(YY.use in ind){
  YY.use<-YY.use+1
  k2<-k2+1
  y<-Y.sim[,holdout.sim[YY.use]]
  true.x<-X.sim[holdout.sim[YY.use],]
  #plot(y,pred.em.bmars(x=true.x,V,n.pc,mm,nc,1,X));abline(a=0,b=1,col=2)
  #aa<-pred.em.bmars(x=true.x,eofs,nc,mcmc.use)
  k3=0
  for(i in mcmc.use){
    k3=k3+1
    temp<-emuSample(as.numeric(true.x[1:6]),which(apply(cats,1,function(x) all(x==true.x[7:11]))))
    aa[,k3]<-temp$eta
    aa2[,k3]<-temp$eta.mean
  }
  #plot(y,rowMeans(aa));abline(a=0,b=1,col=2)
  
  #matplot(t(matrix(y,nrow=137)),t(matrix(rowMeans(aa),nrow=137)),type='l');abline(a=0,b=1,col=2)
  
  #matplot(t(dat[,,tests[YY.use]]),type='l')
  
  pred[,,k2]<-array(y-rowMeans(aa),dim=c(137,34))
  qq<-apply(aa,1,quantile,probs=c(.025,.975))
  mm<-apply(aa2,1,mean)
  
  k<-0
  for(ll in locs){
    k<-k+1
    
    plot(Y.sim[mat[ll,],holdout.sim[YY.use]],type='n',ylim=c(1,5),xaxt='n',yaxt='n')
    #lines(pred[ll,,YY.use],col=4,lty=3,lwd=2)
    kk<-1:length(mat[ll,])
    polygon(c(kk,rev(kk)),c(qq[1,mat[ll,]],rev(qq[2,mat[ll,]])),col='lightgrey',border = F)
    #matplot(aa[mat[ll,],],type='n',add=T,col=2)
    lines(Y.sim[mat[ll,],holdout.sim[YY.use]],col=4,lwd=2,lty=3)
    lines(mm[mat[ll,]],col='darkgrey',lwd=1,lty=5)
    
    #lines(qq[1,mat[ll,]],col=2,lwd=1,lty=2)
    #lines(qq[2,mat[ll,]],col=2,lwd=1,lty=2)
    #matplot(pred[ll,,]-dat[ll,,tests[ind]],type='l',ylim=c(-1.75,1.75),xaxt='n',yaxt='n')
    text(15,4.5,site.names[k])
    if(k %in% c(1,6,11))
      axis(2,at=1:5,labels = c(labels,NA),las=1)
    if(k %in% 11:15)
      axis(1,at = seq(1,34,4),labels = times.disp[seq(1,34,4)],las=2)
    if(k==1)
      legend(c(1,3),c('simulator','surrogate'),col=c(4,'darkgrey'),lty=c(3,5),lwd=c(2,1),bty='n')
  }
  mtext('time',1,outer = T,line=4)
  mtext('concentration',2,outer = T,line=4)
}

# dev.off()



k<-0
for(ll in locs){
  k<-k+1
  #plot(dat[ll,,tests[YY.use]],type='l',ylim=c(1,5),xaxt='n',yaxt='n')
  matplot(pred[ll,,],type='l',ylim=c(-1.75,1.75),xaxt='n',yaxt='n')
  text(15,1.5,site.names[k])
  if(k %in% c(1,6,11))
    axis(2)#,at=1:5,labels = c(labels,NA),las=1)
  if(k %in% 11:15)
    axis(1,at = seq(1,34,4),labels = times.disp[seq(1,34,4)],las=2)
}
mtext('time',1,outer = T,line=4)
mtext('residual',2,outer = T,line=4)


dev.off()
