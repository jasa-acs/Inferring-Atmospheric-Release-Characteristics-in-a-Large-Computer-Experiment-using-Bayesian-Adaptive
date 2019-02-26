set.seed(0)

library(tgp)
library(BASS)
library(BART)
library(randomForest)
library(parallel)


load('flex_reducedDat.rda')
X.sim<-X.sim[-holdout.sim,] # these are left out of the EOFs
X.sim.factor<-X.sim
for(i in 7:11)
  X.sim.factor[,i]<-as.factor(X.sim[,i])

X.sim.bin<-cbind(X.sim[,1:6],model.matrix(lm(X.sim[,1]~.,data=X.sim.factor[,7:11]))[,-1]) # turn factors into binary design matrix

# BASS and randomForest use factors, tgp and BART require factors be broken into binaries

rmse<-function(x,y) sqrt(mean((x-y)^2))

compare<-function(pc,long=F){
  train<-sample(1:nrow(X.sim),size=ntrain,replace=F) # take a random subset (ntrain,ntest,reduced.y.sim,X.sim,... passed from global environment)
  test<-sample((1:nrow(X.sim))[-train],size=ntest,replace=F) # taking a different random subset
  y <- reduced.y.sim[pc,]
  
  ################################################################################################
  ## TGP
  
  if(!long){
    ptm <- proc.time()
    mod.tgp<-btgpllm(X = X.sim.bin[train,], Z = y[train], XX = X.sim.bin[test,],basemax=6) # tried different splitmins and tree priors, this was best
    time.tgp<-(proc.time() - ptm)[3] # time taken
    cover.tgp<-mean(mod.tgp$ZZ.q1 < y[test] & mod.tgp$ZZ.q2 > y[test]) # empirical coverage
    length.tgp<-mean(mod.tgp$ZZ.q2-mod.tgp$ZZ.q1) # average interval length
    pred.mean.tgp<-mod.tgp$ZZ.mean # posterior predictive means
  }
  
  ################################################################################################
  ## BMARS
  
  ptm <- proc.time()
  mod.bass<-bass(X.sim.factor[train,],y[train])
  time.bass<-(proc.time() - ptm)[3] # time taken
  pred.bass<-predict(mod.bass,X.sim.factor[test,]) # posterior predictive samples (without epsilon): rows are mcmc samples, cols are data points
  eps<-matrix(rnorm(nrow(pred.bass)*ncol(pred.bass),0,sqrt(mod.bass$s2)),nrow=nrow(pred.bass)) # epsilon samples that mimic shape of pred.bass; for example, see matrix(rnorm(5*5,1:5,.0001),nrow=5)
  quants<-apply(pred.bass+eps,2,quantile,probs=c(.05,.95)) # posterior predictive quantiles
  cover.bass<-mean(quants[1,] < y[test] & quants[2,] > y[test]) # empirical coverage
  length.bass<-mean(quants[2,]-quants[1,]) # average interval length
  pred.mean.bass<-colMeans(pred.bass) # posterior predictive means
  
  ################################################################################################
  ## BART
  
  ptm <- proc.time()
  mod.bart<-wbart(X.sim.bin[train,],y[train])
  time.bart<-(proc.time() - ptm)[3] # time taken
  pred.bart<-predict(mod.bart,X.sim.bin[test,]) # posterior predictive samples (without epsilon): rows are mcmc samples, cols are data points
  eps<-matrix(rnorm(nrow(pred.bart)*ncol(pred.bart),0,mod.bart$sigma),nrow=nrow(pred.bart)) # epsilon samples that mimic shape of pred.bart
  quants<-apply(pred.bart+eps,2,quantile,probs=c(.05,.95)) # posterior predictive quantiles
  cover.bart<-mean(quants[1,] < y[test] & quants[2,] > y[test]) # empirical coverage
  length.bart<-mean(quants[2,]-quants[1,]) # average interval length
  pred.mean.bart<-colMeans(pred.bart) # posterior predictive means

  ################################################################################################
  ## RF
  
  ptm <- proc.time()
  mod.rf<-randomForest(X.sim.factor[train,],y[train])
  time.rf<-(proc.time() - ptm)[3]
  pred.mean.rf<-predict(mod.rf,X.sim.factor[test,]) # predictions
  
  ################################################################################################
  
  if(long){
    return(list(rmses=c(NA,rmse(pred.mean.bass,y[test]),rmse(pred.mean.bart,y[test]),rmse(pred.mean.rf,y[test])),times=c(NA,time.bass,time.bart,time.rf),cover=c(NA,cover.bass,cover.bart,NA),length=c(NA,length.bass,length.bart,NA)))
  } else{
    return(list(rmses=c(rmse(pred.mean.tgp,y[test]),rmse(pred.mean.bass,y[test]),rmse(pred.mean.bart,y[test]),rmse(pred.mean.rf,y[test])),times=c(time.tgp,time.bass,time.bart,time.rf),cover=c(cover.tgp,cover.bass,cover.bart,NA),length=c(length.tgp,length.bass,length.bart,NA)))
  }
  

}

process.compare<-function(cc){
  rmse<-time<-cover<-length<-matrix(nrow=length(cc),ncol=4)
  for(i in 1:length(cc)){
    if(class(cc[[i]])=='try-error')
      next
    rmse[i,]<-cc[[i]]$rmses
    time[i,]<-cc[[i]]$time
    cover[i,]<-cc[[i]]$cover
    length[i,]<-cc[[i]]$length
  }
  return(list(rmse=rmse,time=time,cover=cover,length=length))
}



nsim<-max(ncores,40)
pc<-1

ntrain<-50
ntest<-1000
out.small<-mclapply(rep(pc,nsim),compare,mc.cores=ncores)
res.small<-process.compare(out.small)

ntrain<-500
ntest<-1000
out.med<-mclapply(rep(pc,nsim),compare,mc.cores=ncores)
res.med<-process.compare(out.med)

ntrain<-5000
ntest<-1000
out.big<-mclapply(rep(pc,nsim),compare,long=T,mc.cores=ncores)
res.big<-process.compare(out.big)



pdf('emu_pc1_sim.pdf',height=8,width=8)

rr.rmse<-log10(range(res.small$rmse,res.med$rmse,res.big$rmse,na.rm = T))
rr.time<-log10(range(res.small$time,res.med$time,res.big$time,na.rm = T))
rr.cover<-range(res.small$cover,res.med$cover,res.big$cover,.9,na.rm = T)
rr.length<-range(res.small$length,res.med$length,res.big$length,na.rm = T)
par(mfrow=c(4,3),oma=c(0,5,0,5),mar=c(2,0,2.5,0))

boxplot(log10(res.small$rmse),ylim=rr.rmse,xaxt='n')
abline(h=pretty(rr.rmse),col='lightgrey')
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART','RF'))
mtext(expression(log[10]~'RMSE'),side = 2,line=3)
mtext('training n = 50',side=3,line=1)
boxplot(log10(res.med$rmse),ylim=rr.rmse,yaxt='n',xaxt='n')
abline(h=pretty(rr.rmse),col='lightgrey')
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART','RF'))
mtext('n = 500',side=3,line=1)
boxplot(log10(res.big$rmse),ylim=rr.rmse,yaxt='n',xaxt='n')
abline(h=pretty(rr.rmse),col='lightgrey')
axis(side = 1,at = 1:4,labels = c('','BMARS','BART','RF'))
mtext('n = 5000',side=3,line=1)

boxplot(res.small$cover,ylim=rr.cover,xaxt='n')
abline(h=pretty(rr.cover),col='lightgrey')
abline(h=.9,col=2,lty=2)
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART',''))
mtext('coverage',side = 2,line=3)
boxplot(res.med$cover,ylim=rr.cover,yaxt='n',xaxt='n')
abline(h=pretty(rr.cover),col='lightgrey')
abline(h=.9,col=2,lty=2)
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART',''))
boxplot(res.big$cover,ylim=rr.cover,yaxt='n',xaxt='n')
abline(h=pretty(rr.cover),col='lightgrey')
abline(h=.9,col=2,lty=2)
axis(side = 1,at = 1:4,labels = c('','BMARS','BART',''))

boxplot(res.small$length,ylim=rr.length,xaxt='n')
abline(h=pretty(rr.length),col='lightgrey')
abline(h=.9,col=2,lty=2)
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART',''))
mtext('interval length',side = 2,line=3)
boxplot(res.med$length,ylim=rr.length,yaxt='n',xaxt='n')
abline(h=pretty(rr.length),col='lightgrey')
abline(h=.9,col=2,lty=2)
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART',''))
boxplot(res.big$length,ylim=rr.length,yaxt='n',xaxt='n')
abline(h=pretty(rr.length),col='lightgrey')
abline(h=.9,col=2,lty=2)
axis(side = 1,at = 1:4,labels = c('','BMARS','BART',''))

boxplot(log10(res.small$time),ylim=rr.time,xaxt='n')
abline(h=pretty(rr.time),col='lightgrey')
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART',''))
mtext(expression(log[10]~'time (seconds)'),side = 2,line=3)
boxplot(log10(res.med$time),ylim=rr.time,yaxt='n',xaxt='n')
abline(h=pretty(rr.time),col='lightgrey')
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART','RF'))
boxplot(log10(res.big$time),ylim=rr.time,yaxt='n',xaxt='n')
abline(h=pretty(rr.time),col='lightgrey')
axis(side = 1,at = 1:4,labels = c('','BMARS','BART','RF'))

dev.off()

save.image('emu_pc1_sim.rda')









pc<-50

ntrain<-50
ntest<-1000
out.small<-mclapply(rep(pc,nsim),compare,mc.cores=ncores)
res.small<-process.compare(out.small)

ntrain<-500
ntest<-1000
out.med<-mclapply(rep(pc,nsim),compare,mc.cores=ncores)
res.med<-process.compare(out.med)

ntrain<-5000
ntest<-1000
out.big<-mclapply(rep(pc,nsim),compare,long=T,mc.cores=ncores)
res.big<-process.compare(out.big)



pdf('emu_pc50_sim.pdf',height=8,width=8)

rr.rmse<-log10(range(res.small$rmse,res.med$rmse,res.big$rmse,na.rm = T))
rr.time<-log10(range(res.small$time,res.med$time,res.big$time,na.rm = T))
rr.cover<-range(res.small$cover,res.med$cover,res.big$cover,.9,na.rm = T)
rr.length<-range(res.small$length,res.med$length,res.big$length,na.rm = T)
par(mfrow=c(4,3),oma=c(0,5,0,5),mar=c(2,0,2.5,0))

boxplot(log10(res.small$rmse),ylim=rr.rmse,xaxt='n')
abline(h=pretty(rr.rmse),col='lightgrey')
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART','RF'))
mtext(expression(log[10]~'RMSE'),side = 2,line=3)
mtext('training n = 50',side=3,line=1)
boxplot(log10(res.med$rmse),ylim=rr.rmse,yaxt='n',xaxt='n')
abline(h=pretty(rr.rmse),col='lightgrey')
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART','RF'))
mtext('n = 500',side=3,line=1)
boxplot(log10(res.big$rmse),ylim=rr.rmse,yaxt='n',xaxt='n')
abline(h=pretty(rr.rmse),col='lightgrey')
axis(side = 1,at = 1:4,labels = c('','BMARS','BART','RF'))
mtext('n = 5000',side=3,line=1)

boxplot(res.small$cover,ylim=rr.cover,xaxt='n')
abline(h=pretty(rr.cover),col='lightgrey')
abline(h=.9,col=2,lty=2)
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART',''))
mtext('coverage',side = 2,line=3)
boxplot(res.med$cover,ylim=rr.cover,yaxt='n',xaxt='n')
abline(h=pretty(rr.cover),col='lightgrey')
abline(h=.9,col=2,lty=2)
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART',''))
boxplot(res.big$cover,ylim=rr.cover,yaxt='n',xaxt='n')
abline(h=pretty(rr.cover),col='lightgrey')
abline(h=.9,col=2,lty=2)
axis(side = 1,at = 1:4,labels = c('','BMARS','BART',''))

boxplot(res.small$length,ylim=rr.length,xaxt='n')
abline(h=pretty(rr.length),col='lightgrey')
abline(h=.9,col=2,lty=2)
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART',''))
mtext('interval length',side = 2,line=3)
boxplot(res.med$length,ylim=rr.length,yaxt='n',xaxt='n')
abline(h=pretty(rr.length),col='lightgrey')
abline(h=.9,col=2,lty=2)
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART',''))
boxplot(res.big$length,ylim=rr.length,yaxt='n',xaxt='n')
abline(h=pretty(rr.length),col='lightgrey')
abline(h=.9,col=2,lty=2)
axis(side = 1,at = 1:4,labels = c('','BMARS','BART',''))

boxplot(log10(res.small$time),ylim=rr.time,xaxt='n')
abline(h=pretty(rr.time),col='lightgrey')
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART',''))
mtext(expression(log[10]~'time (seconds)'),side = 2,line=3)
boxplot(log10(res.med$time),ylim=rr.time,yaxt='n',xaxt='n')
abline(h=pretty(rr.time),col='lightgrey')
axis(side = 1,at = 1:4,labels = c('TGP','BMARS','BART','RF'))
boxplot(log10(res.big$time),ylim=rr.time,yaxt='n',xaxt='n')
abline(h=pretty(rr.time),col='lightgrey')
axis(side = 1,at = 1:4,labels = c('','BMARS','BART','RF'))

dev.off()

save.image('emu_pc50_sim.rda')
