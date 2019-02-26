load('flex_reducedDat.rda')

library(BASS)

x<-X.sim[-holdout.sim,]
x<-data.frame(x[,1:6],apply(x[,7:11],2,as.factor))
ho<-sample(nrow(x),size=200) # a few more holdouts from the BASS models

n.pc<-100 # how many EOFs to use

library(parallel)
mod<-mclapply(1:n.pc,function(i) bass(x[-ho,],reduced.y.sim[i,-ho],maxInt=4,maxInt.cat=3,nmcmc=200000,nburn=190000,thin=10,h2=1e5,save.yhat=F,small=T),mc.cores = ncores,mc.preschedule=F)

save.image('flex_bassEmu.rda')

