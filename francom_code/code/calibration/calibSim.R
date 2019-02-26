################################################################################
## calibration using pre-fit Bayesian MARS emulator
################################################################################
# can calibrate to obs or one of the model runs

out.name<-'calibBassSim'

nmcmc<-40000 
nburn<-30000

cat('Calibration ');timestamp()

library(parallel)
library(BASS)
library(mnormt)


st<-Sys.time()
allowed.time<-3*24*60*60 # seconds allotted on msub
end.time<-st+allowed.time-10*60 # time to wrap things up (only 1 hour left)

## emulator prediction functions
source('code/emulator/predict_funcs.R')

################################################################################
## load emulator
load('flex_bassEmu.rda')
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
 load('flex_dat.rda')
# Y.sim: output from 18000 model runs

# use one of the holdouts to test calibration

hold.use<-7
y<-c(Y.sim[,holdout.sim[hold.use]])
for(i in 7:11)
  X.sim[,i]<-as.factor(X.sim[,i])
true.x<-X.sim[holdout.sim[hold.use],]
use<-!is.na(y)
Tr<-diag(4658)
yobs<-matrix(y,137)

################################################################################
## load observations
# load('obs.rda')
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


# true.lat<-(35.2111-35.1977)/(35.2250-35.1977)
# true.lon<-(-120.8543--120.8708)/(-120.8384--120.8708)
# true.amt<-0.5822002#136.016/(990)#0.5822002#
# true.dur<-true.st<-.5
# true.z<-1/9
# 
# true.x<-c(true.lat,true.lon,true.z,true.st,true.dur,true.amt)
# 
# y<-c(yobs)
# use<-!is.na(y)
# y<-y[use]

################################################################################
## get all 162 categorical combinations, build every categorical basis used
n.pc<-100 # number of EOFs
nmcmc.emu<-1000 # number of emulator MCMC samples
cats<-expand.grid(0:1,0:2,0:2,0:2,0:2)
cats<-data.frame(apply(cats,2,as.factor))
mb<-max(sapply(mod,function(m) max(m$nbasis)))+1
  # max number of basis functions plus intercept

bassCat<-array(raw(0),dim=c(mb,nrow(cats),n.pc,nmcmc.emu))
for(it in 1:nmcmc.emu){
  bassCat[,,,it]<-array(as.raw(bassCatBasis(mod,cats,1,it)),dim = c(mb,n.pc,162))
  print(it)
}
# save(bassCat,file='bassCatRaw.rda')

# load('bassCatRaw.rda')
# bassCat: all categorical bases used in all mod.  
#   bassCat[i,j,k,l] - ith basis function, jth combination of categories (cats[j,]),
#   kth EOF model, lth MCMC iteration.  Stored as raw to save space.

################################################################################
## useful quantities

dd<-c(mb,n.pc) # dimensions used in some predict functions
K<-eofs.sim[,1:n.pc] # eof basis matrix
mm<-0 # data mean subtracted before SVD
nc<-min(4,ncores) # number of cores to use in prediction, more than 4 doesn't really help
TrK<-Tr%*%K # transformed EOFs to be in same space as obs
TrK<-TrK[use,] # exclude values where obs are missing
Trmm<-0 # c(Tr%*%mm)
cpTrK<-crossprod(TrK)
TrKt<-t(TrK)
s.bass<-sqrt(do.call(cbind,lapply(mod,function(x) x$s2))) # matrix of bass std devs,
  # dimension nmcmc.emu x n.pc
cat.names<-apply(x[,7:11],2,function(x){sort(unique(x))})
n<-length(y)
p<-ncol(x)
p1<-6 # continuous variables
p2<-5 # discrete variables
count<-0
cc<-2.4^2/p1
eps<-1e-8
cov.theta1<-diag(p1)*eps

nt1<-apply(yobs,1,function(x) length(which(!is.na(x))))
nt<-nt1[nt1>0] # for each location, how many time observations are there
ns<-length(nt) # number of locations
space.mat<-kronecker(1:137,t(rep(1,12)))
space.vec<-c(space.mat)[use]
B<-matrix(ncol=ns,nrow=sum(nt)) # transformation from spatial locations to 
                                #  repeated spatial locations to match pattern in y
pp<-0
for(j in (1:137)[nt1>0]){
  pp<-pp+1
  B[,pp]<-as.numeric(space.vec==j) 
  # B[i,j] = 1 if y[i] is at location j.  Hence, colSums(B)==nt and rowSums(B)
  #   is a vector of ones.
}


################################################################################
## get leftover error from using only 100 eofs (spatiotemporal)
trunc.error<- Y.sim[,-holdout.sim][,-ho]-K%*%reduced.y.sim[1:n.pc,-ho] # deviations
trunc.error<-Tr[use,]%*%trunc.error 
  # same as distributing multiplication through previous line
# rm(Y.sim)

truncErrSamp<-function(te.mat){ # a function to quickly sample one of the truncation 
                                #   errors at each space and time, vectorized
  # te.mat is truncation error matrix, ncol is number of sims, nrow is length of EOFs
  #   (space time combinations)
  inds<-sample.int(ncol(te.mat),nrow(te.mat),replace=T) 
  samp<-te.mat[cbind(1:nrow(te.mat),inds)]
  return(samp)
}

# example truncation error distrubution for a certain space and time
i=900
hist(trunc.error[i,],freq=F)
lines(density(trunc.error[i,]))
curve(dnorm(x,0,sd(trunc.error[i,])),add=T,n=1000,col=2) # why we don't use Gaussian

################################################################################
## fast sampling functions

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

theta2Samp<-function(qf.stars){
  lpost<- -.5*qf.stars
  lpost<-lpost-max(lpost)
  probs<-exp(lpost)/sum(exp(lpost))
  return(sample(162,size=1,prob=probs))
}

################################################################################
## obs error

# IID ERROR

# initialize
s2<-NA;s2[1]<-1 
vv<-s2[1]

# priors
a<-0;b<-0 # priors for s2

# LOCATION VARYING ERROR

# # initialize
# s2<-NA;s2[1]<-1 
# s2l<-matrix(nrow=nmcmc,ncol=ns) # locations dependent discrepancy
# s2l[1,]<-1
# vv<-c(B%*%s2l[1,]) # this just gets the ordering and repeats of s2l to match y
# 
# # priors
# a<-1;b<-.1 # priors for s2
# a.s2l<-1000 # hyperparm for location dependent variances - make this large to keep
#             #   these from getting too spread out

################################################################################
## discrepancy

# NONE

delta<-0

# EOF SPACE

# # initialize
# d<-matrix(nrow=nmcmc-nburn,ncol=n.pc); d.curr<-rep(0,n.pc) # discrepancy
# s2d<-NA; s2d[1]<-.00001
# delta<-c(TrK%*%d.curr)
# 
# # priors
# ad<-1/.001;bd<-.00001 # priors for s2d

# MODULARIZED - BASS

# x.locs<-B%*%as.matrix(coordsLL[nt1>0,])
# x.times<-c(kronecker(rep(1,137),t(times.obs)))[use]
# 
# e.prior<-emuSample(rep(.5,6),which(apply(cats,1,function(x) all(x==c(0,1,1,1,0)))))
# d.prior<-y-e.prior$eta.mean
# dmod.prior<-bass(cbind(x.locs,x.times),d.prior)
# dmod.s<-sqrt(mean(dmod.prior$s2))
# dmod.m<-dmod.prior$yhat.mean
# itd<-sample(1000,size=1)
# delta<-dmod.prior$yhat[itd,]#rnorm(n,dmod.m,dmod.s)

# with gamma
gam<-NA
gam[1]<-1
library(truncnorm)


################################################################################
## initializations
theta1<-matrix(nrow=nmcmc,ncol=p1);theta1[1,]<-runif(p1) # continuous parameters
theta2.ind<-rep(integer(0),nmcmc);theta2.ind[1]<-sample(162,size=1) 
# categorical parameters (only 162 possible combinations, keep track of the 
#   combination rather than the vector)

system.time(emuSamp.curr<-emuSample(theta1[1,],theta2.ind[1]))
rr<-y-emuSamp.curr$eta-delta
qf<-sum(rr^2/vv)
  
eta<-matrix(nrow=nmcmc-nburn,ncol=n)







################################################################################
################################################################################
#MCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMC
cat('start MCMC ');timestamp()
for(i in 2:nmcmc){
  
  #MCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMC
  ## get a sample from emulator posterior, to propagate uncertainty (monte carlo)
  emuSamp.curr<-emuSample(theta1[i-1,],theta2.ind[i-1])
  rr<-y-emuSamp.curr$eta-delta*gam[i-1]
  qf<-sum(rr^2/vv)
  
  #MCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMC
  ## update obs error
  
  # IID ERROR
  
  # sample s2
  s2[i]<-1/rgamma(1,n/2+a,b+.5*sum(rr^2))
  vv<-s2[i]
  qf<-sum(rr^2/vv)
  
  # LOCATION VARYING ERROR
  
  # # sample location-varying s2l
  # s2l[i,]<-1/rgamma(ns,nt/2+a.s2l+1,a.s2l*s2[i-1]+c(t(rr^2)%*%B)/2) 
  #   # t(rr^2)%*%B gets sum of squares for each location
  # vv<-c(B%*%s2l[i,])
  # qf<-sum(rr^2/vv)
  # 
  # # sample s2
  # s2[i]<-rgamma(1,ns*(a.s2l+1)+a,rate=1/b+a.s2l*sum(1/s2l[i,]))
  
  #MCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMC
  ## update theta1 - adaptive M-H
  theta1[i,]<-theta1[i-1,]
  if(i>200){
    mi<-max(1,i-300)
    cov.theta1<-cov(theta1[mi:(i-1),])*cc+diag(eps*cc,p1)
  }
  
  theta1.star<-rmnorm(1,theta1[i-1,],cov.theta1) # generate candidate
  if(any(theta1.star<0 | theta1.star>1)){ # constraint
    alpha<- -9999
  } else{
    emuSamp.star<-emuSample(theta1.star,theta2.ind[i-1])
    rr.star<-y-emuSamp.star$eta-delta*gam[i-1]
    qf.star<-sum(rr.star^2/vv)
    alpha<- -.5*(qf.star-qf)
  }
  if(log(runif(1))<alpha){
    theta1[i,]<-theta1.star
    qf<-qf.star
    rr<-rr.star
    emuSamp.curr<-emuSamp.star
    count<-count+1
  }
  
  #MCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMC
  ## update theta2 - brute force sampling
  all.eta.star.mean<-getCatPreds(emuSamp.curr$it,TrK,Trmm,emuSamp.curr$desBasis,
                                 emuSamp.curr$betaMat,nc)
  rr.stars<-y-all.eta.star.mean-c(TrK%*%rnorm(n.pc,0,s.bass[emuSamp.curr$it,]))-
    truncErrSamp(trunc.error)-delta*gam[i-1]
  qf.stars<-colSums(rr.stars^2/vv)
  theta2.ind[i]<-theta2Samp(qf.stars)
  
  emuSamp.curr$eta.mean<-all.eta.star.mean[,theta2.ind[i]]
  emuSamp.curr$eta<-emuSamp.curr$eta.mean+TrK%*%rnorm(n.pc,0,s.bass[emuSamp.curr$it,])+
    truncErrSamp(trunc.error)
  
  #MCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMC
  ## update discrepancy
  
  # None
  
  # EOF SPACE
  
  # tr<-TrKt%*%diag(1/vv)%*%TrK
  # tr[upper.tri(tr)]<-t(tr)[upper.tri(tr)]
  # Sd<-pd.solve(tr+diag(s2d[i-1],n.pc))#pd.solve(cpTrKt/s2[i]+diag(s2d[i-1],n.pc))
  # md<-TrKt%*%diag(1/vv)%*%(y-eta.curr)
  # d.curr<-rmnorm(1,Sd%*%md,Sd)
  # delta<-c(TrK%*%d.curr)
  # 
  # s2d[i]<-1/rgamma(1,n.pc/2+ad,bd+.5*crossprod(d.curr))
  
  # MODULARIZED BASS
  
  delta<-0#dmod.prior$yhat[itd,]#rnorm(n,dmod.m,dmod.s)
  gam[i]<-1#rtruncnorm(1,0,2,sum(delta*(y-emuSamp.curr$eta))/sum(delta^2),s2[i]/sum(delta^2))
  
  
  #MCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMC
  if(i>nburn){
    eta[i-nburn,]<-emuSamp.curr$eta.mean
    #d[i-nburn,]<-d.curr
  }
  
  if(i%%100==0)
    cat('it:',i,'acc:',count,timestamp(quiet=T),'\n')

  if(Sys.time()>end.time) # end loop early
    break
}

#MCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMCMC
################################################################################
################################################################################
cat('finish MCMC ');timestamp()
cat('count: ',count)


save.image(paste(out.name,'.rda',sep=''))
#burn<-1:nburn
burn<-c(1:10000,seq(10001,40000,3)) # include thinning
#X<-x

