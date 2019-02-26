## functions for EOF prediction


predBassEOF<-function(mod.list,x,V,mm,nc=1,mcmc.use=1){
  # mod.list: a list of bass models, order corresponding to order of EOFs
  # x: vector or matrix of inputs (vector for 1 prediction, matrix for multiple predictions)
  # V: matrix with columns as EOFs of centered model runs
  # mm: mean of model runs
  # nc: number of cores
  # mcmc.use: with mcmc iterations to use
  # output is an array or matrix: dimensions are size of data x mcmc iterations x number of predictions - if only one prediction (x a vector) this turns into a matrix
  require(parallel)
  weights<-matrix(unlist(mclapply(mod.list,function(m) predict(m,x,mcmc.use=mcmc.use),mc.cores=nc,mc.preschedule=T)),nrow=n.pc,byrow=T)
  out<-V%*%weights+mm
  out<-array(out,dim=c(nrow(out),length(mcmc.use),nrow(x)))
  return(out)
}

cbind.fill0 <- function(...,n=mb){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  if(is.null(n))
    n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(0, n-nrow(x), ncol(x))))) 
}
bassDesBasis<-function(mod.list,xdes,nc=1,mcmc.use){
  # mod.list: a list of bass models, order corresponding to order of EOFs
  # xdes: vector
  # nc: number of cores
  # mcmc.use: integer
  # output is a matrix (n.pc x max(nbasis)+1) that represents the continuous portion of the basis functions evaluated at xdes for each model in mod.list, fill in 0s
  for(i in 1:length(xdes))
    xdes[i]<-BASS:::scale.range(xdes[i],mod.list[[1]]$range.des[,i])
  ll<-mclapply(mod.list,function(m) BASS:::makeBasisMatrix(m$model.lookup[mcmc.use],m$nbasis[mcmc.use],m$vars.des,m$signs.des,m$knotInd.des,m$degree,matrix(xdes),m$n.int.des,m$xx.des),mc.cores=nc,mc.preschedule=T)
  out<-matrix(0,nrow=mb,ncol=n.pc)
  for(i in 1:n.pc){
    lli<-ll[[i]][,1]
    out[1:length(lli),i]<-lli
  }
  #browser()
  return(out)
  #return(do.call(cbind.fill0,ll))
}

bassCatBasis<-function(mod.list,xcat,nc=1,mcmc.use,mb=mb){ # why not just precompute all these
  # mod.list: a list of bass models, order corresponding to order of EOFs
  # xcat: data.frame of factors
  # nc: number of cores
  # mcmc.use: integer
  # output is a matrix (n.pc x max(nbasis)+1) that represents the categorical portion of the basis functions evaluated at xcat for each model in mod.list, fill in 0s
  ll<-mclapply(mod.list,function(m) BASS:::makeBasisMatrixCat(m$model.lookup[mcmc.use],m$nbasis[mcmc.use],m$vars.cat,xcat,m$n.int.cat,m$sub.list),mc.cores=nc,mc.preschedule=T)
  return(do.call(cbind.fill0,ll))
}
bassBeta<-function(mod.list,mcmc.use,mb=mb){
  temp<-lapply(mod.list,function(m) m$beta[mcmc.use,1:mb])#,mc.cores=nc,mc.preschedule=T)
  ret<-matrix(0,nrow=mb,ncol=n.pc)
  for(ii in 1:n.pc){
    ret[,ii]<-temp[[ii]]
  }
  
  #mat<-do.call(rbind,mclapply(mod.list,function(m) m$beta[mcmc.use,]))
  #mnb<-max(sapply(mod.list,function(m) m$nbasis[mcmc.use]))
  #ret<-mat[,1:mb]
  ret[is.na(ret)]<-0
  #browser()
  return(ret)
}


weightsEOF<-function(desBasis,catBasis,betaMat){
  colSums(desBasis*catBasis*betaMat)
}
pred1BassEOF<-function(mod,x,V,mm,nc,it){
  xdes<-as.numeric(x[1:6])
  xcat<-x[7:11]
  betaMat<-bassBeta(mod,it,mb=mb)
  desBasis<-bassDesBasis(mod,xdes,nc,it)
  #catBasis<-bassCatBasis(mod,xcat,ncores,it)
  cc<-which(c(apply(cats,1,function(xc) all(xc==xcat))))
  catBasis<-raw2num(bassCat[,cc,,it])
  weights<-weightsEOF(desBasis,catBasis,betaMat)
  pred<-V%*%weights+mm
  return(pred)
}

raw2num<-function(x,dd=NULL){ # x can be a matrix
  if(is.null(dd))
    dd<-dim(x)
  x<-as.numeric(x)
  dim(x)<-dd
  return(x)
}

getCatPreds<-function(it,V,mm,desBasis,betaMat,nc){
  pred<-mclapply(1:162,function(i){
    c(V%*%weightsEOF(desBasis,raw2num(bassCat[,i,,it]),betaMat)+mm)
  },mc.cores=nc,mc.preschedule=T)
  return(do.call(cbind,pred))
}