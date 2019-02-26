sob.eof<-function(pc.mod,int.order,pcs,mcmc.use=1,nind=NULL,ncores=1,preschedule=T){
  
  cat('Start',timestamp(quiet = T),'\n')
  p<-pc.mod[[1]]$p
  u.list<-lapply(1:int.order,function(i) combn(1:p,i))
  ncombs.vec<-unlist(lapply(u.list,ncol))
  ncombs<-sum(ncombs.vec)
  nxfunc<-nrow(pcs)
  #sob<-matrix(nrow=nxfunc,ncol=ncombs)
  sob<-ints<-list()
  
  n.pc<-ncol(pcs)
  w0<-unlist(lapply(1:n.pc,function(pc) get.f0(pc.mod,pc,mcmc.use)))
  
  #browser()
  
  f0r2<-(pcs%*%w0)^2
  
  max.nbasis<-max(unlist(lapply(pc.mod,function(x) x$nbasis[mcmc.use])))
  C1.array<-array(dim=c(n.pc,p,max.nbasis))
  for(i in 1:n.pc){
    nb<-pc.mod[[i]]$nbasis[mcmc.use]
    mcmc.mod.usei<-pc.mod[[i]]$model.lookup[mcmc.use]
    for(j in 1:p){
      for(k in 1:nb){
        C1.array[i,j,k]<-C1(pc.mod,j,k,i,mcmc.mod.usei)
      }
    }
    #print(i)
  }
  
  # browser()
  # 
  # C2.array<-array(dim=c(n.pc,n.pc,p,max.nbasis,max.nbasis))
  # for(i1 in 1:n.pc){
  #   nb1<-pc.mod[[i1]]$nbasis[mcmc.use]
  #   mcmc.mod.usei1<-pc.mod[[i1]]$model.lookup[mcmc.use]
  #   for(i2 in 1:n.pc){
  #     nb2<-pc.mod[[i2]]$nbasis[mcmc.use]
  #     mcmc.mod.usei2<-pc.mod[[i2]]$model.lookup[mcmc.use]
  #     for(j in 1:p){
  #       for(k1 in 1:nb1){
  #         for(k2 in 1:nb2){
  #           C2.array[i1,i2,j,k1,k2]<-C2(pc.mod,j,k1,k2,i1,i2,mcmc.mod.usei1,mcmc.mod.usei2) #C2(pc.mod,l,mi,mj,i,j,mcmc.mod.usei,mcmc.mod.usej)
  #         }
  #       }
  #     }
  #   }
  #   print(i1)
  # }
  
  
  #browser()
  u.list1<-list()
  for(i in 1:int.order)
    u.list1<-c(u.list1,split(u.list[[i]], col(u.list[[i]])))
  require(parallel)
  #browser()
  cat('Integrating',timestamp(quiet = T),'\n')
  ints1<-mclapply(u.list1,function(x) func.hat(x,pc.mod,pcs,mcmc.use,f0r2,C1.array),mc.cores=ncores,mc.preschedule = preschedule)
  ints<-list()
  ints[[1]]<-do.call(cbind,ints1[1:ncol(u.list[[1]])])
  if(int.order>1){
    for(i in 2:int.order)
      ints[[i]]<-do.call(cbind,ints1[sum(ncombs.vec[1:(i-1)])+1:ncol(u.list[[i]])])
  }
  
  # for(i in 1:length(u.list))
  #   ints[[i]]<-apply(u.list[[i]],2,function(x) func.hat(x,pc.mod,pcs,mcmc.use,f0r2)) # the heavy lifting
  
  
  sob[[1]]<-ints[[1]]
  # matplot(t(apply(sob[[1]],1,cumsum)),type='l')
  # matplot(t(apply(sens.func$S.var[1,1:5,],2,cumsum)),type='l',add=T)
  
  V.tot<-func.hat(1:p,pc.mod,pcs,mcmc.use,f0r2)
  # plot(V.tot)
  # points(apply(sens.func$S.var[1,,],2,sum),col=2)
  
  cat('Shuffling',timestamp(quiet = T),'\n')
  if(length(u.list)>1){
    for(i in 2:length(u.list)){
      sob[[i]]<-matrix(nrow=nxfunc,ncol=ncol(ints[[i]]))
      for(j in 1:ncol(u.list[[i]])){
        cc<-rep(0,nxfunc)
        for(k in 1:(i-1)){
          ind<-which(apply(u.list[[k]],2,function(x) all(x%in%u.list[[i]][,j])))
          cc<-cc+(-1)^(i-k)*rowSums(ints[[k]][,ind])
        }
        sob[[i]][,j]<-ints[[i]][,j]+cc
      }
    }
  }
  
  
  # sens.func.use<-lapply(strsplit(sens.func$names.ind,'x'),as.numeric)
  # sl<-sapply(sens.func.use,length)
  # ind.list<-list()
  # sob.small<-list()
  # for(i in 1:length(u.list)){
  #   ind.list[[i]]<-NA
  #   k<-0
  #   for(j in which(sl==i)){
  #     k<-k+1
  #     ind.list[[i]][k]<-which(apply(u.list[[i]],2,function(x) all(x==sens.func.use[[j]])))
  #   }
  #   sob.small[[i]]<-sob[[i]][,ind.list[[i]]]
  # }
  # 
  # sob.small<-do.call(cbind,sob.small)
  # matplot(t(apply(sob.small,1,cumsum)),type='l')
  # matplot(t(apply(sens.func$S.var[1,,],2,cumsum)),type='l',add=T)
  
  #browser()
  
  if(is.null(nind))
    nind<-ncombs
  
  
  sob.comb.var<-do.call(cbind,sob)
  
  vv<-colMeans(sob.comb.var)
  ord<-order(vv,decreasing = T)
  cutoff<-vv[ord[nind]]
  if(nind>length(ord))
    cutoff<-min(vv)
  use<-sort(which(vv>=cutoff))
  
  
  V.other<-V.tot-rowSums(sob.comb.var[,use])
  use<-c(use,ncombs+1)
  
  sob.comb.var<-t(cbind(sob.comb.var,V.other))
  sob.comb<-t(t(sob.comb.var)/c(V.tot))
  
  sob.comb.var<-sob.comb.var[use,,drop=F]
  sob.comb<-sob.comb[use,,drop=F]
  
  dim(sob.comb)<-c(1,length(use),nxfunc)
  dim(sob.comb.var)<-c(1,length(use),nxfunc)
  
  
  names.ind<-c(unlist(lapply(u.list,function(x) apply(x,2,paste,collapse='x',sep=''))),'other')
  names.ind<-names.ind[use]
  
  cat('Finish',timestamp(quiet = T),'\n')
  
  #browser()
  
  ret<-list(S=sob.comb,S.var=sob.comb.var,Var.tot=V.tot,names.ind=names.ind,xx=seq(0,1,length.out = nxfunc),func=T)
  class(ret)<-'bassSob'
  plot(ret)
  return(ret)
  
}

################################################################################
## Functions
################################################################################
func.hat<-function(u,pc.mod,pcs,mcmc.use,f0r2,C1.array){ # could speed this up
  #browser()
  res<-rep(0,nrow(pcs))
  n.pc<-length(pc.mod)
  for(i in 1:n.pc){
    res<-res+pcs[,i]^2*Ccross(pc.mod,i,i,u,mcmc.use,C1.array)
    
    if(i<n.pc){
      for(j in (i+1):n.pc){
        res<-res+2*pcs[,i]*pcs[,j]*Ccross(pc.mod,i,j,u,mcmc.use,C1.array)
        #print(c(i,j))
      }
    }
  }
  return(res-f0r2)
}

Ccross<-function(pc.mod,i,j,u,mcmc.use=1,C1.array){ # inner product of main effects from different eof models
  p<-pc.mod[[1]]$p
  mcmc.mod.usei<-pc.mod[[i]]$model.lookup[mcmc.use]
  mcmc.mod.usej<-pc.mod[[j]]$model.lookup[mcmc.use]
  
  Mi<-pc.mod[[i]]$nbasis[mcmc.use]
  Mj<-pc.mod[[j]]$nbasis[mcmc.use]
  mat<-matrix(nrow=Mi,ncol=Mj)
  
  #CC<-C2.temp<-CCu<-matrix(1,nrow=Mi,ncol=Mj)
  
  a0i<-pc.mod[[i]]$beta[mcmc.use,1]
  a0j<-pc.mod[[j]]$beta[mcmc.use,1]
  f0i<-get.f0(pc.mod,i,mcmc.use)
  f0j<-get.f0(pc.mod,j,mcmc.use)
  
  out<- a0i*a0j + a0i*(f0j-a0j) + a0j*(f0i-a0i)
  #browser()
  
  if(Mi>0 & Mj>0){
    ai<-pc.mod[[i]]$beta[mcmc.use,1+1:Mi]
    aj<-pc.mod[[j]]$beta[mcmc.use,1+1:Mj]
    
    for(mi in 1:Mi){
      for(mj in 1:Mj){
        temp1<-ai[mi]*aj[mj]
        temp2<-temp3<-1
        for(l in (1:p)[-u]){
          #temp2<-temp2*C1(pc.mod,l,mi,i,mcmc.mod.usei)*C1(pc.mod,l,mj,j,mcmc.mod.usej) # make a C1 lookup table instead (this is the bottleneck)
          temp2<-temp2*C1.array[i,l,mi]*C1.array[j,l,mj]
          #browser()
        }
        #CC[mi,mj]<-temp2
        for(l in u){
          temp3<-temp3*C2(pc.mod,l,mi,mj,i,j,mcmc.mod.usei,mcmc.mod.usej) # would be nice to use a lookup table here too, but its too big
        }
        #C2.temp[mi,mj]<-temp3
        #CCu[mi,mj]<-temp4
        out<-out+temp1*temp2*temp3#(temp3-1) not -1 since we subtract f0^2 later
        #print(out)
        #mat[mi,mj]<-temp
        #print(c(temp1*temp2*temp3))
      }
    }
  }
  #out<-out+ai%*%(CC*C2.temp/CCu)%*%aj
  if(length(out)==0)
    browser()
  return(out)
}

C1<-function(pc.mod,l,m,pc,mcmc.mod.use){ # l = variable, m = basis function, pc = eof index
  if(l<=pc.mod[[pc]]$pdes){
    int.use.l<-which(pc.mod[[pc]]$vars.des[mcmc.mod.use,m,]==l)
    if(length(int.use.l)==0)
      return(1)
    s<-pc.mod[[pc]]$signs[mcmc.mod.use,m,int.use.l]
    t.ind<-pc.mod[[pc]]$knotInd.des[mcmc.mod.use,m,int.use.l]
    t<-pc.mod[[pc]]$xx.des[t.ind,l]
    q<-pc.mod[[pc]]$degree
    return((1/(q+1)*((s+1)/2-s*t))*s^2)
  } else{
    l.cat<-l-pc.mod[[pc]]$pdes # assumes that des vars come before cat vars, which I think we do internally.
    int.use.l<-which(pc.mod[[pc]]$vars.cat[mcmc.mod.use,m,]==l.cat)
    if(length(int.use.l)==0)
      return(1)
    lD1<-pc.mod[[pc]]$sub.size[mcmc.mod.use,m,int.use.l]
    nlevels<-pc.mod[[pc]]$nlevels[l.cat]
    return(lD1/nlevels)
  }
}
C2<-function(pc.mod,l,m1,m2,pc1,pc2,mcmc.mod.use1,mcmc.mod.use2){
  if(l<=pc.mod[[pc1]]$pdes){ # could do pc1 or pc2, they have the same vars
    int.use.l1<-which(pc.mod[[pc1]]$vars.des[mcmc.mod.use1,m1,]==l)
    int.use.l2<-which(pc.mod[[pc2]]$vars.des[mcmc.mod.use2,m2,]==l)
    if(length(int.use.l1)==0 & length(int.use.l2)==0)
      return(1)
    if(length(int.use.l1)==0)
      return(C1(pc.mod,l,m2,pc2,mcmc.mod.use2))
    if(length(int.use.l2)==0)
      return(C1(pc.mod,l,m1,pc1,mcmc.mod.use1))
    
    q<-pc.mod[[pc1]]$degree
    s1<-pc.mod[[pc1]]$signs[mcmc.mod.use1,m1,int.use.l1]
    s2<-pc.mod[[pc2]]$signs[mcmc.mod.use2,m2,int.use.l2]
    t.ind1<-pc.mod[[pc1]]$knotInd.des[mcmc.mod.use1,m1,int.use.l1]
    t.ind2<-pc.mod[[pc2]]$knotInd.des[mcmc.mod.use2,m2,int.use.l2]
    t1<-pc.mod[[pc1]]$xx.des[t.ind1,l]
    t2<-pc.mod[[pc2]]$xx.des[t.ind2,l]
    if(t2<t1){
      temp<-t1
      t1<-t2
      t2<-temp
      temp<-s1
      s1<-s2
      s2<-temp
    }
    return(C22(t1,t2,s1,s2,q,m1,m2,pc1,pc2))
  } else{
    l.cat<-l-pc.mod[[pc1]]$pdes
    
    int.use.l1<-which(pc.mod[[pc1]]$vars.cat[mcmc.mod.use1,m1,]==l.cat)
    int.use.l2<-which(pc.mod[[pc2]]$vars.cat[mcmc.mod.use2,m2,]==l.cat)
    
    if(length(int.use.l1)==0 & length(int.use.l2)==0)
      return(1)
    if(length(int.use.l1)==0)
      return(C1(pc.mod,l,m2,pc2,mcmc.mod.use2))
    if(length(int.use.l2)==0)
      return(C1(pc.mod,l,m1,pc1,mcmc.mod.use1))
    
    #browser()
    sub1<-pc.mod[[pc1]]$sub.list[[mcmc.mod.use1]][[m1]][[int.use.l1]]
    sub2<-pc.mod[[pc2]]$sub.list[[mcmc.mod.use2]][[m2]][[int.use.l2]]
    if(is.na(sub1[1]) & is.na(sub2[1]))
      browser()
    nlevels<-pc.mod[[pc1]]$nlevels[l.cat]
    return(length(intersect(sub1,sub2))/nlevels)
  }
}


C22<-function(t1,t2,s1,s2,q,m1,m2,pc1,pc2){ # t1<t2
  cc<-BASS:::const(signs=c(s1,s2),knots=c(t1,t2),degree=q)
  if((s1*s2)==0){
    return(0)
  }
  if(m1==m2 & pc1==pc2){ #t1=t2, s1=s2 - NOT TRUE, since these could be different eof models
    return(1/(2*q+1)*((s1+1)/2-s1*t1)^(2*q+1)/cc)
  } else{
    if(s1==1){
      if(s2==1){
        return(BASS:::intabq(t2,1,t1,t2,q)/cc)
      } else{
        return(BASS:::intabq(t1,t2,t1,t2,q)*(-1)^q/cc)
      }
    } else{
      if(s2==1){
        return(0)
      } else{
        return(BASS:::intabq(0,t1,t1,t2,q)/cc)
      }
    }
  }
}

get.f0<-function(pc.mod,pc,mcmc.use){ # mcmc.mod.use is mcmc index not model index
  mcmc.mod.use<-pc.mod[[pc]]$model.lookup[mcmc.use]
  out<-pc.mod[[pc]]$beta[mcmc.use,1] # intercept
  if(pc.mod[[pc]]$nbasis[mcmc.use] > 0){
    for(m in 1:pc.mod[[pc]]$nbasis[mcmc.use]){
      out1<-pc.mod[[pc]]$beta[mcmc.use,1+m]
      for(l in 1:pc.mod[[pc]]$p){
        out1<-out1*C1(pc.mod,l,m,pc,mcmc.mod.use)
      }
      out<-out+out1
    }
  }
  return(out)
}




# Ccross<-function(pc.mod,i,j,u,mcmc.use=1){ # inner product of main effects from different eof models
#   
#     Mi<-pc.mod[[i]]$nbasis[mcmc.use]
#     Mj<-pc.mod[[j]]$nbasis[mcmc.use]
#     if(Mi<Mj){ # so that the nested loop below works, and symmetric
#       temp<-i
#       i<-j
#       j<-temp
#     }
#   
#   mcmc.mod.usei<-pc.mod[[i]]$model.lookup[mcmc.use]
#   mcmc.mod.usej<-pc.mod[[j]]$model.lookup[mcmc.use]
#   
#   Mi<-pc.mod[[i]]$nbasis[mcmc.use]
#   Mj<-pc.mod[[j]]$nbasis[mcmc.use]
#   mat<-matrix(nrow=Mi,ncol=Mj)
#   
#   #CC<-C2.temp<-CCu<-matrix(1,nrow=Mi,ncol=Mj)
#   
#   a0i<-pc.mod[[i]]$beta[mcmc.use,1]
#   a0j<-pc.mod[[j]]$beta[mcmc.use,1]
#   f0i<-get.f0(pc.mod,i,mcmc.use)
#   f0j<-get.f0(pc.mod,j,mcmc.use)
#   
#   out<- a0i*a0j + a0i*(f0j-a0j) + a0j*(f0i-a0i)
#   #browser()
#   ai<-pc.mod[[i]]$beta[mcmc.use,1+1:Mi]
#   aj<-pc.mod[[j]]$beta[mcmc.use,1+1:Mj]
#   
#     #pc.modi<-pc.mod[[i]]
#     #pc.modj<-pc.mod[[j]]
#   
#   for(mi in 1:Mi){
#         temp1<-ai[mi]^2
#         temp2<-temp3<-1
#         for(l in (1:p)[-u]){
#           temp2<-temp2*C1(pc.mod,l,mi,i,mcmc.mod.usei)^2
#         }
#         for(l in u){
#           temp3<-temp3*C2(pc.mod,l,mi,mi,i,i,mcmc.mod.usei,mcmc.mod.usei)
#         }
#         out<-out+temp1*temp2*temp3
# 
#         if(mi<Mj){
#           for(mj in (mi+1):Mj){
#             temp1<-ai[mi]*aj[mj]
#             temp2<-temp3<-1
#             for(l in (1:p)[-u]){
#               temp2<-temp2*C1(pc.mod,l,mi,i,mcmc.mod.usei)*C1(pc.mod,l,mj,j,mcmc.mod.usej)
#             }
#             #CC[mi,mj]<-temp2
#             for(l in u){
#               temp3<-temp3*C2(pc.mod,l,mi,mj,i,j,mcmc.mod.usei,mcmc.mod.usej)
#             }
#             #C2.temp[mi,mj]<-temp3
#             #CCu[mi,mj]<-temp4
#             out<-out+2*temp1*temp2*temp3#(temp3-1) not -1 since we subtract f0^2 later
#             #print(out)
#             #mat[mi,mj]<-temp
#           }
#         }
#       }
#   #out<-out+ai%*%(CC*C2.temp/CCu)%*%aj
#   return(out)
# }
# ################################################################################
# ## Functions
# ################################################################################
# Ccross<-function(pc.mod,i,j,u,mcmc.use=1,w0){ # inner product of main effects from different eof models
#   
#   Mi<-pc.mod[[i]]$nbasis[mcmc.use]
#   Mj<-pc.mod[[j]]$nbasis[mcmc.use]
#   if(Mi<Mj){ # so that the nested loop below works, and symmetric
#     temp<-i
#     i<-j
#     j<-temp
#   }
#   Mi<-pc.mod[[i]]$nbasis[mcmc.use]
#   Mj<-pc.mod[[j]]$nbasis[mcmc.use]
#   
#   mcmc.mod.usei<-pc.mod[[i]]$model.lookup[mcmc.use]
#   mcmc.mod.usej<-pc.mod[[j]]$model.lookup[mcmc.use]
#   
#   
#   mat<-matrix(nrow=Mi,ncol=Mj)
#   
#   #CC<-C2.temp<-CCu<-matrix(1,nrow=Mi,ncol=Mj)
#   
#   a0i<-pc.mod[[i]]$beta[mcmc.use,1]
#   a0j<-pc.mod[[j]]$beta[mcmc.use,1]
#   f0i<-w0[i]#get.f0(pc.mod,i,mcmc.use)
#   f0j<-w0[j]#get.f0(pc.mod,j,mcmc.use)
#   
#   out<- a0i*a0j + a0i*(f0j-a0j) + a0j*(f0i-a0i)
#   #browser()
#   ai<-pc.mod[[i]]$beta[mcmc.use,1+1:Mi]
#   aj<-pc.mod[[j]]$beta[mcmc.use,1+1:Mj]
#   
#   pc.modi<-pc.mod[[i]]
#   pc.modj<-pc.mod[[j]]
#   
#   for(mi in 1:Mi){
#     temp1<-ai[mi]^2
#     temp2<-temp3<-1
#     for(l in (1:p)[-u]){
#       temp2<-temp2*C1(pc.modi,l,mi,mcmc.mod.usei)^2
#     }
#     for(l in u){
#       temp3<-temp3*C2(pc.mod,l,mi,mi,i,i,mcmc.mod.usei,mcmc.mod.usei)
#     }
#     out<-out+temp1*temp2*temp3
#     
#     if(mi<Mj){
#       for(mj in (mi+1):Mj){
#         temp1<-ai[mi]*aj[mj]
#         temp2<-temp3<-1
#         for(l in (1:p)[-u]){
#           temp2<-temp2*C1(pc.modi,l,mi,mcmc.mod.usei)*C1(pc.modj,l,mj,mcmc.mod.usej)
#         }
#         #CC[mi,mj]<-temp2
#         for(l in u){
#           temp3<-temp3*C2(pc.mod,l,mi,mj,i,j,mcmc.mod.usei,mcmc.mod.usej)
#         }
#         #C2.temp[mi,mj]<-temp3
#         #CCu[mi,mj]<-temp4
#         out<-out+2*temp1*temp2*temp3#(temp3-1) not -1 since we subtract f0^2 later
#         #print(out)
#         #mat[mi,mj]<-temp
#       }
#     }
#   }
#   #out<-out+ai%*%(CC*C2.temp/CCu)%*%aj
#   return(out)
# }
# 
