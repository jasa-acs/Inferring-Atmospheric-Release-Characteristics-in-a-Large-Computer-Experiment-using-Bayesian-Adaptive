
library(parallel)
library(BASS)
library(mnormt)




labels <- sapply(1:4,function(i)
  as.expression(bquote(10^ .(i)))
)

sites.plot<-c(325,333,340,343,414,412,423,410,426,225,108,526,427,510,501)
ord<-order(sites.plot)
sites.plot<-sites.plot[ord]
site.names<-c('1 mile from plant','4 miles from plant','Whalers Island','Avila Beach','mouth of See Canyon','midway up See Canyon','1 mile inland along Hwy 101','coastline extent of plume','inland from coastline extent','San Luis Obisbo','Morro Bay','Santa Maria','Pismo Beach','Arroyo Grande','Edna')[ord]



n.eofs<-n.pc
eofs<-eofs.sim[,1:n.eofs]
TrK<-eofs
Trmm<-0


# max number of basis functions plus intercept
dd<-c(mb,n.pc) # dimensions used in some predict functions



# load('flex_dat.rda')
# trunc.error<- Y.sim[,-holdout.sim][,-ho]-K%*%reduced.y.sim[1:n.pc,-ho] # deviations
# rm(Y.sim)

times.disp<-paste(rep(c(paste(0,6:9,sep=''),10:22),each=2),rep(c('00','30'),17),sep=':')

########################################################################
## Dispcrepancy plot
space.mat<-kronecker(1:137,t(rep(1,34)))
space.vec<-c(space.mat)
B2<-matrix(ncol=137,nrow=137*34) # transformation from spatial locations to 
#  repeated spatial locations to match pattern in y
for(j in 1:137){
  B2[,j]<-as.numeric(space.vec==j) 
  # B[i,j] = 1 if y[i] is at location j.  Hence, colSums(B)==nt and rowSums(B)
  #   is a vector of ones.
}
newdat<-cbind(B2%*%as.matrix(coordsLL),c(kronecker(rep(1,137),t(times.modelRuns))))
pdiscrep<-predict(dmod.prior,newdat)
dmm<-colMeans(pdiscrep)
dmm.mat<-t(matrix(dmm,ncol=34))

library(fdakma)
cl.test<-kma(x=t(kronecker(1:34,t(rep(1,137))))[,-1],y0=t(dmm.mat)[,-1],y1=t(apply(dmm.mat,2,diff)),warping.method='NOalignment',similarity.method='d1.L2',n.clust = 5,center.method='k-means')
cols=cl.test$labels

## make calibrated output plots for each cluster?
sites.plot<-sites.modelRuns[cols==2]
ord<-order(sites.plot)
sites.plot<-sites.plot[ord]
site.names<-sites.plot


matplot(times.modelRuns,dmm.mat,type='l',col=cols) # cluster these to try to understand the discrepancy better
abline(v=c(7,19))

pdf(paste(out.name,'_discrep.pdf',sep=''),width=10,height=7)
par(mfrow=c(2,3))
for(i in 1:5){
  matplot(times.modelRuns,dmm.mat[,cols==i],type='l',col=i,ylim=range(dmm.mat),ylab='discrepancy',xlab='time',main=paste('Cluster',i))
  abline(v=c(7,19),col='lightgrey',lty=2)
}
plot(coordsLL,col=cols,pch=cols)
legend('bottomleft',paste('cluster',1:5),col=1:5,pch=1:5)
dev.off()

matplot(times.modelRuns[-1],apply(dmm.mat,2,diff),type='l',col=cols)
########################################################################
## Calibrated prediction plot


locs<-which(sites.modelRuns%in%sites.plot)
#times.disp<-paste(rep(c('05','06','07','08','09',10:21),each=2),':',rep(c('00','30'),17),sep='')

#YY.use<-YY.use+1#140 # choose one of the holdouts to calibrate to
mat<-matrix(1:4658,nrow=137)
mcmc.use<-(1:nmcmc)[-burn][seq(1,20000,20)]
pred<-array(dim=c(137,34,length(mcmc.use)))
aa<-aa2<-aa3<-matrix(nrow=137*34,ncol=length(mcmc.use))





  k3=0
  for(i in mcmc.use){
    k3=k3+1
    temp<-emuSample(as.numeric(theta1[i,]),theta2.ind[i])
    aa[,k3]<-temp$eta
    aa2[,k3]<-temp$eta.mean
    itd<-sample(1000,size=1)
    aa3[,k3]<-gam[i]*pdiscrep[itd,] # don't include the discrepancy error?
    #aa3[,k3]<-rnorm(137*34,pdiscrep[it,],sqrt(dmod.prior$s2[it]))
  }
  
  #pred[,,k2]<-array(y-rowMeans(aa),dim=c(137,34))

  
  #pdf('~/Desktop/flex_calib_fitGam-clust3.pdf',width=10,height=6)
  pdf(paste(out.name,'_discrepPred.pdf',sep=''),width=10,height=6)
  par(mfrow=c(3,4),mar=c(0,0,0,0),oma=c(5.5,5.5,2,2))
  qq<-apply(aa,1,quantile,probs=c(.025,.975))
  mm<-apply(aa,1,mean)
  
  k<-0
  for(ll in locs){
    k<-k+1
    
    plot(1,type='n',ylim=c(1,5),xaxt='n',yaxt='n',xlim=range(times.modelRuns),ylab='')
    polygon(c(times.modelRuns,rev(times.modelRuns)),c(qq[1,mat[ll,]],rev(qq[2,mat[ll,]])),col='lightgrey',border = F)
    lines(times.obs,yobs[ll,],col=4,lwd=2,lty=3,type='b')
    lines(times.modelRuns,mm[mat[ll,]],col='darkgrey',lwd=1,lty=5)
    
    
    
    text(21,4.8,site.names[k])
    #text(15,4.5,site.names[k])
    if(k %in% c(1,5,9))
      axis(2,at=1:5,labels = c(labels,NA),las=1)
    if(k %in% 9:15)
      axis(1,at = times.modelRuns[seq(1,34,4)],labels = times.disp[seq(1,34,4)],las=2)
    #if(k==1)
      #legend(x=6,y=3,legend=c('observation','prediction'),col=c(4,'darkgrey'),lty=c(3,5),lwd=c(2,1),bty='n')
  }
  mtext('time',1,outer = T,line=4)
  mtext('concentration',2,outer = T,line=4)

plot.new()
legend('topright',c('observations','mean predictions'),col=c(4,'lightgrey'),bty='n',lty=c(3,5),pch=c(1,-1))
  qq<-apply(aa+aa3,1,quantile,probs=c(.025,.975))
  mm<-apply(aa+dmm,1,mean)

  
  k<-0
  for(ll in locs){
    k<-k+1
    
    plot(1,type='n',ylim=c(1,5),xaxt='n',yaxt='n',xlim=range(times.modelRuns),ylab='')
    polygon(c(times.modelRuns,rev(times.modelRuns)),c(qq[1,mat[ll,]],rev(qq[2,mat[ll,]])),col='lightgrey',border = F)
    lines(times.obs,yobs[ll,],col=4,lwd=2,lty=3,type='b')
    lines(times.modelRuns,mm[mat[ll,]],col='darkgrey',lwd=1,lty=5)
    
    
    
    text(21,4.8,site.names[k])
    #text(15,4.5,site.names[k])
    if(k %in% c(1,5,9))
      axis(2,at=1:5,labels = c(labels,NA),las=1)
    if(k %in% 9:15)
      axis(1,at = times.modelRuns[seq(1,34,4)],labels = times.disp[seq(1,34,4)],las=2)
    #if(k==1)
      #legend(x=6,y=3,legend=c('observation','prediction'),col=c(4,'darkgrey'),lty=c(3,5),lwd=c(2,1),bty='n')
  }
  mtext('time',1,outer = T,line=4)
  mtext('concentration',2,outer = T,line=4)
  
  
  
dev.off()






par(mfrow=c(1,2))

mm1<-rowMeans(Tr%*%aa2)[use]
mm2<-rowMeans(Tr%*%(aa2+aa3))[use]
rran<-range(c(y,mm1,mm2))
pdf(paste(out.name,'_obsPred.pdf',sep=''),width=9,height=4)
par(mfrow=c(1,2))
plot(y,mm1,ylim=rran,xlim=rran,xlab='observed',ylab='predicted',main='Without Discrepancy',cex=.5);abline(a=0,b=1,col=2)
plot(y,mm2,ylim=rran,xlim=rran,xlab='observed',ylab='predicted',main='With Discrepancy',cex=.5);abline(a=0,b=1,col=2)
dev.off()
