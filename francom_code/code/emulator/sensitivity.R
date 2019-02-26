load('flex_bassEmu.rda')
source('code/emulator/sens_eof.R')



# library(BASS)
n.pc<-100
pcs<-eofs.sim[,1:n.pc]
pc.mod<-mod[1:n.pc]
p<-pc.mod[[1]]$p

######################################################################
## do the heavy lifting
######################################################################
pc.sob.func<-sob.eof(pc.mod,int.order=3,pcs,mcmc.use=1,ncores = ncores, preschedule=F)
save(pc.sob.func,file = 'pcSobFunc.rda')
######################################################################  


dim(pc.sob.func$S)
pc.sob.func$names.ind
# matplot(t(pc.sob.func$S[1,,seq(1,4658,137)+100]),type='l') # indices for 100th location


## go through different locations, look at indices over time
ii<-0
(ii=ii+1)
pc.sob.func1<-pc.sob.func
pc.sob.func1$S<-pc.sob.func1$S[,,seq(1,4658,137)+ii,drop=F]
pc.sob.func1$S.var<-pc.sob.func1$S.var[,,seq(1,4658,137)+ii,drop=F]
pc.sob.func1$Var.tot<-pc.sob.func1$Var.tot[seq(1,4658,137)+ii,,drop=F]
pc.sob.func1$xx<-1:34
plot(pc.sob.func1)

# some of these are exactly the same (eg 47-48), probably sites so close they are in the same grid box for the simulations
plot(coordsLL,type='n')
text(coordsLL[,1],coordsLL[,2],1:137)
text(coordsLL[,1],coordsLL[,2],sites.modelRuns,cex=.5)
points(coordsLL[86:89,],col=2)

points(coordsLL[locs,],col=2)

points(coordsLL[76,1],coordsLL[76,2],col=2)


sites.plot<-c(325,333,340,343,414,412,423,410,426,225,108,526,427,510,501)
ord<-order(sites.plot)
sites.plot<-sites.plot[ord]
which(sites.obs%in%sites.plot)
site.names<-c('1 mile from plant','4 miles from plant','Whalers Island','Avila Beach','mouth of See Canyon','midway up See Canyon','1 mile inland along Hwy 101','coastline extent of plume','inland from coastline extent','San Luis Obisbo','Morro Bay','Santa Maria','Pismo Beach','Arroyo Grande','Edna')[ord]




poly.between<-function(y1,y2,x=NULL,...){
  if(is.null(x))
    x<-seq_len(length(y1))
  polygon(c(x,rev(x)), c(y1,rev(y2)), ...)
}
order.mat<-function(mat){
  do.call(order,lapply(seq_len(ncol(mat)), function(i) mat[,i]))
}

devtools::install_github("karthik/wesanderson")

pdf('sensPlots.pdf',width=10,height=6)
par(mfrow=c(3,5),mar=c(0,0,0,0),oma=c(5.5,5.5,2,2))
locs<-which(sites.modelRuns%in%sites.plot)
times.disp<-paste(rep(c('05','06','07','08','09',10:21),each=2),':',rep(c('00','30'),17),sep='')

locs=locs-1 #Y.sim leaves out 50th location

library(wesanderson)

nuse<-5 
nn<-list()
k<-0
for(ll in locs){
  k<-k+1
  
  x.mean.var <- apply(pc.sob.func$S.var[,,seq(1,4658,137)+ll,drop=F], 2:3, mean)
  x.mean.var2<-x.mean.var[-nrow(x.mean.var),]
  rm<-rowMeans(x.mean.var2)
  use<-rm > .007#sort(order(rm,decreasing = T)[1:nuse])
  nn[[k]]<-pc.sob.func$names.ind[-nrow(x.mean.var)][use]
}

effects<-c(unique(unlist(nn)),'other')
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
cols<-c(wes_palette("Darjeeling1"),wes_palette("Moonrise2"),wes_palette("Darjeeling2"))[1:(length(effects)+1)]#sample(color)[1:(length(effects)+1)]#sample(rainbow(length(effects)))
names.X<-c('latitude','longitude','altitude','start time','duration','amount','init time','PBL','LSM','nudge','reanalysis')

effects<-effects[-length(effects)]
l1<-strsplit(effects,'x')
l1.len<-sapply(l1,length)
nint<-max(l1.len)
# effects<-c(effects[l1.len==1],effects[l1.len==2],effects[l1.len==3])
ord<-list()
ord[[1]]<-which(l1.len==1)[order(as.numeric(unlist(l1[l1.len==1])))]
for(i in 2:nint){
  l1.m<-do.call(rbind,l1[l1.len==i])
  mode(l1.m)<-'numeric'
  ord[[i]]<-which(l1.len==i)[order.mat(l1.m)]#+max(ord[[i-1]])
}
effects<-c(effects[unlist(ord)],'other (2nd/3rd)','other (>3rd)')


k<-0
for(ll in locs){
  k<-k+1
  
  x.mean.var <- apply(pc.sob.func$S.var[,,seq(1,4658,137)+ll,drop=F], 2:3, mean)
  x.mean.var2<-x.mean.var[-nrow(x.mean.var),]
  
  rm<-rowMeans(x.mean.var2)
  use<-pc.sob.func$names.ind[-nrow(x.mean.var)]%in%effects#rm > .005#sort(order(rm,decreasing = T)[1:nuse])
  nn<-c(pc.sob.func$names.ind[-nrow(x.mean.var)][use],'other (2nd/3rd)','other (>3rd)')
  
  #use<-sort(order(rm,decreasing = T)[1:nuse])
  #nn<-c(pc.sob.func$names.ind[use],'other')
  x.mean.var2use<-x.mean.var2[use,,drop=F]
  x.mean.var2<-rbind(x.mean.var2use,
                     colSums(x.mean.var2)-colSums(x.mean.var2use), # other 2/3 way interactions
                     pc.sob.func$Var.tot[seq(1,4658,137)+ll,]-colSums(x.mean.var2)) # >3 way interactions
  plot(1:34,type='n', xlab = "", ylab = "", main = "",ylim=c(0,.75),xaxt='n',yaxt='n')
  abline(v=seq(1,34,4),col='lightgrey')
  abline(h=seq(0,.7,.1),col='lightgrey')
  x.mean.var2<-rbind(0,x.mean.var2)
  cs<-t(apply(x.mean.var2, 2, cumsum))
  #matplot(cs, type = "l",add=T)
  #cs<-cbind(0,cs)
  for(i in 1:length(nn)){
    coli<-cols[nn[i]==effects]
    poly.between(cs[,i],cs[,i+1],col=coli,border=NA)
  }
  
  text(23,.7,site.names[k])
  if(k %in% c(1,6,11))
    axis(2)
  if(k %in% 11:15)
    axis(1,at = seq(1,34,4),labels = times.disp[seq(1,34,4)],las=2)
  if(k==1){
    legend('bottomleft',rev(nn),fill=rev(cols),cex=1.1,y.intersp=.7)
    legend('bottomright',rev(paste(1:length(names.X),names.X,sep=' - ')))
  }
  
}
mtext('time',1,outer = T,line=4)
mtext('partitioned variance',2,outer = T,line=4)
effects


dev.off()






