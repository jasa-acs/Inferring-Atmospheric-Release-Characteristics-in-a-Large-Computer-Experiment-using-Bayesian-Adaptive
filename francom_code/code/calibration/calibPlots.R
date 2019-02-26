library(MASS)
library(dismo)

#lon<-theta1[-burn,1]
#lat<-theta1[-burn,2]
# true.lat<-(35.2098-35.1977)/(35.2250-35.1977)
# true.lon<-(-120.8542--120.8708)/(-120.8384--120.8708)
# true.amt<-0.5822002#136.016/(990)#0.5822002#
# true.dur<-true.st<-.5
# true.z<-1/9
# true.x<-c(true.lat,true.lon,true.z,true.st,true.dur,true.amt)
ss.bw<-.05
kde.bw<-.15
pdf(paste(out.name,'.pdf',sep=''),width=10,height=10)
par(mfrow=c(6,6),oma=c(6,4,4,4),mar=c(0,0,0,0))
lab<-c('latitude','longitude','altitude (m)','start time','duration (h)','amount (log10 kg)')
k=0
rr.true.x<-matrix(c(35.1977,35.2250,-120.8708,-120.8384,1,10,7,9,6,10,1,3),ncol=2,byrow=T)
lo<-c(4,4,4,5,5,5)
ro<-c(3,3,1,1,1,1)
par()$mfrow
for(i in 1:6){
  for(j in 1:6){
    if(i==j){
      col<-blues9[4]
      print(par()$mfrow)
      hist(theta1[-burn,i],xlim=c(0,1),main='',freq=F,xlab='',ylab='',yaxt='n',xaxt='n',col=col)
      #lines(density(theta1[-burn,i]),xlim=c(0,1),main='',xlab='',ylab='',yaxt='n',col=blues9[9])
      mtext(lab[i],1,line=4,cex=1.5)
      axis(1,at = seq(0,1,length.out=lo[i]),labels = round(BASS:::unscale.range(seq(0,1,length.out=lo[i]),rr.true.x[i,]),ro[i]),cex.axis=1.5)
      box()
      abline(v=true.x[i],col=2)
    } else if(i<j){
      smoothScatter(theta1[-burn,c(j,i)],xlim=c(0,1),ylim=c(0,1),bandwidth=ss.bw,xlab='',ylab='',xaxt='n',yaxt='n',nbin=300)
      
      dens <- kde2d(theta1[-burn,j], theta1[-burn,i], n=100,h=kde.bw)
      prob <- c(.95)
      dx <- diff(dens$x[1:2])
      dy <- diff(dens$y[1:2])
      sz <- sort(dens$z)
      c1 <- cumsum(sz) * dx * dy
      levels <- sapply(prob, function(x) {
        approx(c1, sz, xout = 1 - x)$y
      })
      contour(dens, levels=levels, labels=prob, add=T,drawlabels = F)

      points(true.x[j],true.x[i],col=2,pch='x',cex=1.5)
    } else{
      plot.new()
    }
  }
}

par()$mfrow
lab2<-c('init time','PBL','LSM','nudge','reanalysis')#names(X[1,7:11])

par(mfrow=c(6,5),oma=c(4,4,4,4),mar=c(0,0,0,0))
for(i in 1:6){
  for(j in 1:5){
    plot(density(theta1[-burn,i]),main='',xlab='',ylab='',xaxt='n',yaxt='n',lwd=3)
    if(j==1)
      mtext(lab[i],2,line=2.5,cex=.5)
    if(i==6)
      mtext(lab2[j],1,line=2.5,cex=.5)
    p<-1
    for(kk in cat.names[[j]]){
      p<-p+1
      th<-theta1[-burn,i][cats[theta2.ind[-burn],j]==kk]
      if(length(th)>1){
        dd<-density(th)
        lines(dd$x,dd$y*sum(cats[theta2.ind[-burn],j]==kk)/(nmcmc-length(burn)),col=p,lwd=3)
      } else{
        abline(h=0,col=p,lwd=3)
      }
    }
  }
}

cat.names.disp2<-list(c('a9','b15'),c('aYSU','bMYJ','cMYNN'),c('aTD','bNoah','cRUC'),c('anone','blow','chigh'),c('aNARR','bECMWF','cCFSR'))
cat.names.disp<-list(c('9','15'),c('YSU','MYJ','MYNN'),c('TD','Noah','RUC'),c('none','low','high'),c('NARR','ECMWF','CFSR'))

par(mfrow=c(5,5),oma=c(4,4,4,4),mar=c(0,0,0,0))
for(i in 1:5){
  thi<-cat.names.disp2[[i]][as.numeric(as.character(cats[theta2.ind[-burn],i]))+1]#as.numeric(as.character(cats[theta2.ind[-burn],i]))
  for(j in 1:5){
    ni<-length(cat.names[[i]])
    if(i==j){
      tab<-table(c(cat.names.disp2[[i]],thi))
      names(tab)<-cat.names.disp[[i]]
      barplot((tab-1)/(nmcmc-length(burn)),col=1:ni+1,main='',yaxt='n')
      mtext(lab2[i],1,line=2.5,cex=.5)
    } else if(i<j){
      thj<-cat.names.disp2[[j]][as.numeric(as.character(cats[theta2.ind[-burn],j]))+1]#a
      nj<-length(cat.names[[j]])
      gg<-expand.grid(1:nj,1:ni)
      plot(gg,xlab='',ylab='',xaxt='n',yaxt='n',type='n',xlim=range(gg[,1])+c(-.5,.5),ylim=range(gg[,2])+c(-.5,.5))
      for(ki in 1:ni){
        for(kj in 1:nj){
          text(kj,ki,round(mean(thi==cat.names.disp2[[i]][[ki]] & thj==cat.names.disp2[[j]][[kj]]),2))
        }
      }
    } else{
      plot.new()
    }
  }
}

par(mfrow=c(1,5),oma=c(4,5,4,4),mar=c(0,0,0,0))
for(i in 1:5){
  thi<-cat.names.disp2[[i]][as.numeric(as.character(cats[theta2.ind[-burn],i]))+1]
  ni<-length(cat.names[[i]])
  tab<-table(c(cat.names.disp2[[i]],thi))
  names(tab)<-cat.names.disp[[i]]
  plot(1,type='n',ylim=c(0,1),xlim=c(-.2,.5),xaxt='n',yaxt='n',bty='n')
  abline(h=seq(0,1,.2),col='lightgrey')
  ttt<-barplot((tab-1)/(nmcmc-length(burn)),col='cornflowerblue',main='',yaxt='n',ylim=c(0,1),width=.1,xlim=c(-.2,.5),add=T,las=2,border=NA,names.arg='')
  barplot((tab-1)/(nmcmc-length(burn)),col='cornflowerblue',main='',yaxt='n',ylim=c(0,1),width=.1,xlim=c(-.2,.5),add=T,las=2,border=NA,names.arg='')
  mtext(lab2[i],3,line=0,cex=1)
  mtext(names(tab),1,at = c(ttt),las=2,line=0,cex=.7)
  if(i==1){
    axis(2)
    mtext('probability',side = 2,line=3.5)
  } else
    axis(2,at = c(-1,2),labels=c('',''))
}



par(mfrow=c(1,1))

if(F){
  
library(RgoogleMaps)
lat.std<-c(35.1977, 35.2098, 35.2250)
lon.std<-c(-120.8708, -120.8542, -120.8384)

mcmc.use<-(1:nmcmc)[-burn]

tlat<-true.x[1]*(lat.std[3]-lat.std[1])+lat.std[1]
tlon<-true.x[2]*(lon.std[3]-lon.std[1])+lon.std[1]
lat<-theta1[mcmc.use,1]*(lat.std[3]-lat.std[1])+lat.std[1]
lon<-theta1[mcmc.use,2]*(lon.std[3]-lon.std[1])+lon.std[1]

MyMap <- GetMap(center=c(lat.std[2],lon.std[2]), zoom=min(MaxZoom(range(lat.std), range(lon.std))), destfile = "MyTile1.png",maptype='satellite',NEWMAP = F)

PlotOnStaticMap(MyMap, lat = lat,lon = lon,cex=.3)
PlotOnStaticMap(MyMap, lat = c(unlist(tlat)),lon = c(unlist(tlon)),cex=1.5,col='red',pch='x',add=T)

#PlotOnStaticMap(MyMap, lat = coordsLL[,2],lon = coordsLL[,1],cex=1,add=T,col=2)

PlotArrowsOnStaticMap(MyMap,lat.std[1],lon.std[1],lat1=lat.std[3],add=T,FUN=segments,col='blue')
PlotArrowsOnStaticMap(MyMap,lat.std[1],lon.std[1],lon1=lon.std[3],add=T,FUN=segments,col='blue')
PlotArrowsOnStaticMap(MyMap,lat.std[3],lon.std[3],lat1=lat.std[1],add=T,FUN=segments,col='blue')
PlotArrowsOnStaticMap(MyMap,lat.std[3],lon.std[3],lon1=lon.std[1],add=T,FUN=segments,col='blue')


cLL<-data.frame(x=lon.std,y=lat.std)
g <- gmap(cLL, type='terrain', interpolate=T, lonlat=T,style = 'element:labels|visibility:off')#https://developers.google.com/maps/documentation/static-maps/styling
#g <- gmap(cLL, type='satellite', interpolate=T, lonlat=T,zoom=16)
dens <- kde2d(lon, lat, n=100,h=kde.bw/50)
prob <- seq(.05,.95,.1)#c(.1,.3,.5,.7,.95)
dx <- diff(dens$x[1:2])
dy <- diff(dens$y[1:2])
sz <- sort(dens$z)
c1 <- cumsum(sz) * dx * dy
levels <- sapply(prob, function(x) {
  approx(c1, sz, xout = 1 - x)$y
})




dev.off()
pdf(paste(out.name,'3.pdf',sep=''),width=5,height=5)

plot(g,interpolate=TRUE)
contour(dens,add=T,drawlabels = T,levels=seq(0,max(dens$z),4000))#pretty(range(dens$z),5))
points(tlon,tlat,col='darkgreen',pch=3,cex=1,lwd=2)
points(-120.8542,35.2098,col='darkred',pch=8,cex=1,lwd=1)
#points(-120.854355,35.211125,col=3,pch='x',cex=1.5)
points(coordsLL[,1],coordsLL[,2],col=2,cex=.5)

segments(lon.std[1],lat.std[1],lon.std[3],lat.std[1],lty=3)
segments(lon.std[1],lat.std[1],lon.std[1],lat.std[3],lty=3)
segments(lon.std[1],lat.std[3],lon.std[3],lat.std[3],lty=3)
segments(lon.std[3],lat.std[3],lon.std[3],lat.std[1],lty=3)

legend(-120.845,y=35.234,c('reported release location','corrected release location','measurement sites','prior range','posterior contours'),col=c('darkred','darkgreen',2,1,1),pt.cex=c(1,1,.5),pch=c(8,3,1,NA,NA),lty=c(-1,-1,-1,3,1),xpd=T,cex=.5,lwd=c(1,2,1,1,1),bg='white')
}

dev.off()

plot(2)
dev.off()

pdf(paste(out.name,'2.pdf',sep=''),width=8,height=2.5)

par(mfrow=c(1,5),oma=c(4,5,4,4),mar=c(0,0,0,0))
for(i in 1:5){
  thi<-cat.names.disp2[[i]][as.numeric(as.character(cats[theta2.ind[-burn],i]))+1]
  ni<-length(cat.names[[i]])
  tab<-table(c(cat.names.disp2[[i]],thi))
  names(tab)<-cat.names.disp[[i]]
  plot(1,type='n',ylim=c(0,1),xlim=c(-.2,.5),xaxt='n',yaxt='n',bty='n')
  abline(h=seq(0,1,.2),col='lightgrey')
  ttt<-barplot((tab-1)/(nmcmc-length(burn)),col='cornflowerblue',main='',yaxt='n',ylim=c(0,1),width=.1,xlim=c(-.2,.5),add=T,las=2,border=NA,names.arg='')
  barplot((tab-1)/(nmcmc-length(burn)),col='cornflowerblue',main='',yaxt='n',ylim=c(0,1),width=.1,xlim=c(-.2,.5),add=T,las=2,border=NA,names.arg='')
  mtext(lab2[i],3,line=0,cex=1)
  mtext(names(tab),1,at = c(ttt),las=2,line=0,cex=.7)
  if(i==1){
    axis(2)
    mtext('probability',side = 2,line=3.5)
  } else
    axis(2,at = c(-1,2),labels=c('',''))
}

dev.off()
