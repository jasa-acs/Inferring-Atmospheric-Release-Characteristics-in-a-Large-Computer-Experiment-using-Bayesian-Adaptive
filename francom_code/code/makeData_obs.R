## Get observational data for one release - remove background

# get model run characteristics (sites and times)
library(ncdf4)
con<-nc_open('data/flxout-dopptex-timeseries-lhs04.nc')
#print(con)
# 2 variables (excluding dimension variables):
#   int DOPPTEX_Sites[dopptex_site_id,record]   
# float CONC[dopptex_site_id,time,record]   
# 
# 3 dimensions:
#   record  Size:18000   *** is unlimited ***
#   dopptex_site_id  Size:137
# time  Size:34

tt<-ncvar_get(con,varid='DOPPTEX_Sites')
sites.modelRuns<-tt[,1]
times.modelRuns<-seq(6,22.5,.5) # times (PST) of model measurements (averages over the following .5 hour).  So the first measurement is the average from 6.0-6.5 PST and so on
nc_close(con)

# read in observations
dd<-list()
for(i in 0:7){ # eight experimental releases, #3 is the one we analyze
  read.table('data/DOPPTEX/FILE1.DAT',skip=i*1801,nrows=1) # description line
  dd[[i+1]]<-read.table('data/DOPPTEX/FILE1.DAT',skip=i*1801+1,nrows=1800) # actual data, 150 sites 12 times => 1800 obs
}
dd2<-dd
for(i in 1:8){ # recode missing
  nas<-which(dd2[[i]][,6]==999999) # these are missing, NOT ZEROS
  dd2[[i]][nas,6]<-NA
}


## get site coordinates, translate to lat/lon
# library(PBSmapping)
# coordsUTM<-dd[[1]][1:150,3:4]
# attr(coordsUTM,'zone')<-10
# attr(coordsUTM,'projection')<-'UTM'
# names(coordsUTM)<-c('X','Y')
# coordsLL<-convUL(coordsUTM,km = F)

## site coordinates according to NARAC - more accurate because of transformation
cdat<-read.csv('data/DOPTEX_sensorLocationsLatLonWGS84_fromUTMNAD27.csv')
coordsLL2<-cdat[,c(2,7:8)]
coordsLL<-coordsLL2[,c(3,2)]

## not all 150 sites are used in model runs (because some are really far away, didn't want to make simulator grid extend so far)
sites.obs<-dd2[[3]][1:150,1]
sites.obs.use<-which(sites.obs%in%sites.modelRuns) # index of the 137 locations used in the model runs
all(sites.obs[sites.obs.use]==sites.modelRuns)

times.obs<-7+0:11 # times for obs (averages over the following hour).  So the first measurement is the average from 7-8
# want to average model runs so that times match this grid


tTr<-cbind(0,0,kronecker(diag(12),t(matrix(c(1,1)))),0,0,0,0,0,0,0,0)/2 # the time transformation matrix, so if K is 137x34 model run output, K%*%t(tTr) transforms it to be on the same time scale as the observations

# if K is 137x34 matrix of data, then c(K), which stacks columns, is the right order
Tr<-kronecker(tTr,diag(137))

# background level
backg<-dd2[[3]][sites.obs.use,6]
elevation<-dd2[[3]][sites.obs.use,5]
ind<-which(is.na(backg))
# backg[ind]<-median(backg,na.rm = T) # impute missing with median
# impute missing spatially
library(fields)
coordsLL<-coordsLL[sites.obs.use,]
kk<-Krig(coordsLL[-ind,],log10(backg[-ind]),sigma2=0,rho=1,theta=.01,m=1) # length-scale of .01 is roughly 1 km
plot(kk)
surface(kk,zlim = c(0,5.5),type='I')
surface(kk,type='I')
kk.pred<-10^(predict(kk,cbind(coordsLL)[ind,]))
backg.krig<-backg
backg.krig[ind]<-kk.pred


pos<-function(x){
  x[x<0]<-0
  x
}
dat<-matrix(nrow=137,ncol=12)
count<-count.all<-count.vneg<-0
for(i in 1:12){
  new.val<-dd2[[3]][150*(i-1)+sites.obs.use,6]-backg.krig
  #print(min(new.val,na.rm = T))
  #plot(coordsLL)
  #points(coordsLL[which(new.val<0),],col=2)
  #quilt.plot(coordsLL,new.val)
  count<-count+length(which(new.val<0))
  count.all<-count.all+length(which(!is.na(new.val)))
  count.vneg<-count.vneg+length(which(new.val< -50))
  #hist(new.val)
  dat[,i]<-log10(50+pos(new.val))
}
#for(i in 1:12){ # power transform
#  dat[,i]<-(20+pos(dd2[[3]][150*(i-1)+sites.obs.use,6]-backg))^(.025)
#}

count/count.all
count.vneg/count.all

yobs<-dat
save(yobs,Tr,coordsLL,tTr,times.obs,sites.obs,sites.obs.use,sites.modelRuns,times.modelRuns,file='obs.rda')




if(F){
  
  
############################# ############################# #############################
## some plots
############################# ############################# #############################

dat2<-matrix(nrow=137,ncol=12)
for(i in 1:12){
  dat2[,i]<-log10(20+pos(dd2[[3]][150*(i-1)+sites.obs.use,6]))
}
matplot(t(dat2),type='l')



ylim<-range(1,dat2,na.rm = T)
aty <- 1:5#axTicks(2)
labels <- sapply(aty,function(i)
  as.expression(bquote(10^ .(i)))
)

# pdf('time_series.pdf',height=5,width=3)
sites<-c(60,85,121)
sites.obs[sites.obs.use][sites]
par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(5,5,5,5))
plot(times.obs,dat2[sites[1],],type='o',pch=16,ylim=ylim,xlab='local time',ylab='',yaxt='n',xaxt='n')
axis(2,at=aty,labels=labels,las=2); mtext('Time Series',3,2)
plot(times.obs,dat2[sites[2],],type='o',pch=16,ylim=ylim,yaxt='n',xlab='local time',ylab='',xaxt='n')
axis(2,at=aty,labels=labels,las=2); mtext(expression(SF[6]~~'concentration'~~(ng/m^3)),2,2.75)
plot(times.obs,dat2[sites[3],],type='o',pch=16,ylim=ylim,yaxt='n',xlab='',ylab='')
axis(2,at=aty,labels=labels,las=2); mtext('local time',1,3)
# dev.off()

# pdf('time_series_sim.pdf',height=5,width=3)
sites<-c(60,85,121)
sites.obs[sites.obs.use][sites]
par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(5,5,5,5))
plot(times.modelRuns,log10(10^dat.sim[sites[1],,1]),type='o',pch=16,ylim=ylim,xlab='local time',ylab='',yaxt='n',xaxt='n')
axis(2,at=aty,labels=labels,las=2); mtext('Time Series',3,2)
plot(times.modelRuns,log10(10^dat.sim[sites[2],,1]),type='o',pch=16,ylim=ylim,xlab='local time',ylab='',yaxt='n',xaxt='n')
axis(2,at=aty,labels=labels,las=2); mtext(expression(SF[6]~~'concentration'~~(ng/m^3)),2,2.75)
plot(times.modelRuns,log10(10^dat.sim[sites[3],,1]),type='o',pch=16,ylim=ylim,xlab='local time',ylab='',yaxt='n')
axis(2,at=aty,labels=labels,las=2); mtext('local time',1,3)
# dev.off()


library(RgoogleMaps)
lat = coordsLL[,2]
lon = coordsLL[,1]
center = c(35.15,-120.64)
zoom <- min(MaxZoom(range(lat), range(lon)))
MyMap <- GetMap(center=center, zoom=zoom, destfile = "MyTile1.png",maptype='terrain') #markers=...for markings
# pdf('obs_loc.pdf')
tmp <- PlotOnStaticMap(MyMap, lat = lat,lon = lon,destfile = "MyTile1.png",cex=.5,col='red')
# dev.off()


library(ggmap)
library(ggplot2)
lat<-coordsLL[,2]
lon<-coordsLL[,1]
map<-get_map(location=c(lon=-120.64,lat=35.15),zoom=MaxZoom(range(lat),range(lon)),maptype='terrain')
zlim<-range(log10(na.omit(dd2[[3]][,6])))
xlim<-range(lon)+c(-.01,.01)
ylim<-range(lat)+c(-.01,.01)

i=1
dato<-log10(dd2[[3]][150*(i-1)+(1:150),6])
dato<-dato[sites.obs.use]
ord<-order(dato)
df<-data.frame(x=lon,y=lat,z=dato)[ord,]
#df<-df[!is.na(dato[ord]),]
quiltplot <- ggmap(map) +
  scale_x_continuous(limits=xlim,expand=c(0,0)) +
  scale_y_continuous(limits=ylim,expand=c(0,0)) +
  geom_point(data=df[is.na(dato[ord]),],aes(x=x,y=y),color='black',size=2,alpha=7/10,inherit.aes=F) +
  geom_point(data=df[!is.na(dato[ord]),],aes(x=x,y=y,color=z),size=3,alpha=7/10,inherit.aes=F) + 
  #scale_colour_gradientn(limits=zlim,colors=c('blue','green','yellow','red'),breaks=1:4, labels=sapply(1:4,function(i) as.expression(bquote(10^ .(i))))) +
  scale_colour_gradientn(limits=zlim,colors=c('darkblue','blue','darkturquoise','green','yellow','orange','red','darkred'),breaks=1:5, labels=sapply(1:5,function(i) as.expression(bquote(10^ .(i))))) +
  labs(x='longitude',y='latitude') + 
  ggtitle(expression('Background'~SF[6])) +
  theme(legend.title=element_blank())
print(quiltplot)

dato<-log10(backg.krig)
ord<-order(dato)
df<-data.frame(x=lon,y=lat,z=dato)[ord,]
quiltplot2 <- ggmap(map) +
  scale_x_continuous(limits=xlim,expand=c(0,0)) +
  scale_y_continuous(limits=ylim,expand=c(0,0)) +
  #geom_point(data=df[is.na(dato[ord]),],aes(x=x,y=y),color='black',size=2,alpha=7/10,inherit.aes=F) +
  geom_point(data=df[!is.na(dato[ord]),],aes(x=x,y=y,color=z),size=3,alpha=7/10,inherit.aes=F) + 
  #scale_colour_gradientn(limits=zlim,colors=c('blue','green','yellow','red'),breaks=1:4, labels=sapply(1:4,function(i) as.expression(bquote(10^ .(i))))) +
  scale_colour_gradientn(limits=zlim,colors=c('darkblue','blue','darkturquoise','green','yellow','orange','red','darkred'),breaks=1:5, labels=sapply(1:5,function(i) as.expression(bquote(10^ .(i))))) +
  labs(x='longitude',y='latitude') + 
  ggtitle(expression('Background'~SF[6])) +
  theme(legend.title=element_blank())
print(quiltplot2)


print(quiltplot)
print(quiltplot2)


}