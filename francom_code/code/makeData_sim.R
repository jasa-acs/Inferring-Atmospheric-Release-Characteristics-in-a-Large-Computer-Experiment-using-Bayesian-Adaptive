## make simulation data, get EOFs

## get inputs
inputs<-read.table('data/flexuq_lhs04-wrfdims.tab',header = T)[1:18000,]

## get FLEXPART outputs
library(ncdf4)
con<-nc_open('data/flxout-dopptex-timeseries-lhs04.nc')
#print(con)
dopsites.sim<-ncvar_get(con,varid='DOPPTEX_Sites',start=c(1,1),count=c(137,1))

# all data
dat.sim<-ncvar_get(con,start = c(1,1,1),count=c(137,34,18000)) # 137 spatial locations, 34 time points, 18000 model runs
dat.sim<-log10(dat.sim+50) # 50 background
nc_close(con)

# reshape
Y.sim<-matrix(nrow=137*34,ncol=18000)
for(i in 1:18000){
  Y.sim[,i]<-c(dat.sim[,,i])
}

# holdout 10%
holdout.sim<-sample(1:18000,1800) #ch

# get EOFs using svd (time consuming)
temp<-svd(Y.sim[,-holdout.sim]) # non-centered on purpose
U<-temp$u
V<-temp$v
D<-diag(temp$d)
n.eof<-200
rm(temp)

# rename and save
X.sim<-inputs[,-c(1,2)]
reduced.y.sim<-D[1:n.eof,1:n.eof]%*%t(V[,1:n.eof])
eofs.sim<-U[,1:n.eof]

save(X.sim,holdout.sim,reduced.y.sim,eofs.sim,file='flex_reducedDat.rda')
save(Y.sim,file='flex_dat.rda')
# 
# rm(list=ls())
# gc()