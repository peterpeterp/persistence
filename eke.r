source("write.r")
source("load.r")

#dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
#print(dat$lon[488])
#print(dat$lat[488])


nc=open.ncdf("../data/sonstiges/eke/EKE_ERA_Interim_1979-2014_calendar96x72.nc")
print(nc)
pressure_level=get.var.ncdf(nc,"levelist")
print(pressure_level)
time=get.var.ncdf(nc,"time")
monate=1979+time/12

u2syn=get.var.ncdf(nc,"u2syn")
v2syn=get.var.ncdf(nc,"v2syn")

eke_850=(u2syn[,,1,]+v2syn[,,1,])/2

if (1==2){
	pdf(file="../plots/sonstiges/eke")
	print(length(monate))
	print(dim(eke_850))
	plot(monate,eke_850[1,18,])
	abline(lm(eke_850[1,18,]~monate),col="red")
	print(summary(lm(eke_850[1,18,]~monate)))
}


eke_850_season=array(NA,dim=c(96,72,4,36))
index=3
sea=1
jahr=1
for (i in 1:(length(monate)/3-1)){
	if (sea==5){
		sea=1
		jahr=jahr+1
	}
	#cat(i,jahr,sea,index,mean(eke_850[,,index:(index+2)]),"\n")
	eke_850_season[,,sea,jahr]=mean(eke_850[,,index:(index+2)])
	index=index+3
	sea=sea+1
}

detrended <- function(y,x){
	lr=summary(lm(y~x))$coefficients
	detrend=y-lr[1,1]-x*lr[2,1]
	return(list(detrended=detrend,lr=lr))
}





pdf(file="../plots/sonstiges/eke")
x=(1979:2014)
plot(1979:2014,eke_850_season[1,18,1,],lty=1,xlim=c(1980,2014))
abline(lm(eke_850_season[1,18,1,]~x),col="red")

nc=open.ncdf("../data/91_5/2_states/markov/91_5_markov2s.nc")
markov_summer=get.var.ncdf(nc,"markov_summer")[488,4,]
x=get.var.ncdf(nc,"year")+1949
par(new=TRUE)
plot(x,markov_summer,col="green",lty=1,xlim=c(1980,2014))
abline(lm(markov_summer~x),col="violet")

eke_detrend=detrended(eke_850_season[1,18,1,],(1979:2014))$detrended
x=get.var.ncdf(nc,"year")+1949
mar_detrend=detrended(markov_summer,x)$detrended

print(eke_detrend)
print(mar_detrend	)
x=(1979:2014)
plot(1979:2014,eke_detrend,lty=1,xlim=c(1980,2014))

x=get.var.ncdf(nc,"year")+1949
par(new=TRUE)
plot(x,mar_detrend,col="green",lty=1,xlim=c(1980,2014))
