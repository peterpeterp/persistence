source("write.r")
source("load.r")


plot_eke_mar <- function(){
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


	plot(eke_detrend,mar_detrend[(length(mar_detrend)-length(eke_detrend)+1):length(mar_detrend)])
	abline(lm(mar_detrend[(length(mar_detrend)-length(eke_detrend)+1):length(mar_detrend)]~eke_detrend))
}


detrend <- function(y,x){
	lr=summary(lm(y~x))$coefficients
	detrend=y-lr[1,1]-x*lr[2,1]
	return(list(detrended=detrend,lr=lr))
}

create_eke <- function(){
	dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
	nc=open.ncdf("../data/sonstiges/eke/EKE_ERA_Interim_1979-2014_calendar96x72.nc")
	pressure_level=get.var.ncdf(nc,"levelist")

	time=get.var.ncdf(nc,"time")
	monate=1979+time/12

	u2syn=get.var.ncdf(nc,"u2syn")
	v2syn=get.var.ncdf(nc,"v2syn")

	eke_850_all=(u2syn[,,1,]+v2syn[,,1,])/2
	eke_500_all=(u2syn[,,3,]+v2syn[,,3,])/2


	ntot=1319
	dat$lon[dat$lon<0]=dat$lon[dat$lon<0]+360
	eke_850=array(NA,dim=c(ntot,dim(eke_850_all)[3]))
	eke_500=array(NA,dim=c(ntot,dim(eke_500_all)[3]))

	for (q in 1:ntot){
		eke_850[q,]=eke_850_all[(dat$lon[q]/3.75+1),((dat$lat[q]+90)/2.5+1),]	
		eke_500[q,]=eke_500_all[(dat$lon[q]/3.75+1),((dat$lat[q]+90)/2.5+1),]	
	}

	eke_850_sea=array(NA,dim=c(ntot,4,36))
	eke_500_sea=array(NA,dim=c(ntot,4,36))

	index=3
	sea=1
	year=1

	while (index<=432-3){
		if (sea==5){
			sea=1
			year=year+1
		}
		for (q in 1:ntot){
			eke_850_sea[q,sea,year]=mean(eke_850[q,(index):(index+2)],na.rm=TRUE)
			eke_500_sea[q,sea,year]=mean(eke_500[q,(index):(index+2)],na.rm=TRUE)
		}
		
		index=index+3
		sea=sea+1
	}

	print(eke_850[488,])
	print(eke_850_sea[488,,])
	print(eke_500[488,])
	print(eke_500_sea[488,,])

    year <- dim.def.ncdf("year",units="year",vals=1979:2014, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    season <- dim.def.ncdf("season",units="uu",vals=1:4,unlim=FALSE)

    varlon <- var.def.ncdf(name="lon",units="bla",dim=ID, missval=-9999.0)
    varlat <- var.def.ncdf(name="lat",units="bla",dim=ID, missval=-9999.0)
    vareke850 <- var.def.ncdf(name="eke_850",units="bla",dim=list(ID,season,year), missval=-9999.0)
    vareke500 <- var.def.ncdf(name="eke_500",units="bla",dim=list(ID,season,year), missval=-9999.0)
    vars=list(varlon,varlat,vareke850,vareke500)
   
    nc = create.ncdf("../data/sonstiges/eke_ID.nc",vars)

	put.var.ncdf(nc,varlon,dat$lon)
	put.var.ncdf(nc,varlat,dat$lat)
	put.var.ncdf(nc,vareke850,eke_850_sea)
	put.var.ncdf(nc,vareke500,eke_500_sea)

    close.ncdf(nc) 
}   

#create_eke()

eke_markov_correl <- function(){
	ntot=1319

	nc=open.ncdf("../data/91_5/2_states/markov/91_5_markov2s.nc")
	markov_spring=get.var.ncdf(nc,"markov_spring")
	markov_summer=get.var.ncdf(nc,"markov_summer")
	markov_autumn=get.var.ncdf(nc,"markov_autumn")
	markov_winter=get.var.ncdf(nc,"markov_winter")

	markov=array(NA,dim=c(ntot,4,4,36))
	markov[1:ntot,1,1:4,1:36]=markov_spring[1:ntot,1:4,30:65]
	markov[1:ntot,2,1:4,1:36]=markov_summer[1:ntot,1:4,30:65]
	markov[1:ntot,3,1:4,1:36]=markov_autumn[1:ntot,1:4,30:65]
	markov[1:ntot,4,1:4,1:36]=markov_winter[1:ntot,1:4,30:65]

	nc=open.ncdf("../data/eke_ID.nc")
	eke_500=get.var.ncdf(nc,"eke_500")
	eke_850=get.var.ncdf(nc,"eke_850")

	cor=array(NA,dim=c(ntot,4,4,2))
	cor_sig=array(NA,dim=c(ntot,4,4,2))
	x=seq(1,36,1)
	for (q in 1:ntot){
		print(q)
		for (sea in 1:4){
			for (trans in 1:4){
				if (length(which(!is.na(markov[q,sea,trans,])))>10 & length(which(!is.na(eke_500[q,sea,])))>10 & length(which(!is.na(eke_850[q,sea,])))>10){

					y=detrend(eke_500[q,sea,],x)$detrended
					z=detrend(markov[q,sea,trans,],x)$detrended

					lr=summary(lm(y~z))
					cor[q,sea,trans,1]=lr$coefficients[2,1]
					cor_sig[q,sea,trans,1]=lr$coefficients[2,4]

					y=detrend(eke_850[q,sea,],x)$detrended
					lr=summary(lm(y~z))
					cor[q,sea,trans,2]=lr$coefficients[2,1]
					cor_sig[q,sea,trans,2]=lr$coefficients[2,4]
				}
			}
		}
	}

	print(cor[1:ntot,2,4,2])

    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    season <- dim.def.ncdf("season",units="uu",vals=1:4,unlim=FALSE)
    transition <- dim.def.ncdf("transition",units="uu",vals=1:4,unlim=FALSE)
    eke <- dim.def.ncdf("eke",units="1 - 500 , 2 - 850",vals=1:2,unlim=FALSE)

    correl <- var.def.ncdf(name="correl",units="bla",dim=list(ID,season,transition,eke), missval=-9999.0)
    correl_sig <- var.def.ncdf(name="correl_sig",units="bla",dim=list(ID,season,transition,eke), missval=-9999.0)
    vars=list(correl,correl_sig)
   
    nc = create.ncdf("../data/sonstiges/eke_markov_correl.nc",vars)


	put.var.ncdf(nc,correl,cor)
	put.var.ncdf(nc,correl_sig,cor_sig)
}

eke_markov_correl_check <- function(){
	ntot=1319

	nc=open.ncdf("../data/91_5/2_states/markov/91_5_markov2s.nc")
	markov_spring=get.var.ncdf(nc,"markov_spring")
	markov_summer=get.var.ncdf(nc,"markov_summer")

	x=seq(1,36,1)

	pdf(file="../plots/adasdas.pdf")
	plot(NA,xlim=c(1,36),ylim=c(0,1))
	lines(x,markov_summer[488,1,30:65],col="blue")
	lines(x,markov_summer[488,3,30:65],col="red")
	plot(NA,xlim=c(1,36),ylim=c(0,1))
	lines(x,markov_summer[488,4,30:65],col="blue")
	lines(x,markov_summer[488,2,30:65],col="red")

	markov_autumn=get.var.ncdf(nc,"markov_autumn")
	markov_winter=get.var.ncdf(nc,"markov_winter")

	markov=array(NA,dim=c(ntot,4,4,36))
	markov[1:ntot,1,1:4,1:36]=markov_spring[1:ntot,1:4,30:65]
	markov[1:ntot,2,1:4,1:36]=markov_summer[1:ntot,1:4,30:65]
	markov[1:ntot,3,1:4,1:36]=markov_autumn[1:ntot,1:4,30:65]
	markov[1:ntot,4,1:4,1:36]=markov_winter[1:ntot,1:4,30:65]

	nc=open.ncdf("../data/eke_ID.nc")
	eke_500=get.var.ncdf(nc,"eke_500")
	eke_850=get.var.ncdf(nc,"eke_850")

	cor=array(NA,dim=c(ntot,4,4,2))
	cor_sig=array(NA,dim=c(ntot,4,4,2))
	x=seq(1,36,1)
	pdf(file="../plots/teetet.pdf")
	for (q in 488){
		print(q)
		for (sea in 1:4){
			for (trans in 1:4){
				if (length(which(!is.na(markov[q,sea,trans,])))>10 & length(which(!is.na(eke_500[q,sea,])))>10 & length(which(!is.na(eke_850[q,sea,])))>10){


					y=detrend(eke_500[q,sea,],x)$detrended
					z=detrend(markov[q,sea,trans,],x)$detrended
					plot(x,markov[q,sea,trans,],main=paste(sea,trans))
					lines(x,markov[q,sea,trans,])
					#plot(x,z,main=paste(sea,trans))
					#plot(z,y,main=paste(sea,trans))

					lr=summary(lm(y~z))
					#abline(lm(y~z))

					cor[q,sea,trans,1]=lr$coefficients[2,1]
					cor_sig[q,sea,trans,1]=lr$coefficients[2,4]

					y=detrend(eke_850[q,sea,],x)$detrended
					#plot(z,y,main=paste(sea,trans))

					lr=summary(lm(y~z))
					#abline(lm(y~z))

					cor[q,sea,trans,2]=lr$coefficients[2,1]
					cor_sig[q,sea,trans,2]=lr$coefficients[2,4]
				}
			}
		}
	}

}

eke_markov_correl_check()
#eke_markov_correl()