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

	lvls=length(pressure_level)

	time=get.var.ncdf(nc,"time")
	monate=1979+time/12

	u2syn=get.var.ncdf(nc,"u2syn")
	v2syn=get.var.ncdf(nc,"v2syn")

	eke_grided=(u2syn[,,,]+v2syn[,,,])/2

	ntot=1319
	dat$lon[dat$lon<0]=dat$lon[dat$lon<0]+360
	eke_ID=array(NA,dim=c(ntot,lvls,dim(eke_grided)[4]))

	for (q in 1:ntot){
		eke_ID[q,,]=eke_grided[(dat$lon[q]/3.75+1),((dat$lat[q]+90)/2.5+1),,]	
	}

	eke_ID_sea=array(NA,dim=c(ntot,lvls,4,36))

	index=3
	sea=1
	year=1

	while (index<=432-3){
		if (sea==5){
			sea=1
			year=year+1
		}
		for (q in 1:ntot){
			for (lvl in 1:lvls)
				eke_ID_sea[q,lvl,sea,year]=mean(eke_ID[q,lvl,(index):(index+2)],na.rm=TRUE)
		}
		
		index=index+3
		sea=sea+1
	}


    year <- dim.def.ncdf("year",units="year",vals=1979:2014, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    season <- dim.def.ncdf("season",units="uu",vals=1:4,unlim=FALSE)
    levelist <- dim.def.ncdf("levelist",units="mbar",vals=c(850,700,500,400,300,250),unlim=FALSE)

    varlon <- var.def.ncdf(name="lon",units="bla",dim=ID, missval=-9999.0)
    varlat <- var.def.ncdf(name="lat",units="bla",dim=ID, missval=-9999.0)
    vareke <- var.def.ncdf(name="eke",units="(m/s)2",dim=list(ID,levelist,season,year), missval=-9999.0)
    vars=list(varlon,varlat,vareke)
   
    nc = create.ncdf("../data/eke_ID.nc",vars)

	put.var.ncdf(nc,varlon,dat$lon)
	put.var.ncdf(nc,varlat,dat$lat)
	put.var.ncdf(nc,vareke,eke_ID_sea)

    close.ncdf(nc) 
}   

eke_markov_correl <- function(trendID){
	ntot=1319

	nc=open.ncdf(paste("../data/",trendID,"/2_states/markov/",trendID,"_markov_2states.nc",sep=""))
	markov=get.var.ncdf(nc,"markov")

	nc=open.ncdf("../data/eke_ID.nc")
	eke=get.var.ncdf(nc,"eke")

	cor=array(NA,dim=c(ntot,6,4,4))
	cor_sig=array(NA,dim=c(ntot,6,4,4))
	x=seq(1,36,1)
	for (q in 1:ntot){
		print(q)
		for (lvl in 1:6){
			for (sea in 1:4){
				for (trans in 1:4){
					if (length(which(!is.na(markov[q,sea,trans,30:65])))>10 & length(which(!is.na(eke[q,lvl,sea,])))>10){

						#y=detrend(eke_500[q,sea,],x)$detrended
						#z=detrend(markov[q,sea,trans,],x)$detrended

						lr=summary(lm(eke[q,lvl,sea,]~markov[q,sea,trans,30:65]))
						cor[q,lvl,sea,trans]=lr$coefficients[2,1]
						cor_sig[q,lvl,sea,trans]=lr$coefficients[2,4]
					}
				}
			}
		}
	}

    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    season <- dim.def.ncdf("season",units="uu",vals=1:4,unlim=FALSE)
    transition <- dim.def.ncdf("transition",units="uu",vals=1:4,unlim=FALSE)
    levelist <- dim.def.ncdf("levelist",units="mbar",vals=c(850,700,500,400,300,250),unlim=FALSE)

    correl <- var.def.ncdf(name="correl",units="bla",dim=list(ID,levelist,season,transition), missval=-9999.0)
    correl_sig <- var.def.ncdf(name="correl_sig",units="bla",dim=list(ID,levelist,season,transition), missval=-9999.0)
    vars=list(correl,correl_sig)
   
    nc = create.ncdf(paste("../data/",trendID,"/",trendID,"_eke_markov_correl.nc",sep=""),vars)


	put.var.ncdf(nc,correl,cor)
	put.var.ncdf(nc,correl_sig,cor_sig)
}

#create_eke()
eke_markov_correl("91_5")
eke_markov_correl("91_3")
