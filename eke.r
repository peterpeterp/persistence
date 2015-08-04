#!/home/pepflei/R/bin/Rscript

source("write.r")
source("load.r")


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

	eke_ID_year=array(eke_ID,dim=c(ntot,lvls,12,36))

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
    monate <- dim.def.ncdf("monate",units="uu",vals=1:12,unlim=FALSE)
    levelist <- dim.def.ncdf("levelist",units="mbar",vals=c(850,700,500,400,300,250),unlim=FALSE)

    varlon <- var.def.ncdf(name="lon",units="bla",dim=ID, missval=-9999.0)
    varlat <- var.def.ncdf(name="lat",units="bla",dim=ID, missval=-9999.0)
    vareke_sea <- var.def.ncdf(name="eke_sea",units="(m/s)2",dim=list(ID,levelist,season,year), missval=-9999.0)
    vareke_year <- var.def.ncdf(name="eke_year",units="(m/s)2",dim=list(ID,levelist,monate,year), missval=-9999.0)
    vars=list(varlon,varlat,vareke_sea,vareke_year)
   
    nc = create.ncdf("../data/eke_ID.nc",vars)

	put.var.ncdf(nc,varlon,dat$lon)
	put.var.ncdf(nc,varlat,dat$lat)
	put.var.ncdf(nc,vareke_sea,eke_ID_sea)
	put.var.ncdf(nc,vareke_year,eke_ID_year)

    close.ncdf(nc) 
}   

eke_markov_correl <- function(trendID,states,transition_names){
	ntot=1319
	transNumb=states*states

	nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/markov/",trendID,"_markov_",states,"states.nc",sep=""))
	markov=get.var.ncdf(nc,"markov")

	nc=open.ncdf("../data/eke_ID.nc")
	eke=get.var.ncdf(nc,"eke")

	correlation=array(NA,dim=c(ntot,6,4,transNumb))
	corSlope=array(NA,dim=c(ntot,6,4,transNumb))
	corSlope_sig=array(NA,dim=c(ntot,6,4,transNumb))
	x=seq(1,36,1)
	for (q in 488){
		print(q)
		for (lvl in 1:6){
			for (sea in 1:4){
				for (trans in 1:transNumb){
					if (length(which(!is.na(markov[q,sea,trans,30:65])))>10 & length(which(!is.na(eke[q,lvl,sea,])))>10){

						#y=detrend(eke_500[q,sea,],x)$detrended
						#z=detrend(markov[q,sea,trans,],x)$detrended
						if (sd(as.vector(eke[q,lvl,sea,]),na.rm=TRUE)!=0 & sd(as.vector(markov[q,sea,trans,30:65]),na.rm=TRUE)!=0){
							correlation[q,lvl,sea,trans]=cor(x=as.vector(eke[q,lvl,sea,]),y=as.vector(markov[q,sea,trans,30:65]),use="complete")
						}
						lr=summary(lm(eke[q,lvl,sea,]~markov[q,sea,trans,30:65]))
						corSlope[q,lvl,sea,trans]=lr$coefficients[2,1]
						corSlope_sig[q,lvl,sea,trans]=lr$coefficients[2,4]
					}
				}
			}
		}
	}

    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    season <- dim.def.ncdf("season",units="uu",vals=1:4,unlim=FALSE)
    transition <- dim.def.ncdf("transition",units=transition_names,vals=1:transNumb,unlim=FALSE)
    levelist <- dim.def.ncdf("levelist",units="mbar",vals=c(850,700,500,400,300,250),unlim=FALSE)


    varCorrelation <- var.def.ncdf(name="correlation",units="bla",dim=list(ID,levelist,season,transition), missval=-9999.0)
    varCorSlope <- var.def.ncdf(name="corSlope",longname="linear regression like correlation",units="bla",dim=list(ID,levelist,season,transition), missval=-9999.0)
    varCorSlope_sig <- var.def.ncdf(name="corSlope_sig",longname="linear regression like correlation_sig",units="bla",dim=list(ID,levelist,season,transition), missval=-9999.0)
    vars=list(varCorrelation,varCorSlope,varCorSlope_sig)
   
    nc = create.ncdf(paste("../data/",trendID,"/",states,"_states/",trendID,"_eke_markov_",states,"states_correl.nc",sep=""),vars)

	put.var.ncdf(nc,varCorrelation,correlation)
	put.var.ncdf(nc,varCorSlope,corSlope)
	put.var.ncdf(nc,varCorSlope_sig,corSlope_sig)

	close.ncdf(nc)
}

eke_duration_correl <- function(trendID,states,seasons=c("spring","summer","autumn","winter"),taus=c(0.25,0.5,0.75,0.9,0.95,0.98),skips=c(304,329,923,961)){
	library(quantreg)
	ntot=1319
	monatStart=c(59,91,120,151,181,212,243,273,304,334,366,396,427)
	monatStart=monatStart/365

	for (sea in 1:length(seasons)){
		monatSelection=monatStart[((sea-1)*3):(sea*3+1)]

		nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/duration/",trendID,"_duration_",states,"s_",seasons[sea],".nc",sep=""))
		duration=get.var.ncdf(nc,"dur")
		duration_mid=get.var.ncdf(nc,"dur_mid")



		nc=open.ncdf("../data/eke_ID.nc")
		eke=get.var.ncdf(nc,"eke_year")

		cor=array(NA,dim=c(ntot,6,4,2,length(taus)))
		cor_sig=array(NA,dim=c(ntot,6,4,2,length(taus)))
		x=seq(1,36,1)
		for (q in 1:ntot){
			if (length(which(skips==q))<2){
				print(paste(sea,q))
				for (state in 1:2){
					size=length(which(duration_mid[q,state,]>1979 & duration_mid[q,state,]<21014))
					if (size>30){    #length(which(!is.na(eke[q,1:6,sea,])))>10
						dur=array(NA,dim=c(size))
						dur_mid=array(NA,dim=c(size))
						eke_ext=array(NA,dim=c(6,size))
						count=0
						for (year in 30:64){
							year=year+1949
							for (month in 1:3){
								if (count!=size){
									inside=which(duration_mid[q,state,]>(year+monatSelection[month]) & duration_mid[q,state,]<(year+monatSelection[(month+1)]))
									if (length(inside)>0){
										dur[(count+1):(count+length(inside))]=duration[q,state,inside]
										dur_mid[(count+1):(count+length(inside))]=duration_mid[q,state,inside]
										eke_ext[1:6,(count+1):(count+length(inside))]=eke[q,1:6,month,(year-1978)]
										count=count+length(inside)
									}
								}
							}
						}

						for (lvl in 1:6){
							#print("------------------------------")
							#print(lvl)
							#print(as.vector(eke_ext[lvl,]))
							#print(as.vector(dur))
							#pdf(file=paste("../plots/lvl",lvl,".pdf",sep=""))
							#plot(as.vector(eke_ext[lvl,]),as.vector(dur))
							#abline(rq(as.vector(dur)~as.vector(eke_ext[lvl,]),0.95),col="blue")
							#abline(rq(as.vector(dur)~as.vector(eke_ext[lvl,]),0.98),col="red")
							#abline(rq(as.vector(dur)~as.vector(eke_ext[lvl,]),0.5),col="green")
							#graphics.off()
							#print(ecdf(dur))
							#print(rq(as.vector(dur)~as.vector(eke_ext[lvl,]),taus))
							#print(dur)
							#return(list(dur=dur,eke_ext=eke_ext,taus=taus,lvl=lvl))
							#eke_ext[lvl,]=eke_ext[lvl,]#+seq(-0.0001,0.0001,length=3)
							sf=try(summary(rq(as.vector(dur)~as.vector(eke_ext[lvl,]),taus),se="nid"))
							#print(sf)
				            if (class(sf)=="try-error"){
				                #print(sf)
				                for (i in 1:length(taus)){
				                    sf=try(summary(rq(as.vector(dur)~as.vector(eke_ext[lvl,]),taus[i]),se="nid"))
				                    if (class(sf)=="try-error"){
				                        cor[q,lvl,sea,state,i]=NA
				                        cor_sig[q,lvl,sea,state,i]=NA
				                    }
				                    else {
				                        cor[q,lvl,sea,state,i]=sf$coefficients[2,1]
				                        cor_sig[q,lvl,sea,state,i]=sf$coefficients[2,4]
				                    }
				                }
				            }
				            else {
				                slope=sapply(sf, function(x) c(tau=x$tau, x$coefficients[-1,]))

				                cor[q,lvl,sea,state,1:length(taus)]=slope[2,1:length(taus)]
				                cor_sig[q,lvl,sea,state,1:length(taus)]=slope[5,1:length(taus)]

				            }
			            #print(cor[q,lvl,sea,state,1:length(taus)])
						}

					}
				}
			}
		}

	}

	ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
	season <- dim.def.ncdf("season",units="uu",vals=1:4,unlim=FALSE)
	states <- dim.def.ncdf("states",units="uu",vals=1:2,unlim=FALSE)
	levelist <- dim.def.ncdf("levelist",units="mbar",vals=c(850,700,500,400,300,250),unlim=FALSE)
	quantiles <- dim.def.ncdf("quantiles",units="0-1",vals=taus,unlim=FALSE)

	correl <- var.def.ncdf(name="quantile regression like correlation",units="bla",dim=list(ID,levelist,season,states,quantiles), missval=-9999.0)
	correl_sig <- var.def.ncdf(name="quantile regression like correlation_sig",units="bla",dim=list(ID,levelist,season,states,quantiles), missval=-9999.0)
	vars=list(correl,correl_sig)
	   
	nc = create.ncdf(paste("../data/",trendID,"/",states,"_states/",trendID,"_eke_duration_cor_",states,"states_",paste(seasons,sep="-"),".nc",sep=""),vars)


	put.var.ncdf(nc,correl,cor)
	put.var.ncdf(nc,correl_sig,cor_sig)
}

#create_eke()

#eke_markov_correl("91_5",states=3,transition_names="cc nc wc cn nn wn cw nw ww")
#eke_markov_correl("91_5",states=2,transition_names="cc wc cw ww")
#eke_markov_correl("91_3")
out=eke_duration_correl("91_5",states=2)
#lvl=out$lvl
#dur=out$dur
#taus=out$taus
##eke_ext=out$eke_ext
#eke_ext[lvl,]=eke_ext[lvl,]+seq(-0.00001,0.00001,length=2)
