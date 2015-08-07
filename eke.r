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

	eke_ID_sea=array(NA,dim=c(ntot,lvls,5,36))

	index=3
	sea=1
	year=1

	while (index<=432-3){
		if (sea==5){
			sea=1
			year=year+1
		}
		for (q in 1:ntot){
			for (lvl in 1:lvls){
				eke_ID_sea[q,lvl,sea,year]=mean(eke_ID[q,lvl,(index):(index+2)],na.rm=TRUE)
			}
		}
		
		index=index+3
		sea=sea+1
	}
	for (q in 1:ntot){
		for (year in 1:36){
			for (lvl in 1:lvls){
				eke_ID_sea[q,lvl,5,year]=mean(eke_ID[q,lvl,((year-1)*12+1):(year*12)],na.rm=TRUE)
			}
		}
	}

    year <- dim.def.ncdf("year",units="year",vals=1979:2014, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    season <- dim.def.ncdf("season",units="uu",vals=1:5,unlim=FALSE)
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

analyse_eke <- function(yearPeriod=c(1979,2014),yearshift=1978){
    # analyse persistence 2 states
    source("functions_support.r")
    nc_orig=open.ncdf(paste("../data/eke_ID.nc",sep=""))
    ntot=1319
    seasons=c("spring","summer","autumn","winter","year")
    for (sea in 1:length(seasons)){
    	season=seasons[sea]
        print(season)
        eke=get.var.ncdf(nc_orig,"eke_sea")[1:ntot,1:6,sea,]
        tmp=global_analysis(toAna=eke,yearPeriod=yearPeriod,yearshift=yearshift)

	    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
	    levelist <- dim.def.ncdf("levelist",units="mbar",vals=c(850,700,500,400,300,250),unlim=FALSE)

	    mean <- var.def.ncdf(name="mean",units="bla",longname=paste("mean",season),dim=list(ID,levelist), missval=-9999.0)
	    std <- var.def.ncdf(name="std",units="bla",longname=paste("standard deviation",season),dim=list(ID,levelist), missval=-9999.0)

	    MK <- var.def.ncdf(name="MK",units="bla",longname=paste("MK",season),dim=list(ID,levelist), missval=-9999.0)
	    MK_sig <- var.def.ncdf(name="MK_sig",units="bla",longname=paste("MK_sig",season),dim=list(ID,levelist), missval=-9999.0)
	    LR <- var.def.ncdf(name="LR",units="bla",longname=paste("LR",season),dim=list(ID,levelist), missval=-9999.0)
	    LR_sig <- var.def.ncdf(name="LR_sig",units="bla",longname=paste("LR_sig",season),dim=list(ID,levelist), missval=-9999.0)

	    vars=list(mean,std,MK,MK_sig,LR,LR_sig)
   
	    nc = create.ncdf(paste("../data/eke_ID_",season,".nc",sep=""),vars)
	    print(dim(tmp))
	    for (i in 1:6){
	        put.var.ncdf(nc,vars[[i]],tmp[1:ntot,1:dim(tmp)[2],i])  
	    }

	    close.ncdf(nc) 
    }
}

eke_markov_correl <- function(trendID,states,transition_names,stations=seq(1,1319,1),plot=FALSE){
	ntot=1319
	transNumb=states*states

	nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/markov/",trendID,"_markov_",states,"states.nc",sep=""))
	markov=get.var.ncdf(nc,"markov")

	nc=open.ncdf("../data/eke_ID.nc")
	eke=get.var.ncdf(nc,"eke_sea")

	# needed for the plots
	if (plot==TRUE){
		pdf(file=paste("../plots/",trendID,"/",states,"_states/station/mar_eke_",stations,".pdf",sep=""))
		pressure_level=get.var.ncdf(nc,"levelist")
		seasons=c("spring","summer","autumn","winter")
		if (states==2){
			transLongName=c("cold to cold","warm to cold","cold to warm","warm to warm")
			color=c("blue","violet","green","red")
		}
		if (states==3){
			transLongName=c("cold to cold","normal to cold","warm to cold","cold to normal","normal to normal","warm to normal","cold to warm","normal to warm","warm to warm")
		}
	}	

	correlation=array(NA,dim=c(ntot,6,4,transNumb))
	corSlope=array(NA,dim=c(ntot,6,4,transNumb))
	corSlope_sig=array(NA,dim=c(ntot,6,4,transNumb))
	x=seq(1,36,1)
	for (q in stations){
		print(q)
		for (lvl in 1:6){
			for (sea in 1:4){
				if (plot==TRUE){
					plot(NA,ylim=c(0,1.1),xlim=c(min(eke[q,lvl,sea,],na.rm=TRUE),max(eke[q,lvl,sea,],na.rm=TRUE)),
						,xlab="eddy kinetic energy",ylab="transition probability",
						main=paste(seasons[sea],"- pressure:",pressure_level[lvl],"mbar"))
					title=c()
					order=order(eke[q,lvl,sea,1:33])
					year=1979:2014
					year=year[order]
				}
				for (trans in 1:transNumb){
					if (length(which(!is.na(markov[q,sea,trans,30:65])))>10 & length(which(!is.na(eke[q,lvl,sea,])))>10){

						if (sd(as.vector(eke[q,lvl,sea,]),na.rm=TRUE)!=0 & sd(as.vector(markov[q,sea,trans,30:65]),na.rm=TRUE)!=0){
							correlation[q,lvl,sea,trans]=cor(x=as.vector(eke[q,lvl,sea,]),y=as.vector(markov[q,sea,trans,30:65]),use="complete")
						}
						lr=summary(lm(eke[q,lvl,sea,]~markov[q,sea,trans,30:65]))
						corSlope[q,lvl,sea,trans]=lr$coefficients[2,1]
						corSlope_sig[q,lvl,sea,trans]=lr$coefficients[2,4]

						if (plot==TRUE){
							lines(eke[q,lvl,sea,order],markov[q,sea,trans,(order+29)],col=color[trans])
							abline(lm(markov[q,sea,trans,30:65]~eke[q,lvl,sea,]),col=color[trans],lty=2)
							title[trans]=paste(transLongName[trans],"transition - EKE correlation =",correlation[q,lvl,sea,trans])
						}
					}
				}
				if (plot==TRUE){
					for (y in 1:36){
						abline(v=as.vector(eke[q,lvl,sea,order])[y],col="grey",lty=2)
						text(as.vector(eke[q,lvl,sea,order])[y],(0.5+0.05*(-1)^y),label=year[y],srt=90,cex=0.5)
					}
					legend("topright",legend=title,lty=array(1,transNumb),col=color)

				}
			}
		}
	}

	if (plot!=TRUE){
	    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
	    season <- dim.def.ncdf("season",units="uu",vals=1:4,unlim=FALSE)
	    transition <- dim.def.ncdf("transition",units=transition_names,vals=1:transNumb,unlim=FALSE)
	    levelist <- dim.def.ncdf("levelist",units="mbar",vals=c(850,700,500,400,300,250),unlim=FALSE)

	    varCorrelation <- var.def.ncdf(name="correlation",units="bla",dim=list(ID,levelist,season,transition), missval=-9999.0)
	    varCorSlope <- var.def.ncdf(name="corSlope",longname="linear regression like correlation",units="bla",dim=list(ID,levelist,season,transition), missval=-9999.0)
	    varCorSlope_sig <- var.def.ncdf(name="corSlope_sig",longname="linear regression like correlation_sig",units="bla",dim=list(ID,levelist,season,transition), missval=-9999.0)
	    
	    vars=list(varCorrelation,varCorSlope,varCorSlope_sig)  
	    nc = create.ncdf(paste("../data/",trendID,"/",states,"_states/",trendID,"_eke_markov_cor_",states,"states.nc",sep=""),vars)

		put.var.ncdf(nc,varCorrelation,correlation)
		put.var.ncdf(nc,varCorSlope,corSlope)
		put.var.ncdf(nc,varCorSlope_sig,corSlope_sig)

		close.ncdf(nc)
	}
}

eke_duration_correl <- function(trendID,states,seasons=c("spring","summer","autumn","winter"),taus=c(0.25,0.5,0.75,0.9,0.95,0.98),stations=seq(1,1319,1),plot=FALSE){
	library(quantreg)
	ntot=1319
	correlation=array(NA,dim=c(ntot,6,4,states,length(taus)))
	nc=open.ncdf("../data/eke_ID.nc")
	eke=get.var.ncdf(nc,"eke_sea")

	# needed for the plots
	if (plot==TRUE){
		pdf(file=paste("../plots/",trendID,"/",states,"_states/stations/dur_eke_",stations,".pdf",sep=""))
		pressure_level=get.var.ncdf(nc,"levelist")
		color=c("lightblue","blue","green","yellow","red","violet")
		if (states==2){
			state_names=c("cold","warm")
		}
		if (states==3){
			state_names=c("cold","normal","warm")
		}
	}

	for (sea in 1:length(seasons)){
		nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/duration/",trendID,"_duration_",states,"s_",seasons[sea],".nc",sep=""))
		duration=get.var.ncdf(nc,"dur")
		duration_mid=get.var.ncdf(nc,"dur_mid")
		print(nc)

		x=seq(1,36,1)
		for (q in stations){
			cat(paste("--",q))
			for (state in 1:states){
				size=length(which(duration_mid[q,state,]>1979 & duration_mid[q,state,]<2014))
				if (size>30 & length(which(!is.na(eke[q,1:6,sea,])))>10){
					durQu=array(NA,dim=c(36,length(taus)))
					eke_ext=array(NA,dim=c(6,size))
					dur_ext=array(NA,size)
					count=0
					for (i in 1:35){
						if (count!=size){
							year=i+29+1949
							# find durations in the selected year (and season)
							# determine quantile "postions" durQu
							inside=which(duration_mid[q,state,]>year & duration_mid[q,state,]<(year+1))
							if (length(inside)>0){
								durQu[i,]=quantile(duration[q,state,inside],taus,na.rm=TRUE)
								# this is just needed for plots
								eke_ext[1:6,(count+1):(count+length(inside))]=eke[q,1:6,sea,(year-1978)]
								dur_ext[(count+1):(count+length(inside))]=duration[q,state,inside]
								count=count+length(inside)
							}
						}
					}
					# calculate correlation between quantile "positions and eke"
					for (lvl in 1:6){
						for (qu in 1:length(taus)){
							if (sd(as.vector(eke[q,lvl,sea,]),na.rm=TRUE)!=0 & sd(as.vector(durQu[1:36,qu]),na.rm=TRUE)!=0){
								correlation[q,lvl,sea,state,qu]=cor(x=as.vector(eke[q,lvl,sea,]),y=as.vector(durQu[1:36,qu]),use="complete")
							}
						}
					}
					# make little explanatory plot ---------------------------------
					if (plot==TRUE){
						for (lvl in 1:6){
							order=order(eke[q,lvl,sea,1:33])
							year=1979:2014
							year=year[order]
							title=c()
							plot(as.vector(eke_ext[lvl,]),dur_ext,xlab="eddy kinetic energy",ylab="duration length",ylim=c(-2,max(dur_ext,na.rm=TRUE)),
								main=paste(seasons[sea],"- pressure:",pressure_level[lvl],"mbar -",state_names[state],"durations"))
							for (y in 1:36){
								abline(v=as.vector(eke[q,lvl,sea,order])[y],col="grey",lty=2)
								text(as.vector(eke[q,lvl,sea,order])[y],(-1.5+0.7*(-1)^y),label=year[y],srt=90,cex=0.5)
							}
							for (qu in 1:length(taus)){
								title[qu]=paste(taus[qu],"quantile cor =",correlation[q,lvl,sea,state,qu])
								lines(as.vector(eke[q,lvl,sea,order]),as.vector(durQu[order,qu]),col=color[qu],
									main=paste(seasons[sea],"- pressure:",pressure_level[lvl],"mbar -",state_names[state],"duration",taus[qu],"quantile"))
							}
							legend("topright",col=color,lty=array(1,length(taus)),legend=title)
						}
					}
					# ------------------------------------------------------------
				}
			}
		}
	}
	if (plot!=TRUE){
		ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
		season <- dim.def.ncdf("season",units="uu",vals=1:4,unlim=FALSE)
		varstates <- dim.def.ncdf("states",units="uu",vals=1:states,unlim=FALSE)
		levelist <- dim.def.ncdf("levelist",units="mbar",vals=c(850,700,500,400,300,250),unlim=FALSE)
		quantiles <- dim.def.ncdf("quantiles",units="0-1",vals=taus,unlim=FALSE)

		varCorrelation <- var.def.ncdf(name="correlation",longname="correaltion between quantile values and eke values in season",units="bla",dim=list(ID,levelist,season,varstates,quantiles), missval=-9999.0)
		vars=list(varCorrelation)
			 
		nc = create.ncdf(paste("../data/",trendID,"/",states,"_states/",trendID,"_eke_duration_cor_",states,"states.nc",sep=""),vars)

		put.var.ncdf(nc,varCorrelation,correlation)
	}
	graphics.off()

}

#create_eke()
analyse_eke()

#eke_markov_correl("91_5",states=3,transition_names="cc nc wc cn nn wn cw nw ww")
#eke_markov_correl("91_5",states=2,transition_names="cc wc cw ww",stations=163,plot=TRUE)


#eke_duration_correl("91_5",states=3,stations=488,plot=TRUE)
#eke_duration_correl("91_3",states=2)

#eke_duration_correl("91_3",states=3)
