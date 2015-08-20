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



#create_eke()
#analyse_eke()

#eke_markov_correl("91_5",states=3,transition_names="cc nc wc cn nn wn cw nw ww")
eke_markov_correl("91_5",states=2,transition_names="cc wc cw ww",stations=488,plot=TRUE)


#eke_duration_correl("91_5",states=3,stations=488,plot=TRUE)
#eke_duration_correl("91_3",states=2)

#eke_duration_correl("91_3",states=3)
