library(ncdf)
source("load.r")

land_percentage_bigger_than <- function(var,threshold=0,lat_ausschnitt=c(35,60),lon_ausschnitt=c(-180,180),period="1950-2014"){
	inside_region=which(dat$lon>lon_ausschnitt[1] & dat$lon<lon_ausschnitt[2] & dat$lat>lat_ausschnitt[1] & dat$lat<lat_ausschnitt[2])
	bigger=length(which(var[inside_region]>threshold))
	smaller=length(which(var[inside_region]<threshold))
	equal=length(which(var[inside_region]==threshold))
	print(paste("mean:",mean(var[inside_region],na.rm=TRUE),"sd:",sd(var[inside_region],na.rm=TRUE)))
	print(paste("bigger:",bigger,"smaller:",smaller,"equal:",equal))
	print(paste("percentage bigger:",bigger/(bigger+smaller+equal)))
}

ntot=1319
dataset="_TMean"
trendID="91_5"
additional_style=""

dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))


if (1==2){
	period="1980-2014"
	season="4seasons"

	nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",period,"/",trendID,dataset,"_duration_analysis_",season,".nc",sep=""))
	y=get.var.ncdf(nc,"dur_ana_full")

	for (state in 1:2){
		land_percentage_bigger_than(y[1:ntot,state,8,3],threshold=,lat_ausschnitt=c(60,90))#,lon_ausschnitt=c(0,60))
	}
}

if (1==1){
	period="1950-2014"
	season="4seasons"

	nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",period,"/",trendID,dataset,"_duration_analysis_",season,".nc",sep=""))
	ground=get.var.ncdf(nc,"dur_ana_full")

	season="JJA"
	nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",period,"/",trendID,dataset,"_duration_analysis_",season,".nc",sep=""))
	y=get.var.ncdf(nc,"dur_ana_full")

	y=y-ground

	for (state in 1:2){
		land_percentage_bigger_than(y[1:ntot,state,8,3],threshold=0.0,lat_ausschnitt=c(20,90))#,lon_ausschnitt=c(0,60))
	}
}
