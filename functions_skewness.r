#!/home/pepflei/R/bin/Rscript
library(moments)
source("load.r")

calc_runskew <- function(nday=91,trendID="91_5",trend_style="_mean",dataset="_TX"){
    # calculates running mean for each grid point
    # can choose between the c script and r function
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    trend=trend_load(paste("../data/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    detrended=dat$tas-trend

    ntot = 1319
    runskew=array(NA,dim=c(ntot,365))
	for (day in 1:365){
		cat("-")
		for (q in 1:ntot){
			if ((day-(nday-1)/2)>0 & (day+(nday-1)/2)<366){
				z=c(detrended[q,(day-(nday-1)/2):(day+(nday-1)/2),])
			}
			if ((day-(nday-1)/2)<1){
				z=c(detrended[q,(365+day-(nday-1)/2):365,],detrended[q,1:(day+(nday-1)/2),])
			}
			if ((day+(nday-1)/2)>365){
				z=c(detrended[q,(day-(nday-1)/2):365,],detrended[q,1:((day+(nday-1)/2)-365),])
			}
			if (which(!is.na(z)>100)){
				runskew[q,day]=skewness(z,na.rm=TRUE)
			}
		}
	}
	write.table(runskew,paste("../data/",trendID,"/",trendID,trend_style,dataset,"_daily_skewness.txt",sep=""))
}

plot_daily_skew_cycle <- function(stations=c(488,52,274,388,467),trendID="91_5",trend_style="_mean",dataset="_TX"){
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
	runskew=read.table(paste("../data/",trendID,"/",trendID,trend_style,dataset,"_daily_skewness.txt",sep=""))
	pdf(file="../plots/zwischenzeugs/skewness/daily_cycle.pdf")
	plot(NA,NA,xlim=c(1,365),ylim=c(-1,1),xlab="day in year",ylab="skewness")
	color=c("red","blue","green","orange","violet")
	x=1:365
	y=x*NA
	data=array(NA,dim=c(1319,365))
	for (day in x){
		data[,day]=runskew[1:1319,day]
	}
	for (sea in c(60,152,244,335)){
		abline(v=sea,col=rgb(0.5,0.5,0.5,0.5))
	}
	for (q in 1:length(stations)){
		lines(x,data[stations[q],],col=color[q])
		#abline(h=mean(data[stations[q],],na.rm=TRUE),col=color[q],lty=2)
	}
	legend("topright",col=color,legend=c(paste("station ID",stations)),lty=c(1,1,1,1))

}

plot_skewness <- function(dat,grid=FALSE,ausschnitt=c(-80,80)){
	reihen=array(NA,dim=c(3,ntot))
	reihen_sig=array(NA,dim=c(3,ntot))
	titel=c("skewness of TX","skewness of TN","skewness(TX)-skewness(TN)")

	for (q in 1:1319){
		cat("-")
		y=as.vector(dat$tas[q,,])
		nona=which(!is.na(y))
		y=y[nona]
		reihen[1,q]=skewness(y)
	}

	#for TN
	dat=dat_load("../data/HadGHCND_TN_data3D.day1-365.1950-2014.nc")

	for (q in 1:1319){
		cat("-")
		y=as.vector(dat$tas[q,,])
		nona=which(!is.na(y))
		y=y[nona]
		reihen[2,q]=skewness(y)
	}

	reihen[3,]=reihen[1,]-reihen[2,]

	map_allgemein(dat=dat,reihen=reihen,titel=titel,farb_mitte="0",
		farb_palette="lila-gruen",
		filename_plot=paste("../plots/zwischenzeugs/skewness/skewness.pdf",sep=""),worldmap=worldmap,ausschnitt=ausschnitt)
}

plot_seasonal_skewness <- function(dat,grid=FALSE,ausschnitt=c(-80,80)){
	# seasonal skewness
	seasonStart=c(60,151,244,335,1)
	seasonStop=c(151,243,334,424,365)
	season_names=c("spring","summer","autumn","winter","year")

	reihen=array(NA,dim=c(9,ntot))
	reihen_sig=array(NA,dim=c(9,ntot))
	titel=c()
	
	for (sea in 1:5){
		for (q in 1:1319){
			if (sea==4){
				z=c(dat$tas[q,1:(seasonStop[sea]-365),],dat$tas[q,seasonStart[sea]:365,])
			}
			else {z=dat$tas[q,seasonStart[sea]:seasonStop[sea],]}

			y=as.vector(z)
			nona=which(!is.na(y))
			y=y[nona]
			reihen[sea,q]=skewness(y)
		}
		titel[sea]=paste("skewness in",season_names[sea])
	}

	for (sea in 1:4){
		reihen[(5+sea),]=reihen[sea,]-reihen[5,]
		titel[(5+sea)]=paste("skewness anomaly in",season_names[sea])
	}

	map_allgemein(dat=dat,reihen=reihen,titel=titel,farb_mitte="0",
		farb_palette="lila-gruen",
		filename_plot=paste("../plots/zwischenzeugs/skewness/skewness_seasonal.pdf",sep=""),worldmap=worldmap,ausschnitt=ausschnitt)
}

plot_detrended_seasonal_skewness <- function(dat,grid=FALSE,ausschnitt=c(-80,80),trendID="91_5",trend_style="_mean",dataset="_TX",additional_style="_median"){
	# detrended seasonal skewness
	seasonStart=c(60,151,244,335,1)
	seasonStop=c(151,243,334,424,365)
	season_names=c("spring","summer","autumn","winter","year")

	reihen=array(NA,dim=c(9,ntot))
	reihen_sig=array(NA,dim=c(9,ntot))
	titel=c()

    trend=trend_load("../data/91_5/91_5_trend_mean_TX.nc")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")

    detrended=dat$tas-trend

	
	for (sea in 1:5){
		for (q in 1:1319){
			if (sea==4){
				z=c(detrended[q,1:(seasonStop[sea]-365),],detrended[q,seasonStart[sea]:365,])
			}
			else {z=detrended[q,seasonStart[sea]:seasonStop[sea],]}

			y=as.vector(z)
			nona=which(!is.na(y))
			y=y[nona]
			reihen[sea,q]=skewness(y)
		}
		titel[sea]=paste("skewness in",season_names[sea])
	}

	for (sea in 1:4){
		reihen[(5+sea),]=reihen[sea,]-reihen[5,]
		titel[(5+sea)]=paste("skewness anomaly in",season_names[sea])
	}

	map_allgemein(dat=dat,reihen=reihen,titel=titel,farb_mitte="0",
		farb_palette="lila-gruen",
		filename_plot=paste("../plots/zwischenzeugs/skewness/skewness_seasonal_detrended_",trendID,trend_style,dataset,additional_style,".pdf",sep=""),worldmap=worldmap,ausschnitt=ausschnitt)
}

#calc_runskew()
plot_daily_skew_cycle()

