


meridional_average <- function(x,y,color,lineStyle,abschnitt){
	yAv=x*NA
	for (k in 1:length(x)){
		drauf=which(dat$lon==x[k] & dat$lat<abschnitt[2] & dat$lat>abschnitt[1] & !is.na(y))
		yAv[k]=mean(y[drauf])
	}

	noNas=which(!is.na(yAv))
	x[x<0]=x[x<0]+360
	lines(x[noNas],yAv[noNas],col=color,lty=lineStyle)
	#points(x,yAv,col=color,pch=1,cex=0.5)
}

plot_meridonal_average <- function(trendID,states,period,abschnitt){
	ntot=1319
	if (states==2){
		state_names=c("cold","warm")
		transition_names=c("cold to cold",NA,NA,"warm to warm")
		state_auswahl=c(1,4)
	}
	if (states==3){
		state_names=c("cold","normal","warm")
		state_auswahl=c(1,5,9)
	}
    seasons=c("spring","summer","autumn","winter","year")	
    colors=c(rgb(0.3,0.9,0.4),rgb(1,0.5,0.5),rgb(0.5,0.3,0.5),rgb(0.3,0.5,0.9),rgb(0,0,0))	

    pdf(file=paste("../plots/",trendID,"/",states,"_states/longitudinal/",trendID,"_",states,"_",abschnitt[1],"N_",abschnitt[2],"N.pdf",sep=""))

    mar=array(NA,dim=c(length(seasons),ntot,states*states))
    for (sea in 1:length(seasons)){
		nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/markov/",period,"/",trendID,"_mar",states,"s_trend_",seasons[sea],".nc",sep=""))
		mar[sea,,]=get.var.ncdf(nc,"mean")
	}

	par(mar=c(5, 4, 4, 6) + 0.1)
	x=c(seq(0,180,3.75),seq(-180,-3.75,3.75))

	for (state in state_auswahl){
		plot(NA,xlim=c(0,360),ylim=c(min(mar[,,state],na.rm=TRUE),max(mar[,,state],na.rm=TRUE)),
			xlab="longitude",ylab="transition probability",
			main=paste("meridional average between ",abschnitt[1],"N and ",abschnitt[2],"N \n for markov transition ",transition_names[state],sep=""))
		for (sea in 1:length(seasons)){
			meridional_average(x,mar[sea,,state],colors[sea],1,abschnitt)
		}
		legend("bottomleft",lty=array(1,length(seasons)),col=colors,legend=seasons)
	}

    mar=array(NA,dim=c(length(seasons),ntot,states*states))
    for (sea in 1:length(seasons)){
		nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/markov/",period,"/",trendID,"_mar",states,"s_trend_",seasons[sea],".nc",sep=""))
		mar[sea,,]=get.var.ncdf(nc,"MK")
	}


	for (state in state_auswahl){
		plot(NA,xlim=c(0,360),ylim=c(min(mar[,,state],na.rm=TRUE),max(mar[,,state],na.rm=TRUE)),
			xlab="longitude",ylab="Mann Kendall Trend Test",
			main=paste("meridional average between ",abschnitt[1],"N and ",abschnitt[2],"N \n for markov transition ",transition_names[state],sep=""))
		for (sea in 1:length(seasons)){
			meridional_average(x,mar[sea,,state],colors[sea],1,abschnitt)
		}
		legend("bottomleft",lty=array(1,length(seasons)),col=colors,legend=seasons)
	}

    dur=array(NA,dim=c(length(seasons),ntot,states,8,5))
    for (sea in 1:length(seasons)){
		nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/duration/",period,"/",trendID,"_duration_",states,"s_analysis_",seasons[sea],".nc",sep=""))
		dur[sea,,,,]=get.var.ncdf(nc,"dur_ana_full")
	}


	for (state in 1:states){
		plot(NA,xlim=c(0,360),ylim=c(min(dur[,,state,5,1],na.rm=TRUE),max(dur[,,state,5,1],na.rm=TRUE)),
			xlab="longitude",ylab="95 quantile slope",
			main=paste("meridional average between ",abschnitt[1],"N and ",abschnitt[2],"N \n for 95 quantile regression ",state_names[state],sep=""))
		for (sea in 1:length(seasons)){
			meridional_average(x,dur[sea,,state,5,1],colors[sea],1,abschnitt)
		}
		legend("bottomleft",lty=array(1,length(seasons)),col=colors,legend=seasons)
	}
}

source("load.r")
dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
plot_meridonal_average("91_5",2,"1950-2014",c(60,70))