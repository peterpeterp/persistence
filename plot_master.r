# teste teste
source("write.r")
source("load.r")

plot_numbWarm <- function(){
	library(SDMTools)
	source("region_average.r")
	source("map_plot.r")
	source("trend_control.r")
	library(rworldmap)
	library(fields)
	worldmap = getMap(resolution = "low")
	nday = 91
	nyr = 5
	ntot=1319
	dat=dat_load("../data/dat_regional.nc",reg=1)
	per=markov_load(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr))

	numbWarm=trend_control_warm_days(dat,per,c(1),c(365),filename=sprintf("../data/warmTage_tremds_year.nc"))
	reihen=array(numbWarm,dim=c(1,ntot))
	map_allgemein(dat=dat,reihen=reihen,titel=c("yearly increase in 'warm days' from 1950 to 2011"),farbe_mitte="0",
		filename_plot=sprintf("../plots/maps/%s_%s_warm_days.pdf",nday,nyr),worldmap=worldmap,ausschnitt=c(35,66))
	map_regional(dat=dat,toPlot=reihen,titles=c("yearly increase in 'warm days' from 1950 to 2011"),worldmap=worldmap,
		filename_plot=sprintf("../plots/regions/%s_%s_warm_days_regio.pdf",nday,nyr))
}

plot_nas <- function(){
	library(SDMTools)
	source("map_plot.r")
	library(rworldmap)
	library(fields)
	worldmap = getMap(resolution = "low")
	nday = 91
	nyr = 5
	ntot=1319
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc",reg=1)
    missings=read.table("../data/number_of_NA_per_station.txt")
    missings2=read.table("../data/number_of_NA_per_station_2011.txt")
    print(missings2[1:ntot,1])
	reihen=array(NA,dim=c(2,ntot))
	reihen[1,]=missings[1:ntot,2]
	reihen[2,]=missings2[1:ntot,2]
	map_allgemein(dat=dat,reihen=reihen,titel=c("missing values to 2014","missing values to 2011"),farbe_mitte="mean",
		filename_plot=sprintf("../plots/%s_%s_missings.pdf",nday,nyr),worldmap=worldmap,ausschnitt=c(-80,80))

}


plot_markov <- function(states,vars,vars_sig,farbe_mitte,name_zusatz,period,region=NA,seasons=c("spring","summer","autumn","winter","year")){
	# markov xx states
	if (states==2){
		state_names=c("cold","warm")
	}
	if (states==3){
		state_names=c("cold","normal","warm")
	}		

    for (season in seasons){
		nc=open.ncdf(paste("../data/91_5/",states,"_states/markov/",period,"/91_5_mar",states,"s_trend_",season,".nc",sep=""))
		reihen=array(NA,dim=c(states*states,ntot))
		reihen_sig=array(NA,dim=c(states*states,ntot))
		titel=c()

		y=array(get.var.ncdf(nc,vars),dim=c(ntot,states,states))
		if (!is.na(vars_sig)){y_sig=array(get.var.ncdf(nc,vars_sig),dim=c(ntot,states,states))}
		for (from in 1:states){
			for (to in 1:states){
				reihen[((from-1)*states+to),]=y[1:ntot,from,to]
				if (!is.na(vars_sig)){reihen_sig[((from-1)*states+to),]=y_sig[1:ntot,from,to]}
				titel[((from-1)*states+to)]=paste(name_zusatz,"for transition from",state_names[from],"to",state_names[to],"in",season,"in",period)
			}
		}
		map_allgemein(dat=dat,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farbe_mitte=farbe_mitte,
			filename_plot=paste("../plots/",states,"_states/maps/markov/",period,"/91_5_mar",states,"s_",vars,"_",season,"_",period,".pdf",sep=""),
			worldmap=worldmap,ausschnitt=c(-80,80),region=region,regionColor="black")		
	}
}

plot_duration_vergleich <- function(states,period,auswahl=c(8,1,2,3,4,5,6),
	titel_zusatz=c("mean","0.25 quantile","0.5 quantile","0.75 quantile","0.9 quantile","0.95 quantile","0.98 quantile"),
	name_zusatz="quaReg"){


	# duration vergleich
    seasons=c("spring","summer","autumn","winter","year")
	if (states==2){
		state_names=c("cold","warm")
	}
	if (states==3){
		state_names=c("cold","normal","warm")
	}		
	
	vars=c("dur_ana_full")

	for (season in seasons){
		nc=open.ncdf(paste("../data/91_5/",states,"_states/duration/",period,"/91_5_duration_",states,"s_analysis_",season,".nc",sep=""))
		titel=c()
		reihen=array(NA,dim=c(length(auswahl)*states,ntot))
		reihen_sig=array(NA,dim=c(length(auswahl)*states,ntot))
		for (i in 1:length(auswahl)){
			for (state in 1:states){
			    reihen[((i-1)*states+state),]=get.var.ncdf(nc,vars[1])[1:ntot,state,auswahl[i],1]
			    reihen_sig[((i-1)*states+state),]=get.var.ncdf(nc,vars[1])[1:ntot,state,auswahl[i],2]
			    titel[((i-1)*states+state)]=paste("trend for",titel_zusatz[i],"of",state_names[state],"period duration in",season,"in",period)
			}
		}
		map_allgemein(dat=dat,
			filename_plot=paste("../plots/",states,"_states/maps/duration/",period,"/91_5_duration_",season,"_",name_zusatz,"_",period,".pdf",sep=""),
			worldmap=worldmap,ausschnitt=c(-80,80),reihen=reihen,reihen_sig=reihen_sig,titel=titel,farbe_mitte="0")
	}
}

plot_duration_climatology <- function(states,period){
	# duration climatology
    seasons=c("spring","summer","autumn","winter","year")
	if (states==2){
		state_names=c("cold","warm")
	}
	if (states==3){
		state_names=c("cold","normal","warm")
	}	
	titel_zusatz=c("0.25","0.5","0.75","0.9","0.95","0.98")
	vars=c("dur_ana_full")
	auswahl=c(1,2,3,4,5,6)

	for (season in seasons){
		nc=open.ncdf(paste("../data/91_5/",states,"_states/duration/",period,"/91_5_duration_",states,"s_analysis_",season,".nc",sep=""))
		titel=c()
		reihen=array(NA,dim=c(length(auswahl)*2,ntot))
		for (i in 1:length(auswahl)){
			for (state in 1:states){
			   	reihen[((i-1)*2+state),]=get.var.ncdf(nc,vars[1])[1:ntot,state,auswahl[i],3]
			    titel[((i-1)*2+state)]=paste("mean value of",titel_zusatz[i],"quantile of",state_names[state],"period duration in",season)
			}
		}			
		map_allgemein(dat=dat,
			filename_plot=paste("../plots/",states,"_states/maps/duration/",period,"/91_5_duration_",season,"_climatology.pdf",sep=""),
		worldmap=worldmap,ausschnitt=c(-80,80),reihen=reihen,titel=titel,farbe_mitte="mean")
	}
}


plot_regional_average <- function(auswahl,titel_zusatz,period,name_zusatz,region_name){
	# regional trend
    seasons=c("spring","summer","autumn","winter")

	state_names=c("cold","warm")
	states=2

    for (sea in 1:length(seasons)){
		nc=open.ncdf(paste("../data/91_5/2_states/regional/",period,"/91_5_",seasons[sea],"_",region_name,".nc",sep=""))
		poli=get.var.ncdf(nc,"region_coordinates")
		reihen=array(NA,dim=c(length(auswahl)*states,dim(poli)[1]))
		reihen_sig=array(NA,dim=c(length(auswahl)*states,dim(poli)[1]))
		titel=c()

		val=get.var.ncdf(nc,"values")
		val_sig=get.var.ncdf(nc,"values_sig")
		for (k in 1:length(auswahl)){
			for (state in 1:2){
				reihen[((k-1)*2+state),]=val[state,auswahl[k],]
				reihen_sig[((k-1)*2+state),]=val_sig[state,auswahl[k],]
				titel[((k-1)*2+state)]=paste(state_names[state],"period duration",titel_zusatz[k],"quantile in",seasons[sea])
			}
		}
		regions_color(reihen=reihen,reihen_sig=reihen_sig,titles=titel,worldmap=worldmap,poli=poli,
			filename_plot=paste("../plots/2_states/regions/",period,"/91_5_",states,"s_",seasons[sea],"_",name_zusatz,"_",region_name,".pdf",sep=""))
	}	
}

plot_meridonal_average <- function(states=2,vars1,vars2,farbe_mitte,period,name_zusatz){
	# markov xx states
	source("meridional_average.r")
	if (states==2){
		state_names=c("cold","warm")
	}
	if (states==3){
		state_names=c("cold","normal","warm")
	}
    seasons=c("spring","summer","autumn","winter","year")		

    for (season in seasons){
		nc=open.ncdf(paste("../data/91_5/2_states/markov/1980-2014/91_5_mar",states,"s_trend_",season,".nc",sep=""))
		climatology=array(NA,dim=c(states*states,ntot))
		change=array(NA,dim=c(states*states,ntot))
		titel=c()

		y1=array(get.var.ncdf(nc,vars1),dim=c(ntot,states,states))
		y2=array(get.var.ncdf(nc,vars2),dim=c(ntot,states,states))
		for (from in 1:states){
			for (to in 1:states){
				climatology[((from-1)*states+to),]=y1[1:ntot,from,to]
				change[((from-1)*states+to),]=y2[1:ntot,from,to]
				titel[((from-1)*states+to)]=paste(vars1,"for transition from",state_names[from],"to",state_names[to],"in",season)
			}
		}
		meridional_average(dat=dat,climatology=climatology,change=change,titel=titel,filename_plot=paste("../plots/2_states/regional/2_states/regional/",period,"/91_5_mar",states,"s_",vars1,"_",season,"_",name_zusatz,".pdf",sep=""))
		
	}
}

#init
library(SDMTools)
source("functions_regional.r")
source("map_plot.r")
library(rworldmap)
library(fields)
worldmap = getMap(resolution = "low")
nday = 91
nyr = 5
ntot=1319
dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")

#wave_region(7,78,c(40,70))
#location_view()




#commands
#plot_markov(vars="mean",vars_sig=NA,farbe_mitte="mean",name_zusatz="climatology",period="1950-2014")
#plot_duration_climatology(period="1950-2014")

#for (yearperiod in c("1950-1980","1950-2014","1980-2014")){
#for (yearperiod in c("1950-2014","1980-2014")){
for (yearperiod in c("1980-2014")){
	#plot_regional_average(auswahl=c(1,2),titel_zusatz=c("MK","LR"),period=yearperiod,name_zusatz="markov",region_name="7wave")
	#plot_regional_average(auswahl=c(6),titel_zusatz=c("95 quantile regression"),period=yearperiod,name_zusatz="95_quantReg")

	#plot_duration_vergleich(period=yearperiod,auswahl=c(5),titel_zusatz=c("95 quantile"),name_zusatz=c("95_quantReg"))
	#plot_duration_vergleich(states=3,period=yearperiod)

	#plot_meridonal_average(vars="mean",vars2="LR",farbe_mitte="mean",name_zusatz="bansfsdfdas")

	plot_markov(states=2,vars="MK",vars_sig="MK_sig",farbe_mitte="0",name_zusatz="MannKendall test",period=yearperiod,
		region="7wave",seasons=c("summer"))

	#plot_markov(states=3,vars="LR",vars_sig="LR_sig",farbe_mitte="0",name_zusatz="Linear regression",period=yearperiod)
}