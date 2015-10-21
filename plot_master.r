#!/home/pepflei/R/bin/Rscript

# teste teste
source("write.r")
source("load.r")

plot_numbWarm <- function(trendID="91_5",grid=FALSE,ausschnitt=c(-80,80),trend_style="_mean",dataset="_TX",additional_style="_median"){
	library(SDMTools)
	source("map_plot.r")
	source("trend_view.r")
	library(rworldmap)
	library(fields)
	worldmap = getMap(resolution = "low")

	data=read.table(paste("../data/",trendID,"/",additional_style,"/sonstiges/",trendID,trend_style,dataset,"_5seasons_warme_tage_control",additional_style,".txt",sep=""))

	seasons=c("spring","summer","autumn","winter","year")
	reihen=array(NA,dim=c(5,ntot))
	reihen_sig=array(NA,dim=c(5,ntot))
	titel=c()
	for (sea in 1:5){
		reihen[sea,]=data[1:ntot,sea]
		reihen_sig[sea,]=data[1:ntot,(5+sea)]
		titel[sea]=paste("yearly increase in 'cold days' from 1950 to 2014 in",seasons[sea])
	}

	map_allgemein(dat=dat,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte="0",
		farb_palette="lila-gruen",
		filename_plot=paste("../plots/",trendID,"/sonstiges/",trendID,trend_style,dataset,"_cold_days_trend",additional_style,".pdf",sep=""),worldmap=worldmap,ausschnitt=ausschnitt)

	for (sea in 1:5){
		reihen[sea,]=1-data[1:ntot,(10+sea)]
		titel[sea]=paste("percentage of 'cold days' from 1950 to 2014 in",seasons[sea])
	}

	map_allgemein(dat=dat,reihen=reihen,titel=titel,farb_mitte=c(0.45,0.55),
		farb_palette="lila-gruen",
		filename_plot=paste("../plots/",trendID,"/sonstiges/",trendID,trend_style,dataset,"_cold_days_percentage",additional_style,".pdf",sep=""),worldmap=worldmap,ausschnitt=ausschnitt)	

	for (sea in 1:5){
		reihen[sea,]=1-data[1:ntot,(15+sea)]
		titel[sea]=paste("variance of percentage of 'cold days' from 1950 to 2014 in",seasons[sea])
	}

	map_allgemein(dat=dat,reihen=reihen,titel=titel,farb_mitte="mean",
		farb_palette="lila-gruen",
		filename_plot=paste("../plots/",trendID,"/sonstiges/",trendID,trend_style,dataset,"_cold_days_percentage_variance",additional_style,".pdf",sep=""),worldmap=worldmap,ausschnitt=ausschnitt)
}

plot_skewness <- function(dat,grid=FALSE,ausschnitt=c(-80,80)){
	library(moments)

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
		filename_plot=paste("../plots/zwischenzeugs/skewness.pdf",sep=""),worldmap=worldmap,ausschnitt=ausschnitt)
}

plot_nas <- function(grid=FALSE,trend_style=""){
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
	map_allgemein(dat=dat,reihen=reihen,titel=c("missing values to 2014","missing values to 2011"),farb_mitte="mean",
		filename_plot=sprintf("../plots/%s_%s_missings",trend_style,".pdf",nday,nyr),worldmap=worldmap,ausschnitt=c(-80,80),grid=grid)

}


plot_markov <- function(trendID,states,vars,vars_sig,farb_mitte,farb_palette="lila-gruen",
	name_zusatz=vars,titel_zusatz,period,ausschnitt=c(-80,80),region=NA,
	seasons=c("spring","summer","autumn","winter","year"),grid=FALSE,trend_style=""){
	# markov xx states
	if (states==2){
		state_names=c("cold","warm")
	}
	if (states==3){
		state_names=c("cold","normal","warm")
	}		

    for (season in seasons){
		nc=open.ncdf(paste("../data/",trendID,"/",states,"_states",trend_style,"/markov/",period,"/",trendID,"_mar",states,"s_trend_",season,".nc",sep=""))
		reihen=array(NA,dim=c(states*states,ntot))
		reihen_sig=array(NA,dim=c(states*states,ntot))
		titel=c()

		y=array(get.var.ncdf(nc,vars),dim=c(ntot,states,states))
		if (!is.na(vars_sig)){y_sig=array(get.var.ncdf(nc,vars_sig),dim=c(ntot,states,states))}
		for (from in 1:states){
			for (to in 1:states){
				reihen[((from-1)*states+to),]=y[1:ntot,from,to]
				if (!is.na(vars_sig)){reihen_sig[((from-1)*states+to),]=y_sig[1:ntot,from,to]}
				titel[((from-1)*states+to)]=paste(titel_zusatz,"for transition from",state_names[from],"to",state_names[to],"in",season,"in",period)
			}
		}
		map_allgemein(dat=dat,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette=farb_palette,
			filename_plot=paste("../plots/",trendID,"/",states,"_states",trend_style,"/maps/markov/",period,"/",trendID,"_mar",states,"s_",name_zusatz,"_",season,"_",period,".pdf",sep=""),
			worldmap=worldmap,ausschnitt=ausschnitt,region=region,regionColor="black",grid=grid)		
	}
}


plot_duration_vergleich <- function(trendID,states,period,trend_style="",ausschnitt=c(-80,80),auswahl=c(8,1,2,3,4,5,6),
	titel_zusatz=c("mean","0.25 quantile","0.5 quantile","0.75 quantile","0.9 quantile","0.95 quantile","0.98 quantile"),
	name_zusatz="quaReg",seasons=c("spring","summer","autumn","winter","year"),farb_mitte="0",farb_palette="lila-gruen",grid=FALSE){

	if (states==2){
		state_names=c("cold","warm")
	}
	if (states==3){
		state_names=c("cold","normal","warm")
	}		
	
	vars=c("dur_ana_full")

	for (season in seasons){
		nc=open.ncdf(paste("../data/",trendID,"/",states,"_states",trend_style,"/duration/",period,"/",trendID,"_duration_",states,"s_analysis_",season,".nc",sep=""))
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
			filename_plot=paste("../plots/",trendID,"/",states,"_states",trend_style,"/maps/duration/",period,"/",trendID,"_duration_",season,"_",name_zusatz,"_",period,".pdf",sep=""),
			worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette=farb_palette,grid=grid,ausschnitt=ausschnitt)
	}
}

plot_duration_climatology <- function(trendID,states,period,trend_style="",ausschnitt=c(-80,80),
	seasons=c("spring","summer","autumn","winter","year"),farb_mitte="mean",farb_palette="regenbogen",grid=FALSE){
	# duration climatology
    
	if (states==2){
		state_names=c("cold","warm")
	}
	if (states==3){
		state_names=c("cold","normal","warm")
	}	
	titel_zusatz=c("0.25","0.5","0.75","0.9","0.95","0.98","mean")
	vars=c("dur_ana_full")
	auswahl=c(1,2,3,4,5,6,8)

	for (season in seasons){
		nc=open.ncdf(paste("../data/",trendID,"/",states,"_states",trend_style,"/duration/",period,"/",trendID,"_duration_",states,"s_analysis_",season,".nc",sep=""))
		titel=c()
		reihen=array(NA,dim=c(length(auswahl)*states,ntot))
		reihen_sig=array(NA,dim=c(length(auswahl)*states,ntot))
		for (i in 1:length(auswahl)){
			for (state in 1:states){
			   	reihen[((i-1)*states+state),]=get.var.ncdf(nc,vars[1])[1:ntot,state,auswahl[i],3]
			    titel[((i-1)*states+state)]=paste("mean value of",titel_zusatz[i],"quantile of",state_names[state],"period duration in",season)
			}
		}			
		map_allgemein(dat=dat,
			filename_plot=paste("../plots/",trendID,"/",states,"_states",trend_style,"/maps/duration/",period,"/",trendID,"_duration_",season,"_climatology.pdf",sep=""),
		worldmap=worldmap,ausschnitt=ausschnitt,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette=farb_palette,grid=grid)
	}
}

plot_duration_distribution <- function(trendID,states,period,ausschnitt=c(-80,80),grid=FALSE,trend_style=""){
	# duration climatology
    seasons=c("spring","summer","autumn","winter","year")
	if (states==2){
		state_names=c("cold","warm")
	}
	if (states==3){
		state_names=c("cold","normal","warm")
	}	
	titel_zusatz=c("lifetime 1/b","chi squared of exp fit")
	vars=c("distr_ana")
	auswahl=c(3,4)

	for (season in seasons){
		nc=open.ncdf(paste("../data/",trendID,"/",states,"_states",trend_style,"/duration/",period,"/",trendID,"_duration_",states,"s_distribution_",season,".nc",sep=""))
		titel=c()
		reihen=array(NA,dim=c(length(auswahl)*states,ntot))
		for (i in 1:length(auswahl)){
			for (state in 1:states){
				print(dim(get.var.ncdf(nc,vars[1])))
			   	reihen[((i-1)*states+state),]=get.var.ncdf(nc,vars[1])[1:ntot,state,auswahl[i]]
			    titel[((i-1)*states+state)]=paste(titel_zusatz[i],"of",state_names[state],"period duration in",season)
			}
		}			
		map_allgemein(dat=dat,
			filename_plot=paste("../plots/",trendID,"/",states,"_states",trend_style,"/maps/duration/",period,"/",trendID,"_duration_",season,"_distr.pdf",sep=""),
		worldmap=worldmap,ausschnitt,reihen=reihen,titel=titel,farb_mitte="mean",farb_palette="regenbogen",grid=grid)
	}
}

plot_regional_average <- function(auswahl,titel_zusatz,period,name_zusatz,region_name,grid=FALSE,trend_style=""){
	# regional trend
    seasons=c("spring","summer","autumn","winter")

	state_names=c("cold","warm")
	states=2

    for (sea in 1:length(seasons)){
		nc=open.ncdf(paste("../data/91_5/2_states",trend_style,"/regional/",period,"/91_5_",seasons[sea],"_",region_name,".nc",sep=""))
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
			filename_plot=paste("../plots/2_states",trend_style,"/regions/",period,"/91_5_",states,"s_",seasons[sea],"_",name_zusatz,"_",region_name,".pdf",sep=""),grid=grid)
	}	
}

plot_eke <- function(vars,vars_sig,farb_mitte,name_zusatz,period,ausschnitt=c(-80,80),region=NA,
	seasons=c("spring","summer","autumn","winter","year"),grid=FALSE,trend_style=""){
	# eke analysis

    for (season in seasons){
		nc=open.ncdf(paste("../data/eke/eke_analysis_",season,".nc",sep=""))
		pressure=get.var.ncdf(nc,"levelist")
		reihen=array(NA,dim=c(6,ntot))
		reihen_sig=array(NA,dim=c(6,ntot))
		titel=c()

		y=get.var.ncdf(nc,vars)
		if (!is.na(vars_sig)){y_sig=get.var.ncdf(nc,vars_sig)}
		for (lvl in 1:6){
			reihen[lvl,]=y[1:ntot,lvl]
			if (!is.na(vars_sig)){reihen_sig[lvl,]=y_sig[1:ntot,lvl]}
			titel[lvl]=paste("Eddy kinetic Energy",name_zusatz,"in",season,"in",period,"at",pressure[lvl],"mbar")
		}
		map_allgemein(dat=dat,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette="regenbogen",
			filename_plot=paste("../plots/eke/eke_",vars,"_",season,"_",period,trend_style,".pdf",sep=""),
			worldmap=worldmap,ausschnitt,region=region,regionColor="black",grid=grid)		
	}
}

full_plot <- function(trendID,states,ausschnitt=c(-80,80),trend_style=""){
	plot_markov(trendID,states=states,vars="mean",vars_sig=NA,farb_mitte="mean",farb_palette="regenbogen",titel_zusatz="climatology",period="1950-2014",ausschnitt=ausschnitt,trend_style=trend_style)
	plot_duration_climatology(trendID,states=states,period="1950-2014",ausschnitt=ausschnitt,trend_style=trend_style)
	for (yearperiod in c("1950-1980","1950-2014","1980-2014")){
		print(yearperiod)
		plot_markov(trendID,states=states,vars="MK",vars_sig="MK_sig",farb_mitte="0",farb_palette="lila-gruen",titel_zusatz="MannKendall test",period=yearperiod,ausschnitt=ausschnitt,trend_style=trend_style)
		plot_markov(trendID,states=states,vars="LR",vars_sig="LR_sig",farb_mitte="0",farb_palette="lila-gruen",titel_zusatz="Linear Regression",period=yearperiod,ausschnitt=ausschnitt,trend_style=trend_style)
		plot_duration_vergleich(trendID,states=states,period=yearperiod,ausschnitt=ausschnitt,trend_style=trend_style)
	}
}



#init
library(SDMTools)
source("functions_regional.r")
source("map_plot.r")
library(rworldmap)
library(fields)
worldmap = getMap(resolution = "low")
#trendID="91_5"
ntot=1319
dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")

#plot_skewness(dat=dat)

#commands
if (1==2){
	trendID="91_5"
	yearPeriod=c(1950,2014)
	region_name="7rect"
	plot_regional_boxplots(trendID,dat,yearPeriod,region_name)
	#plot_regional_distributions(trendID,dat,yearPeriod,region_name)

}


#wave_region(7,78,c(40,70))
#location_view()



#plot_correl_markov("91_5",states=2,vars="correlation",vars_sig=NA,farb_mitte="0")
#plot_correl_markov("91_3",states=2,vars="correlation",vars_sig=NA,farb_mitte="0")

#plot_correl_duration("91_5",states=2,vars="correlation",vars_sig=NA,farb_mitte="0")
#plot_eke(vars="MK",vars_sig="MK_sig",farb_mitte="0",name_zusatz="MannKendall test",period="1979-2014",ausschnitt=c(-80,80))
#plot_eke(vars="mean",vars_sig=NA,farb_mitte="mean",name_zusatz="climatology",period="1979-2014",ausschnitt=c(-80,80))

#plot_duration_distribution("91_5",2,"1950-2014")





#plot_duration_climatology(trendID,states=states,period="1950-2014",farb_mitte=c(-10,30),seasons=c("summer"))
#plot_duration_vergleich("91_5",states=2,period="1980-2014",farb_mitte=c(0.6,0.6),seasons=c("summer"),name_zusatz="quaReg_colore")
#plot_markov(trendID,states=states,vars="MK",vars_sig="MK_sig",farb_mitte="0",titel_zusatz="MannKendall test",name_zusatz="MK_grid",period=yearperiod,grid=TRUE)


#full_plot("91_5",2,trend_style="_mean_TX")
#full_plot("91_5",3)
#full_plot("91_3",2)
#full_plot("91_3",3)

#plot_numbWarm(trendID="91_5",states=2,trend_style="_median_TX")
#plot_numbWarm(trendID="91_5",states=2,trend_style="_mode_TX")
#plot_numbWarm(trendID="91_5",states=2,trend_style="_mean_TX")
plot_numbWarm(trendID="91_5",trend_style="_mean",dataset="_TX",additional_style="_median")
