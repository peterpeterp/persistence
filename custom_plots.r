source("load.r")


plot_multiple_markov <- function(trendID="91_5",states=2,period="1950-2014",var1="mean",var1_sig=NA,
	farb_mitte="mean",farb_palette="regenbogen",name_zusatz="bla",region=NA,seasons=c("spring","summer","autumn","winter","year"),
	grid=FALSE,season_auswahl=c(1,2,3,4),col_row=c(5,2),ausschnitt=c(30,80),cex=cex,subIndex=c("a","b")){

	if (states==2){
		state_names=c("cold","warm")
		trans_names=c("cold to cold","na","na","warm to warm")
		trans_auswahl=c(1,4)
	}
	if (states==3){
		state_names=c("cold","normal","warm")
		trans_auswahl=c(1,5,9)
	}	

	paper=c(12,((col_row[1]-1)*7/5+1))

	reihen=array(NA,dim=c(states*length(season_auswahl),ntot))
	reihen_sig=array(NA,dim=c(states*length(season_auswahl),ntot))
	titel=c()
    for (i in 1:length(season_auswahl)){
    	sea=season_auswahl[i]
    	season=seasons[sea]

		nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/markov/",period,"/",trendID,"_mar",states,"s_trend_",season,".nc",sep=""))
		y=get.var.ncdf(nc,var1)
		if (!is.na(var1_sig)){y_sig=get.var.ncdf(nc,var1_sig)}
		for (trans in 1:length(trans_auswahl)){
			#index=(trans-1)*length(seasons)+sea
			index=(i-1)*states+trans
			reihen[index,]=y[1:ntot,trans_auswahl[trans]]
			if (!is.na(var1_sig)){reihen_sig[index,]=y_sig[1:ntot,trans_auswahl[trans]]}
			titel[index]=paste(state_names[trans],"to",state_names[trans],"-",season)
		}
	}

	map_allgemein(dat=dat,
		filename_plot=paste("../plots/",trendID,"/",states,"_states/maps/",trendID,"_vergleich_",name_zusatz,"_",period,".pdf",sep=""),
		worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette=farb_palette,grid=grid,
		ausschnitt=ausschnitt,col_row=col_row,paper=paper,cex=cex,subIndex=subIndex)
}

plot_seasonal_anomalie <- function(trendID="91_5",states=2,period="1950-2014",var1="mean",var1_sig=NA,
	farb_mitte="mean",farb_palette="regenbogen",name_zusatz="bla",region=NA,seasons=c("spring","summer","autumn","winter"),
	grid=FALSE,season_auswahl=c(1,2,3,4),col_row=c(5,2),ausschnitt=c(30,80),cex=cex,subIndex=c("a","b")){

	if (states==2){
		state_names=c("cold","warm")
		trans_names=c("cold to cold","na","na","warm to warm")
		trans_auswahl=c(1,4)
	}
	if (states==3){
		state_names=c("cold","normal","warm")
		trans_auswahl=c(1,5,9)
	}	

	paper=c(12,((col_row[1]-1)*7/5+1))

	nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/markov/",period,"/",trendID,"_mar",states,"s_trend_year.nc",sep=""))
	ground=get.var.ncdf(nc,var1)

	reihen=array(NA,dim=c(states*length(season_auswahl),ntot))
	reihen_sig=array(NA,dim=c(states*length(season_auswahl),ntot))
	titel=c()
    for (i in 1:length(season_auswahl)){
    	sea=season_auswahl[i]
    	season=seasons[sea]

		nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/markov/",period,"/",trendID,"_mar",states,"s_trend_",season,".nc",sep=""))
		y=get.var.ncdf(nc,var1)
		for (trans in 1:length(trans_auswahl)){
			index=(i-1)*states+trans
			reihen[index,]=y[1:ntot,trans_auswahl[trans]]-ground[1:ntot,trans_auswahl[trans]]
			titel[index]=paste(state_names[trans],"to",state_names[trans],"-",season)
		}
	}

	map_allgemein(dat=dat,
		filename_plot=paste("../plots/",trendID,"/",states,"_states/maps/",trendID,"_vergleich_",name_zusatz,"_",period,".pdf",sep=""),
		worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette=farb_palette,grid=grid,
		ausschnitt=ausschnitt,col_row=col_row,paper=paper,cex=cex,subIndex=subIndex)
}

#init
library(SDMTools)
source("functions_regional.r")
source("map_plot.r")
library(rworldmap)
library(fields)
worldmap = getMap(resolution = "low")
ntot=1319
dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
trendID="91_5"
states=2
yearperiod="1950-2014"

#plot_duration_climatology(trendID,states=states,period="1950-2014",ausschnitt=c(-80,80))

plot_multiple_markov(col_row=c(2,2),cex=0.6,name_zusatz="annual",season_auswahl=c(5),var1="mean",subIndex=c("a","b"))

plot_seasonal_anomalie(farb_mitte=c(-0.08,0.08),farb_palette="lila-gruen",name_zusatz="seas_anom",var1="mean",
	subIndex=c("c","d","e","f","g","h","i","j"),season_auswahl=c(1,2,3,4),col_row=c(5,2),cex=0.6)