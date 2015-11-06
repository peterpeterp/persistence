#!/home/pepflei/R/bin/Rscript

# teste teste
source("write.r")
source("load.r")

plot_numbWarm <- function(trendID="91_5",trend_style="_mean",grid=FALSE,ausschnitt=c(-80,80),dataset="_TX",additional_style="",version="_4seasons"){
	library(SDMTools)
	source("map_plot.r")
	source("trend_view.r")
	library(rworldmap)
	library(fields)
	worldmap = getMap(resolution = "low")

	data=read.table(paste("../data/",trendID,"/",dataset,additional_style,"/sonstiges/",trendID,trend_style,dataset,"_5seasons_warme_tage_control",version,".txt",sep=""))

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
		filename_plot=paste("../plots/",trendID,"/sonstiges/numberColdDays/",trendID,dataset,"_cold_days_trend",version,".pdf",sep=""),worldmap=worldmap,ausschnitt=ausschnitt)

	for (sea in 1:5){
		reihen[sea,]=1-data[1:ntot,(10+sea)]
		titel[sea]=paste("percentage of 'cold days' from 1950 to 2014 in",seasons[sea])
	}

	map_allgemein(dat=dat,reihen=reihen,titel=titel,farb_mitte=c(0.45,0.55),
		farb_palette="lila-gruen",
		filename_plot=paste("../plots/",trendID,"/sonstiges/numberColdDays/",trendID,dataset,"_cold_days_percentage",version,".pdf",sep=""),worldmap=worldmap,ausschnitt=ausschnitt)	

	for (sea in 1:5){
		reihen[sea,]=data[1:ntot,(15+sea)]
		titel[sea]=paste("variance of percentage of 'cold days' from 1950 to 2014 in",seasons[sea])
	}

	map_allgemein(dat=dat,reihen=reihen,titel=titel,farb_mitte="mean",
		farb_palette="lila-gruen",
		filename_plot=paste("../plots/",trendID,"/sonstiges/numberColdDays/",trendID,dataset,"_cold_days_percentage_variance",version,".pdf",sep=""),worldmap=worldmap,ausschnitt=ausschnitt)
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

plot_duration_vergleich <- function(period,trendID="91_5",dataset="_TX",additional_style="_seasonal_median",ausschnitt=c(-80,80),auswahl=c(8,1,2,3,4,5,6),titel_zusatz=c("mean","0.25 quantile","0.5 quantile","0.75 quantile","0.9 quantile","0.95 quantile","0.98 quantile"),name_zusatz="quaReg",seasons=c("MAM","JJA","SON","DJF","year","4seasons"),farb_mitte="0",farb_palette="lila-gruen",grid=FALSE){


	state_names=c("cold","warm")
	vars=c("dur_ana_full")

	for (season in seasons){
        nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",period,"/",trendID,dataset,"_duration_analysis_",season,".nc",sep=""))
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
		map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/duration/",period,"/",season,"/","duration_trend_",trendID,"_",season,"_",name_zusatz,"_",period,additional_style,".pdf",sep=""),worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette=farb_palette,grid=grid,ausschnitt=ausschnitt)
	}
}

plot_duration_climatology <- function(period,trendID="91_5",dataset="_TX",additional_style="",ausschnitt=c(-80,80),seasons=c("MAM","JJA","SON","DJF","year"),farb_mitte="mean",farb_palette="regenbogen",grid=FALSE,auswahl=c(1,2,3,4,5,6,8),titel_zusatz=c("0.25","0.5","0.75","0.9","0.95","0.98","mean"),name_zusatz="",region=NA,color_lab="# days"){
	# duration climatology
    
	state_names=c("cold","warm")

	
	vars=c("dur_ana_full")
	

	for (season in seasons){
        nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",period,"/",trendID,dataset,"_duration_analysis_",season,".nc",sep=""))
		titel=c()
		reihen=array(NA,dim=c(length(auswahl)*states,ntot))
		reihen_sig=array(NA,dim=c(length(auswahl)*states,ntot))
		for (i in 1:length(auswahl)){
			for (state in 1:states){
			   	reihen[((i-1)*states+state),]=get.var.ncdf(nc,vars[1])[1:ntot,state,auswahl[i],3]
			    titel[((i-1)*states+state)]=paste("mean value of",titel_zusatz[i],"quantile of",state_names[state],"period duration in",season)
			}
		}	
		if (titel_zusatz[1]==""){titel=c("")}
		map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/duration/",period,"/",season,"/","duration_climatology_",trendID,"_",season,"_",additional_style,"_",name_zusatz,".pdf",sep=""),worldmap=worldmap,ausschnitt=ausschnitt,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette=farb_palette,grid=grid,region=region,color_lab=color_lab)

		if (length(auswahl)==1){
			reihen2=array(NA,dim=c(2,ntot))
			reihen_sig2=array(NA,dim=c(2,ntot))
			titel=c()
			
			reihen2[1,]=(reihen[1,]+reihen[2,])/2
			titel[1]=paste("mean(cold,warm) duration length",titel_zusatz[i],"quantile in",season)
			reihen2[2,]=reihen[2,]-reihen[1,]
			titel[2]=paste("difference between warm and cold duration lengths of",titel_zusatz[i],"quantile in",season)

			if (titel_zusatz[1]==""){titel=c("")}
			map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/duration/",period,"/",season,"/","duration_climatology_waco_anom_",trendID,"_",season,"_",additional_style,"_",name_zusatz,".pdf",sep=""),worldmap=worldmap,ausschnitt=ausschnitt,reihen=reihen2,reihen_sig=reihen_sig2,titel=titel,farb_mitte=c("individual","mean","0"),farb_palette=c("individual","regenbogen","lila-gruen"),grid=grid,region=region,color_lab=color_lab)
		}
	}
}

plot_duration_distribution <- function(period="1950-2015",trendID="91_5",dataset="_TX",additional_style="_seasonal_median",ausschnitt=c(-80,80),auswahl=c(7,8),titel_zusatz=c("lifetime 1/b","r squared"),name_zusatz="b_R2",seasons=c("MAM","JJA","SON","DJF","year","4seasons"),farb_mitte="mean",farb_palette="lila-gruen",grid=FALSE,region=NA,average=FALSE){
	# duration climatology
	state_names=c("cold","warm")

	vars=c("distr_ana")

	for (season in seasons){
		nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",period,"/",trendID,dataset,"_duration_distribution_ana_",season,".nc",sep=""))
		titel=c()
		reihen=array(NA,dim=c(length(auswahl)*states,ntot))
		reihen_sig=array(NA,dim=c(length(auswahl)*states,ntot))
		for (i in 1:length(auswahl)){
			toplo=auswahl[i]
			for (state in 1:states){
			   	reihen[((i-1)*states+state),]=get.var.ncdf(nc,vars[1])[1:ntot,state,toplo]
			    titel[((i-1)*states+state)]=paste(titel_zusatz[i],"of",state_names[state],"period duration in",season)
			}
		}	
		print(mean(reihen[1,],na.rm=TRUE))	
		print(mean(reihen[2,],na.rm=TRUE))	
		print(sd(reihen[1,],na.rm=TRUE))	
		print(sd(reihen[2,],na.rm=TRUE))	
		map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/duration/",period,"/",season,"/duration_distr_",trendID,"_",season,"_",additional_style,"_",name_zusatz,".pdf",sep=""),worldmap=worldmap,ausschnitt=ausschnitt,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette="regenbogen",region=region,grid=grid,average=average)
	}
}

plot_eke <- function(vars,vars_sig,farb_mitte,name_zusatz,period,ausschnitt=c(-80,80),region=NA,
	seasons=c("MAM","JJA","SON","DJF","year"),grid=FALSE){
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
			filename_plot=paste("../plots/eke/eke_",vars,"_",season,"_",period,".pdf",sep=""),
			worldmap=worldmap,ausschnitt,region=region,regionColor="black",grid=grid)		
	}
}

full_plot <- function(trendID="91_5",dataset="_TX",additional_style="_seasonal_median",states=2,ausschnitt=c(-80,80)){
	#plot_duration_climatology(trendID,states=states,period="1950-2014",ausschnitt=ausschnitt,dataset=dataset,additional_style=additional_style)
	titel_zusatz=c("25","50","75","90","95","98",NA,"mean")
	for (qua in c(3,4,5,6,8)){
		print(qua)
		plot_duration_climatology(trendID,period="1950-2014",ausschnitt=ausschnitt,dataset=dataset,additional_style="",name_zusatz=titel_zusatz[qua],auswahl=c(qua),titel_zusatz=c(titel_zusatz[qua]),farb_mitte="gemeinsam mean")
		plot_duration_climatology(trendID,period="1950-2014",ausschnitt=ausschnitt,dataset=dataset,additional_style="",name_zusatz=titel_zusatz[qua],auswahl=c(qua),titel_zusatz=c(titel_zusatz[qua]),farb_mitte="gemeinsam mean",seasons=c("4seasons"))
	}
	for (yearperiod in c("1950-2014","1980-2014","1950-1980")){
		print(yearperiod)
		for (qua in c(3,4,5,6,8)){
			plot_duration_vergleich(trendID,period=yearperiod,ausschnitt=ausschnitt,dataset=dataset,additional_style=additional_style,auswahl=c(qua),farb_mitte="gemeinsam 0",name_zusatz=titel_zusatz[qua],titel_zusatz=titel_zusatz[qua])			
			plot_duration_vergleich(trendID,period=yearperiod,ausschnitt=ausschnitt,dataset=dataset,additional_style=additional_style,auswahl=c(qua),farb_mitte="gemeinsam 0",name_zusatz=titel_zusatz[qua],titel_zusatz=titel_zusatz[qua],seasons=c("4seasons"))
		}
	}
}



#init
library(SDMTools)
source("functions_regional.r")
source("map_plot.r")
library(rworldmap)
library(fields)
worldmap = getMap(resolution = "low")

ntot=1319

states=2
trendID="91_5"
dataset="_TX"
additional_style=""

dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))


#plot_duration_climatology(trendID,states=states,period="1950-2014",dataset=dataset,additional_style="",name_zusatz="mean_7rect",auswahl=c(8),titel_zusatz=c(""),farb_mitte="gemeinsam mean",seasons=c("4seasons"),region="7rect",grid=TRUE)

titel=c("mean","relative sd","skewness","1/b","R2","A")
name=c("mean","rel_sd","skew","t","R2","A")
farb=c(3,7.5,0.8,1.2,2,4,3,12,0.7,1,0.02,0.2)
auswahl=c(1,3,4,7,8,5)
#for (i in c(1)){
for (i in c(1,2,3,4,5,6)){
	for (period in c("1950-2014","1950-1980","1980-2014")){
		if (is.na(farb[(i-1)*2+2])){farb_mitte="mean"}
		else {farb_mitte=farb[((i-1)*2+1):((i-1)*2+2)]}
		plot_duration_distribution(trendID=trendID,period=period,dataset=dataset,additional_style="",auswahl=c(auswahl[i]),farb_mitte=farb_mitte,titel_zusatz=c(titel[i]),name_zusatz=name[i],seasons=c("4seasons"),average=TRUE,ausschnitt=c(35,60),region="7rect")
	}
}



#full_plot(trendID=trendID,states=states,dataset=dataset,additional_style=additional_style)

#commands
if (1==2){
	trendID="91_5"
	yearPeriod=c(1950,2014)
	region_name="7rect"
	plot_regional_boxplots(trendID,dat,yearPeriod,region_name)
	#plot_regional_distributions(trendID,dat,yearPeriod,region_name)
}


#plot_numbWarm(trendID="91_5",dataset="_TX",additional_style=additional_style,version="_4seasons")
#plot_numbWarm(trendID="91_5",dataset="_TX",additional_style=additional_style,version="_year")


#wave_region(7,78,c(40,70))
#location_view()




#plot_correl_duration("91_5",states=2,vars="correlation",vars_sig=NA,farb_mitte="0")
#plot_eke(vars="MK",vars_sig="MK_sig",farb_mitte="0",name_zusatz="MannKendall test",period="1979-2014",ausschnitt=c(-80,80))
#plot_eke(vars="mean",vars_sig=NA,farb_mitte="mean",name_zusatz="climatology",period="1979-2014",ausschnitt=c(-80,80))

#plot_duration_distribution("91_5",2,"1950-2014")





#plot_duration_climatology(trendID,states=states,period="1950-2014",farb_mitte=c(-10,30),seasons=c("summer"))
#plot_duration_vergleich("91_5",states=2,period="1980-2014",farb_mitte=c(0.6,0.6),seasons=c("summer"),name_zusatz="quaReg_colore")
#plot_markov(trendID,states=states,vars="MK",vars_sig="MK_sig",farb_mitte="0",titel_zusatz="MannKendall test",name_zusatz="MK_grid",period=yearperiod,grid=TRUE)






#plot_detrended_seasonal_skewness(dat,grid=FALSE,ausschnitt=c(-80,80),trendID="91_5",trend_style="_mean",dataset="_TX",additional_style="_median")
