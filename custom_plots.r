source("load.r")

plot_temp_anom_trend <-function(trendID="91_5",dataset="_TMean",additional_style="",yearPeriod=c(1950,2014),var1="dur_ana_full",var1_sig=NA,var_qua=5,farb_mitte="mean",farb_palette="regenbogen",name_zusatz="bla",titel_zusatz="sd",region=NA,seasons=c("MAM","JJA","SON","DJF"),paper=c(12,8),grid=FALSE,season_auswahl=c(1,2,3,4),baseline="4seasons",col_row=c(5,2),ausschnitt=c(30,80),cex=cex,color_lab="mean duration length anomaly to annual mean",cex_axis=1){

	period=paste(yearPeriod[1],"-",yearPeriod[2],sep="")
	if (1==2){
		inside=which(dat$time>yearPeriod[1] & dat$time<=yearPeriod[2])

		reihen=array(NA,dim=c(1,ntot))
		reihen_sig=reihen*NA
		for (q in 1:ntot){
			cat("-")
			if (length(as.vector(dat$tas[q,,])[inside])>1000){
				lr=summary(lm(as.vector(dat$tas[q,,])[inside]~dat$time[inside],na.rm=TRUE))
				#print(lr)
				reihen[1,q]=lr$coefficients[2]
				reihen_sig[1,q]=lr$coefficients[8]
			}
		}
		table=array(NA,dim=c(2,ntot))
		table[1,]=reihen[1,]
		table[2,]=reihen_sig[1,]
		write.table(table,paste("../data/sonstiges/temp_anomaly_trend_",period,dataset,".txt",sep=""))
	}
	if (1==1){
		table=read.table(paste("../data/sonstiges/temp_anomaly_trend_",period,dataset,".txt",sep=""))
		reihen=array(NA,dim=c(1,ntot))
		reihen_sig=reihen*NA
		for (q in 1:ntot){
			reihen[1,q]=table[1,q]
			#reihen_sig[1,q]=table[2,q]
		}
	}

	filename_plot=paste("../plots/zwischenzeugs/temp_trend_",period,dataset,".pdf",sep="")
	map_allgemein(dat=dat,filename_plot=filename_plot,worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=c(""),farb_mitte="0",farb_palette="blau-rot",grid=grid,ausschnitt=c(30,80),col_row=c(1,1),cex_axis=1,cex=1,paper=c(8,3))
}

plot_warm_cold_diff_duration <- function(trendID="91_5",dataset="_TX",additional_style="_seasonal_median",states=2,period="1950-2014",var1="dur_ana_full",var1_sig=NA,	farb_mitte="mean",farb_palette="regenbogen",titel_zusatz="95_quan_diff"	,region=NA,seasons=c("spring","summer","autumn","winter","year"),grid=FALSE,season_auswahl=c(5),col_row=c(5,2),ausschnitt=c(30,80),cex=cex,subIndex=c("a","b")){

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



	reihen=array(NA,dim=c(2*length(season_auswahl),ntot))
	reihen_sig=array(NA,dim=c(2*length(season_auswahl),ntot))
	titel=c()
    for (i in 1:length(season_auswahl)){
    	sea=season_auswahl[i]
    	season=seasons[sea]
    	print(season)

	nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",period,"/",trendID,dataset,"_duration_analysis_",season,".nc",sep=""))
		y=get.var.ncdf(nc,var1)
		index=(i-1)*2+1
		reihen[index,]=y[1:ntot,1,5,3]-y[1:ntot,2,5,3]
		titel[index]=paste("cold - warm for 95. quantile",season)
		index=(i-1)*2+2
		reihen[index,]=y[1:ntot,1,8,3]-y[1:ntot,2,8,3]
		titel[index]=paste("cold - warm for mean duration length",season)
	}

	filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/",trendID,dataset,"_cold_warm_vergleich_",additional_style,"_",period,".pdf",sep="")
	map_allgemein(dat=dat,		filename_plot=filename_plot,worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette=farb_palette,grid=grid,ausschnitt=ausschnitt,col_row=col_row,cex=cex,subIndex=subIndex)
}

plot_diff_trend_duration <- function(trendID="91_5",states=2,period="1950-2014",var="dur_ana_full",version1="_mean_TX",version2="_mean_TN",	farb_mitte="0",farb_palette="lila-gruen",name_zusatz="mean_dur_len",titel_zusatz="95_quan_diff"
	,region=NA,seasons=c("spring","summer","autumn","winter","year"),verglichenee_ergebnis=3,
	grid=FALSE,season_auswahl=c(1,2,3,4,5),quantiles=c(1,2,3,4,5,6,8),col_row=c(5,2),ausschnitt=c(30,80),cex=1,subIndex=c("a","b")){

	if (states==2){
		state_names=c("cold","warm")
		trans_names=c("cold to cold","na","na","warm to warm")
		trans_auswahl=c(1,4)
	}


	paper=c(12,((col_row[1]-1)*7/5+1))

	
	quantile_names=c("0.25 quantile","0.5 quantile","0.75 quantile","0.9 quantile","0.95 quantile","0.98 quantile","na","mean duration length")


	reihen=array(NA,dim=c(2*length(season_auswahl)*length(quantiles),ntot))
	reihen_sig=array(NA,dim=c(2*length(season_auswahl)*length(quantiles),ntot))
	titel=c()
	index=0
    for (i in 1:length(season_auswahl)){
    	sea=season_auswahl[i]
    	season=seasons[sea]
    	print(season)

		nc1=open.ncdf(paste("../data/",trendID,"/",version1,"/duration/",period,"/",trendID,"_duration_",states,"s_analysis_",season,".nc",sep=""))
		nc2=open.ncdf(paste("../data/",trendID,"/",version2,"/duration/",period,"/",trendID,"_duration_",states,"s_analysis_",season,".nc",sep=""))
		y1=get.var.ncdf(nc1,var)
		y2=get.var.ncdf(nc2,var)

		for (quAn in quantiles){
			print(quAn)
			print(quantile_names[quAn])
			for (state in 1:2){
				print(state_names[state])
				index=index+1
				reihen[index,]=y1[1:ntot,state,quAn,verglichenee_ergebnis]-y2[1:ntot,state,quAn,verglichenee_ergebnis]
				titel[index]=paste(version1,"-",version2,quantile_names[quAn],name_zusatz,state_names[state],season)
			}
		}
	}

	filename_plot=paste("../plots/",trendID,"/sonstiges/trend_vergleich/",trendID,"_vergleich_",version1,"_",version2,"_",name_zusatz,".pdf",sep="")

	map_allgemein(dat=dat,
		filename_plot=filename_plot,
		worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette=farb_palette,grid=grid,
		ausschnitt=ausschnitt,col_row=col_row,cex=cex,subIndex=subIndex)
}



plot_seasonal_anomalie_duration <- function(trendID="91_5",dataset="_TX",additional_style="",period="1950-2014",var1="dur_ana_full",var1_sig=NA,var_qua=5,farb_mitte="mean",farb_palette="regenbogen",name_zusatz="bla",titel_zusatz="sd",region=NA,seasons=c("MAM","JJA","SON","DJF"),paper=c(12,8),grid=FALSE,season_auswahl=c(1,2,3,4),baseline="4seasons",col_row=c(5,2),ausschnitt=c(30,80),cex=cex,color_lab="mean duration length anomaly to annual mean",cex_axis=1){

	state_names=c("cold","warm")

	if (col_row[1]>1){
		paper=c(12,((col_row[1]-1)*7/5+1))
		cex_axis=1.5
	}

	nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",period,"/",trendID,dataset,"_duration_analysis_",baseline,".nc",sep=""))
	ground=get.var.ncdf(nc,var1)

	reihen=array(NA,dim=c(states*length(season_auswahl),ntot))
	reihen_sig=array(NA,dim=c(states*length(season_auswahl),ntot))
	titel=c()
	titel[1]=c("")
    for (i in 1:length(season_auswahl)){
    	sea=season_auswahl[i]
    	season=seasons[sea]

		nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",period,"/",trendID,dataset,"_duration_analysis_",season,".nc",sep=""))
		y=get.var.ncdf(nc,var1)
		for (state in 1:2){
			index=(i-1)*states+state
			reihen[index,]=y[1:ntot,state,var_qua,3]-ground[1:ntot,state,var_qua,3]
			if (col_row[1]==1){titel[index]=paste(season,titel_zusatz,state_names[state],"period anomalie to annual mean")}
		}
	}
	filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/",trendID,"_vergleich_",name_zusatz,"_",period,"_",titel_zusatz,"_kompakt.pdf",sep="")
	#map_allgemein(dat=dat,filename_plot=filename_plot,worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette=farb_palette,grid=grid,ausschnitt=ausschnitt,col_row=col_row,paper=paper,cex=cex,color_lab=color_lab,cex_axis=cex_axis)

	# wa co anomaly
	#nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",period,"/",trendID,dataset,"_duration_analysis_",baseline,".nc",sep=""))
	#ground=get.var.ncdf(nc,var1)

	reihen=array(NA,dim=c(9,ntot))
	reihen_sig=array(NA,dim=c(9,ntot))
	titel=c()
	farb_mitte=c()
	farb_palette=c()
	farb_mitte[1]="individual"
	farb_palette[1]="individual"
    for (i in 1:length(season_auswahl)){
    	sea=season_auswahl[i]
    	season=seasons[sea]

		nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",period,"/",trendID,dataset,"_duration_analysis_",season,".nc",sep=""))
		y=get.var.ncdf(nc,var1)
		

		index=(i-1)*states+1
		zwi=(y[,1,var_qua,3]+y[,2,var_qua,3])/2-(ground[,1,var_qua,3]+ground[,2,var_qua,3])/2
		reihen[index,]=zwi
		titel[index]=paste(season,titel_zusatz,"period anomalie to annual mean")
		farb_mitte[index+1]="0"
		farb_palette[index+1]="lila-gruen"
		index=(i-1)*states+2
		zwi=(y[,2,var_qua,3] - ground[,2,var_qua,3])-(y[,1,var_qua,3] - ground[,1,var_qua,3])
		reihen[index,]=zwi	
		titel[index]=paste(season,titel_zusatz,"period anomalie to annual mean difference warm - cold")
		farb_mitte[index+1]="0"
		farb_palette[index+1]="lila-gruen"
	}
	filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/",trendID,"_vergleich_",name_zusatz,"_",period,"_",titel_zusatz,"_kompakt_waco.pdf",sep="")
	#map_allgemein(dat=dat,filename_plot=filename_plot,worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette=farb_palette,grid=grid,ausschnitt=ausschnitt,col_row=col_row,paper=paper,cex=cex,cex_axis=cex_axis)
	
	for (q in 1:ntot){
		reihen[9,q]=sd(reihen[c(1,3,5,7),q])
	}
	reihen2=array(reihen[9,],dim=c(1,ntot))

	filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/",trendID,"_vergleich_",name_zusatz,"_",period,"_",titel_zusatz,"_kompakt_mean(waco).pdf",sep="")
	map_allgemein(dat=dat,filename_plot=filename_plot,worldmap=worldmap,reihen=reihen[c(1,3,5,7),],reihen_sig=reihen_sig[c(1,3,5,7),],titel=c(""),farb_mitte=c(-1.6,1.6),farb_palette=farb_palette,grid=grid,ausschnitt=ausschnitt,col_row=c(4,1),paper=c(8,7),cex=0.8,cex_axis=1.5,color_lab=color_lab)
	
	for (q in 1:ntot){
		reihen[9,q]=sd(reihen[c(1,3,5,7),q])
	}
	reihen2=array(reihen[9,],dim=c(1,ntot))
	col_row=c(4,1)
	legende=c(5,6,6,6)
	mat=c()
	for (row in 1:(col_row[1])){
		mat[((row-1)*7+1):(row*7)]=c(row,row,row,row,row,row,legende[row])
	}
	filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/",trendID,"_vergleich_",name_zusatz,"_",period,"_",titel_zusatz,"_kompakt_sd(seas_anom).pdf",sep="")
	map_allgemein(dat=dat,filename_plot=filename_plot,worldmap=worldmap,reihen=reihen2,titel=c(""),farb_mitte=c(0,1.2),farb_palette="regenbogen",grid=grid,ausschnitt=c(30,80),col_row=c(1,1),cex_axis=1,cex=1,paper=c(8,3),color_lab="interseasonal variability")

}

#init
library(SDMTools)
source("functions_regional.r")
source("map_plot.r")
library(rworldmap)
library(fields)
worldmap = getMap(resolution = "low")
ntot=1319
dat=dat_load("../data/HadGHCND_TMean_data3D.day1-365.1950-2014.nc")
trendID="91_5"
states=2
yearperiod="1950-2014"
dataset="_TMean"


#plot_warm_cold_diff_duration(col_row=c(1,1),cex=0.6,dataset="_TX",additional_style="",season_auswahl=c(1,2,3,4,5),var1="dur_ana_full",ausschnitt=c(-80,80),farb_palette="lila-gruen",farb_mitte="0")

#plot_diff_trend_duration(version1="_mean_TX",version2="_mode_TX",quantiles=c(8),col_row=c(1,1),cex=0.6,ausschnitt=c(-80,80),farb_mitte=c(-3,3))

plot_temp_anom_trend(yearPeriod=c(1950,2014))
plot_temp_anom_trend(yearPeriod=c(1950,1980))
plot_temp_anom_trend(yearPeriod=c(1980,2014))

if (1==2){
	#plot_seasonal_anomalie_duration(trendID="91_5",dataset=dataset,additional_style="",farb_mitte="mean",farb_palette="lila-gruen",name_zusatz="seas_anom_to_DJF",titel_zusatz="95_quantile",var1="dur_ana_full",var_qua=5,season_auswahl=c(1,2,3),baseline="DJF",col_row=c(1,1),cex=1,ausschnitt=c(-80,80))
	#plot_seasonal_anomalie_duration(trendID="91_5",dataset=dataset,additional_style="",farb_mitte="mean",farb_palette="lila-gruen",name_zusatz="seas_anom_to_4sea",titel_zusatz="95_quantile",var1="dur_ana_full",var_qua=5,season_auswahl=c(1,2,3,4),baseline="4seasons",col_row=c(1,1),cex=1,ausschnitt=c(-80,80))
	#plot_seasonal_anomalie_duration(trendID="91_5",dataset=dataset,additional_style="",farb_mitte=c(-12,12),farb_palette="lila-gruen",name_zusatz="seas_anom_klein",titel_zusatz="95_quantile",var1="dur_ana_full",var_qua=5,season_auswahl=c(1,2,3,4),col_row=c(5,2),cex=0.6)
	#plot_seasonal_anomalie_duration(trendID="91_5",dataset=dataset,additional_style="",farb_mitte=c(-2,2),farb_palette="lila-gruen",name_zusatz="seas_anom_klein",titel_zusatz="mean",var1="dur_ana_full",var_qua=8,season_auswahl=c(1,2,3,4),col_row=c(5,2),cex=0.6,cex_axis=1)
	plot_seasonal_anomalie_duration(trendID="91_5",dataset=dataset,additional_style="",farb_mitte=c(-2,2),farb_palette="lila-gruen",name_zusatz="seas_anom",titel_zusatz="mean",var1="dur_ana_full",var_qua=8,season_auswahl=c(1,2,3,4),col_row=c(5,2),cex=0.6,cex_axis=1)

	#plot_seasonal_anomalie_duration(trendID="91_5",dataset="_TX",additional_style="",farb_mitte=c(-2,2),farb_palette="lila-gruen",name_zusatz="seas_anom",titel_zusatz="mean",var1="dur_ana_full",var_qua=8,season_auswahl=c(1,2,3,4),col_row=c(1,1),cex=1,ausschnitt=c(-80,80))

	#plot_seasonal_anomalie_duration(farb_mitte=c(-12,12),farb_palette="lila-gruen",name_zusatz="seas_anom",		titel_zusatz="95 quantile anomalie",var1="dur_ana_full",subIndex=c("c","d","e","f","g","h","i","j"),season_auswahl=c(1,2,3,4),col_row=c(5,2),cex=0.6)
}

