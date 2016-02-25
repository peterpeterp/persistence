

plot_maps <- function(trendID="91_5",dataset="_TMean",additional_style="",period="1950-2014",file="_others",var="other_stuff",
	ausschnitt=c(-80,80),season_auswahl=c(1,2,3,4,5),sub_auswahl=c(NA,NA),value_auswahl=c(1,2),sig_auswahl=c(NA,NA),value_zusatz=c("mean","sd"),sub_zusatz=c("",""),name_zusatz="mean_sd",farb_mitte="mean",farb_palette="lila-gruen",grid=FALSE,region=NA,col_row=c(1,1),mat=NA,paper=c(12,8),pointsize=1.2,subIndex=c("a","b"),sig_style=c(NA),signi_level=0.05){

	season_names=c("MAM","JJA","SON","DJF","4seasons")
	state_names=c("cold","warm")

	nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	values=var.get.nc(nc,var)

	values<<-values
	adas

	if (!is.na(sig_style[1])){
		nc_2=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,sig_style[2],".nc",sep=""))
		values2=var.get.nc(nc_2,var)
	}
	reihen=array(NA,dim=c(length(season_auswahl)*length(value_auswahl)*length(sub_auswahl)*states,ntot))
	reihen_sig=array(NA,dim=c(length(season_auswahl)*length(value_auswahl)*length(sub_auswahl)*states,ntot))
	if (length(farb_mitte)>1){farb_mitte_end=c(999)}
	titel=c("")	
	index=0
	for (sea in season_auswahl){
		season=season_names[sea]
		if (is.na(sub_auswahl[1])){
			for (val in 1:length(value_auswahl)){
				for (state in 1:states){
					index<-index+1
				    if (value_auswahl[val] %in% c(2,6,8)){reihen[index,]=round(exp(-values[sea,,state,value_auswahl[val]])*100,01)}
				    else{reihen[index,]=values[sea,,state,value_auswahl[val]]}
				    if (!is.na(sig_style[1])){reihen_sig[index,]=values[sea,,state,20]-values2[sea,,state,20]}
				    if (!is.na(sig_auswahl[val])){reihen_sig[index,]=values[sea,,state,sig_auswahl[val]]}
				    if (value_zusatz[1]!=""){titel[index]=paste(value_zusatz[val],"of",state_names[state],"period duration in",season,"in",period)}
					if (length(farb_mitte)>1){farb_mitte_end[(2+2*(index-1)):(3+2*(index-1))]=farb_mitte[((val-1)*2+1):((val-1)*2+2)]}

				}
			}
		}
		# for quantile plots for example
		if (!is.na(sub_auswahl[1])){
			for (sub in 1:length(sub_auswahl)){
				for (val in 1:length(value_auswahl)){
					for (state in 1:states){
						index<-index+1
					    reihen[index,]=values[sea,,state,sub_auswahl[sub],value_auswahl[val]]
					    if (!is.na(sig_auswahl[val])){reihen_sig[index,]=values[sea,,state,sub_auswahl[sub],sig_auswahl[val]]}
					    if (value_zusatz[1]!=""){titel[index]=paste(value_zusatz[val],"of the",sub_zusatz[sub],"percentile of",state_names[state],"period duration in",season,"in",period)}
					}
				}
			}
		}
	}
	if (length(farb_mitte)==1){farb_mitte_end=farb_mitte}
	topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/","duration_trend_",trendID,"_",season,"_",name_zusatz,"_",period,additional_style,".pdf",sep=""),reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte_end,farb_palette=farb_palette,grid=grid,region=region,pointsize=0.85,ausschnitt=c(-90,90),paper=c(7,5),signi_level=signi_level) #,reihen_sig=attribution_changes[,]

	#map_allgemein(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/","duration_trend_",trendID,"_",season,"_",name_zusatz,"_",period,additional_style,".pdf",sep=""),reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte_end,farb_palette=farb_palette,grid=grid,ausschnitt=ausschnitt,region=region,col_row=col_row,mat=mat,paper=paper,pointsize=pointsize,subIndex=subIndex,signi_level=signi_level)
	#map_allgemein(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/","duration_trend_",trendID,"_",season,"_",name_zusatz,"_",period,additional_style,".pdf",sep=""),reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte_end,farb_palette=farb_palette,grid=grid,ausschnitt=ausschnitt,region=region,col_row=col_row,mat=mat,paper=paper,pointsize=pointsize,subIndex=subIndex,signi_level=signi_level)
}

plot_diff_maps <- function(trendID="91_5",dataset="_TMean",additional_style="",period="1950-2014",file="_fit_2expo",var="fit_stuff",
	ausschnitt=c(-80,80),season_auswahl=c(1,2,3,4,5),value1_auswahl=c(2),value2_auswahl=c(4),value_zusatz=c("b1-b2"),name_zusatz="diffB",farb_mitte="mean",farb_palette="lila-gruen",grid=FALSE,region=NA,col_row=c(1,1),mat=NA,paper=c(12,8),pointsize=1.2,subIndex=c("a","b")){

	season_names=c("MAM","JJA","SON","DJF","4seasons")
	state_names=c("cold","warm")
	vars=c("dur_ana_full")
	print(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	values=var.get.nc(nc,var)

	reihen=array(NA,dim=c(length(season_auswahl)*length(value1_auswahl)*states,ntot))
	reihen_sig=array(NA,dim=c(length(season_auswahl)*length(value1_auswahl)*states,ntot))
	titel=c("")	
	index=0
	for (sea in season_auswahl){
		season=season_names[sea]
		for (val in 1:length(value1_auswahl)){
			for (state in 1:states){
				index<-index+1
				reihen[index,]=values[sea,,state,value1_auswahl[val]]-values[sea,,state,value2_auswahl[val]]
				if (value_zusatz[1]!=""){titel[index]=paste(value_zusatz[val],"of",state_names[state],"period duration in",season,"in",period)}
			}
		}
	}
	map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/duration/",period,"/","duration_trend_",trendID,"_",season,"_",name_zusatz,"_",period,additional_style,".pdf",sep=""),worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette=farb_palette,grid=grid,ausschnitt=ausschnitt,region=region,col_row=col_row,mat=mat,paper=paper,pointsize=pointsize,subIndex=subIndex)
}

plot_fit_diff_maps <- function(trendID="91_5",dataset="_TMean",additional_style="",period1="1950-2014",period2="1950-2014",file1="_fit_2expo",file2="_fit_2expo_restrict",var="fit_stuff",ausschnitt=c(-80,80),season_auswahl=c(1,2,3,4,5),value_auswahl=c(20,19,2),value_zusatz=c("BIC","R2","b1"),name_zusatz="diff_exp_2exp_restricted",farb_mitte="0",farb_palette="lila-gruen",grid=FALSE,region=NA,col_row=c(1,1),mat=NA,paper=c(12,8),pointsize=1.2,subIndex=c("a","b")){

	season_names=c("MAM","JJA","SON","DJF","4seasons")
	state_names=c("cold","warm")
	vars=c("dur_ana_full")
	print(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period1,"/",trendID,"_",dataset,"_",period1,file1,".nc",sep=""))
	nc1=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period1,"/",trendID,"_",dataset,"_",period1,file1,".nc",sep=""))
	nc2=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period2,"/",trendID,"_",dataset,"_",period2,file2,".nc",sep=""))
	values1=var.get.nc(nc1,var)
	values2=var.get.nc(nc2,var)

	reihen=array(NA,dim=c(length(season_auswahl)*length(value_auswahl)*states,ntot))
	reihen_sig=array(NA,dim=c(length(season_auswahl)*length(value_auswahl)*states,ntot))
	titel=c("")	
	if (length(farb_mitte)>1){farb_mitte_end=c(999)}
	index=0
	for (sea in season_auswahl){
		season=season_names[sea]
		for (val in 1:length(value_auswahl)){
			for (state in 1:states){
				index<-index+1
				reihen[index,]=values1[sea,,state,value_auswahl[val]]-values2[sea,,state,value_auswahl[val]]
				if (value_zusatz[1]!=""){titel[index]=paste("difference in",value_zusatz[val],"of",state_names[state],"period duration in",season,"in",period)}
				if (length(farb_mitte)>1){farb_mitte_end[(2+2*(index-1)):(3+2*(index-1))]=farb_mitte[((val-1)*2+1):((val-1)*2+2)]}
			}
		}
	}
	if (length(farb_mitte)==1){farb_mitte_end=farb_mitte}
	if (period1==period2){period=paste(period1,"/",sep="")}
	if (period1!=period2){period=paste(period1,"-",period2,"_",sep="")}
	map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/duration/",period,"duration_trend_",trendID,"_",season,"_",name_zusatz,"_",additional_style,".pdf",sep=""),worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte_end,farb_palette=farb_palette,grid=grid,ausschnitt=ausschnitt,region=region,col_row=col_row,mat=mat,paper=paper,pointsize=pointsize,subIndex=subIndex)
}


source("load.r")
library(SDMTools)
source("map_plot.r")
library(rworldmap)
library(fields)
worldmap = getMap(resolution = "low")

ntot=1319
states=2
trendID="91_5"
dataset="_TMean"
additional_style=""

dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
#plot_maps()
#plot_maps(file="_quantiles",var="quantile_stuff",sub_auswahl=c(5,7),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("quantile"),sub_zusatz=c("95th","100th"),name_zusatz="quantile",farb_mitte="mean",farb_palette="regenbogen")
#plot_maps(file="_quantiles",var="quantile_stuff",sub_auswahl=c(5,7),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c("gr slope"),sub_zusatz=c("95th","100th"),name_zusatz="qr_slope",farb_mitte="0")
#plot_maps(file="_fit_2expo_thresh_5-15",var="fit_stuff",sub_auswahl=c(NA),value_auswahl=c(2,4,5,18,19,20),sig_auswahl=c(NA,NA,NA,NA,NA),value_zusatz=c("b1","b2","thresh","distr_size","R2","BIC"),sub_zusatz=c(NA),name_zusatz="expo")

period=c("1950-2014","1980-2014")
for (i in c(1)){
	print(period[i])
	plot_maps(file="_fit_2expo_4:100",var="fit_stuff",sub_auswahl=c(NA),value_auswahl=c(6,8,9,14),sig_auswahl=c(17,17,17,17),value_zusatz=c("P1","P2","threshold","distr_size"),sub_zusatz=c(NA),name_zusatz="fit_2expo_4:100",period=period[i],signi_level=0,farb_mitte=c(70,90,70,90,5,15,20,50),farb_palette="spacy")
	plot_maps(file="_quantiles",var="quantile_stuff",sub_auswahl=c(5,7),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("quantile"),sub_zusatz=c("95th"),name_zusatz="quantile",farb_mitte="mean",farb_palette="regenbogen")
	plot_maps(file="_quantiles",var="quantile_stuff",sub_auswahl=c(5,7),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c("qr slope"),sub_zusatz=c("95th","100th"),name_zusatz="qr_slope",farb_mitte="0",signi_level=0.05)
	plot_maps(file="_others",var="other_stuff",sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),sub_zusatz=c(NA),name_zusatz="mean",period=period[i],signi_level=0.05,farb_mitte="mean",farb_palette="regenbogen")
	plot_maps(file="_others",var="other_stuff",sub_auswahl=c(NA),value_auswahl=c(4),sig_auswahl=c(10),value_zusatz=c("linear regression"),sub_zusatz=c(NA),name_zusatz="lm",period=period[i],signi_level=0.05,farb_mitte="0",farb_palette="lila-gruen")
	
}


#location_view(regions=TRUE)