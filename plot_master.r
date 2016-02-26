

plot_maps <- function(period="1950-2014",file="_others",var="other_stuff",ausschnitt=c(-80,80),season_auswahl=c(1,2,3,4,5),sub_auswahl=c(NA,NA),value_auswahl=c(1,2),sig_auswahl=c(NA,NA),value_zusatz=c("mean","sd"),sub_zusatz=c("",""),name_zusatz="mean_sd",farb_mitte="mean",farb_palette="lila-gruen",grid=FALSE,region=NA,col_row=c(1,1),mat=NA,paper=c(12,8),pointsize=1.2,subIndex=c("a","b"),sig_style=c(NA),signi_level=0.05,ntot=1319){

	print(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	values<-var.get.nc(nc,var)

	if (!is.na(sig_style[1])){
		nc_2<-open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,sig_style[2],".nc",sep=""))
		values2<-var.get.nc(nc_2,var)
	}
	reihen<-array(NA,dim=c(length(season_auswahl)*length(value_auswahl)*length(sub_auswahl)*2,ntot))
	reihen_sig<-array(NA,dim=c(length(season_auswahl)*length(value_auswahl)*length(sub_auswahl)*2,ntot))
	if (length(farb_mitte)>1){farb_mitte_end<-c(999)}
	titel<-c("")	
	index<-0
	for (sea in season_auswahl){
		season<-season_names[sea]
		if (is.na(sub_auswahl[1])){
			for (val in 1:length(value_auswahl)){
				for (state in 1:2){
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
					for (state in 1:2){
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



plot_reg_map_old <- function(region_name="ward22",fit_style,region_names,regNumb,ID_select,period){
    # plots worldmap and colored regions on it
    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    print(paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style,".nc",sep=""))
    nc = open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style,".nc",sep=""))
    fit_stuff=var.get.nc(nc,"fit_stuff")

    filename<-paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",region_name,"_",period,"_seasons.pdf",sep="")

    reihen1<-array(NA,c(10,1319))
    reihen2<-array(NA,c(10,1319))
    reihen3<-array(NA,c(10,1319))
    reihen4<-array(NA,c(10,1319))
    index<-0
    for (sea in 1:5){
        for (state in 1:2){
            index<-index+1
            for (reg in 1:regNumb){
                reihen1[index,which(attribution==reg)]=fit_stuff[sea,reg,state,9]
                reihen2[index,which(attribution==reg)]=fit_stuff[sea,reg,state,17]
                reihen3[index,which(attribution==reg)]=fit_stuff[sea,reg,state,6]-fit_stuff[sea,reg,state,8]
                reihen4[index,which(attribution==reg)]=fit_stuff[sea,reg,state,2]
            }
                      
        }
    }
    topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",region_name,"_",period,"_thresh.pdf",sep=""),reihen=reihen1,farb_mitte=c(4,12),farb_palette="regenbogen",pointsize=0.8,ausschnitt=c(-90,90),paper=c(6.4,3.7),region_name=region_name,regNumb=regNumb,land_col=NA)
    topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",region_name,"_",period,"_dBIC.pdf",sep=""),reihen=reihen2,farb_mitte=c(-100,0),farb_palette="regenbogen",pointsize=0.8,ausschnitt=c(-90,90),paper=c(6.4,3.7),region_name=region_name,regNumb=regNumb,land_col=rgb(0,0,0,0))
    topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",region_name,"_",period,"_b1-b2.pdf",sep=""),reihen=reihen3,farb_mitte=c(-0.3,0.3),farb_palette="lila-gruen",pointsize=0.8,ausschnitt=c(-90,90),paper=c(6.4,3.7),region_name=region_name,regNumb=regNumb,land_col=rgb(0,0,0,0))
    topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",region_name,"_",period,"_b.pdf",sep=""),reihen=reihen4,farb_mitte=c(0.15,0.35),farb_palette="regenbogen",pointsize=0.8,ausschnitt=c(-90,90),paper=c(6.4,3.7),region_name=region_name,regNumb=regNumb,land_col=rgb(0,0,0,0))
    return()
}

plot_reg_maps <- function(region_name="ward23",period="1950-2014",file="_others",var="other_stuff",ausschnitt=c(-80,80),season_auswahl=c(1,2,3,4,5),sub_auswahl=c(NA,NA),value_auswahl=c(1,2),sig_auswahl=c(NA,NA),value_zusatz=c("mean","sd"),sub_zusatz=c("",""),name_zusatz="mean_sd",farb_mitte="mean",farb_palette="lila-gruen",grid=FALSE,region=NA,col_row=c(1,1),mat=NA,paper=c(12,8),pointsize=1.2,subIndex=c("a","b"),sig_style=c(NA),signi_level=0.05,ntot=1319){

    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    regNumb<-length(unique(attribution[!is.na(attribution)]))

    print("../data/_TMean/91_5/regional/7rect/1950-2014/91_5__TMean_7rect_1950-2014_others.nc")
	print(paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,file,".nc",sep=""))
	nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,file,".nc",sep=""))
	values<-var.get.nc(nc,var)

	if (!is.na(sig_style[1])){
		nc_2<-open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,sig_style[2],".nc",sep=""))
		values2<-var.get.nc(nc_2,var)
	}
	reihen<-array(NA,dim=c(length(season_auswahl)*length(value_auswahl)*length(sub_auswahl)*2,ntot))
	reihen_sig<-array(NA,dim=c(length(season_auswahl)*length(value_auswahl)*length(sub_auswahl)*2,ntot))
	if (length(farb_mitte)>1){farb_mitte_end<-c(999)}
	titel<-c("")	
	index<-0
	for (sea in season_auswahl){
		season<-season_names[sea]
		if (is.na(sub_auswahl[1])){
			for (val in 1:length(value_auswahl)){
				for (state in 1:2){
					index<-index+1
		            for (reg in 1:regNumb){
		                reihen[index,which(attribution==reg)]=values[sea,reg,state,value_auswahl[val]]
		            }
				    if (!is.na(sig_auswahl[val])){
				    	for (reg in 1:regNumb){
		                	reihen_sig[index,which(attribution==reg)]=values[sea,reg,state,sig_auswahl[val]]
		            	}
		            }
				    if (value_zusatz[1]!=""){titel[index]=paste(value_zusatz[val],"of",state_names[state],"period duration in",season,"in",period)}
					if (length(farb_mitte)>1){farb_mitte_end[(2+2*(index-1)):(3+2*(index-1))]=farb_mitte[((val-1)*2+1):((val-1)*2+2)]}
				}
			}
		}
		# for quantile plots for example
		if (!is.na(sub_auswahl[1])){
			for (sub in 1:length(sub_auswahl)){
				for (val in 1:length(value_auswahl)){
					for (state in 1:2){
						index<-index+1
			            for (reg in 1:regNumb){
			                reihen[index,which(attribution==reg)]=values[sea,reg,state,sub_auswahl[sub],value_auswahl[val]]
			            }
					    if (!is.na(sig_auswahl[val])){
					    	for (reg in 1:regNumb){
			                	reihen_sig[index,which(attribution==reg)]=values[sea,reg,state,sub_auswahl[sub],sig_auswahl[val]]
			            	}
			            }
					    if (value_zusatz[1]!=""){titel[index]=paste(value_zusatz[val],"of the",sub_zusatz[sub],"percentile of",state_names[state],"period duration in",season,"in",period)}
						if (length(farb_mitte)>1){farb_mitte_end[(2+2*(index-1)):(3+2*(index-1))]=farb_mitte[((val-1)*2+1):((val-1)*2+2)]}
					}
				}
			}
		}
	}
	if (length(farb_mitte)==1){farb_mitte_end=farb_mitte}
	topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/","duration_trend_",trendID,"_",region_name,"_",season,"_",name_zusatz,"_",period,additional_style,".pdf",sep=""),reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte_end,farb_palette=farb_palette,grid=grid,region=region,pointsize=0.85,ausschnitt=c(-90,90),paper=c(7,5),signi_level=signi_level) #,reihen_sig=attribution_changes[,]
}