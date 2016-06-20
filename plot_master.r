

plot_maps <- function(period="1950-2014",file="_others",var="other_stuff",sub_auswahl=c(NA,NA),value_auswahl=c(1,2),sig_auswahl=c(NA,NA),value_zusatz=c("mean","sd"),name_zusatz="mean_sd",farb_mitte="mean",farb_palette="lila-gruen",sig_style=c(NA),signi_level=0.05,ntot=length(dat$ID),ID_select=1:ntot){
	# plot ordenary map
	# filenename of the data file is composed of different strings specifying method to calculate trends, detrending params, period ....
	# var is the name of variable in netcdf file
	# value_auswahl has to be checked for the netcdf file

	print(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	values<-var.get.nc(nc,var)

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
				    #if (value_auswahl[val] %in% c(2,6,8)){reihen[index,]=round(exp(-values[sea,,state,value_auswahl[val]])*100,01)}
				    reihen[index,]=values[sea,,state,value_auswahl[val]]
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
					    if (value_zusatz[1]!=""){titel[index]=paste(value_zusatz[val],"of the",sub_zusatz[sub_auswahl[sub]],"percentile of",state_names[state],"period duration in",season,"in",period)}
						if (length(farb_mitte)>1){farb_mitte_end[(2+2*(index-1)):(3+2*(index-1))]=farb_mitte[((val-1)*2+1):((val-1)*2+2)]}
					}
				}
			}
		}
	}
	if (length(farb_mitte)==1){farb_mitte_end=farb_mitte}
	if (is.na(pch_points[1])){pch_points=array(15,dim(reihen))}
	filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",name_zusatz,name_reg_zusatz,"_",period,"_",name_style,".pdf",sep="")
	print(filename_plot)
	topo_map_plot(filename_plot=filename_plot,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte_end,farb_palette=farb_palette,signi_level=signi_level,ID_select=ID_select) #,reihen_sig=attribution_changes[,]
}

plot_diff_maps <- function(trendID="91_5",dataset="_TMean",additional_style="",period="1950-2014",file="_fit_2expo",var="fit_stuff",value1_auswahl=c(2),value2_auswahl=c(4),value_zusatz=c("b1-b2"),name_zusatz="diffB",farb_mitte="mean",farb_palette="lila-gruen",ID_select=1:ntot){
	# plot diff map
	# filenename of the data file is composed of different strings specifying method to calculate trends, detrending params, period ....
	# var is the name of variable in netcdf file
	# values from value2_auswahl will be subtracted from values of value1_auswahl

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
	map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/duration/",period,"/","duration_trend_",trendID,"_",season,"_",name_zusatz,"_",period,additional_style,".pdf",sep=""),worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette=farb_palette,grid=grid,yAusschnitt=yAusschnitt,region=region,col_row=col_row,mat=mat,paper=paper,pointsize=pointsize,subIndex=subIndex,ID_select=ID_select)
}



value <-function(x){return(x)}
inverse <-function(x){return(1/x)}

plot_reg_maps <- function(region_name="ward23",file="_others",var="other_stuff",sub_auswahl=c(NA,NA),value_auswahl=c(1,2),value2_auswahl=c(NA),sig_auswahl=c(NA,NA),value_zusatz=c("mean","sd"),name_zusatz="mean_sd",farb_mitte="mean",farb_palette="lila-gruen",sig_style=c(NA),signi_level=0.05,ntot=1319,reg_select=1:24,operation=value){
	# plot regional map
	# filenename of the data file is composed of different strings specifying method to calculate trends, detrending params, period ....
	# var is the name of variable in netcdf file
	# value_auswahl has to be checked for the netcdf file
	# if value2_auswahl!=NA than it is a diff map between value2_auswahl and value_auswahl
	# inverse of value can be plotted with operation=inverse

    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    regNumb<-length(unique(attribution[!is.na(attribution)]))

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
		            for (reg in reg_select){
		                if (is.na(value2_auswahl[1])){reihen[index,which(attribution==reg)]=operation(values[sea,reg,state,value_auswahl[val]])}
		                if (!is.na(value2_auswahl[1])){reihen[index,which(attribution==reg)]=operation(values[sea,reg,state,value_auswahl[val]])-operation(values[sea,reg,state,value2_auswahl[val]])}
		            }
				    if (!is.na(sig_auswahl[val])){
				    	for (reg in reg_select){
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
			            for (reg in reg_select){
			                if (is.na(value2_auswahl[1])){reihen[index,which(attribution==reg)]=operation(values[sea,reg,state,sub_auswahl[sub],value_auswahl[val]])}
			                if (!is.na(value2_auswahl[1])){reihen[index,which(attribution==reg)]=operation(values[sea,reg,state,sub_auswahl[sub],value_auswahl[val]])-operation(values[sea,reg,state,sub_auswahl[sub],value2_auswahl[val]])}
			            }
					    if (!is.na(sig_auswahl[val])){
					    	for (reg in reg_select){
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
	filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/","duration_trend_",trendID,"_",region_name,"_",name_zusatz,name_reg_zusatz,"_",period,"_",name_style,".pdf",sep="")
	print(filename_plot)
	greyLand<<-FALSE

	reihen_sig<- -reihen_sig
	topo_map_plot(filename_plot=filename_plot,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte_end,farb_palette=farb_palette,signi_level=signi_level,region_name=region_name) ; greyLand<<-TRUE
}

plot_state_mean_maps <- function(period="1950-2014",file="_others",var="other_stuff",sub_auswahl=c(NA,NA),value_auswahl=c(1,2),sig_auswahl=c(NA,NA),value_zusatz=c("mean","sd"),name_zusatz="mean_sd",farb_mitte="mean",farb_palette="lila-gruen",signi_level=0.05,ntot=length(dat$ID),ID_select=1:ntot){

	print(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	values<-var.get.nc(nc,var)

	reihen<-array(NA,dim=c(length(season_auswahl)*length(value_auswahl)*length(sub_auswahl),ntot))
	reihen_sig<-array(NA,dim=c(length(season_auswahl)*length(value_auswahl)*length(sub_auswahl),ntot))
	if (length(farb_mitte)>1){farb_mitte_end<-c(999)}
	titel<-c("")	
	index<-0
	for (sea in season_auswahl){
		season<-season_names[sea]
		if (is.na(sub_auswahl[1])){
			for (val in 1:length(value_auswahl)){
 				index<-index+1
				reihen[index,]=(values[sea,,1,value_auswahl[val]]+values[sea,,2,value_auswahl[val]])/2
				if (length(farb_mitte)>1){farb_mitte_end[(2+2*(index-1)):(3+2*(index-1))]=farb_mitte[((val-1)*2+1):((val-1)*2+2)]}
			}
		}
		# for quantile plots for example
		if (!is.na(sub_auswahl[1])){
			for (sub in 1:length(sub_auswahl)){
				for (val in 1:length(value_auswahl)){
					index<-index+1
					reihen[index,]=(values[sea,,1,sub_auswahl[sub],value_auswahl[val]]+values[sea,,2,sub_auswahl[sub],value_auswahl[val]])/2
					if (length(farb_mitte)>1){farb_mitte_end[(2+2*(index-1)):(3+2*(index-1))]=farb_mitte[((val-1)*2+1):((val-1)*2+2)]}
				}
			}
		}
	}
	if (length(farb_mitte)==1){farb_mitte_end=farb_mitte}
	filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/","duration_",trendID,"_",name_zusatz,name_reg_zusatz,"_",period,additional_style,".pdf",sep="") ; print(filename_plot)
	topo_map_plot(filename_plot=filename_plot,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte_end,farb_palette=farb_palette,signi_level=signi_level,ID_select=ID_select)
}

plot_state_diff_maps <- function(period="1950-2014",file="_others",var="other_stuff",sub_auswahl=c(NA,NA),value_auswahl=c(1,2),sig_auswahl=c(NA,NA),value_zusatz=c("mean","sd"),name_zusatz="mean_sd",farb_mitte="mean",farb_palette="lila-gruen",signi_level=0.05,ntot=length(dat$ID),ID_select=1:ntot){

	print(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	values<-var.get.nc(nc,var)

	reihen<-array(NA,dim=c(length(season_auswahl)*length(value_auswahl)*length(sub_auswahl),ntot))
	reihen_sig<-array(NA,dim=c(length(season_auswahl)*length(value_auswahl)*length(sub_auswahl),ntot))
	if (length(farb_mitte)>1){farb_mitte_end<-c(999)}
	titel<-c("")	
	index<-0
	for (sea in season_auswahl){
		season<-season_names[sea]
		if (is.na(sub_auswahl[1])){
			for (val in 1:length(value_auswahl)){
 				index<-index+1
				reihen[index,]=values[sea,,1,value_auswahl[val]]-values[sea,,2,value_auswahl[val]]
				if (length(farb_mitte)>1){farb_mitte_end[(2+2*(index-1)):(3+2*(index-1))]=farb_mitte[((val-1)*2+1):((val-1)*2+2)]}

			}
		}
		# for quantile plots for example
		if (!is.na(sub_auswahl[1])){
			for (sub in 1:length(sub_auswahl)){
				for (val in 1:length(value_auswahl)){
					index<-index+1
					reihen[index,]=values[sea,,1,sub_auswahl[sub],value_auswahl[val]]-values[sea,,2,sub_auswahl[sub],value_auswahl[val]]
					if (length(farb_mitte)>1){farb_mitte_end[(2+2*(index-1)):(3+2*(index-1))]=farb_mitte[((val-1)*2+1):((val-1)*2+2)]}
				}
			}
		}
	}
	if (length(farb_mitte)==1){farb_mitte_end=farb_mitte}
	filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/","duration_",trendID,"_",name_zusatz,name_reg_zusatz,"_",period,additional_style,".pdf",sep="") ; print(filename_plot)
	topo_map_plot(filename_plot=filename_plot,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte_end,farb_palette=farb_palette,signi_level=signi_level,ID_select=ID_select)
}


plot_seasonal_anomaly_maps <- function(period="1950-2014",file="_others",var="other_stuff",sub_auswahl=c(NA,NA),value_auswahl=c(1,2),sig_auswahl=c(NA,NA),value_zusatz=c("mean","sd"),name_zusatz="mean_sd",farb_mitte="mean",farb_palette="lila-gruen",signi_level=0.05,ntot=length(dat$ID),ID_select=1:ntot){

	print(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	values<-var.get.nc(nc,var)

	reihen<-array(NA,dim=c(4*length(value_auswahl)*length(sub_auswahl),ntot))
	reihen_sig<-array(NA,dim=c(length(season_auswahl)*length(value_auswahl)*length(sub_auswahl),ntot))
	if (length(farb_mitte)>1){farb_mitte_end<-c(999)}
	titel<-c("")	
	index<-0
	for (sea in 1:4){
		season<-season_names[sea]
		if (is.na(sub_auswahl[1])){
			for (val in 1:length(value_auswahl)){
 				index<-index+1
				reihen[index,]=(values[sea,,1,value_auswahl[val]]+values[sea,,2,value_auswahl[val]])/2 - (values[5,,1,value_auswahl[val]]+values[5,,2,value_auswahl[val]])/2
				if (length(farb_mitte)>1){farb_mitte_end[(2+2*(index-1)):(3+2*(index-1))]=farb_mitte[((val-1)*2+1):((val-1)*2+2)]}

			}
		}
		# for quantile plots for example
		if (!is.na(sub_auswahl[1])){
			for (sub in 1:length(sub_auswahl)){
				for (val in 1:length(value_auswahl)){
					index<-index+1
					reihen[index,]=(values[sea,,1,sub_auswahl[sub],value_auswahl[val]]+values[sea,,2,sub_auswahl[sub],value_auswahl[val]])/2 - (values[5,,1,sub_auswahl[sub],value_auswahl[val]]+values[5,,2,sub_auswahl[sub],value_auswahl[val]])/2
					if (length(farb_mitte)>1){farb_mitte_end[(2+2*(index-1)):(3+2*(index-1))]=farb_mitte[((val-1)*2+1):((val-1)*2+2)]}
				}
			}
		}
	}
	if (length(farb_mitte)==1){farb_mitte_end=farb_mitte}
	filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/","duration_",trendID,"_",name_zusatz,name_reg_zusatz,"_",period,additional_style,".pdf",sep="") ; print(filename_plot)
	topo_map_plot(filename_plot=filename_plot,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte_end,farb_palette=farb_palette,signi_level=signi_level,ID_select=ID_select)
}


plot_boxplot <- function(quans,x,width,color="white",border="black",density=NA){

    links=x-width/2
    rechts=x+width/2
    polygon(x=c(rechts,links,links,rechts),y=c(quans[2],quans[2],quans[4],quans[4]),col=color,border=border,density=density)    
    #polygon(x=c(rechts,links,links,rechts),y=c(quans[2],quans[2],quans[4],quans[4]),col=color,border=border,density=density)    
    for (qu in quans[1:5]){
        lines(c(links,rechts),c(qu,qu),col=border)
    }
    lines(c(links,links),quans[c(2,4)],col=border)
    lines(c(rechts,rechts),quans[c(2,4)],col=border)
    lines(c(x,x),quans[c(1,2)],lty=2,col=border)
    lines(c(x,x),quans[c(4,5)],lty=2,col=border)
}

plot_reg_boxplots <- function(region_name="ward24",file="_quantiles",var="quantile_stuff",name_zusatz="quans",ID_select=1:24,hlines=c(30),ID_length=length(ID_select)){

    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    regNumb<-length(unique(attribution[!is.na(attribution)]))

	print(paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,file,".nc",sep=""))
	nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,file,".nc",sep=""))
	values<-var.get.nc(nc,var)

	filename_plot<-paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/","duration_trend_",trendID,"_",region_name,"_",name_zusatz,name_reg_zusatz,"_",period,additional_style,"_boxplot.pdf",sep="") ; print(filename_plot)
	pdf(file=filename_plot,width=6,height=3)
	par(mar=c(1,3,1,1))

	colors<-c("blue","red")

		
	for (sea in season_auswahl){
		plot(NA,xlim=c(0.5,ID_length),ylim=c(-1,25),frame.plot=FALSE,axes=FALSE,xlab="",ylab="")
		axis(2)
		for (i in 1:ID_length){
			for (state in 1:2){
				plot_boxplot(values[sea,ID_select[i],state,,1],i+0.2*(-1)^state,0.3,color=colors[state])
			}
			text(i,0,ID_select[i])
		}
		for (i in 1:length(ID_select)){if (ID_select[i] %in% hlines){abline(v=i+0.5,col="gray",lty=2,lwd=2)}}
	}
	graphics.off()
}

zonaly_averaged_plot <- function(period="1950-2014",file="_others",var="other_stuff",sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),name_zusatz="mean_sd",farb_mitte="mean",farb_palette="lila-gruen",sig_style=c(NA),signi_level=0.05,ntot=length(dat$ID),ID_select=1:ntot){

	nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,file,".nc",sep=""))
	values<-var.get.nc(nc,var)


	pdf(file=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",name_zusatz,name_reg_zusatz,"_",period,"_",name_style,"zonally.pdf",sep=""))
	plot(NA,xlim=c(3,8),ylim=c(0,90),xlab="mean persistence",ylab="latitude")
	#plot(NA,ylim=c(3,8),xlim=c(0,90))
	color=c("green","orange","brown","blue")
	for (sea in 1:4){
		zoMea<-array(NA,36)
		index<-0
		for (lat in seq(0,90,2.5)){
			inside<-which(dat$lat==lat)
			zoMea[index<-index+1]=mean(values[sea,inside,,1],na.rm=TRUE)
		}
		nona<-which(!is.na(zoMea))
		#lines(smooth.spline(x=seq(0,90,2.5)[nona],y=zoMea[nona]),col=color[sea])
		#lines(seq(0,90,2.5)[nona],zoMea[nona],col=color[sea])
		lines(zoMea[nona],seq(0,90,2.5)[nona],col=color[sea])
		print(zoMea)
	}
	graphics.off()

}

