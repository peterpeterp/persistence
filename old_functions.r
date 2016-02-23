#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ master.r ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

master_state_attribution_2_trends <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style="_run_median"){
    # calculate persistence 2 states
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend1=trend_load(paste("../data/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    trend2=trend_load(paste("../data/",trendID,"/",trendID,trend_style,dataset,"_run_median.nc",sep=""))
    detrended=dat$tas-trend1-trend2
    per=state_attribution(dat,detrended,nday,nyr,
        filename=paste("../data/",trendID,"/",trendID,trend_style,dataset,"_state_ind",additional_style,".nc",sep=""))
}


master_daily_median_on_detrended <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style="_daily_median"){
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    trend=trend_load(paste("../data/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    detrended=dat$tas-trend

    ntot = 1319
    dailymedian=array(NA,dim=c(ntot,365))
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
            if (which(length(!is.na(z))>100)){
                dailymedian[q,day]=median(z,na.rm=TRUE)
            }
        }
    }

    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)

    daily_med <- var.def.ncdf(name="_daily_median",units="medain value for each day",dim=list(ID,day), missval=-9999.0)
    filename=paste("../data/",trendID,"/",trendID,trend_style,dataset,additional_style,".nc",sep="")
    nc = create.ncdf(filename,daily_med)
    put.var.ncdf(nc,daily_med,dailymedian)
    close.ncdf(nc)
}

master_runmedian_on_detrended <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style="_run_median"){
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    trend=trend_load(paste("../data/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    detrended=dat$tas-trend

    trash = ((nyr-1)/2*365+(nday-1))
    ntot = length(dat$ID)
    runmedian=dat$tas*NA
    for (q in 1:ntot) {
        temp = r_calc_runmedian_2D(detrended[q,,],nday=nday,nyr=nyr)
        temp[1:trash]=NA
        temp[(length(dat$time)-trash):length(dat$time)]=NA
        runmedian[q,,]=temp
    }
    trend_write(paste("../data/",trendID,"/",trendID,trend_style,dataset,additional_style,".nc",sep=""),dat,runmedian)
}

master_duration_distribution <- function(yearPeriod,trendID,states,seasons=c("spring","summer","autumn","winter","year")
    ,trend_style="_mean",additional_style=""){
    # analyse duration periods 2 states
    for (season in seasons){
        print(season)
        nc=open.ncdf(paste("../data/",trendID,"/",additional_style,"/duration/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_duration_",states,"s_analysis_",season,".nc",sep=""))
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")
        duration_distribution(dur,dur_mid,filename=paste("../data/",trendID,"/",additional_style,"/duration/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_duration_",states,"s_distribution_",season,".nc",sep=""),
            season=season,yearPeriod)#,stations=seq(487,489,1))#,stations=134:136)
    }
}

master_regional_trend <- function(yearPeriod,region_name,trendID,trend_style="_mean",additional_style=""){
    source("functions_regional.r")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    regional_analysis(dat=dat,yearPeriod,filepath=paste("../data/",trendID,"/2_states",trend_style,additional_style,"/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",sep=""),region_name=region_name)
}


master_endaussage <- function(){
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    points=c(1950,2014,1950,1980,1980,2014)
    for (i in 1:3){
        period=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        print(paste(period,sep="-"))
        end_aussage(dat=dat,yearPeriod=paste(period[1],period[2],sep="-"),states=2,trendID="91_5")
    }
}

#master_runmedian_on_detrended(nday=nday,nyr=nyr,trendID=trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)
#master_daily_median_on_detrended(nday=nday,nyr=nyr,trendID=trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)


#master_state_attribution_2_trends(nday,nyr,trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ master.r ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#++++++++++++++++++++++++++++++++++++++ functions_regional.r +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
wave_region <- function(wavenumber,trough,lats){
    # creates a table with coordinates of rectangles representing wave ridges and troughs
    wellen=array(NA,dim=c(40,13))
    brei=360/wavenumber/4
    lon=trough
    lon0=trough
    i=1
    wellen[,7:8]=lats[1]
    wellen[,9:10]=lats[2]
    while ((round(x=lon,digits=3)==round(x=lon0,digits=3) & i>1)==FALSE){
        if (lon>(180-brei) & lon<(180+brei)){
            wellen[i,1:4]=c(180,lon-brei,lon-brei,180)
            wellen[i,13]=i  
            i=i+1   
            wellen[i,1:4]=c(lon+brei-360,-180,-180,lon+brei-360)
            wellen[i,13]=i  
            i=i+1   
            lon=lon+2*brei
        }
        if (lon<=(180-brei)){
            wellen[i,1:4]=c(lon+brei,lon-brei,lon-brei,lon+brei)
            wellen[i,13]=i  
            i=i+1   
            lon=lon+2*brei
        }
        if (lon>=(180+brei)){
            lon=lon-360
            wellen[i,1:4]=c(lon+brei,lon-brei,lon-brei,lon+brei)
            wellen[i,13]=i  
            i=i+1   
            lon=lon+2*brei
        }
        
    }
    write.table(wellen[1:length(which(!is.na(wellen[,13]))),],paste("../data/region_poligons/",wavenumber,"wave.txt",sep=""))
}


regions_color <- function(reihen,reihen_sig,worldmap,titles,poli,filename_plot){
    # plots worldmap and colored regions on it
    jet.colors <- colorRampPalette( c(rgb(0.2,0.6,0.2),rgb(0.5,1,0.5), rgb(0.98,0.98,0.98) ,rgb(1,0.5,1),rgb(0.6,0.2,0.6)))
    
    nbcol <- 101
    color <- jet.colors(nbcol)

    pdf(file = filename_plot,width=12,height=8)

    for (rei in 1:dim(reihen)[1]){            
        y=c()
        index=c()
        signi=c()
        j=0
        for (i in 1:dim(poli)[1]){
            poliLabel=i
            if (!is.na(reihen[rei,poliLabel])){
                j=j+1
                y[j]=reihen[rei,poliLabel]
                index[j]=i 
                if (abs(reihen[rei,poliLabel])>0.0001){
                    signi[j]=sprintf("%.04f",reihen_sig[rei,poliLabel])
                }         
            }
        }
        aushol=max(c(abs(max(y)),abs(min(y))))
        y[j+1]=-aushol
        y[j+2]=aushol
        facetcol <- cut(y,nbcol)  

        print(titles[rei])
        plot(worldmap,main=titles[rei])

        for (i in 1:j){
            lon=poli[index[i],1:6]
            lat=poli[index[i],7:12]
            lon=lon[!is.na(lon)]
            lat=lat[!is.na(lat)]
            polygon(x=lon,y=lat,col=color[facetcol[i]],border="green")
            text(mean(lon),mean(lat),label=signi[i],cex=0.7,col="black")
        }
        image.plot(legend.only=T, zlim=range(y), col=color)
    }
    graphics.off()
    return()
}

plot_regional_boxplots <- function(trendID,dat,yearPeriod,region_name,additional_style,dataset){
    nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",region_name,"_",yearPeriod[1],"-",yearPeriod[2],"_distributions.nc",sep=""))
    regions=get.var.ncdf(nc,"region")
    quantiles=array(get.var.ncdf(nc,"quantile"),dim=c(5,2,length(regions),10))

    color=c(rgb(0,0,1,0.6),rgb(1,0,0,0.4))

    plot_names=c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s")

    region_names=c("western \n N. America","central \n N. America","eastern \n N. America","Europe","western \n Asia","central \n Asia","eastern \n Asia")
    region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")

    season_names=c("spring","summer","autumn","winter","year")
    season_names=c("MAM","JJA","SON","DJF","year")

    # different order of seasonal_auswahl might result in strange things, dangerous
    season_auswahl=c(1,2,3,4,5)

    #regional focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_boxplots_reg.pdf",sep=""),width=4,height=12)
    par(mfrow=c(5,1))
    par(mar=c(1,5,4,3))

    for (sea in season_auswahl){
        plot(NA,xlim=c(0.5,length(regions)*0.65),ylim=c(-3,24),frame.plot=FALSE,axes=FALSE,ylab="# days",xlab="",main=season_names[sea])
        axis(2,ylim=c(-3,24))
        for (i in seq(0,25,5)){
            abline(h=i,col=rgb(0.8,0.8,0.8,0.6))
        }
        for (reg in regions){
            for (state in 1:2){
                quAn=quantiles[sea,state,reg,]
                #print(quAn)
                mitte=reg*0.6-0.1+0.2*(state-1)
                links=mitte-0.1
                rechts=mitte+0.1
                polygon(x=c(rechts,links,links,rechts),y=c(quAn[2],quAn[2],quAn[4],quAn[4]),col=color[state])
                lines(c(mitte,mitte),c(quAn[1],quAn[5]))
                lines(c(links,rechts),c(quAn[1],quAn[1]))
                lines(c(links,rechts),c(quAn[5],quAn[5]))
                lines(c(links,rechts),c(quAn[3],quAn[3]))
                points(mitte,quAn[10],pch=4)
            }
            text(reg*0.6,-3,region_names[reg])
        }
    }
    graphics.off()


    # seasonal focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_boxplots_sea.pdf",sep=""),width=4,height=12)
    par(mfrow=c(7,1))
    par(mar=c(1,5,4,3))

    for (reg in regions){
        plot(NA,xlim=c(0.5,length(season_auswahl)*0.65),ylim=c(-3,24),frame.plot=FALSE,axes=FALSE,ylab="# days",xlab="",main=region_names[reg])
        axis(2,ylim=c(-3,24))
        for (i in seq(0,25,5)){
            abline(h=i,col=rgb(0.8,0.8,0.8,0.6))
        }
        for (sea in season_auswahl){
            for (state in 1:2){
                quAn=quantiles[sea,state,reg,]
                #print(quAn)
                mitte=sea*0.6-0.1+0.2*(state-1)
                links=mitte-0.1
                rechts=mitte+0.1
                polygon(x=c(rechts,links,links,rechts),y=c(quAn[2],quAn[2],quAn[4],quAn[4]),col=color[state])
                lines(c(mitte,mitte),c(quAn[1],quAn[5]))
                lines(c(links,rechts),c(quAn[1],quAn[1]))
                lines(c(links,rechts),c(quAn[5],quAn[5]))
                lines(c(links,rechts),c(quAn[3],quAn[3]))
                points(mitte,quAn[10],pch=4)
            }
            text(sea*0.6,-3,season_names[sea])
        }
    }
    graphics.off()

    #regional focus normalized -> are the behaviours of different quantiles different?
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_normalized_reg.pdf",sep=""),width=8,height=12)
    par(mfrow=c(6,2))
    par(mar=c(1,5,4,3))

    qua_selection=c(1,2,3,4,5,10)
    qua_names=c("5 percentile","25 percentile","median","75 percentile","95 percentile","mean")
    qua_colors=c(rgb(0.1,0.1,0.6),rgb(0.2,0.4,0.7),rgb(0,0,0),rgb(0.7,0.2,0.4),rgb(1,0,0),rgb(0,1,0))
    qua_colors=c("lightblue","blue","black","orange","red","green")
    qua_ltys=c(1,2,1,2,1,1)
    qua_pchs=c(13,1,15,1,13,4)
    state_names=c("warm","cold")

    x=regions

    for (sea in season_auswahl){
        for (state in 1:2){
            plot(NA,xlim=c(0,length(regions)),ylim=c(-0.1,1),frame.plot=FALSE,axes=FALSE,ylab="normalized # days",xlab="",main=paste(season_names[sea],"  ",state_names[state]))
            axis(2,ylim=c(0,1))

            for (reg in regions){
                text(reg,-0.1,region_names[reg])
            }
            for (qua_index in 1:length(qua_selection)){
                qua=qua_selection[qua_index]
                y=quantiles[sea,state,,qua]
                y_norm=(y-min(y,na.rm=TRUE))/(max(y,na.rm=TRUE)-min(y,na.rm=TRUE))
                lines(x,y_norm,col=qua_colors[qua_index],lty=qua_ltys[qua_index])
                points(x,y_norm,col=qua_colors[qua_index],pch=qua_pchs[qua_index])
            }
        }
        
    }
    plot(NA,xlim=c(0,length(regions)),ylim=c(-0.1,1),frame.plot=FALSE,axes=FALSE,ylab="",xlab="",main="")
    legend("topleft",col=qua_colors,lty=qua_ltys,legend=qua_names,pch=qua_pchs)
    graphics.off()

    #seasonal focus normalized -> are the behaviours of different quantiles different?
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_normalized_sea.pdf",sep=""),width=8,height=12)
    par(mfrow=c(7,2))
    par(mar=c(1,5,4,3))

    qua_selection=c(1,2,3,4,5,10)
    qua_names=c("5 percentile","25 percentile","median","75 percentile","95 percentile","mean")
    qua_colors=c(rgb(0.1,0.1,0.6),rgb(0.2,0.4,0.7),rgb(0,0,0),rgb(0.7,0.2,0.4),rgb(1,0,0),rgb(0,1,0))
    qua_colors=c("lightblue","blue","black","orange","red","green")
    qua_ltys=c(1,2,1,2,1,1)
    qua_pchs=c(13,1,15,1,13,4)
    state_names=c("warm","cold")

    x=season_auswahl

    for (reg in regions){
        for (state in 1:2){
            plot(NA,xlim=c(0,length(regions)),ylim=c(-0.1,1),frame.plot=FALSE,axes=FALSE,ylab="normalized # days",xlab="",main=paste(region_names[reg],"  ",state_names[state]))
            axis(2,ylim=c(0,1))

            for (sea in season_auswahl){
                text(sea,-0.1,season_names[sea])
            }
            for (qua_index in 1:length(qua_selection)){
                qua=qua_selection[qua_index]
                y=quantiles[,state,reg,qua]
                y_norm=(y-min(y,na.rm=TRUE))/(max(y,na.rm=TRUE)-min(y,na.rm=TRUE))
                lines(x,y_norm,col=qua_colors[qua_index],lty=qua_ltys[qua_index])
                points(x,y_norm,col=qua_colors[qua_index],pch=qua_pchs[qua_index])
            }
        }
        if (2==1){
            plot(NA,xlim=c(0,length(regions)),ylim=c(-0.1,1),frame.plot=FALSE,axes=FALSE,ylab="",xlab="",main="")
            legend("topleft",col=qua_colors,lty=qua_ltys,legend=qua_names,pch=qua_pchs)
        }
        
    }
    graphics.off()
}


in plot_regional_boxplots

        if (1==2){
            tmp=boxplot(dists,at=at_,col=color,ylim=c(0,(max(quantiles[sea,,,])+4)),boxwex=0.3,names=c(region_names,1:regNumb*NA),outline=FALSE,frame.plot=FALSE,axes=FALSE,main=season)
            axis(2)
            for (reg in 1:regNumb){
                text(reg,max(quantiles[sea,,,])+2,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))
            }
            for (quA in 1:length(taus)){
                points(at_[1:regNumb],quantiles[sea,,1,quA],col="blue",pch=quA)
                lines(at_[1:regNumb],quantiles[sea,,1,quA],col="blue",lty=quA)
                points(at_[(regNumb+1):(regNumb*2)],quantiles[sea,,2,quA],col="red",pch=quA)
                lines(at_[(regNumb+1):(regNumb*2)],quantiles[sea,,2,quA],col="red",lty=quA)
            }
        }

#++++++++++++++++++++++++++++++++++++++ functions_regional.r +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#====================================== plot_master.r ==================================================================

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

#====================================== plot_master.r ==================================================================

#______________________________________ functions_correl.r _____________________________________________________________

markov_correl <- function(trendID,states,toCor,toCor_name,toCor_short,toCor_shortZu,additional_style,
	toCor_startYear=1950,transition_names,stations=seq(1,1319,1),plot=FALSE){
	# toCor is array with dim=c(ntot,seasons,years) or dim=c(seasons,years) 
	# function will follow different procedures depending on dim(toCor)

	# correlation is cov/(var*var)
	# also linear regression with corSlope and corSlope_sig

	# using toCor_startYear and toCor_years the matching years from markov are determined
	# from these years noNa are the years that will be plotted

	ntot=1319
	transNumb=states*states
	toCor_years=dim(toCor)[length(dim(toCor))]

	nc=open.ncdf(paste("../data/",trendID,"/",states,"_states",additional_style,"/markov/",trendID,"_markov_",states,"states.nc",sep=""))
	markov=get.var.ncdf(nc,"markov")

	# needed for the plots only --------------------------------
	if (plot==TRUE){
		pdf(file=paste("../plots/",trendID,"/",states,"_states",additional_style,"/stations/mar_",toCor_short,"_",toCor_shortZu,"_",stations,".pdf",sep=""))
		seasons=c("spring","summer","autumn","winter")
		if (states==2){
			transLongName=c("cold to cold","warm to cold","cold to warm","warm to warm")
			color=c("blue","violet","green","red")
		}
		if (states==3){
			transLongName=c("cold to cold","normal to cold","warm to cold","cold to normal","normal to normal","warm to normal","cold to warm","normal to warm","warm to warm")
		}
	}
	# ---------------------------------------------------------	

	correlation=array(NA,dim=c(ntot,5,transNumb))
	corSlope=array(NA,dim=c(ntot,5,transNumb))
	corSlope_sig=array(NA,dim=c(ntot,5,transNumb))
	x=seq(1,toCor_years,1)
	marSeq=(toCor_startYear-1950+1):(toCor_startYear+toCor_years-1950)
	for (q in stations){
		print(q)
		if (length(dim(toCor))==2){
			toCor_loc=toCor
		}
		if (length(dim(toCor))==3){
			toCor_loc=toCor[q,,]
		}
		for (sea in 1:5){	
			# needed for the plots only --------------------------------						
			if (plot==TRUE){
				plot(NA,ylim=c(0,1.1),xlim=c(min(toCor_loc[sea,],na.rm=TRUE),max(toCor_loc[sea,],na.rm=TRUE)),
					,xlab=toCor_name,ylab="transition probability",
					main=paste(seasons[sea]))
				title=c()

				noNa=which(!is.na(markov[q,sea,1,marSeq]))
				order=order(toCor_loc[sea,noNa])
				year=toCor_startYear:(toCor_startYear+toCor_years)
				year=year[noNa][order]
			}
			# ---------------------------------------------------------
			for (trans in 1:transNumb){
				if (length(which(!is.na(markov[q,sea,trans,marSeq])))>10 & length(which(!is.na(toCor_loc[sea,])))>10){

					if (sd(as.vector(toCor_loc[sea,]),na.rm=TRUE)!=0 & sd(as.vector(markov[q,sea,trans,marSeq]),na.rm=TRUE)!=0){
						correlation[q,sea,trans]=cor(x=as.vector(toCor_loc[sea,]),y=as.vector(markov[q,sea,trans,marSeq]),use="complete")
					}
					lr=summary(lm(toCor_loc[sea,]~markov[q,sea,trans,marSeq]))
					corSlope[q,sea,trans]=lr$coefficients[2,1]
					corSlope_sig[q,sea,trans]=lr$coefficients[2,4]

					if (plot==TRUE){
						lines(toCor_loc[sea,noNa][order],markov[q,sea,trans,marSeq][noNa][order],col=color[trans])
						abline(lm(markov[q,sea,trans,marSeq]~toCor_loc[sea,]),col=color[trans],lty=2)
						title[trans]=paste(transLongName[trans],"transition -",toCor_name ,"=",correlation[q,sea,trans])
					}
				}
			}
			# needed for the plots only --------------------------------
			if (plot==TRUE){
				shift=-1
				for (y in 1:toCor_years){
					if(shift==2){shift=-1}
					abline(v=as.vector(toCor_loc[sea,noNa][order])[y],col="grey",lty=2)
					text(as.vector(toCor_loc[sea,noNa][order])[y],(0.5+0.05*shift),label=year[y],srt=90,cex=0.5)
					shift=shift+1
				}
				legend("topright",legend=title,lty=array(1,transNumb),col=color)
			}
			# ---------------------------------------------------------
		}
	}

	if (plot!=TRUE){
	    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
	    season <- dim.def.ncdf("season",units="uu",vals=1:5,unlim=FALSE)
	    transition <- dim.def.ncdf("transition",units=transition_names,vals=1:transNumb,unlim=FALSE)

	    varCorrelation <- var.def.ncdf(name="correlation",units="bla",dim=list(ID,season,transition), missval=-9999.0)
	    varCorSlope <- var.def.ncdf(name="corSlope",longname="linear regression like correlation",units="bla",dim=list(ID,season,transition), missval=-9999.0)
	    varCorSlope_sig <- var.def.ncdf(name="corSlope_sig",longname="linear regression like correlation_sig",units="bla",dim=list(ID,season,transition), missval=-9999.0)
	    
	    vars=list(varCorrelation,varCorSlope,varCorSlope_sig)  
	    nc = create.ncdf(paste("../data/",trendID,"/",states,"_states",additional_style,"/correlations/",trendID,"_",toCor_short,"_",toCor_shortZu,"_markov_cor_",states,"states.nc",sep=""),vars)

		put.var.ncdf(nc,varCorrelation,correlation)
		put.var.ncdf(nc,varCorSlope,corSlope)
		put.var.ncdf(nc,varCorSlope_sig,corSlope_sig)

		close.ncdf(nc)
	}
}


mar_correlation_plot <- function(additional_style,trendID="91_5",states=2,toCor_short="nao",toCor_name="NAO",toCor_shortZu="",
	transition_names=c("cc","wc","cw","ww"),seasons=c("spring","summer","autumn","winter","year"),
	dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc"),worldmap = getMap(resolution = "low"),ntot=1319){

    nc=open.ncdf(paste("../data/",trendID,"/",states,"_states",additional_style,"/correlations/",trendID,"_",toCor_short,"_",toCor_shortZu,"_markov_cor_",states,"states.nc",sep=""))
    correlation=get.var.ncdf(nc,"correlation")
    corSlope=get.var.ncdf(nc,"corSlope")
    corSlope_sig=get.var.ncdf(nc,"corSlope_sig")
	reihen=array(NA,dim=c(20,ntot))
	reihen_sig=array(NA,dim=c(20,ntot))
	titel=c()

	for (sea in 1:5){
		for (trans in 1:4){
			reihen[((sea-1)*4+trans),]=corSlope[,sea,trans]
			reihen_sig[((sea-1)*4+trans),]=corSlope_sig[,sea,trans]
			titel[((sea-1)*4+trans)]=paste("correlation between",toCor_name,"and transition probability",transition_names[trans],"in",seasons[sea])
		}
	}
    map_allgemein(dat,filename_plot=paste("../plots/",trendID,"/",states,"_states",additional_style,"/maps/mar_cor/",trendID,"_",toCor_short,"_",toCor_shortZu,"_markov_",states,"states.pdf",sep=""),
    	worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte="0",farb_palette="gold-blau")
}



eke_mar_correl <- function(additional_style,trendID="91_5",states=2,level=1,stations=seq(1,1319,1),plot=FALSE,transition_names="cc wc cw ww"){
	nc=open.ncdf("../data/eke/eke_ID.nc")
	eke=get.var.ncdf(nc,"eke_sea")
	pressure=get.var.ncdf(nc,"levelist")

	markov_correl(additional_style=additional_style,trendID=trendID,states=states,stations=stations,plot=plot,toCor=eke[,level,,],toCor_name="eddy kynetic energy",toCor_short="eke",toCor_shortZu=paste(pressure[1],"mbar",sep=""),toCor_startYear=1979,transition_names=transition_names)
}


index_mar_correl <- function(additional_style,trendID="91_5",states=2,toCor_name="NAO",toCor_short="nao",stations=seq(1,1319,1),plot=FALSE,transition_names="cc wc cw ww"){
	nc=open.ncdf(paste("../data/roh_data/",toCor_name,"_sea.nc",sep=""))
	index_sea=get.var.ncdf(nc,toCor_short)

	markov_correl(additional_style=additional_style,trendID=trendID,states=states,stations=stations,plot=plot,toCor=index_sea,toCor_name=toCor_name,toCor_short=toCor_short,toCor_shortZu="",toCor_startYear=1950,transition_names=transition_names)
}

#______________________________________ functions_correl.r _____________________________________________________________

#************************************** functions_support.r *************************************************************

calc_per <- function(dat,trend,nday,nyr,model,states,transition_names,filename){
    source("functions_markov.r")

    ## User parameters 
    #trash is the number of data point which are wasted by detrending
    trash = ((nyr-1)/2*365+(nday-1))
    ntot = length(dat$ID)
    laenge_zeit = length(dat$time)
    transitions=states*states

    # Calculate persistence information
    #cat("Calculating persistence... ")

    markov_per = list(ind=dat$tas*NA,markov=array(NA,dim=c(ntot,5,transitions,65)),markov_conf=array(NA,dim=c(ntot,5,65)))

    for (q in 1:ntot) { 
        cat("-")
        if (length(which(is.na(dat$tas[q,,])))<(2*trash+365*20)){

            # Calculate persistence vector
            y = dat$tas[q,,]
            per_ind = y*NA 

            if (states==3){
                detrended = y-trend[q,,]
                threshold = sd(detrended,na.rm=TRUE)*0.5

                per_ind[detrended>threshold]  =  1
                per_ind[detrended<(-threshold)]  = -1 
                per_ind[detrended<threshold & detrended>(-threshold)] = 0 
            }

            if (states==2){
                per_ind[y > trend[q,,]]=1
                per_ind[y < trend[q,,]]=-1
            # the >= was somehow problematic, since it affects roughly 5% of the datapoints
            # now the datapoints sitting on the trend are randomly attributed to warm or cold
                per_ind[y == trend[q,,]]=0
                per_ind[per_ind==0]=sample(c(-1,1),1)
            }

            markov_per$ind[q,,] = per_ind
            # Go through vector and calculate persistent events 
            # Perform this step on a 1D vector to avoid artificial cutoffs 
            # at day 1 and day 365 of the year 
            per_ind1D = as.vector(per_ind) 

                
            tmp=seasonal_matrix_out(per_ind1D,model,states,array(c(60,150,151,241,242,333,334,424),dim=c(2,4)))
            for (i in 1:4){
                markov_per$markov[q,i,,]=tmp$out[i,,]
                markov_per$markov_conf[q,i,]=tmp$out_conf[i,]
            }

            tmp=seasonal_matrix_out(per_ind1D,model,states,array(c(1,365),dim=c(2,1)))
            markov_per$markov[q,5,,]=tmp$out[1,,]
            markov_per$markov_conf[q,5,]=tmp$out_conf[1,]

        } 
        else {
            cat(sprintf("> ID %s lon %s  lat %s <",dat$ID[q],dat$lon[q],dat$lat[q]))
        }
     
    }
    markov_write(filename,dat,markov_per,transitions,transition_names) 

    cat("done.\n")
    return(0)#list(markov_per=markov_per,shock_per=shock_per))
}

state_attribution_old <- function(dat,detrended,nday,nyr,filename){
    ## User parameters 
    #trash is the number of data point which are wasted by detrending
    trash = ((nyr-1)/2*365+(nday-1))
    ntot = length(dat$ID)
    laenge_zeit = length(dat$time)

    # Calculate persistence information
    #cat("Calculating persistence... ")

    state_ind=dat$tas*NA

    for (q in 1:ntot) { 
        cat("-")
        if (length(which(is.na(dat$tas[q,,])))<(2*trash+365*20)){

            # Calculate persistence vector
            y = dat$tas[q,,]
            per_ind = y*NA 

            threshold = median(detrended[q,,],na.rm=TRUE)
            per_ind[detrended[q,,] < threshold]=-1
            per_ind[detrended[q,,] > threshold]=1
            # the >= was somehow problematic, since it affects roughly 5% of the datapoints
            # now the datapoints sitting on the trend are randomly attributed to warm or cold
            per_ind[detrended[q,,] == threshold]=1
            #per_ind[per_ind==0]=sample(c(-1,1),1)

            state_ind[q,,]=per_ind

        } 
        else {
            cat(sprintf("> ID %s lon %s  lat %s <",dat$ID[q],dat$lon[q],dat$lat[q]))
        }
     
    }

    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    year <- dim.def.ncdf("year",units="year",vals=1:65, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)

    ind <- var.def.ncdf(name="ind",units="cold = -1 , warm = 1",dim=list(ID,day,year), missval=-9999.0)
    nc = create.ncdf(filename,ind)
    put.var.ncdf(nc,ind,state_ind)
    close.ncdf(nc)

    cat("done.\n")
    return(0)#list(markov_per=markov_per,shock_per=shock_per))
}

trend_analysis <- function(x,y){
    library(Kendall)
    if (length(which(is.na(y)))>7 | length(which(y==0))>40){
        return(list(slope=NA,slope_sig=NA,MK=NA,MK_sig=NA))
    }
    lm.r=lm(y~x)
    slope=summary(lm.r)$coefficients[2]
    slope_sig=summary(lm.r)$coefficients[8]
    out=MannKendall(y)
    MK=out[1]$tau
    MK_sig=out[2]$sl
    return(list(slope=slope,slope_sig=slope_sig,MK=MK,MK_sig=MK_sig))
}

global_analysis <- function(toAna,yearPeriod,yearshift=1949){
    #toANA of dim(ntot, something, year)
    yearPeriod=yearPeriod-yearshift
    t=seq(yearPeriod[1],yearPeriod[2],1)
    ntot=1319
    series=dim(toAna)[2]
    analysis=array(NA,dim=c(ntot,series,6))
    for (i in 1:series){
        for (q in 1:ntot){
            tmp=trend_analysis(t,toAna[q,i,yearPeriod[1]:yearPeriod[2]])
            analysis[q,i,1]=mean(toAna[q,i,yearPeriod[1]:yearPeriod[2]],na.rm=TRUE)
            analysis[q,i,2]=sd(toAna[q,i,yearPeriod[1]:yearPeriod[2]],na.rm=TRUE)
            analysis[q,i,3]=tmp$MK
            analysis[q,i,4]=tmp$MK_sig
            analysis[q,i,5]=tmp$slope
            analysis[q,i,6]=tmp$slope_sig
        }
    }
    return(analysis)
}

end_aussage <- function(dat,yearPeriod,trendID,states,seasons=c("spring","summer","autumn","winter","year"),region=c(-180,180,30,60)){
    # calculates percentage of grid points in region having positive trend
    # calculates big average
    print(yearPeriod)
    for (season in seasons){
        print(season)
        nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/duration/",yearPeriod,"/",trendID,"_duration_",states,"s_analysis_",season,".nc",sep=""))
        dur_ana_full=get.var.ncdf(nc,"dur_ana_full")
        #duration_analysis_write(paste("../data/",trendID,"/",states,"_states/duration/",yearPeriod,"/",trendID,"_duration_",states,"s_analysis_",season,".nc",sep=""),dur_ana_full,season)
        quantiles=get.var.ncdf(nc,"quantiles")
        outs=get.var.ncdf(nc,"outs")
        percentage=array(NA,dim=c(states,length(quantiles)))
        for (state in 1:states){
            for (quan in quantiles){
                slope=dur_ana_full[1:1319,state,quan,1]
                inside_region=which(dat$lon>region[1] & dat$lon<region[2] & dat$lat>region[3] & dat$lat< region[4])
                noNa=which(!is.na(slope[inside_region]))
                percentage[state,quan]=length(which(slope[inside_region[noNa]]>0))/length(noNa)
            }
        }
        print(percentage)
    }
}

#************************************** functions_support.r *************************************************************


# ---------------------- rq() function ------------------------------------------------------
function (formula, tau = 0.5, data, subset, weights, na.action, 
    method = "br", model = TRUE, contrasts = NULL, ...) 
{                                                                                                                                                                                   
    call <- match.call()                                                                                                                                                            
    mf <- match.call(expand.dots = FALSE)                                                                                                                                           
    m <- match(c("formula", "data", "subset", "weights", "na.action"),                                                                                                              
        names(mf), 0)                                                                                                                                                               
    mf <- mf[c(1, m)]                                                                                                                                                               
    mf$drop.unused.levels <- TRUE                                                                                                                                                   
    mf[[1]] <- as.name("model.frame")                                                                                                                                               
    mf <- eval.parent(mf)                                                                                                                                                           
    if (method == "model.frame")                                                                                                                                                    
        return(mf)                                                                                                                                                                  
    mt <- attr(mf, "terms")                                                                                                                                                         
    weights <- as.vector(model.weights(mf))                                                                                                                                         
    Y <- model.response(mf)                                                                                                                                                         
    X <- model.matrix(mt, mf, contrasts)                                                                                                                                            
    eps <- .Machine$double.eps^(2/3)                                                                                                                                                
    Rho <- function(u, tau) u * (tau - (u < 0))                                                                                                                                     
    if (length(tau) > 1) {                                                                                                                                                          
        if (any(tau < 0) || any(tau > 1))                                                                                                                                           
            stop("invalid tau:  taus should be >= 0 and <= 1")                                                                                                                      
        if (any(tau == 0))                                                                                                                                                          
            tau[tau == 0] <- eps
        if (any(tau == 1)) 
            tau[tau == 1] <- 1 - eps
        coef <- matrix(0, ncol(X), length(tau))
        rho <- rep(0, length(tau))
        fitted <- resid <- matrix(0, nrow(X), length(tau))
        for (i in 1:length(tau)) {
            z <- {
                if (length(weights)) 
                  rq.wfit(X, Y, tau = tau[i], weights, method, 
                    ...)
                else rq.fit(X, Y, tau = tau[i], method, ...)
            }
            coef[, i] <- z$coefficients
            resid[, i] <- z$residuals
            rho[i] <- sum(Rho(z$residuals, tau[i]))
            fitted[, i] <- Y - z$residuals
        }
        taulabs <- paste("tau=", format(round(tau, 3)))
        dimnames(coef) <- list(dimnames(X)[[2]], taulabs)
        dimnames(resid) <- list(dimnames(X)[[1]], taulabs)
        fit <- z
        fit$coefficients <- coef
        fit$residuals <- resid
        fit$fitted.values <- fitted
        if (method == "lasso") 
            class(fit) <- c("lassorqs", "rqs")
        else if (method == "scad") 
            class(fit) <- c("scadrqs", "rqs")
        else class(fit) <- "rqs"
    }
    else {
        process <- (tau < 0 || tau > 1)
        if (tau == 0) 
            tau <- eps
        if (tau == 1) 
            tau <- 1 - eps
        fit <- {
            if (length(weights)) 
                rq.wfit(X, Y, tau = tau, weights, method, ...)
            else rq.fit(X, Y, tau = tau, method, ...)
        }
        if (process) 
            rho <- list(x = fit$sol[1, ], y = fit$sol[3, ])
        else {
            dimnames(fit$residuals) <- list(dimnames(X)[[1]], 
                NULL)
            rho <- sum(Rho(fit$residuals, tau))
        }
        if (method == "lasso") 
            class(fit) <- c("lassorq", "rq")
        else if (method == "scad") 
            class(fit) <- c("scadrq", "rq")
        else class(fit) <- ifelse(process, "rq.process", "rq")
    }
    fit$na.action <- attr(mf, "na.action")
    fit$formula <- formula
    fit$terms <- mt
    fit$xlevels <- .getXlevels(mt, mf)
    fit$call <- call
    fit$tau <- tau
    fit$weights <- weights
    fit$residuals <- drop(fit$residuals)
    fit$rho <- rho
    fit$method <- method
    fit$fitted.values <- drop(fit$fitted.values)
    attr(fit, "na.message") <- attr(m, "na.message")
    if (model) 
        fit$model <- mf
    fit
}

# ------------------------------ rq.fit() ------------------------------------
Error in as.matrix(x) : 
  Fehler bei der Auswertung des Argumentes 'x' bei der Methodenauswahl
für Funktion 'as.matrix': Fehler: Argument "x" fehlt (ohne Standardwert)
> rq.fit
function (x, y, tau = 0.5, method = "br", ...) 
{
    fit <- switch(method, fn = rq.fit.fnb(x, y, tau = tau, ...), 
        fnb = rq.fit.fnb(x, y, tau = tau, ...), fnc = rq.fit.fnc(x, 
            y, tau = tau, ...), pfn = rq.fit.pfn(x, y, tau = tau, 
            ...), br = rq.fit.br(x, y, tau = tau, ...), lasso = rq.fit.lasso(x, 
            y, tau = tau, ...),   = rq.fit.scad(x, y, tau = tau, 
            ...), {
            what <- paste("rq.fit.", method, sep = "")
            if (exists(what, mode = "function")) (get(what, mode = "function"))(x, 
                y, ...) else stop(paste("unimplemented method:", 
                method))
        })
    fit$fitted.values <- y - fit$residuals
    fit$contrasts <- attr(x, "contrasts")
    fit
}

# ---------------------- rq() function ------------------------------------------------------
