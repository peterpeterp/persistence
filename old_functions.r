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
