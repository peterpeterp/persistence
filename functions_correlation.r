#!/home/pepflei/R/bin/Rscript

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



duration_correl <- function(trendID,states,toCor,toCor_name,toCor_short,toCor_shortZu,additional_style,
	toCor_startYear=1950,seasons=c("spring","summer","autumn","winter","year"),taus=c(0.25,0.5,0.75,0.9,0.95,0.98),stations=seq(1,1319,1),plot=FALSE){
	# toCor is array with dim=c(ntot,seasons,years) or dim=c(seasons,years) 
	# function will follow different procedures depending on dim(toCor)

	# correlation is calculated by determining quantiles and than cov/(var*var)
	# the qu() function had some problem, no idea why?

	# using toCor_startYear and toCor_years the matching years from markov are determined
	# from these years noNa are the years that will be plotted

	toCor_years=dim(toCor)[3]

	library(quantreg)
	ntot=1319
	correlation=array(NA,dim=c(ntot,4,states,length(taus)))

	# needed for the plots ------------------------------------
	if (plot==TRUE){
		pdf(file=paste("../plots/",trendID,"/",states,"_states",additional_style,"/stations/dur_",toCor_short,"_",toCor_shortZu,"_",stations,".pdf",sep=""))
		color=c("lightblue","blue","green","yellow","red","violet")
		if (states==2){
			state_names=c("cold","warm")
		}
		if (states==3){
			state_names=c("cold","normal","warm")
		}
	}
	# ---------------------------------------------------------

	for (sea in 1:length(seasons)){
		nc=open.ncdf(paste("../data/",trendID,"/",states,"_states",additional_style,"/duration/",trendID,"_duration_",states,"s_",seasons[sea],".nc",sep=""))
		duration=get.var.ncdf(nc,"dur")
		duration_mid=get.var.ncdf(nc,"dur_mid")

		x=seq(1,toCor_years,1)
		for (q in stations){
			if (length(dim(toCor))==2){
				toCor_loc=toCor
			}
			if (length(dim(toCor))==3){
				toCor_loc=toCor[q,,]
			}
			cat(paste("--",q))
			for (state in 1:states){
				size=length(which(duration_mid[q,state,]>toCor_startYear & duration_mid[q,state,]<2014))
				if (size>30 & length(which(!is.na(toCor_loc[sea,])))>10){
					durQu=array(NA,dim=c(toCor_years,length(taus)))
					toCor_ext=array(NA,dim=c(size))
					dur_ext=array(NA,size)
					count=0
					for (i in 1:(toCor_years-1)){
						if (count!=size){
							year=i+toCor_startYear-1
							# find durations in the selected year (and season)
							# determine quantile "postions" durQu
							inside=which(duration_mid[q,state,]>year & duration_mid[q,state,]<(year+1))
							if (length(inside)>0){
								durQu[i,]=quantile(duration[q,state,inside],taus,na.rm=TRUE)
								# this is just needed for plots
								toCor_ext[(count+1):(count+length(inside))]=toCor_loc[sea,i]
								dur_ext[(count+1):(count+length(inside))]=duration[q,state,inside]
								count=count+length(inside)
							}
						}
					}
					# calculate correlation between quantile "positions and toCor"
					for (qu in 1:length(taus)){
						if (sd(as.vector(toCor_loc[sea,]),na.rm=TRUE)!=0 & sd(as.vector(durQu[1:toCor_years,qu]),na.rm=TRUE)!=0){
							correlation[q,sea,state,qu]=cor(x=as.vector(toCor_loc[sea,]),y=as.vector(durQu[1:toCor_years,qu]),use="complete")
						}
					}
					# make little explanatory plot ---------------------------------
					if (plot==TRUE){
						noNa=which(!is.na(durQu[,4]))
						order=order(toCor_loc[sea,noNa])
						year=toCor_startYear:2014
						year=year[order]
						title=c()
						plot(as.vector(toCor_ext),dur_ext,xlab=toCor_name,ylab="duration length",ylim=c(-2,max(dur_ext,na.rm=TRUE)),
							main=paste(seasons[sea],state_names[state],"durations"))
						for (y in 1:toCor_years){
							abline(v=as.vector(toCor_loc[sea,order])[y],col="grey",lty=2)
							text(as.vector(toCor_loc[sea,order])[y],(-1.5+0.7*(-1)^y),label=year[y],srt=90,cex=0.5)
						}
						for (qu in 1:length(taus)){
							title[qu]=paste(taus[qu],"quantile cor -",toCor_name,"=",correlation[q,sea,state,qu])
							lines(as.vector(toCor_loc[sea,noNa][order]),as.vector(durQu[noNa,qu][order]),col=color[qu])
						}
						legend("topright",col=color,lty=array(1,length(taus)),legend=title)
					}
					# ------------------------------------------------------------
				}
			}
		}
	}
	if (plot!=TRUE){
		ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
		season <- dim.def.ncdf("season",units="uu",vals=1:5,unlim=FALSE)
		varstates <- dim.def.ncdf("states",units="uu",vals=1:states,unlim=FALSE)
		quantiles <- dim.def.ncdf("quantiles",units="0-1",vals=taus,unlim=FALSE)

		varCorrelation <- var.def.ncdf(name="correlation",longname=paste("correaltion between quantile values and",toCor_name,"values in season"),units="bla",dim=list(ID,season,varstates,quantiles), missval=-9999.0)
		vars=list(varCorrelation)
			 
		nc = create.ncdf(paste("../data/",trendID,"/",states,"_states",additional_style,"/correlations/",trendID,"_",toCor_short,"_",toCor_shortZu,"_duration_cor_",states,"states.nc",sep=""),vars)

		put.var.ncdf(nc,varCorrelation,correlation)
	}
	graphics.off()

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

dur_correlation_plot <- function(additional_style,trendID="91_5",states=2,toCor_short="nao",toCor_name="NAO",toCor_shortZu="",quA=0.95,
	state_names=c("cold","warm"),seasons=c("spring","summer","autumn","winter","year"),
	dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc"),worldmap = getMap(resolution = "low"),ntot=1319){

    nc=open.ncdf(paste("../data/",trendID,"/",states,"_states",additional_style,"/correlations/",trendID,"_",toCor_short,"_",toCor_shortZu,"_duration_cor_",states,"states.nc",sep=""))
    correlation=get.var.ncdf(nc,"correlation")
	reihen=array(NA,dim=c(8,ntot))
	titel=c()

	taus=get.var.ncdf(nc,"quantiles")
	qu=which(taus==quA)

	for (sea in 1:5){
		for (state in 1:states){
			reihen[((sea-1)*states+state),]=correlation[,sea,state,qu]
			titel[((sea-1)*states+state)]=paste("correlation between",toCor_name,"and",quA,"percentile of",state_names[state],"period duration in",seasons[sea])
		}
	}
    map_allgemein(dat,filename_plot=paste("../plots/",trendID,"/",states,"_states",additional_style,"/maps/dur_cor/",trendID,"_",toCor_short,"_",toCor_shortZu,"_duration_",quA,"_",states,"states.pdf",sep=""),
    	worldmap=worldmap,reihen=reihen,titel=titel,farb_mitte="0")
}



eke_mar_correl <- function(additional_style,trendID="91_5",states=2,level=1,stations=seq(1,1319,1),plot=FALSE,transition_names="cc wc cw ww"){
	nc=open.ncdf("../data/eke/eke_ID.nc")
	eke=get.var.ncdf(nc,"eke_sea")
	pressure=get.var.ncdf(nc,"levelist")

	markov_correl(additional_style=additional_style,trendID=trendID,states=states,stations=stations,plot=plot,toCor=eke[,level,,],toCor_name="eddy kynetic energy",toCor_short="eke",toCor_shortZu=paste(pressure[1],"mbar",sep=""),toCor_startYear=1979,transition_names=transition_names)
}

eke_dur_correl <- function(additional_style,trendID="91_5",states=2,level=1,stations=seq(1,1319,1),plot=FALSE){
	nc=open.ncdf("../data/eke/eke_ID.nc")
	eke=get.var.ncdf(nc,"eke_sea")
	pressure=get.var.ncdf(nc,"levelist")

	duration_correl(additional_style=additional_style,trendID=trendID,states=states,stations=stations,plot=plot,toCor=eke[,level,,],toCor_name="eddy kynetic energy",toCor_short="eke",toCor_shortZu=paste(pressure[1],"mbar",sep=""),toCor_startYear=1979)
}

index_mar_correl <- function(additional_style,trendID="91_5",states=2,toCor_name="NAO",toCor_short="nao",stations=seq(1,1319,1),plot=FALSE,transition_names="cc wc cw ww"){
	nc=open.ncdf(paste("../data/roh_data/",toCor_name,"_sea.nc",sep=""))
	index_sea=get.var.ncdf(nc,toCor_short)

	markov_correl(additional_style=additional_style,trendID=trendID,states=states,stations=stations,plot=plot,toCor=index_sea,toCor_name=toCor_name,toCor_short=toCor_short,toCor_shortZu="",toCor_startYear=1950,transition_names=transition_names)
}

index_dur_correl <- function(additional_style,trendID="91_5",states=2,toCor_name="NAO",toCor_short="nao",stations=seq(1,1319,1),plot=FALSE){
	nc=open.ncdf(paste("../data/roh_data/",toCor_name,"_sea.nc",sep=""))
	index_sea=get.var.ncdf(nc,toCor_short)

	duration_correl(additional_style=additional_style,trendID=trendID,states=states,stations=stations,plot=plot,toCor=index_sea,toCor_name=toCor_name,toCor_short=toCor_short,toCor_shortZu="",toCor_startYear=1950)
}



if (1==1){
	source("load.r")

	library(SDMTools)
	source("functions_regional.r")
	source("map_plot.r")
	library(rworldmap)
	library(fields)
	worldmap = getMap(resolution = "low")
	ntot=1319
	dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
}

#mar_correlation_plot(toCor_short="eke",toCor_name="EKE",toCor_shortZu="850mbar")
#mar_correlation_plot(toCor_short="mei",toCor_name="MEI",toCor_shortZu="")
#dur_correlation_plot(toCor_short="mei",toCor_name="MEI",toCor_shortZu="",quA=0.95)

if (1==2){
	source("load.r")

	#eke_mar_correl(level=1)
	indicesGr=c("NAO","AO","MEI","PNA","EAWR","PE")
	indicesKl=c("nao","ao","mei","pna","eawr","pe")
	for (i in 1:length(indicesGr)){
		indexGr=indicesGr[i]
		indexKl=indicesKl[i]
		#index_dur_correl(toCor_name=indexGr,toCor_short=indexKl)
		#index_mar_correl(toCor_name=indexGr,toCor_short=indexKl)

		#dur_correlation_plot(toCor_short=indexKl,toCor_name=indexGr,toCor_shortZu="",quA=0.95)
		mar_correlation_plot(toCor_short=indexKl,toCor_name=indexGr,toCor_shortZu="")

	}
}

#index_mar_correl(toCor_short="nao",toCor_name="NAO",plot=TRUE,stations=488)

#eke_mar_correl(additional_style="_mean_TX_not_random")
#mar_correlation_plot(additional_style="_mean_TX_not_random",toCor_short="eke",toCor_name="EKE",toCor_shortZu="850mbar")

eke_dur_correl(additional_style="_mean_TX_not_random")
dur_correlation_plot(additional_style="_mean_TX_not_random",toCor_short="eke",toCor_name="EKE",toCor_shortZu="",quA=0.95)
