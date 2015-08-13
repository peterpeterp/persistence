
markov_correl <- function(trendID,states,toCor,toCor_name,toCor_short,toCor_shortZu,toCor_startYear=1950,transition_names,stations=seq(1,1319,1),plot=FALSE){
	ntot=1319
	transNumb=states*states
	toCor_years=dim(toCor)[length(dim(toCor))]

	nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/markov/",trendID,"_markov_",states,"states.nc",sep=""))
	markov=get.var.ncdf(nc,"markov")

	# needed for the plots
	if (plot==TRUE){
		pdf(file=paste("../plots/",trendID,"/",states,"_states/stations/mar_",toCor_short,"_",toCor_shortZu,"_",stations,".pdf",sep=""))
		seasons=c("spring","summer","autumn","winter")
		if (states==2){
			transLongName=c("cold to cold","warm to cold","cold to warm","warm to warm")
			color=c("blue","violet","green","red")
		}
		if (states==3){
			transLongName=c("cold to cold","normal to cold","warm to cold","cold to normal","normal to normal","warm to normal","cold to warm","normal to warm","warm to warm")
		}
	}	

	correlation=array(NA,dim=c(ntot,4,transNumb))
	corSlope=array(NA,dim=c(ntot,4,transNumb))
	corSlope_sig=array(NA,dim=c(ntot,4,transNumb))
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
		for (sea in 1:4){							
			if (plot==TRUE){
				plot(NA,ylim=c(0,1.1),xlim=c(min(toCor_loc[sea,],na.rm=TRUE),max(toCor_loc[sea,],na.rm=TRUE)),
					,xlab=toCor_name,ylab="transition probability",
					main=paste(seasons[sea]))
				title=c()

				noNa=which(!is.na(markov[q,sea,1,marSeq]))
				order=order(toCor_loc[sea,noNa])
				year=toCor_startYear:(toCor_startYear+toCor_years)
				year=year[order]
			}
			for (trans in 1:transNumb){
				if (length(which(!is.na(markov[q,sea,trans,marSeq])))>10 & length(which(!is.na(toCor_loc[sea,])))>10){

					if (sd(as.vector(toCor_loc[sea,]),na.rm=TRUE)!=0 & sd(as.vector(markov[q,sea,trans,marSeq]),na.rm=TRUE)!=0){
						correlation[q,sea,trans]=cor(x=as.vector(toCor_loc[sea,]),y=as.vector(markov[q,sea,trans,marSeq]),use="complete")
					}
					lr=summary(lm(toCor_loc[sea,]~markov[q,sea,trans,marSeq]))
					corSlope[q,sea,trans]=lr$coefficients[2,1]
					corSlope_sig[q,sea,trans]=lr$coefficients[2,4]

					if (plot==TRUE){
						lines(toCor_loc[sea,order],markov[q,sea,trans,(order+marSeq[1]-1)],col=color[trans])
						abline(lm(markov[q,sea,trans,marSeq]~toCor_loc[sea,]),col=color[trans],lty=2)
						title[trans]=paste(transLongName[trans],"transition - EKE correlation =",correlation[q,sea,trans])
					}
				}
			}
			if (plot==TRUE){
				for (y in 1:36){
					abline(v=as.vector(toCor_loc[sea,order])[y],col="grey",lty=2)
					text(as.vector(toCor_loc[sea,order])[y],(0.5+0.05*(-1)^y),label=year[y],srt=90,cex=0.5)
				}
				legend("topright",legend=title,lty=array(1,transNumb),col=color)
			}
		}
	}

	if (plot!=TRUE){
	    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
	    season <- dim.def.ncdf("season",units="uu",vals=1:4,unlim=FALSE)
	    transition <- dim.def.ncdf("transition",units=transition_names,vals=1:transNumb,unlim=FALSE)

	    varCorrelation <- var.def.ncdf(name="correlation",units="bla",dim=list(ID,season,transition), missval=-9999.0)
	    varCorSlope <- var.def.ncdf(name="corSlope",longname="linear regression like correlation",units="bla",dim=list(ID,season,transition), missval=-9999.0)
	    varCorSlope_sig <- var.def.ncdf(name="corSlope_sig",longname="linear regression like correlation_sig",units="bla",dim=list(ID,season,transition), missval=-9999.0)
	    
	    vars=list(varCorrelation,varCorSlope,varCorSlope_sig)  
	    nc = create.ncdf(paste("../data/",trendID,"/",states,"_states/",trendID,"_eke_markov_cor_",states,"states.nc",sep=""),vars)

		put.var.ncdf(nc,varCorrelation,correlation)
		put.var.ncdf(nc,varCorSlope,corSlope)
		put.var.ncdf(nc,varCorSlope_sig,corSlope_sig)

		close.ncdf(nc)
	}
}



duration_correl <- function(trendID,states,toCor,toCor_name,toCor_short,toCor_shortZu,toCor_startYear=1950,seasons=c("spring","summer","autumn","winter"),taus=c(0.25,0.5,0.75,0.9,0.95,0.98),stations=seq(1,1319,1),plot=FALSE){
	# toCor is array with dim=c(ntot,seasons,years)
	toCor_years=dim(toCor)[-1]

	library(quantreg)
	ntot=1319
	correlation=array(NA,dim=c(ntot,4,states,length(taus)))

	# needed for the plots
	if (plot==TRUE){
		pdf(file=paste("../plots/",trendID,"/",states,"_states/stations/dur_",toCor_short,"_",stations,".pdf",sep=""))
		color=c("lightblue","blue","green","yellow","red","violet")
		if (states==2){
			state_names=c("cold","warm")
		}
		if (states==3){
			state_names=c("cold","normal","warm")
		}
	}

	for (sea in 1:length(seasons)){
		nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/duration/",trendID,"_duration_",states,"s_",seasons[sea],".nc",sep=""))
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
							title[qu]=paste(taus[qu],"quantile cor =",correlation[q,sea,state,qu])
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
		season <- dim.def.ncdf("season",units="uu",vals=1:4,unlim=FALSE)
		varstates <- dim.def.ncdf("states",units="uu",vals=1:states,unlim=FALSE)
		quantiles <- dim.def.ncdf("quantiles",units="0-1",vals=taus,unlim=FALSE)

		varCorrelation <- var.def.ncdf(name="correlation",longname=paste("correaltion between quantile values and",toCor_name,"values in season"),units="bla",dim=list(ID,season,varstates,quantiles), missval=-9999.0)
		vars=list(varCorrelation)
			 
		nc = create.ncdf(paste("../data/",trendID,"/",states,"_states/",trendID,"_",toCor_short,"_",toCor_shortZu,"_duration_cor_",states,"states.nc",sep=""),vars)

		put.var.ncdf(nc,varCorrelation,correlation)
	}
	graphics.off()

}


source("load.r")

eke_mar_correl <- function(trendID="91_5",states=2){
	nc=open.ncdf("../data/eke/eke_ID.nc")
	eke=get.var.ncdf(nc,"eke_sea")

	markov_correl(trendID=trendID,states=states,stations=488,plot=TRUE,toCor=eke[,1,,],toCor_name="eddy kynetic energy",toCor_short="eke",toCor_shortZu="850mbar",toCor_startYear=1979)
}

eke_dur_correl <- function(trendID="91_5",states=2){
	nc=open.ncdf("../data/eke/eke_ID.nc")
	eke=get.var.ncdf(nc,"eke_sea")

	duration_correl(trendID=trendID,states=states,stations=488,plot=TRUE,toCor=eke[,1,,],toCor_name="eddy kynetic energy",toCor_short="eke",toCor_shortZu="850mbar",toCor_startYear=1979)
}


nao_dur_correl <- function(trendID="91_5",states=2){
	nc=open.ncdf("../data/roh_data/NAO_sea.nc")
	index_sea=get.var.ncdf(nc,"nao")
	print(index_sea)


	duration_correl(trendID=trendID,states=states,stations=488,plot=TRUE,toCor=index_sea,toCor_name="NOA",toCor_short="noa",toCor_shortZu="",toCor_startYear=1950)
}

eke_mar_correl()
#nao_dur_correl()