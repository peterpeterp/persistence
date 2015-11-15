#!/home/pepflei/R/bin/Rscript

duration_correl <- function(dat,trendID,states,toCor,toCor_name,toCor_short,toCor_shortZu,dataset,additional_style,
	toCor_startYear=1950,season_names=c("MAM","JJA","SON","DJF","year"),taus=c(0.25,0.5,0.75,0.9,0.95,0.98),stations=seq(1,1319,1),plot=FALSE){
	# toCor is array with dim=c(ntot,seasons,years) or dim=c(seasons,years) 
	# function will follow different procedures depending on dim(toCor)

	# correlation is calculated by determining quantiles and than cov/(var*var)
	# the qu() function had some problem, no idea why?

	# using toCor_startYear and toCor_years the matching years from markov are determined
	# from these years noNa are the years that will be plotted

	toCor_years=dim(toCor)[length(dim(toCor))]

	library(quantreg)
	ntot=1319
	correlation=array(NA,dim=c(ntot,5,2,7,3))

	# needed for the plots ------------------------------------
	if (plot==TRUE){
		pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/stations/dur_",toCor_short,"_",toCor_shortZu,"_",stations,".pdf",sep=""))
		color=c("lightblue","blue","green","yellow","red","violet")
		state_names=c("cold","warm")
	}
	# ---------------------------------------------------------

	for (sea in 1:length(season_names)){
		print(season_names[sea])
        nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",trendID,dataset,"_duration_",season_names[sea],".nc",sep=""))
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
			cat(paste("-",q))
			for (state in 1:states){
				size=length(which(duration_mid[q,state,]>toCor_startYear & duration_mid[q,state,]<2014))
				if (size>30 & length(which(!is.na(toCor_loc[sea,])))>10){
					durQu=array(NA,dim=c(toCor_years,length(taus)))
					durMean=array(NA,dim=c(toCor_years))
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
								durMean[i]=mean(duration[q,state,inside],na.rm=TRUE)
								# this is just needed for plots
								toCor_ext[(count+1):(count+length(inside))]=toCor_loc[sea,i]
								dur_ext[(count+1):(count+length(inside))]=duration[q,state,inside]
								count=count+length(inside)
							}
						}
					}
					# calculate correlation between quantile "positions and toCor"

					if (1==2){
						for (qu in 1:length(taus)){
							if (sd(as.vector(toCor_loc[sea,]),na.rm=TRUE)!=0 & sd(as.vector(durQu[1:toCor_years,qu]),na.rm=TRUE)!=0){
								correlation[q,sea,state,qu,3]=cor(x=as.vector(toCor_loc[sea,]),y=as.vector(durQu[1:toCor_years,qu]),use="complete")
							}
						}
					}
					# quantile regression ------------------------------------------------------------
					#there is a really strange bug here! unclear what the problem is. maybe if there are to many points associated to the same point.....
					if (1==2){
						sf=try(summary(rq(as.vector(dur_ext)~as.vector(toCor_ext),taus),se="nid"))
		                if (class(sf)=="try-error"){
		                    for (qu in 1:length(taus)){
		                        sf=try(summary(rq(as.vector(dur_ext)~as.vector(toCor_ext),taus[qu]),se="nid"))
		                        if (class(sf)=="try-error"){
		                            correlation[q,sea,state,qu,1:2]=array(NA,2)
		                        }
		                        else {
		                            correlation[q,sea,state,qu,1]=sf$coefficients[2,1]
		                            correlation[q,sea,state,qu,2]=sf$coefficients[2,4]
		                        }
		                    }
		                }
		                else {
		                    slope=sapply(sf, function(x) c(tau=x$tau, x$coefficients[-1,]))
		                    correlation[q,sea,state,1:length(taus),1]=slope[2,1:length(taus)]
		                    correlation[q,sea,state,1:length(taus),2]=slope[5,1:length(taus)]
		                }
		            }
	               	# quantile regression ------------------------------------------------------------


					# home made quantile regression ------------------------------------------------------------
					#there is a really strange bug here! unclear what the problem is. maybe if there are to many points associated to the same point.....
					if (1==1){
						for (qu in 1:length(taus)){
							if (sd(as.vector(toCor_loc[sea,]),na.rm=TRUE)!=0 & sd(as.vector(durQu[1:toCor_years,qu]),na.rm=TRUE)!=0){
								lr=summary(lm(as.vector(durQu[1:toCor_years,qu])~as.vector(toCor_loc[sea,])))
								correlation[q,sea,state,qu,1]=lr$coefficients[2,1]
								correlation[q,sea,state,qu,2]=lr$coefficients[2,4]
								correlation[q,sea,state,qu,3]=cor(x=as.vector(toCor_loc[sea,]),y=as.vector(durMean[1:toCor_years]),use="complete")
							}
						}
		            }
	               	# quantile regression ------------------------------------------------------------


					# calculate correlation between mean duration and toCor
					if (sd(as.vector(toCor_loc[sea,]),na.rm=TRUE)!=0 & sd(as.vector(durMean[1:toCor_years]),na.rm=TRUE)!=0){
						lr=summary(lm(durMean[1:toCor_years]~toCor_loc[sea,]))
						correlation[q,sea,state,7,1]=lr$coefficients[2,1]
						correlation[q,sea,state,7,2]=lr$coefficients[2,4]
						correlation[q,sea,state,7,3]=cor(x=as.vector(toCor_loc[sea,]),y=as.vector(durMean[1:toCor_years]),use="complete")
					}

					# make little explanatory plot ---------------------------------
					if (plot==TRUE){
						noNa=which(!is.na(durQu[,4]))
						order=order(toCor_loc[sea,noNa])
						year=toCor_startYear:2014
						year=year[order]
						title=c()
						plot(as.vector(toCor_ext),dur_ext,xlab=toCor_name,ylab="duration length",ylim=c(-2,max(dur_ext,na.rm=TRUE)),
							main=paste(season_names[sea],state_names[state],"durations"))
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
		varstates <- dim.def.ncdf("states",units="uu",vals=1:2,unlim=FALSE)
		quantiles <- dim.def.ncdf("quantiles",units="quantiles: 0.25, 0.5, 0.75, 0.9, 0.95, 0.98 - mean",vals=c(taus,7),unlim=FALSE)
		outs <- dim.def.ncdf("outs",units="reg_value - sig - cor()",vals=1:3,unlim=FALSE)

		varCorrelation <- var.def.ncdf(name="correlation",longname=paste("correaltion between quantile values and",toCor_name,"values in season"),units="bla",dim=list(ID,season,varstates,quantiles,outs), missval=-9999.0)
		vars=list(varCorrelation)
			 
		nc = create.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/correlations/",trendID,"_",toCor_short,"_",toCor_shortZu,"_duration_cor_",states,"states.nc",sep=""),vars)
		put.var.ncdf(nc,varCorrelation,correlation)
	}
	graphics.off()

}

dur_correlation_plot <- function(dat,trendID="91_5",dataset="_TX",additional_style="",states=2,toCor_short="nao",toCor_name="NAO",toCor_shortZu="",quA=0.95,state_names=c("cold","warm"),seasons=c("spring","summer","autumn","winter","year"),worldmap = getMap(resolution = "low"),ntot=1319,farb_mitte=farb_mitte){

    nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/correlations/",trendID,"_",toCor_short,"_",toCor_shortZu,"_duration_cor_",states,"states.nc",sep=""))
    correlation=get.var.ncdf(nc,"correlation")
	reihen=array(NA,dim=c(10,ntot))
	reihen_sig=array(NA,dim=c(10,ntot))
	titel=c()

	taus=get.var.ncdf(nc,"quantiles")
	qu=which(taus==quA)

	for (sea in 1:5){
		for (state in 1:states){
			reihen[((sea-1)*states+state),]=correlation[,sea,state,qu,1]
			reihen_sig[((sea-1)*states+state),]=correlation[,sea,state,qu,2]
			titel[((sea-1)*states+state)]=paste("correlation between",toCor_name,"and",quA,"percentile of",state_names[state],"period duration in",seasons[sea])
		}
	}
	#print(paste("../plots/",trendID,"/",dataset,additional_style,"/maps/dur_cor/",trendID,"_",toCor_short,"_",toCor_shortZu,"_duration_",quA,"_",states,"states.pdf",sep=""))
    map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/dur_cor/",trendID,"_",toCor_short,"_",toCor_shortZu,"_duration_",quA,"_",states,"states.pdf",sep=""),
    	worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette="gold-blau")

    # mean duration length
	reihen=array(NA,dim=c(10,ntot))
	reihen_sig=array(NA,dim=c(10,ntot))
	titel=c()

	for (sea in 1:5){
		for (state in 1:states){
			reihen[((sea-1)*states+state),]=correlation[,sea,state,7,1]
			reihen_sig[((sea-1)*states+state),]=correlation[,sea,state,7,2]
			titel[((sea-1)*states+state)]=paste("correlation between",toCor_name,"and","mean",state_names[state],"period duration in",seasons[sea])
		}
	}
    map_allgemein(dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/maps/dur_cor/",trendID,"_",toCor_short,"_",toCor_shortZu,"_duration_","mean","_",states,"states.pdf",sep=""),
    	worldmap=worldmap,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=farb_mitte,farb_palette="gold-blau")
}


eke_dur_correl <- function(dat,trendID="91_5",dataset="_TX",additional_style="",states=2,level=1,stations=seq(1,1319,1),plot=FALSE){
	nc=open.ncdf("../data/eke/eke_ID.nc")
	eke=get.var.ncdf(nc,"eke_sea")
	pressure=get.var.ncdf(nc,"levelist")

	duration_correl(dat=dat,dataset=dataset,additional_style=additional_style,trendID=trendID,states=states,stations=stations,plot=plot,toCor=eke[,level,,],toCor_name="eddy kynetic energy",toCor_short="eke",toCor_shortZu=paste(pressure[1],"mbar",sep=""),toCor_startYear=1979)
}


index_dur_correl <- function(dat,trendID="91_5",dataset="_TX",additional_style="",states=2,toCor_name="NAO",toCor_short="nao",stations=seq(1,1319,1),plot=FALSE){
	nc=open.ncdf(paste("../data/climatological_indices/",toCor_name,"_sea.nc",sep=""))
	index_sea=get.var.ncdf(nc,toCor_short)

	duration_correl(dat=dat,dataset=dataset,additional_style=additional_style,trendID=trendID,states=states,stations=stations,plot=plot,toCor=index_sea,toCor_name=toCor_name,toCor_short=toCor_short,toCor_shortZu="",toCor_startYear=1950)
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
	dataset="_TMean"
	additional_style=""
	dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
}


if (1==1){
	source("load.r")

	#eke_mar_correl(level=1)
	indicesGr=c("NAO","AO","MEI","PNA","EAWR","PE")
	indicesKl=c("nao","ao","mei","pna","eawr","pe")
	for (i in 1:length(indicesGr)){
		indexGr=indicesGr[i]
		indexKl=indicesKl[i]

		index_dur_correl(dat,dataset=dataset,additional_style=additional_style,toCor_name=indexGr,toCor_short=indexKl)
		dur_correlation_plot(dat,dataset=dataset,additional_style=additional_style,toCor_short=indexKl,toCor_name=indexGr,toCor_shortZu="",quA=0.95,farb_mitte="0")
	}
}


eke_dur_correl(dat,dataset=dataset,additional_style=additional_style)
dur_correlation_plot(dat,dataset=dataset,additional_style=additional_style,toCor_short="eke",toCor_name="EKE",toCor_shortZu="850mbar",quA=0.95,farb_mitte=c(-2,2))

