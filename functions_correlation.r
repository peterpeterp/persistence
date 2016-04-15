#!/home/pepflei/R/bin/Rscript

duration_correl <- function(toCor,toCor_name,toCor_short,toCor_shortZu,plot=FALSE,toCor_startYear=1950){
	# toCor is array with dim=c(ntot,seasons,years) or dim=c(seasons,years) 
	# function will follow different procedures depending on dim(toCor)

	# correlation is calculated by determining quantiles and than cov/(var*var)
	# the qu() function had some problem, no idea why?

	# using toCor_startYear and toCor_years the matching years from markov are determined
	# from these years noNa are the years that will be plotted

	toCor_years=dim(toCor)[length(dim(toCor))]

	correlation=array(NA,dim=c(5,ntot,2,(length(taus)+2),4))

	# needed for the plots ------------------------------------
	if (plot==TRUE){
		pdf(file=paste("../plots/",dataset,additional_style,"/",trendID,"/gridpoints/dur_",toCor_short,"_",toCor_shortZu,"_",ID_select,".pdf",sep=""),width=4,height=4)
	    par(mar=c(3, 3, 3, 3) + 0.1)
		plot(NA,xlim=c(1,10),ylim=c(1,10),frame.plot=FALSE,axes=FALSE,xlab="",ylab="")
		

		plot(NA,xlab="",ylab="",ylim=c(-10,50),main="",frame.plot=FALSE,axes=FALSE)
	    at_=axis(1,labels=FALSE,col="black")
	    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
	    axis(1,at=at_)

		plot(NA,xlab="",ylab="",ylim=c(-10,50),main="",frame.plot=FALSE,axes=FALSE)
	    at_=axis(2,labels=FALSE,col="black")
	    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
	    axis(2,at=at_)

		color=c("blue","green","red","violet",NA,"blue")
		state_names=c("cold","warm")
	}
	# ---------------------------------------------------------

	for (sea in 1:5){
		season<-season_names[sea]
        nc_dur=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_duration_",season,".nc",sep=""))
        dur=var.get.nc(nc_dur,"dur")
        dur_mid=var.get.nc(nc_dur,"dur_mid")

		x=seq(1,toCor_years,1)

        percentage<-0
        cat(paste("\n0 -> -> -> -> -> 100\n"))
        for (q in ID_select){
            if (q/ntot*100 > percentage){
                cat("-")
                percentage<-percentage+5
            }
			if (length(dim(toCor))==2){
				toCor_loc=toCor
			}
			if (length(dim(toCor))==3){
				toCor_loc=toCor[q,sea,]
				order=order(toCor_loc)
			}



			for (state in 1:2){
				size<-length(which(dur_mid[q,state,]>toCor_startYear & dur_mid[q,state,]<2014))
				if (size>30 & length(which(!is.na(toCor_loc)))>10){
					durQu<-array(NA,dim=c(toCor_years,(length(taus)+2)))
					toCor_ext<-array(NA,dim=c(size))
					dur_ext<-array(NA,size)
					count=0
					for (i in 1:(toCor_years-1)){
						if (count!=size){
							year=i+toCor_startYear-1
							# find durations in the selected year (and season)
							# determine quantile "postions" durQu
							inside=which(dur_mid[q,state,]>year & dur_mid[q,state,]<(year+1))
							if (length(inside)>0){
								durQu[i,1:length(taus)]=quantile(dur[q,state,inside],taus,na.rm=TRUE)
								durQu[i,(length(taus)+2)]=mean(dur[q,state,inside],na.rm=TRUE)
								# this is just needed for plots
								toCor_ext[(count+1):(count+length(inside))]=toCor_loc[i]
								dur_ext[(count+1):(count+length(inside))]=dur[q,state,inside]
								count<-count+length(inside)
							}
						}
					}
					# for regressions lm and rq order of the series makes no difference

					# quantile regression ------------------------------------------------------------
					#there was a really strange bug here! problem solved with dither
					sf=try(summary(rq(dither(as.vector(dur_ext),value=noise_level)~dither(as.vector(toCor_ext),value=noise_level),taus),se="boot"))
		               if (class(sf)!="try-error"){
		                slope=sapply(sf, function(x) c(tau=x$tau, x$coefficients[-1,]))
		                correlation[sea,q,state,1:length(taus),1]=slope[2,1:length(taus)]
		                correlation[sea,q,state,1:length(taus),2]=slope[5,1:length(taus)]
		                correlation[sea,q,state,1:length(taus),3]=sapply(sf, function(x) c(tau=x$tau, x$coefficients[1,1]))[2,]
		            }					

					# calculate correlation between mean duration and toCor
					if (sd(as.vector(toCor_loc),na.rm=TRUE)!=0 & sd(durQu[,(length(taus)+2)],na.rm=TRUE)!=0){
						lr<-summary(lm(as.vector(dur_ext)~dither(as.vector(toCor_ext),value=noise_level)))
						correlation[sea,q,state,(length(taus)+2),1]=lr$coefficients[2,1]
						correlation[sea,q,state,(length(taus)+2),2]=lr$coefficients[2,4]
						correlation[sea,q,state,(length(taus)+2),3]=lr$coefficients[1,1]
						correlation[sea,q,state,(length(taus)+2),4]=cor(x=as.vector(toCor_ext),y=as.vector(dur_ext),use="complete")
					}

					# make little explanatory plot ---------------------------------
					if (plot==TRUE){
						noNa=which(!is.na(durQu[,4]))
						order=order(toCor_loc[noNa])
						year=toCor_startYear:2014
						year=year[order]
						title=c()
						plot(NA,xlim=c(1,10),ylim=c(1,10),frame.plot=FALSE,axes=FALSE,xlab="",ylab="")
						text(5,5,paste(season_names[sea],state_names[state],"durations"))

						plot(as.vector(toCor_ext),dur_ext,xlab="",ylab="",ylim=c(-10,50),main="",pch=20,col="gray",frame.plot=TRUE,axes=FALSE)
						for (y in 1:toCor_years){
							abline(v=as.vector(toCor_loc[order])[y],col="grey",lty=2)
							text(as.vector(toCor_loc[order])[y],(-6+3*(-1)^y),label=year[y],srt=90,cex=0.5)
						}
						for (qu in c(3,6)){
							#title[qu]=paste(taus[qu],"quantile cor -",toCor_name,"=",correlation[sea,q,state,qu,1])
							#lines(as.vector(toCor_loc[noNa][order]),as.vector(durQu[noNa,qu][order]),col=color[qu])
							lines(as.vector(toCor_loc[noNa][order]),(as.vector(toCor_loc[noNa][order])*correlation[sea,q,state,qu,1]+correlation[sea,q,state,qu,3]),col=color[qu])

						}
						#legend("topright",col=color[c(3,6)],lty=array(1,2),legend=c("95th percentile","mean"))
					}
					# ------------------------------------------------------------
				}
			}
		}
	}
	graphics.off()

	if (plot!=TRUE){
		filename<-paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_eke_",toCor_shortZu,".nc",sep="") ; print(filename)
	    nc_out <- create.nc(filename)

	    print(dim(correlation))
	    dim.def.nc(nc_out,"seasons",dimlength=5,unlim=FALSE)
	    dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
	    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
	  
	    dim.def.nc(nc_out,"quants",dimlength=(length(taus)+2),unlim=FALSE)
	    dim.def.nc(nc_out,"outs",dimlength=4,unlim=FALSE)

	    var.def.nc(nc_out,"correlation","NC_DOUBLE",c(0,1,2,3,4))
	    att.put.nc(nc_out, "correlation", "missing_value", "NC_DOUBLE", -99999.9)
	    att.put.nc(nc_out, "correlation", "dim_explanation", "NC_CHAR", "season-ID-state-(taus,NA,mean)-outs")
	    att.put.nc(nc_out, "correlation", "explanation", "NC_CHAR", "1-correl slope, 2 correl-slope-sig, 3-intercept, 4-additional-correl")
	        
	    var.put.nc(nc_out,"correlation",correlation)      
	 
	    close.nc(nc_out) 
	}
	graphics.off()

}

dur_correlation_plot <- function(toCor_short="nao",toCor_name="NAO",toCor_shortZu="",val=1,val_zusatz="_mean",farb_mitte=farb_mitte){

    nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_eke_",toCor_shortZu,".nc",sep=""))
    correlation=var.get.nc(nc,"correlation")
	reihen=array(NA,dim=c(10,ntot))
	reihen_sig=array(NA,dim=c(10,ntot))
	titel=c()

	for (sea in 1:5){
		for (state in 1:2){
			reihen[((sea-1)*2+state),]=correlation[sea,,state,val,1]
			reihen_sig[((sea-1)*2+state),]=correlation[sea,,state,val,2]
		}
	}
	filename_plot<-paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_correl_",toCor_name,"_",toCor_shortZu,val_zusatz,".pdf",sep="")

	topo_map_plot(filename=filename_plot,reihen=reihen,reihen_sig=reihen_sig,farb_mitte=farb_mitte,farb_palette="gold-blau",titel=c(""))
}


eke_dur_correl <- function(level=1,plot=FALSE){
	nc=open.nc("../data/eke/eke_ID.nc")
	eke=var.get.nc(nc,"eke_sea")
	pressure=var.get.nc(nc,"levelist")

	duration_correl(plot=plot,toCor=eke[,level,,],toCor_name="eddy kynetic energy",toCor_short="eke",toCor_shortZu=paste(pressure[1],"mbar",sep=""),toCor_startYear=1979)
}


index_dur_correl <- function(dat,trendID="91_5",dataset="_TX",additional_style="",states=2,toCor_name="NAO",toCor_short="nao",stations=seq(1,1319,1),plot=FALSE){
	nc=open.ncdf(paste("../data/climatological_indices/",toCor_name,"_sea.nc",sep=""))
	index_sea=get.var.ncdf(nc,toCor_short)

	duration_correl(dat=dat,dataset=dataset,additional_style=additional_style,trendID=trendID,states=states,stations=stations,plot=plot,toCor=index_sea,toCor_name=toCor_name,toCor_short=toCor_short,toCor_shortZu="",toCor_startYear=1950)
}



correl_init <- function(){
	source("load.r")
	source("map_plot.r")
    library(quantreg)
    library(RNetCDF)
    library(SDMTools)
    library(fields)

	trendID<<-"91_7"
	dataset<<-"_TMean"
	additional_style<<-""
    #dat<<-dat_load(paste("../data/",dataset,"/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
	ntot<<-1319

    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    taus<<-c(0.75,0.9,0.95,0.99)
    ID_select<<-1:ntot
    noise_level<<-0.0000001
}


if (1==2){
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

correl_init()
ID_select<<-c(488,631,332)
eke_dur_correl(plot=TRUE)

adas
plot_init_Had_multiple()
dur_correlation_plot(toCor_short="eke",toCor_name="EKE",toCor_shortZu="850mbar",val=3,val_zusatz="_95",farb_mitte=c(-1,1))

