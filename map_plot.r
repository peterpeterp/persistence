# teste teste
source("write.r")
source("load.r")
library(SDMTools)

location_finder <- function(station=0,lon=0,lat=0){
	dat=dat_load("../data/mid_lat.nc")
	library(rworldmap)
	library(fields)
	if (station!=0){
		q=which(dat$ID==station)
	}
	worldmap = getMap(resolution = "low")
	pdf(file="../plots/location.pdf")
	plot(worldmap,xlim=c(dat$lon[q]-10,dat$lon[q]+10),ylim=c(dat$lat[q]-10,dat$lat[q]+10))
	points(dat$lon[q],dat$lat[q],pch=15,col="red")
}

location_view <- function(station=0,lon=0,lat=0){
	dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2011.nc")
	library(rworldmap)
	library(fields)
	if (station!=0){
		q=which(dat$ID==station)
	}
	worldmap = getMap(resolution = "low")
	pdf(file="../plots/ID_region_map.pdf")
	plot(worldmap)#,xlim=c(-180,-5),ylim=c(35,60), asp = 3.5)
	regions_to_map()
	for (i in 1:length(dat$ID)){
		text(dat$lon[i],dat$lat[i],label=dat$ID[i],col="red",cex=0.125)
	}
}

map_allgemein <- function(dat,filename_plot,worldmap,ausschnitt,reihen,titel,farbe_mitte,reihen_sig=reihen*NA){
	#dat data form data_load()
	#filename_plot str - where to save plot
	#worldmap background of lon lat plot
	#ausschnitt c(lat_min,lat_max)
	#reihen array(... dim=c(anzahl der plots, anzahl der stationen))
	#titel liste von strings, plot-titles
	#farbe_mitte mid point of color range (white) at 0 for "0" or at the mean for "mean"

	jet.colors <- colorRampPalette( c( "violet","blue","white","yellow","red") )
	nbcol <- 100
	color <- jet.colors(nbcol)	

	pdf(file = filename_plot,width=12,height=8)
	par(mfrow=c(2,1))

	mid_lat = which(dat$lat >= ausschnitt[1] & dat$lat <= ausschnitt[2])
	for (i in 1:dim(reihen)[1]){
		print(titel[i])
		size=length(mid_lat)
		regio=array(NA,size)
		lon=array(NA,size)
		lat=array(NA,size)
		y1=array(NA,size)
		sig=array(NA,size)
		m=0

		for (k in 1:length(dat$ID)){
			if (k %in% mid_lat){
				if (is.na(reihen[i,k])==FALSE){
					m<-m+1
					lon[m]=dat$lon[k]
					lat[m]=dat$lat[k]
					y1[m]=reihen[i,k]
					if (reihen_sig[i,k]<0.05){
						sig[m]=4
					}					
				}
			}
		}
		y=array(NA,(m+2))
		notna=which(is.na(y1)==FALSE)
		for (n in notna){
			y[n+2]=y1[n]
		}
		if (farbe_mitte=="0"){
			aushol=max(c(abs(max(y,na.rm=TRUE)),abs(min(y,na.rm=TRUE))))
			y[1]=-aushol
			y[2]=aushol
		}
		if (farbe_mitte=="mean"){
			mi=mean(y,na.rm=TRUE)
			aushol=max(c(max(y,na.rm=TRUE)-mi,mi-min(y,na.rm=TRUE)))
			y[1]=mi-aushol
			y[2]=mi+aushol
		}			
		facetcol <- cut(y,nbcol)
		plot(worldmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5, main=titel[i])

		points(lon,lat,pch=15,col=color[facetcol[3:(size+2)]],cex=1.2)
		points(lon,lat,pch=sig)

		image.plot(legend.only=T, zlim=range(y), col=color)
	}
    graphics.off()
}

if (1==1){
	source("region_average.r")
	library(rworldmap)
	library(fields)
	worldmap = getMap(resolution = "low")
	nday = 91
	nyr = 5
	ntot=1319
	dat=dat_load("../data/dat_regional.nc",reg=1)

	if (1==1){
		# markov summer
		vars=c("mar_s_w_lr","mar_s_k_lr","mar_s_w_mk","mar_s_k_mk")
        vars_sig=c("mar_s_w_lr_sig","mar_s_k_lr_sig","mar_s_w_mk_sig","mar_s_k_mk_sig")
		nc=open.ncdf(sprintf("../data/%s_%s/%s_%s_markov_trend.nc",nday,nyr,nday,nyr))
		reihen=array(NA,dim=c(4,ntot))
		reihen_sig=array(NA,dim=c(4,ntot))
		titel=c()

		for (i in 1:4){
			reihen[i,]=get.var.ncdf(nc,vars[i])
			reihen_sig[i,]=get.var.ncdf(nc,vars_sig[i])
			for (k in 1:(length(nc$var))){
				if (nc$var[[k]]$name==vars[i]){
					titel[i]=nc$var[[k]]$longname
				}
			}	
		}
		map_allgemein(dat=dat,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farbe_mitte="0",
			filename_plot=sprintf("../plots/maps/%s_%s_markov_trend_summer.pdf",nday,nyr),
			worldmap=worldmap,ausschnitt=c(35,66))		
		map_regional(dat=dat,toPlot=reihen,titles=titel,filename_plot=sprintf("../plots/regions/%s_%s_markov_summer.pdf",nday,nyr))
	}

	if (1==2){
		# summer vergleich
		waka=c("warm","cold")
		titel_zusatz=c("mean","a","a_err","b","b_err","0.05 percentile","0.10 percentile")
		vars=c("dur_ana_warm_before","dur_ana_cold_before",
		    	"dur_ana_warm_after","dur_ana_cold_after")

		nc=open.ncdf(sprintf("../data/%s_%s/%s_%s_duration_analysis_summer.nc",nday,nyr,nday,nyr))

		reihen=array(NA,dim=c(4,ntot))
		titel=c()
		auswahl=c(1,6)
		for (i in 1:length(auswahl)){
		    for (j in 1:2){
		    	reihen[((i-1)*2+j),]=get.var.ncdf(nc,vars[j+2])[1:ntot,auswahl[i]]-get.var.ncdf(nc,vars[j])[1:ntot,auswahl[i]]
		    	titel[((i-1)*2+j)]=paste("difference in",waka[j],"period duration",titel_zusatz[auswahl[i]],"before and after 1980")
		    }
		}
		
		map_regional(dat=dat,toPlot=reihen,titles=titel,filename_plot=sprintf("../plots/regions/%s_%s_duration_summer.pdf",nday,nyr))
		map_allgemein(dat=dat,
			filename_plot=sprintf("../plots/maps/%s_%s_duration_summer_analysis.pdf",nday,nyr),
			worldmap=worldmap,ausschnitt=c(35,66),reihen=reihen,titel=titel,farbe_mitte="0")
	}	
	if (1==2){
	# summer all
	    titel_zusatz=c("mean","a","a_err","b","b_err","0.05 percentile","0.10 percentile")
	    vars=c("dur_ana_warm_full","dur_ana_cold_full")

	    nc=open.ncdf(sprintf("../data/%s_%s/%s_%s_duration_analysis_summer.nc",nday,nyr,nday,nyr))

	    reihen=array(NA,dim=c(14,ntot))
	    titel=c()
	    for (i in 1:7){
	    	for (j in 1:2){
	    		reihen[((i-1)*2+j),]=get.var.ncdf(nc,vars[j])[1:ntot,i]
	    		titel[((i-1)*2+j)]=paste(nc$var[[j]]$longname,titel_zusatz[i])
	    	}
	    }
		map_allgemein(dat=dat,
			filename_plot=sprintf("../plots/maps/%s_%s_duration_summer_diff_1980.pdf",nday,nyr),
			worldmap=worldmap,ausschnitt=c(35,66),reihen=reihen,titel=titel,farbe_mitte="mean")
	}        	        

	
}

