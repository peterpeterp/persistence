# teste teste
source("write_load.r")

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
	dat=dat_load("../data/mid_lat.nc")
	library(rworldmap)
	library(fields)
	if (station!=0){
		q=which(dat$ID==station)
	}
	worldmap = getMap(resolution = "low")
	pdf(file="../plots/location.pdf")
	plot(worldmap,xlim=c(-180,-5),ylim=c(35,60), asp = 3.5)
	for (i in 1:length(dat$ID)){
		text(dat$lon[i],dat$lat[i],label=dat$ID[i],col="red",cex=0.25)
	}
}

map_plot <- function(dat,warm_trend,cold_trend,warm_trend_sig,cold_trend_sig,filename,newmap,ausschnitt){
	mid_lat = which(dat$lat >= ausschnitt[1] & dat$lat <= ausschnitt[2] & ((is.na(warm_trend))==FALSE) & ((is.na(cold_trend))==FALSE))
	size=length(mid_lat)
	lon=array(NA,size)
	lat=array(NA,size)
	warm=array(NA,size)
	cold=array(NA,size)	
	warm_sig=array(NA,size)
	cold_sig=array(NA,size)	
	j=0
	for (i in 1:length(dat$ID)){
		#if (dat$ID[i] %in% mid_lat){
		if (i %in% mid_lat){
			j<-j+1
			lon[j]=dat$lon[i]
			lat[j]=dat$lat[i]
			warm[j]=warm_trend[i]
			cold[j]=cold_trend[i]
			if (warm_trend_sig[i]<0.05){
				warm_sig[j]=4
			}
			if (cold_trend_sig[i]<0.05){
				cold_sig[j]=4
			}
		}
	}

	jet.colors <- colorRampPalette( c("blue","white", "red") )
	nbcol <- 100
	color <- jet.colors(nbcol)

	pdf(file = filename,width=12,height=8)
	par(mfrow=c(2,1))

	y=array(NA,(size+2))
	y[3:(size+2)]=warm
	y[abs(y)>100*mean(y)]=NA


	aushol=max(c(abs(max(warm)),abs(min(warm))))
	y[1]=-aushol
	y[2]=aushol
	facetcol <- cut(y,nbcol)
	plot(newmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5)
	points(lon,lat,pch=15,col=color[facetcol[3:(size+2)]],cex=1.2)
	points(lon,lat,pch=warm_sig,col="black",cex=1.2)
	image.plot(legend.only=T, zlim=range(y), col=color)

	y=array(NA,(size+2))
	y[3:(size+2)]=cold
	y[abs(y)>100*mean(y)]=NA

	aushol=max(c(abs(max(cold)),abs(min(cold))))
	y[1]=-aushol
	y[2]=aushol
	facetcol <- cut(y,nbcol)
	plot(newmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5)
	points(lon,lat,pch=15,col=color[facetcol[3:(size+2)]],cex=1.2)
	points(lon,lat,pch=cold_sig,col="black",cex=1.2)
	image.plot(legend.only=T, zlim=range(y), col=color)

    graphics.off()

}


if (1==1){
	dyn.load("persistence_tools.so")
	library(rworldmap)
	library(fields)
	worldmap = getMap(resolution = "low")
	ndays = c(121,91,61,45)
	nyrs = c(7,5,3,1)

	ndays = c(91)
	nyrs = c(5)

	for (nday in ndays){
	    for (nyr in nyrs){
	        dat=dat_load("../data/mid_lat.nc")
			trend=trend_load(sprintf("../data/%s_%s_trend.nc",nday,nyr))
			per=per_load(sprintf("../data/%s_%s_per_shock_ma_3_neu.nc",nday,nyr))

			map_plot(dat,per$shock_mk[,2],per$shock_mk[,3],per$shock_mk_sig[,2],per$shock_mk_sig[,3],sprintf("../plots/%s_%s_ma_3_mk_neu.pdf",nday,nyr),worldmap,c(35,66))
			map_plot(dat,per$shock_mk[,1],per$shock_lr[,1],per$shock_mk_sig[,1],per$shock_lr_sig[,1],sprintf("../plots/%s_%s_ma_3_year_neu.pdf",nday,nyr),worldmap,c(35,66))			
			map_plot(dat,per$markov_mk[,3],per$markov_mk[,4],per$markov_mk_sig[,3],per$markov_mk_sig[,4],sprintf("../plots/%s_%s_markov_s_mk.pdf",nday,nyr),worldmap,c(35,66))
			map_plot(dat,per$markov_lr[,3],per$markov_lr[,4],per$markov_lr_sig[,3],per$markov_lr_sig[,4],sprintf("../plots/%s_%s_markov_s_lr.pdf",nday,nyr),worldmap,c(35,66))
			#map_plot(dat,per$markov_mk[,1],per$markov_lr[,3],per$markov_mk_sig[,3],per$markov_lr_sig[,3],sprintf("../plots/%s_%s_markov_year.pdf",nday,nyr),worldmap,c(35,66))
		}
	}

	
}

#location_view(572)