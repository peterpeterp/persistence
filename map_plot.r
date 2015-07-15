
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
	nbcol <- 101
	color <- jet.colors(nbcol)	

	pdf(file = filename_plot,width=12,height=8)
	par(mfrow=c(2,1))


	mid_lat = which(dat$lat >= ausschnitt[1] & dat$lat <= ausschnitt[2])
	if (farbe_mitte=="gemeinsam 0"){
		aushol=max(c(abs(max(reihen[1:dim(reihen)[1],mid_lat],na.rm=TRUE)),abs(min(reihen[1:dim(reihen)[1],mid_lat],na.rm=TRUE))))
	}
	if (farbe_mitte=="gemeinsam mean"){	
		mi=mean(reihen[1:dim(reihen)[1],mid_lat],na.rm=TRUE)
		aushol=max(c(abs(max(reihen[1:dim(reihen)[1],mid_lat],na.rm=TRUE))-mi,mi-abs(min(reihen[1:dim(reihen)[1],mid_lat],na.rm=TRUE))))
	}

	for (i in 1:dim(reihen)[1]){
		print(titel[i])
		size=length(mid_lat)
		regio=array(NA,size)
		lon=array(NA,size)
		lat=array(NA,size)
		y1=array(NA,size)
		sig=array(NA,size)
		nas=array(NA,size)
		m=0

		for (k in 1:length(dat$ID)){
			if (k %in% mid_lat){
				if (is.na(reihen[i,k])==FALSE){
					m<-m+1
					lon[m]=dat$lon[k]
					lat[m]=dat$lat[k]
					y1[m]=reihen[i,k]
					if (reihen_sig[i,k]<0.05 & !is.na(reihen_sig[i,k])){
						sig[m]=4
					}					
				}
			}
		}
		nas[which(dat$lat >= ausschnitt[1] & dat$lat <= ausschnitt[2] & is.na(reihen[i,]))]=13

		y=array(NA,(m+2))
		notna=which(is.na(y1)==FALSE)
		for (n in notna){
			y[n+2]=y1[n]
		}
		if (farbe_mitte=="gemeinsam 0"){
			y[1]=-aushol
			y[2]=aushol			
		}
		if (farbe_mitte=="gemeinsam mean"){
			y[1]=mi-aushol
			y[2]=mi+aushol			
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
		points(dat$lon,dat$lat,pch=nas)

		image.plot(legend.only=T, zlim=range(y), col=color)
	}
    graphics.off()
}



