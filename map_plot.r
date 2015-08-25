
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
	dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
	library(rworldmap)
	library(fields)
	if (station!=0){
		q=which(dat$ID==station)
	}
	worldmap = getMap(resolution = "low")
	pdf(file="../plots/ID_region_map.pdf")
	plot(worldmap)#,xlim=c(-180,-5),ylim=c(35,60), asp = 3.5)
	for (i in 1:length(dat$ID)){
		text(dat$lon[i],dat$lat[i],label=dat$ID[i],col="red",cex=0.125)
	}
	region_names=c("srex","7rect","6wave","7wave","8wave")
	color=c("blue","green","red","orange","black")
    for (k in 1:length(region_names)){
    	add_region(region_names[k],color[k])
    }
}

add_region <- function(region_name,farbe){
    poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))
    for (i in 1:dim(poli)[1]){
        lon=poli[i,1:6]
        lat=poli[i,7:12]
        lon=lon[!is.na(lon)]
        lat=lat[!is.na(lat)]
        polygon(x=lon,y=lat,border=farbe)
    }
}

map_allgemein <- function(dat,filename_plot,worldmap,reihen,reihen_sig=reihen*NA,titel,
	farb_mitte="mean",farb_palette="regenbogen",region=NA,regionColor="black",grid=FALSE,ausschnitt=c(-80,80),col_row=c(1,1),paper=c(12,8)){
	#dat data form data_load()
	#filename_plot str - where to save plot
	#worldmap background of lon lat plot
	#ausschnitt c(lat_min,lat_max)
	#reihen array(... dim=c(anzahl der plots, anzahl der stationen))
	#titel liste von strings, plot-titles
	#farb_mitte mid point of color range (white) at 0 for "0" or at the mean for "mean"

	if (farb_palette=="gold-blau"){
		jet.colors <- colorRampPalette( c(rgb(0.2,0.6,0.6),rgb(0.5,1,1), rgb(0.98,0.98,0.98) ,rgb(1,1,0),rgb(0.6,0.6,0)))
	}
	if (farb_palette=="lila-gruen"){
		jet.colors <- colorRampPalette( c(rgb(0.2,0.6,0.2),rgb(0.5,1,0.5), rgb(0.98,0.98,0.98) ,rgb(1,0.5,1),rgb(0.6,0.2,0.6)))
	}
	if (farb_palette=="regenbogen"){
		jet.colors <- colorRampPalette( c( "blue","green","yellow","red") )
	}

	nbcol <- 101
	color <- jet.colors(nbcol)	

	pdf(file = filename_plot,width=paper[1],height=paper[2])
    par(mar=c(1,1,2,4))
	par(mfrow=col_row)


	mid_lat = which(dat$lat >= ausschnitt[1] & dat$lat <= ausschnitt[2])
	if (farb_mitte=="gemeinsam 0"){
		aushol=max(c(abs(max(reihen[1:dim(reihen)[1],mid_lat],na.rm=TRUE)),abs(min(reihen[1:dim(reihen)[1],mid_lat],na.rm=TRUE))))
	}
	if (farb_mitte=="gemeinsam mean"){	
		mi=mean(reihen[1:dim(reihen)[1],mid_lat],na.rm=TRUE)
		aushol=max(c(abs(max(reihen[1:dim(reihen)[1],mid_lat],na.rm=TRUE))-mi,mi-abs(min(reihen[1:dim(reihen)[1],mid_lat],na.rm=TRUE))))
	}

	subCount=0
	for (i in 1:dim(reihen)[1]){
		subCount=subCount+1
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

		# depending on color-cheme --------------------------------
		if (farb_mitte=="cut"){
			mi=mean(y,na.rm=TRUE)
			sd=sd(y,na.rm=TRUE)
			high=mi+1*sd
			low=mi-1*sd
			y[y>high]=high
			y[y<low]=low
			y[1]=low
			y[2]=high
		}

		if (farb_mitte=="gemeinsam 0"){
			y[1]=-aushol
			y[2]=aushol			
		}
		if (farb_mitte=="gemeinsam mean"){
			y[1]=mi-aushol
			y[2]=mi+aushol			
		}
		if (farb_mitte=="0"){
			aushol=max(c(abs(max(y,na.rm=TRUE)),abs(min(y,na.rm=TRUE))))
			y[1]=-aushol
			y[2]=aushol
		}
		if (farb_mitte=="mean"){
			mi=mean(y,na.rm=TRUE)
			y[1]=mi
			y[2]=mi
		}	
		if ((length(farb_mitte)==2 & farb_mitte[1]!=farb_mitte[2])){
			y[1]=farb_mitte[1]
			y[2]=farb_mitte[2]
			y[y>farb_mitte[2]]=farb_mitte[2]
			y[y<farb_mitte[1]]=farb_mitte[1]
		}	
		if ((length(farb_mitte)==2 & farb_mitte[1]==farb_mitte[2])){
			y[1]=farb_mitte[1]
			y[2]=farb_mitte[1]
			y[y>farb_mitte[1]]=farb_mitte[1]
			y[y<(-farb_mitte[1])]=(-farb_mitte[1])
		}	
		if (farb_mitte=="nichts"){
			y[1]=mean(y,na.rm=TRUE)
			y[2]=mean(y,na.rm=TRUE)
		}		
		# -----------------------------------------------------------------

		facetcol <- cut(y,nbcol)
		plot(worldmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5, main=titel[i])

		points(lon,lat,pch=15,col=color[facetcol[3:(size+2)]],cex=1.2)
		points(lon,lat,pch=sig)
		points(dat$lon,dat$lat,pch=nas)

		if (grid==TRUE){
			for (longi in seq(-180,180,30)){
				abline(v=longi,col="grey")
				text(longi,-80,label=longi)
			}
			for (lati in seq(-60,60,10)){
				abline(h=lati,col="grey")
				text(-190,lati,label=lati)
			}
		}
		if (!is.na(region)){
			add_region(region,regionColor)
		}
		if (subCount==(col_row[1]-1)){
			image.plot(legend.only=T,horizontal=TRUE, zlim=range(y), col=color,add=TRUE,legend.mar=5.1)
			subCount=0
			plot(NA,xlim=c(0,1),ylim=c(1,0),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
		}
		
	}
    graphics.off()
}

