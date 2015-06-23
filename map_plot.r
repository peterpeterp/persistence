# teste teste
source("write.r")
source("load.r")

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
		text(dat$lon[i],dat$lat[i],label=dat$ID[i],col="green",cex=0.25)
	}
}

regions_to_map <- function(filename="../data/SREX_regions_all.csv"){
	srex <- read.csv(file=filename, header=TRUE, sep=",")
	latpos=c(3,5,7,9,11,13)
	lonpos=c(4,6,8,10,12,14)
	for (i in 1:30){
		if (srex[i,2]<27){
			lat=c()
			lon=c()
			k=0
			for (j in latpos){
				if (srex[i,j]!=9999){
					k=k+1
					lat[k]=srex[i,j]
				}
			}
			k=0
			for (j in lonpos){
				if (srex[i,j]!=9999){
					k=k+1
					lon[k]=srex[i,j]
				}
			}
			polygon(x=lon,y=lat,col=rgb(1,0,1,0.0),border="green")
		}
	}
}

trend_plot <- function(dat,filename_plot,newmap,ausschnitt,filename_markov=99,filename_shock=99){
	
	jet.colors <- colorRampPalette( c( "violet","blue","white","yellow","red") )
	nbcol <- 100
	color <- jet.colors(nbcol)	

	pdf(file = filename_plot,width=12,height=8)
	par(mfrow=c(2,1))

	if (filename_markov!=99 & filename_shock==99){
		str=c("mar_s_w_lr","mar_s_k_lr","mar_w_w_lr","mar_w_k_lr","mar_y_w_lr","mar_y_k_lr",
        "mar_s_w_mk","mar_s_k_mk","mar_w_w_mk","mar_w_k_mk","mar_y_w_mk","mar_y_k_mk")
        str_sig=c("mar_s_w_lr_sig","mar_s_k_lr_sig","mar_w_w_lr_sig","mar_w_k_lr_sig","mar_y_w_lr_sig","mar_y_k_lr_sig",
        	"mar_s_w_mk_sig","mar_s_k_mk_sig","mar_w_w_mk_sig","mar_w_k_mk_sig","mar_y_w_mk_sig","mar_y_k_mk_sig")
		nc=open.ncdf(filename_markov)
	}
	if (filename_shock!=99 & filename_markov==99){
    	str=c("sho_s_lr","sho_w_lr","sho_y_lr",
        	"sho_s_mk","sho_w_mk","sho_y_mk")
    	str_sig=c("sho_s_lr_sig","sho_w_lr_sig","sho_y_lr_sig",
			"sho_s_mk_sig","sho_w_mk_sig","sho_y_mk_sig")
		nc=open.ncdf(filename_shock)
	}
	for (i in 1:2){
		for (j in 1:(length(str)/2)){
			var1<-nc$var[[j+(i-1)*(length(str)/2)*2]]
			trend_tmp=get.var.ncdf(nc,str[j+(i-1)*(length(str)/2)])
			trend_sig_tmp=get.var.ncdf(nc,str_sig[j+(i-1)*(length(str)/2)])
			mid_lat = which(dat$lat >= ausschnitt[1] & dat$lat <= ausschnitt[2] & is.na(trend_tmp)==FALSE)
			size=length(mid_lat)
			lon=array(NA,size)
			lat=array(NA,size)
			trend=array(NA,size)
			trend_sig=array(NA,size)
			m=0
			for (k in 1:length(dat$ID)){
				if (k %in% mid_lat){
					m<-m+1
					lon[m]=dat$lon[k]
					lat[m]=dat$lat[k]
					trend[m]=trend_tmp[k]
					if (trend_sig_tmp[k]<0.05){
						trend_sig[m]=4
					}
				}
			}
			y=array(NA,(size+2))
			y[3:(size+2)]=trend
			y[abs(y)>100*mean(y)]=NA
			aushol=max(c(abs(max(trend)),abs(min(trend))))
			y[1]=-aushol
			y[2]=aushol
			facetcol <- cut(y,nbcol)
			plot(newmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5,main=var1$longname)
			points(lon,lat,pch=15,col=color[facetcol[3:(size+2)]],cex=1.2)
			points(lon,lat,pch=trend_sig,col="black",cex=1.2)
			image.plot(legend.only=T, zlim=range(y), col=color)
			regions_to_map()

		}
	}
    graphics.off()

}

climatology_markov <- function(dat,filename,filename_plot,newmap,ausschnitt){
	
	jet.colors <- colorRampPalette( c( "violet","blue","white","yellow","red") )
	nbcol <- 100
	color <- jet.colors(nbcol)	

	pdf(file = filename_plot,width=12,height=8)
	par(mfrow=c(2,1))
		
	str=c("markov_s_w","markov_s_k","markov_w_w","markov_w_k","markov_y_w","markov_y_k")
	nc=open.ncdf(filename)

	mid_lat = which(dat$lat >= ausschnitt[1] & dat$lat <= ausschnitt[2])
	for (i in 1:6){
		var1<-nc$var[[i+1]]
		tmp=get.var.ncdf(nc,str[i])
		size=length(mid_lat)
		lon=array(NA,size)
		lat=array(NA,size)
		y1=array(NA,size)
		y2=array(NA,size)
		m=0
		for (k in 1:length(dat$ID)){
			if (k %in% mid_lat){
				if (is.na(mean(tmp[k,],na.rm=TRUE))==FALSE){
					m<-m+1
					lon[m]=dat$lon[k]
					lat[m]=dat$lat[k]
					y1[m]=mean(tmp[k,],na.rm=TRUE)
					y2[m]=sd(tmp[k,],na.rm=TRUE)
				}
			}
		}
		y=array(NA,(m))
		notna=which(is.na(y1)==FALSE)
		for (i in notna){
			y[i]=y1[i]
		}
		facetcol <- cut(y,nbcol)
		plot(newmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5,main=paste(var1$longname,"mean"))
		regions_to_map()
		points(lon,lat,pch=15,col=color[facetcol[3:(size+2)]],cex=1.2)
		image.plot(legend.only=T, zlim=range(y), col=color)
		regions_to_map()

		y=array(NA,(m))
		notna=which(is.na(y2)==FALSE)
		for (i in notna){
			y[i]=y2[i]
		}
		facetcol <- cut(y,nbcol)
		plot(newmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5,main=paste(var1$longname,"standart deviation"))
		points(lon,lat,pch=15,col=color[facetcol[3:(size+2)]],cex=1.2)
		image.plot(legend.only=T, zlim=range(y), col=color)
		regions_to_map()

	}
    graphics.off()
}

climatology_duration <- function(dat,filename,filename_plot,newmap,ausschnitt){
	
	jet.colors <- colorRampPalette( c( "violet","blue","white","yellow","red") )
	nbcol <- 100
	color <- jet.colors(nbcol)	

	pdf(file = filename_plot,width=12,height=8)
	par(mfrow=c(2,1))
		
	str=c("dur_mean_w","dur_mean_c","dur_mean_d_w","dur_mean_d_c","dur_X_d_w","dur_X_d_c")
	nc=open.ncdf(filename)

	mid_lat = which(dat$lat >= ausschnitt[1] & dat$lat <= ausschnitt[2])
	for (i in 1:6){
		var1<-nc$var[[i]]
		tmp=get.var.ncdf(nc,str[i])
		size=length(mid_lat)
		lon=array(NA,size)
		lat=array(NA,size)
		y1=array(NA,size)
		y2=array(NA,size)
		m=0
		for (k in 1:length(dat$ID)){
			if (k %in% mid_lat){
				if (is.na(tmp[k])==FALSE){
					m<-m+1
					lon[m]=dat$lon[k]
					lat[m]=dat$lat[k]
					y1[m]=tmp[k]
				}
			}
		}
		y=array(NA,(m+2))
		notna=which(is.na(y1)==FALSE)
		for (i in notna){
			y[i+2]=y1[i]
		}
		aushol=max(c(abs(max(y,na.rm=TRUE)),abs(min(y,na.rm=TRUE))))
		y[1]=-aushol
		y[2]=aushol
		facetcol <- cut(y,nbcol)
		plot(newmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5,main=paste(var1$longname))
		points(lon,lat,pch=15,col=color[facetcol[3:(size+2)]],cex=1.2)
		image.plot(legend.only=T, zlim=range(y), col=color)
		regions_to_map()

	}
    graphics.off()
}



if (1==1){
	dyn.load("persistence_tools.so")
	library(rworldmap)
	library(fields)
	worldmap = getMap(resolution = "low")
	ndays = c(121,91,61)
	nyrs = c(7,5,3)

	ndays = c(91)
	nyrs = c(5)

	for (nday in ndays){
	    for (nyr in nyrs){
	        dat=dat_load("../data/mid_lat.nc")
	        if (1==1){
				trend_plot(dat,sprintf("../plots/maps/%s_%s_markov_trend.pdf",nday,nyr),worldmap,c(35,66),
					filename_markov=sprintf("../data/%s_%s/%s_%s_markov_trend.nc",nday,nyr,nday,nyr))	        	
	        }
	        if (1==1){
				climatology_markov(dat=dat,filename=sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr),
					filename_plot=sprintf("../plots/maps/%s_%s_markov_climatology.pdf",nday,nyr),
					worldmap,c(35,66))
   	
	        }
	        if (1==1){
				climatology_duration(dat=dat,filename=sprintf("../data/%s_%s/%s_%s_duration_summer.nc",nday,nyr,nday,nyr),
					filename_plot=sprintf("../plots/maps/%s_%s_duration_summer_climatology.pdf",nday,nyr),
					worldmap,c(35,66))
   	
	        }
	        if (1==1){
				climatology_duration(dat=dat,filename=sprintf("../data/%s_%s/%s_%s_duration_winter.nc",nday,nyr,nday,nyr),
					filename_plot=sprintf("../plots/maps/%s_%s_duration_winter_climatology.pdf",nday,nyr),
					worldmap,c(35,66))
   	
	        }
	        if (1==1){
				climatology_duration(dat=dat,filename=sprintf("../data/%s_%s/%s_%s_duration_year.nc",nday,nyr,nday,nyr),
					filename_plot=sprintf("../plots/maps/%s_%s_duration_year_climatology.pdf",nday,nyr),
					worldmap,c(35,66))
   	
	        }	        	        
		}
	}

	
}

