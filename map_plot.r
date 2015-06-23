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
		text(dat$lon[i],dat$lat[i],label=dat$ID[i],col="red",cex=0.25)
	}
}

trend_plot <- function(dat,filename_plot,newmap,ausschnitt,filename_markov=99,filename_shock=99){
	
	jet.colors <- colorRampPalette( c("blue","white", "red") )
	nbcol <- 100
	color <- jet.colors(nbcol)	

	pdf(file = filename_plot,width=12,height=8)
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

		}
	}
    graphics.off()

}

climatology_plot2 <- function(dat,str,filename,filename_plot,newmap,ausschnitt){
	
	jet.colors <- colorRampPalette( c( "violet","blue","green","yellow","red") )
	nbcol <- 100
	color <- jet.colors(nbcol)	

	pdf(file = filename_plot,width=12,height=8)
	par(mfrow=c(2,1))
		
	#str=c("markov_s_w","markov_s_k","markov_w_w","markov_w_k","markov_y_w","markov_y_k")
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
					if (opt[i,2]=="mean"){
						y1[m]=mean(tmp[k,],na.rm=TRUE)
					}
					
				}
			}
		}
		y=array(NA,(m))
		notna=which(is.na(mean)==FALSE)
		for (i in notna){
			y[i]=y1[i]
		}
		print(y)
		facetcol <- cut(y,nbcol)
		plot(newmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5,main=paste(var1$longname))
		points(lon,lat,pch=15,col=color[facetcol[3:(size+2)]],cex=1.2)
		image.plot(legend.only=T, zlim=range(y), col=color)

		y=array(NA,(m))
		notna=which(is.na(y2)==FALSE)
		for (i in notna){
			y[i]=y2[i]
		}
		print(y)
		facetcol <- cut(y,nbcol)
		plot(newmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5,main=paste(var1$longname))
		points(lon,lat,pch=15,col=color[facetcol[3:(size+2)]],cex=1.2)
		image.plot(legend.only=T, zlim=range(y), col=color)
	}
    graphics.off()
}

climatology_plot <- function(dat,opt,filenames,filename_plot,newmap,ausschnitt){
	
	jet.colors <- colorRampPalette( c( "violet","blue","green","yellow","red") )
	nbcol <- 100
	color <- jet.colors(nbcol)	

	pdf(file = filename_plot,width=12,height=8)
	par(mfrow=c(2,1))
		
	#str=c("markov_s_w","markov_s_k","markov_w_w","markov_w_k","markov_y_w","markov_y_k")
	nc1=open.ncdf(filenames[1])
	nc2=open.ncdf(filenames[2])


	mid_lat = which(dat$lat >= ausschnitt[1] & dat$lat <= ausschnitt[2])
	for (i in 1:2){
		if (opt[(i-1)*4+1]==1){
			var1<-nc1$var[[opt[(i-1)*4+4]]]
			tmp=get.var.ncdf(nc1,opt[(i-1)*4+2])
			#print(tmp)
		}
		if (opt[(i-1)*4+1]==2){
			var1<-nc2$var[[opt[(i-1)*4+4]]]
			tmp=get.var.ncdf(nc2,opt[(i-1)*4+2])
			#print(tmp)

		}
		size=length(mid_lat)
		lon=array(NA,size)
		lat=array(NA,size)
		y1=array(NA,size)
		y2=array(NA,size)
		m=0
		if (opt[(i-1)*4+3]=="mean"){
			for (k in 1:length(dat$ID)){
				if (k %in% mid_lat){
					if (is.na(mean(tmp[k,],na.rm=TRUE))==FALSE){
						m<-m+1
						lon[m]=dat$lon[k]
						lat[m]=dat$lat[k]
						y1[m]=mean(tmp[k,],na.rm=TRUE)
					}
				}
			}
		}
		if (opt[(i-1)*4+3]=="value"){
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
		}
		y=array(NA,(m))
		notna=which(is.na(mean)==FALSE)
		for (i in notna){
			y[i]=y1[i]
		}
		print(paste(var1$longname))
		print(y)
		facetcol <- cut(y,nbcol)
		plot(newmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5,main=paste(var1$longname))
		points(lon,lat,pch=15,col=color[facetcol[3:(size+2)]],cex=1.2)
		image.plot(legend.only=T, zlim=range(y), col=color)
		print("done")
		graphics.off()
		dfgdf
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
	        if (1==2){
				trend_plot(dat,sprintf("../plots/maps/%s_%s_markov_trend.pdf",nday,nyr),worldmap,c(35,66),
					filename_markov=sprintf("../data/%s_%s/%s_%s_markov_trend.nc",nday,nyr,nday,nyr))	        	
	        }
	        if (1==1){
	        	filenames=c(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr),sprintf("../data/%s_%s/%s_%s_duration.nc",nday,nyr,nday,nyr))
	        	opt=c(1,"markov_s_w","mean",2,2,"dur_mean_w","value",1)
				climatology_plot(dat=dat,opt=opt,filenames=filenames,
					filename_plot=sprintf("../plots/maps/%s_%s_markov_climatology.pdf",nday,nyr),
					worldmap,c(35,66))
   	
	        }
		}
	}

	
}

#location_view(572)