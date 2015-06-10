# teste teste
source("write_load.r")




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
	warnings()
	plot(newmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5)
	points(lon,lat,pch=15,col=color[facetcol[3:(size+2)]],cex=1.2)
	points(lon,lat,pch=warm_sig,col="black",cex=1.2)
	image.plot(legend.only=T, zlim=range(y), col=color)

	y=array(NA,(size+2))
	y[3:(size+2)]=cold
	print(y[y>100])
	print(y[abs(y)>100])
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
	library(RNetCDF)
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
			per=per_load(sprintf("../data/%s_%s_per_shock_first_test.nc",nday,nyr))
			#map_plot(dat,per$year_warm_lr,per$year_cold_lr,per$year_warm_lr_sig,per$year_cold_lr_sig,sprintf("../plots/%s_%s_year_lr.pdf",nday,nyr),worldmap,c(35,66))
			#map_plot(dat,per$win_warm_lr,per$win_cold_lr,per$win_warm_lr_sig,per$win_cold_lr_sig,sprintf("../plots/%s_%s_winter_lr.pdf",nday,nyr),worldmap,c(35,66))
			#map_plot(dat,per$sum_warm_lr,per$sum_cold_lr,per$sum_warm_lr_sig,per$sum_cold_lr_sig,sprintf("../plots/%s_%s_summer_lr.pdf",nday,nyr),worldmap,c(35,66))
			#map_plot(dat,per$year_warm_mk,per$year_cold_mk,per$year_warm_mk_sig,per$year_cold_mk_sig,sprintf("../plots/%s_%s_year_mk.pdf",nday,nyr),worldmap,c(35,66))
			#map_plot(dat,per$win_warm_mk,per$win_cold_mk,per$win_warm_mk_sig,per$win_cold_mk_sig,sprintf("../plots/%s_%s_winter_mk.pdf",nday,nyr),worldmap,c(35,66))
			#map_plot(dat,per$sum_warm_mk,per$sum_cold_mk,per$sum_warm_mk_sig,per$sum_cold_mk_sig,sprintf("../plots/%s_%s_summer_mk.pdf",nday,nyr),worldmap,c(35,66))
			
			#map_plot(dat,per$shock_mk[,1],per$shock_mk[,2],per$shock_mk_sig[,1],per$shock_mk_sig[,2],sprintf("../plots/%s_%s_shock_mk.pdf",nday,nyr),worldmap,c(35,66))
			map_plot(dat,per$shock_mk[,3],per$shock_lr[,3],per$shock_mk_sig[,3],per$shock_lr_sig[,3],sprintf("../plots/%s_%s_shock_year.pdf",nday,nyr),worldmap,c(35,66))
		}
	}

	
}

