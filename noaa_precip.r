
remap_check <- function(){
	#cdo -f nc sellonlatbox,230,305,20,50 -random,r720x360 grid_0p5_same_area.nc
	#cdo -f nc remapbil,grid_0p5_same_area.nc precip.V1.0.2003.nc precip.V1.0.2003_0p5.nc

	pdf("../plots/_noaa/remap_check.pdf")

	add_name<-c("2003","2003_0p5")
	pointsize<-c(0.3,0.6)
	for (i in 1:2){
		print(paste("../data/raw_data/precip/precip.V1.0.",add_name[i],".nc",sep=""))
		nc <- open.nc(paste("../data/raw_data/precip/precip.V1.0.",add_name[i],".nc",sep=""))
		pp<-var.get.nc(nc,"precip")
		lon<-var.get.nc(nc,"lon")
		lat<-var.get.nc(nc,"lat")

		
		XX<-c()
		YY<-c()
		PP<-c()
		day<-300
		nas<-c()
		XXXX<-c()
		YYYY<-c()
		index<-0
		nas_index<-0

		for (x in 1:length(lon)){
			for (y in 1:length(lat)){
				if (!is.na(pp[x,y,day])){
					index<-index+1
					XX[index]=lon[x]
					YY[index]=lat[y]
					PP[index]=pp[x,y,day]
				}
				if (is.na(pp[x,y,day])){
					nas_index<-nas_index+1
					nas[nas_index]=4
					XXXX[nas_index]=lon[x]
					YYYY[nas_index]=lat[y]
				}
			}
		}
		plot(topoWorld,xlim=c(230,305),ylim=c(20,50), asp = 1,location="none")
		PP[PP>40]=40
		PP[index+1:index+2]=c(0,40)
		jet.colors <- colorRampPalette( c( "orange","yellow","green","blue") )
		color <- jet.colors(101)	
		facetcol <- cut(PP,101)
		points(XX,YY,col=color[facetcol[1:index]],pch=15,cex=pointsize[i])
		image.plot(legend.only=T,horizontal=TRUE, zlim=range(PP[1:index]), col=color,add=TRUE,fill=TRUE)

		#show nas
		plot(topoWorld,xlim=c(230,305),ylim=c(20,50), asp = 1,location="none")
		points(XXXX,YYYY,pch=15,cex=pointsize[i])
	}
	graphics.off()
	
}

merge_files <- function(){

	#nas raussuchen, lon lat klaeren etc
	nc <- open.nc(paste("../data/raw_data/precip/precip.V1.0.2003_0p5.nc",sep=""))
	pp<-var.get.nc(nc,"precip")
	lon<-var.get.nc(nc,"lon")
	lat<-var.get.nc(nc,"lat")

	lon_merged<-c()
	lat_merged<-c()
	XX<-c()
	YY<-c()
	ntot<-0
	day<-300
	for (x in 1:length(lon)){
		for (y in 1:length(lat)){
			if (!is.na(pp[x,y,day])){
				ntot<-ntot+1
				lon_merged[ntot]=lon[x]
				lat_merged[ntot]=lat[y]
				XX[ntot]=x
				YY[ntot]=y
			}
		}
	}	

	pp_merged<-array(NA,c(ntot,365,59))
	for (yr in 1:59){
		print(paste("../data/raw_data/precip/precip.V1.0.",(yr+1947),"_0p5.nc",sep=""))
		nc <- open.nc(paste("../data/raw_data/precip/precip.V1.0.",(yr+1947),"_0p5.nc",sep=""))
		pp<-var.get.nc(nc,"precip")
		for (q in 1:ntot){
			pp_merged[q,,yr]=pp[XX[q],YY[q],1:365]
		}
	}

	print(paste("../data/raw_data/precip/precip.V1.0.1948-2006_0p5.nc",sep=""))
	nc_out <- create.nc(paste("../data/raw_data/precip/precip.V1.0.1948-2006_0p5.nc",sep=""))
	            
	dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
	dim.def.nc(nc_out,"days",dimlength=365,unlim=FALSE)
	dim.def.nc(nc_out,"years",dimlength=59,unlim=FALSE)

	var.def.nc(nc_out,"day","NC_SHORT",c(1))
	var.def.nc(nc_out,"year","NC_SHORT",c(2))
	var.def.nc(nc_out,"ID","NC_SHORT",c(0))
	var.def.nc(nc_out,"lon","NC_FLOAT",c(0))
	var.def.nc(nc_out,"lat","NC_FLOAT",c(0))


	var.def.nc(nc_out,"pp","NC_DOUBLE",c(0,1,2))
	att.put.nc(nc_out, "pp", "missing_value", "NC_DOUBLE", -9999)

	var.put.nc(nc_out,"ID",1:ntot)    
	var.put.nc(nc_out,"lon",lon_merged)    
	var.put.nc(nc_out,"lat",lat_merged) 
	var.put.nc(nc_out,"day",1:365)      
	var.put.nc(nc_out,"year",1:59)      
	var.put.nc(nc_out,"pp",pp_merged)              
	close.nc(nc_out) 
}

check_merged <- function(){
	nc <- open.nc(paste("../data/raw_data/precip/precip.V1.0.1948-2006_0p5.nc",sep=""))
	pp<-var.get.nc(nc,"pp")
	lon<-var.get.nc(nc,"lon")
	lat<-var.get.nc(nc,"lat")	
	ntot<-length(lat)

	pdf("../plots/_noaa/merge_check.pdf")
	
	plot(topoWorld,xlim=c(230,305),ylim=c(20,50), asp = 1,location="none")

	PP<-pp[,300,56]

	PP[PP>40]=40
	PP<-c(PP,c(0,40))
	jet.colors <- colorRampPalette( c( "orange","yellow","green","blue") )
	color <- jet.colors(101)	
	facetcol <- cut(PP,101)
	points(lon,lat,col=color[facetcol[1:ntot]],pch=15,cex=0.6)
	image.plot(legend.only=T,horizontal=TRUE, zlim=range(PP[1:ntot]), col=color,add=TRUE,fill=TRUE)
	graphics.off()
}


library(RNetCDF)
library(SDMTools)
library(fields)
source("load.r")
source("map_plot.r")

#remap_check()
#merge_files()
#check_merged()