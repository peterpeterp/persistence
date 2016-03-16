
library(RNetCDF)

split_original <- function(i=1){
	nc<-open.nc("../data/raw_data/eobs-rr/rr_0.50deg_reg_v12.0.nc")
	rr<-var.get.nc(nc,"rr",start=c((58*(i-1)+1),1,1),count=c(58,101,23922))
	lon<-var.get.nc(nc,"longitude",start=c((58*(i-1)+1)),count=c(58))
	lat<-var.get.nc(nc,"latitude")
	time<-var.get.nc(nc,"time")
	print(dim(rr))

	close.nc(nc) 

	nc_out <- create.nc(paste("../data/raw_data/eobs-rr/rr_0.50deg_reg_v12.0_",i,".nc",sep=""))
	            
	dim.def.nc(nc_out,"lon",dimlength=58, unlim=FALSE)
	dim.def.nc(nc_out,"lat",dimlength=101,unlim=FALSE)
	dim.def.nc(nc_out,"time",dimlength=23922,unlim=FALSE)

	var.def.nc(nc_out,"lon","NC_FLOAT",c(0))
	#att.put.nc(nc_out, "lon", "missing_value", "NC_DOUBLE", 99999)

	var.def.nc(nc_out,"lat","NC_FLOAT",c(1))
	#att.put.nc(nc_out, "lat", "missing_value", "NC_DOUBLE", 99999)

	var.def.nc(nc_out,"time","NC_DOUBLE",c(2))
	#att.put.nc(nc_out, "time", "missing_value", "NC_DOUBLE", 99999)

	var.def.nc(nc_out,"rr","NC_SHORT",c(0,1,2))
	att.put.nc(nc_out, "rr", "missing_value", "NC_SHORT", -99)

	var.put.nc(nc_out,"lon",lon)              
	var.put.nc(nc_out,"lat",lat)              
	var.put.nc(nc_out,"time",time)              
	var.put.nc(nc_out,"rr",rr)              

	close.nc(nc_out) 
}

keep_not_empty <-function(i=1){
	nc<-open.nc(paste("../data/raw_data/eobs-rr/rr_0.50deg_reg_v12.0_",i,".nc",sep=""))
	rr<<-var.get.nc(nc,"rr")
	lon<<-var.get.nc(nc,"lon")
	lat<<-var.get.nc(nc,"lat")
	time_vec<<-var.get.nc(nc,"time")
	print(dim(rr))

	lon_neu<-array(NA,58*101)
	lat_neu<-array(NA,58*101)
	ID_neu<-array(NA,58*101)
	precip<-array(NA,c(58*101,23922))
	index<-0
	for (x in 1:58){
		for (y in 1:101){
			if (length(which(!is.na(rr[x,y,])))>0){
				index<-index+1
				lon_neu[index]=lon[x]
				lat_neu[index]=lat[y]
				ID_neu[index]=index
				precip[index,]=rr[x,y,]
			}
		}
	}
	nc_out <- create.nc(paste("../data/raw_data/eobs-rr/rr_0.50deg_reg_v12.0_reduced",i,".nc",sep=""))
	            
	dim.def.nc(nc_out,"ID",dimlength=index, unlim=FALSE)
	dim.def.nc(nc_out,"time",dimlength=23922,unlim=FALSE)

	var.def.nc(nc_out,"ID","NC_SHORT",c(0))

	var.def.nc(nc_out,"lon","NC_FLOAT",c(0))
	#att.put.nc(nc_out, "lon", "missing_value", "NC_DOUBLE", 99999)

	var.def.nc(nc_out,"lat","NC_FLOAT",c(0))
	#att.put.nc(nc_out, "lat", "missing_value", "NC_DOUBLE", 99999)

	var.def.nc(nc_out,"time","NC_DOUBLE",c(1))
	#att.put.nc(nc_out, "time", "missing_value", "NC_DOUBLE", 99999)

	var.def.nc(nc_out,"precip","NC_SHORT",c(0,1))
	att.put.nc(nc_out, "precip", "missing_value", "NC_SHORT", -99)

	var.put.nc(nc_out,"ID",ID_neu[1:index])    
	var.put.nc(nc_out,"lon",lon_neu[1:index])    
	var.put.nc(nc_out,"lat",lat_neu[1:index]) 
	var.put.nc(nc_out,"time",time_vec)      
	var.put.nc(nc_out,"precip",precip[1:index,])              
	close.nc(nc_out) 

}


merge_reduced <- function(files=1:4){

	pp<-array(NA,c(232*101,23922))
	ID<-array(NA,232*101)
	lon<-array(NA,232*101)
	lat<-array(NA,232*101)
	index<-0
	for (i in files){
		nc <- open.nc(paste("../data/raw_data/eobs-rr/rr_0.50deg_reg_v12.0_reduced",i,".nc",sep=""))
		tmp_rr<-var.get.nc(nc,"precip")
		tmp_ID<-var.get.nc(nc,"ID")
		tmp_lon<-var.get.nc(nc,"lon")
		tmp_lat<-var.get.nc(nc,"lat")
		time_vec<-var.get.nc(nc,"time")

		len<-length(tmp_ID)

		print(dim(pp[(index+1):(index+len),]))
		print(dim(tmp_rr))
		pp[(index+1):(index+len),]=tmp_rr
		ID[(index+1):(index+len)]=tmp_ID+index
		lon[(index+1):(index+len)]=tmp_lon
		lat[(index+1):(index+len)]=tmp_lat

		index<- index+len


	}

	nc_out <- create.nc(paste("../data/raw_data/eobs-rr/rr_0.50deg_reg_v12.0_merged.nc",sep=""))
	            
	dim.def.nc(nc_out,"ID",dimlength=index, unlim=FALSE)
	dim.def.nc(nc_out,"time",dimlength=23922,unlim=FALSE)

	var.def.nc(nc_out,"ID","NC_SHORT",c(0))
	var.def.nc(nc_out,"lon","NC_FLOAT",c(0))
	var.def.nc(nc_out,"lat","NC_FLOAT",c(0))
	var.def.nc(nc_out,"time","NC_DOUBLE",c(1))

	var.def.nc(nc_out,"pp","NC_SHORT",c(0,1))
	att.put.nc(nc_out, "pp", "missing_value", "NC_SHORT", -99)

	var.put.nc(nc_out,"ID",ID[1:index])    
	var.put.nc(nc_out,"lon",lon[1:index])    
	var.put.nc(nc_out,"lat",lat[1:index]) 
	var.put.nc(nc_out,"time",time_vec)      
	var.put.nc(nc_out,"pp",pp[1:index,])              
	close.nc(nc_out) 
}

put_in_dat_format <- function(){
	nc <- open.nc(paste("../data/raw_data/eobs-rr/rr_0.50deg_reg_v12.0_merged.nc",sep=""))
	pp<<-var.get.nc(nc,"pp")
	lon<<-var.get.nc(nc,"lon")
	lat<<-var.get.nc(nc,"lat")
	ID<<-var.get.nc(nc,"ID")
	ntot<-length(ID)

	pp_neu<-array(NA,c(ntot,365,66))
	for (yr in 1:65){
		pp_neu[,,yr]=pp[,(365*(yr-1)+1):(365*yr)]
	}	
	pp_neu[,1:197,66]=pp[,23726:23922]

	nc_out <- create.nc(paste("../data/raw_data/eobs-rr/rr_0.50deg_reg_v12.0_end.nc",sep=""))
	            
	dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
	dim.def.nc(nc_out,"days",dimlength=365,unlim=FALSE)
	dim.def.nc(nc_out,"years",dimlength=66,unlim=FALSE)

	var.def.nc(nc_out,"day","NC_SHORT",c(1))
	var.def.nc(nc_out,"year","NC_SHORT",c(2))
	var.def.nc(nc_out,"ID","NC_SHORT",c(0))
	var.def.nc(nc_out,"lon","NC_FLOAT",c(0))
	var.def.nc(nc_out,"lat","NC_FLOAT",c(0))


	var.def.nc(nc_out,"pp","NC_SHORT",c(0,1,2))
	att.put.nc(nc_out, "pp", "missing_value", "NC_SHORT", -99)

	var.put.nc(nc_out,"ID",ID)    
	var.put.nc(nc_out,"lon",lon)    
	var.put.nc(nc_out,"lat",lat) 
	var.put.nc(nc_out,"day",1:365)      
	var.put.nc(nc_out,"year",1950:2015)      
	var.put.nc(nc_out,"pp",pp_neu)              
	close.nc(nc_out) 
}



data_view <- function(){
	dat <<- dat_load_eobs(paste("../data/raw_data/eobs-rr/rr_0.50deg_reg_v12.0_1950-2015.nc",sep=""))
	pp<-dat$pp
	ntot<-length(dat$ID)


	reihen<-array(NA,c(3,ntot))
	for (q in 1:ntot){
		cat("-")
		reihen[1,q]=mean(pp[q,,],na.rm=TRUE)
		reihen[2,q]=sd(pp[q,,],na.rm=TRUE)
	}
	reihen[3,]=dat$ID

    topo_map_plot(filename_plot=paste("../plots/eobs/EU_PP_mean_sd.pdf",sep=""),reihen=reihen,farb_palette="regenbogen",pointsize=1.0,ausschnitt=c(10,80),xAusschnitt=c(-30,50),paper=c(5,5),ID_select=1:ntot)
}

location_view_eobs <- function(station=0,lon=0,lat=0,regions=NA){
	if (station!=0){
		q=which(dat$ID==station)
	}
	worldmap = getMap(resolution = "low")
	pdf(file="../plots/eobs/ID_region_map.pdf",width=8,height=8)
	plot(worldmap,xlim=c(-30,70),ylim=c(20,80), asp = 1.5)
	for (i in 1:length(dat$ID)){
		text(dat$lon[i],dat$lat[i],label=dat$ID[i],col="red",cex=0.05)
	}
    data(topoWorld)
	plot(topoWorld,xlim=c(-30,70),ylim=c(20,80), asp = 1.5,location="none")
	for (i in 1:length(dat$ID)){
		text(dat$lon[i],dat$lat[i],label=dat$ID[i],col="red",cex=0.05)
	}
	graphics.off()
}

pdf_view <- function(q){
	dat <<- dat_load_eobs(paste("../data/raw_data/eobs-rr/rr_0.50deg_reg_v12.0_1950-2015.nc",sep=""))
	pp<-dat$pp	
	histo(pp[q,,],plot=TRUE)
}

################################

#for (i in 1:4){
#	split_original(i)
#	keep_not_empty(i)
#}
#merge_reduced()
put_in_dat_format()


library(SDMTools)
library(fields)
source("load.r")
source("map_plot.r")

#data_view()

#location_view_eobs()