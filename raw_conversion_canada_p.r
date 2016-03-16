

merge_years <- function(){
	ntot=2182
	pp_merged=array(NA,c(ntot,365,63))
	for (yr in 1950:2012){
		cat(paste(yr,"-"))
		print(dim(read.table(paste("../data/raw_data/canda_P/P_canada",yr,".txt",sep=""))))
		pp_merged[,,(yr-1949)]=read.table(paste("../data/raw_data/canda_P/P_canada",yr,".txt",sep=""))[,c(6:64,66:371)]
	}

	example_data=read.table(paste("../data/raw_data/canda_P/Prnational","2000",".txt",sep=""))

	nc_out=create.nc("../data/canda_P/canada_p.nc")

	dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
	dim.def.nc(nc_out,"days",dimlength=365,unlim=FALSE)
	dim.def.nc(nc_out,"years",dimlength=63,unlim=FALSE)

	var.def.nc(nc_out,"day","NC_SHORT",c(1))
	var.def.nc(nc_out,"year","NC_SHORT",c(2))
	var.def.nc(nc_out,"ID","NC_SHORT",c(0))
	var.def.nc(nc_out,"lon","NC_FLOAT",c(0))
	var.def.nc(nc_out,"lat","NC_FLOAT",c(0))


	var.def.nc(nc_out,"pp","NC_DOUBLE",c(0,1,2))
	att.put.nc(nc_out, "pp", "missing_value", "NC_DOUBLE", -9999)

	var.put.nc(nc_out,"ID",1:ntot)    
	var.put.nc(nc_out,"lon",example_data[,3])    
	var.put.nc(nc_out,"lat",example_data[,4]) 
	var.put.nc(nc_out,"day",1:365)      
	var.put.nc(nc_out,"year",1950:2012)      
	var.put.nc(nc_out,"pp",pp_merged)              
	close.nc(nc_out) 	
}