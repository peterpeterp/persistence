# I used that to make EKE plots

analyse_whole_EKE <- function(variable='spi3'){
	nc<-open.nc('/home/peter/Dokumente/pik/backuped/data/spi/spi_cru_ts_3_1949-2012_96x73.nc4')
	var<-var.get.nc(nc,variable)
	
	lon<-var.get.nc(nc,"lon")
	lat<-var.get.nc(nc,"lat")

	time<-var.get.nc(nc,"time")

	loc<-array(NA,c(96*73,2))

	season_months<-array(NA,c(5,12))
	season_months[1,]=c(3,4,5)
	season_months[2,]=c(6,7,8)
	season_months[3,]=c(9,10,11)
	season_months[4,]=c(12,1,2)
	season_months[5,]=1:12

	season_indices<-array(NA,c(5,432))

	for (sea in 1:5){
		for (i in 1:36){
			for (j in 1:12){
				if (j %in% season_months[sea,]){season_indices[sea,(i-1)*12+j]=1}
			}
		}
	}

	seasonal_eke<-array(NA,c(5,96*73,432))

	index<-0
	for (x in 1:96){
		cat(x)
		for (y in 1:73){
			index<-index+1
			
			loc[index,]=c(lon[x],lat[y])
			for (sea in 1:5){
				seasonal_eke[sea,index,]=eke[x,y,]
				seasonal_eke[sea,index,which(is.na(season_indices[sea,]))]=NA
			}
		}
	}

	print(index)

	nc_out<-create.nc(paste("../data/eke/eke_ana_96x73_",levelist[lvl],"mbar.nc",sep=""))

	dim.def.nc(nc_out,"seasons",dimlength=5,unlim=FALSE)
	dim.def.nc(nc_out,"ID",dimlength=index,unlim=FALSE)
    dim.def.nc(nc_out,"locs",dimlength=2,unlim=FALSE)
    dim.def.nc(nc_out,"time",dimlength=432,unlim=FALSE)


    var.def.nc(nc_out,"loc","NC_DOUBLE",c(1,3))
    att.put.nc(nc_out, "loc", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "loc", "dim_explanation", "NC_CHAR", "ID-locs")
    att.put.nc(nc_out, "loc", "val_explanation", "NC_CHAR", "lon, lat")

    var.def.nc(nc_out,"seasonal_eke","NC_DOUBLE",c(0,1,4))
    att.put.nc(nc_out, "seasonal_eke", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "seasonal_eke", "dim_explanation", "NC_CHAR", "sea-ID-time")
    att.put.nc(nc_out, "seasonal_eke", "val_explanation", "NC_CHAR", "eke")

    var.put.nc(nc_out,"loc",loc)     
    var.put.nc(nc_out,"seasonal_eke",seasonal_eke)     
    close.nc(nc_out)
}

plot_eke <- function(add_name="full"){
	nc<-open.nc(paste("../data/eke/eke_ana_96x73_500mbar.nc",sep=""))
	eke_ana<<-var.get.nc(nc,"eke_ana")
	loc<<-var.get.nc(nc,"loc")
	ntot<<-96*73

	dat<-list(ID=1:ntot,lon=loc[,1],lat=loc[,2])
	dat$lon[dat$lon>180]=dat$lon[dat$lon>180]-360
	dat<<-dat

	topo_map_plot(filename=paste("../plots/eke/eke_1979-2014_",add_name,"_LR.pdf",sep=""),reihen=eke_ana[,,2],reihen_sig=eke_ana[,,3],farb_mitte=c(-0.01,0.01),farb_palette="lila-gruen",titel=c(""))
	topo_map_plot(filename=paste("../plots/eke/eke_1979-2014_",add_name,"_MN.pdf",sep=""),reihen=eke_ana[,,1],farb_mitte="mean",farb_palette="regenbogen",titel=c(""))

}

master_init <- function(){
    source("map_plot.r")
    source("inits_plot.r")

    library(RNetCDF)
    library(SDMTools)
    library(fields)

    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
}

#analyse_whole_EKE(lvl=3)


master_init()


plot_init_Had()

plot_eke(add_name="500_mb")

plot_init_multi_SH()
plot_eke(add_name="500_mb_SH")