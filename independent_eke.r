



analyse_whole_EKE <- function(lvl=1){
	nc<-open.nc("../data/sonstiges/eke_roh/EKE_ERA_Interim_1979-2014_calendar_96x73.nc")
	u2<-var.get.nc(nc,"u2syn")
	v2<-var.get.nc(nc,"v2syn")
	
	lon<-var.get.nc(nc,"lon")
	lat<-var.get.nc(nc,"lat")
	levelist<-var.get.nc(nc,"levelist")

	time<-var.get.nc(nc,"time")+1


	eke<-(u2[,,lvl,]+v2[,,lvl,])/2

	eke_ana<-array(NA,c(5,96*73,5))
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

				eke_ana[sea,index,1]=mean(seasonal_eke[sea,index,],na.rm=TRUE)
				eke_ana[sea,index,2:3]=summary(lm(seasonal_eke[sea,index,which(!is.na(season_indices[sea,]))]~time[which(!is.na(season_indices[sea,]))],na.rm=TRUE))$coef[c(2,8)]
			}
		}
	}

	print(index)

	nc_out<-create.nc(paste("../data/eke/eke_ana_96x73_",levelist[lvl],"mbar.nc",sep=""))

	dim.def.nc(nc_out,"seasons",dimlength=5,unlim=FALSE)
	dim.def.nc(nc_out,"ID",dimlength=index,unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=5,unlim=FALSE)
    dim.def.nc(nc_out,"locs",dimlength=2,unlim=FALSE)
    dim.def.nc(nc_out,"time",dimlength=432,unlim=FALSE)


    var.def.nc(nc_out,"loc","NC_DOUBLE",c(1,3))
    att.put.nc(nc_out, "loc", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "loc", "dim_explanation", "NC_CHAR", "ID-locs")
    att.put.nc(nc_out, "loc", "val_explanation", "NC_CHAR", "lon, lat")

    var.def.nc(nc_out,"eke_ana","NC_DOUBLE",c(0,1,2))
    att.put.nc(nc_out, "eke_ana", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "eke_ana", "dim_explanation", "NC_CHAR", "sea-ID-outs")
    att.put.nc(nc_out, "eke_ana", "val_explanation", "NC_CHAR", "mean, slope ,sig")

    var.def.nc(nc_out,"seasonal_eke","NC_DOUBLE",c(0,1,4))
    att.put.nc(nc_out, "seasonal_eke", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "seasonal_eke", "dim_explanation", "NC_CHAR", "sea-ID-time")
    att.put.nc(nc_out, "seasonal_eke", "val_explanation", "NC_CHAR", "eke")

    var.put.nc(nc_out,"loc",loc)     
    var.put.nc(nc_out,"eke_ana",eke_ana)     
    var.put.nc(nc_out,"seasonal_eke",seasonal_eke)     
    close.nc(nc_out)
}

plot_eke <- function(){
	nc<-open.nc(paste("../data/eke/eke_ana_96x73_850mbar.nc",sep=""))
	eke_ana<<-var.get.nc(nc,"eke_ana")
	loc<<-var.get.nc(nc,"loc")
	ntot<<-96*73

	dat<-list(ID=1:ntot,lon=loc[,1],lat=loc[,2])
	dat$lon[dat$lon>180]=dat$lon[dat$lon>180]-360
	dat<<-dat

	topo_map_plot(filename=paste("../plots/eke/eke_1979-2014_full_LR.pdf",sep=""),reihen=eke_ana[,,2],reihen_sig=eke_ana[,,3],farb_mitte=c(-0.01,0.01),farb_palette="lila-gruen",titel=c(""))

}

master_init <- function(){
    source("map_plot.r")
    source("inits_plot.r")

    library(RNetCDF)
    library(SDMTools)
    library(fields)

    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
}

#analyse_whole_EKE()


master_init()
plot_init_multi_SH()
plot_eke()