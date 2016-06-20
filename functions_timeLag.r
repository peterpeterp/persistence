# calculates the time lag between max(abs(anomaly)) and mid-point of period or beginning of period

timeLag_analysis <- function(period="1950-2014",ID_select=1:length(dat$ID),add_name="quant_other",folder="/gridded/",ID_name="",ID_names=1:length(dat$ID)){    

	filename <- paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_quantiles.nc",sep="") ; print(filename)
	print(dim(var.get.nc(open.nc(filename),"quantile_stuff")))
    taus<-var.get.nc(open.nc(filename),"taus")
    tau<-which(taus==0.95)
	quantile95<-var.get.nc(open.nc(filename),"quantile_stuff")[1:5,,,tau,1]

	lag_ana<-array(NA,c(5,ntot,2,2))
    x<-(1:(365*66))

    for (sea in 1:5){ 
        season<-season_names[sea]
        cat(paste("\n",season))  

        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,ID_name,"_duration_",season,".nc",sep="") ; print(filename)
        dur<-var.get.nc(open.nc(filename),"dur")
        dur_mid<-var.get.nc(open.nc(filename),"dur_mid")

        percentage<-0
        cat(paste("\n0 -> -> -> -> -> 100\n"))
        for (q in 1:ntot){
            if (q/ntot*100 > percentage){
                cat("-")
                percentage<-percentage+5
            }
            #cat(q)
            for (state in 1:2){
            	extremes<-which(!is.na(dur[q,state,]) & !is.na(quantile95[sea,q,state]) & dur[q,state,]>=quantile95[sea,q,state])
                duration<-dur[q,state,extremes]
                duration_mid<-dur_mid[q,state,extremes]
                lags_toMid<-duration_mid*NA
                lags_toStart<-duration_mid*NA
                if (length(extremes)>2){
	                for (i in 1:length(extremes)){
	                	#event<-round((duration_mid[i]-1950)*365-duration[i]/2) : round((duration_mid[i]-1950)*365+duration[i]/2-1)
                        
                        eventcenter<-which.min(abs( x - (duration_mid[i]-1950)*365))
                        if ((-1)^duration[i]>0){
                            if (eventcenter<(duration_mid[i]-1950)*365){event<-(eventcenter-duration[i]/2+1):(eventcenter+duration[i]/2)}
                            if (eventcenter>(duration_mid[i]-1950)*365){event<-(eventcenter-duration[i]/2):(eventcenter+duration[i]/2-1)}
                        }
                        if ((-1)^duration[i]<0){event<-(eventcenter-(duration[i]-1)/2):(eventcenter+(duration[i]-1)/2)}


	                	if (length(which(is.na(dat$tas[q,,][event])))==0){
			                lags_toMid[i]<-which.max(dat$tas[q,,][event]*(-1)^state)-duration[i]/2
			                lags_toStart[i]<-which.max(dat$tas[q,,][event]*(-1)^state)
			            }

			            if (length(which(is.na(dat$tas[q,,][event])))!=0){
                            print(paste(sea,q,state,duration_mid[i],duration[i]))
                        }
		            }
		            if (length(which(!is.na(lags_toMid)))>3){
			            lag_ana[sea,q,state,1]=mean(lags_toMid,na.rm=TRUE)
			            lag_ana[sea,q,state,2]=mean(lags_toStart,na.rm=TRUE)
			        }
			    }
	        }
	    }
	}

    filename <- paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_timeLags.nc",sep="") ; print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)

    dim.def.nc(nc_out,"seasons",dimlength=5,unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)

    dim.def.nc(nc_out,"outs",dimlength=2,unlim=FALSE)

    var.def.nc(nc_out,"lag_ana","NC_DOUBLE",c(0,1,2,3))
    att.put.nc(nc_out, "lag_ana", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "lag_ana", "dim_explanation", "NC_CHAR", "season-ID-state-..")
    att.put.nc(nc_out, "lag_ana", "explanation", "NC_CHAR", "1-lag between max anomaly and mid of period length, 2-lag between max anomaly and begin of period length")
        
    var.put.nc(nc_out,"lag_ana",lag_ana)      
 
    close.nc(nc_out) 	
}


