

timeLag_analysis <- function(period="1950-2014",ID_select=1:length(dat$ID),add_name="quant_other",folder="/gridded/",ID_name="",ID_names=1:length(dat$ID)){    

	filename <- paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_quantiles.nc",sep="") ; print(filename)
	quantile95<-var.get.nc(open.nc(filename),"quantile_stuff")[1:5,,,2,1]

	lag_ana<-array(NA,c(5,ntot,2,2))

	#lags<-array(NA,c(5,ntot,2,length(dat$year)*100,5))
	#print(33)
	#lens<-array(NA,c(5))

    for (sea in 4:5){ 
        season<-season_names[sea]
        cat(paste("\n",season))  

        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,ID_name,"_duration_",season,".nc",sep="") ; print(filename)
        dur<-var.get.nc(open.nc(filename),"dur")
        dur_mid<-var.get.nc(open.nc(filename),"dur_mid")
        #lens[sea]<-dim(dur_mid)[3]

        percentage<-0
        cat(paste("\n0 -> -> -> -> -> 100\n"))
        for (q in 1:3){
        	cat("q")
            if (q/ntot*100 > percentage){
                cat("-")
                percentage<-percentage+5
            }
            for (state in 1:2){
            	extremes<-which(!is.na(dur[q,state,]) & dur[q,state,]>=quantile95[sea,q,state])
                duration<-dur[q,state,extremes]
                duration_mid<-dur_mid[q,state,extremes]
                lags_toMid<-duration_mid*NA
                lags_toStart<-duration_mid*NA
                for (i in 1:length(extremes)){
                	cat(i)
                	cat("-")
	                event<-round((duration_mid[i]-1950)*365-duration[i]/2+1) : round((duration_mid[i]-1950)*365+duration[i]/2)
	                lags_toMid[i]<-which.max(dat$tas[q,,][event]*(-1)^state)-duration[i]/2
	                lags_toStart[i]<-which.max(dat$tas[q,,][event]*(-1)^state)
	                
	                #lags[sea,q,state,i,1]<-which.max(dat$tas[q,,][event]*(-1)^state)-duration[i]/2
	                #lags[sea,q,state,i,2]<-which.max(dat$tas[q,,][event]*(-1)^state)
	                #lags[sea,q,state,i,4]<-duration[i]
	                #lags[sea,q,state,i,5]<-duration_mid[i]
	            }
	            lag_ana[sea,q,state,1]=mean(lags_toMid,na.rm=TRUE)
	            lag_ana[sea,q,state,2]=mean(lags_toStart,na.rm=TRUE)
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

    #dim.def.nc(nc_out,"events",dimlength=max(lens,na.rm=TRUE),unlim=FALSE)
    #dim.def.nc(nc_out,"outs",dimlength=5,unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=2,unlim=FALSE)

    var.def.nc(nc_out,"lag_ana","NC_DOUBLE",c(0,1,2,3))
    att.put.nc(nc_out, "lag_ana", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "lag_ana", "dim_explanation", "NC_CHAR", "season-ID-state-..")
    att.put.nc(nc_out, "lag_ana", "explanation", "NC_CHAR", "1-lag between max anomaly and mid of period length, 2-lag between max anomaly and begin of period length")
        
    var.put.nc(nc_out,"lag_ana",lag_ana)      
 
    close.nc(nc_out) 	
}


