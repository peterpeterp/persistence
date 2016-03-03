
state_ana <- function(seasonStart=c(60,152,244,335,1),seasonStop=c(151,243,334,424,365)){
    ntot<-length(dat$ID)
    stateIndeces<-c(-1,1)
    print(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,dataset,"_state_ind.nc",sep=""))
    nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,dataset,"_state_ind.nc",sep=""))
    ind<-var.get.nc(nc,"ind")

    ind_ana<-array(NA,dim=c(5,3,ntot,7))
    percentage<-0
    cat(paste("\n0 -> -> -> -> -> 100\n"))
    for (q in 1:ntot){
        if (q/ntot*100 > percentage){
            cat("-")
            percentage<-percentage+5
        }        
        for (sea in 1:length(seasonStart)){
            stateCount<-array(NA,c(3,length(dat$year)))
            for (state in 1:2){
	            for (yr in (yearPeriod[1]-dat$year[1]+1):(yearPeriod[2]-dat$year[1]+1)){
	                if (seasonStop[sea]>365 & yr<length(dat$year)){
	                    stateCount[state,yr]=length(which(ind[q,seasonStart[sea]:365,yr]==stateIndeces[state]))                 	

	                    stateCount[state,yr]=stateCount[state,yr]+length(which(ind[q,1:(seasonStop[sea]-365),(yr+1)]==stateIndeces[state])) 
	                }
	                if (seasonStop[sea]<=365){
	                    stateCount[state,yr]=length(which(ind[q,seasonStart[sea]:seasonStop[sea],yr]==stateIndeces[state])) 
	                }               
	            }
	        }
	        if (sum(stateCount,na.rm=TRUE)>((seasonStop[sea]-seasonStart[sea])*10)){
	        	for (state in 1:2){
		            lm.r<-lm(stateCount[state,]~dat$year)
		            ind_ana[sea,state,q,1]=summary(lm.r)$coefficients[2]
		            ind_ana[sea,state,q,2]=summary(lm.r)$coefficients[4]

		            percentage_array<-stateCount[state,]/(stateCount[1,]+stateCount[2,])
		           	lm.r<-lm(percentage_array~dat$year)
		    		ind_ana[sea,state,q,3]=summary(lm.r)$coefficients[2]
		            ind_ana[sea,state,q,4]=summary(lm.r)$coefficients[4]

		            ind_ana[sea,state,q,5]=sum(stateCount[state,],na.rm=TRUE)
		            ind_ana[sea,state,q,6]=sum(stateCount[state,],na.rm=TRUE)/sum(c(stateCount[1:2,]),na.rm=TRUE)
		            ind_ana[sea,state,q,7]=sd(percentage_array,na.rm=TRUE)
		       	}
	        }
	        #same for nas
	       	for (yr in 1:length(dat$year)){
	            if (seasonStop[sea]>365 & yr<length(dat$year)){
	                stateCount[3,yr]=length(which(is.na(ind[q,seasonStart[sea]:365,yr])))                 	

	                stateCount[3,yr]=stateCount[3,yr]+length(which(is.na(ind[q,1:(seasonStop[sea]-365),(yr+1)]))) 
	            }
	            if (seasonStop[sea]<=365){
	                stateCount[3,yr]=length(which(is.na(ind[q,seasonStart[sea]:seasonStop[sea],yr])))
	            }               
	        }
	        if (length(which(is.na(stateCount[3,])))>10){
	        	if (TRUE){
		            lm.r<-lm(stateCount[3,]~dat$year)
		            ind_ana[sea,3,q,1]=summary(lm.r)$coefficients[2]
		            ind_ana[sea,3,q,2]=summary(lm.r)$coefficients[4]

		            percentage_array<-stateCount[3,]/(stateCount[1,]+stateCount[2,])
		           	lm.r<-lm(percentage_array~dat$year)
		    		ind_ana[sea,3,q,3]=summary(lm.r)$coefficients[2]
		            ind_ana[sea,3,q,4]=summary(lm.r)$coefficients[4]

		            ind_ana[sea,3,q,5]=sum(stateCount[3,],na.rm=TRUE)
		            ind_ana[sea,3,q,6]=sum(stateCount[3,],na.rm=TRUE)/sum(c(stateCount[1:2,]),na.rm=TRUE)
		            ind_ana[sea,3,q,7]=sd(percentage_array,na.rm=TRUE)
		        }
		    }
        }
    }
    print(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",period,"_state_percentage.nc",sep=""))
    nc_out<-create.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",period,"_state_percentage.nc",sep=""))
    
    dim.def.nc(nc_out,"sea",dimlength=5, unlim=FALSE)       
    dim.def.nc(nc_out,"state",dimlength=3, unlim=FALSE)   
    dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
    dim.def.nc(nc_out,"out",dimlength=7, unlim=FALSE)
    dim.def.nc(nc_out,"years",dimlength=length(dat$year), unlim=FALSE)

    var.def.nc(nc_out,"ind_ana","NC_DOUBLE",c(0,1,2,3))
    att.put.nc(nc_out, "ind_ana", "missing_value", "NC_DOUBLE", -99999)
    att.put.nc(nc_out, "ind_ana", "dim_explanation", "NC_CHAR", "sea-state-ID-...")
    att.put.nc(nc_out, "ind_ana", "explanation", "NC_CHAR", "increase in dry (wet) days, singificance of increase, increase of percentage, signi of increase, total count of days, total percentage of dry (wet) days, sd of percentage")

    var.put.nc(nc_out,"ind_ana",ind_ana)      
}

state_ana_view <- function(yAusschnitt=c(20,50),xAusschnitt=c(220,310),asp=1,paper=c(8,3.5),pointsize=0.58){
    print(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",period,"_state_percentage.nc",sep=""))
    nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",period,"_state_percentage.nc",sep=""))	
    ind_ana<<-var.get.nc(nc,"ind_ana")

    print(ind_ana[2,2,1:100,2])
    print(ind_ana[2,2,1:100,1])
    print(ind_ana[2,2,1:100,3])
    print(ind_ana[2,2,1:100,1]/ind_ana[2,2,1:100,3])
    asdas

    print(dim(ind_ana))
    seas<-5
    states<-3
    # increase
    reihen=array(NA,c(seas*states,ntot))
    reihen_sig=array(NA,c(seas*states,ntot))
    titel=c()
    index=0
    for (sea in 1:seas){
    	for (state in 1:states){
    		index<-index+1
    		reihen[index,]=ind_ana[sea,state,,1]
    		reihen_sig[index,]=ind_ana[sea,state,,2]
    		titel[index]=paste(season_names[sea],"increase of",state_names[state],"days")
    	}
    }
	topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",period,"_state_incrCount.pdf",sep=""),reihen=reihen,reihen_sig=reihen_sig,titel=titel,signi_level=0.05,farb_mitte="0",farb_palette="lila-gruen",yAusschnitt=yAusschnitt,xAusschnitt=xAusschnitt,asp=asp,paper=paper,pointsize=pointsize)

    # increase
    reihen=array(NA,c(seas*states,ntot))
    reihen_sig=array(NA,c(seas*states,ntot))
    titel=c()
    index=0
    for (sea in 1:seas){
    	for (state in 1:states){
    		index<-index+1
    		reihen[index,]=ind_ana[sea,state,,3]
    		reihen_sig[index,]=ind_ana[sea,state,,4]
    		titel[index]=paste(season_names[sea],"increase of",state_names[state],"percentage")
    	}
    }
	topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",period,"_state_incrPerc.pdf",sep=""),reihen=reihen,reihen_sig=reihen_sig,titel=titel,signi_level=0.05,farb_mitte="0",farb_palette="lila-gruen",yAusschnitt=yAusschnitt,xAusschnitt=xAusschnitt,asp=asp,paper=paper,pointsize=pointsize)

	# percentage
    reihen=array(NA,c(seas*states,ntot))
    reihen_sig=array(NA,c(seas*states,ntot))
    titel=c()
    index=0
    for (sea in 1:seas){
    	for (state in 1:states){
    		index<-index+1
    		reihen[index,]=ind_ana[sea,state,,5]
    		titel[index]=paste(season_names[sea],"count of",state_names[state],"days")
    	}
    }
	topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",period,"_state_count.pdf",sep=""),reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=c(0,1),farb_palette="lila-gruen",yAusschnitt=yAusschnitt,xAusschnitt=xAusschnitt,asp=asp,paper=paper,pointsize=pointsize)

	# percentage
    reihen=array(NA,c(seas*states,ntot))
    reihen_sig=array(NA,c(seas*states,ntot))
    titel=c()
    index=0
    for (sea in 1:seas){
    	for (state in 1:states){
    		index<-index+1
    		reihen[index,]=ind_ana[sea,state,,6]
    		titel[index]=paste(season_names[sea],"percentage of",state_names[state],"days")
    	}
    }
	topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",period,"_state_perc.pdf",sep=""),reihen=reihen,reihen_sig=reihen_sig,titel=titel,farb_mitte=c(0,1),farb_palette="lila-gruen",yAusschnitt=yAusschnitt,xAusschnitt=xAusschnitt,asp=asp,paper=paper,pointsize=pointsize)
}

mean_dur_no_memory <- function(p){
	meanNoMem<-0
	for (i in 1:100){meanNoMem<-meanNoMem+i*p^i}
	return(meanNoMem)
}

mean_dur_no_memory_change <- function(p,dp){
	meanNoMem<-0
	for (i in 1:100){meanNoMem<-meanNoMem+(i^2)*p^(i-1)*dp}
	return(meanNoMem)
}

memory_test <- function(){
    print(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",period,"_state_percentage.nc",sep=""))
    nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",period,"_state_percentage.nc",sep=""))	
    ind_ana<-var.get.nc(nc,"ind_ana")

	print(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_others",".nc",sep=""))
	nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_others",".nc",sep=""))
	values<-var.get.nc(nc,"other_stuff")	

	for (sea in 1:5){
		for (state in 1:2){
			plot(NA,xlim=c(-0.1,0.1),ylim=c(-0.1,0.1))
		}
	}
}