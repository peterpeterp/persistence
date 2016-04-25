
state_ana <- function(seasonStart=c(60,152,244,335,1),seasonStop=c(151,243,334,424,365)){
    ntot<-length(dat$ID)
    stateIndeces<-c(-1,1)
    print(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,dataset,"_state_ind.nc",sep=""))
    nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,dataset,"_state_ind.nc",sep=""))
    ind<<-var.get.nc(nc,"ind")

    ind_ana<-array(NA,dim=c(5,2,ntot,7))
    percentage<-0
    cat(paste("\n0 -> -> -> -> -> 100\n"))
    for (q in 1:ntot){
        if (q/ntot*100 > percentage){
            cat("-")
            percentage<-percentage+5
        }        
        for (sea in 1:length(seasonStart)){
            stateCount<-array(NA,c(2,length(dat$year)))
	        for (yr in (yearPeriod[1]-dat$year[1]+1):(yearPeriod[2]-dat$year[1]+1)){
	            for (state in 1:2){
	                if (seasonStop[sea]>365 & yr<length(dat$year)){
	                    stateCount[state,yr]=length(which(ind[q,seasonStart[sea]:365,yr]==stateIndeces[state]))                 	

	                    stateCount[state,yr]=stateCount[state,yr]+length(which(ind[q,1:(seasonStop[sea]-365),(yr+1)]==stateIndeces[state])) 
	                }
	                if (seasonStop[sea]<=365){
	                    stateCount[state,yr]=length(which(ind[q,seasonStart[sea]:seasonStop[sea],yr]==stateIndeces[state])) 
	                }   
	            }
	            if (sum(stateCount[,yr],na.rm=TRUE)<50){stateCount[,yr]=NA}
	        }
	        if (sum(stateCount,na.rm=TRUE)>((seasonStop[sea]-seasonStart[sea])*10)){
	        	for (state in 1:2){
		            lm.r<-lm(stateCount[state,]~dat$year)
		            ind_ana[sea,state,q,1]=summary(lm.r)$coef[2]
		            ind_ana[sea,state,q,2]=summary(lm.r)$coef[8]

		            percentage_array<-stateCount[state,]/(stateCount[1,]+stateCount[2,])
		           	lm.r<-lm(percentage_array~dat$year)
		    		ind_ana[sea,state,q,3]=summary(lm.r)$coef[2]
		            ind_ana[sea,state,q,4]=summary(lm.r)$coef[8]

		            ind_ana[sea,state,q,5]=sum(stateCount[state,],na.rm=TRUE)
		            ind_ana[sea,state,q,6]=sum(stateCount[state,],na.rm=TRUE)/sum(c(stateCount[1:2,]),na.rm=TRUE)
		            ind_ana[sea,state,q,7]=sd(percentage_array,na.rm=TRUE)
		       	}
	        }
        }
    }
    iii<<-ind_ana
    print(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",period,"_state_percentage.nc",sep=""))
    nc_out<-create.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",period,"_state_percentage.nc",sep=""))
    
    dim.def.nc(nc_out,"sea",dimlength=5, unlim=FALSE)       
    dim.def.nc(nc_out,"state",dimlength=2, unlim=FALSE)   
    dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
    dim.def.nc(nc_out,"out",dimlength=7, unlim=FALSE)

    var.def.nc(nc_out,"ind_ana","NC_DOUBLE",c(0,1,2,3))
    att.put.nc(nc_out, "ind_ana", "missing_value", "NC_DOUBLE", -99999)
    att.put.nc(nc_out, "ind_ana", "dim_explanation", "NC_CHAR", "sea-state-ID-...")
    att.put.nc(nc_out, "ind_ana", "explanation", "NC_CHAR", "increase in dry (wet) days, singificance of increase, increase of percentage, signi of increase, total count of days, total percentage of dry (wet) days, sd of percentage")

    var.put.nc(nc_out,"ind_ana",ind_ana)      
}

state_ana_view <- function(){
    print(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",period,"_state_percentage.nc",sep=""))
    nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",period,"_state_percentage.nc",sep=""))	
    ind_ana<-var.get.nc(nc,"ind_ana")

    print(dim(ind_ana))
    seas<-5
    states<-2
    # increase
    reihen1=array(NA,c(seas*states,ntot))
    reihen2=array(NA,c(seas*states,ntot))
    reihen3=array(NA,c(seas*states,ntot))
    reihen4=array(NA,c(seas*states,ntot))
    reihen_sig1=array(NA,c(seas*states,ntot))
    reihen_sig2=array(NA,c(seas*states,ntot))
    reihen_sig3=array(NA,c(seas*states,ntot))
    reihen_sig4=array(NA,c(seas*states,ntot))
    titel1=c()
    titel2=c()
    titel3=c()
    titel4=c()
    index=0

    for (sea in 1:seas){
    	for (state in 1:states){
    		index<-index+1
            reihen1[index,]=ind_ana[sea,state,,1]
            reihen2[index,]=ind_ana[sea,state,,3]
            reihen3[index,]=ind_ana[sea,state,,5]
            reihen4[index,]=ind_ana[sea,state,,6]
            reihen_sig1[index,]=ind_ana[sea,state,,2]
            reihen_sig2[index,]=ind_ana[sea,state,,4]
    		titel1[index]=paste(season_names[sea],"increase of",state_names[state],"days")
            titel2[index]=paste(season_names[sea],"increase of",state_names[state],"percentage")
            titel3[index]=paste(season_names[sea],"count of",state_names[state],"days")
            titel4[index]=paste(season_names[sea],"percentage of",state_names[state],"days")

    	}
    }
	topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",period,"_state_incrCount.pdf",sep=""),reihen=reihen1,reihen_sig=reihen_sig1,titel=titel1,signi_level=0.05,farb_mitte="0",farb_palette="lila-gruen")
    topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",period,"_state_incrPerc.pdf",sep=""),reihen=reihen2,reihen_sig=reihen_sig2,titel=titel2,signi_level=0.05,farb_mitte="0",farb_palette="lila-gruen")
    topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",period,"_state_count.pdf",sep=""),reihen=reihen3,reihen_sig=reihen_sig3,titel=titel3,farb_mitte=c(0,1),farb_palette="lila-gruen")
    topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",period,"_state_perc.pdf",sep=""),reihen=reihen4,reihen_sig=reihen_sig4,titel=titel4,farb_mitte=c(0,1),farb_palette="lila-gruen")
}

mean_dur_no_memory <- function(p){
	meanNoMem<-0
	for (i in 1:100){meanNoMem<-meanNoMem+i*p^i}
	return(meanNoMem)
}

mean_dur_no_memory_change <- function(p,dp){
	meanNoMem<-p*0
	for (i in 1:100){meanNoMem<-meanNoMem+(i^2)*p^(i-1)*dp}
	return(meanNoMem)
}

memory_test <- function(yAusschnitt=c(20,80),xAusschnitt=c(-30,80),asp=1,paper=c(8,5),pointsize=0.44){
    print(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",period,"_state_percentage.nc",sep=""))
    nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",period,"_state_percentage.nc",sep=""))	
    ind_ana<-var.get.nc(nc,"ind_ana")

	print(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_others",".nc",sep=""))
	nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_others",".nc",sep=""))
	values<-var.get.nc(nc,"other_stuff")	

    reihen1=array(NA,c(10,ntot))
    reihen2=array(NA,c(10,ntot))
    reihen_sig1=array(NA,c(10,ntot))
    reihen_sig2=array(NA,c(10,ntot))
    titel1=c()
    titel2=c()
    index=0
    pdf(paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",period,"_memory.pdf",sep=""))
	for (sea in 1:5){
		for (state in 1:2){
            index<-index+1
            reihen1[index,]=values[sea,,state,4]-mean_dur_no_memory_change(ind_ana[sea,state,,6],ind_ana[sea,state,,3])
            reihen2[index,]=mean_dur_no_memory_change(ind_ana[sea,state,,6],ind_ana[sea,state,,3])
            titel1[index]=paste(season_names[sea],"difference between measured and no memory change in mean dur for",state_names[state],"periods")
            titel2[index]=paste(season_names[sea],"expected no memory change in mean dur for",state_names[state],"periods")

            plot(NA,xlim=c(-0.1,0.1),ylim=c(-0.1,0.1),main=paste(season_names[sea],state_names[state]),xlab="measured increase in mean(dur)",ylab="expected increase in mean(dur) without memory")
            points(values[sea,,state,4],mean_dur_no_memory_change(ind_ana[sea,state,,6],ind_ana[sea,state,,3]),cex=0.5,pch=20)
		}
	}
    graphics.off()
    topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",period,"_diffToNoMemModel.pdf",sep=""),reihen=reihen1,reihen_sig=reihen_sig1,titel=titel1,farb_mitte=c(-0.1,0.1),farb_palette="lila-gruen",yAusschnitt=yAusschnitt,xAusschnitt=xAusschnitt,asp=asp,paper=paper,pointsize=pointsize)
    topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",period,"/",trendID,"_",period,"_noMemModel.pdf",sep=""),reihen=reihen2,reihen_sig=reihen_sig2,titel=titel2,farb_mitte=c(-0.1,0.1),farb_palette="lila-gruen",yAusschnitt=yAusschnitt,xAusschnitt=xAusschnitt,asp=asp,paper=paper,pointsize=pointsize)
}