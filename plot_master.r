# teste teste
source("write.r")
source("load.r")

if (1==2){
	library(SDMTools)
	source("region_average.r")
	source("map_plot.r")
	source("trend_control.r")
	library(rworldmap)
	library(fields)
	worldmap = getMap(resolution = "low")
	nday = 91
	nyr = 5
	ntot=1319
	dat=dat_load("../data/dat_regional.nc",reg=1)
	per=markov_load(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr))

	numbWarm=trend_control_warm_days(dat,per,c(1),c(365),filename=sprintf("../data/warmTage_tremds_year.nc"))
	reihen=array(numbWarm,dim=c(1,ntot))
	map_allgemein(dat=dat,reihen=reihen,titel=c("yearly increase in 'warm days' from 1950 to 2011"),farbe_mitte="0",
		filename_plot=sprintf("../plots/maps/%s_%s_warm_days.pdf",nday,nyr),worldmap=worldmap,ausschnitt=c(35,66))
	map_regional(dat=dat,toPlot=reihen,titles=c("yearly increase in 'warm days' from 1950 to 2011"),worldmap=worldmap,
		filename_plot=sprintf("../plots/regions/%s_%s_warm_days_regio.pdf",nday,nyr))

}


if (1==1){
	library(SDMTools)
	source("region_average.r")
	source("map_plot.r")
	library(rworldmap)
	library(fields)
	worldmap = getMap(resolution = "low")
	nday = 91
	nyr = 5
	ntot=1319
	dat=dat_load("../data/dat_regional.nc",reg=1)

	if (1==1){
		# markov summer 3states
		vars=c("MK","LR")
        vars_sig=c("MK_sig","LR_sig")
       	seasons=c("spring","summer","autumn","winter","year")
       	for (season in seasons){
			nc=open.ncdf(paste("../data/91_5/91_5_mar3s_trend_",season,".nc",sep=""))
			reihen=array(NA,dim=c(9,ntot))
			reihen_sig=array(NA,dim=c(9,ntot))
			titel=c()

			states=c("cold","normal","warm")

			MK=array(get.var.ncdf(nc,vars[1]),dim=c(ntot,3,3))
			MK_sig=array(get.var.ncdf(nc,vars_sig[1]),dim=c(ntot,3,3))
			for (from in 1:3){
				for (to in 1:3){
					reihen[((from-1)*3+to),]=MK[1:ntot,from,to]
					reihen_sig[((from-1)*3+to),]=MK_sig[1:ntot,from,to]
					titel[((from-1)*3+to)]=paste("Mann-Kendall test for transition from",states[from],"to",states[to],"in",season)
				}
			}

			map_allgemein(dat=dat,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farbe_mitte="gemeinsam 0",
				filename_plot=paste("../plots/maps/91_5_mar3s_trend_",season,".pdf",sep=""),
				worldmap=worldmap,ausschnitt=c(35,66))		
			map_regional(dat=dat,toPlot=reihen,titles=titel,worldmap=worldmap,
				filename_plot=paste("../plots/regions/91_5_mar3s_trend_",season,".pdf",sep=""))
		}
	}


	if (1==2){
		# markov summer
		vars=c("mar_s_w_lr","mar_s_k_lr","mar_s_w_mk","mar_s_k_mk")
        vars_sig=c("mar_s_w_lr_sig","mar_s_k_lr_sig","mar_s_w_mk_sig","mar_s_k_mk_sig")
		nc=open.ncdf(sprintf("../data/%s_%s/%s_%s_markov_trend.nc",nday,nyr,nday,nyr))
		reihen=array(NA,dim=c(4,ntot))
		reihen_sig=array(NA,dim=c(4,ntot))
		titel=c()

		for (i in 1:4){
			reihen[i,]=get.var.ncdf(nc,vars[i])
			reihen_sig[i,]=get.var.ncdf(nc,vars_sig[i])
			for (k in 1:(length(nc$var))){
				if (nc$var[[k]]$name==vars[i]){
					titel[i]=nc$var[[k]]$longname
				}
			}	
		}
		map_allgemein(dat=dat,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farbe_mitte="0",
			filename_plot=sprintf("../plots/maps/%s_%s_markov_trend_summer.pdf",nday,nyr),
			worldmap=worldmap,ausschnitt=c(35,66))		
		map_regional(dat=dat,toPlot=reihen,titles=titel,worldmap=worldmap,
			filename_plot=sprintf("../plots/regions/%s_%s_markov_summer.pdf",nday,nyr))
	}

	if (1==2){
		# summer vergleich
		waka=c("warm","cold")
		titel_zusatz=c("mean","a","a_err","b","b_err","0.02 percentile","0.05 percentile","0.10 percentile")
		vars=c("dur_ana_warm_before","dur_ana_cold_before",
		    	"dur_ana_warm_after","dur_ana_cold_after")

		nc=open.ncdf(sprintf("../data/%s_%s/%s_%s_duration_analysis_summer.nc",nday,nyr,nday,nyr))

		titel=c()
		auswahl=c(1,6,7)
		reihen=array(NA,dim=c(length(auswahl)*2,ntot))

		for (i in 1:length(auswahl)){
		    for (j in 1:2){
		    	reihen[((i-1)*2+j),]=get.var.ncdf(nc,vars[j+2])[1:ntot,auswahl[i]]-get.var.ncdf(nc,vars[j])[1:ntot,auswahl[i]]
		    	titel[((i-1)*2+j)]=paste("difference in",waka[j],"period duration",titel_zusatz[auswahl[i]],"before and after 1980")
		    }
		}
		
		map_regional(dat=dat,toPlot=reihen,titles=titel,worldmap=worldmap,
			filename_plot=sprintf("../plots/regions/%s_%s_duration_summer.pdf",nday,nyr))
		map_allgemein(dat=dat,
			filename_plot=sprintf("../plots/maps/%s_%s_duration_summer_analysis.pdf",nday,nyr),
			worldmap=worldmap,ausschnitt=c(35,66),reihen=reihen,titel=titel,farbe_mitte="0")
	}	
	if (1==2){
	# summer all
	    titel_zusatz=c("mean","a","a_err","b","b_err","0.05 percentile","0.10 percentile")
	    vars=c("dur_ana_warm_full","dur_ana_cold_full")

	    nc=open.ncdf(sprintf("../data/%s_%s/%s_%s_duration_analysis_summer.nc",nday,nyr,nday,nyr))

	    reihen=array(NA,dim=c(14,ntot))
	    titel=c()
	    for (i in 1:7){
	    	for (j in 1:2){
	    		reihen[((i-1)*2+j),]=get.var.ncdf(nc,vars[j])[1:ntot,i]
	    		titel[((i-1)*2+j)]=paste(nc$var[[j]]$longname,titel_zusatz[i])
	    	}
	    }
		map_allgemein(dat=dat,
			filename_plot=sprintf("../plots/maps/%s_%s_duration_summer_diff_1980.pdf",nday,nyr),
			worldmap=worldmap,ausschnitt=c(35,66),reihen=reihen,titel=titel,farbe_mitte="mean")
	}        	        

	
}