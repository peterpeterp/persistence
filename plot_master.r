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
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc",reg=1)

	if (1==1){
		# markov xx states
		state_names=c("cold","warm")
		#state_names=c("cold","normal","warm")
       	seasons=c("spring","summer","autumn","winter","year")
		vars="MK"
        vars_sig="MK_sig"
        farbe_mitte="gemeinsam 0"
        name_zusatz="grob"
		
		states=length(state_names)

       	for (season in seasons){
			nc=open.ncdf(paste("../data/91_5/91_5_mar",states,"s_trend_",season,".nc",sep=""))
			reihen=array(NA,dim=c(states*states,ntot))
			reihen_sig=array(NA,dim=c(states*states,ntot))
			titel=c()


			y=array(get.var.ncdf(nc,vars),dim=c(ntot,states,states))
			if (!is.na(vars_sig)){y_sig=array(get.var.ncdf(nc,vars_sig),dim=c(ntot,states,states))}
			for (from in 1:states){
				for (to in 1:states){
					reihen[((from-1)*states+to),]=y[1:ntot,from,to]
					if (!is.na(vars_sig)){reihen_sig[((from-1)*states+to),]=y_sig[1:ntot,from,to]}
					titel[((from-1)*states+to)]=paste(vars,"for transition from",state_names[from],"to",state_names[to],"in",season)
				}
			}

			map_allgemein(dat=dat,reihen=reihen,reihen_sig=reihen_sig,titel=titel,farbe_mitte=farbe_mitte,
				filename_plot=paste("../plots/maps/91_5_mar",states,"s_",vars,"_",season,"_",name_zusatz,".pdf",sep=""),
				worldmap=worldmap,ausschnitt=c(35,66))		
			map_regional(dat=dat,toPlot=reihen,titles=titel,worldmap=worldmap,
				filename_plot=paste("../plots/regions/91_5_mar",states,"s_",vars,"_",season,"_",name_zusatz,".pdf",sep=""))
		}
	}



	if (1==2){
		# duration vergleich
       	seasons=c("spring","summer","autumn","winter")
	    state_zusatz=c("cold","warm")
		titel_zusatz=c("mean","a","a_err","b","b_err","0.02 percentile","0.05 percentile","0.10 percentile")
		vars=c("dur_ana_before","dur_ana_after")
		auswahl=c(1,6,7)

		for (season in seasons){
			nc=open.ncdf(paste("../data/",nday,"_",nyr,"/",nday,"_",nyr,"_duration_2s_analysis_",season,".nc",sep=""))
			titel=c()
			reihen=array(NA,dim=c(length(auswahl)*2,ntot))
			for (i in 1:length(auswahl)){
			    for (state in 1:2){
			    	reihen[((i-1)*2+state),]=get.var.ncdf(nc,vars[1+1])[1:ntot,state,auswahl[i]]-get.var.ncdf(nc,vars[1])[1:ntot,state,auswahl[i]]
			    	titel[((i-1)*2+state)]=paste("difference in",state_zusatz[state],"period duration",titel_zusatz[auswahl[i]],"before and after 1980 in",season)
			    }
			}
			
			map_regional(dat=dat,toPlot=reihen,titles=titel,worldmap=worldmap,
				filename_plot=paste("../plots/regions/",nday,"_",nyr,"_duration_",season,"_diff_1980.pdf",sep=""))
			map_allgemein(dat=dat,
				filename_plot=paste("../plots/maps/",nday,"_",nyr,"_duration_",season,"_diff_1980.pdf",sep=""),
				worldmap=worldmap,ausschnitt=c(35,66),reihen=reihen,titel=titel,farbe_mitte="0")
		}
	}

	if (1==2){
		# duration climatology
       	seasons=c("spring","summer","autumn","winter")
	    state_zusatz=c("cold","warm")
		titel_zusatz=c("mean","a","a_err","b","b_err","0.02 percentile","0.05 percentile","0.10 percentile")
		vars=c("dur_ana_full")
		auswahl=c(1,6,7)

		for (season in seasons){
			nc=open.ncdf(paste("../data/",nday,"_",nyr,"/",nday,"_",nyr,"_duration_2s_analysis_",season,".nc",sep=""))
			titel=c()
			reihen=array(NA,dim=c(length(auswahl)*2,ntot))
			for (i in 1:length(auswahl)){
			    for (state in 1:2){
			    	reihen[((i-1)*2+state),]=get.var.ncdf(nc,vars[1])[1:ntot,state,auswahl[i]]
			    	titel[((i-1)*2+state)]=paste(state_zusatz[state],"period duration",titel_zusatz[auswahl[i]],"from 1950-2014 in",season)
			    }
			}
			
			map_regional(dat=dat,toPlot=reihen,titles=titel,worldmap=worldmap,
				filename_plot=paste("../plots/regions/",nday,"_",nyr,"_duration_",season,"_climatology.pdf",sep=""))
			map_allgemein(dat=dat,
				filename_plot=paste("../plots/maps/",nday,"_",nyr,"_duration_",season,"_climatology.pdf",sep=""),
				worldmap=worldmap,ausschnitt=c(35,66),reihen=reihen,titel=titel,farbe_mitte="mean")
		}
	}
}