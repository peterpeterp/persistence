#!/home/pepflei/R/bin/Rscript
# Load useful functions 
#dyn.load("persistence_tools.so")
source("functions_support.r")
source("functions_duration.r")
source("write.r")
source("load.r")

master_nas <- function(){
    # count nas
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2011.nc")
    find_nas(dat)    
}

master_trend <- function(nday,nyr,trendID){
    # calculate trend
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend=calc_trend(dat,paste("../data/",trendID,"/",trendID,"_trend.nc",sep=""),nday,nyr)
}

master_trend_control <- function(nday,nyr,trendID){
    # trend control
    source("trend_control.r")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    nc=open.ncdf(paste("../data/",trendID,"/",trendID,"_markov2s.nc",sep=""))
    ind=get.var.ncdf(nc,"ind")
    trend_control_warm_days(dat,ind)   
}

master_markov <- function(nday,nyr,trendID,states,transition_names){
    # calculate persistence 2 states
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend=trend_load(paste("../data/",trendID,"/",trendID,"_trend.nc",sep=""))
    if (states==2){
        per=calc_per(dat,trend,nday,nyr,model=markov_2states,states=states,transition_names=transition_names,
            filename=paste("../data/",trendID,"/",states,"_states/markov/",trendID,"_markov_",states,"states.nc",sep=""))
    }
    if (states==3){
        per=calc_per(dat,trend,nday,nyr,model=markov_3states,states=states,transition_names=transition_names,
            filename=paste("../data/",trendID,"/",states,"_states/markov/",trendID,"_markov_",states,"states.nc",sep=""))
    }
}


master_analyse_markov <- function(yearPeriod,trendID,states,transition_names){
    # analyse persistence 2 states
    nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/markov/",trendID,"_markov_",states,"states.nc",sep=""))
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    seasons=c("spring","summer","autumn","winter","year")
    for (sea in 1:length(seasons)){
        print(seasons[sea])
        per=get.var.ncdf(nc,"markov")[1:length(dat$ID),sea,,]
        tmp=global_analysis(toAna=per,yearPeriod=yearPeriod)
        markov_analysis_write(filename=paste("../data/",trendID,"/",states,"_states/markov/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_mar",states,"s_trend_",seasons[sea],".nc",sep=""),
            analysis=tmp,season=seasons[sea],transition_names=transition_names)
    }
}

master_duration <- function(nday,nyr,trendID,states){
    # calculate duration periods 2 states
    trash=((nyr-1)/2*365+(nday-1))

    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")

    if (states==2){
        stateIndeces=c(-1,1)
    }
    if (states==3){
        stateIndeces=c(-1,0,1)
    }

    nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/markov/",trendID,"_markov_",states,"states.nc",sep=""))        
    ind=get.var.ncdf(nc,"ind")

    calc_global_dur(dat=dat,ind=ind,trash=trash,filename=paste("../data/",trendID,"/",states,"_states/duration/",trendID,"_duration_",states,"s_year.nc",sep=""),states=stateIndeces)

    nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/duration/",trendID,"_duration_",states,"s_year.nc",sep=""))
    dur=get.var.ncdf(nc,"dur")
    dur_mid=get.var.ncdf(nc,"dur_mid")
    
    duration_seasons(dur,dur_mid,season=c(59,151),filename=paste("../data/",trendID,"/",states,"_states/duration/",trendID,"_duration_",states,"s_spring.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(151,243),filename=paste("../data/",trendID,"/",states,"_states/duration/",trendID,"_duration_",states,"s_summer.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(243,335),filename=paste("../data/",trendID,"/",states,"_states/duration/",trendID,"_duration_",states,"s_autumn.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(335,425),filename=paste("../data/",trendID,"/",states,"_states/duration/",trendID,"_duration_",states,"s_winter.nc",sep=""))
}

master_analyse_duration <- function(yearPeriod,trendID,states,seasons=c("spring","summer","autumn","winter","year")){
    # analyse duration periods 2 states
    for (season in seasons){
        print(season)
        nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/duration/",trendID,"_duration_",states,"s_",season,".nc",sep=""))
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")
        duration_analysis(dur,dur_mid,filename=paste("../data/",trendID,"/",states,"_states/duration/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_duration_",states,"s_analysis_",season,".nc",sep=""),
            season=season,yearPeriod)#,stations=seq(487,489,1))#,stations=134:136)
    }
}

master_duration_distribution <- function(yearPeriod,trendID,states,seasons=c("spring","summer","autumn","winter","year")){
    # analyse duration periods 2 states
    for (season in seasons){
        print(season)
        nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/duration/",trendID,"_duration_",states,"s_",season,".nc",sep=""))
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")
        duration_distribution(dur,dur_mid,filename=paste("../data/",trendID,"/",states,"_states/duration/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_duration_",states,"s_distribution_",season,".nc",sep=""),
            season=season,yearPeriod)#,stations=seq(487,489,1))#,stations=134:136)
    }
}

master_regional_trend <- function(yearPeriod,region_name,trendID){
    source("functions_regional.r")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    regional_analysis(dat=dat,yearPeriod,filepath=paste("../data/",trendID,"/2_states/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",sep=""),region_name=region_name)
}

master_regional_climatology <- function(yearPeriod,region_name,trendID){
    source("functions_regional.r")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    regional_climatology(trendID=trendID,dat=dat,yearPeriod=yearPeriod,region_name=region_name)
    plot_regional_distributions(trendID,dat,yearPeriod,region_name)
    plot_regional_boxplots(trendID,dat,yearPeriod,region_name)
}




full_2states <- function(nday,nyr){
    #complete 2 states analysis 
    trendID=paste(nday,"_",nyr,sep="")

    #master_trend(nday,nyr,trendID)

    #master_markov(nday,nyr,trendID,states=2,transition_names=c("cc wc cw ww"))

    points=c(1950,2014,1950,1980,1980,2014)
    for (i in 1:3){
        period=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        master_analyse_markov(yearPeriod=period,trendID,states=2,transition_names=c("cc wc cw ww"))
    }

    #master_duration(nday,nyr,trendID,2)
    for (i in 1:3){
        period=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        #master_analyse_duration(yearPeriod=period,trendID,states=2)
        #master_duration_distribution(yearPeriod=period,trendID,states=2)
    }
}

full_3states <- function(nday,nyr){
    #complete 3 states analysis 
    trendID=paste(nday,"_",nyr,sep="")

    #master_trend(nday,nyr,trendID)

    #master_markov(nday,nyr,trendID,states=3,transition_names=c("cc nc wc cn nn nw cw nw ww"))

    points=c(1950,2014,1950,1980,1980,2014)
    for (i in 1:3){
        period=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        #master_analyse_markov(yearPeriod=period,trendID,states=3,transition_names=c("cc nc wc cn nn nw cw nw ww"))
    }

    master_duration(nday,nyr,trendID,states=3)
    for (i in 1:3){
        period=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        master_analyse_duration(yearPeriod=period,trendID,states=3)
    }
}

#full_3states(91,5)
#full_3states(91,3)
#full_2states(91,5)

#master_regional_climatology(yearPeriod=c(1950,2014),region_name="7rect",trendID="91_5")
master_regional_climatology(yearPeriod=c(1950,2014),region_name="midlat",trendID="91_5")
#master_regional_climatology(yearPeriod=c(1950,1980),region_name="7rect",trendID="91_5")
#master_regional_climatology(yearPeriod=c(1980,2014),region_name="7rect",trendID="91_5")