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

master_trend <- function(nday,nyr,trendID,trend_style="_mean",additional_style="_TX"){
    # calculate trend
    # choice between mean, median and estimated mode is possible
    if (trend_style=="_mean"){procedure=r_calc_runmean_2D}
    if (trend_style=="_median"){procedure=r_calc_runmedian_2D}
    if (trend_style=="_mode"){procedure=r_calc_runmode_2D}
    if (additional_style=="_TX"){dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")}
    if (additional_style=="_TN"){dat=dat_load("../data/HadGHCND_TN_data3D.day1-365.1950-2014.nc")}    
    trend=calc_trend(dat,paste("../data/",trendID,"/",trendID,"_trend",trend_style,additional_style,".nc",sep=""),nday,nyr,procedure=procedure)
}

master_trend_control <- function(trendID,trend_style="_mean",dataset="_TX",additional_style="_median"){
    # trend control
    source("trend_view.r")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    nc=open.ncdf(paste("../data/",trendID,"/",trendID,trend_style,dataset,"_state_ind",additional_style,".nc",sep=""))
    ind=get.var.ncdf(nc,"ind")
    trend_control_warm_days(dat,ind,filename=paste("../data/",trendID,"/",additional_style,"/sonstiges/",trendID,trend_style,dataset,"_5seasons_warme_tage_control",additional_style,".txt",sep=""))   
}

master_state_attribution <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style="_median"){
    # calculate persistence 2 states
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend=trend_load(paste("../data/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    per=state_attribution(dat,trend,nday,nyr,
            filename=paste("../data/",trendID,"/",trendID,trend_style,dataset,"_state_ind",additional_style,".nc",sep=""))

}

master_markov <- function(nday,nyr,trendID,states,transition_names,trend_style="_mean",additional_style=""){
    # calculate persistence 2 states
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend=trend_load(paste("../data/",trendID,"/",trendID,"_trend",trend_style,additional_style,".nc",sep=""))
    if (states==2){
        per=calc_per(dat,trend,nday,nyr,model=markov_2states,states=states,transition_names=transition_names,
            filename=paste("../data/",trendID,"/",states,"_states",trend_style,additional_style,"/markov/",trendID,"_markov_",states,"states.nc",sep=""))
    }
    if (states==3){
        per=calc_per(dat,trend,nday,nyr,model=markov_3states,states=states,transition_names=transition_names,
            filename=paste("../data/",trendID,"/",states,"_states",trend_style,additional_style,"/markov/",trendID,"_markov_",states,"states.nc",sep=""))
    }
}


master_analyse_markov <- function(yearPeriod,trendID,states,transition_names,trend_style="_mean",additional_style=""){
    # analyse persistence 2 states
    nc=open.ncdf(paste("../data/",trendID,"/",states,"_states",trend_style,additional_style,"/markov/",trendID,"_markov_",states,"states.nc",sep=""))
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    seasons=c("spring","summer","autumn","winter","year")
    for (sea in 1:length(seasons)){
        print(seasons[sea])
        per=get.var.ncdf(nc,"markov")[1:length(dat$ID),sea,,]
        tmp=global_analysis(toAna=per,yearPeriod=yearPeriod)
        markov_analysis_write(filename=paste("../data/",trendID,"/",states,"_states",trend_style,additional_style,"/markov/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_mar",states,"s_trend_",seasons[sea],".nc",sep=""),
            analysis=tmp,season=seasons[sea],transition_names=transition_names)
    }
}

master_duration <- function(nday,nyr,trendID,states,trend_style="_mean",dataset="_TX",additional_style="_median"){
    # calculate duration periods 2 states
    trash=((nyr-1)/2*365+(nday-1))

    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")

    if (states==2){
        stateIndeces=c(-1,1)
    }

    nc=open.ncdf(paste("../data/",trendID,"/",trendID,trend_style,dataset,"_state_ind",additional_style,".nc",sep=""))        
    ind=get.var.ncdf(nc,"ind")

    calc_global_dur(dat=dat,ind=ind,trash=trash,filename=paste("../data/",trendID,"/",additional_style,"/duration/",trendID,trend_style,dataset,additional_style,"_duration_year.nc",sep=""),states=stateIndeces)

    nc=open.ncdf(paste("../data/",trendID,"/",additional_style,"/duration/",trendID,trend_style,dataset,additional_style,"_duration_year.nc",sep=""))
    dur=get.var.ncdf(nc,"dur")
    dur_mid=get.var.ncdf(nc,"dur_mid")
    
    duration_seasons(dur,dur_mid,season=c(60,151),filename=paste("../data/",trendID,"/",additional_style,"/duration/",trendID,trend_style,dataset,additional_style,"_duration_spring.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(152,243),filename=paste("../data/",trendID,"/",additional_style,"/duration/",trendID,trend_style,dataset,additional_style,"_duration_summer.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(244,334),filename=paste("../data/",trendID,"/",additional_style,"/duration/",trendID,trend_style,dataset,additional_style,"_duration_autumn.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(335,424),filename=paste("../data/",trendID,"/",additional_style,"/duration/",trendID,trend_style,dataset,additional_style,"_duration_winter.nc",sep=""))
}

master_analyse_duration <- function(yearPeriod,trendID,states,seasons=c("spring","summer","autumn","winter","year")
    ,trend_style="_mean",additional_style=""){
    # analyse duration periods 2 states
    for (season in seasons){
        print(season)
        nc=open.ncdf(paste("../data/",trendID,"/",states,"_states",trend_style,additional_style,"/duration/",trendID,"_duration_",states,"s_",season,".nc",sep=""))
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")
        duration_analysis(dur,dur_mid,filename=paste("../data/",trendID,"/",states,"_states",trend_style,additional_style,"/duration/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_duration_",states,"s_analysis_",season,".nc",sep=""),
            season=season,yearPeriod)#,stations=seq(487,489,1))#,stations=134:136)
    }
}

master_duration_distribution <- function(yearPeriod,trendID,states,seasons=c("spring","summer","autumn","winter","year")
    ,trend_style="_mean",additional_style=""){
    # analyse duration periods 2 states
    for (season in seasons){
        print(season)
        nc=open.ncdf(paste("../data/",trendID,"/",states,"_states",trend_style,additional_style,"/duration/",trendID,"_duration_",states,"s_",season,".nc",sep=""))
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")
        duration_distribution(dur,dur_mid,filename=paste("../data/",trendID,"/",states,"_states",trend_style,additional_style,"/duration/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_duration_",states,"s_distribution_",season,".nc",sep=""),
            season=season,yearPeriod)#,stations=seq(487,489,1))#,stations=134:136)
    }
}

master_regional_trend <- function(yearPeriod,region_name,trendID,trend_style="_mean",additional_style=""){
    source("functions_regional.r")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    regional_analysis(dat=dat,yearPeriod,filepath=paste("../data/",trendID,"/2_states",trend_style,additional_style,"/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",sep=""),region_name=region_name)
}

master_regional_climatology <- function(yearPeriod,region_name,trendID,trend_style="_mean",additional_style=""){
    source("functions_regional.r")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    #regional_climatology(trendID=trendID,dat=dat,yearPeriod=yearPeriod,region_name=region_name)
    plot_regional_distributions(trendID,dat,yearPeriod,region_name)
    plot_regional_boxplots(trendID,dat,yearPeriod,region_name)
}




full_2states <- function(nday,nyr,trendID,trend_style="_mean",additional_style="_median",dataset="_TX"){
    #complete 2 states analysis 
    trendID=paste(nday,"_",nyr,sep="")

    #master_trend(nday,nyr,trendID,trend_style=trend_style,dataset=dataset)

    #master_state_attribution(nday,nyr,trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)

    #master_markov(nday,nyr,trendID,states=2,transition_names=c("cc wc cw ww"),trend_style=trend_style,additional_style=additional_style)

    points=c(1950,2014,1950,1980,1980,2014)
    for (i in 1:3){
        period=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        #master_analyse_markov(yearPeriod=period,trendID,states=2,transition_names=c("cc wc cw ww"),trend_style=trend_style,additional_style=additional_style)
    }

    #master_duration(nday,nyr,trendID,2,trend_style=trend_style,additional_style=additional_style)
    for (i in 1:3){
        period=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        #master_analyse_duration(yearPeriod=period,trendID,states=2,trend_style=trend_style,additional_style=additional_style)
        #master_duration_distribution(yearPeriod=period,trendID,states=2,trend_style=trend_style,additional_style=additional_style)
    }
}


#full_3states(91,5)
#full_3states(91,3)
#full_2states(91,5,trend_style="_mode",additional_style="_TX")
#full_2states(91,5,trend_style="_mean",dataset="_TX",additional_style="_median")



#master_trend(trendID="61_5",nday=61,nyr=5)
#master_markov(nday=61,nyr=5,trendID="61_5",states=2,transition_names=c("cc wc cw ww"))

#master_trend_control(trendID="61_5",states=2)
#master_trend_control(trendID="91_5",states=2,trend_style="_median_TX")
#master_trend_control(trendID="91_5",states=2,trend_style="_mode_TX")
#master_trend_control(trendID="91_5",states=2,trend_style="_mean_TX")
master_trend_control(trendID="91_5",trend_style="_mean",dataset="_TX",additional_style="_median")

#master_regional_climatology(yearPeriod=c(1950,2014),region_name="7rect",trendID="91_5")
#master_regional_climatology(yearPeriod=c(1950,2014),region_name="midlat",trendID="91_5")
#master_regional_climatology(yearPeriod=c(1950,1980),region_name="7rect",trendID="91_5")
#master_regional_climatology(yearPeriod=c(1980,2014),region_name="7rect",trendID="91_5")

if (1==2){
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    points=c(1950,2014,1950,1980,1980,2014)
    for (i in 1:3){
        period=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        print(paste(period,sep="-"))
        end_aussage(dat=dat,yearPeriod=paste(period[1],period[2],sep="-"),states=2,trendID="91_5")
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