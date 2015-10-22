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

master_trend_control <- function(trendID,trend_style="_mean",dataset="_TX",additional_style="_median"){
    # trend control
    source("trend_view.r")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    nc=open.ncdf(paste("../data/",trendID,"/",trendID,trend_style,dataset,"_state_ind",additional_style,".nc",sep=""))
    ind=get.var.ncdf(nc,"ind")
    trend_control_warm_days(dat,ind,filename=paste("../data/",trendID,"/",additional_style,"/sonstiges/",trendID,trend_style,dataset,"_5seasons_warme_tage_control",additional_style,".txt",sep=""))   
}

master_trend <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX"){
    # calculate trend
    # choice between mean, median and estimated mode is possible
    if (trend_style=="_mean"){procedure=r_calc_runmean_2D}
    if (trend_style=="_median"){procedure=r_calc_runmedian_2D}
    if (trend_style=="_mode"){procedure=r_calc_runmode_2D}
    if (additional_style=="_TX"){dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")}
    if (additional_style=="_TN"){dat=dat_load("../data/HadGHCND_TN_data3D.day1-365.1950-2014.nc")}    
    trend=calc_trend(dat,paste("../data/",trendID,"/",trendID,"_trend",trend_style,additional_style,".nc",sep=""),nday,nyr,procedure=procedure)
}

master_runmedian_on_detrended <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style="_run_median"){
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    trend=trend_load(paste("../data/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    detrended=dat$tas-trend

    trash = ((nyr-1)/2*365+(nday-1))
    ntot = length(dat$ID)
    runmedian=dat$tas*NA
    for (q in 1:ntot) {
        temp = r_calc_runmedian_2D(detrended[q,,],nday=nday,nyr=nyr)
        temp[1:trash]=NA
        temp[(length(dat$time)-trash):length(dat$time)]=NA
        runmedian[q,,]=temp
    }
    trend_write(paste("../data/",trendID,"/",trendID,trend_style,dataset,additional_style,".nc",sep=""),dat,runmedian)
}

master_daily_median_on_detrended <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style="_daily_median"){
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    trend=trend_load(paste("../data/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    detrended=dat$tas-trend

    ntot = 1319
    dailymedian=array(NA,dim=c(ntot,365))
    for (day in 1:365){
        cat("-")
        for (q in 1:ntot){
            if ((day-(nday-1)/2)>0 & (day+(nday-1)/2)<366){
                z=c(detrended[q,(day-(nday-1)/2):(day+(nday-1)/2),])
            }
            if ((day-(nday-1)/2)<1){
                z=c(detrended[q,(365+day-(nday-1)/2):365,],detrended[q,1:(day+(nday-1)/2),])
            }
            if ((day+(nday-1)/2)>365){
                z=c(detrended[q,(day-(nday-1)/2):365,],detrended[q,1:((day+(nday-1)/2)-365),])
            }
            if (which(length(!is.na(z))>100)){
                dailymedian[q,day]=median(z,na.rm=TRUE)
            }
        }
    }

    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)

    daily_med <- var.def.ncdf(name="_daily_median",units="medain value for each day",dim=list(ID,day), missval=-9999.0)
    filename=paste("../data/",trendID,"/",trendID,trend_style,dataset,additional_style,".nc",sep="")
    nc = create.ncdf(filename,daily_med)
    put.var.ncdf(nc,daily_med,dailymedian)
    close.ncdf(nc)
}

master_seasonal_median_on_detrended <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style="_seasonal_median"){
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    trend=trend_load(paste("../data/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    detrended=dat$tas-trend

    ntot = 1319
    seasonalmedian=array(0,dim=c(ntot,365))

    seasonStart=c(60,152,244,335,1)
    seasonStop=c(151,243,334,424,365)
    for (sea in 1:4){
        cat("-")
        for (q in 1:ntot){
            if (sea==4){
                z=c(detrended[q,1:(seasonStop[sea]-365),],detrended[q,seasonStart[sea]:365,])
                seasonalmedian[q,1:(seasonStop[sea]-365)]=array(median(z,na.rm=TRUE),(seasonStop[sea]-365))
                seasonalmedian[q,seasonStart[sea]:365]=array(median(z,na.rm=TRUE),(365-seasonStart[sea]+1))

            }
            else {
                z=detrended[q,seasonStart[sea]:seasonStop[sea],]
                seasonalmedian[q,seasonStart[sea]:seasonStop[sea]]=array(median(z,na.rm=TRUE),(seasonStop[sea]-seasonStart[sea]+1))
            }   
        }
    }

    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)

    seasonal_med <- var.def.ncdf(name="_seasonal_median",units="medain of season value for each day",dim=list(ID,day), missval=-9999.0)
    filename=paste("../data/",trendID,"/",trendID,trend_style,dataset,additional_style,".nc",sep="")
    nc = create.ncdf(filename,seasonal_med)
    put.var.ncdf(nc,seasonal_med,seasonalmedian)
    close.ncdf(nc)
}

master_state_attribution <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style="_median"){
    # calculate persistence 2 states
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend=trend_load(paste("../data/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    detrended=dat$tas-trend
    detrended=detrended-median(detrended,na.rm=TRUE)
    per=state_attribution(dat,detrended,nday,nyr,
        filename=paste("../data/",trendID,"/",trendID,trend_style,dataset,"_state_ind",additional_style,".nc",sep=""))
}

master_state_attribution_2_trends <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style="_run_median"){
    # calculate persistence 2 states
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend1=trend_load(paste("../data/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    trend2=trend_load(paste("../data/",trendID,"/",trendID,trend_style,dataset,"_run_median.nc",sep=""))
    detrended=dat$tas-trend1-trend2
    per=state_attribution(dat,detrended,nday,nyr,
        filename=paste("../data/",trendID,"/",trendID,trend_style,dataset,"_state_ind",additional_style,".nc",sep=""))
}

master_state_attribution_daily_median <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style="_daily_median"){
    # calculate persistence 2 states
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend=trend_load(paste("../data/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    nc=open.ncdf(paste("../data/",trendID,"/",trendID,trend_style,dataset,additional_style,".nc",sep=""))
    daily_median=get.var.ncdf(nc,additional_style)
    detrended=dat$tas-trend
    for (year in 1:65){
        detrended[,,year]=detrended[,,year]-daily_median
    }
    per=state_attribution(dat,detrended,nday,nyr,
        filename=paste("../data/",trendID,"/",trendID,trend_style,dataset,"_state_ind",additional_style,".nc",sep=""))
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
    library(quantreg)
    # analyse duration periods 2 states
    for (season in seasons){
        print(season)
        nc=open.ncdf(paste("../data/",trendID,"/",additional_style,"/duration/",trendID,trend_style,dataset,additional_style,"_duration_",season,".nc",sep=""))
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")
        duration_analysis(dur,dur_mid,filename=paste("../data/",trendID,"/",additional_style,"/duration/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,trend_style,dataset,additional_style,"_duration_analysis_",season,".nc",sep=""),
            season=season,yearPeriod)#,stations=seq(487,489,1))#,stations=134:136)
    }
}

master_duration_distribution <- function(yearPeriod,trendID,states,seasons=c("spring","summer","autumn","winter","year")
    ,trend_style="_mean",additional_style=""){
    # analyse duration periods 2 states
    for (season in seasons){
        print(season)
        nc=open.ncdf(paste("../data/",trendID,"/",additional_style,"/duration/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_duration_",states,"s_analysis_",season,".nc",sep=""))
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")
        duration_distribution(dur,dur_mid,filename=paste("../data/",trendID,"/",additional_style,"/duration/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_duration_",states,"s_distribution_",season,".nc",sep=""),
            season=season,yearPeriod)#,stations=seq(487,489,1))#,stations=134:136)
    }
}

master_regional_trend <- function(yearPeriod,region_name,trendID,trend_style="_mean",additional_style=""){
    source("functions_regional.r")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    regional_analysis(dat=dat,yearPeriod,filepath=paste("../data/",trendID,"/2_states",trend_style,additional_style,"/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",sep=""),region_name=region_name)
}

master_regional_climatology <- function(yearPeriod,region_name,trendID,trend_style="_mean",dataset="_TX",additional_style="_seasonal_median"){
    source("functions_regional.r")
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    regional_climatology(trendID=trendID,additional_style=additional_style,dat=dat,yearPeriod=yearPeriod,region_name=region_name)
    #plot_regional_distributions(trendID,dat,yearPeriod,region_name,additional_style=additional_style)
    plot_regional_boxplots(trendID,dat,yearPeriod,region_name,additional_style=additional_style)
}

master_endaussage <- function(){
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    points=c(1950,2014,1950,1980,1980,2014)
    for (i in 1:3){
        period=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        print(paste(period,sep="-"))
        end_aussage(dat=dat,yearPeriod=paste(period[1],period[2],sep="-"),states=2,trendID="91_5")
    }
}


full_2states <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style="_seasonal_median"){
    #complete 2 states analysis 
    trendID=paste(nday,"_",nyr,sep="")

    #master_duration(nday,nyr,trendID,2,trend_style=trend_style,additional_style=additional_style)
    points=c(1950,2014,1950,1980,1980,2014)
    for (i in 1:3){
        period=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        #master_analyse_duration(yearPeriod=period,trendID=trendID,states=2,trend_style=trend_style,additional_style=additional_style)
        #master_duration_distribution(yearPeriod=period,trendID,states=2,trend_style=trend_style,additional_style=additional_style)
        master_regional_climatology(yearPeriod=period,region_name="7rect",trendID=trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)

    }
}

# init
nday=91
nyr=5
trendID="91_5"
trend_style="_mean"
dataset="_TX"

additional_style="_seasonal_median"

#master_trend(nday,nyr,trendID,trend_style=trend_style,dataset=dataset)
#master_runmedian_on_detrended(nday=nday,nyr=nyr,trendID=trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)
#master_daily_median_on_detrended(nday=nday,nyr=nyr,trendID=trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)
#master_seasonal_median_on_detrended(nday=nday,nyr=nyr,trendID=trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)

#master_state_attribution_2_trends(nday,nyr,trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)
#master_state_attribution(nday,nyr,trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)
#master_state_attribution_daily_median(nday,nyr,trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)


full_2states(91,5,trend_style="_mean",dataset="_TX",additional_style="_seasonal_median")



#master_trend_control(trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)



