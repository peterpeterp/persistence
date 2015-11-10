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

master_trend_control <- function(trendID,trend_style="_mean",dataset="_TX",additional_style=""){
    # trend control
    source("trend_view.r")
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_state_ind","_median",".nc",sep=""))
    ind=get.var.ncdf(nc,"ind")
    trend_control_warm_days(dat,ind,filename=paste("../data/",trendID,"/",dataset,additional_style,"/sonstiges/",trendID,trend_style,dataset,"_5seasons_warme_tage_control",additional_style,"_year",".txt",sep=""))   
    nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_state_ind","_seasonal_median",".nc",sep=""))
    ind=get.var.ncdf(nc,"ind")
    trend_control_warm_days(dat,ind,filename=paste("../data/",trendID,"/",dataset,additional_style,"/sonstiges/",trendID,trend_style,dataset,"_5seasons_warme_tage_control",additional_style,"_4seasons",".txt",sep=""))  
}

master_trend <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style=""){
    # calculate trend
    # choice between mean, median and estimated mode is possible
    if (trend_style=="_mean"){procedure=r_calc_runmean_2D}
    if (trend_style=="_median"){procedure=r_calc_runmedian_2D}
    if (trend_style=="_mode"){procedure=r_calc_runmode_2D}
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    trend=calc_trend(dat,paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""),nday,nyr,procedure=procedure)
}

master_seasonal_median_on_detrended <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style=""){
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    trend=trend_load(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
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
    filename=paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_seasonal_median",".nc",sep="")
    nc = create.ncdf(filename,seasonal_med)
    put.var.ncdf(nc,seasonal_med,seasonalmedian)
    close.ncdf(nc)
}

master_state_attribution <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style=""){
    # calculate persistence 2 states
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    trend=trend_load(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    detrended=dat$tas-trend
    for (q in 1:1319){
       detrended[q,,]=detrended[q,,]-median(detrended[q,,],na.rm=TRUE) 
    }
    
    per=state_attribution(dat,detrended,nday,nyr,
        filename=paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_state_ind","_median",".nc",sep=""))
}


master_state_attribution_daily_median <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style=""){
    # calculate persistence 2 states
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    trend=trend_load(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_seasonal_median",".nc",sep=""))
    daily_median=get.var.ncdf(nc,"_seasonal_median")
    detrended=dat$tas-trend
    for (year in 1:65){
        detrended[,,year]=detrended[,,year]-daily_median
    }
    per=state_attribution(dat,detrended,nday,nyr,
        filename=paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_state_ind","_seasonal_median",".nc",sep=""))
}


master_duration <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style=""){
    # calculate duration periods 2 states
    trash=((nyr-1)/2*365+(nday-1))
    stateIndeces=c(-1,1)

    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))

    # calculate annual duration
    nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_state_ind","_median",".nc",sep=""))        
    ind=get.var.ncdf(nc,"ind")
    calc_global_dur(dat=dat,ind=ind,trash=trash,filename=paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_year.nc",sep=""),states=stateIndeces)

    # calculate seasonal duration
    nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_state_ind","_seasonal_median",".nc",sep=""))        
    ind=get.var.ncdf(nc,"ind")
    calc_global_dur(dat=dat,ind=ind,trash=trash,filename=paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_4seasons.nc",sep=""),states=stateIndeces)

    # open seasonal duration and seperate into individual files
    nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_4seasons.nc",sep=""))
    dur=get.var.ncdf(nc,"dur")
    dur_mid=get.var.ncdf(nc,"dur_mid")
    
    duration_seasons(dur,dur_mid,season=c(60,151),filename=paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_MAM.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(152,243),filename=paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_JJA.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(244,334),filename=paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_SON.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(335,424),filename=paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_DJF.nc",sep=""))
}

master_analyse_duration <- function(yearPeriod,trendID,seasons=c("MAM","JJA","SON","DJF","year","4seasons"),dataset="_TX",additional_style=""){
    library(quantreg)
    # analyse duration periods 2 states
    for (season in seasons){
        nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_",season,".nc",sep=""))
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")
        duration_analysis(dur,dur_mid,filename=paste("../data/",trendID,"/",dataset,additional_style,"/duration/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,dataset,"_duration_analysis_",season,".nc",sep=""),season=season,yearPeriod)
    }
}

master_duration_distribution <- function(yearPeriod,trendID,seasons=c("MAM","JJA","SON","DJF","year","4seasons"),dataset="_TX",additional_style=""){
    library(moments)
    # analyse duration periods 2 states
    for (season in seasons){
        print(season)
        nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_",season,".nc",sep=""))
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")
        duration_distribution(dur,dur_mid,filename=paste("../data/",trendID,"/",dataset,additional_style,"/duration/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,dataset,"_duration_distribution_ana_",season,".nc",sep=""),
            season=season,yearPeriod)#,stations=seq(487,489,1))#,stations=134:136)
    }
}

master_regional_climatology <- function(region_name,trendID,dataset="_TX",additional_style=""){
    library(quantreg)
    library(moments)
    source("functions_regional.r")
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))

    points=c(1950,2014,1950,1980,1980,2014)
    #points=c(1950,1980,1980,2014)
    for (i in 1:3){
        yearPeriod=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        print(yearPeriod)
        #regional_climatology(dat=dat,yearPeriod=yearPeriod,region_name=region_name,trendID=trendID,additional_style=additional_style,dataset=dataset)
        #plot_regional_distributions(dat,yearPeriod,region_name,trendID,dataset=dataset,additional_style=additional_style)
        #regional_quantiles_fits(dat=dat,yearPeriod=yearPeriod,region_name=region_name,trendID=trendID,dataset=dataset,additional_style=additional_style)
        plot_regional_fit_parameters(dat=dat,yearPeriod=yearPeriod,region_name=region_name,trendID=trendID,dataset=dataset,additional_style=additional_style)
        plot_regional_boxplots(dat=dat,yearPeriod=yearPeriod,region_name=region_name,trendID=trendID,dataset=dataset,additional_style=additional_style)
    }

    #plot_regional_boxplots_vergleich(dat=dat,yearPeriod1=c(1950,1980),yearPeriod2=c(1980,2014),region_name=region_name,trendID=trendID,dataset=dataset,additional_style=additional_style)
}


full <- function(nday,nyr,trend_style="_mean",dataset="_TX",additional_style=""){
    #complete 2 states analysis 
    trendID=paste(nday,"_",nyr,sep="")

    #master_duration(nday,nyr,trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)
    points=c(1950,2014,1950,1980,1980,2014)
    points=c(1950,1980,1980,2014)
    for (i in 1:3){
        period=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        #master_analyse_duration(yearPeriod=period,trendID=trendID,dataset=dataset,additional_style=additional_style)
        #master_analyse_duration(yearPeriod=period,trendID=trendID,dataset=dataset,additional_style=additional_style,seasons=c("year"))
        master_duration_distribution(yearPeriod=period,trendID=trendID,dataset=dataset,additional_style=additional_style)
    }
}

# init
nday=91
nyr=5
trendID=paste(nday,"_",nyr,sep="")
dataset="_TX"
trend_style="_mean"
additional_style=""

#master_trend(nday,nyr,trendID,trend_style=trend_style,dataset=dataset)


#master_seasonal_median_on_detrended(nday=nday,nyr=nyr,trendID=trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)

#master_state_attribution(nday,nyr,trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)
#master_state_attribution_daily_median(nday,nyr,trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)


#full(nday,nyr,trend_style="_mean",dataset="_TX",additional_style="")



#master_trend_control(trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)
master_regional_climatology(region_name="7rect",trendID=trendID,dataset=dataset,additional_style=additional_style)


