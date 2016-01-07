#!/home/pepflei/R/bin/Rscript

###################################################################
# Persistence in Temperature time series
# Author: Peter Pfleiderer
# Institution: Potsdam Institute for Climate Impact Research (PIK)
# Year: 2015
###################################################################

master_nas <- function(){
    # count nas
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2011.nc")
    find_nas(dat)    
}

master_trend_control <- function(trendID,trend_style="_mean",dataset="_TX",additional_style=""){
    # trend control
    source("trend_view.r")
    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_state_ind","_median",".nc",sep=""))
    ind=var.get.nc(nc,"ind")
    trend_control_warm_days(dat,ind,filename=paste("../data/",trendID,"/",dataset,additional_style,"/sonstiges/",trendID,trend_style,dataset,"_5seasons_warme_tage_control",additional_style,"_year",".txt",sep=""))   
    nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_state_ind","_seasonal_median",".nc",sep=""))
    ind=var.get.nc(nc,"ind")
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

    day <- dim.def.nc("day", units="d",vals=1:365, unlim=FALSE)
    ID <- dim.def.nc("ID",units="ID",vals=1:ntot, unlim=FALSE)

    seasonal_med <- var.def.nc(name="_seasonal_median",units="medain of season value for each day",dim=list(ID,day), missval=-9999.0)
    filename=paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_seasonal_median",".nc",sep="")
    nc = create.nc(filename,seasonal_med)
    put.var.nc(nc,seasonal_med,seasonalmedian)
    close.nc(nc)
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

master_duration <- function(nday,nyr,trendID,trend_style="_mean",dataset="_TX",additional_style=""){
    # calculate duration periods 2 states
    trash=((nyr-1)/2*365+(nday-1))
    stateIndeces=c(-1,1)

    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))

    # calculate seasonal duration
    nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/",trendID,trend_style,dataset,"_state_ind","_seasonal_median",".nc",sep=""))
    ind=var.get.nc(nc,"ind")
    cat("\nidentifying persistent periods:")
    calc_global_dur(dat=dat,ind=ind,trash=trash,filename=paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",trendID,dataset,"_duration_4seasons.nc",sep=""),states=stateIndeces)

    # open seasonal duration and seperate into individual files
    nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_4seasons.nc",sep=""))
    dur=var.get.nc(nc,"dur")
    dur_mid=var.get.nc(nc,"dur_mid")

    cat("\nattributing periods to seasons:")    
    duration_seasons(dur,dur_mid,season=c(60,151),filename=paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",trendID,dataset,"_duration_MAM.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(152,243),filename=paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",trendID,dataset,"_duration_JJA.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(244,334),filename=paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",trendID,dataset,"_duration_SON.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(335,424),filename=paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",trendID,dataset,"_duration_DJF.nc",sep=""))
}

master_duration_analysis <- function(trendID,dataset="_TMean",ID_select=1:1319){
    library(quantreg)
    source("analysis_tools.r")
    #dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))

    points=c(1950,2014,1950,1980,1980,2014)
    for (i in 1:3){
        yearPeriod=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        print(yearPeriod)
        print("others")
        #duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(1,0,0,0),add_name="_testin")
        print("quant")
        #duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,1,0,0,0,0),add_name="_testin",noise_level=c(0,0.000001))
        
        #duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,1,0,0,0,0),plot_select=c(488,661,553),add_name="expo")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,0,1,0,0,0),plot_select=c(NA),ID_select=ID_select,add_name="2expo_thresh_5-15",write=TRUE)
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,0,0,0,0,1),plot_select=c(NA),ID_select=ID_select,add_name="2expo_overlap",write=TRUE)
    }
}

master_regional_climatology <- function(trendID,ID_select=1:7,plot_select=1:7,dataset="_TMean",additional_style="",region_name="7rect",region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){
    library(quantreg)
    library(moments)
    
    #dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    #regional_distributions(dat=dat,region_name=region_name,trendID=trendID,dataset=dataset,additional_style=additional_style)
    
    ID_name=paste("_",region_name,sep="")
    points=c(1950,2014,1950,1980,1980,2014)
    #points=c(1950,2014)
    #points=c(1980,1997,1997,2014)
    for (i in 1:3){
        yearPeriod=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        period=paste(yearPeriod[1],"-",yearPeriod[2],sep="")
        print(yearPeriod)
        
        #duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(1,0,0,0,0,0,0),add_name="_testin",folder="/regional/",ID_name=ID_name,ID_select=1:26,plot_select=1:26)
        #duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,1,0,0,0,0),add_name="_testin",folder="/regional/",ID_name=ID_name,ID_length=7,noise_level=c(0,0001,0.001))
        #duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,1,0,0,0,0),add_name="_expo",folder="/regional/",ID_name=ID_name,ID_select=ID_select,plot_select=plot_select,ID_names=region_names,ID_length=26)
        #duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,0,1,0,0),add_name="combi_expo",folder="/regional/",ID_name="_7rect",ID_select=1:7,plot_select=1:7,ID_names=region_names,ID_length=7)
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,0,1,0,0,0),add_name="2expo_thresh_5-15",folder="/regional/",ID_name=ID_name,ID_select=ID_select,plot_select=plot_select,ID_names=region_names,ID_length=26)
        #duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,0,0,0,0,1),add_name="2expo_overlap",folder="/regional/",ID_name=ID_name,ID_select=ID_select,plot_select=plot_select,ID_names=region_names,ID_length=26)
        #duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,1,0,0,1),add_name="_expo_fixed_2expo",folder="/regional/",ID_name="_7rect",ID_select=1:7,plot_select=1:7,ID_names=region_names,ID_length=7)
        

        #plot_regional_boxplots(period=paste(yearPeriod[1],"-",yearPeriod[2],sep=""),region_name=region_name,region_names=region_names,trendID=trendID,dataset=dataset,additional_style=additional_style)
        #plot_regional_fit_parameters(period=period,trendID=trendID,additional_style=additional_style,dataset=dataset,region_name=region_name,fit_style="_fit_2expo_b1>b2_5-10")
        fit_info_to_map(trendID=trendID,region_name=region_name,region_names=region_names,ID_select=ID_select,fit_style1="2expo_b1>b2_5-15",fit_style2="2expo_overlap",period=period)
        #write_regional_fit_table(trendID=trendID,region_name=region_name,region_names=region_names,ID_select=ID_select,fit_style1="2expo_b1>b2_5-15",fit_style2="2expo_overlap",period=period)

    }
    #plot_regional_boxplots_vergleich(period1="1950-1980",period2="1980-2014",region_name=region_name,trendID=trendID,additional_style=additional_style,dataset=dataset)
    #plot_regional_fit_vergleich(period1="1950-1980",period2="1980-2014",fit_style="_fit_2expo_b1>b2_5-15",region_name=region_name,region_names=region_names,trendID=trendID,additional_style=additional_style,dataset=dataset)
    #plot_regional_fit_vergleich(period1="1950-1980",period2="1980-1997",fit_style="_fit_2expo_b1>b2_5-15",region_name=region_name,region_names=region_names,trendID=trendID,additional_style=additional_style,dataset=dataset)
    #plot_regional_fit_vergleich(period1="1950-1980",period2="1997-2014",fit_style="_fit_2expo_b1>b2_5-15",region_name=region_name,region_names=region_names,trendID=trendID,additional_style=additional_style,dataset=dataset)
}


full <- function(nday,nyr,trend_style="_mean",dataset="_TX",additional_style=""){
    #complete 2 states analysis 
    trendID=paste(nday,"_",nyr,sep="")

    master_duration(nday,nyr,trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)
    points=c(1950,2014,1950,1980,1980,2014)
    #points=c(1950,1980,1980,2014)
    for (i in 1:3){
        period=c(points[(2*(i-1)+1)],points[(2*(i-1)+2)])
        master_analyse_duration(yearPeriod=period,trendID=trendID,dataset=dataset,additional_style=additional_style)
        master_duration_distribution(yearPeriod=period,trendID=trendID,dataset=dataset,additional_style=additional_style)
    }
}

###################################################################
# init: loading sources, setting variables ....
###################################################################
source("functions_support.r")
source("functions_duration.r")
source("functions_regional.r")
source("analysis_tools.r")
source("write.r")
source("load.r")

nday=91
nyr=5
trendID=paste(nday,"_",nyr,sep="")
dataset="_TMean"
trend_style="_mean"
additional_style=""

###################################################################
# use some of the functions above
###################################################################

#master_trend(nday,nyr,trendID,trend_style=trend_style,dataset=dataset)


#master_seasonal_median_on_detrended(nday=nday,nyr=nyr,trendID=trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)

#master_state_attribution(nday,nyr,trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)
#master_state_attribution_daily_median(nday,nyr,trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)


#full(nday,nyr,trend_style=trend_style,dataset=dataset,additional_style=additional_style)


#regional_quantiles_fits(dat=dat,yearPeriod=c(1950.2014),region_name="mid_lat_belt",trendID=trendID,dataset=dataset,additional_style=additional_style,plot=TRUE,season_auswahl=c(2),write=FALSE,add_name="_4sea",region_names=c("661","488"),q=c(661,488))


#master_trend_control(trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)
#master_regional_climatology(region_name="7rect",trendID=trendID,dataset=dataset,additional_style=additional_style)
#master_regional_climatology(region_name="mid_lat_belt",region_names=c("mid-lat","polar"),trendID=trendID,dataset=dataset,additional_style=additional_style)

ID_select=c(1,2,3,4,5,6,9,10,11,12,13,14,16,17,18,19,20,21,22,23,24,25,26)
#master_regional_climatology(trendID=trendID,region_name="srex",region_names=1:26,ID_select=ID_select,plot_select=ID_select)

#duration_analysis(yearPeriod=c(1980,2014),trendID=trendID,dataset=dataset,option=c(0,0,0,1,1,0),add_name="_testin2",ID_select=c(460,466,554),plot_select=c(460,466,554),ID_length=1319,write=FALSE)
#duration_analysis(yearPeriod=c(1950,1980),trendID=trendID,dataset=dataset,option=c(0,0,1,1,1,0),add_name="_testin2",ID_select=c(460,466,554),plot_select=c(460,466,554),ID_length=1319,write=FALSE)
#duration_analysis(yearPeriod=c(1950,2014),trendID=trendID,dataset=dataset,option=c(0,0,0,1,0,0),add_name="_testin3",ID_select=c(2,4,30,52,374,506,31),plot_select=c(2,4,30,52,374,506,31),ID_length=1319,write=FALSE)

#master_duration_analysis(trendID=trendID)

#master_duration(nday,nyr,trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)

if (1==1){
    reg=13
    #dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    #IDregions=points_to_regions(dat,"srex")
    #ID_select=which(IDregions==reg)
    #print(ID_select)
    #print(length(ID_select))

    #duration_analysis(yearPeriod=c(1950,2014),trendID=trendID,dataset=dataset,option=c(0,0,0,1,0,0,0),add_name="2expo_thresh_5-15_reg_13",ID_select=ID_select,plot_select=ID_select,ID_length=1319)
    #plot_fits_for_region(period="1950-2014",trendID=trendID,dataset=dataset,fit_style="2expo_thresh_5-15",reg=reg,region_name="srex",ID_select=ID_select)
    distr_nearest_neighbours(period="1950-2014",trendID=trendID,dataset=dataset,fit_style="2expo_thresh_5-15",reg=reg,region_name="srex",ID_select=ID_select)

}