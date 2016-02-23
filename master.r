#!/home/pepflei/R/bin/Rscript

###################################################################
# Persistence in Temperature time series
# Author: Peter Pfleiderer
# Institution: Potsdam Institute for Climate Impact Research (PIK)
# Year: 2015
###################################################################

master_nas <- function(){
    # count nas
    find_nas(dat)    
}

master_trend_control <- function(){
    # trend control
    source("trend_view.r")
    nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_state_ind","_median",".nc",sep=""))
    ind=var.get.nc(nc,"ind")
    trend_control_warm_days(dat,ind,filename=paste("../data/",dataset,additional_style,"/",trendID,"/sonstiges/",trendID,trend_style,dataset,"_5seasons_warme_tage_control",additional_style,"_year",".txt",sep=""))   
    nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_state_ind","_seasonal_median",".nc",sep=""))
    ind=var.get.nc(nc,"ind")
    trend_control_warm_days(dat,ind,filename=paste("../data/",dataset,additional_style,"/",dataset,additional_style,"/sonstiges/",trendID,trend_style,dataset,"_5seasons_warme_tage_control",additional_style,"_4seasons",".txt",sep=""))  
}

master_trend <- function(){
    # calculate trend
    # choice between mean, median and estimated mode is possible
    if (trend_style=="_mean"){procedure=r_calc_runmean_2D}
    if (trend_style=="_median"){procedure=r_calc_runmedian_2D}
    if (trend_style=="_mode"){procedure=r_calc_runmode_2D}
    print(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    trend=calc_trend(dat,paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""),nday,nyr,procedure=procedure)
}

master_seasonal_median_on_detrended <- function(){
    print(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    trend=trend_load(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    detrended=dat$tas-trend

    ntot = 1319
    seasonal_median=array(0,dim=c(ntot,365))

    seasonStart=c(60,152,244,335,1)
    seasonStop=c(151,243,334,424,365)
    for (sea in 1:4){
        cat("-")
        for (q in 1:ntot){
            if (sea==4){
                z=c(detrended[q,1:(seasonStop[sea]-365),],detrended[q,seasonStart[sea]:365,])
                seasonal_median[q,1:(seasonStop[sea]-365)]=array(median(z,na.rm=TRUE),(seasonStop[sea]-365))
                seasonal_median[q,seasonStart[sea]:365]=array(median(z,na.rm=TRUE),(365-seasonStart[sea]+1))

            }
            else {
                z=detrended[q,seasonStart[sea]:seasonStop[sea],]
                seasonal_median[q,seasonStart[sea]:seasonStop[sea]]=array(median(z,na.rm=TRUE),(seasonStop[sea]-seasonStart[sea]+1))
            }   
        }
    }
    print(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_seasonal_median",".nc",sep=""))
    nc_out <- create.nc(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_seasonal_median",".nc",sep=""))

    dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
    dim.def.nc(nc_out,"day",dimlength=365,unlim=FALSE)

    var.def.nc(nc_out,"seasonal_median","NC_DOUBLE",c(0,1))
    att.put.nc(nc_out, "seasonal_median", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "seasonal_median", "dim_explanation", "NC_CHAR", "ID-day")
    att.put.nc(nc_out, "seasonal_median", "explanation", "NC_CHAR","for each season and each gridpoint the median of the detrended timeseries is computed. this seasonal median is stored for each day")

    var.put.nc(nc_out,"seasonal_median",seasonal_median)
    close.nc(nc_out)  
}

master_state_attribution <- function(){
    # calculate persistence 2 states
    print(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    trend=trend_load(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    print(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_seasonal_median",".nc",sep=""))
    nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_seasonal_median",".nc",sep=""))
    seasonal_median=var.get.nc(nc,"seasonal_median")
    detrended=dat$tas-trend
    for (year in 1:65){
        detrended[,,year]=detrended[,,year]-seasonal_median
    }
    per=state_attribution(dat,detrended,nday,nyr,filename=paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_state_ind.nc",sep=""))
}

master_duration <- function(){
    # calculate duration periods 2 states
    trash=((nyr-1)/2*365+(nday-1))
    stateIndeces=c(-1,1)

    # calculate annual period durations
    nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_state_ind.nc",sep=""))
    ind=var.get.nc(nc,"ind")
    cat("\nidentifying persistent periods:")
    calc_global_dur(ind=ind,trash=trash,filename=paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_duration_4seasons.nc",sep=""),states=stateIndeces)

    # open annual duration and seperate into individual files
    nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_duration_4seasons.nc",sep=""))
    dur=var.get.nc(nc,"dur")
    dur_mid=var.get.nc(nc,"dur_mid")

    cat("\nattributing periods to seasons:")    
    duration_seasons(dur,dur_mid,season=c(60,151),filename=paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_duration_MAM.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(152,243),filename=paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_duration_JJA.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(244,334),filename=paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_duration_SON.nc",sep=""))
    duration_seasons(dur,dur_mid,season=c(335,424),filename=paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_duration_DJF.nc",sep=""))
}

master_duration_analysis <- function(ID_select=1:1319){
    yearLimits=c(1950,2014,1980,2014)
    for (i in 1:2){
        yearPeriod=c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        print(yearPeriod)

        print("others")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(1,0,0,0,0,0,0,0))

        print("quant")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,1,0,0,0,0,0,0),noise_level=c(0,0.000001))
        
        print("fit")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,0,1,0,0,0,0),plot_select=c(NA),ID_select=ID_select,add_name="2expo_4:100",xStart=4,write=TRUE)
    }
}

master_regional_climatology <- function(region_name="7rect",region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA"),ID_select=1:7,plot_select=1:7,ID_length=7){

    
    #regional_attribution(region_name=region_name,trendID=trendID,dataset=dataset,additional_style=additional_style)
    
    ID_name=paste("_",region_name,sep="")
    yearLimits=c(1950,2014,1980,2014)
    for (i in 1:2){
        yearPeriod=c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        period=paste(yearPeriod[1],"-",yearPeriod[2],sep="")
        print(yearPeriod)

        print("others")
        #duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(1,0,0,0,0,0,0,0),ID_name=ID_name,ID_select=ID_select,ID_names=region_names,ID_length=ID_length,folder=paste("/regional/",region_name,"/",sep=""))

        print("quant")
        #duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,1,0,0,0,0,0,0),noise_level=c(0,0.000001),ID_name=ID_name,ID_select=ID_select,ID_names=region_names,ID_length=ID_length,folder=paste("/regional/",region_name,"/",sep=""))
        
        print("fit")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,0,1,0,0,0,0),add_name="2expo_4:100",xStart=4,ID_name=ID_name,ID_select=ID_select,plot_select=plot_select,ID_names=region_names,ID_length=ID_length,folder=paste("/regional/",region_name,"/",sep=""))        

        #plot_regional_boxplots(period=paste(yearPeriod[1],"-",yearPeriod[2],sep=""),region_name=region_name,region_names=region_names,trendID=trendID,dataset=dataset,additional_style=additional_style)
        #plot_regional_fit_parameters(period=period,trendID=trendID,additional_style=additional_style,dataset=dataset,region_name=region_name,fit_style="_fit_2expo_b1>b2_5-10")
        write_regional_fit_table(trendID=trendID,region_name=region_name,region_names=region_names,ID_select=ID_select,fit_style="2expo_4:100",period=period)

    }
    #plot_regional_boxplots_vergleich(period1="1950-1980",period2="1980-2014",region_name=region_name,trendID=trendID,additional_style=additional_style,dataset=dataset)
    #plot_regional_fit_vergleich(period1="1950-1980",period2="1997-2014",fit_style="_fit_2expo_b1>b2_5-15",region_name=region_name,region_names=region_names,trendID=trendID,additional_style=additional_style,dataset=dataset)
}


###################################################################
# init: loading sources, setting variables ....
###################################################################

master_init <- function(){
    source("functions_support.r")
    source("functions_duration.r")
    source("functions_regional.r")
    source("analysis_tools.r")
    source("write.r")
    source("load.r")

    library(moments)
    library(quantreg)
    library(stats4)
    library(RNetCDF)
    library(SDMTools)


    nday<<-91
    nyr<<-5
    trendID<<-paste(nday,"_",nyr,sep="")
    dataset<<-"_TMean"
    trend_style<<-"_mean"
    additional_style<<-""
    dat<<-dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))


    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    state_names<<-c("cold","warm")
}

###################################################################
# basic analysis
###################################################################
master_init()

#master_trend()
#master_seasonal_median_on_detrended()
#master_state_attribution()
#master_duration()

###################################################################
# fits, quantiles etc
###################################################################

#master_duration_analysis()


###################################################################
# regional commands
###################################################################

master_regional_climatology(region_name="ward20",region_names=1:20,ID_select=1:20,plot_select=1:20,ID_length=20)
asasd

taus=c(0,0.05,0.25,0.5,0.75,0.95,1)
trend_analysis()
adsas

#master_regional_climatology(region_name="mid_lat_belt",region_names=c("mid-lat","polar"),trendID=trendID,dataset=dataset,additional_style=additional_style)

#ID_select=c(1,2,3,4,5,6,9,10,11,12,13,14,16,17,18,19,20,21,22,23,24,25,26)
#master_regional_climatology(trendID=trendID,region_name="srex",region_names=1:26,ID_select=ID_select,plot_select=ID_select,ID_length=26)

#plot_fits_for_region(period="1950-2014",trendID=trendID,dataset=dataset,fit_style="2expo_thresh_5-15",reg=13,region_name="7rect",ID_select=ID_select)


#regional_quantiles_fits(dat=dat,yearPeriod=c(1950.2014),region_name="mid_lat_belt",trendID=trendID,dataset=dataset,additional_style=additional_style,plot=TRUE,season_auswahl=c(2),write=FALSE,add_name="_4sea",region_names=c("661","488"),q=c(661,488))

#master_regional_climatology(region_name="kmeans_ma5_grou6",trendID=trendID,dataset=dataset,additional_style=additional_style)

#dat<<-dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
#IDregions=points_to_regions(dat,region_name)
#ID_select=which(IDregions==reg)

period="1950-2014"
nGroup=7
add_name="2Distr_3Red"
region_name=paste(add_name,"_",nGroup,sep="")
region_name="7rect"

#duration_analysis(yearPeriod=c(1950,2014),trendID=trendID,dataset=dataset,option=c(0,0,0,1,0,0,0,0),add_name="2expo_thresh_5-15",folder="/regional/",ID_name=paste("_",region_name,sep=""),ID_select=1:nGroup,plot_select=1:nGroup,ID_names=1:nGroup,ID_length=nGroup)



#nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/nearest_neighbors/",period,"/",trendID,"_",period,"_",region_name,".nc",sep=""))
#IDregions=var.get.nc(nc,"attribution")


IDregions=points_to_regions(dat,region_name)
#ID_select=which(IDregions==reg)


for (reg in 1:7){
    plot_fits_for_region(reg=reg,IDregions=c("from polygons"),region_name=region_name)
}



#master_duration_analysis <- function(trendID){
#full(nday,nyr,trend_style=trend_style,dataset=dataset,additional_style=additional_style)

#master_trend_control(trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)

#sdsdsdsds



#duration_analysis(yearPeriod=c(1980,2014),trendID=trendID,dataset=dataset,option=c(0,0,0,1,1,0),add_name="_testin2",ID_select=c(460,466,554),plot_select=c(460,466,554),ID_length=1319,write=FALSE)
#duration_analysis(yearPeriod=c(1950,1980),trendID=trendID,dataset=dataset,option=c(0,0,1,1,1,0),add_name="_testin2",ID_select=c(460,466,554),plot_select=c(460,466,554),ID_length=1319,write=FALSE)
#duration_analysis(yearPeriod=c(1950,2014),trendID=trendID,dataset=dataset,option=c(0,0,0,1,0,0,0,0),plot_select=c(NA),ID_select=1:1319,add_name="2expo_4:100",xStart=4,write=TRUE)

#adsas

#master_duration_analysis(trendID=trendID)

#master_duration(nday,nyr,trendID,trend_style=trend_style,dataset=dataset,additional_style=additional_style)