#!/home/pepflei/R/bin/Rscript

###################################################################
# Persistence in Temperature time series
# Author: Peter Pfleiderer
# Institution: Potsdam Institute for Climate Impact Research (PIK)
# Year: 2015
###################################################################

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
    if (trend_style=="_mean"){procedure=r_calc_runmean_2D}
    print(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    trend=calc_trend(dat,paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""),nday,nyr,procedure=procedure)
}

master_seasonal_median_on_detrended <- function(){
    # calculates the median for each season and each grid-point
    # this seasonal median is stored in aarray of 365 with the correpsonding seasonal median for each day
    print(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    trend=trend_load(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    detrended=dat$tas-trend

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
    # from temp anomalies running mean and seasonal median is subtracted
    # than negative days are cold (stored as -1), positive days are warm (stored as +1)
    print(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    trend=trend_load(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    print(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_seasonal_median",".nc",sep=""))
    nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_seasonal_median",".nc",sep=""))
    seasonal_median=var.get.nc(nc,"seasonal_median")
    detrended=dat$tas-trend
    for (year in 1:length(dat$year)){
        detrended[,,year]=detrended[,,year]-seasonal_median
    }
    per=state_attribution(detrended=detrended,nday=nday,nyr=nyr,filename=paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_state_ind.nc",sep=""))
}

master_duration <- function(){
    # calculate duration periods 2 states
    stateIndeces=c(-1,1)

    # calculate annual period durations
    nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_state_ind.nc",sep=""))
    ind=var.get.nc(nc,"ind")
    cat("\nidentifying persistent periods:")
    calc_global_dur(ind=ind,filename=paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_duration_4seasons.nc",sep=""),states=stateIndeces)

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

master_gridded_analysis <- function(ID_select=1:ntot){
    # analyzes period records for different time periods
    for (i in 1:(length(yearLimits)/2)){
        yearPeriod=c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        period<<-paste(yearLimits[(2*(i-1)+1)],"-",yearLimits[(2*(i-1)+2)],sep="")
        print(yearPeriod)

        # mean persistence, linear regression
        print("others")
        duration_analysis(yearPeriod=yearPeriod,option=c(1,0,0,0,0,0,0,0))

        # persistence quantiles, quantile regressions
        print("quant")
        duration_analysis(yearPeriod=yearPeriod,option=c(0,1,0,0,0,0,0,0),noise_level=0.00001)
        
        # fits
        print("fit")
        duration_analysis(yearPeriod=yearPeriod,option=c(0,0,0,1,0,0,0,0),plot_select=c(NA),ID_select=ID_select,add_name="2expo_4:100",xStart=1,write=TRUE)
    }

}

master_gridded_plots <- function(){
    # plots the results of gridded analysis
    for (i in 1:(length(yearLimits)/2)){
        period<-paste(yearLimits[(2*(i-1)+1)],"-",yearLimits[(2*(i-1)+2)],sep="")
        print(period)

        #others
        plot_maps(file="_others",var="other_stuff",period=period,sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c(""),name_zusatz="mean",signi_level=0.05,farb_mitte=c(2,9),farb_palette="regenbogen")

        #indexTopRight<<-c("")
        #indexBottomLeft<<-c("MAM\ncold","MAM\nwarm","JJA\ncold","JJA\nwarm","SON\ncold","SON\nwarm","DJF\ncold","DJF\nwarm","Annual\ncold","Annual\nwarm")
        plot_maps(file="_others",var="other_stuff",sub_auswahl=c(NA),value_auswahl=c(4),sig_auswahl=c(10),value_zusatz=c(""),name_zusatz="lr",period=period,signi_level=0.05,farb_mitte=c(-0.1,0.1),farb_palette="lila-gruen")

        # quantiles
        plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(2),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c(""),name_zusatz="qu_95",farb_mitte=c(8,26),farb_palette="regenbogen")
        plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(2),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c(""),name_zusatz="qr_sl_95",farb_mitte=c(-0.4,0.4),signi_level=0.05)
        plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(1),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c(""),name_zusatz="qr_sl_75",farb_mitte=c(-0.2,0.2),signi_level=0.05)
        
        #fits
        plot_maps(file="_fit_2expo_4:100",var="fit_stuff",sub_auswahl=c(NA),value_auswahl=c(10,12,13,20),sig_auswahl=c(17,17,17,17),value_zusatz=c("b1","b2","threshold","distr_size"),name_zusatz="fit_2expo_4:100",period=period,signi_level=0,farb_mitte=c(0.1,0.3,0.1,0.3,5,15,20,50),farb_palette="lila-gruen")
    }
}

master_regional_analysis <- function(region_name="7rect",ID_length=7,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA"),ID_select=1:ID_length,plot_select=1:ID_length){
    # peristent periods from grid points are grouped into regional persistence records
    # for shuffling a different data format is needed
    regional_attribution(region_name=region_name,toDo=c(TRUE,TRUE,TRUE))
    # for Mann-Kendall test yearly averages are required
    duration_yearly_values(folder=paste("/regional/",region_name,"/",sep=""),ID_name=region_name,ID_select=ID_select,ID_length=ID_length)
    ID_name=paste("_",region_name,sep="")
    for (i in 1:(length(yearLimits)/2)){
        yearPeriod=c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        period<<-paste(yearPeriod[1],"-",yearPeriod[2],sep="")
        print(yearPeriod)

        # mean persistence, linear regression
        print("others")
        duration_analysis(yearPeriod=yearPeriod,option=c(1,0,0,0,0,0,0,0),ID_name=ID_name,ID_select=ID_select,ID_names=region_names,ID_length=ID_length,folder=paste("/regional/",region_name,"/",sep=""))

        # persistence quantiles, quantile regressions
        print("quant")
        duration_analysis(yearPeriod=yearPeriod,option=c(0,0,1,0,0,0,0,0),noise_level=c(0,0.000001),ID_name=ID_name,ID_select=ID_select,ID_names=region_names,ID_length=ID_length,folder=paste("/regional/",region_name,"/",sep=""))
       
        # fits 
        print("fit")
        duration_analysis(yearPeriod=yearPeriod,option=c(0,0,0,1,0,0,0,0),add_name="2expo_4:100",xStart=4,ID_name=ID_name,ID_select=ID_select,plot_select=plot_select,ID_names=region_names,ID_length=ID_length,folder=paste("/regional/",region_name,"/",sep=""))

    } 
}

master_regional_plots <- function(region_name="7rect",ID_length=7,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA"),ID_select=1:ID_length,plot_select=1:ID_length,hlines=c(30)){

    season_names<<-c("MAM","JJA","SON","DJF","Annual") 
    ID_name=paste("_",region_name,sep="")
    for (i in 1:(length(yearLimits)/2)){
        yearPeriod=c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        period<<-paste(yearPeriod[1],"-",yearPeriod[2],sep="")
        print(yearPeriod)

        #others
        plot_reg_maps(region_name=region_name,file="_others",var="other_stuff",sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),name_zusatz="mean",signi_level=0.05,farb_mitte=c(2,9),farb_palette="regenbogen")
        plot_reg_maps(region_name=region_name,file="_others",var="other_stuff",sub_auswahl=c(NA),value_auswahl=c(4),sig_auswahl=c(NA),value_zusatz=c("linear regression"),name_zusatz="lm",signi_level=0.05,farb_mitte=c(-0.07,0.07),farb_palette="lila-gruen")  
        plot_reg_maps(region_name=region_name,file="_shuffQuant",var="original_slopes",sub_auswahl=c(5),value_auswahl=c(1),sig_auswahl=c(2),value_zusatz=c("linear regression"),name_zusatz="lmSig",signi_level=0.05,farb_mitte=c(-0.07,0.07),farb_palette="lila-gruen")  

        #quants
        plot_reg_maps(region_name=region_name,file="_quantiles",var="quantile_stuff",sub_auswahl=c(2),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("quantile"),name_zusatz="quantile",farb_mitte=c(8,26),farb_palette="regenbogen")
        plot_reg_maps(region_name=region_name,file="_quantiles",var="quantile_stuff",sub_auswahl=c(2),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c("qr slope"),name_zusatz="qr_slope",farb_mitte=c(-0.35,0.35),signi_level=0.05)
        plot_reg_maps(region_name=region_name,file="_shuffQuant",var="original_slopes",sub_auswahl=c(3),value_auswahl=c(1),sig_auswahl=c(2),value_zusatz=c("qr slope"),name_zusatz="qr_slopeSig",farb_mitte=c(-0.35,0.35),signi_level=0.05)

        plot_reg_boxplots(region_name=region_name,file="_quantiles",var="quantile_stuff",name_zusatz="quants",ID_select=ID_select,hlines=hlines)

        #fits
        plot_reg_maps(region_name=region_name,file="_fit_2expo_4:100",var="fit_stuff",sub_auswahl=c(NA),value_auswahl=c(10,12,13,20),sig_auswahl=c(17,17,17,17),value_zusatz=c("b1","b2","threshold","distr_size"),name_zusatz="fit_2expo_4:100",signi_level=0,farb_mitte=c(0.1,0.3,0.1,0.3,5,15,20,50),farb_palette="regenbogen")
   

        plot_reg_fit_table(region_name=region_name,file="_fit_2expo_4:100",var="fit_stuff",name_zusatz="slopeTab",value_auswahl=c(12,14),val_names=c("b1","b2"),colorRange=c(0.05,0.35),ID_select=ID_select,hlines=hlines)
        plot_reg_fit_table(region_name=region_name,file="_fit_2expo_4:100",var="fit_stuff",name_zusatz="threshTab",value_auswahl=c(15),val_names=c(""),colorRange=c(4,14),ID_select=ID_select,hlines=hlines)

        print("MannKendall")
        duration_MannKendall(yearPeriod=yearPeriod,folder=paste("/regional/",region_name,"/",sep=""),ID_name=region_name,ID_select=ID_select,ID_length=ID_length,hlines=hlines)


    }
}

master_special_plots <- function(){
    # just some adapted plot routines 
    # results climatology
    color_legend<<-"right" ; margins<<-c(0,0,0,5)
    indexTopRight<<-c("")
    indexBottomLeft<<-c("")
    
    indexBottomLeft<<-c("MAM","JJA","SON","DJF","Annual")
    plot_state_mean_maps(file="_others",var="other_stuff",period="1950-2014",sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),name_zusatz="mean_cowa_fixed_axis",signi_level=0.05,farb_mitte=c(3,9),farb_palette="regenbogen")
    #color_legend<<-"seperate" ; margins<<-c(0,0,0,0)

    # seas
    indexBottomLeft<<-c("MAM","JJA","SON","DJF","Annual")
    plot_state_mean_maps(file="_others",var="other_stuff",period="1950-2014",sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),name_zusatz="mean_cowa_multi",signi_level=0.05,farb_mitte=c(3,9),farb_palette="regenbogen",ID_select=ID_select)
    plot_state_diff_maps(file="_others",var="other_stuff",period="1950-2014",sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),name_zusatz="diff_cowa_multi",signi_level=0.05,farb_mitte=c(-0.2,0.2),farb_palette="lila-gruen")
    indexTopRight<<-c("e","f","g","h")
    plot_seasonal_anomaly_maps(file="_others",var="other_stuff",period="1950-2014",sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),name_zusatz="seas_anom_cowa_mean",signi_level=0.05,farb_mitte=c(-2,2),farb_palette="lila-gruen",ID_select=ID_select)
    indexTopRight<<-c("a","b","c","d")

    # quant
    plot_state_mean_maps(file="_quantiles",var="quantile_stuff",period="1950-2014",sub_auswahl=c(2),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c(""),name_zusatz="qu95_cowa_multi",signi_level=0.05,farb_mitte=c(8,26),farb_palette="regenbogen")
    plot_state_diff_maps(file="_quantiles",var="quantile_stuff",period="1950-2014",sub_auswahl=c(2),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c(""),name_zusatz="diff_qu95_cowa_multi",signi_level=0.05,farb_mitte=c(-3,3),farb_palette="lila-gruen")
}

master_correlation <- function(){
    # correlation between eke and persistence
    # also possible for other clim. indices as nao, enso...
    
    correl_init()
    # eke 
    eke_dur_correl(level=3,plot=FALSE,detrending=TRUE)
    eke_dur_correl(plot=TRUE,detrending=TRUE,ID_select=c(488))


    indexTopRight<<-c("","","","")
    indexBottomLeft<<-c("MAM\ncold","MAM\nwarm","JJA\ncold","JJA\nwarm","SON\ncold","SON\nwarm","DJF\ncold","DJF\nwarm","Annual\ncold","Annual\nwarm")
    dur_correlation_plot(toCor_short="eke",toCor_name="EKE",toCor_shortZu="_500mbar",val=5,val_zusatz="_mean",farb_mitte=c(-1,1))
    dur_correlation_plot(toCor_short="eke",toCor_name="EKE",toCor_shortZu="_500mbar",val=3,val_zusatz="_95",farb_mitte=c(-5,5))
    
    #plot_init_EU_Had()
    #dur_correlation_plot(toCor_short="eke",toCor_name="EKE",toCor_shortZu="_500mbar",val=1,val_zusatz="_mean_EU",farb_mitte=c(-3,3))
    #dur_correlation_plot(toCor_short="eke",toCor_name="EKE",toCor_shortZu="_500mbar",val=3,val_zusatz="_95_EU",farb_mitte=c(-5,5))
    
    # noa
    #index_dur_correl(toCor_name="NAO",toCor_short="nao")
    #dur_correlation_plot(toCor_short="nao",toCor_name="NAO",toCor_shortZu="",val=1,val_zusatz="_mean",farb_mitte=c(-3,3))
    #dur_correlation_plot(toCor_short="nao",toCor_name="NAO",toCor_shortZu="",val=1,val_zusatz="_95",farb_mitte=c(-5,5))
}

master_sensitivity <- function(){
    # sensitivity tests
    general_trenID_sensitivity(trendIDs=c("91_5","91_7","91_9"),add_name="7",legend=c("5-year","9-year"))
    general_trenID_sensitivity(trendIDs=c("61_7","91_7","121_7"),add_name="91",legend=c("61-day","121-day"),yAxis=FALSE)

    indexTopRight<<-c("","","","")
    sens_gridded(trendIDs=c("61_7","91_7","121_7"),file="_others",var="other_stuff",sub_auswahl=NA,value_auswahl=1,name_zusatz="mean",farb_mitte=c(-5,5))
    
    sens_gridded(trendIDs=c("91_5","91_7","91_9"),file="_others",var="other_stuff",sub_auswahl=NA,value_auswahl=1,name_zusatz="mean",farb_mitte=c(-5,5))

    sens_gridded(period="1979-2011",file="_others",var="other_stuff",sub_auswahl=NA,value_auswahl=4,name_zusatz="lr",farb_mitte=c(-100,100))
    sens_gridded(file="_quantiles",var="quantile_stuff",sub_auswahl=0.95,value_auswahl=1,name_zusatz="qu95",farb_mitte=c(-5,5))


    season_names<<-c("MAM","JJA","SON","DJF","Annual")
    sens_regional_fits()
    sens_regional_trends(period="1979-2011")
    sens_regional_trends(period="1979-2011",region_name="overReg",regPos=26:30,regNumb=5)
    season_names<<-c("MAM","JJA","SON","DJF","4seasons")

}

master_timeLag <- function(){
    # tests between max anomaly and midpoint of period
    timeLag_analysis()
    plot_maps(file="_timeLags",var="lag_ana",period="1950-2014",sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c(""),name_zusatz="timeLagMid",signi_level=0.05,farb_mitte=c(-3,3),farb_palette="lila-gruen")
    plot_maps(file="_timeLags",var="lag_ana",period="1950-2014",sub_auswahl=c(NA),value_auswahl=c(2),sig_auswahl=c(NA),value_zusatz=c(""),name_zusatz="timeLagBeg",signi_level=0.05,farb_mitte=c(10,15),farb_palette="regenbogen")

}

###################################################################
# init: loading sources, setting variables ....
###################################################################

master_init <- function(id){
    print(id)
    source("functions_support.r")           # for 2d running mean, state attribution, missing value analysis

    source("functions_duration.r")          # creation persistence records, linear regressions, quantile regressions, fits
    source("functions_analysis_tools.r")    # small functions needed for fitting, qunatile analysis, ks-goodness of fit

    source("functions_regional.r")          # creation of regional persistence records, creation of regional files
    source("functions_MannKendall.r")       # regional mann-kendall test
    source("functions_correlation.r")       # correlation analysis
    source("functions_sensitivity.r")       # sensitivity tests
    source("functions_timeLag.r")           # time lags between max anomaly and midpoint of period

    source("write.r")                       # netcdf write functions
    source("load.r")                        # load function

    source("map_plot.r")                    # functions to nicely display gridpoints, region borders, ...
    source("plot_master.r")                 # functions to create maps
    source("plot_tables.r")                 # functions to plot tables of trends or fit parameters
    source("inits_plot.r")                  # inits for map resolutions, areas, labels etc...


    library(moments)
    library(quantreg)
    library(stats4)
    library(RNetCDF)
    library(SDMTools)
    library(fields)

    nday<<-91
    nyr<<-id
    trendID<<-paste(nday,"_",nyr,sep="")
    dataset<<-"_TMean"
    trend_style<<-"_mean"
    additional_style<<-""
    dat<<-dat_load(paste("../data/",dataset,"/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    ntot<<-length(dat$ID)

    # borders of periods to be studied
    yearLimits<<-c(1979,2011,1979,1995,1995,2011)
    yearLimits<<-c(1950,2014)

    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    seasonal_boundaries<<-array(c(60,152,244,335,1,151,243,334,424,365),c(5,2))
    state_names<<-c("cold","warm")

    taus<<-c(0.75,0.95,0.99)
}

###################################################################
#parallel for different trends
id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(id)
if (!is.na(id)){id<-id*2+3}
if (is.na(id)){id<-7}

#if (id!=7){dasdas}


###################################################################
# basic analysis
###################################################################
master_init(id)
#plot_init_Had_multiple()
plot_init_Had_multiple_noAA()
#plot_init_multi_SH()

# this could be used to consider only grid-points with less than x% missing
#naRa<-read.table(paste("../data/",dataset,"/naRatio.txt",sep=""))
#ID_select<-which(naRa[4,]==1)

# I used all grid-points
ID_select<-1:ntot

master_trend()
master_seasonal_median_on_detrended()
master_state_attribution()
master_duration()

###################################################################
# fits, quantiles etc on grid-point level
###################################################################

master_gridded_analysis()
master_gridded_plots()


###################################################################
# regional commands
###################################################################

# 24 cluster-region analysis
master_regional_analysis(region_name="ward24",ID_length=24,region_names=1:24)
master_regional_plots(region_name="ward24",ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),ID_length=24,region_names=1:24,hlines=c(23,20,22,8))

# climate zone analysis
master_regional_analysis(region_name="overReg",ID_length=5,region_names=1:5)
duration_MannKendall(yearPeriod=c(1979,2011),folder=paste("/regional/","overReg","/",sep=""),ID_name="overReg",ID_select=1:5,ID_length=5,hlines=c(30),colorbar=TRUE,header=FALSE,regLabel=c("NHpo","NHml","NHst","Tro","SHml"))

###################################################################
# special stuff
###################################################################

master_special_plots()

plot_all_changes_table(ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8))

master_nas()

master_correlation()

master_sensitivity()

master_timeLag()

zonaly_averaged_plot()