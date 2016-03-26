#!/home/pepflei/R/bin/Rscript

###################################################################
# Persistence in Temperature time series
# Author: Peter Pfleiderer
# Institution: Potsdam Institute for Climate Impact Research (PIK)
# Year: 2015
###################################################################

master_nas <- function(){
    # count nas
    naRatio=array(NA,c(5,ntot))
    for (q in 1:ntot){naRatio[1,q]<-length(which(is.na(dat$tas[q,,])))/(365*length(dat$year))}
    naRatio[2,which(naRatio[1,]<0.3)]=naRatio[1,which(naRatio[1,]<0.3)]
    naRatio[3,which(naRatio[1,]<0.2)]=naRatio[1,which(naRatio[1,]<0.2)]
    naRatio[4,which(naRatio[1,]<0.1)]=naRatio[1,which(naRatio[1,]<0.1)]
    naRatio[5,which(naRatio[1,]<0.01)]=naRatio[1,which(naRatio[1,]<0.01)]

    topo_map_plot(filename=paste("../plots/",dataset,"/na_ratio.pdf",sep=""),reihen=naRatio,farb_mitte="mean",farb_palette="regenbogen",titel=c("na ratio","less than 30% nas","less than 20% nas","less than 10% nas","less than 1% nas"))
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

master_gridded_analysis <- function(ID_select=1:1319){
    for (i in 1:(length(yearLimits)/2)){
        yearPeriod=c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        print(yearPeriod)

        print("others")
        #duration_analysis(yearPeriod=yearPeriod,option=c(1,0,0,0,0,0,0,0))

        print("quant")
        #duration_analysis(yearPeriod=yearPeriod,option=c(0,1,0,0,0,0,0,0),noise_level=0.00001)
        
        print("fit")
        duration_analysis(yearPeriod=yearPeriod,option=c(0,0,0,1,0,0,0,0),plot_select=c(NA),ID_select=ID_select,add_name="2expo_4:100",xStart=1,write=TRUE)
    }

}

master_gridded_plots <- function(){
    for (i in 1:(length(yearLimits)/2)){
        period<-paste(yearLimits[(2*(i-1)+1)],"-",yearLimits[(2*(i-1)+2)],sep="")
        print(period)

        #others
        #plot_maps(file="_others",var="other_stuff",period=period,sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),name_zusatz="mean",signi_level=0.05,farb_mitte=c(2,9),farb_palette="regenbogen")
        #plot_maps(file="_others",var="other_stuff",sub_auswahl=c(NA),value_auswahl=c(4),sig_auswahl=c(10),value_zusatz=c("linear regression"),name_zusatz="lm",period=period,signi_level=0.05,farb_mitte=c(-0.07,0.07),farb_palette="lila-gruen")
        # quantiles
        #plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(2),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("quantile"),name_zusatz="qu_95",farb_mitte=c(8,28),farb_palette="regenbogen")
        #plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(2),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c("qr slope"),name_zusatz="qr_sl_95",farb_mitte=c(-0.35,0.35),signi_level=0.05)
        #plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(1),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c("qr slope"),name_zusatz="qr_sl_75",farb_mitte=c(-0.2,0.2),signi_level=0.05)
        
        #fits
        #plot_maps(file="_fit_2expo_4:100",var="fit_stuff",sub_auswahl=c(NA),value_auswahl=c(10,12,13,20),sig_auswahl=c(17,17,17,17),value_zusatz=c("b1","b2","threshold","distr_size"),name_zusatz="fit_2expo_4:100",period=period,signi_level=0,farb_mitte=c(0.1,0.3,0.1,0.3,5,15,20,50),farb_palette="lila-gruen")
    }
}

master_regional_analysis <- function(region_name="7rect",ID_length=7,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA"),ID_select=1:ID_length,plot_select=1:ID_length){
    regional_attribution(region_name=region_name)
    duration_yearly_values(folder=paste("/regional/",region_name,"/",sep=""),ID_name=region_name,ID_select=ID_select,ID_length=ID_length)
    ID_name=paste("_",region_name,sep="")
    for (i in 1:(length(yearLimits)/2)){
        yearPeriod=c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        period<<-paste(yearPeriod[1],"-",yearPeriod[2],sep="")
        print(yearPeriod)

        print("others")
        duration_analysis(yearPeriod=yearPeriod,option=c(1,0,0,0,0,0,0,0),ID_name=ID_name,ID_select=ID_select,ID_names=region_names,ID_length=ID_length,folder=paste("/regional/",region_name,"/",sep=""))

        print("quant")
        duration_analysis(yearPeriod=yearPeriod,option=c(0,0,1,0,0,0,0,0),noise_level=c(0,0.000001),ID_name=ID_name,ID_select=ID_select,ID_names=region_names,ID_length=ID_length,folder=paste("/regional/",region_name,"/",sep=""))
        
        print("fit")
        duration_analysis(yearPeriod=yearPeriod,option=c(0,0,0,1,0,0,0,0),add_name="2expo_4:100",xStart=4,ID_name=ID_name,ID_select=ID_select,plot_select=plot_select,ID_names=region_names,ID_length=ID_length,folder=paste("/regional/",region_name,"/",sep=""))

    } 
}

master_regional_plots <- function(region_name="7rect",ID_length=7,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA"),ID_select=1:ID_length,plot_select=1:ID_length){
    
    ID_name=paste("_",region_name,sep="")
    for (i in 1:(length(yearLimits)/2)){
        yearPeriod=c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        period<<-paste(yearPeriod[1],"-",yearPeriod[2],sep="")
        print(yearPeriod)

        #others
        plot_reg_maps(region_name=region_name,file="_others",var="other_stuff",sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),name_zusatz="mean",signi_level=0.05,farb_mitte=c(2,9),farb_palette="regenbogen")
        plot_reg_maps(region_name=region_name,file="_others",var="other_stuff",sub_auswahl=c(NA),value_auswahl=c(4),sig_auswahl=c(NA),value_zusatz=c("linear regression"),name_zusatz="lm",signi_level=0.05,farb_mitte=c(-0.07,0.07),farb_palette="lila-gruen")  
        #plot_reg_maps(region_name=region_name,file="_shuffQuant",var="original_slopes",sub_auswahl=c(5),value_auswahl=c(1),sig_auswahl=c(2),value_zusatz=c("linear regression"),name_zusatz="lmSig",signi_level=0.05,farb_mitte=c(-0.07,0.07),farb_palette="lila-gruen")  

        #quants
        plot_reg_maps(region_name=region_name,file="_quantiles",var="quantile_stuff",sub_auswahl=c(2),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("quantile"),name_zusatz="quantile",farb_mitte=c(8,28),farb_palette="regenbogen")
        plot_reg_maps(region_name=region_name,file="_quantiles",var="quantile_stuff",sub_auswahl=c(2),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c("qr slope"),name_zusatz="qr_slope",farb_mitte=c(-0.35,0.35),signi_level=0.05)
        #plot_reg_maps(region_name=region_name,file="_shuffQuant",var="original_slopes",sub_auswahl=c(3),value_auswahl=c(1),sig_auswahl=c(2),value_zusatz=c("qr slope"),name_zusatz="qr_slopeSig",farb_mitte=c(-0.35,0.35),signi_level=0.05)

        #fits
        plot_reg_maps(region_name=region_name,file="_fit_2expo_4:100",var="fit_stuff",sub_auswahl=c(NA),value_auswahl=c(10,12,13,20),sig_auswahl=c(17,17,17,17),value_zusatz=c("b1","b2","threshold","distr_size"),name_zusatz="fit_2expo_4:100",signi_level=0,farb_mitte=c(0.1,0.3,0.1,0.3,5,15,20,50),farb_palette="regenbogen")
   

        write_regional_fit_table(region_name=region_name,region_names=region_names,ID_select=ID_select,fit_style="2expo_4:100",ID_length=ID_length)

        print("MannKendall")
        duration_MannKendall(yearPeriod=yearPeriod,folder=paste("/regional/",region_name,"/",sep=""),ID_name=region_name,ID_select=ID_select,ID_length=ID_length)


    }
}

master_special_plots <- function(){


    # results climatology
    plot_init_Had()
    plot_state_mean_maps(file="_others",var="other_stuff",period="1950-2014",sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),name_zusatz="mean_cowa",signi_level=0.05,farb_mitte=c(2,9),farb_palette="regenbogen")

    plot_init_multi_midlat()
    plot_state_mean_maps(file="_others",var="other_stuff",period="1950-2014",sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),name_zusatz="mean_cowa_multi",signi_level=0.05,farb_mitte=c(2,9),farb_palette="regenbogen")

}

###################################################################
# init: loading sources, setting variables ....
###################################################################

master_init <- function(id){
    print(id)
    source("functions_support.r")
    source("functions_duration.r")
    source("functions_regional.r")
    source("functions_MannKendall.r")
    source("analysis_tools.r")
    source("write.r")
    source("load.r")
    source("plot_master.r")
    source("map_plot.r")
    source("inits_plot.r")


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
    #yearLimits<<-c(1950,1980)
    yearLimits<<-c(1980,2014,1950,2014,1950,1980)


    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    state_names<<-c("cold","warm")

    taus<<-c(0.75,0.95,0.99)
}

###################################################################
#parallel for different trends
id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(id)
if (!is.na(id)){id<-id*2+3}
if (is.na(id)){id<-7}


###################################################################
# basic analysis
###################################################################
master_init(id)
plot_init_Had()

#master_trend()
#master_seasonal_median_on_detrended()
#master_state_attribution()
#master_duration()

###################################################################
# fits, quantiles etc
###################################################################

#master_gridded_analysis()
#master_gridded_plots()


###################################################################
# regional commands
###################################################################

#master_regional_analysis(region_name="ward24",ID_length=24,region_names=1:24)
#master_regional_plots(region_name="ward24",ID_select=c(3,4,5,7,11,12,14,16,18,20,22),ID_length=24,region_names=1:24)

###################################################################
# special stuff
###################################################################

master_special_plots()