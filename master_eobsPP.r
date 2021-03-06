#!/home/pepflei/R/bin/Rscript

###################################################################
# Persistence in Precipitation time series
# Author: Peter Pfleiderer
# Institution: Potsdam Institute for Climate Impact Research (PIK)
# Year: 2015
###################################################################


master_state_attribution <- function(threshold=0.5){
    ## User parameters 
    ntot <- length(dat$ID)

    state_ind<-dat$pp*NA

    percentage<-0
    cat(paste("\n0 -> -> -> -> -> 100\n"))
    for (q in 1:ntot){
        if (q/ntot*100 > percentage){
            cat("-")
            percentage<-percentage+5
        }
        if (length(which(is.na(dat$pp[q,,])))<(365*20)){
            # Calculate state vector
            per_ind <- dat$pp[q,,]*NA 

            per_ind[dat$pp[q,,] < threshold]=-1
            per_ind[dat$pp[q,,] >= threshold]=1

            state_ind[q,,]=per_ind
        }     
    }

    print(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,dataset,"_state_ind.nc",sep=""))
    nc_out<-create.nc(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,dataset,"_state_ind.nc",sep=""))

    print(dim(state_ind))

    att.put.nc(nc_out, "NC_GLOBAL", "method", "NC_CHAR", paste("precipitation below (above)",threshold,"mm are dry (wet)."))
    att.put.nc(nc_out, "NC_GLOBAL", "explanation", "NC_CHAR", "state indices -1 (dry), 1 (wet)")
    
    dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)   
    dim.def.nc(nc_out,"day",dimlength=365, unlim=FALSE)   
    dim.def.nc(nc_out,"year",dimlength=length(dat$year), unlim=FALSE)   

    var.def.nc(nc_out,"ind","NC_SHORT",c(0,1,2))
    att.put.nc(nc_out, "ind", "missing_value", "NC_SHORT", -99)
    att.put.nc(nc_out, "ind", "dim_explanation", "NC_CHAR", "ID-day-year")
    att.put.nc(nc_out, "ind", "explanation", "NC_CHAR", "state attribution for each day at each gridpoint")

    var.put.nc(nc_out,"ind",state_ind)             
    cat("done.\n")
}

master_duration <- function(){
    # calculate duration periods 2 states
    stateIndeces=c(-1,1)

    # calculate annual period durations
    nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,dataset,"_state_ind.nc",sep=""))
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
        period<<-paste(yearLimits[(2*(i-1)+1)],"-",yearLimits[(2*(i-1)+2)],sep="")
        print(yearPeriod)

        print("others")
        duration_analysis(yearPeriod=yearPeriod,option=c(1,0,0,0,0,0,0,0))

        print("quant")
        duration_analysis(yearPeriod=yearPeriod,option=c(0,1,0,0,0,0,0,0),noise_level=0.00001)
        
        print("fit")
        #duration_analysis(yearPeriod=yearPeriod,option=c(0,0,0,1,0,0,0,0),plot_select=c(NA),ID_select=ID_select,add_name="2expo_4:100",xStart=1,write=TRUE)
    }

}

master_gridded_plots <- function(){
    for (i in 1:(length(yearLimits)/2)){
        period<<-paste(yearLimits[(2*(i-1)+1)],"-",yearLimits[(2*(i-1)+2)],sep="")
        print(period)

        #others
        #plot_maps(file="_others",var="other_stuff",period=period,sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c(""),name_zusatz="mean",signi_level=0.1,farb_mitte=c(2,9),farb_palette="regenbogen")
        #plot_maps(file="_others",var="other_stuff",sub_auswahl=c(NA),value_auswahl=c(4),sig_auswahl=c(10),value_zusatz=c(""),name_zusatz="lr",period=period,signi_level=0.05,farb_mitte=c(-0.02,0.02),farb_palette="lila-gruen")
        # quantiles
        #plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(2),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c(""),name_zusatz="qu_95",farb_mitte=c(4,23),farb_palette="regenbogen")
        #plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(2),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c(""),name_zusatz="qr_sl_95",farb_mitte=c(-0.15,0.15),signi_level=0.1,farb_palette="for_dim")  #,farb_mitte=c(-0.08,0.08)
        plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(2),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c(""),name_zusatz="qr_sl_95_JJA_dry",farb_mitte=c(-0.15,0.15),signi_level=0.1,farb_palette="for_dim")  #,farb_mitte=c(-0.08,0.08)
        #plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(1),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c(""),name_zusatz="qr_sl_75",farb_mitte=c(-0.1,0.1),signi_level=0.1)
        
        #fits
       # plot_maps(file="_fit_2expo_4:100",var="fit_stuff",sub_auswahl=c(NA),value_auswahl=c(6,8,9,14),sig_auswahl=c(17,17,17,17),value_zusatz=c("b1","b2","threshold","distr_size"),name_zusatz="fit_2expo_4:100",period=period,signi_level=0,farb_mitte=c(0.1,0.3,0.1,0.3,5,15,20,50),farb_palette="lila-gruen")
    }
}

master_support_analysis <- function(){
    source("precip-tools.r")
    yearLimits=c(1950,2015,1980,2015)
    for (i in 1:2){
        yearPeriod<<-c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        period<<-paste(yearLimits[(2*(i-1)+1)],"-",yearLimits[(2*(i-1)+2)],sep="")
        #state_ana()
        state_ana_view()
        #memory_test()
    }
}

master_correlation <- function(){
    # not working due to missing eke file / wrong grid

    #plot_init_Had_multiple()

    #eke_dur_correl(plot=FALSE,detrending=TRUE,ID_select<-1:ntot)


    #eke_dur_correl(plot=TRUE,detrending=TRUE,ID_select<-c(488))

    #dur_correlation_plot(toCor_short="eke",toCor_name="EKE",toCor_shortZu="850mbar",val=1,val_zusatz="_mean",farb_mitte=c(-3,3),ID_select<-1:ntot)
    #dur_correlation_plot(toCor_short="eke",toCor_name="EKE",toCor_shortZu="850mbar",val=3,val_zusatz="_95",farb_mitte=c(-5,5),ID_select<-1:ntot)

    #plot_init_EU_Had()
    #dur_correlation_plot(toCor_short="eke",toCor_name="EKE",toCor_shortZu="850mbar",val=1,val_zusatz="_mean_EU",farb_mitte=c(-3,3),ID_select<-1:ntot)
    #dur_correlation_plot(toCor_short="eke",toCor_name="EKE",toCor_shortZu="850mbar",val=3,val_zusatz="_95_EU",farb_mitte=c(-5,5),ID_select<-1:ntot)
}

###################################################################
# init: loading sources, setting variables ....
###################################################################

master_init <- function(){
    source("functions_support.r")
    source("functions_duration.r")
    source("functions_regional.r")
    source("functions_MannKendall.r")
    source("functions_correlation.r")
    source("functions_analysis_tools.r")
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


    trendID<<-"0p5"
    dataset<<-"_eobsPP"
    additional_style<<-""
    #dat <<- dat_load_precipitation(paste("../data/",dataset,"/rr_0.50deg_reg_v12.0_1950-2015.nc",sep=""))
    ntot<<-length(dat$ID)
    yearLimits<<-c(1980,2015)
    #yearLimits<<-c(1950,2015)


    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    state_names<<-c("dry","wet")

    taus<<-c(0.75,0.95,0.99)
}

###################################################################
# this is just a copy of master_HadGHCND.r
# for comments see master_HadGHCND.r

###################################################################
# basic analysis
###################################################################
master_init()
plot_init_EU()

#master_state_attribution()
#master_duration()

###################################################################
# fits, quantiles etc
###################################################################

#master_gridded_analysis()
master_gridded_plots()

#master_support_analysis()

#master_nas()

#master_correlation()