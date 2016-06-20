#!/home/pepflei/R/bin/Rscript

###################################################################
# Persistence in Precipitation time series
# Author: Peter Pfleiderer
# Institution: Potsdam Institute for Climate Impact Research (PIK)
# Year: 2015
###################################################################

master_nas <- function(){
    # count nas
    find_nas(dat)    
}

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
    print(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,dataset,"_state_ind.nc",sep=""))
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

master_gridded_analysis <- function(){
    for (i in 1:(length(yearLimits)/2)){
        yearPeriod=c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        print(yearPeriod)

        print("others")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(1,0,0,0,0,0,0,0))

        print("quant")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,1,0,0,0,0,0,0),noise_level=0.00001)
        
        print("fit")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,0,1,0,0,0,0),plot_select=c(NA),add_name="2expo_4:100",xStart=4,write=TRUE)
    }
}

master_gridded_plots <- function(){
    for (i in 1:(length(yearLimits)/2)){
        period<-paste(yearLimits[(2*(i-1)+1)],"-",yearLimits[(2*(i-1)+2)],sep="")
        print(period)

        #others
        plot_maps(file="_others",var="other_stuff",period=period,sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),name_zusatz="mean",signi_level=0.05,farb_mitte=c(2,9),farb_palette="regenbogen")
        plot_maps(file="_others",var="other_stuff",sub_auswahl=c(NA),value_auswahl=c(4),sig_auswahl=c(10),value_zusatz=c("linear regression"),name_zusatz="lm",period=period,signi_level=0.05,farb_mitte=c(-0.07,0.07),farb_palette="lila-gruen")
        # quantiles
        plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(2),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("quantile"),name_zusatz="qu_95",farb_mitte=c(8,28),farb_palette="regenbogen")
        plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(2),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c("qr slope"),name_zusatz="qr_sl_95",farb_mitte=c(-0.35,0.35),signi_level=0.05)
        plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(1),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c("qr slope"),name_zusatz="qr_sl_75",farb_mitte=c(-0.2,0.2),signi_level=0.05)
        
        #fits
        plot_maps(file="_fit_2expo_4:100",var="fit_stuff",sub_auswahl=c(NA),value_auswahl=c(6,8,9,14),sig_auswahl=c(17,17,17,17),value_zusatz=c("b1","b2","threshold","distr_size"),name_zusatz="fit_2expo_4:100",period=period,signi_level=0,farb_mitte=c(0.1,0.3,0.1,0.3,5,15,20,50),farb_palette="lila-gruen")
    }
}

master_support_analysis <- function(){
    source("precip-tools.r")
    yearLimits=c(1948,2006,1980,2006)
    for (i in 1:2){
        yearPeriod<<-c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        period<<-paste(yearLimits[(2*(i-1)+1)],"-",yearLimits[(2*(i-1)+2)],sep="")
        print(period)
        #state_ana()
        #state_ana_view(yAusschnitt=c(20,50),xAusschnitt=c(220,310),asp=1,paper=c(8,3.5),pointsize=0.58)
        memory_test(yAusschnitt=c(20,50),xAusschnitt=c(220,310),asp=1,paper=c(8,3.5),pointsize=0.58)
    }
}

master_regional_analysis <- function(region_name="all",ID_length=2,region_names="all_merged",ID_select=c(1),plot_select=c(1)){
    regional_attribution(region_name=region_name,IDregions=array(1,c(ntot,5)),regNumb=ID_length)
    ID_name=paste("_",region_name,sep="")
    for (i in 1:(length(yearLimits)/2)){
        yearPeriod=c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        period=paste(yearPeriod[1],"-",yearPeriod[2],sep="")
        print(yearPeriod)

        print("others")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(1,0,0,0,0,0,0,0),ID_name=ID_name,ID_select=ID_select,ID_names=region_names,ID_length=ID_length,folder=paste("/regional/",region_name,"/",sep=""))

        print("quant")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,1,0,0,0,0,0),noise_level=c(0,0.000001),ID_name=ID_name,ID_select=ID_select,ID_names=region_names,ID_length=ID_length,folder=paste("/regional/",region_name,"/",sep=""))
        
        print("fit")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,0,1,0,0,0,0),add_name="2expo_4:100",xStart=4,ID_name=ID_name,ID_select=ID_select,plot_select=plot_select,ID_names=region_names,ID_length=ID_length,folder=paste("/regional/",region_name,"/",sep=""))
    } 
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
    source("plot_master.r")
    source("map_plot.r")


    library(moments)
    library(quantreg)
    library(stats4)
    library(RNetCDF)
    library(SDMTools)
    library(fields)


    trendID<<-"0p5"
    dataset<<-"_ghcn7"
    additional_style<<-""
    filename<-paste("../data/_ghcn7/ghcn_pp_1950-2015_reg7_<10.nc",sep="") ; print(filename)
    dat <<- dat_load_precipitation(filename)
    ntot<<-length(dat$ID)
    yearLimits<<-c(1980,2015,1950,2015,1950,1980)


    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    state_names<<-c("dry","wet")
}

plot_init <- function(){
    paper<<-c(8,5)
    yAusschnitt<<-c(35,60)
    xAusschnitt<<-c(-100,-50)
    asp<<-1
    pointsize<<-0.44
    pch_points<<-c(1,NA,0.5,0.5)

    pch_sig<<-4
    col_sig<<-rgb(0.1,0.1,0.1,0.6)
    cex_sig<<-0.03

    region<<-NA

    season_auswahl<<-c(1,2,3,4,5)
    sub_zusatz<<-c("75th","95th","99th")
    name_reg_zusatz<<-""

    col_row<<-c(1,1)
    mat<<-NA
    layout_mat<<-c(NA)
    subIndex<<-c("a","b")
}

###################################################################
# basic analysis
###################################################################
master_init()
plot_init()

#master_state_attribution()
#master_duration()

###################################################################
# fits, quantiles etc
###################################################################

#master_gridded_analysis()
master_gridded_plots()

#master_regional_analysis()

#master_support_analysis()

