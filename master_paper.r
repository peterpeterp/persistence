#!/home/pepflei/R/bin/Rscript

###################################################################
# Persistence in Temperature time series
# Author: Peter Pfleiderer
# Institution: Potsdam Institute for Climate Impact Research (PIK)
# Year: 2015
###################################################################

master_regional_analysis <- function(region_name="7rect",ID_length=7,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA"),ID_select=1:ID_length,plot_select=1:ID_length){
    # peristent periods from grid points are grouped into regional persistence records
    # for shuffling a different data format is needed
    # for Mann-Kendall test yearly averages are required
    #duration_yearly_values(folder=paste("/regional/",region_name,"/",sep=""),ID_name=region_name,ID_select=ID_select,ID_length=ID_length)
    ID_name=paste("_",region_name,sep="")
    for (i in 1:(length(yearLimits)/2)){
        yearPeriod=c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        period<<-paste(yearPeriod[1],"-",yearPeriod[2],sep="")
        print(yearPeriod)

        # fits 
        print("fit")
        duration_analysis(yearPeriod=yearPeriod,option=c(0,0,0,1,0,0,0,0),add_name="2expo_4:100_prob",xStart=4,ID_name=ID_name,ID_select=ID_select,plot_select=plot_select,ID_names=region_names,ID_length=ID_length,folder=paste("/regional/",region_name,"/",sep=""))

    } 
}

master_regional_plots_paper <- function(region_name="7rect",ID_length=7,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA"),ID_select=1:ID_length,plot_select=1:ID_length,hlines=c(30)){

    season_names<<-c("MAM","JJA","SON","DJF","Annual") 
    ID_name=paste("_",region_name,sep="")
    for (i in 1:(length(yearLimits)/2)){
        yearPeriod=c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        period<<-paste(yearPeriod[1],"-",yearPeriod[2],sep="")
        print(yearPeriod)

        #fits
        indexBottomLeft<<-c("cold\nMAM","warm\nMAM","cold\nJJA","warm\nJJA","cold\nSON","warm\nSON","cold\nDJF","warm\nDJF","cold\nAnnual","warm\nAnnual")
        #plot_reg_maps(region_name=region_name,file="_fit_2expo_4:100",var="fit_stuff",sub_auswahl=c(NA),value_auswahl=c(12),sig_auswahl=c(24),value_zusatz=c(""),name_zusatz="fit_2expo_4:100_b1",signi_level=0.99,farb_mitte=c(3.5,6),farb_palette="gruen-rot",reg_select=ID_select,operation=inverse)
        #plot_reg_maps(region_name=region_name,file="_fit_2expo_4:100",var="fit_stuff",sub_auswahl=c(NA),value_auswahl=c(14),sig_auswahl=c(24),value_zusatz=c(""),name_zusatz="fit_2expo_4:100_b2_prob",signi_level=0,farb_mitte="mean",farb_palette="gruen-rot",reg_select=ID_select,operation=prob)
        #plot_reg_maps(region_name=region_name,file="_fit_2expo_4:100",var="fit_stuff",sub_auswahl=c(NA),value_auswahl=c(14),value2_auswahl=c(12),sig_auswahl=c(24),value_zusatz=c(""),name_zusatz="fit_2expo_4:100_b2-b1_prob",signi_level=0,farb_mitte=c(-9,9),farb_palette="blau-lila",reg_select=ID_select,operation=prob)
        plot_reg_maps(region_name=region_name,file="_fit_2expo_4:100",var="fit_stuff",sub_auswahl=c(NA),value_auswahl=c(15),sig_auswahl=c(24),value_zusatz=c(""),name_zusatz="fit_2expo_4:100_tr",signi_level=0,farb_mitte=c(4,14),farb_palette="cyan-orange",reg_select=ID_select)
   

        #plot_reg_fit_table(region_name=region_name,file="_fit_2expo_4:100",var="fit_stuff",name_zusatz="slopeTab_nh",value_auswahl=c(12,14),val_names=c("b1","b2"),colorRange=c(0.05,0.35),ID_select=ID_select,regLabel=region_names,hlines=hlines)
        #plot_reg_fit_table(region_name=region_name,file="_fit_2expo_4:100",var="fit_stuff",name_zusatz="threshTab_nh",value_auswahl=c(15),val_names=c(""),colorRange=c(4,14),ID_select=ID_select,regLabel=region_names,hlines=hlines)
    }
}

master_special_plots <- function(){
    indexBottomLeft<<-c("MAM","JJA","SON","DJF","Annual")
    plot_state_mean_maps(file="_others",var="other_stuff",period="1950-2014",sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),name_zusatz="mean_cowa_multi_nh",signi_level=0.05,farb_mitte=c(3,9),farb_palette="regenbogen",ID_select=ID_select)
    plot_seasonal_anomaly_maps(file="_others",var="other_stuff",period="1950-2014",sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),name_zusatz="seas_anom_cowa_mean_nh",signi_level=0.05,farb_mitte=c(-2,2),farb_palette="lila-gruen",ID_select=ID_select)

}

master_correlation <- function(){
    correl_init()
    indexTopRight<<-c("","","","")
    indexBottomLeft<<-c("MAM\ncold","MAM\nwarm","JJA\ncold","JJA\nwarm","SON\ncold","SON\nwarm","DJF\ncold","DJF\nwarm","Annual\ncold","Annual\nwarm")
    dur_correlation_plot(toCor_short="eke",toCor_name="EKE",toCor_shortZu="_500mbar",val=1,val_zusatz="_mean_nh",farb_mitte=c(-1,1))
    dur_correlation_plot(toCor_short="eke",toCor_name="EKE",toCor_shortZu="_500mbar",val=3,val_zusatz="_95_nh",farb_mitte=c(-5,5))
}

master_plot_monthly_correlation <- function(toCor_name="spi3",toCor_shortZu=""){
    filename<-paste('../data/',dataset,'/',trendID,'/gridded/',trendID,'',dataset,'_duration_',toCor_name,'_cor.nc',sep="") ; print(filename)
    correlation=var.get.nc(open.nc(filename),"correlation")
    reihen=array(NA,dim=c(8,ntot))
    reihen_sig=array(NA,dim=c(8,ntot))
    for (sea in 1:4){
        for (state in 1:2){
            reihen[((sea-1)*2+state),]=correlation[1,state,sea,]
            reihen_sig[((sea-1)*2+state),]=correlation[2,state,sea,]
        }
    }

    filename_plot<-paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_correl_",toCor_name,"_",toCor_shortZu,'pearson',".pdf",sep="")
    topo_map_plot(filename=filename_plot,reihen=reihen,reihen_sig=reihen_sig,farb_mitte=c(-1,1),farb_palette="gold-blau",titel=c(""),ID_select=ID_select)

    correlation_qu=var.get.nc(open.nc(filename),"correlation_qu")
    reihen=array(NA,dim=c(8,ntot))
    reihen_sig=array(NA,dim=c(8,ntot))
    for (sea in 1:4){
        for (state in 1:2){
            reihen[((sea-1)*2+state),]=correlation_qu[2,5,state,sea,]
            reihen_sig[((sea-1)*2+state),]=correlation_qu[3,5,state,sea,]
        }
    }

    filename_plot<-paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_correl_",toCor_name,"_",'95qu',".pdf",sep="")
    topo_map_plot(filename=filename_plot,reihen=reihen,reihen_sig=reihen_sig,farb_mitte=c(-3,3),farb_palette="gold-blau",titel=c(""),ID_select=ID_select)

    reihen=array(NA,dim=c(8,ntot))
    reihen_sig=array(NA,dim=c(8,ntot))
    for (sea in 1:4){
        for (state in 1:2){
            reihen[((sea-1)*2+state),]=correlation_qu[2,4,state,sea,]
            reihen_sig[((sea-1)*2+state),]=correlation_qu[3,4,state,sea,]
        }
    }

    filename_plot<-paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_correl_",toCor_name,"_",'75qu',".pdf",sep="")
    topo_map_plot(filename=filename_plot,reihen=reihen,reihen_sig=reihen_sig,farb_mitte=c(-3,3),farb_palette="gold-blau",titel=c(""),ID_select=ID_select)

}

master_plot_monthly_correlation()

###################################################################
# init: loading sources, setting variables ....
###################################################################

master_init <- function(id){
    print(id)
    source("functions_support.r")
    source("functions_duration.r")
    source("functions_regional.r")
    source("functions_MannKendall.r")
    source("functions_correlation.r")
    source("functions_sensitivity.r")
    source("functions_timeLag.r")
    source("functions_analysis_tools.r")
    source("write.r")
    source("load.r")
    source("plot_master.r")
    source("plot_tables.r")
    source("map_plot.r")
    source("inits_plot.r")


    library(latex2exp)
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
    yearLimits<<-c(1980,2014,1950,2014,1950,1980)
    yearLimits<<-c(1979,2011,1979,1995,1995,2011)


    yearLimits<<-c(1979,1995,1995,2011)
    yearLimits<<-c(1950,2014)
    #yearLimits<<-c(1979,2011,1950,2014)

    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    seasonal_boundaries<<-array(c(60,152,244,335,1,151,243,334,424,365),c(5,2))
    state_names<<-c("cold","warm")

    taus<<-c(0.75,0.95,0.99)
}

###################################################################
# this is just a copy of master_HadGHCND.r
# here only the plots needed for the paper are left
# for comments see master_HadGHCND.r

###################################################################
# basic analysis
###################################################################
master_init(7)

plot_init_Had_multiple_NH()

naRa<-read.table(paste("../data/",dataset,"/naRatio.txt",sep=""))
ID_select<-which(naRa[4,]==1)




###################################################################
# special stuff
###################################################################
#master_regional_analysis(region_name="ward24",ID_length=24,region_names=1:24)

#master_regional_plots_paper(region_name="ward24",ID_select=c(2,1,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22),ID_length=19,region_names=1:19,hlines=c(23,20))

#master_special_plots()

#master_correlation()

master_plot_monthly_correlation()