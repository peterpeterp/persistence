
master_init <- function(id){
    print(id)
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

    nday<<-91
    nyr<<-id
    trendID<<-paste(nday,"_",nyr,sep="")
    dataset<<-"_TMean"
    trend_style<<-"_mean"
    additional_style<<-""
    dat<<-dat_load(paste("../data/",dataset,"/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    ntot<<-length(dat$ID)
    #yearLimits<<-c(1950,1980)
    yearLimits<<-c(1980,2014)


    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    state_names<<-c("cold","warm")
}

plot_init_EU <- function(){
    paper<<-c(8,5)
    yAusschnitt<<-c(20,80)
    xAusschnitt<<-c(-30,80)
    asp<<-1
    pointsize<<-1
    pch_points<<-c(1,NA,1.875,1.25)

    pch_sig<<-4
    col_sig<<-rgb(0.1,0.1,0.1,0.6)
    cex_sig<<-0.1

    region<<-NA

    season_auswahl<<-1:5
    sub_zusatz<<-c("75th","95th","99th")
    name_reg_zusatz<<-"_EU"

    col_row<<-c(1,1)
    mat<<-NA
    layout_mat<<-c(NA)
    subIndex<<-c("a","b")
}


plot_init_USA <- function(){
    paper<<-c(8,3.5)
    yAusschnitt<<-c(20,50)
    xAusschnitt<<-c(-140,-50)
    asp<<-1
    pointsize<<-1
    pch_points<<-c(1,NA,1.875,1.25)

    pch_sig<<-4
    col_sig<<-rgb(0.1,0.1,0.1,0.6)
    cex_sig<<-0.1

    region<<-NA

    season_auswahl<<-1:5
    sub_zusatz<<-c("75th","95th","99th")
    name_reg_zusatz<<-"_USA"

    col_row<<-c(1,1)
    mat<<-NA
    layout_mat<<-c(NA)
    subIndex<<-c("a","b")
}

plot_init_reg7 <- function(){
    paper<<-c(8,5)
    yAusschnitt<<-c(35,60)
    xAusschnitt<<-c(-100,-50)
    asp<<-1
    pointsize<<-0.44
    pch_points<<-c(1,NA,1.875,1.25)

    pch_sig<<-4
    col_sig<<-rgb(0.1,0.1,0.1,0.6)
    cex_sig<<-0.03

    region<<-NA

    season_auswahl<<-c(1,2,3,4,5)
    sub_zusatz<<-c("75th","95th","99th")
    name_reg_zusatz<<-"_reg7"

    col_row<<-c(1,1)
    mat<<-NA
    layout_mat<<-c(NA)
    subIndex<<-c("a","b")
}

master_init(7)

#plot_init_EU()
#master_gridded_plots()

#plot_init_USA()
#master_gridded_plots()

plot_init_reg7()
master_gridded_plots(ID_select=which())