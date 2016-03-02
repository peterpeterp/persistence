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
    dim.def.nc(nc_out,"year",dimlength=66, unlim=FALSE)   

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
    yearLimits=c(1950,2015,1980,2015)
    for (i in 1:2){
        yearPeriod=c(yearLimits[(2*(i-1)+1)],yearLimits[(2*(i-1)+2)])
        print(yearPeriod)

        print("others")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(1,0,0,0,0,0,0,0))

        print("quant")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,1,0,0,0,0,0,0))
        
        print("fit")
        duration_analysis(yearPeriod=yearPeriod,trendID=trendID,dataset=dataset,option=c(0,0,0,1,0,0,0,0),plot_select=c(NA),ID_select=ID_select,add_name="2expo_4:100",xStart=4,write=TRUE)
    }

}

master_gridded_plots <- function(ID_select=1:1319){
    yearLimits=c(1950,2015,1980,2015)
    for (i in 1:2){
        period<-paste(yearLimits[(2*(i-1)+1)],"-",yearLimits[(2*(i-1)+2)],sep="")
        print(period)

        #others
        plot_maps(file="_others",var="other_stuff",period=period,sub_auswahl=c(NA),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("mean period length"),sub_zusatz=c(NA),name_zusatz="mean",signi_level=0.05,farb_mitte=c(2,9),farb_palette="regenbogen")
        plot_maps(file="_others",var="other_stuff",sub_auswahl=c(NA),value_auswahl=c(4),sig_auswahl=c(10),value_zusatz=c("linear regression"),sub_zusatz=c(NA),name_zusatz="lm",period=period,signi_level=0.05,farb_mitte=c(-0.07,0.07),farb_palette="lila-gruen")

        # quantiles
        plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(2),value_auswahl=c(1),sig_auswahl=c(NA),value_zusatz=c("quantile"),sub_zusatz=c("95th"),name_zusatz="quantile",farb_mitte=c(8,28),farb_palette="regenbogen")
        plot_maps(file="_quantiles",var="quantile_stuff",period=period,sub_auswahl=c(2),value_auswahl=c(2),sig_auswahl=c(3),value_zusatz=c("qr slope"),sub_zusatz=c("95th","100th"),name_zusatz="qr_slope",farb_mitte=c(-0.35,0.35),signi_level=0.05)
        
        #fits
        plot_maps(file="_fit_2expo_4:100",var="fit_stuff",sub_auswahl=c(NA),value_auswahl=c(6,8,9,14),sig_auswahl=c(17,17,17,17),value_zusatz=c("b1","b2","threshold","distr_size"),sub_zusatz=c(NA),name_zusatz="fit_2expo_4:100",period=period,signi_level=0,farb_mitte=c(0.1,0.3,0.1,0.3,5,15,20,50),farb_palette="lila-gruen")

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
    dataset<<-"_eobs"
    additional_style<<-""
    dat <<- dat_load_precipitation(paste("../data/",dataset,"/rr_0.50deg_reg_v12.0_1950-2015.nc",sep=""))


    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    state_names<<-c("dry","wet")
}


###################################################################
# basic analysis
###################################################################
master_init()

#master_state_attribution()
#master_duration()

###################################################################
# fits, quantiles etc
###################################################################

master_gridded_analysis()
master_gridded_plots()

#-1--2319-2-2219-3-2219-4-2021-5-2220-6-2219-7-2121-8-2025-9-1824-10-2120-11-2222-12-2027-13-2023-14-2123-15-1928-16-1926-17-2024-18-2123-19-1627-20-1927-21-2025-22-2123-23-1923-24-1826-25-1928-26-2324-27-2023-28-1923-29-1826-30-2028-31-2326-32-1824-33-1924-34-1827-35-2028-36-2027-37-1824-38-2023-39-1727-40-1927-41-2027-42-1924-43-2022-44-1827-45-1926-46-2122-47-2020-48-2020-49-2122-50-2221-51-2121-52-2020-53-1923-54-2222-55-2418-56-2418-57-2221-58-2517-59-2418-60-2018-61-2117-62-2517-63-2517-64-2018-65-1820-66-1918-67-2018-68-1821-69-1918-70--71-1820-72-1821-73-1918-74--75--76--77-1821-78-1821-79--80--81--82-1918-83--84--85--86-7-87--88--89-437-90-427-91--92-407-93-418-94--95--96--97--98--99--100--101--102--103--104--105--106--107--108--109-2035-110-1936-111-1934-112-1931-113--114--115--116--117--118--119--120--121--122-3517-123-3318-124-2032-125-1934-126-2030-127-2030-128-2031-129-1932-130--131--132--133--134--135--136--137--138--139-4014-140-4216-141-3617-142-3317-143-3317-144-3420-145-3421-146-3423-147-3127-148-2827-149-2624-150-2325-151-2130-152-2126-153-2131-154-2133-155-1933-156--157--158--159--160--161--162--163--164--165-4314-166-4014-167-4115-168-4116-169-3616-170-3616-171-3420-172-3321-173-3325-174-3423-175-3323-176-3027-177-2725-178-2621-179-2128-180-2227-181-2328-182-2130-183-1931-184-1936-185--186--187--188--189--190--191--192--193--194-4313-195-4212-196-4114-197-4016-198-3916-199-3416-200-3618-201-3721-202-3423-203-3421-204-3124-205-3021-206-2822-207-2420-208-2621-209-2325-210-2229-211-2327-212-2125-213-2031-214-1733-215--216--217--218--219--220--221--222--223-4112-224-4213-225-4412-226-4113-227-4215-228-4015-229-4015-230-3715-231-3817-232-3519-233-3321-234-3522-235-3321-236-3021-237-2820-238-2520-239-2521-240-2428-241-2127-242-2127-243-2025-244-1729-245--246--247--248--249--250--251--252--253-3713-254-4013-255-4411-256-4411-257-4114-258-4113-259-3915-260-3915-261-3817
