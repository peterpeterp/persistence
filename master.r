#!/home/pepflei/R/bin/Rscript
# Load useful functions 
#dyn.load("persistence_tools.so")
source("functions_support.r")
source("functions_duration.r")
source("write.r")
source("load.r")

master_nas <- function(){
    # count nas
    nday=91
    nyr=5
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2011.nc")
    find_nas(dat)    
}

master_trend <- function(){
    # calculate trend
    nday=91
    nyr=5
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend=calc_trend(dat,sprintf("../data/%s_%s/%s_%s_trend.nc",nday,nyr,nday,nyr),nday,nyr)
}

master_trend_control <- function(){
    # trend control
    source("trend_control.r")
    nday=91
    nyr=5
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    nc=open.ncdf("../data/91_5/91_5_markov2s.nc")
    ind=get.var.ncdf(nc,"ind")
    trend_control_warm_days(dat,ind)   
}

master_markov_2states <- function(){
    # calculate persistence 2 states
    nday=91
    nyr=5
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend.nc",nday,nyr,nday,nyr))
    per=calc_per(dat,trend,nday,nyr,model=markov_2states,states=2,transition_names="cc wc cw ww",
        filename=sprintf("../data/%s_%s/%s_%s_markov2s.nc",nday,nyr,nday,nyr))
}

master_analyse_markov_2states <- function(yearPeriod){
    # analyse persistence 2 states
    nc=open.ncdf("../data/91_5/2_states/markov/91_5_markov2s.nc")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    seasons=c("spring","summer","autumn","winter","year")
    for (sea in 1:length(seasons)){
        print(season)
        per=get.var.ncdf(nc,"markov")[1:length(dat$ID),sea,]
        tmp=global_analysis(per=per,filename_neu=paste("../data/91_5/2_states/markov/",yearPeriod[1],"-",yearPeriod[2],"/91_5_mar2s_trend_",seasons[sea],".nc",sep=""),
            season=seasons[sea],transition_names="cc wc cw ww",yearPeriod)
    }
}

master_duration_2states <- function(){
    # calculate duration periods 2 states
    nday=91
    nyr=5
    trash=((nyr-1)/2*365+(nday-1))

    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")

    nc=open.ncdf("../data/91_5/91_5_markov2s.nc")        
    ind=get.var.ncdf(nc,"ind")

    calc_global_dur(dat=dat,ind=ind,trash=trash,filename=sprintf("../data/%s_%s/%s_%s_duration_2s.nc",nday,nyr,nday,nyr),states=c(-1,1))

    nc=open.ncdf(sprintf("../data/%s_%s/%s_%s_duration_2s.nc",nday,nyr,nday,nyr))
    dur=get.var.ncdf(nc,"dur")
    dur_mid=get.var.ncdf(nc,"dur_mid")
    
    duration_seasons(dur,dur_mid,season=c(59,151),filename=sprintf("../data/%s_%s/%s_%s_duration_2s_spring.nc",nday,nyr,nday,nyr))
    duration_seasons(dur,dur_mid,season=c(151,243),filename=sprintf("../data/%s_%s/%s_%s_duration_2s_summer.nc",nday,nyr,nday,nyr))
    duration_seasons(dur,dur_mid,season=c(243,335),filename=sprintf("../data/%s_%s/%s_%s_duration_2s_autumn.nc",nday,nyr,nday,nyr))
    duration_seasons(dur,dur_mid,season=c(335,425),filename=sprintf("../data/%s_%s/%s_%s_duration_2s_winter.nc",nday,nyr,nday,nyr))
}

master_analyse_duration_2states <- function(yearPeriod){
    # analyse duration periods 2 states
    nday=91
    nyr=5
    for (season in c("spring","summer","autumn","winter","year")){
        print(season)
        nc=open.ncdf(paste("../data/",nday,"_",nyr,"/2_states/duration/",nday,"_",nyr,"_duration_2s_",season,".nc",sep=""))
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")
        duration_analysis(dur,dur_mid,filename=paste("../data/91_5/2_states/duration/",yearPeriod[1],"-",yearPeriod[2],"/91_5_duration_2s_analysis_",season,".nc",sep=""),
            season=season,yearPeriod)#,stations=seq(487,489,1))#,stations=134:136)
    }
}


master_markov_3states <- function(){
    # calculate persistence 3 states
    nday=91
    nyr=5
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend.nc",nday,nyr,nday,nyr))
    per=calc_per(dat,trend,nday,nyr,model=markov_3states,states=3,transition_names="cc nc wc cn nn wn cw nw ww",
        filename=sprintf("../data/%s_%s/%s_%s_markov3s.nc",nday,nyr,nday,nyr))
}

master_analyse_markov_3states <- function(yearPeriod){
    # analyse persistence 3 states
    nc=open.ncdf("../data/91_5/3_states/markov/91_5_markov3s.nc")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    seasons=c("spring","summer","autumn","winter","year")
    for (season in seasons){
        print(season)
        per=get.var.ncdf(nc,paste("markov_",season,sep=""))
        tmp=global_analysis(per=per,filename_neu=paste("../data/91_5/3_states/markov/",yearPeriod[1],"-",yearPeriod[2],"/91_5_mar3s_trend_",season,".nc",sep=""),
            season=season,transition_names="cc nc wc cn nn wn cw nw ww",yearPeriod)
    }
}

master_duration_3states <- function(){
    # calculate duration periods 3 states
    nday=91
    nyr=5
    trash=((nyr-1)/2*365+(nday-1))

    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")

    nc=open.ncdf("../data/91_5/91_5_markov3s.nc")        
    ind=get.var.ncdf(nc,"ind")

    calc_global_dur(dat=dat,ind=ind,trash=trash,filename=sprintf("../data/%s_%s/%s_%s_duration_3s.nc",nday,nyr,nday,nyr),states=c(-1,0,1))

    nc=open.ncdf(sprintf("../data/%s_%s/%s_%s_duration_3s.nc",nday,nyr,nday,nyr))
    dur=get.var.ncdf(nc,"dur")
    dur_mid=get.var.ncdf(nc,"dur_mid")
    
    duration_seasons(dur,dur_mid,season=c(59,151),filename=sprintf("../data/%s_%s/%s_%s_duration_3s_spring.nc",nday,nyr,nday,nyr))
    duration_seasons(dur,dur_mid,season=c(151,243),filename=sprintf("../data/%s_%s/%s_%s_duration_3s_summer.nc",nday,nyr,nday,nyr))
    duration_seasons(dur,dur_mid,season=c(243,335),filename=sprintf("../data/%s_%s/%s_%s_duration_3s_autumn.nc",nday,nyr,nday,nyr))
    duration_seasons(dur,dur_mid,season=c(335,425),filename=sprintf("../data/%s_%s/%s_%s_duration_3s_winter.nc",nday,nyr,nday,nyr))
}

master_analyse_duration_3states <- function(yearPeriod){
    # analyse duration periods 3 states
    nday=91
    nyr=5
    for (season in c("spring","summer","autumn","winter","year")){
        print(season)
        nc=open.ncdf(paste("../data/91_5/3_states/duration/91_5_duration_3s_",season,".nc",sep=""))
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")
        duration_analysis(dur,dur_mid,filename=paste("../data/91_5/3_states/duration/",yearPeriod[1],"-",yearPeriod[2],"/91_5_duration_3s_analysis_",season,".nc",sep=""),
            season=season,yearPeriod)
    }
}

master_regional_average <- function(yearPeriod,region_name){
    source("functions_regional.r")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    regional_analysis(dat=dat,yearPeriod,filepath=paste("../data/91_5/2_states/regional/",yearPeriod[1],"-",yearPeriod[2],"/91_5_",sep=""),region_name=region_name)

}

#master_markov_2states()

#master_analyse_markov_2states(yearPeriod=c(1950,1980))
#master_analyse_duration_2states(yearPeriod=c(1950,1980))
#master_analyse_duration_2states(yearPeriod=c(1980,2014))
#master_analyse_duration_2states(yearPeriod=c(1950,2014))
#master_regional_average(yearPeriod=c(1950,1980))
#master_regional_average(yearPeriod=c(1950,2014))
#master_regional_average(yearPeriod=c(1980,2014),region_name="7wave")


#master_analyse_duration_3states(yearPeriod=c(1950,2014))
#master_analyse_markov_3states(yearPeriod=c(1950,2014))
#master_analyse_markov_3states(yearPeriod=c(1980,2014))
#master_analyse_markov_3states(yearPeriod=c(1950,1980))

dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")

ntot=1319

if (1==2){
    nc=open.ncdf("../data/91_5/2_states/markov/91_5_markov2s.nc")
    markov_spring=get.var.ncdf(nc,"markov_spring")
    markov_summer=get.var.ncdf(nc,"markov_summer")
    markov_autumn=get.var.ncdf(nc,"markov_autumn")
    markov_winter=get.var.ncdf(nc,"markov_winter")
    markov_year=get.var.ncdf(nc,"markov_year")

    ind=get.var.ncdf(nc,"ind")

    markov=array(NA,dim=c(ntot,5,4,dim(markov_summer)[3]))
    markov[1:ntot,1,1:4,]=markov_spring[1:ntot,1:4,]
    markov[1:ntot,2,1:4,]=markov_summer[1:ntot,1:4,]
    markov[1:ntot,3,1:4,]=markov_autumn[1:ntot,1:4,]
    markov[1:ntot,4,1:4,]=markov_winter[1:ntot,1:4,]
    markov[1:ntot,5,1:4,]=markov_year[1:ntot,1:4,]

    markov_write("../data/91_5/2_states/markov/91_5_markov_2states.nc",dat,list(ind=ind,markov=markov),4,"cc wc cw ww")
}

if (1==1){
    nc=open.ncdf("../data/91_5/3_states/markov/91_5_markov3s.nc")
    markov_spring=get.var.ncdf(nc,"markov_spring")
    markov_summer=get.var.ncdf(nc,"markov_summer")
    markov_autumn=get.var.ncdf(nc,"markov_autumn")
    markov_winter=get.var.ncdf(nc,"markov_winter")
    markov_year=get.var.ncdf(nc,"markov_year")

    ind=get.var.ncdf(nc,"ind")

    markov=array(NA,dim=c(ntot,5,9,dim(markov_summer)[3]))
    markov[1:ntot,1,1:9,]=markov_spring[1:ntot,1:9,]
    markov[1:ntot,2,1:9,]=markov_summer[1:ntot,1:9,]
    markov[1:ntot,3,1:9,]=markov_autumn[1:ntot,1:9,]
    markov[1:ntot,4,1:9,]=markov_winter[1:ntot,1:9,]
    markov[1:ntot,5,1:9,]=markov_year[1:ntot,1:9,]

    markov_write("../data/91_5/3_states/markov/91_5_markov_3states.nc",dat,list(ind=ind,markov=markov),9,"cc nc wc cn nn wn cw nw ww")
}

for (period in c("1980-2014","1950-2014","1950-1980")){
    for (season in c("spring","summer","autumn","winter","year")){
        nc=open.ncdf(paste("../data/91_5/2_states/markov/",period,"/91_5_mar2s_trend_",season,".nc",sep=""),write=TRUE)
        att.put.ncdf(nc,0,"Description",paste("trend analysis for markov probabilities in",season, "for time interval:",period))
        att.put.ncdf(nc,"transitions","long_name","cc wc cw ww")
        att.put.ncdf(nc,"transitions","Longname","c -> cold w -> warm")
        att.put.ncdf( nc, "mean", "long_name", "mean value")
        att.put.ncdf( nc, "std", "long_name", "standard deviation")
        att.put.ncdf( nc, "MK", "long_name", "Mann Kendall test")
        att.put.ncdf( nc, "MK_sig", "long_name", "significance of Mann Kendall test")
        att.put.ncdf( nc, "LR", "long_name", "linear regression")
        att.put.ncdf( nc, "LR_sig", "long_name", "linear regression significance")
        close.ncdf(nc)
    }
}

for (period in c("1980-2014","1950-2014","1950-1980")){
    for (season in c("spring","summer","autumn","winter","year")){
        nc=open.ncdf(paste("../data/91_5/3_states/markov/",period,"/91_5_mar3s_trend_",season,".nc",sep=""),write=TRUE)
        att.put.ncdf(nc,0,"Description",paste("trend analysis for markov probabilities in",season, "for time interval:",period))
        att.put.ncdf(nc,"transitions","long_name","cc nc wc cn nn wn cw nw ww")
        att.put.ncdf(nc,"transitions","Longname","c -> cold n -> normal w -> warm")
        att.put.ncdf( nc, "mean", "long_name", "mean value")
        att.put.ncdf( nc, "std", "long_name", "standard deviation")
        att.put.ncdf( nc, "MK", "long_name", "Mann Kendall test")
        att.put.ncdf( nc, "MK_sig", "long_name", "significance of Mann Kendall test")
        att.put.ncdf( nc, "LR", "long_name", "linear regression")
        att.put.ncdf( nc, "LR_sig", "long_name", "linear regression significance")
        close.ncdf(nc)
    }
}