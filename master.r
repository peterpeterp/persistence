#!/home/pepflei/R/bin/Rscript
# Load useful functions 
#dyn.load("persistence_tools.so")
source("functions_support.r")
source("functions_duration.r")
source("write.r")
source("load.r")
dyn.load("persistence_tools.so")

if (1==2){
    # calculate trend
    nday=91
    nyr=5
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend=calc_trend(dat,sprintf("../data/%s_%s/%s_%s_trend.nc",nday,nyr,nday,nyr),nday,nyr)
}

if (1==2){
    # trend control
    source("trend_control.r")
    nday=91
    nyr=5
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    nc=open.ncdf("../data/91_5/91_5_markov2s.nc")
    ind=get.var.ncdf(nc,"ind")
    trend_control_warm_days(dat,ind)   
}

if (1==2){
    # calculate persistence 2 states
    nday=91
    nyr=5
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend.nc",nday,nyr,nday,nyr))
    per=calc_per(dat,trend,nday,nyr,model=markov_2states,states=2,transition_names="cc cw wc ww",
        filename=sprintf("../data/%s_%s/%s_%s_markov2s.nc",nday,nyr,nday,nyr))
}

if (1==2){
    # analyse persistence 2 states
    nc=open.ncdf("../data/91_5/91_5_markov2s.nc")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    seasons=c("spring","summer","autumn","winter","year")
    for (season in seasons){
        per=get.var.ncdf(nc,paste("markov_",season,sep=""))
        tmp=global_analysis(per=per,filename_neu=paste("../data/91_5/91_5_mar2s_trend_",season,".nc",sep=""),season=season,transition_names="cc cw wc ww")
    }
}

if (1==2){
    # calculate persistence 3 states
    nday=91
    nyr=5
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend.nc",nday,nyr,nday,nyr))
    per=calc_per(dat,trend,nday,nyr,model=markov_3states,states=3,transition_names="cc cn cw nc nn nw wc wn ww",
        filename=sprintf("../data/%s_%s/%s_%s_markov3s.nc",nday,nyr,nday,nyr))
}

if (1==2){
    # analyse persistence 3 states
    nc=open.ncdf("../data/91_5/91_5_markov3s.nc")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    seasons=c("spring","summer","autumn","winter","year")
    for (season in seasons){
        per=get.var.ncdf(nc,paste("markov_",season,sep=""))
        tmp=global_analysis(per=per,filename_neu=paste("../data/91_5/91_5_mar3s_trend_",season,".nc",sep=""),season=season,transition_names="cc cn cw nc nn nw wc wn ww")
    }
}

if (1==1){
    # calculate duration periods
    nday=91
    nyr=5
    trash=((nyr-1)/2*365+(nday-1))
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    per=markov_load(sprintf("../data/%s_%s/%s_%s_markov2s.nc",nday,nyr,nday,nyr),4)

    calc_global_dur(dat=dat,per=per,trash=trash,filename=sprintf("../data/%s_%s/%s_%s_duration_test.nc",nday,nyr,nday,nyr),states=c(-1,1))
    sfsd
    dur=duration_load(filename=sprintf("../data/%s_%s/%s_%s_duration.nc",nday,nyr,nday,nyr))
    
    duration_seasons(dur,season=c(59,151),filename=sprintf("../data/%s_%s/%s_%s_duration_spring.nc",nday,nyr,nday,nyr))
    duration_seasons(dur,season=c(151,243),filename=sprintf("../data/%s_%s/%s_%s_duration_summer.nc",nday,nyr,nday,nyr))
    duration_seasons(dur,season=c(243,335),filename=sprintf("../data/%s_%s/%s_%s_duration_autumn.nc",nday,nyr,nday,nyr))
    duration_seasons(dur,season=c(335,425),filename=sprintf("../data/%s_%s/%s_%s_duration_winter.nc",nday,nyr,nday,nyr))
}

if (1==1){
    # analyse duration periods
    nday=91
    nyr=5
    for (season in c("spring","summer","autumn","winter")){
        print(season)
        dur=duration_load(filename=paste("../data/",nday,"_",nyr,"/",nday,"_",nyr,"_duration_",season,".nc",sep=""))
        duration_analysis(dur,filename=paste("../data/",nday,"_",nyr,"/",nday,"_",nyr,"_duration_analysis_",season,".nc",sep=""),
            season=season,trenn=1980)#,stations=134:136)
    }
}


