library(ncdf)

dat_write <- function(filename,data3D)
{
    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    year <- dim.def.ncdf("year",units="year",vals=1950:2011, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:length(data3D$ID), unlim=FALSE)

    
    reihen=list(data3D$lon,data3D$lat,data3D$region,data3D$tas)
    names=c("lon","lat","region","tas")

    varlon <- var.def.ncdf(name=names[1],units="bla",dim=ID, missval=-9999.0)
    varlat <- var.def.ncdf(name=names[2],units="bla",dim=ID, missval=-9999.0)
    varregion <- var.def.ncdf(name=names[3],units="bla",longname="SREX-Region",dim=ID, missval=-9999.0)
    vartas <- var.def.ncdf(name=names[4],units="bla",dim=list(ID,day,year), missval=-9999.0)
    vars=list(varlon,varlat,varregion,vartas)
   
    nc = create.ncdf(filename,vars)

    for (i in 1:4){
        put.var.ncdf(nc,vars[[i]],reihen[[i]])
    }

    close.ncdf(nc)    
}


trend_write <- function(filename,data3D,trend)
{
    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    year <- dim.def.ncdf("year",units="year",vals=1:62, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:length(data3D$ID), unlim=FALSE)

    vartrend <- var.def.ncdf(name="trend",units="bla",dim=list(ID,day,year), missval=-9999.0)
    nc = create.ncdf(filename,vartrend)
    put.var.ncdf(nc,vartrend,trend)

    close.ncdf(nc)  
}

markov_write <- function(filename,data3D,per,transitions,transition_names)
{
    ntot=length(data3D$ID)
    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    year <- dim.def.ncdf("year",units="year",vals=1:62, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    transition <- dim.def.ncdf("transition",units="transition",vals=1:transitions, unlim=FALSE)


    ind <- var.def.ncdf(name="ind",units="1 or -1",dim=list(ID,day,year), missval=-9999.0)
    markov_spring <- var.def.ncdf(name="markov_spring",units="0-1",longname=paste("spring markov",transition_names),dim=list(ID,transition,year), missval=-9999.0)
    markov_summer <- var.def.ncdf(name="markov_summer",units="0-1",longname=paste("summer markov",transition_names),dim=list(ID,transition,year), missval=-9999.0)
    markov_autumn <- var.def.ncdf(name="markov_autumn",units="0-1",longname=paste("autumn markov",transition_names),dim=list(ID,transition,year), missval=-9999.0)
    markov_winter <- var.def.ncdf(name="markov_winter",units="0-1",longname=paste("winter markov",transition_names),dim=list(ID,transition,year), missval=-9999.0)
    markov_year <- var.def.ncdf(name="markov_year",units="0-1",longname=paste("year markov",transition_names),dim=list(ID,transition,year), missval=-9999.0)

    markov_spring_conf <- var.def.ncdf(name="markov_spring_conf",units="0-1",longname="spring markov confidenceLevel",dim=list(ID,year), missval=-9999.0)
    markov_summer_conf <- var.def.ncdf(name="markov_summer_conf",units="0-1",longname="summer markov confidenceLevel",dim=list(ID,year), missval=-9999.0)
    markov_autumn_conf <- var.def.ncdf(name="markov_autumn_conf",units="0-1",longname="autumn markov confidenceLevel",dim=list(ID,year), missval=-9999.0)
    markov_winter_conf <- var.def.ncdf(name="markov_winter_conf",units="0-1",longname="winter markov confidenceLevel",dim=list(ID,year), missval=-9999.0)
    markov_year_conf <- var.def.ncdf(name="markov_year_conf",units="0-1",longname="year markov confidenceLevel",dim=list(ID,year), missval=-9999.0)

    vars=list(ind,
        markov_spring,markov_summer,markov_autumn,markov_winter,markov_year,
        markov_spring_conf,markov_summer_conf,markov_autumn_conf,markov_winter_conf,markov_year_conf)
   
    nc = create.ncdf(filename,vars)

    put.var.ncdf(nc,vars[[1]],per$ind)
    for (i in 1:5){
        put.var.ncdf(nc,vars[[i+1]],per$markov[1:ntot,i,,])
    }
    for (i in 1:5){
        put.var.ncdf(nc,vars[[i+6]],per$markov_conf[1:ntot,i,])
    }

    close.ncdf(nc) 
}


markov_analysis_write <- function(filename,analysis,season,transition_names)
{
    ntot=1319
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    transitions <- dim.def.ncdf("transitions",units="transisitons",vals=1:dim(analysis)[2],unlim=FALSE)

    mean <- var.def.ncdf(name="mean",units="bla",longname=paste("mean",season,transition_names),dim=list(ID,transitions), missval=-9999.0)
    std <- var.def.ncdf(name="std",units="bla",longname=paste("standard deviation",season,transition_names),dim=list(ID,transitions), missval=-9999.0)

    MK <- var.def.ncdf(name="MK",units="bla",longname=paste("MK",season,transition_names),dim=list(ID,transitions), missval=-9999.0)
    MK_sig <- var.def.ncdf(name="MK_sig",units="bla",longname=paste("MK_sig",season,transition_names),dim=list(ID,transitions), missval=-9999.0)
    LR <- var.def.ncdf(name="LR",units="bla",longname=paste("LR",season,transition_names),dim=list(ID,transitions), missval=-9999.0)
    LR_sig <- var.def.ncdf(name="LR_sig",units="bla",longname=paste("LR_sig",season,transition_names),dim=list(ID,transitions), missval=-9999.0)

    vars=list(mean,std,MK,MK_sig,LR,LR_sig)
   
    nc = create.ncdf(filename,vars)

    for (i in 1:6){
        put.var.ncdf(nc,vars[[i]],analysis[1:ntot,1:dim(analysis)[2],i])  
    }

    close.ncdf(nc) 
}


duration_write <- function(filename,dur,len)
{
    ntot=1319
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    periods <- dim.def.ncdf("periods",units="periods",vals=1:len, unlim=FALSE)

    dur_warm <- var.def.ncdf(name="dur_warm",units="days",longname="duration of warm periods",dim=list(ID,periods), missval=-9999.0)
    dur_cold <- var.def.ncdf(name="dur_cold",units="days",longname="duration of cold periods",dim=list(ID,periods), missval=-9999.0)
    dur_warm_mid <- var.def.ncdf(name="dur_warm_mid",units="days",longname="midpoint of duration of warm periods",dim=list(ID,periods), missval=-9999.0)
    dur_cold_mid <- var.def.ncdf(name="dur_cold_mid",units="days",longname="midpoint of duration of cold periods",dim=list(ID,periods), missval=-9999.0)


    vars=list(dur_warm,dur_cold,dur_warm_mid,dur_cold_mid)
   
    nc = create.ncdf(filename,vars)
    print(dim(dur))
    for (i in 1:4){
        put.var.ncdf(nc,vars[[i]],dur[1:ntot,i,])              
    }

    close.ncdf(nc) 
}


duration_analysis_write <- function(filename,dur,season,trenn)
{
    ntot=1319
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    output <- dim.def.ncdf("out ID",units="out ID",vals=1:8,unlim=FALSE)

    dur_ana_warm_full <- var.def.ncdf(name="dur_ana_warm_full",units="values",longname=paste("analysis of duration of warm periods in",season,"from 1950 to 2011"),dim=list(ID,output), missval=-9999.0)
    dur_ana_cold_full <- var.def.ncdf(name="dur_ana_cold_full",units="values",longname=paste("analysis of duration of cold periods in",season,"from 1950 to 2011"),dim=list(ID,output), missval=-9999.0)
    dur_ana_warm_before <- var.def.ncdf(name="dur_ana_warm_before",units="values",longname=paste("analysis of duration of warm periods in",season,"from 1950 to",trenn),dim=list(ID,output), missval=-9999.0)
    dur_ana_cold_before <- var.def.ncdf(name="dur_ana_cold_before",units="values",longname=paste("analysis of duration of cold periods in",season,"from 1950 to",trenn),dim=list(ID,output), missval=-9999.0)
    dur_ana_warm_after <- var.def.ncdf(name="dur_ana_warm_after",units="values",longname=paste("analysis of duration of warm periods in",season,"from",trenn,"to 2011"),dim=list(ID,output), missval=-9999.0)
    dur_ana_cold_after <- var.def.ncdf(name="dur_ana_cold_after",units="values",longname=paste("analysis of duration of cold periods in",season,"from",trenn,"to 2011"),dim=list(ID,output), missval=-9999.0)
    

    vars=list(dur_ana_warm_full,dur_ana_cold_full,
        dur_ana_warm_before,dur_ana_cold_before,
        dur_ana_warm_after,dur_ana_cold_after)
   
    nc = create.ncdf(filename,vars)

    put.var.ncdf(nc,vars[[1]],dur[1:ntot,1,1,])      
    put.var.ncdf(nc,vars[[2]],dur[1:ntot,1,2,])      
    put.var.ncdf(nc,vars[[3]],dur[1:ntot,2,1,])      
    put.var.ncdf(nc,vars[[4]],dur[1:ntot,2,2,])      
    put.var.ncdf(nc,vars[[5]],dur[1:ntot,3,1,])      
    put.var.ncdf(nc,vars[[6]],dur[1:ntot,3,2,])      


    close.ncdf(nc) 
}

