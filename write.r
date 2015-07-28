library(ncdf)

dat_write <- function(filename,data3D)
{
    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    year <- dim.def.ncdf("year",units="year",vals=1950:2014, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:length(data3D$ID), unlim=FALSE)

    
    reihen=list(data3D$lon,data3D$lat,data3D$tas)
    names=c("lon","lat","tas")

    varlon <- var.def.ncdf(name=names[1],units="bla",dim=ID, missval=-9999.0)
    varlat <- var.def.ncdf(name=names[2],units="bla",dim=ID, missval=-9999.0)
    vartas <- var.def.ncdf(name=names[3],units="bla",dim=list(ID,day,year), missval=-9999.0)
    vars=list(varlon,varlat,vartas)
   
    nc = create.ncdf(filename,vars)

    for (i in 1:3){
        put.var.ncdf(nc,vars[[i]],reihen[[i]])
    }

    close.ncdf(nc)    
}


trend_write <- function(filename,data3D,trend)
{
    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    year <- dim.def.ncdf("year",units="year",vals=1:65, unlim=FALSE)
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
    year <- dim.def.ncdf("year",units="year",vals=1:65, unlim=FALSE)
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



duration_write <- function(filename,dur,dur_mid,len)
{
    states=dim(dur)[2]
    ntot=1319
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    varstates <- dim.def.ncdf("states",units="states",vals=1:(states),unlim=FALSE)
    periods <- dim.def.ncdf("periods",units="periods",vals=1:len, unlim=FALSE)


    vardur <- var.def.ncdf(name="dur",units="days",longname=paste("duration of periods of same state, states beeing:",states),dim=list(ID,varstates,periods), missval=-9999.0)
    vardur_mid <- var.def.ncdf(name="dur_mid",units="days",longname=paste("midpoints of periods of same state, states beeing:",states),dim=list(ID,varstates,periods), missval=-9999.0)


    vars=list(vardur,vardur_mid)
   
    nc = create.ncdf(filename,vars)
    put.var.ncdf(nc,vardur,dur[1:ntot,1:states,])              
    put.var.ncdf(nc,vardur_mid,dur_mid[1:ntot,1:states,])              

    close.ncdf(nc) 
}

duration_analysis_write <- function(filename,dur,season,trenn){
    print(filename)
    states=dim(dur)[2]
    ntot=1319
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    varstates <- dim.def.ncdf("states",units="states",vals=1:(states),unlim=FALSE)
    quantiles <- dim.def.ncdf("0.25 0.5 0.75 0.9 0.95 0.98 NA lr",units="0-1",vals=1:8,unlim=FALSE)
    outs <- dim.def.ncdf("quantile_slope quantile_slope_sig quantile_mean quantile_mean_sd quantile_intercept",units="0-1",vals=1:5,unlim=FALSE)

    dur_ana_full <- var.def.ncdf(name="dur_ana_full",units="values",longname=paste("analysis of duration of periods in",season,"from 1950 to 2014"),dim=list(ID,varstates,quantiles,outs), missval=-9999.0)
    

    vars=list(dur_ana_full)
   
    nc = create.ncdf(filename,vars)
    put.var.ncdf(nc,vars[[1]],dur[1:ntot,1:states,1:8,1:5])      

    close.ncdf(nc) 
}

regional_analysis_write <- function(filename,y,y_sig,poli)
{
    region <- dim.def.ncdf("region",units="region",vals=1:dim(poli)[1], unlim=FALSE)
    varstates <- dim.def.ncdf("states",units="states",vals=1:2,unlim=FALSE)
    outs <- dim.def.ncdf("analysis",units="slope MK 0.25 0.5 0.75 0.95 0.9 0.99",vals=1:8,unlim=FALSE)

    poli_points <- dim.def.ncdf("poli_points",units="id",vals=1:12,unlim=FALSE)

    region_coordinates <- var.def.ncdf(name="region_coordinates",units="deg",longname="1:6 lon - 7:12 lat",dim=list(region,poli_points),missval=-9999.0)

    values <- var.def.ncdf(name="values",units="values",longname=paste("analysis of: LR MK 0.9 0.95 0.98 0.99"),dim=list(varstates,outs,region), missval=-9999.0)
    values_sig <- var.def.ncdf(name="values_sig",units="values_sig",longname=paste("significance of: LR MK 0.9 0.95 0.98 0.99"),dim=list(varstates,outs,region), missval=-9999.0)
    
    vars=list(values,values_sig,region_coordinates)
   
    nc = create.ncdf(filename,vars)
    put.var.ncdf(nc,vars[[1]],y[1:2,1:8,1:dim(poli)[1]])      
    put.var.ncdf(nc,vars[[2]],y_sig[1:2,1:8,1:dim(poli)[1]]) 

    pol_poi=array(NA,c(dim(poli)[1],12))
    for (i in 1:dim(poli)[1]){
        for (j in 1:12){
            
            if (is.numeric(poli[i,j])){
                pol_poi[i,j]=poli[i,j]
            }
        }
    }
    put.var.ncdf(nc,region_coordinates,pol_poi)      

    close.ncdf(nc) 
}