library(ncdf)

dat_load <- function(filename){
    nc  = open.ncdf(filename)
    dat = list()     
    dat$day  = get.var.ncdf(nc,"day")
    dat$year = get.var.ncdf(nc,"year")
    dat$ID   = get.var.ncdf(nc,"ID")
    dat$lon   = get.var.ncdf(nc,"lon")
    dat$lat   = get.var.ncdf(nc,"lat")
    dat$tas   = get.var.ncdf(nc,"tas")

    # Add additional time dimension (days are decimal points dt=1/365)
    dat$time = c(1:(length(dat$day)*length(dat$year)))/365 - 0.5*1/365 + min(dat$year)

    return(dat)
}

trend_load <- function(filename){
    nc  = open.nc(filename)  
    trend  = get.var.ncdf(nc,"trend")

    return(trend)
}


per_load <- function(filename){
    nc=open.nc(filename)
    ind   = get.var.ncdf(nc,"ind")
    ntot=dim(ind)[1]

    per = list(ind=ind,markov=array(NA,dim=c(ntot,6,62)),markov_mk=array(NA,dim=c(ntot,6)),markov_lr=array(NA,dim=c(ntot,6)),markov_mk_sig=array(NA,dim=c(ntot,6)),markov_lr_sig=array(NA,dim=c(ntot,6)),
            shock=array(NA,dim=c(ntot,3,62)),shock_mk=array(NA,dim=c(ntot,3)),shock_lr=array(NA,dim=c(ntot,3)),shock_mk_sig=array(NA,dim=c(ntot,3)),shock_lr_sig=array(NA,dim=c(ntot,3)))
    str_markov=c("year_warm","year_cold","sum_warm","sum_cold","win_warm","win_cold")
    for (i in 1:6){
        per$markov[1:ntot,i,]=get.var.ncdf(nc,str_markov[i])
    }
    str_markov_mk=c("year_warm_mk","year_cold_mk","sum_warm_mk","sum_cold_mk","win_warm_mk","win_cold_mk")
    for (i in 1:6){
        per$markov_mk[1:ntot,i]=get.var.ncdf(nc,str_markov_mk[i])
    }
    str_markov_mk_sig=c("year_warm_mk_sig","year_cold_mk_sig","sum_warm_mk_sig","sum_cold_mk_sig","win_warm_mk_sig","win_cold_mk_sig")
    for (i in 1:6){
        per$markov_mk_sig[1:ntot,i]=get.var.ncdf(nc,str_markov_mk_sig[i])
    }
    str_markov_lr=c("year_warm_lr","year_cold_lr","sum_warm_lr","sum_cold_lr","win_warm_lr","win_cold_lr")
    for (i in 1:6){
        per$markov_lr[1:ntot,i]=get.var.ncdf(nc,str_markov_lr[i])
    }
    str_markov_lr_sig=c("year_warm_lr_sig","year_cold_lr_sig","sum_warm_lr_sig","sum_cold_lr_sig","win_warm_lr_sig","win_cold_lr_sig")
    for (i in 1:6){
        per$markov_lr_sig[1:ntot,i]=get.var.ncdf(nc,str_markov_lr_sig[i])
    }

    str_shock=c("shock_y","shock_s","shock_w")
    for (i in 1:3){
        per$shock[1:ntot,i,]=get.var.ncdf(nc,str_shock[i])
    }
    str_shock_mk=c("shock_y_mk","shock_s_mk","shock_w_mk")
    for (i in 1:3){
        per$shock_mk[1:ntot,i]=get.var.ncdf(nc,str_shock_mk[i])
    }
    str_shock_mk_sig=c("shock_y_mk_sig","shock_s_mk_sig","shock_w_mk_sig")
    for (i in 1:3){
        per$shock_mk_sig[1:ntot,i]=get.var.ncdf(nc,str_shock_mk_sig[i])
    }
    str_shock_lr=c("shock_y_lr","shock_s_lr","shock_w_lr")
    for (i in 1:3){
        per$shock_lr[1:ntot,i]=get.var.ncdf(nc,str_shock_lr[i])
    }
    str_shock_lr_sig=c("shock_y_lr_sig","shock_s_lr_sig","shock_w_lr_sig")
    for (i in 1:3){
        per$shock_lr_sig[1:ntot,i]=get.var.ncdf(nc,str_shock_lr_sig[i])
    }


    per$ind <- get.var.ncdf(nc,"ind")

    return(per)
}


dat_write <- function(filename,data3D)
{


    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    year <- dim.def.ncdf("year",units="year",vals=1:62, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:length(data3D$ID), unlim=FALSE)

    
    reihen=c(data3D$lon,data3D$lat,data3D$tas)
    names=c("lon","lat","tas")
    vars=c(NA,NA,NA)
    for (i in 1:2){
        varsi <- var.def.ncdf(name=names[i],units="bla",dim=ID, missval=-9999.0)
        print(varsi)


    }
    for (i in 3:3){
        vars[i] <- var.def.ncdf(name=names[i],units="bla",dim=list(ID,day,year), missval=-9999.0)
        print(vars[i])
    }    
    nc = create.ncdf(filename,vars)
    for (i in 1:3){
        print(i)
        put.var.ncdf(nc,vars[i],reihen[i])
    }


    close.ncdf(nc)    
}


trend_write <- function(filename,data3D,trend)
{
    nc = create.nc(filename)
    dim.def.ncdf(nc, "day",  dimlength=length(data3D$day),  unlim=FALSE)
    dim.def.ncdf(nc, "year", dimlength=length(data3D$year), unlim=FALSE)
    dim.def.ncdf(nc, "ID",   dimlength=length(data3D$ID), unlim=FALSE)

    var.def.ncdf(nc,"trend","NC_FLOAT", c("ID","day","year"))
    att.put.ncdf(nc, "trend", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "trend", trend)

    close.nc(nc)
}

per_write <- function(filename,data3D,per)
{
    nc = create.nc(filename)
    dim.def.ncdf(nc, "day",  dimlength=length(data3D$day),  unlim=FALSE)
    dim.def.ncdf(nc, "year", dimlength=length(data3D$year), unlim=FALSE)
    dim.def.ncdf(nc, "ID",   dimlength=length(data3D$ID), unlim=FALSE)

    ntot=length(data3D$ID)

    var.def.ncdf(nc, "ind", "NC_FLOAT", c("ID","day","year"))
    att.put.ncdf(nc, "ind", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "ind",  per$ind)

# -----------shock persistence

#------------- bayesian info criterion
    if (length(per$bic) != 0){
        var.def.ncdf(nc,"bic_s","NC_FLOAT", c("ID","year"))
        att.put.ncdf(nc, "bic_s", "missing_value", "NC_FLOAT", -9999.0)

        var.def.ncdf(nc,"bic_w","NC_FLOAT", c("ID","year"))
        att.put.ncdf(nc, "bic_w", "missing_value", "NC_FLOAT", -9999.0)

        var.def.ncdf(nc,"bic_y","NC_FLOAT", c("ID","year"))
        att.put.ncdf(nc, "bic_y", "missing_value", "NC_FLOAT", -9999.0)

        var.put.nc(nc, "bic_y", per$bic[1:ntot,1,]) 
        var.put.nc(nc, "bic_s", per$bic[1:ntot,2,])
        var.put.nc(nc, "bic_w", per$bic[1:ntot,3,])
    }

#------------- persistence

    var.def.ncdf(nc,"shock_s","NC_FLOAT", c("ID","year"))
    att.put.ncdf(nc, "shock_s", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"shock_w","NC_FLOAT", c("ID","year"))
    att.put.ncdf(nc, "shock_w", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"shock_y","NC_FLOAT", c("ID","year"))
    att.put.ncdf(nc, "shock_y", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "shock_y", per$shock[1:ntot,1,]) 
    var.put.nc(nc, "shock_s", per$shock[1:ntot,2,])
    var.put.nc(nc, "shock_w", per$shock[1:ntot,3,])

    var.def.ncdf(nc,"shock_s_mk","NC_FLOAT", c("ID"))
    att.put.ncdf(nc, "shock_s_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"shock_w_mk","NC_FLOAT", c("ID"))
    att.put.ncdf(nc, "shock_w_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"shock_y_mk","NC_FLOAT", c("ID"))
    att.put.ncdf(nc, "shock_y_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "shock_y_mk", per$shock_mk[1:ntot,1]) 
    var.put.nc(nc, "shock_s_mk", per$shock_mk[1:ntot,2])
    var.put.nc(nc, "shock_w_mk", per$shock_mk[1:ntot,3])

    var.def.ncdf(nc,"shock_s_mk_sig","NC_FLOAT", c("ID"))
    att.put.ncdf(nc, "shock_s_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"shock_w_mk_sig","NC_FLOAT", c("ID"))
    att.put.ncdf(nc, "shock_w_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"shock_y_mk_sig","NC_FLOAT", c("ID"))
    att.put.ncdf(nc, "shock_y_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "shock_y_mk_sig", per$shock_mk_sig[1:ntot,1]) 
    var.put.nc(nc, "shock_s_mk_sig", per$shock_mk_sig[1:ntot,2])
    var.put.nc(nc, "shock_w_mk_sig", per$shock_mk_sig[1:ntot,3])

    var.def.ncdf(nc,"shock_s_lr","NC_FLOAT", c("ID"))
    att.put.ncdf(nc, "shock_s_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"shock_w_lr","NC_FLOAT", c("ID"))
    att.put.ncdf(nc, "shock_w_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"shock_y_lr","NC_FLOAT", c("ID"))
    att.put.ncdf(nc, "shock_y_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "shock_y_lr", per$shock_lr[1:ntot,1]) 
    var.put.nc(nc, "shock_s_lr", per$shock_lr[1:ntot,2])
    var.put.nc(nc, "shock_w_lr", per$shock_lr[1:ntot,3])

    var.def.ncdf(nc,"shock_s_lr_sig","NC_FLOAT", c("ID"))
    att.put.ncdf(nc, "shock_s_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"shock_w_lr_sig","NC_FLOAT", c("ID"))
    att.put.ncdf(nc, "shock_w_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"shock_y_lr_sig","NC_FLOAT", c("ID"))
    att.put.ncdf(nc, "shock_y_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "shock_y_lr_sig", per$shock_lr_sig[1:ntot,1]) 
    var.put.nc(nc, "shock_s_lr_sig", per$shock_lr_sig[1:ntot,2])
    var.put.nc(nc, "shock_w_lr_sig", per$shock_lr_sig[1:ntot,3])

#---------------- markov persistence

    var.def.ncdf(nc,"year_warm","NC_FLOAT", c("ID","year"))
    att.put.ncdf(nc, "year_warm", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"year_cold","NC_FLOAT", c("ID","year"))
    att.put.ncdf(nc, "year_cold", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"sum_warm","NC_FLOAT", c("ID","year"))
    att.put.ncdf(nc, "sum_warm", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"sum_cold","NC_FLOAT", c("ID","year"))
    att.put.ncdf(nc, "sum_cold", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"win_warm","NC_FLOAT", c("ID","year"))
    att.put.ncdf(nc, "win_warm", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc,"win_cold","NC_FLOAT", c("ID","year"))
    att.put.ncdf(nc, "win_cold", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "year_warm", per$markov[1:ntot,1,])
    var.put.nc(nc, "year_cold", per$markov[1:ntot,2,])

    var.put.nc(nc, "sum_warm", per$markov[1:ntot,3,])
    var.put.nc(nc, "sum_cold", per$markov[1:ntot,4,])

    var.put.nc(nc, "win_warm", per$markov[1:ntot,5,])
    var.put.nc(nc, "win_cold", per$markov[1:ntot,6,])


    var.def.ncdf(nc, "sum_warm_mk", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "sum_warm_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "sum_cold_mk", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "sum_cold_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "win_warm_mk", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "win_warm_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "win_cold_mk", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "win_cold_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "year_warm_mk", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "year_warm_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "year_cold_mk", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "year_cold_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "year_warm_mk",  per$markov_mk[1:ntot,1])
    var.put.nc(nc, "year_cold_mk",  per$markov_mk[1:ntot,2])

    var.put.nc(nc, "sum_warm_mk",  per$markov_mk[1:ntot,3])
    var.put.nc(nc, "sum_cold_mk",  per$markov_mk[1:ntot,4])

    var.put.nc(nc, "win_warm_mk",  per$markov_mk[1:ntot,5])
    var.put.nc(nc, "win_cold_mk",  per$markov_mk[1:ntot,6])


    var.def.ncdf(nc, "sum_warm_mk_sig", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "sum_warm_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "sum_cold_mk_sig", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "sum_cold_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "win_warm_mk_sig", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "win_warm_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "win_cold_mk_sig", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "win_cold_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "year_warm_mk_sig", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "year_warm_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "year_cold_mk_sig", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "year_cold_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "year_warm_mk_sig",  per$markov_mk_sig[1:ntot,1])
    var.put.nc(nc, "year_cold_mk_sig",  per$markov_mk_sig[1:ntot,2])

    var.put.nc(nc, "sum_warm_mk_sig",  per$markov_mk_sig[1:ntot,3])
    var.put.nc(nc, "sum_cold_mk_sig",  per$markov_mk_sig[1:ntot,4])

    var.put.nc(nc, "win_warm_mk_sig",  per$markov_mk_sig[1:ntot,5])
    var.put.nc(nc, "win_cold_mk_sig",  per$markov_mk_sig[1:ntot,6])


    var.def.ncdf(nc, "sum_warm_lr", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "sum_warm_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "sum_cold_lr", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "sum_cold_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "win_warm_lr", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "win_warm_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "win_cold_lr", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "win_cold_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "year_warm_lr", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "year_warm_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "year_cold_lr", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "year_cold_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "year_warm_lr",  per$markov_lr[1:ntot,1])
    var.put.nc(nc, "year_cold_lr",  per$markov_lr[1:ntot,2])

    var.put.nc(nc, "sum_warm_lr",  per$markov_lr[1:ntot,3])
    var.put.nc(nc, "sum_cold_lr",  per$markov_lr[1:ntot,4])

    var.put.nc(nc, "win_warm_lr",  per$markov_lr[1:ntot,5])
    var.put.nc(nc, "win_cold_lr",  per$markov_lr[1:ntot,6])


    var.def.ncdf(nc, "sum_warm_lr_sig", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "sum_warm_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "sum_cold_lr_sig", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "sum_cold_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "win_warm_lr_sig", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "win_warm_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "win_cold_lr_sig", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "win_cold_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "year_warm_lr_sig", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "year_warm_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.ncdf(nc, "year_cold_lr_sig", "NC_FLOAT", "ID")
    att.put.ncdf(nc, "year_cold_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "year_warm_lr_sig",  per$markov_lr_sig[1:ntot,1])
    var.put.nc(nc, "year_cold_lr_sig",  per$markov_lr_sig[1:ntot,2])

    var.put.nc(nc, "sum_warm_lr_sig",  per$markov_lr_sig[1:ntot,3])
    var.put.nc(nc, "sum_cold_lr_sig",  per$markov_lr_sig[1:ntot,4])

    var.put.nc(nc, "win_warm_lr_sig",  per$markov_lr_sig[1:ntot,5])
    var.put.nc(nc, "win_cold_lr_sig",  per$markov_lr_sig[1:ntot,6])



    close.nc(nc)

}

dat_write_part <-function(dat,qq,filename){
    size=length(qq)
    dat_klein = list(ID=array(NA,size),day=array(NA,365),year=array(NA,62),lon=array(NA,size),lat=array(NA,size),tas= array(NA,dim=c(size,365,62)))
    dat_klein$day[] = dat$day[]
    dat_klein$year[] = dat$year[]
    for (i in 1:length(qq)){
        q=qq[i]
        dat_klein$ID[i] = q
        dat_klein$lat[i] = dat$lat[q]
        dat_klein$lon[i] = dat$lon[q]
        dat_klein$tas[i,,] = dat$tas[q,,]
    }
    dat_write(filename,dat_klein)
    return(dat_klein)
}

trend_write_part <-function(trend,dat_klein,qq,filename){
    size=length(qq)
    trend_klein=array(NA,dim=c(size,365,62))  
    for (i in 1:length(qq)){
        q=qq[i]
        trend_klein[i,,] = trend[q,,]
    }
    trend_write(filename,dat_klein,trend_klein)
    return(trend_klein)
}

per_write_part <-function(per,dat_klein,qq,filename){
    size=length(qq)
    per_klein = list(ind=array(NA,dim=c(size,365,62)),markov=array(NA,dim=c(size,6,62)),markov_mk=array(NA,dim=c(size,6)),markov_lr=array(NA,dim=c(size,6)),markov_mk_sig=array(NA,dim=c(size,6)),markov_lr_sig=array(NA,dim=c(size,6)),
        shock=array(NA,dim=c(size,3,62)),shock_mk=array(NA,dim=c(size,3)),shock_lr=array(NA,dim=c(size,3)),shock_mk_sig=array(NA,dim=c(size,3)),shock_lr_sig=array(NA,dim=c(size,3)))
    for (i in 1:length(qq)){
        q=qq[i]
        per_klein$ind[i,,] = per$ind[q,,]
        per_klein$markov[i,,] = per$markov[q,,]
        per_klein$markov_mk[i,] = per$markov_mk[q,]
        per_klein$markov_mk_sig[i,] = per$markov_mk_sig[q,]
        per_klein$markov_lr[i,] = per$markov_lr[q,]
        per_klein$markov_lr_sig[i,] = per$markov_lr_sig[q,]

        per_klein$shock[i,,] = per$shock[q,,]
        per_klein$shock_mk[i,] = per$shock_mk[q,]
        per_klein$shock_mk_sig[i,] = per$shock_mk_sig[q,]
        per_klein$shock_lr[i,] = per$shock_lr[q,]
        per_klein$shock_lr_sig[i,] = per$shock_lr_sig[q,]


    }
    per_write(filename,dat_klein,per_klein)
    return(per_klein)
}