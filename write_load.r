library(RNetCDF)






dat_load <- function(filename){
    nc  = open.nc(filename)
    dat = list()     
    dat$day  = var.get.nc(nc,"day")
    dat$year = var.get.nc(nc,"year")
    dat$ID   = var.get.nc(nc,"ID")
    dat$lon   = var.get.nc(nc,"lon")
    dat$lat   = var.get.nc(nc,"lat")
    dat$tas   = var.get.nc(nc,"tas")

    # Add additional time dimension (days are decimal points dt=1/365)
    dat$time = c(1:(length(dat$day)*length(dat$year)))/365 - 0.5*1/365 + min(dat$year)

    return(dat)
}

trend_load <- function(filename){
    nc  = open.nc(filename)  
    trend  = var.get.nc(nc,"trend")

    return(trend)
}


per_load <- function(filename){
    nc=open.nc(filename)
    ind   = var.get.nc(nc,"ind")
    ntot=dim(ind)[1]

    per = list(ind=ind,markov=array(NA,dim=c(ntot,6,62)),markov_mk=array(NA,dim=c(ntot,6)),markov_lr=array(NA,dim=c(ntot,6)),markov_mk_sig=array(NA,dim=c(ntot,6)),markov_lr_sig=array(NA,dim=c(ntot,6)),
            shock=array(NA,dim=c(ntot,3,62)),shock_mk=array(NA,dim=c(ntot,3)),shock_lr=array(NA,dim=c(ntot,3)),shock_mk_sig=array(NA,dim=c(ntot,3)),shock_lr_sig=array(NA,dim=c(ntot,3)))
    str_markov=c("year_warm","year_cold","win_warm","win_cold","sum_warm","sum_cold")
    for (i in 1:6){
        per$markov[1:ntot,i,]=var.get.nc(nc,str_markov[i])
    }
    str_markov_mk=c("year_warm_mk","year_cold_mk","win_warm_mk","win_cold_mk","sum_warm_mk","sum_cold_mk")
    for (i in 1:6){
        per$markov_mk[1:ntot,i]=var.get.nc(nc,str_markov_mk[i])
    }
    str_markov_mk_sig=c("year_warm_mk_sig","year_cold_mk_sig","win_warm_mk_sig","win_cold_mk_sig","sum_warm_mk_sig","sum_cold_mk_sig")
    for (i in 1:6){
        per$markov_mk_sig[1:ntot,i]=var.get.nc(nc,str_markov_mk_sig[i])
    }
    str_markov_lr=c("year_warm_lr","year_cold_lr","win_warm_lr","win_cold_lr","sum_warm_lr","sum_cold_lr")
    for (i in 1:6){
        per$markov_lr[1:ntot,i]=var.get.nc(nc,str_markov_lr[i])
    }
    str_markov_lr_sig=c("year_warm_lr_sig","year_cold_lr_sig","win_warm_lr_sig","win_cold_lr_sig","sum_warm_lr_sig","sum_cold_lr_sig")
    for (i in 1:6){
        per$markov_lr_sig[1:ntot,i]=var.get.nc(nc,str_markov_lr_sig[i])
    }

    str_shock=c("shock_s","shock_w","shock_y")
    for (i in 1:3){
        per$shock[1:ntot,i,]=var.get.nc(nc,str_shock[i])
    }
    str_shock_mk=c("shock_s_mk","shock_w_mk","shock_y_mk")
    for (i in 1:3){
        per$shock_mk[1:ntot,i]=var.get.nc(nc,str_shock_mk[i])
    }
    str_shock_mk_sig=c("shock_s_mk_sig","shock_w_mk_sig","shock_y_mk_sig")
    for (i in 1:3){
        per$shock_mk_sig[1:ntot,i]=var.get.nc(nc,str_shock_mk_sig[i])
    }
    str_shock_lr=c("shock_s_lr","shock_w_lr","shock_y_lr")
    for (i in 1:3){
        per$shock_lr[1:ntot,i]=var.get.nc(nc,str_shock_lr[i])
    }
    str_shock_lr_sig=c("shock_s_lr_sig","shock_w_lr_sig","shock_y_lr_sig")
    for (i in 1:3){
        per$shock_lr_sig[1:ntot,i]=var.get.nc(nc,str_shock_lr_sig[i])
    }


    per$ind <- var.get.nc(nc,"ind")

    return(per)
}


dat_write <- function(filename,data3D)
{
    nc = create.nc(filename)
    dim.def.nc(nc, "day",  dimlength=length(data3D$day),  unlim=FALSE)
    dim.def.nc(nc, "year", dimlength=length(data3D$year), unlim=FALSE)
    dim.def.nc(nc, "ID",   dimlength=length(data3D$ID), unlim=FALSE)

    var.def.nc(nc, "day", "NC_INT", "day")
    var.def.nc(nc, "year", "NC_INT", "year")
    var.def.nc(nc, "ID", "NC_INT", "ID")

    var.def.nc(nc, "lon", "NC_FLOAT", "ID")
    var.def.nc(nc, "lat", "NC_FLOAT", "ID")
 
    var.def.nc(nc, "tas", "NC_FLOAT", c("ID","day","year"))
    att.put.nc(nc, "tas", "long_name", "NC_CHAR", "Near-surface air temperature anomaly")
    att.put.nc(nc, "tas", "units", "NC_CHAR", "degrees Celcius")
    att.put.nc(nc, "tas", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "day",  data3D$day)
    var.put.nc(nc, "year", data3D$year)
    var.put.nc(nc, "ID",   data3D$ID)

    var.put.nc(nc, "lon",  data3D$lon)
    var.put.nc(nc, "lat",  data3D$lat)

    var.put.nc(nc, "tas",  data3D$tas)

    close.nc(nc)
}

trend_write <- function(filename,data3D,trend)
{
    nc = create.nc(filename)
    dim.def.nc(nc, "day",  dimlength=length(data3D$day),  unlim=FALSE)
    dim.def.nc(nc, "year", dimlength=length(data3D$year), unlim=FALSE)
    dim.def.nc(nc, "ID",   dimlength=length(data3D$ID), unlim=FALSE)

    var.def.nc(nc,"trend","NC_FLOAT", c("ID","day","year"))
    att.put.nc(nc, "trend", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "trend", trend)

    close.nc(nc)
}

per_write <- function(filename,data3D,per)
{
    nc = create.nc(filename)
    dim.def.nc(nc, "day",  dimlength=length(data3D$day),  unlim=FALSE)
    dim.def.nc(nc, "year", dimlength=length(data3D$year), unlim=FALSE)
    dim.def.nc(nc, "ID",   dimlength=length(data3D$ID), unlim=FALSE)

    ntot=length(data3D$ID)

    var.def.nc(nc, "ind", "NC_FLOAT", c("ID","day","year"))
    att.put.nc(nc, "ind", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "ind",  per$ind)

# -----------shock persistence

#------------- bayesian info criterion
    if (length(per$bic) != 0){
        var.def.nc(nc,"bic_s","NC_FLOAT", c("ID","year"))
        att.put.nc(nc, "bic_s", "missing_value", "NC_FLOAT", -9999.0)

        var.def.nc(nc,"bic_w","NC_FLOAT", c("ID","year"))
        att.put.nc(nc, "bic_w", "missing_value", "NC_FLOAT", -9999.0)

        var.def.nc(nc,"bic_y","NC_FLOAT", c("ID","year"))
        att.put.nc(nc, "bic_y", "missing_value", "NC_FLOAT", -9999.0)

        var.put.nc(nc, "bic_s", per$bic[1:ntot,1,]) 
        var.put.nc(nc, "bic_w", per$bic[1:ntot,2,])
        var.put.nc(nc, "bic_y", per$bic[1:ntot,3,])
    }

#------------- persistence

    var.def.nc(nc,"shock_s","NC_FLOAT", c("ID","year"))
    att.put.nc(nc, "shock_s", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"shock_w","NC_FLOAT", c("ID","year"))
    att.put.nc(nc, "shock_w", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"shock_y","NC_FLOAT", c("ID","year"))
    att.put.nc(nc, "shock_y", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "shock_s", per$shock[1:ntot,1,]) 
    var.put.nc(nc, "shock_w", per$shock[1:ntot,2,])
    var.put.nc(nc, "shock_y", per$shock[1:ntot,3,])

    var.def.nc(nc,"shock_s_mk","NC_FLOAT", c("ID"))
    att.put.nc(nc, "shock_s_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"shock_w_mk","NC_FLOAT", c("ID"))
    att.put.nc(nc, "shock_w_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"shock_y_mk","NC_FLOAT", c("ID"))
    att.put.nc(nc, "shock_y_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "shock_s_mk", per$shock_mk[1:ntot,1]) 
    var.put.nc(nc, "shock_w_mk", per$shock_mk[1:ntot,2])
    var.put.nc(nc, "shock_y_mk", per$shock_mk[1:ntot,3])

    var.def.nc(nc,"shock_s_mk_sig","NC_FLOAT", c("ID"))
    att.put.nc(nc, "shock_s_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"shock_w_mk_sig","NC_FLOAT", c("ID"))
    att.put.nc(nc, "shock_w_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"shock_y_mk_sig","NC_FLOAT", c("ID"))
    att.put.nc(nc, "shock_y_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "shock_s_mk_sig", per$shock_mk_sig[1:ntot,1]) 
    var.put.nc(nc, "shock_w_mk_sig", per$shock_mk_sig[1:ntot,2])
    var.put.nc(nc, "shock_y_mk_sig", per$shock_mk_sig[1:ntot,3])

    var.def.nc(nc,"shock_s_lr","NC_FLOAT", c("ID"))
    att.put.nc(nc, "shock_s_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"shock_w_lr","NC_FLOAT", c("ID"))
    att.put.nc(nc, "shock_w_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"shock_y_lr","NC_FLOAT", c("ID"))
    att.put.nc(nc, "shock_y_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "shock_s_lr", per$shock_lr[1:ntot,1]) 
    var.put.nc(nc, "shock_w_lr", per$shock_lr[1:ntot,2])
    var.put.nc(nc, "shock_y_lr", per$shock_lr[1:ntot,3])

    var.def.nc(nc,"shock_s_lr_sig","NC_FLOAT", c("ID"))
    att.put.nc(nc, "shock_s_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"shock_w_lr_sig","NC_FLOAT", c("ID"))
    att.put.nc(nc, "shock_w_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"shock_y_lr_sig","NC_FLOAT", c("ID"))
    att.put.nc(nc, "shock_y_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "shock_s_lr_sig", per$shock_lr_sig[1:ntot,1]) 
    var.put.nc(nc, "shock_w_lr_sig", per$shock_lr_sig[1:ntot,2])
    var.put.nc(nc, "shock_y_lr_sig", per$shock_lr_sig[1:ntot,3])

#---------------- markov persistence

    var.def.nc(nc,"year_warm","NC_FLOAT", c("ID","year"))
    att.put.nc(nc, "year_warm", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"year_cold","NC_FLOAT", c("ID","year"))
    att.put.nc(nc, "year_cold", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"sum_warm","NC_FLOAT", c("ID","year"))
    att.put.nc(nc, "sum_warm", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"sum_cold","NC_FLOAT", c("ID","year"))
    att.put.nc(nc, "sum_cold", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"win_warm","NC_FLOAT", c("ID","year"))
    att.put.nc(nc, "win_warm", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"win_cold","NC_FLOAT", c("ID","year"))
    att.put.nc(nc, "win_cold", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "year_warm", per$markov[1:ntot,1,])
    var.put.nc(nc, "year_cold", per$markov[1:ntot,2,])

    var.put.nc(nc, "win_warm", per$markov[1:ntot,3,])
    var.put.nc(nc, "win_cold", per$markov[1:ntot,4,])

    var.put.nc(nc, "sum_warm", per$markov[1:ntot,5,])
    var.put.nc(nc, "sum_cold", per$markov[1:ntot,6,])


    var.def.nc(nc, "sum_warm_mk", "NC_FLOAT", "ID")
    att.put.nc(nc, "sum_warm_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "sum_cold_mk", "NC_FLOAT", "ID")
    att.put.nc(nc, "sum_cold_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "win_warm_mk", "NC_FLOAT", "ID")
    att.put.nc(nc, "win_warm_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "win_cold_mk", "NC_FLOAT", "ID")
    att.put.nc(nc, "win_cold_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "year_warm_mk", "NC_FLOAT", "ID")
    att.put.nc(nc, "year_warm_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "year_cold_mk", "NC_FLOAT", "ID")
    att.put.nc(nc, "year_cold_mk", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "year_warm_mk",  per$markov_mk[1:ntot,1])
    var.put.nc(nc, "year_cold_mk",  per$markov_mk[1:ntot,2])

    var.put.nc(nc, "win_warm_mk",  per$markov_mk[1:ntot,3])
    var.put.nc(nc, "win_cold_mk",  per$markov_mk[1:ntot,4])

    var.put.nc(nc, "sum_warm_mk",  per$markov_mk[1:ntot,5])
    var.put.nc(nc, "sum_cold_mk",  per$markov_mk[1:ntot,6])


    var.def.nc(nc, "sum_warm_mk_sig", "NC_FLOAT", "ID")
    att.put.nc(nc, "sum_warm_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "sum_cold_mk_sig", "NC_FLOAT", "ID")
    att.put.nc(nc, "sum_cold_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "win_warm_mk_sig", "NC_FLOAT", "ID")
    att.put.nc(nc, "win_warm_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "win_cold_mk_sig", "NC_FLOAT", "ID")
    att.put.nc(nc, "win_cold_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "year_warm_mk_sig", "NC_FLOAT", "ID")
    att.put.nc(nc, "year_warm_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "year_cold_mk_sig", "NC_FLOAT", "ID")
    att.put.nc(nc, "year_cold_mk_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "year_warm_mk_sig",  per$markov_mk_sig[1:ntot,1])
    var.put.nc(nc, "year_cold_mk_sig",  per$markov_mk_sig[1:ntot,2])

    var.put.nc(nc, "win_warm_mk_sig",  per$markov_mk_sig[1:ntot,3])
    var.put.nc(nc, "win_cold_mk_sig",  per$markov_mk_sig[1:ntot,4])

    var.put.nc(nc, "sum_warm_mk_sig",  per$markov_mk_sig[1:ntot,5])
    var.put.nc(nc, "sum_cold_mk_sig",  per$markov_mk_sig[1:ntot,6])


    var.def.nc(nc, "sum_warm_lr", "NC_FLOAT", "ID")
    att.put.nc(nc, "sum_warm_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "sum_cold_lr", "NC_FLOAT", "ID")
    att.put.nc(nc, "sum_cold_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "win_warm_lr", "NC_FLOAT", "ID")
    att.put.nc(nc, "win_warm_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "win_cold_lr", "NC_FLOAT", "ID")
    att.put.nc(nc, "win_cold_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "year_warm_lr", "NC_FLOAT", "ID")
    att.put.nc(nc, "year_warm_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "year_cold_lr", "NC_FLOAT", "ID")
    att.put.nc(nc, "year_cold_lr", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "year_warm_lr",  per$markov_lr[1:ntot,1])
    var.put.nc(nc, "year_cold_lr",  per$markov_lr[1:ntot,2])

    var.put.nc(nc, "win_warm_lr",  per$markov_lr[1:ntot,3])
    var.put.nc(nc, "win_cold_lr",  per$markov_lr[1:ntot,4])

    var.put.nc(nc, "sum_warm_lr",  per$markov_lr[1:ntot,5])
    var.put.nc(nc, "sum_cold_lr",  per$markov_lr[1:ntot,6])


    var.def.nc(nc, "sum_warm_lr_sig", "NC_FLOAT", "ID")
    att.put.nc(nc, "sum_warm_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "sum_cold_lr_sig", "NC_FLOAT", "ID")
    att.put.nc(nc, "sum_cold_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "win_warm_lr_sig", "NC_FLOAT", "ID")
    att.put.nc(nc, "win_warm_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "win_cold_lr_sig", "NC_FLOAT", "ID")
    att.put.nc(nc, "win_cold_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "year_warm_lr_sig", "NC_FLOAT", "ID")
    att.put.nc(nc, "year_warm_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "year_cold_lr_sig", "NC_FLOAT", "ID")
    att.put.nc(nc, "year_cold_lr_sig", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "year_warm_lr_sig",  per$markov_lr_sig[1:ntot,1])
    var.put.nc(nc, "year_cold_lr_sig",  per$markov_lr_sig[1:ntot,2])

    var.put.nc(nc, "win_warm_lr_sig",  per$markov_lr_sig[1:ntot,3])
    var.put.nc(nc, "win_cold_lr_sig",  per$markov_lr_sig[1:ntot,4])

    var.put.nc(nc, "sum_warm_lr_sig",  per$markov_lr_sig[1:ntot,5])
    var.put.nc(nc, "sum_cold_lr_sig",  per$markov_lr_sig[1:ntot,6])



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