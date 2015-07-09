library(ncdf)

dat_load <- function(filename,reg=0){
    nc  = open.ncdf(filename)
    dat = list()     
    dat$day  = get.var.ncdf(nc,"day")
    dat$year = get.var.ncdf(nc,"year")
    dat$ID   = get.var.ncdf(nc,"ID")
    dat$lon   = get.var.ncdf(nc,"lon")
    dat$lat   = get.var.ncdf(nc,"lat")
    if (reg==1){ 
        dat$region   = get.var.ncdf(nc,"region")
    }
    dat$tas   = get.var.ncdf(nc,"tas")

    # Add additional time dimension (days are decimal points dt=1/365)
    dat$time = c(1:(length(dat$day)*length(dat$year)))/365 - 0.5*1/365 + min(dat$year)
    dat$time_2D = array(dat$time,dim=c(365,62))

    return(dat)
}

trend_load <- function(filename){
    nc  = open.ncdf(filename)  
    trend  = get.var.ncdf(nc,"trend")

    return(trend)
}

duration_load <- function(filename){
    nc  = open.ncdf(filename)
    dur = list() 

    dur$dur_warm  = get.var.ncdf(nc,"dur_warm")
    dur$dur_cold  = get.var.ncdf(nc,"dur_cold")
    dur$dur_warm_mid  = get.var.ncdf(nc,"dur_warm_mid")
    dur$dur_cold_mid  = get.var.ncdf(nc,"dur_cold_mid")
    return(dur)
}
