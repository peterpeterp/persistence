
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
    dat$time_2D = array(dat$time,dim=c(365,65))

    return(dat)
}

trend_load <- function(filename){
    nc  = open.nc(filename)  
    trend  = var.get.nc(nc,"trend")
    return(trend)
}

dat_load_precipitation <- function(filename){
    nc  = open.nc(filename)
    dat = list()     
    dat$day  = var.get.nc(nc,"day")
    dat$year = var.get.nc(nc,"year")
    dat$ID   = var.get.nc(nc,"ID")
    dat$lon   = var.get.nc(nc,"lon")
    dat$lat   = var.get.nc(nc,"lat")
    dat$pp   = var.get.nc(nc,"pp")

    # Add additional time dimension (days are decimal points dt=1/365)
    dat$time = c(1:(length(dat$day)*length(dat$year)))/365 - 0.5*1/365 + min(dat$year)
    dat$time_2D = array(dat$time,dim=c(365,length(dat$year)))

    return(dat)
}
