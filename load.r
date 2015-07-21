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
    dat$time_2D = array(dat$time,dim=c(365,65))

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

markov_load <- function(filename,transitions){
    nc=open.ncdf(filename)
    ntot=1319

    markov_per = list(ind=array(NA,dim=c(ntot,365,65)),markov=array(NA,dim=c(ntot,5,transitions,65)),markov_err=array(NA,dim=c(ntot,5,65)))
    str_markov=c("markov_spring","markov_summer","markov_autumn","markov_winter","markov_year")
        
    for (i in 1:5){
        markov_per$markov[1:ntot,i,,]=get.var.ncdf(nc,str_markov[i])
    }
    str_markov_conf=c("markov_spring_conf","markov_summer_conf","markov_autumn_conf","markov_winter_conf","markov_year_conf")
    for (i in 1:5){
        markov_per$markov_err[1:ntot,i,]=get.var.ncdf(nc,str_markov_conf[i])
    }

    markov_per$ind <- get.var.ncdf(nc,"ind")

    return(markov_per)
}