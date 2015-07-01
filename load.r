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

markov_jjay_load <- function(filename){
    nc=open.ncdf(filename)
    ntot=1319

    markov_per = list(ind=array(NA,dim=c(ntot,365,62)),markov=array(NA,dim=c(ntot,8,62)),markov_err=array(NA,dim=c(ntot,4,62)))

    str_markov=c("markov_jn_w","markov_jl_w","markov_ag_w","markov_yr_w","markov_jn_c","markov_jl_c","markov_ag_c","markov_yr_c")
    for (i in 1:8){
        markov_per$markov[1:ntot,i,]=get.var.ncdf(nc,str_markov[i])
    }
    str_markov_err=c("markov_jn_err","markov_jl_err","markov_ag_err","markov_yr_err")
    for (i in 1:4){
        markov_per$markov_err[1:ntot,i,]=get.var.ncdf(nc,str_markov_err[i])
    }

    markov_per$ind <- get.var.ncdf(nc,"ind")

    return(markov_per)
}

markov_load <- function(filename){
    nc=open.ncdf(filename)
    ntot=1319

    markov_per = list(ind=array(NA,dim=c(ntot,365,62)),markov=array(NA,dim=c(ntot,6,62)),markov_err=array(NA,dim=c(ntot,3,62)))
    str_markov=c("markov_s_w","markov_w_w","markov_y_w",
        "markov_s_k","markov_w_k","markov_y_k")
    for (i in 1:6){
        markov_per$markov[1:ntot,i,]=get.var.ncdf(nc,str_markov[i])
    }
    str_markov_err=c("markov_s_err","markov_w_err","markov_y_err")
    for (i in 1:3){
        markov_per$markov_err[1:ntot,i,]=get.var.ncdf(nc,str_markov_err[i])
    }

    markov_per$ind <- get.var.ncdf(nc,"ind")

    return(markov_per)
}

shock_load <- function(filename){
    nc=open.ncdf(filename)
    ntot=1319

    shock_per = list(shock=array(NA,dim=c(ntot,6,62)),shock_bic=array(NA,dim=c(ntot,3,62)))
    str_shock=c("shock_s","shock_w","shock_y")
    for (i in 1:3){
        shock_per$shock[1:ntot,i,]=get.var.ncdf(nc,str_shock[i])
    }
    str_bic=c("bic_s","bic_w","bic_y")
    for (i in 1:3){
        shock_per$shock_bic[1:ntot,i,]=get.var.ncdf(nc,str_bic[i])
    }

    return(shock_per)
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
