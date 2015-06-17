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
    nc  = open.ncdf(filename)  
    trend  = get.var.ncdf(nc,"trend")

    return(trend)
}


per_load <- function(filename){
    nc=open.ncdf(filename)
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
    year <- dim.def.ncdf("year",units="year",vals=1:62, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:length(data3D$ID), unlim=FALSE)

    vartrend <- var.def.ncdf(name="trend",units="bla",dim=list(ID,day,year), missval=-9999.0)
    nc = create.ncdf(filename,vartrend)
    put.var.ncdf(nc,vartrend,trend)

    close.ncdf(nc)  
}

markov_write <- function(filename,data3D,per)
{
    ntot=length(data3D$ID)
    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    year <- dim.def.ncdf("year",units="year",vals=1:62, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)


    ind <- var.def.ncdf(name="ind",units="bla",dim=list(ID,day,year), missval=-9999.0)
    markov_s_w <- var.def.ncdf(name="markov_s_w",units="bla",dim=list(ID,year), missval=-9999.0)
    markov_s_k <- var.def.ncdf(name="markov_s_k",units="bla",dim=list(ID,year), missval=-9999.0)
    markov_w_w <- var.def.ncdf(name="markov_w_w",units="bla",dim=list(ID,year), missval=-9999.0)
    markov_w_k <- var.def.ncdf(name="markov_w_k",units="bla",dim=list(ID,year), missval=-9999.0)
    markov_y_w <- var.def.ncdf(name="markov_y_w",units="bla",dim=list(ID,year), missval=-9999.0)
    markov_y_k <- var.def.ncdf(name="markov_y_k",units="bla",dim=list(ID,year), missval=-9999.0)

    markov_s_err <- var.def.ncdf(name="markov_s_err",units="bla",dim=list(ID,year), missval=-9999.0)
    markov_w_err <- var.def.ncdf(name="markov_w_err",units="bla",dim=list(ID,year), missval=-9999.0)
    markov_y_err <- var.def.ncdf(name="markov_y_err",units="bla",dim=list(ID,year), missval=-9999.0)

    vars=list(ind,markov_s_w,markov_s_k,markov_w_w,markov_w_k,markov_y_w,markov_y_k,markov_s_err,markov_w_err,markov_y_err)
   
    nc = create.ncdf(filename,vars)

    put.var.ncdf(nc,vars[[1]],per$ind)
    for (i in 1:6){
        put.var.ncdf(nc,vars[[i+1]],per$markov[1:ntot,i,])
    }
    for (i in 1:3){
        put.var.ncdf(nc,vars[[i+7]],per$markov_err[1:ntot,i,])
    }

    close.ncdf(nc) 
}

shock_write <- function(filename,data3D,per)
{
    ntot=length(data3D$ID)
    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    year <- dim.def.ncdf("year",units="year",vals=1:62, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)


    shock_s <- var.def.ncdf(name="shock_s",units="bla",dim=list(ID,year), missval=-9999.0)
    shock_w <- var.def.ncdf(name="shock_w",units="bla",dim=list(ID,year), missval=-9999.0)
    shock_y <- var.def.ncdf(name="shock_y",units="bla",dim=list(ID,year), missval=-9999.0)

    bic_s <- var.def.ncdf(name="bic_s",units="bla",dim=list(ID,year), missval=-9999.0)
    bic_w <- var.def.ncdf(name="bic_w",units="bla",dim=list(ID,year), missval=-9999.0)
    bic_y <- var.def.ncdf(name="bic_y",units="bla",dim=list(ID,year), missval=-9999.0)
    vars=list(shock_s,shock_w,shock_y,bic_s,bic_w,bic_y)
   
    nc = create.ncdf(filename,vars)

    for (i in 1:3){
        put.var.ncdf(nc,vars[[i]],per$shock[1:ntot,i,])
    }
    for (i in 1:3){
        put.var.ncdf(nc,vars[[i+3]],per$shock_bic[1:ntot,i,])
    }
    close.ncdf(nc) 
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