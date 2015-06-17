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


markov_load <- function(filename){
    nc=open.ncdf(filename)
    ind   = get.var.ncdf(nc,"ind")
    ntot=dim(ind)[1]

    markov_per = list(ind=array(NA,dim=c(ntot,365,62)),markov=array(NA,dim=c(ntot,6,62)),markov_err=array(NA,dim=c(ntot,3,62)))
    str_markov=c("markov_s_w","markov_s_k","markov_w_w","markov_w_k","markov_y_w","markov_y_k")
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
    ntot=819

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


markov_trend_write <- function(filename,per)
{
    ntot=819
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)

    mar_s_w_lr <- var.def.ncdf(name="mar_s_w_lr",units="bla",dim=list(ID), missval=-9999.0)
    mar_s_k_lr <- var.def.ncdf(name="mar_s_k_lr",units="bla",dim=list(ID), missval=-9999.0)
    mar_w_w_lr <- var.def.ncdf(name="mar_w_w_lr",units="bla",dim=list(ID), missval=-9999.0)
    mar_w_k_lr <- var.def.ncdf(name="mar_w_k_lr",units="bla",dim=list(ID), missval=-9999.0)
    mar_y_w_lr <- var.def.ncdf(name="mar_y_w_lr",units="bla",dim=list(ID), missval=-9999.0)
    mar_y_k_lr <- var.def.ncdf(name="mar_y_k_lr",units="bla",dim=list(ID), missval=-9999.0)

    mar_s_w_lr_sig <- var.def.ncdf(name="mar_s_w_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_s_k_lr_sig <- var.def.ncdf(name="mar_s_k_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_w_w_lr_sig <- var.def.ncdf(name="mar_w_w_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_w_k_lr_sig <- var.def.ncdf(name="mar_w_k_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_y_w_lr_sig <- var.def.ncdf(name="mar_y_w_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_y_k_lr_sig <- var.def.ncdf(name="mar_y_k_lr_sig",units="bla",dim=list(ID), missval=-9999.0)

    mar_s_w_mk <- var.def.ncdf(name="mar_s_w_mk",units="bla",dim=list(ID), missval=-9999.0)
    mar_s_k_mk <- var.def.ncdf(name="mar_s_k_mk",units="bla",dim=list(ID), missval=-9999.0)
    mar_w_w_mk <- var.def.ncdf(name="mar_w_w_mk",units="bla",dim=list(ID), missval=-9999.0)
    mar_w_k_mk <- var.def.ncdf(name="mar_w_k_mk",units="bla",dim=list(ID), missval=-9999.0)
    mar_y_w_mk <- var.def.ncdf(name="mar_y_w_mk",units="bla",dim=list(ID), missval=-9999.0)
    mar_y_k_mk <- var.def.ncdf(name="mar_y_k_mk",units="bla",dim=list(ID), missval=-9999.0)

    mar_s_w_mk_sig <- var.def.ncdf(name="mar_s_w_mk_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_s_k_mk_sig <- var.def.ncdf(name="mar_s_k_mk_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_w_w_mk_sig <- var.def.ncdf(name="mar_w_w_mk_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_w_k_mk_sig <- var.def.ncdf(name="mar_w_k_mk_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_y_w_mk_sig <- var.def.ncdf(name="mar_y_w_mk_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_y_k_mk_sig <- var.def.ncdf(name="mar_y_k_mk_sig",units="bla",dim=list(ID), missval=-9999.0)

    vars=list(mar_s_w_lr,mar_s_k_lr,mar_w_w_lr,mar_w_k_lr,mar_y_w_lr,mar_y_k_lr,
        mar_s_w_lr_sig,mar_s_k_lr_sig,mar_w_w_lr_sig,mar_w_k_lr_sig,mar_y_w_lr_sig,mar_y_k_lr_sig,
        mar_s_w_mk,mar_s_k_mk,mar_w_w_mk,mar_w_k_mk,mar_y_w_mk,mar_y_k_mk,
        mar_s_w_mk_sig,mar_s_k_mk_sig,mar_w_w_mk_sig,mar_w_k_mk_sig,mar_y_w_mk_sig,mar_y_k_mk_sig)
   
    nc = create.ncdf(filename,vars)

    for (j in 1:4){
        for (i in 1:6){
            put.var.ncdf(nc,vars[[i+(j-1)*6]],per[1:ntot,i,j])
        }        
    }

    close.ncdf(nc) 
}

shock_trend_write <- function(filename,per)
{
    ntot=819
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)

    sho_s_lr <- var.def.ncdf(name="sho_s_lr",units="bla",dim=list(ID), missval=-9999.0)
    sho_w_lr <- var.def.ncdf(name="sho_w_lr",units="bla",dim=list(ID), missval=-9999.0)
    sho_y_lr <- var.def.ncdf(name="sho_y_lr",units="bla",dim=list(ID), missval=-9999.0)

    sho_s_lr_sig <- var.def.ncdf(name="sho_s_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    sho_w_lr_sig <- var.def.ncdf(name="sho_w_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    sho_y_lr_sig <- var.def.ncdf(name="sho_y_lr_sig",units="bla",dim=list(ID), missval=-9999.0)

    sho_s_mk <- var.def.ncdf(name="sho_s_mk",units="bla",dim=list(ID), missval=-9999.0)
    sho_w_mk <- var.def.ncdf(name="sho_w_mk",units="bla",dim=list(ID), missval=-9999.0)
    sho_y_mk <- var.def.ncdf(name="sho_y_mk",units="bla",dim=list(ID), missval=-9999.0)

    sho_s_mk_sig <- var.def.ncdf(name="sho_s_mk_sig",units="bla",dim=list(ID), missval=-9999.0)
    sho_w_mk_sig <- var.def.ncdf(name="sho_w_mk_sig",units="bla",dim=list(ID), missval=-9999.0)
    sho_y_mk_sig <- var.def.ncdf(name="sho_y_mk_sig",units="bla",dim=list(ID), missval=-9999.0)

    vars=list(sho_s_lr,sho_w_lr,sho_y_lr,
        sho_s_lr_sig,sho_w_lr_sig,sho_y_lr_sig,
        sho_s_mk,sho_w_mk,sho_y_mk,
        sho_s_mk_sig,sho_w_mk_sig,sho_y_mk_sig)
   
    nc = create.ncdf(filename,vars)

    for (j in 1:4){
        for (i in 1:3){
            put.var.ncdf(nc,vars[[i+(j-1)*3]],per[1:ntot,i,j])
        }        
    }

    close.ncdf(nc) 
}