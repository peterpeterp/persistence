library(ncdf)

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


    ind <- var.def.ncdf(name="ind",units="1 or -1",dim=list(ID,day,year), missval=-9999.0)
    markov_s_w <- var.def.ncdf(name="markov_s_w",units="0-1",longname="summer markov warm persistence",dim=list(ID,year), missval=-9999.0)
    markov_w_w <- var.def.ncdf(name="markov_w_w",units="0-1",longname="winter markov warm persistence",dim=list(ID,year), missval=-9999.0)
    markov_y_w <- var.def.ncdf(name="markov_y_w",units="0-1",longname="year markov warm persistence",dim=list(ID,year), missval=-9999.0)

    markov_s_k <- var.def.ncdf(name="markov_s_k",units="0-1",longname="summer markov cold persistence",dim=list(ID,year), missval=-9999.0)
    markov_w_k <- var.def.ncdf(name="markov_w_k",units="0-1",longname="winter markov cold persistence",dim=list(ID,year), missval=-9999.0)
    markov_y_k <- var.def.ncdf(name="markov_y_k",units="0-1",longname="year markov cold persistence",dim=list(ID,year), missval=-9999.0)

    markov_s_err <- var.def.ncdf(name="markov_s_err",units="0-1",dim=list(ID,year), missval=-9999.0)
    markov_w_err <- var.def.ncdf(name="markov_w_err",units="0-1",dim=list(ID,year), missval=-9999.0)
    markov_y_err <- var.def.ncdf(name="markov_y_err",units="0-1",dim=list(ID,year), missval=-9999.0)

    vars=list(ind,markov_s_w,markov_w_w,markov_y_w,
        markov_s_k,markov_w_k,markov_y_k,
        markov_s_err,markov_w_err,markov_y_err)
   
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

markov_jjay_write <- function(filename,data3D,per)
{
    ntot=length(data3D$ID)
    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    year <- dim.def.ncdf("year",units="year",vals=1:62, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)


    ind <- var.def.ncdf(name="ind",units="1 or -1",dim=list(ID,day,year), missval=-9999.0)
    markov_jn_w <- var.def.ncdf(name="markov_jn_w",units="0-1",longname="june markov warm persistence",dim=list(ID,year), missval=-9999.0)
    markov_jl_w <- var.def.ncdf(name="markov_jl_w",units="0-1",longname="july markov warm persistence",dim=list(ID,year), missval=-9999.0)
    markov_ag_w <- var.def.ncdf(name="markov_ag_w",units="0-1",longname="august markov warm persistence",dim=list(ID,year), missval=-9999.0)
    markov_yr_w <- var.def.ncdf(name="markov_yr_w",units="0-1",longname="year markov warm persistence",dim=list(ID,year), missval=-9999.0)
 
    markov_jn_c <- var.def.ncdf(name="markov_jn_c",units="0-1",longname="june markov cold persistence",dim=list(ID,year), missval=-9999.0)
    markov_jl_c <- var.def.ncdf(name="markov_jl_c",units="0-1",longname="july markov cold persistence",dim=list(ID,year), missval=-9999.0)
    markov_ag_c <- var.def.ncdf(name="markov_ag_c",units="0-1",longname="august markov cold persistence",dim=list(ID,year), missval=-9999.0)
    markov_yr_c <- var.def.ncdf(name="markov_yr_c",units="0-1",longname="year markov cold persistence",dim=list(ID,year), missval=-9999.0)

    markov_jn_err <- var.def.ncdf(name="markov_jn_err",units="0-1",dim=list(ID,year), missval=-9999.0)
    markov_jl_err <- var.def.ncdf(name="markov_jl_err",units="0-1",dim=list(ID,year), missval=-9999.0)
    markov_ag_err <- var.def.ncdf(name="markov_ag_err",units="0-1",dim=list(ID,year), missval=-9999.0)
    markov_yr_err <- var.def.ncdf(name="markov_yr_err",units="0-1",dim=list(ID,year), missval=-9999.0)

    vars=list(ind,markov_jn_w,markov_jl_w,markov_ag_w,markov_yr_w,markov_jn_c,markov_jl_c,markov_ag_c,markov_yr_c,markov_jn_err,markov_jl_err,markov_ag_err,markov_yr_err)
   
    nc = create.ncdf(filename,vars)

    put.var.ncdf(nc,vars[[1]],per$ind)
    for (i in 1:8){
        put.var.ncdf(nc,vars[[i+1]],per$markov[1:ntot,i,])
    }
    for (i in 1:4){
        put.var.ncdf(nc,vars[[i+9]],per$markov_err[1:ntot,i,])
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

    mar_s_w_lr <- var.def.ncdf(name="mar_s_w_lr",units="bla",longname="summer warm markov persistence linear regression",dim=list(ID), missval=-9999.0)
    mar_s_k_lr <- var.def.ncdf(name="mar_s_k_lr",units="bla",longname="summer cold markov persistence linear regression",dim=list(ID), missval=-9999.0)
    mar_w_w_lr <- var.def.ncdf(name="mar_w_w_lr",units="bla",longname="winter warm markov persistence linear regression",dim=list(ID), missval=-9999.0)
    mar_w_k_lr <- var.def.ncdf(name="mar_w_k_lr",units="bla",longname="winter cold markov persistence linear regression",dim=list(ID), missval=-9999.0)
    mar_y_w_lr <- var.def.ncdf(name="mar_y_w_lr",units="bla",longname="year warm markov persistence linear regression",dim=list(ID), missval=-9999.0)
    mar_y_k_lr <- var.def.ncdf(name="mar_y_k_lr",units="bla",longname="year cold markov persistence linear regression",dim=list(ID), missval=-9999.0)

    mar_s_w_lr_sig <- var.def.ncdf(name="mar_s_w_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_s_k_lr_sig <- var.def.ncdf(name="mar_s_k_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_w_w_lr_sig <- var.def.ncdf(name="mar_w_w_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_w_k_lr_sig <- var.def.ncdf(name="mar_w_k_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_y_w_lr_sig <- var.def.ncdf(name="mar_y_w_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    mar_y_k_lr_sig <- var.def.ncdf(name="mar_y_k_lr_sig",units="bla",dim=list(ID), missval=-9999.0)

    mar_s_w_mk <- var.def.ncdf(name="mar_s_w_mk",units="bla",longname="summer warm markov persistence Mann Kendall test",dim=list(ID), missval=-9999.0)
    mar_s_k_mk <- var.def.ncdf(name="mar_s_k_mk",units="bla",longname="summer cold markov persistence Mann Kendall test",dim=list(ID), missval=-9999.0)
    mar_w_w_mk <- var.def.ncdf(name="mar_w_w_mk",units="bla",longname="winter warm markov persistence Mann Kendall test",dim=list(ID), missval=-9999.0)
    mar_w_k_mk <- var.def.ncdf(name="mar_w_k_mk",units="bla",longname="winter cold markov persistence Mann Kendall test",dim=list(ID), missval=-9999.0)
    mar_y_w_mk <- var.def.ncdf(name="mar_y_w_mk",units="bla",longname="year warm markov persistence Mann Kendall test",dim=list(ID), missval=-9999.0)
    mar_y_k_mk <- var.def.ncdf(name="mar_y_k_mk",units="bla",longname="year cold markov persistence Mann Kendall test",dim=list(ID), missval=-9999.0)

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

    sho_s_lr <- var.def.ncdf(name="sho_s_lr",units="bla",longname="summer shock persistence linear regression",dim=list(ID), missval=-9999.0)
    sho_w_lr <- var.def.ncdf(name="sho_w_lr",units="bla",longname="winter shock persistence linear regression",dim=list(ID), missval=-9999.0)
    sho_y_lr <- var.def.ncdf(name="sho_y_lr",units="bla",longname="year shock persistence linear regression",dim=list(ID), missval=-9999.0)

    sho_s_lr_sig <- var.def.ncdf(name="sho_s_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    sho_w_lr_sig <- var.def.ncdf(name="sho_w_lr_sig",units="bla",dim=list(ID), missval=-9999.0)
    sho_y_lr_sig <- var.def.ncdf(name="sho_y_lr_sig",units="bla",dim=list(ID), missval=-9999.0)

    sho_s_mk <- var.def.ncdf(name="sho_s_mk",units="bla",longname="summer shock persistence Mann Kendall test",dim=list(ID), missval=-9999.0)
    sho_w_mk <- var.def.ncdf(name="sho_w_mk",units="bla",longname="winter shock persistence Mann Kendall test",dim=list(ID), missval=-9999.0)
    sho_y_mk <- var.def.ncdf(name="sho_y_mk",units="bla",longname="year shock persistence Mann Kendall test",dim=list(ID), missval=-9999.0)

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



duration_write <- function(filename,dur,len)
{
    ntot=819
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


duration_analysis_write <- function(filename,dur,season)
{
    ntot=819
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)

    dur_mean_w <- var.def.ncdf(name="dur_mean_w",units="days",longname=paste("mean duration of warm periods in",season),dim=list(ID), missval=-9999.0)
    dur_mean_c <- var.def.ncdf(name="dur_mean_c",units="days",longname=paste("mean duration of cold periods in",season),dim=list(ID), missval=-9999.0)
    
    dur_mean_d_w <- var.def.ncdf(name="dur_mean_d_w",units="days",longname=paste("difference between mean duration of warm periods in",season,"before and after 1980"),dim=list(ID), missval=-9999.0)
    dur_mean_d_c <- var.def.ncdf(name="dur_mean_d_c",units="days",longname=paste("difference between mean duration of cold periods in",season,"before and after 1980"),dim=list(ID), missval=-9999.0)

    dur_X_d_w <- var.def.ncdf(name="dur_X_d_w",units="#",longname=paste("difference in the amount of extremely long warm periods in",season,"before and after 1980"),dim=list(ID), missval=-9999.0)
    dur_X_d_c <- var.def.ncdf(name="dur_X_d_c",units="#",longname=paste("difference in the amount of extremely long cold periods in",season,"before and after 1980"),dim=list(ID), missval=-9999.0)

    vars=list(dur_mean_w,dur_mean_c,dur_mean_d_w,dur_mean_d_c,dur_X_d_w,dur_X_d_c)
   
    nc = create.ncdf(paste(filename,season,".nc",sep=""),vars)

    for (i in 1:6){
        put.var.ncdf(nc,vars[[i]],dur[1:ntot,i])      
    }

    close.ncdf(nc) 
}