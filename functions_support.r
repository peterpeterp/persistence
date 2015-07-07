
r_calc_runmean_2D <- function(y2D,nday,nyr){
    nbuff1 = (nday-1)/2 
    nbuff2 = (nyr-1)/2 
    buffer = c(nbuff1,nbuff2)

    dims0 = dim(y2D)
    dims  = dim(y2D)
    dims[1] = dims[1]+2*nbuff1
    dims[2] = dims[2]+2*nbuff2

    y2Dex = array(NA,dim=c(dims[1],dims[2]))
    y2Dex[(nbuff1+1):(dims[1]-nbuff1),(nbuff2+1):(dims[2]-nbuff2)] = y2D[,]
    y2Dex[1:nbuff1,(nbuff2+2):(dims[2]-nbuff2+1)] = y2D[(dims0[1]-nbuff1+1):dims0[1],]   
    y2Dex[(dims[1]-nbuff1+1):dims[1],(nbuff2):(dims[2]-nbuff2-1)] = y2D[1:nbuff1,] 


    trend=y2D*NA

    for (i in 1:365){
        for (j in 1:62){
            trend[i,j]=mean(y2Dex[(i+nbuff1-nbuff1):(i+nbuff1+nbuff1),(j+nbuff2-nbuff2):(j+nbuff2+nbuff2)],na.rm=TRUE)
        }
    }
    cat("-")
    return(trend)
}

c_calc_runmean_2D <- function(y2D,nday,nyr)
{   # Input 2D array: y[day,year],
    # Calculate the running mean given a window of nyrs and ndays 
    # around each point in time
    nbuff1 = (nday-1)/2 
    nbuff2 = (nyr-1)/2 
    buffer = c(nbuff1,nbuff2)

    dims0 = dim(y2D)
    dims  = dim(y2D)
    dims[1] = dims[1]+2*nbuff1
    dims[2] = dims[2]+2*nbuff2

    # Fill in buffered array 
    y2D_list = array(y2D,dim=(dims0[1]*dims0[2]))
    y2D_list[is.na(y2D_list)]=0
    y2Dex_list = array(-99,dim=(dims[1]*dims[2])) 
    trend_list = array(1,dim=(dims0[1]*dims0[2]))       

    tempi = .C("c_run_mean2d",daten=as.numeric(y2D_list),datenex=as.numeric(y2Dex_list),temp_trend=as.numeric(trend_list),size=as.integer(dims0),tag=as.integer(buffer))
    trend = array(tempi$temp_trend,dim=c(dims0[1],dims0[2]))

    return(trend)
}

seasonal_1D_out <- function(dat,seasons=array(c(151,242,334,425),dim=c(2,2)),model=markov_calc,order=0,shift=0,interval=365){
    if (seasons[length(seasons)]>365){
        shift=seasons[length(seasons)]-365
        seasons[,]=seasons[,]-shift
    }
    size=length(dat)
    x=seq(1, size, 1)
    i=shift
    j=1
    out1=array(NA,dim=c(dim(seasons)[2],62))
    out2=array(NA,dim=c(dim(seasons)[2],62))
    out_err=array(NA,dim=c(dim(seasons)[2],62))
    out_add=array(NA,dim=c(dim(seasons)[2],62))

    while ((i+interval)<size){
        if ((is.na(dat[i+1])==FALSE) & (is.na(dat[i+interval])==FALSE)){
            for (sea in 1:length(seasons[1,])){
                x=dat[(seasons[1,sea]+i):(seasons[2,sea]+i)]
                tmp=model(x,order)

                out1[sea,j]=tmp$P_w
                out2[sea,j]=tmp$P_c
                out_err[sea,j]=tmp$error
                out_add[sea,j]=tmp$bic

            }
        }
        j=j+1
        i=i+interval
    }
    return(list(out1=out1,out2=out2,out_err=out_err,out_add=out_add))
}  

seasonal_matrix_out <- function(dat,seasons=array(c(151,242,334,425),dim=c(2,2)),shift=0,interval=365){
    if (seasons[length(seasons)]>365){
        shift=seasons[length(seasons)]-365
        seasons[,]=seasons[,]-shift
    }
    size=length(dat)
    x=seq(1, size, 1)
    i=shift
    j=1
    out=array(NA,dim=c(dim(seasons)[2],9,62))
    out_conf=array(NA,dim=c(dim(seasons)[2],62))

    while ((i+interval)<size){
        if ((is.na(dat[i+1])==FALSE) & (is.na(dat[i+interval])==FALSE)){
            for (sea in 1:length(seasons[1,])){
                x=dat[(seasons[1,sea]+i):(seasons[2,sea]+i)]
                tmp=markov_3states(x)

                out[sea,1:9,j]=tmp$transMat
                out_conf[sea,j]=tmp$confidence

            }
        }
        j=j+1
        i=i+interval
    }
    return(list(out=out,out_conf=out_conf))
}   



trend_analysis <- function(x,y){
    library(Kendall)
    if (length(which(is.na(y)))>7){
        return(list(slope=NA,slope_sig=NA,MK=NA,MK_sig=NA))
    }
    lm.r=lm(y~x)
    slope=summary(lm.r)$coefficients[2]
    slope_sig=summary(lm.r)$coefficients[8]
    out=MannKendall(y)
    MK=out[1]$tau
    MK_sig=out[2]$sl
    return(list(slope=slope,slope_sig=slope_sig,MK=MK,MK_sig=MK_sig))
}

global_trend <- function(filename_markov=99,filename_markov_neu=99,filename_shock=99,filename_shock_neu=99){
    t=seq(1,62,1)
    ntot=1319
    trend_markov=array(NA,dim=c(ntot,6,4))
    if (filename_markov!=99){   
        per=markov_load(filename_markov)
        for (i in 1:6){
            for (q in 1:ntot){
                tmp=trend_analysis(t,per$markov[q,i,])
                trend_markov[q,i,1]=tmp$slope
                trend_markov[q,i,2]=tmp$slope_sig
                trend_markov[q,i,3]=tmp$MK
                trend_markov[q,i,4]=tmp$MK_sig
            }
        }
        markov_trend_write(filename_markov_neu,trend_markov)
    }
    trend_shock=array(NA,dim=c(ntot,3,4))
    if (filename_shock!=99){      
        per=shock_load(filename_shock)
        for (i in 1:3){
            for (q in 1:ntot){
                tmp=trend_analysis(t,per$shock[q,i,])
                trend_shock[q,i,1]=tmp$slope
                trend_shock[q,i,2]=tmp$slope_sig
                trend_shock[q,i,3]=tmp$MK
                trend_shock[q,i,4]=tmp$MK_sig
            }
        }
        shock_trend_write(filename_shock_neu,trend_shock)
    }
    return(list(trend_markov=trend_markov,trend_shock=trend_shock))
}





