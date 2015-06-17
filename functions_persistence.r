library(stats)
library(markovchain)

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

bic_selective <- function(x,order){
    shock=array(NA,7)
    bic=array(NA,7)
#    for (i in 1:5){
#        tmp=shock_ar(as.vector(x),i)
#        shock[i]=tmp$P_w
#        bic[i]=tmp$bic
#    }   
    for (ma in 1:7){
        tmp=shock_ma(as.vector(x),ma)
        shock[ma]=tmp$P_w
        bic[ma]=tmp$bic
    }
    return(list(P_w=shock[which(bic==min(bic,na.rm=TRUE))],P_c=NA,error=0,bic=which(bic==min(bic,na.rm=TRUE))))

}

seasonal <- function(dat,seasons=array(c(151,242,334,425),dim=c(2,2)),model=markov_calc,order=0,shift=0,interval=365){
    if (seasons[length(seasons)]>365){
        shift=seasons[length(seasons)]-365
        seasons[,]=seasons[,]-shift
    }
    size=length(dat)
    x=seq(1, size, 1)
    i=shift
    j=1
    summer_w=array(NA,62)
    summer_c=array(NA,62)
    winter_w=array(NA,62)
    winter_c=array(NA,62)
    error_s=array(NA,62)
    error_w=array(NA,62)
    bic_s=array(NA,62)
    bic_w=array(NA,62)
    while ((i+interval)<size){
        if ((is.na(dat[i+1])==FALSE) & (is.na(dat[i+interval])==FALSE)){
            for (sea in 1:length(seasons[1,])){
                x=dat[(seasons[1,sea]+i):(seasons[2,sea]+i)]
                tmp=model(x,order)
                if (sea==1){
                    summer_w[j]=tmp$P_w
                    summer_c[j]=tmp$P_c
                }
                if (sea==2){
                    winter_w[j+1]=tmp$P_w
                    winter_c[j+1]=tmp$P_c
                }
                if (tmp$bic!=0){
                    if (sea==1){
                        bic_s[j]=tmp$bic
                    }
                    if (sea==2){
                        bic_w[j+1]=tmp$bic
                    }             
                }
                if (tmp$error!=0){
                    if (sea==1){
                        error_s[j]=tmp$error
                    }
                    if (sea==2){
                        error_w[j+1]=tmp$error
                    }             
                }
            }
        }
        j=j+1
        i=i+interval
    }
    return(list(summer_w=summer_w,summer_c=summer_c,winter_w=winter_w,winter_c=winter_c,error_s=error_s,error_w=error_w,bic_s=bic_s,bic_w=bic_w))
}   


markov_calc <- function(x,order){
    tmp=.C("c_markov",data=as.integer(x),warm=as.numeric(0.0),cold=as.numeric(0.0),size=as.integer(length(x)))
    if (tmp$warm == 99){
        tmp$warm <- NA
    }
    if (tmp$cold == 99){
        tmp$cold <- NA
    }
    return(list(P_w=tmp$warm,P_c=tmp$cold,error=0,bic=0))
}

markov_chft <- function(x,order){
    tmp=markovchainFit(data=x)
    return(list(P_w=tmp$estimate[2][2],P_c=tmp$estimate[1][1],error=tmp$confidenceInterval$confidenceLevel,bic=0))
}


shock_ar_1 <- function(x,order){
    arma_x=arima(x,order=c(1,1,0),method="ML")
    ar_x=arma_x$coef[1]
    Px=1
    psi=array(0,100)
    for (i in 1:100){
        Px=Px+ar_x^i
    }
    return(list(P_w=Px,P_c=NA,error=0,bic=-2*arma_x$loglik+(1*log(length(x)))))
}

shock_ar <- function(x,ar_order){
    if (ar_order==1){
        return(shock_ar_1(x))
    }
    arma_x=arima(x,order=c(ar_order,1,0),method="ML")

    A <- matrix(0,nrow=ar_order,ncol=ar_order)
    A[1,]=arma_x$coef[1:ar_order]
    for (i in 2:ar_order){
        A[i,(i-1)]=1
    }
    Px=1+A[1]
    amul=A
    for (i in 1:100){
        amul=A %*% amul
        Px=Px+amul[1]
    }
    return(list(P_w=Px,P_c=NA,error=0,bic=-2*arma_x$loglik+ar_order*log(length(x))))
}


shock_ma <- function(x,ma_order){
    arma_x=arima(x,order=c(0,1,ma_order),method="ML")
    Px=1
    for (i in 1:ma_order){
        Px=Px+arma_x$coef[i]
    }
    return(list(P_w=Px,P_c=NA,error=0,bic=-2*arma_x$loglik+ma_order*log(length(x))))
}



calc_markov_per <- function(y,time,trend=NULL,nday=91,nyr=5,trash=(365*2+61))
{   # Given 2D inputs of y[nday,nyr] and the smooth trend[nday,nyr],
    # calculate the persistence index (-1,0,1) and the persistence
    # vector per[nevents,nyr] and the event_day[nevents,nyr]
    # time should be a decimal-based vector 

    if (length(y) != length(time)) {
        cat("calc_persistence:: error: length of time vector does not match length of y.\n")
        return(NULL)
    }

    # Calculate persistence vector
    per_ind = y*0 


    per_ind[y >= trend]  =  1
    per_ind[y <  trend]  = -1 
    per_ind[is.na(per_ind)] = 0  # To avoid NA values

    # Go through vector and calculate persistent events 
    # Perform this step on a 1D vector to avoid artificial cutoffs 
    # at day 1 and day 365 of the year 
    per_ind1D = as.vector(per_ind) 

    size=length(per_ind1D)

    markov=array(99,dim=c(6,62))
    markov_err=array(99,dim=c(3,62))

    tmp=seasonal(per_ind1D,array(c(151,242,334,425),dim=c(2,2)),markov_chft)
    markov[1,]=tmp$summer_w
    markov[2,]=tmp$summer_c
    markov[3,]=tmp$winter_w
    markov[4,]=tmp$winter_c
    markov_err[1,]=tmp$error_s
    markov_err[2,]=tmp$error_s

    tmp=seasonal(per_ind1D,array(c(1,365),dim=c(2,1)),markov_chft)
    markov[5,]=tmp$summer_w
    markov[6,]=tmp$summer_c
    markov_err[3,]=tmp$error_s

    return(list(ind=per_ind,markov=markov,markov_err=markov_err))

}

calc_shock_per <- function(y,time,trend=NULL,nday=91,nyr=5,trash=(365*2+61))
{   # Given 2D inputs of y[nday,nyr] and the smooth trend[nday,nyr],
    # calculate the persistence index (-1,0,1) and the persistence
    # vector per[nevents,nyr] and the event_day[nevents,nyr]
    # time should be a decimal-based vector 

    if (length(y) != length(time)) {
        cat("calc_persistence:: error: length of time vector does not match length of y.\n")
        return(NULL)
    }

    shock=array(NA,dim=c(3,62))
    shock_bic=array(NA,dim=c(3,62))

    tmp=seasonal((y-trend),array(c(151,242,334,425),dim=c(2,2)),shock_ma,3)
    shock[1,]=tmp$summer_w
    shock_bic[1,]=tmp$bic_s
    shock[2,]=tmp$winter_w
    shock_bic[2,]=tmp$bic_w

    tmp=seasonal((y-trend),array(c(1,365),dim=c(2,1)),shock_ma,3)
    shock[3,]=tmp$summer_w
    shock_bic[3,]=tmp$bic_s

    return(list(ind=per_ind,markov=markov,markov_err=markov_err,shock=shock,shock_bic=shock_bic))
}






