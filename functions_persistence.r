library(Kendall)
library(stats)

c_calc_runmean_2D = function(y2D,nday,nyr)
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

bic_selective <- function(x){
    bic_tmp=c(bayesian_i_c(x,c(1,1,0)),bayesian_i_c(x,c(2,1,0)),bayesian_i_c(x,c(0,1,1)))
    if (which(bic_tmp==min(bic_tmp))==1){
        return(list(P_w=shock_ar_1(x)$P_w,P_c=NA,bic=1))
    }
    if (which(bic_tmp==min(bic_tmp))==2){
        return(list(P_w=shock_ar_2(x)$P_w,P_c=NA,bic=2))
    }
    if (which(bic_tmp==min(bic_tmp))==3){
        return(list(P_w=shock_ma_1(x)$P_w,P_c=NA,bic=3))
    }
    return(list(P_w=NA,P_c=NA,bic=NA))

}

seasonal <- function(dat,seasons=array(c(151,242,334,425),dim=c(2,2)),model=markov_calc,shift=0,interval=365){
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
    bic_s=array(NA,62)
    bic_w=array(NA,62)
    while ((i+interval)<size){
        if ((is.na(dat[i+1])==FALSE) & (is.na(dat[i+interval])==FALSE)){
            for (sea in 1:length(seasons[1,])){
                x=dat[(seasons[1,sea]+i):(seasons[2,sea]+i)]
                tmp=model(x)
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
            }
        }
        j=j+1
        i=i+interval
    }
    return(list(summer_w=summer_w,summer_c=summer_c,winter_w=winter_w,winter_c=winter_c,bic_s=bic_s,bic_w=bic_w))
}   

bayesian_i_c <- function(x,order){
    arma_x=arima(x,order=order,method="ML")
    bic=-2*arma_x$loglik+(sum(order)-1)*log(length(x))
    return(bic)
}

markov_calc <- function(x){
    tmp=.C("c_markov",data=as.integer(x),warm=as.numeric(0.0),cold=as.numeric(0.0),size=as.integer(length(x)))
    if (tmp$warm == 99){
        tmp$warm <- NA
    }
    if (tmp$cold == 99){
        tmp$cold <- NA
    }
    return(list(P_w=tmp$warm,P_c=tmp$cold,bic=0))
}

shock_ar_1 <- function(x){
    arma_x=arima(x,order=c(1,1,0),method="ML")
    ar_x=arma_x$coef[1]
    Px=1
    psi=array(0,100)
    for (i in 1:100){
        Px=Px+ar_x^i
    }
    return(list(P_w=Px,P_c=NA,bic=-2*arma_x$loglik+(1*log(length(x)))))
}

shock_ar_2 <- function(x){
    arma_x=arima(x,order=c(2,1,0),method="ML")
    ar_x=arma_x$coef[1:2]
    wurzel=sqrt(ar_x[2]^2+4*ar_x[1])
    y1=(-ar_x[2]+wurzel)/(2*ar_x[1])
    y2=(-ar_x[2]-wurzel)/(2*ar_x[1])
    la1=1/y1
    la2=1/y2
    c1=la1/(la1-la2)
    c2=la2/(la2-la1)
    Px=1
    psi=array(0,100)
    for (i in 1:100){
        for (k in 0:i){
            psi[i]=c1*la1^k+c2*la2^k
        }
        Px=Px+psi[i]
    }
    return(list(P_w=Px,P_c=NA,bic=-2*arma_x$loglik+(2*log(length(x)))))
}


shock_ma_1 <- function(x){
    ma_order=1
    arma_x=arima(x,order=c(0,1,ma_order),method="ML")
    Px=1
    for (i in 1:ma_order){
        Px=Px+arma_x$coef[i]
    }
    return(list(P_w=Px,P_c=NA,bic=-2*arma_x$loglik+(3*log(length(x)))))
}

shock_ma_2 <- function(x){
    ma_order=2
    arma_x=arima(x,order=c(0,1,ma_order),method="ML")
    Px=1
    for (i in 1:ma_order){
        Px=Px+arma_x$coef[i]
    }
    return(list(P_w=Px,P_c=NA,bic=-2*arma_x$loglik+(3*log(length(x)))))
}

shock_ma_3 <- function(x){
    ma_order=3
    arma_x=arima(x,order=c(0,1,ma_order),method="ML")
    Px=1
    for (i in 1:ma_order){
        Px=Px+arma_x$coef[i]
    }
    return(list(P_w=Px,P_c=NA,bic=-2*arma_x$loglik+(3*log(length(x)))))
}


calc_persistence = function(y,time,trend=NULL,nday=91,nyr=5,trash=(365*2+61))
{   # Given 2D inputs of y[nday,nyr] and the smooth trend[nday,nyr],
    # calculate the persistence index (-1,0,1) and the persistence
    # vector per[nevents,nyr] and the event_day[nevents,nyr]
    # time should be a decimal-based vector 

    if (length(y) != length(time)) {
        cat("calc_persistence:: error: length of time vector does not match length of y.\n")
        return(NULL)
    }

    time0=proc.time()[1]

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
    markov_mk=array(NA,6)
    markov_mk_sig=array(NA,6)





    tmp=seasonal(per_ind1D,array(c(1,365),dim=c(2,1)),markov_calc)
    markov[1,]=tmp$summer_w
    markov[2,]=tmp$summer_c

    tmp=seasonal(per_ind1D,array(c(151,242,334,425),dim=c(2,2)),markov_calc)
    markov[3,]=tmp$summer_w
    markov[4,]=tmp$summer_c
    markov[5,]=tmp$winter_w
    markov[6,]=tmp$winter_c

    for (i in 1:6){
        markov[i,][markov[i,] == 99]=NA
        out=MannKendall(markov[i,])
        markov_mk[i]=out[1]$tau
        markov_mk_sig[i]=out[2]$sl
    }

    shock=array(NA,dim=c(3,62))
    bic=array(NA,dim=c(3,62))
    shock_mk=array(NA,3)
    shock_mk_sig=array(NA,3)

    tmp=seasonal(y,array(c(1,365),dim=c(2,1)),shock_ma_3)
    shock[1,]=tmp$summer_w
    bic[1,]=tmp$bic_s
    tmp=seasonal(y,array(c(151,242,334,425),dim=c(2,2)),shock_ma_3)
    shock[2,]=tmp$summer_w
    bic[2,]=tmp$bic_s
    shock[3,]=tmp$winter_w
    bic[3,]=tmp$bic_w

    for (i in 1:3) {
        out=MannKendall(shock[i,])
        shock_mk[i]=out[1]$tau
        shock_mk_sig[i]=out[2]$sl
    }


    return(list(ind=per_ind,markov=markov,markov_mk=markov_mk,markov_mk_sig=markov_mk_sig,shock=shock,shock_mk=shock_mk,shock_mk_sig=shock_mk_sig,bic=bic))

}







