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

seasonal <- function(dat,shift=100,interval=365,seasons=c(0,365),model="ar_2"){
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
        if ((is.na(dat[i])==FALSE) & (is.na(dat[i+interval])==FALSE)){
            for (sea in 1:(length(seasons)-1)){
                x=dat[(seasons[sea]+i):(seasons[sea+1]+i)]
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

markov <- function(x){
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

shock_ma_1 <- function(x){
    arma_x=arima(x,order=c(0,1,1),method="ML")
    Px=arma_x$coef[1]+1
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

    temp=.C("c_markov_year",data=as.integer(per_ind1D),warm_y=as.numeric(markov[1,]),cold_y=as.numeric(markov[2,])
            ,start=as.integer(0),interval=as.integer(365),size=as.integer(size))
    markov[1,]=temp$warm_y
    markov[2,]=temp$cold_y

    seasons=c(182,183)
    start=100

    win_warm = array(99,62)
    win_cold = array(99,62)
    sum_warm = array(99,62)
    sum_cold = array(99,62)

    temp=.C("c_markov_season",data=as.integer(per_ind1D),warm_w=as.numeric(markov[3,]),cold_w=as.numeric(markov[4,]),warm_s=as.numeric(markov[5,]),cold_s=as.numeric(markov[6,])
            ,start=as.integer(start),season_length=as.integer(seasons),season_number=as.integer(length(seasons)),size=as.integer(size))

    markov[3,]=temp$warm_w
    markov[4,]=temp$cold_w
    markov[5,]=temp$warm_s
    markov[6,]=temp$cold_s

    markov_mk=array(NA,6)
    markov_mk_sig=array(NA,6)

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
    tmp=seasonal_bic(y,seasons=c(0,182,365),model="ar_2")
    shock[1,]=tmp$shock_s
    shock[2,]=tmp$shock_w
    bic[1,]=tmp$bic_s
    bic[2,]=tmp$bic_w
    tmp=seasonal_bic(y,seasons=c(0,365),model="ar_2")
    shock[3,]=tmp$shock_s
    bic[3,]=tmp$bic_s

    for (i in 1:3) {
        out=MannKendall(shock[i,])
        shock_mk[i]=out[1]$tau
        shock_mk_sig[i]=out[2]$sl
    }


    return(list(ind=per_ind,markov=markov,markov_mk=markov_mk,markov_mk_sig=markov_mk_sig,shock=shock,shock_mk=shock_mk,shock_mk_sig=shock_mk_sig,bic=bic))

}







