library(stats)

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