library(markovchain)

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
    if (sum(x,na.rm=TRUE)==length(x) | sum(-x,na.rm=TRUE)==length(x)){
        if (x[1]>0){
            P_w=1
            P_c=0
        }
        if (x[1]<0){
            P_w=0
            P_c=1
        }        
    }
    else {
        P_w=tmp$estimate[2][2]
        P_c=tmp$estimate[1][1]
    }
    return(list(P_w=P_w,P_c=P_c,error=tmp$confidenceInterval$confidenceLevel,bic=0))
}