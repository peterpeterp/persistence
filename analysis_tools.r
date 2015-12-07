
quantile_pete <- function(dist,taus,na.rm=TRUE,plot=FALSE){
    # calculates quantiles from empirical cumulative distribution function -> gives out decimal numbers
    if (na.rm==TRUE){dist=dist[which(!is.na(dist))]}

    cdf=array(NA,max(dist))
    out=taus*NA

    for (i in 1:max(dist)){
        cdf[i]=length(which(dist<=i))
    }
    cdf=cdf/cdf[length(cdf)]

    for (qu in 1:length(taus)){
        ueb=which(cdf>taus[qu])[1]
        unt=ueb-1
        if (unt<1){out[qu]=ueb}
        else {out[qu]=ueb-(cdf[ueb]-taus[qu])/(cdf[ueb]-cdf[unt])}
    }  

    if (plot==TRUE){    
        #pdf(file="../plots/zwischenzeugs/quantile_test.pdf")
        plot(ecdf(dist))
        lines(cdf,col="red")
        for (qu in c(2,4,5)){
            abline(v=out[qu],col="blue")
            abline(h=taus[qu],col="blue")
            #print(quantile(dist,prob=taus[qu]))
            abline(v=quantile(dist,prob=taus[qu]),col="green",lty=2)
        }
    }
    return(out)
}

quantile_analysis <- function(x,y,taus){
    # what to do ... 
    print("nothing done here")
}

exponential_fit <- function(X,Y,start_guess=c(a=0.1,b=0.1),lower_limit=c(-Inf,-Inf),upper_limit=c(Inf,Inf)){
    xy=data.frame(y=Y,x=X)
    # try fit
    exp_nls=try(nls(y~(a*exp(-b*x)),data=xy,start=start_guess,lower=lower_limit,upper=upper_limit,na.action=na.exclude)) 
    
    # if succes
    if (class(exp_nls)!="try-error"){
        a=summary(exp_nls)$parameters[1]
        b=summary(exp_nls)$parameters[2]

        expfit=a*exp(-X*b)
        R2=1-sum(((Y-expfit)^2),na.rm=TRUE)/sum(((Y-mean(Y,na.rm=TRUE))^2),na.rm=TRUE)
        BIC=BIC(exp_nls)

        return(list(pars=c(a,b),ana=c(R2,BIC),fit=expfit))
    }
    # if fail
    else{
        return(list(pars=c(NA,NA),ana=c(NA,NA),fit=X*NA))
    }
}

two_exp_fit <- function(X,Y,start_guess=c(a=0.1,b=0.1),lower_limit=c(-Inf,-Inf),upper_limit=c(Inf,Inf)){
    xy=data.frame(y=Y,x=X)
    # try fit
    exp_nls=try(nls(y~(a*exp(-b*x)),data=xy,start=start_guess,lower=lower_limit,upper=upper_limit,na.action=na.exclude)) 
    
    # if succes
    if (class(exp_nls)!="try-error"){
        a=summary(exp_nls)$parameters[1]
        b=summary(exp_nls)$parameters[2]

        expfit=a*exp(-X*b)
        R2=1-sum(((Y-expfit)^2),na.rm=TRUE)/sum(((Y-mean(Y,na.rm=TRUE))^2),na.rm=TRUE)
        BIC=BIC(exp_nls)

        return(list(pars=c(a,b),ana=c(R2,BIC),fit=expfit))
    }
    # if fail
    else{
        return(list(pars=c(NA,NA),ana=c(NA,NA),fit=X*NA))
    }
}

