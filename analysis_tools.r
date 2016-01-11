library(stats4)

quantile_pete <- function(dist,taus,na.rm=TRUE,plot=FALSE){
    # calculates quantiles from empirical cumulative distribution function -> gives out decimal numbers
    if (na.rm==TRUE){dist=dist[which(!is.na(dist))]}

    # calculate cdf
    cdf=array(NA,max(dist))
    for (i in 1:max(dist)){
        cdf[i]=length(which(dist<=i))
    }
    cdf=cdf/cdf[length(cdf)]

    # get quantiles
    out=taus*NA
    for (qu in 1:length(taus)){
        if (taus[qu]==1){out[qu]=max(dist)} # special case
        else{
            ueb=which(cdf>taus[qu])[1]
            unt=ueb-1
            if (unt<1){out[qu]=ueb}
            else {out[qu]=ueb-(cdf[ueb]-taus[qu])/(cdf[ueb]-cdf[unt])}
        }
    }  

    # plot for quality check
    if (plot==TRUE){    
        pdf(file="../plots/zwischenzeugs/quantile_test.pdf")
        plot(ecdf(dist))
        lines(cdf,col="red")
        for (qu in c(2,4,5)){
            abline(v=out[qu],col="blue")
            abline(h=taus[qu],col="blue")
            abline(v=quantile(dist,prob=taus[qu]),col="green",lty=2)
        }
    }
    return(out)
}

quantile_analysis <- function(x,y,taus,noise_level=c(0.000001,0.0001)){
    # x=dur_mid, y=dur, taus=percentiles

    # noise is added to x and y to avoid crash
    x=x+rnorm(length(x),mean=0,sd=1)*noise_level[1]
    y=y+rnorm(length(y),mean=0,sd=1)*noise_level[2]
    
    quantiles=quantile_pete(y,taus=taus,na.rm=TRUE)
    slopes=taus*NA
    slope_sigs=taus*NA

    # try quantile regression for each tau individually 
    for (i in 1:length(taus)){
        quant_zwi=try(summary(rq(y~x,taus[i]),se="boot"),silent=TRUE)
        if (class(quant_zwi)!="try-error"){
            slopes[i]=quant_zwi$coefficients[2]
            slope_sigs[i]=quant_zwi$coefficients[8]
        }
    }
    return(list(quantiles=quantiles,slopes=slopes,slope_sigs=slope_sigs))
}

fit_plot_empty <- function(){
    # creates empty plots with x-axis and y-axis
    par(mar=c(3, 3, 3, 3) + 0.1)
    # x
    plot(NA,xlab="",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=FALSE,cex=0.5)
    at_=axis(1,labels=FALSE,col="black")
    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
    axis(1,at=at_)

    # x
    plot(NA,xlab="",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=FALSE,cex=0.5)
    at_=axis(3,labels=FALSE,col="black")
    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
    axis(3,at=at_)

    # y
    plot(NA,xlab="",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=FALSE,cex=0.5)
    at_=axis(2,labels=FALSE,col="black")
    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
    axis(2,at=at_)

    # y
    plot(NA,xlab="",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=FALSE,cex=0.5)
    at_=axis(4,labels=FALSE,col="black")
    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
    axis(4,at=at_)

    # y
    plot(NA,xlab="",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=FALSE,log="y",cex=0.5)
    at_=axis(2,labels=FALSE,col="black")
    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
    axis(2,at=at_)

    # y
    plot(NA,xlab="",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=FALSE,log="y",cex=0.5)
    at_=axis(4,labels=FALSE,col="black")
    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
    axis(4,at=at_)
}

fit_plot_reference <- function(x,y,sea,q,state){
    state_names=c("cold","warm")

    #reference
    par(mar=c(0, 0, 0, 0) + 0.1)
    plot(NA,xlim=c(0,10),ylim=c(0,10),axes=FALSE,frame.plot=FALSE)
    text(5,9.4,paste(sea,q,state_names[state]))

    #extreme events
    plot(NA,xlim=c(0,10),ylim=c(0,10),axes=FALSE,frame.plot=FALSE)
    order_=order(y,decreasing=TRUE)
    for (extreme in 1:10){
        text(3,extreme,round(x[order_][extreme],01))
        text(8,extreme,round(y[order_][extreme],01))
    }
}

fit_plot_combi <- function(X,Y,counts,expfit,fit,fitstuff,fit_style,sea,q,state){
    state_names=c("cold","warm")
    color=c("blue","red",rgb(0.5,0.1,0.5),rgb(0.1,0.7,0.1),rgb(0,0,0))

    par(mar=c(3, 3, 3, 3) + 0.1)  

    # second plot page
    nonull=which(counts>0)              
    plot(X[nonull],counts[nonull],xlab="days",xlim=c(0,70),ylab="counts",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],cex=0.5)
    axis(2)
    legend("topright",legend=c(paste("     ID=",q,sep=""),paste("dBIC=",round(fitstuff[17],01),sep="")),bty="n")

    # second plot page
    par(mar=c(3, 3, 3, 3) + 0.1)                
    plot(X[nonull],Y[nonull],xlab="days",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],cex=0.5)
    if (!is.na(fitstuff[9])){
        abline(v=fitstuff[9],col="grey")
        #text(thresh,0.22,label=round(thresh,02),col="grey")
    }
    lines(expfit,col=color[3],lty=2)
    lines(fit,col=color[4],lty=1)               
 
    text(50,0.18+0.06,paste("BIC=",round(fitstuff[16],01)),pos=1,col=color[3])                 
    text(50,0.155+0.06,paste("BIC=",round(fitstuff[20],01)),pos=1,col=color[4])      

    text(50,0.11+0.06,paste("R2=",round(fitstuff[15],03)),pos=1,col=color[3])                 
    text(50,0.085+0.06,paste("R2=",round(fitstuff[19],03)),pos=1,col=color[4])            

    # third version
    plot(X[nonull],Y[nonull],xlab="days",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],log="y",cex=0.5)
    if (!is.na(fitstuff[9])){
        abline(v=fitstuff[9],col="grey")
        text(fitstuff[9],0.00002,label=round(fitstuff[9],02),col=rgb(0,0,0))
    }    
    lines(expfit,col=color[3],lty=2)
    lines(fit,col=color[4],lty=1)
    text(50,0.22,paste("P=",round(exp(-fitstuff[2])*100,01)),pos=1,col=color[3])                 
    text(50,0.05,paste("P1=",round(exp(-fitstuff[6])*100,01)),pos=1,col=color[4])                 
    text(50,0.02,paste("P2=",round(exp(-fitstuff[8])*100,01)),pos=1,col=color[4])  
               
}

fit_plot <- function(X,Y,fit,legend,fit_style,sea,q,state,thresh=NA){
    state_names=c("cold","warm")
    color=c("blue","red")

    #first plot page
    par(mar=c(3, 3, 3, 3) + 0.1)                
    plot(X,Y,xlab="days",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],cex=0.5)
    if (!is.na(thresh)){
        abline(v=thresh,col="grey")
        #text(thresh,0.22,label=round(thresh,02),col="grey")
    }
    lines(fit,col="black")
    text(50,0.24,legend[1],pos=1)                 

    # second version
    plot(X,Y,xlab="days",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],log="y",cex=0.5)

    lines(fit,col="black")
    if (!is.na(thresh)){
        abline(v=thresh,col="grey")
        text(thresh,0.00002,label=round(thresh,02),col=rgb(0,0,0))
    }    
    text(50,0.24,legend[2],pos=1)                 
}

exponential_fit <- function(X,Y,start_guess=c(a=0.1,b=0.1),lower_limit=c(-Inf,-Inf),upper_limit=c(Inf,Inf)){
    xy=data.frame(y=Y,x=X)
    # try fit
    exp_nls=try(nls(y~(a*exp(-b*x)),data=xy,algorithm="port",start=start_guess,lower=lower_limit,upper=upper_limit,na.action=na.exclude),silent=TRUE) 
    
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

overlap_expo <-function(x,a1,b1,a2,b2){
    a1*exp(-b1*x)+a2*exp(-b2*x)
}

combi_expo <-function(x,a1,b1,b2,thresh){
    y=x*NA
    y[x<=thresh]=a1*exp(-x[x<=thresh]*b1)
    y[x>thresh]=a1*exp((b2-b1)*thresh)*exp(-(x[x>thresh])*b2)
    return(y)
}

combi_expo_restraint <-function(x,a1,b1,b2,thresh){
    if (b2>b1){return(x*0)}
    else {
        y=x*NA
        y[x<=thresh]=a1*exp(-x[x<=thresh]*b1)
        y[x>thresh]=a1*exp((b2-b1)*thresh)*exp(-(x[x>thresh])*b2)
        return(y)
    }
}

to_be_minimized <- function(data, par){
    if (par[3]>par[2]){return(999999)}
    if (par[4]>=25 | par[4]<=5){return(999999)}
    else {return(with(data, sum((combi_expo(x,par[1],par[2],par[3],par[4]) - y)^2)))}
}

log_lik <- function(data,par){
    with(data,-sum(log(combi_expo(x,par[1],par[2],par[3],par[4]))))
}


expo <-function(x,a,b){
    return(a*exp(-x*b))
}

overlap_two_exp_fit <- function(X,Y,a1_guess=0.1,b1_guess=0.1,b2_guess=0.1,thresh_guess=8,thresh_down=5,thresh_up=15){
    xy=data.frame(y=Y[thresh_guess:length(Y)],x=X[thresh_guess:length(Y)])
    exp_nls=try(nls(y~(a*exp(-b*x)),data=xy,start=c(a=0.1,b=0.3),na.action=na.exclude),silent=TRUE) 
    if (class(exp_nls)!="try-error"){
        a2_guess=summary(exp_nls)$parameters[1]        
        b2_guess=summary(exp_nls)$parameters[2]        
        Y_substracted=Y-a2_guess*exp(-b2_guess*X)
    }

    xy=data.frame(y=Y_substracted[1:thresh_guess],x=X[1:thresh_guess])
    exp_nls=try(nls(y~(a*exp(-b*x)),data=xy,start=c(a=0.1,b=0.1),na.action=na.exclude),silent=TRUE) 
    if (class(exp_nls)!="try-error"){
        a1_guess=summary(exp_nls)$parameters[1]
        b1_guess=summary(exp_nls)$parameters[2]
    }

    xy=data.frame(y=Y,x=X)
    combi_nls=try(nls(y~overlap_expo(x,a1,b1,a2,b2),data=xy,start=c(a1=a1_guess,b1=b1_guess,a2=a2_guess,b2=b2_guess),algorithm="port",lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf),na.action=na.exclude,nls.control(maxiter = 10000, tol = 1e-04, minFactor=1/10024)),silent=TRUE)
    if (class(combi_nls)!="try-error"){
        test=try(summary(combi_nls))
        if (class(test)!="try-error"){
        
            a1=summary(combi_nls)$parameters[1]
            b1=summary(combi_nls)$parameters[2]
            a2=summary(combi_nls)$parameters[3] 
            b2=summary(combi_nls)$parameters[4] 
            comb_fit=overlap_expo(X,a1,b1,a2,b2)
            R2=1-sum(((Y-comb_fit)^2),na.rm=TRUE)/sum(((Y-mean(Y,na.rm=TRUE))^2),na.rm=TRUE)
            BIC=try(BIC(combi_nls),silent=TRUE)
            if (class(BIC)=="try-error"){BIC=NA}
            return(list(pars=c(a1,b1,a2,b2,NA),ana=c(R2,BIC),fit=comb_fit))
        }
    }
    return(list(pars=c(NA,NA,NA,NA,NA),ana=c(NA,NA),fit=X*NA))
}

two_exp_fit <- function(X,Y,y,a1_guess=0.1,b1_guess=0.1,b2_guess=0.1,thresh_guess=8,thresh_down=5,thresh_up=15){
    
    xy=data.frame(y=Y[1:thresh_guess],x=X[1:thresh_guess])
    exp_nls=try(nls(y~(a*exp(-b*x)),data=xy,start=c(a=0.1,b=0.1),na.action=na.exclude),silent=TRUE) 
    if (class(exp_nls)!="try-error"){
        a1_guess=summary(exp_nls)$parameters[1]
        b1_guess=summary(exp_nls)$parameters[2]
    }

    xy=data.frame(y=Y[thresh_guess:length(Y)],x=X[thresh_guess:length(Y)])
    exp_nls=try(nls(y~(a1_guess*exp((b-b1_guess)*thresh_guess)*exp(-b*x)),data=xy,start=c(b=0.3),na.action=na.exclude),silent=TRUE) 
    if (class(exp_nls)!="try-error"){
        b2_guess=summary(exp_nls)$parameters[1]        
        a2_guess=a1_guess*exp((b2_guess-b1_guess)*thresh_guess)
    }

    xy=data.frame(y=Y,x=X)
    #combi_opt=try(optim(par=c(a1_guess,b1_guess,b2_guess,thresh_guess),to_be_minimized,data=xy))
    #if (class(combi_opt)!="try-error"){
    #    a1=combi_opt$par[1]
    #    b1=combi_opt$par[2]
    #    b2=combi_opt$par[3]
    #    thresh=combi_opt$par[4]
    #    a2=a1*exp((b2-b1)*thresh)
    #    comb_fit=combi_expo(X,a1,b1,b2,thresh)
    #}
    combi_nls=try(nls(y~combi_expo(x,a1,b1,b2,thresh),data=xy,start=c(a1=a1_guess,b1=b1_guess,b2=b2_guess,thresh=thresh_guess),algorithm="port",lower=c(0,0,0,thresh_down),upper=c(Inf,Inf,Inf,thresh_up),na.action=na.exclude,nls.control(maxiter = 10000, tol = 1e-04, minFactor=1/10024, warnOnly=TRUE)),silent=TRUE)
    if (class(combi_nls)!="try-error"){
        a1=summary(combi_nls)$parameters[1]
        b1=summary(combi_nls)$parameters[2]
        b2=summary(combi_nls)$parameters[3] 
        thresh=summary(combi_nls)$parameters[4] 
        a2=a1*exp((b2-b1)*thresh)
        comb_fit=combi_expo(X,a1,b1,b2,thresh)
        R2=1-sum(((Y-comb_fit)^2),na.rm=TRUE)/sum(((Y-mean(Y,na.rm=TRUE))^2),na.rm=TRUE)
        BIC=try(BIC(combi_nls),silent=TRUE)
        if (class(BIC)=="try-error"){BIC=NA}
        return(list(pars=c(a1,b1,a2,b2,thresh),ana=c(R2,BIC),fit=comb_fit))
    }
    if (class(combi_nls)!="try-error"){return(list(pars=c(NA,NA,NA,NA,NA),ana=c(NA,NA),fit=X*NA))}
}

two_exp_fit_restricted <- function(X,Y,y,a1_guess=0.1,b1_guess=0.1,b2_guess=0.1,thresh_guess=8,thresh_down=5,thresh_up=15){
    
    xy=data.frame(y=Y[1:thresh_guess],x=X[1:thresh_guess])
    exp_nls=try(nls(y~(a*exp(-b*x)),data=xy,start=c(a=0.1,b=0.1),na.action=na.exclude),silent=TRUE) 
    if (class(exp_nls)!="try-error"){
        a1_guess=summary(exp_nls)$parameters[1]
        b1_guess=summary(exp_nls)$parameters[2]
    }

    xy=data.frame(y=Y[thresh_guess:length(Y)],x=X[thresh_guess:length(Y)])
    exp_nls=try(nls(y~(a1_guess*exp((b-b1_guess)*thresh_guess)*exp(-b*x)),data=xy,start=c(b=0.3),na.action=na.exclude),silent=TRUE) 
    if (class(exp_nls)!="try-error"){
        b2_guess=summary(exp_nls)$parameters[1]        
        a2_guess=a1_guess*exp((b2_guess-b1_guess)*thresh_guess)
    }

    xy=data.frame(y=Y,x=X)
    combi_opt=try(optim(par=c(a1_guess,b1_guess,b2_guess,thresh_guess),to_be_minimized,data=xy))
    if (class(combi_opt)!="try-error"){
        a1=combi_opt$par[1]
        b1=combi_opt$par[2]
        b2=combi_opt$par[3]
        thresh=combi_opt$par[4]
        a2=a1*exp((b2-b1)*thresh)
        comb_fit=combi_expo(X,a1,b1,b2,thresh)
    }
    # how to restrict b1 and b2 - 
    #combi_nls=try(nls(y~combi_expo(x,a1,b1,b2,thresh),data=xy,start=c(a1=a1_guess,b1=b1_guess,b2=b2_guess,thresh=thresh_guess),algorithm="port",lower=c(0,b2_guess-0.001,0,thresh_down),upper=c(Inf,Inf,b1_guess+0.001,thresh_up),na.action=na.exclude,nls.control(maxiter = 10000, tol = 1e-04, minFactor=1/10024, warnOnly=TRUE)),silent=TRUE)
    combi_nls=try(nls(y~combi_expo(x,a1,b1,b2,thresh),data=xy,start=c(a1=a1_guess,b1=b1_guess,b2=b2_guess,thresh=thresh_guess),algorithm="port",lower=c(0,0,0,thresh_down),upper=c(Inf,Inf,Inf,thresh_up),na.action=na.exclude,nls.control(maxiter = 10000, tol = 1e-04, minFactor=1/10024, warnOnly=TRUE)),silent=TRUE)
    if (class(combi_nls)!="try-error"){
        a1=summary(combi_nls)$parameters[1]
        b1=summary(combi_nls)$parameters[2]
        b2=summary(combi_nls)$parameters[3] 
        thresh=summary(combi_nls)$parameters[4] 
        a2=a1*exp((b2-b1)*thresh)
        comb_fit=combi_expo(X,a1,b1,b2,thresh)
        R2=1-sum(((Y-comb_fit)^2),na.rm=TRUE)/sum(((Y-mean(Y,na.rm=TRUE))^2),na.rm=TRUE)
        BIC=try(BIC(combi_nls),silent=TRUE)
        if (class(BIC)=="try-error"){BIC=NA}
    }
    if (b1>b2){return(list(pars=c(a1,b1,a2,b2,thresh),ana=c(R2,BIC),fit=comb_fit))}

    if (b1<b2){
        b=mean(c(b1,b2))
        combi_nls_restricted=try(nls(y~combi_expo(x,a1,b1,b2,thresh),data=xy,start=c(a1=a1_guess,b1=b+0.01,b2=b-0.01,thresh=thresh_guess),algorithm="port",lower=c(0,b,0,thresh_down),upper=c(Inf,Inf,b,thresh_up),na.action=na.exclude,nls.control(maxiter = 10000, tol = 1e-04, minFactor=1/10024, warnOnly=TRUE)),silent=TRUE)
        if (class(combi_nls_restricted)!="try-error"){
            summary=try(summary(combi_nls_restricted),silent=TRUE)
            if (class(summary)!="try-error"){
                a1=summary$parameters[1]
                b1=summary$parameters[2]
                b2=summary$parameters[3] 
                thresh=summary$parameters[4] 
                a2=a1*exp((b2-b1)*thresh)
                comb_fit=combi_expo(X,a1,b1,b2,thresh)
                R2=1-sum(((Y-comb_fit)^2),na.rm=TRUE)/sum(((Y-mean(Y,na.rm=TRUE))^2),na.rm=TRUE)
                BIC=BIC(combi_nls_restricted)
                return(list(pars=c(a1,b1,a2,b2,thresh),ana=c(R2,BIC),fit=comb_fit))
            }
        }
    }
    return(list(pars=c(NA,NA,NA,NA,NA),ana=c(NA,NA),fit=X*NA))
}


two_exp_fit_fixed_thresh <- function(X,Y,a1_guess=0.1,b1_guess=0.1,b2_guess=0.1,lower_limit=c(0,-Inf,0,-Inf,5),upper_limit=c(Inf,Inf,Inf,Inf,Inf)){
    xy=data.frame(y=Y,x=X)
    
    thresh=(4:17)*0.5+2
    trials=length(thresh)
    a1=array(NA,dim=c(trials))
    b1=array(NA,dim=c(trials))
    a2=array(NA,dim=c(trials))
    b2=array(NA,dim=c(trials))
    comb_fit=array(NA,dim=c(trials,length(X)))
    R2=array(NA,dim=c(trials))
    BIC=array(NA,dim=c(trials))

    for (i in 1:trials){
        xy=data.frame(y=Y[1:thresh[i]],x=X[1:thresh[i]])
        exp_nls=try(nls(y~(a*exp(-b*x)),data=xy,start=c(a=0.1,b=0.1),na.action=na.exclude),silent=TRUE) 
        if (class(exp_nls)!="try-error"){
            a1_guess=summary(exp_nls)$parameters[1]
            b1_guess=summary(exp_nls)$parameters[2]
        }

        xy=data.frame(y=Y[thresh[i]:length(Y)],x=X[thresh[i]:length(Y)])
        exp_nls=try(nls(y~(a1_guess*exp((b-b1_guess)*thresh[i])*exp(-b*x)),data=xy,start=c(b=0.3),na.action=na.exclude),silent=TRUE) 
        if (class(exp_nls)!="try-error"){
            b2_guess=summary(exp_nls)$parameters[1]        
            a2_guess=a1_guess*exp((b2_guess-b1_guess)*thresh[i])
        }


        xy=data.frame(y=Y,x=X)
        combi_nls=try(nls(y~combi_expo(x,a1,b1,b2,thresh[i]),data=xy,start=c(a1=a1_guess,b1=b1_guess,b2=b2_guess),na.action=na.exclude),silent=TRUE)

        plot(X,Y,xlab="days",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=TRUE,pch=20,col="grey",cex=0.5,log="y")
        if (class(combi_nls)!="try-error"){
            a1[i]=summary(combi_nls)$parameters[1]
            b1[i]=summary(combi_nls)$parameters[2]
            b2[i]=summary(combi_nls)$parameters[3] 
            a2[i]=a1[i]*exp((b2[i]-b1[i])*thresh[i])
            comb_fit[i,]=combi_expo(X,a1[i],b1[i],b2[i],thresh[i])
            R2[i]=1-sum(((Y-comb_fit[i,])^2),na.rm=TRUE)/sum(((Y-mean(Y,na.rm=TRUE))^2),na.rm=TRUE)
            BIC[i]=BIC(combi_nls)

            lines(X,comb_fit[i,],col="green")
            abline(v=thresh[i],col="lightblue") 
            legend("topright",legend=c(paste(round(BIC[i],01),"\n",round(R2[i],03),"\n",thresh[i])),bty="n")
        }   
    }
    if (length(which(!is.na(R2)))<3){        
        return(list(pars=c(NA,NA,NA,NA,NA),ana=c(NA,NA),fit=X*NA))
    }
    if (length(which(!is.na(R2)))>1){   
        best=which(R2==max(R2,na.rm=TRUE))
        return(list(pars=c(a1[best],b1[best],a2[best],b2[best],thresh[best]),ana=c(R2[best],BIC[best]),fit=comb_fit[best,]))
    }
}

general_extreme_values_fit <- function(X,Y,start_guess=c(xi=0.3,mu=1,beta=2),lower_limit=c(-Inf,-Inf,-Inf),upper_limit=c(Inf,Inf,Inf)){



    #really bad!!!!!!



    xy=data.frame(y=Y,x=X)
    # try fit
    gev_nls=try(nls(y~dgev(x=x,xi=xi,mu=mu,beta=beta),data=xy,algorithm="plinear",start=start_guess,lower=lower_limit,upper=upper_limit,na.action=na.exclude)) 
 
    print(gev_nls)
    print(X)

    x<<-X
    y<<-Y

    # if succes
    if (class(gev_nls)!="try-error"){
        xi=summary(gev_nls)$parameters[1]
        mu=summary(gev_nls)$parameters[2]
        beta=summary(gev_nls)$parameters[3]

        gevfit=dgev(x=X,xi=xi,mu=mu,beta=beta)
        R2=1-sum(((Y-gevfit)^2),na.rm=TRUE)/sum(((Y-mean(Y,na.rm=TRUE))^2),na.rm=TRUE)
        BIC=BIC(gev_nls)

        return(list(pars=c(xi,mu,beta),ana=c(R2,BIC),fit=gevfit))
    }
    # if fail



    else{
        return(list(pars=start_guess,ana=c(NA,NA),fit=dgev(x=X,xi=start_guess[1],mu=start_guess[2],beta=start_guess[3])))
    }
}


plot_fits_for_region <- function(period="1950-2014",trendID="91_5",dataset="_TMean",fit_style="2expo_thresh_5-15",reg=13,model="combi_expo",region_name="srex",ID_select=1:1319){

    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style,".nc",sep=""))
    fit_stuff_reg=var.get.nc(nc,"fit_stuff")
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_fit_",fit_style,".nc",sep=""))
    fit_stuff_individual=var.get.nc(nc,"fit_stuff")
    distr_stuff_individual=var.get.nc(nc,"distr_stuff")

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",period,"/",trendID,"_",region_name,"_",period,"_region_",reg,".pdf",sep=""))
    
    color=c(rgb(0.5,1,0.8),rgb(0.3,0.9,0.8),rgb(0.9,0.6,0.8),rgb(0.8,0.2,0.6),rgb(0.5,0.5,0.5,0.2))

    X=(0:30)+0.5
    for (sea in 1:5){
        for (state in 1:2){
            plot(NA,ylim=c(0,10),xlim=c(0,10),axes=FALSE,frame.plot=FALSE)
            text(5,5,paste(sea,state,reg))
            plot(NA,xlab="days",ylab="counts",ylim=c(0,150),xlim=c(0,30),axes=TRUE,frame.plot=TRUE)
            for (q in ID_select){
                points(distr_stuff_individual[sea,q,state,1,],distr_stuff_individual[sea,q,state,3,],pch=16,col=color[5],cex=1.5)
            }

            plot(NA,xlab="days",ylab="probability density",ylim=c(0.00001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE)
            for (q in ID_select){points(distr_stuff_individual[sea,q,state,1,],distr_stuff_individual[sea,q,state,2,],pch=16,col=color[5],cex=1.5)}
            for (q in ID_select){
                abline(v=fit_stuff_individual[sea,q,state,9],col=color[5],lty=2,lwd=0.3)
                lines(X,combi_expo(X,fit_stuff_individual[sea,q,state,5],fit_stuff_individual[sea,q,state,6],fit_stuff_individual[sea,q,state,8],fit_stuff_individual[sea,q,state,9]),col=color[1],lwd=0.3)
                lines(X,expo(X,fit_stuff_individual[sea,q,state,1],fit_stuff_individual[sea,q,state,2]),col=color[3],lwd=0.3,lty=2)
            }
            abline(v=fit_stuff_reg[sea,reg,state,9],col=color[5],lty=2,lwd=2)
            lines(X,combi_expo(X,fit_stuff_reg[sea,reg,state,5],fit_stuff_reg[sea,reg,state,6],fit_stuff_reg[sea,reg,state,8],fit_stuff_reg[sea,reg,state,9]),col=color[2],lwd=2,lty=1)
            lines(X,expo(X,fit_stuff_reg[sea,reg,state,1],fit_stuff_reg[sea,reg,state,2]),col=color[4],lty=1,lwd=2)
            

            plot(NA,xlab="days",ylab="probability density",ylim=c(0.00001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,log="y")
            for (q in ID_select){points(distr_stuff_individual[sea,q,state,1,],distr_stuff_individual[sea,q,state,2,],pch=16,col=color[5],cex=1.5)}
            for (q in ID_select){
                abline(v=fit_stuff_individual[sea,q,state,9],col=color[5],lty=2,lwd=0.3)
                lines(X,combi_expo(X,fit_stuff_individual[sea,q,state,5],fit_stuff_individual[sea,q,state,6],fit_stuff_individual[sea,q,state,8],fit_stuff_individual[sea,q,state,9]),col=color[1],lwd=0.3)
                lines(X,expo(X,fit_stuff_individual[sea,q,state,1],fit_stuff_individual[sea,q,state,2]),col=color[3],lwd=0.3,lty=2)
            }
            abline(v=fit_stuff_reg[sea,reg,state,9],col=color[5],lty=2,lwd=2)
            lines(X,combi_expo(X,fit_stuff_reg[sea,reg,state,5],fit_stuff_reg[sea,reg,state,6],fit_stuff_reg[sea,reg,state,8],fit_stuff_reg[sea,reg,state,9]),col=color[2],lwd=2,lty=1)
            lines(X,expo(X,fit_stuff_reg[sea,reg,state,1],fit_stuff_reg[sea,reg,state,2]),col=color[4],lty=1,lwd=2)

        }
    }
    
    
}

