
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

quantile_analysis <- function(x,y,taus,noise_level=0){
    # x=dur_mid, y=dur, taus=percentiles
    #cat(length(unique(y)))
    quantiles<-quantile_pete(y,taus=taus,na.rm=TRUE)
    slopes<-taus*NA
    slope_sigs<-taus*NA

    # try quantile regression for each tau individually 
    for (i in 1:length(taus)){
        # noise is added via "dither()" to x and y to avoid crash
        print(proc.time()) 
        quant_zwi1<<-try(summary(rq(y~x,taus[i]),cov=FALSE,se="boot"))#,silent=TRUE
        print(proc.time())        
        quant_zwi2<<-try(summary(rq(y~x,taus[i]),cov=FALSE,se="nid"))#,silent=TRUE
        print(proc.time())
        quant_zwi3<<-try(summary(rq(y~x,taus[i]),cov=FALSE,se="iid"))#,silent=TRUE
        print(proc.time())
        quant_zwi4<<-try(summary(rq(y~x,taus[i]),cov=FALSE,se="ker"))#,silent=TRUE
        print(proc.time())
        quant_zwi5<<-try(summary(rq(dither(y,value=noise_level)~x,taus[i]),cov=FALSE,se="boot"))#,silent=TRUE
        print(proc.time())
        quant_zwi6<<-try(summary(rq(dither(y,value=noise_level)~x,taus[i]),cov=FALSE,se="nid"))#,silent=TRUE
        print(proc.time())
        quant_zwi7<<-try(summary(rq(dither(y,value=noise_level)~x,taus[i]),cov=FALSE,se="iid"))#,silent=TRUE
        print(proc.time())
        quant_zwi8<<-try(summary(rq(dither(y,value=noise_level)~x,taus[i]),cov=FALSE,se="ker"))#,silent=TRUE
        print(proc.time())
        quant_zwi<<-try(summary(rq(dither(y,value=noise_level)~x,taus[i]),cov=FALSE))#,silent=TRUE
        print(proc.time())
        adsasd
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

fit_plot_comparison <- function(distr,fits,sea,q,state,wilcox){
    color=c("green","orange",rgb(0.5,0.1,0.5),rgb(0.1,0.7,0.1),rgb(0,0,0))

    par(mar=c(3, 3, 3, 3) + 0.1)  

    # first plot page
    plot(NA,xlab="days",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],cex=0.5)
    text(50,0.18+0.06,paste("wilcox=",round(wilcox,01),"\n ",season_names[sea],q,state_names[state]),pos=1,col=color[3])                 
    for (i in 1:2){
        nonull<-which(distr[i,3,]>0)
        points(distr[i,1,nonull],distr[i,2,nonull],pch=20,col=color[i],cex=0.5)
        lines(distr[i,1,nonull],combi_expo(distr[i,1,nonull],fits[i,5],fits[i,6],fits[i,8],fits[i,9]),col=color[i],lty=1)

    }  

    # second plot page
    plot(NA,xlab="days",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],cex=0.5,log="y")
    text(50,0.22,paste("wilcox=",round(wilcox,01),"\n ",season_names[sea],q,state_names[state]),pos=1,col=color[3])                 
    for (i in 1:2){
        nonull<-which(distr[i,3,]>0)
        points(distr[i,1,nonull],distr[i,2,nonull],pch=20,col=color[i],cex=0.5)
        lines(distr[i,1,nonull],combi_expo(distr[i,1,nonull],fits[i,5],fits[i,6],fits[i,8],fits[i,9]),col=color[i],lty=1)
    }              
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

exponential_fit <- function(X,Y,start_guess=c(a=0.1,b=0.1),lower_limit=c(-Inf,-Inf),upper_limit=c(Inf,Inf),xStart=1,xStop=100){
    xy=data.frame(y=Y[xStart:xStop],x=X[xStart:xStop])
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

two_exp_fit <- function(X,Y,y,a1_guess=0.1,b1_guess=0.1,b2_guess=0.1,thresh_guess=8,thresh_down=5,thresh_up=15,xStart=1,xStop=100){
    xy<-data.frame(y=Y[xStart:thresh_guess],x=X[xStart:thresh_guess])
    exp_nls<-try(nls(y~(a*exp(-b*x)),data=xy,start=c(a=0.1,b=0.1),na.action=na.exclude),silent=TRUE) 
    if (class(exp_nls)!="try-error"){
        a1_guess<-summary(exp_nls)$parameters[1]
        b1_guess<-summary(exp_nls)$parameters[2]
    }

    xy<-data.frame(y=Y[thresh_guess:xStop],x=X[thresh_guess:xStop])
    exp_nls<-try(nls(y~(a1_guess*exp((b-b1_guess)*thresh_guess)*exp(-b*x)),data=xy,start=c(b=0.3),na.action=na.exclude),silent=TRUE) 
    if (class(exp_nls)!="try-error"){
        b2_guess<-summary(exp_nls)$parameters[1]        
        a2_guess<-a1_guess*exp((b2_guess-b1_guess)*thresh_guess)
    }
    xy<-data.frame(y=Y[xStart:xStop],x=X[xStart:xStop])
    combi_nls<-try(nls(y~combi_expo(x,a1,b1,b2,thresh),data=xy,start=c(a1=a1_guess,b1=b1_guess,b2=b2_guess,thresh=thresh_guess),algorithm="port",lower=c(0,0,0,thresh_down),upper=c(Inf,Inf,Inf,thresh_up),na.action=na.exclude,nls.control(maxiter = 10000, tol = 1e-04, minFactor=1/10024, warnOnly=TRUE)),silent=TRUE)
    if (class(combi_nls)!="try-error"){
        summ_nls<-try(summary(combi_nls))
        if (class(summ_nls)!="try-error"){
            a1<-summ_nls$parameters[1]
            b1<-summ_nls$parameters[2]
            b2<-summ_nls$parameters[3] 
            thresh<-summ_nls$parameters[4] 
            a2<-a1*exp((b2-b1)*thresh)
            comb_fit<-combi_expo(X,a1,b1,b2,thresh)
            R2<-1-sum(((Y-comb_fit)^2),na.rm=TRUE)/sum(((Y-mean(Y,na.rm=TRUE))^2),na.rm=TRUE)
            BIC<-try(BIC(combi_nls),silent=TRUE)
            if (class(BIC)=="try-error"){BIC=NA}
            return(list(pars=c(a1,b1,a2,b2,thresh),ana=c(R2,BIC),fit=comb_fit))
        }
        if (class(summ_nls)=="try-error"){return(list(pars=c(NA,NA,NA,NA,NA),ana=c(NA,NA),fit=X*NA))}

    }
    if (class(combi_nls)=="try-error"){return(list(pars=c(NA,NA,NA,NA,NA),ana=c(NA,NA),fit=X*NA))}
}

two_exp_fit_restricted <- function(X,Y,y,a1_guess=0.1,b1_guess=0.1,b2_guess=0.1,thresh_guess=8,thresh_down=4,thresh_up=15){
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




