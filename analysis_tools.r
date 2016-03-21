library(dgof)

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
        quant_zwi<-try(summary(rq(dither(y,value=noise_level)~x,taus[i]),cov=FALSE,se="boot"))#,silent=TRUE
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

exponential_fit <- function(X,Y,y,start_guess=c(a=0.1,b=0.1),lower_limit=c(-Inf,-Inf),upper_limit=c(Inf,Inf),xStart=1,xStop=100){
    Y<-Y[xStart:xStop]
    X<-X[xStart:xStop]
    nona<-which(Y>0)
    Y<-Y[nona]
    X<-X[nona]
    xy=data.frame(y=Y,x=X)
    # try fit
    exp_nls=try(nls(y~(a*exp(-b*x)),data=xy,algorithm="port",start=start_guess,lower=lower_limit,upper=upper_limit,na.action=na.exclude),silent=TRUE) 
    
    # if succes
    if (class(exp_nls)!="try-error"){
        a<-summary(exp_nls)$parameters[1]
        b<-summary(exp_nls)$parameters[2]

        expfit<-a*exp(-X*b)
        R2<-1-sum(((Y-expfit)^2),na.rm=TRUE)/sum(((Y-mean(Y,na.rm=TRUE))^2),na.rm=TRUE)
        BIC<-BIC(exp_nls)
        chi2<-chisq.test(Y[which(!is.na(Y))],p=expfit[which(!is.na(Y))],rescale=TRUE)$p.value
        ks<-ks.test(Y,expfit)$p.value
        return(list(pars=c(a,b),ana=c(chi2,R2,ks,BIC),fit=expfit))
    }
    # if fail
    else{
        return(list(pars=c(NA,NA),ana=c(NA,NA,NA,NA),fit=X*NA))
    }
}

combi_expo <-function(x,a1,b1,b2,thresh){
    y=x*NA
    y[x<=thresh]=a1*exp(-x[x<=thresh]*b1)
    y[x>thresh]=a1*exp((b2-b1)*thresh)*exp(-(x[x>thresh])*b2)
    return(y)
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
            chi2<-chisq.test(Y[which(!is.na(Y))],p=comb_fit[which(!is.na(Y))],rescale=TRUE)$p.value
            ks<-ks.test(Y,comb_fit)$p.value
            if (class(BIC)=="try-error"){BIC=NA}
            return(list(pars=c(a1,b1,a2,b2,thresh),ana=c(chi2,R2,ks,BIC),fit=comb_fit))
        }
        if (class(summ_nls)=="try-error"){return(list(pars=c(NA,NA,NA,NA,NA),ana=c(NA,NA,NA,NA),fit=X*NA))}
    }
    if (class(combi_nls)=="try-error"){return(list(pars=c(NA,NA,NA,NA,NA),ana=c(NA,NA,NA,NA),fit=X*NA))}
}

