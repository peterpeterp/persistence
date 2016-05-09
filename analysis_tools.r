library(dgof)

cdf_onDist <- function(dist){
    dist[is.na(dist)]=0
    cdf<-dist*NA
    for (i in 1:length(dist)){cdf[i]<-sum(dist[1:i])}
    cdf<-cdf/max(cdf)
    return(cdf)
}

cdf_onData <- function(dist,breaks=1:max(dist,na.rm=TRUE)){
    dist<-dist[!is.na(dist)]
    cdf=array(NA,length(breaks))
    for (i in 1:length(breaks)){
        cdf[i]=length(which(dist<=breaks[i]))
    }
    cdf=cdf/cdf[length(cdf)]
}

ks_distribution <- function(t){
    sum<-0
    for (i in 1:100){sum<-sum+(-1)^(i-1)*exp(-2*i^2*t)}
    return(1-2*sum)
}



ks_test_manual <- function(X_toFit,cdf1,cdf2,n=length(X_toFit)){
    diff_cdf<-abs(cdf1-cdf2)
    D_val<-max(diff_cdf)
    D_pos<-X_toFit[which.max(diff_cdf)]
    
    ks<-0
    for (alpha in c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)){
        #if (D_val<sqrt(-0.5*log(alpha/2)*2/n)){ks<-alpha} # statistics book from sonja for 2 distrs
        if (D_val<sqrt(-0.5*log(alpha/2)/n)){ks<-alpha} # table for one distr
    }
    return(list(D_val=D_val,D_pos=D_pos,ks=ks))
}

quantile_pete <- function(dist,taus,na.rm=TRUE,plot=FALSE){
    # calculates quantiles from empirical cumulative distribution function -> gives out decimal numbers
    if (na.rm==TRUE){dist=dist[which(!is.na(dist))]}

    # get cdf
    cdf<-cdf_onData(dist)

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


fit_plot_empty <- function(xUntil=60){
    # creates empty plots with x-axis and y-axis
    par(mar=c(3, 3, 3, 3) + 0.1)
    # x
    plot(NA,xlab="",ylim=c(0.00001,0.25),xlim=c(0,xUntil),ylab="",axes=FALSE,frame.plot=FALSE,cex=0.5)
    at_=axis(1,labels=FALSE,col="black")
    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
    axis(1,at=at_)

    # x
    plot(NA,xlab="",ylim=c(0.00001,0.25),xlim=c(0,xUntil),ylab="",axes=FALSE,frame.plot=FALSE,cex=0.5)
    at_=axis(3,labels=FALSE,col="black")
    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
    axis(3,at=at_)

    # y
    plot(NA,xlab="",ylim=c(0.00001,0.25),xlim=c(0,xUntil),ylab="",axes=FALSE,frame.plot=FALSE,cex=0.5)
    at_=axis(2,labels=FALSE,col="black")
    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
    axis(2,at=at_)

    # y
    plot(NA,xlab="",ylim=c(0.00001,0.25),xlim=c(0,xUntil),ylab="",axes=FALSE,frame.plot=FALSE,cex=0.5)
    at_=axis(4,labels=FALSE,col="black")
    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
    axis(4,at=at_)

    # y
    plot(NA,xlab="",ylim=c(0.00001,0.25),xlim=c(0,xUntil),ylab="",axes=FALSE,frame.plot=FALSE,log="y",cex=0.5)
    at_=axis(2,labels=FALSE,col="black")
    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
    axis(2,at=at_)

    # y
    plot(NA,xlab="",ylim=c(0.00001,0.25),xlim=c(0,xUntil),ylab="",axes=FALSE,frame.plot=FALSE,log="y",cex=0.5)
    at_=axis(4,labels=FALSE,col="black")
    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
    axis(4,at=at_)

    # ecdf y 
    plot(NA,xlab="",ylim=c(0,1),xlim=c(0,xUntil),ylab="",axes=FALSE,frame.plot=FALSE,cex=0.5)
    at_=axis(2,labels=FALSE,col="black")
    if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
    axis(2,at=at_)

    # ecdf y
    plot(NA,xlab="",ylim=c(0,1),xlim=c(0,xUntil),ylab="",axes=FALSE,frame.plot=FALSE,cex=0.5)
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

fit_plot_combi <- function(distrs,fitstuff,fit_style,sea,q,state,xUntil=60){
    color=c("blue","red",rgb(0.5,0.1,0.5),rgb(0.1,0.7,0.1),rgb(0,0,0))

    par(mar=c(3, 3, 3, 3) + 0.1)  

    # second plot page
    nonull=which(distrs[3,]>0)              
    plot(distrs[1,nonull],distrs[3,nonull],xlab="days",xlim=c(0,xUntil),ylab="counts",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],cex=0.5)
    axis(2)
    legend("topright",legend=c(paste("     ID=",q,sep=""),paste("dBIC=",round(fitstuff[24],01),sep="")),bty="n")

    # second plot page
    par(mar=c(3, 3, 3, 3) + 0.1)                
    plot(distrs[1,nonull],distrs[2,nonull],xlab="days",ylim=c(0.00001,0.25),xlim=c(0,xUntil),ylab="",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],cex=0.5)
    if (!is.na(fitstuff[15])){
        abline(v=fitstuff[15],col="grey")
        #text(thresh,0.22,label=round(thresh,02),col="grey")
    }
    lines(distrs[4,],col=color[3],lty=2)
    lines(distrs[8,],col=color[4],lty=1)               
 
    text(40,0.18+0.06,paste("BIC=",round(fitstuff[9],01)),pos=1,col=color[3])                 
    text(40,0.155+0.06,paste("BIC=",round(fitstuff[22],01)),pos=1,col=color[4])      

    text(40,0.11+0.06,paste("KS=",round(fitstuff[8],03)),pos=1,col=color[3])                 
    text(40,0.085+0.06,paste("KS=",round(fitstuff[21],03)),pos=1,col=color[4])            

    # third version
    plot(distrs[1,nonull],distrs[2,nonull],xlab="days",ylim=c(0.00001,0.25),xlim=c(0,xUntil),ylab="",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],log="y",cex=0.5)
    if (!is.na(fitstuff[15])){
        abline(v=fitstuff[15],col="grey")
        text(fitstuff[15],0.00002,label=round(fitstuff[15],02),col=rgb(0,0,0))
    }    
    lines(distrs[4,],col=color[3],lty=2)
    lines(distrs[8,],col=color[4],lty=1)  

    text(50,0.22,paste("b=",round(fitstuff[2],02)),pos=1,col=color[3])                 
    text(50,0.05,paste("b1=",round(fitstuff[12],02)),pos=1,col=color[4])                 
    text(50,0.02,paste("b2=",round(fitstuff[14],02)),pos=1,col=color[4])                

    plot(NA,xlim=c(0,xUntil),ylim=c(0,1),xlab="days",ylab="",axes=FALSE,frame.plot=TRUE)
    points(distrs[5,],col=color[state],pch=20,cex=0.4)

    x<-fitstuff[7] ; lines(c(x,x),c(distrs[6,x],distrs[5,x]),col=color[3],lwd=2)
    x<-fitstuff[20] ; lines(c(x,x),c(distrs[10,x],distrs[9,x]),col=color[4],lwd=2)

    lines(distrs[6,],col=color[3])
    lines(distrs[10,],col=color[4])
    
    text(40,0.4,paste("D=",round(fitstuff[6],03)),col=color[3],pos=1)                 
    text(40,0.3,paste("D=",round(fitstuff[19],03)),col=color[4],pos=1)               
}


exponential_fit <- function(X,Y,start_guess=c(a=0.1,b=0.1),lower_limit=c(-Inf,-Inf),upper_limit=c(Inf,Inf),xStart=1,plot_cdf=FALSE){

    # only data points until longest period 
    longestPeriod<-max(X[Y>0])   
    Y_toFit<-Y[xStart:longestPeriod]
    X_toFit<-X[xStart:longestPeriod]

    # perform fit
    xy<-data.frame(y=Y_toFit,x=X_toFit)
    nls_fit<-try(nls(y~(a*exp(-b*x)),data=xy,algorithm="port",start=start_guess,lower=lower_limit,upper=upper_limit,na.action=na.exclude),silent=TRUE) 
    if (class(nls_fit)!="try-error"){
        summ_nls<-try(summary(nls_fit))
        if (class(summ_nls)!="try-error"){
            a<-summ_nls$parameters[1]
            b<-summ_nls$parameters[2]

            # goodness analysis
            exp_fit<-a*exp(-X_toFit*b)

            chi2<-chisq.test(Y_toFit,p=exp_fit,rescale=TRUE)$p.value
            R2<-1-sum(((Y_toFit-exp_fit)^2),na.rm=TRUE)/sum(((Y_toFit-mean(Y_toFit,na.rm=TRUE))^2),na.rm=TRUE)


            # ks test
            cdf_Data<-cdf_onDist(Y_toFit)
            cdf_Fit<-cdf_onDist(exp_fit)
            tmpks<-ks_test_manual(X_toFit,cdf_Data,cdf_Fit)

            # BIC
            BIC<-try(BIC(nls_fit),silent=TRUE)
            if (class(BIC)=="try-error"){BIC=NA}

            exp_fit<-a*exp(-X*b)
            return(list(pars=c(a,b),ana=c(chi2,R2,tmpks$D_val,tmpks$D_pos,tmpks$ks,BIC),fit=exp_fit,cdf_Data=cdf_Data,cdf_Fit=cdf_Fit))
        }
        if (class(summ_nls)=="try-error"){return(list(pars=c(NA,NA),ana=c(NA,NA,NA,NA,NA,NA),fit=X_toFit*NA,cdf_Data=X_toFit*NA,cdf_Fit=X_toFit*NA))}
    }
    if (class(nls_fit)=="try-error"){return(list(pars=c(NA,NA),ana=c(NA,NA,NA,NA,NA,NA),fit=X_toFit*NA,cdf_Data=X_toFit*NA,cdf_Fit=X_toFit*NA))}
}

combi_expo <-function(x,a1,b1,b2,thresh){
    y=x*NA
    y[x<=thresh]=a1*exp(-x[x<=thresh]*b1)
    y[x>thresh]=a1*exp((b2-b1)*thresh)*exp(-(x[x>thresh])*b2)
    return(y)
}

two_exp_fit <- function(X,Y,y,a1_guess=0.1,b1_guess=0.1,b2_guess=0.1,thresh_guess=8,thresh_down=5,thresh_up=15,xStart=1,xStop=100,plot_cdf=FALSE){
    # only data points until longest period 
    longestPeriod<-max(X[Y>0])   
    Y_toFit<-Y[xStart:longestPeriod]
    X_toFit<-X[xStart:longestPeriod]

    # prepare guesses for the parameters
    xy<-data.frame(y=Y[which(X>=xStart & X<=thresh_guess)],x=X[which(X>=xStart & X<=thresh_guess)])
    exp_nls<-try(nls(y~(a*exp(-b*x)),data=xy,start=c(a=0.1,b=0.1),na.action=na.exclude),silent=TRUE) 
    if (class(exp_nls)!="try-error"){
        a1_guess<-summary(exp_nls)$parameters[1]
        b1_guess<-summary(exp_nls)$parameters[2]
    }

    xy<-data.frame(y=Y[which(X<=xStop & X>=thresh_guess)],x=X[which(X<=xStop & X>=thresh_guess)])
    exp_nls<-try(nls(y~(a1_guess*exp((b-b1_guess)*thresh_guess)*exp(-b*x)),data=xy,start=c(b=0.3),na.action=na.exclude),silent=TRUE) 
    if (class(exp_nls)!="try-error"){
        b2_guess<-summary(exp_nls)$parameters[1]        
        a2_guess<-a1_guess*exp((b2_guess-b1_guess)*thresh_guess)
    }

    # perform actual fit
    xy<-data.frame(y=Y_toFit,x=X_toFit)
    nls_fit<-try(nls(y~combi_expo(x,a1,b1,b2,thresh),data=xy,start=c(a1=a1_guess,b1=b1_guess,b2=b2_guess,thresh=thresh_guess),algorithm="port",lower=c(0,0,0,thresh_down),upper=c(Inf,Inf,Inf,thresh_up),na.action=na.exclude,nls.control(maxiter = 10000, tol = 1e-04, minFactor=1/10024, warnOnly=TRUE)),silent=TRUE)
    if (class(nls_fit)!="try-error"){
        summ_nls<-try(summary(nls_fit))
        if (class(summ_nls)!="try-error"){
            # fit parameters
            a1<-summ_nls$parameters[1]
            b1<-summ_nls$parameters[2]
            b2<-summ_nls$parameters[3] 
            thresh<-summ_nls$parameters[4] 
            a2<-a1*exp((b2-b1)*thresh)

            # goodness analysis
            comb_fit<-combi_expo(X_toFit,a1,b1,b2,thresh)

            chi2<-chisq.test(Y_toFit,p=comb_fit,rescale=TRUE)$p.value
            R2<-1-sum(((Y_toFit-comb_fit)^2),na.rm=TRUE)/sum(((Y_toFit-mean(Y_toFit,na.rm=TRUE))^2),na.rm=TRUE)

            # ks test
            cdf_Data<-cdf_onDist(Y_toFit)
            cdf_Fit<-cdf_onDist(comb_fit)
            tmpks<-ks_test_manual(X_toFit,cdf_Data,cdf_Fit)

            # BIC
            BIC<-try(BIC(nls_fit),silent=TRUE)
            if (class(BIC)=="try-error"){BIC=NA}

            # store the fit for complete X in order to plot it later on
            comb_fit<-combi_expo(X,a1,b1,b2,thresh)
            return(list(pars=c(a1,b1,a2,b2,thresh),ana=c(chi2,R2,tmpks$D_val,tmpks$D_pos,tmpks$ks,BIC),fit=comb_fit,cdf_Data=cdf_Data,cdf_Fit=cdf_Fit))
        }
        if (class(summ_nls)=="try-error"){return(list(pars=c(NA,NA,NA,NA,NA),ana=c(NA,NA,NA,NA,NA,NA),fit=X_toFit*NA,cdf_Data=X_toFit*NA,cdf_Fit=X_toFit*NA))}
    }
    if (class(combi_nls)=="try-error"){return(list(pars=c(NA,NA,NA,NA,NA),ana=c(NA,NA,NA,NA,NA,NA),fit=X_toFit*NA,cdf_Data=X_toFit*NA,cdf_Fit=X_toFit*NA))}
}

