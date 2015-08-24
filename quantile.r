source("load.r")
library(quantreg)

if (1==2){
	dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
	data=array(c(dat$time,as.vector(dat$tas[488,,])),dim=c(length(dat$time),2))
	write.table(data,"../data/488.txt")
}

if (1==2){
	nc=open.ncdf(paste("../data/",91,"_",5,"/",91,"_",5,"_duration_2s_","summer",".nc",sep=""))
	dur=get.var.ncdf(nc,"dur")
	dur_mid=get.var.ncdf(nc,"dur_mid")
	data=array(c(as.vector(dur_mid[488,,]),as.vector(dur[488,,])),dim=c(length(as.vector(as.vector(dur_mid[488,,]))),2))
	write.table(data,"../data/dur_488.txt")
}


if (1==2){
	dur=read.table("../data/sonstiges/dur_488.txt")
	size=dim(dur)[1]
	size=length(which(!is.na(dur[1:size,1])))

    pdf(file="../plots/quantile.pdf")
    plot(dur[1:size,1],dur[1:size,2])
    perc=c(0.5,0.75,0.95,0.98)
    #perc=c(0.5,0.95)
    qu=array(NA,dim=c(64,length(perc)))

    for (year in 1950:2013){
    	inside=which(dur[1:size,1]>year & dur[1:size,1]<(year+1))
    	x=dur[inside,1]
    	y=dur[inside,2]
    	qu[(year-1949),]=quantile(y,probs=perc,na.rm=TRUE)

    	print(sum(y))
    	print(summary(quantile(y,probs=perc,na.rm=TRUE)))
    	print((quantile(y,probs=perc,na.rm=TRUE)))

    	print(quantile(y,probs=0.5))
    	print(y)
    	print(length(which(y>quantile(y,probs=0.5))))
    	print(length(which(y>quantile(y,probs=0.75))))


    	
    }
    #fhf
    tYear=seq(1950.5,2013.5,1)

    for (p in 1:length(perc)){
	    requ=rq(dur[1:size,2]~dur[1:size,1],perc[p])
		quReg=summary(requ)$coefficients
		print("--------------RQ---------------")

		print(summary.rq(requ))
		print("--------------LR---------------")
		color=rgb(perc[p], ((length(perc)-p)/length(perc)), (1-perc[p]))
		lines(tYear,qu[,p],col=color)
		abline(requ,col=color)
		print(qu[,p])
		abline(lm(qu[,p]~tYear),col=color,lty=2)
		print(summary(lm(qu[,p]~tYear)))
	}


}

if (1==2){
	dd=read.table("../data/sonstiges/dur_488.txt")
	size=dim(dd)[1]
	size=length(which(!is.na(dd[1:size,1])))

    pdf(file="../plots/quantile.pdf")
    plot(dd[1:size,1],dd[1:size,2])

    for (perc in c(0.25,0.5,0.9,0.95,0.98,0.99)){
    	y=dd[1:1100,2]
    	x=dd[1:1100,1]
		test=summary(rq(y~x,0.95))$coefficients
		print(test[c(1,7,2,8)])

		y=c(y,array(24,10))
		x=c(x,array(tail(x,n=1),10))

		test=summary(rq(y~x,0.95))$coefficients
		print(test[c(1,7,2,8)])
		adsa
		t=1:11
		x=c()
		y=c()
		z=c()
		for (i in t){
			qu=quantile(dd[((i-1)*100+1):((i)*100),2],probs=c(perc))
			y[i]=qu
			x[i]=mean(dd[((i-1)*100+1):((i)*100),1],na.rm=TRUE)
			z[i]=dd[(i*100),1]
		}
		print(y)
		print(x)
		lm.r=lm(y~x)


	    lines(x,y,col="red")
	    #abline(lm.r)
	    abline(test,col="green")
	}
    for (i in t){
    	abline(v=z[i])
    }
}

if (1==2){
	nc=open.ncdf("../data/91_5/91_5_duration_2s_analysis_summer.nc")
	dur=get.var.ncdf(nc,"dur_ana_full")
	print(dur[487,2,,])
	print(dur[488,2,,])
	print(dur[489,2,,])
}

if (1==1){
	perc=c(0.25,0.5,0.75,0.9)
	x=seq(0,10000,1)
	y=x*NA
	xplot=x*NA
	samp=seq(-10,+10,length=1000)
	qu=array(NA,dim=c(10,length(perc)))
	xqu=array(NA,dim=c(10))
	for (t in 1:10){
		x[((t-1)*1000+1):(t*1000)]=array((t*1000),1000)
		zwi=dnorm(samp,mean=0,sd=(4+t))+0.001*t

		y[((t-1)*1000+1):(t*1000)]=zwi
		qu[t,]=quantile(zwi,perc)
		xqu[t]=t*1000
	}
	pdf(file="../plots/artificial_quantile.pdf")

	plot(x,y)

	for (p in 1:length(perc)){
		#color=rgb(perc[p], ((length(perc)-p)/length(perc)), (1-perc[p]))
		color=rgb(0.8,0.8,1)
		lines(xqu,as.vector(qu[,p]),col=color,lty=3)
		abline(rq(y~x,perc[p]),col=color)
		abline(lm(as.vector(qu[,p])~xqu),col=color,lty=2)
	}
	graphics.off()
}

# ---------------------- rq() function ------------------------------------------------------
function (formula, tau = 0.5, data, subset, weights, na.action, 
    method = "br", model = TRUE, contrasts = NULL, ...) 
{                                                                                                                                                                                   
    call <- match.call()                                                                                                                                                            
    mf <- match.call(expand.dots = FALSE)                                                                                                                                           
    m <- match(c("formula", "data", "subset", "weights", "na.action"),                                                                                                              
        names(mf), 0)                                                                                                                                                               
    mf <- mf[c(1, m)]                                                                                                                                                               
    mf$drop.unused.levels <- TRUE                                                                                                                                                   
    mf[[1]] <- as.name("model.frame")                                                                                                                                               
    mf <- eval.parent(mf)                                                                                                                                                           
    if (method == "model.frame")                                                                                                                                                    
        return(mf)                                                                                                                                                                  
    mt <- attr(mf, "terms")                                                                                                                                                         
    weights <- as.vector(model.weights(mf))                                                                                                                                         
    Y <- model.response(mf)                                                                                                                                                         
    X <- model.matrix(mt, mf, contrasts)                                                                                                                                            
    eps <- .Machine$double.eps^(2/3)                                                                                                                                                
    Rho <- function(u, tau) u * (tau - (u < 0))                                                                                                                                     
    if (length(tau) > 1) {                                                                                                                                                          
        if (any(tau < 0) || any(tau > 1))                                                                                                                                           
            stop("invalid tau:  taus should be >= 0 and <= 1")                                                                                                                      
        if (any(tau == 0))                                                                                                                                                          
            tau[tau == 0] <- eps
        if (any(tau == 1)) 
            tau[tau == 1] <- 1 - eps
        coef <- matrix(0, ncol(X), length(tau))
        rho <- rep(0, length(tau))
        fitted <- resid <- matrix(0, nrow(X), length(tau))
        for (i in 1:length(tau)) {
            z <- {
                if (length(weights)) 
                  rq.wfit(X, Y, tau = tau[i], weights, method, 
                    ...)
                else rq.fit(X, Y, tau = tau[i], method, ...)
            }
            coef[, i] <- z$coefficients
            resid[, i] <- z$residuals
            rho[i] <- sum(Rho(z$residuals, tau[i]))
            fitted[, i] <- Y - z$residuals
        }
        taulabs <- paste("tau=", format(round(tau, 3)))
        dimnames(coef) <- list(dimnames(X)[[2]], taulabs)
        dimnames(resid) <- list(dimnames(X)[[1]], taulabs)
        fit <- z
        fit$coefficients <- coef
        fit$residuals <- resid
        fit$fitted.values <- fitted
        if (method == "lasso") 
            class(fit) <- c("lassorqs", "rqs")
        else if (method == "scad") 
            class(fit) <- c("scadrqs", "rqs")
        else class(fit) <- "rqs"
    }
    else {
        process <- (tau < 0 || tau > 1)
        if (tau == 0) 
            tau <- eps
        if (tau == 1) 
            tau <- 1 - eps
        fit <- {
            if (length(weights)) 
                rq.wfit(X, Y, tau = tau, weights, method, ...)
            else rq.fit(X, Y, tau = tau, method, ...)
        }
        if (process) 
            rho <- list(x = fit$sol[1, ], y = fit$sol[3, ])
        else {
            dimnames(fit$residuals) <- list(dimnames(X)[[1]], 
                NULL)
            rho <- sum(Rho(fit$residuals, tau))
        }
        if (method == "lasso") 
            class(fit) <- c("lassorq", "rq")
        else if (method == "scad") 
            class(fit) <- c("scadrq", "rq")
        else class(fit) <- ifelse(process, "rq.process", "rq")
    }
    fit$na.action <- attr(mf, "na.action")
    fit$formula <- formula
    fit$terms <- mt
    fit$xlevels <- .getXlevels(mt, mf)
    fit$call <- call
    fit$tau <- tau
    fit$weights <- weights
    fit$residuals <- drop(fit$residuals)
    fit$rho <- rho
    fit$method <- method
    fit$fitted.values <- drop(fit$fitted.values)
    attr(fit, "na.message") <- attr(m, "na.message")
    if (model) 
        fit$model <- mf
    fit
}

# ------------------------------ rq.fit() ------------------------------------
Error in as.matrix(x) : 
  Fehler bei der Auswertung des Argumentes 'x' bei der Methodenauswahl
für Funktion 'as.matrix': Fehler: Argument "x" fehlt (ohne Standardwert)
> rq.fit
function (x, y, tau = 0.5, method = "br", ...) 
{
    fit <- switch(method, fn = rq.fit.fnb(x, y, tau = tau, ...), 
        fnb = rq.fit.fnb(x, y, tau = tau, ...), fnc = rq.fit.fnc(x, 
            y, tau = tau, ...), pfn = rq.fit.pfn(x, y, tau = tau, 
            ...), br = rq.fit.br(x, y, tau = tau, ...), lasso = rq.fit.lasso(x, 
            y, tau = tau, ...),   = rq.fit.scad(x, y, tau = tau, 
            ...), {
            what <- paste("rq.fit.", method, sep = "")
            if (exists(what, mode = "function")) (get(what, mode = "function"))(x, 
                y, ...) else stop(paste("unimplemented method:", 
                method))
        })
    fit$fitted.values <- y - fit$residuals
    fit$contrasts <- attr(x, "contrasts")
    fit
}

