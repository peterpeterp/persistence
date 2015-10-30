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

if (1==2){
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

quantile_pete <- function(dist,taus,na.rm=TRUE){
    if (na.rm==TRUE){dist=dist[which(!is.na(dist))]}

    cdf=array(NA,max(dist))
    out=taus*NA

    cum=0
    for (i in 1:max(dist)){
        cum=cum+length(which(dist>i))
        cdf[i]=cum
    }
    cdf=cdf/cdf[length(cdf)]

    for (qu in 1:length(taus)){
        ueb=which(cdf>taus[qu])[1]
        unt=ueb-1
        if (unt<1){out[qu]=ueb}
        else {out[qu]=ueb-(cdf[ueb]-taus[qu])/(cdf[ueb]-cdf[unt])}
    }  
    return(out)
}

if (1==1){
    dist=read.table("../data/sonstiges/test_dist.txt")
    dist=dist[,1]
    print(quantile(dist,prob=0.95))

    print(quantile_pete(dist,taus=c(0.05,0.25,0.5,0.75,0.95,0.98)))
    print(quantile(dist,prob=c(0.05,0.25,0.5,0.75,0.95,0.98)))
}
