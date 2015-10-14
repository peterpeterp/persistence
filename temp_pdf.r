source("load.r")

if (1==2){
	library(MASS)
	library(moments)

	trendID="91_5"
	trend_median=trend_load(paste("../data/",trendID,"/",trendID,"_trend_median.nc",sep=""))
	dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
	q=1011

	y=dat$tas[q,,]
	nona=which(!is.na(y))
	y=y[nona]
	tre=trend_median[q,,]
	tre=tre[nona]


	warm=which(trend_median[q,,]>=dat$tas[q,,])
	kalt=which(trend_median[q,,]<dat$tas[q,,])
	druuf=which(trend_median[q,,]==dat$tas[q,,])

	print(length(warm))
	print(length(kalt))
	print(length(kalt)/(length(warm)+length(kalt)))
	print(length(druuf)/length(dat$tas[q,,]))

	warm=which(tre>=y)
	kalt=which(tre<y)
	druuf=which(tre==y)

	print(y[druuf][1:10])
	print(tre[druuf][1:10])
	print(y[druuf][15])

	print(length(warm))
	print(length(kalt))
	print(length(kalt)/(length(warm)+length(kalt)))
	print(length(druuf)/length(y))
}


estimate_mode <- function(x,na.rm=TRUE) {
	if (na.rm==TRUE){
		x=x[which(!is.na(x))]
	}
	d <- density(x)
	return(d$x[which.max(d$y)])
}

if (1==1){
	library(MASS)
	library(moments)

	dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
	pdf(file="../plots/zwischenzeugs/temp_pdf.pdf")

	br=seq(-30,30,0.3)

	q=152
	temp=hist(dat$tas[q,,],breaks=br,plot=FALSE)
	plot(temp$mids,temp$density)
	abline(v=mean(dat$tas[q,,],na.rm=TRUE),lty=3)
	abline(v=median(dat$tas[q,,],na.rm=TRUE),lty=2)
	abline(v=estimate_mode(dat$tas[q,,],na.rm=TRUE),lty=1)

	y=as.vector(dat$tas[q,,])
	nona=which(!is.na(y))
	y=y[nona]
	fit=fitdistr(y,"normal")

	#lines(temp$mids,(0.1*exp(-(temp$mids-fit$estimate[1])^2/fit$estimate[2]^2)),col="black")
	lines(density(y,lty=2))

	tab <- data.frame(x=temp$mids, r=temp$density)
	res <- nls( r ~ k*exp(-1/2*(x-mu)^2/sigma^2),data=tab, start=c(mu=0.2,sigma=3,k=1))
	v <- summary(res)$parameters[,"Estimate"]
	#lines(temp$mids,(v[3]*exp(-(temp$mids-v[1])^2/v[2]^2)),col="blue",lty=2)



	print(fit)
	print(v)
	print(mean(y))

	schiefe=1/length(y)*sum(((y-fit$estimate[1])/(fit$estimate[2]))^3)
	print(schiefe)
	print(skewness(y))


	q=351
	temp=hist(dat$tas[q,,],breaks=br,plot=FALSE)
	points(temp$mids,temp$density,col="red",pch=4)
	abline(v=mean(dat$tas[q,,],na.rm=TRUE),lty=3,col="red")
	abline(v=median(dat$tas[q,,],na.rm=TRUE),lty=2,col="red")
	abline(v=estimate_mode(dat$tas[q,,],na.rm=TRUE),lty=1,col="red")

	y=as.vector(dat$tas[q,,])
	nona=which(!is.na(y))
	y=y[nona]
	fit=fitdistr(y,"normal")

	#lines(temp$mids,(0.1*exp(-(temp$mids-fit$estimate[1])^2/fit$estimate[2]^2)),col="red")
	lines(density(y,lty=2),col="red")

	tab <- data.frame(x=temp$mids, r=temp$density)
	res <- nls( r ~ k*exp(-1/2*(x-mu)^2/sigma^2),data=tab, start=c(mu=0,sigma=3,k=1))
	v <- summary(res)$parameters[,"Estimate"]
	#lines(temp$mids,(v[3]*exp(-(temp$mids-v[1])^2/v[2]^2)),col="green",lty=2)

	print(fit)
	print(v)
	print(mean(y))

	schiefe=1/length(y)*sum(((y-fit$estimate[1])/(fit$estimate[2]))^3)
	print(schiefe)
	print(skewness(y))
}


