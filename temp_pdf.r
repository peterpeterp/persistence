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
	print(length(druuf)/length(y)*100)
}


estimate_mode <- function(x,na.rm=TRUE) {
	if (na.rm==TRUE){
		x=x[which(!is.na(x))]
	}
	d <- density(x)
	return(d$x[which.max(d$y)])
}

if (1==2){
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

if (1==2){
	q=52
	dat=dat_load("../data/HadGHCND_TMean_data3D.day1-365.1950-2014.nc")
	pdf(file=paste("../plots/zwischenzeugs/temp_pdf_",q,".pdf",sep=""))

	br=seq(-30,30,0.3)

	color=rgb(0.6,0.8,0.5,0.5)
	temp=hist(dat$tas[q,152:243,],breaks=br,plot=TRUE,col=color,border=color)
	color=rgb(0.6,0.8,0.5,1)
	abline(v=median(dat$tas[q,152:243,],na.rm=TRUE),col=color)
	#color=rgb(0.8,0.6,0.5,0.5)
	#temp=hist(dat$tas[q,244:334,],breaks=br,plot=TRUE,col=color,border=color,add=TRUE)
	#color=rgb(0.8,0.4,0.8,0.5)
	#temp=hist(dat$tas[q,60:151,],breaks=br,plot=TRUE,col=color,border=color,add=TRUE)
	color=rgb(0.6,0.5,0.6,0.5)
	temp=hist(dat$tas[q,c(1:60,335:365),],breaks=br,plot=TRUE,col=color,border=color,add=TRUE)
	color=rgb(0.6,0.5,0.6,1)
	abline(v=median(dat$tas[q,c(1:60,335:365),],na.rm=TRUE),col=color)

	#abline(v=0)
}

if (1==2){

	q1=488
	q2=1233
	dat=dat_load("../data/HadGHCND_TMean_data3D.day1-365.1950-2014.nc")
	pdf(file=paste("../plots/zwischenzeugs/temp_pdf/temp_pdf_",q1,"_",q2,".pdf",sep=""))

	br=seq(-30,30,0.3)

	q=q1
	color=rgb(0.6,0.8,0.5,0.5)
	temp=hist(dat$tas[q,,],breaks=br,plot=TRUE,col=color,border=color)
	color=rgb(0.6,0.8,0.5,1)
	abline(v=median(dat$tas[q,,],na.rm=TRUE),col=color)

	q=q2
	color=rgb(0.6,0.5,0.6,0.5)
	temp=hist(dat$tas[q,,],breaks=br,plot=TRUE,col=color,border=color,add=TRUE)
	color=rgb(0.6,0.5,0.6,1)
	abline(v=median(dat$tas[q,,],na.rm=TRUE),col=color)

	#abline(v=0)
}

plot_dist_cold_warm_stupid <-function(q=52,median=FALSE){
	dat=dat_load("../data/HadGHCND_TMean_data3D.day1-365.1950-2014.nc")
	pdf(file=paste("../plots/zwischenzeugs/temp_pdf/temp_pdf_stupid_",q,"_",median,".pdf",sep=""),width=4,height=4)

	br=seq(-30,30,1)

	y=dat$tas[q,c(1:60,335:365),]
	color="white"
	#plot(NA,xlim=c(-20,20),ylim=c(0,600),axes=TRUE,ylab="counts",xlab="temperature anomaly",main="",frame=FALSE)
	temp=hist(y,plot=TRUE,col=color,border=color,axes=TRUE,ylab="counts",xlab="temperature anomaly",main="")
	color=rgb(0.5,0.5,1,0.7)
	if (median==TRUE){temp=hist(y[y<=median(y,na.rm=TRUE)],plot=TRUE,col=color,border=color,add=TRUE)}
	if (median==FALSE){temp=hist(y[y<=0],breaks=br,plot=TRUE,col=color,border=color,add=TRUE)}
	color=rgb(1,0.5,0.5,0.7)
	if (median==TRUE){temp=hist(y[y>median(y,na.rm=TRUE)],plot=TRUE,col=color,border=color,add=TRUE)}
	if (median==FALSE){temp=hist(y[y>0],breaks=br,plot=TRUE,col=color,border=color,add=TRUE)}
	
	text(-15,200,"cold?")
	text(15,200,"warm?")

	#abline(h=0,col="white",cex=1.01)

    par(new=TRUE,plt=c(0.67,0.98,0.6,0.98))
    plot(worldmap,xlim=c(dat$lon[q]-20,dat$lon[q]+20),ylim=c(dat$lat[q]-20,dat$lat[q]+20))
    points(dat$lon[q],dat$lat[q],pch=15,col="violet")
	graphics.off()
}

plot_dur_dist <-function(q=488,median=FALSE,season="4seasons",trendID="91_5",dataset="_TMean",additional_style=""){
	dat=dat_load("../data/HadGHCND_TMean_data3D.day1-365.1950-2014.nc")

	nc_dur=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_",season,".nc",sep=""))
    dur=get.var.ncdf(nc_dur,"dur")
   	dur_mid=get.var.ncdf(nc_dur,"dur_mid")	

   	pdf(file=paste("../plots/zwischenzeugs/temp_pdf/dur_pdf_",q,"_",season,".pdf",sep=""),width=4,height=4)

	br=seq(0,100,1)

	y=dur[q,2,]
	qua=quantile_pete(y,taus=c(0.05,0.25,0.5,0.75,0.95))

	color=rgb(0.5,0.5,0.8,0.5)
	plot(NA,xlim=c(0,60),ylim=c(0,400),axes=TRUE,ylab="counts",xlab="period length [days]",main="",frame=FALSE)
	abline(v=qua,col=rgb(0.6,0.6,0.6,0.5))

	temp=hist(y,plot=TRUE,ylim=c(0,400),xlim=c(0,60),breaks=br,col=color,border=color,add=TRUE)
	#temp=hist(y,plot=TRUE,ylim=c(0,400),xlim=c(0,60),breaks=br,col=color,border=color,axes=TRUE,ylab="counts",xlab="period length [days]",main="")

	#for (qua in quants){abline(v=qua)}

	#abline(v=quants)
	mi=380
	brei=10
	oben=mi+brei
	unten=mi-brei
	polygon(x=c(qua[2],qua[2],qua[4],qua[4]),y=c(unten,oben,oben,unten),col=color,border="black")
	lines(c(qua[1],qua[1]),c(unten,oben))
	lines(c(qua[3],qua[3]),c(unten,oben))
	lines(c(qua[5],qua[5]),c(unten,oben))
	lines(c(qua[1],qua[2]),c(mi,mi))
	lines(c(qua[4],qua[5]),c(mi,mi))

    par(new=TRUE,plt=c(0.67,0.88,0.6,0.88))
    plot(worldmap,xlim=c(dat$lon[q]-10,dat$lon[q]+10),ylim=c(dat$lat[q]-20,dat$lat[q]+20))
    points(dat$lon[q],dat$lat[q],pch=15,col="violet")
	graphics.off()
}

source("functions_regional.r")
library(rworldmap)
library(fields)
worldmap = getMap(resolution = "low")

plot_dur_dist(488)

#plot_dist_cold_warm_stupid(52)
#plot_dist_cold_warm_stupid(52,median=TRUE)
#plot_dist_cold_warm_stupid(461)
#plot_dist_cold_warm_stupid(461,median=TRUE)