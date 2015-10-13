source("load.r")
library(MASS)
library(moments)


dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")

pdf(file="../plots/zwischenzeugs/temp_pdf.pdf")

br=seq(-30,30,1)

q=325
temp=hist(dat$tas[q,,],breaks=br,plot=FALSE)
plot(temp$mids,temp$density)
abline(v=mean(dat$tas[q,,],na.rm=TRUE),lty=2)
abline(v=median(dat$tas[q,,],na.rm=TRUE))

y=as.vector(dat$tas[q,,])
nona=which(!is.na(y))
y=y[nona]
fit=fitdistr(y,"normal")

lines(temp$mids,(0.1*exp(-(temp$mids-fit$estimate[1])^2/fit$estimate[2]^2)),col="black")
lines(density(y,lty=2))

tab <- data.frame(x=temp$mids, r=temp$density)
res <- nls( r ~ k*exp(-1/2*(x-mu)^2/sigma^2),data=tab, start=c(mu=0.2,sigma=3,k=1))
v <- summary(res)$parameters[,"Estimate"]
lines(temp$mids,(v[3]*exp(-(temp$mids-v[1])^2/v[2]^2)),col="blue",lty=2)

print(fit)
print(v)
print(mean(y))

schiefe=1/length(y)*sum(((y-fit$estimate[1])/(fit$estimate[2]))^3)
print(schiefe)
print(skewness(y))


q=152
temp=hist(dat$tas[q,,],breaks=br,plot=FALSE)
points(temp$mids,temp$density,col="red",pch=4)
abline(v=mean(dat$tas[q,,],na.rm=TRUE),col="red",lty=2)
abline(v=median(dat$tas[q,,],na.rm=TRUE),col="red")

y=as.vector(dat$tas[q,,])
nona=which(!is.na(y))
y=y[nona]
fit=fitdistr(y,"normal")

lines(temp$mids,(0.1*exp(-(temp$mids-fit$estimate[1])^2/fit$estimate[2]^2)),col="red")

tab <- data.frame(x=temp$mids, r=temp$density)
res <- nls( r ~ k*exp(-1/2*(x-mu)^2/sigma^2),data=tab, start=c(mu=0,sigma=3,k=1))
v <- summary(res)$parameters[,"Estimate"]
lines(temp$mids,(v[3]*exp(-(temp$mids-v[1])^2/v[2]^2)),col="green",lty=2)

print(fit)
print(v)
print(mean(y))

schiefe=1/length(y)*sum(((y-fit$estimate[1])/(fit$estimate[2]))^3)
print(schiefe)
print(skewness(y))

