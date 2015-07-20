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


dd=read.table("../data/dur_488.txt")
size=dim(dd)[1]
test=rq(dd[1:size,2]~dd[1:size,1],0.5)
print(test)
test=rq(dd[1:size,2]~dd[1:size,1],0.99)
print(test)
print(summary(test))
F=ecdf(dd[1:size,2])
print(F)
print(F(14))
print(knots(F))
print(summary(F))
summary.stepfun(F)
per95=which(F(dd[1:size,2])>0.95)
print(dd[per95,1])
print(lm(dd[per95,2]~dd[per95,1]))
print(rq(dd[1:size,2]~dd[1:size,1],0.95))