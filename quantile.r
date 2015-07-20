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

if (1==1){
	dd=read.table("../data/dur_488.txt")
	size=dim(dd)[1]
	size=length(which(!is.na(dd[1:size,1])))

    pdf(file="../plots/quantile.pdf")
    plot(dd[1:size,1],dd[1:size,2])

    for (perc in c(0.25,0.5,0.9,0.95,0.98,0.99)){
		test=summary(rq(dd[1:1100,2]~dd[1:1100,1],0.95))$coefficients
		print(test[c(1,7,2,8)])
		print(test$coefficients[1])
		print(test$coefficients[7])
		print(test$coefficients[2])
		print(test$coefficients[8])
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