source("write.r")
source("load.r")
source("functions_support.r")

if (1==1){
	nc=open.ncdf("../data/91_5/91_5_markov3s.nc")
	print(nc)
	dfsdf
	dat=dat_load("../data/dat_regional.nc")
	seasons=c("spring","summer","autumn","winter","year")
	for (season in seasons){
		tmp=get.var.ncdf(nc,paste("markov_",season,sep=""))
		tmp=global_trend(per=tmp,filename_neu=paste("../data/91_5/91_5_mar3s_trend_",season,".nc",sep=""),season=season,transition_names="cc cn cw nc nn nw wc wn ww")
	}
}

if (1==2){
	nc=open.ncdf("../data/91_5/91_5_markov3s.nc")
	dat=dat_load("../data/dat_regional.nc")
	q=488
	pdf(file="../plots/3states!!!")
	tmp=get.var.ncdf(nc,"markov_summer")
	summer=array(tmp[q,,],dim=c(3,3,62))

	states=c("cold","normal","warm")
	label=c()

	jet.colors <- colorRampPalette( c( "violet","blue","green","yellow","red") )
	nbcol <- 9
	color <- jet.colors(nbcol)

	plot(NA,xlim=c(1950,2011),ylim=c(0,1))
	for (from in 1:3){
		for (to in 1:3){
			lines(dat$year,summer[from,to,],col=color[(from-1)*3+to])
			#print(summary(lm(summer[from,to,]~dat$year)))
			label[(from-1)*3+to]=paste(states[from],"to",states[to],summary(lm(summer[from,to,]~dat$year))$coefficients[2])
			abline(lm(summer[from,to,]~dat$year),col=color[(from-1)*3+to])
		}
	}
	abline(v=c(2003,2006))

	legend("topleft",lty=array(1,9),col=color,legend=label)

}

if(1==2){
	q=488
	trend=trend_load("../data/91_5/91_5_trend_r.nc")

	detrended=dat$tas[q,,]-trend[q,,]
	threshold=sd(detrended,na.rm=TRUE)*0.5


	pdf(file="../plots/3_state_description")
	plot(dat$time,dat$tas[q,,],xlim=c(2002,2004),pch=20,cex=0.5,main="3 states",ylab="temperature anomaly in deg C",xlab="")
	lines(dat$time,trend[q,,],col="green")
	lines(dat$time,trend[q,,]+threshold,col="grey")
	lines(dat$time,trend[q,,]-threshold,col="grey")
	warm=dat$tas[q,,]
	warm[warm<trend[q,,]+threshold]=NA
	points(dat$time,warm,pch=20,cex=0.5,col="red")
	cold=dat$tas[q,,]
	cold[cold>trend[q,,]-threshold]=NA
	points(dat$time,cold,pch=20,cex=0.5,col="blue")
	legend("bottomright",pch=c(NA,20,20,20),col=c(NA,"red","black","blue"),legend=c("diff between grey lines = 1 sd","warm day","average day","cold day"))
}