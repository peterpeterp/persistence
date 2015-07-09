source("write.r")
source("load.r")
source("functions_support.r")

if (1==2){
	nc=open.ncdf("../data/91_5/91_5_markov3s.nc")
	dat=dat_load("../data/dat_regional.nc")
	tmp=get.var.ncdf(nc,"markov_summer")

	global_trend(per=tmp,filename_neu="../data/91_5/91_5_mar3s_trend.nc",season="summer",transition_names=c("cc","cn","cw","nc","nn","nw","wc","wn","ww"))
}

if (1==1){
	w=seq(1,9,1)
	print(w)
	wm=array(w,dim=c(3,3))
	print(wm)
	print(as.vector(wm))
	for (i in 1:3){
		for (j in 1:3){
			print(wm[i,j])
		}
	}
	wewer


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
			print(summary(lm(summer[from,to,]~dat$year)))
			label[(from-1)*3+to]=paste(states[from],"to",states[to],summary(lm(summer[from,to,]~dat$year))$coefficients[2])
			abline(lm(summer[from,to,]~dat$year),col=color[(from-1)*3+to])
		}
	}

	abline(v=c(2003,2006))

	legend("topleft",lty=array(1,9),col=color,legend=label)

	print(summer)
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