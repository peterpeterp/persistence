# teste teste
source("write_load.r")
library(Kendall)
source("functions_persistence.r")
dyn.load("persistence_tools.so")



persistence_test <- function(station,start=1,stop=62){
	station=which(dat$ID==station)
	data = dat$tas[station,1:365,start:stop]
	size=length(data)
	shift=110
	interval=365

	if(1==2){
		time0=proc.time()[1]
		warm_test=array(0,size)
		cold_test=array(0,size)
		for (i in (interval/2):(size-interval/2)) {
			mcFitMLE<-markovchainFit(data=data[(i-interval/2):(i+interval/2)])
			warm_test[i]=mcFitMLE$estimate[2][2]
			cold_test[i]=mcFitMLE$estimate[1][1]
		}
		print(proc.time()[1]-time0)
	}

	if (7==7){
		time0=proc.time()[1]
		tmp=seasonal(as.vector(data),100,365,c(0,182,365),shock_ar_2)
		ar_2_s=tmp$summer_w
		ar_2_bic=tmp$bic_s
		print(tmp)
		print(proc.time()[1]-time0)

		time0=proc.time()[1]
		tmp=seasonal(as.vector(data),100,365,c(0,182,365),shock_ma_1)
		ma_1_s=tmp$summer_w
		ma_1_bic=tmp$bic_s
		print(tmp)
		print(proc.time()[1]-time0)

		time0=proc.time()[1]
		tmp=seasonal(as.vector(data),100,365,c(0,182,365),shock_ar_1)
		ar_1_s=tmp$summer_w
		ar_1_bic=tmp$bic_s
		print(tmp)
		print(proc.time()[1]-time0)
	}



	if (2==2){
		per_ind1D=as.vector(per$ind[station,1:365,start:stop])
		time0=proc.time()[1]
	    tmp=seasonal(per_ind1D,100,365,c(0,182,365),markov)
	    markov_s_w=tmp$summer_w
	    markov_s_c=tmp$summer_c
		print(tmp)

		print(proc.time()[1]-time0)
	}


	x=seq(1, size, 1)
    time = c(1:(length(dat$day)*(stop-start+1)))/365 - 0.5*1/365 + min(dat$year+start)

    pdf(file=sprintf("../plots/persistence_test/persistence_tests_%s_%s_%s.pdf",dat$lon[station],dat$lat[station],dat$ID[station]))
	
    par(mfrow=c(2,1))
    print(mean(ar_2_s,na.rm=TRUE))
    par(plt=c(0.12,0.95,0.1,0.95))
    plot(NA,xlim=c(1950,2011),ylim=c(-0.5,0.5),ylab="persistence indices")
	lines(dat$year,(ar_2_s-mean(ar_2_s,na.rm=TRUE)),lty=1,pch=15,col="red")
	lines(dat$year,ma_1_s-mean(ma_1_s,na.rm=TRUE),lty=1,pch=15,col="blue")
	lines(dat$year,ar_1_s-mean(ar_1_s,na.rm=TRUE),lty=1,pch=15,col="green")
	lines(dat$year,markov_s_w-mean(markov_s_w,na.rm=TRUE),lty=1,col="black")
	lines(dat$year,markov_s_w+markov_s_c-mean(markov_s_w+markov_s_c,na.rm=TRUE),lty=2,col="black")
	grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
	legend("topright", pch = c(15,15,15,15), col = c("black", "red", "blue","green"), 
		legend = c("markov","ar 2","ma_1","ar_1"))
    abline(v=2006)

    par(new=TRUE,plt=c(0.13,0.35,0.7,0.9))
	
	q=station
	library(rworldmap)
	worldmap = getMap(resolution = "low")
	plot(worldmap,xlim=c(dat$lon[q]-10,dat$lon[q]+10),ylim=c(dat$lat[q]-10,dat$lat[q]+10))
	points(dat$lon[q],dat$lat[q],pch=15,col="red")



	graphics.off()



    if(2==2){
    	pdf(file=sprintf("../plots/persistence_test/bic_analysis/bic_tests_%s_%s_%s.pdf",dat$lon[station],dat$lat[station],dat$ID[station]))

		plot(dat$year,ar_2_bic,pch=15,col="red",cex=0.0)
		lines(dat$year,ar_2_bic,col="red")
		abline(h=mean(ar_2_bic[1:61]),col="red")

		lines(dat$year,ar_1_bic,col="blue")
		abline(h=mean(ar_1_bic[1:61]),col="blue")

		lines(dat$year,ma_1_bic,col="green")
		abline(h=mean(ma_1_bic[1:61]),col="green")
		legend("topright", pch = c(15,15,15), col = c("red", "blue","green"), 
			legend = c(sprintf("ar 2 %f",mean(ar_2_bic[1:61])),sprintf("ar 1 %f",mean(ar_1_bic[1:61])), sprintf("ma 1 %f",mean(ma_1_bic[1:61]))))


	}
}

dat=dat_load("../data/mid_lat.nc")
trend=trend_load("../data/91_5_trend.nc")
per=per_load("../data/91_5_per_shock_first_test.nc")

persistence_test(605)