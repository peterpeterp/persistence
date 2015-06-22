# teste teste
source("write_load.r")
library(Kendall)
source("functions_persistence.r")
dyn.load("persistence_tools.so")


persistence_test <- function(station,heat_waves=c(2003,2006),start=1,stop=62){
	station=which(dat$ID==station)
	size=length(data)


	x=seq(1, size, 1)
    time = c(1:(length(dat$day)*(stop-start+1)))/365 - 0.5*1/365 + min(dat$year+start)

    pdf(file=sprintf("../plots/jjay_watch/persistence_%s_%s_%s_%s_%s.pdf",nday,nyr,dat$lon[station],dat$lat[station],dat$ID[station]))
	
    par(mfrow=c(2,1))
    par(plt=c(0.12,0.95,0.1,0.95))
    plot(NA,xlim=c(1950,2011),ylim=c(-3,3),ylab="trend")
    points(dat$time,dat$tash[station])
    lines(time,trend[station,,])

    print(per$markov[station,1,])
    par(plt=c(0.12,0.95,0.1,0.95))    
    plot(NA,xlim=c(1950,2011),ylim=c(-0.7,0.7),ylab="persistence indix anomalie")
	lines(dat$year,(per$markov[station,1,]-mean(per$markov[station,1,],na.rm=TRUE)),lty=1,pch=15,col="red")
	lines(dat$year,(per$markov[station,2,]-mean(per$markov[station,2,],na.rm=TRUE)),lty=1,pch=15,col="blue")
	lines(dat$year,(per$markov[station,3,]-mean(per$markov[station,3,],na.rm=TRUE)),lty=1,pch=15,col="green")


	#abline(lm((shock_1-mean(shock_1,na.rm=TRUE))~dat$year),col="red")
	#abline(lm((markov_s_w-mean(markov_s_w,na.rm=TRUE))~dat$year),col="black")


	grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
	legend("topright", pch = c(15,15), col = c("black", "red"), 
		legend = c("markov","shock ma 3"))
	for (j in heat_waves){
		abline(v=j)
	}

 
	graphics.off()

}

nday=91
nyr=5

dat=dat_load("../data/mid_lat.nc")
trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend.nc",nday,nyr,nday,nyr))
per=markov_jjay_load(sprintf("../data/%s_%s/%s_%s_markov_jjay.nc",nday,nyr,nday,nyr))

pos=c(488,577,190,489,744,97,270,307)
pos=c(190)
waves=c(2003,2006)
for (i in 1:length(pos)){
	persistence_test(pos[i],waves)
}
