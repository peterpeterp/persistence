# teste teste
source("write_load.r")
library(Kendall)
source("functions_persistence.r")
dyn.load("persistence_tools.so")


persistence_test <- function(station,heat_waves=c(2003,2006),start=1,stop=62){
	station=which(dat$ID==station)
	data = dat$tas[station,1:365,start:stop]-trend[station,1:365,start:stop]
	size=length(data)
	per_ind1D=as.vector(per$ind[station,1:365,start:stop])
	per_ind1D[per_ind1D==0]=NA

	if (7==7){
		time0=proc.time()[1]
		#tmp=seasonal(as.vector(data),array(c(151,242),dim=c(2,1)),bic_selective)
		tmp=seasonal(as.vector(data),array(c(151,242),dim=c(2,1)),shock_ma,3)
		shock_1=tmp$summer_w
		shock_1_bic=tmp$bic_s
		print(tmp)
		print(proc.time()[1]-time0)

		time0=proc.time()[1]
		#tmp=seasonal(as.vector(data),array(c(151,242),dim=c(2,1)),bic_selective)
		tmp=seasonal(as.vector(data),array(c(151,242),dim=c(2,1)),shock_ma,2)
		shock_2=tmp$summer_w
		shock_2_bic=tmp$bic_s
		print(tmp)
		print(proc.time()[1]-time0)

		time0=proc.time()[1]
		#tmp=seasonal(as.vector(data),array(c(151,242),dim=c(2,1)),bic_selective)
		tmp=seasonal(as.vector(data),array(c(151,242),dim=c(2,1)),shock_ar,4)
		shock_3=tmp$summer_w
		shock_3_bic=tmp$bic_s
		print(tmp)
		print(proc.time()[1]-time0)
	}

	if(3==3){
		time0=proc.time()[1]
		#tmp=seasonal(per_ind1D,100,365,c(0,182,365),markov)
		tmp=seasonal(per_ind1D,array(c(151,242),dim=c(2,1)),markov_chft)
	    markov_s_w=tmp$summer_w
	    markov_s_c=tmp$summer_c
		print(tmp)

		print(proc.time()[1]-time0)
		}


	if (2==3){
		time0=proc.time()[1]
	    #tmp=seasonal(per_ind1D,100,365,c(0,182,365),markov)
	    tmp=seasonal(per_ind1D,array(c(151,242),dim=c(2,1)),markov_calc)
	    c_markov_s_w=tmp$summer_w
	    c_markov_s_c=tmp$summer_c
		print(tmp)

		print(proc.time()[1]-time0)
	}


	x=seq(1, size, 1)
    time = c(1:(length(dat$day)*(stop-start+1)))/365 - 0.5*1/365 + min(dat$year+start)

    pdf(file=sprintf("../plots/persistence_test/persistence_%s_%s_%s_%s_%s.pdf",nday,nyr,dat$lon[station],dat$lat[station],dat$ID[station]))
	
    par(mfrow=c(2,1))
    par(plt=c(0.12,0.95,0.1,0.95))
    plot(NA,xlim=c(1950,2011),ylim=c(-3,3),ylab="trend")
    points(dat$time,dat$tash[station])
    lines(time,trend[station,,])
    #lines(time,dat$tas[station,,],col="red")

    par(plt=c(0.12,0.95,0.1,0.95))    
    plot(NA,xlim=c(1950,2011),ylim=c(-0.7,0.7),ylab="persistence indix anomalie")
	lines(dat$year,(shock_1-mean(shock_1,na.rm=TRUE)),lty=1,pch=15,col="red")
	#lines(dat$year,(shock_2-mean(shock_2,na.rm=TRUE)),lty=1,pch=15,col="blue")
	#lines(dat$year,(shock_3-mean(shock_3,na.rm=TRUE)),lty=1,pch=15,col="green")
	lines(dat$year,(markov_s_w-mean(markov_s_w,na.rm=TRUE)),lty=1,col="black")
	#lines(dat$year,markov_s_w-mean(markov_s_w,na.rm=TRUE)+0.02,lty=2,col="black")
	#lines(dat$year,markov_s_w+markov_s_c-mean(markov_s_w+markov_s_c,na.rm=TRUE),lty=2,col="black")
	abline(lm((shock_1-mean(shock_1,na.rm=TRUE))~dat$year),col="red")
	abline(lm((markov_s_w-mean(markov_s_w,na.rm=TRUE))~dat$year),col="black")


	grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
	legend("topright", pch = c(15,15), col = c("black", "red"), 
		legend = c("markov","shock ma 3"))
	for (j in heat_waves){
		abline(v=j)
	}

    par(new=TRUE,plt=c(0.13,0.35,0.7,0.9))
	
	q=station
	library(rworldmap)
	worldmap = getMap(resolution = "low")
	plot(worldmap,xlim=c(dat$lon[q]-10,dat$lon[q]+10),ylim=c(dat$lat[q]-10,dat$lat[q]+10))
	points(dat$lon[q],dat$lat[q],pch=15,col="red")

    par(plt=c(0.12,0.95,0.1,0.95))    
	plot(dat$year,(shock_1),lty=1,pch=15,col="blue",cex=0)
	for (i in 1:length(dat$year)){
		text(dat$year[i],(shock_1)[i],label=3,col="red",cex=0.7)
		text(dat$year[i],(shock_2)[i],label=2,col="green",cex=0.7)
		text(dat$year[i],(shock_3)[i],label=4,col="blue",cex=0.7)

	}		
	graphics.off()



    if(2==3){
    	pdf(file=sprintf("../plots/persistence_test/bic_analysis/bic_tests_%s_%s_%s.pdf",dat$lon[station],dat$lat[station],dat$ID[station]))

		plot(dat$year,shock_1_bic,pch=15,col="red",cex=0.0)
		lines(dat$year,shock_1_bic,col="red")
		abline(h=mean(shock_1_bic[1:61],na.rm=TRUE),col="red")

		lines(dat$year,shock_2_bic,col="blue")
		abline(h=mean(shock_2_bic[1:61],na.rm=TRUE),col="blue")

		lines(dat$year,shock_3_bic,col="green")
		abline(h=mean(shock_3_bic[1:61],na.rm=TRUE),col="green")


		legend("topright", pch = c(15,15,15), col = c("red", "blue","green"), 
			legend = c(sprintf("ma 1 %f",mean(shock_1_bic[1:61],na.rm=TRUE)),sprintf("ma 2 %f",mean(shock_2_bic[1:61],na.rm=TRUE)), sprintf("ma 3 %f",mean(shock_3_bic[1:61],na.rm=TRUE))))

	}
}

nday=45
nyr=1

dat=dat_load("../data/mid_lat.nc")
trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend.nc",nday,nyr,nday,nyr))
per=markov_load(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr))

pos=c(488,577,190,489,744,97,270,307)
pos=c(487)
waves=c(2003,2006)
for (i in 1:length(pos)){
	persistence_test(pos[i],waves)
}
