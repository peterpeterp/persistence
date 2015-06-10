# teste teste
source("write_load.r")
library(Kendall)
source("functions_persistence.r")
dyn.load("persistence_tools.so")




dat=dat_load("../data/dat_205-210.nc")
trend=trend_load("../data/trend_205-210.nc")
per=per_load("../data/per_205-210.nc")

if (TRUE){
	start=1
	stop=62
	station=3
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
		tmp=seasonal_bic(as.vector(data),100,365,c(0,182,365),"ar_2")
		print(tmp)
		print(proc.time()[1]-time0)
		time0=proc.time()[1]

		ar_2_s=tmp$cum_s
		ar_2_bic=tmp$bic
		tmp=seasonal(as.vector(data),100,365,c(0,182,365),"ar_2")

		print(proc.time()[1]-time0)
		ar_1_s=tmp$cum_s
		ar_1_bic=tmp$bic

		time0=proc.time()[1]

		tmp=seasonal(as.vector(data),100,365,c(0,365),"ma_1")
		ma_1_s=tmp$cum_s
		ma_1_bic=tmp$bic
		print(proc.time()[1]-time0)
	}



	if (2==2){
		per_ind1D=as.vector(per$ind[station,1:365,start:stop])
		time0=proc.time()[1]
	    seasons=c(182,183)
	    startzzz=100

    	markov=array(99,dim=c(6,62))

	    win_warm = array(99,62)
	    win_cold = array(99,62)
	    sum_warm = array(99,62)
	    sum_cold = array(99,62)

	    temp=.C("c_markov_season",data=as.integer(per_ind1D),warm_w=as.numeric(markov[3,]),cold_w=as.numeric(markov[4,]),warm_s=as.numeric(markov[5,]),cold_s=as.numeric(markov[6,])
	            ,startzzz=as.integer(start),season_length=as.integer(seasons),season_number=as.integer(length(seasons)),size=as.integer(size))

	    markov[3,]=temp$warm_w
	    markov[4,]=temp$cold_w
	    markov[5,]=temp$warm_s
    	markov[6,]=temp$cold_s
		print(proc.time()[1]-time0)
	}


	x=seq(1, size, 1)
    time = c(1:(length(dat$day)*(stop-start+1)))/365 - 0.5*1/365 + min(dat$year+start)

    pdf(file=sprintf("../plots/persistence_tests_%s_%s_year.pdf",dat$lon[station],dat$lat[station]))

    markov_s=markov[5,]
    markov_s[markov_s==99]=NA

	plot(dat$year,markov_s,pch=15,col="black",ylim=c(0.5,2))
	points(dat$year,ar_2_s,pch=15,col="red")
	points(dat$year,ar_1_s,pch=15,col="blue")
	points(dat$year,ma_1_s,pch=15,col="green")
	abline(lm(markov_s~dat$year),col="black")
	abline(lm(ar_2_s~dat$year),col="red")
	abline(lm(ar_1_s~dat$year),col="blue")
	abline(lm(ma_1_s~dat$year),col="green")

	abline(v=2006)
	coefs=c(1,1,1,1)
	coefs[1]=lm(markov_s~dat$year)$coefficients[2]
	coefs[2]=lm(ar_2_s~dat$year)$coefficients[2]
	coefs[3]=lm(ar_1_s~dat$year)$coefficients[2]
	coefs[4]=lm(ma_1_s~dat$year)$coefficients[2]

	legend("topright", pch = c(15,NA,15,15,15), col = c("green", NA, "red", "blue","black"), 
		legend = c(sprintf("markov %f",coefs[1]), NA, sprintf("ar 2 %f",coefs[2]),sprintf("ar 1 %f",coefs[3]), sprintf("ma 1 %f",coefs[4])))


    if(2==2){
    	print(dat$lon[station])
    	pdf(file=sprintf("../plots/bic_tests_%s_%s_year.pdf",dat$lon[station],dat$lat[station]))
		plot(dat$year,ar_2_bic,pch=15,col="red")
		abline(h=mean(ar_2_bic[1:61]),col="red")

		points(dat$year,ar_1_bic,pch=15,col="blue")
		abline(h=mean(ar_1_bic[1:61]),col="blue")

		points(dat$year,ma_1_bic,pch=15,col="green")
		abline(h=mean(ma_1_bic[1:61]),col="green")
		legend("topright", pch = c(15,15,15), col = c("red", "blue","black"), 
			legend = c(sprintf("ar 2 %f",mean(ar_2_bic[1:61])),sprintf("ar 1 %f",mean(ar_1_bic[1:61])), sprintf("ma 1 %f",mean(ma_1_bic[1:61]))))


	}
}