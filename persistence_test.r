# teste teste
source("write_load.r")
library(Kendall)
source("functions_persistence.r")
dyn.load("persistence_tools.so")

seasonal_test <- function(dat,shift=100,interval=365,seasons=c(0,365),model=shock_ar_1){
    size=length(dat)
    x=seq(1, size, 1)
    i=100
    j=1

    summer=array(NA,62)
    winter=array(NA,62)
    bic_s=array(NA,62)
    bic_w=array(NA,62)
    while ((j*interval+shift)<size){
        i=(j-1)*interval+shift
        if ((is.na(dat[i])==FALSE) & (is.na(dat[i+interval])==FALSE)){
            for (sea in 1:(length(seasons)-1)){
                x=dat[(seasons[sea]+i):(seasons[sea+1]+i)]
				Py=model(x)
                if (sea==1){
                    summer[j]=Py
                }
                if (sea==2){
                    winter[j+1]=Py
                }
            }
        }

        j=j+1
    }
    return(list(shock_s=summer,shock_w=winter,bic_s=bic_s,bic_w=bic_w))
}


dat=dat_load("../data/mid_lat.nc")
trend=trend_load("../data/91_5_trend.nc")
per=per_load("../data/91_5_per_shock_first_test.nc")

if (TRUE){
	start=1
	stop=62
	station=which(dat$ID==572)
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
		shock_s=tmp$shock_s
		shock_bic=tmp$bic
		print(tmp)
		print(proc.time()[1]-time0)

		time0=proc.time()[1]
		#tmp=seasonal_test(as.vector(data),100,365,c(0,182,365),shock_ar_2)
		#shock_s=tmp$cum_s
		#print(tmp)
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

	    for (i in 1:6){
	        markov[i,][markov[i,] == 99]=NA
	    }

		print(proc.time()[1]-time0)
	}


	x=seq(1, size, 1)
    time = c(1:(length(dat$day)*(stop-start+1)))/365 - 0.5*1/365 + min(dat$year+start)

    pdf(file=sprintf("../plots/persistence_tests_%s_%s_572.pdf",dat$lon[station],dat$lat[station]))
	
    par(mfrow=c(2,1))

    par(plt=c(0.12,0.95,0.1,0.95))
    plot(NA,xlim=c(1950,2011),ylim=c(0.5,2),ylab="persistence indices")
	lines(dat$year,shock_s,lty=1,pch=15,col="black")
	lines(dat$year,markov[5,],lty=1,col="red")
	lines(dat$year,markov[6,],lty=1,col="blue")
	legend("topright", pch = c(15,NA,15,15,15), col = c("green", NA, "red", "blue","black"), 
		#legend = c(sprintf("markov %f",coefs[1]), NA, sprintf("ar 2 %f",coefs[2]),sprintf("ar 1 %f",coefs[3]), sprintf("ma 1 %f",coefs[4])))
    

    par(new=TRUE,plt=c(0.13,0.35,0.7,0.9))

	q=station
	library(rworldmap)
	worldmap = getMap(resolution = "low")
	plot(worldmap,xlim=c(dat$lon[q]-10,dat$lon[q]+10),ylim=c(dat$lat[q]-10,dat$lat[q]+10))
	points(dat$lon[q],dat$lat[q],pch=15,col="red")



	graphics.off()



    if(2==3){
    	print(dat$lon[station])
    	pdf(file=sprintf("../plots/bic_tests_%s_%s_572.pdf",dat$lon[station],dat$lat[station]))
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