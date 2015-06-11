# teste teste
source("write_load.r")



trend_plot <- function(dat,trend,per,filename,station,start,stop){
	# Make a simple test plot
    trash = ((nyr-1)/2*365+(nday-1))
	ntot=length(dat$ID)
	qq=which(dat$lon>-10 & dat$lon<30 & dat$lat>40 & dat$lat<60)
	#qq=which(dat$ID==station)    
    xlim  = c(start+1949,stop+1949)

    ylim1 = c(-40,18)
    ylim2 = c(0.6,2)
    ylim3 = c(1.4,1.8)

    #myfigure(fldr,"peter_per_test",asp=1.0,pointsize=16)
    pdf(file = filename)
    par(mfrow=c(2,1))
    par(plt=c(0.1,0.95,0.18,0.95))

    for (ip in qq){
    	if (length(which(is.na(dat$tas[ip,,])))<(2*trash+730)){

		    plot(xlim,ylim1,type="n",xlab="",ylab="Temperature anomaly (°C)")
		    abline(h=0,col="grey40")
   		    lines(dat$time,dat$tas[ip,,],col="black")

		    lines(dat$time,trend[ip,,],lwd=3,col=2)

		    plot(xlim,ylim2,type="n",xlab="",ylab="shock persistence")
		    lines(dat$year[start:stop],per$shock[ip,1,start:stop],type="l",lwd=2,xlab="",ylab=dat$lon[ip],col="green")
		    lines(dat$year[start:stop],per$markov[ip,5,start:stop],type="l",lwd=2,xlab="",ylab=dat$lat[ip],col="red")
		    lines(dat$year[start:stop],per$markov[ip,6,start:stop],type="l",lwd=2,xlab="",ylab=dat$lat[ip],col="blue")
		    abline(lm(per$shock[ip,1,start:stop]~dat$year[start:stop]),col="green")
		    abline(lm(per$markov[ip,5,start:stop]~dat$year[start:stop]),col="red")
		    abline(lm(per$markov[ip,6,start:stop]~dat$year[start:stop]),col="blue")
			
			abline(v=2003)

			legend("topright", col = c(NA,NA,"red", "blue"), 
			legend = c(sprintf("ID %s lon %s  lat %s",dat$ID[ip],dat$lon[ip],dat$lat[ip]),NA,"summer","winter"))


		    #plot(xlim,ylim3,type="n",xlab="",ylab="Tendency")
		    #lines(dat$year[start:stop],(per$sum_cold[ip,start:stop]+per$sum_warm[ip,start:stop]),type="l",lwd=2,xlab="",ylab="Tendency",col="green")
		    #abline(lm((per$sum_cold[ip,start:stop]+per$sum_warm[ip,start:stop])~dat$year[start:stop]),col="green")

		}
		else {
            cat(sprintf("> ID %s lon %s  lat %s <",dat$ID[ip],dat$lon[ip],dat$lat[ip]))
        }
	}

    graphics.off()
}


if (1==1){
	nday=91
	nyr=5
	dat=dat_load("../data/mid_lat.nc")
	trend=trend_load(sprintf("../data/%s_%s_trend.nc",nday,nyr))
	per=per_load(sprintf("../data/91_5_per_shock_first_test.nc",nday,nyr))

	if(2==5){
		qq=which(dat$lon>-15 & dat$lon<30 & dat$lat>35 & dat$lat<60)
		dat_klein=dat_write_part(dat,qq,"../data/dat_west_europe.nc")
		trend_klein=trend_write_part(trend,dat_klein,qq,"../data/trend_west_europe.nc")
		per_klein=per_write_part(per,dat_klein,qq,"../data/per_west_europe.nc")
	}

	trend_plot(dat,trend,per,sprintf("../plots/persistence_%s_%s.pdf",nday,nyr),209,1,62)



	if (1==2){
		for (i in 1:length(dat$ID)){
			if((per$ma_warm_trend[i]+per$ma_cold_trend[i])>0.002){
				print(dat$ID[i])
				print(dat$lon[i])
				print(dat$lat[i])
				print((dat$ma_warm_trend[i]+dat$ma_cold_trend[i]))
			}
		}
	}
}

