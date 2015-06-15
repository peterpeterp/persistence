# teste teste

a <- matrix(c(-0.7634849,1, -0.4833480,0),nrow=2,ncol=2)
print(a)
print(a %*%a)

px=a[1]+1
amul=a
for (i in 1:1000){
	amul=a %*% amul
	px=px+amul[1]
	print(amul[1])
}
print(px)




source("write_load.r")
library(Kendall)
source("functions_persistence.r")
dyn.load("persistence_tools.so")



persistence_test <- function(station,heat_waves=c(2003,2006),start=1,stop=62){
	station=which(dat$ID==station)
	#data = dat$tas[station,1:365,start:stop]
	shift=110
	interval=365

	if (2==2) {
		tmp = read.table("ocean_anomalie.txt",sep="\t")
		data=tmp[,2]
		time=tmp[,1]
		size=length(data)


		first_dif=data*NA
		for (i in 2:size){
			first_dif[i]=data[i]-data[i-1]
		}
		pdf(file="../plots/land")
		plot(NA,xlim=c(1878,2014),ylim=c(-0.7,1.2))
		lines(time,data)
		graphics.off()
		print(data)
		print(first_dif)
		print(diff(data))
		first_dif=diff(data)

	}

	if (8==9) {
		laenge<-250
		t<-seq(1,laenge,1)
		x<-seq(1,laenge,1)
		y<-seq(1,laenge,1)
		ut<-seq(1,laenge,1)

		for (i in 0:laenge){
			ut[i]=runif(1,-1,1)
		}


		for (i in 0:laenge){
			if(i<3){
				y[i]=ut[i]
				x[i]=ut[i]
			}
			if(i>2){
				x[i]=1.7*x[i-1]-0.7*x[i-2]+ut[i]
				y[i]=0.3*y[i-1]+0.7*y[i-2]+ut[i]
			}
			
		}
		first_dif=y

	}




	if (7==7){
		time0=proc.time()[1]
		tmp=shock_ma(as.vector(first_dif),1)
		shock_1=tmp$summer_w
		shock_1_bic=tmp$bic_s
		print(tmp)
		print(proc.time()[1]-time0)

		time0=proc.time()[1]
		tmp=shock_ma(as.vector(first_dif),2)
		shock_1=tmp$summer_w
		shock_1_bic=tmp$bic_s
		print(tmp)
		print(proc.time()[1]-time0)


		time0=proc.time()[1]
		tmp=shock_ma(as.vector(first_dif),3)
		shock_1=tmp$summer_w
		shock_1_bic=tmp$bic_s
		print(tmp)
		print(proc.time()[1]-time0)

		time0=proc.time()[1]
		tmp=shock_ar_2(as.vector(first_dif))
		shock_2=tmp$summer_w
		shock_2_bic=tmp$bic_s
		print(tmp)
		print(proc.time()[1]-time0)

		time0=proc.time()[1]
		tmp=shock_ar_1(as.vector(first_dif))
		shock_3=tmp$summer_w
		shock_3_bic=tmp$bic_s
		print(tmp)
		print(proc.time()[1]-time0)
		sdfsdf
	}



	if (2==2){
		per_ind1D=as.vector(per$ind[station,1:365,start:stop])
		per_ind1D[per_ind1D==0]=NA
		time0=proc.time()[1]
	    #tmp=seasonal(per_ind1D,100,365,c(0,182,365),markov)
	    tmp=seasonal(per_ind1D,array(c(151,242),dim=c(2,1)),markov_calc)
	    c_markov_s_w=tmp$summer_w
	    c_markov_s_c=tmp$summer_c
		print(tmp)

		print(proc.time()[1]-time0)

		if(3==4){
			time0=proc.time()[1]
		    #tmp=seasonal(per_ind1D,100,365,c(0,182,365),markov)
		    tmp=seasonal(per_ind1D,array(c(151,242),dim=c(2,1)),markov_chft)
		    markov_s_w=tmp$summer_w
		    markov_s_c=tmp$summer_c
			print(tmp)

			print(proc.time()[1]-time0)
		}
	}


	x=seq(1, size, 1)
    time = c(1:(length(dat$day)*(stop-start+1)))/365 - 0.5*1/365 + min(dat$year+start)

    pdf(file=sprintf("../plots/persistence_test/persistence_tests_%s_%s_%s.pdf",dat$lon[station],dat$lat[station],dat$ID[station]))
	
    par(mfrow=c(2,1))
    par(plt=c(0.12,0.95,0.1,0.95))
    plot(NA,xlim=c(1950,2011),ylim=c(-3,3),ylab="trend")
    lines(time,trend[station,,])
    #lines(time,dat$tas[station,,],col="red")

    par(plt=c(0.12,0.95,0.1,0.95))    
    plot(NA,xlim=c(1950,2011),ylim=c(-0.5,0.5),ylab="persistence indix anomalie")
	lines(dat$year,(shock_1-mean(shock_1,na.rm=TRUE)),lty=1,pch=15,col="red")
	lines(dat$year,(shock_2-mean(shock_2,na.rm=TRUE)),lty=1,pch=15,col="blue")
	lines(dat$year,(shock_3-mean(shock_3,na.rm=TRUE)),lty=1,pch=15,col="green")
	lines(dat$year,c_markov_s_w-mean(c_markov_s_w,na.rm=TRUE),lty=1,col="black")
	#lines(dat$year,markov_s_w-mean(markov_s_w,na.rm=TRUE)+0.02,lty=2,col="black")
	#lines(dat$year,markov_s_w+markov_s_c-mean(markov_s_w+markov_s_c,na.rm=TRUE),lty=2,col="black")
	grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
	legend("topright", pch = c(15,15,15,15), col = c("black", "red", "blue","green"), 
		legend = c("markov","ma_1","ar_1","ar_2"))
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
	plot(dat$year,(shock_1-mean(shock_1,na.rm=TRUE)),lty=1,pch=15,col="blue")


	graphics.off()



    if(2==2){
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

dat=dat_load("../data/mid_lat.nc")
trend=trend_load("../data/91_5_trend.nc")
per=per_load("../data/91_5_per_shock_first_test.nc")

persistence_test(777,c(2003,2006))