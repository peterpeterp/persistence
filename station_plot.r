source("write.r")
source("load.r")
source("functions_duration.r")

station_plot <- function(dat,trend,per,dur1,dur2,dur3,dur4,q,start,stop,filename){
    pdf(file=filename)
    par(mfrow=c(2,1))
    par(mar=c(2,2, 2, 0.5)) 
    station=which(dat$ID==q)[1]


    warmeTage=array(NA,62)
    kalteTage=array(NA,62)
    for (i in 1:62){
        warmeTage[i]=length(which(per$ind[station,151:243,i]==1))
        kalteTage[i]=length(which(per$ind[station,151:243,i]==-1)) 
    }

    plot(NA,xlim=c(1950,2010),xaxp=c(1950,2010,12),ylim=c(0,50),xlab="",ylab="duration of persistence")#, main="summer warm persistence")
    points(dur2$dur_warm_mid[station,],dur2$dur_warm[station,],pch=4,col="red")
    markov=per$markov[station,1,]*30
    lines(dat$year+0.5,markov,col="green")
    abline(lm(warmeTage~dat$year),col="yellow")
    for (i in seq(1950,2010,5)){
        abline(v=i,col="gray",lty="dotted")
    }
    lines(dat$time[start:stop],as.vector(trend[station,,])[start:stop],col="black")
    legend("topleft", pch = c(NA,4,NA,NA,NA),lty=c(NA,NA,1,1,1), col = c(NA,"red", "green","yellow","black"), 
            legend = c("summer","warm days in a row","markov cold persistence","# warm days","trend"))

    plot(NA,xlim=c(1950,2011),xaxp=c(1950,2010,12),ylim=c(0,50),xlab="",ylab="duration of persistence")#, main="summer cold persistence")
    points(dur2$dur_cold_mid[station,],dur2$dur_cold[station,],pch=5,col="blue")
    markov=per$markov[station,2,]*30
    lines(dat$year,markov,col="green")
    abline(lm(kalteTage~dat$year),col="violet")
    lines(dat$time[start:stop],as.vector(trend[station,,])[start:stop],col="black")
    for (i in seq(1950,2010,5)){
        abline(v=i,col="gray",lty="dotted")
    }
    legend("topleft", pch = c(NA,5,NA,NA,NA),lty=c(NA,NA,1,1,1), col = c(NA,"blue", "green","violet","black"), 
            legend = c("summer","cold days in a row","markov cold persistence","# cold days","trend"))

    par(new=TRUE,plt=c(0.76,0.96,0.67,0.87))
    q=station
    library(rworldmap)
    worldmap = getMap(resolution = "low")
    plot(worldmap,xlim=c(dat$lon[q]-10,dat$lon[q]+10),ylim=c(dat$lat[q]-10,dat$lat[q]+10))
    points(dat$lon[q],dat$lat[q],pch=15,col="red")

    if (1==1){
        trenn=1980
        br=seq(0,100,2)
        ntot=819
        season="summer"
        vor_warm=which(dur2$dur_warm_mid[q,]<trenn)
        vor_cold=which(dur2$dur_cold_mid[q,]<trenn)
        nach_warm=which(dur2$dur_warm_mid[q,]>trenn)
        nach_cold=which(dur2$dur_cold_mid[q,]>trenn)

        hist1=hist(dur2$dur_warm[q,vor_warm],breaks=br,plot=FALSE)
        hist2=hist(dur2$dur_warm[q,nach_warm],breaks=br,plot=FALSE)

        endx=max(which(hist1$counts!=0),which(hist2$counts!=0))+1
        print(endx)

        par(mfrow=c(1,2))
        par(plt=c(0.2,0.85,0.2,0.85))        
        plot(hist1$mids,hist1$density,xlim=c(0,endx*2+2),col=rgb(0,1,0,1/2),pch=15,ylim=c(0,max(c(hist1$density,hist2$density))),main=paste("warm periods in",season),xlab="duration of period")
        points(hist2$mids,hist2$density,col=rgb(1,0,0,1/2),pch=15)

        x=hist1$mids[1:endx]
        xy=data.frame(y=hist1$density[1:endx],x=x)
        fit=nls(y~exp(a+b*x),data=xy,start=list(a=0,b=0))
        yfit=exp(summary(fit)$parameters[1])*exp(x*summary(fit)$parameters[2])
        lines(x,yfit,col="green")
        x5=(log(0.05)-summary(fit)$parameters[1])/summary(fit)$parameters[2]
        abline(v=x5,col="green")

        x=hist2$mids[1:endx]
        xy=data.frame(y=hist2$density[1:endx],x=x)
        fit=nls(y~exp(a+b*x),data=xy,start=list(a=0,b=0))
        yfit=exp(summary(fit)$parameters[1])*exp(x*summary(fit)$parameters[2])
        lines(x,yfit,col="red")
        x5=(log(0.05)-summary(fit)$parameters[1])/summary(fit)$parameters[2]
        abline(v=x5,col="red")

        legend("topright", pch = c(15,15), col = c("green", "red"), 
                legend = c(paste("before",trenn),paste("after",trenn)))
        print(yfit)
        print(x)
        werwe      
        hist1=hist(dur2$dur_cold[q,vor_cold],breaks=br,plot=FALSE)
        hist2=hist(dur2$dur_cold[q,nach_cold],breaks=br,plot=FALSE)

        plot(hist1,xlim=c(0,max(which(hist1$counts!=0))*2+2),col=rgb(0,1,0,1/2),ylim=c(0,max(c(hist1$counts,hist2$counts))),main=paste("cold periods in",season),xlab="duration of period")
        plot(hist2,col=rgb(1,0,0,1/2),add=TRUE)
        legend("topright", pch = c(15,15), col = c("green", "red"), 
                legend = c(paste("before",trenn),paste("after",trenn)))
                    
    }

    graphics.off()
}





dat=dat_load("../data/dat_regional.nc")



ndays = c(91)
nyrs = c(5)

stations=c(488,510,604,744,920,887,251,98,270,281,169,164,353,121,11,39)
#stations=c(251)
for (nday in ndays){
    for (nyr in nyrs){
        for (qq in stations){
            trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend_r.nc",nday,nyr,nday,nyr))
            cat("loading persistence\n") 
            per=markov_load(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr))
            dur1=duration_load(sprintf("../data/%s_%s/%s_%s_duration_spring.nc",nday,nyr,nday,nyr))
            dur2=duration_load(sprintf("../data/%s_%s/%s_%s_duration_summer.nc",nday,nyr,nday,nyr))
            dur3=duration_load(sprintf("../data/%s_%s/%s_%s_duration_autumn.nc",nday,nyr,nday,nyr))
            dur4=duration_load(sprintf("../data/%s_%s/%s_%s_duration_winter.nc",nday,nyr,nday,nyr))

            station_plot(dat=dat,trend=trend,per=per,dur1=dur1,dur2=dur2,dur3=dur3,dur4=dur4,
                q=qq,start=365*(3)+1,stop=365*(58),filename=sprintf("../plots/station/%s_%s_station->%s.pdf",nday,nyr,qq))
        }
    }
}


