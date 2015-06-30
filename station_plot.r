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

        par(mfrow=c(1,2))

        par(plt=c(0.2,0.85,0.2,0.85))        
        plot(hist1,xlim=c(0,max(which(hist1$counts!=0))*2+2),col=rgb(0,1,0,1/2),ylim=c(0,max(c(hist1$counts,hist2$counts))),main=paste("warm periods in",season),xlab="duration of period")
        plot(hist2,col=rgb(1,0,0,1/2),add=TRUE)
        legend("topright", pch = c(15,15), col = c("green", "red"), 
                legend = c(paste("before",trenn),paste("after",trenn)))
        
        hist1=hist(dur2$dur_cold[q,vor_cold],breaks=br,plot=FALSE)
        hist2=hist(dur2$dur_cold[q,nach_cold],breaks=br,plot=FALSE)

        plot(hist1,xlim=c(0,max(which(hist1$counts!=0))*2+2),col=rgb(0,1,0,1/2),ylim=c(0,max(c(hist1$counts,hist2$counts))),main=paste("cold periods in",season),xlab="duration of period")
        plot(hist2,col=rgb(1,0,0,1/2),add=TRUE)
        legend("topright", pch = c(15,15), col = c("green", "red"), 
                legend = c(paste("before",trenn),paste("after",trenn)))
                    
    }

    graphics.off()
}

trend_view_plot <- function(dat,q,start,stop,ndays,nyrs){
    qq=which(dat$ID==q)
    q=qq[1]
    pdf(file=sprintf("../plots/station/trend_%s.pdf",dat$ID[q]))
    for (k in 1:3){
        plot(dat$time,dat$tas[q,,],pch=20,cex=0.2,xlim=c(start[k],stop[k]),ylab="temperature anomaly in deg C",xlab="",main="different trends")
        linestyle=c()
        color=c()
        label=c("days years")
        for (i in 1:length(ndays)){
            for (j in 1:length(nyrs)){
                trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend.nc",ndays[i],nyrs[j],ndays[i],nyrs[j]))
                linestyle[(i-1)*length(nyrs)+j+1]=j
                color[(i-1)*length(nyrs)+j+1]=rgb(i/3,0,0.5)
                label[(i-1)*length(nyrs)+j+1]=paste(ndays[i],nyrs[j],sep="     ")
                lines(dat$time,trend[q,,],col=rgb(i/3,0,0.5),lty=j)
            }
        }
        rect(2003.5-30/365,-30,2003.5+30/365,30,col=rgb(1/3,0,0.5),density=0.0)
        #text(2003.5+30/365,5,labels=61,col=rgb(1/3,0,0.5))
        rect(2003.5-45/365,-30,2003.5+45/365,30,col=rgb(2/3,0,0.5),density=0.0)
        #text(2003.5+45/365,13.5,labels=91,col=rgb(2/3,0,0.5))
        rect(2003.5-60/365,-30,2003.5+60/365,30,col=rgb(3/3,0,0.5),density=0.0)
        #text(2003.5+60/365,15,labels=121,col=rgb(3/3,0,0.5))

        trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend_r.nc",91,5,91,5))
        lines(dat$time,trend[qq[1],,],col=rgb(0,1,0),lty=2)
        linestyle[(i-1)*length(nyrs)+j+2]=2
        color[(i-1)*length(nyrs)+j+2]=rgb(0,1,0)
        label[(i-1)*length(nyrs)+j+2]=paste("91 5 control")
        
        legend("bottomleft",lty=linestyle, col = color, legend=label)
    }
    graphics.off()

}



dat=dat_load("../data/mid_lat.nc")
#trend_view_plot(dat,488,c(2000,2003,2003.3),c(2005,2004,2003.7),c(61,91,121),c(5,1,7))



ndays = c(91)
nyrs = c(5)

stations=c(488,510,604,744,920,887,251,98,270,281,169,164,353,121,11,39)
#stations=c(251)
for (nday in ndays){
    for (nyr in nyrs){
        for (qq in stations){
            trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend.nc",nday,nyr,nday,nyr))
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


