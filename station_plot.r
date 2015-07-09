source("write.r")
source("load.r")
source("functions_duration.r")

station_plot <- function(dat,trend,per,dur,q,start,stop,filename){
    pdf(file=filename)

    station=which(dat$ID==q)[1]

    season_names=c("spring","summer","autumn","winter","year")
    season_starts=c(59,151,242,335,1)
    season_stops=c(151,242,335,425,1)
    for (sea in 1:5){
        par(mfrow=c(2,1))
        par(mar=c(2,2, 2, 0.5)) 
        if (1==2){
            warmeTage=array(NA,62)
            kalteTage=array(NA,62)
            for (i in 1:62){
                if (length(!is.na(per$ind[station,season[1]:season[2],i]))>(season[2]-season[1]-10)){
                    warmeTage[i]=length(which(per$ind[station,season[1]:season[2],i]==1))
                    kalteTage[i]=length(which(per$ind[station,season[1]:season[2],i]==-1)) 
                }
            }
        }

        plot(NA,xlim=c(1950,2010),xaxp=c(1950,2010,12),ylim=c(0,50),xlab="",ylab="duration of persistence")#, main="summer warm persistence")
        points(dur$dur_warm_mid[station,],dur$dur_warm[station,],pch=4,col="red")
        markov=per$markov[station,sea,4,]*30
        lines(dat$year+0.5,markov,col="green")
        #abline(lm(warmeTage~dat$year),col="yellow")
        for (i in seq(1950,2010,5)){
            abline(v=i,col="gray",lty="dotted")
        }
        lines(dat$time[start:stop],as.vector(trend[station,,])[start:stop],col="black")
        #legend("topleft", pch = c(NA,4,NA,NA,NA),lty=c(NA,NA,1,1,1), col = c(NA,"red", "green","yellow","black"), 
        #        legend = c(season_names[sea],"warm days in a row","markov cold persistence",
        #            sprintf("warm day increase per year %.3f",summary(lm(warmeTage~dat$year))$coefficients[2]),"trend"))

        plot(NA,xlim=c(1950,2011),xaxp=c(1950,2010,12),ylim=c(0,50),xlab="",ylab="duration of persistence")#, main="summer cold persistence")
        points(dur$dur_cold_mid[station,],dur$dur_cold[station,],pch=5,col="blue")
        markov=per$markov[station,sea,1,]*30
        lines(dat$year,markov,col="green")
        markov=per$markov[station,sea,2,]*30
        lines(dat$year,markov,col="violet")
        markov=per$markov[station,sea,3,]*30
        lines(dat$year,markov,col="red")
        #abline(lm(kalteTage~dat$year),col="violet")
        lines(dat$time[start:stop],as.vector(trend[station,,])[start:stop],col="black")
        for (i in seq(1950,2010,5)){
            abline(v=i,col="gray",lty="dotted")
        }
        #legend("topleft", pch = c(NA,5,NA,NA,NA),lty=c(NA,NA,1,1,1), col = c(NA,"blue", "green","violet","black"), 
        #        legend = c(season_names[sea],"cold days in a row","markov cold persistence",
        #            sprintf("warm day increase per year %.3f",summary(lm(kalteTage~dat$year))$coefficients[2]),"trend"))

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
            vor_warm=which(dur$dur_warm_mid[q,]<trenn)
            vor_cold=which(dur$dur_cold_mid[q,]<trenn)
            nach_warm=which(dur$dur_warm_mid[q,]>trenn)
            nach_cold=which(dur$dur_cold_mid[q,]>trenn)

            hist1=hist(dur$dur_warm[q,vor_warm],breaks=br,plot=FALSE)
            hist2=hist(dur$dur_warm[q,nach_warm],breaks=br,plot=FALSE)


            par(mfrow=c(1,2))
            par(plt=c(0.2,0.85,0.2,0.85))  

            endx=max(which(hist1$counts!=0),which(hist2$counts!=0))+1

            plot(hist1$mids,hist1$density,xlim=c(0,endx*2+2),col=rgb(0,1,0,1/2),pch=15,ylim=c(0,max(c(hist1$density,hist2$density))),main=paste("warm periods in",season_names[sea]),xlab="duration of period")
            points(hist2$mids,hist2$density,col=rgb(1,0,0,1/2),pch=15)

            x=hist1$mids
            xy=data.frame(y=hist1$density,x=x)
            fit=nls(y~exp(a+b*x),data=xy,start=list(a=0,b=0))
            yfit=exp(summary(fit)$parameters[1])*exp(x*summary(fit)$parameters[2])
            lines(x,yfit,col="green")
            x5=(log(0.05)-summary(fit)$parameters[1])/summary(fit)$parameters[2]
            abline(v=x5,col="green")

            x=hist2$mids
            xy=data.frame(y=hist2$density,x=x)
            fit=nls(y~exp(a+b*x),data=xy,start=list(a=0,b=0))
            yfit=exp(summary(fit)$parameters[1])*exp(x*summary(fit)$parameters[2])
            lines(x,yfit,col="red")
            x5=(log(0.05)-summary(fit)$parameters[1])/summary(fit)$parameters[2]
            abline(v=x5,col="red")

            legend("topright", pch = c(15,15), col = c("green", "red"), 
                    legend = c(paste("before",trenn),paste("after",trenn)))
        
            hist1=hist(dur$dur_cold[q,vor_cold],breaks=br,plot=FALSE)
            hist2=hist(dur$dur_cold[q,nach_cold],breaks=br,plot=FALSE)

            endx=max(which(hist1$counts!=0),which(hist2$counts!=0))+1
          
            plot(hist1$mids,hist1$density,xlim=c(0,endx*2+2),col=rgb(0,1,0,1/2),pch=15,ylim=c(0,max(c(hist1$density,hist2$density))),main=paste("cold periods in",season_names[sea]),xlab="duration of period")
            points(hist2$mids,hist2$density,col=rgb(1,0,0,1/2),pch=15)

            x=hist1$mids
            xy=data.frame(y=hist1$density,x=x)
            fit=nls(y~exp(a+b*x),data=xy,start=list(a=0,b=0))
            yfit=exp(summary(fit)$parameters[1])*exp(x*summary(fit)$parameters[2])
            lines(x,yfit,col="green")
            x5=(log(0.05)-summary(fit)$parameters[1])/summary(fit)$parameters[2]
            abline(v=x5,col="green")

            x=hist2$mids
            xy=data.frame(y=hist2$density,x=x)
            fit=nls(y~exp(a+b*x),data=xy,start=list(a=0,b=0))
            yfit=exp(summary(fit)$parameters[1])*exp(x*summary(fit)$parameters[2])
            lines(x,yfit,col="red")
            x5=(log(0.05)-summary(fit)$parameters[1])/summary(fit)$parameters[2]
            abline(v=x5,col="red")

            legend("topright", pch = c(15,15), col = c("green", "red"), 
                    legend = c(paste("before",trenn),paste("after",trenn)))
                        
        }
    }
    graphics.off()
}





dat=dat_load("../data/dat_regional.nc")



ndays = c(91)
nyrs = c(5)

stations=c(488,510,604,744,920,887,251,98,270,281,169,164,353,121,11,39)
stations=c(510)
for (nday in ndays){
    for (nyr in nyrs){
        for (qq in stations){
            trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend_r.nc",nday,nyr,nday,nyr))
            cat("loading persistence\n") 
            per=markov_load(sprintf("../data/%s_%s/%s_%s_markov2s.nc",nday,nyr,nday,nyr),4)
            dur=duration_load(sprintf("../data/%s_%s/%s_%s_duration_summer.nc",nday,nyr,nday,nyr))
            #dur1=duration_load(sprintf("../data/%s_%s/%s_%s_duration_summer.nc",nday,nyr,nday,nyr))

            station_plot(dat=dat,trend=trend,per=per,dur=dur,
                q=qq,start=365*(3)+1,stop=365*(58),
                filename=sprintf("../plots/station/%s_%s_station_test_%s.pdf",nday,nyr,qq))
        }
    }
}


