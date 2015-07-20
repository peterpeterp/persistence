source("write.r")
source("load.r")
source("functions_duration.r")

station_plot <- function(dat,trend,per,dur_name,warmeTageIncr,q,filename){
    pdf(file=filename,paper='a4')

    station=which(dat$ID==q)[1]

    if (dim(per$markov)[3]==4){
        state_names=c("cold","warm")
        color=c("blue","red")
    }

    if (dim(per$markov)[3]==9){
        state_names=c("cold","normal","warm")
        color=c("blue","violet","red")
    }

    season_names=c("spring","summer","autumn","winter")
    season_starts=c(59,151,242,335,1)
    season_stops=c(151,242,335,425,1)
    for (sea in 1:4){
        par(mfrow=c(length(state_names),1))
        par(mar=c(2,2, 2, 0.5))


        nc=open.ncdf(paste(dur_name,season_names[sea],".nc",sep="")) 
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")

        for (state in 1:length(state_names)){
            plot(NA,xlim=c(1950,2014),xaxp=c(1950,2015,13),ylim=c(0,50),xlab="",ylab="duration of persistence")#, main="summer warm persistence")
            points(dur_mid[station,state,],dur[station,state,],pch=4,col=color[state])

            markov=per$markov[station,sea,(state*state),]*30
            lines(dat$year+0.5,markov,col="green")

            lines(dat$time,as.vector(trend[station,,]),col="black")

            F=ecdf(dur[station,state,])
            perc95=which(F(dur[station,state,])>0.95)
            abline(lm(dur[station,state,perc95]~dur_mid[station,state,perc95]),col="red")

            for (i in seq(1950,2015,5)){
                abline(v=i,col="gray",lty="dotted")
            }
            legend("topleft", pch = c(NA,4,NA,NA,NA),lty=c(NA,NA,1,NA,1), col = c(NA,color[state], "green",NA,"black"), 
                    legend = c(season_names[sea],paste(state_names[state],"days in a row"),paste("markov",state_names[state],"persistence"),
                        paste("warm day increase per year",sprintf("%.3f",warmeTageIncr[q,sea])),"trend"))

        }


        par(new=TRUE,plt=c(0.76,0.96,0.67,0.87))
        q=station
        library(rworldmap)
        worldmap = getMap(resolution = "low")
        plot(worldmap,xlim=c(dat$lon[q]-10,dat$lon[q]+10),ylim=c(dat$lat[q]-10,dat$lat[q]+10))
        points(dat$lon[q],dat$lat[q],pch=15,col="red")

        if (1==1){
            trenn=1980
            br=seq(0,100,2)
            
            par(mfrow=c(1,length(state_names)))
            par(plt=c(0.2,0.85,0.2,0.85)) 

            for (state in 1:length(state_names)){

                vor=which(dur_mid[q,state,]<trenn)
                nach=which(dur_mid[q,state,]>trenn)

                hist1=hist(dur[q,state,vor],breaks=br,plot=FALSE)
                hist2=hist(dur[q,state,nach],breaks=br,plot=FALSE)


 

                endx=max(which(hist1$counts!=0),which(hist2$counts!=0))+1

                plot(hist1$mids,hist1$density,xlim=c(0,endx*2+2),col=rgb(0,1,0,1/2),pch=15,ylim=c(0,max(c(hist1$density,hist2$density))),main=paste(state_names[state],"periods in",season_names[sea]),xlab="duration of period")
                points(hist2$mids,hist2$density,col=rgb(1,0,0,1/2),pch=15)

                x=hist1$mids
                xy=data.frame(y=hist1$density,x=x)
                fit=nls(y~exp(a+b*x),data=xy,start=list(a=0,b=0))
                yfit=exp(summary(fit)$parameters[1])*exp(x*summary(fit)$parameters[2])
                lines(x,yfit,col="green")
                x2=(log(0.02)-summary(fit)$parameters[1])/summary(fit)$parameters[2]
                abline(v=x2,col="green")

                x=hist2$mids
                xy=data.frame(y=hist2$density,x=x)
                fit=nls(y~exp(a+b*x),data=xy,start=list(a=0,b=0))
                yfit=exp(summary(fit)$parameters[1])*exp(x*summary(fit)$parameters[2])
                lines(x,yfit,col="red")
                x2=(log(0.02)-summary(fit)$parameters[1])/summary(fit)$parameters[2]
                abline(v=x2,col="red")

                legend("topright", pch = c(15,15), col = c("green", "red"), 
                        legend = c(paste("before",trenn),paste("after",trenn)))
            
            }
                        
        }
    }
    graphics.off()
}



dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")



ndays = c(91)
nyrs = c(5)

stations=c(488,510,604,744,920,887,251,98,270,281,169,164,353,121,11,39)
stations=c(488)
for (nday in ndays){
    for (nyr in nyrs){
        for (qq in stations){
            trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend.nc",nday,nyr,nday,nyr))
            cat("loading persistence\n") 
            per=markov_load(sprintf("../data/%s_%s/%s_%s_markov2s.nc",nday,nyr,nday,nyr),4)
            #dur=duration_load(sprintf("../data/%s_%s/%s_%s_duration_",nday,nyr,nday,nyr))
            warmeTageIncr=read.table("../data/warmeTage_trends_5seasons.txt")


            station_plot(dat=dat,trend=trend,per=per,dur_name=sprintf("../data/%s_%s/%s_%s_duration_2s_",nday,nyr,nday,nyr),
                warmeTageIncr=warmeTageIncr,q=qq,
                filename=sprintf("../plots/station/%s_%s_a4_2s_%s.pdf",nday,nyr,qq))
        }
    }
}


