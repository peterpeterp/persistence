source("write.r")
source("load.r")
source("functions_duration.r")

duration_plot <- function(dat,trend,per,station,start,stop,filename){
    tmp=per_duration(ind=as.vector(per$ind[station,,])[start:stop],time=as.vector(dat$time_2D)[start:stop])
    pdf(file=filename)
    #plot(dat$time[start:stop],as.vector(dat$tas[station,,])[start:stop],ylim=c(-30,30))
    plot(NA,xlim=c(1950,2011),ylim=c(-30,30))
    points(tmp$dur_warm_mid,tmp$dur_warm,pch=15,col="red")
    points(tmp$dur_cold_mid,-tmp$dur_cold,pch=15,col="blue")
    lines(dat$time[start:stop],as.vector(trend[station,,])[start:stop],col="black")
    mark_sum=per$markov[station,1,]*30
    lines(dat$year,mark_sum,col="green")
    xy=seq(1950,2011,1)   
    xl=xy+0.41
    xr=xy+0.66
    rect(xl,-100,xr,100, col= "lightblue", density=10)
    graphics.off()
}







ndays = c(91)
nyrs = c(5)

dat=dat_load("../data/mid_lat.nc")

stations=c(1249,1233,1234,1214,1248)
stations=c(482,489,499)
for (nday in ndays){
    for (nyr in nyrs){
        for (qq in stations){
            q=which(dat$ID==qq)
            trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend.nc",nday,nyr,nday,nyr))
            cat("loading persistence\n") 
            per=markov_load(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr))
            print(per$ind[q,1:365,54:58])
            print(dat$tas[q,1:365,54:58])
            print(trend[q,1:365,54:58])
            sdfsf
            duration_plot(dat=dat,trend=trend,per=per,station=q,start=365*(3)+1,stop=365*(58),filename=sprintf("../plots/dur/%s_%s_duration_%s.pdf",nday,nyr,qq))
        }
    }
}


