source("write_load.r")
source("functions_persistence.r")

duration_plot <- function(dat,trend,ind,station,start,stop,filename){
    tmp=per_duration(ind=as.vector(per$ind[station,,])[start:stop],time=as.vector(dat$time_2D)[start:stop])
    pdf(file=filename)
    plot(dat$time[start:stop],as.vector(dat$tas[station,,])[start:stop],ylim=c(-30,30))
    points(tmp$dur_warm_mid,tmp$dur_warm,pch=15,col="red")
    points(tmp$dur_cold_mid,-tmp$dur_cold,pch=15,col="blue")
    lines(dat$time[start:stop],as.vector(trend[station,,])[start:stop],col="green")
    #lines(dat$year,per$markov[station,1,],col="lightblue")
    xy=c(2003,2004,2005,2006,2007)
    xl=xy+0.41
    xr=xy+0.66
    rect(xl,-100,xr,100, col= "lightblue", density=10)
    graphics.off()
}

duration_hist <- function(per,dat,station,beg1,beg2,end1,end2){
    br<-seq(0,60,10)

    dur1=per_duration(ind=per$ind[station,151:243,3:31],time=dat$time_2D[151:243,beg1:end1])
    dur2=per_duration(ind=per$ind[station,151:243,30:58],time=dat$time_2D[151:243,beg2:end2])

    hist1=hist(dur1$dur_warm,breaks=br,plot=FALSE)
    hist2=hist(dur2$dur_warm,breaks=br,plot=FALSE)
    warmX_diff=sum(hist2$density[3:6])-sum(hist1$density[3:6])

    hist1=hist(dur1$dur_cold,breaks=br,plot=FALSE)
    hist2=hist(dur2$dur_cold,breaks=br,plot=FALSE)   
    coldX_diff=sum(hist2$density[3:6])-sum(hist1$density[3:6])

    warm_mean_diff=mean(dur2$dur_warm)-mean(dur1$dur_warm)
    cold_mean_diff=mean(dur2$dur_cold)-mean(dur1$dur_cold)

    if (3==3){
        pdf(file="../plots/histogramm.pdf")
        plot(hist1,col=rgb(0,0,1,1/4))
        plot(hist2,col=rgb(1,0,0,1/4),add=TRUE)
    }

    return(list(warm_mean_diff=warm_mean_diff,cold_mean_diff=cold_mean_diff,warmX_diff=warmX_diff,coldX_diff=coldX_diff))
}





ndays = c(61)
nyrs = c(5)

dat=dat_load("../data/mid_lat.nc")


for (nday in ndays){
    for (nyr in nyrs){
        trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend.nc",nday,nyr,nday,nyr))
        cat("loading persistence\n") 
        per=markov_load(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr))

        tmp=duration_hist(per,dat,station=510,23,41,40,58)
        print(tmp)
        #duration_plot(dat=dat,trend=trend,ind=per$ind,station=487,start=365*(53)+1,stop=365*(58),filename=sprintf("../plots/dur/%s_%s_duration.pdf",nday,nyr))

    }
}


