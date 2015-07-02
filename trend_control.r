source("region_average.r")
source("write.r")
source("load.r")

trend_control_warm_days <- function(dat,per,seasonStart=c(59,151,243,335),seasonStop=c(150,242,334,424)){
    ntot=length(dat$ID)
    waTrend=array(NA,dim=c(ntot,8))
    for (q in 1:ntot){
        for (sea in 1:4){
            warmeTage=array(NA,62)
            kalteTage=array(NA,62)
            for (i in 1:61){
                if (seasonStop[sea]>365){
                    warmeTage[i]=length(which(per$ind[q,seasonStart[sea]:365,i]==1))
                    kalteTage[i]=length(which(per$ind[q,seasonStart[sea]:365,i]==-1)) 

                    warmeTage[i]=warmeTage[i]+length(which(per$ind[q,1:(seasonStop[sea]-365),(i+1)]==1))
                    kalteTage[i]=kalteTage[i]+length(which(per$ind[q,1:(seasonStop[sea]-365),(i+1)]==-1)) 
                }
                else {
                    warmeTage[i]=length(which(per$ind[q,seasonStart[sea]:seasonStop[sea],i]==1))
                    kalteTage[i]=length(which(per$ind[q,seasonStart[sea]:seasonStop[sea],i]==-1))  
                }   
                if ((warmeTage[i]+kalteTage[i])<80){
                    warmeTage[i]=NA
                    kalteTage[i]=NA
                }            
            }
            if (length(which(is.na(warmeTage)))<30){ 
                lm.r=lm(warmeTage~dat$year)
                waTrend[q,sea]=summary(lm.r)$coefficients[2]
                waTrend[q,(4+sea)]=summary(lm.r)$coefficients[8]
            }
        }
    }
    write.table(waTrend,"../data/warmeTage_trends_4seasons.txt")
    return(waTrend[1:ntot,1:4])
}


trend_view_diff_trends <- function(dat,q,start,stop,ndays,nyrs){
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
#trend_view_diff_trends(dat,488,c(2000,2003,2003.3),c(2005,2004,2003.7),c(61,91,121),c(5,1,7))

trend_view_detail <-function(dat,trend,year,q){
	y=trend[488,,year]
    maxi=which(diff(y)==max(diff(y[45:320])))
    print(diff(y))
    p0=maxi
    print(p0)
    print(y[(maxi-2):(maxi+2)])
    print(diff(y)[(maxi-2):(maxi+2)])
    window1=dat$tas[488,(p0-45-1):(p0+45-1),(year-2):(year+2)]
    window2=dat$tas[488,(p0-45):(p0+45),(year-2):(year+2)]
    window3=dat$tas[488,(p0-45+1):(p0+45+1),(year-2):(year+2)]

    window_beg=dat$tas[488,(p0-45):(p0-45+2),(year-3):(year+2)]
    window_end=dat$tas[488,(p0+45):(p0+45+2),(year-3):(year+2)]

    print(window_beg)
    print(window_end)
    print(dat$tas[488,(p0-47):(p0+47),year])

    diff1=sum(window_beg[2,]-window_beg[1,])+sum(window_end[2,]-window_end[1,])
    diff2=sum(window_beg[3,]-window_beg[2,])+sum(window_end[3,]-window_end[2,])

    print(paste(diff1,diff1/(91*5)))
    print(paste(diff2,diff2/(91*5)))

    print(paste(mean(window1),mean(window2),mean(window3))) 

    pdf(file="../plots/station/trend_detail.pdf")
    plot(dat$time,dat$tas[488,,],xlim=c((1950+year+0.15),(year+1950+0.5)),ylim=c(-7,10),pch=20,cex=0.4)

    p=1950+year+p0/365-0.5*1/365
    for (i in 1:5){
    	text((p-(44)/365),window_beg[1,i],label=i,col="blue",cex=0.6)
    }
    for (i in 1:5){
    	text((p-(44-1)/365),window_beg[2,i],label=i,col="green",cex=0.6)
    }

    for (i in 1:5){
    	print(paste(i,window_end[1,i]))
    	text((p+(44)/365),window_end[1,i],label=i,col="blue",cex=0.6)
    }
    for (i in 1:5){
    	text((p+(44+1)/365),window_end[2,i],label=i,col="green",cex=0.6)
    }
    points(dat$time+1,dat$tas[488,,],pch=20,cex=0.2,col="green")
    points(dat$time+2,dat$tas[488,,],pch=20,cex=0.2,col="violet")

    points(dat$time,trend[q,,],pch=20,cex=0.2,col="red")	
    points(p,trend[q,p0,year],col="blue",pch=20,cex=0.3)
    points((p+1/365),trend[q,p0+1,year],col="green",pch=20,cex=0.3)


}










if (1==1){
    nday=91
    nyr=5


    trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend_r.nc",nday,nyr,nday,nyr))
    y=trend[488,45:320,53]


	dat=dat_load("../data/dat_regional.nc",reg=1)

	trend_view_detail(dat,trend,53,488)
	fsdf


    z=dat$tas[488,,53]

    mea_trend_diff=mean(abs(diff(y)))
    print(mea_trend_diff)
    print(mea_trend_diff*91*5/10)

    maxi=which(diff(y)==max(diff(y)))
    p=maxi+45
    print(maxi)
    print(p)
    print(max(diff(y)))

    window1=dat$tas[488,(p-45-1):(p+45-1),51:55]
    window2=dat$tas[488,(p-45):(p+45),51:55]
    window3=dat$tas[488,(p-45+1):(p+45+1),51:55]

    window_beg=dat$tas[488,(p-45-1):(p-45+1),51:55]
    window_end=dat$tas[488,(p+45-1):(p+45+1),51:55]

    print(window_beg)
    print(window_end)

    diff1=sum(window_beg[2,]-window_beg[1,])+sum(window_end[2,]-window_end[1,])
    diff2=sum(window_beg[3,]-window_beg[2,])+sum(window_end[3,]-window_end[2,])

    print(paste(diff1,diff1/(91*5)))
    print(paste(diff2,diff2/(91*5)))

    print(paste(mean(window1),mean(window2),mean(window3)))
    print(trend[488,(p-1):(p+1),53])

    trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend_r.nc",nday,nyr,nday,nyr))
    y=trend[488,,]
    print(max(abs(diff(y)),na.rm=TRUE))


    fsdfsd

    per=markov_load(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr))
    #watrends=trend_control_warm_days(dat,per)
    watrends=read.table("../data/warmeTage_trends_4seasons.txt")
    map_regional(dat,watrends[,1:4],c("spring warm day increase","summer warm day increase","autumn warm day increase","winter warm day increase"))
}