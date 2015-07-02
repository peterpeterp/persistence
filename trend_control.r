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
                waTrend[q,(4+sea)]=summary(lm.r)$coefficients[4]
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
	year_index=year-1950
	y=trend[488,,year_index]
    maxi=which(diff(y)==max(diff(y[45:320])))

    d1=maxi
    d2=maxi+1

    d1_y=year+d1/365 - 0.5/365
    d2_y=year+d2/365 - 0.5/365

    pdf(file="../plots/station/trend_detail.pdf")
    plot(dat$time,dat$tas[q,,],xlim=c((year+0.15),(year+0.5)),ylim=c(-7,10),pch=20,cex=0.4,ylab="temperature anomalie in deg C",main="trend ???")

    for (i in -2:2){
    	text((d1_y-45/365),dat$tas[q,(d1-45),(year_index+(i+1))],label=i,col="blue",cex=0.4)
    }
    for (i in -2:2){
    	text((d2_y-45/365),dat$tas[q,(d2-45),(year_index+(i+1))],label=i,col="green",cex=0.4)
    }

    for (i in -2:2){
    	text((d1_y+45/365),dat$tas[q,(d1+45),(year_index+(i+1))],label=i,col="blue",cex=0.4)
    }
    for (i in -2:2){
    	text((d2_y+45/365),dat$tas[q,(d2+45),(year_index+(i+1))],label=i,col="green",cex=0.4)
    }



    points((d1_y-45/365),mean(dat$tas[q,(d1-45),(year_index-1):(year_index+3)]),pch=15,col="blue",cex=0.6)
    points((d2_y-45/365),mean(dat$tas[q,(d2-45),(year_index-1):(year_index+3)]),pch=15,col="green",cex=0.6)

    points((d1_y+45/365),mean(dat$tas[q,(d1+45),(year_index-1):(year_index+3)]),pch=15,col="blue",cex=0.6)
    points((d2_y+45/365),mean(dat$tas[q,(d2+45),(year_index-1):(year_index+3)]),pch=15,col="green",cex=0.6)

    x=seq((year+1/365-0.5/365),(year+365/365-0.5/365),1/365)	
	y=rowMeans(dat$tas[q,1:365,(year_index-1):(year_index+3)],dims=1)
    points(x,y,pch=15,col="violet",cex=0.3)

    points(d1_y,mean(y[(d1-45):(d1+45)]),col="blue")
    points(d2_y,mean(y[(d2-45):(d2+45)]),col="green")

    points(dat$time,trend[q,,],pch=20,cex=0.2,col="red")	
    points(d1_y,trend[q,maxi,year_index],col="blue",pch=20,cex=0.3)
    points(d2_y,trend[q,maxi+1,year_index],col="green",pch=20,cex=0.3)

    print(y)

    YearMeanDiff=mean(dat$tas[q,(d2-45),(year_index-1):(year_index+3)])-mean(dat$tas[q,(d1-45),(year_index-1):(year_index+3)])+mean(dat$tas[q,(d2+45),(year_index-1):(year_index+3)])-mean(dat$tas[q,(d1+45),(year_index-1):(year_index+3)])
    print(YearMeanDiff)
    print(mean(dat$tas[q,(d2-45),(year_index-1):(year_index+3)]))
    print(mean(dat$tas[q,(d1-45),(year_index-1):(year_index+3)]))
    print(mean(dat$tas[q,(d2+45),(year_index-1):(year_index+3)]))
    print(mean(dat$tas[q,(d1+45),(year_index-1):(year_index+3)]))

    trenddiff=mean(y[(d2-45):(d2+45)])-mean(y[(d1-45):(d1+45)])
    print(trenddiff)

    print(paste((d1-45),(d1+45),(d2-45),(d2+45)))
    print(paste(d1,d2))

    su1=sum(y[(d1-45):(d1+45)])
    su2=sum(y[(d2-45):(d2+45)])
    print(su1)
    print(su2)
    print((su2-su1)/91)

    legend("bottomright",pch=c(20,20,15,NA,NA),col=c("black","red","violet",NA,NA),
        legend=c("temp","trend","5 year mean",
            paste("sum 5 year mean blue",su1),
            paste("sum 5 year mean green",su2),
            paste("diff between trend points",trenddiff )))
}










if (1==1){
    nday=91
    nyr=5


    trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend_r.nc",nday,nyr,nday,nyr))
	dat=dat_load("../data/dat_regional.nc",reg=1)

	#trend_view_detail(dat,trend,2003,488)
    #trend_control_warm_days(dat,per)

    per=markov_load(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr))
    #watrends=trend_control_warm_days(dat,per)
    watrends=read.table("../data/warmeTage_trends_4seasons.txt")
    map_regional(dat,watrends[,1:4],c("spring warm day increase","summer warm day increase","autumn warm day increase","winter warm day increase"))
}