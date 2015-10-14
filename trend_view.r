source("write.r")
source("load.r")

trend_control_warm_days <- function(dat,ind,seasonStart=c(60,151,242,334,1),seasonStop=c(150,241,333,424,365),
    filename="../data/warmeTage_trends_5seasons_1950-2014.txt"){
    ntot=length(dat$ID)
    waTrend=array(NA,dim=c(ntot,15))
    for (q in 1:ntot){
        for (sea in 1:length(seasonStart)){
            warmeTage=array(NA,65)
            kalteTage=array(NA,65)
            for (i in 1:61){
                if (seasonStop[sea]>365){
                    warmeTage[i]=length(which(ind[q,seasonStart[sea]:365,i]==1))
                    kalteTage[i]=length(which(ind[q,seasonStart[sea]:365,i]==-1)) 

                    warmeTage[i]=warmeTage[i]+length(which(ind[q,1:(seasonStop[sea]-365),(i+1)]==1))
                    kalteTage[i]=kalteTage[i]+length(which(ind[q,1:(seasonStop[sea]-365),(i+1)]==-1)) 
                }
                else {
                    warmeTage[i]=length(which(ind[q,seasonStart[sea]:seasonStop[sea],i]==1))
                    kalteTage[i]=length(which(ind[q,seasonStart[sea]:seasonStop[sea],i]==-1))  
                }   
                if ((warmeTage[i]+kalteTage[i])<(seasonStop[sea]-seasonStart[sea]-10)){
                    warmeTage[i]=NA
                    kalteTage[i]=NA
                }            
            }
            if (length(which(is.na(warmeTage)))<30){ 
                lm.r=lm(warmeTage~dat$year)
                waTrend[q,sea]=summary(lm.r)$coefficients[2]
                waTrend[q,(5+sea)]=summary(lm.r)$coefficients[4]
                waTrend[q,(10+sea)]=sum(warmeTage,na.rm=TRUE)/sum(c(warmeTage,kalteTage),na.rm=TRUE)
            }
        }
        cat(waTrend[q,15])
        cat("__")
    }
    write.table(waTrend,filename)
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

trend_3states <- function(){
    q=238
    trend=trend_load("../data/91_5/91_5_trend.nc")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")

    detrended=dat$tas[q,,]-trend[q,,]
    threshold=sd(detrended,na.rm=TRUE)*0.5


    pdf(file="../plots/488")
    plot(dat$time,dat$tas[q,,],xlim=c(2000,2014),pch=20,cex=0.5,main="3 states",ylab="temperature anomaly in deg C",xlab="")
    lines(dat$time,trend[q,,],col="green")
    lines(dat$time,trend[q,,]+threshold,col="grey")
    lines(dat$time,trend[q,,]-threshold,col="grey")
    warm=dat$tas[q,,]
    warm[warm<trend[q,,]+threshold]=NA
    points(dat$time,warm,pch=20,cex=0.5,col="red")
    cold=dat$tas[q,,]
    cold[cold>trend[q,,]-threshold]=NA
    points(dat$time,cold,pch=20,cex=0.5,col="blue")
    legend("bottomright",pch=c(NA,20,20,20),col=c(NA,"red","black","blue"),legend=c("diff between grey lines = 1 sd","warm day","average day","cold day"))
}

method_explanation <- function(yearReal=2003,todo=c(1,1,1,1),name="full"){
    year=yearReal-1949
    q=488
    source("load.r")
    trend=trend_load("../data/91_5/91_5_trend.nc")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    nc=open.ncdf("../data/91_5/2_states/duration/91_5_duration_2s_year.nc")
    dur=get.var.ncdf(nc,"dur")
    dur_mid=get.var.ncdf(nc,"dur_mid")

    nc=open.ncdf("../data/91_5/2_states/markov/91_5_markov_2states.nc")
    mar=get.var.ncdf(nc,"markov")



    pdf(file=paste("../plots/explanatory_",name,".pdf",sep=""),width=8,height=7)
    par(mar=c(4.2,4,1,1))
    plot(NA,xlim=c(yearReal+59/365,yearReal+1.5+59/365),ylim=c(-22,25),pch=20,cex=0.5,
        main="",ylab="temperature anomaly in deg C",xlab="",frame.plot=FALSE)


    if (todo[3]==1){
        waIn=which(dur_mid[q,2,]>yearReal & dur_mid[q,2,]<yearReal+1.3)
        mi=13
        maxi=max(dur[q,2,waIn])
        for (p in waIn){
            oben=mi+2.5#+(-1)^p*1.25
            unten=mi#+(-1)^p*1.25    
            beg=dur_mid[q,2,p]-1/2*dur[q,2,p]/365
            end=dur_mid[q,2,p]+1/2*dur[q,2,p]/365
            polygon(x=c(beg,beg,end,end),y=c(unten,oben,oben,unten),border=NA,col=rgb(dur[q,2,p]/maxi/2+0.5,dur[q,2,p]/maxi/4+0.3,dur[q,2,p]/maxi/4+0.3))
            #abline(v=beg,col=rgb(0.9,0.9,0.9))
        }
        text(yearReal+1.5,mi+1.25,label="duration of \n warm periods")

        kaIn=which(dur_mid[q,1,]>yearReal & dur_mid[q,1,]<yearReal+1.3)
        mi=-13
        maxi=max(dur[q,1,waIn])
        for (p in kaIn){
            oben=mi+2.5#+(-1)^p*1.25
            unten=mi#+(-1)^p*1.25    
            beg=dur_mid[q,1,p]-1/2*dur[q,1,p]/365
            end=dur_mid[q,1,p]+1/2*dur[q,1,p]/365
            polygon(x=c(beg,beg,end,end),y=c(unten,oben,oben,unten),border=NA,col=rgb(dur[q,1,p]/maxi/4+0.3,dur[q,1,p]/maxi/4+0.3,dur[q,1,p]/maxi/2+0.5))
            #abline(v=beg,col=rgb(0.9,0.9,0.9))
        }
        text(yearReal+1.5,mi+1.25,label="duration of \n cold periods")
    }

    if (todo[4]==1){
        seasonStart=c(59,151,243,334,424)
        season_names=c("spring","summer","autumn","winter")
        for (sea in 1:4){
            beg=yearReal+seasonStart[sea]/365
            end=yearReal+seasonStart[sea+1]/365

            text((beg+45/365),24,label=season_names[sea])

            polygon(x=c(beg,beg,end,end),y=c(18,22,22,18),col=rgb(0.9,0.6,0.6))
            text((beg+45/365),20,label=paste(sprintf("%.02f",(mar[q,sea,4,year]*100)),"%",sep=""))

            polygon(x=c(beg,beg,end,end),y=c(-16,-20,-20,-16),col=rgb(0.6,0.6,0.9))
            text((beg+45/365),-18,label=paste(sprintf("%.02f",(mar[q,sea,1,year]*100)),"%",sep=""))
        }
        text(yearReal+1.5,20,label="transition probability \n warm to warm")
        #arrows(yearReal+1.25,22,yearReal+1.15,22,code=2)
        text(yearReal+1.5,-18,label="transition probability \n cold to cold")
        #arrows(yearReal+1.25,-18,yearReal+1.15,-18,code=2)
    }

    points(dat$time,dat$tas[q,,],cex=0.5)
    if (todo[1]==1){
       lines(dat$time,trend[q,,],col="green") 
    }

    if (todo[2]==1){  
        time=seq(yearReal,(yearReal+2),length=730)
        time=time+0.5*1/365
        warm=dat$tas[q,,year:(year+1)]
        warm[warm<trend[q,,year:(year+1)]]=NA
        points(time,warm,pch=20,cex=0.5,col="red")
        cold=dat$tas[q,,year:(year+1)]
        cold[cold>trend[q,,year:(year+1)]]=NA
        points(time,cold,pch=20,cex=0.5,col="blue")
    }
    #legend("bottomleft",pch=c(20,20),col=c("red","blue"),legend=c("warm day","cold day"))
    graphics.off()
}

if (1==2){
    method_explanation(2003,c(2,2,2,2),"1")
    method_explanation(2003,c(1,2,2,2),"2")
    method_explanation(2003,c(1,1,2,2),"3")
    method_explanation(2003,c(1,1,1,2),"4")
    method_explanation(2003,c(1,1,1,1),"5")
}

asymmetry_analysis <- function(q=459){
    trend=trend_load("../data/91_5/91_5_trend.nc")
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    y=dat$tas[q,,]-trend[q,,]
    mittel=mean(y,na.rm=TRUE)
    std=sd(y,na.rm=TRUE)

    sd1h=length(which(y>(1*std)))
    sd2h=length(which(y>(2*std)))
    sd3h=length(which(y>(3*std)))

    sd1t=length(which(y<(-1*std)))
    sd2t=length(which(y<(-2*std)))
    sd3t=length(which(y<(-3*std)))

    print(sd1h)
    print(sd2h)
    print(sd3h)
    print(sd1t)
    print(sd2t)
    print(sd3t)

    pdf(file="../plots/zwischenzeugs/asymmetry.pdf")
    plot(dat$time,y,xlim=c(2000,2010))

}

trend_mean_median <- function(q,start,stop){
    trendID="91_5"
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    trend_mean=trend_load(paste("../data/",trendID,"/",trendID,"_trend.nc",sep=""))
    trend_median=trend_load(paste("../data/",trendID,"/",trendID,"_trend_median.nc",sep=""))
    pdf(file=sprintf("../plots/91_5/sonstiges/trend_91_5_mean_median_station%s.pdf",dat$ID[q]))

    selection=which(dat$time > start & dat$time < stop)
    plot(dat$time[selection],as.vector(dat$tas[q,,])[selection])
    lines(dat$time,trend_mean[q,,],col="blue")
    lines(dat$time,trend_median[q,,],col="green")

    y=as.vector(dat$tas[q,,])[which(!is.na(dat$tas[q,,]))]
    between1=which(y>trend_mean[q,,] & y<trend_median[q,,])
    points(dat$time[between1],y[between1],pch=4,col="red")

    between2=which(y<trend_mean[q,,] & y>trend_median[q,,])
    points(dat$time[between2],y[between2],pch=4,col="orange")
        
    legend("bottomright",legend=paste((length(between1)+length(between2))/length(y)*100,"% of days between the trends",sep=""))

    graphics.off()

}

trend_mean_median(488,2000,2010)
trend_mean_median(351,2000,2010)
trend_mean_median(1011,2000,2010)