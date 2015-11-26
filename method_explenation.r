source("load.r")

method_explanation <- function(yearReal=2003,todo=c(1,1,1,1),name="full"){
    year=yearReal-1949
    q=745
    trend=trend_load("../data/91_5/_TMean/91_5_trend_mean_TMean.nc")
    dat=dat_load("../data/HadGHCND_TMean_data3D.day1-365.1950-2014.nc")
    nc=open.ncdf("../data/91_5/_TMean/91_5_mean_TMean_seasonal_median.nc")
    seasonal_median=get.var.ncdf(nc,"_seasonal_median")
    nc=open.ncdf("../data/91_5/_TMean/91_5_mean_TMean_state_ind_seasonal_median.nc")
    ind=get.var.ncdf(nc,"ind")
    nc=open.ncdf("../data/91_5/_TMean/duration/91_5_TMean_duration_4seasons.nc")
    dur=get.var.ncdf(nc,"dur")
    dur_mid=get.var.ncdf(nc,"dur_mid")

    pdf(file=paste("../plots/explanatory_",name,".pdf",sep=""),width=6,height=5)
    par(mar=c(4.2,4,1,1))
    plot(NA,xlim=c(yearReal+59/365,yearReal+1.48+59/365),ylim=c(-20,17),pch=20,cex=0.5,
        main="",ylab="temperature anomaly in deg C",xlab="",frame.plot=FALSE)

    detrended=dat$tas-trend
    
    if (todo[1]==0){
       points(dat$time,dat$tas[q,,],cex=0.5)
       points(dat$time,dat$tas[q,,]-trend[q,,],cex=0.5,col="green",pch=4)
       lines(dat$time,trend[q,,],col="red")
    }

    if (todo[1]==1){
       points(dat$time,detrended[q,,],cex=0.5)
    }

    if (todo[2]==1){
        sea_limit=c(60,152,244,335)
        season_names=c("spring","summer","autumn","winter")
        x=(1:365)/365+yearReal
        lines(x,seasonal_median[q,],col="green")

        x=(1:365)/365+yearReal+1
        lines(x,seasonal_median[q,],col="green")
        for (i in 1:4){
            sea=sea_limit[i]
            #lines(c(yearReal+sea/365,yearReal+sea/365),c(min(detrended[q,,year:(year+1)],na.rm=TRUE)-2.5,max(detrended[q,,year:(year+1)],na.rm=TRUE)+2.5),col=rgb(0.5,0.5,0.5,0.5))
            abline(v=yearReal+sea/365,col=rgb(0.5,0.5,0.5,0.5))
            abline(v=yearReal+1+sea/365,col=rgb(0.5,0.5,0.5,0.5))
            #lines(c(yearReal+1+sea/365,yearReal+1+sea/365),c(min(detrended[q,,year:(year+1)],na.rm=TRUE)-2.5,max(detrended[q,,year:(year+1)],na.rm=TRUE)+2.5),col=rgb(0.5,0.5,0.5,0.5))
            text(yearReal+(sea+45)/365,16,label=season_names[i],col=rgb(0.5,0.5,0.5,0.5))
            text(yearReal+1+(sea+45)/365,16,label=season_names[i],col=rgb(0.5,0.5,0.5,0.5))

        }

    }

    if (todo[3]==1){  
        time=(1:730)/365+yearReal
        warm=detrended[q,,year:(year+1)]
        warm[ind[q,,year:(year+1)]<0]=NA
        points(time,warm,pch=20,cex=0.5,col="red")
        cold=detrended[q,,year:(year+1)]
        cold[ind[q,,year:(year+1)]>0]=NA
        points(time,cold,pch=20,cex=0.5,col="blue")
    }
    #legend("bottomleft",pch=c(20,20),col=c("red","blue"),legend=c("warm day","cold day"))

    if (todo[4]==1){
        waIn=which(dur_mid[q,2,]>yearReal & dur_mid[q,2,]<yearReal+1.3)
        mi=max(detrended[q,,year:(year+1)],na.rm=TRUE)+0.7
        maxi=max(dur[q,2,waIn])
        oben=mi+2.5#+(-1)^p*1.25
        unten=mi#+(-1)^p*1.25
        polygon(x=c(yearReal,yearReal,yearReal+2,yearReal+2),c(unten-0.5,oben+0.5,oben+0.5,unten-0.5),col="white",border="white")
        for (p in waIn){    
            beg=dur_mid[q,2,p]-1/2*dur[q,2,p]/365
            end=dur_mid[q,2,p]+1/2*dur[q,2,p]/365
            polygon(x=c(beg,beg,end,end),y=c(unten,oben,oben,unten),border=NA,col=rgb(dur[q,2,p]/maxi/2+0.5,dur[q,2,p]/maxi/4+0.3,dur[q,2,p]/maxi/4+0.3))
            #abline(v=beg,col=rgb(0.9,0.9,0.9))
        }
        text(yearReal+1.5,mi+1.25,label="duration of \n warm periods",cex=1)

        kaIn=which(dur_mid[q,1,]>yearReal & dur_mid[q,1,]<yearReal+1.3)
        mi=min(detrended[q,,year:(year+1)],na.rm=TRUE)-0.7
        maxi=max(dur[q,1,waIn])
        oben=mi#+(-1)^p*1.25
        unten=mi-2.5#+(-1)^p*1.25 
        polygon(x=c(yearReal,yearReal,yearReal+2,yearReal+2),c(unten-0.5,oben+0.5,oben+0.5,unten-0.5),col="white",border="white")
        for (p in kaIn){
            beg=dur_mid[q,1,p]-1/2*dur[q,1,p]/365
            end=dur_mid[q,1,p]+1/2*dur[q,1,p]/365
            polygon(x=c(beg,beg,end,end),y=c(unten,oben,oben,unten),border=NA,col=rgb(dur[q,1,p]/maxi/4+0.3,dur[q,1,p]/maxi/4+0.3,dur[q,1,p]/maxi/2+0.5))
            #abline(v=beg,col=rgb(0.9,0.9,0.9))
        }
        text(yearReal+1.5,mi-1.25,label="duration of \n cold periods",cex=1)
    }

    if (todo[5]==1){
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

    graphics.off()
}

if (1==1){
    #method_explanation(2003,c(0,2,2,2),"1")
    method_explanation(2010,c(1,2,2,2,2),"1")
    method_explanation(2010,c(1,1,2,2,2),"2")
    method_explanation(2010,c(1,1,1,2,2),"3")
    method_explanation(2010,c(1,1,1,1,2),"4")
    #method_explanation(2003,c(1,1,1,1,1),"5")
}

method_seasonal_median <- function(q=507,yearReal=2003,trendID="91_5",trend_style="_mean",dataset="_TX",additional_style="_seasonal_median"){

    dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
    trend=trend_load(paste("../data/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    nc=open.ncdf(paste("../data/",trendID,"/",trendID,trend_style,dataset,additional_style,".nc",sep=""))
    seasonal_median=get.var.ncdf(nc,additional_style)

    year=yearReal-1949

    pdf(file="../plots/zwischenzeugs/state_attribution.pdf")
    par(mfrow=c(2,1))
    par(mar=c(3,5,4,3))

    
    plot(NA,xlim=c(1950,2014),ylim=c(-22,25),xlab="",ylab="temp deg C",main="detrending")
    points(dat$time,dat$tas[q,,],xlim=c(year,year+1),cex=0.1)
    lines(dat$time,trend[q,,],col="green")
    abline(v=c(2003,2004),col="violet")

    plot(NA,xlim=c(2003,2004),ylim=c(-5,5),xlab="",ylab="temp deg C",main="difference between annual and seasonal median")
    #points(dat$time,dat$tas[q,,]-trend[q,,])
    x=seq(2003,2004,1/365)[1:365]
    lines(x,seasonal_median[q,],col="orange",lty=1)    
    x=seq(2004,20045,1/365)[1:365]
    lines(x,seasonal_median[q,],col="orange",lty=1)

    abline(h=median(dat$tas[q,,]-trend[q,,],na.rm=TRUE),col="black")

    time=seq(yearReal,(yearReal+1),1/365)[1:365]
    time=time+0.1*1/365
    warm=dat$tas[q,,year]-trend[q,,year]
    warm[warm<median(dat$tas[q,,]-trend[q,,],na.rm=TRUE)]=NA
    points(time,warm,pch=2,cex=0.5,col="red")
    cold=dat$tas[q,,year]-trend[q,,year]
    cold[cold>median(dat$tas[q,,]-trend[q,,],na.rm=TRUE)]=NA
    points(time,cold,pch=2,cex=0.5,col="blue")

    time=seq(yearReal,(yearReal+1),1/365)[1:365]
    time=time+0.9*1/365
    warm=dat$tas[q,,year]-trend[q,,year]
    warm[warm<seasonal_median[q,]]=NA
    points(time,warm,pch=6,cex=0.4,col="red")
    cold=dat$tas[q,,year]-trend[q,,year]
    cold[cold>seasonal_median[q,]]=NA
    points(time,cold,pch=6,cex=0.4,col="blue")

    time=seq(yearReal,(yearReal+1),1/365)[1:365]
    time=time+0.5*1/365
    critical=dat$tas[q,,year]-trend[q,,year]
    critical[critical>seasonal_median[q,] & critical>median(dat$tas[q,,]-trend[q,,],na.rm=TRUE)]=NA
    critical[critical<seasonal_median[q,] & critical<median(dat$tas[q,,]-trend[q,,],na.rm=TRUE)]=NA
    points(time,critical,pch=1,cex=1.5,col="green")



}

#method_seasonal_median()
