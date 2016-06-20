# this only creates explanatory plots

method_explanation <- function(yearReal=2003,todo=c(1,1,1,1),name="full",q=661){
    year=yearReal-1949
    
    trend=trend_load(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,"_trend",trend_style,dataset,".nc",sep=""))
    nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_seasonal_median",".nc",sep=""))
    seasonal_median=var.get.nc(nc,"seasonal_median")
    nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/",trendID,trend_style,dataset,"_state_ind.nc",sep=""))
    ind=var.get.nc(nc,"ind")
    nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_duration_4seasons.nc",sep=""))
    dur=var.get.nc(nc,"dur")
    dur_mid=var.get.nc(nc,"dur_mid")

    par(mar=c(4.2,4,1,1))
    plot(NA,xlim=c(yearReal+59/365,yearReal+1.48+59/365),ylim=c(-23,18),pch=20,cex=0.5,
        main="",ylab="",xlab="",frame.plot=TRUE) #ylab="temperature anomaly in deg C"

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
        lines(x,seasonal_median[q,],col="black",cex=2)

        x=(1:365)/365+yearReal+1
        lines(x,seasonal_median[q,],col="black",cex=2)
        for (i in 1:4){
            sea=sea_limit[i]
            #lines(c(yearReal+sea/365,yearReal+sea/365),c(min(detrended[q,,year:(year+1)],na.rm=TRUE)-2.5,max(detrended[q,,year:(year+1)],na.rm=TRUE)+2.5),col=rgb(0.5,0.5,0.5,0.5))
            abline(v=yearReal+sea/365,col=rgb(0.5,0.5,0.5,1))
            abline(v=yearReal+1+sea/365,col=rgb(0.5,0.5,0.5,1))
            #lines(c(yearReal+1+sea/365,yearReal+1+sea/365),c(min(detrended[q,,year:(year+1)],na.rm=TRUE)-2.5,max(detrended[q,,year:(year+1)],na.rm=TRUE)+2.5),col=rgb(0.5,0.5,0.5,0.5))
            text(yearReal+(sea+45)/365,16,label=season_names[i],col=rgb(0.5,0.5,0.5,1))
            text(yearReal+1+(sea+45)/365,16,label=season_names[i],col=rgb(0.5,0.5,0.5,1))

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
        polygon(x=c(yearReal+0.15,yearReal+0.15,yearReal+1.69,yearReal+1.69),c(unten-0.5,oben+0.5,oben+0.5,unten-0.5),col="white",border="white")
        for (p in waIn){    
            beg=dur_mid[q,2,p]-1/2*dur[q,2,p]/365
            end=dur_mid[q,2,p]+1/2*dur[q,2,p]/365
            polygon(x=c(beg,beg,end,end),y=c(unten,oben,oben,unten),border=NA,col=rgb(dur[q,2,p]/maxi/2+0.5,dur[q,2,p]/maxi/4+0.3,dur[q,2,p]/maxi/4+0.3))
            #if (dur[q,2,p]>20){text(dur_mid[q,2,p],mi+1.25,dur[q,2,p])}
            #abline(v=beg,col=rgb(0.9,0.9,0.9))
        }
        text(yearReal+1.5,mi+1.25,label="duration of \n warm periods",cex=1)

        kaIn=which(dur_mid[q,1,]>yearReal & dur_mid[q,1,]<yearReal+1.3)
        mi=min(detrended[q,,year:(year+1)],na.rm=TRUE)-0.7
        maxi=max(dur[q,1,waIn])
        oben=mi#+(-1)^p*1.25
        unten=mi-2.5#+(-1)^p*1.25 
        polygon(x=c(yearReal+0.15,yearReal+0.15,yearReal+1.69,yearReal+1.69),c(unten-0.5,oben+0.5,oben+0.5,unten-0.5),col="white",border="white")
        for (p in kaIn){
            beg=dur_mid[q,1,p]-1/2*dur[q,1,p]/365
            end=dur_mid[q,1,p]+1/2*dur[q,1,p]/365
            polygon(x=c(beg,beg,end,end),y=c(unten,oben,oben,unten),border=NA,col=rgb(dur[q,1,p]/maxi/4+0.3,dur[q,1,p]/maxi/4+0.3,dur[q,1,p]/maxi/2+0.5))
            #abline(v=beg,col=rgb(0.9,0.9,0.9))
        }
        text(yearReal+1.5,mi-1.25,label="duration of \n cold periods",cex=1)
    }

    if (todo[5]==1){

        waIn=which((dur_mid[q,2,]>(yearReal+152/365) & dur_mid[q,2,]<(yearReal+243/365)) | (dur_mid[q,2,]>(yearReal+1+152/365) & dur_mid[q,2,]<(yearReal+1+243/365)))
        mi=max(detrended[q,,year:(year+1)],na.rm=TRUE)+0.7
        maxi=max(dur[q,2,waIn])
        oben=mi+2.5#+(-1)^p*1.25
        unten=mi#+(-1)^p*1.25
        #polygon(x=c(yearReal+1.2,yearReal+1.2,yearReal+1.9,yearReal+1.9),c(unten-0.9,oben+0.9,oben+0.9,unten-0.9),col="white",border="white")
        polygon(x=c(yearReal+0.15,yearReal+0.15,yearReal+1.69,yearReal+1.69),c(unten-0.5,oben+0.5,oben+0.5,unten-0.5),col="white",border="white")
        for (p in waIn){    
            beg=dur_mid[q,2,p]-1/2*dur[q,2,p]/365
            end=dur_mid[q,2,p]+1/2*dur[q,2,p]/365
            polygon(x=c(beg,beg,end,end),y=c(unten,oben,oben,unten),border=NA,col=rgb(dur[q,2,p]/maxi/2+0.5,dur[q,2,p]/maxi/4+0.3,dur[q,2,p]/maxi/4+0.3))
            #if (dur[q,2,p]>20){text(dur_mid[q,2,p],mi+1.25,dur[q,2,p])}
            #abline(v=beg,col=rgb(0.9,0.9,0.9))
        }
        #text(yearReal+1.5,mi+1.25,label="duration of \n warm periods",cex=1)

        kaIn=which((dur_mid[q,1,]>(yearReal+152/365) & dur_mid[q,1,]<(yearReal+243/365)) | (dur_mid[q,1,]>(yearReal+1+152/365) & dur_mid[q,1,]<(yearReal+1+243/365)))
        mi=min(detrended[q,,year:(year+1)],na.rm=TRUE)-0.7
        maxi=max(dur[q,1,waIn])
        oben=mi#+(-1)^p*1.25
        unten=mi-2.5#+(-1)^p*1.25 
        #polygon(x=c(yearReal+1.2,yearReal+1.2,yearReal+1.9,yearReal+1.9),c(unten-0.9,oben+0.9,oben+0.9,unten-0.9),col="white",border="white")
        polygon(x=c(yearReal+0.15,yearReal+0.15,yearReal+1.69,yearReal+1.69),c(unten-0.5,oben+0.5,oben+0.5,unten-0.5),col="white",border="white")
        for (p in kaIn){
            beg=dur_mid[q,1,p]-1/2*dur[q,1,p]/365
            end=dur_mid[q,1,p]+1/2*dur[q,1,p]/365
            polygon(x=c(beg,beg,end,end),y=c(unten,oben,oben,unten),border=NA,col=rgb(dur[q,1,p]/maxi/4+0.3,dur[q,1,p]/maxi/4+0.3,dur[q,1,p]/maxi/2+0.5))
            #abline(v=beg,col=rgb(0.9,0.9,0.9))
        }
        #text(yearReal+1.5,mi-1.25,label="duration of \n cold periods",cex=1)
    }

}

if (TRUE){
    library(rworldmap)
    library(fields)
    worldmap = getMap(resolution = "low")
    q=661
    pdf(file=paste("../plots/",dataset,"/sonstiges/method_explanation/explanatory.pdf",sep=""),width=6,height=5)      

    plot(worldmap,xlim=c(dat$lon[q]-20,dat$lon[q]+20),ylim=c(dat$lat[q]-20,dat$lat[q]+20))
    polygon(x=c(dat$lon[q]-pch_points[3],dat$lon[q]+pch_points[3],dat$lon[q]+pch_points[3],dat$lon[q]-pch_points[3]),y=c(dat$lat[q]-pch_points[4],dat$lat[q]-pch_points[4],dat$lat[q]+pch_points[4],dat$lat[q]+pch_points[4]),border=rgb(1,1,1,0.0),col="red")
    
    method_explanation(q=q,yearReal=2010,todo=c(0,2,2,2,2))
    method_explanation(q=q,yearReal=2010,todo=c(1,2,2,2,2))
    method_explanation(q=q,yearReal=2010,todo=c(1,1,2,2,2))
    method_explanation(q=q,yearReal=2010,todo=c(1,1,1,2,2))
    method_explanation(q=q,yearReal=2010,todo=c(1,1,1,1,2))
    method_explanation(q=q,yearReal=2010,todo=c(1,1,1,2,1))
    graphics.off()
}

#method_explanation(yearReal=1965,todo=c(1,1,1,1,2),name="test",q=436)

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
