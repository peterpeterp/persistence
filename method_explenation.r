

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