
station_plot <- function(q,yearPeriod=c(1950,2014),state=2,todo=c(1,1,1,1),zusatz=1){
    season_names=c("MAM","JJA","SON","DJF","year","4seasons")
    
    color=c("blue","red")
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/stations/",q,"_",zusatz,".pdf",sep=""),width=5,height=4)
    for (sea in 1:6){
        #nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",period,"/",trendID,dataset,"_duration_analysis_",season_names[sea],".nc",sep=""))        
        nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_",season_names[sea],".nc",sep=""))
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")
        plot(dur_mid[q,state,],dur[q,state,],col=color[state],ylab="days",xlab="",cex=0.4)

        if (todo[1]==1){
            lm=lm(dur[q,state,]~dur_mid[q,state,])
            abline(lm,lty=2)  
        }      

        if (todo[2]==1){
            quantreg=rq(dur[q,state,]~dur_mid[q,state,],c(0.95))
            print(summary(quantreg,se="boot"))
            abline(quantreg)  
        }

        if (todo[3]==1){
            legend("topleft",legend=c(paste("95th percentile slope:",round(quantreg$coefficients[2],02),"\n significance level:",round(summary(quantreg,se="boot")$coef[8],02))),bty="n")
        }

        maxis=order(dur[q,state,],decreasing=TRUE)
        print(season_names[sea])

    }

}

source("load.r")
library(quantreg)

trendID="91_5"
dataset="_TMean"
additional_style=""

dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))

stations=c(661)
for (q in stations){
    #station_plot(q=q,todo=c(2,2,2,2),zusatz=1)
    #station_plot(q=q,todo=c(1,2,2,2),zusatz=2)
    station_plot(q=q,todo=c(1,1,2,2),zusatz=3)
    station_plot(q=q,todo=c(1,1,1,2),zusatz=4)
    #station_plot(q=q,todo=c(1,1,1,2))
}