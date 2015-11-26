
station_plot <- function(q,yearPeriod=c(1950,2014)){
    season_names=c("MAM","JJA","SON","DJF","year","4seasons")
    state=2
    color=c("blue","red")
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/stations/",q,".pdf",sep=""))
    for (sea in 1:6){
        #nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/duration/",period,"/",trendID,dataset,"_duration_analysis_",season_names[sea],".nc",sep=""))        
        nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_",season_names[sea],".nc",sep=""))
        dur=get.var.ncdf(nc,"dur")
        dur_mid=get.var.ncdf(nc,"dur_mid")
        plot(dur_mid[q,state,],dur[q,state,],col=color[state])
        quantreg=rq(dur[q,state,]~dur_mid[q,state,],c(0.5))
        abline(quantreg)        
        quantreg=rq(dur[q,state,]~dur_mid[q,state,],c(0.95))
        abline(quantreg)
        lm=lm(dur[q,state,]~dur_mid[q,state,])
        abline(lm,lty=2)
    }

}

source("load.r")
library(quantreg)

trendID="91_5"
dataset="_TMean"
additional_style=""

dat=dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))

stations=c(488)
for (q in stations){
    station_plot(q=q)
}