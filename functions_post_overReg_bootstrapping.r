

evaluate_bootstraping <- function(region="mid-lat-3-4-7-12-16-20_100",state=2,season="JJA",period="1950-2014",region_name="ward24"){

	trends<-array(NA,c(10000,3))
	for (i in 1:100){
		filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/bootstrap/",trendID,dataset,"_",period,"_linear_trend_boot_",season,"_",state_names[state],"_",region_name,"-",region,".nc",sep="") ; print(filename)
		trends[((i-1)*100+2):(i*100),]<-var.get.nc(open.nc(filename),"trends")[2:100,]
	}
	for (out in c(1,3)){
		print(quantile(trends[2:10000,out],c(0.1,0.9),na.rm=TRUE))
		print(var.get.nc(open.nc(filename),"trends")[1,out])
	}
}



init <- function(){
    library(RNetCDF)
    library(quantreg)
    nday<<-91
    nyr<<-7
    trendID<<-paste(nday,"_",nyr,sep="")
    dataset<<-"_TMean"
    additional_style<<-""
    taus<<-c(0.75,0.95,0.99)
    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    state_names<<-c("cold","warm")

    taus<<-c(0.75,0.95,0.99)
}

init()

region_name<-"ward24"
period<-"1950-2014"

evaluate_bootstraping()
