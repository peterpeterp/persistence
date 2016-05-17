

evaluate_bootstraping <- function(x=1,out=5,region=12,state=2,season="JJA",period="1950-2014",region_name="ward24"){

	filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/bootstrap/",trendID,dataset,"_",period,"_linear_trend_boot_",season,"_",state_names[state],"_",region_name,"-",region,".nc",sep="") ; print(filename)
	trends<-var.get.nc(open.nc(filename),"trends")
	lines(c(x,x),quantile(trends[2:10000,out],c(0.1,0.9)),col="gray")
	points(x,trends[1,out],col="red")
	text(x,0.0008,round(trends[1,out]/quantile(trends[2:10000,out],0.90),03))
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


nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_shuffQuant.nc",sep=""))
original_slopes=var.get.nc(nc,"original_slopes")
confi_quantiles=var.get.nc(nc,"confi_quantiles")

pdf("../plots/test_boot.pdf")
for (out in c(2,5)){
	for (state in 1:2){
		plot(NA,xlim=c(0,10),ylim=c(-0.001,0.001))
		count<-0
		for (reg in c(3,4,7,12,16,20)){
			count<-count+1
			evaluate_bootstraping(x=count,out=out,region=reg,state=state)
		}
	}
}

for (out in c(2,5)){
	for (state in 1:2){
		plot(NA,xlim=c(0,10),ylim=c(-0.1,0.1))
		count<-0
		for (reg in c(3,4,7,12,16,20)){
			count<-count+1
			lines(c(count+0.1,count+0.1),confi_quantiles[2,reg,state,3,c(2,4)],col="lightblue")
			points(count+0.1,original_slopes[2,reg,state,3,1],col="green")
			text(count,0.08,round(original_slopes[2,reg,state,3,1]/confi_quantiles[2,reg,state,3,4],03))
		}
	}
}

graphics.off()