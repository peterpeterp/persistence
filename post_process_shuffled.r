
confidence_interval <- function(seasons=1:4,period="1980-2014",folder="/regional/7rect/",ID_name="7rect",regNumb=7,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){
	confi_quantiles<-array(NA,dim=c(5,regNumb,2,5,5))
    quant_regressions<-array(NA,dim=c(5,regNumb,2,5))
    for (sea in seasons){
        season<-season_names[sea]
        shuffled_mass<-array(NA,dim=c(10000,regNumb,2,5))
        original_mass<-array(NA,dim=c(100,regNumb,2,5))
        for (id in 1:100){
        	#print(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/shuffled/",trendID,dataset,"_",ID_name,"_",period,"_shuffled_trends_",season,"_",id,".nc",sep=""))
        	nc <- try(open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/shuffled/",trendID,dataset,"_",ID_name,"_",period,"_shuffled_trends_",season,"_",id,".nc",sep="")))
        	if (class(nc)!="try-error"){
	        	shuffled<-var.get.nc(nc,"shuffled")
	        	original<-var.get.nc(nc,"original")
	        	shuffled_mass[((id-1)*100+1):(id*100),,,]=shuffled
	        	original_mass[id,,,]=original
	        }
        }
        shuffled_mass<<-shuffled_mass
        original_mass<<-original_mass
        quant_regressions[sea,,,]=original_mass[3,,,]

        #pdf("../plots/test1.pdf")
        for (q in 1:regNumb){
            for (state in 1:2){
                for (out in c(1,2,3,4,5)){
                    confi_quantiles[sea,q,state,out,]=quantile(shuffled_mass[,q,state,out],c(0.025,0.05,0.5,0.95,0.975),na.rm=TRUE)
                }
                #plot(NA,xlim=c(1,5),ylim=c(-0.2,0.2))
                #polygon(x=c(1:5,5:1),y=c(confi_quantiles[q,state,1:5,1],confi_quantiles[q,state,5:1,5]),border="white",col=rgb(0.5,0.5,0.5,0.5))
                #points(1:5,original_mass[1,q,state,1:5],col="red")
                #text(1:5,original_mass[1,q,state,1:5],label=round(original_mass[1,q,state,1:5],03),col="red")
                #abline(h=0)
            }
        }
        #graphics.off()
    }
    pdf(paste("../plots/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,dataset,"_",ID_name,"_",period,"_rq_trends33.pdf",sep=""))
    for (sea in seasons){
        season<-season_names[sea]
        for (state in 1:2){
            for (out in c(2,3,5)){
                plot(NA,xlim=c(1,regNumb),ylim=c(-0.3,0.3),axes=FALSE,ylab="period length per year",xlab="",main=paste(season,state_names[state],out_names[out]))
                axis(2)
                axis(1,1:regNumb,labels=region_names)
                polygon(x=c(1:regNumb,regNumb:1),y=c(confi_quantiles[sea,1:regNumb,state,out,1],confi_quantiles[sea,regNumb:1,state,out,5]),border="white",col=rgb(0.5,0.5,0.5,0.5))
                points(1:regNumb,quant_regressions[sea,1:regNumb,state,out],col="red")
                abline(h=0)
            }
        }
        
    }
    graphics.off()
    quant_regressions<<-quant_regressions

}



init <- function(){
    library(quantreg)
    library(RNetCDF)


    nday<<-91
    nyr<<-5
    trendID<<-paste(nday,"_",nyr,sep="")
    dataset<<-"_TMean"
    additional_style<<-""

    state_names<<-c("cold","warm")
    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    out_names<<-c(0.5,0.75,0.95,0.99,"mean")
}

init()
confidence_interval()