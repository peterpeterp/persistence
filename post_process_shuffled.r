
confidence_interval <- function(seasons=1:4){
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
        quant_regressions[sea,,,]=original_mass[4,,,]

        for (q in 1:regNumb){
            for (state in 1:2){
                for (out in c(1,2,3,4,5)){
                    confi_quantiles[sea,q,state,out,]=quantile(shuffled_mass[,q,state,out],c(0.025,0.05,0.5,0.95,0.975),na.rm=TRUE)
                }
            }
        }
    }
    pdf(paste("../plots/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,dataset,"_",ID_name,"_",period,"_rq_trends.pdf",sep=""))
    for (sea in seasons){
        season<-season_names[sea]
        for (state in 1:2){
            for (out in c(2,3,5)){
                plot(NA,xlim=c(1,plotNumb),ylim=c(-0.3,0.3),axes=FALSE,ylab="period length per year",xlab="",main=paste(season,state_names[state],out_names[out]))
                axis(2)
                axis(1,1:plotNumb,labels=region_names[plot_select])
                polygon(x=c(1:plotNumb,plotNumb:1),y=c(confi_quantiles[sea,plot_select,state,out,2],confi_quantiles[sea,plot_select[plotNumb:1],state,out,4]),border="white",col=rgb(0.5,0.5,0.5,0.5))
                points(1:plotNumb,quant_regressions[sea,plot_select,state,out],col="red")
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

    period<<-"1980-2014"
    ID_name<<-"ward23"
    folder<<-paste("/regional/",ID_name,"/",sep="")
    regNumb<<-23
    region_names<<-1:23

    plot_select<<-c(3,4,5,7,11,12,13,14,16,18,20,21)
    #plot_select<<-c(11,12,16,20)
    plotNumb<<-length(plot_select)
}

init()
confidence_interval()