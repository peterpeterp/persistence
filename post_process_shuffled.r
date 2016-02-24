
confidence_interval <- function(seasons=2,yearPeriod=c(1950,2014),folder="/regional/7rect/",ID_name="7rect",ID_length=7){

	quantile_regression<-array(NA,dim=c(5,ID_length,2,9,3))

    for (sea in seasons){
        season<-season_names[sea]
        shuffled_mass<-array(NA,dim=c(10000,ID_length,2,5))
        original_mass<-array(NA,dim=c(100,ID_length,2,5))
        for (id in 1:60){
        	print(paste("../data/",dataset,additional_style,"/",trendID,folder,"/shuffled/",trendID,dataset,"_",ID_name,"_",yearPeriod[1],"-",yearPeriod[2],"_shuffled_trends_",season,"_",id,".nc",sep=""))
        	nc <- try(open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,"/shuffled/",trendID,dataset,"_",ID_name,"_",yearPeriod[1],"-",yearPeriod[2],"_shuffled_trends_",season,"_",id,".nc",sep="")))
        	if (class(nc)!="try-error"){
	        	shuffled<<-var.get.nc(nc,"shuffled")
	        	original<-var.get.nc(nc,"original")
	        	shuffled_mass[((id-1)*100+1):(id*100),,,]=shuffled
	        	original_mass[id,,,]=original
	        }
        }
        shuffled_mass<<-shuffled_mass
        original_mass<<-original_mass

        pdf("../plots/test1.pdf")
        confi_quantiles<-array(NA,dim=c(ID_length,2,5,5))
        for (q in 1:ID_length){
            for (state in 1:2){
                for (out in c(1,2,3,4,5)){
                    confi_quantiles[q,state,out,]=quantile(shuffled_mass[,q,state,out],c(0.025,0.05,0.5,0.95,0.975),na.rm=TRUE)
                }
                plot(NA,xlim=c(1,5),ylim=c(-0.2,0.2))
                polygon(x=c(1:5,5:1),y=c(confi_quantiles[q,state,1:5,1],confi_quantiles[q,state,5:1,5]),border="white",col=rgb(0.5,0.5,0.5,0.5))
                points(1:5,original_mass[1,q,state,1:5],col="red")
                #text(1:5,original_mass[1,q,state,1:5],label=round(original_mass[1,q,state,1:5],03),col="red")
                abline(h=0)
            }
        }
        graphics.off()
    }

}



init <- function(){
    library(quantreg)
    library(RNetCDF)


    nday<<-91
    nyr<<-5
    trendID<<-paste(nday,"_",nyr,sep="")
    dataset<<-"_TMean"
    additional_style<<-""

    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
}

init()
confidence_interval()