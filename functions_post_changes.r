

evaluate_bootstraping <- function(overRegions=c("NHml","NHpo","NHst","tro","SHml"),period="1979-2011",region_name="ward24"){
	result<-array(NA,c(5,2,30,4,10))
	for (sea in 1:5){
		for (state in 1:2){
			shuffled<-array(NA,c(10001,30,4))
			for (oR in 1:length(overRegions)){
				
				for (i in 1:10){
					filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/bootstrap/",trendID,dataset,"_",period,"_boot_",region_name,"_",overRegions[oR],"_",season_names[sea],"_",state_names[state],"_",i,".nc",sep="") ; print(filename)
					regions<-try(var.get.nc(open.nc(filename),"regions"))
					#print(regions)
					if (regions[1]!="Error : Datei oder Verzeichnis nicht gefunden\n"){
						regions<-regions[1:(length(regions)-1)]

						#  for Regions
						shuffled[((i-1)*1000+2):(i*1000+1),regions,]<-var.get.nc(open.nc(filename),"statistics")[2:1001,1:length(regions),]
						shuffled[1,regions,]<-var.get.nc(open.nc(filename),"statistics")[1,1:length(regions),]

						# for overRegions
						shuffled[((i-1)*1000+2):(i*1000+1),(25+oR),]<-var.get.nc(open.nc(filename),"statistics")[2:1001,(length(regions)+1),]
						shuffled[1,(25+oR),]<-var.get.nc(open.nc(filename),"statistics")[1,(length(regions)+1),]
					}
				}

				for (out in 1:4){
					for (q in c(1:24,26:30)){
						result[sea,state,q,out,1]<-shuffled[1,q,out]
						result[sea,state,q,out,3:8]<-quantile(shuffled[2:10001,q,out],c(0.025,0.05,0.1,0.9,0.95,0.975),na.rm=TRUE)
					}
				}
			}
		}
	}
	statistics<-array(NA,c(5,30,2,4,10))
    for (sea in 1:5){
        for (q in 1:30){
            for (state in 1:2){
                statistics[sea,q,state,,]=result[sea,state,q,,]
            }
        }
    }
	filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_bootstrap.nc",sep="") ; print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "bootstrapping results", "NC_CHAR","see regions")
    att.put.nc(nc_out, "NC_GLOBAL", "replics", "NC_CHAR","10000")
            
    dim.def.nc(nc_out,"seasons",dimlength=5, unlim=FALSE)
    dim.def.nc(nc_out,"regions",dimlength=30,unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2, unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=4,unlim=FALSE)
    dim.def.nc(nc_out,"confidence",dimlength=10,unlim=FALSE)

    var.def.nc(nc_out,"statistics","NC_DOUBLE",c(0,1,2,3,4))
    att.put.nc(nc_out, "statistics", "missing_value", "NC_DOUBLE", -99999)
    att.put.nc(nc_out, "statistics", "dim_explanation", "NC_CHAR", "season-(regions,NA,overRegions)-state-(lr,qr,ks_full,ks_upper)-(original,NA,0.025,0.05,0.1,0.9,0.95,0.975,NA,NA)")
    var.put.nc(nc_out,"statistics",statistics)   
    close.nc(nc_out) 
}

plot_ks_statistic_distrs <- function(overRegions=c("NHml","NHpo","NHst","tro","SHml"),period="1979-2011",periods=c("1979-1995","1995-2011"),region_name="ward24"){
    pdf(paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_dur_ks_test_",periods[1],"_vx_",periods[2],"_KS.pdf",sep=""),width=3,height=3)
    par(mar=c(3, 3, 3, 3) + 0.1)  
    color=c(NA,NA,rgb(0.7,0.2,0.7),rgb(0.5,0.7,0.5))
    ymax<-2000
    plot(NA,xlab="",ylim=c(0,ymax),xlim=c(0,0.2),ylab="",axes=FALSE,frame.plot=TRUE,cex=0.5) ; at_=axis(1,labels=FALSE,col="black") ; if (length(at_)>4){at_=at_[2:(length(at_)-1)]} ; axis(1,at=at_)
    plot(NA,xlab="",ylim=c(0,ymax),xlim=c(0,0.2),ylab="",axes=FALSE,frame.plot=TRUE,cex=0.5) ; at_=axis(3,labels=FALSE,col="black") ; if (length(at_)>4){at_=at_[2:(length(at_)-1)]} ; axis(3,at=at_)
    plot(NA,xlab="",ylim=c(0,ymax),xlim=c(0,0.2),ylab="",axes=FALSE,frame.plot=TRUE,cex=0.5) ; at_=axis(2,labels=FALSE,col="black") ; if (length(at_)>4){at_=at_[2:(length(at_)-1)]} ; axis(2,at=at_)
    plot(NA,xlab="",ylim=c(0,ymax),xlim=c(0,0.2),ylab="",axes=FALSE,frame.plot=TRUE,cex=0.5) ; at_=axis(4,labels=FALSE,col="black") ; if (length(at_)>4){at_=at_[2:(length(at_)-1)]} ; axis(4,at=at_)

    for (sea in 1:4){
        for (state in 1:2){
            shuffled<-array(NA,c(10001,30,4))
            for (oR in 1:length(overRegions)){
                
                for (i in 1:10){
                    filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/bootstrap/",trendID,dataset,"_",period,"_boot_",region_name,"_",overRegions[oR],"_",season_names[sea],"_",state_names[state],"_",i,".nc",sep="") ; print(filename)
                    regions<-try(var.get.nc(open.nc(filename),"regions"))
                    #print(regions)
                    if (regions[1]!="Error : Datei oder Verzeichnis nicht gefunden\n"){
                        regions<-regions[1:(length(regions)-1)]

                        #  for Regions
                        shuffled[((i-1)*1000+2):(i*1000+1),regions,]<-var.get.nc(open.nc(filename),"statistics")[2:1001,1:length(regions),]
                        shuffled[1,regions,]<-var.get.nc(open.nc(filename),"statistics")[1,1:length(regions),]
                    }
                }
                #plot
                for (q in regions){
                    plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=FALSE,frame.plot=FALSE) ; text(0.5,0.5,paste(season_names[sea],q,state_names[state]))
                    for (out in c(3,4)){   
                        if (max(shuffled[2:10001,q,out],na.rm=TRUE)<0.3){     
                            hist(shuffled[2:10001,q,out],br=seq(0,0.3,0.005),ylim=c(0,ymax),xlim=c(0,0.2),main="",ylab="",xlab="",axes=FALSE,col=color[out],border=color[out])
    
                            quks<-quantile(shuffled[2:10001,q,out],c(0.9,0.95),na.rm=TRUE)
                            abline(v=quks[1],lty=3,col="gray")
                            abline(v=quks[2],lty=2,col="gray")
                            abline(v=shuffled[1,q,out],lty=1,col="black")

                            if (out==3){text(0.085,1800,expression('KS'['full']),col="black",pos=4)}
                            if (out==4){text(0.07,1800,expression('KS'['upper']),col="black",pos=4)}
                            if (shuffled[1,q,out]<quks[1]){text(0.135,1800,"<90%",pos=4)}
                            if (shuffled[1,q,out]>quks[1] & shuffled[1,q,out]<quks[2]){text(0.135,1800,">90%",pos=4)}
                            if (shuffled[1,q,out]>quks[2]){text(0.135,1800,">95%",pos=4)}
                        }
                        if (max(shuffled[2:10001,q,out],na.rm=TRUE)>0.3){plot(1:10)}
                    }
                }

            }
        }
    }
    graphics.off() 
}

write_slope_table_control <- function(period="1979-2011",region_name="ward24"){
    filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_quantiles.nc",sep="") ; print(filename)
    qr<-var.get.nc(open.nc(filename),"quantile_stuff")
    filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_others.nc",sep="") ; print(filename)
    mn<-var.get.nc(open.nc(filename),"other_stuff")

    signis<-array(NA,c(5,24,2,5,2))

    values<-array(NA,c(5,24,2,2))
    print(dim(mn))
    print(dim(values))
    values[,,,1]=mn[1:5,,,4]
    values[,,,2]=qr[1:5,,,2,2]

    print(values)
    values[which(is.na(values))]=0
    values<-values

    nbcol<<-101
    plot_reg_table_general(values=values,signis=signis,filename=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,dataset,"_",region_name,"_",period,"_slope_control.pdf",sep=""),val_names=c("mn","95"),region_name="ward24",colorRange=c(-0.1,0.1),farb_palette="lila-gruen",ID_select=reg_order,hlines=hlines)
}

write_slope_table <- function(period="1979-2011",region_name="ward24"){
    filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_bootstrap.nc",sep="") ; print(filename)
    original_slopes<-var.get.nc(open.nc(filename),"statistics")

    signis<-array(NA,c(5,24,2,5,2))
    signis[,,,5,1][which(original_slopes[,1:24,,1,1]>original_slopes[,1:24,,1,6])]=1
    signis[,,,5,1][which(original_slopes[,1:24,,1,1]<original_slopes[,1:24,,1,5])]=1
    signis[,,,4,1][which(original_slopes[,1:24,,1,1]>original_slopes[,1:24,,1,7])]=1
    signis[,,,4,1][which(original_slopes[,1:24,,1,1]<original_slopes[,1:24,,1,4])]=1

    signis[,,,5,2][which(original_slopes[,1:24,,2,1]>original_slopes[,1:24,,2,6])]=1
    signis[,,,5,2][which(original_slopes[,1:24,,2,1]<original_slopes[,1:24,,2,5])]=1
    signis[,,,4,2][which(original_slopes[,1:24,,2,1]>original_slopes[,1:24,,2,7])]=1
    signis[,,,4,2][which(original_slopes[,1:24,,2,1]<original_slopes[,1:24,,2,4])]=1

    values<-array(NA,c(5,24,2,2))
    values[,,,1]=original_slopes[,1:24,,1,1]
    values[,,,2]=original_slopes[,1:24,,2,1]

    values[which(is.na(values))]=0
    values<-values*3650

    plot_reg_table_general(values=values,signis=signis,filename=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,dataset,"_",region_name,"_",period,"_slope_boot.pdf",sep=""),val_names=c("mn","95"),region_name="ward24",colorRange=c(-1,1),farb_palette="lila-gruen",ID_select=reg_order,regLabel=1:24,hlines=hlines,colorbar=FALSE)

    signis<-array(NA,c(5,5,2,5,2))
    #signis[,,,4,1][which(original_slopes[,,,1,1]>original_slopes[,,,1,6])]=1
    #signis[,,,4,1][which(original_slopes[,,,1,1]<original_slopes[,,,1,5])]=1
    signis[,,,5,1][which(original_slopes[,26:30,,1,1]>original_slopes[,26:30,,1,7])]=1
    signis[,,,5,1][which(original_slopes[,26:30,,1,1]<original_slopes[,26:30,,1,4])]=1

    #signis[,,,4,2][which(original_slopes[,,,2,1]>original_slopes[,,,2,6])]=1
    #signis[,,,4,2][which(original_slopes[,,,2,1]<original_slopes[,,,2,5])]=1
    signis[,,,5,2][which(original_slopes[,26:30,,2,1]>original_slopes[,26:30,,2,7])]=1
    signis[,,,5,2][which(original_slopes[,26:30,,2,1]<original_slopes[,26:30,,2,4])]=1


    values<-array(NA,c(5,5,2,2))
    values[,,,1]=original_slopes[,26:30,,1,1]
    values[,,,2]=original_slopes[,26:30,,2,1]

    values[which(is.na(values))]=0
    values<-values*3650

    plot_reg_table_general(values=values,signis=signis,filename=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,dataset,"_",region_name,"_",period,"_slope_boot_overReg.pdf",sep=""),val_names=c("mn","95"),region_name="ward24",colorRange=c(-1,1),farb_palette="lila-gruen",ID_select=1:5,regLabel=c("NHml","NHpo","NHst","Tro","SHml"),paperHeight=6,hlines=c(30),colorbar=TRUE,header=FALSE)
}

plot_distr_compa_table <- function(period="1979-2011",periods=c("1979-1995","1995-2011"),regNumb=24,region_name="ward24"){

    filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_bootstrap.nc",sep="") ; print(filename)
    original_ks<-var.get.nc(open.nc(filename),"statistics")[,1:24,,3:4,]

    fit_params=array(NA,c(2,5,regNumb,2,30))
    quantiles=array(NA,c(2,5,regNumb,2,2,3))
    for (i in 1:length(periods)){
        part_period<-periods[i]
        filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",part_period,"/",trendID,"_",dataset,"_",region_name,"_",part_period,"_fit_","2expo_4:100",".nc",sep=""); print(filename)
        fit_params[i,,,,]=var.get.nc(open.nc(filename),"fit_stuff")[1:5,,,]

        filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",part_period,"/",trendID,"_",dataset,"_",region_name,"_",part_period,"_quantiles",".nc",sep=""); print(filename)
        quantiles[i,,,,2,]=var.get.nc(open.nc(filename),"quantile_stuff")[1:5,,,2,]

        filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",part_period,"/",trendID,"_",dataset,"_",region_name,"_",part_period,"_others",".nc",sep=""); print(filename)
        quantiles[i,,,,1,1]=var.get.nc(open.nc(filename),"other_stuff")[1:5,,,1]
    }

    signis<-array(NA,c(5,regNumb,2,5,2))
    signis[,,,5,1][which(original_ks[,,,1,1]>original_ks[,,,1,6])]=1
    signis[,,,4,1][which(original_ks[,,,1,1]>original_ks[,,,1,7])]=1
    signis[,,,5,2][which(original_ks[,,,2,1]>original_ks[,,,2,6])]=1
    signis[,,,4,2][which(original_ks[,,,2,1]>original_ks[,,,2,7])]=1

    values<-array(NA,c(5,regNumb,2,2))
    values[,,,1]=quantiles[2,,,,1,1]-quantiles[1,,,,1,1]
    values[,,,2]=quantiles[2,,,,2,1]-quantiles[1,,,,2,1]

    plot_reg_table_general(values=values,signis=signis,filename=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,dataset,"_",region_name,"_dur_ks_test_",periods[1],"_vx_",periods[2],"_quantile_diff.pdf",sep=""),val_names=c("mn","95"),region_name="ward24",colorRange=c(-2,2),farb_palette="lila-gruen",ID_select=reg_order,hlines=hlines)



    values<-array(NA,c(5,regNumb,2,2))
    values[,,,1]=fit_params[2,,,,12]-fit_params[1,,,,12]
    values[,,,2]=fit_params[2,,,,14]-fit_params[1,,,,14]

    values[is.na(values)]=0 
    signis[,,,2,][is.na(values)]=1

    signis[,,,2,1][which(fit_params[1,1:5,,,21]<0.99)]=1
    signis[,,,2,2][which(fit_params[1,1:5,,,21]<0.99)]=1
    signis[,,,2,1][which(fit_params[2,1:5,,,21]<0.99)]=1
    signis[,,,2,2][which(fit_params[2,1:5,,,21]<0.99)]=1

    signis[,,,4:5,1:2]=NA
    signis[,,,5,1][which(original_ks[,,,2,1]>original_ks[,,,2,6])]=1
    signis[,,,4,1][which(original_ks[,,,2,1]>original_ks[,,,2,7])]=1
    signis[,,,5,2][which(original_ks[,,,2,1]>original_ks[,,,2,6])]=1
    signis[,,,4,2][which(original_ks[,,,2,1]>original_ks[,,,2,7])]=1


    plot_reg_table_general(values=values,signis=signis,filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,dataset,"_",region_name,"_dur_ks_test_",periods[1],"_vx_",periods[2],"_slopeDiff.pdf",sep=""),val_names=c("b1","b2"),region_name="ward24",colorRange=c(-0.1,0.1),farb_palette="lila-gruen-inv",ID_select=reg_order,hlines=hlines)

    values<-array(NA,c(5,regNumb,2,1))
    values[,,,1]=fit_params[2,,,,15]-fit_params[1,,,,15]   
    values[is.na(values)]=3

    plot_reg_table_general(values=values,signis=signis,filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,dataset,"_",region_name,"_dur_ks_test_",periods[1],"_vx_",periods[2],"_threshDiff.pdf",sep=""),val_names=c(""),region_name="ward24",colorRange=c(-4,4),farb_palette="lila-gruen-inv",ID_select=reg_order,hlines=hlines)

}


init <- function(){
    source("plot_tables.r")
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

    reg_order<<-c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24)

    hlines<<-c(23,20,22,8)
}

init()

#evaluate_bootstraping()

write_slope_table()
write_slope_table_control()
plot_distr_compa_table()
plot_ks_statistic_distrs()
