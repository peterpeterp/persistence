

sens_gridded <- function(trendIDs=c("91_5","91_7","91_9"),period="1950-2014",file="_others",var="other_stuff",sub_auswahl=NA,value_auswahl=1,name_zusatz="mean",farb_mitte=c(-10,10),farb_palette="blau-rot",signi_level=0.05){

	values<-array(NA,c(3,5,ntot,2))
	for (i in 1:length(trendIDs)){
		filename <- paste("../data/",dataset,additional_style,"/",trendIDs[i],"/gridded/",period,"/",trendIDs[i],"_",dataset,"_",period,file,".nc",sep="") ; print(filename)
		if (is.na(sub_auswahl)){values[i,,,]<-var.get.nc(open.nc(filename),var)[1:5,,,value_auswahl]}
		if (!is.na(sub_auswahl)){
			taus<-var.get.nc(open.nc(filename),"taus")
			tau<-which(taus==sub_auswahl)
			values[i,,,]<-var.get.nc(open.nc(filename),var)[1:5,,,tau,value_auswahl]
		}
	}
	
	reihen<-array(NA,dim=c(length(season_auswahl)*2,ntot))
	reihen_sig<-array(NA,dim=c(length(season_auswahl)*2,ntot))
    index<-0
	for (sea in season_auswahl){
		season<-season_names[sea]
		for (state in 1:2){
			index<-index+1
			reihen[index,]=(values[1,sea,,state]-values[2,sea,,state])/values[2,sea,,state]*100
		}
	}


    print(paste("ann >0:",length(which(reihen[9:10,]>0))/length(reihen[9:10,])))
    print(paste("ann abs<5:",length(which(abs(reihen[9:10,])<5))/length(reihen[9:10,])))
    print(paste("ann abs<10:",length(which(abs(reihen[9:10,])<10))/length(reihen[9:10,])))

    print(paste("seas >0:",length(which(reihen[1:8,]>0))/length(reihen[1:8,])))
    print(paste("seas abs<5:",length(which(abs(reihen[1:8,])<5))/length(reihen[1:8,])))
    print(paste("seas abs<10:",length(which(abs(reihen[1:8,])<10))/length(reihen[1:8,])))

    filename_plot<-paste("../plots/",dataset,additional_style,"/sensitivity/","sensitivity_",name_zusatz,"_",period,"_",name_style,"_5-7.pdf",sep="") ; print(filename_plot)
    indexBottomLeft<<-c("5 years\ncold","5 years\nwarm","5 years\ncold","5 years\nwarm","5 years\ncold","5 years\nwarm","5 years\ncold","5 years\nwarm","5 years\ncold","5 years\nwarm")
	topo_map_plot(filename_plot=filename_plot,reihen=reihen,reihen_sig=reihen_sig,farb_mitte=farb_mitte,farb_palette=farb_palette,signi_level=signi_level)

    reihen<-array(NA,dim=c(length(season_auswahl)*2,ntot))
    reihen_sig<-array(NA,dim=c(length(season_auswahl)*2,ntot))
    index<-0
    for (sea in season_auswahl){
        season<-season_names[sea]
        for (state in 1:2){
            index<-index+1
            reihen[index,]=(values[3,sea,,state]-values[2,sea,,state])/values[2,sea,,state]*100
        }
    }

    print(paste("ann >0:",length(which(reihen[9:10,]>0))/length(reihen[9:10,])))
    print(paste("ann abs<5:",length(which(abs(reihen[9:10,])<5))/length(reihen[9:10,])))
    print(paste("ann abs<10:",length(which(abs(reihen[9:10,])<10))/length(reihen[9:10,])))

    print(paste("seas >0:",length(which(reihen[1:8,]>0))/length(reihen[1:8,])))
    print(paste("seas abs<5:",length(which(abs(reihen[1:8,])<5))/length(reihen[1:8,])))
    print(paste("seas abs<10:",length(which(abs(reihen[1:8,])<10))/length(reihen[1:8,])))

    filename_plot<-paste("../plots/",dataset,additional_style,"/sensitivity/","sensitivity_",name_zusatz,"_",period,"_",name_style,"_9-7.pdf",sep="") ; print(filename_plot)
    indexBottomLeft<<-c("9 years\ncold","9 years\nwarm","9 years\ncold","9 years\nwarm","9 years\ncold","9 years\nwarm","9 years\ncold","9 years\nwarm","9 years\ncold","9 years\nwarm")
    topo_map_plot(filename_plot=filename_plot,reihen=reihen,reihen_sig=reihen_sig,farb_mitte=farb_mitte,farb_palette=farb_palette,signi_level=signi_level)

}


sens_regional_fits <- function(region_name="ward24",regNumb=24,trendIDs=c("91_5","91_7","91_9"),period="1950-2014",name_zusatz="fit",farb_mitte="",farb_palette="blau-rot",signi_level=0.05){

	values<-array(NA,c(3,5,regNumb,2,30))
	for (i in 1:length(trendIDs)){
		filename <- paste("../data/",dataset,additional_style,"/",trendIDs[i],"/regional/",region_name,"/",period,"/",trendIDs[i],"_",dataset,"_",region_name,"_",period,"_fit_2expo_4:100.nc",sep="") ; print(filename)
        values[i,,,,]<-var.get.nc(open.nc(filename),"fit_stuff")[1:5,,,]
	}
	
    signis<-array(NA,c(5,regNumb,2,5,2))

    signis[,,,4,1][which(values[1,1:5,,,24]>0)]=1
    signis[,,,4,2][which(values[1,1:5,,,24]>0)]=1
    signis[,,,4,1][which(values[2,1:5,,,24]>0)]=1
    signis[,,,4,2][which(values[2,1:5,,,24]>0)]=1

	reihen<-array(NA,dim=c(5,regNumb,2,2))
    for (sea in season_auswahl){
        season<-season_names[sea]
        for (state in 1:2){
            reihen[sea,,state,1]=(values[1,sea,,state,12]-values[2,sea,,state,12])/values[2,sea,,state,12]*100
            reihen[sea,,state,2]=(values[1,sea,,state,14]-values[2,sea,,state,14])/values[2,sea,,state,14]*100
        }
    }

    BICaccepted<-which(is.na(signis[,,,4,1:2]))
    print(paste(">0:",length(which(reihen[BICaccepted]>0))/length(reihen[BICaccepted])))
    print(paste("abs<5:",length(which(abs(reihen[BICaccepted])<5))/length(reihen[BICaccepted])))
    print(paste("abs<10:",length(which(abs(reihen[BICaccepted])<10))/length(reihen[BICaccepted])))

    reihen[which(is.na(reihen))]=0
    plot_reg_table_general(values=reihen,signis=signis,filename=paste("../plots/",dataset,additional_style,"/sensitivity/",region_name,"_sensitivity_",name_zusatz,"_",period,"_slopes_5-7.pdf",sep=""),val_names=c("b1","b2"),region_name="ward24",colorRange=c(-10,10),farb_palette=farb_palette,ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8),colorbar=FALSE)

    reihen<-array(NA,dim=c(5,regNumb,2,1))
    for (sea in season_auswahl){
        season<-season_names[sea]
        for (state in 1:2){
            reihen[sea,,state,1]=(values[1,sea,,state,15]-values[2,sea,,state,15])/values[2,sea,,state,15]*100
        }
    }

    BICaccepted<-which(is.na(signis[,,,4,1]))
    print(paste(">0:",length(which(reihen[BICaccepted]>0))/length(reihen[BICaccepted])))
    print(paste("abs<5:",length(which(abs(reihen[BICaccepted])<5))/length(reihen[BICaccepted])))
    print(paste("abs<10:",length(which(abs(reihen[BICaccepted])<10))/length(reihen[BICaccepted])))

    plot_reg_table_general(values=reihen,signis=signis,filename=paste("../plots/",dataset,additional_style,"/sensitivity/",region_name,"_sensitivity_",name_zusatz,"_",period,"_tresh_5-7.pdf",sep=""),val_names=c(""),region_name="ward24",colorRange=c(-10,10),farb_palette=farb_palette,ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8),colorbar=FALSE)



    signis<-array(NA,c(5,regNumb,2,5,2))
    signis[,,,4,1][which(values[2,1:5,,,24]>0)]=1
    signis[,,,4,2][which(values[2,1:5,,,24]>0)]=1
    signis[,,,4,1][which(values[3,1:5,,,24]>0)]=1
    signis[,,,4,2][which(values[3,1:5,,,24]>0)]=1

    reihen<-array(NA,dim=c(5,regNumb,2,2))
    for (sea in season_auswahl){
        season<-season_names[sea]
        for (state in 1:2){
            reihen[sea,,state,1]=(values[3,sea,,state,12]-values[2,sea,,state,12])/values[2,sea,,state,12]*100
            reihen[sea,,state,2]=(values[3,sea,,state,14]-values[2,sea,,state,14])/values[2,sea,,state,14]*100
        }
    }

    BICaccepted<-which(is.na(signis[,,,4,1:2]))
    print(paste(">0:",length(which(reihen[BICaccepted]>0))/length(reihen[BICaccepted])))
    print(paste("abs<5:",length(which(abs(reihen[BICaccepted])<5))/length(reihen[BICaccepted])))
    print(paste("abs<10:",length(which(abs(reihen[BICaccepted])<10))/length(reihen[BICaccepted])))

    reihen[which(is.na(reihen))]=0
    plot_reg_table_general(values=reihen,signis=signis,filename=paste("../plots/",dataset,additional_style,"/sensitivity/",region_name,"_sensitivity_",name_zusatz,"_",period,"_slopes_9-7.pdf",sep=""),val_names=c("b1","b2"),region_name="ward24",colorRange=c(-10,10),farb_palette=farb_palette,ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8))

    reihen<-array(NA,dim=c(5,regNumb,2,1))
    for (sea in season_auswahl){
        season<-season_names[sea]
        for (state in 1:2){
            reihen[sea,,state,1]=(values[3,sea,,state,15]-values[2,sea,,state,15])/values[2,sea,,state,15]*100
        }
    }

    BICaccepted<-which(is.na(signis[,,,4,1]))
    print(paste(">0:",length(which(reihen[BICaccepted]>0))/length(reihen[BICaccepted])))
    print(paste("abs<5:",length(which(abs(reihen[BICaccepted])<5))/length(reihen[BICaccepted])))
    print(paste("abs<10:",length(which(abs(reihen[BICaccepted])<10))/length(reihen[BICaccepted])))

    plot_reg_table_general(values=reihen,signis=signis,filename=paste("../plots/",dataset,additional_style,"/sensitivity/",region_name,"_sensitivity_",name_zusatz,"_",period,"_tresh_9-7.pdf",sep=""),val_names=c(""),region_name="ward24",colorRange=c(-10,10),farb_palette=farb_palette,ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8))
}

sens_regional_trends <- function(region_name="ward24",region_name2="ward24",regNumb=24,regPos=1:24,regLabel=1:24,trendIDs=c("91_5","91_7","91_9"),period="1950-2014",name_zusatz="mean",farb_mitte="",farb_palette="weiss-rot",signi_level=0.05){

    filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name2,"/",period,"/",trendID,"_",dataset,"_",region_name2,"_",period,"_bootstrap.nc",sep="") ; print(filename)
    original_slopes<-var.get.nc(open.nc(filename),"statistics")[,regPos,,,]

    filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,dataset,"_",region_name,"_",period,"_MK.nc",sep="") ; print(filename)
    MK<-var.get.nc(open.nc(filename),"MK")[,,,,]


    signis<-array(NA,c(5,regNumb,2,20,2))
    signis[,,,11,1][which(original_slopes[,,,1,1]>original_slopes[,,,1,8])]=1
    signis[,,,11,1][which(original_slopes[,,,1,1]<original_slopes[,,,1,3])]=1
    signis[,,,14,2][which(original_slopes[,,,2,1]>original_slopes[,,,2,8])]=1
    signis[,,,14,2][which(original_slopes[,,,2,1]<original_slopes[,,,2,3])]=1

    signis[,,,12,1][which(original_slopes[,,,3,1]>original_slopes[,,,3,7])]=1
    signis[,,,15,2][which(original_slopes[,,,4,1]>original_slopes[,,,4,7])]=1

    signis[,,,13,1][which(MK[,,,5,2]<=0.05 & MK[,,,5,1]!=0)]=1
    signis[,,,16,2][which(MK[,,,2,2]<=0.05 & MK[,,,2,1]!=0)]=1

    values<-array(NA,c(3,5,regNumb,2,2))
    for (i in 1:length(trendIDs)){
        filename <- paste("../data/",dataset,additional_style,"/",trendIDs[i],"/regional/",region_name,"/",period,"/",trendIDs[i],"_",dataset,"_",region_name,"_",period,"_others.nc",sep="") ; print(filename)
        values[i,,,,1]<-var.get.nc(open.nc(filename),"other_stuff")[1:5,,,4]
        filename <- paste("../data/",dataset,additional_style,"/",trendIDs[i],"/regional/",region_name,"/",period,"/",trendIDs[i],"_",dataset,"_",region_name,"_",period,"_quantiles.nc",sep="") ; print(filename)
        values[i,,,,2]<-var.get.nc(open.nc(filename),"quantile_stuff")[1:5,,,2,2]
    }
    
    reihen<-array(NA,dim=c(5,regNumb,2,2))
    for (sea in season_auswahl){
        season<-season_names[sea]
        for (state in 1:2){
            reihen[sea,,state,1]=(values[1,sea,,state,1]-values[2,sea,,state,1])/values[2,sea,,state,1]*100
            reihen[sea,,state,2]=(values[1,sea,,state,2]-values[2,sea,,state,2])/values[2,sea,,state,2]*100
        }
    }
    reihen[which(is.na(reihen))]=0
    plot_reg_table_general(values=reihen,signis=signis,filename=paste("../plots/",dataset,additional_style,"/sensitivity/",region_name,"_sensitivity_",name_zusatz,"_",period,"_",name_style,"_5-7.pdf",sep=""),val_names=c("mn","95"),region_name="ward24",colorRange=c(-100,100),farb_palette="blau-rot",ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8),colorbar=FALSE)

    reihen<-array(NA,dim=c(5,regNumb,2,2))
    for (sea in season_auswahl){
        season<-season_names[sea]
        for (state in 1:2){
            reihen[sea,,state,1]=(values[3,sea,,state,1]-values[2,sea,,state,1])/values[2,sea,,state,1]*100
            reihen[sea,,state,2]=(values[3,sea,,state,2]-values[2,sea,,state,2])/values[2,sea,,state,2]*100
        }
    }
    reihen[which(is.na(reihen))]=0
    plot_reg_table_general(values=reihen,signis=signis,filename=paste("../plots/",dataset,additional_style,"/sensitivity/",region_name,"_sensitivity_",name_zusatz,"_",period,"_",name_style,"_9-7.pdf",sep=""),val_names=c("mn","95"),region_name="ward24",colorRange=c(-100,100),farb_palette="blau-rot",ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8))
}
