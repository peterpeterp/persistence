

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

    indexTopRight<<-c("a","A","a","A","a","A","a","A","a","A")
    indexBottomLeft<<-c("cold","warm","cold","warm","cold","warm","cold","warm","cold","warm")
	filename_plot<-paste("../plots/",dataset,additional_style,"/sensitivity/","sensitivity_",name_zusatz,"_",period,"_",name_style,"_5-7.pdf",sep="") ; print(filename_plot)
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

    indexTopRight<<-c("b","B","b","B","b","B","b","B","b","B")
    indexBottomLeft<<-c("cold","warm","cold","warm","cold","warm","cold","warm","cold","warm")
    filename_plot<-paste("../plots/",dataset,additional_style,"/sensitivity/","sensitivity_",name_zusatz,"_",period,"_",name_style,"_9-7.pdf",sep="") ; print(filename_plot)
    topo_map_plot(filename_plot=filename_plot,reihen=reihen,reihen_sig=reihen_sig,farb_mitte=farb_mitte,farb_palette=farb_palette,signi_level=signi_level)

    indexTopRight<<-c("a","A","b","B","c","C","d","D","e","E")
    indexBottomLeft<<-c("MAM","MAM","JJA","JJA","SON","SON","DJF","DJF","Annual","Annual")
}


sens_regional_fits <- function(region_name="ward24",regNumb=24,trendIDs=c("91_5","91_7","91_9"),period="1950-2014",name_zusatz="fit",farb_mitte="",farb_palette="weiss-rot",signi_level=0.05){

	values<-array(NA,c(3,5,regNumb,2,30))
	for (i in 1:length(trendIDs)){
		filename <- paste("../data/",dataset,additional_style,"/",trendIDs[i],"/regional/",region_name,"/",period,"/",trendIDs[i],"_",dataset,"_",region_name,"_",period,"_fit_2expo_4:100.nc",sep="") ; print(filename)
        values[i,,,,]<-var.get.nc(open.nc(filename),"fit_stuff")[1:5,,,]
	}
	
    signis<-array(NA,c(5,regNumb,2,5,2))

    signis[,,,4,1][which(values[1,1:5,,,24]>0)]=1
    signis[,,,4,2][which(values[1,1:5,,,24]>0)]=1
    signis[,,,2,1][which(values[2,1:5,,,24]>0)]=1
    signis[,,,2,2][which(values[2,1:5,,,24]>0)]=1
    signis[,,,5,1][which(values[3,1:5,,,24]>0)]=1
    signis[,,,5,2][which(values[3,1:5,,,24]>0)]=1

	reihen<-array(NA,dim=c(5,regNumb,2,2))
    for (sea in season_auswahl){
        season<-season_names[sea]
        for (state in 1:2){
            reihen[sea,,state,1]=(values[1,sea,,state,12]-values[2,sea,,state,12])/values[2,sea,,state,12]*100
            reihen[sea,,state,2]=(values[1,sea,,state,14]-values[2,sea,,state,14])/values[2,sea,,state,14]*100
        }
    }

    nbcol<<-101
    reihen[which(is.na(reihen))]=0
    plot_reg_table_general(values=reihen,signis=signis,filename=paste("../plots/",dataset,additional_style,"/sensitivity/",region_name,"_sensitivity_",name_zusatz,"_",period,"_slopes_5-7.pdf",sep=""),val_names=c("b1","b2"),region_name="ward24",colorRange=c(-10,10),farb_palette="blau-rot",ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8),colorbar=FALSE)

    reihen<-array(NA,dim=c(5,regNumb,2,2))
    for (sea in season_auswahl){
        season<-season_names[sea]
        for (state in 1:2){
            reihen[sea,,state,1]=(values[3,sea,,state,12]-values[2,sea,,state,12])/values[2,sea,,state,12]*100
            reihen[sea,,state,2]=(values[3,sea,,state,14]-values[2,sea,,state,14])/values[2,sea,,state,14]*100
        }
    }

    reihen[which(is.na(reihen))]=0
    plot_reg_table_general(values=reihen,signis=signis,filename=paste("../plots/",dataset,additional_style,"/sensitivity/",region_name,"_sensitivity_",name_zusatz,"_",period,"_slopes_9-7.pdf",sep=""),val_names=c("b1","b2"),region_name="ward24",colorRange=c(-10,10),farb_palette="blau-rot",ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8))

    reihen<-array(NA,dim=c(5,regNumb,2,1))
    for (sea in season_auswahl){
        season<-season_names[sea]
        for (state in 1:2){
            reihen[sea,,state,1]=(values[1,sea,,state,15]-values[2,sea,,state,15])/values[2,sea,,state,15]*100
        }
    }
    plot_reg_table_general(values=reihen,signis=signis,filename=paste("../plots/",dataset,additional_style,"/sensitivity/",region_name,"_sensitivity_",name_zusatz,"_",period,"_tresh_5-7.pdf",sep=""),val_names=c(""),region_name="ward24",colorRange=c(-10,10),farb_palette="blau-rot",ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8),colorbar=FALSE)

    reihen<-array(NA,dim=c(5,regNumb,2,1))
    for (sea in season_auswahl){
        season<-season_names[sea]
        for (state in 1:2){
            reihen[sea,,state,1]=(values[3,sea,,state,15]-values[2,sea,,state,15])/values[2,sea,,state,15]*100
        }
    }
    plot_reg_table_general(values=reihen,signis=signis,filename=paste("../plots/",dataset,additional_style,"/sensitivity/",region_name,"_sensitivity_",name_zusatz,"_",period,"_tresh_9-7.pdf",sep=""),val_names=c(""),region_name="ward24",colorRange=c(-10,10),farb_palette="blau-rot",ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8))
}

sens_regional_trends <- function(region_name="ward24",regNumb=24,trendIDs=c("91_5","91_7","91_9"),period="1950-2014",name_zusatz="mean",farb_mitte="",farb_palette="weiss-rot",signi_level=0.05){

    filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuffQuant.nc",sep="")
    original_slopes<-var.get.nc(open.nc(filename),"original_slopes")

    signis<-array(NA,c(5,regNumb,2,5,2))
    signis[,,,4,1][which(!is.na(original_slopes[,,,5,3]))]=1
    signis[,,,5,1][which(!is.na(original_slopes[,,,5,2]))]=1

    signis[,,,4,2][which(!is.na(original_slopes[,,,3,3]))]=1
    signis[,,,5,2][which(!is.na(original_slopes[,,,3,2]))]=1



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
