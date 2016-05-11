

sens_gridded <- function(trendIDs=c("91_5","91_7","91_9"),period="1950-2014",file="_others",var="other_stuff",sub_auswahl=NA,value_auswahl=1,name_zusatz="mean",farb_mitte="mean",farb_palette="weiss-rot",signi_level=0.05){

	values<-array(NA,c(3,5,ntot,2))
	for (i in 1:length(trendIDs)){
		filename <- paste("../data/",dataset,additional_style,"/",trendIDs[i],"/gridded/",period,"/",trendIDs[i],"_",dataset,"_",period,file,".nc",sep="") ; print(filename)
		#print(dim(values))
		#print(dim(var.get.nc(open.nc(filename),var)))
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
			mn<-abs(colMeans(values[,sea,,state]))
			reihen[index,]=sqrt(1/3*(values[1,sea,,state]-mn)^2+(values[2,sea,,state]-mn)^2+(values[3,sea,,state]-mn)^2)/mn*100
			#reihen[index,]=0
			#reihen[index,][sign(values[1,sea,,state])==sign(values[2,sea,,state]) & sign(values[3,sea,,state])==sign(values[2,sea,,state])]=1
			#reihen[index,2:4]=2:4

		}
	}

	values<<-values
	reihen<<-reihen

	filename_plot=paste("../plots/",dataset,additional_style,"/sensitivity/","sensitivity_",name_zusatz,"_",period,"_",name_style,".pdf",sep="")
	print(filename_plot)
	topo_map_plot(filename_plot=filename_plot,reihen=reihen,reihen_sig=reihen_sig,farb_mitte=farb_mitte,farb_palette=farb_palette,signi_level=signi_level)
}


sens_regional_fits <- function(region_name="ward24",regNumb=24,trendIDs=c("91_5","91_7","91_9"),period="1950-2014",name_zusatz="fit",farb_mitte="",farb_palette="weiss-rot",signi_level=0.05){

	values<-array(NA,c(3,5,regNumb,2,30))
	for (i in 1:length(trendIDs)){
		filename <- paste("../data/",dataset,additional_style,"/",trendIDs[i],"/regional/",region_name,"/",period,"/",trendIDs[i],"_",dataset,"_",region_name,"_",period,"_fit_2expo_4:100.nc",sep="") ; print(filename)
		print(dim(var.get.nc(open.nc(filename),"fit_stuff")))
		print(dim(values))
        values[i,,,,]<-var.get.nc(open.nc(filename),"fit_stuff")[1:5,,,]
	}
	
	reihen<-array(NA,dim=c(5,regNumb,2,2))
	index<-0
	for (sea in season_auswahl){
		season<-season_names[sea]
		for (state in 1:2){
			index<-index+1
			mn<-abs(colMeans(values[,sea,,state,12]))
			reihen[sea,,state,1]=sqrt(1/3*((values[1,sea,,state,12]-mn)^2+(values[2,sea,,state,12]-mn)^2+(values[3,sea,,state,12]-mn)^2))/mn*100
			mn<-abs(colMeans(values[,sea,,state,14]))
			reihen[sea,,state,2]=sqrt(1/3*((values[1,sea,,state,14]-mn)^2+(values[2,sea,,state,14]-mn)^2+(values[3,sea,,state,14]-mn)^2))/mn*100
		}
	}
    signis<-array(NA,c(5,regNumb,2,5,2))

    signis[,,,4,1][which(values[1,1:5,,,24]>0)]=1
    signis[,,,4,2][which(values[1,1:5,,,24]>0)]=1
    signis[,,,2,1][which(values[2,1:5,,,24]>0)]=1
    signis[,,,2,2][which(values[2,1:5,,,24]>0)]=1
    signis[,,,5,1][which(values[3,1:5,,,24]>0)]=1
    signis[,,,5,2][which(values[3,1:5,,,24]>0)]=1

    nbcol<<-101
    reihen[which(is.na(reihen))]=0
    plot_reg_table_general(values=reihen,signis=signis,filename=paste("../plots/",dataset,additional_style,"/sensitivity/",region_name,"_sensitivity_",name_zusatz,"_",period,"_slopes.pdf",sep=""),val_names=c("b1","b2"),region_name="ward24",colorRange=c(0,10),farb_palette="weiss-rot",ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8))

    reihen<-array(NA,dim=c(5,regNumb,2,1))
    index<-0
    for (sea in season_auswahl){
        season<-season_names[sea]
        for (state in 1:2){
            index<-index+1
            mn<-abs(colMeans(values[,sea,,state,15]))
            reihen[sea,,state,1]=sqrt(1/3*((values[1,sea,,state,15]-mn)^2+(values[2,sea,,state,15]-mn)^2+(values[3,sea,,state,15]-mn)^2))/mn*100
        }
    }
    plot_reg_table_general(values=reihen,signis=signis,filename=paste("../plots/",dataset,additional_style,"/sensitivity/",region_name,"_sensitivity_",name_zusatz,"_",period,"_tresh.pdf",sep=""),val_names=c(""),region_name="ward24",colorRange=c(0,10),farb_palette="weiss-rot",ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8))
}

sens_regional_trends <- function(region_name="ward24",regNumb=24,trendIDs=c("91_5","91_7","91_9"),period="1950-2014",name_zusatz="mean",farb_mitte="",farb_palette="weiss-rot",signi_level=0.05){

    values<-array(NA,c(3,5,regNumb,2,2))
    for (i in 1:length(trendIDs)){
        filename <- paste("../data/",dataset,additional_style,"/",trendIDs[i],"/regional/",region_name,"/",period,"/",trendIDs[i],"_",dataset,"_",region_name,"_",period,"_others.nc",sep="") ; print(filename)
        print(dim(var.get.nc(open.nc(filename),"other_stuff")))
        print(dim(values))
        values[i,,,,1]<-var.get.nc(open.nc(filename),"other_stuff")[1:5,,,4]
        filename <- paste("../data/",dataset,additional_style,"/",trendIDs[i],"/regional/",region_name,"/",period,"/",trendIDs[i],"_",dataset,"_",region_name,"_",period,"_quantiles.nc",sep="") ; print(filename)
        print(dim(var.get.nc(open.nc(filename),"quantile_stuff")))
        print(var.get.nc(open.nc(filename),"quantile_stuff")[2,12,2,,])
        values[i,,,,2]<-var.get.nc(open.nc(filename),"quantile_stuff")[1:5,,,2,2]
    }
    
    reihen<-array(NA,dim=c(5,regNumb,2,2))
    index<-0
    for (sea in season_auswahl){
        season<-season_names[sea]
        for (state in 1:2){
            index<-index+1
            mn<-abs(colMeans(values[,sea,,state,1]))
            reihen[sea,,state,1]=sqrt(1/3*((values[1,sea,,state,1]-mn)^2+(values[2,sea,,state,1]-mn)^2+(values[3,sea,,state,1]-mn)^2))/mn*100
            mn<-abs(colMeans(values[,sea,,state,2]))
            reihen[sea,,state,2]=sqrt(1/3*((values[1,sea,,state,2]-mn)^2+(values[2,sea,,state,2]-mn)^2+(values[3,sea,,state,2]-mn)^2))/mn*100
        }
    }
    signis<-array(NA,c(5,regNumb,2,5,2))
    nbcol<<-101
    reihen[which(is.na(reihen))]=0
    plot_reg_table_general(values=reihen,signis=signis,filename=paste("../plots/",dataset,additional_style,"/sensitivity/",region_name,"_sensitivity_",name_zusatz,"_",period,"_",name_style,".pdf",sep=""),val_names=c("mn","95"),region_name="ward24",colorRange=c(0,300),farb_palette="weiss-rot",ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24),hlines=c(23,20,22,8))
}

write_robustness_table <- function(){
    print(paste("../plots/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuffQuant.tex",sep=""))
    table<-file(paste("../plots/",dataset,additional_style,"/slope_robustness_shuffQuant.tex",sep=""))
    options(scipen=100)

    lines=c()
    index=0

    lines[index<-index+1]="\\documentclass[a4paper,12pt]{article}"
    lines[index<-index+1]="\\usepackage{xcolor,colortbl,pgf}"
    lines[index<-index+1]="\\usepackage{makecell}"
    lines[index<-index+1]="\\usepackage{geometry}"
    lines[index<-index+1]="\\geometry{ a4paper, total={190mm,288mm}, left=20mm, top=10mm, }"
    lines[index<-index+1]="\\definecolor{white}{rgb}{1,1,1}"
    lines[index<-index+1]="\\definecolor{red}{rgb}{1,0.3,1}"
    lines[index<-index+1]="\\definecolor{blue}{rgb}{0.3,1,1}"
    lines[index<-index+1]="\\begin{document}"


    slope_mass=array(NA,c(3,5,23,2,5,5))
    slope_robustness=array(NA,c(5,23,2,5,5))
    trendIDs=c("91_5","91_7","91_9")

    for (period in c("1950-2014","1980-2014")){
        for (i in 1:3){
            trendID=trendIDs[i]
            print(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuffQuant.nc",sep=""))
            nc=try(open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuffQuant.nc",sep="")),silent=TRUE)
            if (class(nc)!="try-error"){slope_mass[i,,,,,]=var.get.nc(nc,"original_slopes")}
        }
        for (sea in 1:4){
            for (reg in reg_order){
                for (state in 1:2){
                    for (out in 1:5){
                        slope_robustness[sea,reg,state,out,1]=mean(slope_mass[,sea,reg,state,out,1],na.rm=TRUE)
                        slope_robustness[sea,reg,state,out,2]=sd(slope_mass[,sea,reg,state,out,1],na.rm=TRUE)
                        slope_robustness[sea,reg,state,out,3]=sd(slope_mass[,sea,reg,state,out,1],na.rm=TRUE)/mean(slope_mass[,sea,reg,state,out,1],na.rm=TRUE)*100
                        if (abs(sum(slope_mass[,sea,reg,state,out,1],na.rm=TRUE))<sum(abs(slope_mass[,sea,reg,state,out,1]),na.rm=TRUE)){slope_robustness[sea,reg,state,out,4]=-1}
                        if (abs(sum(slope_mass[,sea,reg,state,out,1],na.rm=TRUE))==sum(abs(slope_mass[,sea,reg,state,out,1]),na.rm=TRUE)){slope_robustness[sea,reg,state,out,4]=1}

                    }
                }
            }
        }

        for (out in c(3,5)){
            lines[index<-index+1]=paste("\\begin{table}[!h]")
            lines[index<-index+1]=paste("\\begin{tabular}{c||cc||cc||cc||cc}")

            lines[index<-index+1]=paste(period,"& \\multicolumn{2}{c}{MAM} & \\multicolumn{2}{c}{JJA} & \\multicolumn{2}{c}{SON} & \\multicolumn{2}{c}{DJF}","\\","\\",sep="")
            lines[index<-index+1]=paste(out_names[out],"& cold & warm & cold & warm & cold & warm & cold & warm","\\","\\",sep="") 
            lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
            lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"

            for (reg in reg_order){
                newline<-paste(reg)
                for (sea in 1:4){
                    for (state in 1:2){
                        if (slope_robustness[sea,reg,state,out,4]>0){color<-"white!70"}
                        if (slope_robustness[sea,reg,state,out,4]<0){color<-"red!100"}               
                        newline<-paste(newline,"&{\\cellcolor{",color,"}{",round(slope_robustness[sea,reg,state,out,3]),"}}",sep="") #{\\cellcolor{",color,"}}
                    }
                }
                lines[index<-index+1]=paste(newline,"\\","\\",sep="")
                if (reg %in% hlines){
                    lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
                    lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
                }
            }
            lines[index<-index+1]=paste("\\end{tabular}")
            lines[index<-index+1]=paste("\\end{table}")
        }
    }

    lines[index<-index+1]="\\end{document}"
    writeLines(lines, table)
    close(table)
}