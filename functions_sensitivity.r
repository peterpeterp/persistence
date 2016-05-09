

sens_climatology <- function(trendIDs=c("91_5","91_7","91_9"),period="1950-2014",file="_others",var="other_stuff",sub_auswahl=NA,value_auswahl=1,name_zusatz="mean",farb_mitte="mean",farb_palette="weiss-rot",signi_level=0.05){

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

	filename_plot=paste("../plots/",dataset,additional_style,"/sensitivity/","sensitivity_",name_zusatz,name_reg_zusatz,"_",period,"_",name_style,".pdf",sep="")
	print(filename_plot)
	topo_map_plot(filename_plot=filename_plot,reihen=reihen,reihen_sig=reihen_sig,farb_mitte=farb_mitte,farb_palette=farb_palette,signi_level=signi_level)

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