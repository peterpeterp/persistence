library(dgof)

ks_wilcoxon_period <- function(ID_name,yearPeriod1,yearPeriod2,periods,folder=paste("/regional/",ID_name,"/",sep=""),ID_select,ID_length=length(ID_select)){

    wilcox_test=array(NA,c(5,ID_length,2,3))
    ks_test=array(NA,c(5,ID_length,2,3))

    pdf(paste("../plots/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_ks_wilcox_",periods[1],"_vs_",periods[2],"_ecdf.pdf",sep=""),width=3,height=3)
    par(mar=c(3, 3, 3, 3) + 0.1)  

    for (sea in 1:5){ 
        season<-season_names[sea]
        cat(paste("\n",season))  

        print(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_duration_",season,".nc",sep=""))
        nc_dur<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_duration_",season,".nc",sep=""))
        dur<-var.get.nc(nc_dur,"dur")
        dur_mid<-var.get.nc(nc_dur,"dur_mid")

        percentage<-0
        cat(paste("\n0 -> -> -> -> -> 100\n"))
        for (q in ID_select){
            if (q/ID_length*100 > percentage){
                cat("-")
                percentage<-percentage+5
            }
            for (state in 1:2){
                ord<-order(dur_mid[q,state,])
                y<-as.vector(dur[q,state,][ord])
                x<-as.vector(dur_mid[q,state,][ord])
                inYearPeriod<-which(x>yearPeriod1[1] & x<yearPeriod1[2])
                y1<-y[inYearPeriod]
                inYearPeriod<-which(x>yearPeriod2[1] & x<yearPeriod2[2])
                y2<-y[inYearPeriod]
                ks_test[sea,q,state,1]=ks.test(y1,y2)$p.value
                wilcox_test[sea,q,state,1]=wilcox.test(x=y1,y=y2,paired=FALSE)$p.value
                plot.ecdf(y1,xlim=c(0,60),ylim=c(0,1),main="",cex=0.1,col="green",ylab="",xlab="",axes=TRUE,frame.plot=TRUE)
                par(new=TRUE)
                plot.ecdf(y2,xlim=c(0,60),ylim=c(0,1),main="",cex=0.1,col="orange",ylab="",xlab="",axes=TRUE,frame.plot=TRUE)
                text(40,0.25,paste("wilcox:",round(wilcox_test[sea,q,state,1],02)),cex=1)
                text(40,0.1,paste("ks:",round(ks_test[sea,q,state,1],02)),cex=1)
                text(40,0.6,paste(season_names[sea],q,state_names[state]),cex=1)
            }
        }
    }
    wilcox_test[,,,2][wilcox_test[,,,1]>0.05]=1
    wilcox_test[,,,2][wilcox_test[,,,1]<0.05]=0
    wilcox_test[,,,3][wilcox_test[,,,1]>0.1]=1
    wilcox_test[,,,3][wilcox_test[,,,1]<0.1]=0

    ks_test[,,,2][ks_test[,,,1]>0.05]=1
    ks_test[,,,2][ks_test[,,,1]<0.05]=0
    ks_test[,,,3][ks_test[,,,1]>0.1]=1
    ks_test[,,,3][ks_test[,,,1]<0.1]=0

    print(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_ks_wilcox_",yearPeriod1[1],"-",yearPeriod1[2],"_vs_",yearPeriod2[1],"-",yearPeriod2[2],".nc",sep=""))
    nc_out <- create.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_ks_wilcox_",yearPeriod1[1],"-",yearPeriod1[2],"_vs_",yearPeriod2[1],"-",yearPeriod2[2],".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "comment", "NC_CHAR", "wilcoxon test to identify changes in distribution")
    
    dim.def.nc(nc_out,"seasons",dimlength=5, unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)

    dim.def.nc(nc_out,"outs",dimlength=3,unlim=FALSE)

    var.def.nc(nc_out,"wilcox_test","NC_DOUBLE",c(0,1,2,3))
    att.put.nc(nc_out, "wilcox_test", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "wilcox_test", "dim_explanation", "NC_CHAR", "season-ID-states-...")
    att.put.nc(nc_out, "wilcox_test", "val_explanation", "NC_CHAR", "wilcoxon test p-value, p-value>0.05 = 1, p-value>0.1 = 1")

    var.def.nc(nc_out,"ks_test","NC_DOUBLE",c(0,1,2,3))
    att.put.nc(nc_out, "ks_test", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "ks_test", "dim_explanation", "NC_CHAR", "season-ID-states-...")
    att.put.nc(nc_out, "ks_test", "val_explanation", "NC_CHAR", "ks test p-value, p-value>0.05 = 1, p-value>0.1 = 1")

    var.put.nc(nc_out,"wilcox_test",wilcox_test)              
    var.put.nc(nc_out,"ks_test",ks_test)              

    close.nc(nc_out) 
    graphics.off()
}

distribution_comparision <- function(ID_name,periods,folder=paste("/regional/",ID_name,"/",sep=""),ID_select,ID_length=length(ID_select),region_names=1:23,hlines=c(30)){

   	filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_wilcox_",periods[1],"_vs_",periods[2],".nc",sep="") ; print(filename)
    wilcox_test <- var.get.nc(open.nc(filename),"wilcox_test")
    ks_test <- var.get.nc(open.nc(filename),"ks_test")

    fit_params=array(NA,c(2,6,ID_length,2,20))
    quantiles=array(NA,c(2,6,ID_length,2,4,3))
    distrs=array(NA,c(2,6,ID_length,2,5,100))
    for (i in 1:length(periods)){
    	period<-periods[i]
    	filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_fit_","2expo_4:100",".nc",sep=""); print(filename)
   	 	fit_params[i,,,,]=var.get.nc(open.nc(filename),"fit_stuff")

        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_quantiles",".nc",sep=""); print(filename)
        quantiles[i,,,,1:3,]=var.get.nc(open.nc(filename),"quantile_stuff")

        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_others",".nc",sep=""); print(filename)
        quantiles[i,,,,4,1]=var.get.nc(open.nc(filename),"other_stuff")[,,,1]
		
		filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_distributions.nc",sep=""); print(filename)
		distrs[i,,,,,]=var.get.nc(open.nc(filename),"distr_stuff")
	}

	pdf(paste("../plots/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_wilcox_",periods[1],"_vs_",periods[2],".pdf",sep=""),width=3,height=3)
    for (sea in 1:5){ 
        season<-season_names[sea]
        cat(paste("\n",season))  
        percentage<-0
        cat(paste("\n0 -> -> -> -> -> 100\n"))
        for (q in ID_select){
            if (q/ID_length*100 > percentage){
                cat("-")
                percentage<-percentage+5
            }
            for (state in 1:2){
            	fit_plot_comparison(distr=distrs[,sea,q,state,,],fits=fit_params[,sea,q,state,],sea,q,state,wilcox=wilcox_test[sea,q,state,1])
            }
        }
    }
    graphics.off()

    table<-file(paste("../plots/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_wilcox_",periods[1],"_vx_",periods[2],"_fit_diff.tex",sep=""))
    options(scipen=100)
    lines=c()
    index=0
    lines[index<-index+1]="\\documentclass[a4paper,12pt]{article}"
    lines[index<-index+1]="\\usepackage{xcolor,colortbl,pgf}"
    lines[index<-index+1]="\\usepackage{makecell}"
    lines[index<-index+1]="\\usepackage{geometry}"
    lines[index<-index+1]="\\geometry{ a4paper, total={190mm,288mm}, left=10mm, top=10mm, }"

    lines[index<-index+1]="\\definecolor{white}{rgb}{1,1,1}"
    lines[index<-index+1]="\\definecolor{green}{rgb}{0.5,1,0.5}"
    lines[index<-index+1]="\\definecolor{turkis}{rgb}{0.5,1,1}"
    lines[index<-index+1]="\\definecolor{violet}{rgb}{1,0.5,1}"

    lines[index<-index+1]="\\begin{document}"
    for (sea in 1:5){
        lines[index<-index+1]=paste("\\begin{table}[!h]")
        lines[index<-index+1]=paste("\\begin{tabular}{c||c||c|c|c|c||c||c||c|c|c|c}")

        lines[index<-index+1]=paste("\\multicolumn{12}{c}{",season_names[sea]," ",period," $",trendID,"$}\\","\\",sep="")
        lines[index<-index+1]=paste(" &\\multicolumn{5}{c}{",state_names[1],"}","& &\\multicolumn{5}{c}{",state_names[2],"}\\","\\",sep="")
        lines[index<-index+1]=paste("reg & $b_{expo}$  & b1 & b2 & thresh & dBIC & & $b_{expo}$ & b1 & b2 & thresh & dBIC  ","\\","\\",sep="")
        for (reg in ID_select){
            newline=paste(region_names[reg],sep="")
            for (state in 1:2){
                if (state==2){newline<-paste(newline," &",sep="")}

                for (i in c(2,6,8,9,17)){
                    if (ks_test[sea,reg,state,1]<=0.1){background<-"red!25"}
                    if (ks_test[sea,reg,state,1]<=0.05){background<-"red!50"}
                    if (ks_test[sea,reg,state,1]>0.1){background<-"white!10"}
                    newline<-paste(newline," &{\\cellcolor{",background,"}{",round(fit_params[2,sea,reg,state,i]-fit_params[1,sea,reg,state,i],02),"}}",sep="")
                }
                
            }
            newline<-paste(newline,paste("\\","\\",sep=""))
            lines[index<-index+1]=newline
            if (reg %in% hlines){
                lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
                lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
            }
        }
        lines[index<-index+1]=paste("\\end{tabular}")
        lines[index<-index+1]=paste("\\end{table}")
        lines[index<-index+1]=paste("\\vspace{0cm}")
      
    }
    lines[index<-index+1]="\\newpage"
    lines[index<-index+1]="\\fcolorbox{red!25}{red!25}{p.value(ks.test)$<$0.1}\\"
    lines[index<-index+1]="\\fcolorbox{red!50}{red!50}{p.value(ks.test)$<$0.05}\\"
    lines[index<-index+1]="\\end{document}"
    writeLines(lines, table)
    close(table)

    table<-file(paste("../plots/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_wilcox_",periods[1],"_vx_",periods[2],"_quantile_diff.tex",sep=""))
    options(scipen=100)
    lines=c()
    index=0
    lines[index<-index+1]="\\documentclass[a4paper,12pt]{article}"
    lines[index<-index+1]="\\usepackage{xcolor,colortbl,pgf}"
    lines[index<-index+1]="\\usepackage{makecell}"
    lines[index<-index+1]="\\usepackage{geometry}"
    lines[index<-index+1]="\\geometry{ a4paper, total={190mm,288mm}, left=10mm, top=10mm, }"

    lines[index<-index+1]="\\definecolor{white}{rgb}{1,1,1}"
    lines[index<-index+1]="\\definecolor{green}{rgb}{0.5,1,0.5}"
    lines[index<-index+1]="\\definecolor{turkis}{rgb}{0.5,1,1}"
    lines[index<-index+1]="\\definecolor{violet}{rgb}{1,0.5,1}"

    lines[index<-index+1]="\\begin{document}"
    for (sea in 1:5){
        lines[index<-index+1]=paste("\\begin{table}[!h]")
        lines[index<-index+1]=paste("\\begin{tabular}{c||c|c|c||c||c||c|c|c||c}")

        lines[index<-index+1]=paste("\\multicolumn{10}{c}{",season_names[sea]," ",period," $",trendID,"$}\\","\\",sep="")
        lines[index<-index+1]=paste(" &\\multicolumn{4}{c}{",state_names[1],"}","& &\\multicolumn{4}{c}{",state_names[2],"}\\","\\",sep="")
        lines[index<-index+1]=paste("reg & 75th & 95th & 99th & mean & & 75th & 95th & 99th & mean &","\\","\\",sep="")
        for (reg in ID_select){
            newline=paste(region_names[reg],sep="")
            for (state in 1:2){
                if (state==2){newline<-paste(newline," &",sep="")}

                for (i in c(1,2,3,4)){
                    val<-quantiles[2,sea,reg,state,i,1]-quantiles[1,sea,reg,state,i,1]
                    if (ks_test[sea,reg,state,1]>0.1 & val>0){background<-"violet!10"}
                    if (ks_test[sea,reg,state,1]<=0.1 & val>0){background<-"violet!25"}
                    if (ks_test[sea,reg,state,1]<=0.05 & val>0){background<-"violet!50"}

                    if (ks_test[sea,reg,state,1]>0.1 & val<0){background<-"turkis!10"}
                    if (ks_test[sea,reg,state,1]<=0.1 & val<0){background<-"turkis!25"}
                    if (ks_test[sea,reg,state,1]<=0.05 & val<0){background<-"turkis!50"}
                    newline<-paste(newline," &{\\cellcolor{",background,"}{",round(val,02),"}}",sep="")
                }
                
            }
            newline<-paste(newline,paste("\\","\\",sep=""))
            lines[index<-index+1]=newline
            if (reg %in% hlines){
                lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
                lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
            }
        }
        lines[index<-index+1]=paste("\\end{tabular}")
        lines[index<-index+1]=paste("\\end{table}")
        lines[index<-index+1]=paste("\\vspace{0cm}")
      
    }
    lines[index<-index+1]="\\newpage"
    lines[index<-index+1]="\\fcolorbox{violet!25}{violet!25}{p.value(ks.test)$>$0.1 and increase}\\"
    lines[index<-index+1]="\\fcolorbox{violet!25}{violet!25}{p.value(ks.test)$<$0.1 and increase}\\"
    lines[index<-index+1]="\\fcolorbox{violet!50}{violet!50}{p.value(ks.test)$<$0.05 and increase}\\"

    lines[index<-index+1]="\\fcolorbox{turkis!25}{turkis!25}{p.value(ks.test)$>$0.1 and decrease}\\"
    lines[index<-index+1]="\\fcolorbox{turkis!25}{turkis!25}{p.value(ks.test)$<$0.1 and decrease}\\"
    lines[index<-index+1]="\\fcolorbox{turkis!50}{turkis!50}{p.value(ks.test)$<$0.05 and decrease}\\"
    lines[index<-index+1]="\\end{document}"
    writeLines(lines, table)
    close(table)
}

ks_wilcoxon_trenID_sensitivity <- function(trendIDs=c("91_5","91_7","91_9")){
    wilcox_test=array(NA,c(2,3))
    ks_test=array(NA,c(2,3))

    pdf(paste("../plots/",dataset,additional_style,"/dur_ks_wilcox_sensibility_ecdf.pdf",sep=""),width=3,height=3)
    par(mar=c(3, 3, 3, 3) + 0.1)  

    filename<-paste("../data/",dataset,additional_style,"/",trendIDs[1],"/gridded/",trendIDs[1],dataset,"_duration_4seasons.nc",sep="") ; print(filename)
    dur5<-var.get.nc(open.nc(filename),"dur")

    filename<-paste("../data/",dataset,additional_style,"/",trendIDs[2],"/gridded/",trendIDs[2],dataset,"_duration_4seasons.nc",sep="") ; print(filename)
    dur7<-var.get.nc(open.nc(filename),"dur")

    filename<-paste("../data/",dataset,additional_style,"/",trendIDs[3],"/gridded/",trendIDs[3],dataset,"_duration_4seasons.nc",sep="") ; print(filename)
    dur9<-var.get.nc(open.nc(filename),"dur")

    colors=c(rgb(1,0.5,0.5,0.5),rgb(0.5,1,0.5,0.5),rgb(0.5,0.5,1,0.5))
    for (state in 1:2){
        y5<-as.vector(dur5[,state,])
        y7<-as.vector(dur7[,state,])
        y9<-as.vector(dur9[,state,])

        br<-0:max(c(max(y5,na.rm=TRUE),max(y7,na.rm=TRUE),max(y9,na.rm=TRUE)),na.rm=TRUE)
        tmp5<<-hist(y5,breaks=br,plot=FALSE)
        tmp7<<-hist(y7,breaks=br,plot=FALSE)
        tmp9<<-hist(y9,breaks=br,plot=FALSE)

        plot(NA,xlim=c(0,50),ylim=c(0,500000),main="",ylab="",xlab="")
        points(tmp5$mids,tmp5$count,col=rgb(1,0.5,0.5,0.5),cex=0.3)
        points(tmp7$mids,tmp7$count,col=rgb(0.5,1,0.5,0.5),cex=0.3)
        points(tmp9$mids,tmp9$count,col=rgb(0.5,0.5,1,0.5),cex=0.3)
        legend("topright",col=colors,legend=trendIDs)

        plot(NA,xlim=c(0,50),ylim=c(100,500000),main="",ylab="",xlab="",log="y")
        points(tmp5$mids[1:50],tmp5$count[1:50],col=rgb(1,0.5,0.5,0.5),cex=0.3)
        points(tmp7$mids[1:50],tmp7$count[1:50],col=rgb(0.5,1,0.5,0.5),cex=0.3)
        points(tmp9$mids[1:50],tmp9$count[1:50],col=rgb(0.5,0.5,1,0.5),cex=0.3)
        legend("topright",col=colors,legend=trendIDs)

        plot(NA,xlim=c(0,50),ylim=c(-500,15000),main="",ylab="",xlab="")
        points(tmp5$mids[1:50],tmp5$count[1:50]-tmp7$count[1:50],col=rgb(1,0.5,0.5,0.5),cex=0.3)
        points(tmp7$mids[1:50],tmp5$count[1:50]-tmp9$count[1:50],col=rgb(0.5,1,0.5,0.5),cex=0.3)
        points(tmp9$mids[1:50],tmp7$count[1:50]-tmp9$count[1:50],col=rgb(0.5,0.5,1,0.5),cex=0.3)
        legend("topright",col=colors,legend=c("5-7","5-9","7-9"))

        plot.ecdf(y5,xlim=c(0,60),ylim=c(0,1),main="",cex=0.1,col=colors[1],ylab="",xlab="",axes=TRUE,frame.plot=TRUE)
        par(new=TRUE)
        plot.ecdf(y7,xlim=c(0,60),ylim=c(0,1),main="",cex=0.1,col=colors[2],ylab="",xlab="",axes=TRUE,frame.plot=TRUE)
        par(new=TRUE)
        plot.ecdf(y9,xlim=c(0,60),ylim=c(0,1),main="",cex=0.1,col=colors[3],ylab="",xlab="",axes=TRUE,frame.plot=TRUE)

        graphics.off()

        adsa
        ks_test[state,1]=ks.test(y5,y7)$p.value
        wilcox_test[state,1]=wilcox.test(y5,y7)$p.value
        ks_test[state,2]=ks.test(y5,y9)$p.value
        wilcox_test[state,2]=wilcox.test(y5,y9)$p.value
        ks_test[state,3]=ks.test(y7,y9)$p.value
        wilcox_test[state,3]=wilcox.test(y5,y7)$p.value
    }


    graphics.off()
}



#ks_wilcoxon_trenID_sensitivity(trendIDs=c("91_5","91_7","91_9"))



#wilcoxon_period(ID_name="ward23",yearPeriod1=c(1950,1980),yearPeriod2=c(1980,2014),periods=c("1950-1980","1980-2014"),ID_select=1:23)

distribution_comparision(ID_name="ward23",periods=c("1950-1980","1980-2014"),ID_select=c(3,4,7,11,12,14,16,18,20,22),ID_length=23)