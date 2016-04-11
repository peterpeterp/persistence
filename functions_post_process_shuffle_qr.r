
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



plot_confi_intervals <- function(){
    print(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuffQuant.nc",sep=""))
    nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuffQuant.nc",sep=""))
    original_slopes=var.get.nc(nc,"original_slopes")
    confi_quantiles=var.get.nc(nc,"confi_quantiles")

    pdf(paste("../plots/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuffQuant.pdf",sep=""),width=19/4,height=29/4)
    par(mfrow=c(4,4),mar=c(2,2,2,0))
    yout=c(1,0.15,0.5,0.7,0.15)
    polygon_col=c(rgb(0.7,0.7,1,0.4),rgb(1,0.7,0.7,0.4))
    state_col=c(rgb(0.4,0.4,1),rgb(1,0.4,0.4))
    for (out in c(2,3,4,5)){
        for (sea in 1:4){
            plot(NA,xlim=c(1,length(reg_order)),ylim=c(-yout[out],yout[out]),main=paste(season_names[sea],out_names[out]),axes=FALSE)
            axis(1,at=1:length(reg_order),label=reg_order,las=2, cex.axis=0.8)
            axis(2)
            for (state in 1:2){
                polygon(x=c(1:length(reg_order),length(reg_order):1),y=c(confi_quantiles[sea,reg_order,state,out,1],confi_quantiles[sea,reg_order[length(reg_order):1],state,out,5]),col=polygon_col[state],border="white")
            }
            for (state in 1:2){
                points(1:length(reg_order),original_slopes[sea,reg_order,state,out,1],pch=20,cex=1,col=state_col[state])
            }
            abline(h=0)
            for (i in 1:length(reg_order)){if (reg_order[i] %in% hlines){abline(v=(i+0.5),lty=2)}}
        }
    }
    graphics.off()

}

confidence_interval <- function(seasons=1:5){
	confi_quantiles<-array(NA,dim=c(5,regNumb,2,5,5))
    original_slopes<-array(NA,dim=c(5,regNumb,2,5,5))
    for (sea in seasons){
        season<-season_names[sea]
        shuffled_mass<-array(NA,dim=c(10000,regNumb,2,5))
        original_mass<-array(NA,dim=c(100,regNumb,2,5))
        for (id in 1:100){
        	nc <- try(open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/shuffled/",trendID,dataset,"_",ID_name,"_",period,"_shuffled_trends_",season,"_",id,".nc",sep="")))
            if (class(nc)=="try-error"){print(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/shuffled/",trendID,dataset,"_",ID_name,"_",period,"_shuffled_trends_",season,"_",id,".nc",sep=""))}
            if (class(nc)!="try-error"){
	        	shuffled<-var.get.nc(nc,"shuffled")
	        	original<-var.get.nc(nc,"original")
	        	shuffled_mass[((id-1)*100+1):(id*100),,,]=shuffled
	        	original_mass[id,,,]=original
	        }
        }
        original_slopes[sea,,,,1]=original_mass[which(!is.na(original_mass[,1,1,1]))[1],,,]

        for (q in 1:regNumb){
            for (state in 1:2){
                for (out in c(1,2,3,4,5)){
                    confi_quantiles[sea,q,state,out,]=quantile(shuffled_mass[,q,state,out],c(0.025,0.05,0.5,0.95,0.975),na.rm=TRUE)
                    if (original_slopes[sea,q,state,out,1]>confi_quantiles[sea,q,state,out,5] | original_slopes[sea,q,state,out,1]<confi_quantiles[sea,q,state,out,1]){original_slopes[sea,q,state,out,2]=-0.05}
                    if (original_slopes[sea,q,state,out,1]>confi_quantiles[sea,q,state,out,4] | original_slopes[sea,q,state,out,1]<confi_quantiles[sea,q,state,out,2]){original_slopes[sea,q,state,out,3]=-0.1}
                }
            }
        }
    }
    print(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuffQuant.nc",sep=""))
    nc_out<-create.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuffQuant.nc",sep=""))

    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)
    
    dim.def.nc(nc_out,"seasons",dimlength=5,unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=regNumb, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)

    dim.def.nc(nc_out,"slopes",dimlength=5,unlim=FALSE)
    dim.def.nc(nc_out,"confi_quants",dimlength=5,unlim=FALSE)
        

    var.def.nc(nc_out,"confi_quantiles","NC_DOUBLE",c(0,1,2,3,4))
    att.put.nc(nc_out, "confi_quantiles", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "confi_quantiles", "dim_explanation", "NC_CHAR", "season-ID-state-...")
    att.put.nc(nc_out, "confi_quantiles", "explanation", "NC_CHAR", "(0.5,0.75,0.95,0.99,mean) x (shuflle_slope quantiles 0.025,0.05,0.5,0.95,0.975)")

    var.def.nc(nc_out,"original_slopes","NC_DOUBLE",c(0,1,2,3,4))
    att.put.nc(nc_out, "original_slopes", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "original_slopes", "dim_explanation", "NC_CHAR", "season-ID-state-...")
    att.put.nc(nc_out, "original_slopes", "explanation", "NC_CHAR", "(0.5,0.75,0.95,0.99,mean) x (slope,0.05 signi, 0.1 signi")
        
    var.put.nc(nc_out,"original_slopes",original_slopes) 
    var.put.nc(nc_out,"confi_quantiles",confi_quantiles) 

    close.nc(nc_out)     
}

write_slope_table <- function(){
    filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuffQuant.nc",sep="")
    original_slopes<-var.get.nc(open.nc(filename),"original_slopes")

    signis<-array(NA,c(5,regNumb,2,5,2))
    #signis[,,,1,1][which(is.na(original_slopes[,,,5,signi]))]=1
    signis[,,,4,1][which(!is.na(original_slopes[,,,5,3]))]=1
    signis[,,,5,1][which(!is.na(original_slopes[,,,5,2]))]=1

    #signis[,,,1,2][which(is.na(original_slopes[,,,3,signi]))]=1
    signis[,,,4,2][which(!is.na(original_slopes[,,,3,3]))]=1
    signis[,,,5,2][which(!is.na(original_slopes[,,,3,2]))]=1

    values<-array(NA,c(5,regNumb,2,2))
    values[,,,1]=original_slopes[,,,5,1]
    values[,,,2]=original_slopes[,,,3,1]

    nbcol<<-101
    plot_reg_table_general(values=values,signis=signis,filename=paste("../plots/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,dataset,"_",ID_name,"_",period,"_slope.pdf",sep=""),val_names=c("mn","95"),region_name="ward24",colorRange=c(-0.1,0.1),farb_palette="lila-gruen",ID_select=reg_order,hlines=hlines)
}


init <- function(){
    library(quantreg)
    library(RNetCDF)

    period<<-"1950-2014"
    trendID<<-"91_7"

    dataset<<-"_TMean"
    additional_style<<-""

    state_names<<-c("cold","warm")
    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    out_names<<-c(0.5,0.75,0.95,0.99,"mean")

    ID_name<<-"ward24"
    folder<<-paste("/regional/",ID_name,"/",sep="")
    regNumb<<-24
    region_names<<-1:24

    plot_select<<-c(3,4,5,7,11,12,13,14,16,18,20,21)
    reg_order<<-c(1,2,6,10,13,19,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,23)
    #reg_order<<-c(3,4,5,7,11,12,14,16,18,20,22)
    #reg_order<<-c(3,4,7,12,16,20)
    hlines<<-c(30)
    hlines<<-c(19,20,22,8)
    #plot_select<<-c(11,12,16,20)
    plotNumb<<-length(plot_select)
}

init()
#confidence_interval()
#write_robustness_table()
write_slope_table()

for (trendID in c("91_7","91_5","91_9")){
    #plot_confi_intervals()
    #write_slope_table(signi=2)
}


