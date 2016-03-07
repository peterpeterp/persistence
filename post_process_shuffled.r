
confidence_interval <- function(seasons=1:4){
	confi_quantiles<-array(NA,dim=c(5,regNumb,2,5,5))
    original_slopes<-array(NA,dim=c(5,regNumb,2,5,5))
    for (sea in seasons){
        season<-season_names[sea]
        shuffled_mass<-array(NA,dim=c(10000,regNumb,2,5))
        original_mass<-array(NA,dim=c(100,regNumb,2,5))
        for (id in 1:100){
        	nc <- try(open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/shuffled/",trendID,dataset,"_",ID_name,"_",period,"_shuffled_trends_",season,"_",id,".nc",sep="")))
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
    print(paste("../plots/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuffQuant.tex",sep=""))
    table<-file(paste("../plots/",dataset,additional_style,"/",trendID,folder,trendID,"_",dataset,"_",ID_name,"_shuffQuant.tex",sep=""))
    options(scipen=100)

    lines=c()
    index=0

    lines[index<-index+1]="\\documentclass[a4paper,12pt]{article}"
    lines[index<-index+1]="\\usepackage{xcolor,colortbl,pgf}"
    lines[index<-index+1]="\\usepackage{makecell}"
    lines[index<-index+1]="\\usepackage{geometry}"
    lines[index<-index+1]="\\geometry{ a4paper, total={190mm,288mm}, left=20mm, top=20mm, }"
    lines[index<-index+1]="\\definecolor{white}{rgb}{1,1,1}"
    lines[index<-index+1]="\\definecolor{red}{rgb}{1,0.3,1}"
    lines[index<-index+1]="\\definecolor{blue}{rgb}{0.3,1,1}"
    lines[index<-index+1]="\\begin{document}"

    for (period in c("1950-2014","1980-2014")){
        print(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuffQuant.nc",sep=""))
        nc=try(open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuffQuant.nc",sep="")),silent=TRUE)
        if (class(nc)!="try-error"){
            original_slopes=var.get.nc(nc,"original_slopes")
            confi_quantiles=var.get.nc(nc,"confi_quantiles")
            lines[index<-index+1]=paste("\\begin{table}[!h]")
            lines[index<-index+1]=paste("\\begin{tabular}{c||cc||cc||cc||cc}")

            lines[index<-index+1]=paste(period,"& \\multicolumn{2}{c}{MAM} & \\multicolumn{2}{c}{JJA} & \\multicolumn{2}{c}{SON} & \\multicolumn{2}{c}{DJF}","\\","\\",sep="")
            lines[index<-index+1]=paste("$",trendID,"$","& cold & warm & cold & warm & cold & warm & cold & warm","\\","\\",sep="") 
            lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
            lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"

            for (reg in reg_order){
                newline<-paste(reg)
                for (sea in 1:4){
                    for (state in 1:2){
                        if (original_slopes[sea,reg,state,3,1]>0 & !is.na(original_slopes[sea,reg,state,3,3])){color<-"red!70"}
                        if (original_slopes[sea,reg,state,3,1]>0 & is.na(original_slopes[sea,reg,state,3,3])){color<-"red!25"}
                        if (original_slopes[sea,reg,state,3,1]<0 & !is.na(original_slopes[sea,reg,state,3,3])){color<-"blue!70"}
                        if (original_slopes[sea,reg,state,3,1]<0 & is.na(original_slopes[sea,reg,state,3,3])){color<-"blue!25"}
                        if (abs(original_slopes[sea,reg,state,3,1])<0.0001){color<-"white!25"}
           
                        newline<-paste(newline,"&{\\cellcolor{",color,"}{",round(original_slopes[sea,reg,state,3,1],03),"}}",sep="") #{\\cellcolor{",color,"}}
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
    par(mfrow=c(4,3),mar=c(2,2,2,0))
    yout=c(1,1,0.5,0.7,0.15)
    polygon_col=c(rgb(0.7,0.7,1,0.4),rgb(1,0.7,0.7,0.4))
    state_col=c(rgb(0.4,0.4,1),rgb(1,0.4,0.4))
    for (sea in 1:4){
        for (out in c(3,4,5)){
        
            plot(NA,xlim=c(1,length(reg_order)),ylim=c(-yout[out],yout[out]),main=paste(season_names[sea],out_names[out]),axes=FALSE)
            axis(1,at=1:length(reg_order),label=reg_order)
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

    period<<-"1950-2014"
    ID_name<<-"ward23"
    folder<<-paste("/regional/",ID_name,"/",sep="")
    regNumb<<-23
    region_names<<-1:23

    plot_select<<-c(3,4,5,7,11,12,13,14,16,18,20,21)
    reg_order<<-c(1,2,6,10,19,3,4,7,12,13,16,20,5,11,14,18,21,22,17,8,9,15,23)
    hlines<<-c(19,20,22,8)
    #plot_select<<-c(11,12,16,20)
    plotNumb<<-length(plot_select)
}

init()
#confidence_interval()

plot_confi_intervals()

write_slope_table()
