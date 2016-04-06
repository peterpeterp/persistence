library(Kendall)

duration_yearly_values <- function(folder="/gridded/",ID_name="",ID_select=1:length(dat$ID),ID_length=length(ID_select)){

    yearly_values=array(NA,c(5,ID_length,2,length(dat$year),7))

    for (sea in 1:5){
    	season<-season_names[sea]
        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_reg_binned_duration_",season,".nc",sep="") ; print(filename)
        nc<-open.nc(filename)
        binned_dur<-var.get.nc(nc,"binned_dur")

        percentage<-0
        cat(paste("\n0 -> -> -> -> -> 100\n"))
        for (q in ID_select){
            if (q/ID_length*100 > percentage){
                cat("-")
                percentage<-percentage+5
            }
        	for (state in 1:2){
        		for (yr in 1:length(dat$year)){
        			yearly_values[sea,q,state,yr,7]=length(which(!is.na(binned_dur[q,state,,yr])))
        			if (yearly_values[sea,q,state,yr,7]>10){
        				#print(binned_dur[q,state,,yr])
        				#print(paste(q,state,yr))
	        			yearly_values[sea,q,state,yr,5]=mean(binned_dur[q,state,,yr],na.rm=TRUE)
	        			yearly_values[sea,q,state,yr,1:length(taus)]=quantile_pete(binned_dur[q,state,,yr],taus=taus,na.rm=TRUE)
	        		}
        		}
        	}
        }
    }

    filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_yearly_values.nc",sep="") ; print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "taus", "NC_CHAR", paste(taus))
            
    dim.def.nc(nc_out,"seasons",dimlength=5, unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
    dim.def.nc(nc_out,"years",dimlength=length(dat$year),unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=7,unlim=FALSE)

    var.def.nc(nc_out,"yearly_values","NC_DOUBLE",c(0,1,2,3,4))
    att.put.nc(nc_out, "yearly_values", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "yearly_values", "dim_explanation", "NC_CHAR", "sea-ID-states-years-outs")
    att.put.nc(nc_out, "yearly_values", "val_explanation", "NC_CHAR", "annual: 7-# of periods in year, 5-mean(dur), 1:3-quantile(dur)")

    var.put.nc(nc_out,"yearly_values",yearly_values)              
    close.nc(nc_out) 
}

duration_MannKendall <- function(yearPeriod,folder="/gridded/",ID_name="",ID_select=1:length(dat$ID),ID_length=length(ID_select),region_names=1:ID_length,hlines=c(30)){
    filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_yearly_values.nc",sep="") ; print(filename)
    yearly_values<-var.get.nc(open.nc(filename),"yearly_values")

    MK=array(NA,c(5,ID_length,2,5,2))
    for (sea in 1:5){
    	for (q in ID_select){
    		for (state in 1:2){
    			for (out in c(1,2,3,5)){
    				tmp<-MannKendall(yearly_values[sea,q,state,(yearPeriod[1]-dat$year[1]+1):(yearPeriod[2]-dat$year[1]+1),out])
    				MK[sea,q,state,out,1]=tmp$tau
    				MK[sea,q,state,out,2]=tmp$sl
    			}
    		}
    	}
    }

    filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,dataset,"_",ID_name,"_",period,"_MK.nc",sep="") ; print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "taus", "NC_CHAR", paste(taus))
            
    dim.def.nc(nc_out,"seasons",dimlength=5, unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=5,unlim=FALSE)
    dim.def.nc(nc_out,"results",dimlength=2,unlim=FALSE)

    var.def.nc(nc_out,"MK","NC_DOUBLE",c(0,1,2,3,4))
    att.put.nc(nc_out, "MK", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "MK", "dim_explanation", "NC_CHAR", "sea-ID-states-(1:3-taus 5-mean)-results")
    att.put.nc(nc_out, "MK", "val_explanation", "NC_CHAR", "1-MK, 2-MK p-value")

    var.put.nc(nc_out,"MK",MK)              
    close.nc(nc_out) 

    filename<-paste("../plots/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,dataset,"_",ID_name,"_",period,"_MK.tex",sep="") ; print(filename)
    table<-file(filename)
    options(scipen=100) ; lines=c() ; index=0
    
    lines[index<-index+1]="\\documentclass[a4paper,12pt]{article}"
    lines[index<-index+1]="\\usepackage{xcolor,colortbl,pgf}"
    lines[index<-index+1]="\\usepackage{makecell}"
    lines[index<-index+1]="\\usepackage{geometry}"
    lines[index<-index+1]="\\geometry{ a4paper, total={190mm,288mm}, left=10mm, top=10mm, }"

    lines[index<-index+1]="\\definecolor{white}{rgb}{1,1,1}"
    lines[index<-index+1]="\\definecolor{red}{rgb}{1,0.3,1}"
    lines[index<-index+1]="\\definecolor{blue}{rgb}{0.3,1,1}"

    lines[index<-index+1]="\\begin{document}"
    lines[index<-index+1]="\\setlength{\\tabcolsep}{4pt}"

    lines[index<-index+1]=paste("\\begin{table}[!h]")
    lines[index<-index+1]=paste("\\begin{tabular}{c||cccc||cccc||cccc||cccc||cccc}")
    lines[index<-index+1]=paste("& \\multicolumn{4}{c}{MAM} & \\multicolumn{4}{c}{JJA} & \\multicolumn{4}{c}{SON} & \\multicolumn{4}{c}{DJF} & \\multicolumn{4}{c}{Annual}","\\\\")
    lines[index<-index+1]=paste("& \\multicolumn{2}{c}{cold} & \\multicolumn{2}{c}{warm}& \\multicolumn{2}{c}{cold} & \\multicolumn{2}{c}{warm} & \\multicolumn{2}{c}{cold} & \\multicolumn{2}{c}{warm}& \\multicolumn{2}{c}{cold} & \\multicolumn{2}{c}{warm} & \\multicolumn{2}{c}{cold} & \\multicolumn{2}{c}{warm}","\\\\")
    lines[index<-index+1]=paste("& m & 95 & m & 95 & m & 95 & m & 95 & m & 95 & m & 95 & m & 95 & m & 95 & m & 95 & m & 95","\\\\") 
    lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
    lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"

    for (reg in ID_select){
        newline<-paste(reg)
        for (sea in 1:5){
            for (state in 1:2){
                for (i in c(5,2)){
                    if (MK[sea,reg,state,i,1]>0){background_color<-"red"}
                    if (MK[sea,reg,state,i,1]<0){background_color<-"blue"}
                    if (MK[sea,reg,state,i,2]>0.1){background_intesity<-"!15"}
                    if (MK[sea,reg,state,i,2]<=0.1){background_intesity<-"!50"}
                    #if (MK[sea,reg,state,i,2]<=0.05){background_intesity<-"!75"}

                    if (MK[sea,reg,state,i,2]>0.1){newline<-paste(newline," &{\\cellcolor{",background_color,"!25}{ }}",sep="")}
                    if (MK[sea,reg,state,i,2]<=0.1){newline<-paste(newline," &{\\cellcolor{",background_color,"!75}{X}}",sep="")}

                    #newline<-paste(newline," &{\\cellcolor{",background_color,background_intesity,"}{",round(MK[sea,reg,state,i,1],03),"}}",sep="")
                }              
            }
        }
        lines[index<-index+1]=paste(newline,"\\\\")
        if (reg %in% hlines){
            lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
            lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
        }
    }
    lines[index<-index+1]=paste("\\end{tabular}")
    lines[index<-index+1]=paste("\\end{table}")
    lines[index<-index+1]=paste("\\vspace{0cm}")

    lines[index<-index+1]="\\newpage"

    lines[index<-index+1]="\\end{document}"
    writeLines(lines, table)
    close(table)
}