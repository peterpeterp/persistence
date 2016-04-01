


write_tex_table<-function(filename,outs,values,signis,header){

    table<-file(filename)
    options(scipen=100)
    lines=c()
    index=0
    lines[index<-index+1]="\\documentclass[a4paper,12pt]{article}"
    lines[index<-index+1]="\\usepackage{xcolor,colortbl,pgf}"
    lines[index<-index+1]="\\usepackage{makecell}"
    lines[index<-index+1]="\\usepackage[top=3cm, bottom=3cm, left=3cm, right=3cm]{geometry}"

    lines[index<-index+1]="\\definecolor{white}{rgb}{1,1,1}"
    lines[index<-index+1]="\\definecolor{green}{rgb}{0.5,1,0.5}"
    lines[index<-index+1]="\\definecolor{blue}{rgb}{0.3,1,1}"
    lines[index<-index+1]="\\definecolor{red}{rgb}{1,0.3,1}"

    lines[index<-index+1]="\\begin{document}"
    lines[index<-index+1]="\\setlength{\\tabcolsep}{4pt}"
    lines[index<-index+1]=paste("\\begin{table}[!h]")
    lines[index<-index+1]=paste("\\begin{tabular}{c||cccc||cccc||cccc||cccc}")
    lines[index<-index+1]=paste("& \\multicolumn{4}{c}{MAM} & \\multicolumn{4}{c}{JJA} & \\multicolumn{4}{c}{SON} & \\multicolumn{4}{c}{DJF}","\\\\")
    lines[index<-index+1]=paste("& \\multicolumn{2}{c}{cold} & \\multicolumn{2}{c}{warm}& \\multicolumn{2}{c}{cold} & \\multicolumn{2}{c}{warm} & \\multicolumn{2}{c}{cold} & \\multicolumn{2}{c}{warm}& \\multicolumn{2}{c}{cold} & \\multicolumn{2}{c}{warm}","\\\\")
    lines[index<-index+1]=header 
    lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
    lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"

    for (reg in ID_select){
        newline<-paste(reg)
        for (sea in 1:4){
            for (state in 1:2){
                for (i in outs){
                    val<-values[sea,reg,state,i]
                    sig<-signis[sea,reg,state,i]

                    if (val>0){background_color<-"red"}
                    if (val<0){background_color<-"blue"}
                    if (sig>0.1){newline<-paste(newline," &{\\cellcolor{",background_color,"!25}{ }}",sep="")}
                    if (sig<=0.1){newline<-paste(newline," &{\\cellcolor{",background_color,"!75}{\\footnotesize{",round(val,01),"}}}",sep="")}
                }
                
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
      
    lines[index<-index+1]="\\end{document}"
    writeLines(lines, table)
    close(table)
}

write_regional_fit_table <- function(region_name="srex",fit_style,region_names,ID_select,ID_length=length(ID_select),hlines=c(30)){
    print(paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style,".nc",sep=""))
    nc = open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style,".nc",sep=""))
    fit_stuff<-var.get.nc(nc,"fit_stuff")


    table<-file(paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",region_name,"_",period,"_fit_",fit_style,".tex",sep=""))
    options(scipen=100)

    lines=c()
    index=0

    lines[index<-index+1]="\\documentclass[a4paper,12pt]{article}"
    lines[index<-index+1]="\\usepackage{xcolor,colortbl,pgf}"
    lines[index<-index+1]="\\usepackage{makecell}"
    lines[index<-index+1]="\\usepackage[top=3cm, bottom=3cm, left=3cm, right=3cm]{geometry}"

    lines[index<-index+1]="\\definecolor{white}{rgb}{1,1,1}"
    lines[index<-index+1]="\\definecolor{green}{rgb}{0.5,1,0.5}"
    lines[index<-index+1]="\\definecolor{blue}{rgb}{0.3,1,1}"
    lines[index<-index+1]="\\definecolor{red}{rgb}{1,0.3,1}"

    lines[index<-index+1]="\\begin{document}"
    lines[index<-index+1]="\\setlength{\\tabcolsep}{4pt}"

    for (sea in 1:5){
        lines[index<-index+1]=paste("\\begin{table}[!h]")
        lines[index<-index+1]=paste("\\begin{tabular}{c||c||c|c|c|c||c||c|c|c|c}")

        lines[index<-index+1]=paste("\\multicolumn{11}{c}{",season_names[sea]," ",period," $",trendID,"$}\\","\\",sep="")
        lines[index<-index+1]=paste(" &\\multicolumn{5}{c}{",state_names[1],"}","&\\multicolumn{5}{c}{",state_names[2],"}\\","\\",sep="")
        lines[index<-index+1]=paste("reg & $b_{expo}$ & b1 & b2 & thresh &b1-b2& $b_{expo}$ & b1 & b2 & thresh & b1-b2","\\","\\",sep="")
        lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
        lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
        for (reg in ID_select){
            newline<-paste(region_names[reg],sep="")
            for (state in 1:2){
                for (i in c(2)){
                    if (fit_stuff[sea,reg,state,24]>0){background_color<-"green"}
                    if (fit_stuff[sea,reg,state,24]<0){background_color<-"white"}
                    if (fit_stuff[sea,reg,state,8]<0.99){newline<-paste(newline," &{\\cellcolor{",background_color,"!25}{ }}",sep="")}
                    if (fit_stuff[sea,reg,state,8]>=0.99){newline<-paste(newline," &{\\cellcolor{",background_color,"!75}{",round(fit_stuff[sea,reg,state,i],02),"}}",sep="")}
                }
                for (i in c(12,14,15)){
                    if (fit_stuff[sea,reg,state,24]>0){background_color<-"white"}
                    if (fit_stuff[sea,reg,state,24]<0){background_color<-"green"}
                    if (fit_stuff[sea,reg,state,21]<0.99){newline<-paste(newline," &{\\cellcolor{",background_color,"!25}{ }}",sep="")}
                    if (fit_stuff[sea,reg,state,21]>=0.99){newline<-paste(newline," &{\\cellcolor{",background_color,"!75}{",round(fit_stuff[sea,reg,state,i],02),"}}",sep="")}
                }
                for (i in 1){
                    if (fit_stuff[sea,reg,state,12]>fit_stuff[sea,reg,state,14]){background_color<-"red"}
                    if (fit_stuff[sea,reg,state,12]<fit_stuff[sea,reg,state,14]){background_color<-"blue"}
                    if (fit_stuff[sea,reg,state,21]<0.99){newline<-paste(newline," &{\\cellcolor{",background_color,"!25}{ }}",sep="")}
                    if (fit_stuff[sea,reg,state,21]>=0.99){newline<-paste(newline," &{\\cellcolor{",background_color,"!75}{",round(fit_stuff[sea,reg,state,12]-fit_stuff[sea,reg,state,14],02),"}}",sep="")}
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
    lines[index<-index+1]="\\end{document}"
    writeLines(lines, table)
    close(table)

    table<-file(paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",region_name,"_",period,"_fit_",fit_style,"_reduced.tex",sep=""))
    options(scipen=100) ;  lines=c() ; index=0
    lines[index<-index+1]="\\documentclass[a4paper,12pt]{article}"
    lines[index<-index+1]="\\usepackage{xcolor,colortbl,pgf}"
    lines[index<-index+1]="\\usepackage{makecell}"
    lines[index<-index+1]="\\usepackage[top=3cm, bottom=3cm, left=3cm, right=3cm]{geometry}"

    lines[index<-index+1]="\\definecolor{white}{rgb}{1,1,1}"
    lines[index<-index+1]="\\definecolor{green}{rgb}{0.5,1,0.5}"
    lines[index<-index+1]="\\definecolor{blue}{rgb}{0.3,1,1}"
    lines[index<-index+1]="\\definecolor{red}{rgb}{1,0.3,1}"

    lines[index<-index+1]="\\begin{document}"
    lines[index<-index+1]=paste("\\begin{figure}[!h]")
    lines[index<-index+1]=paste("\\begin{picture}(500,630)")
    lines[index<-index+1]=paste("\\put(0,0){\\rotatebox{90}{")

    lines[index<-index+1]="\\setlength{\\tabcolsep}{4pt}"
    lines[index<-index+1]=paste("\\begin{tabular}{c||cccccc||cccccc||cccccc||cccccc}")
    lines[index<-index+1]=paste("& \\multicolumn{6}{c}{MAM} & \\multicolumn{6}{c}{JJA} & \\multicolumn{6}{c}{SON} & \\multicolumn{6}{c}{DJF}","\\\\")
    lines[index<-index+1]=paste("& \\multicolumn{3}{c}{cold} & \\multicolumn{3}{c}{warm}& \\multicolumn{3}{c}{cold} & \\multicolumn{3}{c}{warm} & \\multicolumn{3}{c}{cold} & \\multicolumn{3}{c}{warm}& \\multicolumn{3}{c}{cold} & \\multicolumn{3}{c}{warm}","\\\\")
    lines[index<-index+1]="reg&b1&b2&tr&b1&b2&tr&b1&b2&tr&b1&b2&tr&b1&b2&tr&b1&b2&tr&b1&b2&tr&b1&b2&tr\\\\" 
    lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"
    lines[index<-index+1]="\\Xhline{2\\arrayrulewidth}"

    for (reg in ID_select){
        newline<-paste(reg)
        for (sea in 1:4){
            for (state in 1:2){
                for (i in c(12,14,15)){
                    if (fit_stuff[sea,reg,state,24]>0){background_intensity<-"!25"}
                    if (fit_stuff[sea,reg,state,24]<0){background_intensity<-"!75"}
                    if (fit_stuff[sea,reg,state,11]>fit_stuff[sea,reg,state,13]){background_color<-"green"}
                    if (fit_stuff[sea,reg,state,11]<fit_stuff[sea,reg,state,13]){background_color<-"blue"}
                    if (fit_stuff[sea,reg,state,21]<0.99){newline<-paste(newline," &{\\cellcolor{white!10}{ }}",sep="")}
                    if (fit_stuff[sea,reg,state,21]>=0.99 & i!=15){newline<-paste(newline," &{\\cellcolor{",background_color,background_intensity,"}{\\footnotesize{",round(fit_stuff[sea,reg,state,i],02),"}}}",sep="")}
                    if (fit_stuff[sea,reg,state,21]>=0.99 & i==15){newline<-paste(newline," &{\\cellcolor{",background_color,background_intensity,"}{\\footnotesize{",round(fit_stuff[sea,reg,state,i]),"}}}",sep="")}
                }
                
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
    lines[index<-index+1]=paste("}}")
    lines[index<-index+1]=paste("\\end{picture}")
    lines[index<-index+1]=paste("\\end{figure}")
      
    lines[index<-index+1]="\\end{document}"
    writeLines(lines, table)
    close(table)

}


write_multi_plot_tex <- function(filename,ID_select=c(20,14)){
    scale<-0.7
    types<-c(6,4,5)

    table<-file(paste(filename,"multi.tex",sep=""))
    options(scipen=100) ; lines=c() ; index=0
    lines[index<-index+1]="\\documentclass[a4paper,12pt]{article}"
    lines[index<-index+1]="\\usepackage{xcolor,colortbl,pgf}"

    lines[index<-index+1]="\\begin{document}"
    lines[index<-index+1]=paste("\\newcommand{\\source}{",filename,".pdf","}",sep="")

    for (sea in 1:5){
        lines[index<-index+1]=paste("\\begin{figure}[!h]")
        lines[index<-index+1]=paste("\\begin{picture}(500,500)")
        for (i in 1:length(types)){
            type<-types[i]
            count<-0
            for (q in ID_select){
                count<-count+1
                for (state in 1:2){
                    lines[index<-index+1]=paste("\\put(",(count-1)*200+(state-1)*100,",",(length(types)-i)*100,"){\\includegraphics[page=",8+(sea-1)*288+12*(q-1)+6*(state-1)+type,",scale=",scale,"]{\\source}}",sep="")

                    # axes
                    if (count==1 & state==1 & type==4){lines[index<-index+1]=paste("\\put(",(count-1)*200+(state-1)*100,",",(length(types)-i)*100,"){\\includegraphics[page=",3,",scale=",scale,"]{\\source}}",sep="")}
                    if (count==1 & state==1 & type==5){lines[index<-index+1]=paste("\\put(",(count-1)*200+(state-1)*100,",",(length(types)-i)*100,"){\\includegraphics[page=",5,",scale=",scale,"]{\\source}}",sep="")}
                    if (count==1 & state==1 & type==6){lines[index<-index+1]=paste("\\put(",(count-1)*200+(state-1)*100,",",(length(types)-i)*100,"){\\includegraphics[page=",7,",scale=",scale,"]{\\source}}",sep="")}
                    if (i==length(types)){lines[index<-index+1]=paste("\\put(",(count-1)*200+(state-1)*100,",",(length(types)-i)*100,"){\\includegraphics[page=",1,",scale=",scale,"]{\\source}}",sep="")}
                }
            }
        }
        lines[index<-index+1]=paste("\\end{picture}")
        lines[index<-index+1]=paste("\\end{figure}")        
    }
    lines[index<-index+1]="\\end{document}"



    writeLines(lines, table)
    close(table)
}

#write_multi_plot_tex("/home/peter/Dokumente/pik/backuped/plots/_TMean/91_7/regional/ward24/1950-2014/_ward24_dist_diff_fit_plot__TMean_1950-2014_2expo_4:100")