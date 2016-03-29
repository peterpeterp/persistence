


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


write_tex_table_vert<-function(filename,outs,values,signis){

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
    lines[index<-index+1]="\\definecolor{turkis}{rgb}{0.5,1,1}"
    lines[index<-index+1]="\\definecolor{violet}{rgb}{1,0.5,1}"

    lines[index<-index+1]="\\begin{document}"
    lines[index<-index+1]=paste("\\begin{table}[!h]")
    lines[index<-index+1]=paste("\\begin{tabular}{c||cccc||cccc||cccc||cccc}")


    
    for (sea in 1:4){
        for (state in 1:2){
            for (i in outs){
            	newline<-paste(state)
            	for (reg in ID_select){
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
    lines[index<-index+1]=paste("\\vspace{0cm}")
      
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