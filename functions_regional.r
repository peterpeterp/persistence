
wave_region <- function(wavenumber,trough,lats){
    # creates a table with coordinates of rectangles representing wave ridges and troughs
    wellen=array(NA,dim=c(40,13))
    brei=360/wavenumber/4
    lon=trough
    lon0=trough
    i=1
    wellen[,7:8]=lats[1]
    wellen[,9:10]=lats[2]
    while ((round(x=lon,digits=3)==round(x=lon0,digits=3) & i>1)==FALSE){
        if (lon>(180-brei) & lon<(180+brei)){
            wellen[i,1:4]=c(180,lon-brei,lon-brei,180)
            wellen[i,13]=i  
            i=i+1   
            wellen[i,1:4]=c(lon+brei-360,-180,-180,lon+brei-360)
            wellen[i,13]=i  
            i=i+1   
            lon=lon+2*brei
        }
        if (lon<=(180-brei)){
            wellen[i,1:4]=c(lon+brei,lon-brei,lon-brei,lon+brei)
            wellen[i,13]=i  
            i=i+1   
            lon=lon+2*brei
        }
        if (lon>=(180+brei)){
            lon=lon-360
            wellen[i,1:4]=c(lon+brei,lon-brei,lon-brei,lon+brei)
            wellen[i,13]=i  
            i=i+1   
            lon=lon+2*brei
        }
        
    }
    write.table(wellen[1:length(which(!is.na(wellen[,13]))),],paste("../data/region_poligons/",wavenumber,"wave.txt",sep=""))
}

points_to_regions <- function(dat,region_names=c("mid_lat_belt","srex","7rect","6wave","7wave","8wave")){
    # loads region coordinates and writes a file in which grid points are associated to regions
    # outpufile has following columns: ID, regions from region_names
    library(SDMTools)
    ntot=length(dat$ID)
    points=cbind(x=dat$lon,y=dat$lat)

    region = array(NA,dim=c(ntot,(length(region_names)+1)))

    for (j in 1:length(region_names)){
        poli=read.table(paste("../data/region_poligons/",region_names[j],".txt",sep=""))
        for (k in 1:dim(poli)[1]){
            x=c()
            y=c()
            for (i in 1:6){
                if (!is.na(poli[k,i])){
                    x[i]=poli[k,i]
                    y[i]=poli[k,(i+6)]
                }
            }
            poligon=cbind(x=x,y=y)
            inside=pnt.in.poly(points,poligon)$pip
            region[which(inside==1),(j+1)]=poli[k,13]
        }
    }

    if (length(region_names)>1){
        region[1:ntot,1]=dat$ID
        write.table(region,"../data/ID-regions.txt")
    }
    return(region[1:ntot,2])
}


regions_color <- function(reihen,reihen_sig,worldmap,titles,poli,filename_plot){
    # plots worldmap and colored regions on it
    jet.colors <- colorRampPalette( c(rgb(0.2,0.6,0.2),rgb(0.5,1,0.5), rgb(0.98,0.98,0.98) ,rgb(1,0.5,1),rgb(0.6,0.2,0.6)))
    
    nbcol <- 101
    color <- jet.colors(nbcol)

    pdf(file = filename_plot,width=12,height=8)

    for (rei in 1:dim(reihen)[1]){            
        y=c()
        index=c()
        signi=c()
        j=0
        for (i in 1:dim(poli)[1]){
            poliLabel=i
            if (!is.na(reihen[rei,poliLabel])){
                j=j+1
                y[j]=reihen[rei,poliLabel]
                index[j]=i 
                if (abs(reihen[rei,poliLabel])>0.0001){
                    signi[j]=sprintf("%.04f",reihen_sig[rei,poliLabel])
                }         
            }
        }
        aushol=max(c(abs(max(y)),abs(min(y))))
        y[j+1]=-aushol
        y[j+2]=aushol
        facetcol <- cut(y,nbcol)  

        print(titles[rei])
        plot(worldmap,main=titles[rei])

        for (i in 1:j){
            lon=poli[index[i],1:6]
            lat=poli[index[i],7:12]
            lon=lon[!is.na(lon)]
            lat=lat[!is.na(lat)]
             polygon(x=lon,y=lat,col=color[facetcol[i]],border="green")
            text(mean(lon),mean(lat),label=signi[i],cex=0.7,col="black")
        }
        image.plot(legend.only=T, zlim=range(y), col=color)
    }
    graphics.off()
    return()
}

#============================================================================================================================


duration_region <- function(regions,reg,dur,dur_mid){
    # combines all recorded durations of one region to one duration array, same for dur_mid
    inside=which(regions==reg)
    duration=array(NA,dim=c(1000000))
    duration_mid=array(NA,dim=c(1000000))
    count=1
    # combines the recorded periods from all the grid points of one region in one array
    for (i in inside){
        values=length(which(!is.na(dur[i,])))
        duration[count:(count+values)]=dur[i,1:values]
        duration_mid[count:(count+values)]=dur_mid[i,1:values]
        count=count+values
    }
    nona=which(!is.na(duration) & !is.na(duration_mid))
    duration=duration[nona]
    duration_mid=duration_mid[nona]
    return(list(duration=duration,duration_mid=duration_mid))
}



regional_distributions <- function(dat,region_name,trendID,additional_style,dataset){
    # performs the entire regional analysis of markov and duration
    # result will be written in nc file

    # pnly for one trend and 2 states until now
    ntot=length(dat$ID)

    poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))
    regNumb=dim(poli)[1]
    print(regNumb)
    IDregions=points_to_regions(dat,c(region_name))

    season_names=c("MAM","JJA","SON","DJF","4seasons")


    for (sea in 1:5){
        season=season_names[sea]
        print(paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_",season,".nc",sep=""))
        nc_dur=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",trendID,dataset,"_duration_",season,".nc",sep=""))
        dur=var.get.nc(nc_dur,"dur")
        dur_mid=var.get.nc(nc_dur,"dur_mid")   

        reg_dur=array(NA,dim=c(regNumb,2,365*65*100))
        reg_dur_mid=array(NA,dim=c(regNumb,2,365*65*100))
        maxis=array(NA,dim=c(2*regNumb))

        for (state in 1:2){
            for (reg in 1:regNumb){            
                tmp=duration_region(IDregions,reg,dur[1:ntot,state,],dur_mid[1:ntot,state,])
                duration=tmp$duration
                duration_mid=tmp$duration_mid
                if (length(duration)>100){
                    ord=order(duration_mid)
                    maxis[(state-1)*regNumb+reg]=length(duration)
                    reg_dur[reg,state,1:maxis[(state-1)*regNumb+reg]]=duration[ord]
                    reg_dur_mid[reg,state,1:maxis[(state-1)*regNumb+reg]]=duration_mid[ord]
                }
            }
        }
        len=max(maxis,na.rm=TRUE)
        print(paste("../data/",trendID,"/",dataset,additional_style,"/","regional/",trendID,dataset,"_",region_name,"_duration_",season,".nc",sep=""))
        duration_write(filename=paste("../data/",trendID,"/",dataset,additional_style,"/","regional/",trendID,dataset,"_",region_name,"_duration_",season,".nc",sep=""),dur=reg_dur[,,1:len],dur_mid=reg_dur_mid[,,1:len],len=len,ID_length=regNumb,ID_name=region_name)

    }

}


# ------------------------------------------------------------------------------------------------
plot_boxplot <- function(quans,x,width,color="white",border="black",density=NA){

    links=x-width/2
    rechts=x+width/2
    polygon(x=c(rechts,links,links,rechts),y=c(quans[2],quans[2],quans[4],quans[4]),col=color,border=border,density=density)    
    #polygon(x=c(rechts,links,links,rechts),y=c(quans[2],quans[2],quans[4],quans[4]),col=color,border=border,density=density)    
    for (qu in quans[1:5]){
        lines(c(links,rechts),c(qu,qu),col=border)
    }
    lines(c(links,links),quans[c(2,4)],col=border)
    lines(c(rechts,rechts),quans[c(2,4)],col=border)
    lines(c(x,x),quans[c(1,2)],lty=2,col=border)
    lines(c(x,x),quans[c(4,5)],lty=2,col=border)
}


plot_regional_boxplots <- function(period,region_name,trendID,additional_style,dataset,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){
    # performs the entire regional analysis of markov and duration
    # result will be written in nc file

    nc_oth = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_others.nc",sep=""))
    others=var.get.nc(nc_oth,"other_stuff")
    nc_qua = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_quantiles.nc",sep=""))
    quantiles=var.get.nc(nc_qua,"quantile_stuff")

    regNumb=dim(quantiles)[2]
    seaNumb=6
    
    season_names=c("MAM","JJA","SON","DJF","4seasons")

    color=c(rgb(0.5,0.5,1,0.8),rgb(1,0.5,0.5,0.8))
    pos=c(-1,1)*0.15
    width=6
    height=3

    # regional focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",period,"/",region_name,"_",trendID,"_",period,"_boxplots_regional_new.pdf",sep=""),width=width,height=height)
    par(mfrow=c(1,1))
    par(mar=c(1,4,2,3))   
    at_=seq(1, regNumb, 1)
    at_=c(at_-0.15,at_+0.15)
    buchstaben=c("a","b","c","d")
    drueber=2
    for (sea in 1:length(season_names)){   
        season=season_names[sea]

        plot(NA,xlim=c(0,7.5),ylim=c(0,25),frame.plot=FALSE,axes=FALSE,ylab="days")#(max(quantiles[sea,,,1:5])+drueber)
        axis(2)
        for (i in axis(2)){
            abline(h=i,col=rgb(0.5,0.5,0.5,0.5),lty=3)
        }
        legend("topright",legend=c(buchstaben[sea]),bty="n")
        #text(8,23,label=buchstaben[sea],cex=2)#(max(quantiles[sea,,,1:5],na.rm=TRUE)+drueber-0.5)

        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,,1],reg+pos[state],0.3,color[state])
                text(reg,24,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))#max(quantiles[sea,,,1:5])+drueber
            }
        }        
        for (oth in c(1)){
            points(at_[1:regNumb],others[sea,,1,oth],col="blue",pch=oth)
            points(at_[(regNumb+1):(regNumb*2)],others[sea,,2,oth],col="red",pch=oth)
        }
    }
    graphics.off()
}

plot_regional_boxplots_vergleich <- function(period1,period2,region_name,trendID,additional_style,dataset,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){
    # performs the entire regional analysis of markov and duration
    # result will be written in nc file


    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period1,"/",trendID,"_",dataset,"_",region_name,"_",period1,"_quantiles.nc",sep=""))
    quantiles1=var.get.nc(nc,"quantile_stuff")
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period2,"/",trendID,"_",dataset,"_",region_name,"_",period2,"_quantiles.nc",sep=""))
    quantiles2=var.get.nc(nc,"quantile_stuff")
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period1,"/",trendID,"_",dataset,"_",region_name,"_",period1,"_others.nc",sep=""))
    others1=var.get.nc(nc,"other_stuff")
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period2,"/",trendID,"_",dataset,"_",region_name,"_",period2,"_others.nc",sep=""))
    others2=var.get.nc(nc,"other_stuff")

    regions=var.get.nc(nc,"region")
    regNumb=dim(quantiles1)[2]
    seaNumb=5
    
    season_names=c("MAM","JJA","SON","DJF","4seasons")

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/","/",region_name,"_",period1,"_diff_",period2,"_boxplots_regional.pdf",sep=""))
    par(mfrow=c(1,1))
    par(mar=c(1,5,4,3))   
    at_=seq(1, regNumb, 1)
    at_=c(at_-0.15,at_+0.15)
    color=c()
    maxi=c()
    taus=c(0.05,0.25,0.5,0.75,0.95,0.91,0.98)
    color=c(rgb(0.5,0.5,1,0.8),rgb(1,0.5,0.5,0.8))
    pos=c(-1,1)*0.15

    quantiles=quantiles2-quantiles1
    others=others2-others1

    # regional focus
    for (sea in 1:length(season_names)){   
        season=season_names[sea]

        plot(NA,xlim=c(0,8),ylim=c(min(quantiles[sea,,,1:5,1]),(max(quantiles[sea,,,1:5,1])+0.5)),frame.plot=FALSE,axes=FALSE,main=season,ylab="days")
        axis(2)
        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,,1],reg+pos[state],0.3,color[state])
                text(reg,max(quantiles[sea,,,1:5,1])+0.5,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))
            }
        }
        for (oth in c(1)){
            points(at_[1:regNumb],others[sea,,1,oth],col="blue",pch=oth)
            points(at_[(regNumb+1):(regNumb*2)],others[sea,,2,oth],col="red",pch=oth)
        }
        for (quA in c(4,5)){
            points(at_[1:regNumb],quantiles[sea,,1,quA,1],col="black",pch=2)
            points(at_[(regNumb+1):(regNumb*2)],quantiles[sea,,2,quA,1],col="black",pch=2)
        }
    }
    graphics.off()

    # seasonal focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/","/",region_name,"_",period1,"_diff_",period2,"_boxplots_seasonal.pdf",sep=""))
    par(mfrow=c(1,1))
    par(mar=c(1,5,4,3))   
    at_=seq(1, seaNumb, 1)
    at_=c(at_-0.15,at_+0.15)

    for (reg in 1:regNumb){
        plot(NA,xlim=c(0,8),ylim=c(min(quantiles[,reg,,1:5,1]),(max(quantiles[,reg,,1:5,1])+0.5)),frame.plot=FALSE,axes=FALSE,main=region_names[reg],ylab="days")
        axis(2)
        for (sea in 1:seaNumb){   
            season=season_names[sea]     
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,,1],sea+pos[state],0.3,color[state])
                text(sea,max(quantiles[,reg,,1:5,1])+0.5,season_names[sea],col=rgb(0.5,0.5,0.5,0.5))
            }
        }
        for (oth in c(1)){
            points(at_[1:seaNumb],others[,reg,1,oth],col="blue",pch=1)
            points(at_[(seaNumb+1):(seaNumb*2)],others[,reg,2,oth],col="red",pch=1)
        }
        for (quA in c(4,5)){
            points(at_[1:seaNumb],quantiles[,reg,1,quA,1],col="black",pch=2)
            points(at_[(seaNumb+1):(seaNumb*2)],quantiles[,reg,2,quA,1],col="black",pch=2)
        }
    }
    graphics.off()
}

plot_regional_fit_vergleich <- function(period1,period2,region_name,trendID,additional_style,dataset,fit_style,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){
    # performs the entire regional analysis of markov and duration
    # result will be written in nc file


    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period1,"/",trendID,"_",dataset,"_",region_name,"_",period1,fit_style,".nc",sep=""))
    fit_stuff1=var.get.nc(nc,"fit_stuff")
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period2,"/",trendID,"_",dataset,"_",region_name,"_",period2,fit_style,".nc",sep=""))
    fit_stuff2=var.get.nc(nc,"fit_stuff")


    regions=var.get.nc(nc,"region")
    regNumb=dim(fit_stuff2)[2]
    seaNumb=5
    
    season_names=c("MAM","JJA","SON","DJF","4seasons")

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/","/",region_name,"_",period1,"_diff_",period2,fit_style,".pdf",sep=""))
    fit_stuff=fit_stuff1-fit_stuff2
    fit_stuff[,,,2]=1/fit_stuff1[,,,2]-1/fit_stuff2[,,,2]
    fit_stuff[,,,4]=1/fit_stuff1[,,,4]-1/fit_stuff2[,,,4]
    for (sea in 1:seaNumb){
        for (state in 1:2){
            plot(NA,xlim=c(0.5,regNumb+0.5),ylim=c(-2,2),xlab="",ylab="lifetime [days]",axes=FALSE,frame.plot=TRUE,cex=0.5)
            axis(2)
            axis(1, at=1:regNumb, labels=region_names) 
            abline(h=0,col="grey")

            points(fit_stuff[sea,,state,2],col="red")
            lines(fit_stuff[sea,,state,2],col="red")
            points(fit_stuff[sea,,state,4],col="orange")
            lines(fit_stuff[sea,,state,4],col="orange")
            par(new=TRUE)
            plot(NA,xlim=c(0.5,regNumb+0.5),ylim=c(-8,8),axes=FALSE,frame.plot=TRUE,cex=0.5,ylab="",xlab="")
            axis(4)
            points(fit_stuff[sea,,state,5],col="blue")
            lines(fit_stuff[sea,,state,5],col="blue")
        }
    }



    graphics.off()

}



plot_regional_fit_parameters <- function(period,trendID,additional_style,dataset,fit_style,region_name="7rect",region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){
    # performs the entire regional analysis of markov and duration
    # result will be written in nc file

    print(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,fit_style,".nc",sep=""))
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,fit_style,".nc",sep=""))
    fit_stuff=var.get.nc(nc,"fit_stuff")


    regNumb=dim(fit_stuff)[2]
    seaNumb=5
    season_names=c("MAM","JJA","SON","DJF","4seasons")

    color=c(rgb(0.5,0.5,1,0.8),rgb(1,0.5,0.5,0.8))
    color=c("blue","red")
    pos=c(-1,1)*0.15
    width=7
    height=5

    # regional focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",period,"/",trendID,"_",period,fit_style,".pdf",sep=""),width=width,height=height)
    par(mfrow=c(1,1))
    par(mar=c(3,5,0,1))   



    for (sea in 1:seaNumb){
        for (state in 1:2){
            plot(NA,xlim=c(0.5,regNumb+0.5),ylim=c(2,15),frame.plot=TRUE,axes=FALSE,ylab="lifetime [days]",xlab="")
            axis(1, at=1:7, labels=region_names) 
            axis(2)
            lines(1:7,1/fit_stuff[sea,,state,2],col=color[1])
            lines(1:7,1/fit_stuff[sea,,state,4],col=color[2])
            lines(1:7,fit_stuff[sea,,state,5],col="green")
        }
    }

    graphics.off()
}

write_regional_fit_table <- function(trendID="91_5",region_name="srex",period,fit_style1,fit_style2,region_names,ID_select){
    regNumb=length(ID_select)

    print(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style1,".nc",sep=""))
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style1,".nc",sep=""))
    fit_stuff1=var.get.nc(nc,"fit_stuff")
    print(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style2,".nc",sep=""))
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style2,".nc",sep=""))
    fit_stuff2=var.get.nc(nc,"fit_stuff")


    color=c("lightblue","lightred")

    table<-file(paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",period,"/",trendID,"_",region_name,"_",period,"_seasons_latex.tex",sep=""))
    options(scipen=100)

    lines=c()
    index=0

    lines[index+1]="\\documentclass[a3paper,12pt,landscape]{article}"
    lines[index+2]="\\usepackage{xcolor,colortbl,pgf}"
    lines[index+3]="\\definecolor{lightblue}{rgb}{0.75,0.75,1}"
    lines[index+4]="\\definecolor{lightred}{rgb}{1,0.75,0.75}"
    lines[index+5]="\\begin{document}"
    index=index+5

    for (sea in 1:5){
        for (state in 1:2){

            lines[index+1]=paste("\\begin{table}[!h]")
            lines[index+2]=paste("\\begin{tabular}{c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|cc|c|c|c|c|c|c|c|c|c|c|c|c|c}")
            lines[index+3]=paste("\\ reg & R2 & BIC & b & R2 & BIC & b1 & b2 & thresh & R2 & BIC & a1 & b1 & a2 & {b2}","\\","\\",sep="")
            index=index+3
            for (reg in ID_select){
                index=index+1
                BICs=c(fit_stuff1[sea,reg,state,c(16,20)],fit_stuff2[sea,reg,state,20])
                if (length(which(!is.na(BICs)))>0){best=which(BICs==min(BICs,na.rm=TRUE))}
                else{best=0}
                newline=paste("\\ ",region_names[reg],sep="")
                for (i in c(15,16,2)){
                    if (best==1){newline=paste(newline," &{\\cellcolor{",color[state],"}}",round(fit_stuff1[sea,reg,state,i],02),"",sep="")}
                    else{newline=paste(newline," &{",round(fit_stuff1[sea,reg,state,i],02),"}",sep="")}                }
                for (i in c(19,20,6,8,9)){
                    if (best==2){newline=paste(newline," &{\\cellcolor{",color[state],"}}",round(fit_stuff1[sea,reg,state,i],02),"",sep="")}
                    else{newline=paste(newline," &{",round(fit_stuff1[sea,reg,state,i],02),"}",sep="")}                }
                for (i in c(19,20,5,6,7,8)){
                    if (best==3){newline=paste(newline," &{\\cellcolor{",color[state],"}}",round(fit_stuff2[sea,reg,state,i],02),"",sep="")}
                    else{newline=paste(newline," &{",round(fit_stuff2[sea,reg,state,i],02),"}",sep="")}
                }   
                newline=paste(newline,paste("\\","\\",sep=""))
                lines[index]=newline
            }
            lines[index+1]=paste("\\end{tabular}")
            lines[index+2]=paste("\\end{table}")
            index=index+2

        }
        
    }
    lines[index+1]="\\end{document}"
    writeLines(lines, table)
    asdas



    y=fit_stuff[sea,,,18]-fit_stuff[sea,,,16]
    y2=fit_stuff[sea,,,17]-fit_stuff[sea,,,16]
    y[is.na(fit_stuff[sea,,,18])]=y2[is.na(fit_stuff[sea,,,18])]    

    for (state in 1:2){
        if (state==1){newline=paste("\\","multirow{2}{*}{$\\Delta$BIC}",sep="")}
        else{newline=""}
        for (reg in 1:regNumb){
            newline=paste(newline," &{\\cellcolor{",color[state],"}}",sep="")
            if (y[reg,state]<0){newline=paste(newline,"\\textit{\\textbf{",round(y[reg,state]),"}}")}
            else{newline=paste(newline,round(y[reg,state]))}
        }
        newline=paste(newline,paste("\\","\\",sep=""))
        if (state==2){newline=paste(newline,"\n\\hline")}
        lines[(state+2)]=newline
    } 

    y=1/fit_stuff[sea,,,8]
    y2=1/fit_stuff[sea,,,2]
    y[is.na(fit_stuff[sea,,,8])]=y2[is.na(fit_stuff[sea,,,8])]    

    for (state in 1:2){
        if (state==1){newline=paste("\\","multirow{2}{*}{1/b}",sep="")}
        else{newline=""}
        for (reg in 1:regNumb){
            newline=paste(newline," &{\\cellcolor{",color[state],"}}",sep="")
            newline=paste(newline,round(y[reg,state],02))
        }
        newline=paste(newline,paste("\\","\\",sep=""))
        if (state==2){newline=paste(newline,"\n\\hline")}
        lines[(state+4)]=newline
    } 

    y=fit_stuff[sea,,,10]

    for (state in 1:2){
        if (state==1){newline=paste("\\","multirow{2}{*}{A}",sep="")}
        else{newline=""}
        for (reg in 1:regNumb){
            newline=paste(newline," &{\\cellcolor{",color[state],"}}",sep="")
            newline=paste(newline,round(y[reg,state],03))
        }
        newline=paste(newline,paste("\\","\\",sep=""))
        if (state==2){newline=paste(newline,"\n\\hline")}
        lines[(state+6)]=newline
    } 
    y=fit_stuff[sea,,,13]

    for (state in 1:2){
        if (state==1){newline=paste("\\","multirow{2}{*}{R2(exp)}",sep="")}
        else{newline=""}
        for (reg in 1:regNumb){
            newline=paste(newline," &{\\cellcolor{",color[state],"}}",sep="")
            newline=paste(newline,round(y[reg,state],03))
        }
        newline=paste(newline,paste("\\","\\",sep=""))
        if (state==2){newline=paste(newline,"\n\\hline")}
        lines[(state+8)]=newline
    } 
    
    writeLines(lines, table)

    close(table)
}

