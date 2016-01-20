
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



regional_attribution <- function(dat,region_name,trendID,additional_style="",dataset="_TMean",IDregions=c("from polygons"),regNumb=7,comment="polygons"){
    # performs the entire regional analysis of markov and duration
    # result will be written in nc file

    # pnly for one trend and 2 states until now
    ntot=length(dat$ID)


    if (IDregions[1]=="from polygons"){
        poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))
        regNumb=dim(poli)[1]
        IDregions_tmp=points_to_regions(dat,c(region_name))
        IDregions=array(NA,c(ntot,5))
        for (i in 1:5){
            IDregions[,i]=IDregions_tmp
        }
    }

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
                tmp=duration_region(IDregions[,sea],reg,dur[1:ntot,state,],dur_mid[1:ntot,state,])
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
        duration_write(filename=paste("../data/",trendID,"/",dataset,additional_style,"/","regional/",trendID,dataset,"_",region_name,"_duration_",season,".nc",sep=""),dur=reg_dur[,,1:len],dur_mid=reg_dur_mid[,,1:len],len=len,ID_length=regNumb,ID_name=region_name,comment=comment)

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

    season_names=c("MAM","JJA","SON","DJF","4seasons")
    state_names=c("cold","warm")

    table<-file(paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",period,"/",trendID,"_",region_name,"_",period,"_all_latex.tex",sep=""))
    options(scipen=100)

    colors=c("white","groegree","zehngree","funfziggree","hundertgree","turkis","violet")
    lines=c()
    index=0

    lines[index+1]="\\documentclass[a4paper,12pt]{article}"
    lines[index+2]="\\usepackage{xcolor,colortbl,pgf}"
    lines[index+3]="\\definecolor{white}{rgb}{1,1,1}"
    lines[index+4]="\\definecolor{groegree}{rgb}{0.85,1,0.85}"
    lines[index+5]="\\definecolor{zehngree}{rgb}{0.75,1,0.75}"
    lines[index+6]="\\definecolor{funfziggree}{rgb}{0.5,1,0.5}"
    lines[index+7]="\\definecolor{hundertgree}{rgb}{0.3,1,0.3}"
    lines[index+8]="\\definecolor{turkis}{rgb}{0.5,1,1}"
    lines[index+9]="\\definecolor{violet}{rgb}{1,0.5,1}"
    lines[index+10]="\\begin{document}\\vspace{-2cm}"
    lines[index+11]="\\fcolorbox{groegree}{groegree}{dBIC$<$0}\\"
    lines[index+12]="\\fcolorbox{zehngree}{zehngree}{dBIC$<$-10}\\"
    lines[index+13]="\\fcolorbox{funfziggree}{funfziggree}{dBIC$<$-50}\\"
    lines[index+14]="\\fcolorbox{hundertgree}{hundertgree}{dBIC$<$-100}\\"
    lines[index+15]="\\fcolorbox{turkis}{turkis}{P2$>$P1}\\"
    lines[index+16]="\\fcolorbox{violet}{violet}{P2$<$P1}\\"
    index=index+16

    for (sea in 1:5){
        for (state in 1:2){

            lines[index+1]=paste("\\begin{table}[!h]")
            lines[index+2]=paste("\\vspace{0cm}")
            lines[index+3]=paste("\\hspace{0cm}")
            lines[index+4]=paste("\\begin{tabular}{c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|cc|c|c|c|c|c|c|c|c|c|c|c|c|c}")

            lines[index+5]=paste("\\multicolumn{15}{1}{",season_names[sea]," ",state_names[state]," ",period,"}\\","\\",sep="")
            lines[index+6]=paste("\\ reg & R2 & BIC & P & R2 & BIC & thresh & P1 & P2 & P2-P1 & R2 & BIC & a1 & a2 & P1 & P2","\\","\\",sep="")
            #lines[index+6]=paste("\\ reg & R2 & BIC & P & R2 & BIC & thresh & P1 & P2 & P2-P1","\\","\\",sep="")
            index=index+6
            for (reg in ID_select){
                background=c(colors[1],colors[1],colors[1])
                BICs=c(fit_stuff1[sea,reg,state,c(16,20)],fit_stuff2[sea,reg,state,20])
                # whithout overlap
                #BICs=c(fit_stuff1[sea,reg,state,c(16,20)])
                if (length(which(!is.na(BICs)))>0){
                    worst=BICs[which(BICs==max(BICs,na.rm=TRUE))]
                    for (i in 1:3){
                        if (!is.na(BICs[i])){
                            if (BICs[i]<worst){background[i]=colors[2]}
                            if ((BICs[i]+10)<worst){background[i]=colors[3]}
                            if ((BICs[i]+50)<worst){background[i]=colors[4]}
                            if ((BICs[i]+100)<worst){background[i]=colors[5]}
                        }
                    }
                }
                newline=paste("\\ ",region_names[reg],sep="")
                for (i in c(15,16)){
                    newline=paste(newline," &{\\cellcolor{",background[1],"}}",round(fit_stuff1[sea,reg,state,i],02),"",sep="")}
                newline=paste(newline," &{\\cellcolor{",background[1],"}}",round(exp(-fit_stuff1[sea,reg,state,2])*100,01),sep="")
                for (i in c(19,20,9)){
                    newline=paste(newline," &{\\cellcolor{",background[2],"}}",round(fit_stuff1[sea,reg,state,i],02),"",sep="")}
                newline=paste(newline," &{\\cellcolor{",background[2],"}}",round(exp(-fit_stuff1[sea,reg,state,6])*100,01),sep="")
                newline=paste(newline," &{\\cellcolor{",background[2],"}}",round(exp(-fit_stuff1[sea,reg,state,8])*100,01),sep="")
                diffP=round(exp(-fit_stuff1[sea,reg,state,8])*100,03)-round(exp(-fit_stuff1[sea,reg,state,6])*100,03)
                if (diffP>0){farbe=colors[6]}
                else {farbe=colors[7]}
                newline=paste(newline," &{\\cellcolor{",farbe,"}}",round(diffP,01),sep="")

                for (i in c(19,20,5,7)){
                    newline=paste(newline," &{\\cellcolor{",background[3],"}}",round(fit_stuff2[sea,reg,state,i],02),"",sep="")}
                newline=paste(newline," &{\\cellcolor{",background[3],"}}",round(exp(-fit_stuff2[sea,reg,state,6])*100,01),sep="")
                newline=paste(newline," &{\\cellcolor{",background[3],"}}",round(exp(-fit_stuff2[sea,reg,state,8])*100,01),sep="")
                
                newline=paste(newline,paste("\\","\\",sep=""))
                lines[index+1]=newline
                index=index+1
            }
            lines[index+1]=paste("\\end{tabular}")
            #lines[index+2]=paste("\\end{table} \\newpage")
            lines[index+2]=paste("\\end{table}")
            lines[index+3]=paste("\\vspace{0cm}")
            index=index+3

        }
        
    }
    lines[index+1]="\\end{document}"
    writeLines(lines, table)
    close(table)
}

fit_info_to_map <- function(trendID="91_5",region_name="srex",period,fit_style1,fit_style2,region_names,ID_select){
    # plots worldmap and colored regions on it
    library(rworldmap)
    library(fields)
    worldmap = getMap(resolution = "low")

    poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))


    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style1,".nc",sep=""))
    fit_stuff1=var.get.nc(nc,"fit_stuff")

    season_names=c("MAM","JJA","SON","DJF","4seasons")
    state_names=c("cold","warm")

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",period,"/",trendID,"_",region_name,"_",period,"_seasons.pdf",sep=""))

    for (sea in 1:5){
        for (state in 1:2){
            plot(worldmap,main=paste(season_names[sea],state_names[state]))

            for (reg in ID_select){
                BICs=c(fit_stuff1[sea,reg,state,c(16,20)])
                if (is.na(BICs[1])){
                    fit=1
                    farbe="white"#rgb(1,0,1,0.3)
                }
                if (is.na(BICs[2])){
                    fit=2
                    farbe="white"#rgb(0,1,0,0.3)
                }
                if (is.na(BICs[2]) & is.na(BICs[1])){
                    fit=0
                    farbe="white"
                    text=""
                }
                if (length(which(!is.na(BICs)))==2){
                    opacity=abs((BICs[2]-BICs[1])/100)
                    if (opacity>1){opacity=1}
                    if (BICs[1]<BICs[2]){
                        fit=1
                        farbe=rgb(0.9,0.5,0.9,opacity)}
                    if (BICs[1]>BICs[2]){
                        fit=2
                        farbe=rgb(0.5,0.9,0.5,opacity)}
                }
                if (fit==1){
                    text=paste("P=",round(exp(-fit_stuff1[sea,reg,state,2])*100,01))
                }                
                if (fit==2){
                    text=paste("P1=",round(exp(-fit_stuff1[sea,reg,state,6])*100,01),"P2=",round(exp(-fit_stuff1[sea,reg,state,8])*100,01),"\n thresh=",round(fit_stuff1[sea,reg,state,9]))
                }

                lon=poli[reg,1:6]
                lat=poli[reg,7:12]
                lon=lon[!is.na(lon)]
                lat=lat[!is.na(lat)]
                #polygon(x=lon,y=lat,col=rgb(farb,1,farb),border="green")
                polygon(x=lon,y=lat,col="white",border="white")
                polygon(x=lon,y=lat,col=farbe,border="white")
                text(mean(lon),mean(lat),label=text,cex=0.3,col="black")

                
            }
            #image.plot(legend.only=T, zlim=range(y), col=color)
        }
    }

    graphics.off()
    return()
}


plot_fits_for_region <- function(reg,IDregions=c("from polygons"),period="1950-2014",trendID="91_5",dataset="_TMean",fit_style="2expo_thresh_5-15",region_name="srex"){

    if (IDregions[1]=="from polygons"){
        poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))
        regNumb=dim(poli)[1]
        IDregions_tmp=points_to_regions(dat,c(region_name))
        IDregions=array(NA,c(ntot,5))
        for (i in 1:5){
            IDregions[,i]=IDregions_tmp
        }
    }


    print(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style,".nc",sep=""))
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style,".nc",sep=""))
    fit_stuff_reg=var.get.nc(nc,"fit_stuff")
    print(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_fit_",fit_style,".nc",sep=""))
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_fit_",fit_style,".nc",sep=""))
    fit_stuff_individual=var.get.nc(nc,"fit_stuff")
    print(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_distributions.nc",sep=""))
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_distributions.nc",sep=""))
    distr_stuff_individual=var.get.nc(nc,"distr_stuff")

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",period,"/",trendID,"_",region_name,"_",period,"_region_",reg,".pdf",sep=""))
    
    color=c(rgb(0.5,1,0.8),rgb(0.3,0.9,0.8),rgb(0.9,0.6,0.8),rgb(0.8,0.2,0.6),rgb(0.5,0.5,0.5,0.2))

    X=(0:30)+0.5
    for (sea in 1:5){

        ID_select=which(IDregions[,sea]==reg)
        print(ID_select)

        for (state in 1:2){
            plot(NA,ylim=c(0,10),xlim=c(0,10),axes=FALSE,frame.plot=FALSE)
            text(5,5,paste("reg:",reg,"\nseason:",sea,"state:",state,"\npoints:",length(ID_select)))
            plot(NA,xlab="days",ylab="counts",ylim=c(0,150),xlim=c(0,30),axes=TRUE,frame.plot=TRUE)
            for (q in ID_select){
                points(distr_stuff_individual[sea,q,state,1,],distr_stuff_individual[sea,q,state,3,],pch=16,col=color[5],cex=1.5)
            }

            plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE)
            for (q in ID_select){points(distr_stuff_individual[sea,q,state,1,],distr_stuff_individual[sea,q,state,2,],pch=16,col=color[5],cex=1.5)}
            for (q in ID_select){
                abline(v=fit_stuff_individual[sea,q,state,9],col=color[5],lty=2,lwd=0.3)
                lines(X,combi_expo(X,fit_stuff_individual[sea,q,state,5],fit_stuff_individual[sea,q,state,6],fit_stuff_individual[sea,q,state,8],fit_stuff_individual[sea,q,state,9]),col=color[1],lwd=0.3)
                lines(X,expo(X,fit_stuff_individual[sea,q,state,1],fit_stuff_individual[sea,q,state,2]),col=color[3],lwd=0.3,lty=2)
            }
            abline(v=fit_stuff_reg[sea,reg,state,9],col=color[5],lty=2,lwd=2)
            lines(X,combi_expo(X,fit_stuff_reg[sea,reg,state,5],fit_stuff_reg[sea,reg,state,6],fit_stuff_reg[sea,reg,state,8],fit_stuff_reg[sea,reg,state,9]),col=color[2],lwd=2,lty=1)
            lines(X,expo(X,fit_stuff_reg[sea,reg,state,1],fit_stuff_reg[sea,reg,state,2]),col=color[4],lty=1,lwd=2)
            

            plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,log="y")
            for (q in ID_select){points(distr_stuff_individual[sea,q,state,1,],distr_stuff_individual[sea,q,state,2,],pch=16,col=color[5],cex=1.5)}
            for (q in ID_select){
                abline(v=fit_stuff_individual[sea,q,state,9],col=color[5],lty=2,lwd=0.3)
                lines(X,combi_expo(X,fit_stuff_individual[sea,q,state,5],fit_stuff_individual[sea,q,state,6],fit_stuff_individual[sea,q,state,8],fit_stuff_individual[sea,q,state,9]),col=color[1],lwd=0.3)
                lines(X,expo(X,fit_stuff_individual[sea,q,state,1],fit_stuff_individual[sea,q,state,2]),col=color[3],lwd=0.3,lty=2)
            }
            abline(v=fit_stuff_reg[sea,reg,state,9],col=color[5],lty=2,lwd=2)
            lines(X,combi_expo(X,fit_stuff_reg[sea,reg,state,5],fit_stuff_reg[sea,reg,state,6],fit_stuff_reg[sea,reg,state,8],fit_stuff_reg[sea,reg,state,9]),col=color[2],lwd=2,lty=1)
            lines(X,expo(X,fit_stuff_reg[sea,reg,state,1],fit_stuff_reg[sea,reg,state,2]),col=color[4],lty=1,lwd=2)

        }
    }
    
    
}


