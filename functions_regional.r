
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

points_to_regions <- function(dat,region_names=c("srex","7rect","6wave","7wave","8wave")){
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

markov_regional_trend <- function(dat,value,regions,yearPeriod,regNumb){
    # value is an array of dim=c(1319,years)
    # averages values in one region
    # performs linear regression and Mann-Kendall test for averaged values

    yearPeriod=yearPeriod-1949
    years=dim(value)[2]
    valOut=array(NA,dim=c(regNumb,years))
    for (reg in 1:regNumb){
        inside=which(regions==reg)
        for (j in 1:years){
            valOut[reg,j]=mean(value[inside,j],na.rm=TRUE)
        }
    }
    out=array(NA,c(2,regNumb))
    out_sig=array(NA,c(2,regNumb))

    x=seq(yearPeriod[1],yearPeriod[2],1)
    for (reg in 1:regNumb){
        y=valOut[reg,x]
        if (length(which(is.nan(y)==1))<10){
            tmp=MannKendall(y)
            out[1,reg]=tmp[1]$tau
            out_sig[1,reg]=tmp[2]$sl
            lm.r=lm(y~x)
            out[2,reg]=summary(lm.r)$coefficients[2]
            out_sig[2,reg]=summary(lm.r)$coefficients[8]
        }
    }
    return(list(out=out,out_sig=out_sig))
}

duration_region <- function(regions,reg,dur,dur_mid){
    # combines all recorded durations of one region to one duration array, same for dur_mid
    inside=which(regions==reg)
    duration=array(NA,dim=c(100000))
    duration_mid=array(NA,dim=c(100000))
    count=1
    # combines the recorded periods from all the grid points of one region in one array
    for (i in inside){
        values=length(which(!is.na(dur[i,])))
        duration[count:(count+values)]=dur[i,1:values]
        duration_mid[count:(count+values)]=dur_mid[i,1:values]
        count=count+values
    }
    duration=duration[!is.na(duration)]
    duration_mid=duration_mid[!is.na(duration_mid)]
    return(list(duration=duration,duration_mid=duration_mid))
}

duration_regional_quantile <- function(dur,dur_mid,regions,yearPeriod,regNumb){
    # dur and dur_mid is an array of dim=c(1319,# of periods)

    taus=c(0.75,0.9,0.95,0.98)
    out=array(NA,c(length(taus),regNumb))
    out_sig=array(NA,c(length(taus),regNumb))
    for (reg in 1:regNumb){
        tmp=duration_region(regions,reg,dur,dur_mid)
        duration=tmp$duration
        duration_mid=tmp$duration_mid
        ord=order(duration_mid)
        if (length(duration)>1000){
            # performs a quantile regression for each region
            for (i in 1:length(taus)){
                print(paste(reg,taus[i]))
                y=as.vector(duration[ord])
                x=as.vector(duration_mid[ord])
                inYearPeriod=which(x>yearPeriod[1] & x<yearPeriod[2])
                y=y[inYearPeriod]
                x=x[inYearPeriod]

                qu=try(summary(rq(y~x,taus[i]))$coefficients)
                if (class(qu)=="try-error"){qu=array(NA,8)}

                out[i,reg]=qu[2]
                out_sig[i,reg]=qu[8]
            }
        }
    }
    return(list(out=out,out_sig=out_sig))
}

regional_trends <- function(dat,yearPeriod,filepath,region_name){
    # performs the entire regional analysis of markov and duration
    # result will be written in ncdf file

    # pnly for one trend and 2 states until now

    nday=91
    nyr=5
    ntot=length(dat$ID)
    library(Kendall)
    library(quantreg)

    nc_mar=open.ncdf("../data/91_5/2_states/markov/91_5_markov2s.nc")
    IDregions=read.table("../data/ID-regions.txt")

    poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))

    regs=dim(poli)[1]
    IDregions=points_to_regions(dat,c(region_name))
    print(IDregions)

    season_names=c("spring","summer","autumn","winter","year")
    for (season in 1:length(season_names)){   
        print(season_names[season]) 
        result=array(NA,dim=c(2,10,regs))
        sig=array(NA,dim=c(2,10,regs))
        nc_dur=open.ncdf(paste("../data/",nday,"_",nyr,"/2_states/duration/",nday,"_",nyr,"_duration_2s_",season_names[season],".nc",sep=""))
        dur=get.var.ncdf(nc_dur,"dur")
        dur_mid=get.var.ncdf(nc_dur,"dur_mid")
        per=get.var.ncdf(nc_mar,paste("markov_",season_names[season],sep=""))
        for (state in 1:2){
            tmp=markov_regional_trend(dat,per[1:ntot,(state*state),],IDregions,yearPeriod=yearPeriod,regNumb=regs)
            result[state,1:2,1:regs]=tmp$out
            sig[state,1:2,1:regs]=tmp$out_sig
            print(tmp)
            tmp=duration_regional_quantile(dat,dur[1:ntot,state,],dur_mid[1:ntot,state,],IDregions,yearPeriod=yearPeriod,regNumb=regs)
            result[state,3:6,1:regs]=tmp$out
            sig[state,3:6,1:regs]=tmp$out_sig
            print(tmp)

        }
        regional_analysis_write(paste(filepath,season_names[season],"_",region_name,".nc",sep=""),result,sig,poli)
    }


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

duration_regional_distribution <- function(dur,dur_mid,regions,yearPeriod,regNumb,maxDur){
    # dur and dur_mid is an array of dim=c(1319,# of periods)
    breaks=seq(0,maxDur,1)
    density=array(NA,dim=c(regNumb,maxDur))
    quantiles=array(NA,dim=c(regNumb,10))
    for (reg in 1:regNumb){
        tmp=duration_region(regions,reg,dur,dur_mid)
        duration=tmp$duration
        duration_mid=tmp$duration_mid
        ord=order(duration_mid)
        if (length(duration)>1000){
            y=as.vector(duration[ord])
            x=as.vector(duration_mid[ord])
            inYearPeriod=which(x>yearPeriod[1] & x<yearPeriod[2])
            y=y[inYearPeriod]
            x=x[inYearPeriod]
            histo=hist(duration,breaks,plot=FALSE)
            density[reg,]=histo$density
            quantiles[reg,1:6]=quantile(duration,probs=c(0.05,0.25,0.5,0.75,0.95,1))
            quantiles[reg,10]=mean(duration,na.rm=TRUE)
            quantiles[reg,9]=sd(duration,na.rm=TRUE)
        }
    }
    return(list(density=density,quantiles=quantiles))
}

regional_climatology <- function(trendID,dat,yearPeriod,region_name){
    # performs the entire regional analysis of markov and duration
    # result will be written in ncdf file

    # pnly for one trend and 2 states until now
    maxDur=200
    print(seq(0,(maxDur-1),1)+0.5)
    ntot=length(dat$ID)
    library(quantreg)

    IDregions=read.table("../data/ID-regions.txt")

    poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))

    regNumb=dim(poli)[1]
    IDregions=points_to_regions(dat,c(region_name))

    print(IDregions)

    season_names=c("spring","summer","autumn","winter","year")

    distributions=array(NA,dim=c(length(season_names),2,regNumb,maxDur))
    quantiles=array(NA,dim=c(length(season_names),2,regNumb,10))

    for (season in 1:length(season_names)){   
        print(season_names[season]) 

        nc_dur=open.ncdf(paste("../data/",trendID,"/2_states/duration/",trendID,"_duration_2s_",season_names[season],".nc",sep=""))
        dur=get.var.ncdf(nc_dur,"dur")
        dur_mid=get.var.ncdf(nc_dur,"dur_mid")
        for (state in 1:2){
            tmp=duration_regional_distribution(dur[1:ntot,state,],dur_mid[1:ntot,state,],IDregions,yearPeriod=yearPeriod,regNumb=regNumb,maxDur=maxDur)
            distributions[season,state,,]=tmp$density
            quantiles[season,state,,]=tmp$quantiles
        }
    }

    ncRegion <- dim.def.ncdf("region",units="region",vals=1:regNumb, unlim=FALSE)
    ncStates <- dim.def.ncdf("states",units="states",vals=1:2,unlim=FALSE)
    ncSeason <- dim.def.ncdf("seasons",units="seasons",vals=1:5,unlim=FALSE)

    mids=seq(0,(maxDur-1),1)+0.5
    ncMids <- dim.def.ncdf("mids",units="days",vals=mids,unlim=FALSE)
    ncOther <- dim.def.ncdf("other",units="0.05,0.25,0.5,0.75,0.95,1,NA,NA,SD,Mean",vals=1:10,unlim=FALSE)
    
    poli_points <- dim.def.ncdf("poli_points",units="id",vals=1:12,unlim=FALSE)

    region_coordinates <- var.def.ncdf(name="region_coordinates",units="deg",longname="1:6 lon - 7:12 lat",dim=list(ncRegion,poli_points),missval=-9999.0)

    ncDensity <- var.def.ncdf(name="density",units="density 0-1",longname="histogramm density of durations recorded in region",dim=list(ncSeason,ncStates,ncRegion,ncMids), missval=-9999.0)
    ncQuantile <- var.def.ncdf(name="quantile",units="density 0-1",longname="0.05,0.25,0.5,0.75,0.95,1,NA,NA,SD,Mean",dim=list(ncSeason,ncStates,ncRegion,ncOther), missval=-9999.0)
    
    vars=list(ncDensity,ncQuantile,region_coordinates)
   
    nc = create.ncdf(paste("../data/",trendID,"/2_states/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",region_name,"_distributions.nc",sep=""),vars)
    put.var.ncdf(nc,ncDensity,distributions)      
    put.var.ncdf(nc,ncQuantile,quantiles)      

    pol_poi=array(NA,c(dim(poli)[1],12))
    for (i in 1:dim(poli)[1]){
        for (j in 1:12){
            
            if (is.numeric(poli[i,j])){
                pol_poi[i,j]=poli[i,j]
            }
        }
    }
    put.var.ncdf(nc,region_coordinates,pol_poi)      

    close.ncdf(nc) 

}

plot_regional_boxplots <- function(trendID,dat,yearPeriod,region_name){
    nc=open.ncdf(paste("../data/",trendID,"/2_states/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",region_name,"_distributions.nc",sep=""))
    regions=get.var.ncdf(nc,"region")
    quantiles=get.var.ncdf(nc,"quantile")

    color=c(rgb(0,0,1,0.6),rgb(1,0,0,0.4))

    plot_names=c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s")

    region_names=c("western \n N. America","central \n N. America","eastern \n N. America","Europe","western \n Asia","central \n Asia","eastern \n Asia")
    region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")

    season_names=c("spring","summer","autumn","winter","year")
    season_names=c("MAM","JJA","SON","DJF","year")
    season_auswahl=c(5,1,2,3,4)

    pdf(file=paste("../plots/",trendID,"/2_states/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_boxplots.pdf",sep=""),width=8,height=12)
    par(mfrow=c(5,2))

    for (sea in season_auswahl){
        print(sea)
        print(season_names[sea])
        for (state in 1:2){
            plot(NA,xlim=c(0.5,length(regions)),ylim=c(-3,24),frame.plot=FALSE,axes=FALSE,ylab="# days",xlab="",main=season_names[sea])
            axis(2,ylim=c(-3,24))
            for (reg in regions){
                quAn=quantiles[sea,state,reg,]
                mitte=reg
                links=mitte-0.2
                rechts=mitte+0.2
                polygon(x=c(rechts,links,links,rechts),y=c(quAn[2],quAn[2],quAn[4],quAn[4]),col=color[state])
                lines(c(mitte,mitte),c(quAn[1],quAn[5]))
                lines(c(links,rechts),c(quAn[1],quAn[1]))
                lines(c(links,rechts),c(quAn[5],quAn[5]))
                lines(c(links,rechts),c(quAn[3],quAn[3]))
                points(mitte,quAn[10],pch=4)
                text(mitte,-3,region_names[reg])
            }
        }
    }
    graphics.off()
}


plot_regional_distributions <- function(trendID,dat,yearPeriod,region_name){

    nc=open.ncdf(paste("../data/",trendID,"/2_states/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",region_name,"_distributions.nc",sep=""))
    distributions=get.var.ncdf(nc,"density")
    regions=get.var.ncdf(nc,"region")
    mids=get.var.ncdf(nc,"mids")
    quantiles=get.var.ncdf(nc,"quantile")

    color=c(rgb(0,0,1,0.6),rgb(1,0,0,0.4))

    plot_names=c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s")

    region_names=c("western \n N. America","central \n N. America","eastern \n N. America","Europe","western \n Asia","central \n Asia","eastern \n Asia")
    region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")

    season_names=c("spring","summer","autumn","winter","year")
    season_names=c("year","MAM","JJA","SON","DJF")
    season_auswahl=c(5,1,2,3,4)

    state_names=c("warm","cold")

    pdf(file=paste("../plots/",trendID,"/2_states/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_distributions.pdf",sep=""),width=8,height=12)

    par(cex.lab=0.5,cex.axis=0.5)

    plot(NA,xlim=c(0,100),ylim=c(0,100),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)


    breite=1/(length(season_auswahl)+2)
    hoehe=1/(length(regions)+2)



    for (sea in season_auswahl){
        for (reg in regions){
            if (reg==1){
                xPos=breite*sea
                yPos=hoehe*length(regions)
                par(new=TRUE,plt=c(xPos,xPos+breite,yPos,yPos+hoehe))
                plot(NA,xlim=c(0,30),ylim=c(0,0.3),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
                axis(3, col="black",col.axis="black")
                xPos=breite*sea
                yPos=hoehe*length(regions)+hoehe
                par(new=TRUE,plt=c(xPos,xPos+breite,yPos,yPos+hoehe))
                plot(NA,xlim=c(0,30),ylim=c(0,0.3),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
                text(x=15,y=0.1,label=season_names[sea])
            }
            if (sea==1){
                xPos=breite
                yPos=hoehe*reg+hoehe/4
                par(new=TRUE,plt=c(xPos,xPos+breite,yPos,yPos+3/4*hoehe))
                plot(NA,xlim=c(0,30),ylim=c(0,0.3),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
                axis(2, col="black",col.axis="black")
                xPos=0.0
                yPos=hoehe*reg+hoehe/4
                par(new=TRUE,plt=c(xPos,xPos+breite,yPos,yPos+3/4*hoehe))
                plot(NA,xlim=c(0,30),ylim=c(0,0.3),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
                text(x=15,y=0.15,label=region_names[reg])
            }


            xPos=breite*sea
            yPos=hoehe*reg+hoehe/4
            par(new=TRUE,plt=c(xPos,xPos+breite,yPos,yPos+3/4*hoehe))
            plot(NA,xlim=c(0,30),ylim=c(0,0.3),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
            text(x=25,y=0.25,label=((7-reg)*length(season_auswahl)+sea),col=rgb(0.6,0.6,0.6,0.9))

            for (br in 1:100){
                y=distributions[sea,1,reg,br]
                polygon(x=c(br+0.45,br-0.45,br-0.45,br+0.45),y=c(0,0,y,y),border=NA,col=color[1])
            }
            par(new=TRUE)
            for (br in 1:100){
                y=distributions[sea,2,reg,br]
                polygon(x=c(br+0.45,br-0.45,br-0.45,br+0.45),y=c(0,0,y,y),border=NA,col=color[2])
            }                

            xPos=breite*sea
            yPos=hoehe*reg
            par(new=TRUE,plt=c(xPos,xPos+breite,yPos,yPos+1/4*hoehe))
            plot(NA,xlim=c(0,30),ylim=c(0,10),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
            for (state in 1:2){
                quAn=quantiles[sea,state,reg,]
                mitte=4+(state-1)*3
                links=mitte-1
                rechts=mitte+1
                polygon(x=c(quAn[2],quAn[2],quAn[4],quAn[4]),y=c(rechts,links,links,rechts),col=color[state])
                lines(c(quAn[1],quAn[5]),c(mitte,mitte))
                lines(c(quAn[1],quAn[1]),c(links,rechts))
                lines(c(quAn[5],quAn[5]),c(links,rechts))
                lines(c(quAn[3],quAn[3]),c(links,rechts))
                points(quAn[10],mitte,pch=4)
            }  
        }
    }


    graphics.off()

}

