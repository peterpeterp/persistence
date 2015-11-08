
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


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


plot_regional_distributions <- function(trendID,dat,yearPeriod,region_name,additional_style,dataset){

    nc=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",region_name,"_",yearPeriod[1],"-",yearPeriod[2],"_distributions.nc",sep=""))
    regions=get.var.ncdf(nc,"region")
    mids=get.var.ncdf(nc,"mids")
    quantiles=array(get.var.ncdf(nc,"quantile"),dim=c(5,2,length(regions),10))
    distributions=array(get.var.ncdf(nc,"density"),dim=c(5,2,length(regions),length(mids)))


    color=c(rgb(0,0,1,0.6),rgb(1,0,0,0.4))

    plot_names=c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s")

    region_names=c("western \n N. America","central \n N. America","eastern \n N. America","Europe","western \n Asia","central \n Asia","eastern \n Asia")
    region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")

    # differnet order in season_auswahl will lead to wrong plots! dangerous
    season_names=c("spring","summer","autumn","winter","year")
    season_names=c("year","MAM","JJA","SON","DJF")
    season_names=c("MAM","JJA","SON","DJF","year")
    season_auswahl=c(1,2,3,4,5)

    state_names=c("warm","cold")

    # regional focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_distributions_reg.pdf",sep=""),width=8,height=12)
    par(cex.lab=0.5,cex.axis=0.5)
    plot(NA,xlim=c(0,100),ylim=c(0,100),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)

    breite=1/(length(regions)+2)
    hoehe=1/(length(season_auswahl)+2)

    for (reg in regions){
        for (sea in season_auswahl){
            if (sea==1){
                xPos=breite*reg
                yPos=hoehe*length(season_auswahl)
                par(new=TRUE,plt=c(xPos,xPos+breite,yPos,yPos+hoehe))
                plot(NA,xlim=c(0,30),ylim=c(0,0.3),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
                axis(3, col="black",col.axis="black")
                xPos=breite*reg
                yPos=hoehe*length(season_auswahl)+hoehe
                par(new=TRUE,plt=c(xPos,xPos+breite,yPos,yPos+hoehe))
                plot(NA,xlim=c(0,30),ylim=c(0,0.3),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
                text(x=15,y=0.1,label=region_names[reg])
            }
            if (reg==1){
                xPos=breite
                yPos=hoehe*(length(season_auswahl)-sea+1)+hoehe/4
                par(new=TRUE,plt=c(xPos,xPos+breite,yPos,yPos+3/4*hoehe))
                plot(NA,xlim=c(0,30),ylim=c(0,0.3),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
                axis(2, col="black",col.axis="black")
                xPos=0.0
                yPos=hoehe*(length(season_auswahl)-sea+1)+hoehe/4
                par(new=TRUE,plt=c(xPos,xPos+breite,yPos,yPos+3/4*hoehe))
                plot(NA,xlim=c(0,30),ylim=c(0,0.3),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
                text(x=15,y=0.15,label=season_names[sea])
            }


            xPos=breite*reg
            yPos=hoehe*(length(season_auswahl)-sea+1)+hoehe/4
            par(new=TRUE,plt=c(xPos,xPos+breite,yPos,yPos+3/4*hoehe))
            plot(NA,xlim=c(0,30),ylim=c(0,0.3),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
            text(x=25,y=0.25,label=((sea-1)*length(regions)+reg),col=rgb(0.6,0.6,0.6,0.9))
            #text(x=15,y=0.25,label=paste(reg,sea),col=rgb(0.6,0.6,0.6,0.9))

            for (verticals in c(5,10,15,20,25)){
                abline(v=verticals,col=rgb(0.6,0.6,0.6,0.6))
            } 

            for (br in 1:100){
                y=distributions[sea,1,reg,br]
                br=br-0.5
                polygon(x=c(br+0.45,br-0.45,br-0.45,br+0.45),y=c(0,0,y,y),border=NA,col=color[1])
            }
            par(new=TRUE)
            for (br in 1:100){
                y=distributions[sea,2,reg,br]
                br=br-0.5
                polygon(x=c(br+0.45,br-0.45,br-0.45,br+0.45),y=c(0,0,y,y),border=NA,col=color[2])
            }               
            xPos=breite*reg
            yPos=hoehe*(length(season_auswahl)-sea+1)
            par(new=TRUE,plt=c(xPos,xPos+breite,yPos,yPos+1/4*hoehe))
            plot(NA,xlim=c(0,30),ylim=c(0,10),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
            for (verticals in c(5,10,15,20,25)){
                abline(v=verticals,col=rgb(0.6,0.6,0.6,0.6))
            }
            for (state in 1:2){
                quAn=quantiles[sea,state,reg,]
                mitte=4+(state-1)*3
                unten=mitte-1
                oben=mitte+1
                polygon(x=c(quAn[2],quAn[2],quAn[4],quAn[4]),y=c(oben,unten,unten,oben),col=color[state])
                lines(c(quAn[1],quAn[5]),c(mitte,mitte))
                lines(c(quAn[1],quAn[1]),c(unten,oben))
                lines(c(quAn[5],quAn[5]),c(unten,oben))
                lines(c(quAn[3],quAn[3]),c(unten,oben))
                points(quAn[10],mitte,pch=4,cex=0.8)
            }    
        }
    }

    #seasonal focus
    graphics.off()

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_distributions_sea.pdf",sep=""),width=8,height=12)
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
            for (verticals in c(5,10,15,20,25)){
                abline(v=verticals,col=rgb(0.6,0.6,0.6,0.6))
            } 

            for (br in 1:100){
                y=distributions[sea,1,reg,br]
                br=br-0.5
                polygon(x=c(br+0.45,br-0.45,br-0.45,br+0.45),y=c(0,0,y,y),border=NA,col=color[1])
            }
            par(new=TRUE)
            for (br in 1:100){
                y=distributions[sea,2,reg,br]
                br=br-0.5
                polygon(x=c(br+0.45,br-0.45,br-0.45,br+0.45),y=c(0,0,y,y),border=NA,col=color[2])
            }               

            xPos=breite*sea
            yPos=hoehe*reg
            par(new=TRUE,plt=c(xPos,xPos+breite,yPos,yPos+1/4*hoehe))
            plot(NA,xlim=c(0,30),ylim=c(0,10),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
            for (verticals in c(5,10,15,20,25)){
                abline(v=verticals,col=rgb(0.6,0.6,0.6,0.6))
            }
            for (state in 1:2){
                quAn=quantiles[sea,state,reg,]
                mitte=4+(state-1)*3
                unten=mitte-1
                oben=mitte+1
                polygon(x=c(quAn[2],quAn[2],quAn[4],quAn[4]),y=c(oben,unten,unten,oben),col=color[state])
                lines(c(quAn[1],quAn[5]),c(mitte,mitte))
                lines(c(quAn[1],quAn[1]),c(unten,oben))
                lines(c(quAn[5],quAn[5]),c(unten,oben))
                lines(c(quAn[3],quAn[3]),c(unten,oben))
                points(quAn[10],mitte,pch=4,cex=0.8)
            }    
        }
    }


    graphics.off()

}


#============================================================================================================================
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

            # y are now the duration length in selected yearPeriod
            histo=hist(y,breaks,plot=FALSE)
            density[reg,]=histo$density

            quantiles[reg,1:6]=quantile(y,probs=c(0.05,0.25,0.5,0.75,0.95,1))
            quantiles[reg,10]=mean(y,na.rm=TRUE)
            quantiles[reg,9]=sd(y,na.rm=TRUE)
        }

    }
    return(list(density=density,quantiles=quantiles))
}

regional_climatology <- function(dat,yearPeriod,region_name,trendID,additional_style,dataset){
    # performs the entire regional analysis of markov and duration
    # result will be written in ncdf file

    # pnly for one trend and 2 states until now
    maxDur=200
    #print(seq(0,(maxDur-1),1)+0.5)
    ntot=length(dat$ID)
    library(quantreg)

    IDregions=read.table("../data/ID-regions.txt")
    poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))
    regNumb=dim(poli)[1]
    IDregions=points_to_regions(dat,c(region_name))
    season_names=c("MAM","JJA","SON","DJF","year")

    distributions=array(NA,dim=c(length(season_names),2,regNumb,maxDur))
    quantiles=array(NA,dim=c(length(season_names),2,regNumb,10))

    for (season in 1:length(season_names)){   
        print(season_names[season])

        nc_dur=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_",season_names[season],".nc",sep=""))
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
   
    nc = create.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",region_name,"_",yearPeriod[1],"-",yearPeriod[2],"_distributions.nc",sep=""),vars)
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
#==========================================================================================================================

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


plot_regional_boxplots <- function(dat,yearPeriod,region_name,trendID,additional_style,dataset){
    # performs the entire regional analysis of markov and duration
    # result will be written in ncdf file

    # pnly for one trend and 2 states until now
    ntot=length(dat$ID)

    nc = open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",region_name,"_",yearPeriod[1],"-",yearPeriod[2],"_quantiles.nc",sep=""))
    quantiles=get.var.ncdf(nc,"quantiles")
    others=get.var.ncdf(nc,"others")
    regions=get.var.ncdf(nc,"region")
    regNumb=7
    seaNumb=6
    region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")
    season_names=c("MAM","JJA","SON","DJF","year","4seasons")

    color=c(rgb(0.5,0.5,1,0.8),rgb(1,0.5,0.5,0.8))
    pos=c(-1,1)*0.15
    width=6
    height=3

    # regional focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_boxplots_regional.pdf",sep=""),width=width,height=height)
    par(mfrow=c(1,1))
    par(mar=c(1,4,2,3))   
    at_=seq(1, regNumb, 1)
    at_=c(at_-0.15,at_+0.15)
    buchstaben=c("a","b","c","d")
    drueber=2
    for (sea in 1:length(season_names)){   
        season=season_names[sea]

        plot(NA,xlim=c(0,7.5),ylim=c(0,(max(quantiles[sea,,,1:5])+drueber)),frame.plot=FALSE,axes=FALSE,ylab="days")
        axis(2)
        for (i in axis(2)){
            abline(h=i,col=rgb(0.5,0.5,0.5,0.5),lty=3)
        }
        text(8,(max(quantiles[sea,,,1:5],na.rm=TRUE)+drueber-0.5),label=buchstaben[sea],cex=2)

        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,],reg+pos[state],0.3,color[state])
                text(reg,max(quantiles[sea,,,1:5])+drueber,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))
            }
        }        
        #for (oth in c(1,7)){
        #    points(at_[1:regNumb],others[sea,,1,oth],col="blue",pch=oth)
        #    lines(at_[1:regNumb],others[sea,,1,oth],col="blue",lty=3)
        #    points(at_[(regNumb+1):(regNumb*2)],others[sea,,2,oth],col="red",pch=oth)
        #    lines(at_[(regNumb+1):(regNumb*2)],others[sea,,2,oth],col="red",lty=3)
        #}



        #for (i in 1:8){
        #    others[sea,,1,i]=(others[sea,,1,i]-min(others[sea,,1,i],na.rm=TRUE))/(max(others[sea,,1,i],na.rm=TRUE)-min(others[sea,,1,i],na.rm=TRUE))
        #    others[sea,,2,i]=(others[sea,,2,i]-min(others[sea,,2,i],na.rm=TRUE))/(max(others[sea,,2,i],na.rm=TRUE)-min(others[sea,,2,i],na.rm=TRUE))
        #}
        for (oth in c(7,11)){
            par(new=TRUE)
            plot(NA,xlim=c(0,7.5),ylim=c(min(others[sea,,,oth],na.rm=TRUE),max(others[sea,,,oth],na.rm=TRUE)),frame.plot=FALSE,axes=FALSE,ylab="")
            axis(4,col="green",ylab="R2")
            points(at_[1:regNumb],others[sea,,1,oth],col="blue",pch=oth)
            lines(at_[1:regNumb],others[sea,,1,oth],col="blue",lty=3)
            points(at_[(regNumb+1):(regNumb*2)],others[sea,,2,oth],col="red",pch=oth)
            lines(at_[(regNumb+1):(regNumb*2)],others[sea,,2,oth],col="red",lty=3)
        }
    }
    graphics.off()

    # seasonal anomaly
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_boxplots_seasonal_anomaly.pdf",sep=""),width=width,height=height)
    par(mfrow=c(1,1))
    par(mar=c(1,4,2,0))   
    at_=seq(1, regNumb, 1)
    at_=c(at_-0.15,at_+0.15)
    buchstaben=c("e","f","g","h","i","j")
    drueber=0.5

    anomaly=quantiles*NA
    for (sea in 1:4){  
        anomaly[sea,,,]=quantiles[sea,,,]-quantiles[6,,,]
    }
    for (sea in 1:4){   
        season=season_names[sea]
        plot(NA,xlim=c(0,8),ylim=c(min(anomaly[sea,,,1:5],na.rm=TRUE),(max(anomaly[sea,,,1:5],na.rm=TRUE)+drueber)),frame.plot=FALSE,axes=FALSE,ylab="days")
        axis(2)
        for (i in axis(2)){
            abline(h=i,col=rgb(0.5,0.5,0.5,0.5),lty=3)
        }        
        text(8,(max(anomaly[sea,,,1:5],na.rm=TRUE)+drueber-0.5),label=buchstaben[sea],cex=2)

        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(anomaly[sea,reg,state,],reg+pos[state],0.3,color[state])
                text(reg,max(anomaly[sea,,,1:5],na.rm=TRUE)+drueber,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))
            }
        }
        for (quA in c(4,5)){
            points(at_[1:regNumb],anomaly[sea,,1,quA],col="black",pch=2)
            points(at_[(regNumb+1):(regNumb*2)],anomaly[sea,,2,quA],col="black",pch=2)
        }
    }

    pdf(file="../plots/zwischenzeugs/pch_palette.pdf")
    plot(1:30,1:30,pch=1:30)
    graphics.off()    

    # regional focus kind of seasonal anomaly
    color=c(rgb(0.5,0.5,1,0.8),rgb(1,0.5,0.5,0.8),rgb(0.5,0.5,1,0.4),rgb(1,0.5,0.5,0.4))

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_boxplots_seas_kind_anomaly.pdf",sep=""))
    par(mfrow=c(1,1))
    par(mar=c(1,4,2,0))   
    at_=seq(1, regNumb, 1)
    at_=c(at_-0.15,at_+0.15)

    for (sea in 1:4){   
        season=season_names[sea]

        plot(NA,xlim=c(0,8),ylim=c(0,(max(quantiles[sea,,,1:5])+4)),frame.plot=FALSE,axes=FALSE,main=season,ylab="days")
        axis(2)
        for (i in axis(2)){
            abline(h=i,col=rgb(0.5,0.5,0.5,0.5),lty=3)
        }        
        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[6,reg,state,],reg+pos[state],0.3,color=color[2+state],border=rgb(0.5,0.5,0.5,0.5))
            }
        }
        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,],reg+pos[state],0.3,color=color[2+state])
                text(reg,max(quantiles[sea,,,1:5])+2,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))
            }
        }

    }
    graphics.off()

    # regional focus plus seasonal anomaly
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_boxplots_seas_combi_anomaly.pdf",sep=""))
    par(mfrow=c(2,1),cex=1)
    par(mar=c(1,4,2,0))   
    at_=seq(1, regNumb, 1)
    at_=c(at_-0.15,at_+0.15)
    buchstaben=c("a","b","c","d","e","f","g","h","i","j")
    drueber=c(1,1.5)
    index=0
    for (sea in 1:4){   
        season=season_names[sea]
        par(mar=c(1,4,2,0))
        index=index+1
        plot(NA,xlim=c(0,8),ylim=c(min(anomaly[sea,,,1:5],na.rm=TRUE),(max(anomaly[sea,,,1:5],na.rm=TRUE)+drueber[1])),frame.plot=FALSE,axes=FALSE,ylab="days")#,main=season)
        axis(2)
        for (i in axis(2)){
            abline(h=i,col=rgb(0.5,0.5,0.5,0.5),lty=3)
        }        
        text(8,(max(anomaly[sea,,,1:5],na.rm=TRUE)+drueber[1]),label=buchstaben[index],cex=2)
        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(anomaly[sea,reg,state,],reg+pos[state],0.3,color[state])
                #text(reg,max(anomaly[sea,,,1:5],na.rm=TRUE)+0.5,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))
            }
        }

        for (quA in c(4,5)){
            points(at_[1:regNumb],anomaly[sea,,1,quA],col="black",pch=2)
            points(at_[(regNumb+1):(regNumb*2)],anomaly[sea,,2,quA],col="black",pch=2)
        }
        par(mar=c(1,4,1,0))
        index=index+1
        plot(NA,xlim=c(0,8),ylim=c(0,(max(quantiles[sea,,,1:5])+drueber[2])),frame.plot=FALSE,axes=FALSE,ylab="days")
        text(8,(max(quantiles[sea,,,1:5],na.rm=TRUE)+drueber[2]),label=buchstaben[index],cex=2)
        axis(2)
        for (i in axis(2)){
            abline(h=i,col=rgb(0.5,0.5,0.5,0.5),lty=3)
        }   
        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,],reg+pos[state],0.3,color=color[state])
                text(reg,max(quantiles[sea,,,1:5])+drueber[2],region_names[reg],col=rgb(0.5,0.5,0.5,0.5),cex=1.3)
            }
        }


    }
    graphics.off()
}

plot_regional_boxplots_vergleich <- function(dat,yearPeriod1,yearPeriod2,region_name,trendID,additional_style,dataset){
    # performs the entire regional analysis of markov and duration
    # result will be written in ncdf file

    # pnly for one trend and 2 states until now
    ntot=length(dat$ID)

    nc = open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",yearPeriod1[1],"-",yearPeriod1[2],"/",trendID,"_",region_name,"_",yearPeriod1[1],"-",yearPeriod1[2],"_quantiles.nc",sep=""))
    quantiles1=get.var.ncdf(nc,"quantiles")
    nc = open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",yearPeriod2[1],"-",yearPeriod2[2],"/",trendID,"_",region_name,"_",yearPeriod2[1],"-",yearPeriod2[2],"_quantiles.nc",sep=""))
    quantiles2=get.var.ncdf(nc,"quantiles")
    regions=get.var.ncdf(nc,"region")
    regNumb=7
    seaNumb=6
    region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")
    season_names=c("MAM","JJA","SON","DJF","year","4seasons")

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/","/",yearPeriod2[1],"-",yearPeriod2[2],"_diff_",yearPeriod1[1],"-",yearPeriod1[2],"_boxplots_regional.pdf",sep=""))
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

    # regional focus
    for (sea in 1:length(season_names)){   
        season=season_names[sea]

        plot(NA,xlim=c(0,8),ylim=c(min(quantiles[sea,,,1:5]),(max(quantiles[sea,,,1:5])+0.5)),frame.plot=FALSE,axes=FALSE,main=season,ylab="days")
        axis(2)
        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,],reg+pos[state],0.3,color[state])
                text(reg,max(quantiles[sea,,,1:5])+0.5,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))
            }
        }
        for (quA in c(9)){
            points(at_[1:regNumb],quantiles[sea,,1,quA],col="blue",pch=1)
            lines(at_[1:regNumb],quantiles[sea,,1,quA],col="blue",lty=3)
            points(at_[(regNumb+1):(regNumb*2)],quantiles[sea,,2,quA],col="red",pch=1)
            lines(at_[(regNumb+1):(regNumb*2)],quantiles[sea,,2,quA],col="red",lty=3)
        }
        for (quA in c(4,5)){
            points(at_[1:regNumb],quantiles[sea,,1,quA],col="black",pch=2)
            points(at_[(regNumb+1):(regNumb*2)],quantiles[sea,,2,quA],col="black",pch=2)
        }
    }
    graphics.off()

    # seasonal focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/","/",yearPeriod2[1],"-",yearPeriod2[2],"_diff_",yearPeriod1[1],"-",yearPeriod1[2],"_boxplots_seasonal.pdf",sep=""))
    par(mfrow=c(1,1))
    par(mar=c(1,5,4,3))   
    at_=seq(1, seaNumb, 1)
    at_=c(at_-0.15,at_+0.15)

    for (reg in 1:regNumb){
        plot(NA,xlim=c(0,8),ylim=c(min(quantiles[,reg,,1:5]),(max(quantiles[,reg,,1:5])+0.5)),frame.plot=FALSE,axes=FALSE,main=region_names[reg],ylab="days")
        axis(2)
        for (sea in 1:seaNumb){   
            season=season_names[sea]     
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,],sea+pos[state],0.3,color[state])
                text(sea,max(quantiles[,reg,,1:5])+0.5,season_names[sea],col=rgb(0.5,0.5,0.5,0.5))
            }
        }
        for (quA in c(9)){
            points(at_[1:seaNumb],quantiles[,reg,1,quA],col="blue",pch=1)
            lines(at_[1:seaNumb],quantiles[,reg,1,quA],col="blue",lty=3)
            points(at_[(seaNumb+1):(seaNumb*2)],quantiles[,reg,2,quA],col="red",pch=1)
            lines(at_[(seaNumb+1):(seaNumb*2)],quantiles[,reg,2,quA],col="red",lty=3)
        }
    }
    graphics.off()
}


quantile_pete <- function(dist,taus,na.rm=TRUE){
    if (na.rm==TRUE){dist=dist[which(!is.na(dist))]}

    cdf=array(NA,max(dist))
    out=taus*NA

    cum=0
    for (i in 1:max(dist)){
        cum=cum+length(which(dist>i))
        cdf[i]=cum
    }
    cdf=cdf/cdf[length(cdf)]

    for (qu in 1:length(taus)){
        ueb=which(cdf>taus[qu])[1]
        unt=ueb-1
        if (unt<1){out[qu]=ueb}
        else {out[qu]=ueb-(cdf[ueb]-taus[qu])/(cdf[ueb]-cdf[unt])}
    }  
    return(out)
}

regional_quantiles <- function(dat,yearPeriod,region_name,trendID,additional_style,dataset){
    # performs the entire regional analysis of markov and duration
    # result will be written in ncdf file

    # pnly for one trend and 2 states until now
    ntot=length(dat$ID)


    poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))
    regNumb=dim(poli)[1]
    IDregions=points_to_regions(dat,c(region_name))

    region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")
    season_names=c("MAM","JJA","SON","DJF","year","4seasons")

    taus=c(0.05,0.25,0.5,0.75,0.95,0.91,0.98,0.1)
    quantiles=array(NA,dim=c(length(season_names),regNumb,2,length(taus)))
    quantiles_vergleich=array(NA,dim=c(length(season_names),regNumb,2,length(taus)))
    others=array(NA,dim=c(length(season_names),regNumb,2,11))

    pdf(file=paste("../plots/zwischenzeugs/reg_dist_diff_fit_plot_",yearPeriod[1],"-",yearPeriod[2],".pdf",sep=""))

    for (sea in 1:length(season_names)){   
        season=season_names[sea]
        dists=list()

        nc_dur=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_",season,".nc",sep=""))
        dur=get.var.ncdf(nc_dur,"dur")
        dur_mid=get.var.ncdf(nc_dur,"dur_mid")
        
        for (state in 1:2){
            for (reg in 1:regNumb){
            
                tmp=duration_region(IDregions,reg,dur[1:ntot,state,],dur_mid[1:ntot,state,])
                duration=tmp$duration
                duration_mid=tmp$duration_mid
                ord=order(duration_mid)
                if (length(duration)>1000){
                    y=as.vector(duration[ord])
                    x=as.vector(duration_mid[ord])
                    inYearPeriod=which(x>yearPeriod[1] & x<yearPeriod[2])
                    y=y[inYearPeriod]
                    x=x[inYearPeriod]

                    # y are now the duration length in selected yearPeriod
                    # quantile determination
                    dists=append(dists,list(y))
                    quantiles[sea,reg,state,1:length(taus)]=quantile_pete(y,taus=taus,na.rm=TRUE)
                    quantiles_vergleich[sea,reg,state,1:length(taus)]=quantile(y,prob=taus,na.rm=TRUE)

                    # exponential fit + other values
                    br=seq(0,max(y,na.rm=TRUE),1)
                    histo=hist(y,breaks=br,plot=FALSE)

                    # linear log fit
                    if (1==2){
                        dens=histo$density
                        dens[dens==0]=NA
                        fit=summary(lm(log(dens)~histo$mids))
                        A=exp(fit$coefficients[1])
                        b=-fit$coefficients[2]
                        R2=fit$r.squared
                        mean=mean(y,na.rm=TRUE)
                        sd=sd(y,na.rm=TRUE)
                        skew=skewness(y,na.rm=TRUE)

                        others[sea,reg,state,1:8]=c(mean,sd,sd/mean,skew,A,b,1/b,R2)
                        print(others[sea,reg,state,1:8])
                    }

                    Y=histo$density
                    X=histo$mids
                    xy=data.frame(y=Y,x=X)
                    fit=nls(y~(a*exp(-b*x)),data=xy,start=list(a=0.1,b=0.1),na.action=na.exclude) 
                    a=summary(fit)$parameters[1]
                    b=summary(fit)$parameters[2]
                    yfit=a*exp(-X*b)
                    R2=1-sum(((Y-yfit)^2),na.rm=TRUE)/sum(((Y-mean(Y,na.rm=TRUE))^2),na.rm=TRUE)
                    mean=mean(y,na.rm=TRUE)
                    sd=sd(y,na.rm=TRUE)
                    skew=skewness(y,na.rm=TRUE)
                    others[sea,reg,state,1:8]=c(mean,sd,sd/mean,skew,a,b,1/b,R2)

                    # plotstuff
                    par(mar=c(5, 4, 4, 4) + 0.1)
                    plot(histo$density,main=paste(season,region_names[reg],state))
                    lines(yfit,col="red")
                    #text(40,0.03,paste("R^2:",round(R2,03)))

                    z=histo$density-yfit
                    z2=z
                    z2[1:round(mean(y,na.rm=TRUE))]=NA
                    nona=!is.na(z2)
                    xy=data.frame(y=z,x=X)

                    fit=try(nls(y~(a*exp(-(b-x)^2/sigma^2)),data=xy,start=list(a=0.005,b=14,sigma=7),na.action=na.exclude))
                    if (class(fit)!="try-error"){
                        A=summary(fit)$parameters[1]
                        B=summary(fit)$parameters[2]                   
                        sigma=summary(fit)$parameters[3]   
                        print(paste(region_names[reg],a,b,sigma))
                        zfit=(a*exp(-(b-X)^2/sigma^2))      

                    }  
                    else{
                        A=0.005
                        B=14
                        sigma=8
                    }
                    zfit=(a*exp(-b*X)+A*exp(-(B-X)^2/sigma^2))      
                    lines(zfit,col="orange") 
                    A=0.01
                    B=12
                    sigma=7

                    xy=data.frame(y=Y,x=X)
                    nls.control(maxiter = 100, tol = 1e-05)
                    fit=try(nls(y~(a*exp(-b*x)+A*exp(-(B-x)^2/sigma^2)),data=xy,start=list(a=a,b=b,A=A,B=B,sigma=sigma),na.action=na.exclude,nls.control(maxiter = 100, tol = 1e-01, minFactor=1/1024)))
                    if (class(fit)!="try-error"){

                        a=summary(fit)$parameters[1]
                        b=summary(fit)$parameters[2]  
                        A=summary(fit)$parameters[3]
                        B=summary(fit)$parameters[4]                   
                        sigma=summary(fit)$parameters[5]   
                        zfit=(a*exp(-b*X)+A*exp(-(B-X)^2/sigma^2))      
                        lines(zfit,col="blue")   
                    } 





                    par(new=TRUE,plt=c(0.5,0.82,0.5,0.87))
                    z=histo$density-yfit
                    plot(z,ylim=c(-0.02,0.02))
                    abline(v=mean(y,na.rm=TRUE))
                    #text(20,0.017,round(sum(z[round(mean(y,na.rm=TRUE)):length(z)]),02))
                    z2=z
                    z2[1:round(mean(y,na.rm=TRUE))]=NA
                    points(z2,col="green")
                    #text(20,0.014,round(mean(z2,na.rm=TRUE),04))
                    nona=!is.na(z2)
                    z2=z2[nona]
                    perc=length(z[z>0])/length(z)*100
                    #text(30,0.012,round(sum(z,na.rm=TRUE),03))
                   
                    others[sea,reg,state,9:11]=c(sum(z,na.rm=TRUE),sum(z[round(mean(y,na.rm=TRUE)):length(z)]),mean(z2,na.rm=TRUE))
                    print(others[sea,reg,state,8:11])
                    
                    xy=data.frame(y=z,x=X)
                    a=0.01
                    b=7
                    sigma=5
                    #lines(a*exp(-(b-X)^2/sigma^2),col="blue")
                    fit=try(nls(y~(a*exp(-(b-x)^2/sigma^2)),data=xy,start=list(a=0.01,b=12,sigma=7),na.action=na.exclude))
                    if (class(fit)!="try-error"){
                        a=summary(fit)$parameters[1]
                        b=summary(fit)$parameters[2]                   
                        sigma=summary(fit)$parameters[3]   
                        print(paste(region_names[reg],a,b,sigma))
                        zfit=(a*exp(-(b-X)^2/sigma^2))      
                        lines(zfit,col="blue")    

                        gaussR2=1-sum(((z2-zfit)^2),na.rm=TRUE)/sum(((mean(z2,na.rm=TRUE)-zfit)^2),na.rm=TRUE) 
                        text(40,0.015,gaussR2)
                        text(40,0.017,sum(((z2-zfit)^2),na.rm=TRUE))
                        text(40,0.02,sum(((mean(z2,na.rm=TRUE)-zfit)^2),na.rm=TRUE))
                    }  

                }
            }
            adsas
        }
    }

    ncRegion <- dim.def.ncdf("region",units="region",vals=1:regNumb, unlim=FALSE)
    ncStates <- dim.def.ncdf("states",units="states",vals=1:2,unlim=FALSE)
    ncSeason <- dim.def.ncdf("seasons",units="seasons",vals=1:6,unlim=FALSE)

    ncTaus <- dim.def.ncdf("taus",units="dimension for both quantile",vals=1:8,unlim=FALSE)
    ncOuts <- dim.def.ncdf("outs",units="dimension for others",vals=1:11,unlim=FALSE)
   
    poli_points <- dim.def.ncdf("poli_points",units="id",vals=1:12,unlim=FALSE)

    region_coordinates <- var.def.ncdf(name="region_coordinates",units="deg",longname="1:6 lon - 7:12 lat",dim=list(ncRegion,poli_points),missval=-9999.0)

    ncQuantile <- var.def.ncdf(name="quantiles",units="quantile values in days evaluated by quantile_pete",longname="0.05,0.25,0.5,0.75,0.95,0.91,0.98,0.1",dim=list(ncSeason,ncRegion,ncStates,ncTaus), missval=-9999.0)
    ncQuantile_vergleich <- var.def.ncdf(name="quantiles_2",units="quantile values in days evaluated by quantile()",longname="0.05,0.25,0.5,0.75,0.95,0.91,0.98,0.1",dim=list(ncSeason,ncRegion,ncStates,ncTaus), missval=-9999.0)
    ncOthers <- var.def.ncdf(name="others",units="different analysis values",longname="mean sd sd/mean skewness A b 1/b R2 sum(deviation of density from fit) sum(deviation of density from fit[values above mean]) mean(deviation of density from fit[values above mean])",dim=list(ncSeason,ncRegion,ncStates,ncOuts), missval=-9999.0)
    
    vars=list(ncQuantile,ncQuantile_vergleich,ncOthers,region_coordinates)
   
    nc = create.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",region_name,"_",yearPeriod[1],"-",yearPeriod[2],"_quantiles.nc",sep=""),vars)
    put.var.ncdf(nc,ncQuantile,quantiles)      
    put.var.ncdf(nc,ncQuantile_vergleich,quantiles_vergleich)      
    put.var.ncdf(nc,ncOthers,others)      

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


