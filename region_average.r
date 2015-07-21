
points_to_regions <- function(dat,filename="../data/SREX_regions_all.csv"){
    srex <- read.csv(file=filename, header=TRUE, sep=",")
    latpos=c(3,5,7,9,11,13)
    lonpos=c(4,6,8,10,12,14)
    poli=array(NA,dim=c(30,12))
    for (i in 1:30){
        if (srex[i,2] < 27){
            lat=c()
            lon=c()
            k=0
            for (j in latpos){
                if (srex[i,j]!=9999){
                    k=k+1
                    lat[k]=srex[i,j]
                }
            }
            k=0
            for (j in lonpos){
                if (srex[i,j]!=9999){
                    k=k+1
                    lon[k]=srex[i,j]
                }
            }
            poli[i,1:(length(lon))]=lon
            poli[i,7:(6+length(lat))]=lat
        }
    }
    dat$region = dat$ID*NA
    points=cbind(x=dat$lon,y=dat$lat)
    for (k in 1:dim(poli)[1]){
        if (is.na(poli[k,1])==FALSE){
            poligon=cbind(x=poli[k,1:6][which(is.na(poli[k,1:6])==FALSE)],y=poli[k,7:12][which(is.na(poli[k,7:12])==FALSE)])
            print(poligon)
            inside=pnt.in.poly(points,poligon)$pip
            dat$region[which(inside==1)]=srex[k,2]

        }
    }
    dat_write(filename="../data/dat_regional.nc",dat)
    return(poli)
}



regions_color <- function(values,worldmap,title){
    poli=read.table("../data/srex_poligons.txt")
    jet.colors <- colorRampPalette( c( "violet","blue","white","yellow","red") )
    nbcol <- 101
    color <- jet.colors(nbcol)

    y=c()
    index=c()
    j=0
    for (i in 1:26){
        if (!is.na(values[i])){
            j=j+1
            y[j]=values[i]
            index[j]=i           
        }
    }
    aushol=max(c(abs(max(y)),abs(min(y))))
    y[j+1]=-aushol
    y[j+2]=aushol
    facetcol <- cut(y,nbcol)  

    plot(worldmap,main=title)

    for (i in 1:j){
        lon=poli[index[i],1:6]
        lat=poli[index[i],7:12]
        lon=lon[!is.na(lon)]
        lat=lat[!is.na(lat)]
        polygon(x=lon,y=lat,col=color[facetcol[i]],border="green")
        text(mean(lon),mean(lat),label=poli[index[i],13],cex=0.7,col="black")
    }
    image.plot(legend.only=T, zlim=range(y), col=color)
    return()
}

average_regional_simple <- function(dat,value){
	valIn=array(NA,26)
    errIn=array(NA,26)
	for (reg in 1:26){
		inside=which(dat$region==reg)
		valIn[reg]=mean(value[inside],na.rm=TRUE)
        errIn[reg]=sd(value[inside],na.rm=TRUE)

	}
	return(list(mean=valIn,sd=errIn))
}

average_regional_trend <- function(dat,value){
    # value is an array of dim=c(1319,??)
    library(Kendall)
    size=dim(value)[2]
    valOut=array(NA,dim=c(26,size))
    for (reg in 1:26){
        inside=which(dat$region==reg)
        for (j in 1:size){
            valOut[reg,j]=mean(value[inside,j],na.rm=TRUE)
        }
    }
    slope=array(NA,26)
    slope_sig=array(NA,26)
    MK=array(NA,26)
    MK_sig=array(NA,26)
    x=seq(1,65,1)
    for (reg in 1:26){
        y=valOut[reg,]
        if (length(which(is.nan(y)==1))<10){
            print(y)
            print(length(which(is.nan(y)==1)))
            lm.r=lm(y~x)
            slope[reg]=summary(lm.r)$coefficients[2]
            slope_sig[reg]=summary(lm.r)$coefficients[8]
            out=MannKendall(y)
            MK[reg]=out[1]$tau
            MK_sig[reg]=out[2]$sl
        }
    }
    return(list(MK=MK,MK_sig=MK_sig,slope=slope,slope_sig=slope_sig))
}

map_regional <- function(dat,toPlot,titles,worldmap,filename_plot){
    # dat from data_load, loaded with reg=1!
    # toPlot array(... ,dim=c(number of maps, number of stations))
    # titles titles of plots (list)

    ntot=length(dat$ID)
    pdf(file=filename_plot)

    for (i in 1:dim(toPlot)[1]){
        out=average_regional_trend(dat,toPlot[i,1:ntot,])
        print(out)
        regions_color(out$MK,worldmap,paste("MK",titles[i]))
        regions_color(out$MK_sig,worldmap,paste("MK_sig",titles[i]))
        regions_color(out$slope,worldmap,paste("LR",titles[i]))
        regions_color(out$slope_sig,worldmap,paste("LR_sig",titles[i]))
    }
    graphics.off()
}

