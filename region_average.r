source("write.r")
source("load.r")
library(SDMTools)


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
    nbcol <- 20
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

map_regional <- function(dat,toPlot,titles){
    ntot=length(dat$ID)
    library(rworldmap)
    library(fields)
    worldmap = getMap(resolution = "low")
    pdf(file="../plots/regions/region_map.pdf")

    for (i in 1:dim(toPlot)[2]){
        out=average_regional(dat,toPlot[1:ntot,i])
        regions_color(out,worldmap,titles[i])
    }
}

if (1==2){
	dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2011.nc",reg=0)
	points_to_regions(dat)
}


average_regional <- function(dat,value){
	valIn=array(NA,26)
    errIn=array(NA,26)
	for (reg in 1:26){
		inside=which(dat$reg==reg)
		valIn[reg]=mean(value[inside],na.rm=TRUE)
        errIn=sum(value[inside],na.rm=TRUE)

	}
	return(valIn)
}

if (1==1){
    nday=91
    nyr=5
	dat=dat_load("../data/dat_regional.nc",reg=1)

    per=markov_load(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr))
}