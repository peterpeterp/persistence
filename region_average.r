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

trend_control <- function(dat,per,seasonStart=c(59,151,243,335),seasonStop=c(150,242,334,424)){
    ntot=length(dat$ID)
    waTrend=array(NA,dim=c(ntot,8))
    for (q in 1:ntot){
        for (sea in 1:4){
            warmeTage=array(NA,62)
            kalteTage=array(NA,62)
            for (i in 1:61){
                if (seasonStop[sea]>365){
                    warmeTage[i]=length(which(per$ind[q,seasonStart[sea]:365,i]==1))
                    kalteTage[i]=length(which(per$ind[q,seasonStart[sea]:365,i]==-1)) 

                    warmeTage[i]=warmeTage[i]+length(which(per$ind[q,1:(seasonStop[sea]-365),(i+1)]==1))
                    kalteTage[i]=kalteTage[i]+length(which(per$ind[q,1:(seasonStop[sea]-365),(i+1)]==-1)) 
                }
                else {
                    warmeTage[i]=length(which(per$ind[q,seasonStart[sea]:seasonStop[sea],i]==1))
                    kalteTage[i]=length(which(per$ind[q,seasonStart[sea]:seasonStop[sea],i]==-1))  
                }   
                if ((warmeTage[i]+kalteTage[i])<80){
                    warmeTage[i]=NA
                    kalteTage[i]=NA
                }            
            }
            if (length(which(is.na(warmeTage)))<30){ 
                lm.r=lm(warmeTage~dat$year)
                waTrend[q,sea]=summary(lm.r)$coefficients[2]
                waTrend[q,(4+sea)]=summary(lm.r)$coefficients[8]
            }
        }
    }
    write.table(waTrend,"../data/warmeTage_trends_4seasons.txt")
    return(waTrend[1:ntot,1:4])
}

regions_color <- function(values,worldmap){
    poli=read.table("../data/poligons.txt")
    print(poli)

    jet.colors <- colorRampPalette( c( "violet","blue","white","yellow","red") )
    nbcol <- 20
    color <- jet.colors(nbcol)
    print(values)
    values = values[!is.na(values)]
     
    size=length(values)       
    y=array(NA,(size+2))
    y[3:(size+2)]=values

    y[1]=-2*sd(values)
    y[2]=2*sd(values)
    facetcol <- cut(y,nbcol)  

    print(y)

    plot(worldmap)

    for (i in 1:26){
        if (!is.na(values[i])){
            lon=poli[i,1:6]
            lat=poli[i,7:12]
            lon=lon[!is.na(lon)]
            lat=lat[!is.na(lat)]
            polygon(x=lon,y=lat,col=color[facetcol[(2+i)]],border="green")
            text(mean(lon),mean(lat),cex=0.7,col="black")
        }
    }
    image.plot(legend.only=T, zlim=range(values), col=color)
    write.table(poli,"poligons.txt")
    return()
}

map_regional <- function(dat,toPlot){
    ntot=length(dat$ID)
    library(rworldmap)
    library(fields)
    worldmap = getMap(resolution = "low")
    pdf(file="../plots/regions/region_map.pdf")
    for (i in 1:dim(toPlot)[2]){
        out=average_regional(dat,toPlot[1:ntot,i])
        regions_color(out,worldmap)
    }
}

if (1==2){
	dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2011.nc",reg=0)
	points_to_regions(dat)
}


average_regional <- function(dat,value){
	out=array(NA,26)
	for (reg in 1:26){
		inside=which(dat$reg==reg)
		out[reg]=mean(value[inside],na.rm=TRUE)
	}
	return(out)
}

if (1==1){
    nday=91
    nyr=5


	dat=dat_load("../data/dat_regional.nc",reg=1)
    per=markov_load(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr))
    #watrends=trend_control(dat,per)
    watrends=read.table("../data/warmeTage_trends_4seasons.txt")
    map_regional(dat,watrends)
    sdfs
	print(average_regional(dat,dat$lon))
}