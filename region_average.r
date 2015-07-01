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

if (1==2){
	dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2011.nc",reg=0)
	points_to_regions(dat)
}


average_regional <- function(dat,value){
	out=array(NA,30)
	for (reg in 1:30){
		inside=which(dat$reg==reg)
		out[reg]=mean(value[inside])
	}
	return(out)
}

if (1==1){
	dat=dat_load("../data/dat_regional.nc",reg=1)
	print(average_regional(dat,dat$lon))
}