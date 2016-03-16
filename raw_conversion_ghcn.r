
station_to_region_attribution <- function(region_name="ward23"){
    library(ff)

    # loads attribution file and visualizes regions on map
    attribution<-read.table(paste("../data/_TMean/ID_regions/",region_name,".txt",sep=""))[,1]
    mids<-read.table(paste("../data/_TMean/ID_regions/",region_name,"_mids.txt",sep=""))  

    stations<<-read.csv.ffdf(file="../data/raw_data/ghcn/ghcnd-stations_ID.txt",sep="\t",fill=TRUE)    

    station_attribution=array(NA,c(dim(stations)[1],2))
    station_attribution[,1]=1:dim(stations)[1]

    regNumb<-dim(mids)[1]
    for (G in 1:regNumb){
        inside<-which(attribution==G)
        for (q in inside){
            cat("-")
            station_attribution[which(stations[,3]<dat$lon[q]+1.875 & stations[,3]>dat$lon[q]-1.875 & stations[,2]<dat$lat[q]+1.25 & stations[,2]>dat$lat[q]-1.25),2]=G
        }
    }
    write.table(station_attribution,paste("../data/raw_data/ghcn/ghcnd-stations_",region_name,".txt",sep=""))
}

create_region_list <- function(region_name="ward23",G=7){
    station_attribution<<-read.table(paste("../data/raw_data/ghcn/ghcnd-stations_",region_name,".txt",sep=""))  
    inside<<-which(station_attribution[,2]==G)
    write.table(inside,paste("../data/raw_data/ghcn/ghcnd-stations_",region_name,"_",G,".txt",sep=""))      
}

extract_precip <- function(yr=2013,region_name="ward23",G=7){
    library(ff)

    selection<-read.table(paste("../data/raw_data/ghcn/ghcnd-stations_",region_name,"_",G,"_ID.txt",sep=""))
    #tmp<-read.csv.ffdf(file=paste("../data/raw_data/ghcn/",yr,".csv",sep=""))
    tmp<-read.csv(paste("../data/raw_data/ghcn/",yr,".csv",sep=""))
    print("loaded")
    prcp<-which(tmp[,3]=="PRCP" & tmp[,1] %in% selection[,2])


    dat=data.frame(id=tmp[prcp,1],date=tmp[prcp,2],prcp=tmp[prcp,4])

    print(paste("../data/raw_data/ghcn/precip_",region_name,"_",G,"_",yr,".txt",sep=""))
    write.table(dat,paste("../data/raw_data/ghcn/precip_",region_name,"_",G,"_",yr,".txt",sep=""))
    print(paste("../data/raw_data/ghcn/precip_",region_name,"_",G,"_",yr,".txt",sep=""))

    rm(tmp)
    rm(selection)
}

convert_tables <- function(yr=2013,region_name="ward23",G=7){
	dat<<-read.table(paste("../data/raw_data/ghcn/precip_",region_name,"_",G,"_",yr,".txt",sep=""))
	dat[,2]<<-dat[,2]-yr*10000
    selection<-read.table(paste("../data/raw_data/ghcn/ghcnd-stations_",region_name,"_",G,"_ID.txt",sep=""))

	ntot<-dim(selection)[1]
    dat_vec<<-c(1:31+100,1:28+200,1:31+300,1:30+400,1:31+500,1:30+600,1:31+700,1:31+800,1:30+900,1:31+1000,1:30+1100,1:31+1200)

	pp=array(NA,c(ntot,365))

    percentage<-0
    cat(paste("\n0 -> -> -> -> -> 100\n"))
    for (q in 1:ntot){
        #print(proc.time())
        if (q/ntot*100 > percentage){
            cat("-")
            percentage<-percentage+5
        }
		ofStation<<-which(dat[,1]==paste(selection[q,2]))

        if (length(ofStation)>2){
            pp[q,which(dat_vec %in% dat[ofStation,2])]=dat[ofStation,3][which(dat[ofStation,2]!=229)]
        }
	}
    write.table(pp,paste("../data/raw_data/ghcn/toMerge/precip_",region_name,"_",G,"_",yr,"_toMerge.txt",sep=""))
}


# prepare regional information (with python script together)
#station_to_region_attribution()
#create_region_list()


id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(id)
toDos<-c(1956,1960,1964,1968,1972,1976,1980,1984,1992,1996,1999:2015)
toDos<-c(1956,1960,1964,1968,1972,1976,1980,1984,1992,1996)


if (id<=length(toDos)){
    yr<-toDos[id]
    print(yr)
    #extract_precip(yr=yr)
    convert_tables(yr=yr)
}

