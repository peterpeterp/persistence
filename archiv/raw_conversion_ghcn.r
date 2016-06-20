
station_to_region_attribution <- function(region_name="ward23"){
    library(ff)
    # dat from HadGHCND is required for reagion attribution!!!!
    # loads attribution file and visualizes regions on map
    attribution<-read.table(paste("../data/_TMean/ID_regions/",region_name,".txt",sep=""))[,1]
    mids<-read.table(paste("../data/_TMean/ID_regions/",region_name,"_mids.txt",sep=""))  

    stations<<-read.csv.ffdf(file="../data/raw_data/ghcn/ghcnd-stations_ID_lat_lon.txt",sep="\t",fill=TRUE)    

    station_attribution=array(NA,c(dim(stations)[1],4))
    station_attribution[,1]=stations[,1]
    station_attribution[,3]=stations[,2]
    station_attribution[,4]=stations[,3]

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
    station_info<-array(NA,c(length(inside),3))
    station_info[,1]=inside
    station_info[,2]=station_attribution[inside,3]
    station_info[,3]=station_attribution[inside,4]
    write.table(station_info,paste("../data/raw_data/ghcn/ghcnd-stations_",region_name,"_",G,"_lat_lon.txt",sep=""))    
}

extract_precip <- function(yr=2013,region_name="ward23",G=7){
    library(ff)

    selection<-read.table(paste("../data/raw_data/ghcn/ghcnd-stations_",region_name,"_",G,"_ID.txt",sep=""))
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

merge_years <- function(yrs=1950:2015,region_name="ward23",G=7){
    library(RNetCDF)
    ntot=16157
    pp=array(NA,c(ntot,365,66))
    for (yr in yrs){
        cat(paste("-",yr))
        tmp<-read.table(paste("../data/raw_data/ghcn/toMerge/precip_",region_name,"_",G,"_",yr,"_toMerge.txt",sep=""))
        for (day in 1:365){
            pp[,day,(yr-1949)]=tmp[,day]
        }
    }  
    selection<-read.table(paste("../data/raw_data/ghcn/ghcnd-stations_",region_name,"_",G,"_lat_lon.txt",sep=""))

    nc_out=create.nc("../data/raw_data/ghcn/ghcn_pp_1950-2015.nc")

    dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
    dim.def.nc(nc_out,"days",dimlength=365,unlim=FALSE)
    dim.def.nc(nc_out,"years",dimlength=66,unlim=FALSE)

    var.def.nc(nc_out,"day","NC_SHORT",c(1))
    var.def.nc(nc_out,"year","NC_SHORT",c(2))
    var.def.nc(nc_out,"ID","NC_SHORT",c(0))
    var.def.nc(nc_out,"lon","NC_FLOAT",c(0))
    var.def.nc(nc_out,"lat","NC_FLOAT",c(0))


    var.def.nc(nc_out,"pp","NC_SHORT",c(0,1,2))
    att.put.nc(nc_out, "pp", "missing_value", "NC_SHORT", -999)

    var.put.nc(nc_out,"ID",1:ntot)  
    var.put.nc(nc_out,"lon",selection[,3]) 
    var.put.nc(nc_out,"lat",selection[,2]) 
    var.put.nc(nc_out,"day",1:365)  
    var.put.nc(nc_out,"year",1:66)  
    var.put.nc(nc_out,"pp",pp)     
    close.nc(nc_out)    
}

record_length_check <- function(){
    source("load.r")
    source("map_plot.r")
    dat<<-dat_load_precipitation("../data/raw_data/ghcn/ghcn_pp_1950-2015_reg7.nc")
    ntot=length(dat$ID)
    missingRatio<-array(NA,c(ntot))
    asdas
    for (q in 1:ntot){
        cat("-")
        missingRatio[q]=length(which(is.na(dat$pp[q,,])))/(365*66)
    }
    missingRatio<<-missingRatio
    useableStations=which(missingRatio<=0.1)
    ntot<-length(useableStations)

    nc_out=create.nc("../data/raw_data/ghcn/ghcn_pp_1950-2015_reg7_<10.nc")
    att.put.nc(nc_out, "NC_GLOBAL", "explanation", "NC_CHAR", "only stations with less than 10% missing")

    dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
    dim.def.nc(nc_out,"days",dimlength=365,unlim=FALSE)
    dim.def.nc(nc_out,"years",dimlength=66,unlim=FALSE)

    var.def.nc(nc_out,"day","NC_SHORT",c(1))
    var.def.nc(nc_out,"year","NC_SHORT",c(2))
    var.def.nc(nc_out,"ID","NC_SHORT",c(0))
    var.def.nc(nc_out,"lon","NC_FLOAT",c(0))
    var.def.nc(nc_out,"lat","NC_FLOAT",c(0))


    var.def.nc(nc_out,"pp","NC_SHORT",c(0,1,2))
    att.put.nc(nc_out, "pp", "missing_value", "NC_SHORT", -999)

    var.put.nc(nc_out,"ID",1:ntot)  
    var.put.nc(nc_out,"lon",dat$lon[useableStations]) 
    var.put.nc(nc_out,"lat",dat$lat[useableStations]) 
    var.put.nc(nc_out,"day",1:365)  
    var.put.nc(nc_out,"year",1950:2015)  
    var.put.nc(nc_out,"pp",dat$pp[useableStations,,])     
    close.nc(nc_out)  
}

plot_init <- function(){
    paper<<-c(8,5)
    yAusschnitt<<-c(35,60)
    xAusschnitt<<-c(-100,-50)
    asp<<-1
    pointsize<<-0.44
    pch_points<<-c(1,NA,0.1,0.1)

    pch_sig<<-4
    col_sig<<-rgb(0.1,0.1,0.1,0.6)
    cex_sig<<-0.03

    region<<-NA

    season_auswahl<<-c(1,2,3,4,5)
    sub_zusatz<<-c("75th","95th","99th")
    name_reg_zusatz<<-""

    col_row<<-c(1,1)
    mat<<-NA
    layout_mat<<-c(NA)
    subIndex<<-c("a","b")
}



# prepare regional information (with python script together)
#station_to_region_attribution()
#create_region_list()

#id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#print(id)
#toDos<-1950:2015

#if (id<=length(toDos)){
#    yr<-toDos[id]
#    print(yr)
#    extract_precip(yr=yr)
#    convert_tables(yr=yr)
#}

#merge_years(yrs=1950:2015)

plot_init()

#record_length_check()
topo_map_plot(filename_plot=paste("../data/raw_data/ghcn/ghcn_pp_1950-2015_data_coverage_all.pdf",sep=""),reihen=array(1:length(dat$ID),c(1,length(dat$ID))),titel=c("stations in reg 7"),farb_mitte=c(0,1),farb_palette="regenbogen")

adas


reihen_full<-reihen
reihen_full[reihen_full<0.9]=NA
topo_map_plot(filename_plot=paste("../data/raw_data/ghcn/ghcn_pp_1950-2015_data_coverage+90.pdf",sep=""),reihen=reihen_full,titel=c("not missing percentage"),farb_mitte=c(0,1),farb_palette="regenbogen")

#record_length_view()

#topo_map_plot(filename_plot=paste("../data/raw_data/ghcn/ghcn_pp_1950-2015_data_loc.pdf",sep=""),reihen=array(1:length(dat$ID),c(1,length(dat$ID))),titel=c("not missing percentage"),farb_mitte=c(0,1),farb_palette="regenbogen")

#tmp<-read.table("/home/peter/Dokumente/pik/backuped/data/raw_data/ghcn/ghcnd-stations_ward23_7_lat_lon.txt")
#pdf("../data/raw_data/ghcn/ghcnd-stations_ward23_7_lat_lon.pdf")
#plot(topoWorld,xlim=xAusschnitt,ylim=yAusschnitt,asp=asp,location="none",col.land=rgb(0,0,0,0),col.water=rgb(0,0,0,0),mar=c(2,3,2,5))
#plot(tmp[,3],tmp[,2])
#graphics.off()
