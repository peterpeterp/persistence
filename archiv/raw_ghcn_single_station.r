
extract_station <- function(yr=2013,station="GM000003342"){
    tmp<<-read.csv(paste("../data/raw_data/ghcn/",yr,".csv",sep=""))
    print("loaded")

    for (i in 1:length(indices)){
        index<-indices[i]
        vals<-which(tmp[,3]==index & tmp[,1]==station)
        print(vals)
        dat=data.frame(id=tmp[vals,1],date=tmp[vals,2],values=tmp[vals,4])
        filename=paste("../data/raw_data/ghcn/",station,"/",station,"_",yr,"_",index,".txt",sep="") ; print(filename)
        write.table(dat,filename)
    }
    
    #rm(tmp)
}

convert_tables <- function(yr=2013,station="GM000003342"){

    dat_vec<-c(1:31+100,1:28+200,1:31+300,1:30+400,1:31+500,1:30+600,1:31+700,1:31+800,1:30+900,1:31+1000,1:30+1100,1:31+1200)

	pp=array(NA,c(365,length(indices)))

    for (i in 1:length(indices)){
        index<-indices[i]        
        dat<-read.table(paste("../data/raw_data/ghcn/",station,"/",station,"_",yr,"_",index,".txt",sep=""))

        if (dim(dat)[1]>1){
            dat[,2]<-dat[,2]-yr*10000
            pp[which(dat_vec %in% dat[,2]),i]=dat[,3][which(dat[,2]!=229)]
        }
	}
    write.table(pp,paste("../data/raw_data/ghcn/",station,"/",station,"_",yr,"_toMerge.txt",sep=""))
}

merge_years <- function(yrs=1950:2015,station="GM000003342"){
    library(RNetCDF)
    values=array(NA,c(length(indices),365,66))
    for (yr in yrs){
        cat(paste("-",yr))
        tmp<-read.table(paste("../data/raw_data/ghcn/",station,"/",station,"_",yr,"_toMerge.txt",sep=""))
        for (i in 1:length(indices)){
            for (day in 1:365){
                values[i,day,(yr-1949)]=tmp[day,i]
            }
        }
    }  

    nc_out=create.nc(paste("../data/raw_data/ghcn/ghcn_",station,"_1950-2015.nc",sep=""))

    dim.def.nc(nc_out,"indices",dimlength=length(indices), unlim=FALSE)
    dim.def.nc(nc_out,"days",dimlength=365,unlim=FALSE)
    dim.def.nc(nc_out,"years",dimlength=66,unlim=FALSE)


    var.def.nc(nc_out,"values","NC_SHORT",c(0,1,2))
    att.put.nc(nc_out, "values", "missing_value", "NC_SHORT", -999)
    att.put.nc(nc_out, "values", "val_explanation", "NC_CHAR", paste(indices))
 
    var.put.nc(nc_out,"values",values)     
    close.nc(nc_out)    
}

indices<-c("PRCP","SNOW","SNWD","TMAX","TMIN")
merge_years(yrs=1950:2014,station="GM000003342")
adas

id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(id)
toDos<-1950:2015

if (id<=length(toDos)){
    yr<-toDos[id]
    print(yr)
    extract_station(yr=yr,station="GM000003342")
    convert_tables(yr=yr,station="GM000003342")
}

