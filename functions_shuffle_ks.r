
shuffle_mat <- function(durMat,years){
    return(durMat[,,,sample(years,years,replace=FALSE)])
}

ks_analysis <- function(seasons=1,id=2,nShuffle=1000,actual_nShuffle=nShuffle){
    print(id)
    for (sea in seasons){
        season<-season_names[sea]

        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_reg_binned_duration_",season,".nc",sep="")
        nc=open.nc(filename)
        binned_dur<-var.get.nc(nc,"binned_dur")
        periodsInYr<-dim(binned_dur)[3]
        ID_length<<-dim(binned_dur)[1]
        ID_select<<-1:ID_length

        original<-ks_statistics(durMat=binned_dur)

        shuffled<-array(NA,c(nShuffle,ID_length,2,2))
        for (shuff in 1:actual_nShuffle){
        	cat(paste("--",shuff))
        	shuffled[shuff,,,]<-ks_statistics(durMat=shuffle_mat(binned_dur,years))
        }

        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/shuffled_ks/",trendID,dataset,"_",ID_name,"_",period,"_shuffled_ks_",season,"_",id,".nc",sep="")
        nc_out <- create.nc(filename)
        att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
            
        dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
        dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
        dim.def.nc(nc_out,"outs",dimlength=2,unlim=FALSE)
        dim.def.nc(nc_out,"nShuffle",dimlength=nShuffle,unlim=FALSE)

        var.def.nc(nc_out,"original","NC_DOUBLE",c(0,1,2))
        att.put.nc(nc_out, "original", "missing_value", "NC_DOUBLE", 99999)
        att.put.nc(nc_out, "original", "dim_explanation", "NC_CHAR", "ID-states-outs")
        att.put.nc(nc_out, "original", "val_explanation", "NC_CHAR", "ks-statistic for 1-whole distribution, 2-upper tail distribution")

        var.def.nc(nc_out,"shuffled","NC_DOUBLE",c(3,0,1,2))
        att.put.nc(nc_out, "shuffled", "missing_value", "NC_DOUBLE", -99999.9)
        att.put.nc(nc_out, "shuffled", "dim_explanation", "NC_CHAR", "shuflle-ID-states-outs")
        att.put.nc(nc_out, "shuffled", "val_explanation", "NC_CHAR", "ks-statistic for 1-whole distribution, 2-upper tail distribution")

        var.put.nc(nc_out,"original",original)              
        var.put.nc(nc_out,"shuffled",shuffled)              

        close.nc(nc_out) 
        print("---------------original------------")
        print(original)    
        print("-----------------------------------")  
    }
}

ks_statistics <- function(durMat){
    ks_stats<-array(NA,c(ID_length,2,2))
    for (q in ID_select){
        for (state in 1:2){
            y1<-durMat[q,state,,startYear:halfYear]
            y2<-durMat[q,state,,(halfYear+1):endYear]

            ks_stats[q,state,1]<-ks.test(y1,y2)$statistic

            distr_cut<-median(c(y1,y2),na.rm=TRUE)
            y1<-y1[y1>=distr_cut]
            y2<-y2[y2>=distr_cut]
            ks_stats[q,state,2]<-ks.test(y1,y2)$statistic
        }
    }
    return(ks_stats)
}

init <- function(){
    library(quantreg)
    library(RNetCDF)
    source("analysis_tools.r")
    nday<<-91
    nyr<<-7
    trendID<<-paste(nday,"_",nyr,sep="")
    dataset<<-"_TMean"
    additional_style<<-""
    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
}

init()

ID_name<-"ward24"
folder<-paste("/regional/",ID_name,"/",sep="")
period<-"1950-2014"

startYear<-1
halfYear<-32
endYear<-65
years<-65

id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(id)


# MAM, JJA, SON, DJF, 4seasons
if (TRUE){
    if (id<11){ks_analysis(seasons=1,id=(id))}
    if (id<21 & id>10){ks_analysis(seasons=2,id=(id-10))}
    if (id<31 & id>20){ks_analysis(seasons=3,id=(id-20))}
    if (id<41 & id>30){ks_analysis(seasons=4,id=(id-30))}
    if (id<51 & id>40){ks_analysis(seasons=5,id=(id-40))}
}
