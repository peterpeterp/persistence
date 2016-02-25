
shuffle_mat <- function(durMat,years){
    return(durMat[,,,sample(years,years,replace=FALSE)])
}

trend_analysis <- function(seasons=4,id=2,yearPeriod=c(1980,2014),ID_name="7rect",folder="/regional/7rect/"){
    print(id)
    for (sea in seasons){
        season<-season_names[sea]

        print(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_reg_binned_duration_",season,".nc",sep=""))
        nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_reg_binned_duration_",season,".nc",sep=""))
        binned_dur<-var.get.nc(nc,"binned_dur")
        periodsInYr<-dim(binned_dur)[3]
        ID_length<-dim(binned_dur)[1]

        # create time vector for trend analysis
        years<-yearPeriod[2]-yearPeriod[1]+1
        time_vec<-rep(1:years,each=periodsInYr)
        binned_dur<-binned_dur[,,,(yearPeriod[1]-1949):(yearPeriod[2]-1949)]

        original<<-trend_evaluation(durMat=array(binned_dur,c(ID_length,2,periodsInYr*years)),time_vec=time_vec,ID_select=1:ID_length)

        nShuffle<-100
        shuffled<-array(NA,c(nShuffle,ID_length,2,5))
        for (shuff in 1:nShuffle){
        	cat(paste("--",shuff))
        	#shuffled[shuff,,,]<-trend_evaluation(durMat=array(shuffle_mat(binned_dur,years),c(ID_length,2,periodsInYr*years)),time_vec=time_vec,ID_select=1:ID_length)
        }

        period<-paste(yearPeriod[1],"-",yearPeriod[2],sep="")    
        print(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/shuffled/",trendID,dataset,"_",ID_name,"_",period,"_shuffled_trends_",season,"_",id,".nc",sep=""))
        nc_out <- create.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/shuffled/",trendID,dataset,"_",ID_name,"_",period,"_shuffled_trends_",season,"_",id,".nc",sep=""))
        att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
            
        dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
        dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
        dim.def.nc(nc_out,"outs",dimlength=5,unlim=FALSE)
        dim.def.nc(nc_out,"nShuffle",dimlength=nShuffle,unlim=FALSE)

        var.def.nc(nc_out,"original","NC_DOUBLE",c(0,1,2))
        att.put.nc(nc_out, "original", "missing_value", "NC_DOUBLE", 99999)
        att.put.nc(nc_out, "original", "dim_explanation", "NC_CHAR", "ID-states-outs")
        att.put.nc(nc_out, "original", "val_explanation", "NC_CHAR", "trends for quantiles and mean")

        var.def.nc(nc_out,"shuffled","NC_DOUBLE",c(3,0,1,2))
        att.put.nc(nc_out, "shuffled", "missing_value", "NC_DOUBLE", -99999.9)
        att.put.nc(nc_out, "shuffled", "dim_explanation", "NC_CHAR", "shuflle-ID-states-outs")
        att.put.nc(nc_out, "shuffled", "val_explanation", "NC_CHAR", "random trends for quantiles and mean")

        var.put.nc(nc_out,"original",original)              
        var.put.nc(nc_out,"shuffled",shuffled)              

        close.nc(nc_out) 
        print("---------------original------------")
        print(original)    
        print("-----------------------------------")  
    }
}

trend_evaluation <- function(durMat,time_vec,ID_select=1:1319,ID_length=length(ID_select)){
    trends<-array(NA,c(ID_length,2,5))
    for (q in ID_select){
        for (state in 1:2){
            y<-durMat[q,state,]
            if (length(which(!is.na(y)))>100){
                trends[q,state,1:4]=rq(y~time_vec,taus)$coef[2,]
                trends[q,state,5]=lm(y~time_vec)$coef[2]
            }
        }
    }
    return(trends)

}

init <- function(){
    library(quantreg)
    library(RNetCDF)
}



init()


nday<-91
nyr<-5
trendID<-paste(nday,"_",nyr,sep="")
dataset<-"_TMean"
additional_style<-""
taus<-c(0.5,0.75,0.95,0.99)
season_names<-c("MAM","JJA","SON","DJF","4seasons")

id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(id)

name_id<-0
for (i in 1:10){
    if (id>40){
        id<-id-40
        name_id<-name_id+10
    }
}

if (id<11){trend_analysis(seasons=1,id=(id+name_id))}
if (id<21 & id>10){trend_analysis(seasons=2,id=(id-10+name_id))}
if (id<31 & id>20){trend_analysis(seasons=3,id=(id-20+name_id))}
if (id<41 & id>30){trend_analysis(seasons=4,id=(id-30+name_id))}
