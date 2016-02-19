
shuffle_mat <- function(durMat,years){
    return(durMat[,,,sample(years,years,replace=FALSE)])
}

trend_analysis <- function(seasons=2,yearPeriod=c(1950,2014),folder="/regional/",ID_name="7rect"){
    for (sea in seasons){
        season<-season_names[sea]

        print(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_reg_binned_duration_",season,".nc",sep=""))
        nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_reg_binned_duration_",season,".nc",sep=""))
        binned_dur<-var.get.nc(nc,"binned_dur")
        periodsInYr<-dim(binned_dur)[3]
        ID_length<-dim(binned_dur)[1]

        # create time vector for trend analysis
        years<-yearPeriod[2]-yearPeriod[1]
        time_vec<-rep(0:years,each=periodsInYr)
        binned_dur[,,,(yearPeriod[1]-1949+1):(yearPeriod[2]-1949)]

        original<-trend_evaluation(durMat=array(binned_dur,c(ID_length,2,periodsInYr*65)),time_vec=time_vec,ID_select=1:ID_length)

        nShuffle<-100
        shuffled<-array(NA,c(nShuffle,ID_length,2,9))
        for (shuff in 1:nShuffle){
        	cat(paste("--",shuff))
        	shuffled[shuff,,,]<-trend_evaluation(durMat=array(shuffle_mat(binned_dur,years),c(ID_length,2,periodsInYr*65)),time_vec=time_vec,ID_select=1:ID_length)
        }

        nc_out <- create.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,"/shuffled/",trendID,dataset,"_",ID_name,"_",yearPeriod[1],"-",yearPeriod[2],"_shuffled_trends_",season,"_",id,".nc",sep=""))
        att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
            
        dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
        dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
        dim.def.nc(nc_out,"outs",dimlength=9,unlim=FALSE)
        dim.def.nc(nc_out,"nShuffle",dimlength=nShuffle,unlim=FALSE)

        var.def.nc(nc_out,"original","NC_INT",c(0,1,2))
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
    }
}

trend_evaluation <- function(durMat,time_vec,ID_select=1:1319,ID_length=length(ID_select)){    #,noise_level=c(0,0)
    trends<-array(NA,c(ID_length,2,9))
    #time_vec=time_vec+rnorm(length(time_vec),mean=0,sd=1)*noise_level[1]
    for (q in ID_select){
        for (state in 1:2){
            y<-durMat[q,state,]
            if (length(which(!is.na(y)))>100){
                #y=y+rnorm(length(y),mean=0,sd=1)*noise_level[2]
                trends[q,state,1:7]=rq(y~time_vec,taus)$coef[2,]
                trends[q,state,9]=lm(y~time_vec)$coef[2]
            }
        }
    }
    return(trends)

}

init <- function(){
    library(quantreg)
    library(RNetCDF)


    nday<<-91
    nyr<<-5
    trendID<<-paste(nday,"_",nyr,sep="")
    dataset<<-"_TMean"
    additional_style<<-""

    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
}

id   = Sys.getenv("SLURM_ARRAY_TASK_ID")
print(id)
init()
taus<-c(0,0.05,0.25,0.5,0.75,0.95,1)
trend_analysis()