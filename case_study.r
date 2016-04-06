
master_init <- function(){
    source("functions_support.r")
    source("functions_duration.r")
    source("functions_regional.r")
    source("analysis_tools.r")
    source("write.r")
    source("load.r")
    source("plot_master.r")
    source("map_plot.r")
    source("inits_plot.r")


    library(moments)
    library(quantreg)
    library(stats4)
    library(RNetCDF)
    library(SDMTools)
    library(fields)
}

master_init()

match_grids <- function(){
    datPP <<- dat_load_precipitation(paste("../data/_eobsPP/rr_0.50deg_reg_v12.0_1950-2015.nc",sep=""))
    datTG <<- dat_load(filename<-paste("../data/_eobsTG/tg_0.50deg_reg_v12.0_end.nc",sep=""))


    matches<-array(NA,c(80000,2))
    index<-0
    for (q in 1:length(datTG$ID)){
        cat("-")
        for (Q in 1:length(datPP$ID))
            if (datTG$lon[q] %in% datPP$lon[Q] & datTG$lat[q] %in% datPP$lat[Q]){matches[index<-index+1,1:2]=c(q,Q)}
    }

    nc_out <- create.nc(paste("../data/_eobs/matches.nc",sep=""))

    dim.def.nc(nc_out,"ID",dimlength=index, unlim=FALSE)
    dim.def.nc(nc_out,"sets",dimlength=2, unlim=FALSE)

    var.def.nc(nc_out,"matches","NC_SHORT",c(0,1))
    att.put.nc(nc_out, "matches", "missing_value", "NC_SHORT", -999)
    att.put.nc(nc_out, "matches", "dim_explanation", "NC_CHAR", "ID-(TG,PP)")
    att.put.nc(nc_out, "matches", "method", "NC_CHAR","bla")

    var.put.nc(nc_out,"matches",matches[1:index,])

    close.nc(nc_out)  
}



heat_wav_detection <- function(sea=2,stateTemp=2,stateRain=1,filename_plot="/home/peter/Dokumente/pik/backuped/plots/_eobsPP/heatwave_map.pdf"){
    other_temp<-var.get.nc(open.nc("/home/peter/Dokumente/pik/backuped/data/_eobsTG/91_7/gridded/1950-2015/91_7__eobsTG_1950-2015_others.nc"),"other_stuff")
    other_rain<-var.get.nc(open.nc("/home/peter/Dokumente/pik/backuped/data/_eobsPP/0p5/gridded/1950-2015/0p5__eobsPP_1950-2015_others.nc"),"other_stuff")

    matches<-var.get.nc(open.nc(paste("../data/_eobs/matches.nc",sep="")),"matches")

    dat<<-list(lon=datPP$lon[matches[,2]],lat=datPP$lat[matches[,2]])

    simultaneous<-array(0,c(1,dim(matches)[1]))

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,4])==sign(other_rain[sea,matches[,2],stateRain,4]))
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,4][simultan])*(-1)

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,4])==sign(other_rain[sea,matches[,2],stateRain,4]) & (other_temp[sea,matches[,1],stateTemp,10]<0.1 | other_rain[sea,matches[,2],stateRain,10]<0.1))
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,4][simultan])*2*(-1)

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,4])==sign(other_rain[sea,matches[,2],stateRain,4]) & (other_temp[sea,matches[,1],stateTemp,10]<0.1 & other_rain[sea,matches[,2],stateRain,10]<0.1))
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,4][simultan])*3*(-1)

    topo_map_plot(filename_plot=paste(filename_plot,"_lr.pdf",sep=""),reihen=array(simultaneous,c(1,ntot)),farb_mitte="0",farb_palette="dry-wet",ID_select=1:dim(matches)[1])


    other_temp<-var.get.nc(open.nc("/home/peter/Dokumente/pik/backuped/data/_eobsTG/91_7/gridded/1950-2015/91_7__eobsTG_1950-2015_quantiles.nc"),"quantile_stuff")[,,,3,]
    other_rain<-var.get.nc(open.nc("/home/peter/Dokumente/pik/backuped/data/_eobsPP/0p5/gridded/1950-2015/0p5__eobsPP_1950-2015_quantiles.nc"),"quantile_stuff")[,,,3,]

    matches<-var.get.nc(open.nc(paste("../data/_eobs/matches.nc",sep="")),"matches")

    dat<<-list(lon=datPP$lon[matches[,2]],lat=datPP$lat[matches[,2]])

    simultaneous<-array(0,c(1,dim(matches)[1]))

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,2])==sign(other_rain[sea,matches[,2],stateRain,2]))
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,2][simultan])*(-1)

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,2])==sign(other_rain[sea,matches[,2],stateRain,2]) & (other_temp[sea,matches[,1],stateTemp,3]<0.1 | other_rain[sea,matches[,2],stateRain,3]<0.1))
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,2][simultan])*2*(-1)

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,2])==sign(other_rain[sea,matches[,2],stateRain,2]) & (other_temp[sea,matches[,1],stateTemp,3]<0.1 & other_rain[sea,matches[,2],stateRain,3]<0.1))
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,2][simultan])*3*(-1)

    topo_map_plot(filename_plot=paste(filename_plot,"_qu95.pdf",sep=""),reihen=array(simultaneous,c(1,ntot)),farb_mitte="0",farb_palette="dry-wet",ID_select=1:dim(matches)[1])    
}

#match_grids()

#datPP <<- dat_load_precipitation(paste("../data/_eobsPP/rr_0.50deg_reg_v12.0_1950-2015.nc",sep=""))

source("map_plot.r")

plot_init_EU()
nbcol<<-7

heat_wav_detection(sea=2,stateTemp=2,stateRain=1,filename_plot="/home/peter/Dokumente/pik/backuped/plots/_eobsPP/heatwave_map")