
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



heat_wav_detection <- function(sea=2,stateTemp=2,stateRain=1,farb_palette="dry-wet",filename_plot="/home/peter/Dokumente/pik/backuped/plots/_eobsPP/heatwave_map.pdf",stateName=c("dry","warm")){
    other_temp<-var.get.nc(open.nc("/home/peter/Dokumente/pik/backuped/data/_eobsTG/91_7/gridded/1950-2015/91_7__eobsTG_1950-2015_others.nc"),"other_stuff")
    other_rain<-var.get.nc(open.nc("/home/peter/Dokumente/pik/backuped/data/_eobsPP/0p5/gridded/1950-2015/0p5__eobsPP_1950-2015_others.nc"),"other_stuff")

    matches<-var.get.nc(open.nc(paste("../data/_eobs/matches.nc",sep="")),"matches")

    dat<<-list(lon=datPP$lon[matches[,2]],lat=datPP$lat[matches[,2]])
    ntot<<-length(matches[,2])

    simultaneous<-array(0,c(1,dim(matches)[1]))

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,4])==sign(other_rain[sea,matches[,2],stateRain,4]))
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,4][simultan])

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,4])==sign(other_rain[sea,matches[,2],stateRain,4]) & other_rain[sea,matches[,2],stateRain,10]<0.1)
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,4][simultan])*2

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,4])==sign(other_rain[sea,matches[,2],stateRain,4]) & other_temp[sea,matches[,1],stateTemp,10]<0.1)
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,4][simultan])*3

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,4])==sign(other_rain[sea,matches[,2],stateRain,4]) & (other_temp[sea,matches[,1],stateTemp,10]<0.1 & other_rain[sea,matches[,2],stateRain,10]<0.1))
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,4][simultan])*4

    simultaneous[simultaneous<0]=0
    simultaneous[simultaneous==0]=NA

    topo_map_plot(filename_plot=paste(filename_plot,"_lr.pdf",sep=""),reihen=array(simultaneous,c(1,ntot)),farb_mitte=c(0,4),farb_palette=farb_palette,ID_select=1:dim(matches)[1])
    plot(1:10,col=rgb(1,1,1,0),frame.plot=FALSE,axes=FALSE,ylab="",xlab="")
    legend("topleft",pch=c(NA,15,15,15),col=c(color[1],color[3],color[4],color[5]),legend=c("significant increase",paste(stateName[1],"pers."),paste(stateName[2],"pers."),"both"),cex=2,box.lwd=0,box.col=rgb(1,1,1,0))
    graphics.off()

    other_temp<-var.get.nc(open.nc("/home/peter/Dokumente/pik/backuped/data/_eobsTG/91_7/gridded/1950-2015/91_7__eobsTG_1950-2015_quantiles.nc"),"quantile_stuff")[,,,3,]
    other_rain<-var.get.nc(open.nc("/home/peter/Dokumente/pik/backuped/data/_eobsPP/0p5/gridded/1950-2015/0p5__eobsPP_1950-2015_quantiles.nc"),"quantile_stuff")[,,,3,]

    matches<-var.get.nc(open.nc(paste("../data/_eobs/matches.nc",sep="")),"matches")

    dat<<-list(lon=datPP$lon[matches[,2]],lat=datPP$lat[matches[,2]])

    simultaneous<-array(0,c(1,dim(matches)[1]))

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,2])==sign(other_rain[sea,matches[,2],stateRain,2]))
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,2][simultan])

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,2])==sign(other_rain[sea,matches[,2],stateRain,2]) & other_rain[sea,matches[,2],stateRain,3]<0.1)
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,2][simultan])*2

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,2])==sign(other_rain[sea,matches[,2],stateRain,2]) & other_temp[sea,matches[,1],stateTemp,3]<0.1 )
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,2][simultan])*3

    simultan<-which(sign(other_temp[sea,matches[,1],stateTemp,2])==sign(other_rain[sea,matches[,2],stateRain,2]) & (other_temp[sea,matches[,1],stateTemp,3]<0.1 & other_rain[sea,matches[,2],stateRain,3]<0.1))
    simultaneous[1,simultan]=sign(other_rain[sea,matches[,2],stateRain,2][simultan])*4

    simultaneous[simultaneous<0]=0
    simultaneous[simultaneous==0]=NA

    topo_map_plot(filename_plot=paste(filename_plot,"_qu95.pdf",sep=""),reihen=array(simultaneous,c(1,ntot)),farb_mitte=c(0,4),farb_palette=farb_palette,ID_select=1:dim(matches)[1])    
    plot(1:10,col=rgb(1,1,1,0),frame.plot=FALSE,axes=FALSE,ylab="",xlab="")
    legend("topleft",pch=c(NA,15,15,15),col=c(color[1],color[3],color[4],color[5]),legend=c("significant increase",paste(stateName[1],"pers."),paste(stateName[2],"pers."),"both"),cex=2,box.lwd=0,box.col=rgb(1,1,1,0))
    graphics.off()
}

wave_goggle <- function(sea,filename_plot="/home/peter/Dokumente/pik/backuped/plots/_eobsPP/waveGoggle_map",stateName=c("dry","warm")){
    other_temp<-var.get.nc(open.nc("/home/peter/Dokumente/pik/backuped/data/_eobsTG/91_7/gridded/1950-2015/91_7__eobsTG_1950-2015_others.nc"),"other_stuff")
    other_rain<-var.get.nc(open.nc("/home/peter/Dokumente/pik/backuped/data/_eobsPP/0p5/gridded/1950-2015/0p5__eobsPP_1950-2015_others.nc"),"other_stuff")

    matches<-var.get.nc(open.nc(paste("../data/_eobs/matches.nc",sep="")),"matches")

    dat<<-list(lon=datPP$lon[matches[,2]],lat=datPP$lat[matches[,2]])
    ntot<<-length(matches[,2])

    simultaneous<-array(NA,c(1,dim(matches)[1]))

    trends<-array(NA,c(dim(matches)[1],4))
    trends[,1:2]=sign(other_temp[sea,matches[,1],1:2,4])
    trends[,3:4]=sign(other_rain[sea,matches[,2],1:2,4])

    sigs<-array(0,c(dim(matches)[1],4))
    sigs[,1:2]=other_temp[sea,matches[,1],1:2,10]
    sigs[,3:4]=other_rain[sea,matches[,2],1:2,10]

    simultan<-which(trends[,1]>0 & trends[,4]>0 & (trends[,2]<0 | trends[,3]<0))
    simultaneous[1,simultan]=trends[simultan,1]*2

    simultan<-which(trends[,2]>0 & trends[,3]>0 & (trends[,1]<0 | trends[,4]<0))
    simultaneous[1,simultan]=trends[simultan,2]*3

    simultan<-which(trends[,1]>0 & trends[,4]>0 & trends[,2]>0 & trends[,3]>0)
    simultaneous[1,simultan]=trends[simultan,1]*1

    topo_map_plot(filename_plot=paste(filename_plot,"_lr.pdf",sep=""),reihen=array(simultaneous,c(1,ntot)),farb_mitte=c(0,3),farb_palette="wave-goggle",ID_select=1:dim(matches)[1])
    draw.circle(c(6,28.5,51),c(45,51,57),c(10,10,10),lwd=5,lty=2)
    #lines(c(-17.5,72.5),c(40,60),lty=2,cex=5)
    plot(1:10,col=rgb(1,1,1,0),frame.plot=FALSE,axes=FALSE,ylab="",xlab="")
    legend("topleft",pch=c(15,15,15),col=c(color[2],color[3],color[4]),legend=c("general increase in persistence","increase for floodings","increase for heat-waves"),cex=2,box.lwd=0,box.col=rgb(1,1,1,0))

    graphics.off()

    other_temp<-var.get.nc(open.nc("/home/peter/Dokumente/pik/backuped/data/_eobsTG/91_7/gridded/1950-2015/91_7__eobsTG_1950-2015_quantiles.nc"),"quantile_stuff")[,,,3,]
    other_rain<-var.get.nc(open.nc("/home/peter/Dokumente/pik/backuped/data/_eobsPP/0p5/gridded/1950-2015/0p5__eobsPP_1950-2015_quantiles.nc"),"quantile_stuff")[,,,3,]

    simultaneous<-array(NA,c(1,dim(matches)[1]))

    trends<-array(NA,c(dim(matches)[1],4))
    trends[,1:2]=sign(other_temp[sea,matches[,1],1:2,2])
    trends[,3:4]=sign(other_rain[sea,matches[,2],1:2,2])

    sigs<-array(0,c(dim(matches)[1],4))
    sigs[,1:2]=other_temp[sea,matches[,1],1:2,3]
    sigs[,3:4]=other_rain[sea,matches[,2],1:2,3]

    simultan<-which(trends[,1]>0 & trends[,4]>0 & (trends[,2]<0 | trends[,3]<0))
    simultaneous[1,simultan]=trends[simultan,1]*2

    simultan<-which(trends[,2]>0 & trends[,3]>0 & (trends[,1]<0 | trends[,4]<0))
    simultaneous[1,simultan]=trends[simultan,2]*3

    simultan<-which(trends[,1]>0 & trends[,4]>0 & trends[,2]>0 & trends[,3]>0)
    simultaneous[1,simultan]=trends[simultan,1]*1


    topo_map_plot(filename_plot=paste(filename_plot,"_qr95.pdf",sep=""),reihen=array(simultaneous,c(1,ntot)),farb_mitte=c(0,3),farb_palette="wave-goggle",ID_select=1:dim(matches)[1])
    #draw.circle(c(-4,10,24,38),c(46,49,52,55),c(6,6,6,6),lwd=5,lty=2)
    draw.circle(c(-3.75,10,23.75,37.5,51.25,65),c(45,49,53,57,61,65),c(6,6,6,6),lwd=5,lty=2,border="black")
    #lines(c(-17.5,72.5),c(40,60),lty=2,cex=5)
    plot(1:10,col=rgb(1,1,1,0),frame.plot=FALSE,axes=FALSE,ylab="",xlab="")
    legend("topleft",pch=c(15,15,15),col=c(color[2],color[3],color[4]),legend=c("general increase in persistence","increase for floodings","increase for heat-waves"),cex=2,box.lwd=0,box.col=rgb(1,1,1,0))
    graphics.off()
}

#match_grids()

#datPP <<- dat_load_precipitation(paste("../data/_eobsPP/rr_0.50deg_reg_v12.0_1950-2015.nc",sep=""))

source("map_plot.r")
library(grid)
#datPP <<- dat_load_precipitation(paste("../data/_eobsPP/rr_0.50deg_reg_v12.0_1950-2015.nc",sep=""))
#datTG <<- dat_load(filename<-paste("../data/_eobsTG/tg_0.50deg_reg_v12.0_end.nc",sep=""))

plot_init_EU()

nbcol<<-4
closePlot<<-FALSE
color_legend<<-FALSE
wave_goggle(sea=2)



nbcol<<-5
closePlot<<-FALSE
color_legend<<-TRUE

heat_wav_detection(sea=2,stateTemp=2,stateRain=1,farb_palette="norm-hot",filename_plot="/home/peter/Dokumente/pik/backuped/plots/_eobsPP/heatwave_map",stateName=c("dry","warm"))
heat_wav_detection(sea=2,stateTemp=1,stateRain=2,farb_palette="norm-wet",filename_plot="/home/peter/Dokumente/pik/backuped/plots/_eobsPP/flooding_map",stateName=c("wet","cold"))

