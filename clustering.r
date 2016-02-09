
create_regional_distr_out_of_clusters <- function(dataset="_TMean",trendID="91_5",additional_style="",markov_style=5,nGroup=6,add_name="_MarkovDistr",period="1950-2014",region_name="_kmeans"){
    sourceName=paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,".nc",sep="")
    print(sourceName)
    nc=open.nc(sourceName)
    IDregions=var.get.nc(nc,"attribution")
    regional_attribution(dat=dat,region_name=region_name,trendID=trendID,dataset=dataset,additional_style=additional_style,IDregions=IDregions,regNumb=nGroup,comment=sourceName)
}

create_kompakt_map_plot_of_clusters <- function(dataset="_TMean",trendID="91_5",additional_style="",markov_style=5,nGroup=6,add_name="_MarkovDistr",period="1950-2014",region_name="_kmeans"){
    sourceName=paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,".nc",sep="")
    print(sourceName)
    nc=open.nc(sourceName)
    IDregions=var.get.nc(nc,"attribution")[,1:4]
    #IDregions[2:1000,c(1,2,3,4)]=NA
    map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,"_map_kompakt",".pdf",sep=""),worldmap=worldmap,reihen=t(matrix(IDregions,c(1319,4))),farb_palette=c("mixed",nGroup,"groups"),col_row=c(5,4),cex=1,pointsize=0.5,cex_axis=1,paper=c(6.3,4),mat=c(1,1,2,2,1,1,2,2,3,3,4,4,3,3,4,4,5,5,5,5),color_lab="")
}

crossCor <- function(x,y,lagMax){
    tmp<-ccf(x,y,lag.max=lagMax,plot=FALSE,na.action=na.pass)
    return(c(max(tmp$acf),tmp$lag[which.max(tmp$acf)]))
}


dissimilarity_matrix <- function(lagMax=3,ID_select=1:1319,timeRange=4000:11000,add_name="hclust",normalize=FALSE){
    #prepare data
    X=array(dat$tas,c(1319,365*65))
    # range in which data coverage is best
    X=X[ID_select,timeRange]
    # remove nas problematic

    IDlength=dim(X)[1]
    xLen=length(timeRange)
    print(proc.time())

    if (normalize){
        # "normalize" using quantiles
        for (q in ID_select){
            if (length(which(!is.na(X[q,])))>xLen/2){
                X[q,]=X[q,]/(quantile(x=X[q,],probs=0.95,na.rm=TRUE)-quantile(x=X[q,],probs=0.05,na.rm=TRUE))
            }
            else{X[q,]=NA}
        }
        print(proc.time())
    }
    

    # to many nas -> flat 
    X[is.na(X)]=0
    X<<-X
    adasd

    # cross correlation
    dissMat=array(NA,c(IDlength,IDlength))
    choiceMat=array(NA,c(IDlength,IDlength))
    for (i in 1:IDlength){
        print(i)
        print(proc.time())
        for (j in 1:IDlength){
            if (length(which(!is.na(X[i,])))>1000 & length(which(!is.na(X[i,])))>1000){
                tmp=crossCor(X[i,],X[j,],lagMax=lagMax)
                dissMat[i,j]=abs(1-tmp[1])
                choiceMat[i,j]=tmp[2]
            }
        }
    }

    nc_out<-create.nc(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_distance_",lagMax,add_name,".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", "gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)
    att.put.nc(nc_out, "NC_GLOBAL", "explanation", "NC_CHAR", "distance matrix for time series")
    att.put.nc(nc_out, "NC_GLOBAL", "max lags in both directions", "NC_CHAR", paste(lagMax))
    
    dim.def.nc(nc_out,"ID",dimlength=IDlength, unlim=FALSE)   
    dim.def.nc(nc_out,"ID2",dimlength=IDlength, unlim=FALSE)   

    var.def.nc(nc_out,"distanceMat","NC_DOUBLE",c(0,1))
    att.put.nc(nc_out, "distanceMat", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "distanceMat", "dim_explanation", "NC_CHAR", "ID-ID")
    att.put.nc(nc_out, "distanceMat", "explanation", "NC_CHAR", "distences between IDs 1 time lag allowed")

    var.def.nc(nc_out,"choiceMat","NC_DOUBLE",c(0,1))
    att.put.nc(nc_out, "choiceMat", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "choiceMat", "dim_explanation", "NC_CHAR", "ID-ID")
    att.put.nc(nc_out, "choiceMat", "explanation", "NC_CHAR", "which lag is closest")

    var.put.nc(nc_out,"distanceMat",dissMat)      
    var.put.nc(nc_out,"choiceMat",choiceMat)      
 
    close.nc(nc_out)     
}

dissimilarity_view <- function(lagMax=16,nGroup=12,load_name="_Cor",add_name="",method=method){
    print(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_distance_",lagMax,load_name,".nc",sep=""))
    nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_distance_",lagMax,load_name,".nc",sep=""))
    distMat<<-var.get.nc(nc,"distanceMat")
    choiceMat<<-var.get.nc(nc,"choiceMat")
    distMat[is.na(distMat)]=1
    ntot<<-length(ID_select)
    print(proc.time())
    
    #auswahl=c(1250,1171,1148,907,980,692,695,481,479)
    #auswahl=c(333,470,604,890,888)
    #reihen=choiceMat[auswahl,]
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,add_name,"_best_lag.pdf",sep=""),reihen=reihen,farb_mitte=0,farb_palette="lila-gruen",highlight_points=auswahl,highlight_color="red",paper=c(8,10),ausschnitt=c(30,70),layout_mat=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6),nrow=16),pointsize=1.5)
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_groups_",nGroup,"_",method,"_best_lag.pdf",sep=""),reihen=reihen,farb_palette="lila-gruen",highlight_points=auswahl,highlight_color="red")

    #reihen=distMat[auswahl,]
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,add_name,"_distance.pdf",sep=""),reihen=reihen,farb_mitte=c(0,1),farb_palette="regenbogen",highlight_points=auswahl,highlight_color="red",paper=c(8,10),ausschnitt=c(30,70),layout_mat=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6),nrow=16),pointsize=1.5)
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_groups_",nGroup,"_",method,"_distance.pdf",sep=""),reihen=reihen,farb_palette="regenbogen")

    if (method!="kmeans"){
        tmp<<-hclust(as.dist(distMat[ID_select,ID_select]),method=method)
        reihen=array(NA,c(8,1319))
        reihen_sig=array(NA,c(8,1319))
        for (plot in 1:8){
            nGroup=plot+18
            same=array(1,nGroup)
            reihen[plot,ID_select]=cutree(tmp,nGroup)
            if (plot>1){
                for (G in 1:nGroup){
                    GG<<-reihen[plot-1,which(reihen[plot,]==G)[1]]
                    if (length(which(reihen[plot,]==G))==length(which(reihen[plot-1,]==GG))){same[G]=0}
                }
                reihen_sig[plot-1,which(reihen[plot,] %in% (1:nGroup)[which(same==1)])]=-2
            }
        }
    }

    if (method=="kmeans"){
        reihen=array(NA,c(30,1319))
        reihen_sig=array(NA,c(30,1319))
        for (plot in 1:30){
            nGroup=plot
            same=array(1,nGroup)
            reihen[plot,ID_select]=pam(distMat[ID_select,ID_select],nGroup)$clustering
            if (plot>1){
                for (G in 1:nGroup){
                    GG<<-reihen[plot-1,which(reihen[plot,]==G)[1]]
                    if (length(which(reihen[plot,]==G))==length(which(reihen[plot-1,]==GG))){same[G]=0}
                }
                reihen_sig[plot-1,which(reihen[plot,] %in% (1:nGroup)[which(same==1)])]=-2
            }
        }
    }


    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_groups_",nGroup,"_",method,"_map.pdf",sep=""),reihen=reihen,reihen_sig=reihen_sig,farb_palette="viele")
    topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_groups_",nGroup,"_",method,"_map.pdf",sep=""),reihen=reihen,reihen_sig=reihen_sig,pch_style=reihen,farb_palette="viele",paper=c(8,12),ausschnitt=c(30,70),layout_mat=matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9),nrow=17),pointsize=1.5)
}

hcluster_view <- function(lagMax=16,nGroup=12,load_name="_Cor",add_name="",method=method){
    print(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_distance_",lagMax,load_name,".nc",sep=""))
    nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_distance_",lagMax,load_name,".nc",sep=""))
    distMat<<-var.get.nc(nc,"distanceMat")
    choiceMat<<-var.get.nc(nc,"choiceMat")
    distMat[is.na(distMat)]=1
    ntot<<-length(ID_select)
    print(proc.time())
    
    #auswahl=c(1250,1171,1148,907,980,692,695,481,479)
    #auswahl=c(333,470,604,890,888)
    #reihen=choiceMat[auswahl,]
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,add_name,"_best_lag.pdf",sep=""),reihen=reihen,farb_mitte=0,farb_palette="lila-gruen",highlight_points=auswahl,highlight_color="red",paper=c(8,10),ausschnitt=c(30,70),layout_mat=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6),nrow=16),pointsize=1.5)
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_groups_",nGroup,"_",method,"_best_lag.pdf",sep=""),reihen=reihen,farb_palette="lila-gruen",highlight_points=auswahl,highlight_color="red")

    #reihen=distMat[auswahl,]
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,add_name,"_distance.pdf",sep=""),reihen=reihen,farb_mitte=c(0,1),farb_palette="regenbogen",highlight_points=auswahl,highlight_color="red",paper=c(8,10),ausschnitt=c(30,70),layout_mat=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6),nrow=16),pointsize=1.5)
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_groups_",nGroup,"_",method,"_distance.pdf",sep=""),reihen=reihen,farb_palette="regenbogen")

    if (method!="kmeans"){
        tmp<<-hclust(as.dist(distMat[ID_select,ID_select]),method=method)
        reihen=array(NA,c(8,1319))
        reihen_sig=array(NA,c(8,1319))
        for (plot in 1:8){
            nGroup=plot+18
            same=array(1,nGroup)
            reihen[plot,ID_select]=cutree(tmp,nGroup)
            if (plot>1){
                for (G in 1:nGroup){
                    GG<<-reihen[plot-1,which(reihen[plot,]==G)[1]]
                    if (length(which(reihen[plot,]==G))==length(which(reihen[plot-1,]==GG))){same[G]=0}
                }
                reihen_sig[plot-1,which(reihen[plot,] %in% (1:nGroup)[which(same==1)])]=-2
            }
        }
    }

    if (method=="kmeans"){
        reihen=array(NA,c(30,1319))
        reihen_sig=array(NA,c(30,1319))
        for (plot in 1:30){
            nGroup=plot
            same=array(1,nGroup)
            reihen[plot,ID_select]=pam(distMat[ID_select,ID_select],nGroup)$clustering
            if (plot>1){
                for (G in 1:nGroup){
                    GG<<-reihen[plot-1,which(reihen[plot,]==G)[1]]
                    if (length(which(reihen[plot,]==G))==length(which(reihen[plot-1,]==GG))){same[G]=0}
                }
                reihen_sig[plot-1,which(reihen[plot,] %in% (1:nGroup)[which(same==1)])]=-2
            }
        }
    }


    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_groups_",nGroup,"_",method,"_map.pdf",sep=""),reihen=reihen,reihen_sig=reihen_sig,farb_palette="viele")
    topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_groups_",nGroup,"_",method,"_map.pdf",sep=""),reihen=reihen,reihen_sig=reihen_sig,pch_style=reihen,farb_palette="viele",paper=c(8,12),ausschnitt=c(30,70),layout_mat=matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9),nrow=17),pointsize=1.5)
}

init <- function(){
    source("write.r")
    source("load.r")

    trendID<<-"91_5"
    dataset<<-"_TMean"
    additional_style<<-""
    period<<-"1950-2014"

    source("map_plot.r")
    source("functions_regional.r")
    library(rworldmap)
    library(fields)
    library(oce)
    data(topoWorld)
    library(quantreg)
    library(cluster)
    library(clusterCrit)

    worldmap<<-getMap(resolution = "low")
    dat<<-dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
}




ID_select=which(dat$lat>=30 & dat$lat<=70 )
ID_select=1:1319


init()



dissimilarity_matrix(lagMax=15,add_name="_Cor",normalize=FALSE)

dissimilarity_view(lagMax=15,load_name="_Cor",add_name="_ww",method="ward.D2")





