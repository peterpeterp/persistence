
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
    #X<<-X

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

dissimilarity_view <- function(lagMax=16,nGroup=12,load_name="_Cor",add_name="",method=method,auswahl=c(333,470,604,890,888)){
    print(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_distance_",lagMax,load_name,".nc",sep=""))
    nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_distance_",lagMax,load_name,".nc",sep=""))
    distMat<<-var.get.nc(nc,"distanceMat")
    choiceMat<<-var.get.nc(nc,"choiceMat")
    distMat[is.na(distMat)]=1
    ntot<<-length(ID_select)

    reihen=choiceMat[auswahl,]
    topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,add_name,"_best_lag.pdf",sep=""),reihen=reihen,farb_mitte=0,farb_palette="lila-gruen",highlight_points=auswahl,highlight_color="red",paper=c(8,10),ausschnitt=c(30,70),layout_mat=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6),nrow=16),pointsize=1.5)
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_groups_",nGroup,"_",method,"_best_lag.pdf",sep=""),reihen=reihen,farb_palette="lila-gruen",highlight_points=auswahl,highlight_color="red")

    reihen=distMat[auswahl,]
    topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,add_name,"_distance.pdf",sep=""),reihen=reihen,farb_mitte=c(0,1),farb_palette="regenbogen",highlight_points=auswahl,highlight_color="red",paper=c(8,10),ausschnitt=c(30,70),layout_mat=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6),nrow=16),pointsize=1.5)
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_groups_",nGroup,"_",method,"_distance.pdf",sep=""),reihen=reihen,farb_palette="regenbogen")

}

cluster_evaluation <- function(method="ward.D2",untilGroup=30,add_name="",ID_select=1:1319,load_name="_Cor",lagMax=15,timeRange=4000:11000,normalize=FALSE){
    print(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_distance_",lagMax,load_name,".nc",sep=""))
    nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_distance_",lagMax,load_name,".nc",sep=""))
    distMat<<-var.get.nc(nc,"distanceMat")
    distMat[is.na(distMat)]=1
    ntot<<-length(ID_select)

    # X is required for cluster Criteria
    # X is evaluated as in dissimilarity_matrix
    X=array(dat$tas,c(1319,365*65))
    # range in which data coverage is best
    X=X[ID_select,timeRange]
    # remove nas problematic
    xLen=length(timeRange)

    if (normalize){
        # "normalize" using quantiles
        for (q in ID_select){
            if (length(which(!is.na(X[q,])))>xLen/2){
                X[q,]=X[q,]/(quantile(x=X[q,],probs=0.95,na.rm=TRUE)-quantile(x=X[q,],probs=0.05,na.rm=TRUE))
            }
            else{X[q,]=NA}
        }
    }
    # to many nas -> flat 
    X[is.na(X)]=0
    X<<-X

    attribution=array(NA,c(untilGroup,1319))
    attribution_changes=array(NA,c(untilGroup,1319))
    criteria=array(NA,c(untilGroup,100))

    if (method!="kmeans"){
        tmp<<-hclust(as.dist(distMat[ID_select,ID_select]),method=method)
        for (nGroup in 1:untilGroup){
            same=array(1,nGroup)
            attribution[nGroup,ID_select]=cutree(tmp,nGroup)
            critTmp=intCriteria(X[ID_select,],attribution[nGroup,ID_select],"all")
            criteria[nGroup,1:length(critTmp)]=critTmp
            if (nGroup>1){
                for (G in 1:nGroup){
                    GG<<-attribution[nGroup-1,which(attribution[nGroup,]==G)[1]]
                    if (length(which(attribution[nGroup,]==G))==length(which(attribution[nGroup-1,]==GG))){same[G]=0}
                }
                attribution_changes[nGroup-1,which(attribution[nGroup,] %in% (1:nGroup)[which(same==1)])]=-2
            }
        }
    }

    if (method=="kmeans"){
        for (nGroup in 1:30){
            same=array(1,nGroup)
            attribution[nGroup,ID_select]=pam(distMat[ID_select,ID_select],nGroup)$clustering
            intCriteria(X[ID_select,],attribution[nGroup,ID_select],"all")
            if (nGroup>1){
                for (G in 1:nGroup){
                    GG<<-attribution[nGroup-1,which(attribution[nGroup,]==G)[1]]
                    if (length(which(attribution[nGroup,]==G))==length(which(attribution[nGroup-1,]==GG))){same[G]=0}
                }
                attribution_changes[nGroup-1,which(attribution[nGroup,] %in% (1:nGroup)[which(same==1)])]=-2
            }
        }
    }

    nc_out<-create.nc(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_distance_",lagMax,add_name,"_",method,"_",untilGroup,".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", "gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)
    att.put.nc(nc_out, "NC_GLOBAL", "explanation", "NC_CHAR", "different cuttrees with indices etc")
    
    dim.def.nc(nc_out,"ID",dimlength=1319, unlim=FALSE)   
    dim.def.nc(nc_out,"clusters",dimlength=untilGroup, unlim=FALSE)   
    dim.def.nc(nc_out,"indices",dimlength=     , unlim=FALSE)   

    var.def.nc(nc_out,"attribution","NC_DOUBLE",c(1,0))
    att.put.nc(nc_out, "attribution", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "attribution", "dim_explanation", "NC_CHAR", "nClusters-ID")
    att.put.nc(nc_out, "attribution", "explanation", "NC_CHAR", "attributions")

    var.def.nc(nc_out,"attribution_changes","NC_DOUBLE",c(1,0))
    att.put.nc(nc_out, "attribution_changes", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "attribution_changes", "dim_explanation", "NC_CHAR", "nClusters-ID")
    att.put.nc(nc_out, "attribution_changes", "explanation", "NC_CHAR", "changes with respect to previous number of clusters")

    var.def.nc(nc_out,"criteria","NC_DOUBLE",c(1,2))
    att.put.nc(nc_out, "criteria", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "criteria", "dim_explanation", "NC_CHAR", "nClusters-indices")
    att.put.nc(nc_out, "criteria", "explanation", "NC_CHAR", "intCriteria - all - indices")

    var.put.nc(nc_out,"attribution",dissMat)      
    var.put.nc(nc_out,"attribution_changes",dissMat)      
    var.put.nc(nc_out,"criteria",choiceMat)      
 
    close.nc(nc_out)         
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








init()



#dissimilarity_matrix(lagMax=15,add_name="_Cor",normalize=FALSE)

#dissimilarity_view(lagMax=15,load_name="_Cor",add_name="_ww",method="ward.D2")

cluster_evaluation(add_name="_ww",ID_select=1:1319)
#cluster_evaluation(add_name="_ml",ID_select=which(dat$lat>=30 & dat$lat<=70 ))



