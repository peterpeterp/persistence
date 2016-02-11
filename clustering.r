
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


dissimilarity_matrix <- function(lagMax=3,ID_select=1:1319,timeRange=c(4000,11000),add_name="hclust",normalize=FALSE){
    #prepare data
    X=array(dat$tas,c(1319,365*65))
    # range in which data coverage is best
    X=X[ID_select,timeRange[1]:timeRange[2]]
    # remove nas problematic

    IDlength=dim(X)[1]
    xLen=length(timeRange[1]:timeRange[2])
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

    nc_out<-create.nc(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],add_name,"_",lagMax,"_dissimilarity_matrix.nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", "gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", paste(timeRange[1],"-",timeRange[2]))
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

dissimilarity_view <- function(lagMax=15,load_name="_Cor",add_name="",auswahl=c(333,470,604,890,888)){
    print(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],add_name,"_",lagMax,"_dissimilarity_matrix.nc",sep=""))
    nc=open.nc(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],add_name,"_",lagMax,"_dissimilarity_matrix.nc",sep=""))
    distMat<<-var.get.nc(nc,"distanceMat")
    choiceMat<<-var.get.nc(nc,"choiceMat")
    #distMat[is.na(distMat)]=1
    #choiceMat[is.na(choiceMat)]=-15

    reihen=choiceMat[auswahl,]
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,add_name,"_best_lag.pdf",sep=""),reihen=reihen,farb_mitte=0,farb_palette="lila-gruen",highlight_points=auswahl,highlight_color="red",paper=c(8,10),ausschnitt=c(30,70),layout_mat=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6),nrow=16),pointsize=1.5)
    topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_best_lag.pdf",sep=""),reihen=reihen,farb_palette="lila-gruen",highlight_points=auswahl,highlight_color="red",pointsize=1.5)

    reihen=distMat[auswahl,]
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,add_name,"_distance.pdf",sep=""),reihen=reihen,farb_mitte=c(0,1),farb_palette="regenbogen",highlight_points=auswahl,highlight_color="red",paper=c(8,10),ausschnitt=c(30,70),layout_mat=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6),nrow=16),pointsize=1.5)
    topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_distance.pdf",sep=""),reihen=reihen,farb_palette="regenbogen",pointsize=1.5)
}

cluster_evaluation <- function(method="ward.D2",untilGroup=11,offsetGroup=10,add_name="",ID_select=1:1319,load_name="_CorLag",lagMax=15,timeRange=4000:11000,normalize=FALSE){
    print(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_dissimilarity_matrix.nc",sep=""))
    nc=open.nc(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_dissimilarity_matrix.nc",sep=""))
    distMat<<-var.get.nc(nc,"distanceMat")
    distMat[is.na(distMat)]=1

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
    criteria=array(NA,c(untilGroup,42))

    if (method!="kmeans"){
        tmp<-hclust(as.dist(distMat[ID_select,ID_select]),method=method)
        rm(distMat)
        for (eva in 1:untilGroup){
            nGroup<-eva+offsetGroup
            print(paste(eva,nGroup))
            same=array(1,nGroup)
            attribution[eva,ID_select]=cutree(tmp,nGroup)

            if (eva>1){
                critTmp=intCriteria(traj=X[ID_select,],part=attribution[eva,ID_select],crit="all")
                for (index in 1:42){criteria[eva,index]=critTmp[[index]]}
                for (G in 1:nGroup){
                    GG<-attribution[eva-1,which(attribution[eva,]==G)[1]]
                    if (length(which(attribution[eva,]==G))==length(which(attribution[eva-1,]==GG))){same[G]=0}
                }
                attribution_changes[eva-1,which(attribution[eva,] %in% (1:nGroup)[which(same==1)])]=-2
            }
        }
    }
    nc_out<-create.nc(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_clustering",add_name,"_",method,"_",offsetGroup,"-",untilGroup+offsetGroup,".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", "gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)
    att.put.nc(nc_out, "NC_GLOBAL", "explanation", "NC_CHAR", "different cuttrees with indices etc")
    att.put.nc(nc_out, "NC_GLOBAL", "method", "NC_CHAR", method)
    
    dim.def.nc(nc_out,"ID",dimlength=1319, unlim=FALSE)   
    dim.def.nc(nc_out,"clusters",dimlength=untilGroup, unlim=FALSE)   
    dim.def.nc(nc_out,"indices",dimlength=42, unlim=FALSE)   

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

    var.put.nc(nc_out,"attribution",attribution)      
    var.put.nc(nc_out,"attribution_changes",attribution_changes)      
    var.put.nc(nc_out,"criteria",criteria)      
 
    close.nc(nc_out)         
}

cluster_view <- function(lagMax=15,load_name="_Cor",add_name="",untilGroup=11,offsetGroup=10,method="ward.D2"){
    print(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_clustering",add_name,"_",method,"_",offsetGroup,"-",untilGroup+offsetGroup,".nc",sep=""))
    nc=open.nc(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_clustering",add_name,"_",method,"_",offsetGroup,"-",untilGroup+offsetGroup,".nc",sep=""))
    attribution<-var.get.nc(nc,"attribution")
    attribution_changes<-var.get.nc(nc,"attribution_changes")
    criteria<-var.get.nc(nc,"criteria")


    topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_",method,"_",offsetGroup,"-",untilGroup+offsetGroup,"_map.pdf",sep=""),reihen=attribution,reihen_sig=attribution_changes,farb_palette="viele",pointsize=1.5)

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_",method,"_",offsetGroup,"-",untilGroup+offsetGroup,"_criteria.pdf",sep=""))
    plot((1:untilGroup)+offsetGroup,criteria[,1])

}

init <- function(){
    source("load.r")

    trendID<<-"91_5"
    dataset<<-"_TMean"
    additional_style<<-""
    period<<-"1950-2014"

    source("map_plot.r")
    library(fields)
    library(oce)
    data(topoWorld)
    library(cluster)
    library(clusterCrit)

    dat<<-dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
}








init()


dissimilarity_matrix(lagMax=15,timeRange=c(2000,22000),add_name="_CorLag",normalize=FALSE)

#dissimilarity_view(lagMax=15,load_name="_Cor",add_name="_ww")

cluster_evaluation(add_name="_ww",ID_select=1:1319,timeRange=c(2000,22000),method="ward.D2",untilGroup=10,offsetGroup=15)

#cluster_view(add_name="_ww",method="ward.D2",untilGroup=10,offsetGroup=15)



