

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
    return(c(max(abs(tmp$acf)),tmp$lag[which.max(abs(tmp$acf))]))
}


dissimilarity_matrix <- function(lagMax=3,ID_select=1:1319,timeRange=c(4000,11000),load_name="_Cor",normalize=FALSE){
    #prepare data
    X=array(dat$tas,c(1319,365*65))
    # range in which data coverage is best
    X=X[ID_select,timeRange[1]:timeRange[2]]
    # remove nas problematic

    IDlength=dim(X)[1]
    xLen=length(timeRange[1]:timeRange[2])
    print(proc.time())

    
    if (normalize){
        for (q in ID_select){
            if (length(which(!is.na(X[q,])))>xLen/2){
                # "normalize" using quantiles
                #X[q,]=X[q,]/(quantile(x=X[q,],probs=0.95,na.rm=TRUE)-quantile(x=X[q,],probs=0.05,na.rm=TRUE))
                # "normalize" using sd
                X[q,]=X[q,]/sd(X[q,],na.rm=TRUE)

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
                dissMat[i,j]=1-abs(tmp[1])
                choiceMat[i,j]=tmp[2]
            }
        }
    }

    nc_out<-create.nc(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_dissimilarity_matrix.nc",sep=""))
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

dissimilarity_view <- function(lagMax=15,load_name="_Cor",add_name="",timeRange=c(2000,22000),auswahl=c(333,470,604,890,888)){
    print(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_dissimilarity_matrix.nc",sep=""))
    nc=open.nc(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_dissimilarity_matrix.nc",sep=""))
    distMat<<-var.get.nc(nc,"distanceMat")
    choiceMat<<-var.get.nc(nc,"choiceMat")
    #distMat[is.na(distMat)]=1
    #choiceMat[is.na(choiceMat)]=-15

    reihen=choiceMat[auswahl,]
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,add_name,"_best_lag.pdf",sep=""),reihen=reihen,farb_mitte=0,farb_palette="lila-gruen",highlight_points=auswahl,highlight_color="red",paper=c(8,10),ausschnitt=c(30,70),layout_mat=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6),nrow=16),pointsize=1.5)
    topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_best_lag.pdf",sep=""),reihen=reihen,farb_palette="lila-gruen",highlight_points=auswahl,highlight_color="red",pointsize=1.5)

    reihen=distMat[auswahl,]
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,add_name,"_distance.pdf",sep=""),reihen=reihen,farb_mitte=c(0,1),farb_palette="regenbogen",highlight_points=auswahl,highlight_color="red",paper=c(8,10),ausschnitt=c(30,70),layout_mat=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6),nrow=16),pointsize=1.5)
    topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_distance.pdf",sep=""),reihen=reihen,farb_palette="regenbogen",pointsize=1.5)
}

cluster_evaluation <- function(method="ward.D2",untilGroup=11,add_name="",ID_select=1:1319,load_name="_CorLag",lagMax=20,timeRange=4000:11000,normalize=FALSE){
    print(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_dissimilarity_matrix.nc",sep=""))
    nc=open.nc(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_dissimilarity_matrix.nc",sep=""))
    distMat<<-var.get.nc(nc,"distanceMat")
    distMat[is.na(distMat)]=1

    # X is required for cluster Criteria
    # X is evaluated as in dissimilarity_matrix
    X=array(dat$tas,c(1319,365*65))
    # range in which data coverage is best
    X=X[,timeRange[1]:timeRange[2]]
    xLen=length(timeRange[1]:timeRange[2])

    # to many nas -> flat 
    ID_select_not_flat<-c()
    count<-0
    for (q in ID_select){
            if (length(which(is.na(X[q,])))<xLen/2){
            count<-count+1
            ID_select_not_flat[count]=q
        }
    }
    ID_select_not_flat<<-ID_select_not_flat

    attribution=array(NA,c(untilGroup,1319))
    attribution_changes=array(NA,c(untilGroup,1319))
    criteria=array(NA,c(untilGroup,44))

    if (method!="kmeans"){
        pdf(file=paste("../plots/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_",method,"_",1,"-",untilGroup,"_tree.pdf",sep=""),width=10,height=3)
        par(mar=c(0,5,0,0))
        tmp<<-hclust(as.dist(distMat[ID_select_not_flat,ID_select_not_flat]),method=method)
        for (nGroup in 1:untilGroup){
            print(nGroup)
            same=array(1,nGroup)
            attribution[nGroup,ID_select_not_flat]=cutree(tmp,nGroup)
            plot(tmp,main="",xlab="")
            if (nGroup>1){
            rect.hclust(tmp,k=nGroup)
                for (G in 1:nGroup){
                    GG<-attribution[nGroup-1,ID_select_not_flat][which(attribution[nGroup,ID_select_not_flat]==G)[1]]
                    if (length(which(attribution[nGroup,ID_select_not_flat]==G))==length(which(attribution[nGroup-1,ID_select_not_flat]==GG))){same[G]=0}
                }
                attribution_changes[nGroup,which(attribution[nGroup,ID_select_not_flat] %in% (1:nGroup)[which(same==1)])]=-2
            }
        }
        # clustering height
        tmp_length<-length(tmp$height)
        criteria[1:untilGroup,44]=tmp$height[(tmp_length):(tmp_length-untilGroup+1)]
    }
    graphics.off()

    if (method=="kmeans"){
        attribution[1,ID_select_not_flat]=array(1,length(ID_select_not_flat))
        for (nGroup in 2:untilGroup){
            tmp<<-pam(x=as.dist(distMat[ID_select_not_flat,ID_select_not_flat]),k=nGroup,diss=TRUE)
            attribution[nGroup,ID_select_not_flat]=tmp$clustering
        }
    }

    # within sum of dissimilarities
    wss<-array(0,untilGroup)
    for (k in 1:untilGroup){
        for (G in unique(attribution[k,!is.na(attribution[k,])])){
            inGroup<-which(attribution[k,]==G)
            wss[k]<-wss[k]+sum(distMat[inGroup,inGroup]^2,na.rm=TRUE)/2
        }
    }
    criteria[,43]=wss



    print(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_clustering",add_name,"_",method,"_",1,"-",untilGroup,".nc",sep=""))
    nc_out<-create.nc(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_clustering",add_name,"_",method,"_",1,"-",untilGroup,".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", "gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)
    att.put.nc(nc_out, "NC_GLOBAL", "explanation", "NC_CHAR", "different cuttrees with indices etc")
    att.put.nc(nc_out, "NC_GLOBAL", "method", "NC_CHAR", method)
    
    dim.def.nc(nc_out,"ID",dimlength=1319, unlim=FALSE)   
    dim.def.nc(nc_out,"clusters",dimlength=untilGroup, unlim=FALSE)   
    dim.def.nc(nc_out,"indices",dimlength=44, unlim=FALSE)   

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
    att.put.nc(nc_out, "criteria", "explanation", "NC_CHAR", "intCriteria - all - indices ; wss")

    var.put.nc(nc_out,"attribution",attribution)      
    var.put.nc(nc_out,"attribution_changes",attribution_changes)      
    var.put.nc(nc_out,"criteria",criteria)      
 
    close.nc(nc_out)         
}

cluster_view <- function(lagMax=20,load_name="_CorLag",add_name="",timeRange=c(2000,22000),untilGroup=11,method="ward.D2",ID_select=1:1319){
    print(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_clustering",add_name,"_",method,"_",1,"-",untilGroup,".nc",sep=""))
    nc=open.nc(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_clustering",add_name,"_",method,"_",1,"-",untilGroup,".nc",sep=""))
    attribution<<-var.get.nc(nc,"attribution")
    attribution_changes<<-var.get.nc(nc,"attribution_changes")
    criteria<<-var.get.nc(nc,"criteria")


    topo_map_plot(filename_plot=paste("../plots/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_",method,"_",1,"-",untilGroup,"_map.pdf",sep=""),reihen=attribution[,],farb_palette="viele",pointsize=1.0,yAusschnitt=c(-80,80),paper=c(7,4)) #,reihen_sig=attribution_changes[,]
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_",method,"_",1,"-",untilGroup,"_map.pdf",sep=""),reihen=attribution[,],farb_palette="viele",pointsize=1.5,ausschnitt=c(35,75),paper=c(7,2)) #,reihen_sig=attribution_changes[,]


    if (1==1){
        # ellbow criterium
        pdf(file=paste("../plots/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_",method,"_",1,"-",untilGroup,"_ellbow.pdf",sep=""),width=4,height=4)
        for (eva in 1:untilGroup){
            plot(criteria[,43],xlab="number of groups",ylab="whithin group sum of squared dissimilarities")
            points(eva,criteria[eva,43],cex=2,col="red")
        }
        graphics.off()

        # height criterium
        pdf(file=paste("../plots/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_",method,"_",1,"-",untilGroup,"_height.pdf",sep=""),width=4,height=4)
        for (eva in 1:untilGroup){
            plot(criteria[,44],xlab="number of groups",ylab="clustering height")
            points(eva,criteria[eva,44],cex=2,col="red")
        }
        graphics.off()
    }
}

write_cluster_region_files <- function(lagMax=20,load_name="_CorLag",add_name="",timeRange=c(2000,22000),nGroup=22,untilGroup=25,method="ward.D2",region_name="ward22",ID_select=1:1319){
    print(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_clustering",add_name,"_",method,"_",1,"-",untilGroup,".nc",sep=""))
    nc=open.nc(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_clustering",add_name,"_",method,"_",1,"-",untilGroup,".nc",sep=""))
    attribution<<-var.get.nc(nc,"attribution")

    mids<-array(NA,c(nGroup,3))
    for (G in 1:nGroup){
        inside<<-which(attribution[nGroup,]==G)
        mids[G,]=c(G,mean(dat$lon[inside]),mean(dat$lat[inside]))
    }
    write.table(attribution[nGroup,],paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))
    write.table(mids,paste("../data/",dataset,"/ID_regions/",region_name,"_mids_mean.txt",sep=""))
}

cluster_vis_map <- function(lagMax=20,load_name="_CorLag",add_name="",timeRange=c(2000,22000),nGroup=22,untilGroup=25,method="ward.D2",region_name="ward22",ID_select=1:1319){
    print(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_clustering",add_name,"_",method,"_",1,"-",untilGroup,".nc",sep=""))
    nc=open.nc(paste("../data/",dataset,additional_style,"/clustering/",timeRange[1],"-",timeRange[2],load_name,"_",lagMax,"_clustering",add_name,"_",method,"_",1,"-",untilGroup,".nc",sep=""))
    attribution<<-var.get.nc(nc,"attribution")

    pdf(paste("../plots/",dataset,additional_style,"/clustering/lag_",lagMax,load_name,add_name,"_",method,"_",nGroup,"_vis.pdf",sep=""),width=7,height=5)
    plot(topoWorld,xlim=c(-180,180),ylim=c(-90,90),asp=1.5,location="none",col.land=rgb(0,0,0,0),col.water=rgb(0,0,0,0),mar=c(0,0,0,0))

    tmp=put_points(points=attribution[nGroup,],yAusschnitt=c(-90,90),farb_palette="viele",pointsize=0.93,ID_select=ID_select)
    region_border(region_name=region_name,regNumb=nGroup,border_col="black")
    for (i in c(-60,-30,0,30,60)){
        abline(h=i,lty=2,col="grey")
        text(-170,i,label=i)
    }
    #draw over axes
    polygon(x=c(-200,-200,200,200),y=c(-100,-88,-88,-100),col="white",border="white")
    polygon(x=c(-200,-200,200,200),y=c(100,88,88,100),col="white",border="white")
    graphics.off()
}

init <- function(){
    source("load.r")

    trendID<<-"91_5"
    dataset<<-"_TMean"
    additional_style<<-""
    period<<-"1950-2014"

    source("map_plot.r")
    source("write.r")
    source("functions_regional.r")
    library(fields)

    library(cluster)
    library(clusterCrit)
    library(RNetCDF)
    season_names<<-c("MAM","JJA","SON","DJF","4seasons")

    dat<<-dat_load(paste("../data/_TMean/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
}





init()

load_name<-"_CorSdNorm"

#dissimilarity_matrix(lagMax=20,timeRange=c(2000,22000),load_name="_AbsCorSdNorm",normalize=TRUE)

#dissimilarity_view(lagMax=20,timeRange=c(2000,22000),load_name=load_name)

#for (method in c("ward.D2","single","centroid")){

ID_select=which(dat$lat<60 & dat$lat>30)
for (method in c("ward.D2")){
    print(method)
    #cluster_evaluation(add_name="_ml",load_name=load_name,ID_select=ID_select,timeRange=c(2000,22000),method=method,untilGroup=25)
    #cluster_view(add_name="_ml",load_name=load_name,ID_select=ID_select,timeRange=c(2000,22000),method=method,untilGroup=25)
}

write_cluster_region_files(add_name="_ml",load_name=load_name,ID_select=1:1319,timeRange=c(2000,22000),method="ward.D2",nGroup=7,region_name="ml7")
cluster_vis_map(add_name="_ml",load_name=load_name,ID_select=1:1319,timeRange=c(2000,22000),method="ward.D2",nGroup=7,region_name="ml7")
