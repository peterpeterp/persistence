

group_reduction <- function(attribution,start,nGroup,reduce,main_add=""){
    score=array(999,c(nGroup,nGroup))
    for (m in 1:reduce){
        for (p in 1:nGroup){
            for (p2 in 1:nGroup){
                if (p!=p2){
                    #score[p]=sum((start[p,]-toOrder[q,])^2)
                    if (length(which(attribution==p))+length(which(attribution==p2))<1000){score[p,p2]=sum(abs(start[p,]-start[p2,]),na.rm=TRUE)}
                    else {score[p,p2]=999}
                }
            }
        }
        similar=which(score==score[which.min(score)],arr.ind=TRUE)

        start[similar[1],]=colMeans(start[similar[1:2],],na.rm=TRUE)
        start[similar[2],]=NA
        #plot_aktuelles_muster(attribution,start,nGroup)
        attribution[which(attribution==similar[2])]=similar[1]
        start=start[which((1:nGroup)!=similar[2]),]
    }
    if (length(which(attribution==nGroup))!=0){
        empty=which(!((1:nGroup) %in% attribution))
        print(empty)
        print(length(which(attribution==nGroup)))
        attribution[attribution==nGroup]=empty
    }
    return(list(attribution=attribution,start=start))
}

k_clustering <- function(nGroup=7,start_mod="random",runs=30){
    Empty=which(is.na(toOrder[,1]))

    if (start_mod[1]=="random"){attribution=sample(nGroup,ntot,replace=TRUE)}
    if (start_mod[1]!="random"){attribution=start_mod}
    attribution[Empty]=0

    start=array(NA,dim=c(nGroup,dimensionality))
    for (G in 1:nGroup){
        start[G,]=colMeans(toOrder[which(attribution==G),],na.rm=TRUE)
    }
    for (i in 1:runs){
        cat(paste(i," "))
        abbruchCond=0
        for (q in 1:ntot){
            if (!is.na(toOrder[q,1])){
                score=array(NA,nGroup)
                for (G in 1:nGroup){
                    # distance definitions here
                    score[G]=sum(abs(start[G,]-toOrder[q,]),na.rm=TRUE)
                    #score[G]=sum((weight*(start[G,]-toOrder[q,])^2),na.rm=TRUE)
                }
                G=which.min(score)
                G_old=attribution[q]
                if (G_old!=G){
                    attribution[q]=G
                    if (length(which(attribution==G))>1){start[G,]=colMeans(toOrder[which(attribution==G),],na.rm=TRUE)}
                    if (length(which(attribution==G))==1){start[G,]=toOrder[which(attribution==G),]}
                    if (length(which(attribution==G_old))>1){start[G_old,]=colMeans(toOrder[which(attribution==G_old),],na.rm=TRUE)}
                    if (length(which(attribution==G_old))==1){start[G_old,]=toOrder[which(attribution==G_old),]}
                    abbruchCond=1
                }  
            }
        }
        if (abbruchCond==0){break}
    }

    attribution[Empty]=0
    return(list(attribution=attribution,groups=start))    
}



group_structure_evaluation <- function(attribution,groups,nGroup){
    score=0
    for (G in nGroup){
        for (q in which(attribution==G)){
            score=score+sum(abs(groups[G,]-toOrder[q,]),na.rm=TRUE)
        }
    }
    return(score)
}

plot_aktuelles_muster <- function(attribution,start,nGroup,points=FALSE,main_zusatz=""){
    if (dimDistr>1){
        for (logAx in c("","y")){
            plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,log=logAx,main=paste(main_zusatz,"Distribution of",state_names[state],"periods"))
            for (p in 1:nGroup){
                lines(1:dimDistr,start[p,],col=color[p])
            }
        }
    }

    if (points==TRUE){
        for (p in 1:nGroup){
            if (!is.na(start[p,1])){
                for (logAx in c("","y")){
                    plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,log=logAx,main=paste(main_zusatz,"Distribution of",state_names[state],"periods","group:",p,"members:",length(which(attribution==p)))) 
                    for (q in which(attribution==p)){
                         points(1:dimDistr,toOrder[q,],pch=16,col=rgb(0.5,0.5,0.5,0.05),cex=1.5)
                    }
                    lines(1:dimDistr,start[p,],col=color[p])
                }
            }
        }
    }
}



clustering <- function(add_name="_forReal_",seasons=1:5,nGroupS=7:10,versions=30,runsMax=100,ID_select=1:1319,distr_select=1:100,plot=c("testMasseGroups","testMasseMaps","endGroups","endMaps"),method="kmeans"){

    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    state_names<<-c("cold","warm")

    jet.colors <- colorRampPalette(c("black","blue","green","yellow","orange","red",rgb(0.5,0.5,0.5,0.5)))
    color <<- jet.colors(tail(nGroupS,n=1)) 

    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_distributions.nc",sep=""))
    distr_stuff=var.get.nc(nc,"distr_stuff")

    results<<-length(nGroupS)
    ntot<<-1319
    dimDistr<<-length(distr_select)

    regionAttribution=array(NA,dim=c(5,results,ntot))
    GroupDistributions=array(NA,dim=c(5,results,tail(nGroupS,n=1),dimDistr*2))

    for (sea in seasons){
        if (1==1){
            print(paste(season_names[sea]))

            toOrder=array(NA,c(ntot,dimDistr*2))
            toOrder[ID_select,1:dimDistr]=distr_stuff[sea,ID_select,1,2,distr_select]
            toOrder[ID_select,(dimDistr+1):(dimDistr*2)]=distr_stuff[sea,ID_select,2,2,distr_select]

            #re normalize in order to get Karl-Pearson
            for (i in 1:(dimDistr*2)){
                toOrder[,i]=toOrder[,i]/sd(toOrder[,i],na.rm=TRUE)
            }


            toOrder[is.na(toOrder)]=0
            toOrder<<-toOrder

            if (method=="kmeans"){
                for (i in 1:results){
                    nGroup<-nGroupS[i]
                    tmp<<-kmeans(x=toOrder,centers=nGroup,iter.max=100,nstart=20)
                    regionAttribution[sea,i,]=tmp$cluster
                    GroupDistributions[sea,i,1:nGroup,]=tmp$centers
                }
            }

            if (method!="kmeans"){
                tmp<<-agnes(x=toOrder,method=method)
                for (i in 1:results){
                    nGroup<-nGroupS[i]
                    regionAttribution[sea,i,]=cutree(tmp,nGroup)
                    groups<-array(NA,c(nGroup,2*dimDistr))
                    for (G in 1:nGroup){
                        groups[G,]=colMeans(toOrder[which(regionAttribution[sea,i,]==G),],na.rm=TRUE)
                    }
                    groups<<-groups
                    GroupDistributions<<-GroupDistributions
                    GroupDistributions[sea,i,1:nGroup,]=groups
                }
                
            }
        }
    }
    if ("endMaps" %in% plot){
        reihen=array(NA,c(length(seasons)*results,ntot))
        titel=c()
        for (sea in seasons){
            for (i in 1:results){
                index<-(sea-1)*results+i
                titel[index]=paste(season_names[sea],"groups:",nGroupS[i])
                reihen[index,]=regionAttribution[sea,i,]
            }
        }
        topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_",add_name,"_",method,"_",tail(nGroupS,n=1),"_map",".pdf",sep=""),reihen=reihen,farb_palette="viele",paper=c(8,10),ausschnitt=c(30,70),layout_mat=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6),nrow=16),pointsize=1.5)
    }


    nc_out<-create.nc(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_",add_name,"_",method,"_",tail(nGroupS,n=1),".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", "gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)
    
    dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
    dim.def.nc(nc_out,"seasons",dimlength=5,unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
    
    dim.def.nc(nc_out,"nResults",dimlength=results,unlim=FALSE)
    dim.def.nc(nc_out,"nGroup",dimlength=tail(nGroupS,n=1),unlim=FALSE)

    dim.def.nc(nc_out,"dimDistr",dimlength=dimDistr,unlim=FALSE)

    var.def.nc(nc_out,"attribution","NC_DOUBLE",c(1,2,0))
    att.put.nc(nc_out, "attribution", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "attribution", "dim_explanation", "NC_CHAR", "ID-season")
    att.put.nc(nc_out, "attribution", "explanation", "NC_CHAR", "group attribution for eachs season")

    var.def.nc(nc_out,"Distribution","NC_DOUBLE",c(1,2,3,4,5))
    att.put.nc(nc_out, "Distribution", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "Distribution", "dim_explanation", "NC_CHAR", "season,states,groups,distributions")
    att.put.nc(nc_out, "Distribution", "explanation", "NC_CHAR", "group characteristics Distribution")
       
    var.put.nc(nc_out,"attribution",regionAttribution)      
    var.put.nc(nc_out,"Distribution",GroupDistributions)      
 
    close.nc(nc_out) 
}

create_regional_distr_out_of_kmeans <- function(dataset="_TMean",trendID="91_5",additional_style="",markov_style=5,nGroup=6,add_name="_MarkovDistr",period="1950-2014",region_name="_kmeans"){
    sourceName=paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,".nc",sep="")
    print(sourceName)
    nc=open.nc(sourceName)
    IDregions=var.get.nc(nc,"attribution")
    regional_attribution(dat=dat,region_name=region_name,trendID=trendID,dataset=dataset,additional_style=additional_style,IDregions=IDregions,regNumb=nGroup,comment=sourceName)
}

create_kompakt_map_plot_of_kmeans <- function(dataset="_TMean",trendID="91_5",additional_style="",markov_style=5,nGroup=6,add_name="_MarkovDistr",period="1950-2014",region_name="_kmeans"){
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


distance_matrix <- function(lagMax=3,ID_select=1:1319,timeRange=4000:11000,add_name="hclust",normalize=FALSE){
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
    distMat=array(NA,c(IDlength,IDlength))
    choiceMat=array(NA,c(IDlength,IDlength))
    for (i in 1:IDlength){
        print(i)
        print(proc.time())
        for (j in 1:IDlength){
            if (length(which(!is.na(X[i,])))>1000 & length(which(!is.na(X[i,])))>1000){
                tmp=crossCor(X[i,],X[j,],lagMax=lagMax)
                distMat[i,j]=abs(1-tmp[1])
                choiceMat[i,j]=tmp[2]
            }
        }
    }

    nc_out<-create.nc(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_distance_",lagMax,add_name,".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", "gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)
    att.put.nc(nc_out, "NC_GLOBAL", "explanation", "NC_CHAR", "distance matrix for time series")
    att.put.nc(nc_out, "NC_GLOBAL", "max lags in both directions", "NC_CHAR", "3")
    
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

    var.put.nc(nc_out,"distanceMat",distMat)      
    var.put.nc(nc_out,"choiceMat",choiceMat)      
 
    close.nc(nc_out)     
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

kmeans_master <- function(nGroupS=7:10,add_name="default",method="kmeans"){
    clustering(add_name=add_name,seasons=1:5,nGroupS=nGroupS,versions=10,runsMax=100,ID_select=ID_select,plot=c("endMaps"),method=method)
    #create_regional_distr_out_of_kmeans(add_name=add_name,region_name=paste(add_name,"_",nGroup,sep=""),nGroup=nGroup)
}


#create_kompakt_map_plot_of_kmeans(nGroup=7,add_name="KarlPerason_AmpMark")

ID_select=which(dat$lat>=30 & dat$lat<=70 )
ID_select=1:1319

#kmeans_raw(nGroup=7,add_name=paste("raw___",tries,sep=""),starts=10,runsMax=100,ID_select=ID_select,seasons=c(5))

#init()
#kmeans_master(nGroupS=17,add_name="2distr_",method="ward",normalize=TRUE)



#distance_matrix(lagMax=15,add_name="_CorNorm",normalize=TRUE)
distance_matrix(lagMax=15,add_name="_Cor",normalize=FALSE)

#hcluster_view(lagMax=15,add_name="_CorNorm")
#hcluster_view(lagMax=15,load_name="_Cor",add_name="2",method="ward.D2")
#hcluster_view(lagMax=15,load_name="_Cor",add_name="_ww",method="ward.D2")
#hcluster_view(lagMax=4,add_name="_quNorm")
#hcluster_view(lagMax=2,add_name="_quNorm")
#hcluster_view(lagMax=1,add_name="_quNorm")
#hcluster_view(lagMax=0,add_name="_quNorm")





