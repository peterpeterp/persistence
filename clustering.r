

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



clustering <- function(add_name="_forReal_",seasons=1:5,nGroup=7,versions=30,runsMax=100,ID_select=1:1319,distr_select=1:100,plot=c("testMasseGroups","testMasseMaps","endGroups","endMaps")){

    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    state_names<<-c("cold","warm")

    jet.colors <- colorRampPalette(c("black","blue","green","yellow","orange","red",rgb(0.5,0.5,0.5,0.5)))
    color <<- jet.colors(nGroup) 

    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_distributions.nc",sep=""))
    distr_stuff=var.get.nc(nc,"distr_stuff")


    ntot<<-1319
    dimDistr<<-length(distr_select)

    regionAttribution=array(NA,dim=c(6,ntot))
    regionAttribution[6,]=1:ntot
    GroupDistributions=array(NA,dim=c(6,nGroup,dimDistr*2))

    for (sea in seasons){
        if (1==1){
            state<<-state
            print(paste(season_names[sea]))

            toOrder=array(NA,c(ntot,dimDistr*2))
            toOrder[ID_select,1:dimDistr]=distr_stuff[sea,ID_select,1,2,distr_select]
            toOrder[ID_select,(dimDistr+1):(dimDistr*2)]=distr_stuff[sea,ID_select,2,2,distr_select]


            toOrder[is.na(toOrder)]=0
            toOrder<<-toOrder

            tmp<<-kmeans(x=toOrder,centers=nGroup,iter.max=100,nstart=20)

            attribution_unordered=tmp$cluster
            groups_unordered=tmp$centers

            attribution=attribution_unordered*NA
            groups=groups_unordered*NA
            for (G in 1:nGroup){
                Gmax=which.max(groups_unordered[,21])
                attribution[attribution_unordered==Gmax]=G
                groups[G,]=groups_unordered[Gmax,]
                groups_unordered[Gmax,]=-999
            }

            if ("endGroups" %in% plot){
                pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_",season_names[sea],"_",add_name,"_",nGroup,".pdf",sep=""))
                plot_aktuelles_muster(attribution=attribution,start=groups,nGroup=nGroup,points=TRUE)
            graphics.off()
            }

            regionAttribution[sea,]=attribution
            GroupDistributions[sea,,]=groups
        }
    }
    if ("endMaps" %in% plot){
        topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,"_map",".pdf",sep=""),reihen=regionAttribution,farb_palette="spacy")
    }


    nc_out<-create.nc(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", "gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)
    att.put.nc(nc_out, "NC_GLOBAL", "over all distances from points to group mean", "NC_CHAR", paste(group_structure_evaluation(attribution=attribution,groups=groups,nGroup=nGroup)))
    
    dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
    dim.def.nc(nc_out,"seasons",dimlength=6,unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
    
    dim.def.nc(nc_out,"nGroup",dimlength=nGroup,unlim=FALSE)

    dim.def.nc(nc_out,"dimDistr",dimlength=dimDistr,unlim=FALSE)

    var.def.nc(nc_out,"attribution","NC_DOUBLE",c(1,2,0))
    att.put.nc(nc_out, "attribution", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "attribution", "dim_explanation", "NC_CHAR", "ID-season")
    att.put.nc(nc_out, "attribution", "explanation", "NC_CHAR", "group attribution for eachs season")

    var.def.nc(nc_out,"Distribution","NC_DOUBLE",c(1,2,3,5))
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
    #X<<-X

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

hcluster_view <- function(lagMax=16,nGroup=12,add_name="hclust"){
    print(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_distance_",lagMax,add_name,".nc",sep=""))
    nc=open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/clustering/",period,"/",trendID,"_",period,"_distance_",lagMax,add_name,".nc",sep=""))
    distMat<<-var.get.nc(nc,"distanceMat")
    choiceMat<<-var.get.nc(nc,"choiceMat")
    tmp<<-hclust(as.dist(distMat),method="ward.D")
    print(proc.time())
    
    auswahl=c(553,890,132,604,1008)
    reihen=choiceMat[auswahl,]
    #print(dim(reihen))
    #map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,"_best_laf.pdf",sep=""),worldmap=worldmap,reihen=reihen,farb_mitte=0,farb_palette="lila-gruen",paper=c(8,6),highlight_points=auswahl,highlight_color="orange")
    #reihen=distMat[auswahl,]
    #topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,"_distance.pdf",sep=""),reihen=reihen,farb_mitte=c(0,0.006))

    reihen=array(NA,c(20,1319))
    for (nGroup in 1:20){reihen[nGroup,]=cutree(tmp,nGroup)}
    topo_map_plot(filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/clustering/lag_",lagMax,"_groups_",nGroup,"_hcluster_map.pdf",sep=""),reihen=reihen,farb_palette="spacy")

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

    worldmap<<-getMap(resolution = "low")
    dat<<-dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
}

kmeans_master <- function(nGroup=7,add_name="default"){
    clustering(add_name=add_name,seasons=1:5,nGroup=nGroup,versions=10,runsMax=100,ID_select=ID_select,plot=c("endMaps"))
    #create_regional_distr_out_of_kmeans(add_name=add_name,region_name=paste(add_name,"_",nGroup,sep=""),nGroup=nGroup)
}


#create_kompakt_map_plot_of_kmeans(nGroup=7,add_name="KarlPerason_AmpMark")

ID_select=which(dat$lat>=10 & dat$lat<=90 )
#ID_select=1:1319

#kmeans_raw(nGroup=7,add_name=paste("raw___",tries,sep=""),starts=10,runsMax=100,ID_select=ID_select,seasons=c(5))

init()
kmeans_master(nGroup=6,add_name=paste("2distr","neu",sep=""))



#distance_matrix(lagMax=15,add_name="_CorNorm",normalize=TRUE)
#distance_matrix(lagMax=15,add_name="_Cor",normalize=FALSE)

#hcluster_view(lagMax=10,add_name="_Cor")
#hcluster_view(lagMax=4,add_name="_quNorm")
#hcluster_view(lagMax=2,add_name="_quNorm")
#hcluster_view(lagMax=1,add_name="_quNorm")
#hcluster_view(lagMax=0,add_name="_quNorm")





