

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


plot_aktuelles_muster <- function(attribution,start,nGroup,points=FALSE,main_zusatz=""){
    for (state in 1:2){
        for (logAx in c("","y")){
            plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,log=logAx,main=paste(main_zusatz,"Distribution of",state_names[state],"periods"))
            for (p in 1:nGroup){
                lines(1:dimDistr,start[p,startDistr[state]:stopDistr[state]],col=color[p])
            }
        }
    }
    if (dimMarkov>1){
        plot(NA,xlab="days",ylab="probability density",ylim=c(0,1),xlim=c(0,dimMarkov+1),axes=TRUE,frame.plot=TRUE,log="",main=paste(main_zusatz,"MarkovChain"))
        for (p in 1:nGroup){
            lines(1:dimMarkov,start[p,1:dimMarkov],col=color[p])
        }
    }
    if (points==TRUE){
        for (state in 1:2){
            for (p in 1:nGroup){
                if (!is.na(start[p,1])){
                    for (logAx in c("","y")){
                        plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,log=logAx,main=paste(main_zusatz,"Distribution of",state_names[state],"periods","group:",p,"members:",length(which(attribution==p)))) 
                        for (q in which(attribution==p)){
                            points(1:dimDistr,toOrder[q,startDistr[state]:stopDistr[state]],pch=16,col=rgb(0.5,0.5,0.5,0.05),cex=1.5)
                        }
                        lines(1:dimDistr,start[p,startDistr[state]:stopDistr[state]],col=color[p])
                    }
                }
            }
        }
        if (dimMarkov>1){
            for (p in 1:nGroup){
                if (!is.na(start[p,1])){
                    plot(NA,xlab="days",ylab="probability density",ylim=c(0,1),xlim=c(0,dimMarkov+1),axes=TRUE,frame.plot=TRUE,log="",main=paste(main_zusatz,"MarkovChain","group:",p,"members:",length(which(attribution==p)))) 
                    for (q in which(attribution==p)){
                        points(1:dimMarkov,toOrder[q,1:dimMarkov],pch=16,col=rgb(0.5,0.5,0.5,0.05),cex=1.5)
                    }
                    lines(1:dimMarkov,start[p,1:dimMarkov],col=color[p])
                }
            }
        }
    }

}

k_nearest_neighbors <- function(nGroup=7,start_mod="random",runs=30){
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
                    #score[G]=sum(abs(start[G,]-toOrder[q,]),na.rm=TRUE)
                    score[G]=sum((weight*(start[G,]-toOrder[q,])^2),na.rm=TRUE)
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

nearest_neighbors <- function(markov_style=NA,add_name="_forReal_",seasons=1:5,nGroupEnd=7,versions=30,runs=30,plot=c("testMasseGroups","testMasseMaps","endGroups","endMaps")){

    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    state_names<<-c("cold","warm")

    jet.colors <- colorRampPalette(c("black","blue","green","yellow","orange","red",rgb(0.5,0.5,0.5,0.5)))
    color <<- jet.colors(nGroupEnd) 

    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_distributions.nc",sep=""))
    distr_stuff=var.get.nc(nc,"distr_stuff")

    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/markov/",period,"/",trendID,dataset,"_",period,"_markov_order",markov_style,".nc",sep=""))
    eventResult=var.get.nc(nc,"eventResult")

    ntot<<-dim(distr_stuff)[2]
    dimMarkov<<-dim(eventResult)[4]
    dimMarkov<<-0
    dimDistr<<-dim(distr_stuff)[5]
    dimensionality<<-dimMarkov+2*dimDistr

    startDistr<<-c(dimMarkov+1,dimMarkov+dimDistr+1)
    stopDistr<<-c(dimMarkov+dimDistr,dimMarkov+dimDistr+dimDistr)

    regionAttribution=array(NA,dim=c(ntot,6))
    regionAttribution[,6]=1:ntot
    GroupMarkovChain=array(NA,dim=c(6,nGroupEnd,dimMarkov))
    GroupDistributions=array(NA,dim=c(6,2,nGroupEnd,dimDistr))

    for (sea in seasons){

        toOrder=array(NA,dim=c(ntot,dimensionality))
        #toOrder[,1:dimMarkov]=eventResult[sea,,4,]
        toOrder[,startDistr[1]:stopDistr[1]]=distr_stuff[sea,,1,2,]
        toOrder[,startDistr[2]:stopDistr[2]]=distr_stuff[sea,,2,2,]
        
        toOrder[is.na(toOrder)]=0
        toOrder<<-toOrder

        tmp<<-kmeans(x=toOrder,centers=nGroupEnd,iter.max=100,nstart=100)
        attribution=tmp$cluster
        groups=tmp$centers

        if ("endGroups" %in% plot){
            pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbors/",period,"/",trendID,"_",period,"_",season_names[sea],"_",add_name,"_",nGroupEnd,".pdf",sep=""))
            plot_aktuelles_muster(attribution=attribution,start=groups,nGroup=nGroupEnd,points=TRUE)
        graphics.off()
        }

        regionAttribution[,sea]=attribution
        GroupMarkovChain[sea,,]=groups[,1:dimMarkov]
        GroupDistributions[sea,1,,]=groups[,startDistr[1]:stopDistr[1]]
        GroupDistributions[sea,2,,]=groups[,startDistr[2]:stopDistr[2]]
    }
    if ("endMaps" %in% plot){map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbors/",period,"/",trendID,"_",period,"_",add_name,"_",nGroupEnd,"_map",".pdf",sep=""),worldmap=worldmap,reihen=t(matrix(regionAttribution,c(ntot,6))),pointsize=1.5,farb_palette=c("mixed",nGroupEnd,"groups"))}


    nc_out<-create.nc(paste("../data/",trendID,"/",dataset,additional_style,"/nearest_neighbors/",period,"/",trendID,"_",period,"_",add_name,"_",nGroupEnd,".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", "gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)
    att.put.nc(nc_out, "NC_GLOBAL", "over all distances from points to group mean", "NC_CHAR", paste(group_structure_evaluation(attribution=attribution,groups=groups,nGroup=nGroupEnd)))
    
    dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
    dim.def.nc(nc_out,"seasons",dimlength=6,unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
    
    dim.def.nc(nc_out,"nGroup",dimlength=nGroupEnd,unlim=FALSE)

    dim.def.nc(nc_out,"dimMarkov",dimlength=dimMarkov,unlim=FALSE)
    dim.def.nc(nc_out,"dimDistr",dimlength=dimDistr,unlim=FALSE)

    var.def.nc(nc_out,"attribution","NC_DOUBLE",c(0,1))
    att.put.nc(nc_out, "attribution", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "attribution", "dim_explanation", "NC_CHAR", "ID-season")
    att.put.nc(nc_out, "attribution", "explanation", "NC_CHAR", "group attribution for eachs season")

    var.def.nc(nc_out,"MarkovChain","NC_DOUBLE",c(1,3,4))
    att.put.nc(nc_out, "MarkovChain", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "MarkovChain", "dim_explanation", "NC_CHAR", "season,groups,chainparameters")
    att.put.nc(nc_out, "MarkovChain", "explanation", "NC_CHAR", "group characteristics MarkovChain")

    var.def.nc(nc_out,"Distribution","NC_DOUBLE",c(1,2,3,5))
    att.put.nc(nc_out, "Distribution", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "Distribution", "dim_explanation", "NC_CHAR", "season,states,groups,distributions")
    att.put.nc(nc_out, "Distribution", "explanation", "NC_CHAR", "group characteristics Distribution")
       
    var.put.nc(nc_out,"attribution",regionAttribution)      
    var.put.nc(nc_out,"MarkovChain",GroupMarkovChain)      
    var.put.nc(nc_out,"Distribution",GroupDistributions)      
 
    close.nc(nc_out) 
}

create_regional_distr_out_of_kmeans <- function(dataset="_TMean",trendID="91_5",additional_style="",markov_style=5,nGroup=6,add_name="_MarkovDistr",period="1950-2014",region_name="_kmeans"){
    sourceName=paste("../data/",trendID,"/",dataset,additional_style,"/nearest_neighbors/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,".nc",sep="")
    print(sourceName)
    nc=open.nc(sourceName)
    IDregions=var.get.nc(nc,"attribution")
    regional_attribution(dat=dat,region_name=region_name,trendID=trendID,dataset=dataset,additional_style=additional_style,IDregions=IDregions,regNumb=nGroup,comment=sourceName)
}

create_kompakt_map_plot_of_kmeans <- function(dataset="_TMean",trendID="91_5",additional_style="",markov_style=5,nGroup=6,add_name="_MarkovDistr",period="1950-2014",region_name="_kmeans"){
    sourceName=paste("../data/",trendID,"/",dataset,additional_style,"/nearest_neighbors/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,".nc",sep="")
    print(sourceName)
    nc=open.nc(sourceName)
    IDregions=var.get.nc(nc,"attribution")[,1:4]
    #IDregions[2:1000,c(1,2,3,4)]=NA
    map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbors/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,"_map_kompakt",".pdf",sep=""),worldmap=worldmap,reihen=t(matrix(IDregions,c(1319,4))),farb_palette=c("mixed",nGroup,"groups"),col_row=c(5,4),cex=1,pointsize=0.5,cex_axis=1,paper=c(6.3,4),mat=c(1,1,2,2,1,1,2,2,3,3,4,4,3,3,4,4,5,5,5,5),color_lab="")

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
    worldmap<<-getMap(resolution = "low")
    dat<<-dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
}

kmeans_master <- function(nGroupEnd=5,add_name="R"){
    nearest_neighbors(markov_style=5,add_name=add_name,seasons=1:5,nGroupEnd=nGroupEnd,versions=10,runs=30,plot=c("endGroups","endMaps","robustMaps"))
    #create_regional_distr_out_of_kmeans(add_name=add_name,region_name=paste(add_name,"_",nGroupEnd,sep=""),nGroup=nGroupEnd)
}

init()

#create_kompakt_map_plot_of_kmeans(nGroup=7,add_name="KarlPerason_AmpMark")
for (groups in 1:5){
    kmeans_master(nGroupEnd=7,add_name=paste("R",groups,sep=""))
}





