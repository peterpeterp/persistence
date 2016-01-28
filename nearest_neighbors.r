

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
    if (dimDistr>1){
        for (state in 1:2){
            for (logAx in c("","y")){
                plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,log=logAx,main=paste(main_zusatz,"Distribution of",state_names[state],"periods"))
                for (p in 1:nGroup){
                    lines(1:dimDistr,start[p,startDistr[state]:stopDistr[state]],col=color[p])
                }
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
        if (dimDistr>1){
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
        #cat(paste(i," "))
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

nearest_neighbors <- function(markov_style=NA,add_name="_forReal_",seasons=1:5,nGroup=7,versions=30,runsMax=100,ID_select=1:1319,plot=c("testMasseGroups","testMasseMaps","endGroups","endMaps")){

    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    state_names<<-c("cold","warm")

    jet.colors <- colorRampPalette(c("black","blue","green","yellow","orange","red",rgb(0.5,0.5,0.5,0.5)))
    color <<- jet.colors(nGroup) 

    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_distributions.nc",sep=""))
    distr_stuff=var.get.nc(nc,"distr_stuff")

    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/markov/",period,"/",trendID,dataset,"_",period,"_markov_order",markov_style,".nc",sep=""))
    eventResult=var.get.nc(nc,"eventResult")

    ntot<<-dim(eventResult)[2]
    dimMarkov<<-dim(eventResult)[4]
    dimMarkov<<-0
    dimDistr<<-dim(distr_stuff)[5]
    distr_select=3:50
    dimDistr<<-length(distr_select)
    #dimDistr<<-0
    dimensionality<<-dimMarkov+2*dimDistr

    startDistr<<-c(dimMarkov+1,dimMarkov+dimDistr+1)
    stopDistr<<-c(dimMarkov+dimDistr,dimMarkov+dimDistr+dimDistr)

    regionAttribution=array(NA,dim=c(ntot,6))
    regionAttribution[,6]=1:ntot
    GroupMarkovChain=array(NA,dim=c(6,nGroup,dimMarkov))
    GroupDistributions=array(NA,dim=c(6,2,nGroup,dimDistr))

    for (sea in seasons){
        print(season_names[sea])

        toOrder=array(NA,dim=c(ntot,dimensionality))
        #toOrder[ID_select,1:dimMarkov]=eventResult[sea,ID_select,3,]
        toOrder[ID_select,startDistr[1]:stopDistr[1]]=distr_stuff[sea,ID_select,1,2,distr_select]
        toOrder[ID_select,startDistr[2]:stopDistr[2]]=distr_stuff[sea,ID_select,2,2,distr_select]
        
        toOrder[is.na(toOrder)]=0
        toOrder<<-toOrder

        print(proc.time())
        tmp<<-kmeans(x=toOrder,centers=nGroup,iter.max=30,nstart=2)
        print(proc.time())
        for (i in 1:2){tmp2<<-k_nearest_neighbors(nGroup=nGroup,runs=30)}
        
        print(proc.time())
        asfdasfd

        attribution_unordered=tmp$cluster
        groups_unordered=tmp$centers

        attribution=attribution_unordered*NA
        groups=groups_unordered*NA
        for (G in 1:nGroup){
            Gmax=which.max(groups_unordered[,1])
            attribution[attribution_unordered==Gmax]=G
            groups[G,]=groups_unordered[Gmax,]
            groups_unordered[Gmax,]=-999
        }

        if ("endGroups" %in% plot){
            pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbors/",period,"/",trendID,"_",period,"_",season_names[sea],"_",add_name,"_",nGroup,".pdf",sep=""))
            plot_aktuelles_muster(attribution=attribution,start=groups,nGroup=nGroup,points=TRUE)
        graphics.off()
        }

        regionAttribution[,sea]=attribution
        #GroupMarkovChain[sea,,]=groups[,1:dimMarkov]
        GroupDistributions[sea,1,,]=groups[,startDistr[1]:stopDistr[1]]
        GroupDistributions[sea,2,,]=groups[,startDistr[2]:stopDistr[2]]
    }
    if ("endMaps" %in% plot){map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbors/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,"_map",".pdf",sep=""),worldmap=worldmap,reihen=t(matrix(regionAttribution,c(ntot,6))),pointsize=1.5,farb_palette=c("mixed",nGroup,"groups"),region="srex")}


    nc_out<-create.nc(paste("../data/",trendID,"/",dataset,additional_style,"/nearest_neighbors/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", "gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)
    att.put.nc(nc_out, "NC_GLOBAL", "over all distances from points to group mean", "NC_CHAR", paste(group_structure_evaluation(attribution=attribution,groups=groups,nGroup=nGroup)))
    
    dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
    dim.def.nc(nc_out,"seasons",dimlength=6,unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
    
    dim.def.nc(nc_out,"nGroup",dimlength=nGroup,unlim=FALSE)

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



kmeans_raw <- function(add_name="_forReal_",seasons=1:4,nGroup=7,starts=10,runsMax=100,ID_select=1:1319,plot=c("endMaps")){

    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    seasonStart=c(60,152,244,335,1)
    seasonStop=c(151,243,334,424,365)
    state_names<<-c("cold","warm")

    jet.colors <- colorRampPalette(c("black","blue","green","yellow","orange","red",rgb(0.5,0.5,0.5,0.5)))
    color <<- jet.colors(nGroup) 

    ntot<<-length(ID_select)
    year_select<<-10:50

    dimensionality<<-length(year_select)*365*5
    regionAttribution=array(NA,dim=c(1319,length(seasons)))

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbors/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,"_spectra",".pdf",sep=""))

    for (sea in seasons){
        print(season_names[sea])
        print(proc.time())

        # take days in one season
        if (sea<4){
            tempSea=dat$tas[ID_select,,year_select]
            tempSea[,1:(seasonStart[sea]-1),]=-999
            tempSea[,(seasonStop[sea]+1):365,]=-999
        }
        if (sea==4){
            tempSea=dat$tas[ID_select,,year_select]
            tempSea[,(seasonStop[sea]-365+1):(seasonStart[sea]-1),]=-999
        }
        if (sea==5){
            tempSea=dat$tas[ID_select,,year_select]
        }

        tempSea<<-tempSea
        partOfYear<<-which(tempSea[1,,4]!=-999)
        condensated<<-tempSea[,partOfYear,]
        condensated<<-array(condensated,c(ntot,length(partOfYear)*length(year_select)))

        toOrder=array(NA,c(ntot,length(partOfYear)*length(year_select)))
        freqs=array(NA,c(ntot,length(partOfYear)*length(year_select)))
        for (q in 1:ntot){
            if (length(which(!is.na(condensated[q,])))>900){
                tmp=spectrum(condensated[q,],na.action=na.exclude,plot=FALSE)
                toOrder[q,1:length(tmp$spec)]=tmp$spec 
                freqs[q,1:length(tmp$freq)]=tmp$freq
            }
        }
        toOrder<<-toOrder
        toOrder[is.na(toOrder)]=0

        plot(NA,xlim=c(0.0001,0.5),ylim=c(0,500),log="x")
        lines(freqs[295,],toOrder[295,],col="red")
        lines(freqs[488,],toOrder[488,],col="blue")
        lines(freqs[489,],toOrder[489,],col="lightblue")
        lines(freqs[661,],toOrder[661,],col="green")

        graphics.off()
        freqs<<-freqs
        toOrder<<-toOrder

        tmp=kmeans(x=toOrder,centers=nGroup,iter.max=runsMax,nstart=starts)
        tmp<<-tmp
        adasd

        

        attribution=tmp$cluster
        regionAttribution[ID_select,sea]=attribution
        print(proc.time())

    }
    regionAttribution<<-regionAttribution
    if ("endMaps" %in% plot){map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbors/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,"_map",".pdf",sep=""),worldmap=worldmap,reihen=t(matrix(regionAttribution,c(1319,length(seasons)))),pointsize=1.5,farb_palette=c("mixed",nGroup,"groups"),region="srex")}


    nc_out<-create.nc(paste("../data/",trendID,"/",dataset,additional_style,"/nearest_neighbors/",period,"/",trendID,"_",period,"_",add_name,"_",nGroup,".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", "gridpoints")
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)
    att.put.nc(nc_out, "NC_GLOBAL", "over all distances from points to group mean", "NC_CHAR", "dasdasd")
    
    dim.def.nc(nc_out,"ID",dimlength=1319, unlim=FALSE)
    dim.def.nc(nc_out,"seasons",dimlength=length(seasons),unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
    

    var.def.nc(nc_out,"attribution","NC_DOUBLE",c(0,1))
    att.put.nc(nc_out, "attribution", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "attribution", "dim_explanation", "NC_CHAR", "ID-season")
    att.put.nc(nc_out, "attribution", "explanation", "NC_CHAR", "group attribution for eachs season")

    var.put.nc(nc_out,"attribution",regionAttribution)      
 
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
    worldmap<<-getMap(resolution = "low")
    dat<<-dat_load(paste("../data/HadGHCND",dataset,"_data3D.day1-365.1950-2014.nc",sep=""))
}

kmeans_master <- function(nGroup=5,add_name="default"){
    nearest_neighbors(markov_style=5,add_name=add_name,seasons=1:5,nGroup=nGroup,versions=10,runsMax=100,ID_select=ID_select,plot=c("endMaps"))
    #create_regional_distr_out_of_kmeans(add_name=add_name,region_name=paste(add_name,"_",nGroup,sep=""),nGroup=nGroup)
}

init()

#create_kompakt_map_plot_of_kmeans(nGroup=7,add_name="KarlPerason_AmpMark")

ID_select=which(dat$lat>0 & dat$lat<100)
#ID_select=1:1319

for (tries in 1:3){
    kmeans_master(nGroup=5,add_name=paste("distr20",tries,sep=""))
    #kmeans_raw(nGroup=7,add_name=paste("raw___",tries,sep=""),starts=10,runsMax=100,ID_select=ID_select,seasons=c(5))
}




