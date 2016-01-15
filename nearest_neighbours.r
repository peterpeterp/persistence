

group_reduction <- function(attribution,start,nGroup,reduce,main_add=""){

    score=array(999,c(nGroup,nGroup))
    for (m in 1:reduce){
        for (p in 1:nGroup){
            for (p2 in 1:nGroup){
                if (p!=p2){
                    #score[p]=sum((start[p,]-toOrder[q,])^2)
                    if (length(which(attribution==p))+length(which(attribution==p2))<1000){score[p,p2]=sum(X*as.vector((start[p,]-start[p2,])^2))}
                    else {score[p,p2]=999}
                }
            }
        }
        similar=which(score==score[which.min(score)],arr.ind=TRUE)

        start[similar[1],]=colMeans(start[similar[1:2],],na.rm=TRUE)
        start[similar[2],]=NA
        plot_aktuelles_muster(start,nGroup,main=paste(main_add,"  ",similar[1],length(which(attribution==similar[1])),similar[3],length(which(attribution==similar[3]))))
        attribution[which(attribution==similar[2])]=similar[1]
        start=start[which((1:nGroup)!=similar[2]),]
    }
    return(list(attribution=attribution,start=start))
}


plot_aktuelles_muster <- function(start,nGroup,main=""){
    plot(NA,xlab="days",ylab="probability density",ylim=c(0.00001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,main=main)
    for (p in 1:nGroup){
        lines(X,start[p,],col=color[p])
    }
    plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,log="y",main=main)
    for (p in 1:nGroup){
        lines(X,start[p,],col=color[p])
    }
}


k_nearest_neighbours <- function(versions=30,nGroup=7,start_mod="random",runs=30){
    noEmpty=which(toOrder[,1]>0)

    attribution=array(NA,c(versions,ntot))
    groups=array(NA,c(versions,nGroup,distrSize))

    for (version in 1:versions){
        cat(paste("\n version:",version,"\n"))
        if (start_mod[1]=="random"){start=toOrder[c(sample(noEmpty,nGroup)),]}
        if (start_mod[1]!="random"){start=start_mod}

        start[is.na(start)]=0

        attribution_old=array(1:10,ntot)
        for (i in 1:runs){
            cat(paste(i," "))
            for (q in 1:ntot){
                score=array(NA,nGroup)
                for (p in 1:nGroup){
                    # distance definitions here
                    score[p]=sum((start[p,]-toOrder[q,])^2)
                    score[p]=sum(start[p,]-toOrder[q,])
                    #score[p]=sum(X*as.vector((start[p,]-toOrder[q,])^2))
                    #score[p]=sum((X^0.5)*as.vector((start[p,]-toOrder[q,])^2))
                    #score[p]=sum(abs((start[p,]-toOrder[q,])*(0.25-start[p,])^2))

                    #wheight=(start[p,]+toOrder[q,])/2
                    #noZero=which(wheight>0)
                    #score[p]=sum(as.vector((start[p,noZero]-toOrder[q,noZero])^2)/wheight[noZero])
                   

                }
                attribution[version,q]=which.min(score)
            }
            for (p in 1:nGroup){
                if (length(which(attribution[version,]==p))==0){print("no match")}
                if (length(which(attribution[version,]==p))==1){start[p,]=toOrder[which(attribution[version,]==p),]}
                if (length(which(attribution[version,]==p))>1){start[p,]=colMeans(toOrder[which(attribution[version,]==p),],na.rm=TRUE)}

            }
            if (sum(attribution_old-attribution[version,],na.rm=TRUE)==0){break}
            else {attribution_old<-attribution[version,]}
        }
        groups[version,,]=start
    }

    return(list(attribution=attribution,groups=groups))    
}



nearest_neighbours <- function(period="1950-2014",trendID="91_5",dataset="_TMean",fit_style="2expo_thresh_5-15",add_name="_forReal_",seasons=1:5,states=1:2,nGroup=7,nReduce=0,versions=30,runs=30,plot=c("testMasseGroups","testMasseMaps","endGroups","endMaps"),season_names=c("MAM","JJA","SON","DJF","4seasons"),state_names=c("cold","warm")){

    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_fit_",fit_style,".nc",sep=""))
    fit_stuff_individual=var.get.nc(nc,"fit_stuff")
    distr_stuff_individual=var.get.nc(nc,"distr_stuff")

    for (sea in seasons){
        for (state in states){

            toOrder=distr_stuff_individual[sea,,state,2,]
            toOrder[is.na(toOrder)]=0
            toOrder<<-toOrder

            ntot<<-dim(toOrder)[1]
            distrSize<<-dim(toOrder)[2]
            
            X<<-1:distrSize

            tmp=k_nearest_neighbours(versions=versions,nGroup=nGroup,start_mod="random",runs=runs)
            attribution=tmp$attribution
            groups=tmp$groups
            if ("testMasseMaps" %in% plot){map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbours/",period,"/",trendID,"_",period,"_",season_names[sea],"_",state_names[state],"_groups_map",add_name,"testMasse",".pdf",sep=""),worldmap=worldmap,reihen=attribution,pointsize=1.5,farb_palette=c("mixed",nGroup,"groups"))}

            nMass=versions*nGroup

            matches=array(NA,dim=c(nMass,ntot))
            index=0
            for (ver in 1:versions){
                for (p in 1:dim(groups)[2]){
                    index=index+1
                    match=which(attribution[ver,]==p)
                    matches[index,1:length(match)]=match
                }
            }

            overlap=array(0,dim=c(nMass,nMass))
            for (i in 1:nMass){
                for (j in 1:nMass){
                    match_i=as.vector(matches[i,])[!is.na(matches[i,])]
                    match_j=as.vector(matches[j,])[!is.na(matches[j,])]
                    overlap[i,j]=length(match_i[match_i %in% match_j])
                }
            }

            same=array(NA,c(nMass,versions))
            for (i in 1:nMass){
                a=which(overlap[,i]>(2/3*overlap[i,i]))
                same[i,1:length(a)]=a
            }

            contained_in_versions=array(NA,nMass)
            for (i in 1:nMass){
                contained_in_versions[i]=length(which(!is.na(same[i,])))
            }    

            start=array(0,c(nGroup,100))
            for (i in 1:nGroup){
                print(i)
                print(dim(groups))
                target=which.max(contained_in_versions)
                for (p in same[target,!is.na(same[target,])]){
                    vers=1
                    for (k in 1:29){
                        if (p > nGroup){
                            p=p-nGroup
                            vers=vers + 1
                        }
                    }
                    print(paste(vers,p))
                    start[i,]=start[i,]+groups[vers,p,]

                }
                start[i,]=start[i,]/contained_in_versions[target]
                contained_in_versions[same[target,!is.na(same[target,])]]=0
            }
            pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbours/",period,"/",trendID,"_",period,"_",season_names[sea],"_",state_names[state],"_reduction",add_name,".pdf",sep=""))
            jet.colors <- colorRampPalette(c("black",rgb(0.5,1,1),"red", "yellow","green",rgb(1,0.5,1),"orange"))
            color <<- jet.colors(nGroup)  
            for (r in 0:nReduce){
                print(r)
                tmp=k_nearest_neighbours(nGroup=nGroup-r,versions=1,runs=40,start_mod=start)
                attribution=tmp$attribution
                groups=tmp$groups
                if (r < nReduce){
                    tmp=group_reduction(attribution=attribution,start=groups[1,,],nGroup=nGroup-r,reduce=1)
                    attribution=tmp$attribution
                    start=tmp$start
                }
            }
            nGroup=nGroup-nReduce
            if ("endMaps" %in% plot){map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbours/",period,"/",trendID,"_",period,"_",season_names[sea],"_",state_names[state],"_groups_map",add_name,".pdf",sep=""),worldmap=worldmap,reihen=attribution,pointsize=1.5,farb_palette=c("mixed",nGroup,"groups"))}
            if ("endGroups" %in% plot){
                pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbours/",period,"/",trendID,"_",period,"_",season_names[sea],"_",state_names[state],"_groups",add_name,".pdf",sep=""))
                jet.colors <- colorRampPalette(c("black",rgb(0.5,1,1),"red", "yellow","green",rgb(1,0.5,1),"orange"))
                color <<- jet.colors(nGroup)  
                plot_aktuelles_muster(start=groups[1,,],nGroup=nGroup)
                for (p in 1:nGroup){
                    if (!is.na(groups[1,p,1])){
                        plot(NA,xlab="days",ylab="probability density",ylim=c(0.00001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,main=length(which(attribution==p))) 
                        for (q in which(attribution==p)){
                            #print("missing points")
                            points(X,toOrder[q,],pch=16,col=rgb(0.5,0.5,0.5,0.05),cex=1.5)
                        }
                        lines(X,groups[1,p,],col=color[p])

                        plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,log="y",main=length(which(attribution==p))) 
                        for (q in which(attribution==p)){
                            #print("missing points")
                            points(X,toOrder[q,],pch=16,col=rgb(0.5,0.5,0.5,0.05),cex=1.5)
                        }
                        lines(X,groups[1,p,],col=color[p])
                    }
                }
            }

        }
    }
}


