

group_reduction <- function(X,attribution,start,nGroup,reduce,main_add=""){

    score=array(999,c(nGroup,nGroup))
    for (m in 1:reduce){
        for (p in 1:nGroup){
            for (p2 in 1:nGroup){
                if (p!=p2){
                    #score[p]=sum((start[p,]-toOrder[q,])^2)
                    if (length(which(attribution==p))+length(which(attribution==p2))<500){score[p,p2]=sum(X*as.vector((start[p,]-start[p2,])^2))}
                    else {score[p,p2]=999}
                }
            }
        }
        similar=which(score==score[which.min(score)],arr.ind=TRUE)

        #print(score)
        #print(similar)

        start[similar[1],]=colMeans(start[similar[1:2],],na.rm=TRUE)
        start[similar[2],]=NA
        plot_aktuelles_muster(X,start,nGroup,main=paste(main_add,"  ",similar[1],length(which(attribution==similar[1])),similar[3],length(which(attribution==similar[3]))))
        attribution[which(attribution==similar[2])]=similar[1]
    }
    return(list(attribution=attribution,start=start))
}


plot_aktuelles_muster <- function(X,start,nGroup,main=""){
    plot(NA,xlab="days",ylab="probability density",ylim=c(0.00001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,main=main)
    for (p in 1:nGroup){
        lines(X,start[p,],col=color[p])
    }
    plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,log="y",main=main)
    for (p in 1:nGroup){
        lines(X,start[p,],col=color[p])
    }
}

distr_nearest_neighbours <- function(period="1950-2014",trendID="91_5",dataset="_TMean",fit_style="2expo_thresh_5-15",region_name="srex",ID_select=1:1319,state=1,sea=4,nGroup=7,versions=30,runs=30,reduce=0,name_zusatz="_testMasse_",start_mod="random",plot=c(1,1,1),write=TRUE){

    #nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style,".nc",sep=""))
    #fit_stuff_reg=var.get.nc(nc,"fit_stuff")
    #distr_stuff_reg=var.get.nc(nc,"distr_stuff")

    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_fit_",fit_style,".nc",sep=""))
    fit_stuff_individual=var.get.nc(nc,"fit_stuff")
    distr_stuff_individual=var.get.nc(nc,"distr_stuff")

    X=(1:100)

    toOrder=distr_stuff_individual[sea,,state,2,]
    
    toOrder[is.na(toOrder)]=0

    noEmpty=which(distr_stuff_individual[sea,,state,2,1]>0)

    ntot=dim(toOrder)[1]
    distrSize=dim(toOrder)[2]


    if (plot[1]==1){
        color=c("red","blue","green","violet","black","orange","lightblue","grey",rgb(1,0.5,0.6),rgb(0.5,0.7,0.9),rgb(0.5,0.2,0.8),rgb(0.2,0.5,0.6))
        jet.colors <- colorRampPalette( c( "blue","green","yellow","red") )
        nbcol <- nGroup
        color <<- jet.colors(nbcol)  
        reihen=array(NA,dim=c(versions,ntot))
        pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbours/",period,"/",trendID,"_",region_name,"_",period,"nearest_neighbours",name_zusatz,".pdf",sep=""))
    }

    if (versions!=1){
        attribution_speicher=array(NA,c(versions,ntot))
        start_speicher=array(NA,c(versions,nGroup,distrSize))
    }

    for (version in 1:versions){
        cat(paste("\n version:",version,"\n"))
        if (start_mod=="random"){start=distr_stuff_individual[sea,c(sample(noEmpty,nGroup)),state,2,]}
        if (start_mod!="random"){start=start_mod}

        start[is.na(start)]=0

        attribution_old=array(1:10,ntot)
        attribution=array(NA,ntot)
        for (i in 1:runs){
            cat(paste(i," "))
            for (q in 1:ntot){
                score=array(NA,nGroup)
                for (p in 1:nGroup){
                    #score[p]=sum((start[p,]-toOrder[q,])^2)
                    score[p]=sum(X*as.vector((start[p,]-toOrder[q,])^2))
                    #score[p]=sum(abs((start[p,]-toOrder[q,])*(0.25-start[p,])^2))

                }

                attribution[q]=which.min(score)
            }
            for (p in 1:nGroup){
                if (length(which(attribution==p))==0){print("no match")}
                if (length(which(attribution==p))==1){start[p,]=toOrder[which(attribution==p),]}
                if (length(which(attribution==p))>1){start[p,]=colMeans(toOrder[which(attribution==p),],na.rm=TRUE)}

            }
            if (sum(attribution_old-attribution,na.rm=TRUE)==0){break}
            else {attribution_old<-attribution}
            #plot_aktuelles_muster(X,start,nGroup,color)
        }

        if (reduce>0){
            tmp=group_reduction(X=X,attribution=attribution,start=start,nGroup=nGroup,reduce=reduce,main_add=version)
            attribution=tmp$attribution
            start=tmp$start
        }

        if (plot[2]==1){
            for (p in 1:nGroup){
                if (!is.na(start[p,1])){
                    plot(NA,xlab="days",ylab="probability density",ylim=c(0.00001,0.25),xlim=c(0,50),axes=TRUE,frame.plot=TRUE,main=length(which(attribution==p))) 
                    for (q in which(attribution==p)){
                        #print("missing points")
                        points(X,toOrder[q,],pch=16,col=rgb(0.5,0.5,0.5,0.2),cex=1.5)
                    }
                    lines(X,start[p,],col=color[p])

                    plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,50),axes=TRUE,frame.plot=TRUE,log="y",main=length(which(attribution==p))) 
                    for (q in which(attribution==p)){
                        #print("missing points")
                        points(X,toOrder[q,],pch=16,col=rgb(0.5,0.5,0.5,0.2),cex=1.5)
                    }
                    lines(X,start[p,],col=color[p])
                }
            }
        }
 
        if (FALSE==TRUE){
            laggedDiffs=c()
            for (p in 1:nGroup){
                laggedDiffs[p]=sum(diff(start[p,3:10]))
            }
            print(laggedDiffs)
            order=order(laggedDiffs)
            print(order)
            
            order=order(start[,1])

            attri_order=attribution*NA
            for (i in 1:length(order)){
                attri_order[attribution==i]=order[i]
            }
            reihen[version,]=attri_order
        }

        if (plot[3]==1){reihen[1,]=attribution}

        if (versions!=1){
            attribution_speicher[version,]=attribution
            start_speicher[version,,]=start
        }
    }

    if (plot[3]==1){map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbours/",period,"/",trendID,"_",region_name,"_",period,"nearest_neighbours_map",name_zusatz,".pdf",sep=""),worldmap=worldmap,reihen=reihen,pointsize=1.5,farb_palette="regenbogen")}
    
    if (write==TRUE){
        nc_out <- create.nc(paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",period,"/",trendID,"_",region_name,"_",period,"nearest_neighbours_map",name_zusatz,"result.nc",sep=""))
        att.put.nc(nc_out, "NC_GLOBAL", "comment", "NC_CHAR", "bla")

        dim.def.nc(nc_out,"seasons",dimlength=6,unlim=FALSE)
        dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
        dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)

        dim.def.nc(nc_out,"versions",dimlength=versions,unlim=FALSE)
        dim.def.nc(nc_out,"nGroup",dimlength=nGroup,unlim=FALSE)
        dim.def.nc(nc_out,"distrSize",dimlength=distrSize,unlim=FALSE)


        var.def.nc(nc_out,"attribution","NC_DOUBLE",c(3,1))
        att.put.nc(nc_out, "attribution", "missing_value", "NC_DOUBLE", -99999.9)
        att.put.nc(nc_out, "attribution", "dim_explanation", "NC_CHAR", "version-ID")
        att.put.nc(nc_out, "attribution", "explanation", "NC_CHAR", "group numbers of differnt versions")

        var.def.nc(nc_out,"groups","NC_DOUBLE",c(3,4,5))
        att.put.nc(nc_out, "groups", "missing_value", "NC_DOUBLE", -99999.9)
        att.put.nc(nc_out, "groups", "dim_explanation", "NC_CHAR", "version-nGroup-distrSize")
        att.put.nc(nc_out, "groups", "explanation", "NC_CHAR", "end distributions of groups")
            
        var.put.nc(nc_out,"attribution",attribution_speicher)      
        var.put.nc(nc_out,"groups",start_speicher)      
     
        close.nc(nc_out) 
    }
    
}

nearest_neighbours_post_opti <- function(period="1950-2014",trendID="91_5",dataset="_TMean",fit_style="2expo_thresh_5-15",region_name="srex",name_zusatz="_testMasse_"){
    
    print(paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",period,"/",trendID,"_",region_name,"_",period,"nearest_neighbours_map",name_zusatz,"result.nc",sep=""))
    nc=open.nc(paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbours/",period,"/",trendID,"_",region_name,"_",period,"nearest_neighbours_map",name_zusatz,"result.nc",sep=""))
    groups=var.get.nc(nc,"groups")
    attribution=var.get.nc(nc,"attribution")


    nGroup=dim(groups)[2]*dim(groups)[1]
    nStart=7
    versions=dim(groups)[1]
    ntot=1319

    matches=array(NA,dim=c(dim(groups)[1],dim(groups)[2],1319))
    matches=array(NA,dim=c(nGroup,ntot))
    index=0
    for (ver in 1:dim(groups)[1]){
        for (p in 1:dim(groups)[2]){
            index=index+1
            match=which(attribution[ver,]==p)
            matches[index,1:length(match)]=match
        }

    }
    overlap=array(0,dim=c(nGroup,nGroup))
    for (i in 1:nGroup){
        for (j in 1:nGroup){
            match_i=as.vector(matches[i,])[!is.na(matches[i,])]
            match_j=as.vector(matches[j,])[!is.na(matches[j,])]
            overlap[i,j]=length(match_i[match_i %in% match_j])
        }
    }

    same=array(NA,c(nGroup,versions))
    for (i in 1:nGroup){
        a=which(overlap[,i]>(2/3*overlap[i,i]))
        same[i,1:length(a)]=a
    }
    contained_in_versions=array(NA,nGroup)
    for (i in 1:nGroup){
        contained_in_versions[i]=length(which(!is.na(same[i,])))
    }    

    plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,log="y")
    start=array(0,c(nStart,100))
    for (i in 1:nStart){
        target=which.max(contained_in_versions)

        gridmass=c()
        for (p in same[target,!is.na(same[target,])]){
            vers=1
            print("-----")
            print(p)
            for (k in 1:29){
                if (p > 7){
                    p=p-7
                    vers=vers + 1
                }
            }
            print(vers)
            print(p)
            #print(which(attribution[vers,]==p))
            #reihen=array(NA,c(1,ntot))
            #reihen[1,1]=8
            #reihen[1,which(attribution[vers,]==p)]=1
            #reihen[1,1]=3

            #map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbours/",period,"/",trendID,"_",region_name,"_",period,"nearest_neighbours_map_ookokokok",name_zusatz,vers,"_",p,".pdf",sep=""),worldmap=worldmap,reihen=reihen,pointsize=1.5,farb_palette="regenbogen")

            start[i,]=start[i,]+groups[vers,p,]

            #plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,log="y")
            #lines(1:100,groups[vers,p,])
        }
        start[i,]=start[i,]/contained_in_versions[target]
        contained_in_versions[same[target,!is.na(same[target,])]]=0
        print(contained_in_versions[1:10])
        print(same[target,!is.na(same[target,])][1:10])
        print(target)
        print(contained_in_versions[target])
        lines(1:100,start[i,])
    }



    distr_nearest_neighbours(period=period,trendID=trendID,dataset=dataset,fit_style=fit_style,region_name=region_name,state=1,sea=4,nGroup=7,versions=1,runs=40,reduce=0,name_zusatz="_end_",plot=c(1,1,1),write=FALSE,start_mod=start)
}