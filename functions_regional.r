
wave_region <- function(wavenumber,trough,lats){
    # creates a table with coordinates of rectangles representing wave ridges and troughs
    wellen=array(NA,dim=c(40,13))
    brei=360/wavenumber/4
    lon=trough
    lon0=trough
    i=1
    wellen[,7:8]=lats[1]
    wellen[,9:10]=lats[2]
    while ((round(x=lon,digits=3)==round(x=lon0,digits=3) & i>1)==FALSE){
        if (lon>(180-brei) & lon<(180+brei)){
            wellen[i,1:4]=c(180,lon-brei,lon-brei,180)
            wellen[i,13]=i  
            i=i+1   
            wellen[i,1:4]=c(lon+brei-360,-180,-180,lon+brei-360)
            wellen[i,13]=i  
            i=i+1   
            lon=lon+2*brei
        }
        if (lon<=(180-brei)){
            wellen[i,1:4]=c(lon+brei,lon-brei,lon-brei,lon+brei)
            wellen[i,13]=i  
            i=i+1   
            lon=lon+2*brei
        }
        if (lon>=(180+brei)){
            lon=lon-360
            wellen[i,1:4]=c(lon+brei,lon-brei,lon-brei,lon+brei)
            wellen[i,13]=i  
            i=i+1   
            lon=lon+2*brei
        }
        
    }
    write.table(wellen[1:length(which(!is.na(wellen[,13]))),],paste("../data/region_poligons/",wavenumber,"wave.txt",sep=""))
}

points_to_regions <- function(dat,region_names=c("mid_lat_belt","srex","7rect","6wave","7wave","8wave")){
    # loads region coordinates and writes a file in which grid points are associated to regions
    # outpufile has following columns: ID, regions from region_names
    library(SDMTools)
    ntot=length(dat$ID)
    points=cbind(x=dat$lon,y=dat$lat)

    region = array(NA,dim=c(ntot,(length(region_names)+1)))

    for (j in 1:length(region_names)){
        poli=read.table(paste("../data/region_poligons/",region_names[j],".txt",sep=""))
        for (k in 1:dim(poli)[1]){
            x=c()
            y=c()
            for (i in 1:6){
                if (!is.na(poli[k,i])){
                    x[i]=poli[k,i]
                    y[i]=poli[k,(i+6)]
                }
            }
            poligon=cbind(x=x,y=y)
            inside=pnt.in.poly(points,poligon)$pip
            region[which(inside==1),(j+1)]=poli[k,13]
        }
    }

    if (length(region_names)>1){
        region[1:ntot,1]=dat$ID
        write.table(region,"../data/ID-regions.txt")
    }
    return(region[1:ntot,2])
}


regions_color <- function(reihen,reihen_sig,worldmap,titles,poli,filename_plot){
    # plots worldmap and colored regions on it
    jet.colors <- colorRampPalette( c(rgb(0.2,0.6,0.2),rgb(0.5,1,0.5), rgb(0.98,0.98,0.98) ,rgb(1,0.5,1),rgb(0.6,0.2,0.6)))
    
    nbcol <- 101
    color <- jet.colors(nbcol)

    pdf(file = filename_plot,width=12,height=8)

    for (rei in 1:dim(reihen)[1]){            
        y=c()
        index=c()
        signi=c()
        j=0
        for (i in 1:dim(poli)[1]){
            poliLabel=i
            if (!is.na(reihen[rei,poliLabel])){
                j=j+1
                y[j]=reihen[rei,poliLabel]
                index[j]=i 
                if (abs(reihen[rei,poliLabel])>0.0001){
                    signi[j]=sprintf("%.04f",reihen_sig[rei,poliLabel])
                }         
            }
        }
        aushol=max(c(abs(max(y)),abs(min(y))))
        y[j+1]=-aushol
        y[j+2]=aushol
        facetcol <- cut(y,nbcol)  

        print(titles[rei])
        plot(worldmap,main=titles[rei])

        for (i in 1:j){
            lon=poli[index[i],1:6]
            lat=poli[index[i],7:12]
            lon=lon[!is.na(lon)]
            lat=lat[!is.na(lat)]
             polygon(x=lon,y=lat,col=color[facetcol[i]],border="green")
            text(mean(lon),mean(lat),label=signi[i],cex=0.7,col="black")
        }
        image.plot(legend.only=T, zlim=range(y), col=color)
    }
    graphics.off()
    return()
}

#============================================================================================================================


duration_region <- function(regions,reg,dur,dur_mid){
    # combines all recorded durations of one region to one duration array, same for dur_mid
    inside=which(regions==reg)
    duration=array(NA,dim=c(1000000))
    duration_mid=array(NA,dim=c(1000000))
    count=1
    # combines the recorded periods from all the grid points of one region in one array
    for (i in inside){
        values=length(which(!is.na(dur[i,])))
        duration[count:(count+values)]=dur[i,1:values]
        duration_mid[count:(count+values)]=dur_mid[i,1:values]
        count=count+values
    }
    duration=duration[!is.na(duration)]
    duration_mid=duration_mid[!is.na(duration_mid)]
    return(list(duration=duration,duration_mid=duration_mid))
}



regional_quantiles_fits <- function(dat,yearPeriod,region_name,trendID,additional_style,dataset,season_auswahl=c(1,2,3,4,5,6),plot=TRUE,write=TRUE,add_name="",region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){
    # performs the entire regional analysis of markov and duration
    # result will be written in ncdf file

    # pnly for one trend and 2 states until now
    ntot=length(dat$ID)


    poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))
    regNumb=dim(poli)[1]
    print(regNumb)
    IDregions=points_to_regions(dat,c(region_name))

    
    state_names=c("cold","warm")
    season_names=c("MAM","JJA","SON","DJF","year","4seasons")

    taus=c(0.05,0.25,0.5,0.75,0.95,0.91,0.98,0.1)
    quantiles=array(NA,dim=c(length(season_names),regNumb,2,length(taus),3))

    fitstuff=array(NA,dim=c(length(season_names),regNumb,2,5,20))

    if (plot==TRUE){
        pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/",region_name,"_dist_diff_fit_plot_",dataset,"_",yearPeriod[1],"-",yearPeriod[2],add_name,".pdf",sep=""),width=2.2,height=2.2)
    }

    for (sea in season_auswahl){   
        season=season_names[sea]
        dists=list()

        nc_dur=open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/","duration/",trendID,dataset,"_duration_",season,".nc",sep=""))
        dur=get.var.ncdf(nc_dur,"dur")
        dur_mid=get.var.ncdf(nc_dur,"dur_mid")
        
        for (state in 1:2){
            for (reg in 1:regNumb){
            
                tmp=duration_region(IDregions,reg,dur[1:ntot,state,],dur_mid[1:ntot,state,])
                duration=tmp$duration
                duration_mid=tmp$duration_mid
                ord=order(duration_mid)
                if (length(duration)>1000){
                    y=as.vector(duration[ord])
                    x=as.vector(duration_mid[ord])
                    inYearPeriod=which(x>yearPeriod[1] & x<yearPeriod[2])
                    y=y[inYearPeriod]
                    x=x[inYearPeriod]

                    # y are now the duration length in selected yearPeriod
                    # quantile determination
                    dists=append(dists,list(y))
                    quantiles[sea,reg,state,1:length(taus)]=quantile_pete(y,taus=taus,na.rm=TRUE)
                    quantiles_vergleich[sea,reg,state,1:length(taus)]=quantile(y,prob=taus,na.rm=TRUE)

                    #for (i in 1:length(taus)){
                    #    qu=try(summary(rq(y~x,taus[i]))$coefficients)
                    #    if (class(qu)=="try-error"){qu=array(NA,8)}
                    #    quantile_regressions[sea,reg,state,c(i,length(taus)+i)]=qu[c(2,8)]
                    #}

                    # exponential fit + other values
                    br=seq(0,max(y,na.rm=TRUE),1)
                    histo=hist(y,breaks=br,plot=FALSE)

                    # data to be fitted
                    Y=histo$density
                    X=histo$mids

                    tmp=exponential_fit(X,Y)
                    fitstuff[sea,reg,state,1,1:2]=tmp$pars
                    fitstuff[sea,reg,state,1,19:20]=tmp$ana




                    # difference fit with gaussian
                    z=histo$density-expfit
                    xy=data.frame(y=z,x=X)
                    gau_nls=try(nls(y~(A*exp(-(B-x)^2/sigma^2)),algorithm="port",data=xy,start=c(A=0.01,B=12,sigma=7),lower=c(A=0,B=0,sigma=0),na.action=na.exclude))
                    if (class(gau_nls)!="try-error"){
                        A=summary(gau_nls)$parameters[1]
                        B=summary(gau_nls)$parameters[2]                   
                        sigma=summary(gau_nls)$parameters[3]   
                        zfit=(A*exp(-(B-X)^2/sigma^2))      
                    }  
                    else{
                        A=0.005
                        B=14
                        sigma=8
                    }


                    # fit expo + gaussian together
                    # first with big restrainst to get a forced result
                    xy=data.frame(y=Y,x=X)
                    dev=1/1000
                    combi_nls=try(nls(y~(a*exp(-b*x)+A*exp(-(B-x)^2/sigma^2)),algorithm="port",data=xy,start=c(a=a,b=b,A=A,B=B,sigma=sigma),lower=c(a=a-dev*sqrt(a^2),b=b-dev*sqrt(b^2),A=A-dev*sqrt(b^2),B=B-dev*sqrt(b^2),sigma=sigma-dev*sqrt(sigma^2)),upper=c(a=a+dev*sqrt(a^2),b=b+dev*sqrt(b^2),A=A+dev*sqrt(b^2),B=B+dev*sqrt(b^2),sigma=sigma+dev*sqrt(sigma^2)),na.action=na.exclude,nls.control(maxiter = 10000, tol = 1e-04, minFactor=1/10024, warnOnly=TRUE)))
                    #if (class(combi_nls)!="try-error"){
                    if (1==1){
                        a=summary(combi_nls)$parameters[1]
                        b=summary(combi_nls)$parameters[2]  
                        A=summary(combi_nls)$parameters[3]
                        B=summary(combi_nls)$parameters[4]                   
                        sigma=summary(combi_nls)$parameters[5]   
                        combifit=(a*exp(-b*X)+A*exp(-(B-X)^2/sigma^2)) 
                        combiR2=1-sum(((Y-combifit)^2),na.rm=TRUE)/sum(((Y-mean(Y,na.rm=TRUE))^2),na.rm=TRUE)
                    } 
                    else{ 
                        print(reg)
                        print(state)        
                        combifit=x*NA
                        combiR2=NA
                    }
                    fitstuff[sea,reg,state,4:6]=c(A,B,sigma)

                    # now with weaker restrainst to get a good result
                    opti_nls=try(nls(y~(a*exp(-b*x)+A*exp(-(B-x)^2/sigma^2)),algorithm="port",data=xy,start=c(a=a,b=b,A=A,B=B,sigma=sigma),lower=c(a=-Inf,b=-Inf,A=0,B=mean(Y,na.rm=TRUE),sigma=0),upper=c(a=Inf,b=Inf,A=max(z,na.rm=TRUE),B=Inf,sigma=Inf),na.action=na.exclude,nls.control(maxiter = 10000, tol = 1e-04, minFactor=1/10024, warnOnly=FALSE)))
                    if (class(opti_nls)!="try-error"){
                        a=summary(opti_nls)$parameters[1]
                        b=summary(opti_nls)$parameters[2]  
                        A=summary(opti_nls)$parameters[3]
                        B=summary(opti_nls)$parameters[4]                   
                        sigma=summary(opti_nls)$parameters[5]   
                        optifit=(a*exp(-b*X)+A*exp(-(B-X)^2/sigma^2)) 
                        optiR2=1-sum(((Y-optifit)^2),na.rm=TRUE)/sum(((Y-mean(Y,na.rm=TRUE))^2),na.rm=TRUE)
                    } 
                    else {optiR2=NA}
                    fitstuff[sea,reg,state,7:12]=c(a,b,NA,A,B,sigma)
                    fitstuff[sea,reg,state,13:15]=c(expR2,combiR2,optiR2)

                    n=length(X)
                    if (class(exp_nls)!="try-error"){
                        fitstuff[sea,reg,state,16]=BIC(exp_nls)
                        #if (class(gau_nls)!="try-error"){
                        if (1==1){
                            fitstuff[sea,reg,state,17]=BIC(combi_nls)
                            if (class(opti_nls)!="try-error"){
                                fitstuff[sea,reg,state,18]=BIC(opti_nls)
                            }
                        }
                    }
                    fitstuff[sea,reg,state,19]=sum(A*exp(-(B-X)^2/sigma^2))

                    mean=mean(y,na.rm=TRUE)
                    sd=sd(y,na.rm=TRUE)
                    skew=skewness(y,na.rm=TRUE)
                    others[sea,reg,state,1:4]=c(mean,sd,sd/mean,skew)


                    # plotstuff
                    if (plot==TRUE){
                        color=c("blue","red")
                        #reference
                        par(mar=c(0, 0, 0, 0) + 0.1)
                        plot(NA,xlim=c(0,10),ylim=c(0,10),axes=FALSE,frame.plot=FALSE)
                        text(5,9.4,paste(season,region_names[reg],state))
                        order_=order(y,decreasing=TRUE)
                        for (extreme in 1:8){
                            text(3,extreme,round(x[order_][extreme],01))
                            text(8,extreme,round(y[order_][extreme],01))
                        }

                    #first plot page
                        par(mar=c(3, 3, 0, 0) + 0.1)
                        par(mfrow=c(1,1))
                        
                        plot(histo$density,xlab="days",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],cex=0.5)
                        if (1==1){
                            at_=axis(2,labels=FALSE,col="black")
                            if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
                            axis(2,at=at_)
                        }
                        if (1==2){
                            at_=axis(1,labels=FALSE,col="black")
                            if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
                            axis(1,at=at_)
                        }
                        lines(expfit,col="black")
                        legend("topright",legend=c(paste("R2=",round(expR2,03),"\n b=",round(fitstuff[sea,reg,state,2],03),sep="")),bty="n")
                        #legend("bottomleft",legend=c(region_names[reg]),bty="n")   
                        text(12,0.00002,region_names[reg])                 



                    # second version
                        par(mar=c(3, 3, 0, 0) + 0.1)
                        par(mfrow=c(1,1))
                        
                        plot(histo$density,xlab="days",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],log="y",cex=0.5)
                        if ((region_name=="7rect" & 4==4) | (region_name!="7rect" & reg==1)){
                            at_=axis(2,labels=FALSE,col="black")
                            if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
                            axis(2,at=at_)
                        }
                        if (2==2 & 4>3){
                            at_=axis(1,labels=FALSE,col="black")
                            if (length(at_)>4){at_=at_[2:(length(at_)-1)]}
                            axis(1,at=at_)
                        }
                        lines(expfit,col="black")
    

                        #abline(v=11,col="gray",lty=2)
                        legend("topright",legend=c(paste("R2=",round(expR2,03),"\n b=",round(-fitstuff[sea,reg,state,2],03),sep="")),bty="n")
                        #legend("bottomleft",legend=c(region_names[reg]),bty="n")
                        text(12,0.00002,region_names[reg])                 

                    }

                    #print(fitstuff[sea,reg,state,])
                    print(quantile_regressions[sea,reg,state,])
                    print(paste(sea,reg,state))

                }
            }
        }
    }
    graphics.off()

    if (write==TRUE){
        ncRegion <- dim.def.ncdf("region",units="region",vals=1:regNumb, unlim=FALSE)
        ncStates <- dim.def.ncdf("states",units="states",vals=1:2,unlim=FALSE)
        ncSeason <- dim.def.ncdf("seasons",units="seasons",vals=1:6,unlim=FALSE)

        ncTaus <- dim.def.ncdf("taus",units="dimension for both quantile",vals=1:8,unlim=FALSE)
        ncOuts <- dim.def.ncdf("outs",units="dimension for others",vals=1:4,unlim=FALSE)
        ncFits <- dim.def.ncdf("fits",units="dimension for fitstuff",vals=1:19,unlim=FALSE)
       
        poli_points <- dim.def.ncdf("poli_points",units="id",vals=1:12,unlim=FALSE)

        region_coordinates <- var.def.ncdf(name="region_coordinates",units="deg",longname="1:6 lon - 7:12 lat",dim=list(ncRegion,poli_points),missval=-9999.0)

        ncQuantile <- var.def.ncdf(name="quantiles",units="quantile values in days evaluated by quantile_pete",longname="0.05,0.25,0.5,0.75,0.95,0.91,0.98,0.1",dim=list(ncSeason,ncRegion,ncStates,ncTaus), missval=-9999.0)
        ncQuantileReg <- var.def.ncdf(name="quantiles",units="quantile regression with slope and significance of slope",longname="0.05,0.25,0.5,0.75,0.95,0.91,0.98,0.1",dim=list(ncSeason,ncRegion,ncStates,ncTaus), missval=-9999.0)
        ncQuantile_vergleich <- var.def.ncdf(name="quantiles_2",units="quantile values in days evaluated by quantile()",longname="0.05,0.25,0.5,0.75,0.95,0.91,0.98,0.1",dim=list(ncSeason,ncRegion,ncStates,ncTaus), missval=-9999.0)
        ncOthers <- var.def.ncdf(name="others",units="different analysis values",longname="mean sd sd/mean skewness",dim=list(ncSeason,ncRegion,ncStates,ncOuts), missval=-9999.0)
        ncFitstuff <- var.def.ncdf(name="fitstuff",units="different analysis values",longname="Exp: [a,b,NA] combi: [A,B,sigma] opti:[a,b,NA,A,B,sigma] R2: [exp,combi,opti] BIC: [exp,combi,opti] percent of gaussian area",dim=list(ncSeason,ncRegion,ncStates,ncFits), missval=-9999.0)
        
        vars=list(ncQuantile,ncQuantileReg,ncQuantile_vergleich,ncOthers,ncFitstuff,region_coordinates)
       
        nc = create.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",region_name,"_",yearPeriod[1],"-",yearPeriod[2],"_quantiles.nc",sep=""),vars)

        put.var.ncdf(nc,ncQuantile,quantiles)      
        put.var.ncdf(nc,ncQuantileReg,quantile_regressions)      
        put.var.ncdf(nc,ncQuantile_vergleich,quantiles_vergleich)      
        put.var.ncdf(nc,ncOthers,others)      
        put.var.ncdf(nc,ncFitstuff,fitstuff)      

        pol_poi=array(NA,c(dim(poli)[1],12))
        for (i in 1:dim(poli)[1]){
            for (j in 1:12){
                
                if (is.numeric(poli[i,j])){
                    pol_poi[i,j]=poli[i,j]
                }
            }
        }
        put.var.ncdf(nc,region_coordinates,pol_poi)      

        close.ncdf(nc) 
    }
}

#==========================================================================================================================

# ------------------------------------------------------------------------------------------------
plot_boxplot <- function(quans,x,width,color="white",border="black",density=NA){

    links=x-width/2
    rechts=x+width/2
    polygon(x=c(rechts,links,links,rechts),y=c(quans[2],quans[2],quans[4],quans[4]),col=color,border=border,density=density)    
    #polygon(x=c(rechts,links,links,rechts),y=c(quans[2],quans[2],quans[4],quans[4]),col=color,border=border,density=density)    
    for (qu in quans[1:5]){
        lines(c(links,rechts),c(qu,qu),col=border)
    }
    lines(c(links,links),quans[c(2,4)],col=border)
    lines(c(rechts,rechts),quans[c(2,4)],col=border)
    lines(c(x,x),quans[c(1,2)],lty=2,col=border)
    lines(c(x,x),quans[c(4,5)],lty=2,col=border)
}


plot_regional_boxplots <- function(dat,yearPeriod,region_name,trendID,additional_style,dataset,quantile_style="quantiles",region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){
    # performs the entire regional analysis of markov and duration
    # result will be written in ncdf file

    # pnly for one trend and 2 states until now
    ntot=length(dat$ID)

    nc = open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",region_name,"_",yearPeriod[1],"-",yearPeriod[2],"_quantiles.nc",sep=""))
    quantiles=get.var.ncdf(nc,quantile_style)
    others=get.var.ncdf(nc,"others")
    fitstuff=get.var.ncdf(nc,"fitstuff")
    regions=get.var.ncdf(nc,"region")


    regNumb=dim(quantiles)[2]
    seaNumb=6
    
    season_names=c("MAM","JJA","SON","DJF","year","4seasons")

    color=c(rgb(0.5,0.5,1,0.8),rgb(1,0.5,0.5,0.8))
    pos=c(-1,1)*0.15
    width=6
    height=3

    # regional focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",region_name,"_",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_",quantile_style,"_boxplots_regional.pdf",sep=""),width=width,height=height)
    par(mfrow=c(1,1),cex.legend=2)
    par(mar=c(1,4,2,3))   
    at_=seq(1, regNumb, 1)
    at_=c(at_-0.15,at_+0.15)
    buchstaben=c("a","b","c","d")
    drueber=2
    for (sea in 1:length(season_names)){   
        season=season_names[sea]

        plot(NA,xlim=c(0,7.5),ylim=c(0,25),frame.plot=FALSE,axes=FALSE,ylab="days")#(max(quantiles[sea,,,1:5])+drueber)
        axis(2)
        for (i in axis(2)){
            abline(h=i,col=rgb(0.5,0.5,0.5,0.5),lty=3)
        }
        legend("topright",legend=c(buchstaben[sea]),bty="n")
        #text(8,23,label=buchstaben[sea],cex=2)#(max(quantiles[sea,,,1:5],na.rm=TRUE)+drueber-0.5)

        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,],reg+pos[state],0.3,color[state])
                text(reg,24,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))#max(quantiles[sea,,,1:5])+drueber
            }
        }        
        for (oth in c(1)){
            points(at_[1:regNumb],others[sea,,1,oth],col="blue",pch=oth)
            #lines(at_[1:regNumb],others[sea,,1,oth],col="blue",lty=3)
            points(at_[(regNumb+1):(regNumb*2)],others[sea,,2,oth],col="red",pch=oth)
            #lines(at_[(regNumb+1):(regNumb*2)],others[sea,,2,oth],col="red",lty=3)
        }
    }
    graphics.off()

    # seasonal anomaly
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",region_name,"_",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_",quantile_style,"_boxplots_seasonal_anomaly.pdf",sep=""),width=width,height=height)
    par(mfrow=c(1,1))
    par(mar=c(1,4,2,0))   
    at_=seq(1, regNumb, 1)
    at_=c(at_-0.15,at_+0.15)
    buchstaben=c("e","f","g","h","i","j")
    drueber=0.5

    anomaly=quantiles*NA
    other_anom=others*NA
    for (sea in 1:4){  
        anomaly[sea,,,]=quantiles[sea,,,]-quantiles[6,,,]
    }
    for (sea in 1:4){   
        season=season_names[sea]
        plot(NA,xlim=c(0,8),ylim=c(min(anomaly[sea,,,1:5],na.rm=TRUE),(max(anomaly[sea,,,1:5],na.rm=TRUE)+drueber)),frame.plot=FALSE,axes=FALSE,ylab="days")
        axis(2)
        for (i in axis(2)){
            abline(h=i,col=rgb(0.5,0.5,0.5,0.5),lty=3)
        }        
        text(8,(max(anomaly[sea,,,1:5],na.rm=TRUE)+drueber-0.5),label=buchstaben[sea],cex=2)

        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(anomaly[sea,reg,state,],reg+pos[state],0.3,color[state])
                text(reg,max(anomaly[sea,,,1:5],na.rm=TRUE)+drueber,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))
            }
        }
        for (quA in c(4,5)){
            points(at_[1:regNumb],anomaly[sea,,1,quA],col="black",pch=2)
            points(at_[(regNumb+1):(regNumb*2)],anomaly[sea,,2,quA],col="black",pch=2)
        }

        other_anom[sea,,,]=others[sea,,,]-others[6,,,]
        for (oth in c(1)){
            points(at_[1:regNumb],other_anom[sea,,1,oth],col="blue",pch=oth)
            lines(at_[1:regNumb],other_anom[sea,,1,oth],col="blue",lty=3)
            points(at_[(regNumb+1):(regNumb*2)],other_anom[sea,,2,oth],col="red",pch=oth)
            lines(at_[(regNumb+1):(regNumb*2)],other_anom[sea,,2,oth],col="red",lty=3)
        }
    }  

    # regional focus kind of seasonal anomaly
    color=c(rgb(0.5,0.5,1,0.8),rgb(1,0.5,0.5,0.8),rgb(0.5,0.5,1,0.4),rgb(1,0.5,0.5,0.4))

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",region_name,"_",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_",quantile_style,"_boxplots_seas_kind_anomaly.pdf",sep=""))
    par(mfrow=c(1,1))
    par(mar=c(1,4,2,0))   
    at_=seq(1, regNumb, 1)
    at_=c(at_-0.15,at_+0.15)

    for (sea in 1:4){   
        season=season_names[sea]

        plot(NA,xlim=c(0,8),ylim=c(0,(max(quantiles[sea,,,1:5])+4)),frame.plot=FALSE,axes=FALSE,main=season,ylab="days")
        axis(2)
        for (i in axis(2)){
            abline(h=i,col=rgb(0.5,0.5,0.5,0.5),lty=3)
        }        
        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[6,reg,state,],reg+pos[state],0.3,color=color[2+state],border=rgb(0.5,0.5,0.5,0.5))
            }
        }
        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,],reg+pos[state],0.3,color=color[2+state])
                text(reg,max(quantiles[sea,,,1:5])+2,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))
            }
        }

    }
    graphics.off()

    # regional focus plus seasonal anomaly
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",region_name,"_",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_",quantile_style,"_boxplots_seas_combi_anomaly.pdf",sep=""))
    par(mfrow=c(2,1),cex=1)
    par(mar=c(1,4,2,0))   
    at_=seq(1, regNumb, 1)
    at_=c(at_-0.15,at_+0.15)
    buchstaben=c("a","b","c","d","e","f","g","h","i","j")
    drueber=c(1,1.5)
    index=0
    for (sea in 1:4){   
        season=season_names[sea]
        par(mar=c(1,4,2,0))
        index=index+1
        plot(NA,xlim=c(0,8),ylim=c(min(anomaly[sea,,,1:5],na.rm=TRUE),(max(anomaly[sea,,,1:5],na.rm=TRUE)+drueber[1])),frame.plot=FALSE,axes=FALSE,ylab="days")#,main=season)
        axis(2)
        for (i in axis(2)){
            abline(h=i,col=rgb(0.5,0.5,0.5,0.5),lty=3)
        }        
        text(8,(max(anomaly[sea,,,1:5],na.rm=TRUE)+drueber[1]),label=buchstaben[index],cex=2)
        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(anomaly[sea,reg,state,],reg+pos[state],0.3,color[state])
                #text(reg,max(anomaly[sea,,,1:5],na.rm=TRUE)+0.5,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))
            }
        }

        for (quA in c(4,5)){
            points(at_[1:regNumb],anomaly[sea,,1,quA],col="black",pch=2)
            points(at_[(regNumb+1):(regNumb*2)],anomaly[sea,,2,quA],col="black",pch=2)
        }
        par(mar=c(1,4,1,0))
        index=index+1
        plot(NA,xlim=c(0,8),ylim=c(0,(max(quantiles[sea,,,1:5])+drueber[2])),frame.plot=FALSE,axes=FALSE,ylab="days")
        text(8,(max(quantiles[sea,,,1:5],na.rm=TRUE)+drueber[2]),label=buchstaben[index],cex=2)
        axis(2)
        for (i in axis(2)){
            abline(h=i,col=rgb(0.5,0.5,0.5,0.5),lty=3)
        }   
        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,],reg+pos[state],0.3,color=color[state])
                text(reg,max(quantiles[sea,,,1:5])+drueber[2],region_names[reg],col=rgb(0.5,0.5,0.5,0.5),cex=1.3)
            }
        }


    }
    graphics.off()
}

plot_regional_boxplots_vergleich <- function(dat,yearPeriod1,yearPeriod2,region_name,trendID,additional_style,dataset,quantile_style="quantiles",region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){
    # performs the entire regional analysis of markov and duration
    # result will be written in ncdf file

    # pnly for one trend and 2 states until now
    ntot=length(dat$ID)

    nc = open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",yearPeriod1[1],"-",yearPeriod1[2],"/",trendID,"_",region_name,"_",yearPeriod1[1],"-",yearPeriod1[2],"_quantiles.nc",sep=""))
    quantiles1=get.var.ncdf(nc,quantile_style)
    others1=get.var.ncdf(nc,"others")
    nc = open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",yearPeriod2[1],"-",yearPeriod2[2],"/",trendID,"_",region_name,"_",yearPeriod2[1],"-",yearPeriod2[2],"_quantiles.nc",sep=""))
    quantiles2=get.var.ncdf(nc,quantile_style)
    others2=get.var.ncdf(nc,"others")
    regions=get.var.ncdf(nc,"region")
    regNumb=dim(quantiles1)[2]
    seaNumb=6
    
    season_names=c("MAM","JJA","SON","DJF","year","4seasons")

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/","/",region_name,"_",yearPeriod2[1],"-",yearPeriod2[2],"_diff_",yearPeriod1[1],"-",yearPeriod1[2],"_",quantile_style,"_boxplots_regional.pdf",sep=""))
    par(mfrow=c(1,1))
    par(mar=c(1,5,4,3))   
    at_=seq(1, regNumb, 1)
    at_=c(at_-0.15,at_+0.15)
    color=c()
    maxi=c()
    taus=c(0.05,0.25,0.5,0.75,0.95,0.91,0.98)
    color=c(rgb(0.5,0.5,1,0.8),rgb(1,0.5,0.5,0.8))
    pos=c(-1,1)*0.15

    quantiles=quantiles2-quantiles1
    print(quantiles)
    others=others2-others1

    # regional focus
    for (sea in 1:length(season_names)){   
        season=season_names[sea]

        plot(NA,xlim=c(0,8),ylim=c(min(quantiles[sea,,,1:5]),(max(quantiles[sea,,,1:5])+0.5)),frame.plot=FALSE,axes=FALSE,main=season,ylab="days")
        axis(2)
        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,],reg+pos[state],0.3,color[state])
                text(reg,max(quantiles[sea,,,1:5])+0.5,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))
            }
        }
        for (oth in c(1)){
            points(at_[1:regNumb],others[sea,,1,oth],col="blue",pch=oth)
            lines(at_[1:regNumb],others[sea,,1,oth],col="blue",lty=3)
            points(at_[(regNumb+1):(regNumb*2)],others[sea,,2,oth],col="red",pch=oth)
            lines(at_[(regNumb+1):(regNumb*2)],others[sea,,2,oth],col="red",lty=3)
        }
        for (quA in c(4,5)){
            points(at_[1:regNumb],quantiles[sea,,1,quA],col="black",pch=2)
            points(at_[(regNumb+1):(regNumb*2)],quantiles[sea,,2,quA],col="black",pch=2)
        }
    }
    graphics.off()

    # seasonal focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/","/",region_name,"_",yearPeriod2[1],"-",yearPeriod2[2],"_diff_",yearPeriod1[1],"-",yearPeriod1[2],"_",quantile_style,"_boxplots_seasonal.pdf",sep=""))
    par(mfrow=c(1,1))
    par(mar=c(1,5,4,3))   
    at_=seq(1, seaNumb, 1)
    at_=c(at_-0.15,at_+0.15)

    for (reg in 1:regNumb){
        plot(NA,xlim=c(0,8),ylim=c(min(quantiles[,reg,,1:5]),(max(quantiles[,reg,,1:5])+0.5)),frame.plot=FALSE,axes=FALSE,main=region_names[reg],ylab="days")
        axis(2)
        for (sea in 1:seaNumb){   
            season=season_names[sea]     
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,],sea+pos[state],0.3,color[state])
                text(sea,max(quantiles[,reg,,1:5])+0.5,season_names[sea],col=rgb(0.5,0.5,0.5,0.5))
            }
        }
        for (oth in c(1)){
            points(at_[1:seaNumb],others[,reg,1,oth],col="blue",pch=1)
            lines(at_[1:seaNumb],others[,reg,1,oth],col="blue",lty=3)
            points(at_[(seaNumb+1):(seaNumb*2)],others[,reg,2,oth],col="red",pch=1)
            lines(at_[(seaNumb+1):(seaNumb*2)],others[,reg,2,oth],col="red",lty=3)
        }
        for (quA in c(4,5)){
            points(at_[1:seaNumb],quantiles[,reg,1,quA],col="black",pch=2)
            points(at_[(seaNumb+1):(seaNumb*2)],quantiles[,reg,2,quA],col="black",pch=2)
        }
    }
    graphics.off()
}





plot_regional_fit_parameters <- function(dat,yearPeriod,region_name,trendID,additional_style,dataset,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){
    # performs the entire regional analysis of markov and duration
    # result will be written in ncdf file

    # pnly for one trend and 2 states until now
    ntot=length(dat$ID)

    nc = open.ncdf(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",region_name,"_",yearPeriod[1],"-",yearPeriod[2],"_quantiles.nc",sep=""))
    quantiles=get.var.ncdf(nc,"quantiles")
    others=get.var.ncdf(nc,"others")
    fitstuff=get.var.ncdf(nc,"fitstuff")
    regions=get.var.ncdf(nc,"region")

    regNumb=dim(quantiles)[2]
    seaNumb=6
    season_names=c("MAM","JJA","SON","DJF","year","4seasons")

    color=c(rgb(0.5,0.5,1,0.8),rgb(1,0.5,0.5,0.8))
    color=c("blue","red")
    pos=c(-1,1)*0.15
    width=4
    height=3

    # regional focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",yearPeriod[1],"-",yearPeriod[2],"_fit_params_regional.pdf",sep=""),width=width,height=height)
    par(mfrow=c(3,1))
    par(mar=c(1,5,0,1))   



    for (sea in 1:seaNumb){
        plot(NA,xlim=c(0.5,7.5),ylim=c(0,10),frame.plot=FALSE,axes=FALSE,ylab="")
        for (reg in 1:regNumb){
            text(reg,1,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))
        }


        noopti=which(is.na(fitstuff[sea,,,18]))
        nocombi=which((fitstuff[sea,,,17]-fitstuff[sea,,,16])>0)
        # a b NA 

        for (val in c(8)){
            y=1/fitstuff[sea,,,val]
            
            y[noopti]=1/fitstuff[sea,,,2]
            mini=min(y,na.rm=TRUE)
            maxi=max(y,na.rm=TRUE)
            range=maxi-mini
            dev=range/100

            plot(NA,xlim=c(0.5,7.5),ylim=c(mini-dev,maxi+dev),frame.plot=FALSE,axes=FALSE,ylab="lifetime [days]",main="")
            axis(2,col="black",ylab="lifetime [days]")
            text(7.5,maxi-dev,"a",cex=1)
            for (state in 1:2){
                points(y[,state],col=color[state],pch=1)
                lines(y[,state],col=color[state],lty=3)
            }
        }
  

        for (val in c(19)){
           
            y=fitstuff[sea,,,val]
            y[noopti]=0
            y[nocombi]=0

            mini=min(y,na.rm=TRUE)
            maxi=max(y,na.rm=TRUE)
            range=maxi-mini
            dev=range/100

            plot(NA,xlim=c(0.5,7.5),ylim=c(mini-dev,maxi+dev),frame.plot=FALSE,axes=FALSE,ylab="% gaussian",main="")
            axis(2,col="black",ylab="")
            text(7.5,maxi-dev,"b",cex=1)
            for (state in 1:2){
                points(y[,state],col=color[state],pch=1)
                lines(y[,state],col=color[state],lty=3)
            }
        }


    }


    # write column for latex
    bicdiff=fitstuff[,,,18]-fitstuff[,,,16]
    y2=fitstuff[,,,17]-fitstuff[,,,16]
    bicdiff[is.na(fitstuff[,,,18])]=y2[is.na(fitstuff[,,,18])]


    color=c("lightblue","lightred")

    table<-file(paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",yearPeriod[1],"-",yearPeriod[2],"/",trendID,"_",region_name,"_",yearPeriod[1],"-",yearPeriod[2],"_seasons_latex.txt",sep=""))
    options(scipen=100)

    lines=c()
    y=fitstuff[sea,,,15]-fitstuff[sea,,,13]
    y2=fitstuff[sea,,,14]-fitstuff[sea,,,13]
    y[is.na(fitstuff[sea,,,15])]=y2[is.na(fitstuff[sea,,,15])]

    for (state in 1:2){
        if (state==1){newline=paste("\\","multirow{2}{*}{$\\Delta R^2$}",sep="")}
        else{newline=" "}
        for (reg in 1:regNumb){
            newline=paste(newline," &{\\cellcolor{",color[state],"}}",sep="")
            newline=paste(newline,round(y[reg,state],04))
            #}
        }
        newline=paste(newline,paste("\\","\\",sep=""))
        if (state==2){newline=paste(newline,"\n\\hline")}
        lines[state]=newline
    }

    y=fitstuff[sea,,,18]-fitstuff[sea,,,16]
    y2=fitstuff[sea,,,17]-fitstuff[sea,,,16]
    y[is.na(fitstuff[sea,,,18])]=y2[is.na(fitstuff[sea,,,18])]    

    for (state in 1:2){
        if (state==1){newline=paste("\\","multirow{2}{*}{$\\Delta$BIC}",sep="")}
        else{newline=""}
        for (reg in 1:regNumb){
            newline=paste(newline," &{\\cellcolor{",color[state],"}}",sep="")
            if (y[reg,state]<0){newline=paste(newline,"\\textit{\\textbf{",round(y[reg,state]),"}}")}
            else{newline=paste(newline,round(y[reg,state]))}
        }
        newline=paste(newline,paste("\\","\\",sep=""))
        if (state==2){newline=paste(newline,"\n\\hline")}
        lines[(state+2)]=newline
    } 

    y=1/fitstuff[sea,,,8]
    y2=1/fitstuff[sea,,,2]
    y[is.na(fitstuff[sea,,,8])]=y2[is.na(fitstuff[sea,,,8])]    

    for (state in 1:2){
        if (state==1){newline=paste("\\","multirow{2}{*}{1/b}",sep="")}
        else{newline=""}
        for (reg in 1:regNumb){
            newline=paste(newline," &{\\cellcolor{",color[state],"}}",sep="")
            newline=paste(newline,round(y[reg,state],02))
        }
        newline=paste(newline,paste("\\","\\",sep=""))
        if (state==2){newline=paste(newline,"\n\\hline")}
        lines[(state+4)]=newline
    } 

    y=fitstuff[sea,,,10]

    for (state in 1:2){
        if (state==1){newline=paste("\\","multirow{2}{*}{A}",sep="")}
        else{newline=""}
        for (reg in 1:regNumb){
            newline=paste(newline," &{\\cellcolor{",color[state],"}}",sep="")
            newline=paste(newline,round(y[reg,state],03))
        }
        newline=paste(newline,paste("\\","\\",sep=""))
        if (state==2){newline=paste(newline,"\n\\hline")}
        lines[(state+6)]=newline
    } 
    y=fitstuff[sea,,,13]

    for (state in 1:2){
        if (state==1){newline=paste("\\","multirow{2}{*}{R2(exp)}",sep="")}
        else{newline=""}
        for (reg in 1:regNumb){
            newline=paste(newline," &{\\cellcolor{",color[state],"}}",sep="")
            newline=paste(newline,round(y[reg,state],03))
        }
        newline=paste(newline,paste("\\","\\",sep=""))
        if (state==2){newline=paste(newline,"\n\\hline")}
        lines[(state+8)]=newline
    } 
    
    writeLines(lines, table)

    close(table)
}

