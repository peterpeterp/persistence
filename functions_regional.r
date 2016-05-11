
write_midlat_region_file <- function(region_name="midlat"){
    nGroup<-2
    attribution<-array(NA,1319)
    attribution[dat$lat>30 & dat$lat<60]=1
    attribution[dat$lat< -30 & dat$lat> -60]=2
    mids<-array(NA,c(nGroup,3))
    for (G in 1:nGroup){
        inside<<-which(attribution==G)
        mids[G,]=c(G,mean(dat$lon[inside]),mean(dat$lat[inside]))
    }
    write.table(attribution,paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))
    write.table(mids,paste("../data/",dataset,"/ID_regions/",region_name,"_mids_mean.txt",sep=""))
}

points_to_regions <- function(region_name="7rect"){
    # loads region coordinates and writes a file in which grid points are associated to regions
    # outpufile has following columns: ID, regions from region_names
    ntot=length(dat$ID)
    points=cbind(x=dat$lon,y=dat$lat)

    poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))

    region<-array(NA,dim=c(ntot))
    mids<-array(NA,c(dim(poli)[1],3))

    for (k in 1:dim(poli)[1]){
        x<-c()
        y<-c()
        for (i in 1:6){
            if (!is.na(poli[k,i])){
                x[i]=poli[k,i]
                y[i]=poli[k,(i+6)]
            }
        }
        poligon=cbind(x=x,y=y)
        inside=pnt.in.poly(points,poligon)$pip
        region[which(inside==1)]=poli[k,13]
        mids[k,]=c(k,mean(dat$lon[which(inside==1)]),mean(dat$lat[which(inside==1)]))
    }

    write.table(region,paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))
    write.table(mids,paste("../data/",dataset,"/ID_regions/",region_name,"_mids.txt",sep=""))

    return(region)
}

duration_region <- function(regions,reg,dur,dur_mid){
    # combines all recorded durations of one region to one duration array, same for dur_mid
    inside=which(regions==reg)
    if (length(inside)<3){return(list(duration=(1:10)*NA,duration_mid=(1:10)*NA))}

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
    nona=which(!is.na(duration) & !is.na(duration_mid))
    duration=duration[nona]
    duration_mid=duration_mid[nona]
    return(list(duration=duration,duration_mid=duration_mid))
}

regional_attribution <- function(region_name,toDo=c(TRUE,TRUE,TRUE),IDregions=c("from file"),regNumb=7,comment="polygons",years=length(dat$year)){
    # performs the entire regional analysis of markov and duration
    # result will be written in nc file

    ntot=length(dat$ID)

    if (IDregions[1]=="from polygons"){
        poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))
        regNumb=dim(poli)[1]
        IDregions_tmp=points_to_regions(dat,c(region_name))
        IDregions=array(NA,c(ntot,5))
        for (i in 1:5){
            IDregions[,i]=IDregions_tmp
        }
    }

    if (IDregions[1]=="from file"){
        attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
        IDregions<-array(NA,c(ntot,5))
        regNumb<-length(unique(attribution[!is.na(attribution)]))
        for (i in 1:5){
            IDregions[,i]=attribution
        }
    }

    for (sea in 1:5){
        season<-season_names[sea]
        filename<-paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_duration_",season,".nc",sep="")
        nc_dur<-open.nc(filename)
        dur<-var.get.nc(nc_dur,"dur")
        dur_mid<-var.get.nc(nc_dur,"dur_mid")   

        # regional duration series, for general analysis
        if (toDo[1]==TRUE){
            reg_dur<-array(NA,dim=c(regNumb,2,365*years*100))
            reg_dur_mid<-array(NA,dim=c(regNumb,2,365*years*100))
            maxis<-array(NA,dim=c(2*regNumb))

            for (state in 1:2){
                for (reg in 1:regNumb){   
                    tmp<-duration_region(regions=IDregions[,sea],reg=reg,dur=dur[1:ntot,state,],dur_mid=dur_mid[1:ntot,state,])
                    duration<-tmp$duration
                    duration_mid<-tmp$duration_mid
                    if (length(duration)>100){
                        ord=order(duration_mid)
                        maxis[(state-1)*regNumb+reg]=length(duration)
                        reg_dur[reg,state,1:maxis[(state-1)*regNumb+reg]]=duration[ord]
                        reg_dur_mid[reg,state,1:maxis[(state-1)*regNumb+reg]]=duration_mid[ord]
                    }
                }
            }
            len=max(maxis,na.rm=TRUE)

            # write regional durations in form of gridded durations
            filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_duration_",season,".nc",sep="")
            duration_write(filename=filename,dur=reg_dur[,,1:len],dur_mid=reg_dur_mid[,,1:len],len=len,ID_length=regNumb,ID_name=region_name,comment=comment)
        }

        # required for yearly shuffeling where all periods within one season are treated as dependent from each other
        if (toDo[2]==TRUE){
            if (toDo[1]==FALSE){
                filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_duration_",season,".nc",sep="") ; print(filename)
                reg_dur<-var.get.nc(nc.open(filename),"dur")
                reg_dur_mid<-var.get.nc(nc.open(filename),"dur_mid")
            }
            # create binned duration file
            binned_dur<-array(NA,dim=c(regNumb,2,365*100,years))
            periods_in_yr<-array(0,regNumb*2*years)
            index=0
            for (reg in 1:regNumb){
                for (state in 1:2){
                    for (yr in 1:length(dat$year)){
                        index<-index+1
                        inYr<-which(reg_dur_mid[reg,state,]>=(yr+dat$year[1]-1) & reg_dur_mid[reg,state,]<(yr+dat$year[1]-1+1))  
                        if (length(inYr)>0){
                            binned_dur[reg,state,1:length(inYr),yr]=reg_dur[reg,state,inYr]
                            periods_in_yr[index]=length(inYr)
                        }
                    }
                }
            }
            # write regional duration in form of yearly binned durations
            len=max(periods_in_yr,na.rm=TRUE)
            filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_reg_binned_duration_",season,".nc",sep="")
            reg_binned_dur_write(filename=filename,binned_dur=binned_dur[,,1:len,],len=len,ID_length=regNumb,ID_name=region_name,comment=comment)
        }

        # huge daily matrix required for advanced bootstrapping. last significance try. Here following periods are treated as independent form each other
        if (toDo[3]==TRUE){
            # transform actual period midpoints into indices
            dur_mid[which((-1)^dur<0)]=dur_mid[which((-1)^dur<0)]+0.5/365

            # create matrix with days as colomns with period length if period midpoint is on day
            daily_binned_periods<-array(NA,c(2,regNumb,120,365*65))
            for (state in 1:2){
                for (reg in 1:regNumb){
                    cat(paste(reg,"-"))
                    for (q in IDregions[,5]){
                        indices<-(dur_mid[q,state,which(!is.na(dur_mid[q,state,]))]-1950)*365+0.5
                        daily_binned_periods[state,reg,q,indices]=dur[q,state,1:length(indices)]

                    }
                }
            }

            # delete days outside of the season
            seasonal_sequence<-seasonal_boundaries[sea,1]:seasonal_boundaries[sea,2]
            seasonal_days<-seasonal_sequence
            if (sea!=4){endYear<-length(dat$year)-1}
            if (sea==4){endYear<-length(dat$year)-2}
            for (yr in 1:endYear){seasonal_days<-c(seasonal_days,seasonal_sequence+365*yr)}

            filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_reg_daily_binned_duration_",season,".nc",sep="") ; print(filename)
            reg_daily_binned_dur_write(filename=filename,daily_binned_periods=daily_binned_periods[,,,seasonal_days],seasonal_days=seasonal_days)

        }
    }
}

#--------------------------------------------------------------------------

plot_regional_fit_vergleich <- function(period1,period2,region_name,trendID,additional_style,dataset,fit_style,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){
    # performs the entire regional analysis of markov and duration
    # result will be written in nc file


    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",region_name,"/",period1,"/",trendID,"_",dataset,"_",region_name,"_",period1,fit_style,".nc",sep=""))
    fit_stuff1=var.get.nc(nc,"fit_stuff")
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",region_name,"/",period2,"/",trendID,"_",dataset,"_",region_name,"_",period2,fit_style,".nc",sep=""))
    fit_stuff2=var.get.nc(nc,"fit_stuff")


    regions=var.get.nc(nc,"region")
    regNumb=dim(fit_stuff2)[2]
    seaNumb=5
    
    season_names=c("MAM","JJA","SON","DJF","4seasons")

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",region_name,"/",region_name,"_",period1,"_diff_",period2,fit_style,".pdf",sep=""))
    fit_stuff=fit_stuff1-fit_stuff2
    fit_stuff[,,,2]=1/fit_stuff1[,,,2]-1/fit_stuff2[,,,2]
    fit_stuff[,,,4]=1/fit_stuff1[,,,4]-1/fit_stuff2[,,,4]
    for (sea in 1:seaNumb){
        for (state in 1:2){
            plot(NA,xlim=c(0.5,regNumb+0.5),ylim=c(-2,2),xlab="",ylab="lifetime [days]",axes=FALSE,frame.plot=TRUE,cex=0.5)
            axis(2)
            axis(1, at=1:regNumb, labels=region_names) 
            abline(h=0,col="grey")

            points(fit_stuff[sea,,state,2],col="red")
            lines(fit_stuff[sea,,state,2],col="red")
            points(fit_stuff[sea,,state,4],col="orange")
            lines(fit_stuff[sea,,state,4],col="orange")
            par(new=TRUE)
            plot(NA,xlim=c(0.5,regNumb+0.5),ylim=c(-8,8),axes=FALSE,frame.plot=TRUE,cex=0.5,ylab="",xlab="")
            axis(4)
            points(fit_stuff[sea,,state,5],col="blue")
            lines(fit_stuff[sea,,state,5],col="blue")
        }
    }
    graphics.off()
}



plot_regional_fit_parameters <- function(period,trendID,additional_style,dataset,fit_style,region_name="7rect",region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){
    # performs the entire regional analysis of markov and duration
    # result will be written in nc file

    print(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,fit_style,".nc",sep=""))
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,fit_style,".nc",sep=""))
    fit_stuff=var.get.nc(nc,"fit_stuff")


    regNumb=dim(fit_stuff)[2]
    seaNumb=5
    season_names=c("MAM","JJA","SON","DJF","4seasons")

    color=c(rgb(0.5,0.5,1,0.8),rgb(1,0.5,0.5,0.8))
    color=c("blue","red")
    pos=c(-1,1)*0.15
    width=7
    height=5

    # regional focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",region_name,"/",period,"/",trendID,"_",period,fit_style,".pdf",sep=""),width=width,height=height)
    par(mfrow=c(1,1))
    par(mar=c(3,5,0,1))   



    for (sea in 1:seaNumb){
        for (state in 1:2){
            plot(NA,xlim=c(0.5,regNumb+0.5),ylim=c(2,15),frame.plot=TRUE,axes=FALSE,ylab="lifetime [days]",xlab="")
            axis(1, at=1:7, labels=region_names) 
            axis(2)
            lines(1:7,1/fit_stuff[sea,,state,2],col=color[1])
            lines(1:7,1/fit_stuff[sea,,state,4],col=color[2])
            lines(1:7,fit_stuff[sea,,state,5],col="green")
        }
    }

    graphics.off()
}



plot_fits_for_region <- function(reg,IDregions=c("from polygons"),period="1950-2014",trendID="91_5",dataset="_TMean",fit_style="2expo_thresh_5-15",region_name="srex"){
    if (IDregions[1]=="from polygons"){
        ntot=1319
        poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))
        regNumb=dim(poli)[1]
        IDregions_tmp=points_to_regions(dat,c(region_name))
        IDregions=array(NA,c(ntot,5))
        for (i in 1:5){
            IDregions[,i]=IDregions_tmp
        }
    }

    print(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style,".nc",sep=""))
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_fit_",fit_style,".nc",sep=""))
    fit_stuff_reg=var.get.nc(nc,"fit_stuff")
    print(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_fit_",fit_style,".nc",sep=""))
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_fit_",fit_style,".nc",sep=""))
    fit_stuff_individual=var.get.nc(nc,"fit_stuff")
    print(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_distributions.nc",sep=""))
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/gridded/",period,"/",trendID,"_",dataset,"_",period,"_distributions.nc",sep=""))
    distr_stuff_individual=var.get.nc(nc,"distr_stuff")

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",period,"/",trendID,"_",region_name,"_",period,"_region_",reg,".pdf",sep=""))
    
    color=c(rgb(0.5,1,0.8),rgb(0.3,0.9,0.8),rgb(0.9,0.6,0.8),rgb(0.8,0.2,0.6),rgb(0.5,0.5,0.5,0.2))

    X=(0:30)+0.5
    for (sea in 1:5){

        ID_select=which(IDregions[,sea]==reg)
        print(ID_select)

        for (state in 1:2){
            plot(NA,ylim=c(0,10),xlim=c(0,10),axes=FALSE,frame.plot=FALSE)
            text(5,5,paste("reg:",reg,"\nseason:",sea,"state:",state,"\npoints:",length(ID_select)))
            plot(NA,xlab="days",ylab="counts",ylim=c(0,150),xlim=c(0,30),axes=TRUE,frame.plot=TRUE)
            for (q in ID_select){
                points(distr_stuff_individual[sea,q,state,1,],distr_stuff_individual[sea,q,state,3,],pch=16,col=color[5],cex=1.5)
            }

            plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE)
            for (q in ID_select){points(distr_stuff_individual[sea,q,state,1,],distr_stuff_individual[sea,q,state,2,],pch=16,col=color[5],cex=1.5)}
            for (q in ID_select){
                abline(v=fit_stuff_individual[sea,q,state,9],col=color[5],lty=2,lwd=0.3)
                lines(X,combi_expo(X,fit_stuff_individual[sea,q,state,5],fit_stuff_individual[sea,q,state,6],fit_stuff_individual[sea,q,state,8],fit_stuff_individual[sea,q,state,9]),col=color[1],lwd=0.3)
                lines(X,expo(X,fit_stuff_individual[sea,q,state,1],fit_stuff_individual[sea,q,state,2]),col=color[3],lwd=0.3,lty=2)
            }
            abline(v=fit_stuff_reg[sea,reg,state,9],col=color[5],lty=2,lwd=2)
            lines(X,combi_expo(X,fit_stuff_reg[sea,reg,state,5],fit_stuff_reg[sea,reg,state,6],fit_stuff_reg[sea,reg,state,8],fit_stuff_reg[sea,reg,state,9]),col=color[2],lwd=2,lty=1)
            lines(X,expo(X,fit_stuff_reg[sea,reg,state,1],fit_stuff_reg[sea,reg,state,2]),col=color[4],lty=1,lwd=2)
            

            plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,30),axes=TRUE,frame.plot=TRUE,log="y")
            for (q in ID_select){points(distr_stuff_individual[sea,q,state,1,],distr_stuff_individual[sea,q,state,2,],pch=16,col=color[5],cex=1.5)}
            for (q in ID_select){
                abline(v=fit_stuff_individual[sea,q,state,9],col=color[5],lty=2,lwd=0.3)
                lines(X,combi_expo(X,fit_stuff_individual[sea,q,state,5],fit_stuff_individual[sea,q,state,6],fit_stuff_individual[sea,q,state,8],fit_stuff_individual[sea,q,state,9]),col=color[1],lwd=0.3)
                lines(X,expo(X,fit_stuff_individual[sea,q,state,1],fit_stuff_individual[sea,q,state,2]),col=color[3],lwd=0.3,lty=2)
            }
            abline(v=fit_stuff_reg[sea,reg,state,9],col=color[5],lty=2,lwd=2)
            lines(X,combi_expo(X,fit_stuff_reg[sea,reg,state,5],fit_stuff_reg[sea,reg,state,6],fit_stuff_reg[sea,reg,state,8],fit_stuff_reg[sea,reg,state,9]),col=color[2],lwd=2,lty=1)
            lines(X,expo(X,fit_stuff_reg[sea,reg,state,1],fit_stuff_reg[sea,reg,state,2]),col=color[4],lty=1,lwd=2)

        }
    }
    graphics.off()    
}


