
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
    # creates different kinds of regional persistence records depending on the required format
    # for bootstrapping a big matrix is required

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
                    count<-0
                    for (q in which(IDregions[,5]==reg)){
                        #print(q)
                        count<-count+1
                        indices<-(dur_mid[q,state,which(!is.na(dur_mid[q,state,]))]-1950)*365+0.5
                        daily_binned_periods[state,reg,count,indices]=dur[q,state,1:length(indices)]

                    }
                }
            }
            # delete days outside of the season
            if (FALSE){
                seasonal_sequence<-seasonal_boundaries[sea,1]:seasonal_boundaries[sea,2]
                seasonal_days<-seasonal_sequence
                if (sea!=4){endYear<-length(dat$year)-1}
                if (sea==4){endYear<-length(dat$year)-2}
                for (yr in 1:endYear){seasonal_days<-c(seasonal_days,seasonal_sequence+365*yr)}
                filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_reg_daily_binned_duration_",season,".nc",sep="") ; print(filename)
                reg_daily_binned_dur_write(filename=filename,daily_binned_periods=daily_binned_periods[,,,seasonal_days],seasonal_days=seasonal_days)
            }

            filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_reg_daily_binned_dur_",season,".nc",sep="") ; print(filename)
            reg_daily_binned_dur_write(filename=filename,daily_binned_periods=daily_binned_periods,seasonal_days=1:365)

        }
    }
}


