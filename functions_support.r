
estimate_mode <- function(x,na.rm=TRUE) {
    if (na.rm==TRUE){
        x=x[which(!is.na(x))]
    }
    if (length(x)<20){return(NA)}
    else {
        d <- density(x)
        return(d$x[which.max(d$y)])
    }
}


r_calc_runmode_2D <- function(y2D,nday,nyr){
    # Input 2D array: y[day,year],
    # Calculate the running mean given a window of nyrs and ndays 
    # around each point in time
    nbuff1 = (nday-1)/2 
    nbuff2 = (nyr-1)/2 
    buffer = c(nbuff1,nbuff2)

    dims0 = dim(y2D)
    dims  = dim(y2D)
    dims[1] = dims[1]+2*nbuff1
    dims[2] = dims[2]+2*nbuff2

    y2Dex = array(NA,dim=c(dims[1],dims[2]))
    y2Dex[(nbuff1+1):(dims[1]-nbuff1),(nbuff2+1):(dims[2]-nbuff2)] = y2D[,]
    y2Dex[1:nbuff1,(nbuff2+2):(dims[2]-nbuff2+1)] = y2D[(dims0[1]-nbuff1+1):dims0[1],]   
    y2Dex[(dims[1]-nbuff1+1):dims[1],(nbuff2):(dims[2]-nbuff2-1)] = y2D[1:nbuff1,] 


    trend=y2D*NA

    for (i in 1:365){
        for (j in 1:65){
            trend[i,j]=estimate_mode(y2Dex[(i+nbuff1-nbuff1):(i+nbuff1+nbuff1),(j+nbuff2-nbuff2):(j+nbuff2+nbuff2)],na.rm=TRUE)
        }
    }
    cat("-")
    return(trend)
}

r_calc_runmedian_2D <- function(y2D,nday,nyr){
    # Input 2D array: y[day,year],
    # Calculate the running mean given a window of nyrs and ndays 
    # around each point in time
    nbuff1 = (nday-1)/2 
    nbuff2 = (nyr-1)/2 
    buffer = c(nbuff1,nbuff2)

    dims0 = dim(y2D)
    dims  = dim(y2D)
    dims[1] = dims[1]+2*nbuff1
    dims[2] = dims[2]+2*nbuff2

    y2Dex = array(NA,dim=c(dims[1],dims[2]))
    y2Dex[(nbuff1+1):(dims[1]-nbuff1),(nbuff2+1):(dims[2]-nbuff2)] = y2D[,]
    y2Dex[1:nbuff1,(nbuff2+2):(dims[2]-nbuff2+1)] = y2D[(dims0[1]-nbuff1+1):dims0[1],]   
    y2Dex[(dims[1]-nbuff1+1):dims[1],(nbuff2):(dims[2]-nbuff2-1)] = y2D[1:nbuff1,] 


    trend=y2D*NA

    for (i in 1:365){
        for (j in 1:65){
            trend[i,j]=median(y2Dex[(i+nbuff1-nbuff1):(i+nbuff1+nbuff1),(j+nbuff2-nbuff2):(j+nbuff2+nbuff2)],na.rm=TRUE)
        }
    }
    cat("-")
    return(trend)
}

r_calc_runmean_2D <- function(y2D,nday,nyr){
    # Input 2D array: y[day,year],
    # Calculate the running mean given a window of nyrs and ndays 
    # around each point in time
    nbuff1 = (nday-1)/2 
    nbuff2 = (nyr-1)/2 
    buffer = c(nbuff1,nbuff2)

    dims0 = dim(y2D)
    dims  = dim(y2D)
    dims[1] = dims[1]+2*nbuff1
    dims[2] = dims[2]+2*nbuff2

    y2Dex = array(NA,dim=c(dims[1],dims[2]))
    y2Dex[(nbuff1+1):(dims[1]-nbuff1),(nbuff2+1):(dims[2]-nbuff2)] = y2D[,]
    y2Dex[1:nbuff1,(nbuff2+2):(dims[2]-nbuff2+1)] = y2D[(dims0[1]-nbuff1+1):dims0[1],]   
    y2Dex[(dims[1]-nbuff1+1):dims[1],(nbuff2):(dims[2]-nbuff2-1)] = y2D[1:nbuff1,] 


    trend=y2D*NA

    for (i in 1:365){
        for (j in 1:65){
            trend[i,j]=mean(y2Dex[(i+nbuff1-nbuff1):(i+nbuff1+nbuff1),(j+nbuff2-nbuff2):(j+nbuff2+nbuff2)],na.rm=TRUE)
        }
    }
    cat("-")
    return(trend)
}

c_calc_runmean_2D <- function(y2D,nday,nyr){
    # Input 2D array: y[day,year],
    # Calculate the running mean given a window of nyrs and ndays 
    # around each point in time
    nbuff1 = (nday-1)/2 
    nbuff2 = (nyr-1)/2 
    buffer = c(nbuff1,nbuff2)

    dims0 = dim(y2D)
    dims  = dim(y2D)
    dims[1] = dims[1]+2*nbuff1
    dims[2] = dims[2]+2*nbuff2

    # create empty array
    y2D_list = array(y2D,dim=(dims0[1]*dims0[2]))
    y2D_list[is.na(y2D_list)]=0
    y2Dex_list = array(-99,dim=(dims[1]*dims[2])) 
    trend_list = array(1,dim=(dims0[1]*dims0[2]))       

    tempi = .C("c_run_mean2d",daten=as.numeric(y2D_list),datenex=as.numeric(y2Dex_list),temp_trend=as.numeric(trend_list),size=as.integer(dims0),tag=as.integer(buffer))
    trend = array(tempi$temp_trend,dim=c(dims0[1],dims0[2]))
    
    cat("-")
    return(trend)
}

seasonal_matrix_out <- function(input,model=markov_2states,states=2,seasons=array(c(151,242,334,425),dim=c(2,2)),interval=365,shift=0){
    #input array is tas over whole time period
    #this cals the function given in "model" and hands to the function an array containing the
    #input of one season. the step is repeated until the end of the time-line
    #output is stored in out where the first dimension is the season

    #since for winter the end-day-index is lower than the begin-day-index (new year), all the 
    #day-indeces are shifted until winter index ends with day-index=365
    #therefore the starting point i is equal to the shift
    #example: winter season start-day: 334 end-day:425
    #shift=425-365=60  -> start-day:274 end-day=365
    #whole time series ist shifted by 60 days
    if (seasons[length(seasons)]>365){
        shift=seasons[length(seasons)]-365
        seasons[,]=seasons[,]-shift
    }
    size=length(input)
    x=seq(1, size, 1)
    i=shift
    j=1
    transitions=states*states
    out=array(NA,dim=c(dim(seasons)[2],transitions,65))
    out_conf=array(NA,dim=c(dim(seasons)[2],65))

    while ((i+interval)<size){
        #goes through the time series until the end is reached
        if ((is.na(input[i+1])==FALSE) & (is.na(input[i+interval])==FALSE)){
            for (sea in 1:length(seasons[1,])){
                #goes through the seasons, selects parts of input and calculates something
                x=input[(seasons[1,sea]+i):(seasons[2,sea]+i)]
                tmp=model(x)

                out[sea,1:transitions,j]=tmp$transMat
                out_conf[sea,j]=tmp$confidence
            }
        }
        j=j+1
        i=i+interval
    }
    return(list(out=out,out_conf=out_conf))
}   

calc_trend <- function(dat,filename,nday,nyr,procedure){
    # calculates running mean for each grid point
    # can choose between the c script and r function
    trash = ((nyr-1)/2*365+(nday-1))
    ntot = length(dat$ID)
    trend=dat$tas*NA
    for (q in 1:ntot) {
        temp = procedure(dat$tas[q,,],nday=nday,nyr=nyr)
        temp[1:trash]=NA
        temp[(length(dat$time)-trash):length(dat$time)]=NA
        trend[q,,]=temp
    }

    trend_write(filename,dat,trend)
    return(trend)
}

find_nas <- function(dat){
    # counts NAs at each grid-point
    ntot=length(dat$ID)
    nas=array(NA,dim=c(ntot,2))
    for (q in 1:ntot){
        cat("-")
        nas[q,1]=dat$ID[q]
        nas[q,2]=length(which(is.na(dat$tas[q,,])))
    }
    write.table(nas,"../data/number_of_NA_per_station_2011.txt")
}

calc_per <- function(dat,trend,nday,nyr,model,states,transition_names,filename){
    source("functions_markov.r")

    ## User parameters 
    #trash is the number of data point which are wasted by detrending
    trash = ((nyr-1)/2*365+(nday-1))
    ntot = length(dat$ID)
    laenge_zeit = length(dat$time)
    transitions=states*states

    # Calculate persistence information
    #cat("Calculating persistence... ")

    markov_per = list(ind=dat$tas*NA,markov=array(NA,dim=c(ntot,5,transitions,65)),markov_conf=array(NA,dim=c(ntot,5,65)))

    for (q in 1:ntot) { 
        cat("-")
        if (length(which(is.na(dat$tas[q,,])))<(2*trash+365*20)){

            # Calculate persistence vector
            y = dat$tas[q,,]
            per_ind = y*NA 

            if (states==3){
                detrended = y-trend[q,,]
                threshold = sd(detrended,na.rm=TRUE)*0.5

                per_ind[detrended>threshold]  =  1
                per_ind[detrended<(-threshold)]  = -1 
                per_ind[detrended<threshold & detrended>(-threshold)] = 0 
            }

            if (states==2){
                per_ind[y > trend[q,,]]=1
                per_ind[y < trend[q,,]]=-1
            # the >= was somehow problematic, since it affects roughly 5% of the datapoints
            # now the datapoints sitting on the trend are randomly attributed to warm or cold
                per_ind[y == trend[q,,]]=0
                per_ind[per_ind==0]=sample(c(-1,1),1)
            }

            markov_per$ind[q,,] = per_ind
            # Go through vector and calculate persistent events 
            # Perform this step on a 1D vector to avoid artificial cutoffs 
            # at day 1 and day 365 of the year 
            per_ind1D = as.vector(per_ind) 

                
            tmp=seasonal_matrix_out(per_ind1D,model,states,array(c(60,150,151,241,242,333,334,424),dim=c(2,4)))
            for (i in 1:4){
                markov_per$markov[q,i,,]=tmp$out[i,,]
                markov_per$markov_conf[q,i,]=tmp$out_conf[i,]
            }

            tmp=seasonal_matrix_out(per_ind1D,model,states,array(c(1,365),dim=c(2,1)))
            markov_per$markov[q,5,,]=tmp$out[1,,]
            markov_per$markov_conf[q,5,]=tmp$out_conf[1,]

        } 
        else {
            cat(sprintf("> ID %s lon %s  lat %s <",dat$ID[q],dat$lon[q],dat$lat[q]))
        }
     
    }
    markov_write(filename,dat,markov_per,transitions,transition_names) 

    cat("done.\n")
    return(0)#list(markov_per=markov_per,shock_per=shock_per))
}

state_attribution <- function(dat,trend,nday,nyr,filename){
    source("functions_markov.r")

    ## User parameters 
    #trash is the number of data point which are wasted by detrending
    trash = ((nyr-1)/2*365+(nday-1))
    ntot = length(dat$ID)
    laenge_zeit = length(dat$time)
    transitions=states*states

    # Calculate persistence information
    #cat("Calculating persistence... ")

    markov_per = list(ind=dat$tas*NA,markov=array(NA,dim=c(ntot,5,transitions,65)),markov_conf=array(NA,dim=c(ntot,5,65)))

    for (q in 1:ntot) { 
        cat("-")
        if (length(which(is.na(dat$tas[q,,])))<(2*trash+365*20)){

            # Calculate persistence vector
            y = dat$tas[q,,]
            per_ind = y*NA 

            if (states==2){
                detrended = y-trend[q,,]
                threshold = median(detrended,na.rm=TRUE)
                per_ind[detrended < threshold]=-1
                per_ind[detrended > threshold]=1
            # the >= was somehow problematic, since it affects roughly 5% of the datapoints
            # now the datapoints sitting on the trend are randomly attributed to warm or cold
                per_ind[detrended == threshold]=1
                #per_ind[per_ind==0]=sample(c(-1,1),1)
            }

        } 
        else {
            cat(sprintf("> ID %s lon %s  lat %s <",dat$ID[q],dat$lon[q],dat$lat[q]))
        }
     
    }

    day <- dim.def.ncdf("day", units="d",vals=1:365, unlim=FALSE)
    year <- dim.def.ncdf("year",units="year",vals=1:65, unlim=FALSE)
    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    varstates <- dim.def.ncdf("states",units="cold=-1 , warm=1",vals=1:states, unlim=FALSE)


    ind <- var.def.ncdf(name="ind",units="1 or -1",dim=list(ID,day,year), missval=-9999.0)
    nc = create.ncdf(filename,ind)
    put.var.ncdf(nc,vars[[1]],per_ind)
    close.ncdf(nc)

    cat("done.\n")
    return(0)#list(markov_per=markov_per,shock_per=shock_per))
}

trend_analysis <- function(x,y){
    library(Kendall)
    if (length(which(is.na(y)))>7 | length(which(y==0))>40){
        return(list(slope=NA,slope_sig=NA,MK=NA,MK_sig=NA))
    }
    lm.r=lm(y~x)
    slope=summary(lm.r)$coefficients[2]
    slope_sig=summary(lm.r)$coefficients[8]
    out=MannKendall(y)
    MK=out[1]$tau
    MK_sig=out[2]$sl
    return(list(slope=slope,slope_sig=slope_sig,MK=MK,MK_sig=MK_sig))
}

global_analysis <- function(toAna,yearPeriod,yearshift=1949){
    #toANA of dim(ntot, something, year)
    yearPeriod=yearPeriod-yearshift
    t=seq(yearPeriod[1],yearPeriod[2],1)
    ntot=1319
    series=dim(toAna)[2]
    analysis=array(NA,dim=c(ntot,series,6))
    for (i in 1:series){
        for (q in 1:ntot){
            tmp=trend_analysis(t,toAna[q,i,yearPeriod[1]:yearPeriod[2]])
            analysis[q,i,1]=mean(toAna[q,i,yearPeriod[1]:yearPeriod[2]],na.rm=TRUE)
            analysis[q,i,2]=sd(toAna[q,i,yearPeriod[1]:yearPeriod[2]],na.rm=TRUE)
            analysis[q,i,3]=tmp$MK
            analysis[q,i,4]=tmp$MK_sig
            analysis[q,i,5]=tmp$slope
            analysis[q,i,6]=tmp$slope_sig
        }
    }
    return(analysis)
}

end_aussage <- function(dat,yearPeriod,trendID,states,seasons=c("spring","summer","autumn","winter","year"),region=c(-180,180,30,60)){
    # calculates percentage of grid points in region having positive trend
    # calculates big average
    print(yearPeriod)
    for (season in seasons){
        print(season)
        nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/duration/",yearPeriod,"/",trendID,"_duration_",states,"s_analysis_",season,".nc",sep=""))
        dur_ana_full=get.var.ncdf(nc,"dur_ana_full")
        #duration_analysis_write(paste("../data/",trendID,"/",states,"_states/duration/",yearPeriod,"/",trendID,"_duration_",states,"s_analysis_",season,".nc",sep=""),dur_ana_full,season)
        quantiles=get.var.ncdf(nc,"quantiles")
        outs=get.var.ncdf(nc,"outs")
        percentage=array(NA,dim=c(states,length(quantiles)))
        for (state in 1:states){
            for (quan in quantiles){
                slope=dur_ana_full[1:1319,state,quan,1]
                inside_region=which(dat$lon>region[1] & dat$lon<region[2] & dat$lat>region[3] & dat$lat< region[4])
                noNa=which(!is.na(slope[inside_region]))
                percentage[state,quan]=length(which(slope[inside_region[noNa]]>0))/length(noNa)
            }
        }
        print(percentage)
    }
}



