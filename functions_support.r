
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
        for (j in 1:length(dat$year)){
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
        for (j in 1:length(dat$year)){
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
        for (j in 1:length(dat$year)){
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
 

calc_trend <- function(dat,filename,nday,nyr,procedure){
    # calculates running mean for each grid point
    # can choose between the c script and r function
    trash = ((nyr-1)/2*365+(nday-1)/2)
    ntot = length(dat$ID)
    trend=dat$tas*NA
    for (q in 1:ntot) {
        temp = procedure(dat$tas[q,,],nday=nday,nyr=nyr)
        temp[1:trash]=NA
        temp[(length(dat$time)-trash):length(dat$time)]=NA
        trend[q,,]=temp
    }
    trrr<<-trend
    trend_write(filename,trend)
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

state_attribution <- function(detrended,nday,nyr,filename){
    ## User parameters 
    #trash is the number of data point which are wasted by detrending
    ntot<-length(dat$ID)

    # Calculate persistence information
    #cat("Calculating persistence... ")

    state_ind<-dat$tas*NA

    for (q in 1:ntot) { 
        cat("-")
        if (length(which(!is.na(dat$tas[q,,])))>(365*20)){

            # Calculate persistence vector
            y <- detrended[q,,]
            per_ind <- y*NA 

            threshold <- 0
            per_ind[detrended[q,,] < threshold]=-1
            per_ind[detrended[q,,] > threshold]=1
            # the >= was somehow problematic, since it affects roughly 5% of the datapoints
            # now the datapoints sitting on the trend are randomly attributed to warm or cold
            per_ind[detrended[q,,] == threshold]=1
            #per_ind[per_ind==0]=sample(c(-1,1),1)

            state_ind[q,,]=per_ind

        } 
        else {
            cat(sprintf("> ID %s lon %s  lat %s <",dat$ID[q],dat$lon[q],dat$lat[q]))
        }
     
    }
    print(filename)
    nc_out<-create.nc(filename)

    print(dim(state_ind))

    att.put.nc(nc_out, "NC_GLOBAL", "method", "NC_CHAR", "detrended with 2d running mean. than for each season and each grid point the median is calculated. days with temp above (below) median are warm (cold).")
    att.put.nc(nc_out, "NC_GLOBAL", "explanation", "NC_CHAR", "state indices -1 (cold), 1 (warm)")
    
    dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)   
    dim.def.nc(nc_out,"day",dimlength=365, unlim=FALSE)   
    dim.def.nc(nc_out,"year",dimlength=length(dat$year), unlim=FALSE)   

    var.def.nc(nc_out,"ind","NC_DOUBLE",c(0,1,2))
    att.put.nc(nc_out, "ind", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "ind", "dim_explanation", "NC_CHAR", "ID-day-year")
    att.put.nc(nc_out, "ind", "explanation", "NC_CHAR", "state attribution for each day at each gridpoint")

    var.put.nc(nc_out,"ind",state_ind)             

    cat("done.\n")
}




