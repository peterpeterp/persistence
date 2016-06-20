
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
    # calculation is done in c
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


state_attribution <- function(detrended,nday,nyr,filename){
    ## User parameters 
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
            #  with "per_ind[per_ind==0]=sample(c(-1,1),1)"" the datapoints sitting on the trend could be randomly attributed to warm or cold
            per_ind[detrended[q,,] == threshold]=1
            

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

master_nas <- function(){
    # this function is used to analyze the ratio and sampling of NAs
    plot_init_Had_multiple_noAA()

    margins<<-c(0,0,0,5)
    color_legend<<-"right"
    indexTopRight<<-c(NA)
    indexBottomLeft<<-c(NA)

    if (FALSE){
        # count nas and plot ratio
        naRatio=array(0,c(5,ntot))
        for (q in 1:ntot){naRatio[1,q]<-length(which(is.na(dat$tas[q,,])))/(365*length(dat$year))}
        naRatio[2,which(naRatio[1,]<0.3)]=1
        naRatio[3,which(naRatio[1,]<0.2)]=1
        naRatio[4,which(naRatio[1,]<0.1)]=1
        naRatio[5,which(naRatio[1,]<0.01)]=1
        naRatio[,2]=3

        write.table(naRatio,paste("../data/",dataset,"/naRatio.txt",sep=""))

        indexBottomLeft<<-c("a")
        topo_map_plot(filename=paste("../plots/",dataset,"/na_ratio.pdf",sep=""),reihen=array(naRatio[,],c(5,ntot)),farb_mitte=c(0,0.7),farb_palette="weiss-rot",titel=c(""))
    }

    if (FALSE){
        # analyze trends in na ratio
        naseries<-dat$tas*0
        naRatio=array(NA,c(2,ntot))
        x=1:(365*65)
        for (q in 1:ntot){
            naseries[q,,][which(is.na(dat$tas[q,,]))]=1
            y<-as.vector(naseries[q,,])
            naRatio[1:2,q]=summary(lm(y~x))$coef[c(2,8)]
        }
        indexBottomLeft<<-c("b")
        topo_map_plot(filename=paste("../plots/",dataset,"/na_increase.pdf",sep=""),reihen=array(naRatio[1,]*3650,c(1,ntot)),farb_mitte="0",farb_palette="lila-gruen")
    }

    if (FALSE){
        # calculate na persistence just as for warm and cold persistence
        calc_global_dur(ind=naseries,states=c(0,1),filename=paste("../data/",dataset,"/na_duration.nc",sep=""))
    }

    if (FALSE){
        if (FALSE){
            #count missing periods and detect trends in missing periods
            nc=open.nc(paste("../data/",dataset,"/na_duration.nc",sep=""))
            dur<<-var.get.nc(nc,"dur")
            dur_mid<<-var.get.nc(nc,"dur_mid")
            intersects<-array(0,c(ntot,65))
            intersect_increase<-array(0,c(3,ntot))
            x<-1:65
            for (q in 1:ntot){
                cat(q)
                for (yr in 1:65){
                    intersects[q,yr]=length(which(dur_mid[q,2,]>(1949+yr) & dur_mid[q,2,]<(1949+yr+1)))
                }
                if (sum(intersects[q,],na.rm=TRUE)>0){
                    #print(summary(lm(intersects[q,]~x))$coef[c(2,8)])
                    intersect_increase[2:3,q]=summary(lm(intersects[q,]~x))$coef[c(2,8)]
                    intersect_increase[1,q]=sum(intersects[q,],na.rm=TRUE)
                }
            }
            intersects<<-intersects
            intersect_increase<<-intersect_increase
        }

        indexBottomLeft<<-c("b")
        topo_map_plot(filename=paste("../plots/",dataset,"/na_intersect_increase.pdf",sep=""),reihen=array(intersect_increase[2,]*10,c(1,ntot)),farb_mitte="0",farb_palette="lila-gruen")
        indexBottomLeft<<-c("a")
        topo_map_plot(filename=paste("../plots/",dataset,"/na_intersect.pdf",sep=""),reihen=array(intersect_increase[1,],c(1,ntot)),farb_mitte=c(0,1500),farb_palette="weiss-rot")

    }

    # for the period after 1979
    # exactly the same as above
    if (FALSE){
        # count nas and plot ratio
        naRatio=array(NA,c(5,ntot))
        for (q in 1:ntot){naRatio[1,q]<-length(which(is.na(dat$tas[q,,30:62])))/(365*length(dat$year))}  
        indexBottomLeft<<-c("a")     
        topo_map_plot(filename=paste("../plots/",dataset,"/na_ratio_1979.pdf",sep=""),reihen=array(naRatio[1,],c(1,ntot)),farb_mitte=c(0,0.7),farb_palette="weiss-rot",titel=c(""))
    }

    if (FALSE){
        # analyze trends in na ratio
        naseries<-dat$tas*0
        naRatio=array(NA,c(2,ntot))
        x=1:(365*65)
        for (q in 1:ntot){
            naseries[q,,][which(is.na(dat$tas[q,,30:62]))]=1
            y<-as.vector(naseries[q,,])
            naRatio[1:2,q]=summary(lm(y~x))$coef[c(2,8)]
        }
        indexBottomLeft<<-c("b")
        topo_map_plot(filename=paste("../plots/",dataset,"/na_increase_1979.pdf",sep=""),reihen=array(naRatio[1,]*3650,c(1,ntot)),farb_mitte="0",farb_palette="lila-gruen")
    }

    if (FALSE){
        if (FALSE){
            #count missing periods and detect trends in missing periods
            nc=open.nc(paste("../data/",dataset,"/na_duration.nc",sep=""))
            dur<<-var.get.nc(nc,"dur")
            dur_mid<<-var.get.nc(nc,"dur_mid")
            intersects<-array(0,c(ntot,32))
            intersect_increase<-array(0,c(3,ntot))
            x<-1:32
            for (q in 1:ntot){
                print(q)
                for (yr in 1:32){
                    intersects[q,yr]=length(which(dur_mid[q,2,]>(1978+yr) & dur_mid[q,2,]<(1978+yr+1)))
                }
                if (sum(intersects[q,],na.rm=TRUE)>0){
                    print(length(intersects[q,]))
                    #print(summary(lm(intersects[q,]~x))$coef[c(2,8)])
                    intersect_increase[2:3,q]=summary(lm(intersects[q,]~x))$coef[c(2,8)]
                    intersect_increase[1,q]=sum(intersects[q,],na.rm=TRUE)
                }
                else{intersect_increase[,q]=0}
            }
            intersects<<-intersects
            intersect_increase<<-intersect_increase
        }

        indexBottomLeft<<-c("b")        
        topo_map_plot(filename=paste("../plots/",dataset,"/na_intersect_increase_1979.pdf",sep=""),reihen=array(intersect_increase[2,]*10,c(1,ntot)),farb_mitte="0",farb_palette="lila-gruen")
        indexBottomLeft<<-c("a")
        topo_map_plot(filename=paste("../plots/",dataset,"/na_intersect_1979.pdf",sep=""),reihen=array(intersect_increase[1,],c(1,ntot)),farb_mitte=c(0,1000),farb_palette="weiss-rot")

    }

    if (TRUE){
        # plot histogramm of missing period lengths
        nc=open.nc(paste("../data/",dataset,"/na_duration.nc",sep=""))
        dur<<-var.get.nc(nc,"dur")  
        pdf(paste("../plots/",dataset,"/na_pdf.pdf",sep=""),width=4,height=4)
        tmp=hist(dur[,2,],breaks=(0:25000+0.5),plot=FALSE)
        plot(tmp$mids,tmp$density,xlim=c(1,10000),ylim=c(0.00001,1),xlab="",ylab="",log="xy",pch=20)
        abline(v=365,col=rgb(0.4,0.4,0.4,1),lty=2)
        text(365,0.5,365,col=rgb(0.4,0.4,0.4,1))
        graphics.off()
    }

}
