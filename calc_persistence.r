#!/home/pepflei/R/bin/Rscript
# Load useful functions 
#dyn.load("persistence_tools.so")
source("functions_support.r")
source("functions_markov.r")
source("functions_shock.r")
source("functions_duration.r")
source("write.r")
source("load.r")


calc_trend <- function(dat,filename,nday,nyr){
    trash = ((nyr-1)/2*365+(nday-1))
    ntot = length(dat$ID)
    trend=dat$tas*NA
    for (q in 1:ntot) {
        temp = r_calc_runmean_2D(dat$tas[q,,],nday=nday,nyr=nyr)
        temp[1:trash]=NA
        temp[(length(dat$time)-trash):length(dat$time)]=NA
        trend[q,,]=temp
    }

    trend_write(filename,dat,trend)
    return(trend)
}

calc_per <- function(dat,trend,nday,nyr,what,filename_markov){
    ## User parameters 
    trash = ((nyr-1)/2*365+(nday-1))
    ntot = length(dat$ID)
    laenge_zeit = length(dat$time)

    # Calculate persistence information
    #cat("Calculating persistence... ")

    if ((what=="markov")){
        markov_per = list(ind=dat$tas*NA,markov=array(NA,dim=c(ntot,6,62)),markov_err=array(NA,dim=c(ntot,3,62)))

        for (q in 1:ntot ) { 
            cat("-")
            if (length(which(is.na(dat$tas[q,,])))<(2*trash+1000)){

                # Calculate persistence vector
                y = dat$tas[q,,]
                per_ind = y*0 


                per_ind[y >= trend[q,,]]  =  1
                per_ind[y <  trend[q,,]]  = -1 
                per_ind[is.na(per_ind)] = 0  # To avoid NA values

                markov_per$ind[q,,] = per_ind
                # Go through vector and calculate persistent events 
                # Perform this step on a 1D vector to avoid artificial cutoffs 
                # at day 1 and day 365 of the year 
                per_ind1D = as.vector(per_ind) 
                per_ind1D[per_ind1D==0]=NA

                size=length(per_ind1D)

                tmp=seasonal_1D_out(per_ind1D,array(c(151,242,334,425),dim=c(2,2)),markov_chft)
                for (i in 1:2){
                    markov_per$markov[q,i,]=tmp$out1[i,]
                    markov_per$markov[q,(i+3),]=tmp$out2[i,]
                    markov_per$markov_err[q,i,]=tmp$out_err[i,]
                }

                tmp=seasonal(per_ind1D,array(c(1,365),dim=c(2,1)),markov_chft)
                markov_per$markov[q,3,]=tmp$out1[1,]
                markov_per$markov[q,6,]=tmp$out2[1,]
                markov_per$markov_err[q,2,]=tmp$out_err[1,]

            } 
            else {
                cat(dat$ID[q])
            }
     
        }
        markov_write(filename_markov,dat,markov_per) 
    }

    cat("done.\n")
    return(0)#list(markov_per=markov_per,shock_per=shock_per))
}

calc_per_3states <- function(dat,trend,nday,nyr,what,filename_markov){
    ## User parameters 
    trash = ((nyr-1)/2*365+(nday-1))
    ntot = length(dat$ID)
    laenge_zeit = length(dat$time)

    # Calculate persistence information
    #cat("Calculating persistence... ")

    if ((what=="markov")){
        markov_per = list(ind=dat$tas*NA,markov=array(NA,dim=c(ntot,5,9,62)),markov_conf=array(NA,dim=c(ntot,5,62)))

        for (q in 1:ntot) { 
            cat("-")
            if (length(which(is.na(dat$tas[q,,])))<(2*trash+730)){

                # Calculate persistence vector
                y = dat$tas[q,,]
                per_ind = y*NA 

                detrended = y-trend[q,,]
                threshold = sd(detrended,na.rm=TRUE)*0.5

                per_ind[detrended>threshold]  =  1
                per_ind[detrended<(-threshold)]  = -1 
                per_ind[detrended<threshold & detrended>(-threshold)] = 0  # To avoid NA values

                markov_per$ind[q,,] = per_ind
                # Go through vector and calculate persistent events 
                # Perform this step on a 1D vector to avoid artificial cutoffs 
                # at day 1 and day 365 of the year 
                per_ind1D = as.vector(per_ind) 

                
                tmp=seasonal_matrix_out(per_ind1D,array(c(59,151,151,243,243,335,335,425),dim=c(2,4)))
                for (i in 1:4){
                    markov_per$markov[q,i,,]=tmp$out[i,,]
                    markov_per$markov_conf[q,i,]=tmp$out_conf[i,]
                }

                tmp=seasonal_matrix_out(per_ind1D,array(c(1,365),dim=c(2,1)))
                markov_per$markov[q,5,,]=tmp$out[1,,]
                markov_per$markov_conf[q,5,]=tmp$out_conf[1,]

            } 
            else {
                cat(sprintf("> ID %s lon %s  lat %s <",dat$ID[q],dat$lon[q],dat$lat[q]))
            }
     
        }
        markov_3states_write(filename_markov,dat,markov_per) 
    }

    cat("done.\n")
    return(0)#list(markov_per=markov_per,shock_per=shock_per))
}



if (1==2){
    ndays=c(91)
    nyrs=c(5)  

    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2011.nc")

    for (nday in ndays){
        for (nyr in nyrs){
            cat(sprintf("\n%s_%s   ",nday,nyr))
            cat("calculating trend \n")
            #trend=calc_trend(dat,sprintf("../data/%s_%s/%s_%s_trend_r.nc",nday,nyr,nday,nyr),nday,nyr)
            trend=trend_load(sprintf("../data/%s_%s/%s_%s_trend_r.nc",nday,nyr,nday,nyr))
            cat(sprintf("\n%s_%s    ",nday,nyr))
            cat("calculating persistence\n") 
            per=calc_per_3states(dat,trend,nday,nyr,"markov",sprintf("../data/%s_%s/%s_%s_markov3s.nc",nday,nyr,nday,nyr))
            #per=markov_load(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr))
        }
    }
}

if (1==2){
    nday=91
    nyr=5
    trash=((nyr-1)/2*365+(nday-1))
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2011.nc")
    per=markov_load(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr))
    calc_global_dur(dat=dat,per=per,trash=trash,filename=sprintf("../data/%s_%s/%s_%s_duration.nc",nday,nyr,nday,nyr))
    dur=duration_load(filename=sprintf("../data/%s_%s/%s_%s_duration.nc",nday,nyr,nday,nyr))
    
    duration_seasons(dur,season=c(59,151),filename=sprintf("../data/%s_%s/%s_%s_duration_spring.nc",nday,nyr,nday,nyr))
    duration_seasons(dur,season=c(151,243),filename=sprintf("../data/%s_%s/%s_%s_duration_summer.nc",nday,nyr,nday,nyr))
    duration_seasons(dur,season=c(243,335),filename=sprintf("../data/%s_%s/%s_%s_duration_autumn.nc",nday,nyr,nday,nyr))
    duration_seasons(dur,season=c(335,425),filename=sprintf("../data/%s_%s/%s_%s_duration_winter.nc",nday,nyr,nday,nyr))
}

if (1==2){
    nday=91
    nyr=5

    for (season in c("spring","summer","autumn","winter")){
        print(season)
        dur=duration_load(filename=paste("../data/",nday,"_",nyr,"/",nday,"_",nyr,"_duration_",season,".nc",sep=""))
        duration_analysis(dur,filename=paste("../data/",nday,"_",nyr,"/",nday,"_",nyr,"_duration_analysis_",season,".nc",sep=""),
            season=season,trenn=1980)#,stations=134:136)
    }
}

if (1==2){
    nday=91
    nyr=5
    tmp=global_trend(filename_markov=sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr),filename_markov_neu=sprintf("../data/%s_%s/%s_%s_markov_trend.nc",nday,nyr,nday,nyr))

}

