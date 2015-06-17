#!/home/pepflei/R/bin/Rscript
# Load useful functions 
#library(myr)                # Only needed for plotting: myfigure function
#dyn.load("persistence_tools.so")
source("functions_persistence.r")
source("write_load.r")


calc_trend <- function(dat,filename,nday,nyr){
    trash = ((nyr-1)/2*365+(nday-1))
    ntot = length(dat$ID)

    trend=dat$tas*NA
    for (q in 1:ntot) {
        temp = c_calc_runmean_2D(dat$tas[q,,],nday=nday,nyr=nyr)
        temp[1:trash]=NA
        temp[(length(dat$time)-trash):length(dat$time)]=NA
        trend[q,,]=temp
    }

    trend_write(filename,dat,trend)
    return(trend)
}


calc_per <- function(dat,trend,filename_markov,filename_shock,nday,nyr,what){
    ## User parameters 
    trash = ((nyr-1)/2*365+(nday-1))
    ntot = length(dat$ID)
    laenge_zeit = length(dat$time)

    # Calculate persistence information
    #cat("Calculating persistence... ")

    if ((what=="markov") | (what=="both")){
        markov_per = list(ind=dat$tas*NA,markov=array(NA,dim=c(ntot,6,62)),markov_err=array(NA,dim=c(ntot,3,62)))

        for (q in 1:ntot ) { 
            cat("-")
            if (length(which(is.na(dat$tas[q,,])))<(2*trash+730)){

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

                tmp=seasonal(per_ind1D,array(c(151,242,334,425),dim=c(2,2)),markov_chft)
                markov_per$markov[q,1,]=tmp$summer_w
                markov_per$markov[q,2,]=tmp$summer_c
                markov_per$markov[q,3,]=tmp$winter_w
                markov_per$markov[q,4,]=tmp$winter_c
                markov_per$markov_err[q,1,]=tmp$error_s
                markov_per$markov_err[q,2,]=tmp$error_w

                tmp=seasonal(per_ind1D,array(c(1,365),dim=c(2,1)),markov_chft)
                markov_per$markov[q,5,]=tmp$summer_w
                markov_per$markov[q,6,]=tmp$summer_c
                markov_per$markov_err[q,3,]=tmp$error_s
            } 
            else {
                cat(sprintf("> ID %s lon %s  lat %s <",dat$ID[q],dat$lon[q],dat$lat[q]))
            }
     
        }
        markov_write(filename_markov,dat,markov_per) 
    }

    if ((what=="shock") | (what=="both")){
        shock_per = list(shock=array(NA,dim=c(ntot,6,62)),shock_bic=array(NA,dim=c(ntot,3,62)))

        for (q in 1:ntot ) { 
            cat("-")
            if (length(which(is.na(dat$tas[q,,])))<(2*trash+730)){

                tmp=seasonal((dat$tas[q,,]-trend[q,,]),array(c(151,242,334,425),dim=c(2,2)),shock_ma,3)
                shock_per$shock[q,1,]=tmp$summer_w
                shock_per$shock[q,2,]=tmp$winter_w
                shock_per$shock_bic[q,1,]=tmp$bic_s
                shock_per$shock_bic[q,2,]=tmp$bic_s

                tmp=seasonal((dat$tas[q,,]-trend[q,,]),array(c(1,365),dim=c(2,1)),shock_ma,3)
                shock_per$shock[q,3,]=tmp$summer_w
                shock_per$shock_bic[q,3,]=tmp$bic_s
            } 
            else {
                cat(sprintf("> ID %s lon %s  lat %s <",dat$ID[q],dat$lon[q],dat$lat[q]))
            }
     
        }
        shock_write(filename_shock,dat,shock_per)       
    }

    cat("done.\n")
    return(0)#list(markov_per=markov_per,shock_per=shock_per))
}

trend_analysis <- function(x,y){
    library(Kendall)
    if (is.na(min(y))){
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

global_trend <- function(filename_markov=99,filename_markov_neu=99,filename_shock=99,filename_shock_neu=99){
    t=seq(1,62,1)
    ntot=819
    trend_markov=array(NA,dim=c(ntot,6,4))
    if (filename_markov!=99){   
        per=markov_load(filename_markov)
        for (i in 1:6){
            for (q in 1:ntot){
                tmp=trend_analysis(t,per$markov[q,i,])
                trend_markov[q,i,1]=tmp$slope
                trend_markov[q,i,2]=tmp$slope_sig
                trend_markov[q,i,3]=tmp$MK
                trend_markov[q,i,4]=tmp$MK_sig
            }
        }
        markov_trend_write(filename_markov_neu,trend_markov)
    }
    trend_shock=array(NA,dim=c(ntot,3,4))
    if (filename_shock!=99){      
        per=shock_load(filename_shock)
        for (i in 1:3){
            for (q in 1:ntot){
                tmp=trend_analysis(t,per$shock[q,i,])
                trend_shock[q,i,1]=tmp$slope
                trend_shock[q,i,2]=tmp$slope_sig
                trend_shock[q,i,3]=tmp$MK
                trend_shock[q,i,4]=tmp$MK_sig
            }
        }
        shock_trend_write(filename_shock_neu,trend_shock)
    }
    return(list(trend_markov=trend_markov,trend_shock=trend_shock))
}

ndays = c(91,121,61)
nyrs = c(5,7,3)


dat=dat_load("../data/mid_lat.nc")


for (nday in ndays){
    for (nyr in nyrs){
        cat(sprintf("\n%s_%s   ",nday,nyr))
        cat("calculating trend \n")
        trend=calc_trend(dat,sprintf("../data/%s_%s_trend.nc",nday,nyr),nday,nyr)
        #trend=trend_load(sprintf("../data/%s_%s_trend.nc",nday,nyr))
        cat(sprintf("\n%s_%s    ",nday,nyr))
        cat("calculating persistence\n")      
        per=calc_per(dat,trend,sprintf("../data/%s_%s_markov.nc",nday,nyr),sprintf("../data/%s_%s_shock_ma_3_.nc",nday,nyr),nday,nyr,"both")
        tmp=global_trend(filename_markov=sprintf("../data/%s_%s_markov.nc",nday,nyr),filename_markov_neu=sprintf("../data/%s_%s_markov_trend.nc",nday,nyr),
                filename_shock=sprintf("../data/%s_%s_shock_ma_3_.nc",nday,nyr),filename_shock_neu=sprintf("../data/%s_%s_shock_ma_3_trend.nc",nday,nyr))
    }
}




