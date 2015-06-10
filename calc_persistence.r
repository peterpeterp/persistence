
# Load useful functions 
#library(myr)                # Only needed for plotting: myfigure function
dyn.load("persistence_tools.so")
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


calc_per <- function(dat,trend,filename,nday,nyr){
    ## User parameters 
    trash = ((nyr-1)/2*365+(nday-1))
    ntot = length(dat$ID)
    laenge_zeit = length(dat$time)

    # Calculate persistence information
    #cat("Calculating persistence... ")

    #per = list(ind=dat$tas*NA,dayp=dat$tas*NA,perp=dat$tas*NA,
     #                  dayn=dat$tas*NA,pern=dat$tas*NA,year_warm=array(NA,dim=c(ntot,62)),year_cold=array(NA,dim=c(ntot,62)),year_warm_trend=array(NA,ntot),year_cold_trend=array(NA,ntot),
      #                 winter_warm=array(NA,dim=c(ntot,62)),winter_cold=array(NA,dim=c(ntot,62)),summer_warm=array(NA,dim=c(ntot,62)),summer_cold=array(NA,dim=c(ntot,62)),sum_warm_trend=array(NA,ntot),sum_cold_trend=array(NA,ntot),
       #                win_warm_trend=array(NA,ntot),win_cold_trend=array(NA,ntot))

    per = list(ind=dat$tas*NA,markov=array(NA,dim=c(ntot,6,62)),
        markov_mk=array(NA,dim=c(ntot,6)),markov_lr=array(NA,dim=c(ntot,6)),markov_mk_sig=array(NA,dim=c(ntot,6)),markov_lr_sig=array(NA,dim=c(ntot,6)),
        shock=array(NA,dim=c(ntot,6,62)),shock_mk=array(NA,dim=c(ntot,6)),shock_mk_sig=array(NA,dim=c(ntot,6)),shock_lr=array(NA,dim=c(ntot,6)),shock_lr_sig=array(NA,dim=c(ntot,6)),
        bic=array(NA,dim=c(ntot,3,62)))

    qq=which(dat$ID==209)

    for (q in 1:ntot ) { #   1:ntot
        cat("-")
        if (length(which(is.na(dat$tas[q,,])))<(2*trash+730)){
            tmp = calc_persistence(dat$tas[q,,],trend=trend[q,,],time=dat$time,nday=nday,nyr=nyr,trash=trash) 
            
            per$ind[q,,] =tmp$ind
            per$markov[q,,] = tmp$markov
            per$markov_mk[q,] = tmp$markov_mk
            per$markov_mk_sig[q,] = tmp$markov_mk_sig

            per$shock[q,,] = tmp$shock
            per$shock_mk[q,] = tmp$shock_mk
            per$shock_mk_sig[q,] = tmp$shock_mk_sig


            per$bic[q,,]=tmp$bic

            for (i in 1:6){
                lm.r=lm(tmp$markov[i,]~dat$year)
                per$markov_lr[q,i]=summary(lm.r)$coefficients[2]
                per$markov_lr_sig[q,i]=summary(lm.r)$coefficients[8]
            }

            for (i in 1:3){
                lm.r=lm(tmp$shock[i,]~dat$year)
                per$shock_lr[q,i]=summary(lm.r)$coefficients[2]
                per$shock_lr_sig[q,i]=summary(lm.r)$coefficients[8]
            }

        } 
        else {
            cat(sprintf("> ID %s lon %s  lat %s <",dat$ID[q],dat$lon[q],dat$lat[q]))
        }
 
    }
    per_write(filename,dat,per)

    cat("done.\n")
    return(per)
}



ndays = c(121,91,61,45)
nyrs = c(7,5,3,1)

ndays = c(91)
nyrs = c(5)




dat=dat_load("../data/mid_lat.nc")


for (nday in ndays){
    for (nyr in nyrs){
        cat(sprintf("\n%s_%s   ",nday,nyr))
        cat("calculating trend \n")
        #trend=calc_trend(dat,sprintf("../data/%s_%s_trend.nc",nday,nyr),nday,nyr)
        trend=trend_load(sprintf("../data/%s_%s_trend.nc",nday,nyr))
        cat(sprintf("\n%s_%s    ",nday,nyr))
        cat("calculating persistence\n")      
        per=calc_per(dat,trend,sprintf("../data/%s_%s_per_shock--------wow.nc",nday,nyr),nday,nyr)
    }
}




