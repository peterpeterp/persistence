source("write.r")
source("load.r")

per_duration <- function(ind,time){
    state=ind[1]
    warm_period=ind*NA
    cold_period=ind*NA    
    warm_period_mid=ind*NA
    cold_period_mid=ind*NA
    period=1
    j_warm=1
    j_cold=1
    nas=0
    for (i in 2:length(ind)){
        #cat(i,state,ind[i],period,"\n")
        if (is.na(ind[i])){
            nas=nas+1
        }
        else{
            if (state==ind[i]){
                period=period+1
            } else{
                if (state==1){
                    warm_period[j_warm]=period
                    warm_period_mid[j_warm]=time[i]-0.5*period/365
                    j_warm=j_warm+1              
                }
                if (state==-1){
                    cold_period[j_cold]=period
                    cold_period_mid[j_warm]=time[i]-0.5*period/365
                    j_cold=j_cold+1
                }
                period=1
                state=state-2*state
            } 
        }
    }
    if (nas>length(ind)/2){
        return(list(dur_warm=NA,dur_cold=NA,dur_warm_mid=NA,dur_cold_mid=NA))
    }
    else{
        if (state==1){
            warm_period[j_warm]=period
            j_warm=j_warm+1              
        }
        if (state==-1){
            cold_period[j_cold]=period
            j_cold=j_cold+1
        }
        return(list(dur_warm=warm_period[1:(j_warm-1)],dur_cold=cold_period[1:(j_cold-1)],dur_warm_mid=warm_period_mid[1:(j_warm-1)],dur_cold_mid=cold_period_mid[1:(j_cold-1)]))
    }
}

duration_hist <- function(dat,per,station,beg1,beg2,end1,end2){
    br<-seq(0,100,10)

    dur1=per_duration(ind=per$ind[station,151:243,3:31],time=dat$time_2D[151:243,beg1:end1])
    dur2=per_duration(ind=per$ind[station,151:243,30:58],time=dat$time_2D[151:243,beg2:end2])

    if (is.na(dur1)){
        cat(sprintf("> ID %s lon %s  lat %s <",dat$ID[station],dat$lon[station],dat$lat[station]))
        return(c(NA,NA,NA,NA,NA,NA))
    }
    else{

        dur12_warm=c(dur1$dur_warm,dur2$dur_warm)
        warm_mean=mean(dur12_warm)

        dur12_cold=c(dur1$dur_cold,dur2$dur_cold)
        cold_mean=mean(dur12_cold)

        hist1=hist(dur1$dur_warm,breaks=br,plot=FALSE)
        hist2=hist(dur2$dur_warm,breaks=br,plot=FALSE)
        warmX_diff=sum(hist2$density[3:8])-sum(hist1$density[3:8])

        hist1=hist(dur1$dur_cold,breaks=br,plot=FALSE)
        hist2=hist(dur2$dur_cold,breaks=br,plot=FALSE)   
        coldX_diff=sum(hist2$density[3:8])-sum(hist1$density[3:8])

        warm_mean_diff=mean(dur2$dur_warm)-mean(dur1$dur_warm)
        cold_mean_diff=mean(dur2$dur_cold)-mean(dur1$dur_cold)

        if (3==4){
            pdf(file="../plots/histogramm.pdf")
            plot(hist1,col=rgb(0,0,1,1/4))
            plot(hist2,col=rgb(1,0,0,1/4),add=TRUE)
        }  
        return(c(warm_mean,cold_mean,warm_mean_diff,cold_mean_diff,warmX_diff,coldX_diff))
    }
}

calc_global_dur <- function(dat,per,filename){
    ntot=length(dat$ID)

    dur=array(NA,dim=c(ntot,6))
    for (q in 1:ntot){
        cat("-")
        dur[q,]=duration_hist(dat,per,q,3,31,30,58)
    }
    duration_write(filename,dur)
}

nday=91
nyr=5
dat=dat_load("../data/mid_lat.nc")
per=markov_load(sprintf("../data/%s_%s/%s_%s_markov.nc",nday,nyr,nday,nyr))

calc_global_dur(dat,per,sprintf("../data/%s_%s/%s_%s_duration.nc",nday,nyr,nday,nyr))