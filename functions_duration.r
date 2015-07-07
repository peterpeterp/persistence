#!/home/pepflei/R/bin/Rscript
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
        if (is.na(ind[i]) | is.na(state)){
            nas=nas+1
            if (is.na(state)){
                state=ind[i]
            }
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
        return(list(dur_warm=warm_period,dur_cold=cold_period,dur_warm_mid=warm_period_mid,dur_cold_mid=cold_period_mid))
    }
}

seasonal_duration <- function(ind,time,seasons=array(c(59,151,151,243,243,335,335,425),dim=c(2,4)),interval=365){
    size=length(ind)
    x=seq(1, size, 1)
    i=1
    j=1
    out1=array(NA,dim=c(dim(seasons)[2],4,62,100))
    while ((i+seasons[2,4])<size){
        if ((is.na(ind[i+1])==FALSE) & (is.na(ind[i+interval])==FALSE)){
            for (sea in 1:length(seasons[1,])){
                ind_season=ind[(seasons[1,sea]+i+1):(seasons[2,sea]+i)]
                time_season=time[(seasons[1,sea]+i+1):(seasons[2,sea]+i)]
                tmp=per_duration(ind_season,time_season)
                out1[sea,1,j,1:length(tmp$dur_warm)]=tmp$dur_warm
                out1[sea,2,j,1:length(tmp$dur_cold)]=tmp$dur_cold
                out1[sea,3,j,1:length(tmp$dur_warm)]=tmp$dur_warm_mid
                out1[sea,3,j,1:length(tmp$dur_cold)]=tmp$dur_cold_mid
            }
        }
        j=j+1
        i=i+interval
    }
    return(out1)
} 

calc_global_dur <- function(dat,per,trash,filename){
    ntot=length(dat$ID)
    
    dur=array(NA,dim=c(ntot,4,62*365))
    maxis=array(NA,ntot)
    len=array(NA,4)
    for (q in 1:ntot){
        cat("-")
        per$ind[per$ind==0]=NA
        tmp=per_duration(as.vector(per$ind[q,,])[trash:(length(per$ind[q,,])-trash)],dat$time[trash:(length(per$ind[q,,])-trash)])
        dur[q,1,1:length(tmp$dur_warm)]=tmp$dur_warm
        dur[q,2,1:length(tmp$dur_cold)]=tmp$dur_cold
        dur[q,3,1:length(tmp$dur_warm)]=tmp$dur_warm_mid
        dur[q,4,1:length(tmp$dur_cold)]=tmp$dur_cold_mid

        for (i in 1:2){
            len[i]=length(which(is.na(dur[q,i,])==FALSE))
        }
        maxis[q]=max(len,na.rm=TRUE)
    }
    duration_write(filename,dur[1:ntot,1:4,1:max(maxis,na.rm=TRUE)],max(maxis,na.rm=TRUE))
}


duration_seasons <- function(dur,season,filename){
    br<-seq(0,100,10)
    ntot=1319
    start=season[1]/365
    stop=season[2]/365

    dur_neu=array(NA,dim=c(ntot,4,60*90))
    maxis=array(NA,ntot)
    len=array(NA,4)
    for (q in 1:ntot){
        select_warm=c()
        select_cold=c()
        for (year in 1950:2011){
            select_warm=c(select_warm,which(dur$dur_warm_mid[q,]>(start+year) & dur$dur_warm_mid[q,]<(stop+year)))
            select_cold=c(select_cold,which(dur$dur_cold_mid[q,]>(start+year) & dur$dur_cold_mid[q,]<(stop+year)))
        }

        if (length(select_warm)>0){
            dur_neu[q,1,1:length(select_warm)]=dur$dur_warm[q,select_warm]
            dur_neu[q,3,1:length(select_warm)]=dur$dur_warm_mid[q,select_warm]
        }
        if (length(select_cold)>0){
            dur_neu[q,2,1:length(select_cold)]=dur$dur_cold[q,select_cold]
            dur_neu[q,4,1:length(select_cold)]=dur$dur_cold_mid[q,select_cold]
        }
        for (i in 1:2){
            len[i]=length(which(is.na(dur_neu[q,i,])==FALSE))
        }
        maxis[q]=max(len,na.rm=TRUE)
    }
    duration_write(filename,dur_neu[1:ntot,1:4,1:max(maxis,na.rm=TRUE)],max(maxis,na.rm=TRUE))
}
    
duration_analysis <- function(dur,filename,season,trenn=1980,stations=seq(1,1319,1)){
    br=seq(0,200,2)
    ntot=1319
    dur_ana=array(NA,dim=c(ntot,3,2,8))
    start=c(1950,1950,1980)
    stop=c(2011,1980,2011)


    for (q in stations){
        cat("-")
        if (length(which(!is.na(dur$dur_warm[q,])))>50){
            for (i in 1:3){
                for (t in 1:2){
                    if (t==1){
                        mid=dur$dur_warm_mid[q,]
                        duration=dur$dur_warm[q,]
                    }
                    if (t==2){
                        mid=dur$dur_cold_mid[q,]
                        duration=dur$dur_cold[q,]
                    }
                    select=which(mid>start[i] & mid<stop[i])
                    histo=hist(duration[select],breaks=br,plot=FALSE)

                    x=histo$mids
                    xy=data.frame(y=histo$density,x=x)
                    fit=nls(y~exp(a+b*x),data=xy,start=list(a=0,b=0))
                    a=summary(fit)$parameters[1]
                    a_err=summary(fit)$parameters[3]
                    b=summary(fit)$parameters[2]
                    b_err=summary(fit)$parameters[4]

                    perc2=(log(0.02)-a)/b
                    perc5=(log(0.05)-a)/b
                    perc10=(log(0.10)-a)/b

                    dur_ana[q,i,t,1:8]=c(mean(duration[select],na.rm=TRUE),a,a_err,b,b_err,perc2,perc5,perc10)
                }
            }
        } 
        else {
            cat(q)
        }


    }
    if (filename!=FALSE){
        duration_analysis_write(filename,dur_ana,season,trenn)
    }
}


