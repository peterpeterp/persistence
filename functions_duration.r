#!/home/pepflei/R/bin/Rscript
source("write.r")
source("load.r")

per_duration <- function(ind,time,state){
    act_state=ind[1]
    period=ind*NA
    period_mid=ind*NA    
    state_count=1
    period_count=1
    nas=0
    for (i in 2:length(ind)){
        if (is.na(ind[i]) | is.na(act_state)){
            nas=nas+1
            if (is.na(act_state)){
                act_state=ind[i]
            }
        }

        else{
            if (act_state==ind[i] & act_state==state){
                state_count=state_count+1
            } 
            if (act_state!=ind[i] & act_state==state){
                period[period_count]=state_count
                period_mid[period_count]=time[i]-0.5*state_count/365
                period_count=period_count+1
                state_count=1
                act_state=99
            }
            if (act_state!=ind[i] & ind[i]==state){
                act_state=ind[i]
            }
        }
    }
    if (nas>length(ind)/2){
        return(list(period=NA,period_mid=NA))
    }
    else{
        if (act_state==state){
            period[period_count]=state_count
        }
        return(list(period=period,period_mid=period_mid))
    }
}


calc_global_dur <- function(dat,ind,trash,filename,states=c(-1,1)){
    ntot=length(dat$ID)
    
    dur=array(NA,dim=c(ntot,length(states),65*365))
    dur_mid=array(NA,dim=c(ntot,length(states),65*365))

    maxis=array(NA,ntot)
    len=array(NA,length(states)*2)
    for (q in 1:ntot){
        cat("-")
        #per$ind[per$ind==0]=NA
        for (i in 1:length(states)){
            tmp=per_duration(as.vector(ind[q,,])[trash:(length(ind[q,,])-trash)],dat$time[trash:(length(ind[q,,])-trash)],states[i])
            dur[q,i,1:length(tmp$period)]=tmp$period
            dur_mid[q,i,1:length(tmp$period)]=tmp$period_mid
        }

        for (i in 1:length(states)){
            len[i]=length(which(!is.na(dur[q,i,])))
        }
        maxis[q]=max(len,na.rm=TRUE)
    }

    duration_write(filename,dur[1:ntot,1:length(states),1:max(maxis,na.rm=TRUE)],
        dur_mid[1:ntot,1:length(states),1:max(maxis,na.rm=TRUE)],max(maxis,na.rm=TRUE))
}


duration_seasons <- function(dur,dur_mid,season,filename){
    states=dim(dur)[2]
    ntot=1319
    dur_neu=array(NA,dim=c(ntot,states,65*92))
    dur_mid_neu=array(NA,dim=c(ntot,states,65*92))
    maxis=array(NA,ntot)

    start=season[1]/365
    stop=season[2]/365

    len=array(NA,4)
    for (q in 1:ntot){
        for (i in 1:states){
            select=c()
            for (year in 1950:2014){
                select=c(select,which(dur_mid[q,i,]>(start+year) & dur_mid[q,i,]<(stop+year)))
            }   
            if (length(select)>0){
                dur_neu[q,i,1:length(select)]=dur[q,i,select]
                dur_mid_neu[q,i,1:length(select)]=dur_mid[q,i,select]
            }        
        }

        for (i in 1:states){
            len[i]=length(which(!is.na(dur_neu[q,i,])))
        }
        maxis[q]=max(len,na.rm=TRUE)
    }
    duration_write(filename,dur_neu[1:ntot,1:states,1:max(maxis,na.rm=TRUE)],
        dur_mid_neu[1:ntot,1:states,1:max(maxis,na.rm=TRUE)],max(maxis,na.rm=TRUE))
}
    
duration_analysis <- function(dur,dur_mid,filename,season,trenn=1980,stations=seq(1,1319,1)){
    library(quantreg)
    ntot=1319
    states=dim(dur)[2]
    dur_ana=array(NA,dim=c(ntot,states,8,5))

    for (q in stations){
        cat("-")
        for (t in 1:states){
            if (length(which(!is.na(dur[q,t,])))>10){
                duration=dur[q,t,]
                mid=dur_mid[q,t,]

                dur_ana[q,t,1,1]=mean(duration,na.rm=TRUE)
                dur_ana[q,t,1,2]=sd(duration,na.rm=TRUE)

                F=ecdf(duration)
                taus=c(0.25,0.5,0.75,0.9,0.95,0.98,0.99)
                for (i in 1:length(taus)){
                    qu=summary(rq(duration~mid,taus[i]))$coefficients
                    dur_ana[q,t,(1+i),1:4]=qu[c(2,8,1,7)]
                    dur_ana[q,t,(1+i),5]=quantile(duration,probs=c(taus[i]),na.rm=TRUE)
                }
            }
            else {
                cat(q)
            }
        }
    } 

    if (filename!=FALSE){
        duration_analysis_write(filename,dur_ana,season,trenn)
    }
}


