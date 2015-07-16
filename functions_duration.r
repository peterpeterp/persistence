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


calc_global_dur <- function(dat,per,trash,filename,states=c(-1,1)){
    ntot=length(dat$ID)
    

    dur=array(NA,dim=c(ntot,length(states)*2,62*365))
    maxis=array(NA,ntot)
    len=array(NA,length(states)*2)
    for (q in 1:ntot){
        cat("-")
        #per$ind[per$ind==0]=NA
        for (i in 1:length(states)){
            tmp=per_duration(as.vector(per$ind[q,,])[trash:(length(per$ind[q,,])-trash)],dat$time[trash:(length(per$ind[q,,])-trash)],states[i])
            dur[q,i,1:length(tmp$period)]=tmp$period
            dur[q,(length(states)+i),1:length(tmp$period)]=tmp$period_mid
        }

        for (i in 1:length(states)){
            len[i]=length(which(!is.na(dur[q,i,])))
        }
        maxis[q]=max(len,na.rm=TRUE)
    }


    duration_write(filename,dur[1:ntot,1:(length(states)*2),1:max(maxis,na.rm=TRUE)],max(maxis,na.rm=TRUE))
}


duration_seasons <- function(dur,season,states,filename){
    ntot=1319
    dur_neu=array(NA,dim=c(ntot,length(states)*2,65*92))
    maxis=array(NA,ntot)

    start=season[1]/365
    stop=season[2]/365

    len=array(NA,4)
    for (q in 1:ntot){
        for (i in 1:length(states)){
            select=c()
            for (year in 1950:2014){
                select=c(select,which(dur[q,((i-1)*2+2),]>(start+year) & dur[q,((i-1)*2+2),]<(stop+year)))
            }   
            if (length(select)>0){
                dur_neu[q,i,]=dur[q,i,select]
                dur_neu[q,(length(states)+i),]=dur[q,(length(states)+i),select]
            }        
        }

        for (i in 1:(length(states)+i)){
            len[i]=length(which(!is.na(dur_neu[q,i,])))
        }
        maxis[q]=max(len,na.rm=TRUE)
    }
    duration_write(filename,dur_neu[1:ntot,1:(2*length(states)),1:max(maxis,na.rm=TRUE)],max(maxis,na.rm=TRUE))
}
    
duration_analysis <- function(dur,filename,season,trenn=1980,stations=seq(1,1319,1)){
    br=seq(0,300,2)
    ntot=1319
    dur_ana=array(NA,dim=c(ntot,3,2,8))
    start=c(1950,1950,1980)
    stop=c(2011,1980,2011)


    for (q in c(488,238)){
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
                    print(duration[select])
                    print(sum(duration[select]))
                    print(q)

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


