#!/home/pepflei/R/bin/Rscript

###################################################################
# calculate and analyse duration time series
# Author: Peter Pfleiderer
# Institution: Potsdam Institute for Climate Impact Research (PIK)
# Year: 2015
###################################################################

per_duration <- function(ind,time,state){
    # ind: is time serie of states for one grid point
    # finds periods of duration for given state
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
    # dat: is required for time and ID
    # ind: array of states (same dimentions as temp anomalies)
    # trash: number of days that have to be neglected (due to detrending)
    # filename: where to store result
    # states: label of states
    ntot=length(dat$ID)
    
    dur=array(NA,dim=c(ntot,length(states),65*365))
    dur_mid=array(NA,dim=c(ntot,length(states),65*365))

    maxis=array(NA,ntot)
    len=array(NA,length(states)*2)

    percentage=0
    cat(paste("\n0 -> -> -> -> -> 100\n"))
    for (q in 1:ntot){
        if (q/1319*100 > percentage){
            cat("-")
            percentage=percentage+5
        }
        for (state in 1:length(states)){
            tmp=per_duration(as.vector(ind[q,,])[trash:(length(ind[q,,])-trash)],dat$time[trash:(length(ind[q,,])-trash)],states[state])
            dur[q,state,1:length(tmp$period)]=tmp$period
            dur_mid[q,state,1:length(tmp$period)]=tmp$period_mid
        }

        for (state in 1:length(states)){
            len[state]=length(which(!is.na(dur[q,state,])))
        }
        maxis[q]=max(len,na.rm=TRUE)
    }
    duration_write(filename,dur[1:ntot,1:length(states),1:max(maxis,na.rm=TRUE)],
        dur_mid[1:ntot,1:length(states),1:max(maxis,na.rm=TRUE)],max(maxis,na.rm=TRUE))
}


duration_seasons <- function(dur,dur_mid,season,filename){
    # selects periods with midpoint in season
    # creates new duration files
    states=dim(dur)[2]
    ntot=1319
    dur_neu=array(NA,dim=c(ntot,states,65*92))
    dur_mid_neu=array(NA,dim=c(ntot,states,65*92))
    maxis=array(NA,ntot)

    start=season[1]/365
    stop=season[2]/365

    len=array(NA,4)
    for (q in 1:ntot){
        for (state in 1:states){
            select=c()
            for (year in 1950:2014){
                select=c(select,which(dur_mid[q,state,]>(start+year) & dur_mid[q,state,]<(stop+year)))
            }   
            if (length(select)>0){
                dur_neu[q,state,1:length(select)]=dur[q,state,select]
                dur_mid_neu[q,state,1:length(select)]=dur_mid[q,state,select]
            }        
        }

        for (state in 1:states){
            len[state]=length(which(!is.na(dur_neu[q,state,])))
        }
        maxis[q]=max(len,na.rm=TRUE)
    }
    duration_write(filename,dur_neu[1:ntot,1:states,1:max(maxis,na.rm=TRUE)],
        dur_mid_neu[1:ntot,1:states,1:max(maxis,na.rm=TRUE)],max(maxis,na.rm=TRUE))
    cat(paste("\ndays in year:",season[1],"-",season[2]))    

}
    

duration_analysis <- function(yearPeriod,trendID,dataset="_TMean",season_auswahl=c(1,2,3,4,5),option=c(1,0,0,0,0,0,0),ID_select=1:1319,write=TRUE,add_name="quant_other",folder="/gridded/",ID_name="",plot_select=c(488,1232,52,661),ID_names=1:1319,ID_length=length(ID_select),noise_level=c(0,0)){

    

    period=paste(yearPeriod[1],"-",yearPeriod[2],sep="")
    season_names=c("MAM","JJA","SON","DJF","4seasons")
    taus=c(0.05,0.25,0.5,0.75,0.95,0.98,1)

    if (!is.na(plot_select[1])){
        #pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,folder,ID_name,"_dist_diff_fit_plot_",dataset,"_",yearPeriod[1],"-",yearPeriod[2],"_",add_name,".pdf",sep=""),width=8,height=6)
        #par(mfrow=c(4,5))

        pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,folder,ID_name,"_dist_diff_fit_plot_",dataset,"_",yearPeriod[1],"-",yearPeriod[2],"_",add_name,".pdf",sep=""),width=2.2,height=2.2)
        par(mfrow=c(1,1))
        fit_plot_empty()
    }


    quantile_stuff=array(NA,dim=c(length(season_names),ID_length,2,length(taus),3))
    fit_stuff=array(NA,dim=c(length(season_names),ID_length,2,20))
    other_stuff=array(NA,dim=c(length(season_names),ID_length,2,12))

    for (sea in season_auswahl){ 
        season=season_names[sea]
        cat(paste("\n",season))  
        dists=list()

        nc_dur=open.nc(paste("../data/",trendID,"/",dataset,additional_style,folder,trendID,dataset,ID_name,"_duration_",season,".nc",sep=""))
        dur=var.get.nc(nc_dur,"dur")
        dur_mid=var.get.nc(nc_dur,"dur_mid")

        percentage=0
        cat(paste("\n0 -> -> -> -> -> 100\n"))
        for (q in ID_select){
            if (q/ID_length*100 > percentage){
                cat("-")
                percentage=percentage+5
            }
            for (state in 1:2){
                #cat(q)
                duration=dur[q,state,]
                duration_mid=dur_mid[q,state,]
                ord=order(duration_mid)
                if (length(which(!is.na(duration)))>100){
                    y=as.vector(duration[ord])
                    x=as.vector(duration_mid[ord])
                    inYearPeriod=which(x>yearPeriod[1] & x<yearPeriod[2])
                    y=y[inYearPeriod]
                    x=x[inYearPeriod]

                    # other stuff
                    if (option[1]==1){
                        other_stuff[sea,q,state,1]=mean(y,na.rm=TRUE)
                        other_stuff[sea,q,state,2]=sd(y,na.rm=TRUE)
                        other_stuff[sea,q,state,12]=length(!is.na(y))

                        linear=try(lm(y~x),silent=TRUE)
                        if (class(linear)!="try-error"){other_stuff[sea,q,state,3:10]=summary(linear)$coef}
                    }

                    # quantile values and regressions
                    if (option[2]==1){
                        tmp=quantile_analysis(x,y,taus,noise_level=noise_level)
                        quantile_stuff[sea,q,state,,1]=tmp$quantiles
                        quantile_stuff[sea,q,state,,2]=tmp$slopes
                        quantile_stuff[sea,q,state,,3]=tmp$slope_sigs
                    }

                    # data to be fitted
                    if (1==1){
                        br=seq(0,max(y,na.rm=TRUE),1)
                        histo=hist(y,breaks=br,plot=FALSE)
                        
                        Y=histo$density
                        X=histo$mids
                    }

                    if (!is.na(plot_select[1])){fit_plot_reference(x=x,y=y,sea=season_names[sea],q=ID_names[q],state=state)}

                    # exponential fit + other values
                    if (option[3]==1){
                        tmp=exponential_fit(X,Y)
                        fit_stuff[sea,q,state,1:2]=tmp$pars
                        fit_stuff[sea,q,state,19:20]=tmp$ana
                        expifit=tmp$fit
                        if (q %in% plot_select){fit_plot(X=X,Y=Y,fit=expifit,legend=c(paste("expo","\n","R2=",round(fit_stuff[sea,q,state,19],03),"\nBIC=",round(fit_stuff[sea,q,state,20],01),sep=""),paste("b1=",round(fit_stuff[sea,q,state,2],03),"\n ",sep="")),sea=season_names[sea],q=ID_names[q],state=state)}
                    }

                    # combination of 2 exponentials seperated by threshold restricted 
                    if (option[4]==1){
                        if (length(which(!is.na(Y)))>30){
                            tmp=two_exp_fit_restricted(X,Y,y)
                            fit_stuff[sea,q,state,1:5]=tmp$pars
                            fit_stuff[sea,q,state,19:20]=tmp$ana
                            fit=tmp$fit
                            if (q %in% plot_select){fit_plot(X=X,Y=Y,fit=fit,legend=c(paste("2expo_restrict","\n","R2=",round(fit_stuff[sea,q,state,19],03),"\nBIC=",round(fit_stuff[sea,q,state,20],01),sep=""),paste("b1=",round(fit_stuff[sea,q,state,2],03),"\nb2=",round(fit_stuff[sea,q,state,4],03),sep="")),sea=season_names[sea],q=ID_names[q],state=state,thresh=fit_stuff[sea,q,state,5])}
                        }
                    }

                    # combination of 2 exponentials seperated by threshold 
                    if (option[5]==1){
                        if (length(which(!is.na(Y)))>30){
                            tmp=two_exp_fit_gepfuscht(X,Y,y)
                            fit_stuff[sea,q,state,1:5]=tmp$pars
                            fit_stuff[sea,q,state,19:20]=tmp$ana
                            fit=tmp$fit
                            if (q %in% plot_select){fit_plot(X=X,Y=Y,fit=fit,legend=c(paste("2expo","\n","R2=",round(fit_stuff[sea,q,state,19],03),"\nBIC=",round(fit_stuff[sea,q,state,20],01),sep=""),paste("b1=",round(fit_stuff[sea,q,state,2],03),"\nb2=",round(fit_stuff[sea,q,state,4],03),sep="")),sea=season_names[sea],q=ID_names[q],state=state,thresh=fit_stuff[sea,q,state,5])}
                        }
                    }

                    # combination of 2 exponentials seperated by fixed threshold
                    if (option[6]==1){
                        if (length(which(!is.na(Y)))>30){
                            tmp=two_exp_fit_fixed_thresh(X,Y)
                            fit_stuff[sea,q,state,1:5]=tmp$pars
                            fit_stuff[sea,q,state,19:20]=tmp$ana
                            fit=tmp$fit
                            if (q %in% plot_select){fit_plot(X=X,Y=Y,fit=fit,legend=c(paste("2expo_fixed","\n","R2=",round(fit_stuff[sea,q,state,19],03),"\nBIC=",round(fit_stuff[sea,q,state,20],01),sep=""),paste("b1=",round(fit_stuff[sea,q,state,2],03),"\nb2=",round(fit_stuff[sea,q,state,4],03),sep="")),sea=season_names[sea],q=ID_names[q],state=state,thresh=fit_stuff[sea,q,state,5])}
                        }
                    }
                }
            }
        }
        if (option[1]==1){other_write(filename=paste("../data/",trendID,"/",dataset,additional_style,folder,period,"/",trendID,"_",dataset,ID_name,"_",period,"_others.nc",sep=""),ID_length=ID_length,ID_name="grid_points",period=period,other_stuff=other_stuff)}

        if (option[2]==1){quantiles_write(filename=paste("../data/",trendID,"/",dataset,additional_style,folder,period,"/",trendID,"_",dataset,ID_name,"_",period,"_quantiles",add_name,".nc",sep=""),ID_length=ID_length,ID_name="grid_points",period=period,taus=taus,quantile_stuff=quantile_stuff)}
        
        if (option[3]==1 | option[4]==1 | option[5]==1 | option[6]==1){fit_write(filename=paste("../data/",trendID,"/",dataset,additional_style,folder,period,"/",trendID,"_",dataset,ID_name,"_",period,"_fit_",add_name,".nc",sep=""),ID_length=ID_length,ID_name="grid_points",period=period,fit_stuff=fit_stuff)}

    }
    graphics.off()
}

