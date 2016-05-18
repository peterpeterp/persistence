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
    
    if (length(which(is.na(ind)))>length(ind)/2){
        return(list(period=NA,period_mid=NA))
    }
    else{
        ind[which(is.na(ind))]=-99
        act_state=ind[1]
        period=ind*NA
        period_mid=ind*NA    
        state_count=1
        period_count=1
        nas=0
        for (i in 2:length(ind)){
            if (act_state!=ind[i] & act_state==state){
                period[period_count]=state_count
                period_mid[period_count]=time[i]-0.5*state_count/365
                period_count=period_count+1
                state_count=1
                act_state=99
            }
            if (act_state==ind[i] & act_state==state){
                state_count=state_count+1
            } 
            if (act_state!=ind[i] & ind[i]==state){
                act_state=ind[i]
            }
        }
        return(list(period=period,period_mid=period_mid))
    }
}


calc_global_dur <- function(ind,filename,states=c(-1,1),years=length(dat$year)){
    # dat: is required for time and ID
    # ind: array of states (same dimentions as temp anomalies)
    # trash: number of days that have to be neglected (due to detrending)
    # filename: where to store result
    # states: label of states
    
    dur<-array(NA,dim=c(ntot,length(states),years*365))
    dur_mid<-array(NA,dim=c(ntot,length(states),years*365))

    maxis<-array(NA,ntot)
    len<-array(NA,length(states)*2)

    percentage<-0
    cat(paste("\n0 -> -> -> -> -> 100\n"))
    for (q in 1:ntot){
        if (q/ntot*100 > percentage){
            cat("-")
            percentage<-percentage+5
        }
        for (state in 1:length(states)){
            tmp<-per_duration(as.vector(ind[q,,]),dat$time,states[state])
            dur[q,state,1:length(tmp$period)]=tmp$period
            dur_mid[q,state,1:length(tmp$period)]=tmp$period_mid
        }

        for (state in 1:length(states)){
            len[state]=length(which(!is.na(dur[q,state,])))
        }
        maxis[q]=max(len,na.rm=TRUE)
    }
    duration_write(filename,dur[1:ntot,1:length(states),1:max(maxis,na.rm=TRUE)],dur_mid[1:ntot,1:length(states),1:max(maxis,na.rm=TRUE)],max(maxis,na.rm=TRUE))
}

duration_seasons <- function(dur,dur_mid,season,filename,years=length(dat$year)){
    # selects periods with midpoint in season
    # creates new duration files
    states=dim(dur)[2]
    ntot=length(dat$ID)
    dur_neu=array(NA,dim=c(ntot,states,years*92))
    dur_mid_neu=array(NA,dim=c(ntot,states,years*92))
    maxis=array(NA,ntot)

    start=season[1]/365
    stop=season[2]/365

    len=array(NA,4)
    for (q in 1:ntot){
        for (state in 1:states){
            select=c()
            for (year in dat$year){
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
    duration_write(filename=filename,dur=dur_neu[1:ntot,1:states,1:max(maxis,na.rm=TRUE)],dur_mid=dur_mid_neu[1:ntot,1:states,1:max(maxis,na.rm=TRUE)],len=max(maxis,na.rm=TRUE))
    cat(paste("\ndays in year:",season[1],"-",season[2]))    
}
    

duration_analysis <- function(yearPeriod,season_auswahl=c(1,2,3,4,5),option=c(1,0,0,0,0,0,0),ID_select=1:length(dat$ID),write=TRUE,add_name="quant_other",folder="/gridded/",ID_name="",plot_select=c(NA),ID_names=1:length(dat$ID),ID_length=length(ID_select),noise_level=0,xStart=1,xStop=100){    

    if (!is.na(plot_select[1])){
        pdf(file=paste("../plots/",dataset,additional_style,"/",trendID,folder,period,"/",ID_name,"_dist_diff_fit_plot_",dataset,"_",period,"_",add_name,".pdf",sep=""),width=3,height=3)
        par(mfrow=c(1,1))
        fit_plot_empty()
    }

    quantile_stuff=array(NA,dim=c(length(season_names),ID_length,2,length(taus),3))
    fit_stuff=array(NA,dim=c(length(season_names),ID_length,2,30))
    other_stuff=array(NA,dim=c(length(season_names),ID_length,2,12))
    distr_stuff=array(NA,dim=c(length(season_names),ID_length,2,10,100))

    for (sea in season_auswahl){ 
        season=season_names[sea]
        cat(paste("\n",season))  
        dists=list()

        print(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,ID_name,"_duration_",season,".nc",sep=""))
        nc_dur=open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,ID_name,"_duration_",season,".nc",sep=""))
        dur=var.get.nc(nc_dur,"dur")
        dur_mid=var.get.nc(nc_dur,"dur_mid")

        percentage<-0
        cat(paste("\n0 -> -> -> -> -> 100\n"))
        for (q in ID_select){
            if (q/ID_length*100 > percentage){
                cat("-")
                percentage<-percentage+5
            }
            for (state in 1:2){
                duration<<-dur[q,state,]
                duration_mid<<-dur_mid[q,state,]
                ord<-order(duration_mid)
                if (length(which(!is.na(duration)))>100){
                    y<<-as.vector(duration[ord])
                    x<<-as.vector(duration_mid[ord])
                    inYearPeriod<-which(x>yearPeriod[1] & x<yearPeriod[2])
                    y<<-y[inYearPeriod]
                    x<<-x[inYearPeriod]

                    # other stuff
                    if (option[1]==1){
                        other_stuff[sea,q,state,1]=mean(y,na.rm=TRUE)
                        other_stuff[sea,q,state,2]=sd(y,na.rm=TRUE)
                        other_stuff[sea,q,state,12]=length(!is.na(y))

                        linear<-try(lm(y~x),silent=TRUE)
                        if (class(linear)!="try-error"){other_stuff[sea,q,state,3:10]=summary(linear)$coef}
                    }

                    # quantile values and regressions for grid points
                    if (option[2]==1){
                        tmp<-quantile_analysis(x,y,taus,noise_level=noise_level)
                        quantile_stuff[sea,q,state,,1]=tmp$quantiles
                        quantile_stuff[sea,q,state,,2]=tmp$slopes
                        quantile_stuff[sea,q,state,,3]=tmp$slope_sigs
                    }

                    # quantile values and regressions for regions
                    if (option[3]==1){
                        quantile_stuff[sea,q,state,,1]=quantile_pete(y,taus=taus,na.rm=TRUE)
                        quantile_stuff[sea,q,state,,2]=rq(dither(y,value=noise_level)~x,taus)$coef[2,]
                    }

                    # data to be fitted
                    # wie sollte breaks definiert sein ....
                    br<-seq(0.5,max(y,na.rm=TRUE)+0.5,1)
                    histo<-hist(y,breaks=br,plot=FALSE)
                        
                    counts<-histo$counts
                    Y<-histo$density
                    X<-histo$mids

                    # ignore really long periods!!! really????
                    if (length(X)>100){stop<-100}
                    else {stop<-length(X)}

                    distr_stuff[sea,q,state,1,1:stop]=X[1:stop]
                    distr_stuff[sea,q,state,2,1:stop]=Y[1:stop]
                    distr_stuff[sea,q,state,3,1:stop]=counts[1:stop]



                    fit_stuff[sea,q,state,26]=length(which(histo$counts>0))
            
                    if (length(which(!is.na(Y)))>15){ 
                        if (option[4]==1){
                            # exponential fit as starting point
                            tmp_exp=exponential_fit(X,Y,xStart=xStart,plot_cdf=(q %in% plot_select))
                            fit_stuff[sea,q,state,1:2]=tmp_exp$pars
                            fit_stuff[sea,q,state,4:9]=tmp_exp$ana
                            distr_stuff[sea,q,state,4,1:stop]=tmp_exp$fit[1:stop]
                            distr_stuff[sea,q,state,5,xStart:stop]=tmp_exp$cdf_Data[1:(stop-xStart+1)]
                            distr_stuff[sea,q,state,6,xStart:stop]=tmp_exp$cdf_Fit[1:(stop-xStart+1)]

                            # combination of 2 exponentials seperated by threshold (restricted threshold range)
                        
                            tmp=two_exp_fit(X,Y,y,xStart=xStart,xStop=100,plot_cdf=(q %in% plot_select))
                            fit_stuff[sea,q,state,11:15]=tmp$pars
                            fit_stuff[sea,q,state,17:22]=tmp$ana
                            distr_stuff[sea,q,state,8,1:stop]=tmp$fit[1:stop]
                            distr_stuff[sea,q,state,9,xStart:stop]=tmp$cdf_Data[1:(stop-xStart+1)]
                            distr_stuff[sea,q,state,10,xStart:stop]=tmp$cdf_Fit[1:(stop-xStart+1)]

                            distr_stuff[sea,q,state,7,xStart:stop]=xStart:stop
                        }
                        fit_stuff[sea,q,state,24]=fit_stuff[sea,q,state,22]-fit_stuff[sea,q,state,9]
                    }

                }
                if (!is.na(plot_select[1])){
                    if (length(which(!is.na(duration)))<=100){
                        x=(1:10)*NA
                        y=(1:10)*NA
                        X=(1:10)*NA
                        Y=(1:10)*NA
                    }
                    if (q %in% plot_select){
                        if (is.na(fit_stuff[sea,q,state,1])){
                            expfit=X*NA
                            fit=X*NA
                        }
                        fit_plot_reference(x=x,y=y,sea=season_names[sea],q=ID_names[q],state=state)
                        fit_plot_combi(distrs=distr_stuff[sea,q,state,,],fitstuff=fit_stuff[sea,q,state,],sea=season_names[sea],q=ID_names[q],state=state)
                    }
                }
            }
        }
    }
    if (option[1]==1){other_write(filename=paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,ID_name,"_",period,"_others.nc",sep=""),ID_length=ID_length,ID_name="grid_points",period=period,other_stuff=other_stuff)}

    if (option[2]==1 | option[3]==1){quantiles_write(filename=paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,ID_name,"_",period,"_quantiles.nc",sep=""),ID_length=ID_length,ID_name="grid_points",period=period,quantile_stuff=quantile_stuff)}
        
    if (sum(option[4:8],na.rm=TRUE)>0){
        fit_write(filename=paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,ID_name,"_",period,"_fit_",add_name,".nc",sep=""),ID_length=ID_length,ID_name="grid_points",period=period,fit_stuff=fit_stuff)
        distr_write(distr_stuff=distr_stuff,filename=paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,ID_name,"_",period,"_distributions.nc",sep=""),ID_length=ID_length,ID_name="grid_points",period=period)
    }
    graphics.off()
}



