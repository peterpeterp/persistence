
ks_wilcoxon_trenID_sensitivity <- function(trendIDs=c("91_5","91_7","91_9")){
    wilcox_test=array(NA,c(2,3))
    ks_test=array(NA,c(2,3))

    pdf(paste("../plots/",dataset,additional_style,"/dur_ks_wilcox_sensibility_no_ecdf.pdf",sep=""),width=4,height=4)
    #par(mar=c(3, 3, 3, 3) + 0.1)  

    filename<-paste("../data/",dataset,additional_style,"/",trendIDs[1],"/gridded/",trendIDs[1],dataset,"_duration_4seasons.nc",sep="") ; print(filename)
    dur5<-var.get.nc(open.nc(filename),"dur")

    filename<-paste("../data/",dataset,additional_style,"/",trendIDs[2],"/gridded/",trendIDs[2],dataset,"_duration_4seasons.nc",sep="") ; print(filename)
    dur7<-var.get.nc(open.nc(filename),"dur")

    filename<-paste("../data/",dataset,additional_style,"/",trendIDs[3],"/gridded/",trendIDs[3],dataset,"_duration_4seasons.nc",sep="") ; print(filename)
    dur9<-var.get.nc(open.nc(filename),"dur")

    colors=c(rgb(1,0.5,0.5,0.5),rgb(0.5,1,0.5,0.5),rgb(0.5,0.5,1,0.5))
    colors=c("cyan","green","magenta")
    for (state in 1:3){
        if (state<3){
            y5<-as.vector(dur5[,state,])
            y7<-as.vector(dur7[,state,])
            y9<-as.vector(dur9[,state,])
        }

        if (state==3){
            y5<-as.vector(dur5[,,])
            y7<-as.vector(dur7[,,])
            y9<-as.vector(dur9[,,])
        }

        br<-0:max(c(max(y5,na.rm=TRUE),max(y7,na.rm=TRUE),max(y9,na.rm=TRUE)),na.rm=TRUE)
        tmp5<<-hist(y5,breaks=br,plot=FALSE)
        tmp7<<-hist(y7,breaks=br,plot=FALSE)
        tmp9<<-hist(y9,breaks=br,plot=FALSE)

        trend_labels <- c("5 years","7 years","9 years")

        plot(NA,xlim=c(0,50),ylim=c(0,500000),main="",ylab="",xlab="")
        points(tmp5$mids,tmp5$count,col=rgb(1,0.5,0.5,0.5),cex=0.3)
        points(tmp7$mids,tmp7$count,col=rgb(0.5,1,0.5,0.5),cex=0.3)
        points(tmp9$mids,tmp9$count,col=rgb(0.5,0.5,1,0.5),cex=0.3)
        legend("topright",col=colors,pch=c(1,1,1),legend=trend_labels,bty = "n")

        plot(NA,xlim=c(0,50),ylim=c(100,500000),main="",ylab="",xlab="",log="y")
        points(tmp5$mids[1:50],tmp5$count[1:50],col=rgb(1,0.5,0.5,0.5),cex=0.3)
        points(tmp7$mids[1:50],tmp7$count[1:50],col=rgb(0.5,1,0.5,0.5),cex=0.3)
        points(tmp9$mids[1:50],tmp9$count[1:50],col=rgb(0.5,0.5,1,0.5),cex=0.3)
        legend("topright",col=colors,pch=c(1,1,1),legend=trend_labels,bty = "n")

        plot(NA,xlim=c(4,50),ylim=c(-5000,2000),main="",ylab="",xlab="")
        abline(h=0,col="gray",lty=2)
        points(tmp5$mids[1:50],tmp5$count[1:50]-tmp7$count[1:50],col=colors[1],cex=0.3)
        points(tmp7$mids[1:50],tmp9$count[1:50]-tmp7$count[1:50],col=colors[3],cex=0.3)
        legend("bottomright",col=colors[c(1,3)],pch=c(1,1,1),legend=c("5 year","9 year"),bty = "n")

        plot(NA,xlim=c(0,50),ylim=c(-500,15000),main="",ylab="",xlab="")
        points(tmp5$mids[1:50],tmp5$count[1:50]-tmp7$count[1:50],col=rgb(1,0.5,0.5,0.5),cex=0.3)
        points(tmp7$mids[1:50],tmp5$count[1:50]-tmp9$count[1:50],col=rgb(0.5,1,0.5,0.5),cex=0.3)
        points(tmp9$mids[1:50],tmp7$count[1:50]-tmp9$count[1:50],col=rgb(0.5,0.5,1,0.5),cex=0.3)
        legend("topright",col=colors,pch=c(1,1,1),legend=c("5-7","5-9","7-9"))

        if (FALSE){
            plot.ecdf(y5,xlim=c(0,60),ylim=c(0,1),main="",cex=0.1,col=colors[1],ylab="",xlab="",axes=TRUE,frame.plot=TRUE)
            par(new=TRUE)
            plot.ecdf(y7,xlim=c(0,60),ylim=c(0,1),main="",cex=0.1,col=colors[2],ylab="",xlab="",axes=TRUE,frame.plot=TRUE)
            par(new=TRUE)
            plot.ecdf(y9,xlim=c(0,60),ylim=c(0,1),main="",cex=0.1,col=colors[3],ylab="",xlab="",axes=TRUE,frame.plot=TRUE)
            text(40,0.4,)

            ks_test[state,1]=ks.test(y5,y7)$p.value
            wilcox_test[state,1]=wilcox.test(y5,y7)$p.value
            ks_test[state,2]=ks.test(y5,y9)$p.value
            wilcox_test[state,2]=wilcox.test(y5,y9)$p.value
            ks_test[state,3]=ks.test(y7,y9)$p.value
            wilcox_test[state,3]=wilcox.test(y5,y7)$p.value

            text(40,0.5,paste("ks 5~7:",round(ks_test[state,1],02)))
            text(40,0.3,paste("ks 5~9:",round(ks_test[state,2],02)))
            text(40,0.1,paste("ks 7~9:",round(ks_test[state,3],02)))
        }

    }
    graphics.off()
}


distr_comparison_plot <- function(X,cdf1,cdf2,D_val,D_pos,ks,xUntil,season,state,q){
    par(mar=c(3, 3, 3, 3) + 0.1)     
    color=c(rgb(0.2,0.7,0.2),rgb(1,0.5,0))       

    plot(NA,xlim=c(0,xUntil),ylim=c(0,1),xlab="days",ylab="",axes=FALSE,frame.plot=TRUE)
    points(X,cdf1,col=color[1],pch=20,cex=0.4)
    points(X,cdf2,col=color[2],pch=20,cex=0.4)

    lines(c(X[D_pos],X[D_pos]),c(cdf1[D_pos],cdf2[D_pos]),col="black",lwd=2)
    
    text(40,0.4,paste("D =",round(D_val,03)),col="black",pos=1)                 
}

ks_wilcoxon_period <- function(ID_name,yearPeriod1,yearPeriod2,periods,folder=paste("/regional/",ID_name,"/",sep=""),ID_select,ID_length=length(ID_select),distr_cut=4){
    # this old, can only be used for plots

    tests=array(NA,c(5,ID_length,2,10))

    pdf(paste("../plots/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_ks_wilcox_",periods[1],"_vs_",periods[2],"_ecdf.pdf",sep=""),width=3,height=3)
    par(mar=c(3, 3, 3, 3) + 0.1)  

    for (sea in 1:5){ 
        season<-season_names[sea]
        cat(paste("\n",season))  

        print(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_duration_",season,".nc",sep=""))
        nc_dur<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_duration_",season,".nc",sep=""))
        dur<-var.get.nc(nc_dur,"dur")
        dur_mid<-var.get.nc(nc_dur,"dur_mid")

        percentage<-0
        cat(paste("\n0 -> -> -> -> -> 100\n"))
        for (q in ID_select){
            if (q/ID_length*100 > percentage){
                cat("-")
                percentage<-percentage+5
            }
            for (state in 1:2){
                ord<-order(dur_mid[q,state,])
                y<-as.vector(dur[q,state,][ord])
                x<-as.vector(dur_mid[q,state,][ord])
                inYearPeriod<-which(x>yearPeriod1[1] & x<yearPeriod1[2])
                y1<-y[inYearPeriod]
                inYearPeriod<-which(x>yearPeriod2[1] & x<yearPeriod2[2])
                y2<-y[inYearPeriod]



                # for whole distribution
                longest_period<-max(c(y1[!is.na(y1)],y1[!is.na(y1)]))
                cdf1<-cdf_onData(y1,1:longest_period)
                cdf2<-cdf_onData(y2,1:longest_period)

                diff_cdf<-abs(cdf1-cdf2)
                D_val<-max(diff_cdf)
                D_pos<-which.max(diff_cdf)

                tests[sea,q,state,1]=ks.test(y1,y2)$p.value
                tests[sea,q,state,2]=D_val
                tests[sea,q,state,3]=D_pos

                # wilcoxon test
                tests[sea,q,state,4]=wilcox.test(x=y1,y=y2,paired=FALSE)$p.value

                plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=FALSE,frame.plot=FALSE)
                text(0.5,0.5,paste(season,q,state_names[state]))
                distr_comparison_plot(X=1:longest_period,cdf1=cdf1,cdf2=cdf2,D_val=D_val,D_pos=D_pos,ks=tests[sea,q,state,1],xUntil=70,season=season,state=state,q=q)
                #text(20,0.18,expression('KS'['full']),col="black",pos=4)   
                #text(38,0.18,paste("=",round(tests[sea,q,state,1],03)),col="black",pos=4)   


                # for period longer than ...
                distr_cut<-median(c(y1,y2),na.rm=TRUE)
                tests[sea,q,state,10]<-distr_cut

                y1<-y1[y1>=distr_cut]
                y2<-y2[y2>=distr_cut]
                longest_period<-max(c(y1[!is.na(y1)],y1[!is.na(y1)]))
                cdf1<-cdf_onData(y1,distr_cut:longest_period)
                cdf2<-cdf_onData(y2,distr_cut:longest_period)

                diff_cdf<-abs(cdf1-cdf2)
                D_val<-max(diff_cdf)
                D_pos<-which.max(diff_cdf)

                tests[sea,q,state,6]=ks.test(y1,y2)$p.value
                tests[sea,q,state,7]=D_val
                tests[sea,q,state,8]=D_pos

                # wilcoxon test
                tests[sea,q,state,9]=wilcox.test(x=y1,y=y2,paired=FALSE)$p.value

                distr_comparison_plot(X=distr_cut:longest_period,cdf1=cdf1,cdf2=cdf2,D_val=D_val,D_pos=D_pos,ks=tests[sea,q,state,6],xUntil=70,season=season,state=state,q=q)
                abline(v=distr_cut,col="gray",lty=2)
                #text(20,0.18,expression('KS'['upper']),col="black",pos=4)   
                #text(43,0.18,paste("=",round(tests[sea,q,state,6],03)),col="black",pos=4)
            }
        }
    }

    print(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_ks_wilcox_",yearPeriod1[1],"-",yearPeriod1[2],"_vs_",yearPeriod2[1],"-",yearPeriod2[2],".nc",sep=""))
    nc_out <- create.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_ks_wilcox_",yearPeriod1[1],"-",yearPeriod1[2],"_vs_",yearPeriod2[1],"-",yearPeriod2[2],".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "comment", "NC_CHAR", "wilcoxon test to identify changes in distribution")
    
    dim.def.nc(nc_out,"seasons",dimlength=5, unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=10,unlim=FALSE)


    var.def.nc(nc_out,"tests","NC_DOUBLE",c(0,1,2,3))
    att.put.nc(nc_out, "tests", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "tests", "dim_explanation", "NC_CHAR", "season-ID-states-...")
    att.put.nc(nc_out, "tests", "val_explanation", "NC_CHAR", "full distr :ks test: 1-p-value, 2-D_val, 3-D_pos, 5-wilcox-p-value ; part distr :ks test: 6-p-value, 7-D_val, 8-D_pos, 9-wilcox-p-value, 10-median c(y1,y2)")

    var.put.nc(nc_out,"tests",tests)              

    close.nc(nc_out) 
    graphics.off()
}




fit_plot_comparison <- function(distr,fits,sea,q,state,ks){
    color=c(rgb(0.2,0.7,0.2),rgb(1,0.5,0))

    par(mar=c(3, 3, 3, 3) + 0.1)  

    #index page
    plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=FALSE,frame.plot=FALSE)
    text(0.5,0.5,paste(season_names[sea],q,state_names[state]))

    # first plot page
    plot(NA,xlab="days",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],cex=0.5)
    #text(40,0.18+0.06,paste("KS_full=",round(ks[1],03),"\n"," KS_upper=",round(ks[6],03)),pos=1,col="black")   

    if (!is.na(ks[1,1])){
        text(24,0.21,expression('KS'['full']),col="black",pos=4)   
        if (ks[1,1]>=ks[1,7]){text(43,0.21,paste(">95%"),col="black",pos=4)}
        if (ks[1,1]>=ks[1,6] & ks[1,1]<ks[1,7]){text(43,0.21,paste(">90%"),col="black",pos=4)}
        if (ks[1,1]<=ks[1,6]){text(43,0.21,paste("<90%"),col="black",pos=4)}

        text(20,0.18,expression('KS'['upper']),col="black",pos=4)   
        if (ks[2,1]>=ks[2,7]){text(43,0.18,paste(">95%"),col="black",pos=4)}
        if (ks[2,1]>=ks[2,6] & ks[2,1]<ks[2,7]){text(43,0.18,paste(">90%"),col="black",pos=4)}
        if (ks[2,1]<=ks[2,6]){text(43,0.18,paste("<90%"),col="black",pos=4)}
    }

    for (i in 1:2){
        nonull<-which(distr[i,3,]>0)
        points(distr[i,1,nonull],distr[i,2,nonull],pch=20,col=color[i],cex=0.5)
        lines(distr[i,1,nonull],distr[i,8,nonull],col=color[i],lty=1)
    }  
    # second plot page
    plot(NA,xlab="days",ylim=c(0.00001,0.25),xlim=c(0,70),ylab="",axes=FALSE,frame.plot=TRUE,pch=20,col=color[state],cex=0.5,log="y")
    text(50,0.23,paste("b1=",round(fits[1,12],02),"\n\n","b2=",round(fits[1,14],02)),pos=1,col=color[1])                 
    text(50,0.1,paste("b1=",round(fits[2,12],02),"\n\n","b2=",round(fits[2,14],02)),pos=1,col=color[2])                 
    for (i in 1:2){
        nonull<-which(distr[i,3,]>0)
        points(distr[i,1,nonull],distr[i,2,nonull],pch=20,col=color[i],cex=0.5)
        lines(distr[i,1,nonull],distr[i,8,nonull],col=color[i],lty=1)
        abline(v=fits[i,15],col=color[i],lty=2)
    }              
}

distribution_comparision <- function(region_name,periods,period,ID_select,ID_length=length(ID_select),region_names=1:24,hlines=c(30)){

    filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_bootstrap.nc",sep="") ; print(filename)
    original_ks<-var.get.nc(open.nc(filename),"statistics")[,1:24,,3:4,]

    fit_params=array(NA,c(2,5,ID_length,2,30))
    distrs=array(NA,c(2,5,ID_length,2,10,100))
    for (i in 1:length(periods)){
    	part_period<-periods[i]
    	filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",part_period,"/",trendID,"_",dataset,"_",region_name,"_",part_period,"_fit_","2expo_4:100",".nc",sep=""); print(filename)
   	 	fit_params[i,,,,]=var.get.nc(open.nc(filename),"fit_stuff")[1:5,,,]
		
		filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",part_period,"/",trendID,"_",dataset,"_",region_name,"_",part_period,"_distributions.nc",sep=""); print(filename)
		distrs[i,,,,,]=var.get.nc(open.nc(filename),"distr_stuff")[1:5,,,,]
	}

	pdf(paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_dur_wilcox_",periods[1],"_vs_",periods[2],".pdf",sep=""),width=3,height=3)
    for (sea in 1:5){ 
        season<-season_names[sea]
        cat(paste("\n",season))  
        for (q in 1:ID_length){
            for (state in 1:2){
            	fit_plot_comparison(distr=distrs[,sea,q,state,,],fits=fit_params[,sea,q,state,],sea,q,state,ks=original_ks[sea,q,state,,])
            }
        }
    }
    graphics.off()
}



library(dgof)


#ks_wilcoxon_trenID_sensitivity(trendIDs=c("91_5","91_7","91_9"))

# this old, can only be used for plots
#ks_wilcoxon_period(ID_name="ward24",yearPeriod1=c(1979,1995),yearPeriod2=c(1995,2011),periods=c("1979-1995","1995-2011"),ID_select=1:24)
# this old, can only be used for plots

ID_select=c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24)
ID_length=24
hlines=c(23,20,22,8)



distribution_comparision(region_name="ward24",period="1979-2011",periods=c("1979-1995","1995-2011"),ID_select=ID_select,ID_length=ID_length,hlines=hlines)