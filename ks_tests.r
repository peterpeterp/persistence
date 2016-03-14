library(dgof)

wilcoxon_period <- function(ID_name,yearPeriod1,yearPeriod2,folder=paste("/regional/",ID_name,"/",sep=""),ID_select,ID_length=length(ID_select)){

	wilcox_test=array(NA,c(5,ID_length,2,3))
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
                #ks_tests[sea,q,state]=ks.test(y1,y2)$p.value

                noise=0.001

                wilcox_test[sea,q,state,1]=wilcox.test(x=y1,y=y2,paired=FALSE)$p.value
                print(wilcox.test(x=y1,y=y2,paired=FALSE)$p.value)
                print(wilcox.test(x=dither(y1,value=noise),y=dither(y2,value=noise),paired=FALSE)$p.value)
                print(ks.test(x=y1,y=y2)$p.value)
                print(ks.test(x=dither(y1,value=noise),y=dither(y2,value=noise))$p.value)
                print(wilcox.test(x=y1,y=y2,paired=FALSE)$p.value)
                print(wilcox.test(x=dither(y1,value=noise),y=dither(y2,value=noise),paired=FALSE)$p.value)
                print(ks.test(x=y1,y=y2)$p.value)
                print(ks.test(x=y1,y=y1)$p.value)
                print(wilcox.test(x=y1,y=y1,paired=FALSE)$p.value)
                adssa

            }
        }
    }
    wilcox_test[,,,2][wilcox_test[,,,1]>0.05]=1
    wilcox_test[,,,2][wilcox_test[,,,1]<0.05]=0
    wilcox_test[,,,3][wilcox_test[,,,1]>0.1]=1
    wilcox_test[,,,3][wilcox_test[,,,1]<0.1]=0

    print(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_wilcox_",yearPeriod1[1],"-",yearPeriod1[2],"_vs_",yearPeriod2[1],"-",yearPeriod2[2],".nc",sep=""))
    nc_out <- create.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_wilcox_",yearPeriod1[1],"-",yearPeriod1[2],"_vs_",yearPeriod2[1],"-",yearPeriod2[2],".nc",sep=""))
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "comment", "NC_CHAR", "wilcoxon test to identify changes in distribution")
    
    dim.def.nc(nc_out,"seasons",dimlength=5, unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)

    dim.def.nc(nc_out,"outs",dimlength=3,unlim=FALSE)

    var.def.nc(nc_out,"wilcox_test","NC_DOUBLE",c(0,1,2,3))
    att.put.nc(nc_out, "wilcox_test", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "wilcox_test", "dim_explanation", "NC_CHAR", "season-ID-states-...")
    att.put.nc(nc_out, "wilcox_test", "val_explanation", "NC_CHAR", "wilcoxon test p-value, p-value>0.05 = 1, p-value>0.1 = 1")

    var.put.nc(nc_out,"wilcox_test",wilcox_test)              

    close.nc(nc_out) 
}

distribution_comparision <- function(ID_name,periods,folder=paste("/regional/",ID_name,"/",sep=""),ID_select,ID_length=length(ID_select),add_name="2expo_4:100"){

    print(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_wilcox_",periods[1],"_vs_",periods[2],".nc",sep=""))
   	print(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_wilcox_",periods[1],"_vs_",periods[2],".nc",sep=""))
    nc <- open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_wilcox_",periods[1],"_vs_",periods[2],".nc",sep=""))
    wilcox_test <- var.get.nc(nc,"wilcox_test")

    fit_params=array(NA,c(2,6,ID_length,2,20))
    distrs=array(NA,c(2,6,ID_length,2,5,100))
    for (i in 1:length(periods)){
    	period<-periods[i]
    	print(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_fit_",add_name,".nc",sep=""))
   	 	nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_fit_",add_name,".nc",sep=""))
   	 	fit_params[i,,,,]=var.get.nc(nc,"fit_stuff")
		
		print(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_distributions.nc",sep=""))
		nc=open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_distributions.nc",sep=""))
		distrs[i,,,,,]=var.get.nc(nc,"distr_stuff")
	}

	pdf(paste("../plots/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_wilcox_",periods[1],"_vs_",periods[2],".pdf",sep=""),width=3,height=3)
    for (sea in 1:5){ 
        season<-season_names[sea]
        cat(paste("\n",season))  
        percentage<-0
        cat(paste("\n0 -> -> -> -> -> 100\n"))
        for (q in ID_select){
            if (q/ID_length*100 > percentage){
                cat("-")
                percentage<-percentage+5
            }
            for (state in 1:2){
            	fit_plot_comparison(distr=distrs[,sea,q,state,,],fits=fit_params[,sea,q,state,],sea,q,state,wilcox=wilcox_test[sea,q,state,1])
            }
        }
    }
    graphics.off()
}



wilcoxon_period(ID_name="ward23",yearPeriod1=c(1950,1980),yearPeriod2=c(1980,2014),ID_select=1:23)

distribution_comparision(ID_name="ward23",periods=c("1950-1980","1980-2014"),ID_select=1:23,add_name="2expo_4:100")