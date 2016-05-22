library(Kendall)

duration_yearly_values <- function(folder="/gridded/",ID_name="",ID_select=1:length(dat$ID),ID_length=length(ID_select)){

    yearly_values=array(NA,c(5,ID_length,2,length(dat$year),7))

    for (sea in 1:5){
    	season<-season_names[sea]
        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_reg_binned_duration_",season,".nc",sep="") ; print(filename)
        nc<-open.nc(filename)
        binned_dur<-var.get.nc(nc,"binned_dur")

        percentage<-0
        cat(paste("\n0 -> -> -> -> -> 100\n"))
        for (q in ID_select){
            if (q/ID_length*100 > percentage){
                cat("-")
                percentage<-percentage+5
            }
        	for (state in 1:2){
        		for (yr in 1:length(dat$year)){
        			yearly_values[sea,q,state,yr,7]=length(which(!is.na(binned_dur[q,state,,yr])))
        			if (yearly_values[sea,q,state,yr,7]>10){
        				#print(binned_dur[q,state,,yr])
        				#print(paste(q,state,yr))
	        			yearly_values[sea,q,state,yr,5]=mean(binned_dur[q,state,,yr],na.rm=TRUE)
	        			yearly_values[sea,q,state,yr,1:length(taus)]=quantile_pete(binned_dur[q,state,,yr],taus=taus,na.rm=TRUE)
	        		}
        		}
        	}
        }
    }

    filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_yearly_values.nc",sep="") ; print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "taus", "NC_CHAR", paste(taus))
            
    dim.def.nc(nc_out,"seasons",dimlength=5, unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
    dim.def.nc(nc_out,"years",dimlength=length(dat$year),unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=7,unlim=FALSE)

    var.def.nc(nc_out,"yearly_values","NC_DOUBLE",c(0,1,2,3,4))
    att.put.nc(nc_out, "yearly_values", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "yearly_values", "dim_explanation", "NC_CHAR", "sea-ID-states-years-outs")
    att.put.nc(nc_out, "yearly_values", "val_explanation", "NC_CHAR", "annual: 7-# of periods in year, 5-mean(dur), 1:3-quantile(dur)")

    var.put.nc(nc_out,"yearly_values",yearly_values)              
    close.nc(nc_out) 
}

duration_MannKendall <- function(yearPeriod,folder="/gridded/",ID_name="",ID_select=1:length(dat$ID),ID_length=length(ID_select),region_names=1:ID_length,hlines=c(30),colorbar=FALSE,header=TRUE,regLabel=1:length(dat$ID)){
    filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_yearly_values.nc",sep="") ; print(filename)
    yearly_values<-var.get.nc(open.nc(filename),"yearly_values")
    period<-paste(yearPeriod[1],"-",yearPeriod[2],sep="")

    pdf(file=paste("../plots/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,dataset,"_",ID_name,"_",period,"_yearlyAverages.pdf",sep=""),width=6,height=4)
    color=c(rgb(0,1,0),rgb(1,0.5,0),rgb(1,0,0),"black","blue")
    MK=array(NA,c(5,ID_length,2,5,2))
    for (sea in 1:5){
    	for (q in ID_select){
    		for (state in 1:2){
                plot(NA,xlim=c(0,65),ylim=c(0,60),xlab="",ylab="",main=paste(sea,q,state))
    			for (out in c(1,2,3,5)){
                    lines(1:65,yearly_values[sea,q,state,,out],col=color[out])
    				tmp<-MannKendall(yearly_values[sea,q,state,(yearPeriod[1]-dat$year[1]+1):(yearPeriod[2]-dat$year[1]+1),out])
    				MK[sea,q,state,out,1]=tmp$tau
    				MK[sea,q,state,out,2]=tmp$sl
    			}
    		}
    	}
    }
    graphics.off()

    filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,dataset,"_",ID_name,"_",period,"_MK.nc",sep="") ; print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "taus", "NC_CHAR", paste(taus))
            
    dim.def.nc(nc_out,"seasons",dimlength=5, unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=5,unlim=FALSE)
    dim.def.nc(nc_out,"results",dimlength=2,unlim=FALSE)

    var.def.nc(nc_out,"MK","NC_DOUBLE",c(0,1,2,3,4))
    att.put.nc(nc_out, "MK", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "MK", "dim_explanation", "NC_CHAR", "sea-ID-states-(1:3-taus 5-mean)-results")
    att.put.nc(nc_out, "MK", "val_explanation", "NC_CHAR", "1-MK, 2-MK p-value")

    var.put.nc(nc_out,"MK",MK)              
    close.nc(nc_out) 

    filename<-paste("../plots/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,dataset,"_",ID_name,"_",period,"_MK.tex",sep="") ; print(filename)
 
    signis<-array(NA,c(5,ID_length,2,5,2))
    #signis[,,,1,1][which(MK[,,,2,2]>0.1)]=1
    signis[,,,5,1][which(MK[,,,2,2]<=0.1)]=1
    signis[,,,4,1][which(MK[,,,2,2]<=0.05)]=1

    #signis[,,,1,2][which(MK[,,,5,2]>0.1)]=1
    signis[,,,5,2][which(MK[,,,5,2]<=0.1)]=1
    signis[,,,4,2][which(MK[,,,5,2]<=0.05)]=1


    values<-array(NA,c(5,ID_length,2,2))
    values[,,,1][which(MK[,,,5,1]>0)]=1
    values[,,,1][which(MK[,,,5,1]<0)]=-1
    values[,,,2][which(MK[,,,3,1]>0)]=1
    values[,,,2][which(MK[,,,3,1]<0)]=-1


    nbcol<<-2
    plot_reg_table_general(values=values,signis=signis,filename=paste("../plots/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,dataset,"_",ID_name,"_",period,"_MK.pdf",sep=""),val_names=c("mn","95"),region_name="ward24",colorRange=c(-0.1,0.1),farb_palette="lila-gruen",ID_select=ID_select,hlines=hlines,style="yn",colorbar=colorbar,header=header,regLabel=regLabel)
}
