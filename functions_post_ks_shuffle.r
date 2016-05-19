
confidence_interval <- function(seasons=1:5){
    confi_quantiles<-array(NA,dim=c(5,regNumb,2,2,5))
    original_ks<-array(NA,dim=c(5,regNumb,2,2,5))

    pdf(paste("../plots/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_ks_test_",periods[1],"_vx_",periods[2],"_KS.pdf",sep=""),width=3,height=3)
    par(mar=c(3, 3, 3, 3) + 0.1)  
    color=c(rgb(0.7,0.2,0.7),rgb(0.5,0.7,0.5))
    ymax<-2000
    plot(NA,xlab="",ylim=c(0,ymax),xlim=c(0,0.2),ylab="",axes=FALSE,frame.plot=TRUE,cex=0.5) ; at_=axis(1,labels=FALSE,col="black") ; if (length(at_)>4){at_=at_[2:(length(at_)-1)]} ; axis(1,at=at_)
    plot(NA,xlab="",ylim=c(0,ymax),xlim=c(0,0.2),ylab="",axes=FALSE,frame.plot=TRUE,cex=0.5) ; at_=axis(3,labels=FALSE,col="black") ; if (length(at_)>4){at_=at_[2:(length(at_)-1)]} ; axis(3,at=at_)
    plot(NA,xlab="",ylim=c(0,ymax),xlim=c(0,0.2),ylab="",axes=FALSE,frame.plot=TRUE,cex=0.5) ; at_=axis(2,labels=FALSE,col="black") ; if (length(at_)>4){at_=at_[2:(length(at_)-1)]} ; axis(2,at=at_)
    plot(NA,xlab="",ylim=c(0,ymax),xlim=c(0,0.2),ylab="",axes=FALSE,frame.plot=TRUE,cex=0.5) ; at_=axis(4,labels=FALSE,col="black") ; if (length(at_)>4){at_=at_[2:(length(at_)-1)]} ; axis(4,at=at_)

    for (sea in seasons){
        season<-season_names[sea]
        shuffled_mass<-array(NA,dim=c(10000,regNumb,2,2))
        original_mass<-array(NA,dim=c(10,regNumb,2,2))
        for (id in 1:10){
            nc <- try(open.nc(paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/shuffled_ks/",trendID,dataset,"_",ID_name,"_",period,"_shuffled_ks_",season,"_",id,".nc",sep="")))
            if (class(nc)!="try-error"){
                shuffled<-var.get.nc(nc,"shuffled")
                original<-var.get.nc(nc,"original")
                shuffled_mass[((id-1)*1000+1):(id*1000),,,]=shuffled
                original_mass[id,,,]=original
            }
        }
        original_ks[sea,,,,1]=original_mass[which(!is.na(original_mass[,1,1,1]))[1],,,]

        for (q in 1:regNumb){
            for (state in 1:2){
                plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=FALSE,frame.plot=FALSE) ; text(0.5,0.5,paste(season_names[sea],q,state_names[state]))
                for (out in c(1,2)){
                    confi_quantiles[sea,q,state,out,]=quantile(shuffled_mass[,q,state,out],c(0.5,0.75,0.9,0.95,0.99),na.rm=TRUE)
                    if (original_ks[sea,q,state,out,1]>confi_quantiles[sea,q,state,out,3]){original_ks[sea,q,state,out,3]=-0.1}
                    if (original_ks[sea,q,state,out,1]>confi_quantiles[sea,q,state,out,4]){original_ks[sea,q,state,out,2]=-0.05}
                    
                    hist(shuffled_mass[,q,state,out],br=seq(0,0.2,0.005),ylim=c(0,ymax),xlim=c(0,0.2),main="",ylab="",xlab="",axes=FALSE,col=color[out],border=color[out])
                    abline(v=confi_quantiles[sea,q,state,out,3],lty=3,col="gray")
                    abline(v=confi_quantiles[sea,q,state,out,4],lty=2,col="gray")
                    abline(v=original_ks[sea,q,state,out,1],lty=1,col="black")

                    if (out==1){text(0.085,1800,expression('KS'['full']),col="black",pos=4)}
                    if (out==2){text(0.07,1800,expression('KS'['upper']),col="black",pos=4)}
                    if (original_ks[sea,q,state,out,1]<confi_quantiles[sea,q,state,out,3]){text(0.135,1800,"<90%",pos=4)}
  
                    if (original_ks[sea,q,state,out,1]>confi_quantiles[sea,q,state,out,3] & original_ks[sea,q,state,out,1]<confi_quantiles[sea,q,state,out,4]){text(0.135,1800,">90%",pos=4)}
                    if (original_ks[sea,q,state,out,1]>confi_quantiles[sea,q,state,out,4]){text(0.135,1800,">95%",pos=4)}
                }
            }
        }
    }
    original_ks<<-original_ks
    confi_quantiles<<-confi_quantiles

    filename <- paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuff_ks.nc",sep="")
    nc_out<-create.nc(filename)

    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)
    
    dim.def.nc(nc_out,"seasons",dimlength=5,unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=regNumb, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)

    dim.def.nc(nc_out,"outs",dimlength=2,unlim=FALSE)
    dim.def.nc(nc_out,"confi_quants",dimlength=5,unlim=FALSE)
        

    var.def.nc(nc_out,"confi_quantiles","NC_DOUBLE",c(0,1,2,3,4))
    att.put.nc(nc_out, "confi_quantiles", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "confi_quantiles", "dim_explanation", "NC_CHAR", "season-ID-state-...")
    att.put.nc(nc_out, "confi_quantiles", "explanation", "NC_CHAR", "(full,upper) x (shuflle_slope quantiles 0.5,0.75,0.9,0.95,0.975)")

    var.def.nc(nc_out,"original_ks","NC_DOUBLE",c(0,1,2,3,4))
    att.put.nc(nc_out, "original_ks", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "original_ks", "dim_explanation", "NC_CHAR", "season-ID-state-...")
    att.put.nc(nc_out, "original_ks", "explanation", "NC_CHAR", "(full,upper) x (ks-statistic, 0.05 signi, 0.1 signi")
        
    var.put.nc(nc_out,"original_ks",original_ks) 
    var.put.nc(nc_out,"confi_quantiles",confi_quantiles) 

    close.nc(nc_out)
    graphics.off()     
}

plot_distr_compa_table <- function(){

    filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_shuff_ks.nc",sep="")
    original_ks<-var.get.nc(open.nc(filename),"original_ks")

    fit_params=array(NA,c(2,6,ID_length,2,30))
    quantiles=array(NA,c(2,6,ID_length,2,length(taus)+1,3))
    distrs=array(NA,c(2,6,ID_length,2,10,100))
    for (i in 1:length(periods)){
        period<-periods[i]
        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_fit_","2expo_4:100",".nc",sep=""); print(filename)
        fit_params[i,,,,]=var.get.nc(open.nc(filename),"fit_stuff")

        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_quantiles",".nc",sep=""); print(filename)
        quantiles[i,,,,1:length(taus),]=var.get.nc(open.nc(filename),"quantile_stuff")

        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_others",".nc",sep=""); print(filename)
        quantiles[i,,,,7,1]=var.get.nc(open.nc(filename),"other_stuff")[,,,1]
        
        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",ID_name,"_",period,"_distributions.nc",sep=""); print(filename)
        distrs[i,,,,,]=var.get.nc(open.nc(filename),"distr_stuff")
    }

    signis<-array(NA,c(5,ID_length,2,5,2))
    #signis[,,,1][which(is.na(original_ks[,,,2,signi]))]=1
    signis[,,,4,2][which(!is.na(original_ks[,,,2,3]))]=1
    signis[,,,5,2][which(!is.na(original_ks[,,,2,2]))]=1

    signis[,,,4,1][which(!is.na(original_ks[,,,1,3]))]=1
    signis[,,,5,1][which(!is.na(original_ks[,,,1,2]))]=1

    print(length(which(signis[,,,4,]==1)))

    values<-array(NA,c(5,ID_length,2,2))
    values[,,,1]=quantiles[2,1:5,,,7,1]-quantiles[1,1:5,,,7,1]
    values[,,,2]=quantiles[2,1:5,,,5,1]-quantiles[1,1:5,,,5,1]

    nbcol<<-101
    plot_reg_table_general(values=values,signis=signis,filename=paste("../plots/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_ks_test_",periods[1],"_vx_",periods[2],"_quantile_diff.pdf",sep=""),val_names=c("mn","95"),region_name="ward24",colorRange=c(-2,2),farb_palette="lila-gruen",ID_select=ID_select,hlines=hlines)

    signis[,,,3,1][which(fit_params[1,1:5,,,21]<0.99)]=1
    signis[,,,3,2][which(fit_params[1,1:5,,,21]<0.99)]=1
    signis[,,,3,1][which(fit_params[2,1:5,,,21]<0.99)]=1
    signis[,,,3,2][which(fit_params[2,1:5,,,21]<0.99)]=1

    signis[,,,2,1][which(fit_params[1,1:5,,,24]>0)]=1
    signis[,,,2,2][which(fit_params[1,1:5,,,24]>0)]=1
    signis[,,,2,1][which(fit_params[2,1:5,,,24]>0)]=1
    signis[,,,2,2][which(fit_params[2,1:5,,,24]>0)]=1

    #nur noch upper
    signis[,,,4,1]=NA
    signis[,,,5,1]=NA
    signis[,,,4,1][which(!is.na(original_ks[,,,2,3]))]=1
    signis[,,,5,1][which(!is.na(original_ks[,,,2,2]))]=1

    values<-array(NA,c(5,ID_length,2,2))
    values[,,,1]=fit_params[2,1:5,,,12]-fit_params[1,1:5,,,12]
    values[,,,2]=fit_params[2,1:5,,,14]-fit_params[1,1:5,,,14]

    fit_params<<-fit_params

    plot_reg_table_general(values=values,signis=array(signis,c(5,ID_length,2,5,2)),filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_ks_test_",periods[1],"_vx_",periods[2],"_slopeDiff.pdf",sep=""),val_names=c("b1","b2"),region_name="ward24",colorRange=c(-0.1,0.1),farb_palette="lila-gruen-inv",ID_select=ID_select,hlines=hlines)

    values<-array(NA,c(5,ID_length,2,1))
    values[,,,1]=fit_params[2,1:5,,,15]-fit_params[1,1:5,,,15]    

    plot_reg_table_general(values=values,signis=array(signis,c(5,ID_length,2,5,1)),filename_plot=paste("../plots/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",ID_name,"_dur_ks_test_",periods[1],"_vx_",periods[2],"_threshDiff.pdf",sep=""),val_names=c(""),region_name="ward24",colorRange=c(-4,4),farb_palette="lila-gruen-inv",ID_select=ID_select,hlines=hlines)

}


init <- function(){
    library(quantreg)
    library(RNetCDF)
    source("plot_tables.r")

    period<<-"1950-2014"
    trendID<<-"91_7"

    dataset<<-"_TMean"
    additional_style<<-""

    state_names<<-c("cold","warm")
    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    out_names<<-c(0.5,0.75,0.95,0.99,"mean")

    ID_name<<-"ward24"
    folder<<-paste("/regional/",ID_name,"/",sep="")
    regNumb<<-24
    ID_length<<-24
    region_names<<-1:24

    plot_select<<-c(3,4,5,7,11,12,13,14,16,18,20,21)
    reg_order<<-c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24)
    ID_select<<-c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24)
    hlines<<-c(23,20,22,8)

    periods<<-c("1950-1980","1980-2014")

    taus<<-c(0.05,0.25,0.5,0.75,0.95,0.99)


}

init()
confidence_interval()
#write_robustness_table()
plot_distr_compa_table()

for (trendID in c("91_7","91_5","91_9")){
    #plot_confi_intervals()
    #write_slope_table(signi=2)
}


