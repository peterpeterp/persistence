

evaluate_bootstraping <- function(overRegions=c("NHml","NHpo","NHst","tro","SHml"),period="1979-2011",region_name="ward24"){

	result<-array(NA,c(5,2,30,4,10))

	for (sea in 1:5){
		for (state in 1:2){
			shuffled<-array(NA,c(10001,30,4))
			for (oR in 1:length(overRegions)){
				
				for (i in 1:10){
					filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/bootstrap/",trendID,dataset,"_",period,"_boot_",region_name,"_",overRegions[oR],"_",season_names[sea],"_",state_names[state],"_",i,".nc",sep="") ; print(filename)
					regions<-try(var.get.nc(open.nc(filename),"regions"))
					#print(regions)
					if (regions[1]!="Error : Datei oder Verzeichnis nicht gefunden\n"){
						regions<-regions[1:(length(regions)-1)]

						#  for Regions
						shuffled[((i-1)*1000+2):(i*1000+1),regions,]<-var.get.nc(open.nc(filename),"statistics")[2:1001,1:length(regions),]
						shuffled[1,regions,]<-var.get.nc(open.nc(filename),"statistics")[1,1:length(regions),]

						# for overRegions
						shuffled[((i-1)*1000+2):(i*1000+1),(25+oR),]<-var.get.nc(open.nc(filename),"statistics")[2:1001,(length(regions)+1),]
						shuffled[1,(25+oR),]<-var.get.nc(open.nc(filename),"statistics")[1,(length(regions)+1),]
					}
				}

				for (out in 1:4){
					for (q in c(1:24,26:30)){
						result[sea,state,q,out,1]<-shuffled[1,q,out]
						result[sea,state,q,out,3:8]<-quantile(shuffled[2:10001,q,out],c(0.025,0.05,0.1,0.9,0.95,0.975),na.rm=TRUE)
					}
				}
			}
		}
	}
	result<<-result

	filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_bootstrap.nc",sep="") ; print(filename)

    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "bootstrapping results", "NC_CHAR","see regions")
    att.put.nc(nc_out, "NC_GLOBAL", "replics", "NC_CHAR","10000")
            
    dim.def.nc(nc_out,"seasons",dimlength=5, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2, unlim=FALSE)
    dim.def.nc(nc_out,"regions",dimlength=30,unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=4,unlim=FALSE)
    dim.def.nc(nc_out,"confidence",dimlength=10,unlim=FALSE)

    var.def.nc(nc_out,"statistics","NC_DOUBLE",c(0,1,2,3,4))
    att.put.nc(nc_out, "statistics", "missing_value", "NC_DOUBLE", -99999)
    att.put.nc(nc_out, "statistics", "dim_explanation", "NC_CHAR", "season-state-(regions,NA,overRegions)-(lr,qr,ks_full,ks_upper)-(original,NA,0.025,0.05,0.1,0.9,0.95,0.975,NA,NA)")

    var.put.nc(nc_out,"statistics",result)   

    close.nc(nc_out) 
}

write_slope_table <- function(period="1979-2011",region_name="ward24"){
    filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_bootstrap.nc",sep="") ; print(filename)
    original_slopes<-var.get.nc(open.nc(filename),"statistics")[,,1:24,,]

    signis<-array(NA,c(5,24,2,5,2))
    signis[,,,4,1][which(original_slopes[,,,1,1]>original_slopes[,,,1,6])]=1
    signis[,,,4,1][which(original_slopes[,,,1,1]<original_slopes[,,,1,5])]=1
    signis[,,,5,1][which(original_slopes[,,,1,1]>original_slopes[,,,1,7])]=1
    signis[,,,5,1][which(original_slopes[,,,1,1]<original_slopes[,,,1,4])]=1

    signis[,,,4,2][which(original_slopes[,,,2,1]>original_slopes[,,,2,6])]=1
    signis[,,,4,2][which(original_slopes[,,,2,1]<original_slopes[,,,2,5])]=1
    signis[,,,5,2][which(original_slopes[,,,2,1]>original_slopes[,,,2,7])]=1
    signis[,,,5,2][which(original_slopes[,,,2,1]<original_slopes[,,,2,4])]=1

    values<-array(NA,c(5,24,2,2))
    values[,,,1]=original_slopes[,,,1,1]
    values[,,,2]=original_slopes[,,,2,1]

    print(values)
    values[which(is.na(values))]=0
    values<-values*365

    nbcol<<-101
    plot_reg_table_general(values=values,signis=signis,filename=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,dataset,"_",region_name,"_",period,"_slope_boot.pdf",sep=""),val_names=c("mn","95"),region_name="ward24",colorRange=c(-0.1,0.1),farb_palette="lila-gruen",ID_select=reg_order,hlines=hlines)
}

write_slope_table_control <- function(period="1979-2011",region_name="ward24"){
    filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_quantiles.nc",sep="") ; print(filename)
    qr<-var.get.nc(open.nc(filename),"quantile_stuff")
    filename<-paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_others.nc",sep="") ; print(filename)
    mn<-var.get.nc(open.nc(filename),"other_stuff")

    signis<-array(NA,c(5,24,2,5,2))

    values<-array(NA,c(5,24,2,2))
    print(dim(mn))
    print(dim(values))
    values[,,,1]=mn[1:5,,,4]
    values[,,,2]=qr[1:5,,,2,2]

    print(values)
    values[which(is.na(values))]=0
    values<-values

    nbcol<<-101
    plot_reg_table_general(values=values,signis=signis,filename=paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,dataset,"_",region_name,"_",period,"_slope_control.pdf",sep=""),val_names=c("mn","95"),region_name="ward24",colorRange=c(-0.1,0.1),farb_palette="lila-gruen",ID_select=reg_order,hlines=hlines)
}


init <- function(){
    library(RNetCDF)
    library(quantreg)
    nday<<-91
    nyr<<-7
    trendID<<-paste(nday,"_",nyr,sep="")
    dataset<<-"_TMean"
    additional_style<<-""
    taus<<-c(0.75,0.95,0.99)
    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    state_names<<-c("cold","warm")

    reg_order<<-c(1,2,6,10,13,19,23,3,4,7,12,16,20,5,11,14,18,21,22,17,8,9,15,24)

    hlines<<-c(23,20,22,8)
}

init()

#evaluate_bootstraping()

write_slope_table()
write_slope_table_control()
adas




pdf("../plots/test_boot.pdf")
for (out in c(2,5)){
	for (state in 1:2){
		plot(NA,xlim=c(0,10),ylim=c(-0.001,0.001))
		count<-0
		for (reg in c(3,4,7,12,16,20)){
			count<-count+1
			evaluate_bootstraping(x=count,out=out,region=reg,state=state)
		}
	}
}

for (out in c(2,5)){
	for (state in 1:2){
		plot(NA,xlim=c(0,10),ylim=c(-0.1,0.1))
		count<-0
		for (reg in c(3,4,7,12,16,20)){
			count<-count+1
			lines(c(count+0.1,count+0.1),confi_quantiles[2,reg,state,3,c(2,4)],col="lightblue")
			points(count+0.1,original_slopes[2,reg,state,3,1],col="green")
			text(count,0.08,round(original_slopes[2,reg,state,3,1]/confi_quantiles[2,reg,state,3,4],03))
		}
	}
}

graphics.off()