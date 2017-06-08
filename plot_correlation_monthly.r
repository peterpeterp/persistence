#!/home/pepflei/R/bin/Rscript


variable='spi3'

filename<-paste('/home/peter/Dokumente/pik/backuped/data/_TMean/91_7/gridded/91_7_TMean_duration_',variable,'_cor.nc',sep="") ; print(filename)
correlation=var.get.nc(open.nc(filename),"correlation")
reihen=array(NA,dim=c(10,ntot))
reihen_sig=array(NA,dim=c(10,ntot))
titel=c()

for (sea in 1:5){
	for (state in 1:2){
		reihen[((sea-1)*2+state),]=correlation[,sea,state,1]
		reihen_sig[((sea-1)*2+state),]=correlation[,sea,state,2]
	}
}
filename_plot<-paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_correl_",toCor_name,"_",toCor_shortZu,val_zusatz,".pdf",sep="")
topo_map_plot(filename=filename_plot,reihen=reihen,reihen_sig=reihen_sig,farb_mitte=farb_mitte,farb_palette="gold-blau",titel=c(""))

