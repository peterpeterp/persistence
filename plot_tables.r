plot_reg_table <- function(region_name="ward24",file="_quantiles",var="quantile_stuff",name_zusatz="quanzs",sub_auswahl=c(NA,NA),value_auswahl=c(1,2),sig_auswahl=c(NA,NA),colorRange=c(2,4.5,4,9,10,22,15,35),ID_select=1:24,hlines=c(30)){
    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    regNumb<-length(unique(attribution[!is.na(attribution)]))

	print(paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,file,".nc",sep=""))
	nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,file,".nc",sep=""))
	values<-var.get.nc(nc,var)

	filename_plot<-paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/","duration_trend_",trendID,"_",region_name,"_",name_zusatz,name_reg_zusatz,"_",period,additional_style,"_table.pdf",sep="") ; print(filename_plot)
	pdf(file=filename_plot,width=6,height=6)
	par(mar=c(3,0,0,0))

	ID_select<-ID_select[length(ID_select):1]

	jet.colors <- colorRampPalette( c( "blue","green","yellow","red") )
	color <- jet.colors(101)	

	for (s in 1:length(sub_auswahl)){
		sub<-sub_auswahl[s]
		plot(NA,xlim=c(0,11),ylim=c(0,29),frame.plot=FALSE,axes=FALSE,xlab="",ylab="")
		for (i in 1:length(ID_select)){text(x=0.5,y=i+1.5,label=ID_select[i])}
		for (sea in season_auswahl){
			text(x=(sea-1)*2+2,y=length(ID_select)+3.5,label=season_names[sea])
			if (TRUE){
				for (val in 1:length(value_auswahl)){
					for (state in 1:2){
						text(x=(sea-1)*2+state+0.5,y=length(ID_select)+2.5,label=state_names[state])
						y<-c(values[sea,,state,sub,val],colorRange[(s-1)*2+1],colorRange[(s-1)*2+2])
						y[y>colorRange[(s-1)*2+2]]=colorRange[(s-1)*2+2]
						y[y<colorRange[(s-1)*2+1]]=colorRange[(s-1)*2+1]
						facetcol <- cut(y,101)
						farben<-color[facetcol]
			            for (i in 1:length(ID_select)){
			            	xPos<-(sea-1)*2+state
			            	yPos<-1+i
			                polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col=farben[ID_select[i]])
			            }
					}
				}
			}
		}
		for (i in 1:length(ID_select)){if (ID_select[i] %in% hlines){lines(c(1,11),c(i+1,i+1),lwd=2)}}
		lines(c(1,11),c(length(ID_select)+2,length(ID_select)+2),lwd=2)
		lines(c(1,11),c(2,2),lwd=2)
		lines(c(1,1),c(2,length(ID_select)+2),lwd=2)
		lines(c(11,11),c(2,length(ID_select)+2),lwd=2)
		par(new=TRUE)
		plot(NA,xlim=c(0,1),ylim=c(1,0),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
		#image.plot(legend.only=T,horizontal=TRUE, zlim=range(y), col=color,add=TRUE,fill=TRUE,smallplot=c(0.1,0.9,0.5,0.9))
		image.plot(legend.only=T,horizontal=TRUE, zlim=range(y[1:length(y)]), col=color,add=FALSE,fill=TRUE,smallplot=c(0.15,0.93,0.1,0.15))
	}
	graphics.off()
}

plot_reg_fit_table <- function(region_name="ward24",file="_quantiles",var="quantile_stuff",name_zusatz="quanzs",value_auswahl=c(12,14,15),val_names=c("1",2),colorRange=c(0,0.4),ID_select=1:24,hlines=c(30)){
    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    regNumb<-length(unique(attribution[!is.na(attribution)]))
	valNumb<-length(value_auswahl)

	print(paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,file,".nc",sep=""))
	nc<-open.nc(paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,file,".nc",sep=""))
	values<-var.get.nc(nc,var)
	print(dim(values))

	filename_plot<-paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",region_name,"_",name_zusatz,name_reg_zusatz,"_",period,"_",file,"_table.pdf",sep="") ; print(filename_plot)
	pdf(file=filename_plot,width=valNumb*3,height=6)
	par(mar=c(3,0,0,0))

	ID_select<-ID_select[length(ID_select):1]

	jet.colors <- colorRampPalette( c("red","yellow","green","blue") )
	color <- jet.colors(101)	

	if (TRUE){
		plot(NA,xlim=c(0,valNumb*10+1),ylim=c(0,29),frame.plot=FALSE,axes=FALSE,xlab="",ylab="")
		for (i in 1:length(ID_select)){text(x=0.5,y=i+1.5,label=ID_select[i])}
		for (sea in season_auswahl){
			if (valNumb>1){text(x=(sea-1)*2*valNumb+3,y=length(ID_select)+4.5,label=season_names[sea])}
			if (valNumb==1){text(x=(sea-1)*2*valNumb+2,y=length(ID_select)+3.5,label=season_names[sea])}
			if (TRUE){
				for (state in 1:2){
					if (valNumb>1){text(x=(sea-1)*2*valNumb+(state-1)*valNumb+2,y=length(ID_select)+3.5,label=state_names[state])}
					if (valNumb==1){text(x=(sea-1)*2*valNumb+(state-1)*valNumb+1.5,y=length(ID_select)+2.5,label=c("c","w")[state])}
					for (v in 1:valNumb){
						text(x=(sea-1)*2*valNumb+(state-1)*valNumb+v+0.5,y=length(ID_select)+2.5,label=val_names[v])
						val<-value_auswahl[v]
						y<-c(values[sea,,state,val],colorRange[1],colorRange[2])
						y[y>colorRange[2]]=colorRange[2]
						y[y<colorRange[1]]=colorRange[1]
						facetcol <- cut(y,101)
						farben<-color[facetcol]
			            for (i in 1:length(ID_select)){
			            	xPos<-(sea-1)*2*valNumb+(state-1)*valNumb+v
			            	yPos<-1+i
			                polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col=farben[ID_select[i]])
			                if(values[sea,ID_select[i],state,24]>0){polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col="black",density=30)}
			                if(values[sea,ID_select[i],state,12]<values[sea,ID_select[i],state,14]){points(xPos+0.5,yPos+0.5,pch=17,col="white",cex=1.5)}	     
			                if(values[sea,ID_select[i],state,21]<0.99){polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border="white",col="white")}			                
			            }
					}
				}
			}
		}

		# border lines
		for (sea in 1:4){lines(c(sea*2*valNumb+1,sea*2*valNumb+1),c(2,length(ID_select)+2),lwd=2)}
		for (i in 1:length(ID_select)){if (ID_select[i] %in% hlines){lines(c(1,valNumb*10+1),c(i+1,i+1),lwd=2)}}
		lines(c(1,valNumb*10+1),c(length(ID_select)+2,length(ID_select)+2),lwd=2)
		lines(c(1,valNumb*10+1),c(2,2),lwd=2)
		lines(c(1,1),c(2,length(ID_select)+2),lwd=2)
		lines(c(valNumb*10+1,valNumb*10+1),c(2,length(ID_select)+2),lwd=2)
		par(new=TRUE)
		plot(NA,xlim=c(0,1),ylim=c(1,0),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
		#image.plot(legend.only=T,horizontal=TRUE, zlim=range(y), col=color,add=TRUE,fill=TRUE,smallplot=c(0.1,0.9,0.5,0.9))
		image.plot(legend.only=T,horizontal=TRUE, zlim=range(y[1:length(y)]), col=color,add=FALSE,fill=TRUE,smallplot=c(0.15,0.93,0.1,0.15))
	}
	graphics.off()
}


plot_reg_table_general <- function(values,signis,notAccepted,filename_plot,val_names,region_name="ward24",colorRange=c(0,0.4),farb_palette="lila-gruen",ID_select=1:24,hlines=c(30)){
    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    regNumb<-length(unique(attribution[!is.na(attribution)]))
	valNumb<-dim(values)[4]


	pdf(file=filename_plot,width=valNumb*3,height=6)
	par(mar=c(3,0,0,0))

	ID_select<-ID_select[length(ID_select):1]

	if (farb_palette=="lila-gruen"){jet.colors <- colorRampPalette( c(rgb(0.1,0.2,0.4),rgb(0.5,1,1),rgb(0.5,1,0.5), rgb(1,1,1),rgb(1,0.7,0.7),rgb(1,0.5,1),rgb(0.4,0.1,0.4)))}
	if (farb_palette=="lila-gruen-inv"){jet.colors <- colorRampPalette( c(rgb(0.1,0.2,0.4),rgb(0.5,1,1),rgb(0.5,1,0.5), rgb(1,1,1),rgb(1,0.7,0.7),rgb(1,0.5,1),rgb(0.4,0.1,0.4))[7:1])}
	color <- jet.colors(101)	

	if (TRUE){
		plot(NA,xlim=c(0,valNumb*10+1),ylim=c(0,29),frame.plot=FALSE,axes=FALSE,xlab="",ylab="")
		for (i in 1:length(ID_select)){text(x=0.5,y=i+1.5,label=ID_select[i])}
		for (sea in season_auswahl){
			if (valNumb>1){text(x=(sea-1)*2*valNumb+3,y=length(ID_select)+4.5,label=season_names[sea])}
			if (valNumb==1){text(x=(sea-1)*2*valNumb+2,y=length(ID_select)+3.5,label=season_names[sea])}
			if (TRUE){
				for (state in 1:2){
					if (valNumb>1){text(x=(sea-1)*2*valNumb+(state-1)*valNumb+2,y=length(ID_select)+3.5,label=state_names[state])}
					if (valNumb==1){text(x=(sea-1)*2*valNumb+(state-1)*valNumb+1.5,y=length(ID_select)+2.5,label=c("c","w")[state])}
					for (v in 1:valNumb){
						text(x=(sea-1)*2*valNumb+(state-1)*valNumb+v+0.5,y=length(ID_select)+2.5,label=val_names[v])
						val<-v
						y<-c(values[sea,,state,val],colorRange[1],colorRange[2])
						y[y>colorRange[2]]=colorRange[2]
						y[y<colorRange[1]]=colorRange[1]
						facetcol <- cut(y,101)
						farben<-color[facetcol]
			            for (i in 1:length(ID_select)){
			            	xPos<-(sea-1)*2*valNumb+(state-1)*valNumb+v
			            	yPos<-1+i
			                polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col=farben[ID_select[i]])
			                if(signis[sea,ID_select[i],state,1]>0.1){polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col="white",density=50)}	     
			                if(!is.na(notAccepted[sea,ID_select[i],state,1])){
			                	polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col="white")
			                	points(xPos+0.5,yPos+0.5,pch=8,cex=1.5)
			                }	     
			            }
					}
				}
			}
		}

		# border lines
		for (sea in 1:4){lines(c(sea*2*valNumb+1,sea*2*valNumb+1),c(2,length(ID_select)+2),lwd=2)}
		for (i in 1:length(ID_select)){if (ID_select[i] %in% hlines){lines(c(1,valNumb*10+1),c(i+1,i+1),lwd=2)}}
		lines(c(1,valNumb*10+1),c(length(ID_select)+2,length(ID_select)+2),lwd=2)
		lines(c(1,valNumb*10+1),c(2,2),lwd=2)
		lines(c(1,1),c(2,length(ID_select)+2),lwd=2)
		lines(c(valNumb*10+1,valNumb*10+1),c(2,length(ID_select)+2),lwd=2)
		par(new=TRUE)
		plot(NA,xlim=c(0,1),ylim=c(1,0),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
		#image.plot(legend.only=T,horizontal=TRUE, zlim=range(y), col=color,add=TRUE,fill=TRUE,smallplot=c(0.1,0.9,0.5,0.9))
		image.plot(legend.only=T,horizontal=TRUE, zlim=range(y[1:length(y)]), col=color,add=FALSE,fill=TRUE,smallplot=c(0.15,0.93,0.1,0.15))
	}
	graphics.off()
}


plot_all_changes_table <- function(region_name="ward24",period="1950-2014",partPeriods=c("1950-1980","1980-2014"),ID_select=1:24,hlines=c(30),folder=paste("/regional/",region_name,"/",sep="")){
    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    regNumb<-length(unique(attribution[!is.na(attribution)]))

	# shuff
    filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_shuffQuant.nc",sep="") ; print(filename)
	LR<-var.get.nc(nc<-open.nc(filename),"original_slopes")

	# MK
	filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/",trendID,dataset,"_",region_name,"_",period,"_MK.nc",sep="") ; print(filename)
	MK<-var.get.nc(nc<-open.nc(filename),"MK")

	# distr comparison
   	filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,trendID,dataset,"_",region_name,"_dur_ks_wilcox_",partPeriods[1],"_vs_",partPeriods[2],".nc",sep="") ; print(filename)
    ks_test <- var.get.nc(open.nc(filename),"tests")
    fit_params=array(NA,c(2,6,regNumb,2,30))
    quantiles=array(NA,c(2,6,regNumb,2,length(taus)+1,3))
    for (i in 1:length(partPeriods)){
    	partPeriod<-partPeriods[i]
    	filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,partPeriod,"/",trendID,"_",dataset,"_",region_name,"_",partPeriod,"_fit_","2expo_4:100",".nc",sep=""); print(filename)
   	 	fit_params[i,,,,]=var.get.nc(open.nc(filename),"fit_stuff")

        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,partPeriod,"/",trendID,"_",dataset,"_",region_name,"_",partPeriod,"_quantiles",".nc",sep=""); print(filename)
        quantiles[i,,,,1:length(taus),]=var.get.nc(open.nc(filename),"quantile_stuff")

        filename<-paste("../data/",dataset,additional_style,"/",trendID,folder,partPeriod,"/",trendID,"_",dataset,"_",region_name,"_",partPeriod,"_others",".nc",sep=""); print(filename)
        quantiles[i,,,,7,1]=var.get.nc(open.nc(filename),"other_stuff")[,,,1]
	}

	changes<-array(NA,c(5,regNumb,2,40,5))

	changes[,,,1:5,1]=LR[,,,,1]
	changes[,,,1:5,2][which(is.na(LR[,,,,3]))]=1

	changes[,,,11:15,1]=MK[,,,,1]
	changes[,,,11:15,2][which(MK[,,,,2]>0.1)]=1

	changes[,,,21:27,1]=quantiles[2,1:5,,,1:7,1]-quantiles[1,1:5,,,1:7,1]
	for(i in 1:7){changes[,,,(20+i),3]=ks_test[,,,1]}
	for(i in 1:7){changes[,,,(20+i),4]=ks_test[,,,6]}
	changes[,,,20:27,2][which(changes[,,,20:27,3]>0.1)]=1
	changes[,,,20:27,2][which(changes[,,,20:27,4]>0.1)]=1

	changes[,,,31:33,1]=(fit_params[2,1:5,,,c(12,14,15)]-fit_params[1,1:5,,,c(12,14,15)])
	for(i in 1:3){changes[,,,(30+i),3]=ks_test[,,,1]}
	changes[,,,31:33,2][which(changes[,,,31:33,3]>0.1)]=1
	for(i in 1:3){changes[,,,(30+i),3]=ks_test[,,,6]}
	changes[,,,31:33,2][which(changes[,,,31:33,3]>0.1)]=1
	for(i in 1:3){changes[,,,(30+i),3]=fit_params[1,1:5,,,21]}
	changes[,,,31:33,2][which(changes[,,,31:33,3]<0.99)]=1
	for(i in 1:3){changes[,,,(30+i),3]=fit_params[2,1:5,,,21]}
	changes[,,,31:33,2][which(changes[,,,31:33,3]<0.99)]=1

	method_names=c("Regression","Mann-Kendall","Distr. Comparison")
	method_pos=c(2,4,7)
	val_names=c("lr","95","mn","95","mn","95","b1","b2","tr")
	value_auswahl=c(5,3,15,13,27,25,31,32,33)
	valNumb=length(value_auswahl)

	filename_plot<-paste("../plots/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/",trendID,"_",region_name,"_",period,"_allChanges.pdf",sep="") ; print(filename_plot)
	pdf(file=filename_plot,width=valNumb,height=6)
	par(mar=c(3,0,0,0))

	ID_select<-ID_select[length(ID_select):1]

	jet.colors <- colorRampPalette( c("red","yellow","green","blue") )
	farben=c(rgb(1,0.3,1),rgb(1,0.3,1,0.3),rgb(0.3,1,1),rgb(0.3,1,1,0.3))
	farben_fit=c(rgb(0.3,1,0),rgb(0.3,1,0,0.3),rgb(1,0.3,0),rgb(1,0.3,0,0.3))

	for (sea in season_auswahl){
		plot(NA,xlim=c(0,valNumb*2+1),ylim=c(-1,29),frame.plot=FALSE,axes=FALSE,xlab="",ylab="")
		for (i in 1:length(ID_select)){text(x=0.5,y=i+1.5,label=ID_select[i])}
		for (state in 1:2){
			text(x=(state-1)*valNumb+valNumb/2+1,y=length(ID_select)+4.5,label=state_names[state])
			for (met in 1:4){text(x=method_pos[met]+(state-1)*valNumb,y=length(ID_select)+3.5,label=method_names[met],cex=0.8)}
			for (v in 1:valNumb){
				text(x=(state-1)*valNumb+v+0.5,y=length(ID_select)+2.5,label=val_names[v])
			    for (i in 1:length(ID_select)){
			    	q<-ID_select[i]
			    	val<-changes[sea,q,state,value_auswahl[v],1]
			    	sig<-changes[sea,q,state,value_auswahl[v],2]
			    	xPos<-(state-1)*valNumb+v
			        yPos<-1+i
			        if (value_auswahl[v]<30){color<-farben}
			        if (value_auswahl[v]>30){color<-farben_fit}
			    	if (is.na(sig) & val>0){polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col=color[1])}
			    	if (!is.na(sig) & val>0){polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col=color[2])}
			    	if (is.na(sig) & val<0){polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col=color[3])}
			    	if (!is.na(sig) & val<0){polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col=color[4])}
				}
			}
		}

		# border lines
		for (i in 1:length(ID_select)){if (ID_select[i] %in% hlines){lines(c(1,valNumb*2+1),c(i+1,i+1),lwd=2)}}
		lines(c(1,valNumb*2+1),c(length(ID_select)+2,length(ID_select)+2),lwd=2)
		lines(c(1,valNumb*2+1),c(2,2),lwd=2)
		lines(c(1,1),c(2,length(ID_select)+2),lwd=2)
		lines(c(valNumb*2+1,valNumb*2+1),c(2,length(ID_select)+2),lwd=2)
		for (i in c(3,5,7,10,12,14,16)){lines(c(i,i),c(2,length(ID_select)+2),lwd=2)}

		polygon(x=c(1,1,4,4),y=c(-0.5,0.5,0.5,-0.5),border=rgb(1,1,1,0.0),col=farben[1])
		text(x=2.5,y=0,label="sign. increase")
		polygon(x=c(5,5,8,8),y=c(-0.5,0.5,0.5,-0.5),border=rgb(1,1,1,0.0),col=farben[3])
		text(x=6.5,y=0,label="sign. decrease")
		polygon(x=c(1,1,4,4),y=c(-2,-1,-1,-2),border=rgb(1,1,1,0.0),col=farben_fit[1])
		text(x=2.5,y=-1.5,label="sign. increase")
		polygon(x=c(5,5,8,8),y=c(-2,-1,-1,-2),border=rgb(1,1,1,0.0),col=farben_fit[3])
		text(x=6.5,y=-1.5,label="sign. decrease")

		polygon(x=c(11,11,14,14),y=c(-0.5,0.5,0.5,-0.5),border=rgb(1,1,1,0.0),col=farben[2])
		text(x=12.5,y=0,label="increase")
		polygon(x=c(15,15,18,18),y=c(-0.5,0.5,0.5,-0.5),border=rgb(1,1,1,0.0),col=farben[4])
		text(x=16.5,y=0,label="decrease")
		polygon(x=c(11,11,14,14),y=c(-2,-1,-1,-2),border=rgb(1,1,1,0.0),col=farben_fit[2])
		text(x=12.5,y=-1.5,label="increase")
		polygon(x=c(15,15,18,18),y=c(-2,-1,-1,-2),border=rgb(1,1,1,0.0),col=farben_fit[4])
		text(x=16.5,y=-1.5,label="decrease")
	}
	graphics.off()
}
