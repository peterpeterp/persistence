
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


plot_reg_table_general <- function(values,signis,filename_plot,val_names,region_name="ward24",colorRange=c(0,0.4),farb_palette="lila-gruen",ID_select=1:24,hlines=c(30),season_auswahl=c(1:5),style="range"){
    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    regNumb<-length(unique(attribution[!is.na(attribution)]))
	valNumb<-dim(values)[4]


	pdf(file=filename_plot,width=valNumb*3,height=6)
	par(mar=c(3,0,0,0))

	ID_select<-ID_select[length(ID_select):1]

	if (farb_palette=="lila-gruen"){jet.colors <- colorRampPalette( c(rgb(0.5,1,1),rgb(0.5,1,0.5), rgb(1,1,1),rgb(1,0.7,0.7),rgb(1,0.5,1)))}
	if (farb_palette=="lila-gruen-inv"){jet.colors <- colorRampPalette( c(rgb(0.5,1,1),rgb(0.5,1,0.5), rgb(1,1,1),rgb(1,0.7,0.7),rgb(1,0.5,1))[5:1])}
	color <- jet.colors(nbcol)	

	plot(NA,xlim=c(0,valNumb*10+1),ylim=c(0,29),frame.plot=FALSE,axes=FALSE,xlab="",ylab="")
	for (i in 1:length(ID_select)){text(x=0.5,y=i+1.5,label=ID_select[i])}
	for (sea in season_auswahl){
		if (valNumb>1){text(x=(sea-1)*2*valNumb+3,y=length(ID_select)+4.5,label=season_names[sea])}
		if (valNumb==1){text(x=(sea-1)*2*valNumb+2,y=length(ID_select)+3.5,label=season_names[sea])}
		for (state in 1:2){
			if (valNumb>1){text(x=(sea-1)*2*valNumb+(state-1)*valNumb+2,y=length(ID_select)+3.5,label=state_names[state])}
			if (valNumb==1){text(x=(sea-1)*2*valNumb+(state-1)*valNumb+1.5,y=length(ID_select)+2.5,label=c("c","w")[state])}
			for (v in 1:valNumb){
				text(x=(sea-1)*2*valNumb+(state-1)*valNumb+v+0.5,y=length(ID_select)+2.5,label=val_names[v])
				val<-v
				y<-c(values[sea,,state,val],colorRange[1],colorRange[2])
				y[y>colorRange[2]]=colorRange[2]
				y[y<colorRange[1]]=colorRange[1]
				facetcol <- cut(y,nbcol)
				#print(facetcol)
				farben<-color[facetcol]
			    for (i in 1:length(ID_select)){
			        xPos<-(sea-1)*2*valNumb+(state-1)*valNumb+v
			        yPos<-1+i
			        polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col=farben[ID_select[i]])
			        if(!is.na(signis[sea,ID_select[i],state,1,val])){polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col="white",density=30)}	     
			        if(!is.na(signis[sea,ID_select[i],state,2,val])){
			            #polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col="white")
			            polygon(x=c(xPos,xPos+1,xPos+1,xPos),y=c(yPos,yPos,yPos+1,yPos+1),border=rgb(1,1,1,0.0),col="black",density=30)
			            #points(xPos+0.5,yPos+0.5,pch=8,cex=1.5)	     
					}
			        if(!is.na(signis[sea,ID_select[i],state,3,val])){
			            points(xPos+0.5,yPos+0.5,pch=17,cex=1.5,col="white")	     
					}
			        if(!is.na(signis[sea,ID_select[i],state,4,val])){
			            points(xPos+0.5,yPos+0.5,pch=20,cex=1.5,col="black")	     
					}			        
					if(!is.na(signis[sea,ID_select[i],state,5,val])){
			            points(xPos+0.5,yPos+0.5,pch=6,cex=1.5,col="black")	     
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

	if (nbcol>5){image.plot(legend.only=T,horizontal=TRUE, zlim=range(y[1:length(y)]), col=color,add=FALSE,fill=TRUE,smallplot=c(0.15,0.93,0.1,0.15))}

	graphics.off()
}
