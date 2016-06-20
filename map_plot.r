
location_view <- function(station=0,lon=0,lat=0,regions=NA){
	# shows region around grid-point
	if (station!=0){
		q=which(dat$ID==station)
	}
	worldmap = getMap(resolution = "low")
	pdf(file="../plots/ID_region_map.pdf",width=12,height=8)
	plot(worldmap)#,xlim=c(-180,-5),ylim=c(35,60), asp = 3.5)
	for (i in 1:length(dat$ID)){
		text(dat$lon[i],dat$lat[i],label=dat$ID[i],col="red",cex=0.125)
	}
    data(topoWorld)
	plot(topoWorld)
	for (i in 1:length(dat$ID)){
		text(dat$lon[i],dat$lat[i],label=dat$ID[i],col="red",cex=0.125)
	}
	if (!is.na(regions)){
		region_names=c("srex","7rect")#,"6wave","7wave","8wave")
		color=c("blue","green","red","orange","black")
	    for (k in 1:length(region_names)){
	    	plot(worldmap)
	    	for (i in 1:length(dat$ID)){
				text(dat$lon[i],dat$lat[i],label=dat$ID[i],col="red",cex=0.125)
			}
	    	add_region(region_names[k],color[k])
	    }
	}
	graphics.off()
}

add_region <- function(region_name,farbe){
	# adds polygon regions to map
    poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))
    for (i in 1:dim(poli)[1]){
        lon=poli[i,1:6]
        lat=poli[i,7:12]
        lon=lon[!is.na(lon)]
        lat=lat[!is.na(lat)]
        polygon(x=lon,y=lat,border=farbe)
    }
    reg_name=read.table(paste("../data/region_poligons/",region_name,"_labels.txt",sep=""))
    for (i in 1:dim(poli)[1]){
        lon=reg_name[i,1]
        lat=reg_name[i,2]
        text(label=reg_name[i,3],x=lon,y=lat,col=farbe,cex=1.5)
    }
}

region_border <- function(ID_select=1:1319,region_name="ward22",reg_select=1:24,reg_names=1:24,border_col="white",text_col=border_col){
    # loads attribution file and visualizes regions on map
    # can be used for regions obtained from clustering
    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    mids<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,"_mids.txt",sep=""))    

    lon<-dat$lon
    lat<-dat$lat
	if (xAusschnitt[2]>180){lon[lon<xAusschnitt[1]]=lon[lon<xAusschnitt[1]]+360}

    regNumb<-dim(mids)[1]
    for (g in 1:length(reg_select)){
    	G<-reg_select[g]
        inside<-which(attribution==G)
        for (q in inside){
            if (!((lon[q]+3.75) %in% lon[inside[which(lat[inside]==lat[q])]])){lines(c(lon[q]+3.75/2,lon[q]+3.75/2),c(lat[q]-1.25,lat[q]+1.25),col=border_col)}
            if (!((lon[q]-3.75) %in% lon[inside[which(lat[inside]==lat[q])]])){lines(c(lon[q]-3.75/2,lon[q]-3.75/2),c(lat[q]-1.25,lat[q]+1.25),col=border_col)}
            if (!((lat[q]+2.5) %in% lat[inside[which(lon[inside]==lon[q])]])){lines(c(lon[q]-3.75/2,lon[q]+3.75/2),c(lat[q]+1.25,lat[q]+1.25),col=border_col)}
            if (!((lat[q]-2.5) %in% lat[inside[which(lon[inside]==lon[q])]])){lines(c(lon[q]-3.75/2,lon[q]+3.75/2),c(lat[q]-1.25,lat[q]-1.25),col=border_col)}
        }
        text(mids[G,2],mids[G,3],label=reg_names[g],col=text_col)
    }
}


put_points <- function(points,points_sig=points*NA,farb_mitte="mean",farb_palette="regenbogen",signi_level=0,i=1,ID_select=1:1319){
	# adds colored grid-points to the map

	if (length(which(!is.na(points[ID_select])))<3){return(list(y=c(NA),color=c(NA)))}
	if (length(unique(points[ID_select]))<3){return(list(y=c(NA),color=c(NA)))}

	y1=points[ID_select]
	sig=points_sig[ID_select]
	pch=pch_points[ID_select]
	sig[(!is.na(sig) & sig>signi_level)]=NA
	sig[(!is.na(sig) & sig<signi_level)]=pch_sig

	notna=which(!is.na(y1))

	y=c(NA,NA,y1[notna])
	sig=sig[notna]

	# depending on color-cheme --------------------------------
	if (length(farb_mitte)>=3){
		if (is.na(farb_mitte[i+2])){farb_mitte_loc=farb_mitte[(i+1)]}
		if (!is.na(farb_mitte[i+2])){farb_mitte_loc=farb_mitte[((i-1)*2+2):((i-1)*2+3)]}
	}
	else {farb_mitte_loc=farb_mitte}

	if ((length(farb_mitte_loc)==2 & farb_mitte_loc[1]!=farb_mitte_loc[2])){
		y[1]=farb_mitte_loc[1]
		y[2]=farb_mitte_loc[2]
		y[y>farb_mitte_loc[2]]=farb_mitte_loc[2]
		y[y<farb_mitte_loc[1]]=farb_mitte_loc[1]
	}	
	if ((length(farb_mitte_loc)==2 & farb_mitte_loc[1]==farb_mitte_loc[2])){
		y[1]=farb_mitte_loc[1]
		y[2]=farb_mitte_loc[1]
		y[y>farb_mitte_loc[1]]=farb_mitte_loc[1]
		y[y<(-farb_mitte_loc[1])]=(-farb_mitte_loc[1])
	}	

	if (length(farb_mitte_loc)==1){
		if (farb_mitte_loc=="cut"){
			mi=mean(y,na.rm=TRUE)
			sd=sd(y,na.rm=TRUE)
			high=mi+1*sd
			low=mi-1*sd
			y[y>high]=high
			y[y<low]=low
			y[1]=low
			y[2]=high
		}

		if (farb_mitte_loc=="gemeinsam 0"){
			y[1]=-aushol
			y[2]=aushol			
		}
		if (farb_mitte_loc=="gemeinsam mean"){
			y[1]=mi-aushol
			y[2]=mi+aushol			
		}
		if (farb_mitte_loc=="0"){
			aushol=max(c(abs(max(y,na.rm=TRUE)),abs(min(y,na.rm=TRUE))))
			y[1]=-aushol
			y[2]=aushol
		}
		if (farb_mitte_loc=="mean"){
			mi=mean(y,na.rm=TRUE)
			y[1]=mi
			y[2]=mi
		}	
		if (farb_mitte_loc=="nichts"){
			y[1]=mean(y,na.rm=TRUE)
			y[2]=mean(y,na.rm=TRUE)
		}
	}


	if ((length(farb_palette)>1 & farb_palette[1]=="mixed")){
		farb_palette_loc=farb_palette[3]
		nbcol=strtoi(farb_palette[2])
	}

	if ((length(farb_palette)>1 & farb_palette[1]=="individual")){
		farb_palette_loc=farb_palette[(i+1)]
	}
	if (length(farb_palette)==1){farb_palette_loc=farb_palette}	


	if (farb_palette_loc=="gold-blau"){
		#jet.colors <- colorRampPalette( c(rgb(0,0.5,1),rgb(0,1,1), rgb(1,1,1) ,rgb(1,1,0),rgb(1,0.5,0)))
		jet.colors <- colorRampPalette( c("blue",rgb(1,1,1),"red"))
		jet.colors <- colorRampPalette( c(rgb(0.2,0.6,0.6),rgb(0.5,1,1), rgb(0.98,0.98,0.98) ,rgb(1,1,0),rgb(0.6,0.6,0)))
		jet.colors <- colorRampPalette( c("blue","cyan",rgb(1,1,1),"magenta","red"))
	}		
	if (farb_palette_loc=="blau-rot"){
		jet.colors <- colorRampPalette( c("blue","white","red"))
	}
	if (farb_palette_loc=="cyan-magenta"){
		jet.colors <- colorRampPalette(c("cyan","white","magenta"))
	}

	if (farb_palette_loc=="drei"){
		jet.colors <- colorRampPalette(c("cyan","green","magenta"))
	}

	if (farb_palette_loc=="lila-gruen"){
		jet.colors <- colorRampPalette( c(rgb(0.1,0.2,0.4),rgb(0.5,1,1),rgb(0.5,1,0.5), rgb(1,1,1),rgb(1,0.7,0.7),rgb(1,0.5,1),rgb(0.4,0.1,0.4)))
	}
	if (farb_palette_loc=="gruen-lila"){
		jet.colors <- colorRampPalette(c(rgb(0.1,0.2,0.4),rgb(0.5,1,1),rgb(0.5,1,0.5), rgb(1,1,1),rgb(1,0.7,0.7),rgb(1,0.5,1),rgb(0.4,0.1,0.4))[7:1])
	}
	if (farb_palette_loc=="rot-gruen"){
		jet.colors <- colorRampPalette(c("green","yellow","red")[3:1])
	}	
	if (farb_palette_loc=="gruen-rot"){
		jet.colors <- colorRampPalette(c("green","yellow","red"))
	}
	if (farb_palette_loc=="cyan-orange"){
		jet.colors <- colorRampPalette(c("cyan","yellow","orange"))
	}
	if (farb_palette_loc=="lila-blau"){
		jet.colors <- colorRampPalette(c(rgb(0.1,0.2,0.4),rgb(0.5,1,1), rgb(1,1,1),rgb(1,0.5,1),rgb(0.4,0.1,0.4))[5:1])
	}	
	if (farb_palette_loc=="blau-lila"){
		jet.colors <- colorRampPalette(c(rgb(0.1,0.2,0.4),rgb(0.5,1,1), rgb(1,1,1),rgb(1,0.5,1),rgb(0.4,0.1,0.4)))
	}	
	if (farb_palette_loc=="regenbogen"){
		jet.colors <- colorRampPalette( c( "blue","green","yellow","red","violet") )
	}
	if (farb_palette_loc=="weiss-rot"){
		jet.colors <- colorRampPalette( c( "white","yellow","red") )
	}
	if (farb_palette_loc=="spacy"){
		jet.colors <- colorRampPalette(c("blue","red","green","yellow","orange","violet"))
	}
	if (farb_palette_loc=="niederschlag"){
		jet.colors <- colorRampPalette( c( "orange","yellow","green","blue") )
	}

	if (farb_palette_loc=="wave-goggle"){
		jet.colors <- colorRampPalette( c(rgb(1,1,1,0),rgb(0.8,1,0.1),rgb(0,1,1),rgb(1,0.5,0)))
	}

	if (farb_palette_loc=="norm-hot"){
		jet.colors <- colorRampPalette( c(rgb(1,1,1,0),rgb(1,1,0),rgb(1,0.5,0),rgb(1,0,0),rgb(1,0,0.5)))
	}

	if (farb_palette_loc=="norm-wet"){
		jet.colors <- colorRampPalette( c(rgb(1,1,1,0),"lightblue",rgb(0,1,0.5),rgb(0,1,1),rgb(0,0,1)) )
	}

	if (farb_palette_loc=="dry-wet"){
		jet.colors <- colorRampPalette( c( "orange","yellow","white","green","blue") )
	}

	if (farb_palette_loc=="viele"){
		jet.colors <- colorRampPalette(c(rgb(0.5,0.5,0.5),rgb(1,0.5,0.5,0.5),"black",rgb(0.8,0.5,1),rgb(1,0,0),rgb(0.2,1,0.5),"yellow",rgb(0.5,0.5,1),rgb(1,0.6,1),rgb(0.5,1,1),rgb(0.2,0.5,1)))
		jet.colors <- colorRampPalette(brewer.pal(n=9,name="Set1"))
	}
	if (farb_palette_loc=="groups"){
		palette=array(c("black","blue","green","yellow","orange","red",rgb(0.5,0.5,0.5,0.5)),nbcol)
		jet.colors <- colorRampPalette( palette)
	}
	color <- jet.colors(nbcol)	
	facetcol <- cut(y,nbcol)

	farben=color[facetcol[3:(length(notna)+2)]]

	lon=dat$lon[ID_select][notna]
	lat=dat$lat[ID_select][notna]


	if (xAusschnitt[2]>180){lon[lon<xAusschnitt[1]]=lon[lon<xAusschnitt[1]]+360}

	#delete out of yAusschnitt
	#ID_select_local=which(lat >= yAusschnitt[1] & lat <= yAusschnitt[2] & lon >= xAusschnitt[1] & lon <= xAusschnitt[2])
	ID_select_local<-1:length(lat)

	if (!is.na(pch_points[2])){
		points(lon[ID_select_local],lat[ID_select_local],pch=pch[ID_select_local],col=farben[ID_select_local],cex=pointsize)#
		points(lon[ID_select_local],lat[ID_select_local],pch=sig[ID_select_local],cex=pointsize,col=col_sig)
	}

	if (is.na(pch_points[2])){
		for (q in 1:length(ID_select_local)){
			polygon(x=c(lon[q]-pch_points[3],lon[q]+pch_points[3],lon[q]+pch_points[3],lon[q]-pch_points[3]),y=c(lat[q]-pch_points[4],lat[q]-pch_points[4],lat[q]+pch_points[4],lat[q]+pch_points[4]),border=rgb(1,1,1,0.0),col=farben[q])
			if (!is.na(sig[q])){
				points(lon[q],lat[q],pch=3,cex=pointsize)
				#lines(c(lon[q]-pch_points[3],lon[q]),c(lat[q]-pch_points[4],lat[q]+pch_points[4]),col=col_sig,cex=cex_sig)
				#lines(c(lon[q],lon[q]+pch_points[3]),c(lat[q]-pch_points[4],lat[q]+pch_points[4]),col=col_sig,cex=cex_sig)
			}
		}
		#points(lon[ID_select],lat[ID_select],pch=sig[ID_select],cex=pointsize,col=col_sig)
	}

	return(list(y=y,color=color))
}

topo_map_plot <- function(filename_plot=filename_plot,reihen=reihen,reihen_sig=reihen*NA,titel=c(""),signi_level=0.05,farb_mitte="mean",farb_palette="regenbogen",regionColor="black",average=FALSE,cex=1,color_lab="",cex_axis=1,highlight_points=c(NA),highlight_color=c(NA),main="",ID_select=1:length(dat$ID),region_name=NA,regNumb=NA){
	# create a nice map on which can be plotted
	
	pdf(file=filename_plot,width=paper[1],height=paper[2])	

	for (i in 1:dim(reihen)[1]){
		par(mar=margins)
		plot(NA,xlim=xAusschnitt,ylim=yAusschnitt,axes=FALSE,frame.plot=FALSE,asp=asp)
		

	    if (titel[1]!=""){
	    	main<-titel[i]
	    	print(titel[i])
	    	#text(0,-85,main)
	    }

	    # color land in grey for no coverage
	    if (greyLand==TRUE){
		    lon<-allPoints$Var1
		    lat<-allPoints$Var2
			for (q in 1:7081){
				polygon(x=c(lon[q]-pch_points[3],lon[q]+pch_points[3],lon[q]+pch_points[3],lon[q]-pch_points[3]),y=c(lat[q]-pch_points[4],lat[q]-pch_points[4],lat[q]+pch_points[4],lat[q]+pch_points[4]),border=rgb(1,1,1,0.0),col=c("white","gray")[land[q]])
			}
		}

	    #data points
	    tmp=put_points(points=reihen[i,],points_sig=reihen_sig[i,],signi_level=signi_level,i=i,farb_mitte=farb_mitte,farb_palette=farb_palette,ID_select=ID_select)

	    #caostline
	    worldCoast<-map(interior=FALSE,resolution=0,plot=FALSE,xlim=xAusschnitt,ylim=yAusschnitt)
	    lines(worldCoast)

	    #highlight points
		for (rad in c(1,1.5,2,2.5)){
			points(dat$lon[highlight_points[i]],dat$lat[highlight_points[i]],col=highlight_color,pch=1,cex=(pointsize*rad))
		}

		#region contours
		if (!is.na(region_name)){
			tmp=put_points(points=reihen[i,],points_sig=reihen_sig[i,],signi_level=signi_level,i=i,farb_mitte=farb_mitte,farb_palette=farb_palette,ID_select=ID_select)
			region_border(region_name=region_name,border_col="black",reg_select=c(2,6,10,13,19,23,1,3,4,7,12,16,20,5,11,14,18,21,22),reg_names=c(""))
		}
		if (!is.na(region)){region_border(region_name=region,border_col=border_col)}

		# clear islands for index
		polygon(x=c(-180,-150,-150,-180),y=c(-30,-30,30,30),border="white",col="white")
		#index topright
		if (length(indexTopRight>=dim(reihen)[1]) & !is.na(indexTopRight[i])){text(posTopRight[1],posTopRight[2],indexTopRight[i],pos=4,cex=cexIndex)}
		#index lowleft
		if (length(indexBottomLeft>=dim(reihen)[1]) & !is.na(indexBottomLeft[i])){text(posBottomLeft[1],posBottomLeft[2],indexBottomLeft[i],pos=4,cex=cexIndex)}

		# delete outer coastline
		if (!is.na(outer_cut)){polygon(x=c(xAusschnitt[1]-outer_cut,xAusschnitt[2]+outer_cut,xAusschnitt[2]+outer_cut,xAusschnitt[1]-outer_cut),y=c(yAusschnitt[1]-outer_cut,yAusschnitt[1]-outer_cut,yAusschnitt[2]+outer_cut,yAusschnitt[2]+outer_cut),col=rgb(1,1,1,0),border="white",lwd=inner_cut)}

		#color bar
		color<<-tmp$color
		y=tmp$y[]
	    if (!is.na(color[1]) & color_legend=="right"){image.plot(legend.only=T,horizontal=FALSE, zlim=range(y), col=color,add=FALSE,legend.lab=color_lab)}
	}
	if (color_legend=="seperate"){
		plot(NA,xlim=c(0,1),ylim=c(1,0),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
		image.plot(legend.only=T,horizontal=TRUE, zlim=range(y), col=color,add=TRUE,fill=TRUE,smallplot=c(0.1,0.9,0.8,0.9),cex=3)
	}
	if (!is.na(color[1])){
		#par(new=TRUE)
		plot(NA,xlim=c(0,1),ylim=c(1,0),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
		image.plot(legend.only=T,horizontal=TRUE, zlim=range(y), col=color,add=TRUE,fill=TRUE,smallplot=c(0.1,0.9,0.5,0.9))
	}
	if (closePlot==TRUE){graphics.off()}
}

library(RColorBrewer)
library(rworldmap)
library(fields)
require(maps)
library(maptools)
data(wrld_simpl)


allPoints <- expand.grid(seq(-180,180,3.75), seq(-90,90,2.5))
allPoints <- SpatialPoints(allPoints, proj4string=CRS(proj4string(wrld_simpl)))
land<-array(1,7081)
land[which(!is.na(over(allPoints, wrld_simpl)$FIPS))]=2

