
location_finder <- function(dat,station=0,lon=0,lat=0){
	library(rworldmap)
	library(fields)
	worldmap = getMap(resolution = "low")
	pdf(file="../plots/location.pdf")
	if (station!=0){
		q=which(dat$ID==station)
		plot(worldmap,xlim=c(dat$lon[q]-10,dat$lon[q]+10),ylim=c(dat$lat[q]-10,dat$lat[q]+10))
		points(dat$lon[q],dat$lat[q],pch=15,col="red")
	}
	else {
		
		plot(worldmap,xlim=c(lon-20,lon+20),ylim=c(lat-20,lat+20))
		points(lon,lat,pch=15,col="red")
		for (q in 1:1319){
			text(dat$lon[q],dat$lat[q],dat$ID[q],cex=0.7)
		}
	}
}

location_view <- function(station=0,lon=0,lat=0,regions=NA){
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

region_border <- function(ID_select=1:1319,region_name="ward22",border_col="white",text_col=border_col){
    # loads attribution file and visualizes regions on map
    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    mids<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,"_mids.txt",sep=""))    

    regNumb<-dim(mids)[1]
    for (G in 1:regNumb){
        inside<-which(attribution==G)
        for (q in inside){
            if (!((dat$lon[q]+3.75) %in% dat$lon[inside[which(dat$lat[inside]==dat$lat[q])]])){lines(c(dat$lon[q]+3.75/2,dat$lon[q]+3.75/2),c(dat$lat[q]-1.25,dat$lat[q]+1.25),col=border_col)}
            if (!((dat$lon[q]-3.75) %in% dat$lon[inside[which(dat$lat[inside]==dat$lat[q])]])){lines(c(dat$lon[q]-3.75/2,dat$lon[q]-3.75/2),c(dat$lat[q]-1.25,dat$lat[q]+1.25),col=border_col)}
            if (!((dat$lat[q]+2.5) %in% dat$lat[inside[which(dat$lon[inside]==dat$lon[q])]])){lines(c(dat$lon[q]-3.75/2,dat$lon[q]+3.75/2),c(dat$lat[q]+1.25,dat$lat[q]+1.25),col=border_col)}
            if (!((dat$lat[q]-2.5) %in% dat$lat[inside[which(dat$lon[inside]==dat$lon[q])]])){lines(c(dat$lon[q]-3.75/2,dat$lon[q]+3.75/2),c(dat$lat[q]-1.25,dat$lat[q]-1.25),col=border_col)}
        }
        text(mids[G,2],mids[G,3],label=G,col=text_col)
    }
}


put_points <- function(points,points_sig=points*NA,farb_mitte="mean",farb_palette="regenbogen",signi_level=0,i=1,ID_select=1:1319){
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
		jet.colors <- colorRampPalette( c(rgb(0.2,0.6,0.6),rgb(0.5,1,1), rgb(0.98,0.98,0.98) ,rgb(1,1,0),rgb(0.6,0.6,0)))
	}		
	if (farb_palette_loc=="blau-rot"){
		jet.colors <- colorRampPalette( c("blue","white","red"))
	}
	if (farb_palette_loc=="lila-gruen"){
		jet.colors <- colorRampPalette( c(rgb(0.5,1,0.5),rgb(0.2,0.6,0.2), rgb(0.0,0.0,0.0),rgb(0.6,0.2,0.6),rgb(1,0.5,1)))
		jet.colors <- colorRampPalette( c(rgb(0.5,1,1),rgb(0.5,1,0.5), rgb(0.0,0.0,0.0),rgb(1,0.5,1),rgb(1,0.7,0.7)))
		jet.colors <- colorRampPalette( c(rgb(0.1,0.2,0.4),rgb(0.5,1,1),rgb(0.5,1,0.5), rgb(1,1,1),rgb(1,0.7,0.7),rgb(1,0.5,1),rgb(0.4,0.1,0.4)))
	}
	if (farb_palette_loc=="regenbogen"){
		jet.colors <- colorRampPalette( c( "blue","green","yellow","red") )
	}
	if (farb_palette_loc=="spacy"){
		#jet.colors <- colorRampPalette(c("blue",rgb(0.1,0.2,0.4),"green",rgb(0.5,1,1),rgb(0.5,1,0.5), "yellow",rgb(1,0.7,0.7),rgb(1,0.5,1),"orange",rgb(0.4,0.1,0.4),"red"))
		jet.colors <- colorRampPalette(c("blue","red","green","yellow","orange","violet"))
	}
	if (farb_palette_loc=="niederschlag"){
		jet.colors <- colorRampPalette( c( "orange","yellow","green","blue") )
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
	print(y)
	color <- jet.colors(nbcol)	
	facetcol <- cut(y,nbcol)

	farben=color[facetcol[3:(length(notna)+2)]]

	lon=dat$lon[notna]
	lat=dat$lat[notna]
	
	if (xAusschnitt[2]>180){lon[lon<xAusschnitt[1]]=lon[lon<xAusschnitt[1]]+360}

	#delete out of yAusschnitt
	ID_select=which(lat >= yAusschnitt[1] & lat <= yAusschnitt[2] & lon >= xAusschnitt[1] & lon <= xAusschnitt[2])

	if (!is.na(pch_points[2])){
		points(lon[ID_select],lat[ID_select],pch=pch[ID_select],col=farben[ID_select],cex=pointsize)#
		points(lon[ID_select],lat[ID_select],pch=sig[ID_select],cex=pointsize,col=col_sig)
	}

	if (is.na(pch_points[2])){
		for (q in ID_select){
			polygon(x=c(lon[q]-pch_points[3],lon[q]+pch_points[3],lon[q]+pch_points[3],lon[q]-pch_points[3]),y=c(lat[q]-pch_points[4],lat[q]-pch_points[4],lat[q]+pch_points[4],lat[q]+pch_points[4]),border=rgb(1,1,1,0.0),col=farben[q])
			if (!is.na(sig[q])){
				lines(c(lon[q]-pch_points[3],lon[q]),c(lat[q]-pch_points[4],lat[q]+pch_points[4]),col=col_sig,cex=cex_sig)
				lines(c(lon[q],lon[q]+pch_points[3]),c(lat[q]-pch_points[4],lat[q]+pch_points[4]),col=col_sig,cex=cex_sig)
			}
		}
		#points(lon[ID_select],lat[ID_select],pch=sig[ID_select],cex=pointsize,col=col_sig)
	}

	return(list(y=y,color=color))
}

map_allgemein <- function(filename_plot=filename_plot,reihen=reihen,reihen_sig=reihen*NA,titel=c(""),signi_level=0.05,
	farb_mitte="mean",farb_palette="regenbogen",region=NA,regionColor="black",average=FALSE,pointsize=1.2,grid=FALSE,ausschnitt=c(-80,80),col_row=c(1,1),paper=c(12,8),cex=1,color_lab="",cex_axis=1,highlight_points=c(NA),highlight_color=c(NA),mat=c(NA),subIndex=c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"),layout_mat=c(NA),topo=FALSE){
	#dat data form data_load()
	#filename_plot str - where to save plot
	#worldmap background of lon lat plot
	#ausschnitt c(lat_min,lat_max)
	#reihen array(... dim=c(anzahl der plots, anzahl der stationen))
	#titel liste von strings, plot-titles
	#farb_mitte mid point of color range (white) at 0 for "0" or at the mean for "mean"

	if (ausschnitt[1]!=-80 & col_row[1]==1){
		paper[2]=paper[2]*(ausschnitt[2]-ausschnitt[1])/160+1
	}

		#if (col_row[1]>1){
		#	paper=c(12,((col_row[1]-1)*7/5+1))
		#	cex_axis=1.5
		#}
	pdf(file = filename_plot,width=paper[1],height=paper[2])
	par(mar=c(1,1,2,4))
	par(cex.axis=cex_axis,cex.lab=cex_axis)


	# create layout matrix for multiplots in style c(1,2,1,2,1,2,3,4,3,4,3,4 ..)
	if (col_row[1]>1 & col_row[2]>1 & is.na(mat[1])){
		par(cex=cex)
		pointsize=1
		mat=c()
		index=0
		for (row in 1:(col_row[1]-1)){
			index=index+1
			mat[((row-1)*6+1):(row*6)]=c(index,index+1,index,index+1,index,index+1)
			index=index+1
		}

		mat[((col_row[1]*col_row[2]-2)*3+1):((col_row[1]*col_row[2]-2)*3+4)]=c(index+1,index+1,index+1,index+1)#,index+1,index+1)
		layout(matrix(mat,length(mat)/2,2, byrow = TRUE))
	}

	if (col_row[1]>1 & col_row[2]==1 & is.na(mat[1])){
		par(cex=cex)
		pointsize=1
		mat=c()
		for (row in 1:(col_row[1])){
			mat[((row-1)*7+1):(row*7)]=c(row,row,row,row,row,row,col_row[1]+1)
		}
		layout(matrix(mat,col_row[1],length(mat)/col_row[1], byrow = TRUE))
	}

	if (col_row[1]==1 & dim(reihen)[1]==1 & is.na(mat[1]) & !topo){
		par(cex=cex)
		pointsize=pointsize
		layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,2,3),1,13, byrow = TRUE))
	}
	if (col_row[1]==1 & dim(reihen)[1]>1){
		pointsize=pointsize
		par(mfrow=c(1,1))
	}

	   # use given mat
	if (!is.na(mat[1])){
		par(cex=cex)
		par(mar=c(1,1,1,1))
		pointsize=pointsize
		print(length(mat))
		layout(matrix(mat,col_row[1],col_row[2], byrow = TRUE), heights=c(2,2,2,2,1))#, heights=c(2,2,2,2,1)
	}

	ID_select=which(dat$lat >= ausschnitt[1] & dat$lat <= ausschnitt[2])

	# if same color range scheme for all plots
	if (farb_mitte[1]=="gemeinsam 0"){
		aushol=max(c(abs(max(reihen[1:dim(reihen)[1],ID_select],na.rm=TRUE)),abs(min(reihen[1:dim(reihen)[1],ID_select],na.rm=TRUE))))
	}

	if (farb_mitte[1]=="gemeinsam mean"){	
		mi=mean(reihen[1:dim(reihen)[1],ID_select],na.rm=TRUE)
		aushol=max(c(abs(max(reihen[1:dim(reihen)[1],ID_select],na.rm=TRUE))-mi,mi-abs(min(reihen[1:dim(reihen)[1],ID_select],na.rm=TRUE))))
	}

	subCount=0
	for (i in 1:dim(reihen)[1]){
		subCount=subCount+1
		if (titel[1]!=""){print(titel[i])}
		if (!topo){
			if (titel[1]==""){plot(worldmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5)}
			else{plot(worldmap,ylim=c(ausschnitt[1],ausschnitt[2]), asp = 1.5, main=titel[i])}
		}
		if (topo){
			if (titel[1]==""){plot(topoWorld,location="none",ylim=c(ausschnitt[1],ausschnitt[2]),xlim=c(0,360),col.land="white",col.water="white",frame.plot=FALSE)}
			else{plot(topoWorld,location="none",ylim=c(ausschnitt[1],ausschnitt[2]),xlim=c(-180,180),col.land="white",main=titel[i])}
		}
		tmp=put_points(points=reihen[i,],points_sig=reihen_sig[i,],ausschnitt=ausschnitt,farb_mitte=farb_mitte,farb_palette=farb_palette,signi_level=signi_level,i=i,pointsize=pointsize)
		color=tmp$color
		y=tmp$y

		for (rad in c(1,1.5)){
			points(dat$lon[highlight_points[i]],dat$lat[highlight_points[i]],col=highlight_color,pch=1,cex=(pointsize*rad))
		}
		#box("figure", col="blue") 

		if (average==TRUE){
			text(-165,ausschnitt[1]+10,paste("mean:",round(mean(y,na.rm=TRUE),02)))
			text(-165,ausschnitt[1]+5,paste("sd:",round(sd(y,na.rm=TRUE),02)))
		}
		#mark nas
		#points(dat$lon,dat$lat,pch=nas,cex=pointsize)

		if (grid==TRUE){
			for (longi in seq(-180,180,30)){
				abline(v=longi,col="grey")
				text(longi,-80,label=longi)
			}
			for (lati in seq(-60,60,10)){
				abline(h=lati,col="grey")
				text(-190,lati,label=lati)
			}
		}
		if (!is.na(region)){
			add_region(region,regionColor)
		}
		#print(y)

		if (col_row[1]>1 & col_row[2]>1){
			legend("topright",legend=c(subIndex[i]),bty="n",cex=cex_axis)
			if (subCount==dim(reihen)[1]){
				plot(NA,xlim=c(0,1),ylim=c(1,0),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
				image.plot(legend.only=T,horizontal=TRUE, zlim=range(y), col=color,add=TRUE,fill=TRUE,legend.mar=8,smallplot=c(0.1,0.9,0.85,0.95),legend.lab=color_lab)
			}
		}
		if (col_row[1]>1 & col_row[2]==1){
			legend("topright",legend=c(subIndex[i]),bty="n",cex=cex_axis)
			if (subCount==dim(reihen)[1]){
				plot(NA,xlim=c(0,1),ylim=c(1,0),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
				image.plot(legend.only=T,horizontal=TRUE, zlim=range(y), col=color,add=FALSE,fill=TRUE,smallplot=c(0.1,0.9,0.6,0.80))
			}
		}
		if (col_row[1]==1 & dim(reihen)[1]==1){
			par(mar=c(1,0,1,0))
			legend("topright",legend=c(subIndex[i]),bty="n",cex=cex_axis)
			plot(NA,xlim=c(0,1),ylim=c(1,0),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
			image.plot(legend.only=T,horizontal=FALSE, zlim=range(y), col=color,add=TRUE,fill=TRUE,smallplot=c(0.1,0.2,0.1,0.90))
			plot(NA,xlim=c(0,1),ylim=c(0,1),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
			text(0.3,0.5,label=color_lab,cex=1,srt=90)
		}
		if (col_row[1]==1 & dim(reihen)[1]>1){
			#image.plot(legend.only=T, zlim=range(y), col=color,add=TRUE,smallplot=c(0.97,0.99,0.1,0.90),legend.lab=color_lab)
			image.plot(legend.only=T, zlim=range(y), col=color,add=TRUE,legend.lab=color_lab)
		}
	}
    graphics.off()
}


topo_map_plot <- function(filename_plot=filename_plot,reihen=reihen,reihen_sig=reihen*NA,titel=c(""),signi_level=0.05,farb_mitte="mean",farb_palette="regenbogen",regionColor="black",average=FALSE,cex=1,color_lab="",cex_axis=1,highlight_points=c(NA),highlight_color=c(NA),main="",ID_select=1:length(dat$ID),region_name=NA,regNumb=NA,border_col="black",land_col=rgb(0.5,0.5,0.5,0.5),water_col=rgb(0,0.2,0.8,0.0)){
	
	pdf(file=filename_plot,width=paper[1],height=paper[2])	

	# use given mat
	if (!is.na(layout_mat[1])){
		par(cex=cex)
		print(length(layout_mat))
		layout(matrix(layout_mat,col_row[1],col_row[2], byrow = TRUE), heights=c(2,2,2,2,1))#, heights=c(2,2,2,2,1)
	}

	if (is.na(layout_mat[1])){par(mfrow=c(1,1))}
	for (i in 1:dim(reihen)[1]){
	    #mapPlot(coastlineWorld,latitudelim=yAusschnitt,longitudelim=xAusschnitt,proj="stereographic",axes=FALSE, fill='white',mar=c(2,3,2,5),main=titel[i])
	    plot(topoWorld,xlim=xAusschnitt,ylim=yAusschnitt,asp=asp,location="none",col.land=rgb(0,0,0,0),col.water=rgb(0,0,0,0),mar=margins,main=titel[i])
	    if (titel[1]!=""){
	    	main<-titel[i]
	    	print(titel[i])
	    	#text(0,-85,main)
	    }
	    axis(2,at=seq(yAusschnitt[1],yAusschnitt[2],10))

	    #data points
	    tmp=put_points(points=reihen[i,],points_sig=reihen_sig[i,],signi_level=signi_level,i=i,farb_mitte=farb_mitte,farb_palette=farb_palette,ID_select=ID_select)
	    #highlight points
		for (rad in c(1,1.5,2,2.5)){
			points(dat$lon[highlight_points[i]],dat$lat[highlight_points[i]],col=highlight_color,pch=1,cex=(pointsize*rad))
		}

		#region contours
		if (!is.na(region_name)){region_border(region_name=region_name,border_col=border_col)}
		if (!is.na(region)){region_border(region_name=region,border_col=border_col)}

		# contour lines topography
		if (!is.na(land_col)){
		    par(new=TRUE)
		    try(plot(topoWorld,xlim=xAusschnitt,ylim=yAusschnitt,asp=asp,location="none",col.land=land_col,col.water=water_col,mar=margins),silent=TRUE)
		}
		#draw over axes

		polygon(x=c(-200,-200,200,200),y=c(yAusschnitt[1]-outer_cut,yAusschnitt[1]+inner_cut,yAusschnitt[1]+inner_cut,yAusschnitt[1]-outer_cut),col="white",border="white")
		polygon(x=c(-200,-200,200,200),y=c(yAusschnitt[2]-inner_cut,yAusschnitt[2]+outer_cut,yAusschnitt[2]+outer_cut,yAusschnitt[2]-inner_cut),col="white",border="white")
		polygon(x=c(xAusschnitt[1]-outer_cut,xAusschnitt[1]-outer_cut,xAusschnitt[1]+inner_cut,xAusschnitt[1]+inner_cut),y=c(yAusschnitt[1],yAusschnitt[2],yAusschnitt[2],yAusschnitt[1]),col="white",border="white")
		polygon(x=c(xAusschnitt[2]-inner_cut,xAusschnitt[2]-inner_cut,xAusschnitt[2]+outer_cut,xAusschnitt[2]+outer_cut),y=c(yAusschnitt[1],yAusschnitt[2],yAusschnitt[2],yAusschnitt[1]),col="white",border="white")


		#color bar
		color=tmp$color
		y=tmp$y[]
	    if (is.na(layout_mat[1]) & !is.na(color[1]) & color_legend=="right"){image.plot(legend.only=T,horizontal=FALSE, zlim=range(y), col=color,add=FALSE,legend.lab=color_lab)}
	}
	if (color_legend=="seperate"){
		plot(NA,xlim=c(0,1),ylim=c(1,0),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
		image.plot(legend.only=T,horizontal=TRUE, zlim=range(y), col=color,add=TRUE,fill=TRUE,smallplot=c(0.1,0.9,0.8,0.9),cex=3)
	}
	if (!is.na(layout_mat[1]) & !is.na(color[1])){
		#par(new=TRUE)
		plot(NA,xlim=c(0,1),ylim=c(1,0),ylab="",xlab="",frame.plot=FALSE,axes=FALSE)
		image.plot(legend.only=T,horizontal=TRUE, zlim=range(y), col=color,add=TRUE,fill=TRUE,smallplot=c(0.1,0.9,0.5,0.9))
	}
	graphics.off()
}

library(oce)
data(topoWorld)
data(coastlineWorld)
library(RColorBrewer)
library(rworldmap)
library(fields)