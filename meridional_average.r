


meridional_average <- function(dat,climatology,change,titel,filename_plot){
	pdf(file=filename_plot)
	par(mar=c(5, 4, 4, 6) + 0.1)
	for (i in 1:dim(climatology)[1]){
		x=c(seq(0,180,3.75),seq(-180,0,3.75))
		y1=x*NA
		y2=x*NA
		for (k in 1:length(x)){
			print(k)
			drauf=which(dat$lon==x[k] & dat$lat<70 & dat$lat>50 & !is.na(climatology[i,]))
			print(drauf)
			print(climatology[i,drauf])
			y1[k]=mean(climatology[i,drauf])
			y2[k]=mean(change[i,drauf])

		}
		x[50:98]=x[50:98]+360
		plot(x,y1-mean(y1,na.rm=TRUE),main=titel[i],axes=FALSE,xlab="",ylab="")
		axis(1,col="black")
		mtext("Longitude",side=1,col="black",line=2.5)  
		axis(2,col="black",las=1)
		mtext("mean",side=2,line=2.5)
		par(new = TRUE)
		y2=y2-mean(y2,na.rm=TRUE)
		toPlot=which(!is.na(y2))
		plot(x[toPlot], y2[toPlot],col="red", type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")
		mtext("change",side=4,col="red",line=4) 
		axis(4, col="red",col.axis="red",las=1)
	}
	graphics.off()
}