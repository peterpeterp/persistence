source("write.r")
source("load.r")
source("functions_support.r")

if (1==1){
	nc=open.ncdf("../data/91_5/91_5_markov3s.nc")
	dat=dat_load("../data/dat_regional.nc")
	q=554	
	pdf(file=paste("../plots/3states_",q,".pdf",sep=""))
	tmp=get.var.ncdf(nc,"markov_summer")
	summer=array(tmp[q,,],dim=c(3,3,62))

	states=c("cold","normal","warm")
	label=c()

	jet.colors <- colorRampPalette( c( "violet","blue","green","yellow","red") )
	nbcol <- 9
	color <- jet.colors(nbcol)

	plot(NA,xlim=c(1950,2011),ylim=c(0,1))
	for (from in 1:3){
		for (to in 1:3){
			lines(dat$year,summer[from,to,],col=color[(from-1)*3+to])
			#print(summary(lm(summer[from,to,]~dat$year)))
			label[(from-1)*3+to]=paste(states[from],"to",states[to],summary(lm(summer[from,to,]~dat$year))$coefficients[2])
			abline(lm(summer[from,to,]~dat$year),col=color[(from-1)*3+to])
		}
	}
	abline(v=c(2003,2006))

	legend("topleft",lty=array(1,9),col=color,legend=label)

}

