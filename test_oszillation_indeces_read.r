
create_index <- function(nameGr="NAO",nameKl="nao"){
	source("load.r")

	print(paste("../data/roh_data/",nameGr,".txt",sep=""))
	index=read.table(paste("../data/roh_data/",nameGr,".txt",sep=""))

	index=array(index[,3],dim=c(12,66))
	index[8:12,66]=NA

	index_sea=array(NA,dim=c(5,65))
	for (year in 1:65){
		index_sea[1,year]=mean(index[3:5,year],na.rm=TRUE)
		index_sea[2,year]=mean(index[6:8,year],na.rm=TRUE)
		index_sea[3,year]=mean(index[9:11,year],na.rm=TRUE)
		index_sea[4,year]=mean(c(index[12,year],index[1:2,(year+1)]),na.rm=TRUE)
		index_sea[5,year]=mean(index[1:12,year],na.rm=TRUE)
	}
	#write.table(index_sea,paste("../data/roh_data/",nameGr,"_sea.txt",sep=""))

	years <- dim.def.ncdf("years",units="1",vals=1:65,unlim=FALSE)
	seasons <- dim.def.ncdf("seasons",units="spri,sum,aut,win,yea",vals=1:5, unlim=FALSE)

	nao <- var.def.ncdf(name=nameKl,units="bla",longname=paste(nameGr),dim=list(seasons,years), missval=-9999.0)

	vars=list(nao)
	   
	nc = create.ncdf(paste("../data/roh_data/",nameGr,"_sea.nc",sep=""),vars)

	put.var.ncdf(nc,nameKl,index_sea)  


	close.ncdf(nc) 
}


if (1==2){
	pdf(file="../plots/NAO.pdf")
	x=1:(66*12)
	plot(x,as.vector(index),cex=0.5,xlim=c(50,100))
	x=seq(4,(64*12+4),12)
	points(x,index_sea[1,],col="green")
	x=seq(7,(64*12+7),12)
	points(x,index_sea[2,],col="red")
	x=seq(10,(64*12+10),12)
	points(x,index_sea[3,],col="violet")
	x=seq(13,(64*12+13),12)
	points(x,index_sea[4,],col="blue")
}

if (1==2){
	create_index("NAO","nao")
	create_index("AO","ao")
	create_index("PNA","pna")
	create_index("EAWR","eawr")
	create_index("PE","pe")
}

if (1==1){
	source("load.r")

	index=read.table(paste("../data/roh_data/","MEI",".txt",sep=""),sep="\t",fill=TRUE,skip=1)
	index=array(as.matrix(index[,2:13]),dim=c(66,12))
	print(index)
	index_sea=array(NA,dim=c(5,65))

	for (year in 1:65){
		print(index[year,])
		index_sea[1,year]=sum(c(index[year,3]*0.5,index[year,4:5],index[year,6]*0.5))/3
		index_sea[2,year]=sum(c(index[year,6]*0.5,index[year,7:8],index[year,9]*0.5))/3
		index_sea[3,year]=sum(c(index[year,9]*0.5,index[year,10:11],index[year,12]*0.5))/3
		index_sea[4,year]=sum(c(index[year,12]*0.5,index[(year+1),1:2],index[(year+1),3]*0.5))/3
		index_sea[5,year]=sum(c(index[year,1]*0.5,index[year,2:11],index[year,12]*0.5))/12
	}

	print(index_sea)

	years <- dim.def.ncdf("years",units="1",vals=1:65,unlim=FALSE)
	seasons <- dim.def.ncdf("seasons",units="spri,sum,aut,win,yea",vals=1:5, unlim=FALSE)

	mei <- var.def.ncdf(name="mei",units="bla",longname=paste("MEI"),dim=list(seasons,years), missval=-9999.0)

	vars=list(mei)
	   
	nc = create.ncdf(paste("../data/roh_data/","MEI","_sea.nc",sep=""),vars)

	put.var.ncdf(nc,"mei",index_sea)  


	close.ncdf(nc) 
}

