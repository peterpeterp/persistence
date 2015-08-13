
source("load.r")

index=read.table("../data/roh_data/NAO.txt")

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
write.table(index_sea,"../data/roh_data/NAO_sea.txt")

years <- dim.def.ncdf("years",units="1",vals=1:65,unlim=FALSE)
seasons <- dim.def.ncdf("seasons",units="spri,sum,aut,win,yea",vals=1:5, unlim=FALSE)

nao <- var.def.ncdf(name="nao",units="bla",longname=paste("NAO"),dim=list(seasons,years), missval=-9999.0)

vars=list(nao)
   
nc = create.ncdf("../data/roh_data/NAO_sea.nc",vars)

put.var.ncdf(nc,"nao",index_sea)  


close.ncdf(nc) 


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
	ntot=1319
	trendID="91_5"
	states=2
	transNumb=states*states
	transition_names=c("cc wc cw ww")
    nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/markov/",trendID,"_markov_",states,"states.nc",sep=""))
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")
    seasons=c("spring","summer","autumn","winter","year")

    correlation=array(NA,dim=c(ntot,5,transNumb))
    for (sea in 1:length(seasons)){
        markov=get.var.ncdf(nc,"markov")[1:length(dat$ID),sea,,]
        for (q in 1:ntot){
        	for (trans in 1:transNumb){
        		print(dim(markov))
				if (length(which(!is.na(markov[q,trans,])))>10){
					correlation[q,sea,trans]=cor(x=as.vector(index_sea[sea,]),y=as.vector(markov[q,trans,]),use="complete")
					print(correlation[q,sea,trans])
				}
			}
        }
    }
    print(correlation)

    close.ncdf(nc)

    ID <- dim.def.ncdf("ID",units="ID",vals=1:ntot, unlim=FALSE)
    transitions <- dim.def.ncdf("transitions",units=transition_names,vals=1:transNumb,unlim=FALSE)
    seasons <- dim.def.ncdf("seasons",units="spri,sum,aut,win,yea",vals=1:5, unlim=FALSE)

    cor <- var.def.ncdf(name="cor_nao",units="bla",longname=paste("NAO markov cor",transition_names),dim=list(ID,seasons,transitions), missval=-9999.0)


    vars=list(cor)
   
    nc = create.ncdf(paste("../data/",trendID,"/",states,"_states/oszi/NAO.nc",sep=""),vars)

    put.var.ncdf(nc,"cor_nao",correlation)  


    close.ncdf(nc) 

}

if (1==1){
	trendID="91_5"
	states=2
	transition_names=c("cc","wc","cw","ww")
    seasons=c("spring","summer","autumn","winter","year")

    library(SDMTools)
	source("functions_regional.r")
	source("map_plot.r")
	library(rworldmap)
	library(fields)
	worldmap = getMap(resolution = "low")
	ntot=1319
    dat=dat_load("../data/HadGHCND_TX_data3D.day1-365.1950-2014.nc")

    nc=open.ncdf(paste("../data/",trendID,"/",states,"_states/oszi/NAO.nc",sep=""))
    correlation=get.var.ncdf(nc,"cor_nao")
	reihen=array(NA,dim=c(20,ntot))
	titel=c()

	for (sea in 1:5){
		for (trans in 1:4){
			reihen[((sea-1)*4+trans),]=correlation[,sea,trans]
			titel[((sea-1)*4+trans)]=paste("correlation between NAO and transition probability",transition_names[trans],"in",seasons[sea])
		}
	}


    map_allgemein(dat,filename_plot="../plots/cor_nao.pdf",worldmap=worldmap,reihen=reihen,titel=titel,farbe_mitte="0")
}


