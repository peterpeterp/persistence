
linear_trend <- function(y,t){
	return(c(rq(dither(y,value=noise_level)~t,taus)$coef[2,],NA,lm(y~t)$coef[2]))
}

bootstrap <- function(x,l,R){
	dimX<-dim(x)
	needed<-0
	for ( i in 1:l){
		if (round((dimX[2]+needed)/l)!=(dimX[2]+needed)/l){needed<-needed+1}
	}
	X<-array(NA,c(dimX[1],dimX[2]+needed))
	X[1:dimX[1],1:dimX[2]]=x

	data_length<-dim(X)[2]
	original_order<-array(1:data_length,c(l,data_length/l))

	t<-array(NA,dim(X))
	for (i in 1:dim(X)[1]){t[i,]=1:dim(X)[2]}
	t<-as.vector(t)

	original_order<<-original_order

	statistics<-array(NA,c(R,5))
	statistics[1,1:5]=linear_trend(y=as.vector(X),t=t)

	print("original slope")
	print(statistics[1,1:5])
    percentage<-0
    cat(paste("\n0 -> -> -> -> -> 100\n"))
    for (r in 2:R){
        if (r/R*100 > percentage){
            cat("-")
            percentage<-percentage+5
        }
		statistics[r,1:5]=linear_trend(y=as.vector(X[,as.vector(original_order[,sample(1:(data_length/l),data_length/l,replace=FALSE)])]),t=t)
	}
	return(statistics)
}

trend_bootstrap <- function(state=1,region=1,replics=9999){
	filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_reg_daily_binned_duration_",season,".nc",sep="")
	durMat<-var.get.nc(open.nc(filename),"daily_binned_periods")
    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    gridPoints<-length(which(attribution==region))

    print(paste(season,state_names[state],region_name,"region",region,"grid points in region",gridPoints))
	statistics<-bootstrap(durMat[state,region,1:gridPoints,],7,replics)

	filename <- paste("../data/",dataset,additional_style,"/",trendID,folder,"/",period,"/bootstrap/",trendID,dataset,"_",period,"linear_trend_boot_",season,"_",state_names[state],"_",region_name,"-",region,".nc",sep="") ; print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", paste(region_name,region))
            
    dim.def.nc(nc_out,"replics",dimlength=10000, unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=5,unlim=FALSE)

    var.def.nc(nc_out,"trends","NC_DOUBLE",c(0,1))
    att.put.nc(nc_out, "trends", "missing_value", "NC_DOUBLE", -99999)
    att.put.nc(nc_out, "trends", "dim_explanation", "NC_CHAR", "(1-original, 9999-random) - (qr 0.75,0.95,0.99, NA, mean")

    var.put.nc(nc_out,"trends",statistics)              

    close.nc(nc_out) 
}


bootstrap_combining_regions <- function(Matrix,gridPoints,l,R){


	dimX2<-dim(Matrix)[3]
	needed<-0
	for ( i in 1:l){
		if (round((dimX2+needed)/l)!=(dimX2+needed)/l){needed<-needed+1}
	}
	
	X<-array(NA,c(sum(gridPoints),dimX2+needed))
	pos<-0
	for (i in 1:length(gridPoints)){
		X[(pos+1):(pos+gridPoints[i]),1:dimX2]=Matrix[i,1:gridPoints[i],]
		pos<-pos+gridPoints[i]
	}

	data_length<-dim(X)[2]
	original_order<-array(1:data_length,c(l,data_length/l))

	t<-array(NA,dim(X))
	for (i in 1:dim(X)[1]){t[i,]=1:dim(X)[2]}
	t<-as.vector(t)

	statistics<-array(NA,c(R,5))
	statistics[1,1:5]=linear_trend(y=as.vector(X),t=t)

	print("original slope")
	print(statistics[1,1:5])
    percentage<-0
    cat(paste("\n0 -> -> -> -> -> 100\n"))
    for (r in 2:R){
        if (r/R*100 > percentage){
            cat("-")
            percentage<-percentage+5
        }
        shuffled_X<-X*NA
        pos<-0
       	for (i in 1:length(gridPoints)){
       		print(length(original_order))
       		print(length(as.vector(original_order[,sample(1:(data_length/l),data_length/l,replace=FALSE)])))
       		print(length(1:data_length))
			shuffled_X[(pos+1):(pos+gridPoints[i]),1:data_length]=X[(pos+1):(pos+gridPoints[i]),as.vector(original_order[,sample(1:(data_length/l),data_length/l,replace=FALSE)])]
			pos<-pos+gridPoints[i]
		}
        shuffled_X<<-shuffled_X

        adas
		statistics[r,1:5]=linear_trend(y=as.vector(X[,as.vector(original_order[,sample(1:(data_length/l),data_length/l,replace=FALSE)])]),t=t)
	}
	return(statistics)
}

trend_bootstrap_for_large_region <- function(state=1,regions=c(3,4,7,12,16,20),replics=9999){
	filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_reg_daily_binned_duration_",season,".nc",sep="")
	durMat<-var.get.nc(open.nc(filename),"daily_binned_periods")

    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    for (i in 1:length(regions)){gridPoints[i]=length(which(attribution==regions[i]))}

    Matrix<-array(NA,c(length(regions),120,dim(durMat)[4]))
    for (i in 1:length(regions)){Matrix[i,,]=durMat[state,regions[i],,]}

	statistics<-bootstrap(Matrix,gridPoints,7,replics)







	filename <- paste("../data/",dataset,additional_style,"/",trendID,folder,"/",period,"/bootstrap/",trendID,dataset,"_",period,"linear_trend_boot_",season,"_",state_names[state],"_",region_name,"-",region,".nc",sep="") ; print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", paste(region_name,region))
            
    dim.def.nc(nc_out,"replics",dimlength=10000, unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=5,unlim=FALSE)

    var.def.nc(nc_out,"trends","NC_DOUBLE",c(0,1))
    att.put.nc(nc_out, "trends", "missing_value", "NC_DOUBLE", -99999)
    att.put.nc(nc_out, "trends", "dim_explanation", "NC_CHAR", "(1-original, 9999-random) - (qr 0.75,0.95,0.99, NA, mean")

    var.put.nc(nc_out,"trends",statistics)              

    close.nc(nc_out) 
}



init <- function(){
	library(boot)
    library(quantreg)
    library(RNetCDF)
    nday<<-91
    nyr<<-7
    trendID<<-paste(nday,"_",nyr,sep="")
    dataset<<-"_TMean"
    additional_style<<-""
    taus<<-c(0.75,0.95,0.99)
    season_names<<-c("MAM","JJA","SON","DJF","4seasons")
    state_names<<-c("cold","warm")

    taus<<-c(0.75,0.95,0.99)

    #noise level to avoid crash of rq()
    noise_level<<-0.00001
}

init()

region_name<-"ward24"
folder<-paste("/regional/",region_name,"/",sep="")
period<-"1950-2014"
season<-"JJA"


a=array(NA,c(4,10,20))
for (i in 1:4){
	for (j in 1:10){
		a[i,j,]=(1:20)*10^(i-3)
	}
}

#print(a)


bootstrap_combining_regions(a,c(3,4,2,3),7,2)




adsasd

id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(id)

regions_to_analyse<-c(3,4,7,12,16,20)
regNumb<-length(regions_to_analyse)
regions_to_analyse<-c(regions_to_analyse,regions_to_analyse)

if (id<=regNumb){trend_bootstrap(state=1,region=regions_to_analyse[id])}
if (id>regNumb){trend_bootstrap(state=2,region=regions_to_analyse[id])}
