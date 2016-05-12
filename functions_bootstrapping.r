library(boot)




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

trend_bootstrap <- function(state=1,region=1,region_name="ward24",season="JJA",period="1950-2014",replics=9999){
	filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_reg_daily_binned_duration_",season,".nc",sep="")
	#durMat<-var.get.nc(open.nc(filename),"daily_binned_periods")
    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    gridPoints<-length(which(attribution==region))

    print(paste(season,state_names[state],region_name,"region",region,"grid points in region",gridPoints))
	statistics<-bootstrap(durMat[state,region,1:gridPoints,],7,replics)

	filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",period,"/bootstrap/",trendID,dataset,"_",period,"linear_trend_boot_",season,"_",state_names[state],"_",region_name,"-",region,".nc",sep="") ; print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", paste(region_name,region))
            
    dim.def.nc(nc_out,"replics",dimlength=10000, unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=5,unlim=FALSE)

    var.def.nc(nc_out,"trends","NC_DOUBLE",c(0,1))
    att.put.nc(nc_out, "trends", "missing_value", "NC_DOUBLE", 99999)
    att.put.nc(nc_out, "trends", "dim_explanation", "NC_CHAR", "(1-original, 9999-random) - (qr 0.75,0.95,0.99, NA, mean")

    var.put.nc(nc_out,"trends",statistics)              

    close.nc(nc_out) 

}

noise_level<-0.000001
trend_bootstrap(state=2,region=12,region_name="ward24",season="JJA",replics=3)


if (FALSE){
	a<-array(NA,c(10,20))
	for (i in 1:10){a[i,]=1:20}
	print(bootstrap(a,7,5))
}



#file test
if(FALSE){
	filename<-paste("../data/",dataset,additional_style,"/",trendID,"/gridded/",trendID,dataset,"_duration_","JJA",".nc",sep="")
	nc_dur<-open.nc(filename)
	dur<-var.get.nc(nc_dur,"dur")
	dur_mid<-var.get.nc(nc_dur,"dur_mid") 

	print(dur[458,2,1:10])
	print(dur_mid[458,2,1:10])
	print(durMat[2,12,1,276:360])
	print(durMat[1,12,1,276:360])
}

