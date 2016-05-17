
linear_trend <- function(y,t){
	return(c(rq(dither(y,value=noise_level)~t,0.95)$coef[2],NA,lm(y~t)$coef[2]))
	#return(c(rq(dither(y,value=noise_level)~t,taus,na.action="na.omit")$coef[2,],NA,lm(y~t,na.action="na.omit")$coef[2]))
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

	statistics<-array(NA,c(R,3))
	statistics[1,1:3]=linear_trend(y=as.vector(X),t=t)

	print("original slope")
	print(statistics[1,1:3])
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
			shuffled_X[(pos+1):(pos+gridPoints[i]),1:data_length]=X[(pos+1):(pos+gridPoints[i]),as.vector(original_order[,sample(1:(data_length/l),data_length/l,replace=FALSE)])]
			pos<-pos+gridPoints[i]
		}

		statistics[r,1:3]=linear_trend(y=as.vector(shuffled_X),t=t)
	}
	return(statistics)
}

trend_bootstrap_for_large_region <- function(state=1,regions=c(3,4,7,12,16,20),over_region_name="mid-lat-3-4-7-12-16-20",replics=100,add_name=1){
	filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_reg_daily_binned_duration_",season,".nc",sep="")
	durMat<-var.get.nc(open.nc(filename),"daily_binned_periods")

    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    gridPoints<-array(NA,length(regions))
    for (i in 1:length(regions)){gridPoints[i]=length(which(attribution==regions[i]))}

    Matrix<-array(NA,c(length(regions),120,dim(durMat)[4]))
    for (i in 1:length(regions)){Matrix[i,,]=durMat[state,regions[i],,]}

    print(paste(season,state_names[state],region_name,over_region_name,"(contains several regions)"))
	statistics<-bootstrap_combining_regions(Matrix,gridPoints,7,replics)

	filename <- paste("../data/",dataset,additional_style,"/",trendID,folder,"/",period,"/bootstrap/",trendID,dataset,"_",period,"_linear_trend_boot_",season,"_",state_names[state],"_",region_name,"-",over_region_name,"_",add_name,".nc",sep="") ; print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR",over_region_name)
            
    dim.def.nc(nc_out,"replics",dimlength=100, unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=3,unlim=FALSE)

    var.def.nc(nc_out,"trends","NC_DOUBLE",c(0,1))
    att.put.nc(nc_out, "trends", "missing_value", "NC_DOUBLE", -99999)
    att.put.nc(nc_out, "trends", "dim_explanation", "NC_CHAR", "(1-original, 99-random) - (qr 0.95, NA, mean")

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


id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(id)

trend_bootstrap_for_large_region(state=2,add_name=id)

