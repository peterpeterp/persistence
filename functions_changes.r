
linear_trend <- function(y,t){
    return(c(lm(y~t)$coef[2],rq(dither(y,value=noise_level)~t,0.95)$coef[2]))
}

ks_statistics <- function(y){
    L_half<-round(length(y)/2)
    ks_full<-ks.test(y[1:L_half],y[(L_half+1):length(y)])$statistic

    distr_cut<-median(y,na.rm=TRUE)
    y[y<=distr_cut]=NA
    ks_upper<-ks.test(y[1:L_half],y[(L_half+1):length(y)])$statistic
    return(c(ks_full,ks_upper))
}


bootstrap_combining_regions <- function(Matrix,gridPoints,l,R){

    # create nice matrix that can be handled for shuffling
	dimX2<-dim(Matrix)[3]
    regNumb<-length(gridPoints)

    # add NAs so that the day length is a multiple of l
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

    # create an index matrix that can be shuffled
	original_order<-array(1:data_length,c(l,data_length/l))

    # create a time vector for regressions
	t<-array(NA,dim(X))
	for (i in 1:dim(X)[1]){t[i,]=1:dim(X)[2]}
	t<-as.vector(t)

	statistics<-array(NA,c(R,(regNumb+1),4))

    # calculate original statistics
    pos<-0
    for (q in 1:regNumb){
        y<-as.vector(X[(pos+1):(pos+gridPoints[q]),])

        tReg<-array(NA,dim(X[(pos+1):(pos+gridPoints[q]),]))
        for (i in 1:dim(tReg)[1]){tReg[i,]=1:dim(tReg)[2]}
        tReg<-as.vector(tReg)

        nona<-which(!is.na(y))
        statistics[1,q,1:2]=linear_trend(y=y[nona],t=tReg[nona])
        statistics[1,q,3:4]=ks_statistics(y=y)
        pos<-pos+gridPoints[q]
    }

    y<-as.vector(X)
    nona<-which(!is.na(y))
    statistics[1,(regNumb+1),1:2]=linear_trend(y=y[nona],t=t[nona])
    statistics[1,(regNumb+1),3:4]=ks_statistics(y=y)

	print("original slope")
    print(statistics[1,,])
    print("original slope [day/year]")
    print(statistics[1,,]*365)

    # start shuffling
    percentage<-0
    cat(paste("\n0 -> -> -> -> -> 100\n"))
    for (r in 2:R){
        if (r/R*100 > percentage){
            cat("-")
            percentage<-percentage+5
        }
        shuffled_X<-X*NA
        pos<-0
       	for (q in 1:regNumb){
            # shuffle one region and evaluate trends for that shuffled region
            singleRegion<-X[(pos+1):(pos+gridPoints[q]),as.vector(original_order[,sample(1:(data_length/l),data_length/l,replace=FALSE)])]
            y<-as.vector(singleRegion)
            tReg<-array(NA,dim(X[(pos+1):(pos+gridPoints[q]),]))
            for (i in 1:dim(tReg)[1]){tReg[i,]=1:dim(tReg)[2]}
            nona<-which(!is.na(y))
            statistics[r,q,1:2]=linear_trend(y=y[nona],t=tReg[nona])
            statistics[r,q,3:4]=ks_statistics(y=y)

            # put shuffled region block in "shuffled overregion" matrix
			shuffled_X[(pos+1):(pos+gridPoints[q]),1:data_length]=singleRegion
			pos<-pos+gridPoints[q]
		}

        #evaluate for overRegion
        y<-as.vector(shuffled_X)
        nona<-which(!is.na(y))
        statistics[r,(regNumb+1),1:2]=linear_trend(y=y[nona],t=tReg[nona])
        statistics[r,(regNumb+1),3:4]=ks_statistics(y=y)
	}
	return(statistics)
}

trend_bootstrap_for_large_region <- function(add_name=1){
    filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_reg_daily_binned_dur_",season,".nc",sep="")
    durMat<<-var.get.nc(open.nc(filename),"daily_binned_periods")[,,,((yearPeriod[1]-1950)*365+1):((yearPeriod[2]-1950)*365)]
    print(dim(durMat))

    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    gridPoints<-array(NA,length(regions))
    for (i in 1:length(regions)){gridPoints[i]=length(which(attribution==regions[i]))}

    Matrix<-array(NA,c(length(regions),120,dim(durMat)[4]))
    for (i in 1:length(regions)){Matrix[i,,]=durMat[state,regions[i],,]}

    print(paste(season,state_names[state],region_name,over_region_name,"(contains several regions)"))
    print(regions)

	statistics<-bootstrap_combining_regions(Matrix,gridPoints,7,replics)

	filename <- paste("../data/",dataset,additional_style,"/",trendID,folder,period,"/bootstrap/",trendID,dataset,"_",period,"_boot_",region_name,"_",over_region_name,"_",season,"_",state_names[state],"_",add_name,".nc",sep="") ; print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "included regions:", "NC_CHAR","see regions")
            
    dim.def.nc(nc_out,"replics",dimlength=replics, unlim=FALSE)
    dim.def.nc(nc_out,"regs",dimlength=(length(regions)+1), unlim=FALSE)
    dim.def.nc(nc_out,"outs",dimlength=4,unlim=FALSE)

    var.def.nc(nc_out,"regions","NC_SHORT",c(1))

    var.def.nc(nc_out,"statistics","NC_DOUBLE",c(0,1,2))
    att.put.nc(nc_out, "statistics", "missing_value", "NC_DOUBLE", -99999)
    att.put.nc(nc_out, "statistics", "dim_explanation", "NC_CHAR", "(Regional values, big Region value) x (1-original, 1000-random) x (lr, qr 0.95, ks_full, ks_upper")

    var.put.nc(nc_out,"statistics",statistics)   
    var.put.nc(nc_out,"regions",c(regions,999))           

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
yearPeriod<-c(1979,2011)
period<-"1979-2011"
replics<-1001


#regions<-c(9,15,24) ; over_region_name<-"SHml"

#regions<-c(5,11,14,18,21,22) ; over_region_name<-"NHst"

#regions<-c(1,2,6,10,13,19,23) ; over_region_name<-"NHpo"

#regions<-c(8,17) ; over_region_name<-"tro"

regions<-c(3,4,7,12,16,20) ; over_region_name<-"NHml"


if (FALSE){
    season<-"DJF"
    state<-1
    trend_bootstrap_for_large_region(add_name=1)
}

id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(id)

season<-"4seasons"
state<-2
trend_bootstrap_for_large_region(add_name=id)
adsas

for (season in season_names){
    for (state in 1:2){
        trend_bootstrap_for_large_region(add_name=id)
    }
}



