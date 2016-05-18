
linear_trend <- function(y,t){
    #print(lm(y[1:length(y)]~t[1:length(t)],na.rm=TRUE))
    return(c(lm(y~t)$coef[2],rq(dither(y,value=noise_level)~t,0.95)$coef[2]))
    #return(c(lm(y[1:length(y)]~t[1:length(t)])$coef[2],rq(dither(y[1:length(y)],value=noise_level)~t[1:length(t)],0.95)$coef[2]))
}

ks_statistics <- function(y){

    L_half<-round(length(y)/2)
    ks_full<-ks.test(y[1:L_half],y[(L_half+1):length(y)])$statistic
    print("ks full ")
    print(tact-proc.time())
    tact<<-proc.time()

    distr_cut<-median(y,na.rm=TRUE)
    y[y<=distr_cut]=NA
    ks_upper<-ks.test(y[1:L_half],y[(L_half+1):length(y)])$statistic

    print("ks upper ")
    print(tact-proc.time())
    tact<<-proc.time()

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



	statistics<-array(NA,c(R,(regNumb+1),5))

    # calculate original statistics
    pos<-0
    for (q in 1:regNumb){
        print(q)

        tact<<-proc.time()
        tstart<<-proc.time()

        y<-as.vector(X[(pos+1):(pos+gridPoints[q]),])

        print("y ")
        print(tact-proc.time())
        tact<<-proc.time()

        tReg<-array(NA,dim(X[(pos+1):(pos+gridPoints[q]),]))
        for (i in 1:dim(tReg)[1]){tReg[i,]=1:dim(tReg)[2]}

        print("tReg ")
        print(tact-proc.time())
        tact<<-proc.time()

        nona<-which(!is.na(y))
        statistics[1,q,1:2]=linear_trend(y=y[nona],t=tReg[nona])

        print("linear trend ")
        print(tact-proc.time())
        tact<<-proc.time()

        statistics[1,q,3:4]=ks_statistics(y=y)
        #statistics[1,q,5]=MannKendall(y[which(!is.na(y))])$tau

        print("MK ")
        print(tact-proc.time())
        tact<<-proc.time()

        print(statistics[1,,])
    }
    print(statistics[1,,])

    y<-as.vector(X)
    nona<-which(!is.na(y))
    statistics[1,(regNumb+1),1:2]=linear_trend(y=y[nona],t=t[nona])
    statistics[1,(regNumb+1),3:4]=ks_statistics(y=y)
    #statistics[1,(regNumb+1),5]=MannKendall(y)$tau


	print("original slope")
	print(statistics[1,,])
    print(tstart-proc.time())

    adsas



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
            statistics[r,q,1:2]=linear_trend( y=y,t=t)
            statistics[r,q,3:4]=ks_statistics(y=y)
            statistics[r,q,5]=MannKendall(y=y)$tau

            # put shuffled region block in "shuffled overregion" matrix
			shuffled_X[(pos+1):(pos+gridPoints[q]),1:data_length]=singleRegion
			pos<-pos+gridPoints[q]
		}

        #evaluate for overRegion
        y<-as.vector(shuffled_X)
        statistics[1,q,1:2]=linear_trend( y=y,t=t)
        statistics[1,q,3:4]=ks_statistics(y=y)
        statistics[1,q,5]=MannKendall(y=y)$tau

        print(statistics[2,,])
        adasd
	}
	return(statistics)
}

trend_bootstrap_for_large_region <- function(state=1,regions=c(3,4,7,12,16,20),over_region_name="mid-lat-3-4-7-12-16-20",replics=100,add_name=1){
    filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_reg_daily_binned_dur_",season,".nc",sep="")
    #filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_reg_daily_binned_duration_",season,".nc",sep="")
    #print(dim(var.get.nc(open.nc(filename),"daily_binned_periods")))
    #durMatact<<-var.get.nc(open.nc(filename),"daily_binned_periods")[,,,((yearPeriod[1]-1950)*365+1):((yearPeriod[2]-1950)*365)]
    print(dim(durMat))


    attribution<-read.table(paste("../data/",dataset,"/ID_regions/",region_name,".txt",sep=""))[,1]
    gridPoints<-array(NA,length(regions))
    for (i in 1:length(regions)){gridPoints[i]=length(which(attribution==regions[i]))}

    Matrix<<-array(NA,c(length(regions),120,dim(durMat)[4]))
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
    datasetact<<-"_TMean"
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
season<-"JJA"


trend_bootstrap_for_large_region(state=2,add_name=id)

adsads

id<-as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(id)

trend_bootstrap_for_large_region(state=2,add_name=id)

