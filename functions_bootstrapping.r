library(boot)

trend_bootstrap <- function(season){
	filename <- paste("../data/",dataset,additional_style,"/",trendID,"/regional/",region_name,"/",trendID,dataset,"_",region_name,"_reg_daily_binned_duration_",season,".nc",sep="")
	durMat<<-var.get.nc(open.nc(filename),"daily_binned_periods")

}




#trend_bootstrap("JJA")



a=array(NA,c(10,20))
for (i in 1:10){a[i,]=sample(1:100,20,replace=FALSE)}
#for (i in 1:10){a[i,]=1:20}

linear_trend <- function(y,t){
	return(c(rq(dither(y,value=noise_level)~t,taus)$coef[2,],NA,lm(y~t)$coef[2]))
}

bootstrap <- function(x,l,R){
	dimX<-dim(x)
	if (round(dimX[2]/l)!=dimX[2]/l){
		needed<-abs(round(dimX[2]/l)-dimX[2]/l)*l
		X<-array(NA,c(dimX[1],dimX[2]+needed))
		X[1:dimX[1],1:dimX[2]]=x
	}
	data_length<-dim(X)[2]
	original_order<-array(1:data_length,c(l,data_length/l))

	t<-array(NA,dim(X))
	for (i in 1:dim(X)[1]){t[i,]=1:dim(X)[2]}
	t<-as.vector(t)

	statistics<-array(NA,c(R,5))
	for (r in 1:R){
		statistics[r,1:5]=linear_trend(y=as.vector(X[,as.vector(original_order[,sample(1:(data_length/l),data_length/l,replace=FALSE)])]),t=t)
	}
	print(statistics)
}


noise_level<-0.000001
bootstrap(a,7,50)