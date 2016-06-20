

multi <- function(x1,x2,x3,a,b,c){
	return(x1*(abs(a)+abs(b)+abs(c))+x2*(abs(b)+abs(c))+x3*abs(c))
}

y=c(8.13,9.11,10.88,8.87,9.51,11.73)    
x1=c(12.53,12.17,12.82,11.98,12.17,10.46)
x2=c(45.62,45.56,52.38,49.09,52.02 ,56.24)
x3=c(90.59,86.95,97.97,94.07,94.54,99.68)
xy<-data.frame(y=y,x1=x1,x2=x2,x3=x3)

data=read.table("/home/peter/Downloads/demandRegion.txt")
xy<-data.frame(y=data[,5],x1=data[,2],x2=data[,3],x3=data[,4])

tmp1<-nls(y~multi(x1,x2,x3,a,b,c),data=xy,algorithm="port",start=c(a=0.01,b=0.01,c=0.01),lower=c(0.0,0,0),upper=c(0.1,0.1,0.1))
param1<-summary(tmp1)$parameters
tmp2<-nls(y~multi(x1,x2,x3,a,b,c),data=xy,algorithm="port",start=c(a=param1[1],b=param[2],c=param[3]),lower=c(0.0,0,0),upper=c(Inf,Inf,Inf))
param2<-summary(tmp1)$parameters
tmp3<-nls(y~multi(x1,x2,x3,a,b,c),data=xy,start=c(a=0.01,b=0.01,c=0.01),nls.control(warnOnly=TRUE))
param3<-summary(tmp3)$parameters
print(sum(abs(param3[1:3])))
print(sum(abs(param3[2:3])))
print(abs(param3[3]))
