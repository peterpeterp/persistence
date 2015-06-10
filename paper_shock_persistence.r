#load 
library(RNetCDF)
library(stats)

laenge<-250
t<-seq(1,laenge,1)
x<-seq(1,laenge,1)
y<-seq(1,laenge,1)
ut<-seq(1,laenge,1)

for (i in 0:laenge){
	ut[i]=runif(1,-1,1)
}


for (i in 0:laenge){
	if(i<3){
		y[i]=ut[i]
		x[i]=ut[i]
	}
	if(i>2){
		x[i]=1.7*x[i-1]-0.7*x[i-2]+ut[i]
		y[i]=0.3*y[i-1]+0.7*y[i-2]+ut[i]
	}
	
}


Px=1
Py=1
for (i in 1:1000){
	Px=Px+0.7^i
	Py=Py+(-0.7)^i
}

print(Px)
print(Py)

bayesian_i_c <- function(x,order){
	arma_x=arima(x,order=order,method="ML")
	bic=-2*arma_x$loglik+(sum(order)-1)*log(length(x))
	return(bic)
}

shock_ar_1 <- function(x){
	arma_x=arima(x,order=c(1,1,0),method="ML")
	ar_x=arma_x$coef[1]
	Px=1
	psi=array(0,100)
	for (i in 1:100){
		Px=Px+ar_x^i
	}
	return(Px)
}

shock_ar_2 <- function(x){
	arma_x=arima(x,order=c(2,1,0),method="ML")
	ar_x=arma_x$coef[1:2]
	a=ar_x[1]
	b=ar_x[2]
	y1=(-b+sqrt(b^2+4*a))/(2*a)
	y2=(-b-sqrt(b^2+4*a))/(2*a)
	la1=1/y1
	la2=1/y2
	c1=la1/(la1-la2)
	c2=la2/(la2-la1)
	Px=1
	psi=array(0,100)
	for (i in 1:100){
		for (k in 0:i){
			psi[i]=c1*la1^k+c2*la2^k
		}
		Px=Px+psi[i]
	}
	return(Px)
}

print(bayesian_i_c(x,c(1,1,0)))
print(bayesian_i_c(x,c(2,1,0)))
print(bayesian_i_c(x,c(3,1,0)))
print(bayesian_i_c(x,c(0,1,1)))
print(bayesian_i_c(x,c(0,1,2)))
print(bayesian_i_c(x,c(1,1,1)))

print(shock_ar_1(x))
print(shock_ar_2(x))