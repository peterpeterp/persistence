library(stats4)


N<-100
x<-rnorm(N,1,2)

ll<-function(mu,sigma){-sum(log(dnorm(x,mu,sigma)))}

mle(ll,start=list(mu=1,sigma=1))