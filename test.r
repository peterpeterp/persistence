
x=seq(1,10,1)
y=c(100,50,60,30,30,10,0,20,10,5)

xy=data.frame(y=y,x=x)
fit=nls(y~exp(a+b*x),data=xy,start=list(a=0,b=0))
print(summary(fit))

pdf(file="../plots/test.pdf")
plot(x,y)
print(exp(summary(fit)$parameters[1]))
yfit=exp(summary(fit)$parameters[1])*exp(x*summary(fit)$parameters[2])
print(yfit)
lines(x,yfit)