

recurrence_plot<- function(q=269){
	library(nonlinearTseries)
   	pdf(paste("../plots/",dataset,additional_style,"/",trendID,"/gridded/",trendID,"_",q,"_","reccurence.pdf",sep=""))
   	recurrencePlot(time.series=dat$tas[q,151:252,1:30],embedding.dim=1,time.lag=0,radius=1)
   	recurrencePlot(time.series=dat$tas[q,151:252,31:60],embedding.dim=1,time.lag=0,radius=1)
   	tmp1<<-rqa(time.series=dat$tas[q,151:252,1:30],embedding.dim=1,time.lag=0,radius=1)
   	tmp2<<-rqa(time.series=dat$tas[q,151:252,31:60],embedding.dim=1,time.lag=0,radius=1)

   	graphics.off()

}