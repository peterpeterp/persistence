# ------------------------------------------------------------------------------------------------
plot_boxplot <- function(quans,x,width,color="white",border="black",density=NA){
    links=x-width/2
    rechts=x+width/2
    polygon(x=c(rechts,links,links,rechts),y=c(quans[2],quans[2],quans[4],quans[4]),col=color,border=border,density=density)    
    #polygon(x=c(rechts,links,links,rechts),y=c(quans[2],quans[2],quans[4],quans[4]),col=color,border=border,density=density)    
    for (qu in quans[1:5]){
        lines(c(links,rechts),c(qu,qu),col=border)
    }
    lines(c(links,links),quans[c(2,4)],col=border)
    lines(c(rechts,rechts),quans[c(2,4)],col=border)
    lines(c(x,x),quans[c(1,2)],lty=2,col=border)
    lines(c(x,x),quans[c(4,5)],lty=2,col=border)
}


plot_regional_boxplots <- function(period,region_name,trendID,additional_style,dataset,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){

    nc_oth = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_others.nc",sep=""))
    others=var.get.nc(nc_oth,"other_stuff")
    nc_qua = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",region_name,"/",period,"/",trendID,"_",dataset,"_",region_name,"_",period,"_quantiles.nc",sep=""))
    quantiles=var.get.nc(nc_qua,"quantile_stuff")

    regNumb=dim(quantiles)[2]
    
    color=c(rgb(0.5,0.5,1,0.8),rgb(1,0.5,0.5,0.8))
    pos=c(-1,1)*0.15
    width=6
    height=3

    # regional focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",region_name,"/",period,"/",region_name,"_",trendID,"_",period,"_boxplots_regional_new.pdf",sep=""),width=width,height=height)
    par(mfrow=c(1,1))
    par(mar=c(1,4,2,3))   
    at_=seq(1, regNumb, 1)
    at_=c(at_-0.15,at_+0.15)
    buchstaben=c("a","b","c","d")
    drueber=2
    for (sea in 1:length(season_names)){   
        season=season_names[sea]

        plot(NA,xlim=c(0,7.5),ylim=c(0,25),frame.plot=FALSE,axes=FALSE,ylab="days")#(max(quantiles[sea,,,1:5])+drueber)
        axis(2)
        for (i in axis(2)){
            abline(h=i,col=rgb(0.5,0.5,0.5,0.5),lty=3)
        }
        legend("topright",legend=c(buchstaben[sea]),bty="n")
        #text(8,23,label=buchstaben[sea],cex=2)#(max(quantiles[sea,,,1:5],na.rm=TRUE)+drueber-0.5)

        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,,1],reg+pos[state],0.3,color[state])
                text(reg,24,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))#max(quantiles[sea,,,1:5])+drueber
            }
        }        
        for (oth in c(1)){
            points(at_[1:regNumb],others[sea,,1,oth],col="blue",pch=oth)
            points(at_[(regNumb+1):(regNumb*2)],others[sea,,2,oth],col="red",pch=oth)
        }
    }
    graphics.off()
}


plot_regional_boxplots_vergleich <- function(period1,period2,region_name,trendID,additional_style,dataset,region_names=c("wNA","cNA","eNA","Eu","wA","cA","eA")){
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",region_name,"/",period1,"/",trendID,"_",dataset,"_",region_name,"_",period1,"_quantiles.nc",sep=""))
    quantiles1=var.get.nc(nc,"quantile_stuff")
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",region_name,"/",period2,"/",trendID,"_",dataset,"_",region_name,"_",period2,"_quantiles.nc",sep=""))
    quantiles2=var.get.nc(nc,"quantile_stuff")
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",region_name,"/",period1,"/",trendID,"_",dataset,"_",region_name,"_",period1,"_others.nc",sep=""))
    others1=var.get.nc(nc,"other_stuff")
    nc = open.nc(paste("../data/",trendID,"/",dataset,additional_style,"/regional/",region_name,"/",period2,"/",trendID,"_",dataset,"_",region_name,"_",period2,"_others.nc",sep=""))
    others2=var.get.nc(nc,"other_stuff")

    regions=var.get.nc(nc,"region")
    regNumb=dim(quantiles1)[2]
    seaNumb=5
    
    season_names=c("MAM","JJA","SON","DJF","4seasons")

    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",region_name,"/",region_name,"_",period1,"_diff_",period2,"_boxplots_regional.pdf",sep=""))
    par(mfrow=c(1,1))
    par(mar=c(1,5,4,3))   
    at_=seq(1, regNumb, 1)
    at_=c(at_-0.15,at_+0.15)
    color=c()
    maxi=c()
    taus=c(0.05,0.25,0.5,0.75,0.95,0.91,0.98)
    color=c(rgb(0.5,0.5,1,0.8),rgb(1,0.5,0.5,0.8))
    pos=c(-1,1)*0.15

    quantiles=quantiles2-quantiles1
    others=others2-others1

    # regional focus
    for (sea in 1:length(season_names)){   
        season=season_names[sea]

        plot(NA,xlim=c(0,8),ylim=c(min(quantiles[sea,,,1:5,1]),(max(quantiles[sea,,,1:5,1])+0.5)),frame.plot=FALSE,axes=FALSE,main=season,ylab="days")
        axis(2)
        for (reg in 1:regNumb){
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,,1],reg+pos[state],0.3,color[state])
                text(reg,max(quantiles[sea,,,1:5,1])+0.5,region_names[reg],col=rgb(0.5,0.5,0.5,0.5))
            }
        }
        for (oth in c(1)){
            points(at_[1:regNumb],others[sea,,1,oth],col="blue",pch=oth)
            points(at_[(regNumb+1):(regNumb*2)],others[sea,,2,oth],col="red",pch=oth)
        }
        for (quA in c(4,5)){
            points(at_[1:regNumb],quantiles[sea,,1,quA,1],col="black",pch=2)
            points(at_[(regNumb+1):(regNumb*2)],quantiles[sea,,2,quA,1],col="black",pch=2)
        }
    }
    graphics.off()

    # seasonal focus
    pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",region_name,"/",region_name,"_",period1,"_diff_",period2,"_boxplots_seasonal.pdf",sep=""))
    par(mfrow=c(1,1))
    par(mar=c(1,5,4,3))   
    at_=seq(1, seaNumb, 1)
    at_=c(at_-0.15,at_+0.15)

    for (reg in 1:regNumb){
        plot(NA,xlim=c(0,8),ylim=c(min(quantiles[,reg,,1:5,1]),(max(quantiles[,reg,,1:5,1])+0.5)),frame.plot=FALSE,axes=FALSE,main=region_names[reg],ylab="days")
        axis(2)
        for (sea in 1:seaNumb){   
            season=season_names[sea]     
            for (state in 1:2){
                plot_boxplot(quantiles[sea,reg,state,,1],sea+pos[state],0.3,color[state])
                text(sea,max(quantiles[,reg,,1:5,1])+0.5,season_names[sea],col=rgb(0.5,0.5,0.5,0.5))
            }
        }
        for (oth in c(1)){
            points(at_[1:seaNumb],others[,reg,1,oth],col="blue",pch=1)
            points(at_[(seaNumb+1):(seaNumb*2)],others[,reg,2,oth],col="red",pch=1)
        }
        for (quA in c(4,5)){
            points(at_[1:seaNumb],quantiles[,reg,1,quA,1],col="black",pch=2)
            points(at_[(seaNumb+1):(seaNumb*2)],quantiles[,reg,2,quA,1],col="black",pch=2)
        }
    }
    graphics.off()
}