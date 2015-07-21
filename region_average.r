
points_to_regions <- function(dat){
    ntot=length(dat$ID)
    region = array(NA,dim=c(ntot,3))
    points=cbind(x=dat$lon,y=dat$lat)

    poli=read.table("../data/srex.txt")
    for (k in 1:dim(poli)[1]){
        x=c()
        y=c()
        for (i in 1:6){
            if (!is.na(poli[k,i])){
                x[i]=poli[k,i]
                y[i]=poli[k,(i+6)]
            }
        }
        poligon=cbind(x=x,y=y)
        print(poligon)
        inside=pnt.in.poly(points,poligon)$pip
        region[which(inside==1),1]=poli[k,13]
    }

    poli=read.table("../data/7rect.txt")
    for (k in 1:dim(poli)[1]){
        x=c()
        y=c()
        for (i in 1:6){
            if (!is.na(poli[k,i])){
                x[i]=poli[k,i]
                y[i]=poli[k,(i+6)]
            }
        }
        poligon=cbind(x=x,y=y)
        print(poligon)
        inside=pnt.in.poly(points,poligon)$pip
        region[which(inside==1),2]=poli[k,13]
    }

    region[1:ntot,3]=dat$ID
    write.table(region,"../data/ID-regions.txt")
    return(0)
}



regions_color <- function(values,sigs,worldmap,title,poli){
    jet.colors <- colorRampPalette( c( "violet","blue","white","yellow","red") )
    nbcol <- 101
    color <- jet.colors(nbcol)

    y=c()
    index=c()
    signi=c()
    j=0
    for (i in 1:dim(poli)[1]){
        if (!is.na(values[i])){
            j=j+1
            y[j]=values[i]
            index[j]=i 
            signi[j]=sprintf("%.02f",sigs[i])         
        }
    }
    aushol=max(c(abs(max(y)),abs(min(y))))
    y[j+1]=-aushol
    y[j+2]=aushol
    facetcol <- cut(y,nbcol)  

    plot(worldmap,main=title)

    for (i in 1:j){
        lon=poli[index[i],1:6]
        lat=poli[index[i],7:12]
        lon=lon[!is.na(lon)]
        lat=lat[!is.na(lat)]
        polygon(x=lon,y=lat,col=color[facetcol[i]],border="green")
        text(mean(lon),mean(lat),label=signi[i],cex=0.7,col="black")
    }
    image.plot(legend.only=T, zlim=range(y), col=color)
    return()
}

average_regional_trend <- function(dat,value,regions,regNumb){
    # value is an array of dim=c(1319,??)
    size=dim(value)[2]
    valOut=array(NA,dim=c(regNumb,size))
    for (reg in 1:regNumb){
        inside=which(regions==reg)
        for (j in 1:size){
            valOut[reg,j]=mean(value[inside,j],na.rm=TRUE)
        }
    }
    out=array(NA,c(2,regNumb))
    out_sig=array(NA,c(2,regNumb))

    x=seq(1,65,1)
    for (reg in 1:regNumb){
        y=valOut[reg,]
        if (length(which(is.nan(y)==1))<10){
            tmp=MannKendall(y)
            out[1,reg]=tmp[1]$tau
            out_sig[1,reg]=tmp[2]$sl
            lm.r=lm(y~x)
            out[2,reg]=summary(lm.r)$coefficients[2]
            out_sig[2,reg]=summary(lm.r)$coefficients[8]
        }
    }
    return(list(out=out,out_sig=out_sig))
}

combine_regional_quantile <- function(dat,dur,dur_mid,regions,regNumb){
    # value is an array of dim=c(1319,??)
    taus=c(0.9,0.95,0.98,0.99)
    out=array(NA,c(length(taus),regNumb))
    out_sig=array(NA,c(length(taus),regNumb))
    for (reg in 1:regNumb){
        inside=which(regions==reg)
        duration=array(NA,dim=c(100000))
        duration_mid=array(NA,dim=c(100000))
        count=1
        for (i in inside){
            values=length(which(!is.na(dur[i,])))
            duration[count:(count+values)]=dur[i,1:values]
            duration_mid[count:(count+values)]=dur_mid[i,1:values]
            count=count+values
        }
        duration=duration[!is.na(duration)]
        duration_mid=duration_mid[!is.na(duration_mid)]
        ord=order(duration_mid)
        if (length(duration)>1000){
            for (i in 1:length(taus)){
                y=as.vector(duration[ord])
                x=as.vector(duration_mid[ord])
                qu=summary(rq(y~x,taus[i]))$coefficients
                out[i,reg]=qu[2]
                out_sig[i,reg]=qu[8]
            }
        }
    }
    return(list(out=out,out_sig=out_sig))
}

regional_analysis <- function(dat){
    nday=91
    nyr=5
    ntot=length(dat$ID)
    library(Kendall)
    library(quantreg)

    nc_mar=open.ncdf("../data/91_5/91_5_markov2s.nc")
    IDregions=read.table("../data/ID-regions.txt")
    region_names=c("srex","7rect")
    regs=c(26,7)

    result=array(NA,dim=c(4,2,10,(26+7)))
    sig=array(NA,dim=c(4,2,10,(26+7)))

    season_names=c("spring","summer","autumn","winter")
    for (season in 1:4){    
        nc_dur=open.ncdf(paste("../data/",nday,"_",nyr,"/",nday,"_",nyr,"_duration_2s_",season_names[season],".nc",sep=""))
        dur=get.var.ncdf(nc_dur,"dur")
        dur_mid=get.var.ncdf(nc_dur,"dur_mid")
        per=get.var.ncdf(nc_mar,paste("markov_",season_names[season],sep=""))
        for (state in 1:2){
            tmp=average_regional_trend(dat,per[1:ntot,(state*state),],IDregions[1:ntot,1],regNumb=regs[1])
            result[season,state,1:2,1:26]=tmp$out
            sig[season,state,1:2,1:26]=tmp$out_sig
            tmp=average_regional_trend(dat,per[1:ntot,(state*state),],IDregions[1:ntot,2],regNumb=regs[2])
            result[season,state,1:2,27:33]=tmp$out
            sig[season,state,1:2,27:33]=tmp$out_sig

            tmp=combine_regional_quantile(dat,dur[1:ntot,state,],dur_mid[1:ntot,state,],IDregions[1:ntot,1],regNumb=regs[1])
            result[season,state,3:6,1:26]=tmp$out
            sig[season,state,3:6,1:26]=tmp$out_sig
            tmp=combine_regional_quantile(dat,dur[1:ntot,state,],dur_mid[1:ntot,state,],IDregions[1:ntot,2],regNumb=regs[2])
            result[season,state,3:6,27:33]=tmp$out
            sig[season,state,3:6,27:33]=tmp$out_sig
            print(tmp)
        }
    }
    regional_analysis_write("../data/test.nc",result,sig)

}






map_regional <- function(dat,toPlot1,toPlot2,titles,worldmap,filename_plot){
    # dat from data_load, loaded with reg=1!
    # toPlot array(... ,dim=c(number of maps, number of stations))
    # titles titles of plots (list)

    ntot=length(dat$ID)
    library(Kendall)
    IDregions=read.table("../data/ID-regions.txt")
    region_names=c("srex","7rect")
    for (k in 2:3){
        pdf(file = paste(filename_plot,"_",region_names[(k-1)],".pdf",sep=""),width=12,height=8)

        for (i in 1:dim(toPlot1)[1]){
            poli=read.table(paste("../data/",region_names[(k-1)],".txt",sep=""))
            tmp=average_regional_trend(dat,toPlot1[i,1:ntot,],IDregions[1:ntot,k],regNumb=dim(poli)[1])
            print(paste("MK",titles[i]))
            regions_color(tmp$out[1,],tmp$out_sig[1,],worldmap,paste("MK",titles[i]),poli)
            print(paste("LR",titles[i]))
            regions_color(tmp$out[2,],tmp$out_sig[2,],worldmap,paste("LR",titles[i]),poli)
        }
        graphics.off()
    }
}

