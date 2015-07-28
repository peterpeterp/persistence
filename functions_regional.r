
wave_region <- function(wavenumber,trough){
    wellen=array(NA,dim=c(40,13))
    brei=360/wavenumber/4
    lon=trough
    lon0=trough
    i=1
    wellen[,7:8]=35
    wellen[,9:10]=60
    while ((round(x=lon,digits=3)==round(x=lon0,digits=3) & i>1)==FALSE){
        if (lon>(180-brei) & lon<(180+brei)){
            wellen[i,1:4]=c(180,lon-brei,lon-brei,180)
            wellen[i,13]=i  
            i=i+1   
            wellen[i,1:4]=c(lon+brei-360,-180,-180,lon+brei-360)
            wellen[i,13]=i  
            i=i+1   
            lon=lon+2*brei
        }
        if (lon<=(180-brei)){
            wellen[i,1:4]=c(lon+brei,lon-brei,lon-brei,lon+brei)
            wellen[i,13]=i  
            i=i+1   
            lon=lon+2*brei
        }
        if (lon>=(180+brei)){
            lon=lon-360
            wellen[i,1:4]=c(lon+brei,lon-brei,lon-brei,lon+brei)
            wellen[i,13]=i  
            i=i+1   
            lon=lon+2*brei
        }
        
    }
    write.table(wellen[1:length(which(!is.na(wellen[,13]))),],paste("../data/",wavenumber,"wave.txt",sep=""))
}

points_to_regions <- function(dat,region_names=c("srex","7rect","6wave","7wave","8wave")){
    library(SDMTools)
    ntot=length(dat$ID)
    points=cbind(x=dat$lon,y=dat$lat)

    region = array(NA,dim=c(ntot,(length(region_names)+1)))

    for (j in 1:length(region_names)){
        poli=read.table(paste("../data/region_poligons/",region_names[j],".txt",sep=""))
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
            inside=pnt.in.poly(points,poligon)$pip
            region[which(inside==1),(j+1)]=poli[k,13]
        }
    }

    if (length(region_names)>1){
        region[1:ntot,1]=dat$ID
        write.table(region,"../data/ID-regions.txt")
    }
    return(region[1:ntot,2])
}

markov_regional_trend <- function(dat,value,regions,yearPeriod,regNumb){
    # value is an array of dim=c(1319,??)
    yearPeriod=yearPeriod-1949
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

    x=seq(yearPeriod[1],yearPeriod[2],1)
    for (reg in 1:regNumb){
        y=valOut[reg,x]
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

duration_regional_quantile <- function(dat,dur,dur_mid,regions,yearPeriod,regNumb){
    # value is an array of dim=c(1319,??)
    taus=c(0.75,0.9,0.95,0.98)
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
                print(paste(reg,taus[i]))
                y=as.vector(duration[ord])
                x=as.vector(duration_mid[ord])
                inYearPeriod=which(x>yearPeriod[1] & x<yearPeriod[2])
                y=y[inYearPeriod]
                x=x[inYearPeriod]

                qu=try(summary(rq(y~x,taus[i]))$coefficients)
                if (class(qu)=="try-error"){qu=array(NA,8)}

                out[i,reg]=qu[2]
                out_sig[i,reg]=qu[8]
            }
        }
    }
    return(list(out=out,out_sig=out_sig))
}

regional_analysis <- function(dat,yearPeriod,filepath,region_name){
    nday=91
    nyr=5
    ntot=length(dat$ID)
    library(Kendall)
    library(quantreg)

    nc_mar=open.ncdf("../data/91_5/2_states/markov/91_5_markov2s.nc")
    IDregions=read.table("../data/ID-regions.txt")

    poli=read.table(paste("../data/region_poligons/",region_name,".txt",sep=""))

    regs=dim(poli)[1]
    IDregions=points_to_regions(dat,c(region_name))
    print(IDregions)

    season_names=c("spring","summer","autumn","winter","year")
    for (season in 1:length(season_names)){   
        print(season_names[season]) 
        result=array(NA,dim=c(2,10,regs))
        sig=array(NA,dim=c(2,10,regs))
        nc_dur=open.ncdf(paste("../data/",nday,"_",nyr,"/2_states/duration/",nday,"_",nyr,"_duration_2s_",season_names[season],".nc",sep=""))
        dur=get.var.ncdf(nc_dur,"dur")
        dur_mid=get.var.ncdf(nc_dur,"dur_mid")
        per=get.var.ncdf(nc_mar,paste("markov_",season_names[season],sep=""))
        for (state in 1:2){
            tmp=markov_regional_trend(dat,per[1:ntot,(state*state),],IDregions,yearPeriod=yearPeriod,regNumb=regs)
            result[state,1:2,1:regs]=tmp$out
            sig[state,1:2,1:regs]=tmp$out_sig
            print(tmp)
            tmp=duration_regional_quantile(dat,dur[1:ntot,state,],dur_mid[1:ntot,state,],IDregions,yearPeriod=yearPeriod,regNumb=regs)
            result[state,3:6,1:regs]=tmp$out
            sig[state,3:6,1:regs]=tmp$out_sig
            print(tmp)

        }
        regional_analysis_write(paste(filepath,season_names[season],"_",region_name,".nc",sep=""),result,sig,poli)
    }


}

regions_color <- function(reihen,reihen_sig,worldmap,titles,filename_plot){
    jet.colors <- colorRampPalette( c( "violet","blue","white","yellow","red") )
    nbcol <- 101
    color <- jet.colors(nbcol)

    IDregions=read.table("../data/ID-regions.txt")
    region_names=c("srex","7rect")
    for (k in 1:2){
        pdf(file = paste(filename_plot,"_",region_names[k],".pdf",sep=""),width=12,height=8)
        poli=read.table(paste("../data/",region_names[k],".txt",sep=""))

        for (rei in 1:dim(reihen)[1]){            
            y=c()
            index=c()
            signi=c()
            j=0
            for (i in 1:dim(poli)[1]){
                poliLabel=poli[i,13]
                if (!is.na(reihen[rei,poliLabel])){
                    j=j+1
                    y[j]=reihen[rei,poliLabel]
                    index[j]=i 
                    if (abs(reihen[rei,poliLabel])>0.0001){
                        signi[j]=sprintf("%.04f",reihen_sig[rei,poliLabel])
                    }         
                }
            }
            aushol=max(c(abs(max(y)),abs(min(y))))
            y[j+1]=-aushol
            y[j+2]=aushol
            facetcol <- cut(y,nbcol)  

            print(titles[rei])
            plot(worldmap,main=titles[rei])

            for (i in 1:j){
                lon=poli[index[i],1:6]
                lat=poli[index[i],7:12]
                lon=lon[!is.na(lon)]
                lat=lat[!is.na(lat)]
                polygon(x=lon,y=lat,col=color[facetcol[i]],border="green")
                text(mean(lon),mean(lat),label=signi[i],cex=0.7,col="black")
            }
            image.plot(legend.only=T, zlim=range(y), col=color)
        }
    graphics.off()
    }
    return()
}



