
if(make_box_wisker){
  # box-wisker plot for all regions and all seasons
#  load("data3D")
  in_celcius = FALSE
  per.decade = TRUE

  # these should be stored in a list
  qr_list <- list()

  ii = ANN_ii
  experiment_tag = paste(".day", ii[1],"-",ii[length(ii)],".",years[1],"-",years[length(years)],".",sep="")
  file = paste("HadGHCND_TX_qr_data3D",experiment_tag,"RData",sep="")
  load(paste(dir,file,sep=""))
  qr_list <- lappend(qr_list, qr_data)

  for(s in 1:4){
    ii = seasonal_ii[[s]]
    experiment_tag = paste(".day", ii[1],"-",ii[length(ii)],".",years[1],"-",years[length(years)],".",sep="")
    file = paste("HadGHCND_TX_qr_data3D",experiment_tag,"RData",sep="")
    load(paste(dir,file,sep=""))
    qr_list <- lappend(qr_list, qr_data)
  }
# load(paste(dir,"/winter/","HadGHCND_TX_qr_data3D.day335-59.1950-2011.RData",sep=""))
#   qr_list <- lappend(qr_list, qr_data)
# load(paste(dir,"/spring/","HadGHCND_TX_qr_data3D.day60-151.1950-2011.RData",sep=""))
#   qr_list <- lappend(qr_list, qr_data)
# load(paste(dir,"/summer/","HadGHCND_TX_qr_data3D.day152-243.1950-2011.RData",sep=""))
#   qr_list <- lappend(qr_list, qr_data)
# load(paste(dir,"/autumn/","HadGHCND_TX_qr_data3D.day244-334.1950-2011.RData",sep=""))
#   qr_list <- lappend(qr_list, qr_data)
  # names for labels
  lbls = c("ANN","DJF","MAM","JJA","SON")

  # revert
  qr_list <- rev(qr_list)
  lbls    <- rev(lbls)

  # open new pdf
pdf(file=paste(dir,"/","box.and.wisker.3regions.in.perdecade.pdf",sep=""),height=3,width=10)
#  par(mfrow=c(4,2))
  par(mfrow=c(1,3))
  par(xaxs = "i")
  par(yaxs = "i")
  ylim_ <- c(-0.05,0.05)
  if(in_celcius) ylim_ <- c(-3,3)
  if(per.decade) ylim_ <- c(-1,1)
#  allregions2 = allregions[2:9] # exclude "global"
  allregions2 = c(allregions[2],allregions[8],allregions[9]) # mid-lats, america, eurasia

  # loop over all regions: 1 plot per region
  for(reg in allregions2){
    # for each region, plot for annual and each season (5 box-plots)
    # the upper-tail, lower-tail and full distribution (3 box-plots)
    #    plot.new()
    # 1) Full distribution
    # loop over annual + seasons and add to list "widths"
    widths <- list()
    for(i in 1:length(qr_list)){
      region <- selectRegion3D(data3D$grid, qr_list[[i]], region=reg)
      widening <- region$q_trend[5,] - region$q_trend[1,]
      if(in_celcius) widening <- widening*length(years)
      if(per.decade) widening <- widening*10.
      widths <- lappend(widths, widening)
    }
    # make a box-wisker plot
    at_=seq(1,5,1) # center
    boxplot(widths,names=lbls,range = 1.5,boxwex=0.15,outline=FALSE,
            border="black",col="grey", horizontal=TRUE,axes=TRUE,
            las=1,xaxt="n",add=FALSE, at=at_, ylim=ylim_)
    axis(side=1) # add x-axis

    # 2) Upper quartile
    # loop over annual + seasons and add to list "widths"
    widths <- list()
    for(i in 1:length(qr_list)){
      region <- selectRegion3D(data3D$grid, qr_list[[i]], region=reg)
      widening <- region$q_trend[5,] - region$q_trend[3,]
      if(in_celcius) widening <- widening*length(years)
      if(per.decade) widening <- widening*10.
      widths <- lappend(widths, widening)
    }
    # add a box-wisker plot
    at_= at_ +0.2 # above the full distribution
    boxplot(widths,names=FALSE,range = 1.5,boxwex=0.15,outline=FALSE,
            border="black",col="red", horizontal=TRUE,axes=FALSE,
            add=TRUE, at=at_)

    # 3) Lower quartile
    # loop over annual + seasons and add to list "widths"
    widths <- list()
    for(i in 1:length(qr_list)){
      region <- selectRegion3D(data3D$grid, qr_list[[i]], region=reg)
      widening <- region$q_trend[3,] - region$q_trend[1,]
      if(in_celcius) widening <- widening*length(years)
      if(per.decade) widening <- widening*10.
      widths <- lappend(widths, widening)
    }
    # add a box-wisker plot
    at_=seq(1,5,1) # center
    at_= at_ - 0.2 # above the full distribution
    boxplot(widths,names=FALSE,range = 1.5,boxwex=0.15,outline=FALSE,
            border="black",col="blue", horizontal=TRUE,axes=FALSE,
            add=TRUE, at=at_)
    abline(h=seq(1.5,4.5,1),lty=3,lwd=0.5)
    mtext(reg)
    at_=seq(1,5,1) # center
    x0 = array(0,dim=4)
    y0 = seq(1.375,4.375,1)
    y1 = y0+0.25
    segments(x0=x0,y0=y0,x1=x0,y1=y1)

  }
  dev.off()
} # end if-loop "box-and-wisker plots"









#---------------------------nearest neighbours

if (1==2){
    if (plot[1]==1){
        color=c("red","blue","green","violet","black","orange","lightblue","grey",rgb(1,0.5,0.6),rgb(0.5,0.7,0.9),rgb(0.5,0.2,0.8),rgb(0.2,0.5,0.6))
        jet.colors <- colorRampPalette( c( "blue","green","yellow","red") )
        nbcol <- nGroup
        color <<- jet.colors(nbcol)  
        reihen=array(NA,dim=c(versions,ntot))
        pdf(file=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbours/",period,"/",trendID,"_",region_name,"_",period,"nearest_neighbours",name_zusatz,".pdf",sep=""))
    }


        if (reduce>0){
            tmp=group_reduction(X=X,attribution=attribution,start=start,nGroup=nGroup,reduce=reduce,main_add=version)
            attribution=tmp$attribution
            start=tmp$start
        }

        if (plot[2]==1){
            for (p in 1:nGroup){
                if (!is.na(start[p,1])){
                    plot(NA,xlab="days",ylab="probability density",ylim=c(0.00001,0.25),xlim=c(0,50),axes=TRUE,frame.plot=TRUE,main=length(which(attribution==p))) 
                    for (q in which(attribution==p)){
                        #print("missing points")
                        points(X,toOrder[q,],pch=16,col=rgb(0.5,0.5,0.5,0.2),cex=1.5)
                    }
                    lines(X,start[p,],col=color[p])

                    plot(NA,xlab="days",ylab="probability density",ylim=c(0.001,0.25),xlim=c(0,50),axes=TRUE,frame.plot=TRUE,log="y",main=length(which(attribution==p))) 
                    for (q in which(attribution==p)){
                        #print("missing points")
                        points(X,toOrder[q,],pch=16,col=rgb(0.5,0.5,0.5,0.2),cex=1.5)
                    }
                    lines(X,start[p,],col=color[p])
                }
            }
        }

    if (plot[3]==1){map_allgemein(dat=dat,filename_plot=paste("../plots/",trendID,"/",dataset,additional_style,"/nearest_neighbours/",period,"/",trendID,"_",region_name,"_",period,"nearest_neighbours_map",name_zusatz,".pdf",sep=""),worldmap=worldmap,reihen=reihen,pointsize=1.5,farb_palette="regenbogen")}
    
    if (write==TRUE){
        nc_out <- create.nc(paste("../plots/",trendID,"/",dataset,additional_style,"/regions/",period,"/",trendID,"_",region_name,"_",period,"nearest_neighbours_map",name_zusatz,"result.nc",sep=""))
        att.put.nc(nc_out, "NC_GLOBAL", "comment", "NC_CHAR", "bla")

        dim.def.nc(nc_out,"seasons",dimlength=6,unlim=FALSE)
        dim.def.nc(nc_out,"ID",dimlength=ntot, unlim=FALSE)
        dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)

        dim.def.nc(nc_out,"versions",dimlength=versions,unlim=FALSE)
        dim.def.nc(nc_out,"nGroup",dimlength=nGroup,unlim=FALSE)
        dim.def.nc(nc_out,"distrSize",dimlength=distrSize,unlim=FALSE)


        var.def.nc(nc_out,"attribution","NC_DOUBLE",c(3,1))
        att.put.nc(nc_out, "attribution", "missing_value", "NC_DOUBLE", -99999.9)
        att.put.nc(nc_out, "attribution", "dim_explanation", "NC_CHAR", "version-ID")
        att.put.nc(nc_out, "attribution", "explanation", "NC_CHAR", "group numbers of differnt versions")

        var.def.nc(nc_out,"groups","NC_DOUBLE",c(3,4,5))
        att.put.nc(nc_out, "groups", "missing_value", "NC_DOUBLE", -99999.9)
        att.put.nc(nc_out, "groups", "dim_explanation", "NC_CHAR", "version-nGroup-distrSize")
        att.put.nc(nc_out, "groups", "explanation", "NC_CHAR", "end distributions of groups")
            
        var.put.nc(nc_out,"attribution",attribution_speicher)      
        var.put.nc(nc_out,"groups",start_speicher)      
     
        close.nc(nc_out) 
    }

        if (FALSE==TRUE){
            laggedDiffs=c()
            for (p in 1:nGroup){
                laggedDiffs[p]=sum(diff(start[p,3:10]))
            }
            print(laggedDiffs)
            order=order(laggedDiffs)
            print(order)
            
            order=order(start[,1])

            attri_order=attribution*NA
            for (i in 1:length(order)){
                attri_order[attribution==i]=order[i]
            }
            reihen[version,]=attri_order
        }

        if (plot[3]==1){reihen[1,]=attribution}

        if (versions!=1){
            attribution_speicher[version,]=attribution
            start_speicher[version,,]=start
        }
}
