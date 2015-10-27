
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




