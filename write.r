library(RNetCDF)


trend_write <- function(filename,trend,ID_length=1319,method="2D running mean")
{
    nc_out <- create.nc(filename)

    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"day",dimlength=365,unlim=FALSE)
    dim.def.nc(nc_out,"year",dimlength=65,unlim=FALSE)

    var.def.nc(nc_out,"trend","NC_DOUBLE",c(0,1,2))
    att.put.nc(nc_out, "trend", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "trend", "dim_explanation", "NC_CHAR", "ID-day-year")
    att.put.nc(nc_out, "trend", "method", "NC_CHAR",method)

    var.put.nc(nc_out,"trend",trend)

    close.nc(nc_out)  
}





duration_write <- function(filename,dur,dur_mid,len,ID_length=1319)
{
    nc_out <- create.nc(filename)
    
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)

    dim.def.nc(nc_out,"periods",dimlength=len,unlim=FALSE)

    var.def.nc(nc_out,"dur","NC_INT",c(0,1,2))
    att.put.nc(nc_out, "dur", "missing_value", "NC_INT", 99999)
    att.put.nc(nc_out, "dur", "dim_explanation", "NC_CHAR", "ID-seasons-states-...")
    att.put.nc(nc_out, "dur", "val_explanation", "NC_CHAR", "length of persistent period")

    var.def.nc(nc_out,"dur_mid","NC_DOUBLE",c(0,1,2))
    att.put.nc(nc_out, "dur_mid", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "dur_mid", "dim_explanation", "NC_CHAR", "ID-seasons-states-...")
    att.put.nc(nc_out, "dur_mid", "val_explanation", "NC_CHAR", "mid-point of persistent period")

    var.put.nc(nc_out,"dur",dur)              
    var.put.nc(nc_out,"dur_mid",dur_mid)              

    close.nc(nc_out) 
}



quantiles_write <- function(filename,ID_length,ID_name,period,taus,quantile_stuff,comment="quantile_analysis"){

    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "comment", "NC_CHAR", comment)
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)
    
    dim.def.nc(nc_out,"seasons",dimlength=6,unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)

    dim.def.nc(nc_out,"taus",dimlength=length(taus),unlim=FALSE)
    dim.def.nc(nc_out,"quant_outs",dimlength=3,unlim=FALSE)
        

    var.def.nc(nc_out,"quantile_stuff","NC_DOUBLE",c(0,1,2,3,4))
    att.put.nc(nc_out, "quantile_stuff", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "quantile_stuff", "dim_explanation", "NC_CHAR", "season-ID-state-...")
    att.put.nc(nc_out, "quantile_stuff", "explanation", "NC_CHAR", "(0.05,0.25,0.5,0.75,0.95,0.98,1) x (quantiles, slopes, slope_sigs)")
        
    var.put.nc(nc_out,"quantile_stuff",quantile_stuff) 

    close.nc(nc_out) 
}


other_write <- function(filename,ID_length,ID_name,period,other_stuff,comment="other_analysis"){

    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "comment", "NC_CHAR", comment)
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)

    dim.def.nc(nc_out,"seasons",dimlength=6,unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)
  
    dim.def.nc(nc_out,"outs",dimlength=12,unlim=FALSE)

    var.def.nc(nc_out,"other_stuff","NC_DOUBLE",c(0,1,2,3))
    att.put.nc(nc_out, "other_stuff", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "other_stuff", "dim_explanation", "NC_CHAR", "season-ID-state-...")
    att.put.nc(nc_out, "other_stuff", "explanation", "NC_CHAR", "mean,sd,summary(lm)$coef")
        
    var.put.nc(nc_out,"other_stuff",other_stuff)      
 
    close.nc(nc_out) 
}

fit_write <- function(filename,ID_length,ID_name,period,fit_stuff,comment="distribution_fits"){

    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "comment", "NC_CHAR", comment)
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)

    dim.def.nc(nc_out,"seasons",dimlength=6,unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)

    dim.def.nc(nc_out,"fit_outs",dimlength=20,unlim=FALSE)

    var.def.nc(nc_out,"fit_stuff","NC_DOUBLE",c(0,1,2,3))
    att.put.nc(nc_out, "fit_stuff", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "fit_stuff", "dim_explanation", "NC_CHAR", "season-ID-state-...")
    att.put.nc(nc_out, "fit_stuff", "explanation", "NC_CHAR", "first values: parameters, 19=R^2, 20=BIC")
    att.put.nc(nc_out, "fit_stuff", "explanation", "NC_CHAR", "first values: parameters, 19=R^2, 20=BIC")
        
    var.put.nc(nc_out,"fit_stuff",fit_stuff)      
 
    close.nc(nc_out) 
}