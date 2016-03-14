
trend_write <- function(filename,trend,ID_length=length(dat$ID),method="2D running mean")
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


duration_write <- function(filename,dur,dur_mid,len,ID_length=length(dat$ID),ID_name="grid_points",comment="no comment")
{
    print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "comment", "NC_CHAR", comment)
    
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)

    dim.def.nc(nc_out,"periods",dimlength=len,unlim=FALSE)

    var.def.nc(nc_out,"dur","NC_INT",c(0,1,2))
    att.put.nc(nc_out, "dur", "missing_value", "NC_INT", 99999)
    att.put.nc(nc_out, "dur", "dim_explanation", "NC_CHAR", "ID-states-...")
    att.put.nc(nc_out, "dur", "val_explanation", "NC_CHAR", "length of persistent period")

    var.def.nc(nc_out,"dur_mid","NC_DOUBLE",c(0,1,2))
    att.put.nc(nc_out, "dur_mid", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "dur_mid", "dim_explanation", "NC_CHAR", "ID-states-...")
    att.put.nc(nc_out, "dur_mid", "val_explanation", "NC_CHAR", "mid-point of persistent period")

    var.put.nc(nc_out,"dur",dur)              
    var.put.nc(nc_out,"dur_mid",dur_mid)              

    close.nc(nc_out) 
}

reg_binned_dur_write <- function(filename,binned_dur,len,ID_length=length(dat$ID),ID_name="grid_points",comment="no comment"){
    print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "comment", "NC_CHAR", comment)
    
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)

    dim.def.nc(nc_out,"periods",dimlength=len,unlim=FALSE)

    dim.def.nc(nc_out,"years",dimlength=length(dat$year),unlim=FALSE)

    var.def.nc(nc_out,"binned_dur","NC_INT",c(0,1,2,3))
    att.put.nc(nc_out, "binned_dur", "missing_value", "NC_INT", 99999)
    att.put.nc(nc_out, "binned_dur", "dim_explanation", "NC_CHAR", "ID-states-years-...")
    att.put.nc(nc_out, "binned_dur", "val_explanation", "NC_CHAR", "length of persistent period as set for each year")

    var.put.nc(nc_out,"binned_dur",binned_dur)              

    close.nc(nc_out) 
}

quantiles_write <- function(filename,ID_length,ID_name,period,taus,quantile_stuff,comment="quantile_analysis"){
    print(filename)
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
    att.put.nc(nc_out, "quantile_stuff", "explanation", "NC_CHAR", "(0.75,0.95,0.99) x (quantiles, slopes, slope_sigs)")
        
    var.put.nc(nc_out,"quantile_stuff",quantile_stuff) 

    close.nc(nc_out) 
}


other_write <- function(filename,ID_length,ID_name,period,other_stuff,comment="other_analysis"){

    print(filename)
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

    print(filename)
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
    att.put.nc(nc_out, "fit_stuff", "explanation", "NC_CHAR", "expo: 1-a, 2-b, 15-R2, 16-BIC special fit: 5-a1, 6-b1, 7-a2, 8-b2, 9-thresh 19=R^2, 20=BIC  special: 14-distr_length, 17-BIc_diff")
        
    var.put.nc(nc_out,"fit_stuff",fit_stuff)      
 
    close.nc(nc_out) 
}

distr_write <- function(distr_stuff,filename,ID_length,ID_name,period,comment="distribution_fits"){

    print(filename)
    nc_out <- create.nc(filename)
    att.put.nc(nc_out, "NC_GLOBAL", "ID_explanation", "NC_CHAR", ID_name)
    att.put.nc(nc_out, "NC_GLOBAL", "comment", "NC_CHAR", comment)
    att.put.nc(nc_out, "NC_GLOBAL", "period", "NC_CHAR", period)

    dim.def.nc(nc_out,"seasons",dimlength=6,unlim=FALSE)
    dim.def.nc(nc_out,"ID",dimlength=ID_length, unlim=FALSE)
    dim.def.nc(nc_out,"states",dimlength=2,unlim=FALSE)

    dim.def.nc(nc_out,"distr_outs",dimlength=5,unlim=FALSE)
    dim.def.nc(nc_out,"distr_length",dimlength=100,unlim=FALSE)

    var.def.nc(nc_out,"distr_stuff","NC_DOUBLE",c(0,1,2,3,4))
    att.put.nc(nc_out, "distr_stuff", "missing_value", "NC_DOUBLE", -99999.9)
    att.put.nc(nc_out, "distr_stuff", "dim_explanation", "NC_CHAR", "season-ID-state-(X,Y,counts,expifit,specialfit)-...")
    att.put.nc(nc_out, "distr_stuff", "explanation", "NC_CHAR", "distribution values: X, Y, counts, exponential fit, special fit")
        
    var.put.nc(nc_out,"distr_stuff",distr_stuff)      
 
    close.nc(nc_out) 
}