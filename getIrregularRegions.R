getIrregularRegions <- function( region.table, lon, lat, region.name, nregions, missval.reg ){
################################################################################
# This function extract the region information from the region file.
################################################################################

### create uniform lat lon grid
xlon = lon
xlat = lat
nlon <- length(xlon)
nlat <- length(xlat)

# define longitude value where boxes should be drawn
if(dat == "HadEX2"){
    start.lon = -11.25
}else if(dat == "HadEX2.Rx5"){
    start.lon = 176.25
}else if(dat == "GPCC_test"){
    start.lon = -12.5
}else{
    stop("error in getIrregularRegions: start.lon needs to be defined for this data set")
}
ind.lon = which(lon==start.lon)

######## find points inside the irregular polygon
reg.index = list()
for (i in 1:nregions){
    reg                   = region.name[i]
    reg.index[[reg]]      = list()
    reg.index[[reg]]$name = reg
    poly.x                = c()
    poly.y                = c()

    for( j in c(3, 5, 7, 9, 11, 13, 15) ){
	if( as.numeric(region.table[i,j]) == missval.reg ){
	    next
	}else{
	    poly.x = c(poly.x, as.numeric(region.table[i,j+1]))
	    poly.y = c(poly.y, as.numeric(region.table[i,j]))
	}
    }

    #### get the lat lon datapoints belonging to the region
    lat.mat <- t(matrix(xlat, nlat, nlon))
    lon.mat <- t(matrix(xlon, nlat, nlon, byrow=T))
    pp <- matrix( point.in.polygon(lon.mat, lat.mat, poly.x, poly.y,
                  mode.checked=FALSE), ncol = nlat, byrow = F )
    pp[pp==0] = NA
    pp[pp!=0] = 1
    reg.index[[reg]]$pp  = pp
    reg.index[[reg]]$lon = poly.x
    reg.index[[reg]]$lat = poly.y

    ## for full-lon-averaged regions also create small area that can be plotted in worldmaps
    poly.y = c(-90, -62.5, -62.5, -90)
    if( any(reg %in% c("TRO", "NST", "SST", "NET")) ){
	if(reg == "NET") poly.x = c(lon[ind.lon+9], lon[ind.lon+9], lon[ind.lon+15], lon[ind.lon+15])
	if(reg == "NST") poly.x = c(lon[ind.lon+18], lon[ind.lon+18], lon[ind.lon+24], lon[ind.lon+24])
	if(reg == "TRO") poly.x = c(lon[ind.lon+27], lon[ind.lon+27], lon[ind.lon+33], lon[ind.lon+33])
	if(reg == "SST") poly.x = c(lon[ind.lon+36], lon[ind.lon+36], lon[ind.lon+42], lon[ind.lon+42])
# 	poly.x = c(-178.125, -178.125, -174.375, -174.375)
	pp <- matrix( point.in.polygon(lon.mat, lat.mat, poly.x, poly.y,
		      mode.checked=FALSE), ncol = nlat, byrow = F )
	pp[pp==0] = NA
	pp[pp!=0] = 1
	reg.index[[reg]]$extra.pp  = pp
	reg.index[[reg]]$extra.lon = poly.x
	reg.index[[reg]]$extra.lat = poly.y
    }

}

## add region = "global"
reg.index[["global"]]$lat  = c(lat[1], lat[length(lat)], lat[length(lat)], lat[1])
reg.index[["global"]]$lon  = c(lon[1], lon[1], lon[length(lon)], lon[length(lon)])
reg.index[["global"]]$name = "global"
reg.index[["global"]]$pp   = matrix(1, length(lon), length(lat))
# poly.x = c(-163.125, -163.125, -148.125, -148.125)
# poly.y = c(-6.75, 6.75, 6.75, -6.75)
poly.x = c(start.lon, start.lon, lon[ind.lon+6], lon[ind.lon+6])
poly.y = c(-90, -62.5, -62.5, -90)
pp <- matrix( point.in.polygon(lon.mat, lat.mat, poly.x, poly.y,
	      mode.checked=FALSE), ncol = nlat, byrow = F )
pp[pp==0] = NA
pp[pp!=0] = 1
reg.index[["global"]]$extra.pp   = pp
reg.index[["global"]]$extra.lon  = poly.x
reg.index[["global"]]$extra.lat  = poly.y

return( reg.index )
}
