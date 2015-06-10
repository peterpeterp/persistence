
dat_write = function(filename,data3D)
{
    nc = create.nc(filename)
    dim.def.nc(nc, "day",  dimlength=length(data3D$day),  unlim=FALSE)
    dim.def.nc(nc, "year", dimlength=length(data3D$year), unlim=FALSE)
    dim.def.nc(nc, "ID",   dimlength=length(data3D$ID), unlim=FALSE)

    var.def.nc(nc, "day", "NC_INT", "day")
    var.def.nc(nc, "year", "NC_INT", "year")
    var.def.nc(nc, "ID", "NC_INT", "ID")

    var.def.nc(nc, "lon", "NC_FLOAT", "ID")
    var.def.nc(nc, "lat", "NC_FLOAT", "ID")
 
    var.def.nc(nc, "tas", "NC_FLOAT", c("ID","day","year"))
    att.put.nc(nc, "tas", "long_name", "NC_CHAR", "Near-surface air temperature anomaly")
    att.put.nc(nc, "tas", "units", "NC_CHAR", "degrees Celcius")
    att.put.nc(nc, "tas", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "day",  data3D$day)
    var.put.nc(nc, "year", data3D$year)
    var.put.nc(nc, "ID",   data3D$ID)

    var.put.nc(nc, "lon",  data3D$lon)
    var.put.nc(nc, "lat",  data3D$lat)

    var.put.nc(nc, "tas",  data3D$tas)

    close.nc(nc)
}

trend_write = function(filename,data3D,trend)
{
    nc = create.nc(filename)
    dim.def.nc(nc, "day",  dimlength=length(data3D$day),  unlim=FALSE)
    dim.def.nc(nc, "year", dimlength=length(data3D$year), unlim=FALSE)
    dim.def.nc(nc, "ID",   dimlength=length(data3D$ID), unlim=FALSE)

    var.def.nc(nc,"trend","NC_FLOAT", c("ID","day","year"))
    att.put.nc(nc, "trend", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "trend", trend)

    close.nc(nc)
}

per_write = function(filename,data3D,per)
{
    nc = create.nc(filename)
    dim.def.nc(nc, "day",  dimlength=length(data3D$day),  unlim=FALSE)
    dim.def.nc(nc, "year", dimlength=length(data3D$year), unlim=FALSE)
    dim.def.nc(nc, "ID",   dimlength=length(data3D$ID), unlim=FALSE)

    var.def.nc(nc, "ma_warm_trend", "NC_FLOAT", "ID")
    att.put.nc(nc, "ma_warm_trend", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc, "ma_cold_trend", "NC_FLOAT", "ID")
    att.put.nc(nc, "ma_cold_trend", "missing_value", "NC_FLOAT", -9999.0)


    var.def.nc(nc,"ind","NC_FLOAT", c("ID","day","year"))
    att.put.nc(nc, "ind", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"perp","NC_FLOAT", c("ID","day","year"))
    att.put.nc(nc, "perp", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"dayp","NC_FLOAT", c("ID","day","year"))
    att.put.nc(nc, "dayp", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"pern","NC_FLOAT", c("ID","day","year"))
    att.put.nc(nc, "pern", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"dayn","NC_FLOAT", c("ID","day","year"))
    att.put.nc(nc, "dayn", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"ma_warm","NC_FLOAT", c("ID","day","year"))
    att.put.nc(nc, "ma_warm", "missing_value", "NC_FLOAT", -9999.0)

    var.def.nc(nc,"ma_cold","NC_FLOAT", c("ID","day","year"))
    att.put.nc(nc, "ma_cold", "missing_value", "NC_FLOAT", -9999.0)

    var.put.nc(nc, "ma_warm_trend",  per$ma_cold_trend)
    var.put.nc(nc, "ma_cold_trend",  per$ma_cold_trend)

    var.put.nc(nc, "ind", per$ind)
    var.put.nc(nc, "perp", per$perp)
    var.put.nc(nc, "dayp", per$dayp)
    var.put.nc(nc, "pern", per$pern)
    var.put.nc(nc, "dayn", per$dayn)
    var.put.nc(nc, "ma_warm", per$ma_warm)
    var.put.nc(nc, "ma_cold", per$ma_cold)

    close.nc(nc)

}

