import os,sys,glob
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date




def grid_to_ID(ID,lat,lon,in_file,out_file,variables,time,year,month):	
	nc_in=Dataset(in_file,'r')

	lat_=nc_in.variables['lat'][:]   ;   ny=len(lat)
	lon_=nc_in.variables['lon'][:]   ;   nx=len(lon)

	for variable in variables:
		tas=nc_in.variables[variable][:,:,:]
		mask=np.ma.getmask(tas)
		tas=np.ma.getdata(tas)
		tas[mask]=np.nan



		tas_new=np.zeros([nt,len(ID)])*np.nan
		count=0
		for x in range(len(lon_)):
			for y in range(len(lat_)):
				i = np.where((lon==lon_[x]) & (lat==lat_[y]))[0]
				if len(i)==1:
					tas_new[:,i[0]]=tas[:,y,x]


		
		os.system('rm '+out_file)
		nc_out=Dataset(out_file,'w')
		nc_out.createDimension('time',nt)
		nc_out.createDimension('ID',len(ID))



		outVar = nc_out.createVariable('ID', 'i2', 'ID')    ;   outVar[:] = ID
		outVar = nc_out.createVariable('lat', 'f', 'ID')    ;   outVar[:] = lat
		outVar = nc_out.createVariable('lon', 'f', 'ID')    ;   outVar[:] = lon

		outVar = nc_out.createVariable('time', 'i2', 'time')    ;   outVar[:] = time
		outVar = nc_out.createVariable('year', 'i2', 'time')    ;   outVar[:] = year
		outVar = nc_out.createVariable('month', 'i2', 'time')    ;   outVar[:] = month

		outVar = nc_out.createVariable(variable, 'f', ('time','ID',))    ;   outVar[:] = tas_new[:,:]

		nc_out.close()


	nc_in.close()


# duration meta
in_file='/home/peter/Dokumente/pik/backuped/data/_TMean/HadGHCND_TMean_data3D.day1-365.1950-2014.nc'
nc_raw=Dataset(in_file,'r')

ID=nc_raw.variables['ID'][:]
lon=nc_raw.variables['lon'][:]
lat=nc_raw.variables['lat'][:]

# in_file='/home/peter/Dokumente/pik/backuped/data/spi/spi_cru_ts_3_1949-2012_96x73.nc4'
# out_file=in_file.replace('.nc','_ID.nc').replace('spi_',variable+'_')
# nc_in=Dataset(in_file,'r')
# # handle time
# time=nc_in.variables['time'][:]
# nt=len(time)
# datevar = []
# time_unit=nc_in.variables['time'].units
# try:    
# cal_temps = nc_in.variables['time'].calendar
# datevar.append(num2date(time,units = time_unit,calendar = cal_temps))
# except:
# datevar.append(num2date(time,units = time_unit))
# # create index variable
# year=np.array([int(str(date).split("-")[0])    for date in datevar[0][:]])
# month=np.array([int(str(date).split("-")[1])    for date in datevar[0][:]])
# grid_to_ID(ID,lat,lon,in_file,out_file,variables=['spi3','spi6','spi12'],time=time,year=year,month=month)

in_file='/home/peter/Dokumente/pik/backuped/data/sonstiges/eke_roh/eke850_1979-2014_calendar_96x73.nc'
out_file=in_file.replace('.nc','_ID.nc')
nc_in=Dataset(in_file,'r')
time=nc_in.variables['time'][:]
year=time/12+1979
month=(time/12.-time/12)*12+1

grid_to_ID(ID,lat,lon,in_file,out_file,variables=['eke850'],time=time,year=year,month=month)





# # get eke
# in_file='/home/peter/Dokumente/pik/backuped/data/sonstiges/eke_roh/EKE_ERA_Interim_1979-2014_calendar_96x73.nc'

# nc_in=Dataset(in_file,'r')
# # handle time
# time=nc_in.variables['time'][:]

# lat_=nc_in.variables['lat'][:]
# lon_=nc_in.variables['lon'][:]

# levelist=nc_in.variables['levelist'][:]

# u2=nc_in.variables['u2syn'][:,0,:,:]
# v2=nc_in.variables['v2syn'][:,0,:,:]

# eke850=0.5*(u2+v2)

# out_file=in_file.replace('EKE_ERA_Interim_','eke850_')
# os.system('rm '+out_file)
# nc_out=Dataset(out_file,'w')

# nc_out.createDimension('time',len(time))
# nc_out.createDimension('lat',len(lat_))
# nc_out.createDimension('lon',len(lon_))

# outVar = nc_out.createVariable('lat', 'f', 'lat')    ;   outVar[:] = lat_
# outVar = nc_out.createVariable('lon', 'f', 'lon')    ;   outVar[:] = lon_

# outVar = nc_out.createVariable('time', 'i2', 'time')    ;   outVar[:] = time
# outVar.units="months since 1979-01-01"

# outVar = nc_out.createVariable('eke850', 'f', ('time','lat','lon',))    ;   outVar[:] = eke850[:,:,:]
# outVar.level='850mbar'

# nc_out.close()
