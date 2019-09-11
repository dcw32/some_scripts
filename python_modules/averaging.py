from netCDF4 import Dataset
from gridbox_areas import gbox_areas
import numpy as np
def time_sa(file,var,area):
	field=file.variables[var][:].squeeze()
	field=np.mean(field,axis=0)
	field_av=np.average(field,weights=area)
	return field_av
def monthly_to_annual(field,lats,lons,tims,axix_no):
	area=gbox_areas(len(lats),len(lons))
	# Do weighting
	field_av=np.zeros(len(tims))
	for i in range(len(tims)):
	        field_av[i]=np.average(field[i,:,:],weights=area)
	nyrs=len(tims)/12
	print nyrs
	yrs=np.arange(nyrs)
	field_ann=np.zeros([nyrs])
	for i in range(0,nyrs):
	        field_ann[i]=np.average(field_av[12*i:12*i+12])
	return field_ann
def monthly_to_xannual(field,lats,lons,tims,axix_no,x):
        area=gbox_areas(len(lats),len(lons))
        # Do weighting
        field_av=np.zeros(len(tims))
        for i in range(len(tims)):
                field_av[i]=np.average(field[i,:,:],weights=area)
        nyrs=len(tims)/(12*x)
        print nyrs
        yrs=np.arange(nyrs)
        field_ann=np.zeros([nyrs])
        for i in range(0,nyrs):
                field_ann[i]=np.average(field_av[12*x*i:12*x*i+12*x])
                field_ann[i]=np.average(field_av[12*x*i:12*x*i+12*x])
        return field_ann
def EXTR(loc,field,jobid,var):
	file=Dataset(loc)
	field[jobid+var]=file.variables[var][:].squeeze()
	lats=file.variables['latitude'][:]
	lons=file.variables['longitude'][:]
	hts=file.variables['atmosphere_hybrid_height_coordinate'][:]
	hts=hts*87
	return field,lats,lons,hts
