from netCDF4 import Dataset
import numpy as N
def plant(jobid):
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_plant_prod.nc')
	secf1=file.variables['gross_primary_productivity'][:]
        secf2=file.variables['net_primary_productivity'][:]
        secf3=file.variables['plant_respiration'][:]
	return (secf1,secf2,secf3)
