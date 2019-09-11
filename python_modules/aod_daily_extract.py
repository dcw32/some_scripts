from netCDF4 import Dataset
import numpy as N
def aodk(jobid):
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_dailyacc_solaod.nc')
        aod_3=file.variables['aod_accum_sol'][:]
        file.close()
	return aod_3
