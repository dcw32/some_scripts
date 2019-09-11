from netCDF4 import Dataset
import numpy as N
def aodk(jobid):
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_coa_ins_aod.nc')
        aod_0=file.variables['aod_coarse_insol'][:]
        file.close()
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_coa_sol_aod.nc')
        aod_1=file.variables['aod_coarse_sol'][:]
        file.close()
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_acc_ins_aod.nc')
        aod_2=file.variables['aod_accum_insol'][:]
        file.close()
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_acc_sol_aod.nc')
        aod_3=file.variables['aod_accum_sol'][:]
        file.close()
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_ait_ins_aod.nc')
        aod_4=file.variables['unknown'][:]
        file.close()
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_ait_sol_aod.nc')
        aod_5=file.variables['aod_aitken_sol'][:]
        file.close()
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_dust_aod.nc')
        aod_6=file.variables['atmosphere_absorption_optical_thickness_due_to_dust_ambient_aerosol'][:]
        file.close()
	aodkclim=aod_0+aod_1+aod_2+aod_3+aod_4+aod_5+aod_6
	return aodkclim,aodkstdv
