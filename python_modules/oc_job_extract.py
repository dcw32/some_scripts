from netCDF4 import Dataset
import numpy as N
def oc_mmrk(jobid):
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_coa_sol_oc_mmr.nc')
        oc_mmr_1=file.variables['coa_sol_oc_mmr'][:]
        file.close()
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_acc_sol_oc_mmr.nc')
        oc_mmr_3=file.variables['acc_sol_oc_mmr'][:]
        file.close()
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_ait_insol_oc_mmr.nc')
        oc_mmr_4=file.variables['ait_insol_oc_mmr'][:]
        file.close()
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_ait_sol_oc_mmr.nc')
        oc_mmr_5=file.variables['ait_sol_oc_mmr'][:]
        file.close()
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_nuc_sol_oc_mmr.nc')
        oc_mmr_6=file.variables['nuc_sol_oc_mmr'][:]
        file.close()
	oc_mmrk=oc_mmr_1+oc_mmr_3+oc_mmr_4+oc_mmr_5+oc_mmr_6
	return oc_mmrk
