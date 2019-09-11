from netCDF4 import Dataset
import numpy as N
def secf(jobid):
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_cond_secorg_oc.nc')
	secf1=file.variables['cond_oc_nuc_sol'][:]
        secf2=file.variables['cond_oc_ait_sol'][:]
        secf3=file.variables['cond_oc_acc_sol'][:]
        secf4=file.variables['cond_oc_coa_sol'][:]
        secf5=file.variables['cond_oc_ait_ins'][:]
	secf=secf1+secf2+secf3+secf4+secf5
	return secf
