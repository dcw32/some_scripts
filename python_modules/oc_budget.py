from netCDF4 import Dataset
import numpy as N
def ocbu(jobid):
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_cond_secorg_oc.nc')
	secf1=file.variables['cond_oc_nuc_sol'][:]
        secf2=file.variables['cond_oc_ait_sol'][:]
        secf3=file.variables['cond_oc_acc_sol'][:]
        secf4=file.variables['cond_oc_coa_sol'][:]
        secf5=file.variables['cond_oc_ait_ins'][:]
	file.close()
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_prim_oc.nc')
	prim1=file.variables['prim_oc_to_ait_sol'][:]
	prim2=file.variables['prim_oc_to_ait_ins'][:]
	file.close()	
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_ddep_oc.nc')
	ddep1=file.variables['ddep_oc_ait_ins'][:]
	ddep2=file.variables['ddep_oc_ait_sol'][:]
	ddep3=file.variables['ddep_oc_acc_sol'][:]
	ddep4=file.variables['ddep_oc_ait_sol'][:]
	ddep5=file.variables['ddep_oc_nuc_sol'][:]
	file.close()
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_scav_imp_oc.nc')
	scavi1=file.variables['scav_imp_oc_ait_ins'][:]
	scavi2=file.variables['scav_imp_oc_coa_sol'][:]
	scavi3=file.variables['scav_imp_oc_acc_sol'][:]
	scavi4=file.variables['scav_imp_oc_ait_sol'][:]
	scavi5=file.variables['scav_imp_oc_nuc_sol'][:]
        file.close()
	file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_scav_nuc_oc.nc')
	scavn1=file.variables['scav_nuc_oc_ait_ins'][:]
	scavn2=file.variables['scav_nuc_oc_coa_sol'][:]
	scavn3=file.variables['scav_nuc_oc_acc_sol'][:]
	scavn4=file.variables['scav_nuc_oc_ait_sol'][:]
	scavn5=file.variables['scav_nuc_oc_nuc_sol'][:]
	file.close()
	secf=secf1+secf2+secf3+secf4+secf5
	prim=prim1+prim2
	ddep=ddep1+ddep2+ddep3+ddep4+ddep5
	scavi=scavi1+scavi2+scavi3+scavi4+scavi5
	scavn=scavn1+scavn2+scavn3+scavn4+scavn5
	return (secf, prim, ddep, scavi, scavn)
	return secf
