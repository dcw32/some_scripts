from netCDF4 import Dataset
import numpy as N
def orw(jobid,nyrs):
        file=Dataset('/scratch/dcw32/netscratch/um/'+jobid+'/'+jobid+'_toa_fluxes.nc')
	or1=file.variables['outgoing_sw_toa_flux'][:]
        or2=file.variables['outgoing_lw_toa_flux'][:]
	or3=file.variables['clearsky_sw_toa_flux'][:]
	or4=file.variables['clearsky_lw_toa_flux'][:]
	allsky=or1+or2-or3-or4
	alls=N.empty([nyrs,145,192])
        for i in range (0,nyrs):
                q=12*i
                allt=allsky[q:q+12,:,:]
                alls[i,:,:]=N.mean(allt,axis=0)
        allstd=N.std(alls,axis=0)
	clearsky=or3+or4
        alls=N.empty([nyrs,145,192])
        for i in range (0,nyrs):
                q=12*i
                allt=clearsky[q:q+12,:,:]
                alls[i,:,:]=N.mean(allt,axis=0)
        clearstd=N.std(alls,axis=0)
	return (allsky, allstd, clearsky, clearstd)
