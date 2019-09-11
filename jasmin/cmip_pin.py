from netCDF4 import Dataset
import numpy as np
from gridbox_areas import gbox_areas
import pylab as plt
import time
import sys
import csv
import pandas as pd
import os
#This plots the Northern Hemisphere SAT over the 1257 volcanic eruption
#
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/CMIP5/historical/'
models=['BCC','CCSM4','GISS','IPSL','MIROC','MPI']
MMM=np.zeros(len(models))
save=False
#
fig=plt.figure()
ax = plt.axes()
for model in range(len(models)):
	file=Dataset(loc+models[model]+'_regrid.nc')
	lons=file.variables['lon'][:]
	lats=file.variables['lat'][:]
	tas=file.variables['tas'][:]
	time=file.variables['time'][:]
	time=(time-time[0])/365.24+1850
	nyrs=len(time)/12
	yrs=np.zeros(nyrs)
	gmst=np.zeros(nyrs)
	areas=gbox_areas(len(lats),len(lons))
        annual_data=np.zeros([nyrs,len(lats),len(lons)])
	for i in range(nyrs):
        	annual_data[i,:,:]=np.mean(tas[12*i:12*i+12,:,:],axis=0)
        	yrs[i]=np.mean(time[12*i:12*i+12])
        	gmst[i]=np.average(annual_data[i,:,:],weights=areas)
	yrs=np.floor(yrs)
	gmst=gmst[np.logical_and(yrs>1976,yrs<2008)]
	yrs=yrs[np.logical_and(yrs>1976,yrs<2008)]
	gmst_fun=np.polyfit(yrs,gmst,3)
        gmst_detrend=gmst-np.polyval(gmst_fun,yrs)
	yrP=np.where(yrs==1992)[0][0]
	print models[model]
	print gmst_detrend[yrP]-np.mean(gmst_detrend[yrP-5:yrP+5])
	print gmst_detrend[yrP]-gmst_detrend[yrP-1]
	MMM[model]=gmst_detrend[yrP]-gmst_detrend[yrP-1]
	ax.plot(yrs,gmst_detrend)
	ax.set_xlim(1990,1995)
#if save==True:
#                print "SAVING "+outdir+outname+'.'+format
#                if not os.path.exists(outdir):
#                        os.makedirs(outdir)
#                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight",transparent=True)
print np.mean(MMM)
plt.show()
