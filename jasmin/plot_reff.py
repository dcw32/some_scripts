from netCDF4 import Dataset
import numpy as np
import calendar as calendar
from gridbox_areas import gbox_areas
from plot_1D import *
import time
######################### Input ####################################
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
jobid2='xnkls'
jobid='xnklq'
jobid2='xnofn'
jobid='xnsvt'
#jobid2='xnoce'
#jobid='xntpi'
#jobid2='xnoce'
#jobid='xnocf'
print jobid
var='field1878_2'
var2='pad_2'
print var
# Plotting Arguments
levs=[50,100,150,200,250,300,350,400,450,500]
colmap='Spectral'
cbar_label='Ozone Column / DU'
outdir='/home/users/dcw32/figures/agu/'
outname=jobid+'_'+jobid2+'_reff_'+time.strftime("%Y%m%d")
save=True
format='png'
extend='both'
####################################################################
#
# Extract
#file=Dataset(loc+jobid+'/netcdf/'+jobid+'_atmos.nc')
#file=Dataset(loc+jobid+'/netcdf/'+jobid+'_pm.nc')
file=Dataset(loc+jobid+'/netcdf/'+jobid+'_aer.nc')
file2=Dataset(loc+jobid2+'/netcdf/'+jobid2+'_aer.nc')
#file2=Dataset(loc+jobid2+'/netcdf/'+jobid2+'_pm.nc')
lats=file.variables['latitude'][:]
lons=file.variables['longitude'][:]
time=file.variables['t'][:]
time2=file2.variables['t'][:]
field1_1=file.variables[var][:].squeeze()#*1.327/2.0
field1_2=file2.variables[var][:].squeeze()
field2_1=file.variables[var2][:].squeeze()#*1.327/2.0
field2_2=file2.variables[var2][:].squeeze()
#
reff1=np.zeros([len(time)])
reff2=np.zeros([len(time2)])
for i in range(len(time)):
	reff1[i]=np.average(field1_1[i,:,:,:],weights=field2_1[i,:,:,:])
for j in range(len(time2)):
	reff2[j]=np.average(field1_2[j,:,:,:],weights=field2_2[j,:,:,:])
reff1=reff1*1.327/2.0
reff2=reff2*1.327/2.0
reff1=reff1*1.0E6
reff2=reff2*1.0E6
time2=(time2-time2[0])/30.
time=(time-time[0])/30.
print np.max(reff1)
print np.max(reff2)
plt.plot(time,reff1,color='#3f2a56',linewidth=3)
plt.plot(time2,reff2,color='#af95a6',linewidth=3)
plt.xlim([0,36])
plt.ylabel(u'R$\mathregular{_{eff}}$ / $\mathregular{\u00B5}$m')
plt.xlabel('Months After Eruption')
if save==True:
                print "SAVING "+outdir+outname+'.'+format
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
plt.show()
