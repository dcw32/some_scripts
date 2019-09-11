from netCDF4 import Dataset
import numpy as np
import calendar as calendar
from gridbox_areas import gbox_areas
from plot_1D import *
import time
######################### Input ####################################
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
#jobid='xnoce'
jobid='xnklp'
print jobid
var='field1878_2'
var2='pad_2'
print var
# Plotting Arguments
levs=[50,100,150,200,250,300,350,400,450,500]
colmap='Spectral'
cbar_label='Ozone Column / DU'
outdir='/home/users/dcw32/figures/agu/'
outname=jobid+'_reff_contour_'+time.strftime("%Y%m%d")
save=True
format='png'
extend='both'
####################################################################
#
# Extract

cols=['#663398','#330065','#333365','#330065','#010066','#003466','#3266cb','#339967','#006632','#336601','#669934','#ffff66','#cc6601','#ff3334','#fe3200','#cd3301','#cd3301','#993400']
levs=[0.0,0.0005,0.001,0.005,0.01,0.05,0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.5,3.0,5.0]
#cols=['midnightblue','royalblue','green','darkgreen','forestgreen','darkolivegreen','gold','saddlebrown','firebrick','indianred','tomato']
#file=Dataset(loc+jobid+'/netcdf/'+jobid+'_atmos.nc')
file=Dataset(loc+jobid+'/netcdf/'+jobid+'_pm.nc')
lats=file.variables['latitude'][:]
lons=file.variables['longitude'][:]
time=file.variables['t'][:]
field1_1=file.variables[var][:].squeeze()
field2_1=file.variables[var2][:].squeeze()
#
print field1_1.shape
reff1=np.zeros([len(time),len(lats)])
for i in range(len(time)):
 for j in range(len(lats)):
	#reff1[i,j]=np.average(field1_1[i,35,j,:],weights=field2_1[i,35,j,:])
	reff1[i,j]=np.average(field1_1[i,:,j,:],weights=field2_1[i,:,j,:])
reff1=reff1*1.327/2.0
reff1=reff1*1.0E6
time=(time-time[0])/30.
print np.max(reff1)
print reff1[len(time)-1,:]
plt.contourf(time,lats,np.transpose(reff1),levs,colors=cols)
cbar=plt.colorbar()
cbar.set_label(r'R$\mathregular{_{eff}}$ / $\mathregular{\mu}$m')
plt.xlim([0,36])
#plt.ylabel(r'R$_{eff}$ / $\mu$m')

plt.xlabel('Months After Eruption')
#plt.text(20,1.2,'Feedbacks On',color='r',fontsize=16)
#plt.text(20,1.30,'Feedbacks Off',color='k',fontsize=16)
if save==True:
                print "SAVING "+outdir+outname+'.'+format
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
plt.show()
