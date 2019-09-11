from netCDF4 import Dataset
import numpy as np
import calendar as calendar
from gridbox_areas import gbox_areas
from plot_1D import *
import time
import scipy.ndimage as ndimage
######################### Input ####################################
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
jobid=['xntpl','xntpm','xntpn','xntpo','xntpp']
cols=['#af95a6','#E89CAE','#0072CE','#BE4D00','#3f2a56']
labell=['0.5x','1x','2x','5x','10x']
#jobid2='xnoce'
#jobid='xntpi'
#jobid2='xnoce'
#jobid='xnocf'
print jobid
var='field1878_6'
var2='pad_2'
print var
# Plotting Arguments
levs=[50,100,150,200,250,300,350,400,450,500]
colmap='Spectral'
cbar_label='Ozone Column / DU'
outdir='/home/users/dcw32/figures/thesis/snowball/volc/'
outname=jobid[0]+'_reff_'+time.strftime("%Y%m%d")
save=True
format='png'
extend='both'
####################################################################
#
# Extract
#file=Dataset(loc+jobid+'/netcdf/'+jobid+'_atmos.nc')
for job in range(len(jobid)):
	file=Dataset(loc+jobid[job]+'/netcdf/'+jobid[job]+'_pm.nc')
	lats=file.variables['latitude'][:]
	lons=file.variables['longitude'][:]
	time=file.variables['t'][:]
	field1_1=file.variables[var][:].squeeze()
	field2_1=file.variables[var2][:].squeeze()
	#
	reff1=np.zeros([len(time)])
	for i in range(len(time)):
		reff1[i]=np.average(field1_1[i,:,:,:],weights=field2_1[i,:,:,:])
	reff1=reff1*1.327/2.0
	reff1=reff1*1.0E6
	time=(time-time[0])/360.
	print np.max(reff1)
	plt.plot(time,reff1,color=cols[job],linewidth=2,alpha=0.2)
        plt.plot(time,ndimage.filters.gaussian_filter(reff1,12.),color=cols[job],linewidth=2,label=labell[job])
plt.ylim(0.21,1.5)
plt.legend(loc="upper left",ncol=5,fancybox=True,columnspacing=1,handletextpad=0)
plt.xticks([0,5,10,15])
plt.ylabel(u'R$\mathregular{_{eff}}$ / $\mathregular{\u00B5}$m')
plt.xlabel('Years After Eruptions Begin')
if save==True:
                print "SAVING "+outdir+outname+'.'+format
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
plt.show()
