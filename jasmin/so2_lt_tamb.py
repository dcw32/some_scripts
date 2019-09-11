from netCDF4 import Dataset
import numpy as np
import calendar as calendar
from gridbox_areas import gbox_areas
from plot_1D import *
import time
######################### Input ####################################
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
jobid='xnwdn'
print jobid
var='field571'
var2='field643'
print var
# Plotting Arguments
levs=[50,100,150,200,250,300,350,400,450,500]
colmap='Spectral'
cbar_label='Ozone Column / DU'
outdir='/home/users/dcw32/figures/thesis/samalas/'
outname=jobid+'_so2_lifetime_'+time.strftime("%Y%m%d")
save=False
format='png'
extend='both'
####################################################################
#
# Extract
def extract(jobid,col,lab):
	file=Dataset(loc+jobid+'/netcdf/'+jobid+'_pd.nc')
	lats=file.variables['latitude'][:]
	lons=file.variables['longitude'][:]
	time=file.variables['t'][:]
	field=file.variables[var][:].squeeze()
	field2=file.variables[var2][:].squeeze()
#
	file.close()
	mass_boxes=field*field2
	mass=np.sum(mass_boxes,axis=3)
	mass=np.sum(mass,axis=2)
	mass=np.sum(mass,axis=1)
	mass=mass*1E-9
	y=np.where(mass==np.max(mass))[0][0]
	time=time-time[y]
	plt.plot(time,mass,color=col,linewidth=2,label=lab)
	return mass
mass=extract('xnwdn','#D50032','Tambora')
plt.axhline(0.36788*np.max(mass),c='k')
plt.legend(fancybox=True)
#plt.xlim([0,360])
plt.xlabel('Days After Eruption')
plt.ylabel(r'SO$\mathregular{_{2}}$ Burden / Tg')
if save==True:
                print "SAVING "+outdir+outname+'.'+format
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight",dpi=200)
plt.show()
