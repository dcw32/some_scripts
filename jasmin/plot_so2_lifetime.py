from netCDF4 import Dataset
import numpy as np
import calendar as calendar
from gridbox_areas import gbox_areas
from plot_1D import *
from scipy.interpolate import CubicSpline
import time
######################### Input ####################################
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
jobid='xnofh'
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
save=True
format='png'
extend='both'
####################################################################
#
# Extract
def extract(jobid,col,lab):
	store=np.zeros([len(jobid),361])
	for i in range(len(jobid)):
		file=Dataset(loc+jobid[i]+'/netcdf/'+jobid[i]+'_pd.nc')
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
		mass=mass/np.max(mass)
		y=np.where(mass==np.max(mass))[0][0]
		time=time-time[y]
		plt.plot(time,mass,color=col,linewidth=2,linestyle=':')
		cs=CubicSpline(time,mass)
		xs=np.arange(0,361,1)
		store[i,:]=cs(xs)
		#plt.plot(xs,cs(xs),color=col,linewidth=2)
	store=np.mean(store,axis=0)
	plt.plot(np.arange(0,361,1),store,color=col,linewidth=2,label=lab)
	return mass
#mass=extract(['xnofa','xnofc','xnofd','xnofe','xnoff'],'#E87722','HI-HAL')
mass=extract(['xnofa'],'#E87722','HI-HAL')
#mass=extract(['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl'],'#4E5B31','LO-HAL')
mass=extract(['xnofg'],'#4E5B31','LO-HAL')
#mass=extract(['xnofn','xnofo','xnofp','xnofq','xnofr','xnofs'],'#D50032','HI-SO2')
mass=extract(['xnofn'],'#D50032','HI-SO2')
#mass=extract(['xnoft','xnofu','xnofv','xnofw','xnofx','xnofy'],'#EF3340','LO-SO2')
mass=extract(['xnofu'],'#EF3340','LO-SO2')
#mass=extract(['xnofc'],'#D50032','HI-SO2')
plt.axhline(0.36788,c='k')
plt.legend(fancybox=True)
plt.xlim([0,360])
plt.xlabel('Days After Eruption')
plt.ylabel(r'SO$\mathregular{_{2}}$ Burden / Tg')
if save==True:
                print "SAVING "+outdir+outname+'.'+format
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight",dpi=200)
plt.show()
