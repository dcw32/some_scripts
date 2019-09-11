from netCDF4 import Dataset
from matplotlib import gridspec
import numpy as np
import calendar as calendar
from gridbox_areas import gbox_areas
from plot_1D import *
import time
######################### Input ####################################
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
jobids1=['xnofa','xnofb','xnofc','xnofd','xnofe','xnoff']
jobids2=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl']
jobids3=['xnofn','xnofo','xnofp','xnofq','xnofr','xnofs','xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
#jobids=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl']
var='salinity'
varn='samalas_ozone_col_stdev'
print var
# Plotting Arguments
#levs=[50,100,150,200,250,300,350,400,450,500]
colmap='Reds'
#cbar_label='Ozone Column / DU'
cmap=plt.cm.get_cmap(colmap)
levs=[0,5,10,20,40,80]
#levs2=np.insert(levs,0,levs[0]-1)
levs2=np.append(levs,levs[len(levs)-1]+1)
norm=mpl.colors.BoundaryNorm(levs2, ncolors=cmap.N, clip=True)
outdir='/home/users/dcw32/figures/thesis_volcano/samalas/'
outname=varn+'_surftimese_plot'+time.strftime("%Y%m%d")
save=True
format='png'
extend='both'
mm1=240
mm2=120
mm3=120
#mm=240
####################################################################
#
# Extract
fig=plt.figure()
gs1 = gridspec.GridSpec(2, 2)
gs1.update(wspace=0.0)
ax1 = plt.subplot(gs1[0,:])
ax2 = plt.subplot(gs1[1,0])
ax3 = plt.subplot(gs1[1,1])
def extract_col(axis,jobids,mm,l_right):
	njob=len(jobids)
	store1=np.zeros([njob,mm,73])
	for i in range(njob):
		file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_oz.nc')
		lats=file.variables['latitude'][:]
		lons=file.variables['longitude'][:]
		time=file.variables['t'][:]
		field=file.variables[var][:].squeeze()
		field=field[:,0,:,:]
		#TimexLat
		field=np.mean(field,axis=2)
		store1[i,:,:]=field
	time=time/360.
	time=time-time[0]
	store=np.std(store1,axis=0)
	# Start Plotting
	cf1=axis.contourf(time,lats,np.transpose(store),levs,cmap='Reds',extend='max',norm=norm)
	#cf1=axis.contourf(time,lats,np.transpose(store),extend='max')
	axis.contour(time,lats,np.transpose(store),[220],colors='k',linewidths=4)
        axis.set_yticks([-90,-45,0,45,90])
        axis.set_yticklabels(['90S','45S','EQ','45N','90N'])
	if l_right==True:
                axis.yaxis.tick_right()
                axis.yaxis.set_ticks_position('both')
	return cf1
cf1=extract_col(ax1,jobids1,mm1,False)
ax1.text(0.02,1.02,'A',transform=ax1.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='14',color='k')
cf2=extract_col(ax2,jobids2,mm2,False)
ax2.text(0.04,1.02,'B',transform=ax2.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='14',color='k')
cf3=extract_col(ax3,jobids3,mm3,True)
ax3.text(0.04,1.02,'C',transform=ax3.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='14',color='k')
cbaxes=fig.add_axes([0.0, 0.20, 0.02, 0.6])
cbar=plt.colorbar(cf1,cax=cbaxes,orientation="vertical")
cbar.set_label('Ozone Column / Dobson Units')
cbaxes.yaxis.set_label_position('left')
ax3.text(0.00,-0.12,'Years After Eruption',transform=ax3.transAxes\
        ,verticalalignment='top',horizontalalignment='center'\
        ,fontsize='14',color='k')
if save==True:
                print "SAVING "+outdir+outname+'.'+format
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")

#plt.savefig(outdir+outname+'.'+format)
plt.show()
