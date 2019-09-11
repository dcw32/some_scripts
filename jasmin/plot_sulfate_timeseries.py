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
jobids3=['xnofn','xnofo','xnofp','xnofq','xnofr','xnofs']
jobids4=['xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
#jobids=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl']
var='field608'
#var_1='temp_2' #AIT
#var_2='field614' #COA
var2='field643'
#
varn='samalas_sulfate'
print var
# Plotting Arguments
#levs=[50,100,150,200,250,300,350,400,450,500]
colmap='Spectral'
#cbar_label='Ozone Column / DU'
levs=[0.1,0.2,0.5,1.0,2.0,5.0]
cmap=plt.cm.get_cmap(colmap)
levs2=np.insert(levs,0,levs[0]-1)
levs2=np.append(levs2,levs2[len(levs2)-1]+1)
norm=mpl.colors.BoundaryNorm(levs2, ncolors=cmap.N, clip=True)
#levs=[0,60,120,180,240,300,360,420]
outdir='/home/users/dcw32/figures/thesis_volcano/samalas/'
outname=varn+'_surftimese_plot'+time.strftime("%Y%m%d")
save=True
format='png'
extend='both'
mm1=24
mm2=24
mm3=24
mm4=24
#mm=240
####################################################################
#
# Extract
fig=plt.figure()
gs1 = gridspec.GridSpec(2, 2)
gs1.update(wspace=0.0)
ax1 = plt.subplot(gs1[0,0])
ax2 = plt.subplot(gs1[0,1])
ax3 = plt.subplot(gs1[1,0])
ax4 = plt.subplot(gs1[1,1])
def extract_col(axis,jobids,mm,l_right):
	njob=len(jobids)
	store1=np.zeros([njob,mm,73])
	for i in range(njob):
		file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_aer.nc')
		file2=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_oz.nc')
		lats=file.variables['latitude'][:]
		lons=file.variables['longitude'][:]
		time=file.variables['t'][:]
		#field=file.variables[var][:].squeeze()+file.variables[var_1][:].squeeze()+file.variables[var_2][:].squeeze()
		field=file.variables[var][:mm,:,:,:].squeeze()
		file.close()
		field2=file2.variables[var2][:mm,:,:,:].squeeze()
		mass_boxes=field*field2
		mass=np.sum(mass_boxes,axis=3)
		mass=np.sum(mass,axis=1)
		#TimexLat
		#field=np.mean(field,axis=2)
		store1[i,:,:]=mass[:mm,:]*10**-9
	time=time/360.
	time=time-time[0]
	store=np.mean(store1,axis=0)
	# Start Plotting
	cf1=axis.contourf(time[:mm],lats,np.transpose(store),levs,cmap=cmap,norm=norm,extend='both')
	#axis.contour(time,lats,np.transpose(store),[220],colors='k',linewidths=4)
        axis.set_yticks([-90,-45,0,45,90])
        axis.set_yticklabels(['90S','45S','EQ','45N','90N'])
	if l_right==True:
                axis.yaxis.tick_right()
                axis.yaxis.set_ticks_position('both')
	return cf1
cf1=extract_col(ax1,jobids1,mm1,False)
ax1.text(0.02,1.02,'A: HI-HAL',transform=ax1.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='14',color='k')
cf2=extract_col(ax2,jobids2,mm2,True)
ax2.text(0.02,1.02,'B: LO-HAL',transform=ax2.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='14',color='k')
cf3=extract_col(ax3,jobids3,mm3,False)
ax3.text(0.02,1.02,'C: HI-SO2',transform=ax3.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='14',color='k')
cf4=extract_col(ax4,jobids4,mm4,True)
ax4.text(0.02,1.02,'D: LO-SO2',transform=ax4.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='14',color='k')
cbaxes=fig.add_axes([0.0, 0.20, 0.02, 0.6])
cbar=plt.colorbar(cf3,cax=cbaxes,orientation="vertical")
cbar.set_label('Sulfate Column / Tg')
cbaxes.yaxis.set_label_position('left')
ax4.text(0.00,-0.12,'Years After Eruption',transform=ax4.transAxes\
        ,verticalalignment='top',horizontalalignment='center'\
        ,fontsize='14',color='k')
if save==True:
                print "SAVING "+outdir+outname+'.'+format
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")

#plt.savefig(outdir+outname+'.'+format)
plt.show()
