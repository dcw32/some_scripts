from netCDF4 import Dataset
from matplotlib import gridspec
import numpy as np
import calendar as calendar
from gridbox_areas import gbox_areas
from plot_1D import *
from calc_solzen import solar_zen
import time
######################### Input ####################################
loc='/shared/netscratch/dcw32/um/'
#loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
jobids1=['xnofa','xnofb','xnofc','xnofd','xnofe','xnoff']
jobids2=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl']
jobids3=['xnofn','xnofo','xnofp','xnofq','xnofr','xnofs','xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
jobid_base='xnfiy'
#jobids=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl']
var='salinity'
varn='samalas_uvi'
print var
# Plotting Arguments
levs=[0.0,2.5,5.5,7.5,10.5,100.5]
#levs=[-1600,-400,-100,-50,-20,0,20,50,100,400,1600]
cols=['#389e29','#fff200','#f08200','#e32d12','#ad5e9b']
#colmap='RdGy_r'
#cmap=plt.cm.get_cmap(colmap)
#levs2=np.insert(levs,0,levs[0]-1)
#levs2=np.append(levs2,levs2[len(levs2)-1]+1)
#norm=mpl.colors.BoundaryNorm(levs2, ncolors=cmap.N, clip=True)
#cbar_label='Ozone Column / DU'
#levs=[0,4,8,12,16,20,24,28]
outdir='./../figures/'
outname='fig4'
save=True
format='pdf'
extend='neither'
mm1=240
mm2=120
mm3=120
#mm=240
####################################################################
#
# Extract
#
file=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_oz.nc')
file2=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_atmos.nc')
lats=file.variables['latitude'][:]
lons=file.variables['latitude'][:]
time=file.variables['t'][:]
ozcol=file.variables[var][:,0,:,:]
clrsk=file2.variables['field208'][:].squeeze()
allsk=file2.variables['field203'][:].squeeze()
CRF=allsk/clrsk
#ozcol=ozcol[:,0,:,:]
yrbas=time/12
clima=np.zeros([12,len(lats),len(lons)])
for mon in range(12):
	clima[mon,:,:]=np.mean(ozcol[mon::12])
solzen=np.zeros([360,len(lats),len(lons)])
for day in range(360):
	for lat in range(len(lats)):
		for lon in range(len(lons)):
			solzen[day,lat,lon]=solar_zen(day,lons[lon],lats[lat])
solzen_clima=np.zeros([12,len(lats),len(lons)])
for mon in range(12):
        solzen_clima[mon,:,:]=np.max(solzen[mon*30:mon*30+30,:,:],axis=0)
clima=np.mean(clima,axis=2)
clima_proc=(clima/300.)**-1.23
solzen_clima=np.mean(solzen_clima,axis=2)
clima_proc=np.roll(clima_proc,5,axis=0)
solzen_clima=np.roll(solzen_clima,5,axis=0)
uvi_clim=4.2*solzen_clima*clima_proc
#plt.contourf(np.arange(12),lats,np.transpose(uvi_clim))
#plt.colorbar()
#plt.show()
#
fig=plt.figure()
gs1 = gridspec.GridSpec(2, 2)
gs1.update(wspace=0.0)
ax1 = plt.subplot(gs1[0,:])
ax2 = plt.subplot(gs1[1,0])
ax3 = plt.subplot(gs1[1,1])
def extract_col(axis,jobids,mm,l_right,clima_proc,solzen_clima):
	njob=len(jobids)
	store1=np.zeros([njob,mm,73])
	clima_now=np.tile(clima_proc,(mm/12,1))
	solz_now=np.tile(solzen_clima,(mm/12,1))
	for i in range(njob):
		file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_oz.nc')
		file2=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_atmos.nc')
		lats=file.variables['latitude'][:]
		lons=file.variables['longitude'][:]
		time=file.variables['t'][:]
		field=file.variables[var][:].squeeze()
		field=field[:,0,:,:]
		clrsk=file2.variables['field208'][:].squeeze()
		allsk=file2.variables['field203'][:].squeeze()
		CRF=allsk/clrsk
		file.close()
		file2.close()
		#TimexLat
		field=np.mean(field,axis=2)
		field=(field/300.)**-1.23
		store1[i,:,:]=field
	time=time/360.
	time=time-time[0]
	store=np.mean(store1,axis=0)
	#UVI index
	store=4.2*solz_now*store
	#store=100*(store-clima_now)/clima_now
	# Start Plotting
	cf1=axis.contourf(time,lats,np.transpose(store),levs,colors=cols)
	#axis.contour(time,lats,np.transpose(store),[220],colors='k',linewidths=4)
        axis.set_yticks([-90,-45,0,45,90])
        axis.set_yticklabels(['90S','45S','EQ','45N','90N'])
	if l_right==True:
                axis.yaxis.tick_right()
                axis.yaxis.set_ticks_position('both')
	return cf1
cf1=extract_col(ax1,['xnofa'],mm1,False,clima_proc,solzen_clima)
ax1.text(0.02,1.02,'A: HI-HAL',transform=ax1.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='14',color='k')
cf2=extract_col(ax2,['xnofg'],mm2,False,clima_proc,solzen_clima)
ax2.text(0.04,1.02,'B: LO-HAL',transform=ax2.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='14',color='k')
cf3=extract_col(ax3,['xnofn'],mm3,True,clima_proc,solzen_clima)
ax3.text(0.04,1.02,'C: BOTH-SO2',transform=ax3.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='14',color='k')
cbaxes=fig.add_axes([-0.06, 0.20, 0.02, 0.6])
cbar=plt.colorbar(cf3,cax=cbaxes,orientation="vertical",ticks=[1.25,3.5,6.5,9.,50.5])
cbar.ax.set_yticklabels(['Low','Medium','High','Very High','Extreme'])
cbar.set_label('UV Index')
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
