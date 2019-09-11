from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np
import pylab as plt
import time
import os
import sys
from gridbox_areas import *
from plot_1D import *
from matplotlib import gridspec
from scipy.stats import ks_2samp
#Data locations & jobids
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
jobid_base='xnwdg'
#jobids=['xnofa','xnofb','xnofc','xnofd','xnofe','xnoff','xnofg','xnofh','xnofi','xnofj','xnofk','xnofl','xnofn','xnofo','xnofp','xnofq','xnofr','xnofs','xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
jobids=['xnwgt','xnwgu','xnwgv','xnwgw','xnwgx','xnwgy']
jobids=['xnwgt','xnwgv','xnwgw','xnwgx']
var='temp_1' # Surface Air Temperature
njobs=len(jobids)
#Indices for seasons
nplots=4
ind_1=[6,12,18,24]
ind_2=[9,15,21,27]
index=3
(guilletmean,guilletmin,guilletmax)=(-0.7145,-1.3906,-0.0867)
#(guilletmean,guilletmin,guilletmax)=(-1.1900,-1.7372,-0.5608)
llab='1258 JJA '
outname='ALLRUNS_nhreg_satdelta_1258JJA_'+time.strftime("%Y%m%d")
labs=['HI-HAL','LO-HAL','HI-SO2','LO-SO2']
conf=0.2
#Plotting arguments
cols=('#673387','#923293','#0083c8','#60cbf4','#bcdef1','#d3e9ca','#f9e7bd','#f6b123','#d9372a','#ac2e2c')
levs=[-2.0,-1.2,-0.8,-0.4,-0.2,0.0,0.3,0.6,1.0]
outdir='/home/users/dcw32/figures/samalas/'
save=False
format='pdf'
METHOD=1
################################################################
# BEGIN
################################################################
#Extract baseline
baseline_file=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_atmos.nc')
baseline_var=baseline_file.variables[var][:,:,:].squeeze()
baseline_t=baseline_file.variables['t'][:]
lons=baseline_file.variables['longitude'][:]
lats=baseline_file.variables['latitude'][:]
lsmfile=Dataset(loc+'/pi_gem.nc')
lsm=lsmfile.variables['lsm'][:].squeeze()
lsmlats=lsmfile.variables['latitude'][:].squeeze()
lsmlons=lsmfile.variables['longitude'][:].squeeze()
areas=gbox_areas(len(lats),len(lons))
for lat in range(len(lats)):
	if lats[lat]<40.:
		areas[lat,:]=0.
areas=areas*lsm
store=np.zeros([27,njobs,len(lats),len(lons)])
store_base=np.zeros([27,njobs,len(lats),len(lons)])
baseline_yr=len(baseline_t)/12
base_clim=np.zeros([12,len(lats),len(lons)])
base_sat_reshape=np.zeros([baseline_yr,12,len(lats),len(lons)])
for lat in range(len(lats)):
	for lon in range(len(lons)):
		base_sat_reshape[:,:,lat,lon]=np.reshape(baseline_var[:,lat,lon],[baseline_yr,12])
for month in range(12):
	base_clim[month,:,:]=np.mean(base_sat_reshape[:,month,:,:],axis=0)
#sys.exit()
store_anom=np.zeros([27,njobs,len(lats),len(lons)])
fig=plt.figure(figsize=(8,11))
#fig=plt.figure()
#For each jobid
for i in range(njobs):
	#Extract the data
	print jobids[i]
	anomaly_file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_atmos.nc')
	anomaly_var=anomaly_file.variables[var][:].squeeze()
	anomaly_t=anomaly_file.variables['t'][:]
	nts=len(anomaly_t)
	#Find the location of the corresponding initial time in the base array
	if METHOD==1:
		base_loc=np.where(baseline_t==anomaly_t[0])[0][0]
		#Subtract the anomaly from the corresponding locations in the base
		store_base[:,i,:,:]=baseline_var[base_loc:base_loc+27,:,:]
		store_anom[:,i,:,:]=anomaly_var[:27,:,:]
		delta=anomaly_var-baseline_var[base_loc:base_loc+nts,:,:]
	elif METHOD==2:
		delta=np.zeros([27,len(lats),len(lons)])
		for l in range(27):
			val=int(((anomaly_t[l]-baseline_t[0])/30.)%12)
			delta[l,:,:]=anomaly_var[l,:,:]-base_clim[val,:,:]
	else:
		sys.exit("METHOD NOT RECOGNISED")
	#Plot time in months
	t_plot=(anomaly_t-anomaly_t[0])/30.
	store[:,i,:,:]=delta[:27,:,:]
#enddo(njobs)
gs=gridspec.GridSpec(6,4)
gs.update(wspace=0.025, hspace=0.05)
for i in range(nplots):
	ax=fig.add_subplot(gs[i])
	im0 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')
	store_1=np.mean(store[ind_1[index]:ind_2[index],:,:,:],axis=0)
	store_base1=np.mean(store_base[ind_1[index]:ind_2[index],:,:,:],axis=0)
	store_anom1=np.mean(store_anom[ind_1[index]:ind_2[index],:,:,:],axis=0)
	#Now take the ensemble mean
	store_1=store_1[i,:,:]
	nhbar=np.average(store_1,weights=areas)
	lonsi,latsi,store_1=latlon_shift(lons,lats,store_1)
	a,b=im0(lonsi,latsi)
	cay=im0.contourf(a,b,store_1,levs,colors=cols,extend='both')
	im0.drawcoastlines()
	im0.drawparallels(np.arange(-90.,120.,30.))
	im0.drawmeridians(np.arange(0.,420.,60.))
        im0.drawmapboundary()
        if nhbar<guilletmin:
            colpl='#D50032'
        elif nhbar>guilletmax:
            colpl='#D50032'
        else:
            colpl='#64A70B'
	if i<4:
		ax.set_title(labs[i],fontsize=12)
        ax.text(0.87,0.95,"{0:+.2f}".format(nhbar),transform=ax.transAxes\
             ,verticalalignment='center',horizontalalignment='left'\
             ,fontsize='12',color=colpl)
#enddo(nplots)
#plt.figtext(0.15,0.75,'1258')
#plt.figtext(0.55,0.75,'1259')
#fig.suptitle(stitle)
cbaxes=fig.add_axes([0.2, 0.065, 0.6, 0.015])
cbar=plt.colorbar(cay,cax=cbaxes,orientation="horizontal")
cbar.set_label(str(llab)+r'Temperature Anomaly / $\degree$C')
#cbar=fig.add_colorbar(cay,cax=cbaxes,orientation="horizontal")
#plt.tight_layout()
if save==True:
 if not os.path.exists(outdir):
  os.makedirs(outdir)
 plt.savefig(outdir+outname+'.'+format,dpi=200,bbox_inches='tight')
 print "SAVING "+outdir+outname+'.'+format
plt.show()
