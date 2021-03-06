from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np
import pylab as plt
import time
import os
from gridbox_areas import *
from plot_1D import *
from matplotlib import gridspec
from scipy.stats import ks_2samp
#Data locations & jobids
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
jobid_base='xnfiy'
#jobids=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl','xnofa','xnofb','xnofc','xnofd','xnofe','xnoff','xnofn','xnofo','xnofp','xnofq','xnofr','xnofs','xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
#jobids=['xnofn','xnofo','xnofp','xnofq','xnofr','xnofs','xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
jobids=['xnwdn','xnwdo','xnwdp','xnwdq','xnwdr','xnwds']
var='temp_1' # Surface Air Temperature
njobs=len(jobids)
#Indices for seasons
nplots=4
ind_1=[90,96,102,108]
ind_2=[93,99,105,111]
labs=['1816 DJF','1816 JJA','1817 DJF','1817 JJA']
stitle='Whole Ensemble'
conf=0.10
#Plotting arguments
#cols=('#673387','#923293','#0083c8','#60cbf4','#bcdef1','#d3e9ca','#f9e7bd','#f6b123','#d9372a','#ac2e2c')
levs=[-2.0,-1.2,-0.8,-0.4,-0.2,0.2,0.4,0.8,1.2,2.0]
#levs=[-50,-30,-10,-5,5,10,30,50]
outdir='/home/users/dcw32/figures/samalas/'
outname='ALLRUNS_nhreg_precipdelta_REG_WIN_'+time.strftime("%Y%m%d")
save=False
format='png'
################################################################
# BEGIN
################################################################
#Extract baseline
baseline_file=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_atmos.nc')
baseline_var=baseline_file.variables[var][:].squeeze()
baseline_t=baseline_file.variables['t'][:]
lons=baseline_file.variables['longitude'][:]
lats=baseline_file.variables['latitude'][:]
areas=gbox_areas(len(lats),len(lons))
store=np.zeros([120,njobs,len(lats),len(lons)])
store_base=np.zeros([120,njobs,len(lats),len(lons)])
store_anom=np.zeros([120,njobs,len(lats),len(lons)])
fig=plt.figure(figsize=(8,8))
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
	base_loc=np.where(baseline_t==anomaly_t[0])[0][0]
	#Subtract the anomaly from the corresponding locations in the base
	store_base[:,i,:,:]=baseline_var[base_loc:base_loc+120,:,:]
	store_anom[:,i,:,:]=anomaly_var[:120,:,:]
	#delta=100*(anomaly_var-baseline_var[base_loc:base_loc+nts,:,:])/baseline_var[base_loc:base_loc+nts,:,:]
	delta=(anomaly_var-baseline_var[base_loc:base_loc+nts,:,:])
	#delta=anomaly_var-baseline_var[base_loc:base_loc+nts,:,:]
	#delta=100*delta/baseline_var[base_loc:base_loc+nts,:,:]
	#Plot time in months
	t_plot=(anomaly_t-anomaly_t[0])/30.
	store[:,i,:,:]=delta[:120,:,:]
#enddo(njobs)
gs=gridspec.GridSpec(nplots/2,nplots/2)
for i in range(nplots):
	ax=fig.add_subplot(gs[i])
	#im0 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')
	im0=Basemap(projection='robin',lon_0=0)
	store_1=np.mean(store[ind_1[i]:ind_2[i],:,:,:],axis=0)
	store_base1=np.mean(store_base[ind_1[i]:ind_2[i],:,:,:],axis=0)
	store_anom1=np.mean(store_anom[ind_1[i]:ind_2[i],:,:,:],axis=0)
	store_statmask1=np.zeros([len(lats),len(lons)])
	# If the difference between the baseline and anomaly ensembles is statistically significant
	# at the 95% confidence level according to a K-S test don't stipple
	for a in range(len(lats)):
		for b in range(len(lons)):
			(stat1,pval1)=ks_2samp(store_base1[:,a,b],store_anom1[:,a,b])
			if pval1<conf:
				store_statmask1[a,b]=1.0
		#enddo(lons)
	#enddo(lats)
	#Now take the ensemble mean
	store_1=np.mean(store_1,axis=0)
	lonsi,latsi,store_1=latlon_shift(lons,lats,store_1)
	lonsi,latsi,store_statmask1=latlon_shift(lons,lats,store_statmask1)
	a,b=im0(lonsi,latsi)
	cay=im0.contourf(a,b,store_1,levs,cmap=plt.cm.get_cmap('RdBu_r'),extend='both')
	im0.contourf(a,b,store_statmask1,hatches=['...',None],cmap=plt.get_cmap('gray'),alpha=0.0,extend='both',rasterized=True)
	im0.drawcoastlines()
	im0.drawparallels(np.arange(-90.,120.,30.))
	im0.drawmeridians(np.arange(0.,420.,60.))
        im0.drawmapboundary()
	ax.set_title(labs[i],fontsize=12)
#enddo(nplots)
#plt.figtext(0.15,0.75,'1258')
#plt.figtext(0.55,0.75,'1259')
fig.suptitle(stitle)
#cbaxes=fig.add_axes([0.2, 0.05, 0.6, 0.03])
#cbar=plt.colorbar(cay,cax=cbaxes,orientation="horizontal")
#cbar=fig.add_colorbar(cay,cax=cbaxes,orientation="horizontal")
plt.tight_layout()
if save==True:
 if not os.path.exists(outdir):
  os.makedirs(outdir)
 plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
 print "SAVING "+outdir+outname+'.'+format
plt.show()
