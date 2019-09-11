from netCDF4 import Dataset
import numpy.ma as ma
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
jobid_base='xnwdg'
jobids=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl','xnofa','xnofb','xnofc','xnofd','xnofe','xnoff','xnofn','xnofo','xnofp','xnofq','xnofr','xnofs','xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
jobids=['xnywg','xnywh','xnywi','xnywj','xnywk','xnywl','xnywa','xnywb','xnywc','xnywd','xnywe','xnywf','xnywn','xnywo','xnywp','xnywq','xnywr','xnyws','xnywt','xnywu','xnywv','xnyww','xnywx','xnywy']
var='temp' # Surface Air Temperature
njobs=len(jobids)
#Indices for seasons
nplots=2
ind_1=[12,24]
ind_2=[15,27]
stitle='Whole Ensemble'
conf=0.2
#Plotting arguments
labs=['1258 JJA','1259 JJA']
cols=('#673387','#923293','#0083c8','#60cbf4','#bcdef1','#d3e9ca','#f9e7bd','#f6b123','#d9372a','#ac2e2c')
levs=[-2.0,-1.2,-0.8,-0.4,-0.2,0.0,0.3,0.6,1.0]
outdir='/home/users/dcw32/figures/thesis/samalas/'
outname='ALLRUNS_nhreg_satdelta_guillet_'+time.strftime("%Y%m%d")
save=True
format='pdf'
################################################################
# BEGIN
################################################################
#Extract baseline
#baseline_file=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_atmos.nc')
baseline_file=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_pm.nc')
baseline_var=baseline_file.variables[var][:].squeeze()
baseline_t=baseline_file.variables['t'][:]
lons=baseline_file.variables['longitude'][:]
lats=baseline_file.variables['latitude'][:]
areas=gbox_areas(len(lats),len(lons))
store=np.zeros([48,njobs,len(lats),len(lons)])
store_base=np.zeros([48,njobs,len(lats),len(lons)])
store_anom=np.zeros([48,njobs,len(lats),len(lons)])
#
#tree_file=Dataset('/group_workspaces/jasmin2/ukca/vol2/dcw32/Obs/samalas/guillet_N48.nc')
tree_file=Dataset('/group_workspaces/jasmin2/ukca/vol2/dcw32/Obs/samalas/guillet2017nhtempgrid.nc')
tree_anomy=tree_file.variables['JJA_temperature_anomaly'][:].squeeze()
tree_lat=tree_file.variables['latitude'][:].squeeze()
tree_lon=tree_file.variables['longitude'][:].squeeze()
tree_anomy=tree_anomy[1:3,:,:]
ce=tree_file.variables['Calibration_validation_statistics'][2,:,:]
#
ncloc='/group_workspaces/jasmin2/ukca/vol2/dcw32/CMIP5/lmil_files/'
models=['BCC','GISS1','GISS4','GISS7','MPI','IPSL','MIROC','CCSM4']
sdate=[850,1201,1201,1201,850,850,850]
nmodels=len(models)
year=1258-1240
fig=plt.figure(figsize=(6,8.5))
mmm=np.zeros([nmodels,180,360])
#CMIP5
for mod in range(nmodels):
        #ax=fig.add_subplot(gs[mod+2])
        print models[mod]
        file=Dataset(ncloc+models[mod]+'_regrid.nc')
        varf=file.variables['tas'][:].squeeze()
        lats_cm=file.variables['lat'][:]
        lons_cm=file.variables['lon'][:]
        time=file.variables['time'][:]
        time=(time-time[0])/360.+sdate[0]
        nyrs=len(time)/12
        jja=np.zeros([nyrs,len(lats_cm),len(lons_cm)])
        for yr in range(nyrs):
                jja[yr,:,:]=np.mean(varf[12*yr+5:12*yr+8,:,:],axis=0)
        jja_avg=np.mean(jja[year-15:year+16,:,:],axis=0)
        jja_plot=jja[year,:,:]-jja_avg
        mmm[mod,:,:]=jja_plot
mmm=np.mean(mmm,axis=0)
mmm2=np.zeros([nmodels,180,360])
#CMIP5
for mod in range(nmodels):
        #ax=fig.add_subplot(gs[mod+2])
        print models[mod]
        file=Dataset(ncloc+models[mod]+'_regrid.nc')
        varf=file.variables['tas'][:].squeeze()
        lats_cm=file.variables['lat'][:]
        lons_cm=file.variables['lon'][:]
        time=file.variables['time'][:]
        time=(time-time[0])/360.+sdate[0]
        nyrs=len(time)/12
        jja=np.zeros([nyrs,len(lats_cm),len(lons_cm)])
        for yr in range(nyrs):
                jja[yr,:,:]=np.mean(varf[12*yr+5:12*yr+8,:,:],axis=0)
        jja_avg=np.mean(jja[year-14:year+17,:,:],axis=0)
        jja_plot=jja[year+1,:,:]-jja_avg
        mmm2[mod,:,:]=jja_plot
mmm2=np.mean(mmm2,axis=0)
#
for i in range(njobs):
	#Extract the data
	print jobids[i]
	#anomaly_file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_atmos.nc')
	anomaly_file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_pm.nc')
	anomaly_var=anomaly_file.variables[var][:].squeeze()
	anomaly_t=anomaly_file.variables['t'][:]
	nts=len(anomaly_t)
	#Find the location of the corresponding initial time in the base array
	base_loc=np.where(baseline_t==anomaly_t[0])[0][0]
	#Subtract the anomaly from the corresponding locations in the base
	store_base[:,i,:,:]=baseline_var[base_loc:base_loc+48,:,:]
	store_anom[:,i,:,:]=anomaly_var[:48,:,:]
	delta=anomaly_var-baseline_var[base_loc:base_loc+nts,:,:]
	#Plot time in months
	t_plot=(anomaly_t-anomaly_t[0])/30.
	store[:,i,:,:]=delta[:48,:,:]
#enddo(njobs)
gs=gridspec.GridSpec(3,2)
for i in range(nplots):
	ax=fig.add_subplot(gs[i])
	#im0 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')
	im0 = Basemap(projection='robin',resolution='l',lon_0=0)
        #im0 = Basemap(projection='npstere',boundinglat=35,lon_0=0,resolution='l')
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
	cay=im0.contourf(a,b,store_1,levs,colors=cols,extend='both')
	im0.contourf(a,b,store_statmask1,hatches=['...',None],cmap=plt.get_cmap('gray'),alpha=0.2,extend='both',rasterized=True)
	im0.drawcoastlines()
	#im0.drawparallels(np.arange(-90.,120.,30.))
	#im0.drawmeridians(np.arange(0.,420.,60.))
        im0.drawmapboundary()
	ax.set_title(labs[i],fontsize=12)
	ax2=fig.add_subplot(gs[i+2])	
	im0 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')
        #im0 = Basemap(projection='npstere',boundinglat=35,lon_0=0,resolution='l')
	store_1=ma.masked_array(tree_anomy[i,:,:],mask=[ce<0.1])
	#lonsi,latsi,store_1=latlon_shift(lons,lats,tree_anomy[i,:,:])
	lonsi,latsi,store_1=latlon_shift(tree_lon,tree_lat,store_1)
	c,d=im0(lonsi,latsi)
	cay=im0.contourf(c,d,store_1,levs,colors=cols,extend='both')
	#cay=im0.contourf(a,b,ma.masked_array(store_1,mask=[ce>0.1]),levs,colors=cols,extend='both')
        im0.drawcoastlines()
        #im0.drawparallels(np.arange(-90.,120.,30.))
        #im0.drawmeridians(np.arange(0.,420.,60.))
        im0.drawmapboundary()
	#CMIP
	ax3=fig.add_subplot(gs[i+4])
	im1 = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')
	if i==0:
		mmm_plot=mmm
	elif i==1:
		mmm_plot=mmm2
	print mmm_plot.shape
	lonsj,latsj,store_cmip=latlon_shift(lons_cm,lats_cm,mmm_plot)
	e,f=im1(lonsj,latsj)
	caz=im1.contourf(e,f,store_cmip,levs,colors=cols,extend='both')
	im1.drawcoastlines()
	im1.drawmapboundary()
#enddo(nplots)
#plt.figtext(0.15,0.75,'1258')
#plt.figtext(0.55,0.75,'1259')
cbaxes=fig.add_axes([0.2, -0.02, 0.6, 0.02])
cbar=plt.colorbar(cay,cax=cbaxes,orientation="horizontal")
cbar.set_label(r'Temperature Anomaly / $\degree$C')
plt.tight_layout()
if save==True:
 if not os.path.exists(outdir):
  os.makedirs(outdir)
 plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
 print "SAVING "+outdir+outname+'.'+format
plt.show()
