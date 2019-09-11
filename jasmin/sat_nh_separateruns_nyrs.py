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
#jobid_base='xnfiy'
jobid_base='xnwdg'
jobids=['xnofa','xnofb','xnofc','xnofd','xnofe','xnoff','xnofg','xnofh','xnofi','xnofj','xnofk','xnofl','xnofn','xnofo','xnofp','xnofq','xnofr','xnofs','xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
jobids=['xnofn','xnofo','xnofp','xnofq','xnofr','xnofs','xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
jobids=['xnofa','xnofb','xnofc','xnofd','xnofe','xnoff']
jobids=['xnywa','xnywb','xnywc','xnywd','xnywf','xnywg','xnywh','xnywi','xnywj','xnywk','xnywl','xnywn','xnywo','xnywp','xnywq','xnywr','xnyws','xnywt','xnywu','xnywv','xnyww','xnywx','xnywy']
jobids=['xnywn','xnywo','xnywp','xnywq','xnywr','xnyws','xnywt','xnywu','xnywv','xnyww','xnywx','xnywy']
jobids=['xnywa','xnywb','xnywc','xnywd','xnywf']
#jobids=['xnywa','xnywg','xnywn','xnywt','xnywb','xnywh','xnywo','xnywu','xnywc','xnywi','xnywp','xnywv','xnywd','xnywj','xnywq','xnyww','xnywe','xnywk','xnywr','xnywx','xnywf','xnywl','xnyws','xnywy']
var='temp' # Surface Air Temperature
njobs=len(jobids)
#Indices for seasons
nplots=24
ind_1=[6,12,18,24]
ind_2=[9,15,21,27]
index=1
(guilletmean,guilletmin,guilletmax)=(-0.7145,-1.3906,-0.0867)
#(guilletmean,guilletmin,guilletmax)=(-1.1900,-1.7372,-0.5608)
llab='1258 JJA '
#llab='1259 JJA '
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
#baseline_file=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_atmos_long.nc')
baseline_file=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_pm.nc')
baseline_var=baseline_file.variables[var][:600,:,:].squeeze()
baseline_t=baseline_file.variables['t'][:600]
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
store=np.zeros([120,njobs,len(lats),len(lons)])
store_base=np.zeros([120,njobs,len(lats),len(lons)])
baseline_yr=len(baseline_t)/12
base_clim=np.zeros([12,len(lats),len(lons)])
base_sat_reshape=np.zeros([baseline_yr,12,len(lats),len(lons)])
for lat in range(len(lats)):
	for lon in range(len(lons)):
		base_sat_reshape[:,:,lat,lon]=np.reshape(baseline_var[:,lat,lon],[baseline_yr,12])
for month in range(12):
	base_clim[month,:,:]=np.mean(base_sat_reshape[:,month,:,:],axis=0)
#sys.exit()
store_anom=np.zeros([120,njobs,len(lats),len(lons)])
#fig=plt.figure()
#For each jobid
for i in range(njobs):
	#Extract the data
	print jobids[i]
	#anomaly_file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_atmos.nc')
	anomaly_file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_atmos.nc')
	anomaly_var=anomaly_file.variables[var][:120,:,:].squeeze()
	anomaly_t=anomaly_file.variables['t'][:120]
	nts=len(anomaly_t)
	#Find the location of the corresponding initial time in the base array
	if METHOD==1:
		base_loc=np.where(baseline_t==anomaly_t[0])[0][0]
		print base_loc
		print base_loc+120
		#Subtract the anomaly from the corresponding locations in the base
		store_base[:,i,:,:]=baseline_var[base_loc:base_loc+120,:,:]
		store_anom[:,i,:,:]=anomaly_var[:120,:,:]
		delta=anomaly_var-baseline_var[base_loc:base_loc+nts,:,:]
	elif METHOD==2:
		delta=np.zeros([120,len(lats),len(lons)])
		for l in range(120):
			val=int(((anomaly_t[l]-baseline_t[0])/30.)%12)
			delta[l,:,:]=anomaly_var[l,:,:]-base_clim[val,:,:]
	else:
		sys.exit("METHOD NOT RECOGNISED")
	#Plot time in months
	t_plot=(anomaly_t-anomaly_t[0])/30.
	store[:,i,:,:]=delta[:120,:,:]
#enddo(njobs)
nyrs=8
minims=np.zeros(nyrs)
maxims=np.zeros(nyrs)
for i in range(nyrs):
	B=np.mean(store[12+12*i:15+12*i,:,:,:],axis=0)
	A=np.zeros(len(jobids))
	for j in range(len(jobids)):
		A[j]=np.average(B[j,:,:],weights=areas)
	minims[i]=np.min(A)
	maxims[i]=np.max(A)
print minims
print maxims
