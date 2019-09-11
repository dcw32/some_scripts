from netCDF4 import Dataset
import numpy as np
import pylab as plt
import time
import os
from gridbox_areas import *
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
jobid_base='xnfiy'
jobid_base='xnwdg'
#jobids=['xnofa','xnofb','xnofc','xnofd','xnofe','xnoff']
#jobids=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl','xnofa','xnofb','xnofc','xnofd','xnofe','xnoff','xnofn','xnofo','xnofp','xnofq','xnofr','xnofs','xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
jobids=['xnsvt','xnsvw','xnsvx','xnsvy']
jobids=['xnyxt','xnyxu','xnyxv','xnyxw','xnyxx','xnyxy']
#jobids=['xnyxa','xnyxb','xnyxc','xnyxd','xnyxe','xnyxf']
nnx=48
#nnx=120
#jobids=['xnofn','xnofo','xnofp','xnofq','xnofr','xnofs']
#jobids=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl']
#jobids=['xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
#jobids=['xnofa','xnofb','xnofc','xnofd','xnofe','xnoff','xnofg','xnofh','xnofi','xnofj','xnofk','xnofl','xnofn','xnofo','xnofp','xnofq','xnofr','xnofs','xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
#jobids=['xnofn','xnofo','xnofp','xnofq','xnofr','xnofs','xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
njobs=len(jobids)
label='HI-SO2+HAL'
var='temp_1'
var='temp'
#levs=[-50.,-30.,-20,-10.,10.,20.,30.,50.]
outdir='/home/users/dcw32/figures/thesis/pinatubo/'
outname=jobids[0]+'_nh_satdelta_timeseries'+time.strftime("%Y%m%d")
save=True
format='png'
################################################################
# BEGIN
################################################################
lsm_file=Dataset(loc+'pi_gem.nc')
lsm=lsm_file.variables['lsm'][:].squeeze()
baseline_file=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_pm.nc')
baseline_var=baseline_file.variables[var][:].squeeze()
baseline_t=baseline_file.variables['t'][:]
lons=baseline_file.variables['longitude'][:]
lats=baseline_file.variables['latitude'][:]
areas=gbox_areas(len(lats),len(lons))
#NH AVG
base_sat=np.zeros(len(baseline_t))
store=np.zeros([nnx,njobs])
for i in range(len(lats)):
	if lats[i]<-90.:
		#if lats[i]>75.:
			areas[i,:]=0.0
	#elif lats[i]==0.:
#		areas[i,:]=0.5*areas[i,:]
#areas=areas*lsm
for l in range(len(baseline_t)):
	base_sat[l]=np.average(baseline_var[l,:,:],weights=areas)
fig=plt.figure()
#For each jobid
for i in range(njobs):
	#Extract the data
	print jobids[i]
	anomaly_file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_pm.nc')
	anomaly_var=anomaly_file.variables[var][:].squeeze()
	anomaly_t=anomaly_file.variables['t'][:]
	nts=len(anomaly_t)
	#NH AVG
	anom_sat=np.zeros(len(anomaly_t))
	for l in range(len(anomaly_t)):
		anom_sat[l]=np.average(anomaly_var[l,:,:],weights=areas)
	#Find the location of the corresponding initial time in the base array
	base_loc=np.where(baseline_t==anomaly_t[0])[0][0]
	#Subtract the anomaly from the corresponding locations in the base
	delta=anom_sat-base_sat[base_loc:base_loc+nts]
	#Plot time in months
	t_plot=(anomaly_t-anomaly_t[0])/30.
	store[:,i]=delta[:nnx]
	plt.plot(t_plot,delta,c='lightgrey')
#store=np.mean(store,axis=1)
#for l in range(120-12):
#	if np.mean(store[l:l+12])>0.:
#		xx=l
#		print str(l)+" MONTH OF FIRST RECOVERY"
#		break
guillet_anoms=np.zeros([njobs,10])
for jobs in range(njobs):
	for years in range(10):
		guillet_anoms[jobs,years]=np.mean(store[12+12*years:15+12*years,jobs])
yrs_out=np.arange(1258,1258+10)
gu_mean=np.average(guillet_anoms,axis=0)
gu_uppe=gu_mean+np.std(guillet_anoms,axis=0)
gu_lowe=gu_mean-np.std(guillet_anoms,axis=0)
#print gu_uppe
#print gu_lowe
#store=np.mean(store,axis=1)
#print str(np.mean(np.mean(store[12:15,:],axis=1),axis=0))+" First Year NHLAND 40-90N dSAT ANOMALY"
print str(np.mean(np.mean(store[7:19,:],axis=1),axis=0))+" First Year GMST"
print str(np.std(np.mean(store[7:19,:],axis=1),axis=0))+" First Year GMST STDEV"
print str(np.min(np.mean(store[7:19,:],axis=1),axis=0))+" First Year GMST Lowest"
print str(np.max(np.mean(store[7:19,:],axis=1),axis=0))+" First Year GMST Highest"
print str(np.mean(np.mean(store[19:31,:],axis=1),axis=0))+" Second Year GMST ANOMALY"
print str(np.std(np.mean(store[19:31,:],axis=1),axis=0))+" Second Year GMST STDEV"
print str(np.min(np.mean(store[19:31,:],axis=1),axis=0))+" Second Year GMST Lowest"
print str(np.max(np.mean(store[19:31,:],axis=1),axis=0))+" Second Year GMST Highest"
store=np.mean(store,axis=1)
#print str(np.average(store[23+12:27+12]))+" Third Year NHLAND 40-90N dSAT ANOMALY"
#print str(np.average(store[23+24:27+24]))+" Fourth Year NHLAND 40-90N dSAT ANOMALY"
plt.plot(t_plot,store,c='k',linewidth=2.)
plt.xlim([0,nnx])
plt.ylim([-2.0,1.2])
if 'xx' in locals():
	plt.axvline(xx,c='r',linestyle='--')
plt.axhline(0.,c='k')
plt.title("SAT DELTA "+label)
if save==True:
 if not os.path.exists(outdir):
  os.makedirs(outdir)
 plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
 print "SAVING "+outdir+outname+'.'+format
plt.show()
