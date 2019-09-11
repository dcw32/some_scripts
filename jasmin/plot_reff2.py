from netCDF4 import Dataset
import numpy as np
import calendar as calendar
from gridbox_areas import gbox_areas
from plot_1D import *
import time
######################### Input ####################################
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
jobid='xnofn'
#jobid2='xnoce'
#jobid='xntpi'
#jobid2='xnoce'
#jobid='xnocf'
print jobid
var='field607'
var2='field608'
print var
# Plotting Arguments
levs=[50,100,150,200,250,300,350,400,450,500]
colmap='Spectral'
cbar_label='Ozone Column / DU'
outdir='/home/users/dcw32/figures/agu/'
save=False
format='png'
extend='both'
####################################################################
#
# Extract
#file=Dataset(loc+jobid+'/netcdf/'+jobid+'_atmos.nc')
#file=Dataset(loc+jobid+'/netcdf/'+jobid+'_pm.nc')
file=Dataset(loc+jobid+'/netcdf/'+jobid+'_aer.nc')
file2=Dataset(loc+jobid+'/netcdf/'+jobid+'_atmos.nc')
lats=file.variables['latitude'][:]
lons=file.variables['longitude'][:]
time=file.variables['t'][:]
nd_nucsol=file.variables[var][:,20:,:,:].squeeze()
nucsol_mmr=field2=file.variables[var2][:,20:,:,:].squeeze()
pressure=file2.variables['p_1'][:,20:,:,:]
theta=file2.variables['theta'][:,20:,:,:]
#nd_nucsol=file.variables[var][:,28:36,42:44,:].squeeze()
#nucsol_mmr=field2=file.variables[var2][:,28:36,42:44,:].squeeze()
#pressure=file2.variables['p_1'][:,28:36,42:44,:]
#theta=file2.variables['theta'][:,28:36,42:44,:]
temp=theta*(pressure/100000.)**(287.05/1005.46)
aird=pressure/(temp*1E6*1.38E-23)
#
sigma=1.4
h2so4_mdnd=nucsol_mmr*aird*0.02897/0.098
nd_nucsol=nd_nucsol*aird
h2so4_mdnd=np.sum(h2so4_mdnd,axis=3)
h2so4_mdnd=np.sum(h2so4_mdnd,axis=2)
h2so4_mdnd=np.sum(h2so4_mdnd,axis=1)
nd_nucsol=np.sum(nd_nucsol,axis=3)
nd_nucsol=np.sum(nd_nucsol,axis=2)
nd_nucsol=np.sum(nd_nucsol,axis=1)
#h2so4_md=np.zeros([len(time),60,73,96])
#for ti in range(len(time)):
#	for hi in range(60):
#		for la in range(73):
#			for lo in range(96):
#				if nd_nucsol[ti,hi,la,lo]>0.:
#					h2so4_md[ti,hi,la,lo]=h2so4_mdnd[ti,hi,la,lo]/nd_nucsol[ti,hi,la,lo]
h2so4_md=h2so4_mdnd/nd_nucsol
h2so4_v=h2so4_md*0.098/(6.022E23*1769.0)
x_nucsol=np.exp(4.5*np.log(sigma)**2)
rbar_nucsol=0.5*(6.0*h2so4_v/np.pi/x_nucsol)**(1./3.)
y_nucsol=np.pi*np.exp(2.0*np.log(sigma)**2)
y2_nucsol=np.pi*np.exp(4.5*np.log(sigma)**2)
saconc=4.0*nd_nucsol*rbar_nucsol**2*y_nucsol
volconc=(4.0/3.0)*nd_nucsol*rbar_nucsol**3*y2_nucsol
reff=volconc/saconc
#reff=np.nanmean(reff,axis=3)
#reff=np.nanmean(reff,axis=2)
#reff=np.nanmean(reff,axis=1)
reff2=rbar_nucsol*np.exp(2.5*np.log(sigma)**2)
#print np.sum(np.isnan(reff2))
#reff2=np.nanmean(reff2,axis=3)
#print np.max(reff2)
#reff2=np.nanmean(reff2,axis=2)
#reff2=np.nanmean(reff2,axis=1)
#print reff2
#
time=(time-time[0])/30.
#plt.plot(time,reff,color='#3f2a56',linewidth=3)
plt.plot(time,reff2,color='#3f2a56',linewidth=3,linestyle=':')
#CS=plt.contourf(time,lats,np.transpose(h2so4_mdnd),color='#3f2a56',linewidth=3,linestyle=':')
#levs=[0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6]
#CS=plt.contourf(time,lats,1.E6*np.transpose(reff2),levs,extend='both')
#plt.clabel(CS, inline=1, fontsize=10)
#plt.colorbar(CS)
#plt.xlim([0,36])
plt.ylabel(u'R$\mathregular{_{eff}}$ / $\mathregular{\u00B5}$m')
plt.xlabel('Months After Eruption')
plt.show()
