from netCDF4 import Dataset
import numpy as np
import pylab as plt
from gridbox_areas import gbox_areas
import time
import os
from scipy import signal
#
from scipy.optimize import leastsq
################################################################
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
jobid='xnfiy'
var='temp'
xs=0.5+np.array([4.,46.,27.,17.,21.,23.])
nens=xs.shape[0]
outdir='/home/users/dcw32/figures/samalas/'
outname=jobid+'_enso_powerspectrum_'+time.strftime("%Y%m%d")
save=True
format='pdf'
file=Dataset(loc+jobid+'/netcdf/'+jobid+'_atmos.nc')
#Vars for ENSO 3.4 Index
dlon=3.75 #Lon box size
#lat_min=-5.0 #5S
#lat_max=5.0 #5N
lon_mini=190.0 #170W
lon_maxi=240.0 #120W
#Extract
lats=file.variables['latitude'][:]
lons=file.variables['longitude'][:]
time=file.variables['t'][:]
time=time/360.0
time=time-time[0]
sst=file.variables[var][:].squeeze()
################################################################
#Formulae
def second_largest(numbers):
    count = 0
    m1 = m2 = float('-inf')
    for x in numbers:
        count += 1
        if x > m2:
            if x >= m1:
                m1, m2 = x, m1            
            else:
                m2 = x
    return m2 if count >= 2 else None
##
#Gridbox weights
lon_maxs=lons+dlon/2.
lon_mins=lons-dlon/2.
enso_w=np.zeros([len(lats),len(lons)])
for i in range(len(lons)):
  if lon_maxs[i]<lon_maxi:
    if lon_mins[i]>lon_mini:
#5N-5S
      enso_w[34:38,i]=1.0
    elif lon_maxs[i]>lon_mini:
      enso_w[34:38,i]=(lon_maxs[i]-lon_mini)/dlon
  elif lon_mins[i]<lon_maxi:
    if lon_mins[i]>lon_mini:
      enso_w[34:38,i]=(lon_maxi-lon_mins[i])/dlon
#Now
area=gbox_areas(len(lats),len(lons))
area=area*enso_w
enso=np.zeros([len(time)])
enso_ra=np.zeros([len(time)])
for k in range(len(time)):
  enso[k]=np.average(sst[k,:,:],weights=area)
for k in range(len(time)):
  if k>1:
    enso_ra[k]=np.mean(enso[k-2:k+1])
  elif k>0:
    enso_ra[k]=np.mean(enso[k-1:k+1])
  else:
    enso_ra[k]=enso[k]
for mons in range(12):
        gm_mon=np.mean(enso_ra[mons::12])
        enso_ra[mons::12]=enso_ra[mons::12]-gm_mon
#
fig,ax=plt.subplots(figsize=(4,4))
f,Pxx_den=signal.periodogram(enso_ra[:480],12.)
# Autocorr lag1
#plt.plot(f, Pxx_den ,color='#EF3340',linewidth=2)
ax.set_xlabel(r'Frequency / yr$^{-1}$')
ax.set_ylabel('Spectral Power Density')
#plt.ylim([9,8.E7])
ax.set_xlim([0,1.0])
print str(round(1.0/(f[np.where(Pxx_den==np.max(Pxx_den))[0][0]]),2))
print str(round(1.0/(f[np.where(Pxx_den==second_largest(Pxx_den))[0][0]]),2))    
#Data
sst_new=np.genfromtxt('/group_workspaces/jasmin2/ukca/vol2/dcw32/Obs/NINO_34.1870-2008.ANOM.txt')
sst_new=np.reshape(sst_new[:,1:],sst_new.shape[0]*12)
print sst_new.shape
enso_ra=np.zeros(sst_new.shape[0])
for k in range(sst_new.shape[0]):
  if k>1:
    enso_ra[k]=np.mean(sst_new[k-2:k+1])
  elif k>0:
    enso_ra[k]=np.mean(sst_new[k-1:k+1])
  else:
    enso_ra[k]=sst_new[k]
for mons in range(12):
        gm_mon=np.mean(enso_ra[mons::12])
        enso_ra[mons::12]=enso_ra[mons::12]-gm_mon
g,Qxx_den=signal.periodogram(enso_ra,12.)
# Autocorr lag1
#plt.plot(np.arange(0,480),enso_detrend[:480],c='k')
#plt.plot(np.arange(0,480),enso_model[:480],c='r')
ax.plot(g, Qxx_den ,color='#E89CAE',linewidth=2)
ax.plot(f, Pxx_den ,color='#D50032',linewidth=2)
print str(round(1.0/(g[np.where(Qxx_den==np.max(Qxx_den))[0][0]]),2))
print str(round(1.0/(g[np.where(Qxx_den==second_largest(Qxx_den))[0][0]]),2))
ax.text(0.9,0.83,"Observations",transform=ax.transAxes\
         ,verticalalignment='center',horizontalalignment='right'\
         ,fontsize='14',color='#E89CAE')
ax.text(0.9,0.9,"Model",transform=ax.transAxes\
         ,verticalalignment='center',horizontalalignment='right'\
         ,fontsize='14',color='#D50032')
ax2=ax.twiny()
ax2.set_xlim([0.,1.])
ax2.set_xticks([0.1,0.2,0.5,1.])
ax2.set_xticklabels(['10','5','2','1'])
ax2.set_xlabel('Period / yr')
if save==True:
 if not os.path.exists(outdir):
  os.makedirs(outdir)
 plt.savefig(outdir+outname+'.'+format,bbox_inches="tight",dpi=400)
 print "SAVING "+outdir+outname+'.'+format
plt.show()
