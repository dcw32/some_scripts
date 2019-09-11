from netCDF4 import Dataset
import numpy as np
import pylab as plt
import time
from scipy import signal
import os
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
jobid='xnfiy'
var='v_1'
xs=0.5+np.array([4.,46.,27.,17.,21.,23.])
xs=0.5+np.array([6.,12.,18.,24.,30.,36.])
nens=xs.shape[0]
levs=np.arange(-30.,30.1,5.0)
#levs=[-50.,-30.,-20,-10.,10.,20.,30.,50.]
outdir='/home/users/dcw32/figures/samalas/'
outname=jobid+'_qbo_psd_at30hPa_timeseries'+time.strftime("%Y%m%d")
save=True
format='png'
################################################################
# BEGIN
################################################################
#Extract Data
file=Dataset(loc+jobid+'/netcdf/'+jobid+'_clim.nc')
lats=file.variables['latitude'][:]
lons=file.variables['longitude'][:]
pres=file.variables['p_1'][:]
time=file.variables['t'][:]
qbo=file.variables['u_1'][:]
#Identify limits of 2.5S-2.5N
latbox_min=np.where(lats==-1.25)[0][0]
latbox_max=np.where(lats==1.25)[0][0]
#Zonally average
qbo=np.mean(qbo,axis=3)
qbo=np.mean(qbo[:,:,latbox_min:latbox_max+1],axis=2)
#Time in days since 2160 --> Years from 0
time=time/360.0
time=time-time[0]
qbo=qbo[:,np.where(pres==30.0)[0][0]]
print qbo.shape
#Plot
fig=plt.figure(figsize=(6,6))
ax=fig.add_subplot(1,1,1)
#f,Pxx_den=signal.periodogram(qbo[:480],3.8580247E-7)
f,Pxx_den=signal.periodogram(qbo[:480],12.)
#f=f*31104000.0
#plt.semilogy(f, Pxx_den,color='k')
#plt.plot(f, Pxx_den,c='#D50032',linewidth=2)
ax.set_xlim([0,1])
#plt.ylim([1E3,1E12])
ax.set_xlabel(r'Frequency / yr$^{-1}$')
ax.set_ylabel(r'Power Spectral Density')
qbo_freq=1.0/(f[np.where(Pxx_den==np.max(Pxx_den))[0][0]])
#plt.axvline(f[np.where(Pxx_den==np.max(Pxx_den))[0][0]],c='k',linestyle='--')
ax.text(0.9,0.9,"Model",transform=ax.transAxes\
         ,verticalalignment='center',horizontalalignment='right'\
         ,fontsize='14',color='#D50032')
# DATA
data=np.genfromtxt('/group_workspaces/jasmin2/ukca/vol2/dcw32/Obs/qbo/qbo.txt',delimiter=' ')
print data.shape
qbo_sing=np.reshape(data,[data.shape[0]*data.shape[1]])*0.1
#g,Qxx_den=signal.periodogram(qbo_sing,3.8580247E-7)
g,Qxx_den=signal.periodogram(qbo_sing,12.)
#g=g*31104000.0
qbo_freq2=1.0/(g[np.where(Qxx_den==np.max(Qxx_den))[0][0]])
ax.text(0.9,0.85,"Observations",transform=ax.transAxes\
         ,verticalalignment='center',horizontalalignment='right'\
         ,fontsize='14',color='#E89CAE')
ax.plot(g,Qxx_den,c='#E89CAE',linewidth=2)
ax.plot(f, Pxx_den,c='#D50032',linewidth=2)
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.yaxis.set_ticks_position('both')
ax2=ax.twiny()
ax2.set_xlim([0.,1.])
ax2.set_xticks([0.1,0.2,0.5,1.])
ax2.set_xticklabels(['10','5','2','1'])
ax2.set_xlabel('Period / yr')
#plt.title(jobid+" PSD for tropicalU at 30hPa")
if save==True:
 if not os.path.exists(outdir):
  os.makedirs(outdir)
 plt.savefig(outdir+outname+'.'+format,bbox_inches="tight",dpi=400)
 print "SAVING "+outdir+outname+'.'+format
plt.show()
