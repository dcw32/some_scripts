from netCDF4 import Dataset
import numpy as np
import pylab as plt
import time
import os
from gridbox_areas import *
from matplotlib import gridspec
loc='/shared/netscratch/dcw32/um/'
jobid_base='xnfiy'
jobids1=['xnofa','xnofb','xnofc','xnofd','xnofe','xnoff']
jobids2=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl']
jobids3=['xnofn','xnofo','xnofp','xnofq','xnofr','xnofs']
#jobids3=['xnsvt','xnsvu','xnsvv','xnsvw','xnsvx','xnsvy']
#jobids3=['xnwdu','xnwdy']
#jobids3=['xnwdn','xnwdo','xnwdp','xnwdq','xnwdr','xnwds']
#jobids3=['xnsvt','xnsvw','xnsvx','xnsvy']
#jobids3=['xnwda','xnwdb','xnwdc','xnwdd','xnwde','xnwdf']
jobids4=['xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
labelm=['HI-HAL','LO-HAL','HI-SO2','LO-SO2']
cols=['#E87722','#4E5B31','#D50032','#93328e','#003C71']
njobs=len(jobids3)
var2='solar' #net sw down
var3='longwave' #net lw down
var4='sh' # sensible
var5='lh' # latent
var0='field201' # SW OUT
var1='olr' # LW OUT
#levs=[-50.,-30.,-20,-10.,10.,20.,30.,50.]
outdir='./../figures/'
outname='fig_s4a'
save=True
format='png'
################################################################
# BEGIN
################################################################
baseline_file=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_atmos_long.nc')
baseline_var0=baseline_file.variables[var0][:].squeeze()
baseline_var1=baseline_file.variables[var1][:].squeeze()
baseline_var2=baseline_file.variables[var2][:].squeeze()
baseline_var3=baseline_file.variables[var3][:].squeeze()
baseline_var4=baseline_file.variables[var4][:].squeeze()
baseline_var5=baseline_file.variables[var5][:].squeeze()
baseline_t=baseline_file.variables['t'][:]
lons=baseline_file.variables['longitude'][:]
lats=baseline_file.variables['latitude'][:]
fig=plt.figure(figsize=(10,5))
areas1=gbox_areas(len(lats),len(lons))
gs=gridspec.GridSpec(1,2)
gs.update(wspace=0.025, hspace=0.02)
ax0=fig.add_subplot(gs[0])
ax1=fig.add_subplot(gs[1])
#ax0.text(0.05,1.05,'A',transform=ax0.transAxes\
#                        ,verticalalignment='bottom',horizontalalignment='left'\
#                       ,fontsize='12',color='k')
#For each jobid
def plot_on_toa(ax,jobids,njobs,loc,areas,col,lab,verbose,mm):
	base_sat0=np.zeros(len(baseline_t))
	base_sat1=np.zeros(len(baseline_t))
	store0=np.zeros([mm,njobs])
	store1=np.zeros([mm,njobs])
	for l in range(len(baseline_t)):
		base_sat0[l]=np.average(baseline_var0[l,:,:],weights=areas)
		base_sat1[l]=np.average(baseline_var1[l,:,:],weights=areas)
	base_sat=0.-base_sat0-base_sat1
	for i in range(njobs):
		#Extract the data
		print jobids[i]
		anomaly_file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_atmos.nc')
		anomaly_var0=anomaly_file.variables[var0][:].squeeze()
		anomaly_var1=anomaly_file.variables[var1][:].squeeze()
		anomaly_t=anomaly_file.variables['t'][:]
		nts=len(anomaly_t)
		anom_sat0=np.zeros(len(anomaly_t))
		anom_sat1=np.zeros(len(anomaly_t))
		for l in range(len(anomaly_t)):
			anom_sat0[l]=np.average(anomaly_var0[l,:,:],weights=areas)
			anom_sat1[l]=np.average(anomaly_var1[l,:,:],weights=areas)
		#Find the location of the corresponding initial time in the base array
		base_loc=np.where(baseline_t==anomaly_t[0])[0][0]
		#Subtract the anomaly from the corresponding locations in the base
		delta0=anom_sat0-base_sat0[base_loc:base_loc+nts]
		delta1=anom_sat1-base_sat1[base_loc:base_loc+nts]
		#Plot time in months
		t_plot=(anomaly_t-anomaly_t[0])/360.
		store0[:,i]=delta0[:mm]
		store1[:,i]=delta1[:mm]
		#ax.plot(t_plot,-delta0,c='#85B09A',alpha=0.5)
		#ax.plot(t_plot,-delta1,c='#85B09A',alpha=0.5)
		#ax.plot(t_plot,0.-delta0-delta1,c='#85B09A',alpha=0.5)
	store0=np.mean(store0,axis=1)
	store1=np.mean(store1,axis=1)
	if verbose==True:
		ax.plot(t_plot[:mm],-store0,c='darkorange',linewidth=2.,label='SW')
		ax.plot(t_plot[:mm],-store1,c='darkred',linewidth=2.,label='LW')
		ax.plot(t_plot[:mm],0.-store0-store1,c='k',linewidth=2.,label='Net')
	else:
                ax.plot(t_plot[:mm],0.-store0-store1,c=col,linewidth=2.,label=lab)
def plot_on_surf(ax,jobids,njobs,loc,areas,col,lab,verbose,mm):
        base_sat2=np.zeros(len(baseline_t))
        base_sat3=np.zeros(len(baseline_t))
        base_sat4=np.zeros(len(baseline_t))
        base_sat5=np.zeros(len(baseline_t))
        store2=np.zeros([mm,njobs])
        store3=np.zeros([mm,njobs])
        store4=np.zeros([mm,njobs])
        store5=np.zeros([mm,njobs])
        for l in range(len(baseline_t)):
                base_sat2[l]=np.average(baseline_var2[l,:,:],weights=areas)
                base_sat3[l]=np.average(baseline_var3[l,:,:],weights=areas)
                base_sat4[l]=np.average(baseline_var4[l,:,:],weights=areas)
                base_sat5[l]=np.average(baseline_var5[l,:,:],weights=areas)
        for i in range(njobs):
                #Extract the data
                print jobids[i]
                anomaly_file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_atmos.nc')
                anomaly_var2=anomaly_file.variables[var2][:].squeeze()
                anomaly_var3=anomaly_file.variables[var3][:].squeeze()
                anomaly_var4=anomaly_file.variables[var4][:].squeeze()
                anomaly_var5=anomaly_file.variables[var5][:].squeeze()
                anomaly_t=anomaly_file.variables['t'][:]
                nts=len(anomaly_t)
                anom_sat2=np.zeros(len(anomaly_t))
                anom_sat3=np.zeros(len(anomaly_t))
                anom_sat4=np.zeros(len(anomaly_t))
                anom_sat5=np.zeros(len(anomaly_t))
                for l in range(len(anomaly_t)):
                        anom_sat2[l]=np.average(anomaly_var2[l,:,:],weights=areas)
                        anom_sat3[l]=np.average(anomaly_var3[l,:,:],weights=areas)
                        anom_sat4[l]=np.average(anomaly_var4[l,:,:],weights=areas)
                        anom_sat5[l]=np.average(anomaly_var5[l,:,:],weights=areas)
                #Find the location of the corresponding initial time in the base array
                base_loc=np.where(baseline_t==anomaly_t[0])[0][0]
                #Subtract the anomaly from the corresponding locations in the base
                delta2=anom_sat2-base_sat2[base_loc:base_loc+nts]
                delta3=anom_sat3-base_sat3[base_loc:base_loc+nts]
                delta4=anom_sat4-base_sat4[base_loc:base_loc+nts]
                delta5=anom_sat5-base_sat5[base_loc:base_loc+nts]
                #Plot time in months
                t_plot=(anomaly_t-anomaly_t[0])/360.
                store2[:,i]=delta2[:mm]
                store3[:,i]=delta3[:mm]
                store4[:,i]=delta4[:mm]
                store5[:,i]=delta5[:mm]
                #ax.plot(t_plot,-delta0,c='#85B09A',alpha=0.5)
                #ax.plot(t_plot,-delta1,c='#85B09A',alpha=0.5)
                #ax.plot(t_plot,0.-delta0-delta1,c='#85B09A',alpha=0.5)
        store2=np.mean(store2,axis=1)
        store3=np.mean(store3,axis=1)
        store4=np.mean(store4,axis=1)
        store5=np.mean(store5,axis=1)
	if verbose==True:
        	ax.plot(t_plot[:mm],store2,c='darkorange',linewidth=2.,label='SW')
        	ax.plot(t_plot[:mm],store3,c='darkred',linewidth=2.,label='LW')
        	ax.plot(t_plot[:mm],-store4,c='darkgreen',linewidth=2.,label='Sensible')
        	ax.plot(t_plot[:mm],-store5,c='darkblue',linewidth=2.,label='Latent')
        	ax.plot(t_plot[:mm],store2+store3-store4-store5,c='k',linewidth=2.,label='Net')
	else:
                ax.plot(t_plot[:mm],store2+store3-store4-store5,c=col,linewidth=2.,label=lab)
plot_on_toa(ax0,jobids3,njobs,loc,areas1,cols[2],labelm[2],True,60)
#plot_on_toa(ax0,jobids2,njobs,loc,areas1,cols[1],labelm[1])
#plot_on_toa(ax0,jobids3,njobs,loc,areas1,cols[2],labelm[2])
#plot_on_toa(ax0,jobids4,njobs,loc,areas1,cols[3],labelm[3])
plot_on_surf(ax1,jobids3,njobs,loc,areas1,cols[2],labelm[2],True,60)
#plot_on(0,baseline_t,baseline_var,jobids2,njobs,var,loc,areas1,1)
#plot_on(0,baseline_t,baseline_var,jobids3,njobs,var,loc,areas1,1)
#plot_on(0,baseline_t,baseline_var,jobids4,njobs,var,loc,areas1,1)
#plot_on(1,baseline_t,baseline_var,jobids2,njobs,var,loc,areas)
ax0.set_xlim([0,5])
ax1.set_xlim([0,5])
ax0.set_ylim([-20,20])
ax1.set_ylim([-20,20])
ax0.set_xticks([1,2,3,4,5])
ax0.set_xticklabels(['1','2','3','4','5'])
ax1.set_xticks([1,2,3,4,5])
ax1.set_xticklabels(['1','2','3','4','5'])
ax0.axhline(0.,c='k')
ax1.axhline(0.,c='k')
ax0.legend(fancybox=True,ncol=3,handletextpad=0.1,columnspacing=0.2)
ax1.legend(fancybox=True,ncol=3,handletextpad=0.1,columnspacing=0.2)
ax1.yaxis.tick_right()
ax1.yaxis.set_ticks_position('both')
ax0.text(0.02,1.02,'A',transform=ax0.transAxes\
                        ,verticalalignment='bottom',horizontalalignment='left'\
                        ,fontsize='14',color='k')
ax0.text(0.95,0.05,'Top of Atmosphere',transform=ax0.transAxes\
                        ,verticalalignment='bottom',horizontalalignment='right'\
                        ,fontsize='14',color='k')
ax1.text(0.95,0.05,'Surface',transform=ax1.transAxes\
                        ,verticalalignment='bottom',horizontalalignment='right'\
                        ,fontsize='14',color='k')
plt.figtext(0.08,0.5,r"Change to Radiative Budget / W m$\mathregular{^{-2}}$"\
         ,verticalalignment='center',horizontalalignment='center'\
         ,fontsize='12',color='k',rotation='90')
#plt.figtext(0.5,0.03,"Years After Eruption"\
#         ,verticalalignment='center',horizontalalignment='center'\
#         ,fontsize='12',color='k')
if save==True:
 if not os.path.exists(outdir):
  os.makedirs(outdir)
 plt.savefig(outdir+outname+'.'+format,bbox_inches="tight",dpi=200)
 print "SAVING "+outdir+outname+'.'+format
plt.show()
