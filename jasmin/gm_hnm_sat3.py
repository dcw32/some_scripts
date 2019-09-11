from netCDF4 import Dataset
import numpy as np
import pylab as plt
import time
import os
from gridbox_areas import *
from matplotlib import gridspec
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
#jobid_base='xnfiy'
jobid_base='xnwdg'
#jobids1=['xnwda','xnwdb','xnwdc','xnwdd','xnwde','xnwdf']
jobids1=['xnywa','xnywb','xnywc','xnywd','xnywf']
#jobids2=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl']
jobids2=['xnywg','xnywh','xnywi','xnywj','xnywk','xnywl']
#jobids3=['xnofn','xnofo','xnofp','xnofq','xnofr','xnofs']
jobids3=['xnywn','xnywo','xnywp','xnywq','xnywr','xnyws']
#jobids4=['xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
jobids4=['xnywt','xnywu','xnywv','xnyww','xnywx','xnywy']
labell=['Global Mean','Northern Hemisphere','Land North of 40$\mathregular{\degree}$N']
labelm=['HI-HAL','LO-HAL','HI-SO2','LO-SO2']
#jobids=['xnofa','xnofb','','xnofc','xnofd','xnofe','xnoff','xnofg','xnofh','xnofi','xnofj','xnofk','xnofl','xnofn','xnofo','xnofp','xnofq','xnofr','xnofs','xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
njobs=len(jobids1)
#label='HI-SO2+HAL'
#label='LO-SO2+HAL'
#label='HI-SO2'
label='LO-SO2'
var='temp'
#levs=[-50.,-30.,-20,-10.,10.,20.,30.,50.]
outdir='/home/users/dcw32/figures/thesis/samalas/'
outname='gl_satdelta_timeseries'+time.strftime("%Y%m%d")
save=False
format='pdf'
################################################################
# BEGIN
################################################################
baseline_file=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_pm.nc')
baseline_var=baseline_file.variables[var][:].squeeze()
baseline_t=baseline_file.variables['t'][:]
lons=baseline_file.variables['longitude'][:]
lats=baseline_file.variables['latitude'][:]
fig=plt.figure(figsize=(8,11))
areas1=gbox_areas(len(lats),len(lons))
areas2=np.copy(areas1)
areas3=np.copy(areas1)
for i in range(len(lats)):
	if lats[i]<0.:
		areas2[i,:]=0.
	if lats[i]==0.:
		areas2[i,:]=0.5*areas2[i,:]
	if lats[i]<40.:
		areas3[i,:]=0.
gs=gridspec.GridSpec(4,3)
gs.update(wspace=0.025, hspace=0.025)
#For each jobid
def plot_on(number,baseline_t,baseline_var,jobids,njobs,var,loc,areas,yy):
	ax=fig.add_subplot(gs[number])
	base_sat=np.zeros(len(baseline_t))
	store=np.zeros([120,njobs])
	for l in range(len(baseline_t)):
		base_sat[l]=np.average(baseline_var[l,:,:],weights=areas)
	for i in range(njobs):
		#Extract the data
		print jobids[i]
		anomaly_file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_atmos.nc')
		anomaly_var=anomaly_file.variables[var][:].squeeze()
		anomaly_t=anomaly_file.variables['t'][:]
		nts=len(anomaly_t)
		anom_sat=np.zeros(len(anomaly_t))
		for l in range(len(anomaly_t)):
			anom_sat[l]=np.average(anomaly_var[l,:,:],weights=areas)
		#Find the location of the corresponding initial time in the base array
		base_loc=np.where(baseline_t==anomaly_t[0])[0][0]
		#Subtract the anomaly from the corresponding locations in the base
		delta=anom_sat-base_sat[base_loc:base_loc+nts]
		#Plot time in months
		t_plot=(anomaly_t-anomaly_t[0])/30.
		store[:,i]=delta[:120]
		if number%3==2:
			#ax.plot(t_plot,delta,c='#E89CAE',alpha=0.5)
			ax.plot(t_plot,delta,c='#6CACE4',alpha=0.5)
		else:
			ax.plot(t_plot,delta,c='#85B09A',alpha=0.5)
	store=np.mean(store,axis=1)
	xx=-100
	for l in range(120-12):
		if np.mean(store[l:l+12])>0.:
			xx=l
			print str(l)+" MONTH OF FIRST RECOVERY"
			break
	print t_plot.shape
	print store.shape
	if number%3==2:
		#ax.plot(t_plot[:120],store,c='#8A1538',linewidth=2.)
		ax.plot(t_plot[:120],store,c='#003C71',linewidth=2.)
	else:
		ax.plot(t_plot[:120],store,c='#095D67',linewidth=2.)
	#ax.set_xlabel('Months After Eruption')
	#ax.set_ylabel('Surface Temperature Change / C')
	ax.set_xlim([0,120])
	if yy==1:
		ax.set_ylim([-1.5,1.0])
	elif yy==2:
		ax.set_ylim([-2.0,1.5])
	ax.set_xticks([24,48,72,96,120])
	if number>8:
		ax.set_xticklabels(['2','4','6','8','10'])
	else:
		ax.set_xticklabels(['']*4)
	#if yy==1:
	#	ax.axvline(xx,c='#EF3340',linestyle='--')
	ax.axhline(0.,c='k')
	if number%3==2:
                ax.set_ylim([-2.5,1.5])
		ax.set_yticks([-2,-1,0,1])
		ax.yaxis.tick_right()
		ax.yaxis.set_ticks_position('both')
	else:
		ax.set_ylim([-1.5,1.0])
		ax.set_yticks([-1,0,1])
	if number%3==1:
		ax.set_yticklabels(['']*3)
	if number<3:
		ax.set_title(labell[number],fontsize=13)
	if number%3==0:
		xx=number/3
		ax.text(0.05,0.95,labelm[xx],transform=ax.transAxes\
			,verticalalignment='top',horizontalalignment='left'\
			,fontsize='12',color='#095D67')
	#plt.title("SAT DELTA "+label)
plot_on(0,baseline_t,baseline_var,jobids1,njobs,var,loc,areas1,1)
plot_on(1,baseline_t,baseline_var,jobids1,njobs,var,loc,areas2,1)
plot_on(2,baseline_t,baseline_var,jobids1,njobs,var,loc,areas3,2)
plot_on(3,baseline_t,baseline_var,jobids2,njobs,var,loc,areas1,1)
plot_on(4,baseline_t,baseline_var,jobids2,njobs,var,loc,areas2,1)
plot_on(5,baseline_t,baseline_var,jobids2,njobs,var,loc,areas3,2)
plot_on(6,baseline_t,baseline_var,jobids3,njobs,var,loc,areas1,1)
plot_on(7,baseline_t,baseline_var,jobids3,njobs,var,loc,areas2,1)
plot_on(8,baseline_t,baseline_var,jobids3,njobs,var,loc,areas3,2)
plot_on(9,baseline_t,baseline_var,jobids4,njobs,var,loc,areas1,1)
plot_on(10,baseline_t,baseline_var,jobids4,njobs,var,loc,areas2,1)
plot_on(11,baseline_t,baseline_var,jobids4,njobs,var,loc,areas3,2)
#plot_on(1,baseline_t,baseline_var,jobids2,njobs,var,loc,areas)
plt.figtext(0.06,0.5,r"Temperature Anomaly / $\mathregular{\degree}$C"\
         ,verticalalignment='center',horizontalalignment='center'\
         ,fontsize='14',color='k',rotation='90')
plt.figtext(0.5,0.06,"Years After Eruption"\
         ,verticalalignment='center',horizontalalignment='center'\
         ,fontsize='14',color='k')
if save==True:
 if not os.path.exists(outdir):
  os.makedirs(outdir)
 plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
 print "SAVING "+outdir+outname+'.'+format
plt.show()
