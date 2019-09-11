from netCDF4 import Dataset
from matplotlib import gridspec
import numpy as np
import calendar as calendar
from gridbox_areas import gbox_areas
from plot_1D import *
import time
######################### Input ####################################
loc='/group_workspaces/jasmin2/ukca/vol2/dcw32/'
jobids1=['xnofa','xnofb','xnofc','xnofd','xnofe','xnoff']
jobids2=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl']
jobids3=['xnofn','xnofo','xnofp','xnofq','xnofr','xnofs']
jobids4=['xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
jobid_base='xnfiy'
#jobids5=['xnsvt','xnsvu','xnsvv','xnsvw','xnsvx','xnsvy']
#jobids1=['xnsvt']
#jobids2=['xnsvu']
#jobids3=['xnsvv']
#jobids4=['xnsvw']
#jobids5=['xnsvx']

cols=['#E87722','#4E5B31','#D50032','#93328e','#003C71']
gbvolf=Dataset('/group_workspaces/jasmin2/ukca/vol1/dcw32/n48_l60_geovol.nc')
gbvol=gbvolf.variables['vol_theta'][:].squeeze()
gbvolf.close()
#jobids=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl']
#var_1='temp_2' #AIT
#var_2='field614' #COA
#
varn='samalas_gm_ozonecol'
# Plotting Arguments
#levs=[50,100,150,200,250,300,350,400,450,500]
colmap='Spectral'
#cbar_label='Ozone Column / DU'
levs=[0.1,0.2,0.5,1.0,2.0,5.0]
cmap=plt.cm.get_cmap(colmap)
levs2=np.insert(levs,0,levs[0]-1)
levs2=np.append(levs2,levs2[len(levs2)-1]+1)
norm=mpl.colors.BoundaryNorm(levs2, ncolors=cmap.N, clip=True)
#levs=[0,60,120,180,240,300,360,420]
areas=gbox_areas(73,96)
outdir='/home/users/dcw32/figures/thesis_volcano/samalas/'
outname=varn+'_timese_plot'+time.strftime("%Y%m%d")
save=True
format='png'
extend='both'
mm1=240
mm2=120
mm3=120
mm4=120
mm5=25+12
#mm=240
###Climatology
file=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_oz.nc')
lats=file.variables['latitude'][:]
time=file.variables['t'][:]
wstar=file.variables['salinity'][:,0,:,:].squeeze()
file.close()
wstar=np.mean(wstar,axis=2)
wstar=np.mean(wstar,axis=1)
clima=np.zeros([12])
for mon in range(12):
        clima[mon]=np.mean(wstar[mon::12])
clima=clima
clima=np.roll(clima,7)
print np.mean(clima)
####################################################################
#
# Extract
fig=plt.figure()
gs1 = gridspec.GridSpec(1, 1)
gs1.update(hspace=0.2)
ax1 = plt.subplot(gs1[0])
#ax2 = plt.subplot(gs1[1])
#ax3 = plt.subplot(gs1[2])
def extract_col(axis,jobids,mm,l_right,col,label):
	njob=len(jobids)
	store1=np.zeros([njob,mm])
        clima_now=np.tile(clima,(mm/12))
	for i in range(njob):
		file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_oz.nc')
		lats=file.variables['latitude'][:]
		lons=file.variables['longitude'][:]
		time=file.variables['t'][:mm]
		#field=file.variables[var][:].squeeze()+file.variables[var_1][:].squeeze()+file.variables[var_2][:].squeeze()
		field=file.variables['salinity'][:mm,0,:,:].squeeze()
		file.close()
		field=np.mean(field,axis=2)
		mass=np.mean(field,axis=1)
		#TimexLat
		#field=np.mean(field,axis=2)
		store1[i,:]=mass[:mm]-clima_now
		time=time/360.
		time=time-time[0]
		#axis.plot(time,mass,linewidth=1,linestyle=':',color=col)
	store1=np.mean(store1,axis=0)
	print label
	print np.max(store1)
	print np.min(store1)
	print np.sum(store1)/12.
	# Start Plotting
	cf1=axis.plot(time,store1,linewidth=2,color=col,label=label)
	#cf1=axis.contourf(time[:mm],lats,np.transpose(store),levs,cmap=cmap,norm=norm,extend='both')
	#axis.contour(time,lats,np.transpose(store),[220],colors='k',linewidths=4)
        #axis.set_yticks([-90,-45,0,45,90])
        #axis.set_yticklabels(['90S','45S','EQ','45N','90N'])
	#if l_right==True:
        #        axis.yaxis.tick_right()
        #        axis.yaxis.set_ticks_position('both')
	return cf1
#Aerosol mass
#cf1=extract_col(ax1,jobids1,mm1,False,cols[0],'HI-HAL')
#cf1=extract_col(ax1,jobids2,mm2,True,cols[1],'LO-HAL')
cf1=extract_col(ax1,['xnofa'],mm1,True,cols[0],'20%')
cf1=extract_col(ax1,['xnwdj'],mm2,True,'forestgreen','10%')
cf1=extract_col(ax1,['xnwdi'],mm2,True,'deeppink','5%')
cf1=extract_col(ax1,['xnofg'],mm2,True,cols[1],'1%')
#cf1=extract_col(ax1,jobids3,mm3,False,cols[2],'HI-SO2')
#cf1=extract_col(ax1,jobids4,mm4,True,cols[3],'LO-SO2')
#cf1=extract_col(ax1,jobids5,mm5,True,cols[4],'PIN20')
ax1.legend(fancybox=True,loc='lower right')
ax1.set_xlim(0,9.9)
#ax1.set_xlabel('Months After Eruption')
ax1.set_ylabel('Ozone Column Change / DU')
#Aerosol number
#        ,verticalalignment='top',horizontalalignment='center'\
#        ,fontsize='14',color='k')
if save==True:
                print "SAVING "+outdir+outname+'.'+format
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")

#plt.savefig(outdir+outname+'.'+format)
plt.show()
