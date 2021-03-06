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
#jobids1=['xnsvt','xnsvu','xnsvv','xnsvw','xnsvx','xnsvy']
jobids2=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl']
jobids3=['xnofn','xnofo','xnofp','xnofq','xnofr','xnofs']
jobids4=['xnoft','xnofu','xnofv','xnofw','xnofx','xnofy']
cols=['#E87722','#4E5B31','#D50032','#EF3340']
gbvolf=Dataset('/group_workspaces/jasmin2/ukca/vol1/dcw32/n48_l60_geovol.nc')
gbvol=gbvolf.variables['vol_theta'][:].squeeze()
gbvolf.close()
#jobids=['xnofg','xnofh','xnofi','xnofj','xnofk','xnofl']
var='field608'
#var_1='temp_2' #AIT
#var_2='field614' #COA
var2='field643'
#
varn='samalas_sulfate_glob'
print var
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
outname=varn+'_surftimese_plot'+time.strftime("%Y%m%d")
save=True
format='png'
extend='both'
mm1=25+12
mm2=25+12
mm3=25+12
mm4=25+12
#mm=240
####################################################################
#
# Extract
fig=plt.figure()
gs1 = gridspec.GridSpec(2, 1)
gs1.update(hspace=0.05)
ax1 = plt.subplot(gs1[0])
ax2 = plt.subplot(gs1[1])
#ax3 = plt.subplot(gs1[2])
def extract_col(axis,jobids,mm,l_right,col,label):
	njob=len(jobids)
	store1=np.zeros([njob,mm])
	for i in range(njob):
		file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_aer.nc')
		file2=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_oz.nc')
		lats=file.variables['latitude'][:]
		lons=file.variables['longitude'][:]
		time=file.variables['t'][:mm]
		#field=file.variables[var][:].squeeze()+file.variables[var_1][:].squeeze()+file.variables[var_2][:].squeeze()
		field=file.variables[var][:mm,:,:,:].squeeze()
		field2=file2.variables[var2][:mm,:,:,:].squeeze()
		file.close()
		file2.close()
		mass_boxes=field*field2
		mass=np.sum(mass_boxes,axis=3)
		mass=np.sum(mass,axis=2)
		mass=np.sum(mass,axis=1)
		mass=mass*10**-9
		#TimexLat
		#field=np.mean(field,axis=2)
		store1[i,:]=mass[:mm]
		time=time/30.
		time=time-time[0]
		axis.plot(time,mass,linewidth=1,linestyle=':',color=col)
	store1=np.mean(store1,axis=0)
	# Start Plotting
	print time.shape
	print store1.shape
	cf1=axis.plot(time,store1,linewidth=2,color=col,label=label)
	#cf1=axis.contourf(time[:mm],lats,np.transpose(store),levs,cmap=cmap,norm=norm,extend='both')
	#axis.contour(time,lats,np.transpose(store),[220],colors='k',linewidths=4)
        #axis.set_yticks([-90,-45,0,45,90])
        #axis.set_yticklabels(['90S','45S','EQ','45N','90N'])
	#if l_right==True:
        #        axis.yaxis.tick_right()
        #        axis.yaxis.set_ticks_position('both')
	return cf1
def extract_num(axis,jobids,mm,l_right,col,label):
        njob=len(jobids)
        store1=np.zeros([njob,mm])
        for i in range(njob):
                file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_aer.nc')
                file2=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_oz.nc')
                mass=file2.variables[var2][:mm,:,:,:].squeeze()
                lats=file.variables['latitude'][:]
                lons=file.variables['longitude'][:]
                time=file.variables['t'][:mm]
                #field=file.variables[var][:].squeeze()+file.variables[var_1][:].squeeze()+file.variables[var_2][:].squeeze()
                field=file.variables[var][:mm,:,:,:].squeeze()
		field=field*mass*6.022E23/0.02897
                #field3=file.variables['field1878_6'][:mm,:,:,:].squeeze()
		#field=field/(field3**3)
		#field=field*300*np.pi
                #mass_boxes=field*field2
                mass_boxes=field
		file.close()
                #mass_boxes=mass_boxes*gbvol*1E6
                mass=np.sum(mass_boxes,axis=3)
                mass=np.sum(mass,axis=2)
                mass=np.sum(mass,axis=1)
                mass=np.divide(mass,np.sum(gbvol*1E6))
                #TimexLat
                #field=np.mean(field,axis=2)
                store1[i,:]=mass[:mm]
                time=time/30.
                time=time-time[0]
                axis.plot(time,mass,linewidth=1,linestyle=':',color=col)
        store1=np.mean(store1,axis=0)
        # Start Plotting
        print time.shape
        print store1.shape
        cf1=axis.plot(time,store1,linewidth=2,color=col,label=label)
        #cf1=axis.contourf(time[:mm],lats,np.transpose(store),levs,cmap=cmap,norm=norm,extend='both')
        #axis.contour(time,lats,np.transpose(store),[220],colors='k',linewidths=4)
        #axis.set_yticks([-90,-45,0,45,90])
        #axis.set_yticklabels(['90S','45S','EQ','45N','90N'])
        #if l_right==True:
        #        axis.yaxis.tick_right()
        #        axis.yaxis.set_ticks_position('both')
        return cf1
def extract_reff(axis,jobids,mm,l_right,col,label):
        njob=len(jobids)
        store1=np.zeros([njob,mm])
        for i in range(njob):
                file=Dataset(loc+jobids[i]+'/netcdf/'+jobids[i]+'_aer.nc')
                lats=file.variables['latitude'][:]
                lons=file.variables['longitude'][:]
                time=file.variables['t'][:mm]
                #field=file.variables[var][:].squeeze()+file.variables[var_1][:].squeeze()+file.variables[var_2][:].squeeze()
                field=file.variables['field1878_6'][:mm,:,:,:].squeeze()
                field2=file.variables['pad_2'][:mm,:,:,:].squeeze()
                file.close()
		reff1=np.zeros([len(time)])
		for ti in range(len(time)):
		        reff1[ti]=np.average(field[ti,:,:,:],weights=field2[ti,:,:,:])
		reff1=reff1*1E6*1.327/2.0
                #mass_boxes=mass_boxes*gbvol*1E6
                store1[i,:]=reff1[:mm]
                time=time/30.
                time=time-time[0]
                axis.plot(time,reff1,linewidth=1,linestyle=':',color=col)
        store1=np.mean(store1,axis=0)
        # Start Plotting
        print time.shape
        print store1.shape
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
cf1=extract_col(ax1,jobids1,mm1,False,cols[0],'HI-HAL')
cf1=extract_col(ax1,jobids2,mm2,True,cols[1],'LO-HAL')
cf1=extract_col(ax1,jobids3,mm3,False,cols[2],'HI-SO2')
cf1=extract_col(ax1,jobids4,mm4,True,cols[3],'LO-SO2')
ax1.legend(fancybox=True)
ax1.set_xlim(0,24)
#ax1.set_xlabel('Months After Eruption')
ax1.set_ylabel('Aerosol Burden / Tg')
#Aerosol number
cf2=extract_reff(ax2,jobids1,mm1,False,cols[0],'HI-HAL')
cf2=extract_reff(ax2,jobids2,mm2,False,cols[1],'LO-HAL')
cf2=extract_reff(ax2,jobids3,mm3,False,cols[2],'HI-SO2')
cf2=extract_reff(ax2,jobids4,mm4,False,cols[3],'LO-SO2')
#ax2.legend(fancybox=True)
ax2.set_xlim(0,24)
ax2.set_xlabel('Months After Eruption')
ax2.set_ylabel('Aerosol Effective Radius / um')

#cbar=plt.colorbar(cf3,cax=cbaxes,orientation="vertical")
#cbar.set_label('Ozone Column / Dobson Units')
#cbaxes.yaxis.set_label_position('left')
#ax1.text(0.50,-0.12,'Years After Eruption',transform=ax1.transAxes\
#        ,verticalalignment='top',horizontalalignment='center'\
#        ,fontsize='14',color='k')
if save==True:
                print "SAVING "+outdir+outname+'.'+format
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")

#plt.savefig(outdir+outname+'.'+format)
plt.show()
