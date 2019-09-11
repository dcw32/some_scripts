from matplotlib import gridspec
from netCDF4 import Dataset
import matplotlib as mpl
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
jobids4=['xnofu','xnofv','xnofw','xnofx','xnofy']
jobid_base='xnfiy'
var='field1076'
print var
# Plotting Arguments
levs=[-0.5,-0.3,-0.1,0.1,0.3,0.5]

cmap=cmap=plt.cm.get_cmap('RdBu_r')
levs2=np.insert(levs,0,levs[0]-1)
levs2=np.append(levs2,levs2[len(levs2)-1]+1)
norm=mpl.colors.BoundaryNorm(levs2, ncolors=cmap.N, clip=True)
outdir='/home/users/dcw32/figures/thesis/samalas/'
outname='wstarbar_diff_plot'+time.strftime("%Y%m%d")
save=True
format='pdf'
extend='both'
####################################################################
#
# Extract w*bar climatology
#
file=Dataset(loc+jobid_base+'/netcdf/'+jobid_base+'_dyn.nc')
lats=file.variables['latitude'][:]
pres=file.variables['p_1'][:]
time=file.variables['t'][:]
wstar=file.variables[var][:,np.where(pres==70.0)[0][0],:].squeeze()
clima=np.zeros([12,len(lats)])
climb=np.zeros([12,len(lats)])
for mon in range(12):
        clima[mon,:]=np.mean(wstar[mon::12,:],axis=0)
        climb[mon,:]=1.96*np.std(wstar[mon::12,:],axis=0)
clima=clima*1000
climb=climb*1000
print np.mean(clima)
clima=np.roll(clima,7,axis=0)
climb=np.roll(climb,7,axis=0)
fig=plt.figure()
gs1 = gridspec.GridSpec(2, 2)
gs1.update(wspace=0.0)
ax1 = plt.subplot(gs1[0,0])
ax2 = plt.subplot(gs1[0,1])
ax3 = plt.subplot(gs1[1,0])
ax4 = plt.subplot(gs1[1,1])
## Function to plot for each run
def wstarb(axis,jobids,mm,clima,l_right):
	njs=len(jobids)
	clima_now=clima
        clima_now=np.tile(clima,(mm/12,1))
        climb_now=np.tile(climb,(mm/12,1))
        store1=np.zeros([njs,mm,72])
	for nj in range(njs):
		print loc+jobids[nj]+'/netcdf/'+jobids[nj]+'_dyn.nc'
		file=Dataset(loc+jobids[nj]+'/netcdf/'+jobids[nj]+'_dyn.nc')
		lats=file.variables['latitude'][:]
		pres=file.variables['p_1'][:]
		time=file.variables['t'][:mm]
		wstar=1000*file.variables[var][:mm,np.where(pres==70.0)[0][0],:].squeeze()
		store1[nj,:,:]=wstar-clima_now
	time=(time-time[0])/360.
	store=np.mean(store1,axis=0)
	stdevmap=np.abs(store)>climb_now
        cf1=axis.contourf(time,lats,np.transpose(store),levs,cmap=cmap,norm=norm,extend=extend)
	cf2=axis.contourf(time,lats,np.transpose(stdevmap),hatches=['xxx',None],cmap=plt.get_cmap('RdBu'),alpha=0.0,extend='both',rasterized=True)
        cf3=axis.contour(time,lats,np.transpose(clima_now),[0],colors='#4E5B31',linewidths=2)
	axis.set_yticks([-90,-45,0,45,90])
	axis.set_yticklabels(['90S','45S','EQ','45N','90N'])
        if l_right==True:
                axis.yaxis.tick_right()
                axis.yaxis.set_ticks_position('both')
	return cf1
cax=wstarb(ax1,jobids1,60,clima,False)
cax=wstarb(ax2,jobids2,60,clima,True)
cax=wstarb(ax3,jobids3,60,clima,False)
cax=wstarb(ax4,jobids4,60,clima,True)
ax1.text(0.04,1.02,'A: HI-HAL',transform=ax1.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='12',color='k')
ax2.text(0.04,1.02,'B: LO-HAL',transform=ax2.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='12',color='k')
ax3.text(0.04,1.02,'C: HI-SO2',transform=ax3.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='12',color='k')
ax4.text(0.04,1.02,'D: LO-SO2',transform=ax4.transAxes\
        ,verticalalignment='bottom',horizontalalignment='left'\
        ,fontsize='12',color='k')
cbaxes=fig.add_axes([0.0, 0.20, 0.02, 0.6])
cbar=plt.colorbar(cax,cax=cbaxes,orientation="vertical")
cbar.set_label(r'Change in $\overline{w}$* at 70hPa / mm s$^{-1}$')
cbaxes.yaxis.set_label_position('left')
plt.figtext(0.5,0.02,'Years After Eruption'\
        ,verticalalignment='bottom',horizontalalignment='center'\
        ,fontsize='12',color='k')
if save==True:
                print "SAVING "+outdir+outname+'.'+format
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight",dpi=200)
plt.show()
