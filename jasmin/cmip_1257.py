from netCDF4 import Dataset
import numpy as np
from gridbox_areas import gbox_areas
import pylab as plt
import time
import sys
import csv
import pandas as pd
import os
from locate import locate2
######################### Input ####################################
loc2=locate2()
#This plots the Northern Hemisphere SAT over the 1257 volcanic eruption
#
def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)
#
#
if loc2=='jasmin':
	outdir='/home/users/dcw32/figures/thesis/samalas/'
	csvloc='/group_workspaces/jasmin2/ukca/vol2/dcw32/CMIP5/lmil_40_75_land/'
	proxyfile='/group_workspaces/jasmin2/ukca/vol2/dcw32/Obs/samalas/ntrend.csv'
	guilletfile='/group_workspaces/jasmin2/ukca/vol2/dcw32/Obs/samalas/guillet2017nhtemp.txt'
	sch15file='/group_workspaces/jasmin2/ukca/vol2/dcw32/Obs/samalas/SCH15.csv'
elif loc2=='cambridge':
        outdir='/home/dcw32/figures/thesis/samalas/'
        csvloc='/home/dcw32/CMIP5/lmil_40_75_land/'
        proxyfile='/home/dcw32/Obs/samalas/ntrend.csv'
        guilletfile='/home/dcw32/Obs/samalas/guillet2017nhtemp.txt'
        sch15file='/home/dcw32/Obs/samalas/SCH15.csv'
else:
	print "ERROR: UNKNOWN PLATFORM"
outname='cmip5_2_'+time.strftime("%Y%m%d")
save=True
format='pdf'
#
cols=['#F1BE48','#85B09A','#6CACE4','#E89CAE']
#models=['BCC','CCSM4','FGOALS','GISS1','GISS4','GISS7','MPI','IPSL','MIROC']
models=['BCC','CCSM4','GISS1','GISS4','GISS7','MPI','IPSL','MIROC']
nmodels=len(models)
#Note 2 of the GISS-ER do not include volcanic forcing, these are not used. An average of the remaining ensemble is used to determine the GISS average. The overall ensemble weights the 5 models equally
#fig=plt.figure()
fig,ax=plt.subplots(figsize=(6,4))
yrs=np.arange(850,1850)
#plt.fill_between(yrs,mmin,mmax,facecolor=cols[0],alpha=0.3)
#plt.plot(yrs,mmm,c=cols[0],linewidth=3,label='CMIP5')
#NTREND
tree_year=[]
tree_anom=[]
tree_usig=[]
tree_lsig=[]
data=pd.read_csv(proxyfile,delimiter=',',dtype=float)
#with open(proxyfile,'r') as csvfile:
#	data=csv.reader(csvfile, delimiter=',',dtype=float)
#	for row in data:
#		tree_year=np.append(tree_year,row[0])
#		tree_anom=np.append(tree_anom,row[1])
#		tree_usig=np.append(tree_usig,row[2])
#		tree_lsig=np.append(tree_lsig,row[3])
#tree_year=tree_year[xs:xs+1000]
tree_yrs=data['Year']
tree_mmm=data['NTREND2015']
tree_low=data['lower2sigma']
tree_high=data['upper2sigma']
print np.array(tree_mmm).shape
tree_rm=running_mean(np.array(tree_mmm),31)
tree_yrs=tree_yrs[15:-15]
tree_mmm=tree_mmm[15:-15]-tree_rm
tree_low=tree_low[15:-15]-tree_rm
tree_high=tree_high[15:-15]-tree_rm
print tree_lsig
plt.fill_between(tree_yrs,tree_low,tree_high,facecolor=cols[1],alpha=0.3)
l1,=plt.plot(tree_yrs,tree_mmm,c=cols[1],linewidth=3)
#
# Extract Schenider data
data_sch=np.genfromtxt(sch15file,skip_header=94,delimiter=' ')
sch_yrs=data_sch[:,0]
sch_mmm=data_sch[:,1]
sch_low=data_sch[:,2]
sch_high=data_sch[:,3]
sch_rm=running_mean(np.array(sch_mmm),31)
sch_yrs=sch_yrs[15:-15]
sch_mmm=sch_mmm[15:-15]-sch_rm
sch_low=sch_low[15:-15]-sch_rm
sch_high=sch_high[15:-15]-sch_rm
plt.fill_between(sch_yrs,sch_low,sch_high,facecolor=cols[2],alpha=0.3)
l2,=plt.plot(sch_yrs,sch_mmm,c=cols[2],linewidth=3)
#plt.xlim([850,1850])
##
#Stoffel values
#
yrs=np.array([1258,1259])
v1=np.array([-1.18,-1.31])
v2=np.array([-1.10,-1.23])
#plt.scatter(yrs,0.5*(v1+v2),marker='o',s=80,c='k',label='Stoffel',zorder=1000)
##
#Guillet values
#
data=np.genfromtxt(guilletfile,skip_header=140)
print data.shape
yrs=data[:,0]
up=data[:,7]
lo=data[:,6]
mn=data[:,5]
##
#for years in range(len(yrs)):
#	if years>0:
#		plt.plot([yrs[years],yrs[years]],[v1[years],v2[years]],c=cols[2],linewidth=3)
#	else:
#                plt.plot([yrs[years],yrs[years]],[v1[years],v2[years]],c=cols[2],linewidth=3,label='Stoffel')
#Guillet plot
plt.fill_between(yrs,lo,up,facecolor=cols[3],alpha=0.3)
l3,=plt.plot(yrs,mn,c=cols[3],linewidth=3)
fig.legend((l1,l2,l3),('NTREND','SCH15','SG17'),loc=9,bbox_to_anchor=(0.5,0.99),bbox_transform=ax.transAxes,ncol=3,handletextpad=0.2,columnspacing=0.8)
#CMIP5
#CMIP5
yrs=np.arange(850,1850)
print "CMIP5 MMM"
anoms=np.zeros([len(yrs),nmodels])
for model in range(nmodels):
        anom=[]
        with open(csvloc+models[model]+'.csv', 'r') as csvfile:
                data=csv.reader(csvfile, delimiter=',')
                for row in data:
                        anom=np.append(anom,row[1])
        anoms[:len(yrs),model]=anom[:len(yrs)]
mmm=np.mean(anoms,axis=1)
mmm_rm=running_mean(mmm,31)
mms=2.*np.std(anoms,axis=1)
mmax=np.max(anoms,axis=1)
mmin=np.min(anoms,axis=1)
yrs=yrs[15:-15]
mms=mms[15:-15]-mmm_rm
mmax=mmax[15:-15]-mmm_rm
mmin=mmin[15:-15]-mmm_rm
mmm=mmm[15:-15]-mmm_rm
#plt.fill_between(yrs,mmm-mms,mmm+mms,facecolor=cols[0],alpha=0.3)
m1=ax.errorbar(yrs-0.15,0.5*(mmax+mmin),yerr=0.5*(mmax-mmin),fmt=None,label='CMIP5-PMIP3',ecolor='#BE4D00',elinewidth=2,capthick=2,zorder=1001)
#Model values
#CESM-LME values
yrs=np.arange(1257,1257+10)
up=np.array([0.6235933, -0.96896905, -2.25005077, -0.82435928, -0.11464864, 0.38096276, 0.40699631, 0.24000263, 0.51365507, 0.39220913])
mid=np.array([0.22051601, -1.32054488, -2.87124004, -1.10222612, -0.36790776, -0.05325477, 0.06939979, -0.05777415,  0.09809921, 0.22515684])
lo=np.array([-0.10546944, -1.72477337, -3.20228151, -1.32906544, -0.77845975, -0.42046459, -0.34098647, -0.55129546, -0.37126957,  0.00436205])
m2=ax.errorbar(yrs-0.05,0.5*(up+lo),yerr=0.5*(up-lo),fmt=None,label='CESM-LME',ecolor='#003C71',elinewidth=2,capthick=2,zorder=1001)
#
#
yrs=np.arange(1258,1258+9)
up=np.array([-0.47525397,-0.08778505,0.48086479,0.34097769,0.52171893,0.41081254,0.61957062,0.73997964,0.80645339])
lo=np.array([-1.12296562,-0.58541222,-0.27965465,-0.20408419,-0.17367071,-0.21379843,0.05188414,-0.31109694,0.18331199 ])
m3=ax.errorbar(yrs+0.05,0.5*(up+lo),yerr=0.5*(up-lo),fmt=None,label='HI-SO2+LO-SO2',ecolor='#93328e',elinewidth=2,capthick=2,zorder=1001)
#
up=np.array([-0.14033673, -0.40392519, -0.16691037, -0.11573025 , 0.0240886 ,  0.01864214 ,0.26670808 , 0.36082187 , 0.72014242 ])
lo=np.array([-0.90095665, -0.6978948 , -0.7790608 , -0.64716563 ,-0.76453361, -0.69714966, -0.4804607 , -0.65716596 ,-0.11894622])
m4=ax.errorbar(yrs+0.15,0.5*(up+lo),yerr=0.5*(up-lo),fmt=None,label='HI-HAL',ecolor='#EF3340',elinewidth=2,capthick=2,zorder=1001)
#
fig.legend((m1,m2,m3,m4),('CMIP5','CESM-LME','BOTH-SO2','HI-HAL'),loc=4,bbox_to_anchor=(0.99,0.01),bbox_transform=ax.transAxes,handletextpad=0.2,columnspacing=0.8)
#
plt.xlim([1256.5,1265.5])	
plt.xticks([1257,1258,1259,1260,1261,1262,1263,1264,1265],['1257','1258','1259','1260','1261','1262','1263','1264','1265'])	
#plt.xlim([850,1900])	
plt.axhline(0.,c='k')
plt.ylim([-4.5,2.])
plt.xlabel('Calendar Year')
#legend = plt.legend(loc='lower right')
plt.ylabel(r'JJA 40-75$\degree$N Land Surface'+'\n'+r'Temperature Anomaly / $\degree$C')
if save==True:
                print "SAVING "+outdir+outname+'.'+format
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight",transparent=True)
plt.show()
