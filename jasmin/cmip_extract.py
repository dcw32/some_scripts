from netCDF4 import Dataset
import numpy as np
from gridbox_areas import gbox_areas
import pylab as plt
import time
import sys
import csv
#This plots the Northern Hemisphere SAT over the 1257 volcanic eruption
#
giss1='/NASA-GISS/GISS-E2-R/past1000/mon/atmos/Amon/r1i1p121/files/tas_20120531/'
giss2='/NASA-GISS/GISS-E2-R/past1000/mon/atmos/Amon/r1i1p124/files/tas_20111208/'
giss3='/NASA-GISS/GISS-E2-R/past1000/mon/atmos/Amon/r1i1p127/files/tas_20120824/'
fgoals='/LASG-IAP/FGOALS-s2/past1000/mon/atmos/Amon/r1i1p1/files/tas_20130315/'
outloc='/group_workspaces/jasmin2/ukca/vol2/dcw32/CMIP5/lmil_40_75_land'
#files=['/BCC/bcc-csm1-1/past1000/mon/atmos/Amon/r1i1p1/files/tas_20120606/tas_Amon_bcc-csm1-1_past1000_r1i1p1_085001-185012.nc']
#files=['/NCAR/CCSM4/past1000/mon/atmos/Amon/r1i1p1/files/tas_20120604/tas_Amon_CCSM4_past1000_r1i1p1_085001-185012.nc']
#files=['/MPI-M/MPI-ESM-P/past1000/mon/atmos/Amon/r1i1p1/files/tas_20111028/tas_Amon_MPI-ESM-P_past1000_r1i1p1_085001-184912.nc']
#files=[['/IPSL/IPSL-CM5A-LR/past1000/mon/atmos/Amon/r1i1p1/files/tas_20120804/tas_Amon_IPSL-CM5A-LR_past1000_r1i1p1_085001-104912.nc','/IPSL/IPSL-CM5A-LR/past1000/mon/atmos/Amon/r1i1p1/files/tas_20120804/tas_Amon_IPSL-CM5A-LR_past1000_r1i1p1_105001-124912.nc','/IPSL/IPSL-CM5A-LR/past1000/mon/atmos/Amon/r1i1p1/files/tas_20120804/tas_Amon_IPSL-CM5A-LR_past1000_r1i1p1_125001-144912.nc','/IPSL/IPSL-CM5A-LR/past1000/mon/atmos/Amon/r1i1p1/files/tas_20120804/tas_Amon_IPSL-CM5A-LR_past1000_r1i1p1_145001-164912.nc','/IPSL/IPSL-CM5A-LR/past1000/mon/atmos/Amon/r1i1p1/files/tas_20120804/tas_Amon_IPSL-CM5A-LR_past1000_r1i1p1_165001-184912.nc']]
#files=[[giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_085001-089912.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_090001-094912.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_095001-099912.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_100001-105012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_105101-110012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_110101-115012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_115101-120012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_120101-125012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_125101-130012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_130101-135012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_135101-140012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_140101-145012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_145101-150012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_150101-155012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_155101-160012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_160101-165012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_165101-170012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_170101-175012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_175101-180012.nc',giss3+'tas_Amon_GISS-E2-R_past1000_r1i1p127_180101-185012.nc']]
#files=['/MIROC/MIROC-ESM/past1000/mon/atmos/Amon/r1i1p1/files/tas_20111228/tas_Amon_MIROC-ESM_past1000_r1i1p1_085001-184912.nc']
files=[[fgoals+'tas_Amon_FGOALS-s2_past1000_r1i1p1_085001-094912.nc',fgoals+'tas_Amon_FGOALS-s2_past1000_r1i1p1_095001-104912.nc',fgoals+'tas_Amon_FGOALS-s2_past1000_r1i1p1_105001-114912.nc',fgoals+'tas_Amon_FGOALS-s2_past1000_r1i1p1_115001-124912.nc',fgoals+'tas_Amon_FGOALS-s2_past1000_r1i1p1_125001-134912.nc',fgoals+'tas_Amon_FGOALS-s2_past1000_r1i1p1_135001-144912.nc',fgoals+'tas_Amon_FGOALS-s2_past1000_r1i1p1_145001-154912.nc',fgoals+'tas_Amon_FGOALS-s2_past1000_r1i1p1_155001-164912.nc',fgoals+'tas_Amon_FGOALS-s2_past1000_r1i1p1_165001-174912.nc',fgoals+'tas_Amon_FGOALS-s2_past1000_r1i1p1_175001-185012.nc']]
#lafs=['/BCC/bcc-csm1-1/piControl/fx/atmos/fx/r0i0p0/files/sftlf_20110101/sftlf_fx_bcc-csm1-1_piControl_r0i0p0.nc']
#lafs=['/NCAR/CCSM4/piControl/fx/atmos/fx/r0i0p0/files/sftlf_20120213/sftlf_fx_CCSM4_piControl_r0i0p0.nc']
#lafs=['/MPI-M/MPI-ESM-P/piControl/fx/atmos/fx/r0i0p0/files/sftlf_20120625/sftlf_fx_MPI-ESM-P_piControl_r0i0p0.nc']
#lafs=['/IPSL/IPSL-CM5A-LR/piControl/fx/atmos/fx/r0i0p0/files/sftlf_20110324/sftlf_fx_IPSL-CM5A-LR_piControl_r0i0p0.nc']
#lafs=['/NASA-GISS/GISS-E2-R/piControl/fx/atmos/fx/r0i0p0/files/sftlf_20130312/sftlf_fx_GISS-E2-R_piControl_r0i0p0.nc']
#lafs=['/MIROC/MIROC-ESM/piControl/fx/atmos/fx/r0i0p0/files/sftlf_20120828/sftlf_fx_MIROC-ESM_piControl_r0i0p0.nc']
lafs=['/LASG-IAP/FGOALS-s2/piControl/fx/atmos/fx/r0i0p0/files/sftlf_4/sftlf_fx_FGOALS-s2_piControl_r0i0p0.nc']
#models=['BCC']
#models=['CCSM4']
#models=['MPI']
#models=['IPSL']
#models=['GISS7']
#models=['MIROC']
models=['FGOALS']
nmodels=len(models)
#Note 2 of the GISS-ER do not include volcanic forcing, these are not used. An average of the remaining ensemble is used to determine the GISS average. The overall ensemble weights the 5 models equally
types=['mul']
for model in range(nmodels):
	print files[model]
	lsm_file=Dataset('/badc/cmip5/data/cmip5/output1'+lafs[model])
	lsm=lsm_file.variables['sftlf'][:].squeeze()
	#for segment in range()
	if types[model]=='sin':
		file=Dataset('/badc/cmip5/data/cmip5/output1'+files[model])
		time=file.variables['time'][:]
		lats=file.variables['lat'][:]
		lons=file.variables['lon'][:]
		var=file.variables['tas'][:].squeeze()
	elif types[model]=='mul':
		file=Dataset('/badc/cmip5/data/cmip5/output1'+files[model][0])
		lats=file.variables['lat'][:]
		lons=file.variables['lon'][:]
		time=[]
		var=np.empty((0,len(lats),len(lons)),float)
		print var.shape
		print len(files[model])
		for segs in range(len(files[model])):
			file=Dataset('/badc/cmip5/data/cmip5/output1'+files[model][segs])
			if models[model]=='FGOALS':
				time=np.append(time,(file.variables['time'][:]/365.)+850.+segs*100.)
			print file.variables['tas'][:].squeeze().shape
			var=np.append(var,file.variables['tas'][:].squeeze(),axis=0)
	print "A"
	#time=(time/365.)+850.
	print time
	nts=len(time)
	areas=gbox_areas(len(lats),len(lons))
	for lat in range(len(lats)):
		if lats[lat]<40.:
			if lats[lat]>75.:
				areas[lat,:]=0.
	print "B"
	areas=areas*lsm
	gmsat=np.zeros(nts)
	first_1800=0
	for ts in range(nts):
		gmsat[ts]=np.average(var[ts,:,:],weights=areas)
		if time[ts]<1800:
			first_1800=first_1800+1
	print first_1800
	print nts
	print "C"
	#norm to 1800-1850
	#deseasonalise
	for mons in range(12):
		gm_mon=np.mean(gmsat[mons::12])
		gmsat[mons::12]=gmsat[mons::12]-gm_mon
	print "D"
	gmsat_norm=gmsat-np.mean(gmsat[first_1800:first_1800+12*49])
	print "E"
	gmsat_mjja=np.zeros(nts/12)
	for yr in range(nts/12):
		gmsat_mjja[yr]=np.average(gmsat_norm[12*yr+4:12*yr+8])
	print "F"
	yrs_out=np.arange(850,850+nts/12)
	ofile=open(outloc+'/'+models[model]+'.csv',"wb")
	writer=csv.writer(ofile,delimiter=',')
	for yr in range(nts/12):
		writer.writerow([yrs_out[yr],gmsat_mjja[yr]])
	ofile.close()
