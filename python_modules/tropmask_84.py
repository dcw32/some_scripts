'''
Tropospheric mask: 1 for troposphere, 0 for statosphere
(Creates a python mask that can be used for any calculation)
Tropospheric height (STASH 30453)
UM: 2K.km-1 lapse rate threshold [WMO, 1057]

INPUT: 
ncfile	= file containing tropopause height
LTROP  	= 1 for troposphere (i.e. mask out stratosphere); 0 otherwise
tmean	= 0 (no time meaning),
	  1 (climatological monthly mean), 
	  2 (full time mean), 
	  in order of most to least expensive!
          5 (no time dimension in variable), 
nyrs	= length of run in years (0 if the entire data is required

OUTPUT:
mask
'''
#*****************************************************************************************************************
# Script to create masks
#*****************************************************************************************************************

import netCDF4 as ncdf
import scipy
import numpy as np
import numpy.ma as ma

#*****************************************************************************************************************
def mask(nctrop,ncpres,LTROP,tmean,nyrs):
    troppres = nctrop.variables['tropopause_pressure'][:]#[-12*nyrs:]   
    pres     = ncpres.variables['air_pressure'][:]#[-12*nyrs:]
#    pres     = np.mean(pres\
#        , axis = 0, dtype=np.float64)
#
    lat  = ncpres.variables['latitude'][:]
    lon  = ncpres.variables['longitude'][:]
    hgt  = range(len(ncpres.variables['level_height'][:]))
    time = ncpres.variables['time']#[-12*nyrs:]
#    time = N.zeros(pres.shape[0])
#
    if tmean==0:
      mask = np.zeros([len(time),len(hgt),len(lat),len(lon)],dtype=np.float64)
      # no time meaning
      for itme in range(len(time)):
          for ihgt in range(len(hgt)):
              for ilat in range(len(lat)):
                  for ilon in range(len(lon)):
                      if pres[itme,ihgt,ilat,ilon] >= troppres[itme,ilat,ilon]:
                          mask[itme,ihgt,ilat,ilon] = LTROP
    elif tmean==5:
      # no time dimension
      mask = np.zeros([len(hgt),len(lat),len(lon)],dtype=np.float64)
      for ilat in range(len(lat)):
         for ilon in range(len(lon)):
            for ihgt in range(len(hgt)):
                if pres[ihgt,ilat,ilon] >= troppres[ilat,ilon]:
                   mask[ihgt,ilat,ilon]  = LTROP
    mmask = np.ma.masked_where(mask==0,mask)
    return mmask




#    if tmean==0:
#    mask = np.zeros([len(time),len(hgt),len(lat),len(lon)],dtype=np.float64)
#      # no time meaning
#      for i in range(len(time)):
#          for j in range(len(hgt)):
#              for k in range(len(lat)):
#                  for l in range(len(lon)):
#                      if modhgt[0,j,k,l] <= trophgt[i,k,l]: 
#                          mask[i,j,k,l] = LTROP
#
#    elif tmean==1:
#      # climatological, monthly mean
#      trophgt_mean = np.empty([12,len(lat),len(lon)],dtype=np.float64)
#      for it in range(12):
#         jt = range(it,nyrs*12,12)
#         trophgt_mean[it,:,:] = np.mean(trophgt[jt,:,:],axis=0)
#      mask = np.empty([12,len(hgt),len(lat),len(lon)],dtype=np.float64)
#      for ilat in range(len(lat)):
#         for ilon in range(len(lon)):
#            for ihgt in range(len(hgt)):
#               for itim in range(12):
#                  if modhgt[0,ihgt,ilat,ilon]<=trophgt_mean[itim,ilat,ilon]:
#                     mask[itim,ihgt,ilat,ilon] = LTROP
#    elif tmean==2:
#      # full time mean
#      trophgt_mean = np.mean(trophgt,axis=0,dtype=np.float64)
#      mask = np.empty([len(hgt),len(lat),len(lon)],dtype=np.float64)
#      for ilat in range(len(lat)):
#         for ilon in range(len(lon)):
#            for ihgt in range(len(hgt)):
#                if modhgt[0,ihgt,ilat,ilon]<=trophgt_mean[ilat,ilon]:
#                   mask[ihgt,ilat,ilon] = LTROP
#    elif tmean==5:
#      # no time dimension
#      mask = np.empty([len(hgt),len(lat),len(lon)],dtype=np.float64)
#      for ilat in range(len(lat)):
#         for ilon in range(len(lon)):
#            for ihgt in range(len(hgt)):
#                if pres[0,ihgt,ilat,ilon]<=troppres[ilat,ilon]:
#                   mask[ihgt,ilat,ilon] = LTROP
#    mmask = np.ma.masked_where(mask==0,mask)
#    return mmask

#*****************************************************************************************************************

