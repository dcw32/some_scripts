from netCDF4 import Dataset
import numpy as np
import math as m
def gbox_areas(x,y):
# lats x lons
	area=np.zeros([x,y])
	R=6.371E6
	for j in range(x):
		area[j,:]=(R**2)*m.radians(360./y)*(m.sin(m.radians(90.-(j-0.5)*180./(x-1)))-m.sin(m.radians(90.-(180./(x-1))*(j+0.5))))
	return area
def sel_lats(x,y,lats,desc):
	box=np.ones([x,y])
	for j in range(x):
		if desc=='poles':
			if lats[j]<-66.4 or lats[j]>66.4:
				box[j,:]=1.
			else:
				box[j,:]=0.
                if desc=='nhp':
                        if lats[j]>66.4:
                                box[j,:]=1.
                        else:
                                box[j,:]=0.
		elif desc=='tropics':
			if abs(lats[j])<1.1:
				box[j,:]=1.
			else:
				box[j,:]=0.
	return box
def air_mass(file,mmr_air,Rgas):
	# mass of air = p*V*Mr(air) / R*T
	pres=file.variables['air_pressure_1'][:].squeeze()
	pres=np.mean(pres,axis=0)
	temp=file.variables['air_temperature_1'][:].squeeze()
	temp=np.mean(temp,axis=0)
	volfile=Dataset('/group_workspaces/jasmin2/ukca/vol1/dcw32/n48_l60_geovol.nc')
	vol=volfile.variables['vol_theta'][:].squeeze()
	mass=(pres*vol*mmr_air)/(Rgas*temp)
	print str(np.sum(mass))+' kg total mass of atmosphere'
	return mass
def air_mass2(file,mmr_air,Rgas):
        # mass of air = p*V*Mr(air) / R*T
        pres=file.variables['air_pressure'][:].squeeze()
        pres=np.mean(pres,axis=0)
        temp=file.variables['air_temperature'][:].squeeze()
        temp=np.mean(temp,axis=0)
        volfile=Dataset('/group_workspaces/jasmin2/ukca/vol1/dcw32/n48_l60_geovol.nc')
        vol=volfile.variables['vol_theta'][:].squeeze()
        mass=(pres*vol*mmr_air)/(Rgas*temp)
        print str(np.sum(mass))+' kg total mass of atmosphere'
        return mass
