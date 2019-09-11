from netCDF4 import Dataset
import numpy as np
import math as m
def area():
        file=Dataset('/scratch/dcw32/netscratch/um/N96_area.nc')
	ar=file.variables['Area'][:]
	return ar
def gbox_areas(x,y):
        area=np.zeros([x,y])
        R=6.371E3
        for j in range(x):
                area[j,:]=(R**2)*m.radians(360./y)*(m.sin(m.radians(90.-(j-0.5)*180./(x-1)))-m.sin(m.radians(90.-(180./(x-1))*(j+0.5))))
        return area
    
