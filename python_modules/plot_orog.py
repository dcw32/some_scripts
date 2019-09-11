from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
import calendar as calendar
import matplotlib as mpl
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math as m
import os
from mpl_toolkits.axes_grid1 import AxesGrid
def orog_plot(field,lons,lats,lat_0,levs,colmap,cbar_label,outdir,outname,save,format,extend):
	fig=plt.figure()
	im0=Basemap(lon_0=0,lat_0=lat_0,projection='ortho')
	a,b=im0(lons,lats)
	im0.drawcoastlines()
	im0.drawmapboundary()
	cmap=plt.cm.get_cmap(colmap)
	norm=mpl.colors.BoundaryNorm(levs, ncolors=cmap.N, clip=True)
	#norm=mpl.colors.BoundaryNorm(levs,len(levs))
	#if extend=='max':
	#	cmap=plt.cm.get_cmap(colmap,len(levs))
	#elif extend=='both':
	#	cmap=plt.cm.get_cmap(colmap,len(levs)+1)
	#else:
	#	cmap=plt.cm.get_cmap(colmap,len(levs)-1)
	cay=im0.contourf(a,b,field,levs,cmap=cmap,extend=extend)
	cbar=im0.colorbar(cay)
        cbar.set_label(cbar_label)
        if save==True:
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
	plt.show()
