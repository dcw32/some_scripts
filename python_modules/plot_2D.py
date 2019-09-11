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
def zonal_plot(hts,lats,field,levs,colmap,ymin,ymax,cbartitle,outdir,outname,save,format,extend):
	fig=plt.figure()
	plt.xlabel('Latitude')
	plt.ylabel('Height / km')
	#levi=np.copy(levs)
	#np.insert(levi,0,0.)
        cmap=plt.cm.get_cmap(colmap)
        levs2=np.insert(levs,0,levs[0]-1)
        levs2=np.append(levs2,levs2[len(levs2)-1]+1)
        norm=mpl.colors.BoundaryNorm(levs2, ncolors=cmap.N, clip=True)
	#cmap=plt.cm.get_cmap(colmap,len(levi))
	plt.contourf(lats,hts,field,levs,cmap=cmap,norm=norm,extend=extend)
	#plt.title(plot_title)
	plt.ylim([ymin,ymax])
	cbar=plt.colorbar()
	cbar.set_label(cbartitle)
        if save==True:
		print "SAVING "+outdir+outname+'.'+format
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
        plt.show()

