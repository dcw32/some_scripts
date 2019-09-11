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
def profile(field,hts,col,xlab,ymin,ymax,xmin,xmax,outdir,outname,format,save):
	fig=plt.figure()
	plt.xlabel(xlab)
	plt.ylabel('Hybrid Height / km')
	plt.plot(field,hts,linewidth=2,color=col)
	axes=plt.gca()
	axes.set_ylim([ymin,ymax])
	axes.set_xlim([xmin,xmax])
        if save==True:
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
        plt.show()	
