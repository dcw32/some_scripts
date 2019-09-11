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
def latlon_shift(lons,lats,data):
	lonsi=lons-360
	data,lonsi=addcyclic(data,lonsi)
	data,lonsi=shiftgrid(-180.,data,lonsi,cyclic=360.)
	lonsi,latsi=np.meshgrid(lonsi,lats)
	return(lonsi,latsi,data)
def latlon_shift2(lons,data):
        lonsi=lons-360
        data,lonsi=addcyclic(data,lonsi)
        data,lonsi=shiftgrid(-180.,data,lonsi,cyclic=360.)
        return(lonsi,data)
def robin_had(data,lsm,lons,lats,levs,colmap,cbar_label,outdir,outname,save,format,extend):
        fig=plt.figure()
        ax1=fig.add_subplot(111)
        im0=Basemap(lon_0=0,projection='robin')
        a,b=im0(lons,lats)
        print a,b
        cmap=plt.cm.get_cmap(colmap)
        levs2=np.insert(levs,0,levs[0]-1)
        levs2=np.append(levs2,levs2[len(levs2)-1]+1)
        norm=mpl.colors.BoundaryNorm(levs2, ncolors=cmap.N, clip=True)
        #norm=mpl.colors.BoundaryNorm(levs,len(levs))
        #if extend=='max':      
#               cmap=plt.cm.get_cmap(colmap,len(levs)+1)
#       else:
#               cmap=plt.cm.get_cmap(colmap,len(levs)+2)
        cay=im0.contourf(a,b,data,levs,cmap=cmap,extend=extend,norm=norm)
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
	cby=im0.contour(a,b,lsm,1,colors='k')
        im0.drawmapboundary()
	cbaxes=fig.add_axes([0.15, 0.20, 0.7, 0.03])
        cbar=plt.colorbar(cay,cax=cbaxes,orientation="horizontal")
        cbar.set_label(cbar_label)
        if save==True:
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
                print "SAVING "+outdir+outname+'.'+format
        plt.show()
def robin_surf(data,lons,lats,levs,colmap,cbar_label,outdir,outname,save,format,extend):
	fig=plt.figure()
	ax1=fig.add_subplot(111)
	im0=Basemap(lon_0=0,projection='robin')
	a,b=im0(lons,lats)
	print a,b
	im0.drawcoastlines()
	im0.drawmapboundary()
	cmap=plt.cm.get_cmap(colmap)
	levs2=np.insert(levs,0,levs[0]-1)
	levs2=np.append(levs2,levs2[len(levs2)-1]+1)
	norm=mpl.colors.BoundaryNorm(levs2, ncolors=cmap.N, clip=True)
	#norm=mpl.colors.BoundaryNorm(levs,len(levs))
	#if extend=='max':	
#		cmap=plt.cm.get_cmap(colmap,len(levs)+1)
#	else:
#		cmap=plt.cm.get_cmap(colmap,len(levs)+2)
	cay=im0.contourf(a,b,data,levs,cmap=cmap,extend=extend,norm=norm)
	cbaxes=fig.add_axes([0.15, 0.20, 0.7, 0.03])
	cbar=plt.colorbar(cay,cax=cbaxes,orientation="horizontal")
	cbar.set_label(cbar_label)
	if save==True:
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
		print "SAVING "+outdir+outname+'.'+format
	plt.show()
def quiverplot(field,field2,lons,lats,colmap,scale,label,leng):
	mp=Basemap(lon_0=0,projection='robin')
	fig1=plt.figure()
	mp.drawcoastlines()
	uproj,vproj,xx,yy=mp.transform_vector(field,field2,lons,lats,31,31,returnxy=True,masked=True)
	lons,lats=np.meshgrid(lons,lats)
	a,b=mp(lons,lats)
	Q=mp.quiver(xx,yy,uproj,vproj,scale=scale,cmap=colmap)
	qk=plt.quiverkey(Q, 0.1, 0.05, leng, label, labelpos='W')
	plt.show()
def streamplot(field,field2,lons,lats,col):
        #mp=Basemap(lon_0=0,projection='robin')
	speed=np.sqrt(field*field+field2*field2)
	lw=5*speed/speed.max()
	mp = Basemap(projection='mill',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
        fig1=plt.figure()
	b,a=np.meshgrid(lons,lats)
	print b.shape
	print field.shape
	b,a=mp(b,a)
        mp.drawcoastlines()
        Q=mp.streamplot(b,a,field,field2,color=col,linewidth=lw)
        plt.show()
def robin_surf2(data,data_2,lons,lats,levs,colmap,cbar_label,outdir,outname,save,format,extend):
        fig=plt.figure()
	ax1=fig.add_subplot(121)
        im0=Basemap(lon_0=0,projection='robin')
	#a,b=np.meshgrid(lons,lats)
        #lons,lats=im0(lons,lats)
        a,b=im0(lons,lats)
        im0.drawcoastlines()
        im0.drawmapboundary()
	plt.text(0.05,1.0,'(a)',transform=ax1.transAxes)
	cmap=plt.cm.get_cmap(colmap)
        levs2=np.insert(levs,0,levs[0]-1)
        levs2=np.append(levs2,levs2[len(levs2)-1]+1)
        norm=mpl.colors.BoundaryNorm(levs2, ncolors=cmap.N, clip=True)
        #norm=mpl.colors.BoundaryNorm(levs,len(levs))
        #if extend=='max':      
#               cmap=plt.cm.get_cmap(colmap,len(levs)+1)
#       else:
#               cmap=plt.cm.get_cmap(colmap,len(levs)+2)
        cay=im0.contourf(a,b,data,levs,cmap=cmap,extend=extend,norm=norm)
	ax2=fig.add_subplot(122)
	plt.text(0.05,1.0,'(b)',transform=ax2.transAxes)
        im1=Basemap(lon_0=0,projection='robin')
        im1.drawcoastlines()
        im1.drawmapboundary()
	caz=im1.contourf(a,b,data_2,levs,cmap=cmap,extend=extend,norm=norm)
        cbaxes=fig.add_axes([0.25, 0.20, 0.5, 0.025])
        cbar=plt.colorbar(cay,cax=cbaxes,orientation="horizontal")
        cbar.set_label(cbar_label)
        if save==True:
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
                print "SAVING "+outdir+outname+'.'+format
        plt.show()
