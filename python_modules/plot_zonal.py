from netCDF4 import Dataset
import numpy as np
import calendar as calendar
from gridbox_areas import gbox_areas
import time
import matplotlib.pyplot as plt
import os
def zonal_pl(lats,field,col,lwid,xlab,xlim,ylab,ylim,aa,bb,outdir,outname,format,save):
	plt.plot(lats,field,color=col,linewidth=lwid)
	plt.xlabel(xlab)
	if aa==True:
		plt.xlim(xlim)
	plt.ylabel(ylab)
	if bb==True:
		plt.ylim(ylim)
        if save==True:
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
		print "SAVING "+outdir+outname+'.'+format
	plt.show()
def zonal_plr(lats,field,col,lwid,xlab,xlim,ylab,ylim,aa,bb,outdir,outname,format,save):
        plt.plot(field,lats,color=col,linewidth=lwid)
        plt.xlabel(ylab)
        if bb==True:
                plt.xlim(ylim)
        plt.ylabel(xlab)
        if aa==True:
                plt.ylim(xlim)
        if save==True:
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
                print "SAVING "+outdir+outname+'.'+format
        plt.show()
def zonal_plr_rhticks(lats,field,col,lwid,xlab,xlim,ylab,ylim,aa,bb,outdir,outname,format,save):
        f=plt.figure()
        ax=f.add_subplot(111)
        ax.yaxis.tick_right()
        plt.plot(field,lats,color=col,linewidth=lwid)
        plt.xlabel(ylab,fontsize=18)
        if bb==True:
                plt.xlim(ylim)
        ylab=plt.ylabel(xlab,fontsize=18)
	ax.yaxis.set_label_position("right")        
	ax.tick_params(axis='both',labelsize=16)
        if aa==True:
                plt.ylim(xlim)
        if save==True:
                if not os.path.exists(outdir):
                        os.makedirs(outdir)
                plt.savefig(outdir+outname+'.'+format,bbox_inches="tight")
                print "SAVING "+outdir+outname+'.'+format
        plt.show()
