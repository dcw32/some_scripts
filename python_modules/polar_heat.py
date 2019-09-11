import math as m
import numpy as np
from scipy import interpolate
# Meridional Heat Transport
# MHT(lat)=-2*pi*R**2 *int|x=sin(lat),1 [ASR(x)-OLR(x)] dx
def mhtr_inter(field1,field2,field4,lats):
	olr=field4
	olr_fun=interpolate.interp1d(lats,olr,kind='cubic')
	asr=field1-field2
	asr_fun=interpolate.interp1d(lats,asr,kind='cubic')
	sens=180./360.
	lats_new=np.arange(-90.,90.+sens,sens)
	olr_new=olr_fun(lats_new)
	asr_new=asr_fun(lats_new)
	dx=np.zeros([len(lats_new)])
	ran=0.5*sens
	for i in range(len(lats_new)):
	        dx[i]=np.sin(np.radians(lats_new[i]+ran))-np.sin(np.radians(lats_new[i]-ran))
	mht=np.zeros([len(lats_new)])
	R=6.371E6 # m
	integrand_running=0
	for i in range(len(lats_new)-1,-1,-1):
	        integrand=dx[i]*(asr_new[i]-olr_new[i])
	        integrand_running=integrand_running+integrand
	        mht[i]=-2*m.pi*R**2*integrand_running
	mht=mht*1E-15
	return (mht,lats_new)
def mhtr(field1,field2,field4,lats):
        olr=field4
        asr=field1-field2
        sens=180./(len(lats)-1)
        dx=np.zeros([len(lats)])
        ran=0.5*sens
        for i in range(len(lats)):
                dx[i]=np.sin(np.radians(lats[i]+ran))-np.sin(np.radians(lats[i]-ran))
        mht=np.zeros([len(lats)])
        mht2=np.zeros([len(lats)])
        R=6.371E6 # m
        integrand_running=0
        for i in range(len(lats)-1,-1,-1):
                integrand=dx[i]*(asr[i]-olr[i])
                integrand_running=integrand_running+integrand
                mht[i]=-2*m.pi*R**2*integrand_running
        mht=mht*1E-15
        integrand_running=0
        for i in range(len(lats)):
                integrand=dx[i]*(asr[i]-olr[i])
                integrand_running=integrand_running+integrand
                mht2[i]=-2*m.pi*R**2*integrand_running
	mht2=mht2*1E-15
	mht=0.5*(mht-mht2)
        return mht
def mhtr_new(asr,olr,lats):
        sens=180./(len(lats)-1)
        dx=np.zeros([len(lats)])
        ran=0.5*sens
        for i in range(len(lats)):
                dx[i]=np.sin(np.radians(lats[i]+ran))-np.sin(np.radians(lats[i]-ran))
        mht=np.zeros([len(lats)])
        mht2=np.zeros([len(lats)])
        R=6.371E6 # m
        integrand_running=0
        for i in range(len(lats)-1,-1,-1):
                integrand=dx[i]*(asr[i]-olr[i])
                integrand_running=integrand_running+integrand
                mht[i]=-2*m.pi*R**2*integrand_running
        mht=mht*1E-15
        for i in range(len(lats)):
                integrand=dx[i]*(asr[i]-olr[i])
                integrand_running=integrand_running+integrand
                mht2[i]=-2*m.pi*R**2*integrand_running
        mht2=mht2*1E-15
        mht=0.5*(mht-mht2)
        return mht

