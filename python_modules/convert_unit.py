''' Convert mixing ratio units '''
#################################################################### 
#################################################################### 
import netCDF4 as ncdf
import numpy as np

import variables_attributes

def convert(unit,var):
    mmrtovmr=1/(variables_attributes.attributes(var)[4])
    molmass=variables_attributes.attributes(var)[-1]
    if unit=='ppmv':
        conversion = 1E6*mmrtovmr
    elif unit=='ppbv':
        conversion = 1E9*mmrtovmr
    elif unit=='pptv':
        conversion = 1E12*mmrtovmr
    elif unit=='h':
        conversion = (1/3600.0)		# from s
    elif unit=='days':
        conversion = (1/(3600.0*24))	# from s
    else:
        conversion=1
    return conversion

def convert2(unit,var,airmass):
    mmrtovmr=1/(variables_attributes.attributes(var)[4])
    molmass=variables_attributes.attributes(var)[-1]
    if unit=='ppmv':
        conversion = 1E6*mmrtovmr
    elif unit=='ppbv':
        conversion = 1E9*mmrtovmr
    elif unit=='pptv':
        conversion = 1E12*mmrtovmr
    elif unit==r'$\mathrm{molecules\ cm^{-3}}$': 
# OH(kg/kg) * mAir(kg) * Na(mol-1) / (Vgridcell(cm-3) * mMol(kg.mol-1) )
        volfile = ncdf.Dataset('/homes/pjt50/plotting/idl/newncfiles/ukca_geovol.nc','r')
        Vcell = np.array(volfile.variables['vol_theta'], dtype=np.float64)[:]
        conversion = airmass*6.022E23/(1E6*Vcell*molmass*1E-3)
    elif unit=='h':
        conversion = (1/3600.0)		# from s
    elif unit=='days':
        conversion = (1/(3600.0*24))	# from s
    else:
        conversion=1
    return conversion


