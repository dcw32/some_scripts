''' Assign information to STASH numbers.
Such as variable names, units, conversion factors.
'''
#################################################################### 
# ASSIGN INFORMATION TO STASH NUMBERS

# Standard names than iris reads can be in file in "python:> help(iris.std_names)": /usr/local/shared/ubuntu-12.04/x86_64/python2.7-iris/1.6.1/local/lib/python2.7/site-packages/Iris-1.6.1-py2.7.egg/iris/std_names.py

#-------------------------------------------------------------------

# conversion factor= cspecies=M(var)/M(air)=M(var)/28.97g.mol-1
#################################################################### 

def attributes(var):
    if var=='air_temperature':
        long_name='T on theta levels in K'
        name=r'Temperature'
        outputname=r'T'
        formula=r'$\mathbf{T}$'
        unit=r'K'
        conversion_factor=1
        molmass=1
    elif var=='air_pressure':
        long_name='P on theta levels in Pa'
        name=r'Pressure'
        outputname=r'P'
        formula=r'$\mathbf{P}$'
        unit='Pa'
        conversion_factor=1
        molmass=1
    elif var=='age_of_air':
        long_name='Age of air in s'
        name=r'Age of air'
        outputname=r'Age_of_air'
        formula=r'$\mathbf{Age\ of\ air}$'
        unit='s'
        conversion_factor=1
        molmass=1
    elif var=='specific_humidity':
        long_name='Specific humidity in kg(H2O)/kg(air)'
        name=r'Specific humidity'
        outputname=r'Spec_humidity'
        formula=r'$\mathbf{H_{2}O}$'
        unit='1'
        conversion_factor=0.621
        molmass=18
    elif var=='tropopause_pressure':
        long_name='P at Tropopause Level in Pa'
        name=r'Pressure at Tropopause'
        outputname=r'P_tropopause'
        formula=r'$\mathbf{P_{Tropopause}}$'
        unit='Pa'
        conversion_factor=1
        molmass=1
    elif var=='tropopause_temperature':
        long_name='T at Tropopause Level'
        name=r'Temperature at Tropopause'
        outputname=r'T_tropopause'
        formula=r'$\mathbf{T_{Tropopause}}$'
        unit='K'
        conversion_factor=1
        molmass=1
    elif var=='tropopause_height':
        long_name='z at Tropopause Level'
        name=r'Height at Tropopause'
        outputname=r'Height_tropopause'
        formula=r'$\mathbf{z_{Tropopause}}$'
        unit='m'
        conversion_factor=1
        molmass=1
#--------------------------------------------------------------
# Tracers UKCA
#--------------------------------------------------------------
    elif var=='mass_fraction_of_ozone_in_air':
        long_name='O3 MMR in kg(O3)/kg(air)'
        name=r'Ozone'
        outputname=r'O3'
        formula=r'$\mathbf{O_{3}}$'
        unit='1'
        conversion_factor=1.657
        molmass=48
    elif var=='mass_fraction_of_nitrogen_monoxide_in_air':
        long_name='NO MMR in kg(NO)/kg(air)'
        name=r'Nitrogen monoxide'
        outputname=r'NO'
        formula=r'$\mathbf{NO}$'
        unit='1'
        conversion_factor=1.036
        molmass=30
    elif var=='mass_fraction_of_NO3_in_air':
        long_name='NO3 MMR in kg(NO3)/kg(air)'
        name=r'NO3'
        outputname=r'NO3'
        formula=r'$\mathbf{NO_{3}}$'
        unit='1'
        conversion_factor=2.14
        molmass=62
    elif var=='mass_fraction_of_methane_in_air':
        long_name='CH4 MMR in kg(CH4)/kg(air)'
        name=r'Methane'
        outputname=r'CH4'
        formula=r'$\mathbf{CH_{4}}$'
        unit='1'
        conversion_factor=0.552
        molmass=16
    elif var=='mass_fraction_of_carbon_monoxide_in_air':
        long_name='CO MMR in kg(CO)/kg(air)'
        name=r'Carbon monoxide'
        outputname=r'CO'
        formula=r'$\mathbf{CO}}$'
        unit='1'
        conversion_factor=0.967
        molmass=28
    elif var=='mass_fraction_of_peroxyacetyl_nitrate_in_air':
        long_name='PAN MMR in kg(PAN)/kg(air)'
        name=r'PAN'
        outputname=r'PAN'
        formula=r'$\mathbf{PAN}}$'
        unit='1'
        conversion_factor=4.177
        molmass=121
    elif var=='mass_fraction_of_chlorine_in_air':
        long_name='Cl MMR in kg(Cl)/kg(air)'
        name=r'Chlorine'
        outputname=r'Cl'
        formula=r'$\mathbf{Cl}$'
        unit='1'
        conversion_factor=1.225
        molmass=35.5
    elif var=='mass_fraction_of_chlorine_monoxide_in_air':
        long_name='ClO MMR in kg(ClO)/kg(air)'
        name=r'Chlorine monoxide'
        outputname=r'ClO'
        formula=r'$\mathbf{ClO}$'
        unit='1'
        conversion_factor=1.778
        molmass=51.5
    elif var=='mass_fraction_of_dichlorine_peroxide_in_air':
        long_name='Cl2O2 MMR in kg(Cl2O2)/kg(air)'
        name=r'Dichlorine peroxide'
        outputname=r'Cl2O2'
        formula=r'$\mathbf{Cl_{2}O_{2}}$'
        unit='1'
        conversion_factor=3.555
        molmass=103
#    elif var=='mass_fraction_of_chlorine_dioxide_in_air':
#        cube.long_name='OClO MMR in kg(OClO)/kg(air)'
#        cube.attributes['name']=r'Chlorine dioxide'
#        cube.attributes['outputname']=r'ClO2'
#        cube.attributes['formula']=r'$\mathbf{OClO}}$'
#        cube.attributes['unit']='1'
#        cube.attributes['conversion_factor']=2.330
#        cube.attributes['molmass']=67.5
    elif var=='mass_fraction_of_nitrous_oxide_in_air':
        long_name='N2O MMR in kg(N2O)/kg(air)'
        name=r'Nitrous oxide'
        outputname=r'N2O'
        formula=r'$\mathbf{N_{2}O}$'
        unit='1'
        conversion_factor=1.519
        molmass=44
    elif var=='mass_fraction_of_hydroxyl_radical_in_air':
        long_name='OH MMR in kg(OH)/kg(air)'
        name=r'Hydroxyl radical'
        outputname=r'OH'
        formula=r'$\mathbf{OH}$'
        unit='1'
        conversion_factor=0.587
        molmass=17
#    elif var=='mass_fraction_of_hydroperoxyl_radical_in_air':
#        cube.long_name='HO2 MMR in kg(HO2)/kg(air)'
#        cube.attributes['name']=r'Hydroperoxyl radical'
#        cube.attributes['outputname']=r'HO2'
#        cube.attributes['formula']=r'$\mathbf{HO_{2}}$'
#        cube.attributes['unit']='1'
#        cube.attributes['conversion_factor']=1.139
#        cube.attributes['molmass']=33
#    elif var=='mass_fraction_of_passive_ozone_in_air':
#        cube.long_name='passive O3 MMR in kg(O3)/kg(air)'
#        cube.attributes['name']=r'Passive ozone'
#        cube.attributes['outputname']=r'O3_passive'
#        cube.attributes['formula']=r'$\mathbf{O_{3}}$'
#        cube.attributes['unit']='1'
#        cube.attributes['conversion_factor']=1.657
#        cube.attributes['molmass']=48
#    elif var=='mass_fraction_of_singlett_oxygen_in_air':
#        cube.long_name='O(1D) MMR in kg(O1D)/kg(air)'
#        cube.attributes['name']=r'Singlett Oxygen'
#        cube.attributes['outputname']=r'O1D'
#        cube.attributes['formula']=r'$\mathbf{O(^{1}D)}$'
#        cube.attributes['unit']='1'
#        cube.attributes['conversion_factor']=0.552
#        cube.attributes['molmass']=16
#    elif var=='mass_fraction_of_nitrogen_dioxide_in_air':
#        cube.long_name='NO2 MMR in kg(NO2)/kg(air)'
#        cube.attributes['name']=r'Nitrogen dioxide'
#        cube.attributes['outputname']=r'NO2'
#        cube.attributes['formula']=r'$\mathbf{NO_{2}}$'
#        cube.attributes['unit']='1'
#        cube.attributes['conversion_factor']=1.588
#        cube.attributes['molmass']=46
#--------------------------------------------------------------
#?# unit Dobson? conversion factor
#    elif var=='Ozone_column':
#        cube.long_name='O3 column in Dobson'
#        cube.attributes['name']=r'Ozone column'
#        cube.attributes['outputname']=r'O3_column'
#        cube.attributes['formula']=r'$\mathbf{O_{3}^{column}}}$'
#        cube.attributes['unit']='Dobson'
#        cube.attributes['conversion_factor']=1
#        cube.attributes['molmass']=48
#--------------------------------------------------------------
# Reaction fluxes
#--------------------------------------------------------------
    elif var=='ch4_oh_rxn_flux':
        long_name='CH4 + OH reaction flux in mol/gridcell/s'
        name=r'CH4+OH'
        outputname=r'CH4+OH'
        formula=r'$\mathbf{CH_{4}\ +\ OH}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
    elif var=='ox_loss_O1D_H2O':
        long_name='O1D + H2O reaction flux in mol/gridcell/s'
        name=r'O1D+H2O'
	outputname=name
        formula=r'$\mathbf{O(1D)\ +\ H_{2}O}$'
	unit='mol gridcell-1 s-1'
	conversion_factor=1
	molmass=1
#    elif var=='STE_ozone':
#        cube.long_name='Stratosphere-Troposphere-Exchange of ozone'
#        cube.attributes['name']=r'STE Ozone'
#        cube.attributes['outputname']=r'O3_STE'
#        cube.attributes['formula']=r'$\mathbf{STE\ O_{3}}$'
#        cube.attributes['unit']='mol gridcell-1 s-1'
#        cube.attributes['conversion_factor']=1
#        cube.attributes['molmass']=48
#    elif var=='air_mass_trop':
#        cube.long_name='Air mass of troposphere'
#        cube.attributes['name']=r'Air mass troposphere'
#        cube.attributes['outputname']=r'Air_mass_tropos'
#        cube.attributes['formula']=r'$\mathbf{Air\ mass^{troposphere}}$'
#        cube.attributes['unit']='kg gridcell-1'
#        cube.attributes['conversion_factor']=1
#        cube.attributes['molmass']=1
#    elif var=='air_mass_atm':
#        cube.long_name='Air mass of whole atmosphere'
#        cube.attributes['name']=r'Air mass atmosphere'
#        cube.attributes['outputname']=r'Air_mass_atmos'
#        cube.attributes['formula']=r'$\mathbf{Air\ mass^{atmosphere}}$'
#        cube.attributes['unit']='kg gridcell-1'
#        cube.attributes['conversion_factor']=1
#        cube.attributes['molmass']=1
#    elif var=='tropospheric_mask':
#        cube.long_name='Tropospheric mask'
#        cube.attributes['name']=r'Tropospheric mask'
#        cube.attributes['outputname']=r'Tropos_Mask'
#        cube.attributes['formula']=r'$\mathbf{Tropospheric mask}$'
#        cube.attributes['unit']='1:tropos,0:other'
#        cube.attributes['conversion_factor']=1
#        cube.attributes['molmass']=1
#--------------------------------------------------------------
# Emissions
#--------------------------------------------------------------
    elif var=='ch4_apparent_ems':
        long_name='CH4 forced surf ems flux in mol/gridcell/s'
        name=r'CH4 surf emissions forced'
        outputname=r'CH4_surfems_forced'
        formula=r'$\mathbf{CH_{4}\ surface\ emissions^{forced}}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=16
    elif var=='ch4_ems':
        long_name='CH4 ems flux in kg/m-2/s'
        name=r'CH4 surf emissions'
        outputname=r'CH4_surfems'
        formula=r'$\mathbf{CH_{4}\ surface\ emissions}$'
        unit='kg m-2 s-1'
        conversion_factor=1
        molmass=16
#--------------------------------------------------------------
# Offline oxidants
#--------------------------------------------------------------
    elif var=='mmr_CH4anthropogenic':
        long_name='CH4 anthropogenic MMR in kg(CH4)/kg(air)'
        name=r'Anthropogenic methane'
        outputname=r'CH4anth'
        formula=r'$\mathbf{CH_{4}\ anthropogenic}$'
        unit='1'
        conversion_factor=0.552
        molmass=16
    elif var=='flux_OHfix_CH4anth':
        long_name='CH4anth + OHfix reaction flux in mol/gridcell/s'
        name=r'CH4anth+OHfix'
        outputname=r'CH4anth+OHfix'
        formula=r'$\mathbf{CH_{4}\ anthropogenic\ +\ OH\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
    elif var=='flux_O1Dfix_CH4anth':
        long_name='CH4anth + O1Dfix reaction flux in mol/gridcell/s'
        name=r'CH4anth+O1Dfix'
        outputname=r'CH4anth+O1Dfix'
        formula=r'$\mathbf{CH_{4}\ anthropogenic\ +\ O^{1}D\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
    elif var=='flux_Clfix_CH4anth':
        long_name='CH4anth + Clfix reaction flux in mol/gridcell/s'
        name=r'CH4anth+Clfix'
        outputname=r'CH4anth+Clfix'
        formula=r'$\mathbf{CH_{4}\ anthropogenic\ +\ Cl\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
#
    elif var=='mmr_CH4wetlands':
        long_name='CH4 wetlands MMR in kg(CH4)/kg(air)'
        name=r'Wetland methane'
        outputname=r'CH4wetl'
        formula=r'$\mathbf{CH_{4}\ wetlands}$'
        unit='1'
        conversion_factor=0.552
        molmass=16
    elif var=='flux_OHfix_CH4wetl':
        long_name='CH4wetl + OHfix reaction flux in mol/gridcell/s'
        name=r'CH4wetl+OHfix'
        outputname=r'CH4wetl+OHfix'
        formula=r'$\mathbf{CH_{4}\ wetlands\ +\ OH\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
    elif var=='flux_O1Dfix_CH4wetl':
        long_name='CH4wetl + O1Dfix reaction flux in mol/gridcell/s'
        name=r'CH4wetl+O1Dfix'
        outputname=r'CH4wetl+O1Dfix'
        formula=r'$\mathbf{CH_{4}\ wetlands\ +\ O^{1}D\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
    elif var=='flux_Clfix_CH4wetl':
        long_name='CH4wetl + Clfix reaction flux in mol/gridcell/s'
        name=r'CH4wetl+Clfix'
        outputname=r'CH4wetl+Clfix'
        formula=r'$\mathbf{CH_{4}\ wetlands\ +\ Cl\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
#
    elif var=='mmr_CH4termites':
        long_name='CH4 termites MMR in kg(CH4)/kg(air)'
        name=r'Termite methane'
        outputname=r'CH4term'
        formula=r'$\mathbf{CH_{4}\ termites}$'
        unit='1'
        conversion_factor=0.552
        molmass=16
    elif var=='flux_OHfix_CH4term':
        long_name='CHterm + OHfix reaction flux in mol/gridcell/s'
        name=r'CH4term+OHfix'
        outputname=r'CH4term+OHfix'
        formula=r'$\mathbf{CH_{4}\ termites\ +\ OH\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
    elif var=='flux_O1Dfix_CH4term':
        long_name='CH4term + O1Dfix reaction flux in mol/gridcell/s'
        name=r'CH4term+O1Dfix'
        outputname=r'CH4term+O1Dfix'
        formula=r'$\mathbf{CH_{4}\ termites\ +\ O^{1}D\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
    elif var=='flux_Clfix_CH4term':
        long_name='CH4term + Clfix reaction flux in mol/gridcell/s'
        name=r'CH4term+Clfix'
        outputname=r'CH4term+Clfix'
        formula=r'$\mathbf{CH_{4}\ termites\ +\ Cl\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
#
    elif var=='mmr_CH4hydrates':
        long_name='CH4 hydr MMR in kg(CH4)/kg(air)'
        name=r'Hydrate methane'
        outputname=r'CH4hydr'
        formula=r'$\mathbf{CH_{4}\ hydrates}$'
        unit='1'
        conversion_factor=0.552
        molmass=16
    elif var=='flux_OHfix_CH4hydr':
        long_name='CH4hydr + OHfix reaction flux in mol/gridcell/s'
        name=r'CH4hydr+OHfix'
        outputname=r'CH4hydr+OHfix'
        formula=r'$\mathbf{CH_{4}\ hydrates\ +\ OH\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
    elif var=='flux_O1Dfix_CH4hydr':
        long_name='CH4hydr + O1Dfix reaction flux in mol/gridcell/s'
        name=r'CH4hydr+O1Dfix'
        outputname=r'CH4hydr+O1Dfix'
        formula=r'$\mathbf{CH_{4}\ hydrates\ +\ O^{1}D\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
    elif var=='flux_Clfix_CH4hydr':
        long_name='CH4hydr + Clfix reaction flux in mol/gridcell/s'
        name=r'CH4hydr+Clfix'
        outputname=r'CH4hydr+Clfix'
        formula=r'$\mathbf{CH_{4}\ hydrates\ +\ Cl\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
#
    elif var=='mmr_CH4soilloss':
        long_name='CH4 soil loss MMR in kg(CH4)/kg(air)'
        name=r'Soil loss methane'
        outputname=r'CH4soil'
        formula=r'$\mathbf{CH_{4}\ soil loss}$'
        unit='1'
        conversion_factor=0.552
        molmass=16
    elif var=='flux_OHfix_CH4soil':
        long_name='CH4soil + OHfix reaction flux in mol/gridcell/s'
        name=r'CH4soil+OHfix'
        outputname=r'CH4soil+OHfix'
        formula=r'$\mathbf{CH_{4}\ soil\ loss\ +\ OH\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
    elif var=='flux_O1Dfix_CH4soil':
        long_name='CH4soil + O1Dfix reaction flux in mol/gridcell/s'
        name=r'CH4soil+O1Dfix'
        outputname=r'CH4soil+O1Dfix'
        formula=r'$\mathbf{CH_{4}\ soil\ loss\ +\ O^{1}D\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
    elif var=='flux_Clfix_CH4soil':
        long_name='CH4soil + Clfix reaction flux in mol/gridcell/s'
        name=r'CH4soil+Clfix'
        outputname=r'CH4soil+Clfix'
        formula=r'$\mathbf{CH_{4}\ soil\ loss\ +\ Cl\ fixed}$'
        unit='mol gridcell-1 s-1'
        conversion_factor=1
        molmass=1
#--------------------------------------------------------------
    return name, outputname, formula, unit, conversion_factor, molmass


