''' Assign information to STASH numbers.
Such as variable names, units, conversion factors.
'''
#################################################################### 
# ASSIGN INFORMATION TO STASH NUMBERS

# Standard names than iris reads can be in file in "python:> help(iris.std_names)": /usr/local/shared/ubuntu-12.04/x86_64/python2.7-iris/1.6.1/local/lib/python2.7/site-packages/Iris-1.6.1-py2.7.egg/iris/std_names.py

#-------------------------------------------------------------------

# conversion factor= cspecies=M(var)/M(air)=M(var)/28.97g.mol-1
#################################################################### 

def UKCA_callback(cube, field, filename):
    if cube.attributes['STASH'] == 'm01s16i004':
        cube.standard_name='air_temperature'

    if cube.attributes['STASH'] == 'm01s00i408':
        cube.standard_name='air_pressure'

    if cube.attributes['STASH'] == 'm01s34i150':
        cube.var_name='age_of_air'

    if cube.attributes['STASH'] == 'm01s00i010':
        cube.standard_name='specific_humidity'

    if cube.attributes['STASH'] == 'm01s30i451':
        cube.var_name='tropopause_pressure'

    if cube.attributes['STASH'] == 'm01s30i452':
        cube.var_name='tropopause_temperature'

    if cube.attributes['STASH'] == 'm01s30i453':
        cube.var_name='tropopause_height'
#--------------------------------------------------------------
# Tracers UKCA
#--------------------------------------------------------------
    if cube.attributes['STASH'] == 'm01s34i001':
        cube.standard_name='mass_fraction_of_ozone_in_air'

    if cube.attributes['STASH'] == 'm01s34i002':
        cube.standard_name='mass_fraction_of_nitrogen_monoxide_in_air'

    if cube.attributes['STASH'] == 'm01s34i003':
        cube.var_name='mass_fraction_of_NO3_in_air'

    if cube.attributes['STASH'] == 'm01s34i009':
        cube.standard_name='mass_fraction_of_methane_in_air'

    if cube.attributes['STASH'] == 'm01s34i010':
        cube.standard_name='mass_fraction_of_carbon_monoxide_in_air'

    if cube.attributes['STASH'] == 'm01s34i017':
        cube.standard_name='mass_fraction_of_peroxyacetyl_nitrate_in_air'

    if cube.attributes['STASH'] == 'm01s34i041':
        cube.var_name='mass_fraction_of_chlorine_in_air'

    if cube.attributes['STASH'] == 'm01s34i042':
        cube.standard_name='mass_fraction_of_chlorine_monoxide_in_air'

    if cube.attributes['STASH'] == 'm01s34i043':
        cube.standard_name='mass_fraction_of_dichlorine_peroxide_in_air'

    if cube.attributes['STASH'] == 'm01s34i044':
        cube.standard_name='mass_fraction_of_chlorine_dioxide_in_air'

    if cube.attributes['STASH'] == 'm01s34i049':
        cube.standard_name='mass_fraction_of_nitrous_oxide_in_air'

    if cube.attributes['STASH'] == 'm01s34i081':
        cube.standard_name='mass_fraction_of_hydroxyl_radical_in_air'

    if cube.attributes['STASH'] == 'm01s34i082':
        cube.standard_name='mass_fraction_of_hydroperoxyl_radical_in_air'

    if cube.attributes['STASH'] == 'm01s34i149':
        cube.var_name='mass_fraction_of_passive_ozone_in_air'

    if cube.attributes['STASH'] == 'm01s34i151':
        cube.var_name='mass_fraction_of_singlett_oxygen_in_air'

    if cube.attributes['STASH'] == 'm01s34i152':
        cube.standard_name='mass_fraction_of_nitrogen_dioxide_in_air'
#--------------------------------------------------------------
#?# unit Dobson? conversion factor
    if cube.attributes['STASH'] == 'm01s34i172':
        cube.var_name='Ozone_column'
#--------------------------------------------------------------
# Reaction fluxes
#--------------------------------------------------------------
    if cube.attributes['STASH'] == 'm01s34i301':
        cube.var_name='ox_prod_HO2_NO'

    if cube.attributes['STASH'] == 'm01s34i302':
        cube.var_name='ox_prod_MeOO_NO'

    if cube.attributes['STASH'] == 'm01s34i303':
        cube.var_name='ox_prod_NO_RO2'

    if cube.attributes['STASH'] == 'm01s34i304':
        cube.var_name='ox_prod_OH_inorgAcid'

    if cube.attributes['STASH'] == 'm01s34i305':
        cube.var_name='ox_prod_OH_orgNitrate'

    if cube.attributes['STASH'] == 'm01s34i306':
        cube.var_name='ox_prod_orgNitrate_photol'

    if cube.attributes['STASH'] == 'm01s34i307':
        cube.var_name='ox_prod_OH_PANrxns'

    if cube.attributes['STASH'] == 'm01s34i311':
        cube.var_name='ox_loss_O1D_H2O'

    if cube.attributes['STASH'] == 'm01s34i312':
        cube.var_name='ox_loss_minor_rxns'

    if cube.attributes['STASH'] == 'm01s34i313':
        cube.var_name='ox_loss_HO2_O3'

    if cube.attributes['STASH'] == 'm01s34i314':
        cube.var_name='ox_loss_OH_O3'

    if cube.attributes['STASH'] == 'm01s34i315':
        cube.var_name='ox_loss_O3_alkene'

    if cube.attributes['STASH'] == 'm01s34i316':
        cube.var_name='ox_loss_N2O5_H2O'

    if cube.attributes['STASH'] == 'm01s34i317':
        cube.var_name='ox_loss_NO3_chemloss'

    if cube.attributes['STASH'] == 'm01s34i321':
        cube.var_name='ozone_dry_dep_3D'

    if cube.attributes['STASH'] == 'm01s34i322':
        cube.var_name='noy_dry_dep_3D'

    if cube.attributes['STASH'] == 'm01s34i331':
        cube.var_name='noy_wet_dep_3D'

    if cube.attributes['STASH'] == 'm01s34i341':
        cube.var_name='ch4_oh_rxn_flux'

    if cube.attributes['STASH'] == 'm01s34i351':
        cube.var_name='STE_ozone'

    if cube.attributes['STASH'] == 'm01s34i352':
        cube.var_name='tendency_ozone_troposphere'

    if cube.attributes['STASH'] == 'm01s34i354':
        cube.var_name='tendency_ozone_atm'
# ?unit
    if cube.attributes['STASH'] == 'm01s34i353':
        cube.var_name='tropos_ozone'

    if cube.attributes['STASH'] == 'm01s34i361':
        cube.var_name='air_mass_trop'

    if cube.attributes['STASH'] == 'm01s34i363':
        cube.var_name='air_mass_atm'

    if cube.attributes['STASH'] == 'm01s34i362':
        cube.var_name='tropospheric_mask'
#--------------------------------------------------------------
# Emissions
#--------------------------------------------------------------
    if cube.attributes['STASH'] == 'm01s34i451':
        cube.var_name='ch4_apparent_ems'

    if cube.attributes['STASH'] == 'm01s00i302':
        cube.var_name='ch4_ems'




