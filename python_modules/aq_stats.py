# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 15:35:36 2015

Funcions for calculating AQ stats,
comaring Model (Mod) to Observations (Obs)
based on equations:
http://www3.epa.gov/scram001/reports/O3%20NAAQS%20Final%20Rule%20RIA%20AQModeling%20Platform%20_03-11-08.pdf

Initial form of functions copied from Zhang et al., (2006).

All functions:
  * Mod = model output
  * Obs = observational data
  * N   = number of non-NaN pairs of data

IMPORTANT: Mod and Obs must be of same size & shape (paired data readings).

# =============================================================================

UPDATE: Jan 18 2015 - new functions to calculate unbiased and symetric AQ
                      evaluation statistics based on Yu et al. (2006).

* B_MNFB  = mean normalised factor bias
* E_MNAFE = mean narmalised absolute factor error
* B_NMBF  = normalised mean bias factor
* E_NMAEF = Normalized Mean Absolute Error Factor

Defined using 'S' for calculations where:
    S_i = (M_i - O_i)/|M_i - O_i|
    S = (M_avg - O_avg)/|M_avg - O_avg|

    * S_i == 1  if M_i >= O_i
    * S_i == -1 if O_i > M_i
'''

@author: scottarcher-nicholls
"""

# Load the modules we will need
import numpy as np

# Define list of statistics that can currently be calculated
AQ_stat_list = [
                'MB', 'NMB', 'RMSE', 'MAE', 'NME',  # Standard metrics
                'FB', 'FAE',                        # Fractional metric
                'NMFB', 'MNAFE', 'NMBF', 'NMAEF'    # Normalised factor metrics
                ]


# =============================================================================


def B_MB(Mod, Obs, check_nans=False):
    '''
    Function to calculate Mean Bias (in units of original)
       MB = SUM(M_i - O_i)/N
    '''
    # Check arrays are ame size
    assert len(Mod) == len(Obs)

    # cross-reference NaNs
    if(check_nans):
        Mod = np.where(np.isnan(Obs), np.nan, Mod)
        Obs = np.where(np.isnan(Mod), np.nan, Obs)

    # N = number of non-NaN entries
    N = np.sum(~np.isnan(Mod))

    return np.divide(np.nansum(np.subtract(Mod, Obs)), N)

# =============================================================================


def B_NMB(Mod, Obs, check_nans=False):
    '''
    Function to calculate Normalised Mean Bias
       NMB = SUM(M_i - O_i)/SUM(O_i)
    '''
    # Check arrays are ame size
    assert len(Mod) == len(Obs)

    # Cross-reference NaNs
    if(check_nans):
        Mod = np.where(np.isnan(Obs), np.nan, Mod)
        Obs = np.where(np.isnan(Mod), np.nan, Obs)

    # Calculate sum of all observations
    O_tot = np.nansum(Obs)

    return np.divide(np.nansum(np.subtract(Mod, Obs)), O_tot)  # *100.

# =============================================================================


def E_RMSE(Mod, Obs, check_nans=False):
    '''
    Function to calculate Root Mean Square Error (in units of original)
    RMSE = sqrt(SUM(M_i-O_i)^2/N)
    '''
    # assume two one-dimentional arrays
    assert len(Mod) == len(Obs)

    # Cross-reference NaNs
    if(check_nans):
        Mod = np.where(np.isnan(Obs), np.nan, Mod)
        Obs = np.where(np.isnan(Mod), np.nan, Obs)

    # N = number of non-NaN entries
    N = np.sum(~np.isnan(Mod))

    return np.sqrt(np.divide(np.nansum(np.square(np.subtract(Mod, Obs))), N))

# =============================================================================


def E_MAE(Mod, Obs, check_nans=False):
    '''
    Function to calculate Mean Absolute Arror (in units of original)
    MAE = SUM(|M_i-O_i|)/N)
    '''

    # assume two one-dimentional arrays
    assert len(Mod) == len(Obs)

    # cross-reference NaNs
    if(check_nans):
        Mod = np.where(np.isnan(Obs), np.nan, Mod)
        Obs = np.where(np.isnan(Mod), np.nan, Obs)


    N = np.sum(~np.isnan(Mod))

    return np.divide(np.nansum(np.abs(np.subtract(Mod, Obs))), N)

# =============================================================================


def E_NME(Mod, Obs, check_nans=False):
    '''
    Function to calculate Normalised Mean Error
    NME = SUM(|M_i-O_i|)/SUM(O_i)
    '''

    # assume two one-dimentional arrays 
    assert len(Mod) == len(Obs)

    # Pass NaN between equivalant array points
    if(check_nans):
        Mod   = np.where(np.isnan(Obs), np.nan, Mod)
        Obs   = np.where(np.isnan(Mod), np.nan, Obs)
    O_tot =  np.nansum(Obs)

    return np.divide(np.nansum(np.abs(np.subtract(Mod, Obs))),O_tot)  # *100.


# =============================================================================
#   Fractional based functions
# =============================================================================


def B_FB(Mod, Obs, check_nans=False):
    '''
    Fractional bias (varies -2 to 2):
    FB = 2/N * SUM( (M_i -O_i) / (M_i + O_i) )
    '''

    assert len(Mod) == len(Obs)

    # Pass NaN between equivalant array points
    if(check_nans):
        Mod = np.where(np.isnan(Obs), np.nan, Mod)
        Obs = np.where(np.isnan(Mod), np.nan, Obs)
    N   = np.sum(~np.isnan(Obs))

    return np.multiply(np.nansum(np.divide(np.subtract(Mod, Obs),
                                           np.add(Mod, Obs))), 2./N)

# =============================================================================


def E_FAE(Mod, Obs, check_nans=False):
    '''
    Fractional Absolute Error (varies 0 to 2):
    FAE = 2/N * SUM( |M_i - O_i| / (M_i + O_i) )
    '''

    # Pass NaN between equivalant array points
    assert len(Mod) == len(Obs)

    # Pass NaN between equivalant array points
    if(check_nans):
        Mod = np.where(np.isnan(Obs), np.nan, Mod)
        Obs = np.where(np.isnan(Mod), np.nan, Obs)
    N   = np.sum(~np.isnan(Obs))

    return np.multiply( np.nansum(np.divide(np.abs(np.subtract(Mod, Obs)),
                                            np.add(Mod, Obs))), 2./N)

# =============================================================================
#   Factor based functions
# =============================================================================


def B_NMFB(Mod, Obs, check_nans=False):
    '''
    Define Mean Normalised Factor Bias:
    NMBF = SUM(S_i * [exp(|ln(M_i/O_i)|)-1] ) / N
    '''

    # Pass NaN between equivalant array points
    assert len(Mod) == len(Obs)

    # Pass NaN between equivalant array points
    if(check_nans):
        Mod = np.where(np.isnan(Obs), np.nan, Mod)
        Obs = np.where(np.isnan(Mod), np.nan, Obs)

    N   = np.sum(~np.isnan(Mod))
    S   = np.where(Mod >= Obs, 1., -1.)

    # Calculate exponant part of function
    expo = np.exp(np.abs(np.log(np.divide(Mod, Obs))))

    return np.divide(np.nansum(np.multiply(S, (expo - 1))), N)

# =============================================================================


def E_MNAFE(Mod, Obs, check_nans=False):
    '''
    Define Mean normalised absolute factor error
    MNAFE = SUM( |exp(|ln(M/O)|) - 1| ) / N
    '''

    assert len(Mod) == len(Obs)

    # Pass NaN between equivalant array points
    # cross-reference NaNs
    if(check_nans):
        Mod = np.where(np.isnan(Obs), np.nan, Mod)
        Obs = np.where(np.isnan(Mod), np.nan, Obs)
    N = np.sum(~np.isnan(Mod))

    # Calculat exponant part of function
    expo = np.exp(np.abs(np.log(np.divide(Mod, Obs))))

    return np.divide(np.nansum(np.abs(expo - 1)), N)

# =============================================================================


def B_NMBF(Mod, Obs, check_nans=False):
    '''
    Define Normalised mean bias factor:
    NMBF = Sbar * [exp(|ln(Mbar/Obar)|) -1]
    '''

    assert len(Mod) == len(Obs)

    # Pass NaN between equivalant array points
    # cross-reference NaNs
    if(check_nans):
        Mod = np.where(np.isnan(Obs), np.nan, Mod)
        Obs = np.where(np.isnan(Mod), np.nan, Obs)

    # Calculate means of M and O
    Mbar = np.nanmean(Mod)
    Obar = np.nanmean(Obs)

    if (Mbar >= Obar):
        Sbar = 1.
    else:
        Sbar = -1.

    return Sbar * (np.exp(abs(np.log(Mbar / Obar))) - 1.)

# =============================================================================


def E_NMAEF(Mod, Obs, check_nans=False):
    '''
    Define Normalised mean absolute error factor:
    NMAEF = SUM(}M_i - O_i|) / (SUM(O_I)^(1+S/2))*(SUM(M_i)^(1-S)/2)
    '''

    assert len(Mod) == len(Obs)

    # Pass NaN between equivalant array points
    Mod = np.where(np.isnan(Obs), np.nan, Mod)
    Obs = np.where(np.isnan(Mod), np.nan, Obs)

    # Calculate means of M and O
    # cross-reference NaNs
    if(check_nans):
        Mtot = np.nansum(Mod)
        Otot = np.nansum(Obs)
    else:
        Mtot = np.sum(Mod)
	Otot = np.sum(Obs)

    if (Mtot >= Otot):
        Sbar = 1.
    else:
        Sbar = -1.

    return np.divide(np.nansum(np.abs(np.subtract(Mod, Obs))),
                     np.multiply(Otot**((1.+Sbar)/2.),
                                 Mtot**((1.-Sbar)/2.))
                     )

# =============================================================================

