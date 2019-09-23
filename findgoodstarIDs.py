#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 13:58:45 2018

Code to generate a list of IDs for the stars that are used in determining PSFs

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots

### get data ###
sdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data

### Extract magnitudes ###
smag = vari_funcs.mag5_stacks(sdata)

### set those outside range to nan ###
smag[smag < 15] = np.nan
smag[smag > 19] = np.nan

### remove rows containing nans ###
mask = ~np.isnan(smag).any(axis=1)
smag = smag[mask]
sID = sdata['DR11_IDs']
goodID = np.copy(sID)
goodID = goodID[mask]

np.save('PSF_IDs_original',goodID)