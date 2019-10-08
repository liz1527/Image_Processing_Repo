#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:11:49 2019

code to investigate stars going into PSFs

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all')


def get_star_data_DR11(sdata, data):
    colname = 'MAG_APER'
    
    # limit magnitude range used in PSF
    mag = sdata[colname][:,4]
    mask1 = mag > 15
    mask2 = mag < 19
    mask3 = mag-sdata[colname][:,1] > -0.7 #pointlike-ness criterion
    mask4 = mag-sdata[colname][:,1] < -0.5 #pointlike-ness criterion
    mask = mask1*mask2*mask3*mask4
    
    tempsdata = sdata[mask]
    
    x = tempsdata[colname][:,4]
    y = tempsdata[colname][:,4] - tempsdata[colname][:,1]
    
    allx = data[colname][:,4]
    ally = data[colname][:,4] - data[colname][:,1]
    
    return x, y, allx, ally, tempsdata
    
sdata = fits.open('UDS_catalogues/DR11_output_stars.fits')[1].data
data = fits.open('UDS_catalogues/DR11_output.fits')[1].data
semesters = ['08B']#, '07B', '08B', '09B', '10B', '11B', '12B']


x, y, allx, ally, psf_star_data = get_star_data_DR11(sdata, data)


plt.figure()
plt.scatter(allx, ally, c='tab:grey',marker='+')
plt.scatter(x,y,c='b')
plt.xlabel('MAG_APER[4]')
plt.ylabel('MAG_APER[4] - MAG_APER[1]')
plt.xlim(xmax=21, xmin=10)
plt.ylim(ymax=-0.4, ymin=-1)
plt.title('DR11 K')
plt.tight_layout()

### Save star data for use in psfs ###
table_save = Table(psf_star_data)
table_save.write('UDS_catalogues/DR11_stars_for_PSFs.fits')