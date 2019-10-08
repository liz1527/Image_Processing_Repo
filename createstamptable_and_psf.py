#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 10:48:02 2018

code to create star stamp table to detemine psf

@author: ppxee
"""

### Import Modules Required ###
from __future__ import print_function, division
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.io import fits
import sys
import numpy as np

sem = str(sys.argv[1]) #1st arguement from call is the semester
print(sem)

def norm(array):
    return array/np.nansum(array)
  
### read in fits_LDAC table and stars table ###
#keyword = {'ignore_missing_end': True}
data = Table.read('SE_stamp_outputs/cleaned_'+sem+'_H_stamp_output.fits', format='fits',hdu=2)
print('Read stamp table')
stars = Table.read('UDS_catalogues/DR11_stars_for_PSFs.fits')
print('Read stars table')

### Define coordinates ###
stampcoord = SkyCoord(data['ALPHA_J2000'], data['DELTA_J2000'])
starscoord = SkyCoord(stars['ALPHA_J2000_1']*u.degree, stars['DELTA_J2000_1']*u.degree)
print('Defined coordinates')

### Match catalogues and create new table ###
idx, d2d , _ = match_coordinates_sky(starscoord, stampcoord) #match these 'good' stars to create table
starstamps = data[idx]

starstamps.write('star_stamps_tables/cleaned_'+sem+'_H_star_stamps_PSF_table.fits')
#print(' finished')

print(len(starstamps))

# extract stamps column and set non source regions to nans
stamps = starstamps['VIGNET']
stamps[stamps<-5e29] = np.nan

# Normalse all of the images - I presume to 1? Can't really see the difference...?
stampsnorm = stamps
for i in range(stamps.shape[0]):
    stampsnorm[i,:,:] = norm(stamps[i,:,:])

# find the median image
stack = np.nanmedian(stampsnorm, axis=0)

# print shape to check its right?
print(stack.shape)

#normalise the median image
stack = norm(stack)

hdu = fits.PrimaryHDU(stack)
hdu.writeto('cleaned_Kstars_'+sem+'_H_PSF.fits', overwrite=True)

# This image is the typical PSF of the stack - can be used to create a matching kernel














