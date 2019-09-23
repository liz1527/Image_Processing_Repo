#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 16:38:35 2018

@author: ppxee
"""

### Import Modules Required ###
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.io import fits
#import sys
#import numpy as np

sem = '05B'#str(sys.argv[1]) #1st arguement from call is the semester
print(sem)

### read in fits_LDAC table and stars table ###
hdul = fits.open('star_stamps_tables/'+sem+'_star_stamps_table.fits', ignore_missing_end=True)
data = hdul[1].data
print('Read stamp table')
stars = Table.read('UDS_catalogues/DR11-secure-stars.fits')
print('Read stars table')

### Define coordinates ###
stampcoord = SkyCoord(data['ALPHA_J2000']*u.degree, data['DELTA_J2000']*u.degree)
starscoord = SkyCoord(stars['RA']*u.degree, stars['DEC']*u.degree)
print('Defined coordinates')

### Match catalogues and create new table ###
#idx, d2d , _ = match_coordinates_sky(starscoord, stampcoord) #match these 'good' stars to create table
idx = [2998]
starstamps = data[idx]

startb = fits.BinTableHDU.from_columns(starstamps)
startb.writeto('star_stamps_tables/test_'+sem+'_star_stamps_table.fits')


### old code ###
#data = Table.read('SE_stamp_outputs/small_'+sem+'_stamp_output.fits', ignore_missing_end=True)
#print('Read stamp table')
#stars = Table.read('UDS_catalogues/DR11-secure-stars.fits')
#print('Read stars table')
#
#### Define coordinates ###
#stampcoord = SkyCoord(data['ALPHA_J2000'], data['DELTA_J2000'])
#starscoord = SkyCoord(stars['RA']*u.degree, stars['DEC']*u.degree)
#print('Defined coordinates')
#
#### Match catalogues and create new table ###
#idx, d2d , _ = match_coordinates_sky(starscoord, stampcoord) #match these 'good' stars to create table
#starstamps = data[idx]
#
#starstamps.write('star_stamps_tables/small_'+sem+'_star_stamps_table.fits')

#-CATALOG_TYPE FITS_LDAC