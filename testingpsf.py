#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 10:44:49 2018

This is Will's PSF notebook put into code so that I can understand what it does

@author: ppxee
"""

from __future__ import print_function, division
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
import sys

sem = str(sys.argv[1]) #1st arguement from call is the semester

def norm(array):
    return array/np.nansum(array)

# reads in file with stamps
data = fits.open(sem+'_star_stamps_table.fits')[1].data

# extract stamps column and set non source regions to nans
stamps = data['VIGNET']
stamps[stamps<-5e29] = np.nan

# Plot some stars- not sure why, just for demonstration?
#plt.figure()
#plt.imshow(np.log(data['VIGNET'][24]))
#plt.figure()
#plt.imshow(np.log(stamps[1000,:,:]))

# Normalse all of the images - I presume to 1? Can't really see the difference...?
stampsnorm = stamps
for i in range(stamps.shape[0]):
    stampsnorm[i,:,:] = norm(stamps[i,:,:])

# find the median image
stack = np.nanmedian(stampsnorm, axis=0)

# print shape to check its right?
print(stack.shape)

#normalise and plot the median image
stack = norm(stack)
plt.figure()
plt.imshow(np.log(stack))
hdu = fits.PrimaryHDU(stack)
hdu.writeto(sem+'_K_PSF.fits', overwrite=True)

# This image is the typical PSF of the stack - can be used to create a matching kernel