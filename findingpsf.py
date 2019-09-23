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

print(sem)

def norm(array):
    return array/np.nansum(array)

def quadrants(basedata, initdata):
    
    bra = basedata['ALPHA_J2000_05B']
    bdec = basedata['DELTA_J2000_05B']
    ira = initdata['ALPHA_J2000_05B']
    idec = initdata['DELTA_J2000_05B']

    ### define bounds of quadrants ###
    maxra = np.max(bra)
    minra = np.min(bra)
    maxdec = np.max(bdec)
    mindec = np.min(bdec)
    
    midra = (maxra+minra)/2
    middec = (maxdec+mindec)/2
    
    ### create masks for quadrant ###
    mask1 = ira < midra
    mask2 = idec < middec
    quad1data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec < middec
    quad2data = initdata[mask1*mask2]
    
    mask1 = ira < midra
    mask2 = idec >= middec
    quad3data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec >= middec
    quad4data = initdata[mask1*mask2]
    
    return quad1data, quad2data, quad3data, quad4data

# reads in file with stamps
data = fits.open('star_stamps_tables/small_'+sem+'_star_stamps_table.fits')[1].data

# limit magnitude range used in PSF
mag = data['MAG_APER'][:,5]
mask1 = mag > 15
mask2 = mag < 20
mask = mask1*mask2

data = data[mask]

print(len(data))

#squad1data, squad2data, squad3data, squad4data = quadrants(data, data)
quaddata = quadrants(data, data)

for n, quad in enumerate(quaddata):
    # extract stamps column and set non source regions to nans
    stamps = quad['VIGNET']
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
    #plt.figure()
    #plt.imshow(np.log(stack))
    
#    hdu = fits.PrimaryHDU(stack)
#    hdu.writeto(sem+str(n)+'_K_PSF.fits', overwrite=True)
    
    # This image is the typical PSF of the stack - can be used to create a matching kernel
