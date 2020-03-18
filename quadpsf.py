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
from matplotlib.patches import Rectangle
plt.close('all')
import sys

sem = '05B'#str(sys.argv[1]) #1st arguement from call is the semester

print(sem)

def norm(array):
    return array/np.nansum(array)

def quadrants(initdata):
    
    ira = initdata['X_IMAGE']
    idec = initdata['Y_IMAGE']

    ### define bounds of quadrants ###
    midra = 12450
    middec = 13310
    
    
    ### create masks for quadrant ###
    mask1 = ira < midra
    mask2 = idec >= middec
    quad1data = initdata[mask1*mask2]
    rect1 = Rectangle([0,middec],midra,middec, color='k')
    
    mask1 = ira >= midra
    mask2 = idec >= middec
    quad2data = initdata[mask1*mask2]
    rect2 = Rectangle([midra,middec],midra,middec, color='m')
    
    mask1 = ira < midra
    mask2 = idec < middec
    quad3data = initdata[mask1*mask2]
    rect3 = Rectangle([0,0],midra,middec, color='b')
    
    mask1 = ira >= midra
    mask2 = idec < middec
    quad4data = initdata[mask1*mask2]
    rect4 = Rectangle([midra,0],midra,middec, color='r')
    
    
#     plot quads to check
    plt.figure(1)
    ax = plt.subplot(111)
    ax.add_patch(rect1)
    ax.add_patch(rect2)
    ax.add_patch(rect3)
    ax.add_patch(rect4)
    plt.xlim(xmin=0,xmax=midra*2)
    plt.ylim(ymin=0,ymax=middec*2)
    
    
    return quad1data, quad2data, quad3data, quad4data

## reads in file with stamps
data = fits.open('star_stamps_tables/small_'+sem+'_star_stamps_table.fits')[1].data
#
## limit magnitude range used in PSF
#mag = data['MAG_APER'][:,5]
#mask1 = mag > 15
#mask2 = mag < 19
#mask = mask1*mask2
#
#data = data[mask]
#
#print(len(data))
#
#squad1data, squad2data, squad3data, squad4data = quadrants(data, data)
quaddata = quadrants(data)
#
### Plot to check ###

plt.scatter(quaddata[0]['X_IMAGE'],quaddata[0]['Y_IMAGE'],c='r',zorder=3)
plt.scatter(quaddata[1]['X_IMAGE'],quaddata[1]['Y_IMAGE'],c='b',zorder=3)
plt.scatter(quaddata[2]['X_IMAGE'],quaddata[2]['Y_IMAGE'],c='m',zorder=3)
plt.scatter(quaddata[3]['X_IMAGE'],quaddata[3]['Y_IMAGE'],c='k',zorder=3)
#
#for n, quad in enumerate(quaddata):
#    # extract stamps column and set non source regions to nans
#    stamps = quad['VIGNET']
#    stamps[stamps<-5e29] = np.nan
#    
#    # Plot some stars- not sure why, just for demonstration?
#    #plt.figure()
#    #plt.imshow(np.log(data['VIGNET'][24]))
#    #plt.figure()
#    #plt.imshow(np.log(stamps[1000,:,:]))
#    
#    # Normalse all of the images - I presume to 1? Can't really see the difference...?
#    stampsnorm = stamps
#    for i in range(stamps.shape[0]):
#        stampsnorm[i,:,:] = norm(stamps[i,:,:])
#    
#    # find the median image
#    stack = np.nanmedian(stampsnorm, axis=0)
#    
#    # print shape to check its right?
#    print(stack.shape)
#    
#    #normalise and plot the median image
#    stack = norm(stack)
##    plt.figure(2)
##    plt.subplot(2,2,n+1)
##    plt.imshow(np.log(stack))
#    
#    hdu = fits.PrimaryHDU(stack)
#    hdu.writeto(sem+'_'+str(n+1)+'_K_PSF.fits', overwrite=True)
##    
    # This image is the typical PSF of the stack - can be used to create a matching kernel
