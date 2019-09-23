#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 14:00:59 2017

Code for degrading images

@author: ppxee
"""
### Import required libraries ###
#import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
import numpy as np #for handling arrays
from astropy.convolution import convolve
from photutils import create_matching_kernel, TopHatWindow, CosineBellWindow
import matplotlib.pyplot as plt
plt.close('all')
    
sdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

colname = 'FWHM_WORLD_'
#data = sem05B[colname][:,1]
sems = ['05B', '06B']#, '07B', '08B', '09B', '10B', '11B', '12B']

#Extract the flux radii and remove negative values
avgFWHM = np.zeros(8)

for n, sem in enumerate(sems):
    colnames = colname+sem
    mag = sdata['MAG_APER_'+sem][:,4]
    mask1 = mag > 15 #removes saturated
    mask2 = mag < 20 #removes very faint stars
    mask = mask1 * mask2
    tempsdata = sdata[mask]
    avgFWHM[n] = np.median(tempsdata[colnames]) #* 3600

   
### Find maximum FWHM as this is what all the others willl become ###
aimind = np.argmax(avgFWHM)
aimsem = sems[aimind]
psf = {}

for semester in sems:
    psf[semester] = fits.open('PSFs/small_'+semester+'_K_PSF.fits')[0].data
    
aimpsf = psf[aimsem]

for semester in sems:
    if semester == aimsem:
#        plt.figure()
#        plt.imshow(np.log(psf[semester]))
        continue
    kernel = create_matching_kernel(psf[semester], aimpsf, window=TopHatWindow(0.5))
    
    plt.figure()
    plt.subplot(121)
#    plt.imshow(kernel)
    plt.imshow(np.log(kernel))
    plt.subplot(122)
#    plt.imshow(psf[semester])
    plt.imshow(np.log(psf[semester]))
#    ### Open image ###
#    im05Bfull = fits.open('UDS_'+semester+'_K.fits', memmap=True)
#    im05B = im05Bfull[0].data
#    hdr = im05Bfull[0].header
#    
#    ### Convolve Image ###
#    newim05B = convolve(im05B, kernel)
#    
#    ### Save the file ###
#    hdu = fits.PrimaryHDU(newim05B, header=hdr)
#    hdulist = fits.HDUList([hdu])
#    hdulist.writeto('newPSF_UDS_'+semester+'_K.fits', overwrite=True)
#    
#    ### CLose the file ###
#    im05Bfull.close()
#    del im05Bfull[0].data
    














