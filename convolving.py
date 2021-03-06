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
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
#from scipy import ndimage
#import math
#from astropy.stats import median_absolute_deviation

def FWHM2sigma(FWHM):
    ''' Function to convert the FWHM of a distribution into a sigma for that
    distribution. It assumes the distribution is gaussian.
    Input:
        FWHM = Full width half maximum of a distriubtution (in my case usually
                of an object from SExtractor)
    Output:
        sigma = standard deviation value of a guassian distribution with the 
                given FWHM. This roughly equates to the psf of the object. '''
    hdr08B = fits.getheader('UDS_08B_K.fits')
    const = hdr08B['CD1_1'] # constant that defines unit conversion for FWHM
    FWHM /= const
    return FWHM/np.sqrt(8*np.log(2))

sem05B = fits.open('SE_outputs_yearstacks/05B_output.fits')[1].data
sem07B = fits.open('SE_outputs_yearstacks/07B_output.fits')[1].data
sem08B = fits.open('SE_outputs_yearstacks/08B_output.fits')[1].data
sem09B = fits.open('SE_outputs_yearstacks/09B_output.fits')[1].data
sem10B = fits.open('SE_outputs_yearstacks/10B_output.fits')[1].data
sem11B = fits.open('SE_outputs_yearstacks/11B_output.fits')[1].data
sem12B = fits.open('SE_outputs_yearstacks/12B_output.fits')[1].data

#sem05B = fits.open('05B_output.fits')[1].data
#sem07B = fits.open('07B_output.fits')[1].data
#sem08B = fits.open('08B_output.fits')[1].data
#sem09B = fits.open('09B_output.fits')[1].data
#sem10B = fits.open('10B_output.fits')[1].data
#sem11B = fits.open('11B_output.fits')[1].data
#sem12B = fits.open('12B_output.fits')[1].data

colname = 'FWHM_WORLD'

#Put data in array
avgFWHM = np.array([np.mean(sem05B[colname]), np.mean(sem07B[colname]), 
                 np.mean(sem08B[colname]), np.mean(sem09B[colname]), 
                 np.mean(sem10B[colname]), np.mean(sem11B[colname]),
                 np.mean(sem12B[colname])])

### Find maximum FWHM as this is what all the others willl become ###
aimind = np.argmax(avgFWHM)

### Convert FWHM into a sigma ###
sigmaold = np.array([FWHM2sigma(fwhm) for fwhm in avgFWHM])
sigmabroad = sigmaold[aimind]

### Find required sigma ###
# sigker^2 = sigbroad^2 - signar^2
sigmakernel = np.array([np.sqrt(sigmabroad**2 - sigma**2) for sigma in sigmaold])

### Open images ###
#im05Bfull = fits.open('UDS_05B_K_bin2x2.fits', memmap=True)
#im05B = im05Bfull[0].data
#hdr = im05Bfull[0].header

### Convolve Images ###
kernel05B = Gaussian2DKernel(sigmakernel[0])
#newim05B = convolve(im05B, kernel05B)
#im05Bfull[0].data = newim05B

### Save the file ###
#hdu = fits.PrimaryHDU(newim05B, header=hdr)
#hdulist = fits.HDUList([hdu])
#hdulist.writeto('new_UDS_05B_K_bin2x2.fits', overwrite=True)
#
#### CLose the file ###
#im05Bfull.close()
#del im05Bfull[0].data















