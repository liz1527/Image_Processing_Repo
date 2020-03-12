#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 14:40:09 2018

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
import numpy as np #for handling arrays
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
#import vari_funcs
#from scipy import ndimage
#import math
#from astropy.stats import median_absolute_deviation
plt.close('all')
import time
start = time.time()
print(start)

def FWHM2sigma(FWHM, const):
    ''' Function to convert the FWHM of a distribution into a sigma for that
    distribution. It assumes the distribution is gaussian.
    Input:
        FWHM = Full width half maximum of a distriubtution (in my case usually
                of an object from SExtractor)
    Output:
        sigma = standard deviation value of a guassian distribution with the 
                given FWHM. This roughly equates to the psf of the object. '''
    FWHM /= const
    return FWHM/np.sqrt(8*np.log(2))
def fluxrad2sigma(fluxrad):
    return fluxrad/np.sqrt(8*np.log(2))
    
stars = fits.open('mag_flux_tables/stars_mag_flux_table.fits')
sdata = stars[1].data
hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

colname = 'FWHM_WORLD_'
#data = sem05B[colname][:,1]
semesters = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']#['05B','10B']

avgFWHM = np.zeros(8)
psf = {}

for n, sem in enumerate(semesters):
    colnames = colname+sem
    mag = sdata['MAG_APER_'+sem][:,4]
    mask1 = mag > 15 #removes saturated
    mask2 = mag < 20 #removes very faint stars
    mask = mask1 * mask2
    tempsdata = sdata[mask]
    avgFWHM[n] = np.median(tempsdata[colnames]) #* 3600
    psf[sem] = fits.open('PSFs/limited_'+sem+'_K_PSF.fits')[0].data
    
   
### Find maximum FWHM as this is what all the others willl become ###
aimind = np.argmax(avgFWHM)
aimsem = semesters[aimind]
aimpsf = psf[aimsem]

#
### Convert FWHM into a sigma ###
sigmaold = np.array([FWHM2sigma(fwhm, const) for fwhm in avgFWHM])
sigmabroad = sigmaold[aimind]

### Find required sigma ###
# sigker^2 = sigbroad^2 - signar^2
sigmakernel = np.array([np.sqrt(sigmabroad**2 - sigma**2) for sigma in sigmaold])

#vari_funcs.avg_lightcurve(avgFWHM)
#plt.ylabel('FWHM (arcsec)')
#plt.title('Median FWHM of stars between 15th and 20th magnitude')

def convolve_image(filename, sigmakernel):
    ## Open image ###
    im05Bfull = fits.open(filename, memmap=True)
    im05B = im05Bfull[0].data
    hdr5 = im05Bfull[0].header
    ## Convolve Image ###
    print('Convolving', filename)
    kernel = Gaussian2DKernel(sigmakernel)
    plt.figure()
    plt.imshow(kernel)
    newim05B = convolve(im05B, kernel, normalize_kernel=True)
    plt.figure()
    plt.subplot(211)
    plt.imshow(np.log(newim05B))
    plt.subplot(212)
    plt.imshow(np.log(aimpsf))
#    ### Save the file ###
    hdu = fits.PrimaryHDU(newim05B, header=hdr5)
    hdulist = fits.HDUList([hdu])
    newfilename = 'conv_limited_' + sem +'_K_PSF.fits'
    hdulist.writeto(newfilename, overwrite=True)
    im05Bfull.close()
    del im05Bfull[0].data
    end = time.time()
    print(end-start)

for n, sem in enumerate(semesters):
    if sem == semesters[aimind]:
        continue
    filename = 'PSFs/limited_' + sem + '_K_PSF.fits'
    convolve_image(filename, sigmakernel[n])

end = time.time()
print(end-start)