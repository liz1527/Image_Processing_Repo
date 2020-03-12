#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 15:39:33 2018

Convolve each quadrant of an image to the worst quadrants seeing

@author: ppxee
"""

### Import required libraries ###
#import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
import numpy as np #for handling arrays
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
#import vari_funcs_no06 #my module to help run code neatly
#from scipy import ndimage
#import math
#from astropy.stats import median_absolute_deviation

def fluxrad2sigma(fluxrad):
    return fluxrad/np.sqrt(8*np.log(2))

def quadrants(initdata,sem):
    
    ira = initdata['X_IMAGE_'+sem]
    idec = initdata['Y_IMAGE_'+sem]

    ### define bounds of quadrants ###
    midra = 12450
    middec = 13310
    
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
hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM
semesters = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
sdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data

for n, sem in enumerate(semesters):
    # limit data
    mag = sdata['MAG_APER_'+sem][:,4]
    mask1 = mag > 15 #removes saturated
    mask2 = mag < 19 #removes very faint stars
    mask = mask1 * mask2
    if sem == '05B':
        ids = sdata['NUMBER_05B'][mask]
    else:
        ids = np.intersect1d(ids, sdata['NUMBER_05B'][mask])
    
mask = np.isin(sdata['NUMBER_05B'], ids)
tempsdata = sdata[mask] 

for sem in semesters:
    colname = 'FWHM_WORLD_'+sem
    
    quad1data, quad2data, quad3data, quad4data = quadrants(tempsdata, sem)
    
    avgFWHM = np.array([np.median(quad1data[colname])*3600,
                        np.median(quad2data[colname])*3600,
                        np.median(quad3data[colname])*3600,
                        np.median(quad4data[colname])*3600])
    print(avgFWHM)

    ### Find maximum FWHM as this is what all the others willl become ###
    aimind = np.argmax(avgFWHM)
    print(aimind)
    #
    ### Convert FWHM into a sigma ###
    sigmaold = np.array([fluxrad2sigma(fwhm) for fwhm in avgFWHM])
    sigmabroad = sigmaold[aimind]
    sigmakernel = np.array([np.sqrt(sigmabroad**2 - sigma**2) for sigma in sigmaold])
    
    
    ### Open image ###
#    filename = 'Images/UDS_' + sem + '_K_bin2x2.fits'
#    im05Bfull = fits.open(filename, memmap=True)
#    im05B = im05Bfull[0].data
#    hdr5 = im05Bfull[0].header
#    
#    quadsv = np.vsplit(im05B, 2)
#    quadsh1 = np.hsplit(quadsv[0], 2)
#    quadsh2 = np.hsplit(quadsv[1], 2)
#    quad1 = quadsh1[0] #TL
#    quad2 = quadsh1[1] #TR
#    quad3 = quadsh2[0] #BL
#    quad4 = quadsh2[1] #BR
#    
#    ### Convolve Image ###
#    if sigmakernel[0] != 0:
#        print('Convolving quad 1')
#        kernel = Gaussian2DKernel(sigmakernel[0])
#        newquad1 = convolve(quad1, kernel, normalize_kernel=True)
#    else:
#        newquad1 = quad1
#    if sigmakernel[1] != 0:
#        print('Convolving quad 2')
#        kernel = Gaussian2DKernel(sigmakernel[1])
#        newquad2 = convolve(quad2, kernel, normalize_kernel=True)
#    else:
#        newquad2 = quad2
#    
#    if sigmakernel[2] != 0:    
#        print('Convolving quad 3')
#        kernel = Gaussian2DKernel(sigmakernel[2])
#        newquad3 = convolve(quad3, kernel, normalize_kernel=True)
#    else:
#        newquad3 = quad1
#    
#    if sigmakernel[3] != 0:
#        print('Convolving quad 4')
#        kernel = Gaussian2DKernel(sigmakernel[3])
#        newquad4 = convolve(quad4, kernel, normalize_kernel=True)
#    else:
#        newquad4 = quad4
#    
#    ## Recombine image
#    newim05Bbot = np.hstack([newquad3, newquad4])
#    newim05Btop = np.hstack([newquad1, newquad2])
#    newim05B = np.vstack([newim05Btop, newim05Bbot])
#    diff = newim05B-im05B
#    print(np.size(diff[diff!=0]))
#    
#    ### Save the file ###
#    hdu = fits.PrimaryHDU(newim05B, header=hdr5)
#    hdulist = fits.HDUList([hdu])
#    newfilename = 'newq_' + filename
#    hdulist.writeto(newfilename, overwrite=True)
#    im05Bfull.close()
#    del im05Bfull[0].data