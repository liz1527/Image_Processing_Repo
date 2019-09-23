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

hdr08B = fits.getheader('UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM
semesters = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
all = fits.open('SE_outputs_yearstacks/05B_output.fits')[1].data
stars = fits.open('stars_mag_flux_table.fits')[1].data

for sem in semesters:
    #stars = fits.open('SE_outputs_yearstacks/newq_12B_stars.fits')[1].data
    #colname = 'FWHM_WORLD'
    colname = 'FLUX_RADIUS_'+sem
    
    #sra = stars['ALPHA_J2000']#_05B']
    #sdec = stars['DELTA_J2000']#_05B']
    
    sra = stars['ALPHA_J2000_05B']
    sdec = stars['DELTA_J2000_05B']
    
    ra = all['ALPHA_J2000']
    dec = all['DELTA_J2000']
    
    ### define bounds of quadrants ###
    maxra = np.max(ra)
    minra = np.min(ra)
    maxdec = np.max(dec)
    mindec = np.min(dec)
    
    midra = (maxra+minra)/2
    middec = (maxdec+mindec)/2
    
    ### separate quadrants ###
    #BR
    mask1 = sra < midra
    mask2 = sdec < middec
    quad1data = stars[mask1*mask2]
    #BL
    mask1 = sra >= midra
    mask2 = sdec < middec
    quad2data = stars[mask1*mask2]
    #TR
    mask1 = sra < midra
    mask2 = sdec >= middec
    quad3data = stars[mask1*mask2]
    #TL
    mask1 = sra >= midra
    mask2 = sdec >= middec
    quad4data = stars[mask1*mask2]
        
    
#    avgFWHM = np.array([np.median(quad1data[colname])/const,
#                        np.median(quad2data[colname])/const,
#                        np.median(quad3data[colname])/const,
#                        np.median(quad4data[colname])/const])
    avgFWHM = np.array([2*np.median(quad1data[colname][:,1]),
                        2*np.median(quad2data[colname][:,1]),
                        2*np.median(quad3data[colname][:,1]),
                        2*np.median(quad4data[colname][:,1])])
    print(avgFWHM)

    ### Find maximum FWHM as this is what all the others willl become ###
    aimind = np.argmax(avgFWHM)
    #
    ### Convert FWHM into a sigma ###
    sigmaold = np.array([fluxrad2sigma(fwhm) for fwhm in avgFWHM])
    sigmabroad = sigmaold[aimind]
    sigmakernel = np.array([np.sqrt(sigmabroad**2 - sigma**2) for sigma in sigmaold])
    
    
    ### Open image ###
#    filename = 'UDS_' + sem + '_K.fits'
#    im05Bfull = fits.open(filename, memmap=True)
#    im05B = im05Bfull[0].data
#    hdr5 = im05Bfull[0].header
#    
#    quadsv = np.vsplit(im05B, 2)
#    quadsh1 = np.hsplit(quadsv[0], 2)
#    quadsh2 = np.hsplit(quadsv[1], 2)
#    quad1 = quadsh2[1] #BR
#    quad2 = quadsh2[0] #BL
#    quad3 = quadsh1[1] #TR
#    quad4 = quadsh1[0] #TL
#    
#    ### Convolve Image ###
#    if sigmakernel[0] != 0:
#        print('Convolving quad 1')
#        kernel = Gaussian2DKernel(sigmakernel[0])
#        newquad3 = convolve(quad3, kernel, normalize_kernel=True)
#    else:
#        newquad3 = quad3
#    if sigmakernel[1] != 0:
#        print('Convolving quad 2')
#        kernel = Gaussian2DKernel(sigmakernel[1])
#        newquad4 = convolve(quad4, kernel, normalize_kernel=True)
#    else:
#        newquad4 = quad4
#    
#    if sigmakernel[2] != 0:    
#        print('Convolving quad 3')
#        kernel = Gaussian2DKernel(sigmakernel[2])
#        newquad1 = convolve(quad1, kernel, normalize_kernel=True)
#    else:
#        newquad1 = quad1
#    
#    if sigmakernel[3] != 0:
#        print('Convolving quad 4')
#        kernel = Gaussian2DKernel(sigmakernel[3])
#        newquad2 = convolve(quad2, kernel, normalize_kernel=True)
#    else:
#        newquad2 = quad2
#    
#    
#    ## Recombine image
#    newim05Bbot = np.hstack([newquad2, newquad1])
#    newim05Btop = np.hstack([newquad4, newquad3])
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