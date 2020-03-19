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
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
#import vari_funcs_no06 #my module to help run code neatly
#from scipy import ndimage
#import math
#from astropy.stats import median_absolute_deviation

def fluxrad2sigma(fluxrad):
    return fluxrad/np.sqrt(8*np.log(2))

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

def quadrants(initdata,sem):
    
    ira = initdata['X_IMAGE_'+sem]
    idec = initdata['Y_IMAGE_'+sem]

    ### define bounds of quadrants ###
    midra = 12450
    middec = 13310
    
    ### create masks for quadrant ###
    mask1 = ira < midra
    mask2 = idec >= middec
    quad1data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec >= middec
    quad2data = initdata[mask1*mask2]
    
    mask1 = ira < midra
    mask2 = idec < middec
    quad3data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec < middec
    quad4data = initdata[mask1*mask2]    
    
    return quad1data, quad2data, quad3data, quad4data

psf_data = fits.open('UDS_catalogues/DR11_stars_for_PSFs.fits')[1].data
sdata = fits.open('mag_flux_tables/K/month/month_stars_mag_flux_table_K_cleaned.fits')[1].data
hdr08B = fits.getheader('Images/UDS-DR11-K.mef.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM
colname = 'FWHM_WORLD_'
extras = np.load('extrascleanedK_quad_month.npy')

months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07',  
          'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09',  
          'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 
          'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 
          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']

for i, mon in enumerate(months):
    print(mon)
    # for new
    colnames = colname+mon
    
    ### Define coordinates ###
    refcoord = SkyCoord(psf_data['ALPHA_J2000_1']*u.degree, 
                        psf_data['DELTA_J2000_1']*u.degree)
    moncoord = SkyCoord(sdata['ALPHA_J2000_'+mon]*u.degree, 
                        sdata['DELTA_J2000_'+mon]*u.degree)
    
    ### Match catalogues and create new table ###
    idx, d2d , _ = match_coordinates_sky(refcoord, moncoord) 
    tempsdata = sdata[idx]      

    quaddata = quadrants(tempsdata, mon)
    
    avgFWHM = np.array([np.median(quaddata[0][colnames]),
                        np.median(quaddata[1][colnames]),
                        np.median(quaddata[2][colnames]),
                        np.median(quaddata[3][colnames])])
#    print(avgFWHM)

    ### Find 0 extra value ###
    aimind = np.array([0,1,2,3])[extras[i,:]==99]
    print(aimind==np.argmax(avgFWHM))
    #
    ### Convert FWHM into a sigma ###
    sigmaold = np.array([FWHM2sigma(fwhm, const) for fwhm in avgFWHM])
    sigmabroad = sigmaold[aimind]
    sigmakernel = np.zeros(len(sigmaold))
    for n, sigma in enumerate(sigmaold):
        if n==aimind:
            continue
        sigmakernel[n] = np.sqrt(sigmabroad**2 - sigma**2) + extras[i,n]
#    print(sigmakernel)
    
#    ## Open image ###
#    filename = 'month/UDS_' + mon[0:3] + '_' + mon[3:5] + '_K_cleaned.fits'
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
#        newquad3 = quad3
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
##    diff = newim05B-im05B
##    print(np.size(diff[diff!=0]))
#    
#    ### Save the file ###
#    hdu = fits.PrimaryHDU(newim05B, header=hdr5)
#    hdulist = fits.HDUList([hdu])
#    newfilename = 'month/UDS_'+mon[0:3]+'_'+mon[3:5]+'_K_quad_cleaned.fits'
#    hdulist.writeto(newfilename, overwrite=True)
#    im05Bfull.close()
#    del im05Bfull[0].data