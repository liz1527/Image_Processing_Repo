#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 11:51:58 2018

Code to test models of PSFs

@author: ppxee
"""


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
plt.close('all')

def radial_profile(data, center):
    y, x = np.indices((data.shape)) #create coordinate grid
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2) #get radius values for grid
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel()) # counts number of times value
                                                # of radius occurs in the psf
                                                # weighted by the data
    nr = np.bincount(r.ravel()) # counts number of radii values in psf
    radialprofile = tbin / nr # as weighted is r*data then get profile by 
                              # dividing by unweighted counts of r values.
    return radialprofile 

def FWHM2sigma_arcsecs(FWHM):
    ''' Function to convert the FWHM of a distribution into a sigma for that
    distribution. It assumes the distribution is gaussian.
    Input:
        FWHM = Full width half maximum of a distriubtution (in my case usually
                of an object from SExtractor)
    Output:
        sigma = standard deviation value of a guassian distribution with the 
                given FWHM. This roughly equates to the psf of the object. '''
    return FWHM/np.sqrt(8*np.log(2))

def normalise(array):
    return array/np.nansum(array)

def getguass(sdata, sem):
    colname = 'FWHM_WORLD_'
    colnames = colname+sem
    mag = sdata['MAG_APER_'+sem][:,4]
    mask1 = mag > 15 #removes saturated
    mask2 = mag < 20 #removes very faint stars
    mask = mask1 * mask2
    tempsdata = sdata[mask]
    fwhm = np.median(tempsdata[colnames]) * 3600
    print(sem+'='+str(fwhm))
    sig = FWHM2sigma_arcsecs(fwhm)
    guas = norm.pdf(r, 0, sig)
#    guas = normalise(guas)
    return guas, sig

def my_gaussian(xdata, mean, sigma):
    guass = (1/(sigma*np.sqrt(2*np.pi))) * np.exp(-0.5 * ((xdata - mean)/sigma)**2)
    return guass

def origin_gaussian(xdata, sigma):
    guass = (1/(sigma*np.sqrt(2*np.pi))) * np.exp(-0.5 * ((xdata - 0)/sigma)**2)
    return guass

def basic_gaussian(xdata, a, b, c):
    guass = a * np.exp(-((xdata - b)**2)/(2*(c**2)))
    return guass

sdata = fits.open('mag_flux_tables/stars_mag_flux_table_new_limited.fits')[1].data
oldsdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

sems = ['10B']#['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
n=1
for sem in sems:
    if sem == '10B':
        psf = fits.getdata('PSFs/limited_'+sem+'_K_PSF.fits')
    else:
        psf = fits.getdata('PSFs/new_limited_'+sem+'_K_PSF.fits')


    # find radial profiles
    radialprofile = radial_profile(psf, [63,63])
#    radialprofile = normalise(radialprofile)
    sqrtrp = np.sqrt(radialprofile)
    
    r = np.arange(0,90,1) * const * 3600 # define radius values

    # find gaussian profile from FWHM
    guas, sig1 = getguass(sdata, sem)
    sqrtguas = np.sqrt(guas)
#    
#    # plot psf (logged so you can see it)
#    plt.figure(1)
#    plt.subplot(121)
#    plt.imshow(np.log(psf))
#    plt.title('new PSF')
#
#    # plot radial profile on same plot with its model
#    plt.subplot(322)
#    plt.plot(r, radialprofile, label=sem)
#    plt.xlabel('Radius (arcsecs)')
#    plt.ylabel('Flux')
#    plt.legend()
#    plt.subplot(324)
#    plt.plot(r, sqrtrp, label=sem)    
#    plt.xlabel('Radius (arcsecs)')
#    plt.ylabel('sqrt(Flux)')
#    plt.legend()
#    plt.subplot(326)
#    plt.plot(r, radialprofile, label=sem)
#    plt.yscale('log')
#    plt.xlabel('Radius (arcsecs)')
#    plt.ylabel('log(Flux)')
#    plt.legend()
#    plt.tight_layout(pad=1)
#    n+=1
#    n += 3


### set up plots wit radial profile ###
plt.figure(1) 
plt.subplot(211)
plt.title('My gaussian function')
plt.plot(r, sqrtrp, label=sem+' real')
plt.xlabel('Radius (arcsecs)')
plt.ylabel('sqrt(Flux)')
plt.xlim(xmin=0, xmax=1.5)
plt.legend()

plt.subplot(212)
plt.plot(r, radialprofile, label=sem+' after')
plt.xlabel('Radius (arcsecs)')
plt.ylabel('Flux')
plt.xlim(xmin=0, xmax=1.5)
plt.legend()
    
## Make test guassians so you can see how it changes ###
    
#sigs = np.linspace(sig1, sig1+1, 10)
#for sig in sigs:
#    
#    gaus2 = my_gaussian(r, 0, sig)
#    sqrtgaus2 = np.sqrt(gaus2)
#    
#    plt.figure(1) 
#    plt.subplot(211)
#    plt.plot(r, sqrtgaus2, '--', label=sem+' Model')
#    plt.legend()
#    
#    plt.subplot(212)
#    plt.plot(r, gaus2, '--', label=sem+' Model')
#    plt.legend()

# Try curve_fit function
init_vals = [0, sig1]
popt, pcov = curve_fit(my_gaussian, r, radialprofile, p0=init_vals)

fitgaus = my_gaussian(r, popt[0], popt[1])
sqrtfitgaus = np.sqrt(fitgaus)
fitfwhm = popt[1] * np.sqrt(8*np.log(2))
plt.figure(1) 
plt.subplot(211)
plt.plot(r, sqrtfitgaus, '--', label=sem+' 2 param fit')
plt.legend()

plt.subplot(212)
plt.plot(r, fitgaus, '--', label=sem+' 2 param fit')
plt.legend()

## Try curve_fit function with fixed mean
#popt, pcov = curve_fit(origin_gaussian, r, radialprofile, p0=sig1)
#
#fitgaus2 = origin_gaussian(r, popt)
#sqrtfitgaus2 = np.sqrt(fitgaus2)
#fitfwhm = popt * np.sqrt(8*np.log(2))
#
#plt.figure(1) 
#plt.subplot(211)
#plt.plot(r, sqrtfitgaus2, '--', label=sem+' 1 param fit')
#plt.legend()
#
#plt.subplot(212)
#plt.plot(r, fitgaus2, '--', label=sem+' 1 param fit')
#plt.legend()

# Try curve_fit function with 3 params
init_vals = [1, 0, sig1]
popt, pcov = curve_fit(basic_gaussian, r, radialprofile, p0=init_vals)

fitgaus2 = basic_gaussian(r, popt[0], popt[1], popt[2])
sqrtfitgaus2 = np.sqrt(fitgaus2)

plt.figure(1) 
plt.subplot(211)
plt.plot(r, sqrtfitgaus2, '--', label=sem+' 3 param fit')
plt.legend()

plt.subplot(212)
plt.plot(r, fitgaus2, '--', label=sem+' 3 param fit')
plt.legend()