#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 14:35:15 2018

Code that plots aperture flux and median flux data for a single set of PSFs

@author: ppxee
"""
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.stats import norm
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from photutils import CircularAperture
from photutils import aperture_photometry
import vari_funcs
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

def convolve_psf(sem, psf, newpsf, sigmakernel):
    ## Open image ###
    im05B = psf[sem]
    
    ## Convolve Image ###
    print('Convolving ', sem)
    kernel = Gaussian2DKernel(sigmakernel)
    newpsf[sem] = convolve(im05B, kernel, normalize_kernel=True) 

def convolve_one_psf(psf, sigmakernel):
#    print(np.sum(psf))
    ## Convolve Image ###
#    print('Convolving ', sem)
    kernel = Gaussian2DKernel(sigmakernel)
    newpsf = convolve(psf, kernel, normalize_kernel=True) 
#    print(np.sum(newpsf))
    return newpsf
    
sdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

colname = 'FWHM_WORLD_'
#data = sem05B[colname][:,1]
semesters = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']#['05B','10B']
centre = [29,29]

avgFWHM = np.zeros(len(semesters))
avgflux = np.zeros(len(semesters))
psf = {}
for n, sem in enumerate(semesters):
    # for new
    colnames = colname+sem
    mag = sdata['MAG_APER_'+sem][:,4]
    mask1 = mag > 15 #removes saturated
    mask2 = mag < 20 #removes very faint stars
    mask = mask1 * mask2
    tempsdata = sdata[mask]
    avgFWHM[n] = np.median(tempsdata[colnames]) #* 3600
    avgflux[n] = np.median(tempsdata['FLUX_APER_'+sem][:,4])
    if sem == '10B':
        psf[sem] = fits.getdata('PSFs/small_'+sem+'_K_PSF.fits')
    else:
        psf[sem] = fits.getdata('PSFs/matched_'+sem+'_K_PSF.fits')
#    psf[sem] = fits.open('PSFs/small_'+sem+'_K_PSF.fits')[0].data
## get flux curve

flux = vari_funcs.flux5_stacks(tempsdata)
flux = vari_funcs.normalise_flux(flux)
   
### Find maximum FWHM as this is what all the others willl become ###
aimind = np.argmax(avgFWHM)
aimsem = semesters[aimind]
aimpsf = psf[aimsem]

### Convert FWHM into a sigma ###
sigmaold = np.array([FWHM2sigma(fwhm, const) for fwhm in avgFWHM])
sigmabroad = sigmaold[aimind]
plt.figure(4)
plt.imshow(np.log(psf['10B']))
#plt.imshow(psf['10B'])

phot = {}
flux = {}

arcsecr = np.array([0.25, 0.5, 0.75, 1, 1.25, 1.5,1.75,2,2.5,3])#[1.5]

for m, r in enumerate(arcsecr):
    pixelr = (r/3600) / const
    aperture = CircularAperture(centre, pixelr)
    plt.figure(4)
    aperture.plot()
    phot[r] = {}
    flux[r] = np.zeros(len(semesters))
    for n, sem in enumerate(semesters):
        ### Determine flux within 3 arcsec apertures ###
        phot[r][sem] = aperture_photometry(psf[sem], aperture)
        flux[r][n] = phot[r][sem]['aperture_sum'][0]


    #set up time variable for plot
    t = np.linspace(1, 8, num=8)
    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
    
    
    plt.figure(1)
    plt.subplot(5,2,m+1)
    plt.plot(t, flux[r], 's',label='r = '+str(r))
#    plt.ylabel('PSF Flux within 3 arcsec')
    #plt.xlim(xmax=7)
#    ax = plt.axes()
    plt.xticks(t, years)
    plt.xlabel('Semester')
    plt.legend()#bbox_to_anchor=(1.1, 1.05)
#    plt.tight_layout()
    
plt.figure(2)
plt.plot(t,avgFWHM*3600,'s')
plt.xticks(t, years)
plt.ylabel('FWHM')
plt.xlabel('Semester')
plt.tight_layout()

plt.figure(3)
plt.plot(t,avgflux, 's')
plt.ylabel('Median Flux of stars')
ax = plt.axes()
plt.xticks(t, years)
plt.xlabel('Semester')
plt.tight_layout()