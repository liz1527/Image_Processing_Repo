#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 15:16:54 2018

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
    
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_extra.fits')[1].data
oldsdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

colname = 'FWHM_WORLD_'
#data = sem05B[colname][:,1]
semesters = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']#['05B','10B']
centre = [63,63]
centresmall = [29,29]

avgFWHM = np.zeros(len(semesters))
oldavgFWHM = np.zeros(len(semesters))
avgflux = np.zeros(len(semesters))
#oldavgflux = np.zeros(len(semesters))
psf = {}
oldpsf = {}
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
    psf[sem] = fits.open('PSFs/limited_'+sem+'_K_PSF.fits')[0].data
#    if sem == '10B':
#        psf[sem] = fits.open('PSFs/limited_'+sem+'_K_PSF.fits')[0].data
#    else:
#        psf[sem] = fits.open('PSFs/extra_'+sem+'_K_PSF.fits')[0].data
#    # for old
    oldmag = sdata['MAG_APER_'+sem][:,4]
    mask1 = mag > 15 #removes saturated
    mask2 = mag < 20 #removes very faint stars
    oldmask = mask1 * mask2
    tempsdata = oldsdata[oldmask]
    oldavgFWHM[n] = np.median(tempsdata[colnames]) #* 3600
#    flux = tempsdata['FLUX_APER_'+sem][:,4]
#    oldavgflux[n] = np.median(flux)
    oldpsf[sem] = fits.open('PSFs/small_'+sem+'_K_PSF.fits')[0].data

## get flux curve

flux = vari_funcs.flux5_stacks(tempsdata)
flux = vari_funcs.normalise_flux(flux)
oldavgflux = np.median(flux, axis=0)
   
### Find maximum FWHM as this is what all the others willl become ###
aimind = np.argmax(oldavgFWHM)
aimsem = semesters[aimind]
aimpsf = oldpsf[aimsem]

### Convert FWHM into a sigma ###
sigmaold = np.array([FWHM2sigma(fwhm, const) for fwhm in oldavgFWHM])
sigmabroad = sigmaold[aimind]

phot = {}
flux = np.zeros(len(semesters))
oldphot = {}
oldflux = np.zeros(len(semesters))

pixelr = (1.5/3600) / const
aperture = CircularAperture(centre, pixelr)
smallaperture = CircularAperture(centresmall, pixelr)
plt.figure(4)
plt.imshow(np.log(psf['12B']))
aperture.plot()

for n, sem in enumerate(semesters):
    ### Determine flux within 3 arcsec apertures ###
    phot[sem] = aperture_photometry(psf[sem], aperture)
    flux[n] = phot[sem]['aperture_sum'][0]
    oldphot[sem] = aperture_photometry(oldpsf[sem], smallaperture)
    oldflux[n] = oldphot[sem]['aperture_sum'][0]

#plt.figure()
#plt.plot(oldflux, 'bs')
#plt.plot(flux,'ro')
#plt.ylabel('Flux within 3 arcsec aperture')
#plt.figure()
#ax = plt.subplot(111)
#plt.plot(oldavgFWHM, 'bs')
#plt.plot(avgFWHM, 'ro')
#plt.ylabel('FWHM')
#ax.invert_yaxis()

#set up time variable for plot
t = np.linspace(1, 8, num=8)
years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']

plt.figure(1)
#plt.plot(t, oldflux, 'bs')
plt.plot(t, flux, 'bs')
plt.ylabel('PSF Flux within 3 arcsec')
#plt.xlim(xmax=7)
ax = plt.axes()
plt.xticks(t, years)
plt.xlabel('Semester')
plt.tight_layout()
#        plt.legend()

plt.figure(3)
plt.plot(t,oldavgflux, 'bs')
plt.ylabel('Median Flux of stars')
ax = plt.axes()
plt.xticks(t, years)
plt.xlabel('Semester')
plt.tight_layout()

### testing the extra factor method ###
tests = np.linspace(1,2.8,500)
r = np.arange(0,90,1) * const * 3600 # define radius values

flux10B = oldflux[3]

newpsf = {}
newphot = {}
newflux = np.zeros(len(semesters))
extras = np.zeros(len(semesters))
#plt.figure()
#for n, sem in enumerate(semesters):
#    if sem == '08B':
#        radialprofile = radial_profile(oldpsf[sem], centre)
#    #    radialprofile = normalise(radialprofile)
#        sqrtrp = np.sqrt(radialprofile)
#        
#        plt.figure(2)
#        plt.plot(r,sqrtrp, label=sem)
#        plt.ylabel('sqrt(Flux)')
#        plt.xlabel('Radius (arcsec)')
#        plt.legend()
#        continue
#    
#    sigma = sigmaold[n]
#    singlephot = {}
#    singleflux = np.zeros(len(tests))
#    
#    for m, extra in enumerate(tests):
##        print(extra)
#        sigmakernel = np.sqrt(sigmabroad**2 - sigma**2) + extra
#
#        new = convolve_one_psf(oldpsf[sem], sigmakernel)    
#        ### Determine flux within 3 arcsec apertures ###
#        singlephot[extra] = aperture_photometry(new, aperture)
#        singleflux[m] = singlephot[extra]['aperture_sum'][0]
#        
##        plt.figure(1)
##        plt.plot(n, singleflux[m],'o', label=extra)
##        plt.legend()
#
#        plt.figure(1)
#        diff = flux10B - singleflux[m]
#        if diff > 0:
#            diff2 = flux10B - singleflux[m-1]
#            if diff2 < diff:
#                print('extra ='+str(tests[m-1]))
#                extras[n] = tests[m-1]
#                plt.plot(t[n], singleflux[m-1],'o', label=extras[n])
#            else:
#                print('extra ='+str(extra))
#                extras[n] = extra
#                plt.plot(t[n], singleflux[m],'o', label=extras[n])
##        else:
##            diffold = diff
#            sigmakernel = np.sqrt(sigmabroad**2 - sigma**2) + extras[n]
#            newpsf[sem] = convolve_one_psf(oldpsf[sem], sigmakernel)
#    #             ##find radial profiles
#            radialprofile = radial_profile(newpsf[sem], centre)
#        #    radialprofile = normalise(radialprofile)
#            sqrtrp = np.sqrt(radialprofile)
#            
#            plt.figure(2)
#            plt.plot(r,sqrtrp,'--', label=sem)
#            plt.ylabel('sqrt(Flux)')
#            plt.xlabel('Radius (arcsec)')
#            plt.xlim(xmin=0, xmax=1.5)
#            plt.legend()
#            break
#
#        