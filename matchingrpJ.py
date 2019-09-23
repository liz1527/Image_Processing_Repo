#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 12:24:26 2018

@author: ppxee
"""
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.stats import norm
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from photutils import CircularAperture, aperture_photometry
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
    newpsf /= np.nansum(newpsf)
#    print(np.sum(newpsf))
    return newpsf
    
#sdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
oldsdata = fits.open('mag_flux_tables/stars_mag_flux_table_J.fits')[1].data
hdr08B = fits.getheader('Images/UDS-DR11-K.mef.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

colname = 'FWHM_WORLD_'
#data = sem05B[colname][:,1]
semesters = ['05B', '06B', '07B', '08B', '10B', '11B', '12B']#['05B','10B']
centre = [29,29]

avgFWHM = np.zeros(len(semesters))
oldavgFWHM = np.zeros(len(semesters))
avgflux = np.zeros(len(semesters))
#oldavgflux = np.zeros(len(semesters))
psf = {}
oldpsf = {}
for n, sem in enumerate(semesters):
    # for new
    colnames = colname+sem
    # for old
    oldmag = oldsdata['MAG_APER_'+sem][:,4]
    mask1 = oldmag > 15 #removes saturated
    mask2 = oldmag < 19 #removes very faint stars
    oldmask = mask1 * mask2
    tempsdata = oldsdata[oldmask]
    oldavgFWHM[n] = np.median(tempsdata[colnames]) #* 3600
#    flux = tempsdata['FLUX_APER_'+sem][:,4]
#    oldavgflux[n] = np.median(flux)
    oldpsf[sem] = fits.open('PSFs/'+sem+'_J_PSF.fits')[0].data
#    if sem == '10B':
#        oldpsf[sem] = fits.open('PSFs/limited_'+sem+'_K_PSF.fits')[0].data
#    else:
#        oldpsf[sem] = fits.open('PSFs/extra_'+sem+'_K_PSF.fits')[0].data

## get flux curve

flux = vari_funcs.jflux4_stacks(tempsdata)
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

### testing the extra factor method ###
tests =np.linspace(-0.5,0.11,1000)
r = np.arange(0,42,1) * const * 3600 # define radius values

#flux10B = oldflux[3]
aimrp = radial_profile(oldpsf[aimsem], centre)
sqrtaimrp = np.sqrt(aimrp)


newpsf = {}
newphot = {}
newflux = np.zeros(len(semesters))
extras = np.zeros(len(semesters))
aperflux = np.empty(len(semesters))
oldaperflux = np.empty(len(semesters))
pixelr = (1.5/3600) / const
aperture = CircularAperture(centre, pixelr)
t = np.linspace(1, 8, num=8)

plt.figure()
plt.plot(r, sqrtaimrp,label=aimsem)   
for n, sem in enumerate(semesters):
    if sem == aimsem:
        phot = aperture_photometry(aimpsf, aperture)
        aperflux[n] = phot['aperture_sum'][0]
        oldaperflux[n] = phot['aperture_sum'][0]
        plt.figure()
        plt.imshow(np.log(aimpsf))
        vari_funcs.no_ticks()
        continue
    
    sigma = sigmaold[n]
    singlephot = {}
    singleflux = np.zeros(len(tests))
    sumdiffold = 10
    for m, extra in enumerate(tests):
#        print(extra)
        sigmakernel = np.sqrt(sigmabroad**2 - sigma**2) + extra
        if sigmakernel <= 0:
            continue
        new = convolve_one_psf(oldpsf[sem], sigmakernel)    
        
        ### Get radial profile ###
        radialprofile = radial_profile(new, centre)
        sqrtrp = np.sqrt(radialprofile)
        
        diff = aimrp[:12] - radialprofile[:12]
        sumdiff = np.nansum(diff)
#        plt.figure(1)
##            plt.subplot(4,2,n+1)
##            plt.plot(r, sqrtaimrp,label='10B')    
#        plt.plot(r,sqrtrp, '--', label=sem+' '+str(extra))
#        plt.ylabel('sqrt(Flux)')
#        plt.xlabel('Radius (arcsec)')
#        plt.legend()
#        print(sumdiff)
        if sumdiff > 0:
            if sumdiffold < sumdiff:
                print('extra ='+str(tests[m-1]))
                extras[n] = tests[m-1]
#                plt.plot(t[n], singleflux[m-1],'o', label=extras[n])
            else:
                print('extra ='+str(extra))
                extras[n] = extra    
            
            plt.figure(1)
#            plt.subplot(4,2,n+1)
#            plt.plot(r, sqrtaimrp,label='10B')    
            plt.plot(r,sqrtrp, '--', label=sem+' '+str(extra))
            plt.ylabel('sqrt(Flux)')
            plt.xlabel('Radius (arcsec)')
#            plt.xlim(xmax=1.5)
            plt.legend()
            plt.figure(2)
            plt.plot(r[:12], diff, label=sem)
            plt.hlines(0,0,1.5)
            plt.xlim(xmin=0, xmax=1.5)
            plt.legend()
            plt.ylabel('Difference from aim rp')
            plt.xlabel('Radius (arcsec)')
            
            phot = aperture_photometry(new, aperture)
            aperflux[n] = phot['aperture_sum'][0]
            oldphot = aperture_photometry(oldpsf[sem], aperture)
            oldaperflux[n] = oldphot['aperture_sum'][0]
            
#            plt.figure()
#            plt.subplot(121)
#            plt.imshow(np.log(oldpsf[sem]))
#            vari_funcs.no_ticks()
#            
#            plt.subplot(122)
#            plt.imshow(np.log(new))
#            vari_funcs.no_ticks()
            break
        else:
            sumdiffold = sumdiff       

x = [1,2,3,4,6,7,8]
years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
plt.figure(figsize=[9,6])
plt.plot(x, aperflux,'o-', label='new')
plt.plot(x, oldaperflux, 'o-', label='old')
plt.xticks(t, years)
plt.xlabel('Semester')
plt.legend()
plt.tight_layout()
#np.save('extrascleanedno06', extras)