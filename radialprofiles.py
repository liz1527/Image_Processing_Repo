#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 10:59:35 2018

@author: ppxee
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.stats import norm
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
    guas = normalise(guas)
    return guas
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_matched.fits')[1].data
oldsdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

r = np.arange(0,42,1) * const * 3600 #define radius values
centre = [29,29]
pixelr = (1.5/3600) / const
aperture = CircularAperture(centre, pixelr)


sems = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
n=1
grid = plt.GridSpec(8, 4, wspace=0.4, hspace=0.05)
phot = {}
flux = np.zeros(len(sems))
oldphot = {}
oldflux = np.zeros(len(sems))
radialprofile = {}
oldradialprofile = {}
for m, sem in enumerate(sems):
    if sem == '10B':
        psf = fits.getdata('PSFs/small_'+sem+'_K_PSF.fits')
    else:
        psf = fits.getdata('PSFs/matched_'+sem+'_K_PSF.fits')
    
#    if sem == '10B':
#        oldpsf = fits.getdata('PSFs/limited_'+sem+'_K_PSF.fits')
#    else:
#        oldpsf = fits.getdata('PSFs/new_limited_'+sem+'_K_PSF.fits')
        
    oldpsf = fits.getdata('PSFs/small_'+sem+'_K_PSF.fits')


    # find radial profiles
    rp = radial_profile(psf, [29,29])
    radialprofile[sem] = normalise(rp)
    sqrtrp = np.sqrt(radialprofile[sem])
    
    oldrp = radial_profile(oldpsf, [29,29])
    oldradialprofile[sem] = normalise(oldrp)
    oldsqrtrp = np.sqrt(oldradialprofile[sem])
    
    ### Determine flux within 3 arcsec apertures ###
    phot[sem] = aperture_photometry(psf, aperture)
    flux[m] = phot[sem]['aperture_sum'][0]
    oldphot[sem] = aperture_photometry(oldpsf, aperture)
    oldflux[m] = oldphot[sem]['aperture_sum'][0]
    
    # find gaussian profile from FWHM
    guas = getguass(sdata, sem)
    guas = normalise(guas)
    sqrtguas = np.sqrt(guas)
    oldguas = getguass(oldsdata, sem)
    oldguas = normalise(oldguas)
    sqrtoldguas = np.sqrt(oldguas)
    
    # plot psf (logged so you can see it)
    plt.figure(1,figsize=(8.27, 11.69), dpi=100)
    plt.subplot(4, 6, n)
#    plt.subplot(grid[n,0])
    plt.imshow(np.log(oldpsf))
    plt.plot(29,29,'+')
    plt.title('PSF '+sem)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    
    plt.subplot(4, 6, n+1)
    plt.imshow(np.log(psf))
    plt.title('new PSF')
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

#     plot radial profile on same plot with its model
    plt.subplot(4, 6, n+2)
#    plt.subplot(4, 2, n)
#    plt.subplot(grid[n,1:])
#    plt.plot(r, radialprofile, label=sem)
#    plt.plot(r, oldradialprofile, '--', label=sem+' before')
#    plt.plot(r, sqrtoldguas, '--', label=sem+' Model')
    plt.plot(r, sqrtrp, label=sem)    
    plt.plot(r, oldsqrtrp, label=sem+' before')
#    plt.yscale('log')
    plt.xlabel('Radius (arcsecs)')
#    plt.ylabel('Flux')
    plt.ylabel('sqrt(Flux)')
#    if sem == '12B':
#        continue
#    plt.tick_params(
#    axis='x',          # changes apply to the x-axis
#    which='both',      # both major and minor ticks are affected
#    bottom=False,      # ticks along the bottom edge are off
#    top=False,         # ticks along the top edge are off
#    labelbottom=False)
    plt.legend()
#    n+=1
    n += 3
    
#     plot radial profile of all semesters
    plt.figure(2)
#    plt.subplot(212)
#    plt.plot(r, radialprofile, label=sem+' after')
    if sem == '10B':
        plt.plot(r, sqrtrp, label=sem+' after')
    else:
        plt.plot(r, sqrtrp, '--', label=sem+' after')
    plt.xlabel('Radius (arcsecs)')
    plt.ylim(ymin=-0.02,ymax=0.57)
    plt.xlim(xmin=0,xmax=1.5)
    plt.ylabel('sqrt(Flux)')
#    plt.yscale('log')
#    plt.ylabel('log(Flux)')
#    plt.ylabel('Flux')
    plt.legend()
#    plt.subplot(211)
##    plt.plot(r, oldradialprofile, '--', label=sem+' before')
#    plt.plot(r, oldsqrtrp, label=sem+' before')
##    plt.yscale('log')
#    plt.xlabel('Radius (arcsecs)')
#    plt.ylabel('Flux')
#    plt.ylim(ymin=-0.02,ymax=0.57)
#    plt.xlim(xmax=1.5)
##    plt.ylabel('log(Flux)')
#    plt.legend()
#    
    ### Plot models ###
#    plt.figure(3)
#    plt.plot(r, guas, '--', label=sem+' Model')
#    plt.plot(r, sqrtguas, label=sem+' Model After')
#    plt.plot(r, sqrtoldguas, '--', label=sem+' Model Before')
#    plt.xlabel('Radius (arcsecs)')
#    plt.ylabel('sqrt(Flux)')
#    plt.legend()
#
#plt.tight_layout()
    
#plt.figure(1)
#plt.subplots_adjust(left=0.02, right=0.95, bottom=0.07, top=0.97)
#plt.xlabel('Radius (arcsecs)')
##    plt.ylabel('Flux')
#plt.ylabel('sqrt(Flux)')

plt.figure()
vari_funcs.avg_lightcurve(oldflux, shape='s', size=10, label='PSF Flux Before')
vari_funcs.avg_lightcurve(flux, label='PSF Flux After')
plt.ylabel('Flux of PSF withing 3 arcsec aperture')
plt.title('')
plt.tight_layout()

#compare new radial profiles
aim = radialprofile['10B']
plt.figure()
for key in radialprofile:
    if key == '10B':
        continue
    diff = aim - radialprofile[key]
    plt.plot(r, diff, '--', label=key)

plt.hlines(0,0,1.5, label='10B')
plt.xlim(xmin=0,xmax=1.5)
plt.ylim(ymin=-0.04, ymax=0.02)
plt.xlabel('Radius (arcsecs)')
plt.ylabel('Flux Difference')
plt.legend()

### Compare old radial profiles ###
aim = oldradialprofile['10B']
plt.figure()
for key in oldradialprofile:
    if key == '10B':
        continue
    diff = aim - oldradialprofile[key]
    plt.plot(r, diff, '--', label=key)

plt.hlines(0,0,1.5, label='10B')
plt.xlim(xmin=0,xmax=1.5)
plt.ylim(ymin=-0.04, ymax=0.02)
plt.xlabel('Radius (arcsecs)')
plt.ylabel('Flux Difference')
plt.legend()