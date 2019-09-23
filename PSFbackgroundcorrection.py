#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 13:29:55 2018

Code to sort out background in PSFs

@author: ppxee
"""

### Import required libraries ###
import matplotlib
matplotlib.use('Agg') #so I can save them remotely
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from photutils import CircularAperture, aperture_photometry
plt.close('all') #close any open plots

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

def get_avg_flux(tbdata):
    flux = vari_funcs.flux5_stacks(tbdata)
    flux = vari_funcs.normalise_flux(flux)
    return np.nanmedian(flux, axis=0)

def psf_and_profile(sem):
    centre = [29,29]
    psf = fits.getdata('PSFs/small_'+sem+'_K_PSF.fits')
#    if sem == '10B':
#        psf = fits.getdata('PSFs/small_'+sem+'_K_PSF.fits')
#    else:
#        psf = fits.getdata('PSFs/matched_'+sem+'_K_PSF.fits')

    rp = vari_funcs.radial_profile(psf, centre)
    return psf, rp, np.sqrt(rp)

def psf_and_profile_mag(upper, lower, sem):
    centre = [29,29]
    if sem == '10B':
        psf = fits.getdata('PSFs/small_'+sem+'_K_PSF.fits')
    else:
        psf = fits.getdata('PSFs/'+str(lower)+str(upper)+'bkgcorr_'+sem+'_K_PSF.fits')
#    psf = fits.getdata('PSFs/mag_changes/small/'+str(lower)+str(upper)+'_'+sem+'_K_PSF.fits')
    rp = vari_funcs.radial_profile(psf, centre)
    return psf, rp, np.sqrt(rp)

def circular_mask(radius, imagesize=[59,59], centre=[29,29]):
    x = np.arange(imagesize[0])
    y = np.arange(imagesize[1])
    xgrid, ygrid = np.meshgrid(x, y)
    dist = np.sqrt((ygrid - centre[1])**2+(xgrid - centre[0])**2)
    mask = np.ones(imagesize)
    mask[dist <= radius] = 0
    mask = mask.astype(bool)
    return mask

semesters = ['05B','06B', '07B', '08B', '09B', '10B', '11B', '12B']
hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM
r = np.arange(0,42,1) * const * 3600 #define radius values
#centre = [63,63]
centre = [29,29]

#sdata = fits.open('mag_flux_tables/stars_mag_flux_table_matched.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_1518bkgcorr.fits')[1].data

#set up time variable for plot
t = np.linspace(1, 8, num=8)

ids = {}
mags = [18]#[17,18,19,20]#16,
aper_sizes = np.zeros(5) +4#[4,3.5,3,2.5,2]
lower = 15

for m, upper in enumerate(mags):
    aper = aper_sizes[m]
    for sem in semesters:
        # limit data
        mag = sdata['MAG_APER_'+sem][:,4]
        mask1 = mag > lower #removes saturated
        mask2 = mag < upper #removes very faint stars
        mask = mask1 * mask2
        if sem == '05B':
            ids = sdata['NUMBER_05B'][mask]
        else:
            ids = np.intersect1d(ids, sdata['NUMBER_05B'][mask])
    
    mask = np.isin(sdata['NUMBER_05B'], ids)
    tempsdata = sdata[mask] 
    print(len(tempsdata['MAG_APER_'+sem][:,4]))
    
    ### get average FWHM ###
    avgFWHM1 = np.zeros(len(semesters))
    for n, sem in enumerate(semesters):
        avgFWHM1[n] = np.nanmedian(tempsdata['FWHM_WORLD_'+sem]) * 3600
    
    ### get average flux ### 
    avgflux1 = get_avg_flux(tempsdata)
    
    ## get psfs, aper flux, and profiles ###
    pixelr = (aper/3600) / const
    rflux = (1.5/3600) / const
    aperture = CircularAperture(centre, rflux)
    smallaperture = CircularAperture(centre, pixelr)
    psf = {}
    rp = {}
    sqrtrp = {}
    aperflux = np.empty(8)
    newaperflux = np.empty(8)
    for n, sem in enumerate(semesters):
        psf[sem], rp[sem], sqrtrp[sem] = psf_and_profile_mag(upper,lower, sem)
#            psf[sem], rp[sem], sqrtrp[sem] = psf_and_profile(sem)
        
        ### Determine flux within 3 arcsec apertures ###
        phot = aperture_photometry(psf[sem], aperture)
        aperflux[n] = phot['aperture_sum'][0]
        
        ### Mask out source ### 
        mask = circular_mask(pixelr)
        bkg = np.copy(psf[sem])
        bkg[~mask] = np.nan

        ### Plot the psfs ###
        plt.figure(m, figsize=[9,6])
        plt.subplot(2,4,n+1)
        plt.imshow(np.log(psf[sem]), vmax=-4.0, vmin=-20)
        vari_funcs.no_ticks()
        plt.title(sem+' '+str(lower)+'<M<'+str(upper))
        smallaperture.plot()
        
#        ### Plot masked PSFs ###
#        plt.figure(m+5, figsize=[9,6])
#        plt.subplot(2,4,n+1)
##        plt.imshow(mask)
#        plt.imshow(np.log(bkg), vmax=-4.0, vmin=-20)
#        vari_funcs.no_ticks()
#        plt.title(sem+' '+str(lower)+'<M<'+str(upper))
#        smallaperture.plot()
        
        ### Find average value of background ###
        avg = np.nanmean(bkg)
        avg2 = np.nanmedian(bkg)
        print(sem+' '+str(lower)+'<M<'+str(upper))
        print('mean = '+str(avg))
        print('median = '+str(avg2))
        
#        bkg = bkg/avg
#        print(np.nanmean(bkg))
        new = np.copy(psf[sem])
        new = (psf[sem]-avg)
        new /= np.nansum(new)
        print(np.nansum(new))
        newbkg = np.copy(new)
        newbkg[~mask] = np.nan
        print('new bkg avg =', np.nanmean(newbkg))
        
        plt.figure(m+10, figsize=[9,6])
        plt.subplot(2,4,n+1)
#        plt.imshow(mask)
        plt.imshow(np.log(new))#, vmax=-4.0, vmin=-20)
        vari_funcs.no_ticks()
        plt.title(sem+' '+str(lower)+'<M<'+str(upper))
        smallaperture.plot()
        
        
        phot = aperture_photometry(new, aperture)
        newaperflux[n] = phot['aperture_sum'][0]
        print(newaperflux[n])
        
        ### Save corrected PSF ###
        hdu = fits.PrimaryHDU(new)
        hdul = fits.HDUList([hdu])
        hdul.writeto('PSFs/mag_changes/small/'+str(lower)+str(upper)+'bkgcorr_'+sem+'_K_PSF_corrected.fits')
        
    ### Plot aper flux curves ###
    plt.figure(20, figsize=[9,6])
    plt.plot(t,newaperflux,'o-',label=str(lower)+'<M<'+str(upper))
    plt.ylabel('Aperture Flux of PSF')
#    plt.ylim(ymax=0.965, ymin=0.944)
    plt.xticks(t, semesters)
    plt.xlabel('Semester')
    plt.tight_layout()
    plt.legend()
    plt.title('3 arcsec aperture')
#        plt.ylim(ymax=1.31, ymin=0.16)
#    lower = upper