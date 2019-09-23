#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 16:05:33 2018

Code to look at psf data of individual quadrants

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
#    if upper == 20:
#        psf = fits.getdata('PSFs/small_'+sem+'_K_PSF.fits')
#    else:
#        psf = fits.getdata('PSFs/matched_'+sem+'_K_PSF.fits')
    psf = fits.getdata('PSFs/mag_changes/small/'+str(lower)+str(upper)+'_'+sem+'_K_PSF_corrected.fits')
    rp = vari_funcs.radial_profile(psf, centre)
    return psf, rp, np.sqrt(rp)

semesters = ['05B','06B', '07B', '08B', '09B', '10B', '11B', '12B']
hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM
r = np.arange(0,42,1) * const * 3600 #define radius values
#centre = [63,63]
centre = [29,29]

#sdata = fits.open('mag_flux_tables/stars_mag_flux_table_matched.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
#set up time variable for plot
t = np.linspace(1, 8, num=8)

ids = {}
mags = [16,17,18,19,20]
aper_sizes = [1.5]#[0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4]

for k, aper in enumerate(aper_sizes):
    lower = 15
    for m, upper in enumerate(mags):
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
        aperture = CircularAperture(centre, pixelr)
        smallaperture = CircularAperture(centre, pixelr)
        psf = {}
        rp = {}
        sqrtrp = {}
        aperflux = np.empty(8)
        for n, sem in enumerate(semesters):
            psf[sem], rp[sem], sqrtrp[sem] = psf_and_profile_mag(upper,lower, sem)
#            psf[sem], rp[sem], sqrtrp[sem] = psf_and_profile(sem)
            ### Determine flux within 3 arcsec apertures ###
            phot = aperture_photometry(psf[sem], aperture)
            aperflux[n] = phot['aperture_sum'][0]
#            plt.figure(m+8, figsize=[9,6])
            if k==0:
                ### Plot the psfs ###
#                plt.subplot(2,4,n+1)
                plt.figure(5, figsize=[9,6])
                plt.subplot(2,4,m+1)
                plt.imshow(np.log(psf[sem]), vmax=-4.0, vmin=-20)
                vari_funcs.no_ticks()
                plt.title(sem+' '+str(lower)+'<M<'+str(upper))
                smallaperture.plot()
            else:
                plt.subplot(2,4,m+1)
                smallaperture.plot()
##            plt.subplot(241)
#            smallaperture.plot()
#        plt.savefig(str(lower)+str(upper)+'cleanedpsfs.png')        
#    
        ### Plot FWHM curves ###
        plt.figure(1, figsize=[9,6])
        plt.plot(t,avgFWHM1,'o-', label=str(lower)+'<M<'+str(upper))
        plt.ylim(ymax=0.90, ymin=0.70)
        plt.xticks(t, semesters)
        plt.ylabel('FWHM')
        plt.xlabel('Semester')
        plt.legend()
        plt.tight_layout()
    #    
        ### Plot median flux curves ###
        plt.figure(2, figsize=[9,6])
        plt.plot(t,avgflux1,'o-', label=str(lower)+'<M<'+str(upper))
        plt.ylabel('Median Flux of stars')
        plt.ylim(ymax=1.015, ymin=0.985)
        plt.xticks(t, semesters)
        plt.xlabel('Semester')
        plt.legend()
        plt.tight_layout()
        
#        ### Plot radial profiles ###
#        plt.figure(3, figsize=[12,9])
##        plt.subplot(3,2,m+1)
#        for sem in sqrtrp:
#            plt.plot(r, rp[sem], label=str(lower)+'<M<'+str(upper))#sem)
#            plt.xlabel('Radius (arcsec)')
#            plt.ylabel('Flux')
##            plt.ylim(ymax=0.2, ymin=0)
#            plt.xlim(xmin=0, xmax=6)
#            plt.legend()
##        plt.title(str(lower)+'<M<'+str(upper))
#        plt.tight_layout()
#        
        ### Plot aper flux curves ###
        plt.figure(4, figsize=[9,6])
        plt.plot(t,aperflux,'o-',label=str(lower)+'<M<'+str(upper))
        plt.ylabel('Aperture Flux of PSF')
    #    plt.ylim(ymax=0.965, ymin=0.944)
        plt.xticks(t, semesters)
        plt.xlabel('Semester')
        plt.tight_layout()
        plt.legend()
        plt.title(str(2*aper)+' arcsec aperture')
#        plt.ylim(ymax=1.31, ymin=0.16)
        lower = upper
    
    #plt.figure(1)
    #plt.savefig('magchangesFWHM.png')
    #plt.figure(2)
    #plt.savefig('magchangesmedianflux.png')
    #plt.figure(3)
    #plt.savefig('magchangessqrtrp.png')
#    plt.figure(k)
#    plt.savefig('Changing_mag_range/set_axes_size_changes/cleanedmagchangesaper'+str(int(20*aper))+'fluxlim.png')
    
