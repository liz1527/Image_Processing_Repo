#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 16:05:33 2018

Code to look at psf data of individual quadrants

@author: ppxee
"""


### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from photutils import CircularAperture, aperture_photometry
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
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

def get_psf(sem):
    return fits.getdata('PSFs/'+sem+'_K_PSF.fits')

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

def get_avg_flux(tbdata):
    flux = vari_funcs.hflux4_stacks(tbdata)
    flux = vari_funcs.normalise_flux(flux)
    return np.nanmedian(flux, axis=0)

def psf_and_profile(quad, sem):
    centre = [29,29]
    psf = fits.getdata('PSFs/H/Quad_PSFs/cleaned_Kstars_'+sem+'_'+str(quad)+'_H_PSF.fits')
#    if sem == '10B':
#        psf = fits.getdata('PSFs/Quad_PSFs/'+sem+'_'+str(quad)+'_K_PSF.fits')
#    else:
#        psf = fits.getdata('PSFs/Quad_PSFs/extraq_'+sem+'_'+str(quad)+'_K_PSF.fits')

    rp = vari_funcs.radial_profile(psf, centre)
    return psf, rp, np.sqrt(rp)

semesters = ['06B', '07B', '08B', '09B', '10B', '11B', '12B']
hdr08B = fits.getheader('Images/UDS-DR11-K.mef.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM
r = np.arange(0,42,1) * const * 3600 #define radius values
centre = [29,29]

psf_data = fits.open('UDS_catalogues/DR11_stars_for_PSFs.fits')[1].data
sdata = fits.open('mag_flux_tables/H/stars_mag_flux_table_H_cleaned.fits')[1].data
#set up time variable for plot
t = np.linspace(2, 8, num=7)

for n, sem in enumerate(semesters):### Define coordinates ###
    
    refcoord = SkyCoord(psf_data['ALPHA_J2000_1']*u.degree, psf_data['DELTA_J2000_1']*u.degree)
    semcoord = SkyCoord(sdata['ALPHA_J2000_'+sem]*u.degree, sdata['DELTA_J2000_'+sem]*u.degree)
    
    ### Match catalogues and create new table ###
    idx, d2d , _ = match_coordinates_sky(refcoord, semcoord) #match these 'good' stars to create table

    if sem == '06B':
        ids = sdata['NUMBER_06B'][idx]
    else:
        ids = np.intersect1d(ids, sdata['NUMBER_06B'][idx])

mask = np.isin(sdata['NUMBER_06B'], ids)
tempsdata = sdata[mask] 
print(len(tempsdata['MAG_APER_'+sem][:,4]))
squad1data, squad2data, squad3data, squad4data = quadrants(tempsdata,'06B')

### get average FWHM ###
avgFWHM1 = np.zeros(len(semesters))
avgFWHM2 = np.zeros(len(semesters))
avgFWHM3 = np.zeros(len(semesters))
avgFWHM4 = np.zeros(len(semesters))
for n, sem in enumerate(semesters):
    avgFWHM1[n] = np.nanmedian(squad1data['FWHM_WORLD_'+sem]) * 3600
    avgFWHM2[n] = np.nanmedian(squad2data['FWHM_WORLD_'+sem]) * 3600
    avgFWHM3[n] = np.nanmedian(squad3data['FWHM_WORLD_'+sem]) * 3600
    avgFWHM4[n] = np.nanmedian(squad4data['FWHM_WORLD_'+sem]) * 3600

### get average flux ### 
avgflux1 = get_avg_flux(squad1data)
avgflux2 = get_avg_flux(squad2data)
avgflux3 = get_avg_flux(squad3data)
avgflux4 = get_avg_flux(squad4data)

## get psfs, aper flux, and profiles ###
pixelr = (1.5/3600) / const
aperture = CircularAperture(centre, pixelr)
smallaperture = CircularAperture(centre, pixelr)
psf = {}
rp = {}
sqrtrp = {}
aperflux = {1:np.empty(len(semesters)),
            2:np.empty(len(semesters)),
            3:np.empty(len(semesters)),
            4:np.empty(len(semesters))}
for m, sem in enumerate(semesters):
    psf[sem] = {}
    rp[sem] = {}
    sqrtrp[sem] ={}
    for n in [1,2,3,4]:
        psf[sem][n], rp[sem][n], sqrtrp[sem][n] = psf_and_profile(n,sem)
        ### Determine flux within 3 arcsec apertures ###
        phot = aperture_photometry(psf[sem][n], aperture)
        aperflux[n][m] = phot['aperture_sum'][0]
        ### Plot the psfs ###
        plt.figure(m+5)
        plt.subplot(2,2,n)
        plt.imshow(np.log(psf[sem][n]), vmax=-4.0, vmin=-20)
        vari_funcs.no_ticks()
        plt.title(sem+str(n))
#        print(np.nanmin(np.log(psf[sem][n])))
        
### Plot FWHM curves ###
plt.figure(1, figsize=[9,6])
plt.subplot(221)
plt.plot(t,avgFWHM1,'o')
plt.ylim(ymax=0.95, ymin=0.780)
plt.xticks(t, semesters)
plt.ylabel('FWHM')
plt.xlabel('Semester')

plt.subplot(222)
plt.plot(t,avgFWHM2,'o')
plt.ylim(ymax=0.95, ymin=0.780)
plt.xticks(t, semesters)
plt.ylabel('FWHM')
plt.xlabel('Semester')

plt.subplot(223)
plt.plot(t,avgFWHM3,'o')
plt.ylim(ymax=0.95, ymin=0.780)
plt.xticks(t, semesters)
plt.ylabel('FWHM')
plt.xlabel('Semester')

plt.subplot(224)
plt.plot(t,avgFWHM4,'o')
plt.ylim(ymax=0.95, ymin=0.780)
plt.xticks(t, semesters)
plt.ylabel('FWHM')
plt.xlabel('Semester')
plt.tight_layout()

### Plot median flux curves ###
plt.figure(2, figsize=[9,6])
plt.subplot(221)
plt.plot(t,avgflux1,'o')
plt.ylabel('Median Flux of stars')
plt.ylim(ymax=1.03, ymin=0.97)
plt.xticks(t, semesters)
plt.xlabel('Semester')

plt.subplot(222)
plt.plot(t,avgflux2,'o')
plt.ylabel('Median Flux of stars')
plt.ylim(ymax=1.03, ymin=0.97)
plt.xticks(t, semesters)
plt.xlabel('Semester')

plt.subplot(223)
plt.plot(t,avgflux3,'o')
plt.ylabel('Median Flux of stars')
plt.ylim(ymax=1.03, ymin=0.97)
plt.xticks(t, semesters)
plt.xlabel('Semester')

plt.subplot(224)
plt.plot(t,avgflux4,'o')
plt.ylabel('Median Flux of stars')
plt.ylim(ymax=1.03, ymin=0.97)
plt.xticks(t, semesters)
plt.xlabel('Semester')
plt.tight_layout()

### Plot radial profiles ###
plt.figure(3, figsize=[12,9])
for sem in sqrtrp:
    for n in sqrtrp[sem]:
        plt.subplot(2,2,n)
        plt.plot(r, sqrtrp[sem][n], label=sem)
        plt.xlabel('Radius (arcsec)')
        plt.ylabel('sqrt(Flux)')
        plt.ylim(ymax=0.16, ymin=0)
        plt.legend()
plt.tight_layout()

### Plot aper flux curves ###
plt.figure(4, figsize=[9,6])
for n in [1,2,3,4]:
    plt.subplot(2,2,n)
    plt.plot(t,aperflux[n],'o')
    plt.ylabel('Aperture Flux of PSF')
    plt.ylim(ymax=0.96, ymin=0.925)
    plt.xticks(t, semesters)
    plt.xlabel('Semester')
    plt.tight_layout()
