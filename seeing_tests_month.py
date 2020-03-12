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
plt.style.use('default')

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

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
    flux = vari_funcs.k_mag_flux.flux5_stacks(tbdata)
    flux = vari_funcs.flux_funcs.normalise_flux(flux)
    return np.nanmedian(flux, axis=0)

def psf_and_profile(sem):
    centre = [29,29]
#    psf = fits.getdata('PSFs/cleaned_'+sem+'_K_PSF.fits')
    if sem == 'dec06':
        psf = fits.getdata('PSFs/K/month/cleaned_'+sem+'_K_PSF.fits')
    else:
        psf = fits.getdata('PSFs/K/month/extra_clean_'+sem+'_K_PSF.fits')

    rp = vari_funcs.correction_funcs.radial_profile(psf, centre)
    return psf, rp, np.sqrt(rp)

def psf_and_profile_old(sem):
    centre = [29,29]
    psf = fits.getdata('PSFs/K/month/cleaned_'+sem+'_K_PSF.fits')
#    if sem == '10B':
#        psf = fits.getdata('PSFs/small_'+sem+'_K_PSF.fits')
#    else:
#        psf = fits.getdata('PSFs/matched_'+sem+'_K_PSF.fits')
    rp = vari_funcs.correction_funcs.radial_profile(psf, centre)
    return psf, rp, np.sqrt(rp)

#semesters = ['05B','07B', '08B', '09B', '10B', '11B', '12B']#'06B', 
months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07',  
          'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09',  
          'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 
          'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 
          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
hdr08B = fits.getheader('Images/UDS-DR11-K.mef.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM
r = np.arange(0,42,1) * const * 3600 #define radius values
centre = [29,29] 

psf_data = fits.open('UDS_catalogues/DR11_stars_for_PSFs.fits')[1].data
sdata = fits.open('mag_flux_tables/K/month/month_stars_mag_flux_table_K_extra_clean.fits')[1].data
sdataold = fits.open('mag_flux_tables/K/month/month_stars_mag_flux_table_K_cleaned.fits')[1].data

### set up month tick details ###
month_info = fits.open('monthly_numbers.fits')[1].data #get month count data
full_months = month_info['Month'] #extract month nanes
tick_inds = np.load('tick_inds_K.npy') #load tick locations
mask = np.zeros(len(full_months)) #set up mask
mask[tick_inds] = 1
mask = mask.astype(bool)
month_ticks = np.copy(full_months)
month_ticks = month_ticks[mask]#retrieve tick details

x = np.arange(0, len(month_info['Frames in v11']))
mask = np.isin(full_months, months)
x_months = x[mask]

for sem in months:
    # limit data
    ### Define coordinates ###
    refcoord = SkyCoord(psf_data['ALPHA_J2000_1']*u.degree, psf_data['DELTA_J2000_1']*u.degree)
    semcoord = SkyCoord(sdata['ALPHA_J2000_'+sem]*u.degree, sdata['DELTA_J2000_'+sem]*u.degree)
    
    ### Match catalogues and create new table ###
    idx, d2d , _ = match_coordinates_sky(refcoord, semcoord) #match these 'good' stars to create table

#    mag = sdata['MAG_APER_'+sem][:,4]
#    mask1 = mag > 15 #removes saturated
#    mask2 = mag < 19 #removes very faint stars
#    mask = mask1 * mask2
    if sem == 'sep05':
        ids = sdata['NUMBER'][idx]
    else:
        ids = np.intersect1d(ids, sdata['NUMBER'][idx])

mask = np.isin(sdata['NUMBER'], ids)
oldmask = np.isin(sdataold['NUMBER'], ids)
tempsdata = sdata[mask] 
tempsdataold = sdataold[mask]
print(len(tempsdata['MAG_APER_'+sem][:,4]))

### get average FWHM ###
avgFWHM1 = np.zeros(len(months))
avgFWHM2 = np.zeros(len(months))
for n, sem in enumerate(months):
    tempsdata['FWHM_WORLD_'+sem][tempsdata['FWHM_WORLD_'+sem]==0] = np.nan
    avgFWHM1[n] = np.nanmedian(tempsdata['FWHM_WORLD_'+sem]) * 3600
    avgFWHM2[n] = np.nanmedian(tempsdataold['FWHM_WORLD_'+sem]) * 3600

#### get average flux ### 
#avgflux1 = get_avg_flux(tempsdata)
#avgflux2 = get_avg_flux(tempsdataold)

### get psfs, aper flux, and profiles ###
pixelr = (1.5/3600) / const
aperture = CircularAperture(centre, pixelr)
smallaperture = CircularAperture(centre, pixelr)
psf = {}
rp = {}
sqrtrp = {}
aperflux = np.empty(len(months))
psfold = {}
rpold = {}
sqrtrpold = {}
aperfluxold = np.empty(len(months))
for m, sem in enumerate(months):
    psf[sem], rp[sem], sqrtrp[sem] = psf_and_profile(sem)
    psfold[sem], rpold[sem], sqrtrpold[sem] = psf_and_profile_old(sem)
    ### Determine flux within 3 arcsec apertures ###
    phot = aperture_photometry(psf[sem], aperture)
    photold = aperture_photometry(psfold[sem], aperture)
    aperflux[m] = phot['aperture_sum'][0]
    aperfluxold[m] = photold['aperture_sum'][0]
    ### Plot the psfs ###
#    plt.figure(5, figsize=[10,2])
#    plt.subplot(1,7,m+1)
#    plt.imshow(np.log(psf[sem]), vmax=-4.0, vmin=-20)
#    vari_funcs.no_ticks()
#    plt.title(sem)
#    plt.tight_layout()
#    plt.figure(7, figsize=[10,2])
#    plt.subplot(1,7,m+1)
#    plt.imshow(np.log(psfold[sem]), vmax=-4.0, vmin=-20)
#    vari_funcs.no_ticks()
#    plt.title(sem)
#    plt.tight_layout()
        
### Plot FWHM curves ###
plt.figure(1, figsize=[9,6])
plt.plot(x_months,avgFWHM2,'s-',label='Before', markersize=8)
plt.plot(x_months,avgFWHM1,'o-',label='After')
#plt.ylim(ymax=0.87, ymin=0.70)
plt.xticks(tick_inds, month_ticks, rotation = 'vertical')
plt.ylabel('FWHM')
plt.xlabel('Month')
plt.legend()
plt.tight_layout()

#### Plot median flux curves ###
#plt.figure(2, figsize=[9,6])
#plt.plot(x_months,avgflux2,'s-',label='old', markersize=8)
#plt.plot(x_months,avgflux1,'o-',label='new')
#plt.ylabel('Median Flux of stars')
##plt.ylim(ymax=1.015, ymin=0.985)
#plt.xticks(tick_inds, month_ticks, rotation = 'vertical')
#plt.xlabel('Month')
#plt.legend()
#plt.tight_layout()
#
### Plot radial profiles ###
plt.figure(3, figsize=[10,7])
base = rp['dec06']
baseold = rpold['dec06']
for sem in sqrtrp:
    plt.figure(3)
    plt.subplot(212)
    plt.plot(r, sqrtrp[sem], label=sem+' after')
    plt.xlabel('Radius (arcsec)')
    plt.ylabel('sqrt(Flux)')
    plt.ylim(ymax=0.17, ymin=0)
    plt.xlim(xmax=2, xmin=0)
#    plt.legend()
    plt.subplot(211)
    plt.plot(r, sqrtrpold[sem],'--', label=sem+' before')
    plt.xlabel('Radius (arcsec)')
    plt.ylabel('sqrt(Flux)')
    plt.ylim(ymax=0.17, ymin=0)
    plt.xlim(xmax=2, xmin=0)
#    plt.legend()
    plt.tight_layout()
    
    ### figure out differences ###
    diff = rp[sem] - base
    diffold = rpold[sem] - baseold
    
    ### Plot differences ###
    plt.figure(4, figsize=[8,7])
    plt.subplot(212)
    plt.plot(r, diff, label=sem+' after')
    plt.xlabel('Radius (arcsec)')
    plt.ylabel('$Flux_{sem} - Flux_{10B}$')
    plt.ylim(ymax=0.011, ymin=-0.002)
    plt.xlim(xmax=2, xmin=0)
#    plt.legend(loc='upper right')
    plt.subplot(211)
    plt.plot(r, diffold,'--', label=sem+' before')
    plt.xlabel('Radius (arcsec)')
    plt.ylabel('$Flux_{sem} - Flux_{10B}$')
    plt.ylim(ymax=0.011, ymin=-0.002)
    plt.xlim(xmax=2, xmin=0)
#    plt.legend(loc='upper right')
    plt.tight_layout()

#
#### Plot aper flux curves ###
#plt.figure(6, figsize=[9,6])
#plt.plot(x,aperfluxold,'s-',label='old', markersize=8)
#plt.plot(x,aperflux,'o-',label='new')
#plt.ylabel('Aperture Flux of PSF')
##plt.ylim(ymax=0.965, ymin=0.944)
#plt.xticks(t, years)
#plt.xlabel('Semester')
#plt.legend()
#plt.tight_layout()
#
