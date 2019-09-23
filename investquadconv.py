#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 16:05:33 2018

Code to compare psf data of individual quadrants before and after convolution

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
    mask2 = idec < middec
    quad1data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec < middec
    quad2data = initdata[mask1*mask2]
    
    mask1 = ira < midra
    mask2 = idec >= middec
    quad3data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec >= middec
    quad4data = initdata[mask1*mask2]
    
    return quad1data, quad2data, quad3data, quad4data

def get_avg_flux(tbdata):
    flux = vari_funcs.flux1_stacks(tbdata)
    flux = vari_funcs.normalise_flux(flux)
    return np.nanmedian(flux, axis=0)

def old_psf_and_profile(quad, sem):
    centre = [29,29]
    psf = fits.getdata('PSFs/Quad_PSFs/'+sem+'_'+str(quad)+'_K_PSF.fits')
#    if sem == '10B':
#        psf = fits.getdata('PSFs/Quad_PSFs/'+sem+'_'+str(quad)+'_K_PSF.fits')
#    else:
#        psf = fits.getdata('PSFs/Quad_PSFs/extraq_'+sem+'_'+str(quad)+'_K_PSF.fits')

    rp = vari_funcs.radial_profile(psf, centre)
    return psf, rp, np.sqrt(rp)

def psf_and_profile(quad, sem):
    centre = [29,29]
    psf = fits.getdata('PSFs/Quad_PSFs/extraq_'+sem+'_'+str(quad)+'_K_PSF.fits')
#    if sem == '10B':
#        psf = fits.getdata('PSFs/Quad_PSFs/'+sem+'_'+str(quad)+'_K_PSF.fits')
#    else:
#        psf = fits.getdata('PSFs/Quad_PSFs/extraq_'+sem+'_'+str(quad)+'_K_PSF.fits')

    rp = vari_funcs.radial_profile(psf, centre)
    return psf, rp, np.sqrt(rp)

semesters = ['05B','06B', '07B', '08B', '09B', '10B', '11B', '12B']
hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM
r = np.arange(0,42,1) * const * 3600 #define radius values
centre = [29,29]

oldsdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_extraq.fits')[1].data
#set up time variable for plot
t = np.linspace(1, 8, num=8)

for n, sem in enumerate(semesters):
    # limit old data
    oldmag = oldsdata['MAG_APER_'+sem][:,4]
    oldmask1 = oldmag > 15 #removes saturated
    oldmask2 = oldmag < 19 #removes very faint stars
    oldmask = oldmask1 * oldmask2
    if sem == '05B':
        oldids = oldsdata['NUMBER_05B'][oldmask]
    else:
        oldids = np.intersect1d(oldids, oldsdata['NUMBER_05B'][oldmask])
    # limit new data
    mag = sdata['MAG_APER_'+sem][:,4]
    mask1 = mag > 15 #removes saturated
    mask2 = mag < 19 #removes very faint stars
    mask = mask1 * mask2
    if sem == '05B':
        ids = sdata['NUMBER_05B'][mask]
    else:
        ids = np.intersect1d(ids, sdata['NUMBER_05B'][mask])

# mask old data
oldmask = np.isin(oldsdata['NUMBER_05B'], oldids)
tempoldsdata = oldsdata[oldmask] 
print(len(tempoldsdata['MAG_APER_'+sem][:,4]))
#oldsquad1data, oldsquad2data, oldsquad3data, oldsquad4data = quadrants(tempoldsdata,'05B')
oldsquaddata = quadrants(tempoldsdata,'05B')

#mask new data
mask = np.isin(sdata['NUMBER_05B'], ids)
tempsdata = sdata[mask] 
print(len(tempsdata['MAG_APER_'+sem][:,4]))
#squad1data, squad2data, squad3data, squad4data = quadrants(tempsdata,'05B')
squaddata = quadrants(tempsdata,'05B')

### get average FWHM ###
oldavgFWHM = np.zeros([4, len(semesters)])
avgFWHM = np.zeros([4, len(semesters)])
oldavgflux = np.zeros([4, len(semesters)])
avgflux = np.zeros([4, len(semesters)])
for m in range(4):
    oldavgflux[m,:] = get_avg_flux(oldsquaddata[m])
    avgflux[m,:] = get_avg_flux(squaddata[m])
    for n, sem in enumerate(semesters):
        oldavgFWHM[m,n] = np.nanmedian(oldsquaddata[m]['FWHM_WORLD_'+sem]) * 3600
        avgFWHM[m,n] = np.nanmedian(squaddata[m]['FWHM_WORLD_'+sem]) * 3600
        

### get psfs, aper flux, and profiles ###
pixelr = (0.25/3600) / const
aperture = CircularAperture(centre, pixelr)
smallaperture = CircularAperture(centre, pixelr)
oldpsf = {}
oldrp = {}
oldsqrtrp = {}
oldaperflux = {1:np.empty(8),
            2:np.empty(8),
            3:np.empty(8),
            4:np.empty(8)}
psf = {}
rp = {}
sqrtrp = {}
aperflux = {1:np.empty(8),
            2:np.empty(8),
            3:np.empty(8),
            4:np.empty(8)}
for m, sem in enumerate(semesters):
    oldpsf[sem] = {}
    oldrp[sem] = {}
    oldsqrtrp[sem] ={}
    psf[sem] = {}
    rp[sem] = {}
    sqrtrp[sem] ={}
    for n in [1,2,3,4]:
        oldpsf[sem][n], oldrp[sem][n], oldsqrtrp[sem][n] = old_psf_and_profile(n,sem)
        psf[sem][n], rp[sem][n], sqrtrp[sem][n] = psf_and_profile(n,sem)
        ### Determine flux within 3 arcsec apertures ###
        oldphot = aperture_photometry(oldpsf[sem][n], aperture)
        oldaperflux[n][m] = oldphot['aperture_sum'][0]
        phot = aperture_photometry(psf[sem][n], aperture)
        aperflux[n][m] = phot['aperture_sum'][0]
#        ### Plot the psfs ###
#        plt.figure(n+4)
#        plt.subplot(2,4,m+1)
#        plt.imshow(np.log(psf[sem][n]), vmax=-4.0, vmin=-20)
#        vari_funcs.no_ticks()
#        plt.title(sem+str(n))
#        print(np.nanmin(np.log(psf[sem][n])))
        
for m in range(4):
    ### Plot FWHM curves ###
    plt.figure(1, figsize=[9,6])
    plt.subplot(2,2,m+1)
    plt.plot(t,oldavgFWHM[m,:],'o')
    plt.plot(t,avgFWHM[m,:],'o')
    #plt.ylim(ymax=0.90, ymin=0.70)
    plt.xticks(t, semesters)
    plt.ylabel('FWHM')
    plt.xlabel('Semester')

    ### Plot median flux curves ###
    plt.figure(2, figsize=[9,6])
    plt.subplot(2,2,m+1)
    plt.plot(t,oldavgflux[m,:],'o')
    plt.plot(t,avgflux[m,:],'o')
    plt.ylabel('Median Flux of stars')
    #plt.ylim(ymax=1.015, ymin=0.985)
    plt.xticks(t, semesters)
    plt.xlabel('Semester')

quad = [1,2,3,4]
for n,sem in enumerate(semesters):
    plt.figure(3, figsize=[9,9])
    plt.subplot(4,2,n+1)
    plt.plot(quad, oldavgFWHM[:,n],'o')
    plt.plot(quad, avgFWHM[:,n],'o')
    plt.ylabel('FWHM')
    plt.xlabel('Quadrant '+sem)
#    plt.ylim(ymin=0.7, ymax=1)
    plt.tight_layout()

    plt.figure(4, figsize=[9,9])
    plt.subplot(4,2,n+1)
    plt.plot(quad, oldavgflux[:,n],'o')
    plt.plot(quad, avgflux[:,n],'o')
    plt.ylabel('0.7" Flux')
    plt.xlabel('Quadrant '+sem)
    plt.tight_layout()
# Plot radial profiles ###
for m, sem in enumerate(sqrtrp):
    for n in sqrtrp[sem]:        
        plt.figure(5, figsize=[12,9])
        plt.subplot(2,4,m+1)
        plt.plot(r, oldsqrtrp[sem][n],'--', label='old '+sem+str(n))
        plt.xlabel('Radius (arcsec)')
        plt.ylabel('sqrt(Flux)')
        plt.ylim(ymax=0.16, ymin=0)
        plt.legend()
        plt.tight_layout()
        
        plt.figure(6, figsize=[12,9])
        plt.subplot(2,4,m+1)
        plt.plot(r, sqrtrp[sem][n], label=sem+str(n))
        plt.xlabel('Radius (arcsec)')
        plt.ylabel('sqrt(Flux)')
        plt.ylim(ymax=0.16, ymin=0)
        plt.legend()
        plt.tight_layout()

### Plot aper flux curves ###
plt.figure(7, figsize=[9,6])
for n in [1,2,3,4]:
    plt.subplot(2,2,n)
    plt.plot(t,oldaperflux[n],'o')
    plt.plot(t,aperflux[n],'o')
    plt.ylabel('Aperture Flux of PSF')
#    plt.ylim(ymax=0.965, ymin=0.944)
    plt.xticks(t, semesters)
    plt.xlabel('Semester')
    plt.tight_layout()
