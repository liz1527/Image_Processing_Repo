#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 12:16:06 2018

Code to investigate the behaviour of a single star through the semesters

Designed to test if too much averaging = mistakes introduced.

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
import sys
plt.close('all')
hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

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

def get_stamp(data, mask):
    stamp = data['VIGNET'][mask]
    stamp[stamp<-5e29] = np.nan
    return np.reshape(stamp, [59,59])

def plot_rp(stamp,sem=''):
    # set up coordinates
    centre = [29,29]
    hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
    const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM
    r = np.arange(0,42,1) * const * 3600 #define radius values
    
    rp = radial_profile(stamp, centre)
    sqrt = np.sqrt(rp)
    plt.plot(r, sqrt, label=sem)
    plt.xlabel('Radius (arcsec)')
    plt.ylabel('sqrt(Flux)')
    return rp

def get_psf(sem):
    return fits.getdata('PSFs/'+sem+'_K_PSF.fits')

### Read in all of the stamp tables and star table ###
data05 = fits.getdata('star_stamps_tables/small_05B_star_stamps_table.fits')
data06 = fits.getdata('star_stamps_tables/small_06B_star_stamps_table.fits')
data07 = fits.getdata('star_stamps_tables/small_07B_star_stamps_table.fits')
data08 = fits.getdata('star_stamps_tables/small_08B_star_stamps_table.fits')
data09 = fits.getdata('star_stamps_tables/small_09B_star_stamps_table.fits')
data10 = fits.getdata('star_stamps_tables/small_10B_star_stamps_table.fits')
data11 = fits.getdata('star_stamps_tables/small_11B_star_stamps_table.fits')
data12 = fits.getdata('star_stamps_tables/small_12B_star_stamps_table.fits')
sdata = fits.getdata('mag_flux_tables/stars_mag_flux_table.fits')
#psf05 = get_psf('05B')
#psf06 = get_psf('06B')
#psf07 = get_psf('07B')
#psf08 = get_psf('08B')
#psf09 = get_psf('09B')
#psf10 = get_psf('10B')
#psf11 = get_psf('11B')
#psf12 = get_psf('12B')

### Pick a star for analysis ###

obnum = 3322 #chosen at random for now
if ~np.isin(obnum, data05['NUMBER']):
    sys.exit('Invalid object number') # quits script if object number invalid
mask = data05['NUMBER'] == obnum

### extract flux curve ###
flux = np.array([data05['FLUX_APER'][:,4][mask],data06['FLUX_APER'][:,4][mask],
                data07['FLUX_APER'][:,4][mask],data08['FLUX_APER'][:,4][mask],
                data09['FLUX_APER'][:,4][mask],data10['FLUX_APER'][:,4][mask], 
                data11['FLUX_APER'][:,4][mask],data12['FLUX_APER'][:,4][mask]])
flux = np.reshape(flux, 8)
fluxerr = np.array([data05['FLUXERR_APER'][:,4][mask],data06['FLUXERR_APER'][:,4][mask],
                    data07['FLUXERR_APER'][:,4][mask],data08['FLUXERR_APER'][:,4][mask],
                    data09['FLUXERR_APER'][:,4][mask],data10['FLUXERR_APER'][:,4][mask], 
                    data11['FLUXERR_APER'][:,4][mask],data12['FLUXERR_APER'][:,4][mask]])
fluxerr = np.reshape(fluxerr, 8)

plt.figure(1)
vari_funcs.avg_lightcurve(flux, fluxerr, label='sextractor flux')
plt.title('')
plt.tight_layout()

### extract FWHM curve ###
fwhm = np.array([data05['FWHM_WORLD'][mask],data06['FWHM_WORLD'][mask],
                data07['FWHM_WORLD'][mask],data08['FWHM_WORLD'][mask],
                data09['FWHM_WORLD'][mask],data10['FWHM_WORLD'][mask], 
                data11['FWHM_WORLD'][mask],data12['FWHM_WORLD'][mask]])*3600
fwhm = np.reshape(fwhm, 8)

plt.figure(2)
vari_funcs.avg_lightcurve(fwhm)
plt.ylabel('FWHM (arcsec)')
plt.title('')
plt.tight_layout()
plt.savefig(str(obnum)+'FWHM.png')

### Get stamps and plot them ###

stamp05 = get_stamp(data05, mask)
stamp06 = get_stamp(data06, mask)
stamp07 = get_stamp(data07, mask)
stamp08 = get_stamp(data08, mask)
stamp09 = get_stamp(data09, mask)
stamp10 = get_stamp(data10, mask)
stamp11 = get_stamp(data11, mask)
stamp12 = get_stamp(data12, mask)

plt.figure(3, figsize=[4,9])
#plt.subplot(441)
plt.subplot(421)
plt.imshow(np.log(stamp05))
vari_funcs.no_ticks()
plt.title('05B')
#plt.subplot(443)
plt.subplot(422)
plt.imshow(np.log(stamp06))
vari_funcs.no_ticks()
plt.title('06B')
#plt.subplot(445)
plt.subplot(423)
plt.imshow(np.log(stamp07))
vari_funcs.no_ticks()
plt.title('07B')
#plt.subplot(447)
plt.subplot(424)
plt.imshow(np.log(stamp08))
vari_funcs.no_ticks()
plt.title('08B')
#plt.subplot(449)
plt.subplot(425)
plt.imshow(np.log(stamp09))
vari_funcs.no_ticks()
plt.title('09B')
#plt.subplot(4,4,11)
plt.subplot(426)
plt.imshow(np.log(stamp10))
vari_funcs.no_ticks()
plt.title('10B')
#plt.subplot(4,4,13)
plt.subplot(427)
plt.imshow(np.log(stamp11))
vari_funcs.no_ticks()
plt.title('11B')
#plt.subplot(4,4,15)
plt.subplot(428)
plt.imshow(np.log(stamp12))
vari_funcs.no_ticks()
plt.title('12B')
plt.savefig(str(obnum)+'PSFs.png')

### Run aperture phot to see if flux curve is the same ###
centre=[29,29]
pixelr = (1.5/3600) / const
aperture = CircularAperture(centre, pixelr)
aper05 = aperture_photometry(stamp05, aperture)
aper06 = aperture_photometry(stamp06, aperture)
aper07 = aperture_photometry(stamp07, aperture)
aper08 = aperture_photometry(stamp08, aperture)
aper09 = aperture_photometry(stamp09, aperture)
aper10 = aperture_photometry(stamp10, aperture)
aper11 = aperture_photometry(stamp11, aperture)
aper12 = aperture_photometry(stamp12, aperture)
aperflux = np.array([aper05['aperture_sum'][0],aper06['aperture_sum'][0],
                aper07['aperture_sum'][0],aper08['aperture_sum'][0],
                aper09['aperture_sum'][0],aper10['aperture_sum'][0], 
                aper11['aperture_sum'][0],aper12['aperture_sum'][0]])
plt.figure(1)
vari_funcs.avg_lightcurve(aperflux, label='python flux')
plt.title('')
plt.tight_layout()
plt.legend()
plt.savefig(str(obnum)+'fluxcurve.png')

### Get and display profiles ###

#plt.figure(3)
#plt.subplot(442)
#plot_rp(stamp05)
#plt.subplot(444)
#plot_rp(stamp06)
#plt.subplot(446)
#plot_rp(stamp07)
#plt.subplot(448)
#plot_rp(stamp08)
#plt.subplot(4,4,10)
#plot_rp(stamp09)
#plt.subplot(4,4,12)
#plot_rp(stamp10)
#plt.subplot(4,4,14)
#plot_rp(stamp11)
#plt.subplot(4,4,16)
#plot_rp(stamp12)

plt.figure(4)
plot_rp(stamp05,'05B')
plot_rp(stamp06,'06B')
plot_rp(stamp07,'07B')
plot_rp(stamp08,'08B')
plot_rp(stamp09,'09B')
plot_rp(stamp10,'10B')
plot_rp(stamp11,'11B')
plot_rp(stamp12,'12B')

#plot_rp(psf05,'05B')
#plot_rp(psf06,'06B')
#plot_rp(psf07,'07B')
#plot_rp(psf08,'08B')
#plot_rp(psf09,'09B')
#plot_rp(psf10,'10B')
#plot_rp(psf11,'11B')
#plot_rp(psf12,'12B')
plt.legend()
plt.savefig(str(obnum)+'profiles.png')









