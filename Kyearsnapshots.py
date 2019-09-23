#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 17:23:43 2018

@author: ppxee
"""


### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
import matplotlib.animation as animation
#plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
from photutils import CircularAperture
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots
hdr08B = fits.getheader('J/UDS_08B_J.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

obnum=202081
### Read in fits files ###
#obdata = fits.getdata('mag_flux_tables/mag_flux_table_best.fits', 1)
obdata = fits.getdata('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits', 1)
obdata = obdata[obdata['NUMBER_05B']==obnum] #restrict to just 252446
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

#flux = np.array([obdata['FLUX_APER_05B'][:,4],#obdata['FLUX_APER_06B'][:,4],
#                obdata['FLUX_APER_07B'][:,4],obdata['FLUX_APER_08B'][:,4],
#                obdata['FLUX_APER_09B'][:,4],obdata['FLUX_APER_10B'][:,4], 
#                obdata['FLUX_APER_11B'][:,4],obdata['FLUX_APER_12B'][:,4]])
flux, fluxerr, obdata = vari_funcs.create_quad_error_array(sigtb, obdata)
flux = np.reshape(flux, 7)
fluxerr = np.reshape(fluxerr, 7)
semesters = ['05B', '07B', '08B', '09B', '10B', '11B', '12B']

n=1
fig = plt.figure(figsize=[10,5])
snaps = []
for sem in semesters:
#    print(sem)
    if sem == '10B':
        imdata = fits.getdata('cleaned_UDS_'+sem+'_K.fits')
    else:
        imdata = fits.getdata('extra_clean_no06_UDS_'+sem+'_K.fits')
#    
#    imdata = fits.getdata('extra_clean_no06_UDS_'+sem+'_K.fits')
    ### Find coordinates of objects ###
    x = obdata['X_IMAGE_05B']#+sem]
    x = x.astype(int)#int(x)
    y = obdata['Y_IMAGE_05B']#'+sem]
    y = y.astype(int)#int(x) 
    
    size = 100 # size of square
    newim = imdata[y[0]-size:y[0]+size,x[0]-size:x[0]+size]
    
    del imdata
    print(newim[size,size])
#    
#    imuppthresh = 50
#    newim[newim>imuppthresh] = imuppthresh
#    imlowthresh = 0 
#    newim[newim<imlowthresh] = imlowthresh
#    
##    plt.subplot(4, 2, n)
#    snap = plt.imshow(newim, animated=True)
##    plt.plot(size, size, 'k+', markersize=10, mfc='none')
#    plt.xticks([])
#    plt.yticks([])
##    plt.title(sem)
#    
#    snaps.append([snap])
#    
#    n += 1

def func(sem):
    if sem == '10B':
        imdata = fits.getdata('cleaned_UDS_'+sem+'_K.fits')
    else:
        imdata = fits.getdata('extra_clean_no06_UDS_'+sem+'_K.fits')
#    imdata = fits.getdata('extra_clean_no06_UDS_'+sem+'_K.fits')
    ### Find coordinates of objects ###
    x = obdata['X_IMAGE_05B']#'+sem]
    x = x.astype(int)#int(x)
    y = obdata['Y_IMAGE_05B']#'+sem]
    y = y.astype(int)#int(x) 
    
    size = 75 # size of square
    newim = imdata[y[0]-size:y[0]+size,x[0]-size:x[0]+size]
    
    del imdata
    
    imuppthresh = 460
    newim[newim>imuppthresh] = imuppthresh
    imlowthresh = 0
    newim[newim<imlowthresh] = 0#imlowthresh
#    
#    ### restrict stars to good range ###
#    mag = sdata['MAG_APER_'+mon][:,4]
#    mask1 = mag > 15 #removes saturated
#    mask2 = mag < 19 #removes very faint stars
#    mask = mask1 * mask2
#    tempsdata = sdata[mask]
#    
#    ### find sigma ### 
#    fwhm = np.nanmean(tempsdata['FWHM_WORLD_'+mon]*3600)
#    sigma = fwhm/np.sqrt(8*np.log(2))
#    ### Convolve Images ###
#    kernel = Gaussian2DKernel(sigma-0.01)
#    newim = convolve(newim, kernel, normalize_kernel=True)
    
    plt.subplot(1, 2, 1)
    ax = fig.gca()
    plt.imshow(newim, animated=True)
    ax.set_xticks([])
    ax.set_yticks([])
#    centre = [size,size]
#    pixelr = (1.5/3600) / const
#    aperture = CircularAperture(centre, pixelr)
    
#    aperture.plot()
    ax.set_title(sem)
    plt.subplot(1, 2, 2)
    t = np.linspace(1, 8, num=8)
    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
    x = [1,3,4,5,6,7,8]
    plt.xticks(t, years)
    plt.xlabel('Semester')
    plt.ylabel('K-band flux')
    plt.xlim(xmin=0, xmax=9)
#    plt.ylim(ymin=12500, ymax=22000)

    if sem == semesters[0]:
        plt.errorbar(x[0], flux[0], yerr=fluxerr[0], fmt = 'ro')
    elif sem == semesters[1]:
        plt.errorbar(x[0:2], flux[0:2], yerr=fluxerr[0:2], fmt = 'ro')
    elif sem == semesters[2]:
        plt.errorbar(x[0:3], flux[0:3], yerr=fluxerr[0:3], fmt = 'ro')
    elif sem == semesters[3]:
        plt.errorbar(x[0:4], flux[0:4], yerr=fluxerr[0:4], fmt = 'ro')
    elif sem == semesters[4]:
        plt.errorbar(x[0:5], flux[0:5], yerr=fluxerr[0:5], fmt = 'ro')
    elif sem == semesters[5]:
        plt.errorbar(x[0:6], flux[0:6], yerr=fluxerr[0:6], fmt = 'ro')
    else:
        plt.errorbar(x, flux, yerr=fluxerr, fmt = 'ro')
#    plt.errorbar(x, flux, yerr=fluxerr, fmt = 'ro', animated=True)

    plt.tight_layout()
#    return ax

#%%
ani = animation.FuncAnimation(fig, func, semesters, interval=1000)
Writer = animation.writers['ffmpeg']
writer = Writer(fps=1, metadata=dict(artist='Me'), bitrate=1800)

#ani.save('62253_K_withcurve.mp4',writer=writer)