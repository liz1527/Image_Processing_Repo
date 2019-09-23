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
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
#import vari_funcs_no06 #my module to help run code neatly
plt.close('all') #close any open plots

obnum=57824
### Read in fits files ###
#obdata = fits.getdata('mag_flux_tables/mag_flux_table_best.fits', 1)
obdata = fits.getdata('mag_flux_tables/mag_flux_table_best_J.fits', 1)
obdata = obdata[obdata['NUMBER_05B']==obnum] #restrict to just 252446

#flux = np.array([obdata['FLUX_APER_05B'][:,4],#obdata['FLUX_APER_06B'][:,4],
#                obdata['FLUX_APER_07B'][:,4],obdata['FLUX_APER_08B'][:,4],
#                obdata['FLUX_APER_09B'][:,4],obdata['FLUX_APER_10B'][:,4], 
#                obdata['FLUX_APER_11B'][:,4],obdata['FLUX_APER_12B'][:,4]])

semesters = ['05B', '06B', '07B', '08B']#, '09B', '10B', '11B', '12B']

n=1
fig = plt.figure()
snaps = []
for sem in semesters:
#    print(sem)
#    if sem == '10B':
#        imdata = fits.getdata('cleaned_UDS_'+sem+'_K.fits')
#    else:
#        imdata = fits.getdata('extra_clean_no06_UDS_'+sem+'_K.fits')
#    
    imdata = fits.getdata('J/UDS_'+sem+'_J.fits')
    ### Find coordinates of objects ###
    x = obdata['X_IMAGE_'+sem]
    x = x.astype(int)#int(x)
    y = obdata['Y_IMAGE_'+sem]
    y = y.astype(int)#int(x) 
    
    size = 50 # size of square
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
    ax = fig.gca()
    imdata = fits.getdata('J/UDS_'+sem+'_J.fits')
    ### Find coordinates of objects ###
    x = obdata['X_IMAGE_'+sem]
    x = x.astype(int)#int(x)
    y = obdata['Y_IMAGE_'+sem]
    y = y.astype(int)#int(x) 
    
    size = 50 # size of square
    newim = imdata[y[0]-size:y[0]+size,x[0]-size:x[0]+size]
    
    del imdata
    
    imuppthresh = 50
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
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(sem)
#    plt.subplot(4, 10, n)
    ax.imshow(newim, animated=True)
#    return ax

#%%
ani = animation.FuncAnimation(fig, func, semesters, interval=500)
Writer = animation.writers['ffmpeg']
writer = Writer(fps=1, metadata=dict(artist='Me'), bitrate=1800)

ani.save('62253_05B-08B_J_small.mp4',writer=writer)