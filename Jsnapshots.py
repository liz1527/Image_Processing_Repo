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
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
#import math
#from astropy.stats import median_absolute_deviation
#import vari_funcs_no06 #my module to help run code neatly
plt.close('all') #close any open plots

obnum=62243
### Read in fits files ###
obdata = fits.getdata('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits', 1)
obdata = obdata[obdata['NUMBER_05B']==obnum] #restrict to just 252446

#flux = np.array([obdata['FLUX_APER_05B'][:,4],#obdata['FLUX_APER_06B'][:,4],
#                obdata['FLUX_APER_07B'][:,4],obdata['FLUX_APER_08B'][:,4],
#                obdata['FLUX_APER_09B'][:,4],obdata['FLUX_APER_10B'][:,4], 
#                obdata['FLUX_APER_11B'][:,4],obdata['FLUX_APER_12B'][:,4]])
semesters = ['sep+oct_1', 'oct_2','nov']#'06B', 
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
    imdata = fits.getdata('J/UDS_'+sem+'_08_J.fits')
    ### Find coordinates of objects ###
    x = obdata['X_IMAGE_08B']
    x = x.astype(int)#int(x)
    y = obdata['Y_IMAGE_08B']
    y = y.astype(int)#int(x) 
    
    size = 100 # size of square
    newim = imdata[y[0]-size:y[0]+size,x[0]-size:x[0]+size]
    
    del imdata
    print(newim[size,size])
    
    imuppthresh = 25
    newim[newim>imuppthresh] = imuppthresh
    imlowthresh = 0 
    newim[newim<imlowthresh] = imlowthresh
    
    ### Convolve Images ###
    kernel = Gaussian2DKernel(0.5)
    newim = convolve(newim, kernel, normalize_kernel=True)
    
#    plt.subplot(4, 2, n)
    snap = plt.imshow(newim, animated=True)
#    plt.plot(size, size, 'k+', markersize=10, mfc='none')
    plt.xticks([])
    plt.yticks([])
#    plt.title(sem)
    
    snaps.append([snap])
    
    n += 1

#%%
ani = animation.ArtistAnimation(fig, snaps, interval=500, blit=True)
Writer = animation.writers['ffmpeg']
writer = Writer(fps=1, metadata=dict(artist='Me'), bitrate=1800)
##
#ani.save('62253_sep08-nov08_J_smoothed.mp4',writer=writer)