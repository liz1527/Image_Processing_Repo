#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 15:59:07 2018

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
import matplotlib.gridspec as gridspec # to nest subplots
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots

### Read in fits files ###
imdata = fits.getdata('UDS-DR11-K.mef.fits') ### SHOULD CHANGE TO THE DR11 IMAGE
varydata = fits.getdata('variable_tables/no06_variables_chi30_DR11data_restframe.fits', 1)

### stick to those above 1e4 ###
flux = vari_funcs.flux5_stacks(varydata)
flux, varydata = vari_funcs.noneg(flux, varydata)
meanflux = np.nanmean(flux, axis=1)
varydata = varydata[meanflux >= 1e4]

#imthresh = 200 #fits.getheader('UDS_05B_K.fits')['SATURATE']
#imdata[imdata>imthresh] = imthresh

### Create mask array ###
maskx = np.zeros(imdata.shape[1]).astype(bool)
masky = np.zeros(imdata.shape[0]).astype(bool)

### Find coordinates of objects ###
x = varydata['X_IMAGE']
x = x.astype(int)
y = varydata['Y_IMAGE']
y = y.astype(int)
size = 25 # size of square

outer = gridspec.GridSpec(1,2, wspace=0.1) # define outer subplots
innerX = gridspec.GridSpecFromSubplotSpec(8, 8, subplot_spec=outer[0])#, 
                                         # wspace=0.1)
innernonX = gridspec.GridSpecFromSubplotSpec(8, 8, subplot_spec=outer[1])#, 
                                         # wspace=0.1)
fig = plt.figure(figsize=[19,10])
n=0
X=0

for ind in range(398):
    if n >= 64 and X>=64:
        continue
    print(ind)
    maskx[x[ind]-size:x[ind]+size] = 1
    masky[y[ind]-size:y[ind]+size] = 1
    
    newim = imdata[masky,:]
    newim = newim[:, maskx]
    
    #Threshold the image
#    maxim = np.max(newim)
#    minim = np.min(newim)
#    imthresh = 0.7*(maxim - minim)
    imthresh = newim[size,size]
    newim[newim>imthresh] = imthresh
    
#    ax1 = plt.subplot(9, 14, ind+1)
#    plt.plot(size, size, 'k+', markersize=10, mfc='none')
    
    if varydata['X-Ray'][ind] == False:
        if n >= 64:
            #reset mask
            maskx = np.zeros(imdata.shape[1]).astype(bool)
            masky = np.zeros(imdata.shape[0]).astype(bool)
            continue
        ax1 = plt.Subplot(fig, innernonX[n])
        ax1.imshow(newim)
        ax1.spines['bottom'].set_color('blue')
        ax1.spines['top'].set_color('blue')
        ax1.spines['left'].set_color('blue')
        ax1.spines['right'].set_color('blue')
        n += 1
    else:
        if X >= 64:
            #reset mask
            maskx = np.zeros(imdata.shape[1]).astype(bool)
            masky = np.zeros(imdata.shape[0]).astype(bool)
            continue
        ax1 = plt.Subplot(fig, innerX[X])
        ax1.imshow(newim)
        ax1.spines['bottom'].set_color('red')
        ax1.spines['top'].set_color('red')
        ax1.spines['left'].set_color('red')
        ax1.spines['right'].set_color('red')
        X += 1

    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['top'].set_linewidth(2)
    ax1.spines['left'].set_linewidth(2)
    ax1.spines['right'].set_linewidth(2)
    fig.add_subplot(ax1)

    plt.xticks([])
    plt.yticks([])
    
    #reset mask
    maskx = np.zeros(imdata.shape[1]).astype(bool)
    masky = np.zeros(imdata.shape[0]).astype(bool)

#fig.tight_layout()
plt.savefig('littleimages_poster.pdf')