#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 16:44:00 2018

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

### Open the fits files and get data ###
varys = fits.open('variable_tables/no06_variables_chi30_2arcsec.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

#varys = vari_funcs.chandra_only(varys)

flux = vari_funcs.flux4_stacks(varys)
flux, varys = vari_funcs.noneg(flux, varys)
flux, fluxerr, newvarys = vari_funcs.create_quad_error_array(sigtb, varys, aper=4)
#flux,fluxerr = vari_funcs.normalise_flux_and_errors(flux, fluxerr)

#set up time variable for plot
t = np.linspace(1, 8, num=8)
years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
xcurve = [1,3,4,5,6,7,8]

### reduce to just single lightcurve ###
obnum = 62243
mask = newvarys['NUMBER_05B'] == obnum
flux = flux[mask].reshape(len(xcurve))
fluxerr = fluxerr[mask].reshape(len(xcurve))
obdata = newvarys[mask]

#### Subtract non flare semester to get flare ###
#flux = flux - flux[-1]

### Plot ###
plt.figure(figsize=[8,6])
if newvarys['X-ray'][mask] == True:
    plt.errorbar(xcurve, flux, yerr=fluxerr, fmt='o', color='r')
else:
    plt.errorbar(xcurve, flux, yerr=fluxerr, fmt='o', color='b')
plt.xlabel('Semester')
plt.ylabel('K-band flux')
plt.title('Lightcurve of Object '+str(obnum))
plt.xticks(t, years)
plt.tight_layout()
#plt.savefig('Chi40Lightcurves/cleaned/no06/mag_'+str(n))#+str(varys['NUMBER_05B'][n])+'_lightcurve.png')
#plt.close('all')

#### Get little image ### 
#semesters = ['05B', '07B', '08B', '09B', '10B', '11B', '12B']#'06B', 
#artists = []
#for n, sem in enumerate(semesters):
#    print(sem)
#    if sem == '10B':
#        imdata = fits.getdata('cleaned_UDS_'+sem+'_K.fits')
#    else:
#        imdata = fits.getdata('extra_clean_no06_UDS_'+sem+'_K.fits')
#    
#    ### Find coordinates of objects ###
#    x = obdata['X_IMAGE_'+sem]
#    x = x.astype(int)#int(x)
#    y = obdata['Y_IMAGE_'+sem]
#    y = y.astype(int)#int(x) 
#    
#    size = 75 # size of half side of square
#    newim = imdata[y[0]-size:y[0]+size,x[0]-size:x[0]+size]
#    
#    del imdata
#    print(newim[size,size])
#    
#    imuppthresh = 200
#    newim[newim>imuppthresh] = imuppthresh
#    imlowthresh = 0 
#    newim[newim<imlowthresh] = imlowthresh
#    
#    ### code from stack exchange
#    ax = plt.gca()
#    im = OffsetImage(newim, zoom=0.75)
#    ab = AnnotationBbox(im, (xcurve[n], flux[n]), xycoords='data', frameon=False)
#    artists.append(ax.add_artist(ab))
#ax.update_datalim(np.column_stack([xcurve, flux]))
#ax.autoscale()
###    plt.subplot(4, 2, n)
##    snap = plt.imshow(newim, extent=(xcurve[n], xcurve[n]+0.2, flux[n], flux[n]+0.2))
###    plt.plot(size, size, 'k+', markersize=10, mfc='none')
##    plt.xticks([])
##    plt.yticks([])
##    plt.title(obnum)
#    
##    snaps.append([snap])
#
