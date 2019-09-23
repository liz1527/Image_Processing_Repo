#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 14:35:24 2019

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from photutils import CircularAperture, aperture_photometry
from photutils import centroid_com, centroid_1dg, centroid_2dg
plt.close('all') #close any open plots

#varys = fits.open('variable_tables/no06_variables_chi40.fits')[1].data
#imdata = fits.getdata('extra_clean_no06_UDS_08B_K.fits')
varys = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
imdata = fits.getdata('extra_clean_no06_UDS_08B_K.fits')
imdata12 = fits.getdata('extra_clean_no06_UDS_12B_K.fits')
hdr08B = fits.getheader('extra_clean_no06_UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM
# 
obnum = 62243#57824
obvarys = varys[varys['NUMBER_05B'] == obnum]
#  
### Find coordinates of objects ###
x = obvarys['X_IMAGE_08B']
x = x.astype(int)#int(x)
y = obvarys['Y_IMAGE_08B']
y = y.astype(int)#int(x) 
x12 = obvarys['X_IMAGE_12B']
x12 = x12.astype(int)#int(x)
y12 = obvarys['Y_IMAGE_12B']
y12 = y12.astype(int)#int(x) 

#### Subtract images before cutting down ###
#newimdata = imdata - imdata12


#
size = 15 # size of half side of square
im = imdata[y[0]-size:y[0]+size,x[0]-size:x[0]+size]
im12 = imdata12[y12[0]-size:y12[0]+size,x12[0]-size:x12[0]+size]
#newim = newimdata[y[0]-size:y[0]+size,x[0]-size:x[0]+size]

#### Find coordinates of object ###
#x = 1606
#y = 5722

#size = 105 # size of half side of square
#newim = imdata[y-size:y+size,x-size:x+size]

del imdata, imdata12
#print(newim[size,size])

newim = im - im12

plt.figure(figsize=[18,8])
plt.subplot(131)
plt.imshow(newim)
vari_funcs.no_ticks()
#plt.plot(size,size,'k+')
plt.title('08B - 12B')
#x1, y1 = centroid_com(newim)
#x2, y2 = centroid_1dg(newim)
x3, y3 = centroid_2dg(newim)
#plt.plot(x1,y1,'r+')
#plt.plot(x2,y2,'m+')
plt.plot(x3,y3,'b+')

imuppthresh = 190
im[im>imuppthresh] = imuppthresh
imlowthresh = 0 
im[im<imlowthresh] = imlowthresh
im12[im12>imuppthresh] = imuppthresh
im12[im12<imlowthresh] = imlowthresh
#newim[newim>imuppthresh] = imuppthresh
#newim[newim<imlowthresh] = imlowthresh

#plt.figure(figsize=[8,8])
plt.subplot(132)
plt.imshow(im)
vari_funcs.no_ticks()
#plt.plot(size,size,'k+')
plt.title('08B')
##x1, y1 = centroid_com(im)
##x2, y2 = centroid_1dg(im)
#x3, y3 = centroid_2dg(im)
##plt.plot(x1,y1,'r+')
##plt.plot(x2,y2,'m+')
#plt.plot(x3,y3,'b+')


#plt.figure(figsize=[8,8])
plt.subplot(133)
plt.imshow(im12)
vari_funcs.no_ticks()
#plt.plot(size,size,'k+')
plt.tight_layout()
plt.title('12B')
#x1, y1 = centroid_com(im12)
#x2, y2 = centroid_1dg(im12)
x3_12, y3_12 = centroid_2dg(im12)
#plt.plot(x1,y1,'r+')
#plt.plot(x2,y2,'m+')
#plt.plot(x3_12,y3_12,'b+')

### find dist between SN and galaxy ###
xdiff = x3_12-x3
ydiff = y3_12-y3
dist_sq = (xdiff)**2 + (ydiff)**2
dist = np.sqrt(dist_sq)
print('offset of flare is '+str(dist*0.1342)+' arcsecs')

err_meas = 0.5
err_xdiff = xdiff * np.sqrt( (err_meas/x3)**2 + (err_meas/x3_12)**2 )
err_ydiff = ydiff * np.sqrt( (err_meas/y3)**2 + (err_meas/y3_12)**2 )
err_r_sq = dist_sq * np.sqrt( (err_xdiff/xdiff)**2 + (err_ydiff/ydiff)**2 )
err_r = dist * (err_r_sq/dist_sq)
err_dist_arcsec = err_r * 0.1342

plt.errorbar(x3,y3,yerr=err_meas, xerr=err_meas, fmt='bx')
plt.errorbar(x3_12,y3_12,yerr=err_meas, xerr=err_meas, fmt='bx')

centre = [size,size]
pixelr = (0.7/3600) / const
aperture = CircularAperture(centre, pixelr)

### Determine flux within 3 arcsec apertures ###
phot = aperture_photometry(newim, aperture)
aperflux = phot['aperture_sum'][0]
#aperture.plot()

### Determine flux off side flare ###
centre = [100,105.5]
pixelr = (0.7/3600) / const
flareaperture = CircularAperture(centre, pixelr)
flarephot = aperture_photometry(newim, flareaperture)
flareflux = flarephot['aperture_sum'][0]
#flareaperture.plot()

mag = 30 - 2.5*np.log10(aperflux)
flaremag = 30 - 2.5*np.log10(flareflux)
