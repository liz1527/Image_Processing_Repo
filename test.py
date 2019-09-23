#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 15:45:02 2017

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
#import vari_funcs_no06 #my module to help run code neatly
plt.close('all') #close any open plots

#imdata = fits.open('UDS_12B_K_bin2x2.fits')[0].data

sem05B = fits.open('SE_outputs_yearstacks/05B_output.fits')[1].data
sem07B = fits.open('SE_outputs_yearstacks/07B_output.fits')[1].data
sem08B = fits.open('SE_outputs_yearstacks/08B_output.fits')[1].data
sem09B = fits.open('SE_outputs_yearstacks/09B_output.fits')[1].data
sem10B = fits.open('SE_outputs_yearstacks/10B_output.fits')[1].data
sem11B = fits.open('SE_outputs_yearstacks/11B_output.fits')[1].data
sem12B = fits.open('SE_outputs_yearstacks/12B_output.fits')[1].data

sem05Bnew = fits.open('SE_outputs_yearstacks/newFR_05B_output.fits')[1].data
sem07Bnew = fits.open('SE_outputs_yearstacks/07B_output.fits')[1].data
sem08Bnew = fits.open('SE_outputs_yearstacks/newFR_08B_output.fits')[1].data
sem09Bnew = fits.open('SE_outputs_yearstacks/newFR_09B_output.fits')[1].data
sem10Bnew = fits.open('SE_outputs_yearstacks/newFR_10B_output.fits')[1].data
sem11Bnew = fits.open('SE_outputs_yearstacks/newFR_11B_output.fits')[1].data
sem12Bnew = fits.open('SE_outputs_yearstacks/newFR_12B_output.fits')[1].data
#sem10B = fits.open('starsfwhm.fits')[1].data

#hdr08B = fits.getheader('UDS_08B_K.fits') # random year (same in all)
#const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

colname = 'MAG_APER'
#
##Extract the flux radii and remove negative values
#fluxrad10 = sem10B[colname][:,1]
#fluxrad10 = fluxrad10[fluxrad10>0]
#fluxrad05 = sem05B[colname][:,1]
#fluxrad05 = fluxrad05[fluxrad05>0]
#fluxrad05new = sem05Bnew[colname][:,1]
#fluxrad05new = fluxrad05new[fluxrad05new>0]
#fluxrad08 = sem07B[colname][:,1]
#fluxrad08 = fluxrad08[fluxrad08>0]
#fluxrad08new = sem08Bnew[colname][:,1]
#fluxrad08new = fluxrad08new[fluxrad08new>0]
#
###Put data in array
avgFRnew = np.array([np.median(sem05Bnew[colname][:,4]),
                    np.median(sem07Bnew[colname][:,4]),
                    np.median(sem08Bnew[colname][:,4]), 
                    np.median(sem09Bnew[colname][:,4]), 
                    np.median(sem10Bnew[colname][:,4]), 
                    np.median(sem11Bnew[colname][:,4]), 
                    np.median(sem12Bnew[colname][:,4])])
avgFR = np.array([np.median(sem05B[colname][:,4]),
                    np.median(sem07B[colname][:,4]),
                    np.median(sem08B[colname][:,4]), 
                    np.median(sem09B[colname][:,4]), 
                    np.median(sem10B[colname][:,4]), 
                    np.median(sem11B[colname][:,4]), 
                    np.median(sem12B[colname][:,4])])

#avgFWHM = np.array([np.median(sem05Bnew['FWHM_WORLD']),
#                    np.median(sem07Bnew['FWHM_WORLD']), 
#                    np.median(sem08Bnew['FWHM_WORLD']), 
#                    np.median(sem09Bnew['FWHM_WORLD']),
#                    np.median(sem10B['FWHM_10B']), 
#                    np.median(sem11Bnew['FWHM_WORLD']),
#                    np.median(sem12Bnew['FWHM_WORLD'])])
t = np.array([1,3,4,5,6,7,8])
years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']

#Plot graph in new figure
plt.figure()
ax = plt.axes()
plt.xticks(np.linspace(1,8,8), years)
ax.plot(t, avgFRnew, 'ro')
plt.xlabel('Semester')
plt.ylabel('0.5 Flux Radius (pixels)')
plt.title('Average FWHM curve for Convolved')
#plt.ylim(2.9, 3.4)

plt.figure()
ax = plt.axes()
plt.xticks(np.linspace(1,8,8), years)
ax.plot(t, avgFR, 'ro')
plt.xlabel('Semester')
plt.ylabel('0.5 Flux Radius (pixels)')
plt.title('Average FWHM curve for unconvolved')
#plt.ylim(2.9, 3.4)
#avgFWHM /= const            
#old05B = sem05B[colname]
#new05B = sem05Bnew[colname]
#diff = new05B - old05B
#
#print(np.size(diff[diff>0]))
#print(np.size(diff[diff==0]))
#print(np.size(diff[diff<0]))



#old = fits.open('UDS_05B_K_bin2x2.fits')[0].data
#new = fits.open('new_UDS_05B_K_bin2x2.fits')[0].data
#print(np.max(old))
#print(np.max(new))
#
#### Read in fits files ###
#imdata = fits.getdata('From_Captain/new_UDS_09B_K.fits')
#imdata2 = fits.getdata('new_UDS_09B_K.fits')
#print(np.sum(sem05Bnew[colname][:,1]-sem05B[colname][:,1]))

#newsem = sem05Bnew[(sem05Bnew[colname]-sem05B[colname])!=0]