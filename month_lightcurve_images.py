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

def month_avg_lightcurve(avgflux, avgfluxerr):
    months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07',  
          'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09',  
          'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 
          'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 
          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
   
       
    #set up time variable for plot
    nums = fits.open('monthly_numbers.fits')[1].data
    t = np.linspace(1, len(nums), num=len(nums))
    tdataind = np.isin(nums['Month'], months)
    tdata = t[tdataind]
    
    ticks = nums['Month']
    mask = np.zeros(len(t))
    inds = [0,4,14,16,23,26,36,38,46,53,59,65,71,77,82,86]
    mask[inds] = 1
    mask = mask.astype(bool)
    ticks[~mask] = ''
    
    #Plot graph in new figure
    plt.figure(figsize=[19,7])
    plt.xticks(t, ticks, rotation='vertical')
    plt.errorbar(tdata, avgflux, yerr=avgfluxerr, fmt = 'bo')
    plt.xlabel('Month')
    plt.ylabel('K-band flux of object')
    plt.title('Average Lightcurve')
    plt.tight_layout()
    return



tbdata = fits.open('mag_flux_tables/month_mag_flux_table_best.fits')[1].data
varys = fits.open('variable_tables/no06_variables_chi40.fits')[1].data
obnum = 62242


#mask = fitsdata['NUMBER_1'] == ob
obdata = tbdata[tbdata['NUMBER'] == obnum]
obvarys = varys[varys['NUMBER_05B'] == 62243]

months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07',  
          'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09',  
          'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 
          'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 
          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']

for month in months:
#    if month == 'sep05':
#        flux = obdata['MAG_APER_'+month][:,4]
#        fluxerr = obdata['MAGERR_APER_'+month][:,4]
#    else:
#        flux = np.append(flux, obdata['MAG_APER_'+month][:,4])
#        fluxerr = np.append(fluxerr, obdata['MAGERR_APER_'+month][:,4])
        
    if month == 'sep05':
        flux = obdata['FLUX_APER_'+month][:,4]
        fluxerr = obdata['FLUXERR_APER_'+month][:,4]
    else:
        flux = np.append(flux, obdata['FLUX_APER_'+month][:,4])
        fluxerr = np.append(fluxerr, obdata['FLUXERR_APER_'+month][:,4])

#mask = flux == 99
mask = flux <= 0
flux[mask] = np.nan
fluxerr[mask] = np.nan
month_avg_lightcurve(flux, fluxerr)
plt.title('Light curve for DR11 ID 62253')#%i' % obnum)

### Get little image ### 
semesters = ['05B', '08B', '12B']#['05B', '07B', '08B', '09B', '10B', '11B', '12B']#'06B', 
xcoords = [15, 50, 80]
ycoords = [20000,20000,20000]
artists = []
for n, sem in enumerate(semesters):
    print(sem)
    if sem == '10B':
        imdata = fits.getdata('cleaned_UDS_'+sem+'_K.fits')
    else:
        imdata = fits.getdata('extra_clean_no06_UDS_'+sem+'_K.fits')
    
    ### Find coordinates of objects ###
    x = obvarys['X_IMAGE_'+sem]
    x = x.astype(int)#int(x)
    y = obvarys['Y_IMAGE_'+sem]
    y = y.astype(int)#int(x) 
    
    size = 30 # size of half side of square
    newim = imdata[y[0]-size:y[0]+size,x[0]-size:x[0]+size]
    
    del imdata
#    print(newim[size,size])
    
    imuppthresh = 200
    newim[newim>imuppthresh] = imuppthresh
    imlowthresh = 0 
    newim[newim<imlowthresh] = imlowthresh
    
    ### code from stack exchange
    ax = plt.gca()
    im = OffsetImage(newim, zoom=2.5)
    ab = AnnotationBbox(im, (xcoords[n], ycoords[n]), xycoords='data', frameon=False)
    artists.append(ax.add_artist(ab))
    
ax.update_datalim(np.column_stack([xcoords, ycoords]))
ax.autoscale()
### put on arrows ###
plt.arrow(12,19000,-7,-3000, head_width=1, 
          head_length=500, color='k') # for first snapshot
plt.arrow(45.5,21000,-5,0, head_width=300, 
          head_length=1, color='k') # for first snapshot
plt.arrow(80,19000,+3,-3500, head_width=1, 
          head_length=500, color='k') # for first snapshot

