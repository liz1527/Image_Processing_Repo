#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:11:49 2019

code to investigate stars going into PSFs

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
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
plt.close('all')

def get_star_data(sdata, psf_data, data, sem):
    colname = 'MAG_APER_'+sem
    
    ### Define coordinates ###
    refcoord = SkyCoord(psf_data['ALPHA_J2000_1']*u.degree, psf_data['DELTA_J2000_1']*u.degree)
    semcoord = SkyCoord(sdata['ALPHA_J2000_'+sem]*u.degree, sdata['DELTA_J2000_'+sem]*u.degree)
    
    ### Match catalogues and create new table ###
    idx, d2d , _ = match_coordinates_sky(refcoord, semcoord) #match these 'good' stars to create table
    tempsdata = sdata[idx]
    
#    # limit magnitude range used in PSF
#    mag = sdata[colname][:,4]
#    mask1 = mag > 15
#    mask2 = mag < 19
#    mask = mask1*mask2
#    
#    tempsdata = sdata[mask]
    
    x = tempsdata[colname][:,4]
    y = tempsdata[colname][:,4] - tempsdata[colname][:,1]
    
    allx = data[colname][:,4]
    ally = data[colname][:,4] - data[colname][:,1]
    
    return x, y, allx, ally

def get_star_data_DR11(sdata, data):
    colname = 'MAG_APER'
    
    # limit magnitude range used in PSF
    mag = sdata[colname][:,4]
    mask1 = mag > 15
    mask2 = mag < 19
    mask3 = mag-sdata[colname][:,1] > -0.7 #pointlike-ness criterion
    mask4 = mag-sdata[colname][:,1] < -0.5 #pointlike-ness criterion
    mask = mask1*mask2*mask3*mask4
    
    tempsdata = sdata[mask]
    
    x = tempsdata[colname][:,4]
    y = tempsdata[colname][:,4] - tempsdata[colname][:,1]
    
    allx = data[colname][:,4]
    ally = data[colname][:,4] - data[colname][:,1]
    
    return x, y, allx, ally, tempsdata
    
Jdata = fits.open('mag_flux_tables/J/mag_flux_table_J_cleaned.fits')[1].data
Jsdata = fits.open('mag_flux_tables/J/stars_mag_flux_table_J_cleaned.fits')[1].data
Hdata = fits.open('mag_flux_tables/H/mag_flux_table_H_cleaned.fits')[1].data
Hsdata = fits.open('mag_flux_tables/H/stars_mag_flux_table_H_cleaned.fits')[1].data
Kdata = fits.open('mag_flux_tables/K/mag_flux_table_cleaned.fits')[1].data
Ksdata = fits.open('mag_flux_tables/K/stars_mag_flux_table_cleaned.fits')[1].data
sdata = fits.open('UDS_catalogues/DR11_output_stars.fits')[1].data
data = fits.open('UDS_catalogues/DR11_output.fits')[1].data
semesters = ['06B', '07B', '08B', '09B', '10B', '11B', '12B']


x, y, allx, ally, psf_star_data = get_star_data_DR11(sdata, data)


plt.figure()
plt.scatter(allx, ally, c='tab:grey',marker='+')
plt.scatter(x,y,c='b')
plt.xlabel('MAG_APER[4]')
plt.ylabel('MAG_APER[4] - MAG_APER[1]')
plt.xlim(xmax=21, xmin=10)
plt.ylim(ymax=-0.4, ymin=-1)
plt.title('DR11 K')
plt.tight_layout()

for sem in semesters: 
    
#    Jx, Jy, allJx, allJy = get_star_data(Jsdata, Jdata, sem)
#    
#    plt.figure()
#    plt.scatter(allJx, allJy, c='tab:grey',marker='+')
#    plt.scatter(Jx,Jy,c='b')
#    plt.xlabel('MAG_APER[4]')
#    plt.ylabel('MAG_APER[4] - MAG_APER[1]')
#    plt.xlim(xmax=21, xmin=10)
#    plt.ylim(ymax=0, ymin=-2.5)
#    plt.title(sem+' J')
    
    Hx, Hy, allHx, allHy = get_star_data(Hsdata, psf_star_data, Hdata, sem)
    
    plt.figure()
    plt.scatter(allHx, allHy, c='tab:grey',marker='+')
    plt.scatter(Hx,Hy,c='b')
    plt.xlabel('MAG_APER[4]')
    plt.ylabel('MAG_APER[4] - MAG_APER[1]')
    plt.xlim(xmax=21, xmin=10)
    plt.ylim(ymax=0, ymin=-2.5)
    plt.title(sem+' H')
    
#    Kx, Ky, allKx, allKy = get_star_data(Ksdata, Kdata, sem)
#    
#    plt.figure()
#    plt.scatter(allKx, allKy, c='tab:grey',marker='+')
#    plt.scatter(Kx,Ky,c='b')
#    plt.xlabel('MAG_APER[4]')
#    plt.ylabel('MAG_APER[4] - MAG_APER[1]')
#    plt.xlim(xmax=21, xmin=10)
#    plt.ylim(ymax=0, ymin=-2.5)
#    plt.title(sem+' K')