#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 12:24:26 2018

@author: ppxee
"""
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.stats import norm
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from photutils import CircularAperture, aperture_photometry
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
#import vari_funcs
plt.close('all')

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

def FWHM2sigma(FWHM, const):
    ''' Function to convert the FWHM of a distribution into a sigma for that
    distribution. It assumes the distribution is gaussian.
    Input:
        FWHM = Full width half maximum of a distriubtution (in my case usually
                of an object from SExtractor)
    Output:
        sigma = standard deviation value of a guassian distribution with the 
                given FWHM. This roughly equates to the psf of the object. '''
    FWHM /= const
    return FWHM/np.sqrt(8*np.log(2))

def fluxrad2sigma(fluxrad):
    return fluxrad/np.sqrt(8*np.log(2))

def convolve_psf(mon, psf, newpsf, sigmakernel):
    ## Open image ###
    im05B = psf[mon]
    
    ## Convolve Image ###
    print('Convolving ', mon)
    kernel = Gaussian2DKernel(sigmakernel)
    newpsf[mon] = convolve(im05B, kernel, normalize_kernel=True) 

def convolve_one_psf(psf, sigmakernel):
#    print(np.sum(psf))
    ## Convolve Image ###
#    print('Convolving ', mon)
    kernel = Gaussian2DKernel(sigmakernel)
    newpsf = convolve(psf, kernel, normalize_kernel=True) 
    newpsf /= np.nansum(newpsf)
#    print(np.sum(newpsf))
    return newpsf

def quadrants(initdata,sem):
    
    ira = initdata['X_IMAGE_'+sem]
    idec = initdata['Y_IMAGE_'+sem]

    ### define bounds of quadrants ###
    midra = 12450
    middec = 13310
    
    ### create masks for quadrant ###
    mask1 = ira < midra
    mask2 = idec >= middec
    quad1data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec >= middec
    quad2data = initdata[mask1*mask2]
    
    mask1 = ira < midra
    mask2 = idec < middec
    quad3data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec < middec
    quad4data = initdata[mask1*mask2]    
    
    return quad1data, quad2data, quad3data, quad4data

#sdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
psf_data = fits.open('UDS_catalogues/DR11_stars_for_PSFs.fits')[1].data
oldsdata = fits.open('mag_flux_tables/K/month/month_stars_mag_flux_table_K_cleaned.fits')[1].data
hdr08B = fits.getheader('Images/UDS-DR11-K.mef.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

### set up month tick details ###
month_info = fits.open('monthly_numbers.fits')[1].data #get month count data
full_months = month_info['Month'] #extract month nanes
tick_inds = np.load('tick_inds_K.npy') #load tick locations
mask = np.zeros(len(full_months)) #set up mask
mask[tick_inds] = 1
mask = mask.astype(bool)
month_ticks = np.copy(full_months)
month_ticks = month_ticks[mask]#retrieve tick details

#### set up month x array ###
#x_month = np.copy(full_months)
#mask = month_info['Frames in V11'] != 0
#x_month = x_month[mask]

colname = 'FWHM_WORLD_'
#data = mon05B[colname][:,1]
#semesters = ['05B', '07B', '08B', '09B', '10B', '11B', '12B']#['05B','10B']
months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07',  
          'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09',  
          'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 
          'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 
          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
centre = [29,29]
quads = np.arange(1,5,1)

x = np.arange(1, len(month_info['Frames in v11'])+1)
mask = np.isin(full_months, months)
x_months = x[mask]

avgFWHM = np.zeros(len(months))
oldavgFWHM = np.zeros(len(months))
avgflux = np.zeros(len(months))
extras = np.zeros([len(months),4]) + 99# 4 columns, 1 per quad, 99 is null
#oldavgflux = np.zeros(len(months))
psf = {}
oldpsf = {}
for i, mon in enumerate(months):
    # for new
    colnames = colname+mon
    
    ### Define coordinates ###
    refcoord = SkyCoord(psf_data['ALPHA_J2000_1']*u.degree, 
                        psf_data['DELTA_J2000_1']*u.degree)
    moncoord = SkyCoord(oldsdata['ALPHA_J2000_'+mon]*u.degree, 
                        oldsdata['DELTA_J2000_'+mon]*u.degree)
    
    ### Match catalogues and create new table ###
    idx, d2d , _ = match_coordinates_sky(refcoord, moncoord) 
    tempsdata = oldsdata[idx]
    
    ### get quaddata
    quaddata = quadrants(tempsdata,mon)
    
    for n, quad in enumerate(quads):
        # for new
        colnames = colname+mon
        # for old
        oldmag = quaddata[n]['MAG_APER_'+mon][:,4]
        oldavgFWHM[n] = np.median(quaddata[n][colnames]) #* 3600
        oldpsf[quad] = fits.open('PSFs/K/month/Quad_PSFs/cleaned_'+
              mon+'_'+str(quad)+'_K_PSF.fits')[0].data
   
    ### Find maximum FWHM as this is what all the others willl become ###
    aimind = np.argmax(oldavgFWHM)
    aimquad = quads[aimind]
    aimpsf = oldpsf[aimquad]
    
    ### Convert FWHM into a sigma ###
    sigmaold = np.array([FWHM2sigma(fwhm, const) for fwhm in oldavgFWHM])
    sigmabroad = sigmaold[aimind]
    
    phot = {}
    flux = np.zeros(len(months))
    oldphot = {}
    oldflux = np.zeros(len(months))
    
    ### testing the extra factor method ###
    tests =np.linspace(-1,0.5,1000)
    r = np.arange(0,42,1) * const * 3600 # define radius values
    
    #flux10B = oldflux[3]
    aimrp = radial_profile(oldpsf[aimquad], centre)
    sqrtaimrp = np.sqrt(aimrp)
    
    
    newpsf = {}
    newphot = {}
    newflux = np.zeros(len(months))
    aperflux = np.empty(len(months))
    oldaperflux = np.empty(len(months))
    pixelr = (1.5/3600) / const
    aperture = CircularAperture(centre, pixelr)
    t = np.linspace(1, 8, num=8)
    
    plt.figure(1, figsize=[9,6])
    plt.plot(r, sqrtaimrp,label=aimquad)   
    
    for n, quad in enumerate(quads):
        if quad == aimquad:
#            phot = aperture_photometry(aimpsf, aperture)
#            aperflux[n] = phot['aperture_sum'][0]
#            oldaperflux[n] = phot['aperture_sum'][0]
    #        plt.figure()
    #        plt.imshow(aimpsf)
    #        vari_funcs.no_ticks()
            continue
        sigma = sigmaold[n]
        singlephot = {}
        singleflux = np.zeros(len(tests))
        sumdiffold = 10
        for m, extra in enumerate(tests):
    #        print(extra)
            sigmakernel = np.sqrt(sigmabroad**2 - sigma**2) + extra
            if sigmakernel <= 0:
                continue
            new = convolve_one_psf(oldpsf[quad], sigmakernel)    
            
            ### Get radial profile ###
            radialprofile = radial_profile(new, centre)
            sqrtrp = np.sqrt(radialprofile)
            
            diff = aimrp[:12] - radialprofile[:12]
            sumdiff = np.nansum(diff)
    #        plt.figure(1)
    ##            plt.subplot(4,2,n+1)
    ##            plt.plot(r, sqrtaimrp,label='10B')    
    #        plt.plot(r,sqrtrp, '--', label=mon+' '+str(extra))
    #        plt.ylabel('sqrt(Flux)')
    #        plt.xlabel('Radius (arcsec)')
    #        plt.legend()
    #        print(sumdiff)
            if sumdiff > 0:
                if sumdiffold < sumdiff:
                    print('extra ='+str(tests[m-1]))
                    extras[i,n] = tests[m-1]
    #                plt.plot(t[n], singleflux[m-1],'o', label=extras[n])
                else:
                    print('extra ='+str(extra))
                    extras[i,n] = extra    
                
#                plt.figure(1, figsize=[9,6])
#    #            plt.subplot(4,2,n+1)
#    #            plt.plot(r, sqrtaimrp,label='10B')    
#                plt.plot(r,sqrtrp, '--', label=str(n)+' '+str(extra))
#                plt.ylabel('sqrt(Flux)')
#                plt.xlabel('Radius (arcsec)')
#                plt.title(mon)
#                plt.xlim(xmax=2.5)
#    #            plt.legend()    
#                plt.tight_layout()
#    
#                plt.figure(2, figsize=[9,6])
#                plt.plot(r[:12], diff, label=str(n))
#                plt.ylim(ymax=0.005, ymin=-0.002)
#                plt.hlines(0,0,1.5)
#                plt.xlim(xmin=0, xmax=1.5)
#    #            plt.legend()
#                plt.ylabel('Difference from aim rp')
#                plt.xlabel('Radius (arcsec)')
#                plt.title(mon)
#                plt.tight_layout()
    
                phot = aperture_photometry(new, aperture)
                aperflux[n] = phot['aperture_sum'][0]
                oldphot = aperture_photometry(oldpsf[quad], aperture)
                oldaperflux[n] = oldphot['aperture_sum'][0]
                
    #            plt.figure()
    #            plt.subplot(121)
    #            plt.imshow(np.log(oldpsf[mon]))
    #            vari_funcs.no_ticks()
    #            
    #            plt.subplot(122)
    #            plt.imshow(np.log(new))
    #            vari_funcs.no_ticks()
                break
            else:
                sumdiffold = sumdiff       
#    plt.figure(1)
#    plt.savefig('month_stacks/quad_rp_matching/sqrtrp/'+mon+'_sqrtrp.png')
#    plt.close(1)
#    
#    plt.figure(2)
#    plt.savefig('month_stacks/quad_rp_matching/rpdiff/'+mon+'_rpdiff.png')
#    plt.close(2)
#
##x = [1,3,4,5,6,7,8]
###years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
#plt.figure(figsize=[9,6])
#plt.plot(x_months, aperflux,'o-', label='new')
#plt.plot(x_months, oldaperflux, 'o-', label='old')
#plt.xticks(tick_inds, month_ticks, rotation = 'vertical')
#plt.xlabel('Month')
#plt.legend()
#plt.tight_layout()
#
#plt.figure(figsize=[9,6])
#plt.plot(x_months, oldavgFWHM, 'o-', label='old')
#plt.xticks(tick_inds, month_ticks, rotation = 'vertical')
#plt.xlabel('Month')
#plt.legend()
#plt.tight_layout()
#
#np.save('extrascleanedK_quad_month', extras)